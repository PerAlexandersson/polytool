//! Stdio tool discovery, schema budgets and handle-based proof workflow.

use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, Command, Stdio};

fn start_server(
    directory: &std::path::Path,
) -> (Child, ChildStdin, BufReader<std::process::ChildStdout>) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool-mcp"))
        .arg("--report-dir")
        .arg(directory)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::null())
        .spawn()
        .expect("spawn polytool-mcp");
    let stdin = child.stdin.take().expect("child stdin");
    let stdout = BufReader::new(child.stdout.take().expect("child stdout"));
    (child, stdin, stdout)
}

fn send(stdin: &mut ChildStdin, value: Value) {
    writeln!(stdin, "{}", serde_json::to_string(&value).unwrap()).expect("write request");
    stdin.flush().expect("flush request");
}

fn read_response(stdout: &mut BufReader<std::process::ChildStdout>, id: i64) -> Value {
    let mut line = String::new();
    loop {
        line.clear();
        let n = stdout.read_line(&mut line).expect("read response");
        assert!(n > 0, "server closed stdout before response {id}");
        let value: Value = serde_json::from_str(line.trim_end()).expect("json response");
        if value.get("id").and_then(Value::as_i64) == Some(id) {
            return value;
        }
    }
}

#[test]
fn lists_tools_and_calls_properties() {
    let directory = std::env::temp_dir().join(format!(
        "polytool-stdio-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    std::fs::create_dir(&directory).unwrap();
    let (mut child, mut stdin, mut stdout) = start_server(&directory);

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 1,
            "method": "initialize",
            "params": {
                "protocolVersion": "2025-03-26",
                "capabilities": {},
                "clientInfo": { "name": "polytool-mcp-test", "version": "0.1.0" }
            }
        }),
    );
    let init = read_response(&mut stdout, 1);
    assert!(init.get("result").is_some(), "initialize failed: {init}");
    let instructions = init["result"]["instructions"].as_str().unwrap();
    assert!(instructions.len() <= 400);
    assert!(instructions.contains("report_id"));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "method": "notifications/initialized"
        }),
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 2,
            "method": "tools/list",
            "params": {}
        }),
    );
    let tools = read_response(&mut stdout, 2);
    let tools_array = tools["result"]["tools"].as_array().expect("tool list");
    let tool_names: Vec<&str> = tools_array
        .iter()
        .filter_map(|tool| tool["name"].as_str())
        .collect();
    assert!(tool_names.contains(&"polynomial_properties"));
    assert!(tool_names.contains(&"coefficient_tests"));
    assert!(tool_names.contains(&"hstar_inequalities"));
    assert!(tool_names.contains(&"cyclic_sieving"));
    assert!(tool_names.contains(&"cyclic_sieving_sequence"));
    assert!(tool_names.contains(&"bicyclic_sieving"));
    assert!(tool_names.contains(&"check_polynomial_family"));
    assert!(tool_names.contains(&"find_recurrence"));
    assert!(tool_names.contains(&"enumerate_real_rooted"));
    assert!(tool_names.contains(&"read_report"));
    assert!(tool_names.contains(&"stability_univariate"));
    assert!(tool_names.contains(&"tool_help"));
    let stability_tool = tools_array
        .iter()
        .find(|t| t["name"] == "stability_univariate")
        .unwrap();
    assert!(stability_tool["description"].as_str().unwrap().len() <= 600);
    assert_eq!(
        stability_tool["inputSchema"]["properties"]["polynomials"]["maxItems"],
        512
    );
    assert_eq!(
        stability_tool["inputSchema"]["$defs"]["gaussian"]["oneOf"][0]["type"],
        "integer"
    );
    let description_bytes: usize = tools_array
        .iter()
        .map(|t| t["description"].as_str().unwrap_or("").len())
        .sum();
    let input_bytes: usize = tools_array
        .iter()
        .map(|t| serde_json::to_vec(&t["inputSchema"]).unwrap().len())
        .sum();
    let largest = tools_array
        .iter()
        .map(|t| t["description"].as_str().unwrap_or("").len())
        .max()
        .unwrap();
    for (name, actual) in [
        ("description_bytes", description_bytes),
        ("input_schema_bytes", input_bytes),
        ("largest_description_bytes", largest),
    ] {
        let limit: usize = include_str!("fixtures/tool-list-budget.tsv")
            .lines()
            .find_map(|line| {
                let (key, value) = line.split_once('\t')?;
                (key == name).then(|| value.parse().unwrap())
            })
            .unwrap();
        assert!(actual <= limit, "{name}: {actual} > {limit}");
    }
    let proof_description = tools_array
        .iter()
        .find(|t| t["name"] == "discover_recurrence_proofs")
        .unwrap()["description"]
        .as_str()
        .unwrap();
    assert!(proof_description.len() <= 600);
    for name in ["check_tnn_network", "recognize_brenti_strip_digraph"] {
        let tool = tools_array
            .iter()
            .find(|tool| tool["name"] == name)
            .unwrap_or_else(|| panic!("missing {name}"));
        assert_eq!(tool["inputSchema"]["type"], "object", "{name}");
        assert_eq!(tool["outputSchema"]["type"], "object", "{name}");
        assert!(tool["description"]
            .as_str()
            .is_some_and(|description| description.contains("NOT")));
    }
    let menu_tool = tools_array
        .iter()
        .find(|tool| tool["name"] == "find_recurrence_representations")
        .expect("representation tool");
    assert_eq!(menu_tool["outputSchema"]["type"], "object");
    let recurrence_schema = &tools_array
        .iter()
        .find(|tool| tool["name"] == "find_recurrence")
        .expect("find_recurrence tool")["inputSchema"];
    assert_eq!(recurrence_schema["type"], "object");
    assert!(recurrence_schema.get("oneOf").is_none());
    let recurrence_properties = recurrence_schema["properties"]
        .as_object()
        .expect("find_recurrence properties");
    for name in [
        "coefficients",
        "polynomials",
        "expressions",
        "text",
        "options",
        "include_code",
        "skip_prefix",
        "min_recurrence_length",
        "max_recurrence_length",
        "min_variable_degree",
        "max_variable_degree",
        "min_index_degree",
        "max_index_degree",
        "min_derivative_order",
        "max_derivative_order",
        "try_inhomogeneous",
        "min_inhomogeneous_variable_degree",
        "max_inhomogeneous_variable_degree",
        "min_inhomogeneous_index_degree",
        "max_inhomogeneous_index_degree",
        "try_denominator",
        "try_alternating_sign",
        "max_denominator_variable_degree",
        "max_denominator_index_degree",
        "min_equation_margin",
        "no_held_out_rows",
        "fit_extra_rows",
        "modular_prefilter",
        "max_candidates",
    ] {
        let property = recurrence_properties
            .get(name)
            .unwrap_or_else(|| panic!("missing {name}"));
        assert!(
            property["description"]
                .as_str()
                .is_some_and(|description| !description.is_empty()),
            "{name} lacks a description: {property}"
        );
    }

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 3,
            "method": "tools/call",
            "params": {
                "name": "polynomial_properties",
                "arguments": {
                    "text": "1, 11, 11, 1"
                }
            }
        }),
    );
    let result = read_response(&mut stdout, 3);
    let structured = &result["result"]["structuredContent"];
    assert_eq!(structured["items"][0]["real_rooted"], true);
    assert_eq!(
        structured["items"][0]["gamma_coefficients"],
        json!(["1", "8"])
    );

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 4,
            "method": "tools/call",
            "params": {
                "name": "check_polynomial_family",
                "arguments": {
                    "sequence": "eulerian",
                    "max_n": 4,
                    "options": {
                        "require_real_rooted": true,
                        "require_gamma_positive": true
                    }
                }
            }
        }),
    );
    let family = read_response(&mut stdout, 4);
    let structured = &family["result"]["structuredContent"];
    assert_eq!(structured["all_required_checks_passed"], true);
    assert!(structured["markdown"]
        .as_str()
        .expect("markdown")
        .contains("Polynomial family check"));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 5,
            "method": "tools/call",
            "params": {
                "name": "find_recurrence",
                "arguments": {
                    "coefficients": [[1], [1], [2], [3], [5], [8]],
                    "min_rec_len": 1,
                    "max_rec_len": 2,
                    "max_var_deg": 0,
                    "max_idx_deg": 0,
                    "max_diff_deg": 0,
                    "include_code": false
                }
            }
        }),
    );
    let recurrence = read_response(&mut stdout, 5);
    let structured = &recurrence["result"]["structuredContent"];
    assert_eq!(structured["found"], true);
    assert_eq!(structured["recurrence"], "P(n) = P(n-1) + P(n-2)");
    assert!(structured["latex"].is_string());
    assert!(structured["candidates_considered"].is_number());
    for absent in ["mathematica", "python", "sage", "recurrence_json"] {
        assert!(structured.get(absent).is_none(), "unexpected {absent}");
    }

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 6,
            "method": "tools/call",
            "params": {
                "name": "find_recurrence",
                "arguments": {
                    "coefficients": [[1], [2], [4], [8], [16]],
                    "include_code": true,
                    "options": {
                        "min_rec_len": 1,
                        "max_rec_len": 1,
                        "max_var_deg": 0,
                        "max_idx_deg": 0,
                        "max_diff_deg": 0
                    }
                }
            }
        }),
    );
    let legacy = read_response(&mut stdout, 6);
    assert_eq!(
        legacy["result"]["structuredContent"]["recurrence"],
        "P(n) = 2 P(n-1)"
    );
    assert!(legacy["result"]["structuredContent"]["python"].is_string());

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0",
            "id": 7,
            "method": "tools/call",
            "params": {
                "name": "find_recurrence",
                "arguments": {
                    "coefficients": [[1], [2], [4]],
                    "expressions": ["1", "2", "4"]
                }
            }
        }),
    );
    let invalid = read_response(&mut stdout, 7);
    assert!(invalid["error"]["message"]
        .as_str()
        .is_some_and(|message| message.contains("expected exactly one")));

    send(
        &mut stdin,
        json!({
            "jsonrpc": "2.0", "id": 8, "method": "tools/call",
            "params": {
                "name": "find_recurrence_representations",
                "arguments": {
                    "coefficients": [[1], [2], [4], [8], [16]],
                    "max_rec_len": 1, "max_diff_deg": 0,
                    "max_var_deg": 0, "max_idx_deg": 0, "max_candidates": 1
                }
            }
        }),
    );
    let menu = read_response(&mut stdout, 8);
    let report = &menu["result"]["structuredContent"];
    assert_eq!(report["schema"], "polytool.recurrence-representations.v1");
    assert_eq!(
        report["representations"][0]["recurrence"],
        "P(n) = 2 P(n-1)"
    );
    assert_eq!(report["representations"][0]["evidence"], "exact_fit");

    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":9,"method":"tools/call","params":{"name":"discover_recurrence_proofs","arguments":{"oeis_id":"A001263"}}}),
    );
    let proof = read_response(&mut stdout, 9);
    let summary = &proof["result"]["structuredContent"];
    assert_eq!(summary["schema"], "polytool.recurrence-proof-summary.v1");
    assert!(serde_json::to_vec(summary).unwrap().len() <= 2048);
    let id = &summary["report"]["id"];
    assert!(id.is_string());
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":10,"method":"tools/call","params":{"name":"read_report","arguments":{"report_id":id,"select":["/status","/attempts/*/theorem"]}}}),
    );
    let read = read_response(&mut stdout, 10);
    assert_eq!(
        read["result"]["structuredContent"]["selected"]["/status"],
        "proved"
    );
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":11,"method":"tools/call","params":{"name":"verify_recurrence_proofs","arguments":{"report_id":id}}}),
    );
    let verified = read_response(&mut stdout, 11);
    assert_eq!(verified["result"]["structuredContent"]["verified"], true);

    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":12,"method":"tools/call","params":{"name":"stability_univariate","arguments":{"polynomials":["(t+1)^3","t-I","0"],"region":"left"}}}),
    );
    let stability = read_response(&mut stdout, 12);
    let summary = &stability["result"]["structuredContent"];
    assert_eq!(
        summary["items"][0]["counts"],
        json!({"inside":3,"boundary":0,"outside":0})
    );
    assert_eq!(summary["items"][2]["counts"], Value::Null);
    assert!(summary["items"][0].get("certificate").is_none());
    assert!(serde_json::to_vec(summary).unwrap().len() <= 1100);
    assert!(serde_json::to_vec(&stability).unwrap().len() <= 2500);
    assert!(summary.get("report_id").is_none());
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":15,"method":"tools/call","params":{"name":"stability_univariate","arguments":{"polynomials":["(t+1)^3","t-I","0"],"region":"left","detail":"full"}}}),
    );
    let full = read_response(&mut stdout, 15);
    let id = &full["result"]["structuredContent"]["report_id"];
    assert!(id.is_string());
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":13,"method":"tools/call","params":{"name":"read_report","arguments":{"report_id":id,"select":["/schema","/items/0/certificate/real_part"]}}}),
    );
    let read = read_response(&mut stdout, 13);
    assert_eq!(
        read["result"]["structuredContent"]["selected"]["/schema"],
        "polytool.stability.v1"
    );
    assert!(
        read["result"]["structuredContent"]["selected"]["/items/0/certificate/real_part"]
            .is_array()
    );
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":14,"method":"tools/call","params":{"name":"tool_help","arguments":{"tool":"stability_univariate"}}}),
    );
    let help = read_response(&mut stdout, 14);
    assert!(help["result"]["structuredContent"]["help"]
        .as_str()
        .unwrap()
        .contains("50-second"));

    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":66,"method":"tools/call","params":{
            "name":"check_polynomial_family","arguments":{
                "polynomials":[{"coefficients":[1,1,1]},{"coefficients":[1,2,1]}],
                "options":{"find_recurrence":false}
            }
        }}),
    );
    let family = read_response(&mut stdout, 66);
    assert_ne!(family["result"]["isError"], true, "{family}");
    assert_eq!(
        family["result"]["structuredContent"]["parse_errors"],
        json!([])
    );

    assert!(tool_names.contains(&"discover_interlacing_refinements"));
    send(
        &mut stdin,
        json!({"jsonrpc":"2.0","id":67,"method":"tools/call","params":{
            "name":"discover_interlacing_refinements","arguments":{"request":{
                "states":[[["1"]],[["2"]],[["4"]],[["8"]]],
                "target_rows":[["1"],["2"],["4"],["8"]],
                "derive_candidates":false,
                "limits":{"index_degree":0,"variable_degree":0},
                "candidates":[{"name":"finite scalar fit"}]
            }}
        }}),
    );
    let refinement = read_response(&mut stdout, 67);
    assert_ne!(refinement["result"]["isError"], true, "{refinement}");
    assert_eq!(
        refinement["result"]["structuredContent"]["status"],
        "finite_fit"
    );

    drop(stdin);
    let _ = child.kill();
    let _ = child.wait();
    std::fs::remove_dir_all(directory).unwrap();
}
