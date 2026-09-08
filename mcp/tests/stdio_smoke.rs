use serde_json::{json, Value};
use std::io::{BufRead, BufReader, Write};
use std::process::{Child, ChildStdin, Command, Stdio};

fn start_server() -> (Child, ChildStdin, BufReader<std::process::ChildStdout>) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool-mcp"))
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
    let (mut child, mut stdin, mut stdout) = start_server();

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
    assert!(tool_names.contains(&"check_polynomial_family"));
    assert!(tool_names.contains(&"find_recurrence"));
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
        "min_rec_len",
        "max_rec_len",
        "min_var_deg",
        "max_var_deg",
        "min_idx_deg",
        "max_idx_deg",
        "min_diff_deg",
        "max_diff_deg",
        "try_inhomogeneous",
        "min_inhomo_var_deg",
        "max_inhomo_var_deg",
        "min_inhomo_idx_deg",
        "max_inhomo_idx_deg",
        "try_denominator",
        "try_alternating_sign",
        "max_denom_var_deg",
        "max_denom_idx_deg",
        "min_margin",
        "no_verify",
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

    drop(stdin);
    let _ = child.kill();
    let _ = child.wait();
}
