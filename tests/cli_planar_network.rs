//! CLI tests for `polytool tnn-network` and `polytool brenti-recognize`.

use serde_json::Value;
use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run_polytool(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn polytool");
    // The CLI may exit on an argument error before reading its input.
    if let Err(error) = child
        .stdin
        .take()
        .expect("polytool stdin")
        .write_all(input.as_bytes())
    {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().expect("wait for polytool")
}

fn json(output: &Output) -> Value {
    serde_json::from_slice(&output.stdout).expect("stdout is JSON")
}

fn stdout(output: &Output) -> String {
    String::from_utf8(output.stdout.clone()).expect("stdout is utf8")
}

fn strings(values: &[&str]) -> Value {
    Value::from(values.iter().map(|v| Value::from(*v)).collect::<Vec<_>>())
}

const PASCAL: &str = "1,1\n1,2,1\n1,3,3,1\n1,4,6,4,1\n";

#[test]
fn tnn_network_certifies_pascal_with_unit_weights_and_chain_polynomials() {
    let output = run_polytool(&["tnn-network", "--chain", "--json"], PASCAL);
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["schema"], "polytool.tnn-network.v1");
    assert_eq!(report["input_mode"], "polynomials");
    assert_eq!(report["input_rows"], 4);
    assert_eq!(report["prepended_p0"], true);
    assert_eq!(
        report["coefficient_matrix"][4],
        strings(&["1", "4", "6", "4", "1"])
    );
    assert_eq!(report["tnn"], true);
    assert_eq!(report["failure"], Value::Null);
    assert_eq!(report["certificate_verified"], true);
    let lambda = report["lambda"].as_array().unwrap();
    assert_eq!(lambda.len(), 4);
    for (n, row) in lambda.iter().enumerate() {
        let row = row.as_array().unwrap();
        assert_eq!(row.len(), n + 1);
        assert!(row.iter().all(|w| w == "1"), "{row:?}");
    }
    let chain = &report["chain"];
    assert_eq!(chain["error"], Value::Null);
    assert_eq!(chain["conclusion_holds"], true);
    assert_eq!(chain["contradicts_theorem"], false);
    // p_n = sum_k k! S(n,k) t^k.
    assert_eq!(
        chain["polynomials"][4]["coefficients"],
        strings(&["0", "1", "14", "36", "24"])
    );
    let polys = chain["polynomials"].as_array().unwrap();
    assert!(polys
        .iter()
        .all(|p| p["real_rooted"] == true && p["roots_in_minus_one_zero"] == true));
    let interlacing = chain["interlacing"].as_array().unwrap();
    assert_eq!(interlacing.len(), 4);
    assert!(interlacing.iter().all(|pair| pair["interlaces"] == true));

    let text = stdout(&run_polytool(&["tnn-network", "--chain"], PASCAL));
    assert!(text.contains("- Verdict: lower-unitriangular and totally nonnegative (TNN)"));
    assert!(text.contains("  n=3: 1, 1, 1, 1"));
    assert!(text.contains("the network's path matrix equals the input"));
    assert!(text.contains("p_4 = t + 14t^2 + 36t^3 + 24t^4: real-rooted, zeros in [-1, 0]"));
    assert!(text.contains("Consecutive interlacing p_n << p_(n+1): all hold"));
}

#[test]
fn tnn_network_guardrail_example_is_tnn_without_row_claims() {
    // [[1,0,0],[1,1,0],[1,1,1]] is TNN although its last row 1 + t + t^2 has
    // no real zeros: the report must not claim anything about the rows.
    for (args, input) in [
        (vec!["tnn-network", "--json"], "1\n1,1\n1,1,1\n"),
        (
            vec!["tnn-network", "--matrix", "--json"],
            "1,0,0\n1,1,0\n1,1,1\n",
        ),
    ] {
        let output = run_polytool(&args, input);
        assert_eq!(output.status.code(), Some(0));
        let report = json(&output);
        assert_eq!(report["tnn"], true);
        assert_eq!(report["prepended_p0"], false);
        assert_eq!(report["certificate_verified"], true);
        assert_eq!(report["lambda"], serde_json::json!([["1"], ["1", "0"]]));
        assert_eq!(report["chain"], Value::Null);
        let note = report["note"].as_str().unwrap();
        assert!(note.contains("does not imply that the row polynomials are real-rooted"));
        let keys: Vec<&String> = report.as_object().unwrap().keys().collect();
        assert!(
            keys.iter().all(|k| !k.contains("real_rooted")),
            "no row real-rootedness field: {keys:?}"
        );
    }
    let text = stdout(&run_polytool(&["tnn-network"], "1\n1,1\n1,1,1\n"));
    assert!(text.contains("totally nonnegative (TNN)"));
    assert!(text.contains("does not imply that the row polynomials are real-rooted"));
}

#[test]
fn tnn_network_reports_a_negative_minor_witness() {
    let output = run_polytool(
        &["tnn-network", "--matrix", "--json"],
        "1,0,0\n1,1,0\n2,0,1\n",
    );
    assert_eq!(output.status.code(), Some(1));
    let report = json(&output);
    assert_eq!(report["tnn"], false);
    assert_eq!(report["lambda"], Value::Null);
    assert_eq!(report["certificate_verified"], Value::Null);
    let failure = &report["failure"];
    assert_eq!(failure["kind"], "not_tnn");
    assert_eq!(failure["witness"]["rows"], serde_json::json!([1, 2]));
    assert_eq!(failure["witness"]["cols"], serde_json::json!([0, 1]));
    assert_eq!(failure["witness"]["value"], "-2");

    let text = stdout(&run_polytool(
        &["tnn-network", "--matrix"],
        "1,0,0\n1,1,0\n2,0,1\n",
    ));
    assert!(text.contains("- Verdict: not lower-unitriangular TNN"));
    assert!(text.contains("rows [1, 2], columns [0, 1] (zero-based) give the minor -2"));
}

#[test]
fn tnn_network_reports_shape_problems_and_input_errors() {
    let non_monic = json(&run_polytool(&["tnn-network", "--json"], "1,2\n1,2,1\n"));
    assert_eq!(non_monic["failure"]["kind"], "non_monic_polynomial");
    assert_eq!(non_monic["failure"]["row"], 1);
    assert_eq!(non_monic["failure"]["value"], "2");

    let above = run_polytool(&["tnn-network", "--matrix", "--json"], "1,2\n1,1\n");
    assert_eq!(above.status.code(), Some(1));
    let above = json(&above);
    assert_eq!(above["failure"]["kind"], "nonzero_above_diagonal");
    assert_eq!(above["failure"]["col"], 1);
    assert_eq!(above["coefficient_matrix"], Value::Null);

    let negative = json(&run_polytool(
        &["tnn-network", "--matrix", "--chain", "--json"],
        "1\n-1,1\n",
    ));
    assert_eq!(negative["failure"]["kind"], "negative_entry");
    // Chain polynomials exist for any unitriangular shape; the theorem
    // simply does not apply.
    assert_eq!(negative["chain"]["error"], Value::Null);
    assert_eq!(negative["chain"]["contradicts_theorem"], false);

    assert_eq!(run_polytool(&["tnn-network"], "").status.code(), Some(2));
    assert_eq!(
        run_polytool(&["tnn-network", "--bogus"], "1\n")
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        run_polytool(&["tnn-network"], "1,x,\n").status.code(),
        Some(2)
    );
}

fn stirling2_input(max_n: usize) -> String {
    let mut rows: Vec<Vec<u64>> = vec![vec![1]];
    for n in 1..=max_n {
        let prev = &rows[n - 1];
        let row = (0..=n)
            .map(|k| {
                let stay = prev.get(k).copied().unwrap_or(0) * k as u64;
                let up = if k > 0 { prev[k - 1] } else { 0 };
                stay + up
            })
            .collect();
        rows.push(row);
    }
    rows[1..]
        .iter()
        .map(|row| {
            row.iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",")
                + "\n"
        })
        .collect()
}

#[test]
fn brenti_recognize_certifies_pascal_and_stirling() {
    let output = run_polytool(&["brenti-recognize", "--json"], PASCAL);
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["schema"], "polytool.brenti-recognize.v1");
    assert_eq!(report["verdict"], "certified");
    assert_eq!(report["certificate_verified"], true);
    assert_eq!(report["row_indices"], serde_json::json!([1, 2, 3, 4]));
    assert_eq!(
        report["digraph"]["outgoing"][0],
        serde_json::json!([{"to": 0, "weight": "1"}, {"to": 1, "weight": "1"}])
    );

    let output = run_polytool(
        &["brenti-recognize", "--weights", "integer", "--json"],
        &stirling2_input(6),
    );
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["verdict"], "certified");
    assert_eq!(report["options"]["weights"], "integer");
    let outgoing = report["digraph"]["outgoing"].as_array().unwrap();
    assert_eq!(outgoing[0], serde_json::json!([{"to": 1, "weight": "1"}]));
    // S(n+1,k) = S(n,k-1) + k S(n,k): y -> y with weight y, y -> y + 1.
    for (y, edges) in outgoing.iter().enumerate().take(6).skip(1) {
        assert_eq!(
            *edges,
            serde_json::json!([
                {"to": y, "weight": y.to_string()},
                {"to": y + 1, "weight": "1"}
            ])
        );
    }

    let text = stdout(&run_polytool(&["brenti-recognize"], PASCAL));
    assert!(text.contains("- Verdict: certified"));
    assert!(text.contains("  0 -> 0: 1, 1: 1"));
    assert!(text.contains("It does not show that the rows are real-rooted"));
}

#[test]
fn brenti_recognize_rejects_eulerian_rows_within_bounds() {
    let eulerian = "1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n";
    let output = run_polytool(&["brenti-recognize", "--json"], eulerian);
    assert_eq!(output.status.code(), Some(1));
    let report = json(&output);
    assert_eq!(report["verdict"], "impossible_within_bounds");
    assert_eq!(report["certified"], false);
    assert_eq!(report["digraph"], Value::Null);
    assert_eq!(report["supports_checked"], 0);
    assert!(report["obstruction"]
        .as_str()
        .unwrap()
        .contains("no admissible weights"));
}

#[test]
fn brenti_recognize_certificate_is_not_a_real_rootedness_claim() {
    let output = run_polytool(&["brenti-recognize", "--json"], "1,1\n1,1,1\n");
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["verdict"], "certified");
    assert!(report["note"]
        .as_str()
        .unwrap()
        .contains("does not show that the rows are real-rooted"));
}

#[test]
fn brenti_recognize_budgets_and_errors() {
    let output = run_polytool(
        &["brenti-recognize", "--max-search-nodes", "1", "--json"],
        "1,1,1\n1,2,3,2,1\n",
    );
    assert_eq!(output.status.code(), Some(3));
    assert_eq!(json(&output)["verdict"], "inconclusive");

    let span = json(&run_polytool(
        &["brenti-recognize", "--max-span", "2", "--json"],
        PASCAL,
    ));
    assert_eq!(span["verdict"], "certified");
    assert_eq!(span["options"]["max_target_span"], 2);

    for args in [
        vec!["brenti-recognize", "--first-row", "2"],
        vec!["brenti-recognize", "--first-row", "0"],
        vec!["brenti-recognize", "--max-span", "0"],
        vec!["brenti-recognize", "--weights", "real"],
        vec!["brenti-recognize", "--max-span"],
    ] {
        let output = run_polytool(&args, PASCAL);
        assert_eq!(output.status.code(), Some(2), "{args:?}");
        assert!(output.stdout.is_empty());
    }
    assert_eq!(
        run_polytool(&["brenti-recognize"], "").status.code(),
        Some(2)
    );
}

#[test]
fn planar_network_commands_have_help() {
    for command in ["tnn-network", "brenti-recognize"] {
        let output = run_polytool(&["help", command], "");
        assert_eq!(output.status.code(), Some(0));
        let text = stdout(&output);
        assert!(text.contains("Guardrails:"), "{command}");
        assert!(text.contains("Exit status:"), "{command}");
    }
    let top = stdout(&run_polytool(&["--help"], ""));
    assert!(top.contains("tnn-network"));
    assert!(top.contains("brenti-recognize"));
}
