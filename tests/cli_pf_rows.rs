//! CLI tests for `polytool pf-network` and `polytool left-production`.

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

const PASCAL: &str = "1,1\n1,2,1\n1,3,3,1\n1,4,6,4,1\n";
const EULERIAN: &str = "1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n";
const STIRLING: &str = "0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n0,1,15,25,10,1\n";

#[test]
fn pf_network_certifies_pascal_and_first_eulerian_rows() {
    let output = run_polytool(&["pf-network", "--json"], PASCAL);
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["schema"], "polytool.pf-network.v1");
    assert_eq!(report["verdict"], "certified");
    assert_eq!(report["certificate_verified"], true);
    assert_eq!(report["network"]["last_column"], 4);
    let first = &report["network"]["steps"][0];
    assert_eq!(first[0]["length"], 1);
    assert_eq!(first[0]["rise"], 0);
    assert_eq!(first[1]["rise"], 1);
    assert_eq!(first[1]["weight"], "1");
    assert!(report["note"]
        .as_str()
        .unwrap()
        .contains("Pólya frequency sequence"));

    // A_1..A_3 = Brenti's pattern; A_3 = (1 + t) A_2 + 6 t A_1.
    let three: String = EULERIAN.lines().take(3).map(|l| format!("{l}\n")).collect();
    let output = run_polytool(&["pf-network", "--json", "--weights", "integer"], &three);
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    let third = report["network"]["steps"][2].as_array().unwrap();
    assert!(third
        .iter()
        .any(|s| s["length"] == 2 && s["rise"] == 1 && s["weight"] == "6"));
}

#[test]
fn pf_network_rejects_non_real_rooted_and_bounded_rows() {
    for input in ["1,1,1\n", "1,1\n1,1,1\n", "1,2\n1,1,1\n"] {
        let output = run_polytool(&["pf-network", "--json"], input);
        assert_eq!(output.status.code(), Some(1), "{input}");
        let report = json(&output);
        assert_eq!(report["verdict"], "impossible_within_bounds");
        assert!(report["obstruction"]
            .as_str()
            .unwrap()
            .contains("not real-rooted"));
    }
    let output = run_polytool(
        &["pf-network", "--max-length", "4", "--max-rise", "4"],
        EULERIAN,
    );
    assert_eq!(output.status.code(), Some(1));
    assert!(stdout(&output).contains("impossible within bounds"));
}

#[test]
fn pf_network_budgets_and_errors() {
    let output = run_polytool(&["pf-network", "--max-search-nodes", "1"], PASCAL);
    assert_eq!(output.status.code(), Some(3));
    assert_eq!(
        run_polytool(&["pf-network", "--max-length", "0"], PASCAL)
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        run_polytool(&["pf-network", "--first-row", "0"], PASCAL)
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        run_polytool(&["pf-network", "--bogus"], PASCAL)
            .status
            .code(),
        Some(2)
    );
}

#[test]
fn left_production_certifies_stirling_and_pascal() {
    let output = run_polytool(&["left-production", "--json"], STIRLING);
    assert_eq!(output.status.code(), Some(0));
    let report = json(&output);
    assert_eq!(report["schema"], "polytool.left-production.v1");
    assert_eq!(report["tnn"], true);
    assert_eq!(report["rows_certified_pf"], true);
    assert_eq!(report["certificate_verified"], true);
    assert_eq!(
        report["production_matrix"][5],
        Value::from(vec!["0", "1", "4", "6", "4", "1"])
    );
    assert!(report["production_network"]["lambda"].is_array());

    let output = run_polytool(&["left-production"], PASCAL);
    assert_eq!(output.status.code(), Some(0));
    assert!(stdout(&output).contains("rows P_0..P_4 are PF"));
}

#[test]
fn left_production_reports_failures_without_row_claims() {
    let output = run_polytool(&["left-production", "--json"], EULERIAN);
    assert_eq!(output.status.code(), Some(1));
    let report = json(&output);
    assert_eq!(report["tnn"], false);
    assert_eq!(report["rows_certified_pf"], false);
    assert_eq!(report["production_matrix"][4][2], "-5");
    assert_eq!(report["negative_minor"]["value"], "-5");
    assert!(report["message"]
        .as_str()
        .unwrap()
        .contains("says nothing about the rows"));

    let output = run_polytool(&["left-production", "--first-row", "0"], "1\n1,1\n1,1,1\n");
    assert_eq!(output.status.code(), Some(1));
    assert!(stdout(&output).contains("Witness"));

    // deg P_n must equal n.
    assert_eq!(
        run_polytool(&["left-production"], "1,1,1\n").status.code(),
        Some(2)
    );
    assert_eq!(
        run_polytool(&["left-production", "--bogus"], PASCAL)
            .status
            .code(),
        Some(2)
    );
}

#[test]
fn pf_row_commands_have_help() {
    for command in ["pf-network", "left-production"] {
        let output = run_polytool(&["help", command], "");
        assert_eq!(output.status.code(), Some(0));
        let text = stdout(&output);
        assert!(text.contains("Guardrails:"), "{command}");
        assert!(text.contains("Exit status:"), "{command}");
    }
    let top = stdout(&run_polytool(&["--help"], ""));
    assert!(top.contains("pf-network"));
    assert!(top.contains("left-production"));
}
