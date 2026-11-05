//! CLI discovery, explicit menu selection, certificate replay and usage failures.

use std::io::Write;
use std::process::{Command, Output, Stdio};

fn run(args: &[&str], input: &str) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Err(error) = child.stdin.take().unwrap().write_all(input.as_bytes()) {
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().unwrap()
}

fn input() -> serde_json::Value {
    serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
            "coeff": {"coeffs": [["1", "1"], ["1"]]}}]}
    })
}

#[test]
fn discover_and_replay_then_reject_tampering() {
    let result = run(
        &["recurrence-prove", "--recurrence", "-", "--json"],
        &input().to_string(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut report: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report["status"], "proved");
    assert_eq!(report["attempts"][0]["theorem"], "product_recurrence");
    assert_eq!(report["lean_replay"], "not_attempted");
    let verified = run(
        &["recurrence-prove", "--verify", "--json"],
        &report.to_string(),
    );
    assert!(verified.status.success());
    report["attempts"][0]["conclusion"]["simple_roots"] = serde_json::json!(true);
    let bad = run(&["recurrence-prove", "--verify"], &report.to_string());
    assert_eq!(bad.status.code(), Some(2));
    assert!(bad.stdout.is_empty());
}

#[test]
fn explicit_menu_selection_and_no_proof_status() {
    let mut bad = input();
    bad["recurrence"]["terms"][0]["coeff"]["coeffs"] = serde_json::json!([["1", "0", "1"]]);
    let menu = serde_json::json!({"schema": "polytool.recurrence-representations.v1",
        "representations": [{"recurrence_json": input()}, {"recurrence_json": bad}]});
    let ambiguous = run(
        &["recurrence-prove", "--recurrence", "-"],
        &menu.to_string(),
    );
    assert_eq!(ambiguous.status.code(), Some(2));
    let yes = run(
        &["recurrence-prove", "--recurrence", "-", "--candidate", "0"],
        &menu.to_string(),
    );
    assert!(yes.status.success());
    assert!(String::from_utf8_lossy(&yes.stdout).contains("ProductRecurrence: Proved"));
    let no = run(
        &[
            "recurrence-prove",
            "--recurrence",
            "-",
            "--candidate",
            "1",
            "--json",
        ],
        &menu.to_string(),
    );
    assert_eq!(no.status.code(), Some(1));
    let report: serde_json::Value = serde_json::from_slice(&no.stdout).unwrap();
    assert_eq!(report["status"], "not_proved");
    assert!(run(&["recurrence-prove", "--verify"], &report.to_string())
        .status
        .success());
}

#[test]
fn help_and_invalid_flags_do_not_consume_input() {
    assert!(run(&["help", "recurrence-prove"], "").status.success());
    for args in [
        vec!["recurrence-prove", "--unknown"],
        vec!["recurrence-prove", "--verify", "--recurrence", "-"],
        vec![
            "recurrence-prove",
            "--max-start-delay",
            "17",
            "--recurrence",
            "-",
        ],
    ] {
        let result = run(&args, "invalid JSON");
        assert_eq!(result.status.code(), Some(2));
        assert!(result.stdout.is_empty());
        assert!(!String::from_utf8_lossy(&result.stderr).contains("parse recurrence JSON"));
    }
}
