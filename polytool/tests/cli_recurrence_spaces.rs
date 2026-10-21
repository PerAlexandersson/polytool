use std::io::Write;
use std::process::{Command, Stdio};

fn call(flag: &str, input: &str) -> std::process::Output {
    let mut p = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .args(["recurrence", flag])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    p.stdin.take().unwrap().write_all(input.as_bytes()).unwrap();
    p.wait_with_output().unwrap()
}

#[test]
fn cli_complete_search_and_certificate_replay() {
    let req = r#"{"source":{"coefficients":[["2"]],"seeds":[["1"]]},"shapes":[{"kind":"box","length":1}]}"#;
    let output = call("--space", req);
    assert!(output.status.success(), "{:?}", output);
    let report: String = String::from_utf8(output.stdout).unwrap();
    assert!(report.contains("conditional_on_source"));
    assert!(call("--verify-space", &report).status.success());
    let mut bad: serde_json::Value = serde_json::from_str(&report).unwrap();
    bad["results"][0]["h"] = serde_json::json!(999);
    assert_eq!(
        call("--verify-space", &bad.to_string()).status.code(),
        Some(2)
    );
}

#[test]
fn cli_resource_exhaustion_and_bad_source_are_not_absence() {
    let req = r#"{"source":{"coefficients":[["2"]],"seeds":[["1"]]},"shapes":[{"kind":"box","length":1}],"limits":{"max_operations":0}}"#;
    let output = call("--space", req);
    assert_eq!(output.status.code(), Some(3));
    assert!(String::from_utf8(output.stdout)
        .unwrap()
        .contains("incomplete_budget"));
    assert_eq!(call("--space", "{}").status.code(), Some(2));
}
