use std::io::Write;
use std::process::{Command, Stdio};

fn call(args: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .arg("recurrence-implicit")
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    if let Err(error) = child.stdin.take().unwrap().write_all(input.as_bytes()) {
        // The CLI may exit before reading its input.
        assert_eq!(error.kind(), std::io::ErrorKind::BrokenPipe);
    }
    child.wait_with_output().unwrap()
}

const WEIGHTED_EULERIAN: &str = r#"{"oeis_id":"A065826","row_count":8,"search":{},
  "transformation":{"kind":"diagonal","weight":{"kind":"inverse_binomial","m":1}}}"#;

#[test]
fn cli_recovers_eulerian_rows_and_replays_the_certificates() {
    let output = call(&["--compact"], WEIGHTED_EULERIAN);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-implicit.v1");
    assert_eq!(
        report["transformed_rows"][3],
        serde_json::json!(["1", "11", "11", "1"])
    );
    assert_eq!(report["assessment"]["status"], "equivalent_description");
    assert_eq!(
        report["independent_search"]["evidence"],
        "certified_by_reconstruction"
    );
    let verified = call(&["--verify"], &report.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = report.clone();
    tampered["transformed_rows"][3][1] = serde_json::json!("12");
    assert_eq!(
        call(&["--verify"], &tampered.to_string()).status.code(),
        Some(2)
    );
}

#[test]
fn cli_assesses_a_description_without_its_constant_condition() {
    let request = r#"{"oeis_id":"A007318","row_count":6,
      "transformation":{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}}"#;
    let output = call(&["--compact"], request);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let mut description = report["description"].clone();
    let full = call(&["--assess"], &description.to_string());
    let full: serde_json::Value = serde_json::from_slice(&full.stdout).unwrap();
    assert_eq!(full["status"], "equivalent_description");
    description["exceptional_index_conditions"] = serde_json::json!([]);
    let weaker = call(&["--assess"], &description.to_string());
    let weaker: serde_json::Value = serde_json::from_slice(&weaker.stdout).unwrap();
    assert_eq!(weaker["status"], "certified_consequence");
    assert_eq!(weaker["uniqueness"]["status"], "not_certified");
}

#[test]
fn cli_rejects_other_maps_and_unknown_options() {
    let other = r#"{"oeis_id":"A007318","transformation":{"kind":"affine","a":"1","b":"1"}}"#;
    let output = call(&[], other);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("transform_recurrence"));
    assert_eq!(call(&["--rows"], "{}").status.code(), Some(2));
    let help = call(&["--help"], "");
    assert!(String::from_utf8_lossy(&help.stdout).contains("derivative_antecedent"));
}
