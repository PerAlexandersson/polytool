use std::io::Write;
use std::process::{Command, Stdio};

fn call(args: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .arg("recurrence-transform")
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

const EULERIAN_F_POLYNOMIAL: &str = r#"{"oeis_id":"A008292","row_count":7,"search":{},
  "transformation":{"kind":"f_polynomial","schedule":{"first_index":1,"first_degree":0,"step":1}}}"#;

#[test]
fn cli_transports_eulerian_rows_and_replays_the_report() {
    let output = call(&["--compact"], EULERIAN_F_POLYNOMIAL);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-transform.v1");
    assert_eq!(report["transported"]["evidence"], "transported_identity");
    assert_eq!(
        report["transported"]["text"],
        "P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)"
    );
    assert_eq!(report["independent_search"]["evidence"], "finite_fit");
    assert_eq!(
        report["transformed_states"][2][0],
        serde_json::json!(["1", "6", "6"])
    );

    let verified = call(&["--verify"], &report.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = report.clone();
    tampered["transformed_states"][2][0][1] = serde_json::json!("7");
    assert_eq!(
        call(&["--verify"], &tampered.to_string()).status.code(),
        Some(2)
    );
}

#[test]
fn cli_rejects_singular_maps_and_unknown_options() {
    let singular = r#"{"rows":[["1"]],"transformation":{"kind":"mobius","a":"1","b":"2","c":"2",
      "e":"4","schedule":{"first_index":0,"first_degree":1}}}"#;
    let output = call(&[], singular);
    assert_eq!(output.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&output.stderr).contains("singular"));
    assert_eq!(call(&["--rows"], "{}").status.code(), Some(2));
    let help = call(&["--help"], "");
    assert!(String::from_utf8_lossy(&help.stdout).contains("symmetric_decomposition"));
}
