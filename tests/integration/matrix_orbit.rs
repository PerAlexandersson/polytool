//! Matrix-orbit CLI, supplied-system and tamper regressions.
use super::common::run_polytool;
use polytool::recurrence::{check_matrix_orbit, MatrixOrbitRequest, RecurrenceJson};
use serde_json::{json, Value};

fn request(id: &str) -> Value {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    json!({"recurrence":RecurrenceJson::from_recurrence_rational(&rule,entry.recurrence_first_index,&seeds,None)})
}
#[test]
fn cli_automatically_recognizes_and_replays() {
    for id in ["A120987", "A304429", "A330509", "A340620", "A245963"] {
        let output = run_polytool(&["matrix-orbit"], &request(id).to_string());
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let report: Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["status"], "proved");
        assert_eq!(report["lean_replay"], "not_attempted");
        let replay = run_polytool(&["matrix-orbit", "--verify"], &report.to_string());
        assert!(replay.status.success());
        let mut changed = report;
        changed["certificate"]["pencils"][0]["zero_side"] = json!(!changed["certificate"]
            ["pencils"][0]["zero_side"]
            .as_bool()
            .unwrap());
        assert_eq!(
            run_polytool(&["matrix-orbit", "--verify"], &changed.to_string())
                .status
                .code(),
            Some(2)
        );
    }
    let compact = run_polytool(
        &["matrix-orbit", "--format", "json", "--compact"],
        &request("A245963").to_string(),
    );
    assert!(compact.status.success());
    assert_eq!(
        String::from_utf8(compact.stdout).unwrap().lines().count(),
        1
    );
}
#[test]
fn discovery_records_real_rootedness_without_temporal_claims() {
    let input = request("A245963")["recurrence"].to_string();
    let output = run_polytool(&["recurrence-prove", "--json"], &input);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    let attempt = report["attempts"]
        .as_array()
        .unwrap()
        .iter()
        .find(|a| a["theorem"] == "matrix_orbit")
        .unwrap();
    assert_eq!(attempt["status"], "proved");
    assert!(attempt["conclusion"].get("interlacing").is_none());
    assert!(run_polytool(
        &["recurrence-prove", "--verify", "--json"],
        &report.to_string()
    )
    .status
    .success());
}
#[test]
fn supplied_system_checks_identity_and_invalid_shapes() {
    let mut input = request("A245963");
    let automatic =
        check_matrix_orbit(&serde_json::from_value::<MatrixOrbitRequest>(input.clone()).unwrap())
            .unwrap();
    input["state_system"] =
        serde_json::to_value(automatic.certificate.unwrap().state_system).unwrap();
    let output = run_polytool(&["matrix-orbit"], &input.to_string());
    assert!(output.status.success());
    input["state_system"]["start_vector"][0][0] = json!("2/2");
    let noncanonical = run_polytool(&["matrix-orbit"], &input.to_string());
    assert!(noncanonical.status.success());
    assert!(run_polytool(
        &["matrix-orbit", "--verify"],
        &String::from_utf8(noncanonical.stdout).unwrap()
    )
    .status
    .success());
    input["state_system"]["start_vector"][0][0] = json!("999");
    assert_eq!(
        run_polytool(&["matrix-orbit"], &input.to_string())
            .status
            .code(),
        Some(1)
    );
    input["state_system"]["projection"] = json!(["1"]);
    assert_eq!(
        run_polytool(&["matrix-orbit"], &input.to_string())
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        run_polytool(&["matrix-orbit"], r#"{"recurrence":{},"unknown":1}"#)
            .status
            .code(),
        Some(2)
    );
}
#[test]
fn rejects_the_path_cube_supplied_realization() {
    let rec = json!({"schema":"polytool.recurrence.v1","first_index":0,"initial_polynomials":[["1"],["1","1"],["1","4","1"],["1","11","11","1"]],"recurrence":{"terms":[
        {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1","3"]]}},
        {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","2","-2"]]}},
        {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","2","-1","-1"]]}},
        {"offset":4,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","0","1","-2","1"]]}}
    ]}});
    let input = json!({"recurrence":rec,"state_system":{
        "matrix":[[["1","1"],["0","2"],["0"],["0"]],[["1"],["0","1"],["2","1"],["0","3"]],[["1"],["0","1"],["0"],["0"]],[["0"],["0"],["1"],["0","1"]]],
        "start_vector":[["1"],["0"],["0"],["0"]],"projection":["1","0","0","0"]}});
    let output = run_polytool(&["matrix-orbit"], &input.to_string());
    assert_eq!(output.status.code(), Some(1));
    let report: Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "not_proved");
    assert!(
        run_polytool(&["matrix-orbit", "--verify"], &report.to_string())
            .status
            .success()
    );
}
