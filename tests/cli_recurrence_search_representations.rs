//! CLI `recurrence-search-representations`: JSON round trip and replay.

use polytool::recurrence::{verify_representation_search, RepresentationSearchReport};
use std::io::Write;
use std::process::{Command, Stdio};

fn call(args: &[&str], input: &str) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
        .arg("recurrence-search-representations")
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

const FACTORIAL_BINOMIAL: &str = r#"{"oeis_id":"A007318","row_count":10,
  "limits":{"max_depth":1},
  "transformations":[{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}]}"#;

#[test]
fn cli_search_round_trips_and_replays() {
    let output = call(&["--compact"], FACTORIAL_BINOMIAL);
    assert!(output.status.success(), "{output:?}");
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(
        value["schema"],
        "polytool.recurrence-representation-search.v1"
    );
    assert_eq!(value["status"], "complete");
    let labels: Vec<&str> = value["representations"][1]["descriptions"]
        .as_array()
        .unwrap()
        .iter()
        .map(|d| d["label"].as_str().unwrap())
        .collect();
    assert!(labels.contains(&"equivalent_generating_description"));
    // The CLI JSON is the crate report.
    let report: RepresentationSearchReport = serde_json::from_value(value.clone()).unwrap();
    verify_representation_search(&report).unwrap();
    assert_eq!(serde_json::to_value(&report).unwrap(), value);

    let verified = call(&["--verify"], &value.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let summary: serde_json::Value = serde_json::from_slice(&verified.stdout).unwrap();
    assert_eq!(summary["verified"], true);
    assert_eq!(summary["representations"], 2);

    let mut tampered = value;
    tampered["representations"][1]["states"][2][0][1] = "3".into();
    let rejected = call(&["--verify"], &tampered.to_string());
    assert_eq!(rejected.status.code(), Some(2));

    let invalid = call(&[], r#"{"transformations": []}"#);
    assert_eq!(invalid.status.code(), Some(2));
    let help = call(&["--help"], "");
    assert!(String::from_utf8_lossy(&help.stdout).contains("duplicate_map"));
}

const EULERIAN_VECTOR: &str = r#"{"oeis_id":"A008292","row_count":8,
  "limits":{"max_depth":2,"max_representations":8},
  "transformations":[
    {"kind":"symmetric_decomposition","schedule":{"first_index":1,"first_degree":1,"step":1}},
    {"kind":"f_polynomial","schedule":{"first_index":1,"first_degree":1,"step":1}}]}"#;

#[test]
fn cli_transforms_vector_representations_and_replays() {
    let output = call(&["--compact"], EULERIAN_VECTOR);
    assert!(output.status.success(), "{output:?}");
    let value: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let report: RepresentationSearchReport = serde_json::from_value(value.clone()).unwrap();
    verify_representation_search(&report).unwrap();
    // Decomposition, then f-polynomials of both components.
    let node = report
        .representations
        .iter()
        .find(|n| n.chain == [0, 1])
        .unwrap();
    assert_eq!(node.components, ["a", "b"]);
    let step = node.steps[0].vector_transport.as_ref().unwrap();
    assert_eq!(step.evidence, "transported_identity");
    assert!(node.descriptions.iter().any(|d| {
        d.kind == "transported_vector_system"
            && d.label == "equivalent_generating_description"
            && d.vector_assessment.as_ref().unwrap().status == "certified"
    }));
    let verified = call(&["--verify"], &value.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = value;
    tampered["representations"][node.id]["steps"][0]["vector_transport"]["evidence"] =
        "checked_supplied_rows".into();
    let rejected = call(&["--verify"], &tampered.to_string());
    assert_eq!(rejected.status.code(), Some(2));
    let help = call(&["--help"], "");
    let help = String::from_utf8_lossy(&help.stdout);
    assert!(help.contains("componentwise") && help.contains("implicit_spaces"));
}
