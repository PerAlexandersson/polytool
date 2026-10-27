//! CLI tests for implicit recurrence transformations.

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

#[test]
fn cli_searches_implicit_relations_directly_and_replays_them() {
    let request = r#"{"oeis_id":"A007318",
      "transformation":{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}},
      "supports":[{"role":"implicit_successor","max_shift":1,"derivative_order":1},
        {"role":"row_equation","derivative_order":2,"index_degree":1,"variable_degree":1}]}"#;
    let output = call(&["--search", "--compact"], request);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-implicit-search.v1");
    let implicit = &report["results"][0]["representatives"][0];
    assert_eq!(implicit["relation"]["text"], "D Q[n+1] - D Q[n] - Q[n] = 0");
    assert_eq!(implicit["assessment"]["status"], "equivalent_description");
    let row = &report["results"][1]["representatives"][0];
    assert_eq!(
        row["relation"]["text"],
        "x D^2 Q[n] + (1 + x) D Q[n] - n Q[n] = 0"
    );
    assert_eq!(row["relation"]["role"], "row_equation");
    let verified = call(&["--verify-search"], &report.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = report.clone();
    tampered["results"][1]["representatives"][0]["relation"]["evidence"] =
        serde_json::json!("finite_fit");
    assert_eq!(
        call(&["--verify-search"], &tampered.to_string())
            .status
            .code(),
        Some(2)
    );
    assert_eq!(
        call(&["--search", "--verify"], request).status.code(),
        Some(2)
    );
}

#[test]
fn cli_attaches_and_replays_implicit_root_domain() {
    let request = r#"{"oeis_id":"A007318","row_count":6,"root_domain":{},
      "transformation":{"kind":"diagonal","weight":{"kind":"factorial_power","power":1}}}"#;
    let output = call(&["--compact"], request);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let root = &report["root_domain"];
    assert_eq!(
        root["implications"][0]["theorem"],
        "polya_schur_multiplier_sequence"
    );
    assert_eq!(root["jensen_fixed_cap"]["evidence"], "fixed_cap");
    let verified = call(&["--verify"], &report.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = report.clone();
    tampered["root_domain"]["claims"] = serde_json::json!([]);
    assert_eq!(
        call(&["--verify"], &tampered.to_string()).status.code(),
        Some(2)
    );
}

#[test]
fn cli_certifies_complete_spaces_after_mobius_maps_and_weights() {
    for (map, theorem) in [
        (
            r#"{"kind":"reversal","schedule":{"first_index":0,"first_degree":0,"step":1}}"#,
            "transported-source theorem",
        ),
        (
            r#"{"kind":"diagonal","weight":{"kind":"binomial","m":1}}"#,
            "weight theorem",
        ),
    ] {
        let request = format!(
            r#"{{"oeis_id":"A007318","transformation":{map},
              "supports":[{{"role":"successor","max_shift":1,"index_degree":1,"variable_degree":2}}]}}"#
        );
        let output = call(&["--search", "--compact"], &request);
        assert!(output.status.success(), "{output:?}");
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert!(report["source"]["completeness_basis"]
            .as_str()
            .unwrap()
            .contains(theorem));
        assert_eq!(report["results"][0]["evidence"], "certified_complete_space");
        let verified = call(&["--verify-search"], &report.to_string());
        assert!(verified.status.success(), "{verified:?}");
        let mut tampered = report;
        tampered["results"][0]["residual_count"] = serde_json::json!(1);
        assert_eq!(
            call(&["--verify-search"], &tampered.to_string())
                .status
                .code(),
            Some(2)
        );
    }
}
