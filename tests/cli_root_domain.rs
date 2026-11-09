//! CLI tests for root-domain analysis of recurrences.

mod common;
use common::run_polytool_subcommand as call;

const EULERIAN_ROOT_DOMAIN: &str = r#"{"oeis_id":"A008292","row_count":6,"root_domain":{},
  "transformation":{"kind":"f_polynomial","schedule":{"first_index":1,"first_degree":0,"step":1}}}"#;

#[test]
fn cli_root_domain_closes_the_eulerian_induction_and_replays() {
    let output = call("recurrence-transform", &["--compact"], EULERIAN_ROOT_DOMAIN);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    let root = &report["root_domain"];
    assert_eq!(root["recurrence_theorem"]["matched"], true);
    assert_eq!(
        root["recurrence_theorem"]["factor"],
        serde_json::json!(["0", "1", "1"])
    );
    let evidence: Vec<&str> = root["claims"]
        .as_array()
        .unwrap()
        .iter()
        .map(|c| c["evidence"].as_str().unwrap())
        .collect();
    assert!(evidence.contains(&"completed_recurrence_induction"));
    assert!(evidence.contains(&"certified_family_invariant"));
    assert!(evidence.contains(&"proposed_interval"));
    assert!(evidence.contains(&"finite_exact_observation"));

    let verified = call("recurrence-transform", &["--verify"], &report.to_string());
    assert!(verified.status.success(), "{verified:?}");
    let mut tampered = report.clone();
    tampered["root_domain"]["claims"][0]["evidence"] =
        serde_json::json!("certified_family_invariant");
    assert_eq!(
        call("recurrence-transform", &["--verify"], &tampered.to_string())
            .status
            .code(),
        Some(2)
    );

    // Without the field the report has no root-domain section.
    let plain = EULERIAN_ROOT_DOMAIN.replace(r#""root_domain":{},"#, "");
    let output = call("recurrence-transform", &["--compact"], &plain);
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert!(report.get("root_domain").is_none());
}

#[test]
fn cli_root_enclosures_reports_multiplicities_and_interval_tests() {
    let request = r#"{"polynomials":[["0","0","1","3","3","1"],["1","1","1"],[]],
      "interval":{"lower":{"value":"-1"},"upper":{"value":"0","closed":false}}}"#;
    let output = call("root-enclosures", &["--compact"], request);
    assert!(output.status.success(), "{output:?}");
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.root-enclosures.v1");
    let items = report["items"].as_array().unwrap();
    assert_eq!(items[0]["valuation"], 2);
    assert_eq!(items[0]["roots"][0]["multiplicity"], 3);
    assert_eq!(items[0]["roots"][1]["lower"], "0");
    // The root 0 sits on the open endpoint.
    assert_eq!(items[0]["real_rooted_in_interval"], false);
    assert_eq!(items[1]["real_rooted"], false);
    assert_eq!(items[1]["real_roots_in_interval"], true);
    assert!(items[2].get("degree").is_none());

    let bad = call(
        "root-enclosures",
        &[],
        r#"{"polynomials":[["1"]],"max_width":"-1"}"#,
    );
    assert_eq!(bad.status.code(), Some(2));
    assert_eq!(
        call("root-enclosures", &["--rows"], "{}").status.code(),
        Some(2)
    );
    let help = call("root-enclosures", &["--help"], "");
    assert!(String::from_utf8_lossy(&help.stdout).contains("multiplicity"));
}
