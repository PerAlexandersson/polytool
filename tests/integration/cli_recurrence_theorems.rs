//! CLI: conditional theorem matches stay separate from finite-fit evidence.

/// Rows of `1/(1 - x R(t x))` with `R = 1 + 3z + 2z^2` (labels 0..9).
const RENEWAL_ROWS: &[u8] = b"1\n1\n1,3\n1,6,2\n1,9,13\n1,12,33,12\n1,15,62,63,4\n\
1,18,100,180,66\n1,21,147,390,321,36\n1,24,203,720,985,360,8\n";

/// Rows of `Q_n = Q_(n-1) + n t Q_(max(0, n-2))`.
const PREDECESSOR_ROWS: &[u8] =
    b"1\n1,1\n1,3\n1,6,3\n1,10,15\n1,15,45,15\n1,21,105,105\n1,28,210,420,105\n";

use crate::common::run_polytool_bytes as run;

const SEARCH: [&str; 7] = [
    "recurrence",
    "--max-var-deg",
    "2",
    "--max-idx-deg",
    "0",
    "--max-rec-len",
    "3",
];

#[test]
fn default_json_is_unchanged_and_match_is_noted_on_stderr() {
    let output = run(&[&SEARCH[..], &["--json"]].concat(), RENEWAL_ROWS);
    assert!(output.status.success());
    let record: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(record["schema"], "polytool.recurrence.v1");
    assert!(record.get("theorem_matches").is_none());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("conditional theorem match weighted_renewal_window"));

    // The supplied record is checked by recurrence-theorems.
    let checked = run(
        &["recurrence-theorems", "--recurrence", "-", "--json"],
        &output.stdout,
    );
    assert!(checked.status.success());
    let report: serde_json::Value = serde_json::from_slice(&checked.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-theorems.v1");
    assert_eq!(report["finite_evidence"]["source"], "supplied_recurrence");
    let theorem = &report["theorem_matches"][0];
    assert_eq!(theorem["theorem"], "weighted_renewal_window");
    assert_eq!(theorem["window"]["max_step"], 1);
    assert_eq!(theorem["boundary"]["status"], "canonical_after_relabel");
}

#[test]
fn theorem_matches_flag_emits_separate_report() {
    let output = run(
        &[&SEARCH[..], &["--json", "--theorem-matches"]].concat(),
        RENEWAL_ROWS,
    );
    assert!(output.status.success());
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["schema"], "polytool.recurrence-theorems.v1");
    let plain = run(&[&SEARCH[..], &["--json"]].concat(), RENEWAL_ROWS);
    let plain: serde_json::Value = serde_json::from_slice(&plain.stdout).unwrap();
    assert_eq!(report["recurrence_json"], plain);
    assert_eq!(report["finite_evidence"]["source"], "adaptive_search_fit");
    assert_eq!(report["finite_evidence"]["status"], "finite_evidence_only");
    assert_eq!(report["finite_evidence"]["rows_compared"], 10);
    let theorem = &report["theorem_matches"][0];
    assert_eq!(theorem["status"], "conditional_theorem_match");
    assert_eq!(theorem["weighted_renewal"]["seed"][2], "2");
    assert_eq!(theorem["proof_status"], "paper_proof_not_lean_checked");

    let output = run(&["recurrence", "--theorem-matches"], PREDECESSOR_ROWS);
    assert!(output.status.success());
    let text = String::from_utf8(output.stdout).unwrap();
    assert!(text.contains("Conditional theorem match: fixed_delay_predecessor_window"));
    assert!(text.contains("Not matched weighted_renewal_window"));
}

#[test]
fn recurrence_theorems_exit_status() {
    let fibonacci = run(&["recurrence", "--json"], b"1\n1\n2\n3\n5\n8\n13\n");
    let none = run(
        &["recurrence-theorems", "--recurrence", "-"],
        &fibonacci.stdout,
    );
    assert_eq!(none.status.code(), Some(1));
    assert!(String::from_utf8(none.stdout)
        .unwrap()
        .contains("No conditional theorem match."));
    let invalid = run(&["recurrence-theorems", "--recurrence", "-"], b"{}");
    assert_eq!(invalid.status.code(), Some(2));
}

#[test]
fn product_zero_seeds_use_the_shared_tail_and_residue_proofs() {
    use polytool::recurrence::{
        discover_recurrence_proofs, RealRootednessTheoremOptions, RecurrenceJson,
    };
    let cases = [
        serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["1"],["0"]],"recurrence":{"terms":[
                {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}}]}}),
        serde_json::json!({"schema":"polytool.recurrence.v1","first_index":0,
            "initial_polynomials":[["0"],["1","1"]],"recurrence":{
                "denominator":{"coeffs":[["-1"],["1"]]},"terms":[
                {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","2"],["2","1"],["1","0"]]}}]}}),
    ];
    for record in cases {
        let record: RecurrenceJson = serde_json::from_value(record).unwrap();
        let proof =
            discover_recurrence_proofs(&record, &RealRootednessTheoremOptions::default()).unwrap();
        let input = serde_json::to_vec(&record).unwrap();
        let output = run(&["recurrence-theorems", "--json"], &input);
        assert_eq!(output.status.code(), Some(1)); // Existing match-based exit convention.
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["trivially_real_rooted_if_rule_holds"], true);
        assert_eq!(report["trivially_real_rooted"], true);
        assert_eq!(
            report.get("real_rooted_from_index").cloned(),
            proof.real_rooted_from_index.map(|n| serde_json::json!(n))
        );
        if let Some(from) = proof.real_rooted_from_index {
            assert_eq!(report["real_rooted_from_index"], from);
        }
        assert_eq!(
            report
                .get("real_rooted_outside_zero_residue_classes")
                .cloned(),
            proof
                .real_rooted_outside_zero_residue_classes
                .map(|r| serde_json::to_value(r).unwrap())
        );
        let text = run(&["recurrence-theorems"], &input);
        assert!(String::from_utf8(text.stdout)
            .unwrap()
            .contains("Product proof if the rule holds:"));
        let mut fitted = serde_json::to_value(&record).unwrap();
        // Reuse a valid fitted record's search metadata; the conditional flag is
        // still available while the supplied-rule flag is absent.
        let fitted_record = run(&["recurrence", "--json"], b"1\n2\n4\n8\n16\n32\n64\n");
        let search: serde_json::Value = serde_json::from_slice(&fitted_record.stdout).unwrap();
        fitted["search"] = search["search"].clone();
        let fitted = run(
            &["recurrence-theorems", "--json"],
            &serde_json::to_vec(&fitted).unwrap(),
        );
        let fitted: serde_json::Value = serde_json::from_slice(&fitted.stdout).unwrap();
        assert!(fitted.get("trivially_real_rooted").is_none());
        assert_eq!(fitted["trivially_real_rooted_if_rule_holds"], true);
    }
}

#[test]
fn product_zero_seeds_do_not_make_invalid_live_rows_trivial() {
    for live in [serde_json::json!(["0"]), serde_json::json!(["1", "0", "1"])] {
        let record = serde_json::json!({"schema":"polytool.recurrence.v1", "first_index":0,
            "initial_polynomials":[live,["0"]],"recurrence":{"terms":[
                {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}}]}});
        let output = run(
            &["recurrence-theorems", "--json"],
            &serde_json::to_vec(&record).unwrap(),
        );
        assert_eq!(output.status.code(), Some(1));
        let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(report["trivially_real_rooted_if_rule_holds"], false);
        assert!(report.get("real_rooted_from_index").is_none());
        assert!(report
            .get("real_rooted_outside_zero_residue_classes")
            .is_none());
    }
}
