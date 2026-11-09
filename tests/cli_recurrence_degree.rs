//! CLI degree laws: a periodic law, certificate replay, tampering, a
//! not-proved family and usage failures.

mod common;
use common::run_polytool as run;

/// `P_n = P_(n-1) + x P_(n-3)`, `P_0 = P_1 = P_2 = 1`: `deg P_n = floor(n/3)`.
fn input() -> serde_json::Value {
    serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"], ["1"], ["1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1"]]}},
            {"offset": 3, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"]]}}
        ]}
    })
}

#[test]
fn periodic_law_replays_and_tampering_is_refused() {
    let result = run(
        &["recurrence-degree", "--recurrence", "-", "--json"],
        &input().to_string(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let mut certificate: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(
        certificate["schema"],
        "polytool.recurrence-leading-coefficient.v1"
    );
    assert_eq!(certificate["status"], "proved");
    assert_eq!(certificate["route"], "cone");
    assert_eq!(certificate["degree_law"]["period"], 3);
    let verified = run(
        &["recurrence-degree", "--verify", "--json"],
        &certificate.to_string(),
    );
    assert!(verified.status.success());
    certificate["degree_law"]["growth_per_period"] = serde_json::json!(2);
    let bad = run(&["recurrence-degree", "--verify"], &certificate.to_string());
    assert_eq!(bad.status.code(), Some(2));
    assert!(bad.stdout.is_empty());

    let text = run(
        &["recurrence-degree", "--recurrence", "-"],
        &input().to_string(),
    );
    assert!(text.status.success());
    let text = String::from_utf8(text.stdout).unwrap();
    assert!(text.contains("Degree law: Proved"), "{text}");
    assert!(text.contains("Related Lean: route cone"), "{text}");
}

#[test]
fn a_vanishing_leader_exits_one_and_bad_flags_exit_two() {
    // P_n = ((n - 50) x + 1) P_(n-1): the leader vanishes at n = 50.
    let input = serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
            "coeff": {"coeffs": [["1", "-50"], ["0", "1"]]}}]}
    });
    let result = run(
        &["recurrence-degree", "--recurrence", "-", "--json"],
        &input.to_string(),
    );
    assert_eq!(result.status.code(), Some(1));
    let certificate: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(certificate["status"], "not_proved");
    for args in [
        &[
            "recurrence-degree",
            "--recurrence",
            "-",
            "--max-period",
            "0",
        ][..],
        &["recurrence-degree", "--verify", "--rows", "3"][..],
        &["recurrence-degree", "--frobnicate"][..],
    ] {
        let result = run(args, &input.to_string());
        assert_eq!(result.status.code(), Some(2), "{args:?}");
        assert!(result.stdout.is_empty());
    }
    // Without --recurrence the record is read from stdin (or --input FILE).
    let bare = run(&["recurrence-degree"], "invalid JSON");
    assert_eq!(bare.status.code(), Some(2));
    assert!(bare.stdout.is_empty());
    let bare = run(&["recurrence-degree", "--json"], &input.to_string());
    assert_eq!(bare.status.code(), Some(1));
    assert_eq!(bare.stdout, result.stdout);
    let help = run(&["recurrence-degree", "--help"], "");
    assert!(help.status.success());
    assert!(String::from_utf8(help.stdout)
        .unwrap()
        .contains("rr_row_leadingCoeff_pos"));
}
