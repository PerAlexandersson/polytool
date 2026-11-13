//! Exact Toeplitz region normalization, conditional status and replay.

use polytool::bkw::certified::{
    certify_recurrence, certify_toeplitz_regions, CertificationOptions, LimitSetStatus,
};
use polytool::recurrence::RecurrenceJson;
use serde_json::json;

fn tails(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| s.to_string()).collect()
}

#[test]
fn rational_invariants_cover_scaling_reflection_and_boundaries() {
    let first = certify_toeplitz_regions(&tails(&["3", "1"]))
        .unwrap()
        .unwrap();
    assert_eq!(first.limit_real, Some(true));
    assert_eq!(first.off_real_tie.as_ref().unwrap().point, tails(&["1/27"]));
    first.verify().unwrap();
    for b in [["3", "-1"], ["12", "8"], ["3", "1"]] {
        let certificate = certify_toeplitz_regions(&tails(&b)).unwrap().unwrap();
        assert_eq!(certificate.limit_real, first.limit_real);
        assert_eq!(
            certificate.off_real_tie.as_ref().unwrap().point,
            first.off_real_tie.as_ref().unwrap().point
        );
    }
    assert_eq!(
        certify_toeplitz_regions(&tails(&["3", "2"]))
            .unwrap()
            .unwrap()
            .limit_real,
        Some(false)
    );
    assert_eq!(
        certify_toeplitz_regions(&tails(&["0", "1"]))
            .unwrap()
            .unwrap()
            .limit_real,
        Some(false)
    );
    assert_eq!(
        certify_toeplitz_regions(&tails(&["-1", "1"]))
            .unwrap()
            .unwrap()
            .limit_real,
        Some(false)
    );
    let lower = certify_toeplitz_regions(&tails(&["3", "1", "0"]))
        .unwrap()
        .unwrap();
    assert_eq!(lower, first);
    assert!(certify_toeplitz_regions(&tails(&["1", "0", "0", "1"]))
        .unwrap()
        .is_none());
    let mut tampered = first;
    tampered.off_real_tie.as_mut().unwrap().point[0] = "1".into();
    assert!(tampered.verify().is_err());
}

#[test]
fn closed_form_diagonal_and_tridiagonal_geometry_is_exact() {
    for (b, real) in [
        (&[][..], true),
        (&["0"][..], true),
        (&["1"][..], true),
        (&["-1"][..], false),
    ] {
        let certificate = certify_toeplitz_regions(&tails(b)).unwrap().unwrap();
        assert_eq!(certificate.limit_real, Some(real));
        certificate.verify().unwrap();
    }
}

#[test]
fn circle_theorem_proves_beyond_the_precomputed_geometry_range() {
    // A genuine bandwidth-four family, with a circle radius supplied.
    let record:RecurrenceJson=serde_json::from_value(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"],["0","1"],["-10","0","1"],["0","-20","0","1"],["100","0","-30","0","1"]],
        "recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}},
            {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-10"]]}},
            {"offset":5,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1"]]}}
        ]}})).unwrap();
    let report = certify_recurrence(
        &record,
        &CertificationOptions {
            maximum_row_index: 12,
            circle_radius: Some("1".into()),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.result, "proved");
    assert_eq!(
        report.proof_route.as_deref(),
        Some("hessenberg_toeplitz_circle")
    );
    assert_eq!(report.circle.as_ref().unwrap().radius, "1");
    assert!(report.conditional_on.is_empty());
    assert!(report.toeplitz_region.is_none());
    assert_eq!(report.limit_set_status, LimitSetStatus::Inconclusive);
    report.verify().unwrap();
}

#[test]
fn quartic_boundary_geometry_uses_the_limit_hypothesis_when_no_rational_radius_is_found() {
    let record: RecurrenceJson =
        serde_json::from_value(json!({"schema":"polytool.recurrence.v1","first_index":0,
        "initial_polynomials":[["1"],["0","1"],["-1","0","1"],["0","-2","0","1"]],
        "recurrence":{"terms":[
            {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","1"]]}},
            {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1"]]}},
            {"offset":4,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["1/12"]]}}
        ]}}))
        .unwrap();
    let report = certify_recurrence(
        &record,
        &CertificationOptions {
            maximum_row_index: 48,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(report.result, "conditional");
    assert_eq!(report.limit_set_status, LimitSetStatus::RealLimitSet);
    assert!(report.circle.is_none());
    assert!(report.conditional_on[0].starts_with("hessenberg_real_limit_spectrum:"));
    report.verify().unwrap();
}
