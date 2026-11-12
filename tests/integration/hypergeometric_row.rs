//! Tests of `recurrence::check_hypergeometric_row_theorem` and its parent
//! route: row-hypergeometric catalog families in the identity coordinates,
//! after the reversal, with a valuation, with seed rows that do not follow
//! the rule and with a non-hypergeometric first coefficient; families outside
//! the multiplier scheme (not real-rooted, alternating, not hypergeometric);
//! and replay with tampered reports and parent records.

use polytool::oeis;
use polytool::recurrence::{
    check_hypergeometric_row_theorem, check_recorded_parent_proof,
    find_recurrence_adaptive_rational, prove_parent_real_rootedness, AdaptiveSearchOptions,
    BigRational, BivarPoly, LiftHypothesisStatus, ParentRecurrenceRecord, RealRootednessTheorem,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, RealRootednessTheoremStatus,
    Recurrence, RecurrenceTerm, VeroneseParentProofRoute, VeroneseParentProofStatus,
};

use crate::common::q;

fn catalog_report(id: &str) -> RealRootednessTheoremReport {
    let family = oeis::by_id(id).expect("catalog family");
    let (rule, seeds) = family.recurrence_parts().expect("catalog recurrence");
    check_hypergeometric_row_theorem(
        &rule,
        family.recurrence_first_index,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .expect("the check runs")
}

/// The proved report of a catalog family, replayed, with its region
/// `(-inf, 0)` or `(-inf, 0]`, and the parent provers proving the family by
/// this route (no earlier route proves these families), with a replayable
/// record.
fn assert_proved(id: &str, closed_at_zero: bool) -> RealRootednessTheoremReport {
    let report = catalog_report(id);
    assert_eq!(
        report.status,
        RealRootednessTheoremStatus::Proved,
        "{id}: {:?}",
        report.open_hypotheses()
    );
    assert_eq!(report.theorem, RealRootednessTheorem::HypergeometricRow);
    report.verify().unwrap_or_else(|e| panic!("{id}: {e}"));
    let region = report
        .conclusion
        .as_ref()
        .and_then(|c| c.root_region.clone())
        .expect("a root region")
        .at(0)
        .unwrap();
    assert!(region.lower.is_none());
    assert_eq!(region.upper, Some(q(0)));
    assert_eq!(region.upper_closed, closed_at_zero, "{id}");
    let family = oeis::by_id(id).unwrap();
    let (rule, seeds) = family.recurrence_parts().unwrap();
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: family.recurrence_first_index,
    };
    let proof =
        prove_parent_real_rootedness(&record, family.recurrence_first_index, &[], true).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved, "{id}");
    // Rows proportional to powers of one polynomial (A124860, A128540) are
    // proved earlier in the cascade by the Liu--Wang three-term route.
    assert!(
        matches!(
            proof.route,
            Some(
                VeroneseParentProofRoute::HypergeometricRow
                    | VeroneseParentProofRoute::LiuWangThreeTerm
            )
        ),
        "{id}: {:?}",
        proof.route
    );
    check_recorded_parent_proof(&record, family.recurrence_first_index, &[], true, &proof).unwrap();
    report
}

#[test]
fn jacobi_interpolation_rows_are_proved() {
    // A330798: sum_k binom(n, k) (2n + 1)_k / (n + 2)_k x^k, the base
    // binom(n, k) with the paired polynomial (k + n + 2) ... (k + 2n), from
    // n = 1 on (the pair has negative length at n = 0, checked exactly).
    let report = assert_proved("A330798", false);
    let names: Vec<&str> = report.hypotheses.iter().map(|h| h.name.as_str()).collect();
    assert!(names.contains(&"row_relation"));
    assert!(names.iter().any(|name| name.starts_with("pair_")));
    assert!(report.normal_form.contains("(k + 2 n)"));
}

#[test]
fn identity_coordinate_families_are_proved() {
    // binom(n, k) (k + 1)(k + 2)(k + 3) / 6: a pair of length 3 (Laguerre).
    assert_proved("A124848", false);
    // binom(n, k) binom(2n, n - k): a reciprocal rising factorial.
    assert_proved("A110608", false);
    assert_proved("A286785", false);
    assert_proved("A123160", false);
    // Jacobsthal--Pascal rows: the first coefficient is not hypergeometric
    // in n, its two-lag scalar recurrence has positive coefficients.
    assert_proved("A124860", false);
    assert_proved("A376827", false);
}

#[test]
fn reversed_and_shifted_families_are_proved() {
    // The reversal at the degree d(n), with the degree bound of the rule.
    assert_proved("A119308", false);
    assert_proved("A357613", false);
    // Rows divisible by x (valuation 1): closed at 0.
    assert_proved("A253283", true);
    assert_proved("A359363", true);
    // Valuation n - 1, a third seed row that does not follow the rule.
    assert_proved("A128540", true);
}

#[test]
fn families_outside_the_scheme_are_not_proved() {
    // binom(n, k)^2 binom(2(n - k), n - k) (a double Narayana transform) and
    // binom(n, k) (n + 2k + 3)! ... : hypergeometric rows whose factors do
    // not match the multiplier menu.
    for id in ["A318397", "A087727"] {
        let report = catalog_report(id);
        assert_eq!(
            report.status,
            RealRootednessTheoremStatus::Unresolved,
            "{id}"
        );
        assert!(report.hypotheses.iter().any(
            |h| h.name == "multiplier_matching" && h.status == LiftHypothesisStatus::Unresolved
        ));
        report.verify().unwrap();
    }
    // The Eulerian numbers: no hypergeometric row relation.
    let report = catalog_report("A008292");
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
}

#[test]
fn a_hypergeometric_family_that_is_not_real_rooted_is_not_proved() {
    // sum_k binom(n, k) binom(2k, k) x^k: 1 + 4x + 6x^2 at n = 2 has no real
    // root; the factor (k - 1/2) has no partner.
    let binomial = |n: i64, k: i64| -> i64 { (0..k).fold(1i64, |acc, i| acc * (n - i) / (i + 1)) };
    let rows: Vec<Vec<BigRational>> = (0..16)
        .map(|n| {
            (0..=n)
                .map(|k| q(binomial(n, k) * binomial(2 * k, k)))
                .collect()
        })
        .collect();
    let fit = find_recurrence_adaptive_rational(&rows, &AdaptiveSearchOptions::default())
        .expect("a recurrence");
    let rule = fit.recurrence;
    let seeds = rows[..rule.max_offset()].to_vec();
    let report = check_hypergeometric_row_theorem(
        &rule,
        0,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    assert!(!report.is_proved(), "{:?}", report.status);
}

#[test]
fn alternating_rows_fail_the_scale() {
    // P_n = (1 - x) P_(n-1): the ratio -(n + 1 - k)/k has the wrong sign.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(-1)]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let report = check_hypergeometric_row_theorem(
        &rule,
        0,
        &[vec![q(1)]],
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    assert!(!report.is_proved());
}

#[test]
fn tampered_reports_and_records_are_rejected() {
    let report = assert_proved("A124848", false);
    let mut tampered = report.clone();
    tampered.normal_form.push_str(" (edited)");
    assert!(tampered.verify().is_err());
    let mut tampered = report.clone();
    tampered.hypotheses[0].evidence = "trust me".into();
    assert!(tampered.verify().is_err());
    let mut tampered = report.clone();
    if let Some(conclusion) = tampered.conclusion.as_mut() {
        conclusion.root_region = None;
    }
    assert!(tampered.verify().is_err());

    // A parent record claiming this route for a family it does not prove.
    let family = oeis::by_id("A318397").unwrap();
    let (rule, seeds) = family.recurrence_parts().unwrap();
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: family.recurrence_first_index,
    };
    let mut proof =
        prove_parent_real_rootedness(&record, family.recurrence_first_index, &[], true).unwrap();
    assert_ne!(proof.status, VeroneseParentProofStatus::Proved);
    proof.status = VeroneseParentProofStatus::Proved;
    proof.route = Some(VeroneseParentProofRoute::HypergeometricRow);
    assert!(
        check_recorded_parent_proof(&record, family.recurrence_first_index, &[], true, &proof)
            .is_err()
    );
}
