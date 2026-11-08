//! The degree-law certificates on catalog families, with the example of
//! `doc/library.md`, and their use by the theorem checks.

use num_bigint::BigInt;
use polytool::recurrence::{
    certify_leading_coefficient, check_hypergeometric_row_theorem,
    check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem, BigRational, BivarPoly,
    LeadingCoefficientOptions, LeadingCoefficientSign, LiftHypothesis, LiftHypothesisStatus,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, RealRootednessTheoremStatus,
    Recurrence, RecurrenceTerm,
};

/// The example of `doc/library.md`.
#[test]
fn library_documentation_example() {
    let q = |v: i64| BigRational::from_integer(BigInt::from(v));
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]),
            ),
            RecurrenceTerm::new(
                3,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let certificate = certify_leading_coefficient(
        &rule,
        0,
        &[vec![q(1)], vec![q(1)], vec![q(1)]],
        &LeadingCoefficientOptions::default(),
    )
    .unwrap();
    assert!(certificate.is_proved());
    assert_eq!(certificate.degree_at(31), Some(10));
    assert_eq!(
        certificate.sign_at(31),
        Some(LeadingCoefficientSign::Positive)
    );
    certificate.verify().unwrap();
}

/// A theorem check of `polytool::recurrence`.
type TheoremCheck =
    fn(
        &Recurrence,
        usize,
        &[Vec<BigRational>],
        &RealRootednessTheoremOptions,
    ) -> Result<RealRootednessTheoremReport, polytool::recurrence::RecurrenceError>;

fn catalog_report(id: &str, check: TheoremCheck) -> RealRootednessTheoremReport {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let report = check(
        &rule,
        entry.recurrence_first_index,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    report.verify().unwrap();
    report
}

fn degree_law_hypothesis(report: &RealRootednessTheoremReport) -> &LiftHypothesis {
    report
        .hypotheses
        .iter()
        .find(|h| h.name == "degree_law")
        .expect("a shared degree-law hypothesis")
}

/// The derivative-lag check takes its degree law from the shared module
/// where a lag multiplier of the top coefficients is negative (A054142,
/// before: `lag_degree_step_nonnegative` failed); the report replays.
#[test]
fn derivative_lag_uses_the_shared_degree_law() {
    for id in ["A054142", "A084534", "A121314"] {
        let report = catalog_report(id, check_liu_wang_derivative_lag_theorem);
        assert!(report.is_proved(), "{id}: {:#?}", report.open_hypotheses());
        let hypothesis = degree_law_hypothesis(&report);
        assert_eq!(hypothesis.status, LiftHypothesisStatus::Verified);
        assert!(hypothesis
            .evidence
            .contains("shared degree-law certificate"));
    }
}

/// The Liu--Wang check admits higher x-degrees when the shared module
/// proves growth one (A008970: the top terms cancel); before, it was not
/// applicable.
#[test]
fn liu_wang_admits_cancelling_top_terms() {
    let report = catalog_report("A008970", check_liu_wang_theorem);
    assert_ne!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert_eq!(
        degree_law_hypothesis(&report).status,
        LiftHypothesisStatus::Verified
    );
}

/// The hypergeometric-row check proves its first coefficient by the shared
/// module where it is neither hypergeometric nor of positive recurrence
/// (A134346: `U_0(n) = 2^(n+1) - 1 = 3 U_0(n-1) - 2 U_0(n-2)`, a growth-ratio
/// invariant); before, `first_coefficient` was unresolved.
#[test]
fn hypergeometric_row_first_coefficient_by_the_shared_module() {
    let report = catalog_report("A134346", check_hypergeometric_row_theorem);
    assert!(report.is_proved(), "{:#?}", report.open_hypotheses());
    let hypothesis = report
        .hypotheses
        .iter()
        .find(|h| h.name == "first_coefficient")
        .unwrap();
    assert!(
        hypothesis.evidence.contains("route recurrence"),
        "{}",
        hypothesis.evidence
    );
}
