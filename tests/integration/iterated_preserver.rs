//! Tests of `recurrence::check_iterated_preserver_theorem` and its parent
//! route: first-order rules `P_n = a_n(x) S_n(P_(n-1))` after seed rows that
//! do not follow the rule (OEIS catalog records), index-dependent
//! multipliers decided by the sign checker, menu operators, a multiplier
//! that is not real-rooted, and replay with tampered certificates.

use polytool::oeis;
use polytool::recurrence::{
    check_iterated_preserver_theorem, prove_parent_real_rootedness, BigRational, BivarPoly,
    IteratedPreserverCertificate, LiftHypothesisStatus, ParentRecurrenceRecord, PreserverOperator,
    RealRootednessTheorem, RealRootednessTheoremOptions, RealRootednessTheoremReport,
    RealRootednessTheoremStatus, Recurrence, RecurrenceTerm, RouteCertificate, RowInterlacing,
    VeroneseParentProofRoute, VeroneseParentProofStatus,
};

use crate::common::{q, rows};

/// The iterated-preserver certificate of a report, for tampering.
fn iterated_preserver_mut(
    report: &mut RealRootednessTheoremReport,
) -> &mut IteratedPreserverCertificate {
    match report.certificate.as_mut() {
        Some(RouteCertificate::IteratedPreserver(certificate)) => certificate,
        _ => panic!("the report has no iterated-preserver certificate"),
    }
}

/// `coeffs[i][j]` = coefficient of `n^i x^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows(coeffs))
}

/// `P_n = sum_d c_d(n, x) D^d P_(n-1)` from `(d, c_d)`.
fn first_order(terms: &[(usize, &[&[i64]])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(order, coeffs)| RecurrenceTerm::new(1, *order, bivar(coeffs)))
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

fn catalog(id: &str) -> (Recurrence, usize, Vec<Vec<BigRational>>) {
    let entry = oeis::by_id(id).expect("catalog entry");
    let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
    (rule, entry.recurrence_first_index, seeds)
}

/// The check with default options, replayed (also through JSON), with a
/// tampered status and a tampered certificate refused.
fn check(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
) -> RealRootednessTheoremReport {
    let report = check_iterated_preserver_theorem(
        rule,
        first_index,
        seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .expect("the check runs");
    assert_eq!(report.theorem, RealRootednessTheorem::IteratedPreserver);
    report.verify().expect("the report replays");
    let json = serde_json::to_string(&report).unwrap();
    let back: RealRootednessTheoremReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    back.verify().unwrap();
    let mut tampered = report.clone();
    tampered.status = if report.is_proved() {
        RealRootednessTheoremStatus::Unresolved
    } else {
        RealRootednessTheoremStatus::Proved
    };
    assert!(tampered.verify().is_err());
    if report.iterated_preserver_certificate().is_some() {
        // A different multiplier: both the recomputation and the replay of
        // the operator identity refuse it.
        let mut tampered = report.clone();
        let certificate = iterated_preserver_mut(&mut tampered);
        certificate.multiplier_fixed_factor.push("1".into());
        let rule_only = certificate.verify_operator_identity(rule);
        assert!(rule_only.is_err(), "the tampered identity replays");
        assert!(tampered.verify().is_err());
        // A proved report without its certificate is refused.
        if report.is_proved() {
            let mut stripped = report.clone();
            stripped.certificate = None;
            assert!(stripped.verify().is_err());
        }
    }
    report
}

fn region_text(report: &RealRootednessTheoremReport) -> Option<String> {
    report
        .conclusion
        .as_ref()
        .and_then(|c| c.root_region.as_ref())
        .map(ToString::to_string)
}

#[test]
fn three_pascal_rows_follow_the_multiplier_from_the_last_seed() {
    // A028262: P_n = (1 + x) P_(n-1) for n >= 3 after 1, 1 + x, 1 + 3x + x^2.
    let (rule, first_index, seeds) = catalog("A028262");
    assert_eq!(seeds.len(), 3);
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.first_index, 0);
    assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"));
    // A linear multiplier: weak interlacing from the rule on, extended over
    // the seeds as far as it is checked exactly.
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
    assert!(conclusion.interlacing_from.unwrap() <= 3);
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert_eq!(certificate.induction_start, 2);
    assert_eq!(certificate.start_degree, 2);
    assert_eq!(certificate.degree_step, 1);
    assert_eq!(certificate.multiplier_fixed_factor, vec!["1", "1"]);
    assert!(certificate.preserver.is_none());
    certificate.verify_operator_identity(&rule).unwrap();
}

#[test]
fn catalog_multipliers_with_two_or_more_seeds() {
    // A029635 (two seeds, multiplier 1 + x) and A024462 (multiplier 1 + 3x).
    for (id, seeds_expected) in [("A029635", 2), ("A024462", 3)] {
        let (rule, first_index, seeds) = catalog(id);
        assert_eq!(seeds.len(), seeds_expected, "{id}");
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
        assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"), "{id}");
    }
    let (rule, first_index, seeds) = catalog("A024462");
    let report = check(&rule, first_index, &seeds);
    assert_eq!(
        report
            .iterated_preserver_certificate()
            .unwrap()
            .multiplier_fixed_factor,
        vec!["1/3", "1"]
    );
}

#[test]
fn index_dependent_multipliers_use_the_sign_checker() {
    // A051150: P_n = (x - 5n - 10) P_(n-1): one positive root -> [0, inf).
    let (rule, first_index, seeds) = catalog("A051150");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(region_text(&report).as_deref(), Some("[0, +inf)"));
    let decided = report
        .hypotheses
        .iter()
        .filter(|h| h.name.starts_with("multiplier_index_factor"))
        .collect::<Vec<_>>();
    assert!(!decided.is_empty());
    assert!(decided.iter().any(|h| h.sign_decision.is_some()));
    // A216118: (n - 1) P_n = (n + 3)(x + n) P_(n-1) from a zero seed on.
    let (rule, first_index, seeds) = catalog("A216118");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"));
    // P_0 = 0: the conclusion starts after the zero seed row.
    assert_eq!(report.conclusion.as_ref().unwrap().first_index, 1);
    assert!(report
        .hypotheses
        .iter()
        .any(|h| h.name == "left_factor_nonzero" && h.status == LiftHypothesisStatus::Verified));
}

#[test]
fn an_index_dependent_root_crossing_zero_gives_the_line_only() {
    // P_n = (x + 3 - n) P_(n-1): the root n - 3 leaves (-inf, 0] at n = 4
    // and is negative for n <= 2.
    let rule = first_order(&[(0, &[&[3, 1], &[-1]])]);
    let report = check(&rule, 0, &rows(&[&[1, 1]]));
    assert!(report.is_proved(), "{report:#?}");
    // Real-rooted on the line only.
    assert_eq!(region_text(&report), None);
    let options = RealRootednessTheoremOptions::default();
    let ray = RealRootednessTheoremOptions {
        root_region: Some(polytool::real_rootedness::RealInterval::at_most(q(0))),
        ..options
    };
    let report = check_iterated_preserver_theorem(&rule, 0, &rows(&[&[1, 1]]), &ray).unwrap();
    assert!(!report.is_proved());
    report.verify().unwrap();
}

#[test]
fn a_theta_polynomial_iterates() {
    // P_n = (theta + 1)(theta + 2) P_(n-1) = (2 + 4 x D + x^2 D^2) P_(n-1).
    let rule = first_order(&[(0, &[&[2]]), (1, &[&[0, 4]]), (2, &[&[0, 0, 1]])]);
    let seeds = rows(&[&[1, 3, 3, 1]]);
    let report = check(&rule, 0, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"));
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert!(matches!(
        certificate.preserver,
        Some(PreserverOperator::EulerPolynomial { .. })
    ));
    assert_eq!(certificate.degree_step, 0);
    assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
    assert!(!report.lean_declarations.is_empty());
}

#[test]
fn a_multiplier_after_a_menu_operator() {
    // P_n = (1 + x)(theta + 1) P_(n-1) = (1 + x) P + (x + x^2) P'.
    let rule = first_order(&[(0, &[&[1, 1]]), (1, &[&[0, 1, 1]])]);
    let report = check(&rule, 0, &rows(&[&[1]]));
    assert!(report.is_proved(), "{report:#?}");
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert_eq!(certificate.multiplier_fixed_factor, vec!["1", "1"]);
    assert!(certificate.preserver.is_some());
    // The central factorial rows A036969: (theta + 1)^2 + x, an Euler step.
    let (rule, first_index, seeds) = catalog("A036969");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert!(matches!(
        report.iterated_preserver_certificate().unwrap().preserver,
        Some(PreserverOperator::EulerStep { .. })
    ));
}

#[test]
fn euler_steps_use_the_degree_local_raising_condition() {
    // A166960, A166961, A166962, A166972: Euler steps with v < 0 and
    // u(n) + v deg P_(n-1) > 0 (the degree-local hypothesis).
    for id in ["A166960", "A166961", "A166962", "A166972"] {
        let (rule, first_index, seeds) = catalog(id);
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
        assert!(matches!(
            report.iterated_preserver_certificate().unwrap().preserver,
            Some(PreserverOperator::EulerStep { .. })
        ));
        assert!(report
            .hypotheses
            .iter()
            .any(|h| h.name == "euler_step_raising_positive_up_to_degree"));
        assert!(report
            .hypotheses
            .iter()
            .all(|h| h.name != "euler_step_raising_slope_nonnegative"));
        assert!(report.reference.contains("generalStep_of_natDegree_le"));
        assert!(report
            .lean_declarations
            .iter()
            .any(|d| d == "RealRooted.EulerBidiagonal.generalRowsDep_spec"));
    }
    // The other catalog Euler steps (n-dependent a = b for A371081).
    for id in [
        "A036969", "A071951", "A080248", "A160562", "A371081", "A371259", "A390433",
    ] {
        let (rule, first_index, seeds) = catalog(id);
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
    }
}

#[test]
fn a_multiplier_that_is_not_real_rooted_gives_no_claim() {
    // P_n = (1 + x^2) P_(n-1).
    let rule = first_order(&[(0, &[&[1, 0, 1]])]);
    let report = check(&rule, 0, &rows(&[&[1, 1]]));
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    assert!(report.conclusion.is_none());
    assert!(report.hypotheses.iter().any(|h| {
        h.name == "multiplier_fixed_factor_in_region"
            && matches!(h.status, LiftHypothesisStatus::Failed { .. })
    }));
    // A seed that is not real-rooted fails with its index.
    let rule = first_order(&[(0, &[&[1, 1]])]);
    let report = check(&rule, 0, &rows(&[&[1, 0, 1]]));
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    // A second-order rule does not apply.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[0, 1]])),
            RecurrenceTerm::new(2, 0, bivar(&[&[1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report = check(&rule, 0, &rows(&[&[1], &[1, 1]]));
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(report.iterated_preserver_certificate().is_none());
}

#[test]
fn a_tampered_region_or_degree_step_is_refused() {
    let (rule, first_index, seeds) = catalog("A028262");
    let report = check(&rule, first_index, &seeds);
    let mut tampered = report.clone();
    iterated_preserver_mut(&mut tampered).degree_step = 2;
    assert!(tampered.verify().is_err());
    assert!(tampered
        .iterated_preserver_certificate()
        .unwrap()
        .verify_operator_identity(&rule)
        .is_err());
    let mut tampered = report.clone();
    let certificate = iterated_preserver_mut(&mut tampered);
    certificate.region.upper = Some(polytool::recurrence::RootEndpoint::rational(&q(-1), true));
    assert!(tampered.verify().is_err());
    // A different preserver for the same rule.
    let mut tampered = report.clone();
    iterated_preserver_mut(&mut tampered).preserver =
        Some(PreserverOperator::DerivativePolynomial {
            coefficients: vec![polytool::recurrence::IndexRationalFunction::constant(&q(1))],
        });
    assert!(tampered.verify().is_err());
}

#[test]
fn the_parent_provers_use_the_iterated_route() {
    // P_n = (1 + x)(theta + 1)^2 P_(n-1), a second-order differential step
    // (outside the Liu--Wang and Ma--Wang shapes).
    let rule = first_order(&[(0, &[&[1, 1]]), (1, &[&[0, 3, 3]]), (2, &[&[0, 0, 1, 1]])]);
    let seeds = rows(&[&[1]]);
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index: 0,
        },
        0,
        &[],
        false,
    )
    .unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(
        proof.route,
        Some(VeroneseParentProofRoute::IteratedPreserver)
    );
    assert!(proof.interlacing.is_none());
    // x P_(n-1) after the seeds 0, 1 (A129184): the root 0 is outside
    // (-inf, 0), the region of Liu--Wang, so the iterated route proves it;
    // the multiplier is linear, so the rows interlace.
    let (rule, first_index, seeds) = catalog("A129184");
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        },
        first_index,
        &[],
        false,
    )
    .unwrap();
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Proved,
        "{proof:#?}"
    );
    assert_eq!(
        proof.route,
        Some(VeroneseParentProofRoute::IteratedPreserver)
    );
    assert!(proof.interlacing.is_some());
    // A region [0, inf) does not give the parent statement (one-signed
    // coefficients).
    let (rule, first_index, seeds) = catalog("A051150");
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        },
        first_index,
        &[],
        false,
    )
    .unwrap();
    assert_ne!(proof.status, VeroneseParentProofStatus::Proved);
    // A second-order parent records no route-8 attempt (the theorem does not
    // apply), so its parent proof is that of the routes before it.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[0, 2]])),
            RecurrenceTerm::new(2, 0, bivar(&[&[0, 0, -1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = rows(&[&[0], &[1]]);
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &rule,
            seeds: &seeds,
            first_index: 0,
        },
        0,
        &[],
        false,
    )
    .unwrap();
    assert!(proof
        .attempts
        .iter()
        .all(|a| a.route != VeroneseParentProofRoute::IteratedPreserver));
}
