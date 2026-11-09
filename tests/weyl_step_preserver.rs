//! Tests of Brändén's Weyl step `A + B D + C D^2` in the iterated-preserver
//! theorem (`recurrence::check_iterated_preserver_theorem`) and of
//! `recurrence::check_weyl_step_preserver`: the catalog rows of Brändén's
//! Corollary 3.5, index-dependent coefficients, single-lag rules
//! `P_n = T(P_(n-2))`, the coordinates (power gauge, translation, reversal),
//! the Corollary 2.5 interlacing, refusals and tampered certificates.

use num_bigint::BigInt;
use polytool::oeis;
use polytool::recurrence::{
    check_iterated_preserver_theorem, check_recorded_parent_proof, check_weyl_step_preserver,
    prove_parent_real_rootedness, BigRational, BivarPoly, IteratedPreserverCertificate,
    IteratedPreserverCoordinates, LiftHypothesisStatus, ParentRecurrenceRecord, PreserverOperator,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, Recurrence, RecurrenceTerm,
    RouteCertificate, RowInterlacing, VeroneseParentProof, VeroneseParentProofRoute,
    VeroneseParentProofStatus, WeylStepPreserverVerdict,
};

mod common;
use common::q;

/// The iterated-preserver certificate of a report, for tampering.
fn iterated_preserver_mut(
    report: &mut RealRootednessTheoremReport,
) -> &mut IteratedPreserverCertificate {
    match report.certificate.as_mut() {
        Some(RouteCertificate::IteratedPreserver(certificate)) => certificate,
        _ => panic!("the report has no iterated-preserver certificate"),
    }
}

fn row(values: &[i64]) -> Vec<BigRational> {
    values.iter().map(|&v| q(v)).collect()
}

/// `coeffs[i][j]` = coefficient of `n^i x^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(coeffs.iter().map(|r| row(r)).collect())
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

/// The check with default options, replayed (also through JSON).
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
    report.verify().expect("the report replays");
    let json = serde_json::to_string(&report).unwrap();
    let back: RealRootednessTheoremReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
    report
}

fn region_text(report: &RealRootednessTheoremReport) -> Option<String> {
    report
        .conclusion
        .as_ref()
        .and_then(|c| c.root_region.as_ref())
        .map(ToString::to_string)
}

fn parent_proof(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
) -> VeroneseParentProof {
    prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule,
            seeds,
            first_index,
        },
        first_index,
        &[],
        false,
    )
    .unwrap()
}

fn is_weyl_step(report: &RealRootednessTheoremReport) -> bool {
    matches!(
        report.iterated_preserver_certificate().unwrap().preserver,
        Some(PreserverOperator::WeylOperator { .. })
    )
}

/// Brändén's Corollary 3.5 rows of the catalog (single lag one or two):
/// proved with every root in (-inf, 0], also as parents.
#[test]
fn the_corollary_three_five_catalog_rows_are_proved() {
    // The first group is matched as a Weyl step; the others by the earlier
    // menu operators (a multiplier after phi(D), the Euler step).
    let weyl = [
        "A048854", "A086915", "A143497", "A151775", "A176230", "A241171", "A286724", "A292219",
        "A304330", "A304334", "A397048", "A305540", "A321434",
    ];
    let others = ["A059110", "A156289", "A271703"];
    for id in weyl.iter().chain(others.iter()) {
        let (rule, first_index, seeds) = catalog(id);
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
        assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"), "{id}");
        assert_eq!(is_weyl_step(&report), weyl.contains(id), "{id}");
        if weyl.contains(id) {
            assert!(report.reference.contains("Corollary 3.5"), "{id}");
            for name in [
                "weyl_step_discriminant_nonnegative",
                "weyl_step_interlacing",
                "weyl_step_leading_coefficients_one_signed",
            ] {
                assert!(
                    report.hypotheses.iter().any(|h| h.name == name),
                    "{id}: {name}"
                );
            }
        }
        let proof = parent_proof(&rule, first_index, &seeds);
        assert_eq!(proof.status, VeroneseParentProofStatus::Proved, "{id}");
        assert_eq!(
            proof.route,
            Some(VeroneseParentProofRoute::IteratedPreserver),
            "{id}"
        );
    }
}

#[test]
fn index_dependent_weyl_steps_use_the_sign_checker() {
    // A397048: (x + n) + (2x + n) D + x D^2, discriminant n^2, the root -n
    // of A left of the root -n/2 of B.
    let (rule, first_index, seeds) = catalog("A397048");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    let decided: Vec<_> = report
        .hypotheses
        .iter()
        .filter(|h| h.name.starts_with("weyl_step") && h.sign_decision.is_some())
        .collect();
    assert!(decided.len() >= 3, "{decided:#?}");
    // A tampered sign certificate is refused.
    let mut tampered = report.clone();
    let hypothesis = tampered
        .hypotheses
        .iter_mut()
        .find(|h| h.name == "weyl_step_interlacing")
        .unwrap();
    hypothesis.evidence.push_str(" (edited)");
    assert!(tampered.verify().is_err());
}

#[test]
fn single_lag_rules_run_in_each_residue_class() {
    for id in ["A305540", "A321434"] {
        let (rule, first_index, seeds) = catalog(id);
        assert_eq!(seeds.len(), 2, "{id}");
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
        let certificate = report.iterated_preserver_certificate().unwrap();
        assert_eq!(certificate.lag, 2, "{id}");
        assert!(certificate.coordinates.is_identity());
        assert!(certificate.verify_operator_identity(&rule).is_ok());
        for class in 0..2 {
            assert!(
                report
                    .hypotheses
                    .iter()
                    .any(|h| h.name == format!("degree_step_class_{class}")
                        && h.status == LiftHypothesisStatus::Verified),
                "{id}"
            );
        }
        // No interlacing is claimed for a lag above one.
        assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
        // A tampered lag is refused.
        let mut tampered = report.clone();
        iterated_preserver_mut(&mut tampered).lag = 1;
        assert!(tampered.verify().is_err());
        assert!(tampered
            .iterated_preserver_certificate()
            .unwrap()
            .verify_operator_identity(&rule)
            .is_err());
    }
    // A lag-two multiplier rule: P_n = (1 + 3x + x^2) P_(n-2) after rows that
    // do not follow it.
    let (rule, first_index, seeds) = catalog("A026386");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(report.iterated_preserver_certificate().unwrap().lag, 2);
}

#[test]
fn a_power_gauge_turns_the_rows_into_an_euler_step() {
    // A269945: x + x D + x^2 D^2 = x^(-1)-conjugate of (theta + 1)^2 + x, on
    // the rows P_n / x from P_2 = x on (the start moved one row).
    let (rule, first_index, seeds) = catalog("A269945");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    assert_eq!(report.start_delay, 1);
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert_eq!(
        certificate.coordinates,
        IteratedPreserverCoordinates::PowerGauge { power: 1 }
    );
    assert!(matches!(
        certificate.preserver,
        Some(PreserverOperator::EulerStep { .. })
    ));
    assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"));
    // Recorded in other coordinates, the identity no longer replays.
    let mut tampered = report.clone();
    iterated_preserver_mut(&mut tampered).coordinates = IteratedPreserverCoordinates::Identity;
    assert!(tampered.verify().is_err());
}

#[test]
fn the_reversal_turns_the_rows_into_an_euler_step() {
    // A008957: 1 + x (theta - n)^2, reversed at the nominal degree n:
    // (theta + 1)^2 + x.  A191935: 1 + x (theta - n)(theta - n - 1), reversed:
    // (theta + 1)(theta + 2) + x.
    for id in ["A008957", "A191935"] {
        let (rule, first_index, seeds) = catalog(id);
        let report = check(&rule, first_index, &seeds);
        assert!(report.is_proved(), "{id}: {report:#?}");
        let certificate = report.iterated_preserver_certificate().unwrap();
        assert_eq!(
            certificate.coordinates,
            IteratedPreserverCoordinates::Reversal {
                nominal_previous_degree: vec!["-1".into(), "1".into()]
            },
            "{id}"
        );
        assert!(matches!(
            certificate.preserver,
            Some(PreserverOperator::EulerStep { .. })
        ));
        assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"), "{id}");
        let proof = parent_proof(&rule, first_index, &seeds);
        assert_eq!(proof.status, VeroneseParentProofStatus::Proved, "{id}");
        // A different nominal degree is refused.
        let mut tampered = report.clone();
        iterated_preserver_mut(&mut tampered).coordinates =
            IteratedPreserverCoordinates::Reversal {
                nominal_previous_degree: vec!["0".into(), "1".into()],
            };
        assert!(tampered.verify().is_err());
    }
}

#[test]
fn a_translation_moves_the_euler_step_to_the_origin() {
    // P_n(x) = R_n(x + 1) for the central factorial rows R_n of A036969,
    // R_n = ((1 + x) + 3x D + x^2 D^2) R_(n-1): the step of P has the
    // coefficients c_d(x + 1), (2 + x) + 3(1 + x) D + (1 + x)^2 D^2, which is
    // not a global preserver (discriminant (1 + x)^2 (1 - 4x)).
    let rule = first_order(&[(0, &[&[2, 1]]), (1, &[&[3, 3]]), (2, &[&[1, 2, 1]])]);
    let seeds = vec![row(&[1])];
    let report = check(&rule, 0, &seeds);
    assert!(report.is_proved(), "{report:#?}");
    let certificate = report.iterated_preserver_certificate().unwrap();
    assert_eq!(
        certificate.coordinates,
        IteratedPreserverCoordinates::Translation { shift: "-1".into() }
    );
    // Q_n = P_n(x - 1) has its roots in (-inf, 0], so P_n in (-inf, -1].
    assert_eq!(region_text(&report).as_deref(), Some("(-inf, 0]"));
    let mut tampered = report.clone();
    iterated_preserver_mut(&mut tampered).coordinates =
        IteratedPreserverCoordinates::Translation { shift: "1".into() };
    assert!(tampered.verify().is_err());
    assert!(tampered
        .iterated_preserver_certificate()
        .unwrap()
        .verify_operator_identity(&rule)
        .is_err());
}

#[test]
fn fixed_weyl_steps_are_decided_exactly() {
    // The unsigned Lah-type step (2 + x) + (8 + 8x) D + 16x D^2 (A048854).
    let verdict = check_weyl_step_preserver(&row(&[2, 1]), &row(&[8, 8]), &row(&[0, 16]));
    assert_eq!(verdict, WeylStepPreserverVerdict::GlobalPreserver);
    // The central factorial step of A036969, (1 + x) + 3x D + x^2 D^2: the
    // discriminant x^2 (5 - 4x) is negative for x > 5/4, so it is not a
    // global preserver (it is a half-line Euler step).
    let verdict = check_weyl_step_preserver(&row(&[1, 1]), &row(&[0, 3]), &row(&[0, 0, 1]));
    match verdict {
        WeylStepPreserverVerdict::NotGlobalPreserver {
            point,
            discriminant,
        } => {
            assert!(point > BigRational::new(5.into(), 4.into()));
            assert!(discriminant < q(0));
        }
        other => panic!("{other:?}"),
    }
    // A nonnegative discriminant, but the roots of A and B do not
    // alternate: x + (x - 1) D is a preserver, x + (x + 1) D has A's root 0
    // right of B's root -1.
    assert_eq!(
        check_weyl_step_preserver(&row(&[0, 1]), &row(&[-1, 1]), &[]),
        WeylStepPreserverVerdict::GlobalPreserver
    );
    assert!(matches!(
        check_weyl_step_preserver(&row(&[0, 1]), &row(&[1, 1]), &[]),
        WeylStepPreserverVerdict::Undecided { .. }
    ));
    // Leading coefficients that change sign with k: 1 - x D maps x^k to
    // (1 - k) x^k.
    assert!(matches!(
        check_weyl_step_preserver(&row(&[1]), &row(&[0, -1]), &[]),
        WeylStepPreserverVerdict::Undecided { .. }
    ));
}

#[test]
fn the_half_line_euler_step_still_wins_where_the_weyl_step_fails() {
    // A036969 is refused as a global preserver (above) and proved by the
    // half-line Euler step, tried before the Weyl step.
    let (rule, first_index, seeds) = catalog("A036969");
    let report = check(&rule, first_index, &seeds);
    assert!(report.is_proved());
    assert!(matches!(
        report.iterated_preserver_certificate().unwrap().preserver,
        Some(PreserverOperator::EulerStep { .. })
    ));
    // A step with no other menu shape and a negative discriminant at x = 0:
    // (1 + x) + D + (1 + x^2) D^2.
    let rule = first_order(&[(0, &[&[1, 1]]), (1, &[&[1]]), (2, &[&[1, 0, 1]])]);
    let report = check(&rule, 0, &[row(&[1])]);
    assert!(!report.is_proved(), "{report:#?}");
    let discriminant = report
        .hypotheses
        .iter()
        .find(|h| h.name == "weyl_step_discriminant_nonnegative")
        .expect("the Weyl-step attempt is reported");
    assert!(matches!(
        discriminant.status,
        LiftHypothesisStatus::Failed { .. }
    ));
    assert!(discriminant
        .evidence
        .contains("not a real-rootedness preserver"));
}

#[test]
fn a_step_free_of_n_gives_the_corollary_two_five_interlacing() {
    for id in ["A048854", "A143497", "A151775"] {
        let (rule, first_index, seeds) = catalog(id);
        let report = check(&rule, first_index, &seeds);
        let conclusion = report.conclusion.as_ref().unwrap();
        assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak), "{id}");
        assert!(report.reference.contains("Corollary 2.5"), "{id}");
        // Evidence: consecutive generated rows interlace.
        let from = conclusion.interlacing_from.unwrap();
        let rows = rule
            .generate_rows_rational(&seeds, first_index, 14)
            .unwrap();
        let integer = |p: &[BigRational]| -> Vec<BigInt> {
            let mut p = p.to_vec();
            while p.last().is_some_and(num_traits::Zero::is_zero) {
                p.pop();
            }
            let lcm = p.iter().fold(BigInt::from(1), |acc, c| {
                num_integer::Integer::lcm(&acc, c.denom())
            });
            p.iter().map(|c| c.numer() * (&lcm / c.denom())).collect()
        };
        for m in from.max(first_index + 1)..first_index + rows.len() {
            let previous = integer(&rows[m - 1 - first_index]);
            let current = integer(&rows[m - first_index]);
            assert_eq!(
                polytool::check_weak_interlacing_bigint_coeffs(&previous, &current),
                Some(true),
                "{id}: P_{} and P_{m}",
                m - 1
            );
        }
    }
    // An index-dependent step claims no interlacing.
    let (rule, first_index, seeds) = catalog("A397048");
    let report = check(&rule, first_index, &seeds);
    assert!(report.conclusion.as_ref().unwrap().interlacing.is_none());
}

#[test]
fn a_tampered_weyl_step_is_refused() {
    let (rule, first_index, seeds) = catalog("A143497");
    let report = check(&rule, first_index, &seeds);
    let mut tampered = report.clone();
    let certificate = iterated_preserver_mut(&mut tampered);
    let Some(PreserverOperator::WeylOperator { terms, .. }) = certificate.preserver.as_mut() else {
        panic!("a Weyl step");
    };
    terms[0].coeff.coeffs[0][0] = "5".into();
    assert!(tampered
        .iterated_preserver_certificate()
        .unwrap()
        .verify_operator_identity(&rule)
        .is_err());
    assert!(tampered.verify().is_err());
    // A proved report without its certificate is refused.
    let mut stripped = report.clone();
    stripped.certificate = None;
    assert!(stripped.verify().is_err());
}

#[test]
fn weaker_recorded_parent_verdicts_keep_replaying() {
    // A304334 was open for the parent provers before the Weyl step: a
    // recorded `not_proved` verdict is still supported, a recorded proof by
    // the Weyl step replays, and a proof claimed by another route is refused.
    let (rule, first_index, seeds) = catalog("A304334");
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let current = prove_parent_real_rootedness(&record, first_index, &[], false).unwrap();
    assert_eq!(current.status, VeroneseParentProofStatus::Proved);
    check_recorded_parent_proof(&record, first_index, &[], false, &current).unwrap();
    let mut open = current.clone();
    open.status = VeroneseParentProofStatus::NotProved;
    open.route = None;
    open.interlacing = None;
    open.statement = "open".into();
    check_recorded_parent_proof(&record, first_index, &[], false, &open).unwrap();
    let mut other = current.clone();
    other.route = Some(VeroneseParentProofRoute::LiuWang);
    assert!(check_recorded_parent_proof(&record, first_index, &[], false, &other).is_err());
    // A026386 was proved by the product induction; the iterated route now
    // settles it last, and a record by the product route replays by it.
    let (rule, first_index, seeds) = catalog("A026386");
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let current = prove_parent_real_rootedness(&record, first_index, &[], false).unwrap();
    assert_eq!(
        current.route,
        Some(VeroneseParentProofRoute::IteratedPreserver)
    );
    let mut product = current.clone();
    product.route = Some(VeroneseParentProofRoute::ProductRecurrence);
    product.interlacing = None;
    check_recorded_parent_proof(&record, first_index, &[], false, &product).unwrap();
}
