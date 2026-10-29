//! Unit tests of the v2 kernel: pseudo-division, gauges and the parent
//! provers (the transport map has its own tests in `section_operators`).

use super::super::implicit::algebra::{Biv, DOp};
use super::parent_proof::{periodic_product_parent, prove_parent, ParentRecord};
use super::section_certify::{
    division_identity_holds, gauge_parent_relation, pseudo_divide_row, residue_row, x_free_leader,
};
use super::section_operators::transport_to_sections;
use super::*;
use num_bigint::BigInt;

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|s| s.to_string()).collect()
}

/// `P_(n+1) = (1 + x) P_n` at base `n`: `S - (1 + x)`.
fn binomial_relation() -> DOp {
    let mut op = DOp::single(1, 0, Biv::one());
    op.add_term((0, 0), Biv::from_y(&[q(-1), q(-1)]));
    op
}

#[test]
fn pseudo_division_of_the_worked_example() {
    // M = S^2 - 2 S - (y - 1) annihilates the even sections of (1 + x)^n.
    let mut input = DOp::single(2, 0, Biv::one());
    input.add_term((1, 0), Biv::constant(q(-2)));
    input.add_term((0, 0), Biv::from_y(&[q(1), q(-1)]));
    let parent = binomial_relation();
    let (k, leader) = x_free_leader(&parent).unwrap();
    let system = transport_to_sections(&parent, 2, 10_000).unwrap();
    let row = residue_row(&input, 2, 0);
    let division = pseudo_divide_row(&row, &system, k, &leader, 10_000).unwrap();
    assert!(division.remainder_is_zero());
    assert!(division_identity_holds(&row, &system, &division));
    assert_eq!(division.multiplier, Biv::one());

    // S^2 - 2 S - y does not annihilate them: a nonzero remainder.
    let mut wrong = DOp::single(2, 0, Biv::one());
    wrong.add_term((1, 0), Biv::constant(q(-2)));
    wrong.add_term((0, 0), Biv::from_y(&[q(0), q(-1)]));
    let row = residue_row(&wrong, 2, 0);
    let division = pseudo_divide_row(&row, &system, k, &leader, 10_000).unwrap();
    assert!(!division.remainder_is_zero());
    assert!(division_identity_holds(&row, &system, &division));
    // The term budget is enforced.
    assert!(matches!(
        pseudo_divide_row(&row, &system, k, &leader, 2),
        Err(RecurrenceError::BudgetExhausted(_))
    ));
}

#[test]
fn gauges_build_the_same_relation_in_both_forms() {
    // lambda(n) = n + 1: ratio (n + 1)/n from n0 = 0 ... but the ratio
    // denominator n vanishes at n = 0 only, which is not > n0.
    let closed = VeroneseSectionGauge::Closed {
        numerator: strings(&["1", "1"]),
        denominator: strings(&["1"]),
    };
    let ratio = VeroneseSectionGauge::Ratio {
        numerator: strings(&["1", "1"]),
        denominator: strings(&["0", "1"]),
        initial_value: "1".into(),
    };
    let from_closed = closed.ratio_form(0).unwrap();
    let from_ratio = ratio.ratio_form(0).unwrap();
    assert_eq!(from_closed.values(0, 5), from_ratio.values(0, 5));
    assert_eq!(from_ratio.values(0, 3), vec![q(1), q(2), q(3)]);
    let parent = binomial_relation();
    let left = gauge_parent_relation(&parent, 1, &from_closed).unwrap();
    let right = gauge_parent_relation(&parent, 1, &from_ratio).unwrap();
    // Both are multiples of (n + 1) P~_(n+1) = (n + 2)(1 + x) P~_n up to a
    // factor in n; compare their canonical forms.
    assert!(ImplicitRelation::from_op(&left).same_relation(&ImplicitRelation::from_op(&right)));

    // Vanishing gauges are rejected.
    let vanishing = VeroneseSectionGauge::Closed {
        numerator: strings(&["-2", "1"]),
        denominator: strings(&["1"]),
    };
    assert!(matches!(
        vanishing.ratio_form(0),
        Err(RecurrenceError::InvalidInput(_))
    ));
    assert!(vanishing.ratio_form(3).is_ok());
    let zero_initial = VeroneseSectionGauge::Ratio {
        numerator: strings(&["1"]),
        denominator: strings(&["1"]),
        initial_value: "0".into(),
    };
    assert!(zero_initial.ratio_form(0).is_err());
}

#[test]
fn periodic_product_route() {
    let rule = |c: &[i64]| Recurrence {
        terms: vec![RecurrenceTerm::new(
            2,
            0,
            BivarPoly::from_coefficients_normalized(vec![c.iter().map(|&v| q(v)).collect()]),
        )],
        denominator: Some(BivarPoly::from_coefficients_normalized(vec![vec![q(2)]])),
        inhomogeneous: None,
    };
    let seeds = vec![vec![q(1)], vec![q(1), q(3)]];
    assert!(periodic_product_parent(&rule(&[2, 2]), &seeds).is_ok());
    assert!(periodic_product_parent(&rule(&[2, -2]), &seeds).is_err());
    assert!(periodic_product_parent(&rule(&[1, 1, 1]), &seeds).is_err());
    let bad_seeds = vec![vec![q(1)], vec![q(1), q(1), q(1)]];
    assert!(periodic_product_parent(&rule(&[2, 2]), &bad_seeds).is_err());
}

/// The Eulerian parent (A008292) is proved by the root-domain route: the
/// f-polynomial map with the fitted schedule `d_n = n - 1` turns it into
/// `Q_(n+1) = (x (1 + x) Q_n)'`, closed by Rolle on `[-1, 0]`.
#[test]
fn eulerian_parent_is_proved_by_the_root_domain_route() {
    let (rule, seeds, first_row) = crate::oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    let first_index = usize::try_from(first_row).unwrap();
    let record = ParentRecord {
        rule: &rule,
        seeds: &seeds,
        first_index,
    };
    let proof = prove_parent(&record, first_index, &[], true).unwrap();
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Proved,
        "{proof:#?}"
    );
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::RootDomain));
    let success = proof.attempts.iter().find(|a| a.success).unwrap();
    assert_eq!(success.variant, "f_polynomial");
    // Without the root-domain route the Liu--Wang route proves it.
    let proof = prove_parent(&record, first_index, &[], false).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::LiuWang));
}

/// The Eulerian rule `A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1)`
/// from `A_0 = 1`: every route fails at the first index (`A_1 = A_0`), and
/// moving the start one row later proves the parent.
#[test]
fn parent_routes_move_the_induction_start() {
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(-1)], vec![q(0), q(1)]]),
            ),
            RecurrenceTerm::new(
                1,
                1,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1), q(-1)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = vec![vec![q(1)]];
    let record = ParentRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: 0,
    };
    let proof = prove_parent(&record, 0, &[], true).unwrap();
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Proved,
        "{proof:#?}"
    );
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::RootDomain));
    let success = proof.attempts.iter().find(|a| a.success).unwrap();
    assert_eq!(
        success.variant,
        "f_polynomial, induction start moved 1 row(s)"
    );
    assert!(success.detail.contains("checked exactly"));
    assert!(proof.statement.contains("n >= 0"));
    // Without the root-domain route, Liu--Wang moves its start itself.
    let proof = prove_parent(&record, 0, &[], false).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::LiuWang));
    let success = proof.attempts.iter().find(|a| a.success).unwrap();
    assert!(success.detail.starts_with("induction start moved 1 row(s)"));
}

fn catalog_record(id: &str) -> (Recurrence, Vec<Vec<BigRational>>, usize) {
    let (rule, seeds, first_row) = crate::oeis::by_id(id)
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    (rule, seeds, usize::try_from(first_row).unwrap())
}

/// The open parents of the Veronese v2 lifts: the Narayana parent (A001263,
/// order two) by the three-term route, the type-B Eulerian parent (A060187)
/// by the Liu--Wang route; the root-domain menu proves neither.
#[test]
fn narayana_and_type_b_eulerian_parents_are_proved() {
    for (id, route) in [
        ("A001263", VeroneseParentProofRoute::ThreeTermRecurrence),
        ("A060187", VeroneseParentProofRoute::LiuWang),
    ] {
        let (rule, seeds, first_index) = catalog_record(id);
        let record = ParentRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof = prove_parent(&record, first_index, &[], true).unwrap();
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::Proved,
            "{id}: {proof:#?}"
        );
        assert_eq!(proof.route, Some(route), "{id}");
        assert!(proof
            .attempts
            .iter()
            .filter(|a| a.route == VeroneseParentProofRoute::RootDomain)
            .all(|a| !a.success));
    }
}

/// The Laguerre family `L_n(-x)` (roots in `(-inf, 0)`) is proved by the
/// three-term route; the Hermite family (roots of both signs) is refuted by
/// the row check before any route runs.
#[test]
fn three_term_route_proves_laguerre_and_refutes_hermite() {
    // Laguerre L_n(-x): n P_n = (2n - 1 + x) P_(n-1) - (n - 1) P_(n-2).
    let laguerre = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(-1), q(1)], vec![q(2)]]),
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(1)], vec![q(-1)]]),
            ),
        ],
        denominator: Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(0)],
            vec![q(1)],
        ])),
        inhomogeneous: None,
    };
    let seeds = vec![vec![q(1)], vec![q(1), q(1)]];
    let record = ParentRecord {
        rule: &laguerre,
        seeds: &seeds,
        first_index: 0,
    };
    let proof = prove_parent(&record, 0, &[], false).unwrap();
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Proved,
        "{proof:#?}"
    );
    assert_eq!(
        proof.route,
        Some(VeroneseParentProofRoute::ThreeTermRecurrence)
    );

    // Hermite: H_n = 2 x H_(n-1) - 2 (n - 1) H_(n-2); H_2 = 4 x^2 - 2.
    let hermite = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(2)]]),
            ),
            RecurrenceTerm::new(
                2,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![q(2)], vec![q(-2)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = vec![vec![q(1)], vec![q(0), q(2)]];
    let record = ParentRecord {
        rule: &hermite,
        seeds: &seeds,
        first_index: 0,
    };
    let proof = prove_parent(&record, 0, &[], false).unwrap();
    assert_eq!(proof.status, VeroneseParentProofStatus::Refuted);
    assert_eq!(proof.failing_row_index, Some(2));
}

/// `f(n, t) P_n = sum c_k(n, t) D^d P_(n-k)` from `(k, d, rows)`, `rows[i][j]`
/// the coefficient of `n^i t^j`.
fn summary_rule(terms: &[(usize, usize, &[&[i64]])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, order, coefficients)| {
                RecurrenceTerm::new(
                    *offset,
                    *order,
                    BivarPoly::from_coefficients_normalized(
                        coefficients
                            .iter()
                            .map(|r| r.iter().map(|&c| q(c)).collect())
                            .collect(),
                    ),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

#[test]
fn lift_summaries_of_proved_and_refuted_certifications() {
    // The even sections of the Eulerian polynomials (A008292), proved through
    // the root-domain parent route.
    let input = summary_rule(&[
        (1, 0, &[&[-1], &[1]]),
        (1, 1, &[&[0, 4]]),
        (2, 0, &[&[2, 4], &[-1, -4], &[0, 1]]),
        (2, 1, &[&[0, 0, 12], &[0, -4, -4]]),
        (2, 2, &[&[0, 0, -4, 4]]),
        (3, 0, &[&[0, 12], &[0, -19], &[0, 8], &[0, -1]]),
        (3, 1, &[&[0, -6, 18], &[0, 6, -22], &[0, 0, 4]]),
        (3, 2, &[&[0, 0, -4, 4], &[0, 0, 4, -4]]),
    ]);
    let (parent_rule, parent_seeds, first_row) = crate::oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    assert_eq!(first_row, 1);
    let even: Vec<Vec<BigRational>> = parent_rule
        .generate_rows_rational(&parent_seeds, 1, 3)
        .unwrap()
        .iter()
        .map(|row| row.iter().step_by(2).cloned().collect())
        .collect();
    let candidate = VeroneseSectionCandidate {
        modulus: 2,
        residue: 0,
        parent: RecurrenceJson::from_recurrence_rational(&parent_rule, 1, &parent_seeds, None),
        index_shift: 0,
        gauge: VeroneseSectionGauge::None,
    };
    let options = VeroneseSectionCertifyOptions::default();
    let certification =
        certify_veronese_section_lift(&input, &even, 1, &candidate, &options).unwrap();
    assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
    let summary = certification.lift_summary.clone().expect("summary present");
    assert_eq!(summary, certification.recompute_lift_summary());
    assert_eq!(summary.verdict, LiftVerdict::Proved);
    assert!(summary.is_consistent());
    assert_eq!(summary.input_first_index, 1);
    assert_eq!(
        summary.identity,
        LiftIdentityStatus::Certified {
            method: "forward_division".into()
        }
    );
    assert!(matches!(
        &summary.parent_proof,
        LiftParentProof::Proved { route, .. } if route == "root_domain"
    ));
    assert_eq!(summary.parent.gauge, Some(VeroneseSectionGauge::None));
    assert_eq!(summary.hypotheses.len(), 1);
    assert_eq!(summary.hypotheses[0].status, LiftHypothesisStatus::Verified);
    certification.verify().unwrap();

    // A tampered summary is refused by the replay.
    let mut tampered = certification.clone();
    if let Some(summary) = tampered.lift_summary.as_mut() {
        summary.verdict = LiftVerdict::Open;
    }
    assert!(tampered.verify().is_err());
    // A record without the field (written before it existed) still replays.
    let mut legacy = serde_json::to_value(&certification).unwrap();
    legacy.as_object_mut().unwrap().remove("lift_summary");
    let legacy: VeroneseSectionCertification = serde_json::from_value(legacy).unwrap();
    assert!(legacy.lift_summary.is_none());
    legacy.verify().unwrap();

    // A wrong seed row is refuted at index 3.
    let mut wrong = even.clone();
    wrong[2] = vec![q(1), q(2)];
    let refuted = certify_veronese_section_lift(&input, &wrong, 1, &candidate, &options).unwrap();
    let summary = refuted.lift_summary.clone().unwrap();
    assert_eq!(summary.verdict, LiftVerdict::RouteFailed);
    assert!(matches!(
        summary.identity,
        LiftIdentityStatus::Refuted { index: Some(3), .. }
    ));
    refuted.verify().unwrap();
}
