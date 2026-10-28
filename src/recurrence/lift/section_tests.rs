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
    // Without the root-domain route it stays open.
    let open = prove_parent(&record, first_index, &[], false).unwrap();
    assert_eq!(open.status, VeroneseParentProofStatus::NotProved);
}
