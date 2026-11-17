//! Boundary checks for the ordered, budgeted rational policy.
use super::*;
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

type P = Polynomial<RationalField>;
fn polynomial(variables: usize, terms: &[(Vec<u32>, i64)]) -> P {
    P::from_terms(
        variables,
        terms
            .iter()
            .map(|(m, c)| (m.clone(), BigRational::from_integer((*c).into())))
            .collect(),
    )
}
fn budget(steps: usize) -> GroebnerBudget {
    GroebnerBudget {
        remaining_steps: steps,
        max_coefficient_bits: 4096,
        max_basis_elements: 400,
    }
}

#[test]
fn ordered_initialization_and_shared_last_step() {
    let x = polynomial(1, &[(vec![1], 1)]);
    let input = [x.clone(), x.clone()];
    let mut limits = budget(0);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&input, &mut limits),
        Err(GroebnerError::StepBudgetExhausted)
    );
    let mut limits = budget(1);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&input, &mut limits),
        Ok(vec![x])
    );
    assert_eq!(limits.remaining_steps, 0);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&input, &mut limits),
        Err(GroebnerError::StepBudgetExhausted)
    );
}

#[test]
fn tails_are_reduced_in_insertion_order() {
    let x_y = polynomial(2, &[(vec![1, 0], 1), (vec![0, 1], 1)]);
    let y = polynomial(2, &[(vec![0, 1], 1)]);
    let x = polynomial(2, &[(vec![1, 0], 1)]);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[x_y.clone(), y.clone()], &mut budget(0)),
        Err(GroebnerError::StepBudgetExhausted)
    );
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[x_y, y.clone()], &mut budget(1)),
        Ok(vec![x, y])
    );
}

#[test]
fn bit_check_follows_subtraction_and_step_debit() {
    let x = polynomial(2, &[(vec![1, 0], 1)]);
    let mut p = x.clone();
    p.add_term(
        vec![0, 1],
        &BigRational::from_integer(BigInt::one() << 4096),
    );
    for (steps, error) in [
        (0, GroebnerError::StepBudgetExhausted),
        (1, GroebnerError::CoefficientBitBudgetExceeded),
    ] {
        let mut limits = budget(steps);
        assert_eq!(
            reduced_lexicographic_groebner_basis_with_budget(&[x.clone(), p.clone()], &mut limits),
            Err(error)
        );
        assert_eq!(limits.remaining_steps, 0);
    }
    // Merely constructing or scaling the same large coefficient is not checked.
    assert!(reduced_lexicographic_groebner_basis_with_budget(&[p], &mut budget(0)).is_ok());
}

#[test]
fn basis_limit_is_only_checked_on_spair_append() {
    let input = [
        polynomial(2, &[(vec![1, 1], 1), (vec![0, 0], -1)]),
        polynomial(2, &[(vec![0, 2], 1), (vec![0, 0], -1)]),
    ];
    let mut limits = budget(100);
    limits.max_basis_elements = 2;
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&input, &mut limits),
        Err(GroebnerError::BasisSizeBudgetExceeded)
    );
    let mut limits = budget(0);
    limits.max_basis_elements = 0;
    let x = polynomial(2, &[(vec![1, 0], 1)]);
    let y = polynomial(2, &[(vec![0, 1], 1)]);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[x.clone(), y.clone()], &mut limits),
        Ok(vec![x, y])
    );
}

#[test]
fn zero_constant_and_invalid_rings() {
    let zero = P::zero(0);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[zero], &mut budget(0)),
        Ok(vec![])
    );
    let one = polynomial(0, &[(vec![], 2)]);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[one], &mut budget(0)),
        Ok(vec![P::monomial(0, vec![], BigRational::one())])
    );
    let mut limits = budget(3);
    assert_eq!(
        reduced_lexicographic_groebner_basis_with_budget(&[P::zero(1), P::zero(2)], &mut limits),
        Err(GroebnerError::IncompatibleVariableCounts)
    );
    assert_eq!(limits.remaining_steps, 3);
    let mut p = P::zero(1);
    p.add_term(vec![1], &BigRational::one());
    p.add_term(vec![1], &-BigRational::one());
    assert!(p.is_zero());
    p.add_term(vec![1], &BigRational::zero());
    assert!(p.is_zero());
}
