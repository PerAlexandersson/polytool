//! Direct checks of the owned sparse representation and generic kernels.

use super::*;
use num_rational::BigRational;

type P = Polynomial<RationalField>;
fn term(e: &[u32], c: i64) -> P {
    P::monomial(e.len(), e.to_vec(), BigRational::from_integer(c.into()))
}

#[test]
fn duplicate_generators_and_equal_leading_terms_preserve_the_ideal() {
    let x = term(&[1, 0], 1);
    let y = term(&[0, 1], 1);
    for order in MonomialOrder::STANDARD_ORDERS {
        for batch_reduce_same_lcm_degree in [false, true] {
            for (input, expected) in [
                (vec![x.clone(), x.clone()], vec![x.clone()]),
                (
                    vec![x.clone(), x.clone() + y.clone()],
                    vec![x.clone(), y.clone()],
                ),
            ] {
                let result = reduced_groebner_basis_with_options(
                    &input,
                    order,
                    GroebnerOptions {
                        batch_reduce_same_lcm_degree,
                        ..GroebnerOptions::default()
                    },
                );
                assert_eq!(result, expected);
                assert!(input
                    .iter()
                    .all(|p| normal_form(p, &result, order).is_zero()));
            }
        }
    }
}

#[test]
fn ordered_division_reconstructs_the_dividend_and_cached_normal_forms() {
    let p = term(&[2, 0], 1) + term(&[0, 1], 2) + term(&[0, 0], 1);
    let divisors = vec![term(&[1, 0], 1) + term(&[0, 0], -1), term(&[0, 1], 1)];
    for order in MonomialOrder::STANDARD_ORDERS {
        let division = divide_by_polynomials(&p, &divisors, order);
        // Both divisors are linear, so the quotients have these known values.
        assert_eq!(
            division.quotients,
            vec![term(&[1, 0], 1) + term(&[0, 0], 1), term(&[0, 0], 2)]
        );
        assert_eq!(division.remainder, term(&[0, 0], 2));
        let leading: Vec<_> = divisors
            .iter()
            .map(|g| leading_term(g, order).unwrap())
            .collect();
        assert_eq!(
            normal_form_with_leading_terms(&p, &divisors, &leading, order),
            division.remainder
        );
        let inputs = vec![p.clone(), p.clone() + term(&[0, 0], 1)];
        assert_eq!(
            matrix_normal_forms_with_leading_terms(&inputs, &divisors, &leading, order),
            inputs
                .iter()
                .map(|p| normal_form(p, &divisors, order))
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn sparse_arithmetic_cancels_terms_and_preserves_zero_ambient_dimensions() {
    let p = term(&[1, 0, 0], 1);
    assert_eq!(p.clone() - p, P::zero(3));
    let basis = GroebnerBasis::new(vec![P::zero(3)], MonomialOrder::Lex);
    assert_eq!(basis.num_vars, 3);
    assert_eq!(basis.normal_form(&term(&[1, 0, 0], 2)), term(&[1, 0, 0], 2));
    assert_eq!(
        reduced_groebner_basis(&[term(&[], 2)], MonomialOrder::Lex),
        vec![term(&[], 1)]
    );
}

#[test]
#[should_panic(expected = "all exponent vectors must have length num_vars")]
fn malformed_exponent_lengths_are_rejected() {
    let _ = P::monomial(2, vec![1], BigRational::from_integer(1.into()));
}

#[test]
#[should_panic(expected = "division by a zero polynomial")]
fn zero_divisors_are_rejected() {
    let _ = normal_form(&term(&[1], 1), &[P::zero(1)], MonomialOrder::Lex);
}

#[test]
#[should_panic(expected = "all generators must have the same number of variables")]
fn inconsistent_generator_rings_are_rejected() {
    let _ = reduced_groebner_basis(&[term(&[1], 1), term(&[1, 0], 1)], MonomialOrder::Lex);
}

#[test]
#[should_panic(expected = "monomial exponent overflow during multiplication")]
fn multiplying_exponents_checks_overflow() {
    let _ = multiply_by_monomial(
        &term(&[u32::MAX], 1),
        &[1],
        BigRational::from_integer(1.into()),
    );
}
