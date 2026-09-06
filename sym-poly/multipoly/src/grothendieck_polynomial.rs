//! Connective-K Grothendieck polynomials and their Lascoux expansions.
//!
//! The deformation parameter is supplied as an element `beta` of the
//! coefficient ring. Use `UnivariatePolynomial<i64>::variable()` to retain a
//! formal parameter, `0` for the Schubert specialization, and `-1` for the
//! usual Grothendieck specialization.

use std::collections::BTreeMap;

use combinatoric_core::{
    compose_permutations, inverse_permutation, longest_permutation, reduced_word,
};
use sym_poly_core::Ring;

use crate::key_polynomial::{dominant_rearrangement, sorting_reduced_word};
use crate::multipoly::MultiPoly;
use crate::operators::{demazure_lascoux_partial_word, demazure_lascoux_pi_word};

/// A Lascoux-basis expansion indexed by weak compositions.
pub type LascouxExpansion<C> = BTreeMap<Vec<u32>, C>;

/// Compute the connective-K Grothendieck polynomial `G_w^(beta)`.
///
/// The permutation uses one-indexed values. The recursion starts with
/// `x^(n-1,n-2,...,0)` at the longest permutation and applies connective-K
/// divided differences toward `w`.
pub fn beta_grothendieck_polynomial<C: Ring>(perm: &[usize], beta: &C) -> MultiPoly<C> {
    let n = perm.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let w0 = longest_permutation(n);
    let inverse = inverse_permutation(perm);
    let quotient = compose_permutations(&inverse, &w0);
    let word = reduced_word(&quotient);
    let rho = (0..n).map(|i| (n - 1 - i) as u32).collect();
    let top = MultiPoly::x_power(n, rho);
    demazure_lascoux_partial_word(&top, &word, beta)
}

/// Compute the ordinary Grothendieck polynomial, the `beta = -1`
/// specialization of [`beta_grothendieck_polynomial`].
pub fn grothendieck_polynomial<C: Ring>(perm: &[usize]) -> MultiPoly<C> {
    beta_grothendieck_polynomial(perm, &C::minus_one())
}

/// Compute the Lascoux polynomial `L_alpha^(beta)` by operators.
///
/// This is independent of the existing K-Kohnert implementation and is
/// useful both as a faster constructor and as a cross-check of conventions.
pub fn lascoux_polynomial_by_operators<C: Ring>(alpha: &[u32], beta: &C) -> MultiPoly<C> {
    let n = alpha.len();
    if n == 0 {
        return MultiPoly::constant(0, C::one());
    }

    let dominant = dominant_rearrangement(alpha);
    let word = sorting_reduced_word(alpha);
    let top = MultiPoly::x_power(n, dominant);
    let operator_word = word.iter().rev().copied().collect::<Vec<_>>();
    demazure_lascoux_pi_word(&top, &operator_word, beta)
}

/// Expand `G_w^(beta)` in the Lascoux basis.
///
/// The result maps a weak composition `alpha` to the coefficient of
/// `L_alpha^(beta)`. The reduction uses the unitriangular leading monomials of
/// Lascoux polynomials.
pub fn beta_grothendieck_to_lascoux<C: Ring>(perm: &[usize], beta: &C) -> LascouxExpansion<C> {
    let grothendieck = beta_grothendieck_polynomial(perm, beta);
    polynomial_to_lascoux(&grothendieck, beta)
}

/// Expand the ordinary Grothendieck polynomial in ordinary Lascoux
/// polynomials at `beta = -1`.
pub fn grothendieck_to_lascoux<C: Ring>(perm: &[usize]) -> LascouxExpansion<C> {
    beta_grothendieck_to_lascoux(perm, &C::minus_one())
}

/// Expand a polynomial in the Lascoux basis for the supplied deformation
/// parameter.
pub fn polynomial_to_lascoux<C: Ring>(polynomial: &MultiPoly<C>, beta: &C) -> LascouxExpansion<C> {
    let mut residual = polynomial.clone();
    let mut expansion = BTreeMap::new();

    while let Some((alpha, coefficient)) = residual
        .terms()
        .iter()
        .next()
        .map(|(alpha, coefficient)| (alpha.clone(), coefficient.clone()))
    {
        let lascoux = lascoux_polynomial_by_operators(&alpha, beta);
        assert_eq!(
            lascoux.coefficient(&alpha),
            C::one(),
            "Lascoux polynomial does not have its expected unit leading term"
        );

        residual = residual - lascoux.scale(&coefficient);
        if let Some(next_alpha) = residual.terms().keys().next() {
            assert!(
                next_alpha > &alpha,
                "Lascoux reduction did not advance the leading monomial"
            );
        }
        let entry = expansion.entry(alpha).or_insert_with(C::zero);
        *entry = entry.clone() + coefficient;
    }

    expansion.retain(|_, coefficient| !coefficient.is_zero());
    expansion
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{lascoux_polynomial, schubert_polynomial};
    use sym_poly_core::UnivariatePolynomial;

    type BetaPolynomial = UnivariatePolynomial<i64>;

    fn beta() -> BetaPolynomial {
        BetaPolynomial::variable()
    }

    fn reconstruct(expansion: &LascouxExpansion<BetaPolynomial>) -> MultiPoly<BetaPolynomial> {
        let num_vars = expansion.keys().next().map(Vec::len).unwrap_or(0);
        expansion
            .iter()
            .fold(MultiPoly::zero(num_vars), |sum, (alpha, coefficient)| {
                sum + lascoux_polynomial_by_operators(alpha, &beta()).scale(coefficient)
            })
    }

    #[test]
    fn beta_zero_is_schubert() {
        for perm in [
            vec![1, 2, 3],
            vec![2, 1, 3],
            vec![1, 3, 2],
            vec![2, 3, 1],
            vec![3, 1, 2],
            vec![2, 1, 4, 3],
        ] {
            assert_eq!(
                beta_grothendieck_polynomial::<i64>(&perm, &0),
                schubert_polynomial::<i64>(&perm),
                "beta=0 mismatch for {perm:?}"
            );
        }
    }

    #[test]
    fn operator_lascoux_matches_k_kohnert() {
        for alpha in [vec![0, 2, 1], vec![1, 0, 2], vec![0, 1, 2]] {
            let by_operators = lascoux_polynomial_by_operators(&alpha, &beta());
            let by_kohnert = lascoux_polynomial(&alpha, &beta(), 10_000).unwrap();
            assert_eq!(by_operators, by_kohnert, "mismatch for {alpha:?}");
        }
    }

    #[test]
    fn grothendieck_lascoux_expansion_reconstructs() {
        for perm in [vec![2, 1, 4, 3], vec![2, 1, 3, 5, 4]] {
            let expansion = beta_grothendieck_to_lascoux(&perm, &beta());
            assert_eq!(
                reconstruct(&expansion),
                beta_grothendieck_polynomial(&perm, &beta()),
                "reconstruction failed for {perm:?}"
            );
            assert!(expansion
                .values()
                .all(|coefficient| { coefficient.coeffs().iter().all(|&integer| integer >= 0) }));
        }
    }

    #[test]
    fn grothendieck_2143_lascoux_expansion() {
        let expansion = beta_grothendieck_to_lascoux(&[2, 1, 4, 3], &beta());
        let expected = BTreeMap::from([
            (vec![1, 0, 1, 0], BetaPolynomial::constant(1)),
            (vec![2, 0, 0, 0], BetaPolynomial::constant(1)),
            (vec![2, 0, 1, 0], BetaPolynomial::variable()),
        ]);
        assert_eq!(expansion, expected);
    }
}
