//! Divided difference operators for multivariate polynomials.
//!
//! - `partial_i(f)` = (f - s_i(f)) / (x_i - x_{i+1})  (simple divided difference)
//! - `pi_i(f)` = partial_i(x_i * f)  (isobaric/Demazure operator)
//! - `theta_i(f)` = pi_i(f) - f      (Demazure atom operator)
//! - `demazure_lascoux_partial_i(f)` = partial_i((1 + beta*x_{i+1}) * f)
//! - `demazure_lascoux_pi_i(f)` = partial_i(x_i * (1 + beta*x_{i+1}) * f)
//!
//! Both always produce polynomials when applied to polynomials.

use crate::multipoly::MultiPoly;
use sym_poly_core::Ring;

/// Simple divided difference ∂_i: (f - s_i(f)) / (x_i - x_{i+1}).
///
/// The i-th simple divided difference (0-indexed). Satisfies:
/// - ∂_i^2 = 0
/// - ∂_i ∂_j = ∂_j ∂_i when |i-j| >= 2
/// - ∂_i ∂_{i+1} ∂_i = ∂_{i+1} ∂_i ∂_{i+1} (braid relation)
pub fn partial_i<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    // Compute f - s_i(f), then divide by (x_i - x_{i+1})
    let si_f = f.swap_vars(i);
    let numerator = f.clone() - si_f;

    // Division by (x_i - x_{i+1}):
    // For each term c * x^e in numerator, we need to express it as
    // (x_i - x_{i+1}) * something.
    //
    // Key identity: x_i^a * x_{i+1}^b - x_i^b * x_{i+1}^a = (x_i - x_{i+1}) * S
    // where S = Σ_{k=0}^{a-b-1} x_i^{a-1-k} * x_{i+1}^{b+k} (for a > b).
    //
    // Since f - s_i(f) is antisymmetric in x_i, x_{i+1}, it's always divisible.
    //
    // Algorithm: process terms greedily. For each monomial where exp[i] > exp[i+1],
    // contribute (exp[i] - exp[i+1]) terms to the quotient.
    divide_by_xi_minus_xi1(&numerator, i)
}

/// Isobaric (Demazure) divided difference π_i = ∂_i(x_i * f).
///
/// Equivalently: π_i(f) = (x_i * f - x_{i+1} * s_i(f)) / (x_i - x_{i+1}).
///
/// Satisfies:
/// - π_i^2 = π_i (idempotent)
/// - π_i π_{i+1} π_i = π_{i+1} π_i π_{i+1} (braid relation)
pub fn pi_i<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    let xi_f = f.mul_var(i);
    partial_i(&xi_f, i)
}

/// Demazure atom operator θ_i = π_i - id.
pub fn theta_i<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    pi_i(f, i) - f.clone()
}

/// Connective-K divided difference
/// `partial_i^beta(f) = partial_i((1 + beta*x_{i+1}) * f)`.
///
/// This specializes to [`partial_i`] at `beta = 0` and satisfies
/// `(partial_i^beta)^2 = -beta * partial_i^beta`.
pub fn demazure_lascoux_partial_i<C: Ring>(f: &MultiPoly<C>, i: usize, beta: &C) -> MultiPoly<C> {
    let deformed = f.clone() + f.mul_var(i + 1).scale(beta);
    partial_i(&deformed, i)
}

/// Connective-K isobaric divided difference
/// `pi_i^beta(f) = partial_i^beta(x_i * f)`.
///
/// This specializes to [`pi_i`] at `beta = 0` and is the operator defining
/// Lascoux polynomials.
pub fn demazure_lascoux_pi_i<C: Ring>(f: &MultiPoly<C>, i: usize, beta: &C) -> MultiPoly<C> {
    demazure_lascoux_partial_i(&f.mul_var(i), i, beta)
}

/// t-deformed Demazure operator:
/// (1 - t) π_i(f) + t s_i(f).
///
/// This is the operator-level ingredient used in the recursive construction
/// of several nonsymmetric Macdonald-type families.
pub fn tpi_i<C: Ring>(f: &MultiPoly<C>, i: usize, t: &C) -> MultiPoly<C> {
    let one_minus_t = C::one() - t.clone();
    pi_i(f, i).scale(&one_minus_t) + f.swap_vars(i).scale(t)
}

/// t-deformed atom operator:
/// (1 - t) θ_i(f) + t s_i(f).
pub fn ttheta_i<C: Ring>(f: &MultiPoly<C>, i: usize, t: &C) -> MultiPoly<C> {
    let one_minus_t = C::one() - t.clone();
    theta_i(f, i).scale(&one_minus_t) + f.swap_vars(i).scale(t)
}

/// Apply simple divided differences in the given order.
pub fn partial_word<C: Ring>(f: &MultiPoly<C>, word: &[usize]) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| partial_i(&acc, i))
}

/// Apply Demazure operators in the given order.
pub fn pi_word<C: Ring>(f: &MultiPoly<C>, word: &[usize]) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| pi_i(&acc, i))
}

/// Apply atom operators in the given order.
pub fn theta_word<C: Ring>(f: &MultiPoly<C>, word: &[usize]) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| theta_i(&acc, i))
}

/// Apply connective-K divided differences in the given order.
pub fn demazure_lascoux_partial_word<C: Ring>(
    f: &MultiPoly<C>,
    word: &[usize],
    beta: &C,
) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| {
        demazure_lascoux_partial_i(&acc, i, beta)
    })
}

/// Apply connective-K isobaric divided differences in the given order.
pub fn demazure_lascoux_pi_word<C: Ring>(
    f: &MultiPoly<C>,
    word: &[usize],
    beta: &C,
) -> MultiPoly<C> {
    word.iter()
        .fold(f.clone(), |acc, &i| demazure_lascoux_pi_i(&acc, i, beta))
}

/// Apply t-deformed Demazure operators in the given order.
pub fn tpi_word<C: Ring>(f: &MultiPoly<C>, word: &[usize], t: &C) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| tpi_i(&acc, i, t))
}

/// Apply t-deformed atom operators in the given order.
pub fn ttheta_word<C: Ring>(f: &MultiPoly<C>, word: &[usize], t: &C) -> MultiPoly<C> {
    word.iter().fold(f.clone(), |acc, &i| ttheta_i(&acc, i, t))
}

/// Divide a polynomial by (x_i - x_{i+1}).
///
/// The polynomial must be antisymmetric in x_i, x_{i+1} (guaranteed when
/// this is called from partial_i on f - s_i(f)).
fn divide_by_xi_minus_xi1<C: Ring>(f: &MultiPoly<C>, i: usize) -> MultiPoly<C> {
    let n = f.num_vars();
    let mut result = MultiPoly::zero(n);

    // For each term with exp[i] > exp[i+1], expand the geometric series
    for (exp, coeff) in f.terms() {
        let a = exp[i];
        let b = exp[i + 1];
        if a <= b {
            // Terms with a <= b should cancel in pairs due to antisymmetry.
            // We only process a > b terms.
            continue;
        }
        // x_i^a * x_{i+1}^b / (x_i - x_{i+1})
        // = x_i^{a-1} * x_{i+1}^b + x_i^{a-2} * x_{i+1}^{b+1} + ... + x_i^b * x_{i+1}^{a-1}
        for k in 0..(a - b) {
            let mut new_exp = exp.clone();
            new_exp[i] = a - 1 - k;
            new_exp[i + 1] = b + k;
            let entry = result.terms.entry(new_exp).or_insert_with(C::zero);
            *entry = entry.clone() + coeff.clone();
        }
    }

    result.strip_zeros();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_partial_i_squared_is_zero() {
        // ∂_0^2 = 0. Test on x_1^3 * x_2.
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let d1 = partial_i(&f, 0);
        let d2 = partial_i(&d1, 0);
        assert!(d2.is_zero(), "∂_0^2 should be 0, got {:?}", d2);
    }

    #[test]
    fn test_partial_i_on_xi() {
        // ∂_0(x_1) = (x_1 - x_2) / (x_1 - x_2) = 1
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let result = partial_i(&x1, 0);
        assert_eq!(result.coefficient(&[0, 0]), 1);
        assert_eq!(result.terms().len(), 1);
    }

    #[test]
    fn test_partial_i_on_x1_squared() {
        // ∂_0(x_1^2) = (x_1^2 - x_2^2)/(x_1 - x_2) = x_1 + x_2
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        let result = partial_i(&f, 0);
        assert_eq!(result.coefficient(&[1, 0]), 1);
        assert_eq!(result.coefficient(&[0, 1]), 1);
        assert_eq!(result.terms().len(), 2);
    }

    #[test]
    fn test_pi_idempotent() {
        // π_0^2 = π_0. Test on x_1^2 * x_2.
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![2, 1, 0]);
        let p1 = pi_i(&f, 0);
        let p2 = pi_i(&p1, 0);
        assert_eq!(p1, p2, "π_0 should be idempotent");
    }

    #[test]
    fn test_pi_on_x1() {
        // π_0(x_1) = ∂_0(x_1^2) = x_1 + x_2
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let result = pi_i(&x1, 0);
        assert_eq!(result.coefficient(&[1, 0]), 1);
        assert_eq!(result.coefficient(&[0, 1]), 1);
    }

    #[test]
    fn test_theta_on_x1_squared() {
        // θ_0(x_1^2) = π_0(x_1^2) - x_1^2 = x_1 x_2 + x_2^2
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        let result = theta_i(&f, 0);
        assert_eq!(result.coefficient(&[2, 0]), 0);
        assert_eq!(result.coefficient(&[1, 1]), 1);
        assert_eq!(result.coefficient(&[0, 2]), 1);
    }

    #[test]
    fn test_demazure_lascoux_specializations() {
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        assert_eq!(demazure_lascoux_partial_i(&f, 0, &0), partial_i(&f, 0));
        assert_eq!(demazure_lascoux_pi_i(&f, 0, &0), pi_i(&f, 0));
    }

    #[test]
    fn test_demazure_lascoux_pi_on_x1() {
        let x1: MultiPoly<i64> = MultiPoly::var(2, 0);
        let result = demazure_lascoux_pi_i(&x1, 0, &2);
        assert_eq!(result.coefficient(&[1, 0]), 1);
        assert_eq!(result.coefficient(&[0, 1]), 1);
        assert_eq!(result.coefficient(&[1, 1]), 2);
        assert_eq!(result.terms().len(), 3);
    }

    #[test]
    fn test_demazure_lascoux_quadratic_relations() {
        let beta = 2;
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let partial_once = demazure_lascoux_partial_i(&f, 0, &beta);
        let partial_twice = demazure_lascoux_partial_i(&partial_once, 0, &beta);
        assert_eq!(partial_twice, partial_once.scale(&(-beta)));

        let pi_once = demazure_lascoux_pi_i(&f, 0, &beta);
        let pi_twice = demazure_lascoux_pi_i(&pi_once, 0, &beta);
        assert_eq!(pi_twice, pi_once);
    }

    #[test]
    fn test_demazure_lascoux_braid_relation() {
        let beta = 2;
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let lhs = demazure_lascoux_pi_i(
            &demazure_lascoux_pi_i(&demazure_lascoux_pi_i(&f, 0, &beta), 1, &beta),
            0,
            &beta,
        );
        let rhs = demazure_lascoux_pi_i(
            &demazure_lascoux_pi_i(&demazure_lascoux_pi_i(&f, 1, &beta), 0, &beta),
            1,
            &beta,
        );
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_tpi_specializations() {
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        assert_eq!(tpi_i(&f, 0, &0), pi_i(&f, 0));
        assert_eq!(tpi_i(&f, 0, &1), f.swap_vars(0));
    }

    #[test]
    fn test_ttheta_specializations() {
        let f: MultiPoly<i64> = MultiPoly::x_power(2, vec![2, 0]);
        assert_eq!(ttheta_i(&f, 0, &0), theta_i(&f, 0));
        assert_eq!(ttheta_i(&f, 0, &1), f.swap_vars(0));
    }

    #[test]
    fn test_braid_relation_partial() {
        // ∂_0 ∂_1 ∂_0 = ∂_1 ∂_0 ∂_1 on x_1^3 * x_2^2 * x_3
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 2, 1]);
        let lhs = partial_i(&partial_i(&partial_i(&f, 0), 1), 0);
        let rhs = partial_i(&partial_i(&partial_i(&f, 1), 0), 1);
        assert_eq!(lhs, rhs, "braid relation for ∂ failed");
    }

    #[test]
    fn test_braid_relation_pi() {
        // π_0 π_1 π_0 = π_1 π_0 π_1 on x_1^3 * x_2
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let lhs = pi_i(&pi_i(&pi_i(&f, 0), 1), 0);
        let rhs = pi_i(&pi_i(&pi_i(&f, 1), 0), 1);
        assert_eq!(lhs, rhs, "braid relation for π failed");
    }

    #[test]
    fn test_commutation_partial() {
        // ∂_0 ∂_2 = ∂_2 ∂_0 (|0-2| >= 2)
        let f: MultiPoly<i64> = MultiPoly::x_power(4, vec![3, 1, 2, 1]);
        let lhs = partial_i(&partial_i(&f, 0), 2);
        let rhs = partial_i(&partial_i(&f, 2), 0);
        assert_eq!(lhs, rhs, "commutation for ∂ failed");
    }

    #[test]
    fn test_word_application() {
        let f: MultiPoly<i64> = MultiPoly::x_power(3, vec![3, 1, 0]);
        let word = vec![0, 1, 0];
        assert_eq!(pi_word(&f, &word), pi_i(&pi_i(&pi_i(&f, 0), 1), 0));
        assert_eq!(
            theta_word(&f, &word),
            theta_i(&theta_i(&theta_i(&f, 0), 1), 0)
        );
    }
}
