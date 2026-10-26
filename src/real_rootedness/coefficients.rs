//! Coefficient-sequence properties: log-concavity, unimodality,
//! ultra-log-concavity, palindromicity, gamma-positivity, the Stapledon
//! decomposition and the Hermite--Biehler parts.
//!
//! These checks look only at the coefficient sequence (ascending degree order)
//! and never at roots.  The real-rootedness entry points live in the sibling
//! `real_rooted` module.

use crate::vec_poly;
use crate::Polynomial;
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

/// Check if a coefficient sequence is log-concave.
///
/// A sequence a_0, a_1, ..., a_d is log-concave if a_k^2 >= a_{k-1} * a_{k+1}
/// for all 1 <= k <= d-1.
pub fn is_log_concave(coeffs: &[i64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    for i in 1..coeffs.len() - 1 {
        let lhs = i128::from(coeffs[i]) * i128::from(coeffs[i]);
        let rhs = i128::from(coeffs[i - 1]) * i128::from(coeffs[i + 1]);
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a `BigInt` coefficient sequence is log-concave.
pub fn is_log_concave_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }
    for i in 1..coeffs.len() - 1 {
        let lhs = &coeffs[i] * &coeffs[i];
        let rhs = &coeffs[i - 1] * &coeffs[i + 1];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a coefficient sequence is unimodal.
///
/// A sequence is unimodal if it is weakly increasing up to some index and
/// weakly decreasing afterwards.  Plateaus are allowed.
pub fn is_unimodal(coeffs: &[i64]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }

    let mut i = 1;
    while i < coeffs.len() && coeffs[i - 1] <= coeffs[i] {
        i += 1;
    }
    while i < coeffs.len() && coeffs[i - 1] >= coeffs[i] {
        i += 1;
    }
    i == coeffs.len()
}

/// Check if a `BigInt` coefficient sequence is unimodal.
pub fn is_unimodal_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    if coeffs.len() <= 2 {
        return true;
    }

    let mut i = 1;
    while i < coeffs.len() && coeffs[i - 1] <= coeffs[i] {
        i += 1;
    }
    while i < coeffs.len() && coeffs[i - 1] >= coeffs[i] {
        i += 1;
    }
    i == coeffs.len()
}

/// Check if a coefficient sequence is ultra-log-concave (satisfies Newton's inequalities).
///
/// A polynomial a_0 + a_1 t + ... + a_d t^d of degree d is ultra-log-concave if
/// a_k / C(d,k) is log-concave, where C(d,k) is the binomial coefficient.
/// This is [`is_ultra_log_concave_of_degree`] with the nominal degree equal to
/// the actual degree.
///
/// This is implied by real-rootedness but is strictly weaker.
pub fn is_ultra_log_concave(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_ultra_log_concave_bigint_coeffs(&coeffs)
}

/// Check if a `BigInt` coefficient sequence is ultra-log-concave.
pub fn is_ultra_log_concave_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = vec_poly::degree(coeffs).unwrap_or(0);
    is_ultra_log_concave_of_degree_bigint_coeffs(coeffs, d)
}

/// Check ultra-log-concavity with respect to a prescribed nominal degree `n`.
///
/// The sequence `a_0, ..., a_n` (with `a_k = 0` past the end of `coeffs`) is
/// *ultra-log-concave of order `n`* (ULC(n)) if `a_k / C(n, k)` is
/// log-concave, i.e.
///
/// ```text
/// a_k^2 · C(n, k-1) · C(n, k+1) ≥ a_(k-1) · a_(k+1) · C(n, k)^2   for 1 ≤ k ≤ n-1.
/// ```
///
/// For `n = deg p` this is [`is_ultra_log_concave`]; the limiting condition
/// as `n → ∞` is [`is_ultra_log_concave_infinite`] (`k a_k^2 ≥ (k+1) a_(k-1) a_(k+1)`).  ULC(n) implies ULC(n')
/// for every `n' ≥ n` for nonnegative sequences, so a nominal degree above the
/// actual degree (for example the rank of a matroid, or the size of a ground
/// set) is a *weaker* condition.  By Newton's inequalities, a real-rooted
/// polynomial with nonnegative coefficients and `deg p ≤ n` is ULC(n).
///
/// Conventions:
///
/// - If `deg p > n` the sequence is not supported on `0..=n`, and the
///   function returns `false`.
/// - The zero polynomial is ULC(n) for every `n` (returns `true`).
/// - Trailing zeros in `coeffs` are ignored.  As in [`is_log_concave`], the
///   inequalities are checked pointwise; the absence of internal zeros and the
///   signs of the coefficients are not checked.
///
/// ```
/// use polytool::{is_ultra_log_concave, is_ultra_log_concave_of_degree};
///
/// // 1 + 2t + t^2 = (1 + t)^2 is ULC(2), hence ULC(n) for all n ≥ 2.
/// assert!(is_ultra_log_concave_of_degree(&[1, 2, 1], 2));
/// assert!(is_ultra_log_concave_of_degree(&[1, 2, 1], 5));
/// assert!(!is_ultra_log_concave_of_degree(&[1, 2, 1], 1));
/// // 1 + 3t + 3t^2 is not ULC at its degree 2 (9 < 3 · 2^2), but it is
/// // ULC(3): 9 · C(3,0) · C(3,2) = 27 ≥ 3 · C(3,1)^2 = 27.
/// assert!(!is_ultra_log_concave(&[1, 3, 3]));
/// assert!(is_ultra_log_concave_of_degree(&[1, 3, 3], 3));
/// ```
pub fn is_ultra_log_concave_of_degree(coeffs: &[i64], n: usize) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_ultra_log_concave_of_degree_bigint_coeffs(&coeffs, n)
}

/// Check infinite ultra-log-concavity, ULC(∞), of a coefficient sequence.
///
/// The sequence `a_0, a_1, ...` is *ultra-log-concave of infinite order* if
/// `a_k k!` is log-concave, i.e.
///
/// ```text
/// k · a_k^2 ≥ (k + 1) · a_(k-1) · a_(k+1)   for every k ≥ 1.
/// ```
///
/// This is the limit `n → ∞` of the degree-based condition
/// [`is_ultra_log_concave_of_degree`] (ULC(n) uses the factors `C(n, k)`, and
/// `a_k / C(n, k)` tends to a multiple of `a_k k!` as `n` grows), and it is the
/// weakest of the ULC(n) family: for nonnegative sequences ULC(n) implies
/// ULC(∞) for every `n`.  The condition does not depend on the nominal degree.  In contrast, [`is_ultra_log_concave`] tests ULC at the
/// nominal degree `deg p`, which depends on the polynomial's length.  For
/// example `1 + 3t + 3t^2` is not ULC at its degree 2, but it is ULC(∞).
///
/// Conventions: trailing zeros are ignored, the zero polynomial returns `true`,
/// and (as for [`is_log_concave`]) the inequalities are checked pointwise, with
/// no check of coefficient signs or internal zeros.
///
/// ```
/// use polytool::{is_ultra_log_concave, is_ultra_log_concave_infinite};
///
/// assert!(!is_ultra_log_concave(&[1, 3, 3]));
/// assert!(is_ultra_log_concave_infinite(&[1, 3, 3]));
/// // 1 + t + t^2 + t^3 fails at k = 1: 1 · a_1^2 = 1 < 2 · a_0 a_2 = 2.
/// assert!(!is_ultra_log_concave_infinite(&[1, 1, 1, 1]));
/// ```
pub fn is_ultra_log_concave_infinite(coeffs: &[i64]) -> bool {
    let coeffs: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_ultra_log_concave_infinite_bigint_coeffs(&coeffs)
}

/// Check ULC(∞) for `BigInt` coefficients; see [`is_ultra_log_concave_infinite`].
pub fn is_ultra_log_concave_infinite_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(d) = vec_poly::degree(coeffs) else {
        return true;
    };
    // For k ≥ d the right-hand side vanishes, so only 1 ≤ k < d needs checking.
    (1..d).all(|k| {
        BigInt::from(k) * &coeffs[k] * &coeffs[k]
            >= BigInt::from(k + 1) * &coeffs[k - 1] * &coeffs[k + 1]
    })
}

/// Check ultra-log-concavity of order `n` for `BigInt` coefficients.
///
/// Same convention as [`is_ultra_log_concave_of_degree`].
pub fn is_ultra_log_concave_of_degree_bigint_coeffs(coeffs: &[BigInt], n: usize) -> bool {
    let Some(d) = vec_poly::degree(coeffs) else {
        return true;
    };
    if d > n {
        return false;
    }
    // For k ≥ d the right-hand side vanishes and the left-hand side is a
    // square times positive binomials, so only 1 ≤ k < d needs checking.
    let mut binom = vec![BigInt::from(1); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(n - k + 1) / BigInt::from(k);
    }
    for k in 1..d {
        let lhs = coeffs[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &coeffs[k - 1] * &coeffs[k + 1] * &binom[k] * &binom[k];
        if lhs < rhs {
            return false;
        }
    }
    true
}

/// Check if a polynomial has palindromic (symmetric) coefficients.
///
/// A polynomial a_0 + a_1 t + ... + a_d t^d is palindromic if a_i = a_{d-i}
/// for all i, equivalently t^d p(1/t) = p(t).
///
/// Returns true for the zero polynomial and constant polynomials.
/// Trailing zeros are trimmed before checking.
pub fn is_palindromic(coeffs: &[i64]) -> bool {
    let d = match coeffs.iter().rposition(|&c| c != 0) {
        Some(d) => d,
        None => return true,
    };
    for i in 0..=d / 2 {
        if coeffs[i] != coeffs[d - i] {
            return false;
        }
    }
    true
}

/// Check if a `BigInt` coefficient polynomial is palindromic.
pub fn is_palindromic_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = match coeffs.iter().rposition(|c| !c.is_zero()) {
        Some(d) => d,
        None => return true,
    };
    for i in 0..=d / 2 {
        if coeffs[i] != coeffs[d - i] {
            return false;
        }
    }
    true
}

/// Compute the gamma coefficients of a palindromic polynomial.
///
/// For a palindromic polynomial p(t) of degree d, there is a unique expansion
///
/// ```text
/// p(t) = sum_{i=0}^{floor(d/2)} gamma_i * t^i * (1+t)^{d-2i}
/// ```
///
/// Returns `None` if the polynomial is not palindromic or if an exact gamma
/// coefficient does not fit in `i64`; use
/// [`gamma_coefficients_bigint_coeffs`] for the no-overflow result.
/// Returns `Some(gamma_coefficients)` where `gamma[i]` is the coefficient
/// of t^i (1+t)^{d-2i} in the expansion.
///
/// The polynomial is gamma-positive iff all returned coefficients are non-negative.
pub fn gamma_coefficients(coeffs: &[i64]) -> Option<Vec<i64>> {
    let coeffs = coeffs.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    gamma_coefficients_bigint_coeffs(&coeffs)?
        .into_iter()
        .map(|gamma| gamma.to_i64())
        .collect()
}

/// Compute gamma coefficients for a palindromic `BigInt` coefficient polynomial.
pub fn gamma_coefficients_bigint_coeffs(coeffs: &[BigInt]) -> Option<Vec<BigInt>> {
    let d = match coeffs.iter().rposition(|c| !c.is_zero()) {
        Some(d) => d,
        None => return Some(vec![]),
    };
    if !is_palindromic_bigint_coeffs(coeffs) {
        return None;
    }
    let half = d / 2;

    let mut binomials = vec![vec![BigInt::from(0); d + 1]; d + 1];
    for n in 0..=d {
        binomials[n][0] = BigInt::from(1);
        for k in 1..=n {
            binomials[n][k] = &binomials[n - 1][k - 1] + &binomials[n - 1][k];
        }
    }

    let mut gamma = vec![BigInt::from(0); half + 1];
    for i in 0..=half {
        let mut val = coeffs[i].clone();
        for j in 0..i {
            val -= &gamma[j] * &binomials[d - 2 * j][i - j];
        }
        gamma[i] = val;
    }

    Some(gamma)
}

/// Check if a polynomial is gamma-positive.
///
/// A palindromic polynomial is gamma-positive if all its gamma coefficients
/// (in the basis {t^i (1+t)^{d-2i}}) are non-negative.
///
/// Returns false for non-palindromic polynomials.
pub fn is_gamma_positive(coeffs: &[i64]) -> bool {
    match gamma_coefficients(coeffs) {
        None => false,
        Some(gammas) => gammas.iter().all(|&g| g >= 0),
    }
}

/// Check if a `BigInt` coefficient polynomial is gamma-positive.
pub fn is_gamma_positive_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    match gamma_coefficients_bigint_coeffs(coeffs) {
        None => false,
        Some(gammas) => gammas.iter().all(|g| g >= &BigInt::zero()),
    }
}

/// Remove all initial zeros, i.e. divide out the largest power of `t`.
pub fn strip_initial_zeros(coeffs: &[i64]) -> &[i64] {
    let start = match coeffs.iter().position(|&c| c != 0) {
        Some(start) => start,
        None => return &[],
    };
    &coeffs[start..]
}

/// Remove all initial zeros from a `BigInt` coefficient vector.
pub fn strip_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> &[BigInt] {
    let start = match coeffs.iter().position(|c| !c.is_zero()) {
        Some(start) => start,
        None => return &[],
    };
    &coeffs[start..]
}

/// Check palindromicity after dividing out all factors of `t`.
pub fn is_palindromic_ignoring_initial_zeros(coeffs: &[i64]) -> bool {
    is_palindromic(strip_initial_zeros(coeffs))
}

/// Check palindromicity after dividing out all factors of `t`.
pub fn is_palindromic_ignoring_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    is_palindromic_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Compute gamma-coefficients after dividing out all factors of `t`.
pub fn gamma_coefficients_ignoring_initial_zeros(coeffs: &[i64]) -> Option<Vec<i64>> {
    gamma_coefficients(strip_initial_zeros(coeffs))
}

/// Compute gamma coefficients after dividing out all factors of `t`.
pub fn gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(
    coeffs: &[BigInt],
) -> Option<Vec<BigInt>> {
    gamma_coefficients_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Check gamma-positivity after dividing out all factors of `t`.
pub fn is_gamma_positive_ignoring_initial_zeros(coeffs: &[i64]) -> bool {
    is_gamma_positive(strip_initial_zeros(coeffs))
}

/// Check gamma-positivity after dividing out all factors of `t`.
pub fn is_gamma_positive_ignoring_initial_zeros_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    is_gamma_positive_bigint_coeffs(strip_initial_zeros_bigint_coeffs(coeffs))
}

/// Compute the Stapledon decomposition with respect to a degree bound `n`.
///
/// For a polynomial `p(x)` of degree at most `n`, there is a unique decomposition
///
/// ```text
/// p(x) = a(x) + x b(x),
/// ```
///
/// where `a(x)` is symmetric with center `n/2` and `b(x)` is symmetric with center
/// `(n-1)/2`. This returns the coefficient vectors of `(a(x), b(x))` in ascending
/// degree order, or `None` if `deg(p) > n`.
///
/// Panics if a coefficient of `a` or `b` does not fit in `i64` (their
/// coefficients are integer combinations of the input's); use
/// [`stapledon_decomposition_bigint_coeffs`] for large values.
pub fn stapledon_decomposition(coeffs: &[i64], n: usize) -> Option<(Vec<i64>, Vec<i64>)> {
    let coeffs = coeffs.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    stapledon_decomposition_bigint_coeffs(&coeffs, n).map(|(a, b)| {
        (
            a.iter()
                .map(|c| i64::try_from(c).expect("Stapledon coefficient too large for i64"))
                .collect(),
            b.iter()
                .map(|c| i64::try_from(c).expect("Stapledon coefficient too large for i64"))
                .collect(),
        )
    })
}

/// Compute the Stapledon decomposition for a `BigInt` coefficient polynomial.
pub fn stapledon_decomposition_bigint_coeffs(
    coeffs: &[BigInt],
    n: usize,
) -> Option<(Vec<BigInt>, Vec<BigInt>)> {
    crate::polynomial::Polynomial::<BigInt>::new(coeffs.to_vec())
        .stapledon_decomposition(n)
        .map(|(a, b)| (a.coeffs().to_vec(), b.coeffs().to_vec()))
}

/// Return the Hermite--Biehler decomposition `(E, O)` where
///
/// ```text
/// p(t) = E(t^2) + t O(t^2).
/// ```
pub fn hermite_biehler_parts(coeffs: &[i64]) -> (Vec<i64>, Vec<i64>) {
    let p = Polynomial::<i64>::from_i64_coeffs(coeffs);
    let (even, odd) = p.hermite_biehler_decomposition();
    (even.coeffs().to_vec(), odd.coeffs().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::real_rootedness::Q;

    #[test]
    fn test_log_concave() {
        assert!(is_log_concave(&[1, 2, 1]));
        assert!(is_log_concave(&[1, 3, 3, 1]));
        assert!(!is_log_concave(&[1, 1, 3]));
    }

    #[test]
    fn ultra_log_concave_infinite_basics_and_relation_to_degree_based() {
        assert!(is_ultra_log_concave_infinite(&[]));
        assert!(is_ultra_log_concave_infinite(&[0, 0]));
        assert!(is_ultra_log_concave_infinite(&[5]));
        assert!(is_ultra_log_concave_infinite(&[1, 3, 3]));
        assert!(!is_ultra_log_concave_infinite(&[1, 1, 1, 1]));
        // Trailing zeros are ignored.
        assert!(is_ultra_log_concave_infinite(&[1, 3, 3, 0, 0]));
        // Poisson weights 1/k! are exactly the boundary case.
        assert!(is_ultra_log_concave_infinite(&[
            720, 720, 360, 120, 30, 6, 1
        ]));
        // Real-rooted nonnegative polynomials of degree <= n are ULC(n), hence ULC(inf):
        // check ULC(n) => ULC(inf) on small nonnegative sequences.
        for a in 0..4i64 {
            for b in 0..4i64 {
                for c in 0..4i64 {
                    for d in 0..4i64 {
                        let p = [a, b, c, d];
                        if is_ultra_log_concave(&p) {
                            assert!(is_ultra_log_concave_infinite(&p), "{p:?}");
                        }
                        if is_ultra_log_concave_of_degree(&p, 9) {
                            assert!(is_ultra_log_concave_infinite(&p), "{p:?}");
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_unimodal() {
        assert!(is_unimodal(&[]));
        assert!(is_unimodal(&[1]));
        assert!(is_unimodal(&[1, 2, 2, 1]));
        assert!(is_unimodal(&[3, 2, 1]));
        assert!(is_unimodal(&[1, 2, 3]));
        assert!(!is_unimodal(&[1, 3, 2, 4]));
    }

    #[test]
    fn test_log_concave_uses_wide_products() {
        assert!(is_log_concave(&[i64::MAX, i64::MAX, i64::MAX]));
        assert!(!is_log_concave(&[i64::MAX, 0, i64::MAX]));
    }

    #[test]
    fn test_ultra_log_concave() {
        assert!(is_ultra_log_concave(&[1, 3, 3, 1])); // (1+t)^3
        assert!(is_ultra_log_concave(&[1, 4, 6, 4, 1])); // (1+t)^4
    }

    #[test]
    fn test_bigint_shape_properties() {
        let scale = BigInt::from(10u64).pow(40);
        let coeffs = [1_i64, 3, 3, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect::<Vec<_>>();
        assert!(is_palindromic_bigint_coeffs(&coeffs));
        assert!(is_unimodal_bigint_coeffs(&coeffs));
        assert!(is_log_concave_bigint_coeffs(&coeffs));
        assert!(is_ultra_log_concave_bigint_coeffs(&coeffs));
    }

    #[test]
    fn test_hermite_biehler_parts() {
        let (even, odd) = hermite_biehler_parts(&[1, 2, 3, 4, 5]);
        assert_eq!(even, vec![1, 3, 5]);
        assert_eq!(odd, vec![2, 4]);
    }

    // -- Palindromic / gamma tests --

    #[test]
    fn test_is_palindromic() {
        assert!(is_palindromic(&[1, 2, 1]));
        assert!(is_palindromic(&[1, 11, 11, 1]));
        assert!(is_palindromic(&[1, 4, 6, 4, 1]));
        assert!(is_palindromic(&[1]));
        assert!(is_palindromic(&[]));
        assert!(!is_palindromic(&[0, 1, 1]));
        assert!(!is_palindromic(&[0, 0, 1, 2, 1]));
        assert!(!is_palindromic(&[1, 2, 3]));
        // Trailing zeros trimmed: [1, 2, 1, 0] is palindromic
        assert!(is_palindromic(&[1, 2, 1, 0]));
        assert!(!is_palindromic(&[0, 1, 2, 1, 0]));
    }

    #[test]
    fn test_gamma_coefficients() {
        // (1+t)^3 = [1,3,3,1]: gamma = [1, 0]
        assert_eq!(gamma_coefficients(&[1, 3, 3, 1]), Some(vec![1, 0]));

        // A_4(t) = [1,11,11,1]: gamma_0=1, gamma_1=11-3=8
        assert_eq!(gamma_coefficients(&[1, 11, 11, 1]), Some(vec![1, 8]));

        // (1+t)^4 = [1,4,6,4,1]: gamma_0=1, gamma_1=4-C(4,1)=0, gamma_2=6-C(4,2)=0
        assert_eq!(gamma_coefficients(&[1, 4, 6, 4, 1]), Some(vec![1, 0, 0]));

        // Not palindromic
        assert_eq!(gamma_coefficients(&[1, 2, 3]), None);

        // Zero/constant
        assert_eq!(gamma_coefficients(&[]), Some(vec![]));
        assert_eq!(gamma_coefficients(&[5]), Some(vec![5]));
        assert_eq!(gamma_coefficients(&[0, 1, 1]), None);
        assert_eq!(gamma_coefficients(&[0, 0, 1, 3, 3, 1]), None);

        // [1,2,2,1] degree 3: gamma_0=1, gamma_1=2-3=-1
        assert_eq!(gamma_coefficients(&[1, 2, 2, 1]), Some(vec![1, -1]));

        // The exact gamma expansion exists, but later coefficients do not fit
        // in i64. The convenience wrapper must fail cleanly without overflowing
        // its former i128 Pascal triangle.
        let mut high_degree = vec![0; 201];
        high_degree[0] = 1;
        high_degree[200] = 1;
        assert_eq!(gamma_coefficients(&high_degree), None);
    }

    #[test]
    fn test_bigint_gamma_coefficients() {
        let scale = BigInt::from(10u64).pow(40);
        let coeffs = [1_i64, 11, 11, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &scale)
            .collect::<Vec<_>>();
        assert_eq!(
            gamma_coefficients_bigint_coeffs(&coeffs),
            Some(vec![scale.clone(), BigInt::from(8) * scale])
        );
    }

    #[test]
    fn test_is_gamma_positive() {
        assert!(is_gamma_positive(&[1, 3, 3, 1]));
        assert!(is_gamma_positive(&[1, 11, 11, 1]));
        assert!(is_gamma_positive(&[1, 4, 6, 4, 1]));
        assert!(!is_gamma_positive(&[0, 1, 1]));
        // [1,2,2,1] has gamma_1 = -1, NOT gamma-positive
        assert!(!is_gamma_positive(&[1, 2, 2, 1]));
        // Non-palindromic: not gamma-positive
        assert!(!is_gamma_positive(&[1, 2, 3]));
    }

    #[test]
    fn test_strip_initial_zero_wrappers() {
        let empty: &[i64] = &[];
        assert_eq!(strip_initial_zeros(empty), empty);
        assert_eq!(strip_initial_zeros(&[0, 0, 1, 2, 1]), &[1, 2, 1]);
        assert!(is_palindromic_ignoring_initial_zeros(&[0, 1, 1]));
        assert!(is_palindromic_ignoring_initial_zeros(&[0, 0, 1, 2, 1]));
        assert_eq!(
            gamma_coefficients_ignoring_initial_zeros(&[0, 1, 1]),
            Some(vec![1])
        );
        assert_eq!(
            gamma_coefficients_ignoring_initial_zeros(&[0, 0, 1, 3, 3, 1]),
            Some(vec![1, 0])
        );
        assert!(is_gamma_positive_ignoring_initial_zeros(&[0, 1, 1]));
    }

    #[test]
    fn test_stapledon_decomposition_basic() {
        let (a, b) = stapledon_decomposition(&[1, 2, 3], 2).unwrap();
        assert_eq!(a, vec![1, 0, 1]);
        assert_eq!(b, vec![2, 2]);
    }

    #[test]
    fn test_stapledon_decomposition_palindromic() {
        let (a, b) = stapledon_decomposition(&[1, 11, 11, 1], 3).unwrap();
        assert_eq!(a, vec![1, 11, 11, 1]);
        assert_eq!(b, Vec::<i64>::new());
    }

    #[test]
    fn test_stapledon_decomposition_with_larger_bound() {
        let (a, b) = stapledon_decomposition(&[1, 4, 1], 4).unwrap();
        assert_eq!(a, vec![1, 5, 6, 5, 1]);
        assert_eq!(b, vec![-1, -5, -5, -1]);
    }

    #[test]
    fn test_stapledon_decomposition_rejects_large_degree() {
        assert_eq!(stapledon_decomposition(&[1, 2, 3], 1), None);
    }

    /// Independent ULC(n) oracle: log-concavity of `a_k / C(n, k)` over ℚ.
    fn ulc_oracle(coeffs: &[i64], n: usize) -> bool {
        if vec_poly::degree(coeffs).is_some_and(|d| d > n) {
            return false;
        }
        let normalized: Vec<Q> = (0..=n)
            .map(|k| {
                Q::new(
                    BigInt::from(vec_poly::coeff(coeffs, k)),
                    crate::numeric::binomial_bigint(n, k),
                )
            })
            .collect();
        (1..n).all(|k| &normalized[k] * &normalized[k] >= &normalized[k - 1] * &normalized[k + 1])
    }

    #[test]
    fn ulc_of_degree_matches_rational_oracle_exhaustively() {
        let values = [-1i64, 0, 1, 2, 3, 5];
        for len in 0..=4u32 {
            for code in 0..values.len().pow(len) {
                let mut rest = code;
                let coeffs: Vec<i64> = (0..len)
                    .map(|_| {
                        let v = values[rest % values.len()];
                        rest /= values.len();
                        v
                    })
                    .collect();
                for n in 0..=6 {
                    let expected = ulc_oracle(&coeffs, n);
                    assert_eq!(
                        is_ultra_log_concave_of_degree(&coeffs, n),
                        expected,
                        "{coeffs:?}, n = {n}"
                    );
                    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
                    assert_eq!(
                        is_ultra_log_concave_of_degree_bigint_coeffs(&big, n),
                        expected
                    );
                }
                let d = vec_poly::degree(&coeffs).unwrap_or(0);
                assert_eq!(
                    is_ultra_log_concave(&coeffs),
                    is_ultra_log_concave_of_degree(&coeffs, d),
                    "{coeffs:?}"
                );
                // ULC(n) implies ULC(n + 1) for nonnegative sequences.
                if coeffs.iter().all(|&c| c >= 0) {
                    for n in d..6 {
                        if is_ultra_log_concave_of_degree(&coeffs, n) {
                            assert!(is_ultra_log_concave_of_degree(&coeffs, n + 1));
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn ulc_of_degree_edge_cases() {
        assert!(is_ultra_log_concave_of_degree(&[], 0));
        assert!(is_ultra_log_concave_of_degree(&[0, 0], 3));
        assert!(is_ultra_log_concave_of_degree(&[4], 0));
        assert!(is_ultra_log_concave_of_degree(&[1, 2, 0, 0], 1));
        assert!(!is_ultra_log_concave_of_degree(&[1, 2, 1], 1));
        assert!(!is_ultra_log_concave(&[1, 3, 3]));
        assert!(is_ultra_log_concave_of_degree(&[1, 3, 3], 3));
        // Large coefficients: (1 + t)^60 is ULC(60) with equality throughout.
        let row: Vec<BigInt> = (0..=60)
            .map(|k| crate::numeric::binomial_bigint(60, k))
            .collect();
        assert!(is_ultra_log_concave_of_degree_bigint_coeffs(&row, 60));
        let mut bumped = row.clone();
        bumped[30] -= 1;
        assert!(!is_ultra_log_concave_of_degree_bigint_coeffs(&bumped, 60));
        assert!(is_ultra_log_concave_of_degree_bigint_coeffs(&row, 61));
    }
}
