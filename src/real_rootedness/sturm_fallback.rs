//! Sturm-chain routes: real-rootedness, simple roots, root approximation and
//! interlacing by explicit root isolation.
//!
//! These are slower than the default PRS and Bézout paths and are kept for
//! cases where root locations are needed; see [`crate::sturm`] for the chain
//! itself.

use super::resultants::resultant;
use super::Q;
use crate::sturm::SturmChain;
use crate::vec_poly;
use crate::BigRational;
use crate::Polynomial;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Check if a polynomial is real-rooted using Sturm chains with exact arithmetic.
///
/// This is the slower rational Sturm-chain method. Prefer
/// [`is_real_rooted`] unless you need the Sturm-chain route explicitly.
///
/// Returns true for the zero polynomial and for constant/linear polynomials.
///
/// [`is_real_rooted`]: super::is_real_rooted
pub fn is_real_rooted_sturm(coeffs: &[i64]) -> bool {
    let Some(start) = coeffs.iter().position(|&c| c != 0) else {
        return true;
    };
    let end = coeffs
        .iter()
        .rposition(|&c| c != 0)
        .expect("nonzero coefficient must have a last position");
    let trimmed = &coeffs[start..=end];
    if trimmed.len() <= 2 {
        return true;
    }

    let sc = SturmChain::from_i64_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    sc.count_real_roots() == sf_degree
}

/// Check if a `BigInt`-coefficient polynomial is real-rooted using Sturm
/// chains with exact rational arithmetic.
pub fn is_real_rooted_sturm_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(start) = coeffs.iter().position(|c| !c.is_zero()) else {
        return true;
    };
    let end = coeffs
        .iter()
        .rposition(|c| !c.is_zero())
        .expect("nonzero coefficient must have a last position");
    let trimmed = &coeffs[start..=end];
    if trimmed.len() <= 2 {
        return true;
    }

    let sc = SturmChain::from_bigint_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    sc.count_real_roots() == sf_degree
}

/// Check whether a polynomial has only simple roots.
///
/// Returns `false` for the zero polynomial and `true` for nonzero constants.
pub fn has_simple_roots(coeffs: &[i64]) -> bool {
    let p = Polynomial::<Q>::new(
        coeffs
            .iter()
            .map(|&c| Q::from_integer(BigInt::from(c)))
            .collect(),
    );
    p.has_simple_roots()
}

/// Check whether a BigInt-coefficient polynomial has only simple roots.
///
/// Returns `false` for the zero polynomial and `true` for nonzero constants.
pub fn has_simple_roots_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let p = Polynomial::<Q>::new(coeffs.iter().map(|c| Q::from_integer(c.clone())).collect());
    p.has_simple_roots()
}

/// Find all real roots of a polynomial as rational interval midpoints
/// (approximations typed as rationals; formerly `real_roots`).
///
/// Returns the **distinct** roots in increasing order (multiplicities are not
/// reported).  A root `0` (from initial zero coefficients) and the root of an
/// input that is linear after removing powers of `t` are exact; all other
/// roots, even rational ones of a linear factor, are the midpoints of
/// isolating intervals of width at most `1e-6`, so they are approximations
/// with error below `1e-6`.
/// Use [`is_real_rooted`] for an exact yes/no decision.
///
/// Returns `None` if the polynomial is not real-rooted.
/// Returns `Some(vec![])` for constant/zero polynomials.
///
/// [`is_real_rooted`]: super::is_real_rooted
pub fn approximate_real_roots(coeffs: &[i64]) -> Option<Vec<BigRational>> {
    let Some(start) = coeffs.iter().position(|&c| c != 0) else {
        return Some(vec![]);
    };
    let end = coeffs
        .iter()
        .rposition(|&c| c != 0)
        .expect("nonzero polynomial has a final nonzero coefficient");
    let trimmed = &coeffs[start..=end];
    let mut roots = if start > 0 {
        vec![Q::zero()]
    } else {
        Vec::new()
    };
    if trimmed.len() <= 1 {
        return Some(roots);
    }
    if trimmed.len() == 2 {
        let a = Q::from_integer(BigInt::from(trimmed[0]));
        let b = Q::from_integer(BigInt::from(trimmed[1]));
        roots.push(-a / b);
        roots.sort();
        return Some(roots);
    }

    let sc = SturmChain::from_i64_coeffs(trimmed);
    let sf_degree = sc.square_free_degree();
    let eps = Q::new(BigInt::from(1), BigInt::from(1_000_000));
    let intervals = sc.isolate_roots(&eps);

    if intervals.len() == sf_degree {
        roots.extend(
            intervals
                .into_iter()
                .map(|(lo, hi)| (&lo + &hi) / Q::from_integer(BigInt::from(2))),
        );
        roots.sort();
        Some(roots)
    } else {
        None
    }
}

/// Deprecated name of [`approximate_real_roots`].
#[deprecated(since = "0.6.0", note = "renamed to `approximate_real_roots`")]
pub fn real_roots(coeffs: &[i64]) -> Option<Vec<BigRational>> {
    approximate_real_roots(coeffs)
}

/// Check if roots of p interlace roots of q, using Sturm chains.
///
/// This is the slower Sturm-chain method that isolates all roots. Prefer
/// [`check_interlacing`] (Bézout-based) for better performance.
///
/// Convention: `check_interlacing_sturm(f, g)` checks f ≪ g (f has smaller degree).
///
/// Returns `None` if either polynomial is not real-rooted.
/// Returns `Some(false)` if deg(g) ≠ deg(f) + 1 for nonzero `f`; `f = 0`
/// gives `Some(true)` for real-rooted `g`, the convention of
/// [`check_interlacing`].
///
/// [`check_interlacing`]: super::check_interlacing
pub fn check_interlacing_sturm(f: &[i64], g: &[i64]) -> Option<bool> {
    if !is_real_rooted_sturm(f) || !is_real_rooted_sturm(g) {
        return None;
    }

    // Same convention as `check_interlacing`: 0 ≪ g for every real-rooted g.
    let Some(df) = vec_poly::degree(f) else {
        return Some(true);
    };
    let Some(dg) = vec_poly::degree(g) else {
        return Some(false);
    };
    if dg != df + 1 {
        return Some(false);
    }
    if !has_simple_roots(f) || !has_simple_roots(g) || resultant(f, g).is_zero() {
        return Some(false);
    }

    let f_chain = SturmChain::from_i64_coeffs(f);
    let g_chain = SturmChain::from_i64_coeffs(g);
    let mut epsilon = Q::one();
    loop {
        let f_intervals = f_chain.isolate_roots(&epsilon);
        let g_intervals = g_chain.isolate_roots(&epsilon);
        debug_assert_eq!(f_intervals.len(), df);
        debug_assert_eq!(g_intervals.len(), dg);

        let certified = f_intervals
            .iter()
            .enumerate()
            .all(|(i, (f_lo, f_hi))| g_intervals[i].1 < *f_lo && *f_hi < g_intervals[i + 1].0);
        if certified {
            return Some(true);
        }

        let definitely_wrong = f_intervals
            .iter()
            .enumerate()
            .any(|(i, (f_lo, f_hi))| *f_hi < g_intervals[i].0 || g_intervals[i + 1].1 < *f_lo);
        if definitely_wrong {
            return Some(false);
        }
        epsilon /= Q::from_integer(BigInt::from(2));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_first_argument_follows_the_check_interlacing_convention() {
        // 0 ≪ g vacuously for every real-rooted g, in both methods.
        for g in [vec![-1, 1], vec![-1, 0, 1], vec![0, -1, 0, 1]] {
            assert_eq!(check_interlacing_sturm(&[0], &g), Some(true), "{g:?}");
            assert_eq!(
                super::super::check_interlacing(&[0], &g),
                Some(true),
                "{g:?}"
            );
        }
        // A non-real-rooted g is rejected by both.
        assert_eq!(
            super::super::check_interlacing(&[0], &[1, 0, 1]),
            Some(false)
        );
        assert_eq!(check_interlacing_sturm(&[0], &[1, 0, 1]), None);
    }

    #[test]
    fn approximate_roots_of_a_quadratic_are_within_the_documented_error() {
        // (3t - 1)(t - 5): not linear, so both roots are interval midpoints.
        let roots = approximate_real_roots(&[5, -16, 3]).unwrap();
        let tolerance = Q::new(BigInt::from(1), BigInt::from(1_000_000));
        let exact = [
            Q::new(BigInt::from(1), BigInt::from(3)),
            Q::from_integer(BigInt::from(5)),
        ];
        for (root, expected) in roots.iter().zip(&exact) {
            assert!(
                num_traits::Signed::abs(&(root - expected)) < tolerance,
                "{root} vs {expected}"
            );
        }
        // A linear input after removing powers of t has an exact root.
        assert_eq!(
            approximate_real_roots(&[0, 0, -2, 3]).unwrap(),
            vec![
                Q::from_integer(BigInt::from(0)),
                Q::new(BigInt::from(2), BigInt::from(3))
            ]
        );
    }

    #[test]
    fn test_sturm_zero_polynomial_edge_cases() {
        assert!(is_real_rooted_sturm(&[]));
        assert!(is_real_rooted_sturm(&[0, 0, 0]));

        let empty: Vec<BigInt> = Vec::new();
        assert!(is_real_rooted_sturm_bigint_coeffs(&empty));
        assert!(is_real_rooted_sturm_bigint_coeffs(&[
            BigInt::from(0),
            BigInt::from(0),
        ]));
    }

    #[test]
    fn test_has_simple_roots_helpers() {
        assert!(!has_simple_roots(&[1, 2, 1])); // (1+t)^2
        assert!(has_simple_roots(&[1, 0, -1])); // 1 - t^2
        assert!(has_simple_roots(&[7])); // constant
        assert!(!has_simple_roots(&[0])); // zero

        let repeated: Vec<BigInt> = [1_i64, 2, 1].into_iter().map(BigInt::from).collect();
        assert!(!has_simple_roots_bigint_coeffs(&repeated));

        let simple: Vec<BigInt> = [1_i64, 0, -1].into_iter().map(BigInt::from).collect();
        assert!(has_simple_roots_bigint_coeffs(&simple));
    }

    #[test]
    fn test_approximate_real_roots() {
        // (t-1)(t-2) = 2 - 3t + t^2
        let roots = approximate_real_roots(&[2, -3, 1]).unwrap();
        assert_eq!(roots.len(), 2);

        // t^2 + 1: not real-rooted
        assert!(approximate_real_roots(&[1, 0, 1]).is_none());

        // Initial zero coefficients represent a root at zero, not coefficients
        // that may simply be discarded from the returned root set.
        assert_eq!(approximate_real_roots(&[0, 1]), Some(vec![Q::zero()]));
        assert_eq!(approximate_real_roots(&[0, 0, 1]), Some(vec![Q::zero()]));
    }

    #[test]
    fn test_a046741_degree_eight_real_roots_regression() {
        let polynomial = [1, 22, 191, 838, 1982, 2504, 1577, 424, 34];
        assert_eq!(
            approximate_real_roots(&polynomial).map(|roots| roots.len()),
            Some(8)
        );
    }

    #[test]
    fn test_interlacing_sturm() {
        // (t-2) ≪ (t-1)(t-3): roots 1 < 2 < 3
        assert_eq!(check_interlacing_sturm(&[-2, 1], &[3, -4, 1]), Some(true));

        // Wrong order: big first → Some(false) since deg diff is -1
        assert_eq!(check_interlacing_sturm(&[3, -4, 1], &[-2, 1]), Some(false));

        // Same-degree: not interlacing (requires deg diff = 1)
        assert_eq!(
            check_interlacing_sturm(&[-1, 0, 1], &[-4, 0, 1]),
            Some(false)
        );

        // Not real-rooted: t^2+1 → None
        assert_eq!(check_interlacing_sturm(&[1, 1], &[1, 0, 1]), None);

        // The linear root 333333334/1000000000 lies just to the right of
        // 1/3, the smaller root of 2 - 7t + 3t^2. Fixed-width midpoint
        // approximations used to reverse this ordering.
        assert_eq!(
            check_interlacing_sturm(&[-333_333_334, 1_000_000_000], &[2, -7, 3]),
            Some(true)
        );
    }
}
