//! Real-rootedness entry points: the default PRS root-counting path and the
//! explicit Bézout and Hermite/Newton-sum matrix criteria, plus squarefree parts.
//!
//! The default `is_real_rooted*` functions delegate to [`crate::root_count`];
//! the `*_bezout*` and `*_hermite*` variants are kept for matrix certificates
//! and comparisons.

use super::bezout::{
    bezout_matrix_bigint, bezout_matrix_bigint_coeffs, is_positive_definite_bigint,
    is_positive_semidefinite_bigint,
};
use super::interlacing::check_interlacing_bigint_coeffs;
use super::Q;
use crate::vec_poly;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Signed, ToPrimitive, Zero};

/// Return the squarefree part `p / gcd(p, p')`, with primitive integer
/// coefficients.
///
/// This is useful for boolean real-rootedness checks: a polynomial is weakly
/// real-rooted if and only if its squarefree part is strictly real-rooted.
/// The computation uses a primitive integer pseudo-remainder GCD, avoiding
/// rational coefficient growth. The zero polynomial is returned as `0`.
pub fn squarefree_part_bigint_coeffs(coeffs: &[BigInt]) -> Vec<BigInt> {
    let Some(degree) = vec_poly::degree(coeffs) else {
        return vec![BigInt::zero()];
    };
    let negative_leading_coefficient = coeffs[degree].is_negative();
    let mut squarefree = crate::root_count::squarefree_part_primitive_bigint_coeffs(coeffs);
    if squarefree.is_empty() {
        return vec![BigInt::zero()];
    }
    if negative_leading_coefficient {
        for coefficient in &mut squarefree {
            *coefficient = -coefficient.clone();
        }
    }
    squarefree
}

fn newton_sums_rational_bigint_coeffs(coeffs: &[BigInt], max_power: usize) -> Option<Vec<Q>> {
    let degree = vec_poly::degree(coeffs)?;
    let lc = coeffs[degree].clone();
    if lc.is_zero() {
        return None;
    }

    let lc_q = Q::from_integer(lc);
    let monic_coeffs: Vec<Q> = coeffs
        .iter()
        .take(degree)
        .map(|c| Q::from_integer(c.clone()) / lc_q.clone())
        .collect();

    let mut sums = vec![Q::from_integer(BigInt::zero()); max_power + 1];
    sums[0] = Q::from_integer(BigInt::from(degree as u64));

    for k in 1..=max_power {
        let mut total = Q::from_integer(BigInt::zero());
        let upper = (k - 1).min(degree);
        for i in 1..=upper {
            total += monic_coeffs[degree - i].clone() * sums[k - i].clone();
        }
        if k <= degree {
            total += Q::from_integer(BigInt::from(k as u64)) * monic_coeffs[degree - k].clone();
        }
        sums[k] = -total;
    }

    Some(sums)
}

/// Build the Hermite/Newton-sum matrix for a polynomial.
///
/// For a degree `d` polynomial with roots `alpha_i`, the rational Hermite
/// matrix has entry `s_{i+j}`, where `s_k=sum_i alpha_i^k`.  This function
/// clears denominators by a single positive scalar and returns an integer
/// matrix with the same definiteness.
pub fn hermite_matrix_bigint_coeffs(coeffs: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let degree = vec_poly::degree(coeffs)?;
    if degree == 0 {
        return Some(Vec::new());
    }

    let sums = newton_sums_rational_bigint_coeffs(coeffs, 2 * degree - 2)?;
    let zero = Q::from_integer(BigInt::zero());
    let denominator_lcm = sums
        .iter()
        .take(2 * degree - 1)
        .filter(|s| **s != zero)
        .fold(BigInt::from(1), |acc, s| acc.lcm(s.denom()));
    let scale = Q::from_integer(denominator_lcm);

    let mut matrix = vec![vec![BigInt::zero(); degree]; degree];
    for i in 0..degree {
        for j in 0..degree {
            let scaled = sums[i + j].clone() * scale.clone();
            debug_assert_eq!(scaled.denom(), &BigInt::from(1));
            matrix[i][j] = scaled.to_integer();
        }
    }
    Some(matrix)
}

/// Check real-rootedness by first taking the squarefree part, then checking
/// strict real-rootedness with the Bézout positive-definite criterion.
///
/// This is kept separate from [`is_real_rooted_bezout_bigint_coeffs`] so the
/// semidefinite Bézout implementation remains available for comparison.
pub fn is_real_rooted_bezout_squarefree_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let squarefree = squarefree_part_bigint_coeffs(coeffs);
    is_strictly_real_rooted_bezout_bigint_coeffs(&squarefree)
}

/// Check strict real-rootedness with the Bézout positive-definite criterion,
/// without removing repeated factors first.
pub fn is_strictly_real_rooted_bezout_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(degree) = vec_poly::degree(coeffs) else {
        return true;
    };
    if degree <= 1 {
        return true;
    }
    let derivative = vec_poly::derivative(coeffs);
    check_interlacing_bigint_coeffs(&derivative, coeffs).unwrap_or(false)
}

/// Check real-rootedness by first taking the squarefree part, then checking
/// positive definiteness of the Hermite/Newton-sum matrix.
pub fn is_real_rooted_hermite_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let squarefree = squarefree_part_bigint_coeffs(coeffs);
    is_strictly_real_rooted_hermite_bigint_coeffs(&squarefree)
}

/// Check strict real-rootedness with the Hermite/Newton-sum positive-definite
/// criterion, without removing repeated factors first.
pub fn is_strictly_real_rooted_hermite_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let Some(degree) = vec_poly::degree(coeffs) else {
        return true;
    };
    if degree <= 1 {
        return true;
    }
    let Some(matrix) = hermite_matrix_bigint_coeffs(coeffs) else {
        return false;
    };
    is_positive_definite_bigint(&matrix)
}

/// Explicit squarefree+Bézout real-rootedness check for `i64` coefficients.
pub fn is_real_rooted_bezout_squarefree(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_bezout_squarefree_bigint_coeffs(&big)
}

/// Explicit squarefree+Hermite/Newton real-rootedness check for `i64`
/// coefficients.
pub fn is_real_rooted_hermite(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_real_rooted_hermite_bigint_coeffs(&big)
}

/// Explicit strict Bézout real-rootedness check for `i64` coefficients.
pub fn is_strictly_real_rooted_bezout(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_strictly_real_rooted_bezout_bigint_coeffs(&big)
}

/// Explicit strict Hermite/Newton real-rootedness check for `i64`
/// coefficients.
pub fn is_strictly_real_rooted_hermite(coeffs: &[i64]) -> bool {
    let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
    is_strictly_real_rooted_hermite_bigint_coeffs(&big)
}

/// Check if a polynomial is real-rooted.
///
/// This is the public default exact path.  It delegates to the primitive
/// integer PRS root-counting implementation in [`crate::root_count`].  Use
/// [`is_real_rooted_bezout`], [`is_real_rooted_bezout_squarefree`], or
/// [`is_real_rooted_hermite`] when an explicit matrix-based comparison or
/// certificate is desired.
pub fn is_real_rooted(coeffs: &[i64]) -> bool {
    crate::root_count::is_real_rooted_fast_i64(coeffs)
}

fn is_real_rooted_bezout_i64_impl(coeffs: &[i64]) -> bool {
    let d = match vec_poly::degree(coeffs) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }

    // Compute f'(t) in BigInt to avoid i64 overflow.
    // f' = Σ (k+1) a_{k+1} t^k, degree d-1.
    let mut fp: Vec<BigInt> = Vec::with_capacity(d);
    let mut overflowed = false;
    for k in 0..d {
        let mult = (k as i64) + 1;
        let ak1 = vec_poly::coeff(coeffs, k + 1);
        match mult.checked_mul(ak1) {
            Some(v) => fp.push(BigInt::from(v)),
            None => {
                fp.push(BigInt::from(mult) * BigInt::from(ak1));
                overflowed = true;
            }
        }
    }

    // Ensure leading coefficients have the same sign.
    let lc_f = BigInt::from(vec_poly::coeff(coeffs, d));
    let lc_fp = fp.last().cloned().unwrap_or_else(|| BigInt::from(0));
    let zero_big = BigInt::from(0);
    if (lc_f > zero_big) != (lc_fp > zero_big) {
        for c in fp.iter_mut() {
            *c = -c.clone();
        }
    }

    if overflowed {
        // Use fully BigInt Bézout path
        let f_big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
        let mat = match bezout_matrix_bigint_coeffs(&f_big, &fp) {
            Some(m) => m,
            None => return false,
        };
        is_positive_semidefinite_bigint(&mat)
    } else {
        // Convert back to i64 for the existing path
        let fp_i64: Vec<i64> = fp
            .iter()
            .map(|b| {
                use num_traits::ToPrimitive;
                b.to_i64().unwrap()
            })
            .collect();
        let mat = match bezout_matrix_bigint(coeffs, &fp_i64) {
            Some(m) => m,
            None => return false,
        };
        is_positive_semidefinite_bigint(&mat)
    }
}

/// Check if a polynomial with `BigInt` coefficients is real-rooted.
///
/// This is the public default exact path.  It uses the primitive integer PRS
/// root-counting implementation in [`crate::root_count`], with a specialized
/// one-signed shortcut when applicable.
///
/// Coefficients are given in ascending degree order.
pub fn is_real_rooted_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    crate::root_count::is_real_rooted_fast_bigint_coeffs(coeffs)
}

/// Check if a polynomial with `BigInt` coefficients is real-rooted using the
/// semidefinite Bézout-matrix criterion.
///
/// Exact, but slower than [`is_real_rooted_bigint_coeffs`]; it is kept public for
/// benchmarking and cases where the matrix certificate itself is desired.
pub fn is_real_rooted_bezout_bigint_coeffs(coeffs: &[BigInt]) -> bool {
    let d = match vec_poly::degree(coeffs) {
        Some(d) => d,
        None => return true,
    };
    if d <= 1 {
        return true;
    }

    let mut gcd = BigInt::zero();
    for c in coeffs.iter().take(d + 1) {
        if !c.is_zero() {
            gcd = if gcd.is_zero() {
                c.abs()
            } else {
                gcd.gcd(&c.abs())
            };
        }
    }
    if !gcd.is_zero() {
        let reduced: Vec<BigInt> = coeffs.iter().take(d + 1).map(|c| c / &gcd).collect();
        if reduced.iter().all(|c| c.to_i64().is_some()) {
            let reduced_i64: Vec<i64> = reduced
                .iter()
                .map(|c| c.to_i64().expect("checked above"))
                .collect();
            return is_real_rooted_bezout(&reduced_i64);
        }
    }

    let mut fp: Vec<BigInt> = Vec::with_capacity(d);
    for k in 0..d {
        fp.push(&coeffs[k + 1] * BigInt::from((k + 1) as u64));
    }

    let lc_f = coeffs[d].clone();
    let lc_fp = fp.last().cloned().unwrap_or_else(|| BigInt::from(0));
    let zero_big = BigInt::from(0);
    if (lc_f > zero_big) != (lc_fp > zero_big) {
        for c in &mut fp {
            *c = -c.clone();
        }
    }

    let mat = match bezout_matrix_bigint_coeffs(coeffs, &fp) {
        Some(m) => m,
        None => return false,
    };
    is_positive_semidefinite_bigint(&mat)
}

/// Explicit Bézout-matrix real-rootedness check for `i64` coefficients.
pub fn is_real_rooted_bezout(coeffs: &[i64]) -> bool {
    is_real_rooted_bezout_i64_impl(coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::real_rootedness::qpoly::{poly_exact_div_q, poly_gcd_q, q_poly_to_primitive_bigint};

    #[test]
    fn test_real_rooted() {
        assert!(is_real_rooted(&[1, 2, 1])); // (1+t)^2
        assert!(is_real_rooted(&[1, 11, 11, 1])); // Eulerian A_4(t)
        assert!(!is_real_rooted(&[1, 43, 196, 168, 23, 1])); // counterexample
        assert!(is_real_rooted(&[0])); // zero
        assert!(is_real_rooted(&[5])); // constant
        assert!(is_real_rooted(&[1, 1])); // linear
    }

    #[test]
    fn test_bezout_real_rooted() {
        assert!(is_real_rooted_bezout(&[1, 2, 1])); // (1+t)^2
        assert!(is_real_rooted_bezout(&[1, 3, 3, 1])); // (1+t)^3
        assert!(is_real_rooted_bezout(&[1, 11, 11, 1])); // Eulerian A_4(t)
        assert!(!is_real_rooted_bezout(&[1, 0, 1])); // t^2 + 1
        assert!(is_real_rooted_bezout(&[1])); // constant
        assert!(is_real_rooted_bezout(&[1, 1])); // linear
    }

    #[test]
    fn test_default_real_rooted_bigint_coeffs() {
        let rr: Vec<BigInt> = [1_i64, 3, 3, 1].into_iter().map(BigInt::from).collect();
        assert!(is_real_rooted_bigint_coeffs(&rr));

        let not_rr: Vec<BigInt> = [1_i64, 0, 1].into_iter().map(BigInt::from).collect();
        assert!(!is_real_rooted_bigint_coeffs(&not_rr));

        let huge_scale = BigInt::from(10u64).pow(40);
        let scaled: Vec<BigInt> = [1_i64, 2, 1]
            .into_iter()
            .map(|c| BigInt::from(c) * &huge_scale)
            .collect();
        assert!(is_real_rooted_bigint_coeffs(&scaled));
    }

    #[test]
    fn test_default_real_rootedness_uses_prs_for_mixed_signs() {
        let cases: Vec<(&[i64], bool)> = vec![
            (&[-6, 11, -6, 1], true),  // (t-1)(t-2)(t-3)
            (&[-2, 1, -2, 1], false),  // (t-2)(t^2+1)
            (&[0, 0, 1, -2, 1], true), // t^2(t-1)^2
            (&[1, 0, 2, 0, 1], false), // (t^2+1)^2
        ];

        for (coeffs, expected) in cases {
            let big: Vec<BigInt> = coeffs.iter().copied().map(BigInt::from).collect();
            assert_eq!(
                is_real_rooted_bigint_coeffs(&big),
                crate::root_count::is_real_rooted_fast_bigint_coeffs(&big),
                "default should delegate to PRS backend for {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted(coeffs),
                expected,
                "i64 default disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_bigint_coeffs(&big),
                expected,
                "BigInt default disagrees on {coeffs:?}"
            );
        }
    }

    #[test]
    fn test_squarefree_part_bigint_coeffs() {
        let repeated: Vec<BigInt> = [0_i64, 0, 1, -2, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            squarefree_part_bigint_coeffs(&repeated),
            vec![BigInt::from(0), BigInt::from(-1), BigInt::from(1)]
        );

        let complex_repeated: Vec<BigInt> =
            [1_i64, 0, 2, 0, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            squarefree_part_bigint_coeffs(&complex_repeated),
            vec![BigInt::from(1), BigInt::from(0), BigInt::from(1)]
        );

        let negative_repeated = repeated
            .iter()
            .map(|coefficient| -coefficient)
            .collect::<Vec<_>>();
        assert_eq!(
            squarefree_part_bigint_coeffs(&negative_repeated),
            vec![BigInt::from(0), BigInt::from(1), BigInt::from(-1)]
        );
        assert_eq!(
            squarefree_part_bigint_coeffs(&[BigInt::from(-6)]),
            vec![BigInt::from(-1)]
        );
        assert_eq!(
            squarefree_part_bigint_coeffs(&[BigInt::zero()]),
            vec![BigInt::zero()]
        );
    }

    #[test]
    fn test_primitive_squarefree_matches_rational_reference_on_small_grid() {
        for degree in 0usize..=4 {
            let total = 5usize.pow((degree + 1) as u32);
            for mut mask in 0..total {
                let mut coeffs = Vec::with_capacity(degree + 1);
                for _ in 0..=degree {
                    coeffs.push(BigInt::from((mask % 5) as i64 - 2));
                    mask /= 5;
                }

                let expected = if let Some(actual_degree) = vec_poly::degree(&coeffs) {
                    let trimmed = coeffs
                        .iter()
                        .take(actual_degree + 1)
                        .cloned()
                        .collect::<Vec<_>>();
                    let p = trimmed
                        .iter()
                        .cloned()
                        .map(Q::from_integer)
                        .collect::<Vec<_>>();
                    if actual_degree <= 1 {
                        q_poly_to_primitive_bigint(&p)
                    } else {
                        let derivative = vec_poly::derivative(&trimmed)
                            .into_iter()
                            .map(Q::from_integer)
                            .collect::<Vec<_>>();
                        let gcd = poly_gcd_q(&p, &derivative);
                        q_poly_to_primitive_bigint(&poly_exact_div_q(&p, &gcd))
                    }
                } else {
                    vec![BigInt::zero()]
                };

                assert_eq!(
                    squarefree_part_bigint_coeffs(&coeffs),
                    expected,
                    "coeffs={coeffs:?}"
                );
            }
        }
    }

    #[test]
    fn test_hermite_matrix_bigint_coeffs() {
        // (t - 1)(t - 2) has Newton sums s_0=2, s_1=3, s_2=5.
        let p: Vec<BigInt> = [2_i64, -3, 1].into_iter().map(BigInt::from).collect();
        assert_eq!(
            hermite_matrix_bigint_coeffs(&p),
            Some(vec![
                vec![BigInt::from(2), BigInt::from(3)],
                vec![BigInt::from(3), BigInt::from(5)],
            ])
        );

        // Non-monic example: roots are still 1 and 2, and denominators are
        // cleared by a positive scalar.
        let scaled: Vec<BigInt> = [4_i64, -6, 2].into_iter().map(BigInt::from).collect();
        assert_eq!(
            hermite_matrix_bigint_coeffs(&scaled),
            Some(vec![
                vec![BigInt::from(2), BigInt::from(3)],
                vec![BigInt::from(3), BigInt::from(5)],
            ])
        );
    }

    #[test]
    fn test_squarefree_real_rootedness_paths() {
        let cases: Vec<(&[i64], bool)> = vec![
            (&[1, 2, 1], true),        // repeated real root
            (&[4, -4, 1], true),       // non-monic repeated real root
            (&[0, 0, 1, -2, 1], true), // t^2(t-1)^2
            (&[-2, 3, -1], true),      // negative scalar times (t-1)(t-2)
            (&[1, 0, 1], false),       // t^2+1
            (&[1, 0, 2, 0, 1], false), // (t^2+1)^2
            (&[1, 43, 196, 168, 23, 1], false),
            (&[1, 11, 11, 1], true),
        ];

        for (coeffs, expected) in cases {
            assert_eq!(
                is_real_rooted_bezout_squarefree(coeffs),
                expected,
                "squarefree Bézout disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_hermite(coeffs),
                expected,
                "Hermite/Newton disagrees on {coeffs:?}"
            );
            assert_eq!(
                is_real_rooted_bezout(coeffs),
                expected,
                "old Bézout disagrees on {coeffs:?}"
            );
        }

        assert!(!is_strictly_real_rooted_bezout(&[1, 2, 1]));
        assert!(!is_strictly_real_rooted_hermite(&[1, 2, 1]));
        assert!(is_strictly_real_rooted_bezout(&[2, -3, 1]));
        assert!(is_strictly_real_rooted_hermite(&[2, -3, 1]));
    }
}
