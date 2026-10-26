//! Exact Hurwitz matrices and the Routh--Hurwitz stability criterion.
//!
//! A real polynomial is *Hurwitz stable* if every complex root has strictly
//! negative real part.  This module decides that property exactly over the
//! integers, with no root approximation.
//!
//! Entry points: [`hurwitz_matrix_bigint`], [`hurwitz_leading_minors_bigint`]
//! and [`is_hurwitz_stable_bigint`].
//!
//! # Conventions
//!
//! - **Input** coefficient vectors are ascending (`coeffs[i]` is the
//!   coefficient of `t^i`), as everywhere in Polytool; trailing zeros are
//!   ignored.  Write `n = deg p` and `a_k = coeffs[n - k]`, so `a_0` is the
//!   leading coefficient and `a_n` the constant term (the descending labels
//!   of the classical literature).
//! - **Matrix layout**: the Hurwitz matrix is the `n x n` matrix with
//!   (0-indexed) entry `H[i][j] = a_(2j - i + 1)`, where `a_k = 0` for
//!   `k < 0` or `k > n`:
//!
//!   ```text
//!   | a_1  a_3  a_5  ... |
//!   | a_0  a_2  a_4  ... |
//!   |  0   a_1  a_3  ... |
//!   |  0   a_0  a_2  ... |
//!   |  ...               |
//!   ```
//!
//!   Its leading principal minors are the Hurwitz determinants
//!   `Δ_1 = a_1, Δ_2 = a_1 a_2 - a_0 a_3, ..., Δ_n = a_n Δ_(n-1)`.
//! - **Signs**: the matrix and minors are built from the coefficients exactly
//!   as given; they are not normalized.  Replacing `p` by `-p` multiplies
//!   `Δ_k` by `(-1)^k`.  [`is_hurwitz_stable_bigint`] multiplies by `-1`
//!   first when the leading coefficient is negative, so stability is a
//!   property of the roots alone.
//! - **Criterion** (Hurwitz 1895): for `a_0 > 0`, `p` is Hurwitz stable iff
//!   `Δ_1, ..., Δ_n` are all positive.
//! - **Degenerate inputs**: the zero polynomial has no Hurwitz matrix (`None`)
//!   and is not stable.  A nonzero constant has the empty matrix, no minors,
//!   and is (vacuously) stable because it has no roots.  A linear polynomial
//!   `c_0 + c_1 t` is stable iff `c_0` and `c_1` are nonzero with equal signs.
//! - **Roots at the origin are not removed.**  A zero constant term is a root
//!   at `0`, which is not in the open left half-plane, so such a polynomial
//!   is not stable (and `Δ_n = 0`).  To test the remaining roots, first divide
//!   out the power of `t`, for example with
//!   [`strip_initial_zeros_bigint_coeffs`](crate::strip_initial_zeros_bigint_coeffs)
//!   or [`vec_poly::exact_div_by_t`].
//!
//! Every real-rooted polynomial whose coefficients are all positive (for
//! example an Eulerian polynomial `A_n`, `n >= 1`) is Hurwitz stable, since
//! its roots are negative reals.
//!
//! Related: the Hermite--Biehler split into even and odd parts is
//! [`hermite_biehler_parts`](crate::hermite_biehler_parts); exact determinants
//! come from [`crate::linalg`].
//!
//! ```
//! use num_bigint::BigInt;
//! use polytool::hurwitz::{hurwitz_leading_minors_bigint, is_hurwitz_stable_bigint};
//!
//! let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
//! // 2 + t + t^2 + t^3: Δ_2 = 1·1 - 1·2 < 0, so not stable.
//! let p = big(&[2, 1, 1, 1]);
//! assert_eq!(hurwitz_leading_minors_bigint(&p), Some(big(&[1, -1, -2])));
//! assert!(!is_hurwitz_stable_bigint(&p));
//! // (1 + t)^3 is stable.
//! assert!(is_hurwitz_stable_bigint(&big(&[1, 3, 3, 1])));
//! ```

use crate::linalg::{bareiss_leading_principal_minors_bigint, determinant};
use crate::vec_poly;
use num_bigint::BigInt;
use num_traits::{Signed, Zero};

/// The `n x n` Hurwitz matrix of a polynomial of degree `n`.
///
/// See the [module documentation](self) for the exact layout: with
/// `a_k = coeffs[n - k]`, entry `[i][j]` (0-indexed) is `a_(2j - i + 1)`,
/// so the first row is `a_1, a_3, ...` and the second `a_0, a_2, ...`.
/// Coefficients are used as given (no sign normalization).
///
/// Returns `None` for the zero polynomial and the empty matrix for a nonzero
/// constant.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::hurwitz::hurwitz_matrix_bigint;
///
/// let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
/// // p = 5 + 4t + 3t^2 + 2t^3 + t^4, so (a_0, ..., a_4) = (1, 2, 3, 4, 5).
/// let h = hurwitz_matrix_bigint(&big(&[5, 4, 3, 2, 1])).unwrap();
/// assert_eq!(
///     h,
///     vec![big(&[2, 4, 0, 0]), big(&[1, 3, 5, 0]), big(&[0, 2, 4, 0]), big(&[0, 1, 3, 5])]
/// );
/// ```
pub fn hurwitz_matrix_bigint(coeffs: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let n = vec_poly::degree(coeffs)?;
    // a_k = coeffs[n - k] for 0 <= k <= n, zero otherwise.
    let a = |k: isize| -> BigInt {
        if (0..=n as isize).contains(&k) {
            coeffs[n - k as usize].clone()
        } else {
            BigInt::zero()
        }
    };
    Some(
        (0..n)
            .map(|i| (0..n).map(|j| a(2 * j as isize - i as isize + 1)).collect())
            .collect(),
    )
}

/// The Hurwitz determinants `Δ_1, ..., Δ_n`: the leading principal minors
/// of [`hurwitz_matrix_bigint`].
///
/// Entry `k - 1` of the result is `Δ_k`, the determinant of the leading
/// `k x k` block.  Minors are exact: fraction-free Bareiss elimination is
/// used when all its pivots are nonzero, and otherwise each minor is computed
/// by a row-swapping exact determinant, so zero minors are reported rather
/// than causing a failure.  Coefficients are used as given; for a negative
/// leading coefficient the criterion applies to the minors of `-p`, which are
/// `(-1)^k Δ_k`.
///
/// Returns `None` for the zero polynomial and `Some(vec![])` for a nonzero
/// constant.  For `n >= 1` the last minor satisfies `Δ_n = a_n Δ_(n-1)`
/// (with `Δ_0 = 1`), where `a_n` is the constant term.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::hurwitz::hurwitz_leading_minors_bigint;
///
/// let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
/// // 1 + t + t^2: Δ_1 = 1, Δ_2 = 1.
/// assert_eq!(hurwitz_leading_minors_bigint(&big(&[1, 1, 1])), Some(big(&[1, 1])));
/// // 1 + t^2 + t^4 has Δ_1 = 0; the remaining minors are still computed.
/// assert_eq!(hurwitz_leading_minors_bigint(&big(&[1, 0, 1, 0, 1])), Some(big(&[0, 0, 0, 0])));
/// ```
pub fn hurwitz_leading_minors_bigint(coeffs: &[BigInt]) -> Option<Vec<BigInt>> {
    let matrix = hurwitz_matrix_bigint(coeffs)?;
    if let Some(minors) = bareiss_leading_principal_minors_bigint(&matrix) {
        return Some(minors);
    }
    Some(
        (1..=matrix.len())
            .map(|k| {
                let block: Vec<Vec<BigInt>> =
                    matrix[..k].iter().map(|row| row[..k].to_vec()).collect();
                determinant(&block)
            })
            .collect(),
    )
}

/// Whether every root of the polynomial has strictly negative real part.
///
/// Decided exactly by the Hurwitz criterion after normalizing the leading
/// coefficient to be positive (see the [module documentation](self) for the
/// conventions).  A coefficient that is zero or of the wrong sign rejects the
/// polynomial before any determinant is computed (positivity of all
/// coefficients is necessary for stability).
///
/// The zero polynomial is not stable; a nonzero constant is stable; a zero
/// constant term (a root at the origin) is not stable.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::hurwitz::is_hurwitz_stable_bigint;
///
/// let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
/// assert!(is_hurwitz_stable_bigint(&big(&[1, 1, 1])));        // t^2 + t + 1
/// assert!(is_hurwitz_stable_bigint(&big(&[-1, -1, -1])));     // sign does not matter
/// assert!(!is_hurwitz_stable_bigint(&big(&[1, 0, 1])));       // roots ±i
/// assert!(!is_hurwitz_stable_bigint(&big(&[0, 1, 1])));       // root 0
/// assert!(is_hurwitz_stable_bigint(&big(&[7])));              // no roots
/// assert!(!is_hurwitz_stable_bigint(&big(&[0])));
/// ```
pub fn is_hurwitz_stable_bigint(coeffs: &[BigInt]) -> bool {
    let Some(n) = vec_poly::degree(coeffs) else {
        return false;
    };
    let coeffs = &coeffs[..=n];
    let normalized: Vec<BigInt> = if coeffs[n].is_negative() {
        coeffs.iter().map(|c| -c).collect()
    } else {
        coeffs.to_vec()
    };
    if !normalized.iter().all(Signed::is_positive) {
        return false;
    }
    if n <= 2 {
        // Degrees 0, 1 and 2: positive coefficients are also sufficient.
        return true;
    }
    hurwitz_leading_minors_bigint(&normalized)
        .expect("nonzero polynomial")
        .iter()
        .all(Signed::is_positive)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;
    use num_traits::One;

    type Q = Ratio<BigInt>;

    fn big(p: &[i64]) -> Vec<BigInt> {
        p.iter().map(|&c| BigInt::from(c)).collect()
    }

    /// Independent oracle: the Routh array over the rationals.  For a positive
    /// leading coefficient, `p` is stable iff all `n + 1` first-column entries
    /// are positive; a zero first-column entry means "not stable".
    fn routh_stable(coeffs: &[i64]) -> bool {
        let Some(n) = vec_poly::degree(coeffs) else {
            return false;
        };
        let sign = if coeffs[n] < 0 { -1 } else { 1 };
        // Descending a_0, ..., a_n.
        let a: Vec<Q> = (0..=n)
            .map(|k| Q::from_integer(BigInt::from(sign * coeffs[n - k])))
            .collect();
        let width = n / 2 + 1;
        let row = |start: usize| -> Vec<Q> {
            (0..width)
                .map(|j| a.get(start + 2 * j).cloned().unwrap_or_else(Q::zero))
                .collect()
        };
        let mut prev = row(0);
        let mut cur = row(1);
        let mut first_column = vec![prev[0].clone()];
        for _ in 0..n {
            first_column.push(cur[0].clone());
            if cur[0].is_zero() {
                return false;
            }
            let next: Vec<Q> = (0..width)
                .map(|j| {
                    let p_next = prev.get(j + 1).cloned().unwrap_or_else(Q::zero);
                    let c_next = cur.get(j + 1).cloned().unwrap_or_else(Q::zero);
                    (&cur[0] * &p_next - &prev[0] * &c_next) / &cur[0]
                })
                .collect();
            prev = cur;
            cur = next;
        }
        first_column.iter().all(|x| x > &Q::zero())
    }

    fn from_factors(linear: &[i64], quadratic: &[(i64, i64)]) -> Vec<i64> {
        let mut p = vec![1i64];
        for &r in linear {
            p = vec_poly::mul(&p, &[r, 1]); // t + r, root -r
        }
        for &(b, c) in quadratic {
            p = vec_poly::mul(&p, &[c, b, 1]); // t^2 + b t + c
        }
        p
    }

    #[test]
    fn known_stable_polynomials() {
        for n in 0..=14 {
            let p = vec_poly::pow(&[BigInt::one(), BigInt::one()], n);
            assert!(is_hurwitz_stable_bigint(&p), "(t + 1)^{n}");
            let minors = hurwitz_leading_minors_bigint(&p).unwrap();
            assert_eq!(minors.len(), n);
            assert!(minors.iter().all(Signed::is_positive));
        }
        assert!(is_hurwitz_stable_bigint(&big(&[1, 1, 1])));
        assert!(is_hurwitz_stable_bigint(&big(&[6, 11, 6, 1]))); // (t+1)(t+2)(t+3)
        assert!(is_hurwitz_stable_bigint(&big(&[2, 3]))); // root -2/3
        assert!(is_hurwitz_stable_bigint(&big(&[-2, -3])));
        assert!(is_hurwitz_stable_bigint(&big(&[-1, -3, -3, -1])));
        assert!(is_hurwitz_stable_bigint(&big(&[5, 0, 0])));
        // Complex roots with negative real part.
        assert!(is_hurwitz_stable_bigint(&big(&from_factors(
            &[1, 2],
            &[(1, 5), (2, 2)]
        ))));
    }

    #[test]
    fn known_unstable_polynomials() {
        assert!(!is_hurwitz_stable_bigint(&big(&[2, 1, 1, 1])));
        assert!(!is_hurwitz_stable_bigint(&big(&[1, 1, 1, 1]))); // (t+1)(t^2+1)
        assert!(!is_hurwitz_stable_bigint(&big(&[1, -1, 1]))); // roots in RHP
        assert!(!is_hurwitz_stable_bigint(&big(&[-2, 3]))); // root 2/3
        assert!(!is_hurwitz_stable_bigint(&big(&[0, 1, 1]))); // root 0
        assert!(!is_hurwitz_stable_bigint(&big(&[0, 1])));
        assert!(!is_hurwitz_stable_bigint(&[]));
        assert!(!is_hurwitz_stable_bigint(&big(&[0, 0])));
        // Positive coefficients but a pair of roots in the right half-plane:
        // (t + 5)^2 (t^2 - t + 4) = 100 + 15t + 19t^2 + 9t^3 + t^4.
        let hidden = from_factors(&[5, 5], &[(-1, 4)]);
        assert!(hidden.iter().all(|&c| c > 0), "{hidden:?}");
        assert!(!is_hurwitz_stable_bigint(&big(&hidden)));
        assert!(hurwitz_leading_minors_bigint(&big(&hidden))
            .unwrap()
            .iter()
            .any(|d| !d.is_positive()));
    }

    #[test]
    fn degenerate_inputs() {
        assert_eq!(hurwitz_matrix_bigint(&[]), None);
        assert_eq!(hurwitz_matrix_bigint(&big(&[0, 0])), None);
        assert_eq!(hurwitz_leading_minors_bigint(&big(&[0])), None);
        assert_eq!(hurwitz_matrix_bigint(&big(&[3, 0])), Some(vec![]));
        assert_eq!(hurwitz_leading_minors_bigint(&big(&[3])), Some(vec![]));
        assert_eq!(hurwitz_matrix_bigint(&big(&[3, 4])), Some(vec![big(&[3])]));
        // Trailing zeros are ignored.
        assert_eq!(
            hurwitz_matrix_bigint(&big(&[1, 2, 3, 0, 0])),
            hurwitz_matrix_bigint(&big(&[1, 2, 3]))
        );
    }

    #[test]
    fn minors_are_not_sign_normalized() {
        let p = big(&[1, 3, 3, 1]);
        let neg = vec_poly::neg(&p);
        let minors = hurwitz_leading_minors_bigint(&p).unwrap();
        let neg_minors = hurwitz_leading_minors_bigint(&neg).unwrap();
        for (k, (d, e)) in minors.iter().zip(&neg_minors).enumerate() {
            let expected = if k % 2 == 0 { -d } else { d.clone() };
            assert_eq!(e, &expected);
        }
    }

    #[test]
    fn last_minor_is_constant_term_times_previous() {
        for p in [
            big(&[5, 4, 3, 2, 1]),
            big(&[2, 1, 1, 1]),
            big(&[1, 0, 1, 0, 1]),
            big(&[0, 1, 2, 3]),
            big(&[7, -3, 0, 2, 1, 4]),
        ] {
            let minors = hurwitz_leading_minors_bigint(&p).unwrap();
            let n = minors.len();
            let previous = if n >= 2 {
                minors[n - 2].clone()
            } else {
                BigInt::one()
            };
            assert_eq!(minors[n - 1], &p[0] * previous, "p = {p:?}");
        }
    }

    #[test]
    fn minors_agree_with_direct_determinants_including_zero_pivots() {
        let mut state = 0x5eed_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..300 {
            let len = (next() % 8) as usize + 1;
            let p: Vec<i64> = (0..len).map(|_| (next() % 5) as i64 - 2).collect();
            let Some(h) = hurwitz_matrix_bigint(&big(&p)) else {
                continue;
            };
            let direct: Vec<BigInt> = (1..=h.len())
                .map(|k| {
                    let block: Vec<Vec<BigInt>> =
                        h[..k].iter().map(|row| row[..k].to_vec()).collect();
                    determinant(&block)
                })
                .collect();
            assert_eq!(hurwitz_leading_minors_bigint(&big(&p)), Some(direct));
        }
    }

    #[test]
    fn real_rooted_positive_coefficients_imply_stable() {
        for (n, row) in crate::sequences::eulerian_polynomials_bigint(16)
            .iter()
            .enumerate()
            .skip(1)
        {
            assert!(crate::is_real_rooted_bigint_coeffs(row));
            assert!(is_hurwitz_stable_bigint(row), "A_{n}");
        }
        for roots in [&[1i64, 1, 1][..], &[1, 2, 3, 4, 5, 6, 7], &[3, 9, 27, 81]] {
            assert!(is_hurwitz_stable_bigint(&big(&from_factors(roots, &[]))));
        }
    }

    #[test]
    fn agrees_with_routh_array_oracle() {
        let mut state = 0xdead_beef_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        let mut stable_count = 0;
        for trial in 0..3000 {
            let len = (next() % 8) as usize + 1;
            // Half of the trials use positive coefficients, where both answers
            // occur often; the rest have arbitrary signs.
            let p: Vec<i64> = (0..len)
                .map(|_| {
                    if trial % 2 == 0 {
                        (next() % 6) as i64 + 1
                    } else {
                        (next() % 7) as i64 - 3
                    }
                })
                .collect();
            let expected = routh_stable(&p);
            stable_count += usize::from(expected);
            assert_eq!(is_hurwitz_stable_bigint(&big(&p)), expected, "p = {p:?}");
        }
        assert!(stable_count > 100, "too few stable samples: {stable_count}");
    }

    #[test]
    fn agrees_with_root_construction() {
        let mut state = 0x0123_4567_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for _ in 0..300 {
            let linear: Vec<i64> = (0..next() % 3).map(|_| (next() % 4) as i64 + 1).collect();
            let quadratic: Vec<(i64, i64)> = (0..next() % 3)
                .map(|_| ((next() % 4) as i64 + 1, (next() % 6) as i64 + 1))
                .collect();
            let stable = from_factors(&linear, &quadratic);
            assert!(is_hurwitz_stable_bigint(&big(&stable)), "{stable:?}");
            // A root at a nonnegative real part: t - r (r >= 0) or
            // t^2 - b t + c (b >= 0, c > 0).
            let bad_linear = vec_poly::mul(&stable, &[-((next() % 3) as i64), 1]);
            assert!(
                !is_hurwitz_stable_bigint(&big(&bad_linear)),
                "{bad_linear:?}"
            );
            let b = (next() % 3) as i64;
            let c = (next() % 4) as i64 + 1;
            let bad_quadratic = vec_poly::mul(&stable, &[c, -b, 1]);
            assert!(
                !is_hurwitz_stable_bigint(&big(&bad_quadratic)),
                "{bad_quadratic:?}"
            );
        }
    }
}
