//! Determinants, fraction-free Bareiss elimination over `Z` and `Z[t]`,
//! characteristic polynomials and closed-walk counts.
//!
//! All results are exact; polynomial entries and characteristic polynomials use
//! ascending degree order.

use crate::vec_poly::strip_trailing_zeros;
use crate::Polynomial;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// Compute the determinant of a BigInt matrix using fraction-free elimination.
///
/// The `0 x 0` matrix has determinant `1`.
///
/// # Panics
///
/// Panics if `mat` is not square (a ragged or rectangular matrix).
pub fn determinant(mat: &[Vec<BigInt>]) -> BigInt {
    assert!(
        is_square_bigint_matrix(mat),
        "determinant requires a square matrix"
    );
    let n = mat.len();
    if n == 0 {
        return BigInt::one();
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut sign_is_negative = false;

    for k in 0..(n - 1) {
        if a[k][k].is_zero() {
            let Some(pivot_row) = ((k + 1)..n).find(|&i| !a[i][k].is_zero()) else {
                return BigInt::zero();
            };
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }

        let pivot = a[k][k].clone();
        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                debug_assert!(
                    (&numerator % &denom).is_zero(),
                    "Bareiss division should be exact"
                );
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    if sign_is_negative {
        -a[n - 1][n - 1].clone()
    } else {
        a[n - 1][n - 1].clone()
    }
}

/// Characteristic polynomial `det(tI - A)` of a square integer matrix.
///
/// Returns the `n + 1` coefficients in ascending degree order (entry `i` is the
/// coefficient of `t^i`), so the result is monic: the last entry is `1`, the
/// entry at index `n - 1` is `-trace(A)`, and entry `0` is `(-1)^n det(A)`.  The
/// `0 x 0` matrix has characteristic polynomial `[1]`.
///
/// The method is the Faddeev--LeVerrier recurrence over `Z`: it uses only
/// matrix products, traces, and divisions by `1, 2, ..., n` that are exact
/// (each quotient is an integer coefficient), so no rational arithmetic or
/// pivoting is involved.  The cost is `n` matrix products, i.e. `O(n^4)`
/// `BigInt` operations.
///
/// # Panics
///
/// Panics if `mat` is not square.
pub fn characteristic_polynomial_bigint(mat: &[Vec<BigInt>]) -> Vec<BigInt> {
    assert!(
        is_square_bigint_matrix(mat),
        "characteristic_polynomial_bigint requires a square matrix"
    );
    let n = mat.len();
    let mut coefficients = vec![BigInt::zero(); n + 1];
    coefficients[n] = BigInt::one();
    // M_1 = I, then M_{k+1} = A M_k + c_{n-k} I with c_{n-k} = -tr(A M_k) / k.
    let mut m: Vec<Vec<BigInt>> = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    if i == j {
                        BigInt::one()
                    } else {
                        BigInt::zero()
                    }
                })
                .collect()
        })
        .collect();
    for k in 1..=n {
        let am = bigint_matrix_product(mat, &m);
        let trace: BigInt = (0..n).map(|i| &am[i][i]).sum();
        let divisor = BigInt::from(k);
        assert!(
            (&trace % &divisor).is_zero(),
            "Faddeev-LeVerrier division should be exact"
        );
        let c = -(trace / divisor);
        m = am;
        for (i, row) in m.iter_mut().enumerate() {
            row[i] += &c;
        }
        coefficients[n - k] = c;
    }
    coefficients
}

/// Closed-walk counts `tr(A^k)` for `k = 0, 1, ..., max_len`.
///
/// For the adjacency matrix of a (di)graph, `tr(A^k)` is the number of closed
/// walks of length `k`; entry `0` is the matrix dimension `n`.  The input may
/// be any square integer matrix.  Returns `max_len + 1` values.
///
/// # Panics
///
/// Panics if `adjacency` is not square.
pub fn closed_walk_counts_bigint(adjacency: &[Vec<BigInt>], max_len: usize) -> Vec<BigInt> {
    assert!(
        is_square_bigint_matrix(adjacency),
        "closed_walk_counts_bigint requires a square matrix"
    );
    let n = adjacency.len();
    let mut counts = Vec::with_capacity(max_len + 1);
    counts.push(BigInt::from(n));
    let mut power = adjacency.to_vec();
    for k in 1..=max_len {
        if k > 1 {
            power = bigint_matrix_product(&power, adjacency);
        }
        counts.push((0..n).map(|i| &power[i][i]).sum());
    }
    counts
}

fn bigint_matrix_product(a: &[Vec<BigInt>], b: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    let n = a.len();
    let mut out = vec![vec![BigInt::zero(); n]; n];
    for i in 0..n {
        for k in 0..n {
            if a[i][k].is_zero() {
                continue;
            }
            for j in 0..n {
                if !b[k][j].is_zero() {
                    out[i][j] += &a[i][k] * &b[k][j];
                }
            }
        }
    }
    out
}

pub(super) fn is_square_bigint_matrix(mat: &[Vec<BigInt>]) -> bool {
    let n = mat.len();
    mat.iter().all(|row| row.len() == n)
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over `BigInt`, without row swaps.
///
/// If every Bareiss pivot is nonzero, the `k`th returned value is the
/// determinant of the leading `(k+1) x (k+1)` principal submatrix.  This is
/// the exact integer version of the fixed-order elimination used in
/// Sylvester's criterion.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact Bareiss division unexpectedly fails.  A nonfinal zero
/// pivot does not imply that the whole matrix is singular; it only means this
/// no-pivot leading-principal computation cannot continue.
pub fn bareiss_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant by no-pivot fraction-free Bareiss elimination.
///
/// This returns `None` in the same cases as
/// [`bareiss_leading_principal_minors_bigint`].  Use [`determinant`] when a
/// row-swapping determinant routine is needed.
pub fn bareiss_determinant_bigint(mat: &[Vec<BigInt>]) -> Option<BigInt> {
    if mat.is_empty() {
        return Some(BigInt::one());
    }
    bareiss_leading_principal_minors_bigint(mat).and_then(|mut minors| minors.pop())
}

/// Compute the leading principal determinants by fraction-free Bareiss
/// elimination over matrices with entries in `Z[t]`.
///
/// The entries are [`Polynomial<BigInt>`] in ascending degree order.  All
/// divisions are checked as exact polynomial divisions over `Z[t]`.  This is
/// useful when a Bézout or Sylvester matrix depends polynomially on a parameter
/// and one wants exact leading principal minors without passing through
/// rational functions.
///
/// Returns `None` if the matrix is not square, a nonfinal required pivot is
/// zero, or an exact polynomial division fails.
pub fn bareiss_leading_principal_minors_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Vec<Polynomial<BigInt>>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = Polynomial::<BigInt>::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = a[i][j].clone() * pivot.clone() - a[i][k].clone() * a[k][j].clone();
                let value = exact_div_polynomial_bigint(&numerator, &denom)?;
                next_entries.push((i, j, value));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}

/// Compute a determinant in `Z[t]` by no-pivot fraction-free Bareiss
/// elimination.
///
/// Use this when the fixed pivot order is part of the certificate.  The return
/// value is `None` under the same conditions as
/// [`bareiss_leading_principal_minors_polynomial_bigint`].
pub fn bareiss_determinant_polynomial_bigint(
    mat: &[Vec<Polynomial<BigInt>>],
) -> Option<Polynomial<BigInt>> {
    if mat.is_empty() {
        return Some(Polynomial::one());
    }
    bareiss_leading_principal_minors_polynomial_bigint(mat).and_then(|mut minors| minors.pop())
}

fn exact_div_polynomial_bigint(
    numerator: &Polynomial<BigInt>,
    denominator: &Polynomial<BigInt>,
) -> Option<Polynomial<BigInt>> {
    if denominator.is_zero() {
        return None;
    }
    if numerator.is_zero() {
        return Some(Polynomial::zero());
    }

    let den = denominator.coeffs();
    let den_degree = den.len() - 1;
    let den_lc = den.last().expect("nonzero denominator");
    let mut rem = numerator.coeffs().to_vec();

    if rem.len() < den.len() {
        return if rem.iter().all(BigInt::is_zero) {
            Some(Polynomial::zero())
        } else {
            None
        };
    }

    let mut quotient = vec![BigInt::zero(); rem.len() - den_degree];
    strip_trailing_zeros(&mut rem);
    while rem.len() >= den.len() && !rem.is_empty() {
        let shift = rem.len() - den.len();
        let rem_lc = rem.last().expect("nonempty remainder").clone();
        if (&rem_lc % den_lc) != BigInt::zero() {
            return None;
        }
        let q = rem_lc / den_lc;
        quotient[shift] = q.clone();
        for (i, c) in den.iter().enumerate() {
            rem[shift + i] -= &q * c;
        }
        strip_trailing_zeros(&mut rem);
    }

    if rem.is_empty() {
        Some(Polynomial::new(quotient))
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::linalg::test_support::{bi, bi_mat};

    fn poly(coeffs: &[i64]) -> Polynomial<BigInt> {
        Polynomial::from_i64_coeffs(coeffs)
    }

    // -----------------------------------------------------------------------
    // Determinant
    // -----------------------------------------------------------------------

    #[test]
    fn test_determinant_empty() {
        assert_eq!(determinant(&[]), BigInt::one());
    }

    #[test]
    fn test_determinant_1x1() {
        assert_eq!(determinant(&bi_mat(&[&[7]])), bi(7));
    }

    #[test]
    fn test_determinant_2x2() {
        // [[1, 2], [3, 4]] -> det = 1*4 - 2*3 = -2
        let m = bi_mat(&[&[1, 2], &[3, 4]]);
        assert_eq!(determinant(&m), bi(-2));
    }

    #[test]
    fn test_determinant_3x3() {
        // [[1, 2, 3], [4, 5, 6], [7, 8, 10]] -> det = 1*(50-48) - 2*(40-42) + 3*(32-35) = 2+4-9 = -3
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6], &[7, 8, 10]]);
        assert_eq!(determinant(&m), bi(-3));
    }

    #[test]
    fn test_determinant_row_pivot_after_first_step() {
        let m = bi_mat(&[&[1, 1, 0], &[1, 1, 1], &[0, 1, 1]]);
        assert_eq!(determinant(&m), bi(-1));
    }

    #[test]
    fn test_determinant_singular() {
        // [[1, 2], [2, 4]] -> det = 0
        let m = bi_mat(&[&[1, 2], &[2, 4]]);
        assert_eq!(determinant(&m), BigInt::zero());
    }

    #[test]
    #[should_panic(expected = "determinant requires a square matrix")]
    fn test_determinant_rejects_nonsquare_matrix() {
        let m = bi_mat(&[&[1, 2, 3], &[4, 5, 6]]);

        let _ = determinant(&m);
    }

    #[test]
    fn test_determinant_identity() {
        let m = bi_mat(&[&[1, 0, 0], &[0, 1, 0], &[0, 0, 1]]);
        assert_eq!(determinant(&m), bi(1));
    }

    #[test]
    fn test_bareiss_leading_principal_minors_bigint() {
        let m = bi_mat(&[&[2, 1], &[1, 2]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(2), bi(3)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(3)));
    }

    #[test]
    fn test_bareiss_bigint_no_pivot_zero_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(determinant(&m), bi(-1));
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
    }

    #[test]
    fn test_bareiss_bigint_final_zero_pivot() {
        let m = bi_mat(&[&[1, 1], &[1, 1]]);
        assert_eq!(
            bareiss_leading_principal_minors_bigint(&m),
            Some(vec![bi(1), bi(0)])
        );
        assert_eq!(bareiss_determinant_bigint(&m), Some(bi(0)));
    }

    #[test]
    fn test_bareiss_polynomial_bigint_2x2() {
        // [[1+z, z], [z, 1+z]] has leading minors 1+z and 1+2z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1])],
            vec![poly(&[0, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 2]))
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_final_zero_pivot() {
        let m = vec![
            vec![poly(&[1, 1]), poly(&[1, 1])],
            vec![poly(&[1, 1]), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), Polynomial::zero()])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(Polynomial::zero())
        );
    }

    #[test]
    fn test_bareiss_polynomial_bigint_nonconstant_exact_division() {
        // The final stage divides by the nonconstant first pivot 1+z.
        let m = vec![
            vec![poly(&[1, 1]), poly(&[0, 1]), Polynomial::zero()],
            vec![poly(&[0, 1]), poly(&[1, 1]), Polynomial::zero()],
            vec![Polynomial::zero(), Polynomial::zero(), poly(&[1, 1])],
        ];
        assert_eq!(
            bareiss_leading_principal_minors_polynomial_bigint(&m),
            Some(vec![poly(&[1, 1]), poly(&[1, 2]), poly(&[1, 3, 2])])
        );
        assert_eq!(
            bareiss_determinant_polynomial_bigint(&m),
            Some(poly(&[1, 3, 2]))
        );
    }

    // Characteristic polynomial and closed walks

    fn lcg_matrix(n: usize, seed: &mut u64, bound: i64) -> Vec<Vec<BigInt>> {
        (0..n)
            .map(|_| {
                (0..n)
                    .map(|_| {
                        *seed = seed
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        bi(((*seed >> 33) as i64) % (2 * bound + 1) - bound)
                    })
                    .collect()
            })
            .collect()
    }

    /// Laplace cofactor expansion along the first row.
    fn cofactor_determinant(mat: &[Vec<BigInt>]) -> BigInt {
        let n = mat.len();
        if n == 0 {
            return BigInt::one();
        }
        let mut total = BigInt::zero();
        for j in 0..n {
            let minor: Vec<Vec<BigInt>> = mat[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter(|&(c, _)| c != j)
                        .map(|(_, v)| v.clone())
                        .collect()
                })
                .collect();
            let term = &mat[0][j] * cofactor_determinant(&minor);
            if j % 2 == 0 {
                total += term;
            } else {
                total -= term;
            }
        }
        total
    }

    fn path_adjacency(n: usize) -> Vec<Vec<BigInt>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| bi(i64::from(i + 1 == j || j + 1 == i)))
                    .collect()
            })
            .collect()
    }

    fn cycle_adjacency(n: usize) -> Vec<Vec<BigInt>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| bi(i64::from((i + 1) % n == j || (j + 1) % n == i)))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn characteristic_polynomial_small_cases() {
        assert_eq!(characteristic_polynomial_bigint(&[]), vec![bi(1)]);
        assert_eq!(
            characteristic_polynomial_bigint(&bi_mat(&[&[5]])),
            vec![bi(-5), bi(1)]
        );
        // det(tI - A) = t^2 - 5t - 2 for [[1,2],[3,4]].
        assert_eq!(
            characteristic_polynomial_bigint(&bi_mat(&[&[1, 2], &[3, 4]])),
            vec![bi(-2), bi(-5), bi(1)]
        );
        // Nilpotent matrix.
        assert_eq!(
            characteristic_polynomial_bigint(&bi_mat(&[&[0, 1, 0], &[0, 0, 1], &[0, 0, 0]])),
            vec![bi(0), bi(0), bi(0), bi(1)]
        );
    }

    #[test]
    fn characteristic_polynomial_matches_cofactor_expansion() {
        let mut seed = 0x9e3779b97f4a7c15u64;
        for n in 1..=5usize {
            for _ in 0..8 {
                let a = lcg_matrix(n, &mut seed, 9);
                let chi = characteristic_polynomial_bigint(&a);
                assert_eq!(chi.len(), n + 1);
                assert_eq!(chi[n], bi(1));
                for t in -3i64..=3 {
                    let shifted: Vec<Vec<BigInt>> = (0..n)
                        .map(|i| {
                            (0..n)
                                .map(|j| {
                                    let diag = if i == j { bi(t) } else { bi(0) };
                                    diag - &a[i][j]
                                })
                                .collect()
                        })
                        .collect();
                    let expected = cofactor_determinant(&shifted);
                    assert_eq!(crate::vec_poly::evaluate(&chi, &bi(t)), expected);
                    assert_eq!(determinant(&shifted), expected);
                }
            }
        }
    }

    #[test]
    fn characteristic_polynomial_of_path_cycle_and_complete_graphs() {
        // Path: p_n = t p_{n-1} - p_{n-2}, i.e. U_n(t/2).
        let mut previous = vec![bi(1)];
        let mut current = vec![bi(0), bi(1)];
        assert_eq!(
            characteristic_polynomial_bigint(&path_adjacency(1)),
            vec![bi(0), bi(1)]
        );
        for n in 2..=9 {
            let next = crate::vec_poly::sub(&crate::vec_poly::shift(&current), &previous);
            assert_eq!(
                characteristic_polynomial_bigint(&path_adjacency(n)),
                next,
                "P_{n}"
            );
            previous = current;
            current = next;
        }
        // Cycle: C_n(t) = 2 T_n(t/2) - 2 = (monic Chebyshev) - 2, via p_n - p_{n-2} - 2.
        let mut paths = vec![vec![bi(1)], vec![bi(0), bi(1)]];
        for n in 2..=9 {
            let next = crate::vec_poly::sub(&crate::vec_poly::shift(&paths[n - 1]), &paths[n - 2]);
            paths.push(next);
        }
        for n in 3..=9 {
            let expected =
                crate::vec_poly::sub(&crate::vec_poly::sub(&paths[n], &paths[n - 2]), &[bi(2)]);
            assert_eq!(
                characteristic_polynomial_bigint(&cycle_adjacency(n)),
                expected,
                "C_{n}"
            );
        }
        // Complete graph: (t - (n-1)) (t + 1)^(n-1).
        for n in 1..=7usize {
            let adjacency: Vec<Vec<BigInt>> = (0..n)
                .map(|i| (0..n).map(|j| bi(i64::from(i != j))).collect())
                .collect();
            let expected = crate::vec_poly::mul(
                &[bi(1 - n as i64), bi(1)],
                &crate::vec_poly::pow(&[bi(1), bi(1)], n - 1),
            );
            assert_eq!(
                characteristic_polynomial_bigint(&adjacency),
                expected,
                "K_{n}"
            );
        }
    }

    #[test]
    fn closed_walk_counts_known_values() {
        // C_5: tr(A^k) = 5, 0, 10, 0, 30, 10 (no triangles; two 5-cycles); k = 0 gives n.
        assert_eq!(
            closed_walk_counts_bigint(&cycle_adjacency(5), 5),
            vec![bi(5), bi(0), bi(10), bi(0), bi(30), bi(10)]
        );
        // K_4: eigenvalues 3, -1, -1, -1.
        let k4: Vec<Vec<BigInt>> = (0..4)
            .map(|i| (0..4).map(|j| bi(i64::from(i != j))).collect())
            .collect();
        let counts = closed_walk_counts_bigint(&k4, 6);
        for (k, count) in counts.iter().enumerate() {
            let expected = 3i64.pow(k as u32) + 3 * (-1i64).pow(k as u32);
            assert_eq!(*count, bi(expected), "k = {k}");
        }
        assert_eq!(closed_walk_counts_bigint(&[], 2), vec![bi(0); 3]);
        assert_eq!(closed_walk_counts_bigint(&k4, 0), vec![bi(4)]);
    }

    #[test]
    fn closed_walks_satisfy_newton_identities_with_characteristic_polynomial() {
        let mut seed = 12345u64;
        for n in 1..=5usize {
            let a = lcg_matrix(n, &mut seed, 4);
            let chi = characteristic_polynomial_bigint(&a);
            let p = closed_walk_counts_bigint(&a, n);
            // Newton: p_k + c_{n-1} p_{k-1} + ... + c_{n-k+1} p_1 + k c_{n-k} = 0.
            for k in 1..=n {
                let mut sum = BigInt::from(k) * &chi[n - k];
                for i in 1..=k {
                    sum += &chi[n - k + i] * &p[i];
                }
                assert!(sum.is_zero(), "n = {n}, k = {k}");
            }
        }
    }

    #[test]
    #[should_panic(expected = "square")]
    fn characteristic_polynomial_rejects_non_square() {
        let _ = characteristic_polynomial_bigint(&bi_mat(&[&[1, 2]]));
    }
}
