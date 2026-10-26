//! Bézout matrices and Bézoutians, and the exact definiteness checks used on them.
//!
//! `bezout_matrix*` builds the degree-difference-one Bézout matrix and
//! `bezoutian_bigint_coeffs` the general Bézoutian; the interlacing and
//! real-rootedness criteria built on them live in the sibling `interlacing` and
//! `real_rooted` modules.

use crate::vec_poly;
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

/// Compute the Bézout matrix B(f, g) where deg(f) = deg(g) + 1.
///
/// This is the degree-difference-one Bézoutian used in the interlacing
/// criterion of Kummer--Naldi--Plaumann, Thm. 2.13, after Krein--Naimark,
/// §2.2.  Same-degree interlacing is handled by [`check_weak_interlacing`]
/// via a degree-extension reduction; see Fisk, §9.21, Cor. 9.145 for the
/// same-degree positive definite form.
///
/// The (i, j) entry (0-indexed) is the coefficient of x^i y^j in
///
/// ```text
/// (f(x)g(y) - f(y)g(x)) / (x - y)
/// ```
///
/// Returns a d x d matrix (as `i64`) where d = deg(f), or `None` if the degree
/// constraint is not satisfied or an exact entry does not fit in `i64`.
///
/// The division is exact because h(x,y) = f(x)g(y) - f(y)g(x) vanishes on x = y.
/// Expanding: for each pair (k, l) with k > l, the contribution is
///
/// ```text
/// (f_k g_l - f_l g_k) * sum_{m=0}^{k-l-1} x^{l+m} y^{k-1-m}
/// ```
///
/// which follows from factoring x^k y^l - x^l y^k = x^l y^l (x^{k-l} - y^{k-l})
/// and using the geometric sum (x^n - y^n)/(x - y) = sum x^i y^{n-1-i}.
///
/// [`check_weak_interlacing`]: super::check_weak_interlacing
pub fn bezout_matrix(f: &[i64], g: &[i64]) -> Option<Vec<Vec<i64>>> {
    bezout_matrix_bigint(f, g)?
        .into_iter()
        .map(|row| row.into_iter().map(|entry| entry.to_i64()).collect())
        .collect()
}

/// Compute the Bézout matrix with BigInt entries (no overflow).
pub(super) fn bezout_matrix_bigint(f: &[i64], g: &[i64]) -> Option<Vec<Vec<BigInt>>> {
    let f_big: Vec<BigInt> = f.iter().map(|&c| BigInt::from(c)).collect();
    let g_big: Vec<BigInt> = g.iter().map(|&c| BigInt::from(c)).collect();
    bezout_matrix_bigint_coeffs(&f_big, &g_big)
}

/// Compute the Bézout matrix from BigInt coefficient vectors.
///
/// The coefficient vectors are in ascending degree order.  This is the exact
/// no-overflow variant of [`bezout_matrix`], useful for benchmarking and for
/// downstream linear-algebra experiments.
pub fn bezout_matrix_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let zero_big = BigInt::from(0);
    let df = {
        let mut d = None;
        for i in (0..f.len()).rev() {
            if f[i] != zero_big {
                d = Some(i);
                break;
            }
        }
        d
    }?;
    let dg = {
        let mut d = None;
        for i in (0..g.len()).rev() {
            if g[i] != zero_big {
                d = Some(i);
                break;
            }
        }
        d
    }?;
    if df != dg + 1 {
        return None;
    }
    let d = df;
    let zero = BigInt::from(0);
    let mut b = vec![vec![zero.clone(); d]; d];

    for k in 0..=df {
        for l in 0..k {
            let fk = if k < f.len() { &f[k] } else { &zero };
            let fl = if l < f.len() { &f[l] } else { &zero };
            let gk = if k < g.len() { &g[k] } else { &zero };
            let gl = if l < g.len() { &g[l] } else { &zero };
            let c = fk * gl - fl * gk;
            if c == zero {
                continue;
            }
            for m in 0..=(k - l - 1) {
                let xi = l + m;
                let yj = k - 1 - m;
                if xi < d && yj < d {
                    b[xi][yj] = &b[xi][yj] + &c;
                }
            }
        }
    }
    Some(b)
}

/// The Bézoutian of two polynomials of arbitrary degrees.
///
/// Let `N = max(deg f, deg g)`, where the zero polynomial does not count
/// (so `N = 0` when both are constant or zero).  The result is the `N x N`
/// matrix `B` whose (0-indexed) entry `B[i][j]` is the coefficient of
/// `x^i y^j` in
///
/// ```text
/// (f(x) g(y) - f(y) g(x)) / (x - y).
/// ```
///
/// This is the same orientation as [`bezout_matrix_bigint_coeffs`], which
/// is the special case `deg f = deg g + 1` (and returns `None` for every other
/// degree pair); here every degree pair is accepted.  Coefficient vectors are
/// ascending and may carry trailing zeros.
///
/// Properties (all exact):
///
/// - `B` is symmetric, bilinear in `(f, g)`, and antisymmetric under swapping
///   the arguments: `B(g, f) = -B(f, g)`.
/// - `rank B = N - deg gcd(f, g)` (with `deg gcd = N` when one argument is
///   zero), so `B` is nonsingular iff `f` and `g` are coprime.
/// - The zero polynomial gives the `N x N` zero matrix, and two constants give
///   the empty matrix.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::bezoutian_bigint_coeffs;
///
/// let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
/// // f = t^2 - 1, g = t: (f(x)y - f(y)x) / (x - y) = xy + 1.
/// let b = bezoutian_bigint_coeffs(&big(&[-1, 0, 1]), &big(&[0, 1]));
/// assert_eq!(b, vec![big(&[1, 0]), big(&[0, 1])]);
/// // Same degree: f = t^2, g = t^2 + t gives xy.
/// let b = bezoutian_bigint_coeffs(&big(&[0, 0, 1]), &big(&[0, 1, 1]));
/// assert_eq!(b, vec![big(&[0, 0]), big(&[0, 1])]);
/// ```
pub fn bezoutian_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Vec<Vec<BigInt>> {
    let n = vec_poly::degree(f)
        .unwrap_or(0)
        .max(vec_poly::degree(g).unwrap_or(0));
    let zero = BigInt::zero();
    let mut b = vec![vec![zero.clone(); n]; n];
    let coeff = |p: &[BigInt], k: usize| -> BigInt { p.get(k).cloned().unwrap_or_default() };
    // For k > l the pair contributes (f_k g_l - f_l g_k) (x^k y^l - x^l y^k),
    // and (x^k y^l - x^l y^k) / (x - y) = sum_{m < k - l} x^(l + m) y^(k - 1 - m).
    for k in 1..=n {
        for l in 0..k {
            let c = coeff(f, k) * coeff(g, l) - coeff(f, l) * coeff(g, k);
            if c.is_zero() {
                continue;
            }
            for m in 0..(k - l) {
                b[l + m][k - 1 - m] += &c;
            }
        }
    }
    b
}

/// Check positive definiteness of a symmetric BigInt matrix via Gaussian elimination.
pub(super) fn is_positive_definite_bigint(mat: &[Vec<BigInt>]) -> bool {
    crate::linalg::is_positive_definite(mat)
}

pub(super) fn is_positive_semidefinite_bigint(mat: &[Vec<BigInt>]) -> bool {
    crate::linalg::is_positive_semidefinite(mat)
}

#[cfg(test)]
mod tests {
    use super::*;

    // -- Bézout matrix tests --

    #[test]
    fn test_bezout_matrix_simple() {
        // f = (t-1)(t-3) = 3 - 4t + t^2, g = t - 2
        // B should be 2x2... wait deg(f)=2, deg(g)=1, so d=2.
        let b = bezout_matrix(&[3, -4, 1], &[-2, 1]).unwrap();
        assert_eq!(b.len(), 2);
        // Should be positive definite since g interlaces f (roots: 1 < 2 < 3).
        let b_big = bezout_matrix_bigint(&[3, -4, 1], &[-2, 1]).unwrap();
        assert!(is_positive_definite_bigint(&b_big));
    }

    #[test]
    fn test_bezout_matrix_checked_narrowing() {
        let max = i64::MAX;
        let min = i64::MIN;
        assert!(bezout_matrix(&[max, min, max], &[min, max]).is_none());
    }
}

#[cfg(test)]
mod bezoutian_tests {
    use super::*;
    use crate::Polynomial;
    use num_rational::Ratio;

    type Q = Ratio<BigInt>;

    fn big(p: &[i64]) -> Vec<BigInt> {
        p.iter().map(|&c| BigInt::from(c)).collect()
    }

    /// Deterministic xorshift generator for reproducible random inputs.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn poly(&mut self, max_len: usize, bound: i64) -> Vec<i64> {
            let len = (self.next() % (max_len as u64 + 1)) as usize;
            (0..len)
                .map(|_| (self.next() % (2 * bound as u64 + 1)) as i64 - bound)
                .collect()
        }
    }

    fn rank_over_q(mat: &[Vec<BigInt>]) -> usize {
        let mut a: Vec<Vec<Q>> = mat
            .iter()
            .map(|row| row.iter().map(|x| Q::from_integer(x.clone())).collect())
            .collect();
        let rows = a.len();
        let cols = a.first().map_or(0, Vec::len);
        let mut rank = 0;
        for col in 0..cols {
            let Some(pivot) = (rank..rows).find(|&r| !a[r][col].is_zero()) else {
                continue;
            };
            a.swap(rank, pivot);
            for r in 0..rows {
                if r != rank && !a[r][col].is_zero() {
                    let factor = &a[r][col] / &a[rank][col];
                    let pivot_row = a[rank].clone();
                    for (entry, p) in a[r].iter_mut().zip(&pivot_row).skip(col) {
                        *entry -= &factor * p;
                    }
                }
            }
            rank += 1;
        }
        rank
    }

    fn gcd_degree(f: &[i64], g: &[i64]) -> Option<usize> {
        let to_q = |p: &[i64]| {
            Polynomial::new(
                p.iter()
                    .map(|&c| Q::from_integer(BigInt::from(c)))
                    .collect::<Vec<_>>(),
            )
        };
        to_q(f).gcd(&to_q(g)).degree()
    }

    /// Evaluate `sum B[i][j] x^i y^j`.
    fn eval_form(b: &[Vec<BigInt>], x: &BigInt, y: &BigInt) -> BigInt {
        let mut total = BigInt::zero();
        for (i, row) in b.iter().enumerate() {
            for (j, entry) in row.iter().enumerate() {
                total += entry * x.pow(i as u32) * y.pow(j as u32);
            }
        }
        total
    }

    #[test]
    fn bezoutian_agrees_with_degree_difference_one_bezout_matrix() {
        let mut rng = Rng(0xbe20_u64);
        let mut checked = 0;
        while checked < 200 {
            let f = rng.poly(7, 6);
            let g = rng.poly(6, 6);
            let (Some(df), Some(dg)) = (vec_poly::degree(&f), vec_poly::degree(&g)) else {
                continue;
            };
            if df != dg + 1 {
                continue;
            }
            assert_eq!(
                Some(bezoutian_bigint_coeffs(&big(&f), &big(&g))),
                bezout_matrix_bigint_coeffs(&big(&f), &big(&g)),
                "f = {f:?}, g = {g:?}"
            );
            checked += 1;
        }
    }

    #[test]
    fn bezoutian_matches_defining_identity_symmetry_and_gcd_rank() {
        let mut rng = Rng(424_242);
        for _ in 0..300 {
            let f = rng.poly(6, 5);
            let g = rng.poly(6, 5);
            // Force common factors in a share of the cases.
            if rng.next() % 3 == 0 {
                let common = rng.poly(4, 3);
                check_pair(&vec_poly::mul(&f, &common), &vec_poly::mul(&g, &common));
            }
            check_pair(&f, &g);
        }
    }

    fn check_pair(f: &[i64], g: &[i64]) {
        let (fb, gb) = (big(f), big(g));
        let b = bezoutian_bigint_coeffs(&fb, &gb);
        let n = vec_poly::degree(f)
            .unwrap_or(0)
            .max(vec_poly::degree(g).unwrap_or(0));
        assert_eq!(b.len(), n);
        assert!(b.iter().all(|row| row.len() == n));
        for (i, row) in b.iter().enumerate() {
            for (j, entry) in row.iter().enumerate() {
                assert_eq!(entry, &b[j][i], "symmetry, f = {f:?}, g = {g:?}");
            }
        }
        // Antisymmetry in the arguments.
        let swapped = bezoutian_bigint_coeffs(&gb, &fb);
        for (row, srow) in b.iter().zip(&swapped) {
            for (x, y) in row.iter().zip(srow) {
                assert_eq!(x, &-y);
            }
        }
        // (x - y) B(x, y) = f(x) g(y) - f(y) g(x) at integer points.
        for x in -3i64..=3 {
            for y in -3i64..=3 {
                let (xb, yb) = (BigInt::from(x), BigInt::from(y));
                let lhs = (&xb - &yb) * eval_form(&b, &xb, &yb);
                let rhs = vec_poly::evaluate(&fb, &xb) * vec_poly::evaluate(&gb, &yb)
                    - vec_poly::evaluate(&fb, &yb) * vec_poly::evaluate(&gb, &xb);
                assert_eq!(lhs, rhs, "f = {f:?}, g = {g:?}, x = {x}, y = {y}");
            }
        }
        // rank = N - deg gcd(f, g); if f = g = 0 the matrix is empty.
        let expected_rank = gcd_degree(f, g).map_or(0, |d| n - d);
        assert_eq!(rank_over_q(&b), expected_rank, "f = {f:?}, g = {g:?}");
    }

    #[test]
    fn bezoutian_edge_cases() {
        assert!(bezoutian_bigint_coeffs(&[], &[]).is_empty());
        assert!(bezoutian_bigint_coeffs(&big(&[3]), &big(&[0, 0])).is_empty());
        let zero_against_cubic = bezoutian_bigint_coeffs(&big(&[0]), &big(&[1, 2, 3, 4]));
        assert_eq!(zero_against_cubic, vec![big(&[0, 0, 0]); 3]);
        // (2 y^2 - 2 x^2) / (x - y) = -2 (x + y)
        let b = bezoutian_bigint_coeffs(&big(&[2]), &big(&[0, 0, 1]));
        assert_eq!(b, vec![big(&[0, -2]), big(&[-2, 0])]);
    }
}
