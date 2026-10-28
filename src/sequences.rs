//! Standard polynomial sequences for combinatorial research and testing.
//!
//! Each function returns a vector of coefficient vectors (ascending degree order).
//!
//! # Indexing conventions
//!
//! Every generator takes a `max_n` and returns the rows of one polynomial
//! family.  One rule applies to all of them:
//!
//! * `rows[n]` is the `n`-th member of the family, for every `n` from `0` to
//!   `max_n` inclusive, wherever the family is defined at `n = 0`.  The result
//!   has `max_n + 1` rows.  This is the offset used by most OEIS triangles
//!   (up to the triangle's own first row).
//! * Where a generator cannot or should not define a row `n = 0`, its rustdoc
//!   states the label of `rows[0]` and the number of rows explicitly.
//!
//! | Generator | `rows[0]` | Rows for `max_n` |
//! |---|---|---|
//! | [`eulerian_polynomials`], [`narayana_polynomials`] (and `_bigint`) | `A_0 = N_0 = 1` | `max_n + 1` |
//! | [`type_b_eulerian_polynomials`], [`chebyshev_polynomials_t`], [`chebyshev_polynomials_u`], [`hermite_polynomials`] (and `_bigint`) | member 0 | `max_n + 1` |
//! | [`shifted_legendre_polynomials_bigint`] | `P~_0 = 1` | `max_n + 1` |
//! | [`derangement_descent_polynomials_bigint`], [`fixed_point_descent_polynomials_bigint`] | permutations of the empty set | `max_n + 1` |
//! | [`binomial_eulerian_polynomials_bigint`] | `B~_0 = 1` | `max_n + 1` |
//! | [`refined_eulerian_polynomials_bigint`] | empty (no first letter) | `max_n + 1`; `rows[n][i]` has first letter `i + 1` |
//! | [`s_eulerian_polynomial_bigint`] | (one polynomial, not a row family) | n/a |
//! | [`literature`]: Delannoy, Delannoy-square, Eulerian-square | row 0 (paper normalization) | `max_n + 1` |
//! | [`literature::hoggatt_polynomials_bigint`] | `H_0 = 0` (the paper's rows start at `H_1`; use `rows[1..]`) | `max_n + 1` |
//! | [`deco::joint_histograms_bigint`] and the deco layer generators | size `N = 1` | `max_size` |
//! | [`deco::deco_polynomials_bigint`] | `P_0` | `max_index + 1` |
//!
//! Recurrence-defined families also have endless lazy counterparts
//! ([`eulerian_polynomials_iter`], [`type_b_eulerian_polynomials_iter`],
//! [`chebyshev_polynomials_t_iter`], [`chebyshev_polynomials_u_iter`],
//! [`hermite_polynomials_iter`], and in [`literature`] the Delannoy families):
//! item `n` is member `n`, they are driven by the shared
//! [`crate::recurrence::RecurrenceRowIterator`], and they keep only the last few
//! members.  Closed-form generators (Narayana, shifted Legendre, Hoggatt) and
//! matrix- or array-valued recurrences have no lazy counterpart.
//!
//! The Eulerian row `A_0 = 1` equals `A_1 = 1`; callers that need a family
//! without repeated rows (for example a consecutive-interlacing check) should
//! use `rows[1..]`.  The OEIS triangles A008292 (Eulerian) and A001263
//! (Narayana) start at `n = 1`, so they correspond to `rows[1..]`.
//!
//! These sequences exercise many properties in [`crate::real_rootedness`]:
//!
//! | Sequence | Palindromic | Gamma-positive | Real-rooted | Consecutive interlace |
//! |---|---|---|---|---|
//! | Eulerian A_n | yes | yes | yes | yes |
//! | Narayana N_n | yes | yes | yes | yes |
//! | Type B Eulerian B_n | yes | yes | yes | yes |
//! | Chebyshev T_n, U_n | no | no | yes | yes |
//! | Hermite He_n | no | no | yes | yes |
//! | Shifted Legendre P~_n | no | no | yes (roots in (0, 1)) | yes |

use crate::recurrence::{
    BigRational, BivarPoly, Recurrence, RecurrenceRowIterator, RecurrenceTerm,
};
use crate::vec_poly;
use num_bigint::BigInt;
use std::iter::FusedIterator;

pub mod deco;
pub mod literature;

/// Lazy, endless iterator over the members of a recurrence-defined polynomial
/// family, as exact integer coefficient vectors (ascending degree order).
///
/// Created by the `iter_*` functions of this module and of [`literature`], for
/// example [`eulerian_polynomials_iter`].  The `n`-th item (counting from
/// zero) is member `n` of the family, so
/// `x_iter().take(max_n + 1).collect::<Vec<_>>()` equals the batch generator's
/// `rows` for `max_n`, and the iterator can be paused and continued without
/// recomputing earlier members (`by_ref().take(k)`, then more calls).
///
/// It is driven by the shared [`RecurrenceRowIterator`] kernel: a member is
/// computed only when requested, and only the last few members the recurrence
/// refers to are retained ([`Self::retained_row_count`]).  The recurrences of
/// these named families have integral values, so, unlike the OEIS adapter,
/// no item is fallible.  The iterator is fused and never ends.  The batch
/// generators keep their dedicated integer implementations, which avoid
/// rational arithmetic and are the independent check for these iterators.
#[derive(Debug, Clone)]
pub struct PolynomialFamilyRowIterator {
    rows: RecurrenceRowIterator,
}

impl PolynomialFamilyRowIterator {
    /// Build from sparse recurrence terms `(offset, derivative_order, c)`,
    /// where `c` lists the monomials `(n_degree, t_degree, coefficient)`, and
    /// integer seed rows for the members `0, 1, ...`.
    pub(crate) fn from_sparse_recurrence(terms: &[SparseTerm<'_>], seeds: &[&[i64]]) -> Self {
        let recurrence = Recurrence {
            terms: terms
                .iter()
                .map(|&(offset, derivative_order, monomials)| {
                    RecurrenceTerm::new(offset, derivative_order, bivariate(monomials))
                })
                .collect(),
            denominator: None,
            inhomogeneous: None,
        };
        let seed_rows = seeds
            .iter()
            .map(|row| {
                row.iter()
                    .map(|&c| BigRational::from_integer(BigInt::from(c)))
                    .collect()
            })
            .collect();
        Self {
            rows: RecurrenceRowIterator::new(recurrence, seed_rows, 0),
        }
    }

    /// Number of members yielded so far; also the index of the next member.
    pub fn rows_yielded(&self) -> usize {
        self.rows.rows_yielded()
    }

    /// Number of members currently held in memory: the recurrence's bounded
    /// history plus unyielded seeds, never the whole sequence.
    pub fn retained_row_count(&self) -> usize {
        self.rows.retained_row_count()
    }
}

impl Iterator for PolynomialFamilyRowIterator {
    type Item = Vec<BigInt>;

    fn next(&mut self) -> Option<Self::Item> {
        let row = self
            .rows
            .next()?
            .expect("the named recurrence families are well defined at every index");
        Some(
            row.into_iter()
                .map(|coefficient| {
                    assert!(
                        coefficient.is_integer(),
                        "the named recurrence families have integral coefficients"
                    );
                    coefficient.to_integer()
                })
                .collect(),
        )
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (usize::MAX, None)
    }
}

impl FusedIterator for PolynomialFamilyRowIterator {}

/// `(offset, derivative_order, monomials)`, each monomial `(n_degree, t_degree, coefficient)`.
type SparseTerm<'a> = (usize, usize, &'a [(usize, usize, i64)]);

/// `sum coefficient * n^n_degree * t^t_degree` as a bivariate polynomial.
fn bivariate(monomials: &[(usize, usize, i64)]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let t_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![BigRational::from_integer(BigInt::from(0)); t_size]; n_size];
    for &(n_degree, t_degree, coefficient) in monomials {
        coeffs[n_degree][t_degree] = BigRational::from_integer(BigInt::from(coefficient));
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

/// Lazily generate the Eulerian polynomials `A_0, A_1, A_2, ...` one at a time.
///
/// Item `n` equals `eulerian_polynomials_bigint(max_n)[n]` for every
/// `max_n >= n` (so item 0 and item 1 are both `[1]`).  The family is driven
/// by the shared recurrence kernel through
/// `A_n = (1 + (n-1) t) A_{n-1} + t (1-t) A'_{n-1}` with seed `A_0 = 1`; see
/// [`PolynomialFamilyRowIterator`].
///
/// ```
/// use polytool::sequences::{eulerian_polynomials_bigint, eulerian_polynomials_iter};
///
/// let mut rows = eulerian_polynomials_iter();
/// let first: Vec<_> = rows.by_ref().take(3).collect();
/// let next = rows.next().unwrap(); // A_3, not recomputing A_0..A_2
/// assert_eq!(first, eulerian_polynomials_bigint(2));
/// assert_eq!(next, eulerian_polynomials_bigint(3)[3]);
/// ```
pub fn eulerian_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[
            (1, 0, &[(0, 0, 1), (1, 1, 1), (0, 1, -1)]),
            (1, 1, &[(0, 1, 1), (0, 2, -1)]),
        ],
        &[&[1]],
    )
}

/// Lazily generate the type-B Eulerian polynomials `B_0, B_1, ...`; item `n`
/// equals `type_b_eulerian_polynomials_bigint(max_n)[n]`.  Uses
/// `B_n = (1 + (2n-1) t) B_{n-1} + 2 t (1-t) B'_{n-1}`; see
/// [`PolynomialFamilyRowIterator`].
pub fn type_b_eulerian_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[
            (1, 0, &[(0, 0, 1), (1, 1, 2), (0, 1, -1)]),
            (1, 1, &[(0, 1, 2), (0, 2, -2)]),
        ],
        &[&[1]],
    )
}

/// Lazily generate the Chebyshev polynomials of the first kind `T_0, T_1, ...`;
/// item `n` equals `chebyshev_polynomials_t_bigint(max_n)[n]`.  Uses
/// `T_n = 2t T_{n-1} - T_{n-2}` with seeds `T_0 = 1`, `T_1 = t`; see
/// [`PolynomialFamilyRowIterator`].
pub fn chebyshev_polynomials_t_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[(1, 0, &[(0, 1, 2)]), (2, 0, &[(0, 0, -1)])],
        &[&[1], &[0, 1]],
    )
}

/// Lazily generate the Chebyshev polynomials of the second kind `U_0, U_1, ...`;
/// item `n` equals `chebyshev_polynomials_u_bigint(max_n)[n]`.  Uses
/// `U_n = 2t U_{n-1} - U_{n-2}` with seeds `U_0 = 1`, `U_1 = 2t`; see
/// [`PolynomialFamilyRowIterator`].
pub fn chebyshev_polynomials_u_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[(1, 0, &[(0, 1, 2)]), (2, 0, &[(0, 0, -1)])],
        &[&[1], &[0, 2]],
    )
}

/// Lazily generate the probabilist's Hermite polynomials `He_0, He_1, ...`;
/// item `n` equals `hermite_polynomials_bigint(max_n)[n]`.  Uses
/// `He_n = t He_{n-1} - (n-1) He_{n-2}` with seeds `He_0 = 1`, `He_1 = t`; see
/// [`PolynomialFamilyRowIterator`].
pub fn hermite_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[(1, 0, &[(0, 1, 1)]), (2, 0, &[(1, 0, -1), (0, 0, 1)])],
        &[&[1], &[0, 1]],
    )
}

/// Compute Eulerian polynomials A_0(t), A_1(t), ..., A_{max_n}(t).
///
/// The Eulerian polynomial A_n(t) = sum_{k=0}^{n-1} A(n,k) t^k where A(n,k) is
/// the number of permutations of {1,...,n} with exactly k descents (positions
/// `i` with `w(i) > w(i+1)` in one-line notation).  By convention `A_0 = 1`
/// (the empty permutation), which equals `A_1 = 1`.
///
/// Recurrence: A_n(t) = (1 + (n-1)t) A_{n-1}(t) + t(1-t) A'_{n-1}(t).
///
/// Properties: palindromic, gamma-positive, real-rooted (all roots negative),
/// and consecutive Eulerian polynomials interlace.
///
/// Returns `polys[n]` = A_n(t) for `0 <= n <= max_n` (`max_n + 1` rows), so
/// `polys[0]` = `polys[1]` = `[1]`.  Use `polys[1..]` for the OEIS A008292 rows.
///
/// **Indexing change in 0.4.0:** earlier versions returned `max_n` rows with
/// `polys[i]` = A_{i+1}; use `eulerian_polynomials(n)[1..]` to reproduce that.
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`eulerian_polynomials_bigint`] for larger `max_n`.
pub fn eulerian_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(eulerian_polynomials_bigint(max_n))
}

/// Compute Eulerian polynomials A_0, ..., A_{max_n} with arbitrary-size
/// integer coefficients; `polys[n]` = A_n (see [`eulerian_polynomials`]).
pub fn eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]); // A_0 = 1
    let t_one_minus_t = [BigInt::from(0), BigInt::from(1), BigInt::from(-1)];

    for n in 1..=max_n {
        let prev = &polys[n - 1];
        // A_n = (1 + (n-1)t) A_{n-1} + t(1-t) A'_{n-1}
        let linear = [BigInt::from(1), usize_to_bigint(n - 1)];
        let term1 = vec_poly::mul(&linear, prev);
        let term2 = vec_poly::mul(&t_one_minus_t, &vec_poly::derivative(prev));
        polys.push(vec_poly::add(&term1, &term2));
    }
    polys
}

fn bigint_polys_to_i64(polys: Vec<Vec<BigInt>>) -> Vec<Vec<i64>> {
    polys
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|c| i64::try_from(&c).expect("sequence coefficient overflow"))
                .collect()
        })
        .collect()
}

fn usize_to_bigint(value: usize) -> BigInt {
    BigInt::from(value)
}

fn binomial_row_bigint(n: usize) -> Vec<BigInt> {
    let mut row = vec![BigInt::from(1); n + 1];
    for j in 1..=n {
        row[j] = &row[j - 1] * usize_to_bigint(n - j + 1) / usize_to_bigint(j);
    }
    row
}

fn exact_div_bigint(numerator: BigInt, denominator: BigInt) -> BigInt {
    debug_assert!(denominator != BigInt::from(0));
    debug_assert_eq!(&numerator % &denominator, BigInt::from(0));
    numerator / denominator
}

/// Compute Narayana polynomials N_0, ..., N_{max_n} with arbitrary-size
/// integer coefficients; `polys[n]` = N_n (see [`narayana_polynomials`]).
pub fn narayana_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]); // N_0 = 1 (row sum C_0 = 1)

    for n in 1..=max_n {
        let binom_n = binomial_row_bigint(n);
        let n_big = usize_to_bigint(n);
        let mut coeffs = Vec::with_capacity(n);
        for k in 0..n {
            let j = k + 1;
            coeffs.push(exact_div_bigint(
                &binom_n[j] * &binom_n[j - 1],
                n_big.clone(),
            ));
        }
        vec_poly::strip_trailing_zeros(&mut coeffs);
        polys.push(coeffs);
    }
    polys
}

/// Compute type-B Eulerian polynomials B_0, ..., B_{max_n} (`polys[n]` = B_n)
/// with arbitrary-size integer coefficients.
pub fn type_b_eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]); // B_0 = 1

    let two_t_one_minus_t = [BigInt::from(0), BigInt::from(2), BigInt::from(-2)];

    for n in 1..=max_n {
        let prev = &polys[n - 1];
        // B_n = (1 + (2n-1)t) B_{n-1} + 2t(1-t) B'_{n-1}
        let linear = [BigInt::from(1), usize_to_bigint(2 * n - 1)];
        let term1 = vec_poly::mul(&linear, prev);
        let term2 = vec_poly::mul(&two_t_one_minus_t, &vec_poly::derivative(prev));
        polys.push(vec_poly::add(&term1, &term2));
    }
    polys
}

/// Compute Chebyshev polynomials T_0, ..., T_{max_n} of the first kind (`polys[n]` = T_n)
/// with arbitrary-size integer coefficients.
pub fn chebyshev_polynomials_t_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(1)]);

    for n in 2..=max_n {
        let term1 = vec_poly::shift(&vec_poly::scale(&polys[n - 1], &BigInt::from(2)));
        polys.push(vec_poly::sub(&term1, &polys[n - 2]));
    }
    polys
}

/// Compute Chebyshev polynomials U_0, ..., U_{max_n} of the second kind (`polys[n]` = U_n)
/// with arbitrary-size integer coefficients.
pub fn chebyshev_polynomials_u_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(2)]);

    for n in 2..=max_n {
        let term1 = vec_poly::shift(&vec_poly::scale(&polys[n - 1], &BigInt::from(2)));
        polys.push(vec_poly::sub(&term1, &polys[n - 2]));
    }
    polys
}

/// Compute probabilist's Hermite polynomials He_0, ..., He_{max_n} (`polys[n]` = He_n)
/// with arbitrary-size integer coefficients.
pub fn hermite_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(1)]);

    for n in 2..=max_n {
        let term1 = vec_poly::shift(&polys[n - 1]);
        let term2 = vec_poly::scale(&polys[n - 2], &usize_to_bigint(n - 1));
        polys.push(vec_poly::sub(&term1, &term2));
    }
    polys
}

/// Shifted Legendre polynomials `P~_0(x), ..., P~_{max_n}(x)` on `[0, 1]`
/// (`polys[n]` = `P~_n`), with exact integer coefficients.
///
/// `P~_n(x) = P_n(2x - 1)`, where `P_n` is the classical Legendre polynomial
/// on `[-1, 1]` normalized by `P_n(1) = 1`.  Explicitly,
///
/// ```text
/// P~_n(x) = sum_{k=0}^{n} (-1)^(n+k) C(n, k) C(n+k, k) x^k
///         = sum_{j=0}^{n} C(n, j)^2 x^j (x - 1)^(n-j),
/// ```
///
/// so all coefficients are integers, the leading coefficient is `C(2n, n)`,
/// and `polys[n][k]` has sign `(-1)^(n+k)`.  This is the standard
/// normalization `P~_n(1) = 1`, `P~_n(0) = (-1)^n` (OEIS A063007 up to the
/// sign `(-1)^(n+k)`), with the three-term recurrence
/// `(n+1) P~_{n+1} = (2n+1)(2x-1) P~_n - n P~_{n-1}`, the orthogonality
/// `∫_0^1 P~_m P~_n dx = δ_{mn} / (2n+1)`, and `n` simple roots in `(0, 1)`.
///
/// Returns `max_n + 1` rows; `polys[0] = [1]`, `polys[1] = [-1, 2]`,
/// `polys[2] = [1, -6, 6]`.
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::sequences::shifted_legendre_polynomials_bigint;
///
/// let rows = shifted_legendre_polynomials_bigint(3);
/// let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
/// assert_eq!(rows[3], big(&[-1, 12, -30, 20]));
/// ```
pub fn shifted_legendre_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    (0..=max_n)
        .map(|n| {
            // C(n, k) C(n+k, k) = C(n+k, k) C(n, k), built incrementally:
            // c_k = c_{k-1} * (n - k + 1)(n + k) / k^2, exact at every step.
            let mut row = Vec::with_capacity(n + 1);
            let mut c = BigInt::from(1);
            for k in 0..=n {
                if k > 0 {
                    c = c * usize_to_bigint(n - k + 1) * usize_to_bigint(n + k)
                        / usize_to_bigint(k * k);
                }
                row.push(if (n + k) % 2 == 0 {
                    c.clone()
                } else {
                    -c.clone()
                });
            }
            row
        })
        .collect()
}

/// Compute Narayana polynomials N_0(t), N_1(t), ..., N_{max_n}(t).
///
/// The Narayana polynomial N_n(t) = sum_{k=0}^{n-1} N(n,k+1) t^k where
/// N(n,k) = (1/n) C(n,k) C(n,k-1) is the Narayana number.
///
/// These are the h-vectors of the associahedra and of the type-A Coxeter complex.
///
/// Properties: palindromic, gamma-positive, real-rooted (all roots negative),
/// and row sums are the Catalan numbers.
///
/// By convention `N_0 = 1` (row sum `C_0 = 1`), equal to `N_1 = 1`.
///
/// Returns `polys[n]` = N_n(t) for `0 <= n <= max_n` (`max_n + 1` rows).  Use
/// `polys[1..]` for the OEIS A001263 rows.
///
/// **Indexing change in 0.4.0:** earlier versions returned `max_n` rows with
/// `polys[i]` = N_{i+1}; use `narayana_polynomials(n)[1..]` to reproduce that.
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`narayana_polynomials_bigint`] for larger `max_n`.
pub fn narayana_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(narayana_polynomials_bigint(max_n))
}

/// Compute type-B Eulerian polynomials B_0(t), B_1(t), ..., B_n(t).
///
/// The type-B Eulerian polynomial B_n(t) counts signed permutations of
/// {±1, ..., ±n} by number of descents (in type B sense).
///
/// Recurrence: B_n(t) = (1 + (2n-1)t) B_{n-1}(t) + 2t(1-t) B'_{n-1}(t).
///
/// Properties: palindromic, gamma-positive, real-rooted.
///
/// Returns `polys[n]` = B_n(t) for `0 <= n <= max_n` (`max_n + 1` rows), so
/// `polys[0]` = B_0 = `[1]`.
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`type_b_eulerian_polynomials_bigint`] for larger `max_n`.
pub fn type_b_eulerian_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(type_b_eulerian_polynomials_bigint(max_n))
}

/// Compute Chebyshev polynomials of the first kind T_0(t), T_1(t), ..., T_n(t).
///
/// Recurrence: T_0 = 1, T_1 = t, T_n = 2t T_{n-1} - T_{n-2}.
///
/// Properties: real-rooted (roots are cos(kπ/n)), consecutive T_n interlace.
///
/// Returns `polys[n]` = T_n(t) for `0 <= n <= max_n` (`max_n + 1` rows).
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`chebyshev_polynomials_t_bigint`] for larger `max_n`.
pub fn chebyshev_polynomials_t(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(chebyshev_polynomials_t_bigint(max_n))
}

/// Compute Chebyshev polynomials of the second kind U_0(t), U_1(t), ..., U_n(t).
///
/// Recurrence: U_0 = 1, U_1 = 2t, U_n = 2t U_{n-1} - U_{n-2}.
///
/// Properties: real-rooted, consecutive U_n interlace,
/// and T_n and U_{n-1} interlace (T'_n = n U_{n-1}).
///
/// Returns `polys[n]` = U_n(t) for `0 <= n <= max_n` (`max_n + 1` rows).
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`chebyshev_polynomials_u_bigint`] for larger `max_n`.
pub fn chebyshev_polynomials_u(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(chebyshev_polynomials_u_bigint(max_n))
}

/// Compute probabilist's Hermite polynomials He_0(t), He_1(t), ..., He_n(t).
///
/// Recurrence: He_0 = 1, He_1 = t, He_n = t He_{n-1} - (n-1) He_{n-2}.
///
/// Properties: real-rooted, consecutive He_n interlace.
/// These are the orthogonal polynomials for the standard Gaussian measure.
///
/// Returns `polys[n]` = He_n(t) for `0 <= n <= max_n` (`max_n + 1` rows).
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`; use
/// [`hermite_polynomials_bigint`] for larger `max_n`.
pub fn hermite_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(hermite_polynomials_bigint(max_n))
}

/// Bivariate fixed-point / descent generating polynomials of the symmetric groups.
///
/// For `0 <= n <= max_n`, `rows[n][j][k]` is the number of permutations `w` of
/// `{1, ..., n}` with exactly `j` fixed points (`w(i) = i`) and exactly `k`
/// descents (positions `i` in `1..n` with `w(i) > w(i+1)`, in one-line
/// notation).  Equivalently `rows[n]` is the coefficient array of
///
/// `F_n(a, t) = sum_w a^{fix(w)} t^{des(w)}`,
///
/// with the outer index the power of `a` and the inner (ascending) vector the
/// power of `t`.  Each `rows[n]` has exactly `n + 1` inner vectors
/// (`j = 0..=n`); every inner vector is canonical (trailing zeros removed,
/// zero polynomial `[0]`), so `rows[n][n - 1]` is `[0]` for `n >= 1`.
///
/// `rows[0] = [[1]]` (the empty permutation) and `rows[1] = [[0], [1]]`.
/// Setting `a = 1` gives the Eulerian polynomial `A_n(t)`, and `a = 0` gives
/// [`derangement_descent_polynomials_bigint`].
///
/// Recurrence (`n >= 1`), with `t(1-t) d/dt` acting on `F_{n-1}`:
///
/// `F_n = (1 + (n-1) t) F_{n-1} + t (1-t) dF_{n-1}/dt + (-1)^{n-1} (a-1) (t-a)^{n-1}`.
pub fn fixed_point_descent_polynomials_bigint(max_n: usize) -> Vec<Vec<Vec<BigInt>>> {
    let mut rows: Vec<Vec<Vec<BigInt>>> = Vec::with_capacity(max_n + 1);
    rows.push(vec![vec![BigInt::from(1)]]);
    let t_one_minus_t = [BigInt::from(0), BigInt::from(1), BigInt::from(-1)];

    for n in 1..=max_n {
        let linear = [BigInt::from(1), usize_to_bigint(n - 1)];
        let previous = &rows[n - 1];
        let mut next: Vec<Vec<BigInt>> = Vec::with_capacity(n + 1);
        for j in 0..=n {
            next.push(match previous.get(j) {
                Some(slice) => vec_poly::add(
                    &vec_poly::mul(&linear, slice),
                    &vec_poly::mul(&t_one_minus_t, &vec_poly::derivative(slice)),
                ),
                None => vec![BigInt::from(0)],
            });
        }
        // Boundary term (-1)^{n-1} (a - 1) (t - a)^{n-1}
        //   = sum_r (-1)^{n-1+r} C(n-1, r) (a^{r+1} - a^r) t^{n-1-r}.
        let binomials = binomial_row_bigint(n - 1);
        for (r, binomial) in binomials.iter().enumerate() {
            let signed = if (n - 1 + r) % 2 == 0 {
                binomial.clone()
            } else {
                -binomial
            };
            let t_degree = n - 1 - r;
            add_at(&mut next[r + 1], t_degree, &signed);
            add_at(&mut next[r], t_degree, &-signed);
        }
        for slice in &mut next {
            *slice = vec_poly::trim(slice);
        }
        rows.push(next);
    }
    rows
}

fn add_at(poly: &mut Vec<BigInt>, index: usize, value: &BigInt) {
    if poly.len() <= index {
        poly.resize(index + 1, BigInt::from(0));
    }
    poly[index] += value;
}

/// Descent polynomials of the derangements of `{1, ..., n}`, `n = 0..=max_n`.
///
/// `rows[n] = sum over derangements w of [n] of t^{des(w)}` in ascending
/// degree order, where a derangement has no fixed point (`w(i) != i`) and
/// `des(w)` counts positions `i` in `1..n` with `w(i) > w(i+1)` in one-line
/// notation (the standard descent statistic, **not** the excedance
/// statistic).  Thus `rows[0] = [1]` (empty permutation), `rows[1] = [0]`
/// (no derangement), `rows[2] = [0, 1]`, `rows[3] = [0, 2]`,
/// `rows[4] = [0, 4, 4, 1]`.  This is the `a^0` slice of
/// [`fixed_point_descent_polynomials_bigint`].
///
/// Matches OEIS A219836 (derangements of `[n]` with `k` descents, `n >= 0`),
/// whose rows are padded with zeros to length `n + 1`; the same triangle
/// satisfies `T(n,k) = (k+1) T(n-1,k) + (n-k) T(n-1,k-1) + (-1)^n [k = n-1]`.
///
/// This is *not* OEIS A046739: that triangle counts derangements by
/// excedances (`w(i) > i`), whose polynomials are `t + t^2` for `n = 3` and
/// `t + 7t^2 + t^3` for `n = 4` and which are palindromic after dividing by
/// the lowest power of `t`.  The descent polynomials here are different (for
/// instance `2t` for `n = 3`).
pub fn derangement_descent_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    fixed_point_descent_polynomials_bigint(max_n)
        .into_iter()
        .map(|mut slices| slices.swap_remove(0))
        .collect()
}

/// Binomial Eulerian polynomials `B~_0(t), ..., B~_{max_n}(t)` (Athanasiadis).
///
/// `rows[n] = 1 + t sum_{m=1}^{n} C(n, m) A_m(t)` where `A_m` is the Eulerian
/// polynomial of [`eulerian_polynomials_bigint`] (descents on `S_m`).  Thus
/// `rows[0] = [1]`, `rows[1] = [1, 1]`, `rows[2] = [1, 3, 1]`,
/// `rows[3] = [1, 7, 7, 1]`.  The rows are palindromic of degree `n`, with
/// `rows[n]` of total weight `sum_m C(n, m) m!` (A000522) at `t = 1`.
pub fn binomial_eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let eulerian = eulerian_polynomials_bigint(max_n);
    (0..=max_n)
        .map(|n| {
            let binomials = binomial_row_bigint(n);
            let mut sum = vec![BigInt::from(0)];
            for m in 1..=n {
                sum = vec_poly::add(&sum, &vec_poly::scale(&eulerian[m], &binomials[m]));
            }
            vec_poly::add(&[BigInt::from(1)], &vec_poly::shift(&sum))
        })
        .collect()
}

/// Descent polynomials of permutations of `[n]` refined by the first letter.
///
/// For `0 <= n <= max_n` and `0 <= i < n`, `rows[n][i]` is the polynomial
/// `sum t^{des(w)}` over permutations `w` of `{1, ..., n}` with `w(1) = i + 1`
/// (one-line notation, descents are positions `k < n` with `w(k) > w(k+1)`).
/// Note the 0-based index `i = w(1) - 1`.  `rows[0]` is empty (the empty
/// permutation has no first letter), and for `n >= 1` the polynomials in
/// `rows[n]` sum to the Eulerian polynomial `A_n(t)`.
///
/// Recurrence: with `R(n, j)` the polynomial for first letter `j`
/// (1-based), `R(1, 1) = 1` and `R(n, j) = sum_{k=1}^{n-1} t^{[k < j]} R(n-1, k)`.
pub fn refined_eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<Vec<BigInt>>> {
    let mut rows: Vec<Vec<Vec<BigInt>>> = Vec::with_capacity(max_n + 1);
    rows.push(Vec::new());
    if max_n == 0 {
        return rows;
    }
    rows.push(vec![vec![BigInt::from(1)]]);
    for n in 2..=max_n {
        let previous = &rows[n - 1];
        let mut next = Vec::with_capacity(n);
        for j in 1..=n {
            let mut total = vec![BigInt::from(0)];
            for (k0, poly) in previous.iter().enumerate() {
                let k = k0 + 1;
                total = if k < j {
                    vec_poly::add(&total, &vec_poly::shift(poly))
                } else {
                    vec_poly::add(&total, poly)
                };
            }
            next.push(total);
        }
        rows.push(next);
    }
    rows
}

/// The `s`-Eulerian polynomial of a vector `s = (s_1, ..., s_m)` of positive integers.
///
/// It is `sum t^{asc(e)}` over all sequences `e = (e_1, ..., e_m)` with
/// `0 <= e_i < s_i`, where `asc(e)` counts the indices `i` in `0..m` with
/// `e_i / s_i < e_{i+1} / s_{i+1}` (set `e_0 = 0`, `s_0 = 1`, so the index
/// `i = 0` counts a positive first entry).  For `s = (1, 2, ..., m)` this is
/// the Eulerian polynomial `A_m(t)` (inversion sequences by ascents), and for
/// `s = (1, ..., 1)` it is `1`.
///
/// An empty `s`, or an `s` containing a zero entry (an empty range of
/// sequences), gives the zero polynomial `[0]`.  The coefficients are
/// computed by dynamic programming over the last entry, in
/// `O(sum_i s_i s_{i+1})` polynomial operations, rather than by enumerating
/// the `prod s_i` sequences.
pub fn s_eulerian_polynomial_bigint(s: &[usize]) -> Vec<BigInt> {
    if s.is_empty() || s.contains(&0) {
        return vec![BigInt::from(0)];
    }
    let one = vec![BigInt::from(1)];
    // dp[e] = polynomial generating the prefix ending with value e.
    let mut dp: Vec<Vec<BigInt>> = (0..s[0])
        .map(|e| {
            if e > 0 {
                vec_poly::shift(&one)
            } else {
                one.clone()
            }
        })
        .collect();
    for i in 1..s.len() {
        let mut next = vec![vec![BigInt::from(0)]; s[i]];
        for (e, poly) in dp.iter().enumerate() {
            for (f, slot) in next.iter_mut().enumerate() {
                let term = if e * s[i] < f * s[i - 1] {
                    vec_poly::shift(poly)
                } else {
                    poly.clone()
                };
                *slot = vec_poly::add(slot, &term);
            }
        }
        dp = next;
    }
    dp.iter()
        .fold(vec![BigInt::from(0)], |acc, p| vec_poly::add(&acc, p))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::real_rootedness::*;

    #[test]
    fn test_eulerian_values() {
        let polys = eulerian_polynomials(5);
        assert_eq!(polys.len(), 6);
        assert_eq!(polys[0], vec![1]); // A_0 = 1
        assert_eq!(polys[1], vec![1]); // A_1 = 1
        assert_eq!(polys[2], vec![1, 1]); // A_2 = 1 + t
        assert_eq!(polys[3], vec![1, 4, 1]); // A_3 = 1 + 4t + t^2
        assert_eq!(polys[4], vec![1, 11, 11, 1]); // A_4
        assert_eq!(polys[5], vec![1, 26, 66, 26, 1]); // A_5
        assert_eq!(eulerian_polynomials(0), vec![vec![1]]);
    }

    #[test]
    fn test_eulerian_properties() {
        let polys = eulerian_polynomials(10);
        for p in &polys {
            assert!(is_palindromic(p), "Eulerian not palindromic: {:?}", p);
            assert!(is_real_rooted(p), "Eulerian not real-rooted: {:?}", p);
            assert!(is_gamma_positive(p), "Eulerian not gamma-positive: {:?}", p);
        }
        // Consecutive interlacing (degree diff = 1); A_0 = A_1 has equal degrees.
        for n in 2..polys.len() {
            let result = check_interlacing(&polys[n - 1], &polys[n]);
            assert_eq!(
                result,
                Some(true),
                "A_{} and A_{} don't interlace",
                n,
                n - 1
            );
        }
    }

    #[test]
    fn test_narayana_values() {
        let polys = narayana_polynomials(5);
        assert_eq!(polys.len(), 6);
        assert_eq!(polys[0], vec![1]); // N_0 = 1
        assert_eq!(polys[1], vec![1]); // N_1 = 1
        assert_eq!(polys[2], vec![1, 1]); // N_2 = 1 + t
        assert_eq!(polys[3], vec![1, 3, 1]); // N_3 = 1 + 3t + t^2
        assert_eq!(polys[4], vec![1, 6, 6, 1]); // N_4
        assert_eq!(polys[5], vec![1, 10, 20, 10, 1]); // N_5
    }

    #[test]
    fn test_narayana_properties() {
        let polys = narayana_polynomials(10);
        for p in &polys {
            assert!(is_palindromic(p), "Narayana not palindromic: {:?}", p);
            assert!(is_real_rooted(p), "Narayana not real-rooted: {:?}", p);
            assert!(is_gamma_positive(p), "Narayana not gamma-positive: {:?}", p);
        }
        // Row sums are Catalan numbers C_n: 1, 1, 2, 5, 14, 42, ...
        let catalan = [1i64, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796];
        for (n, p) in polys.iter().enumerate() {
            let sum: i64 = p.iter().sum();
            assert_eq!(sum, catalan[n], "N_{} row sum wrong", n);
        }
    }

    #[test]
    #[should_panic(expected = "sequence coefficient overflow")]
    fn test_narayana_rejects_i64_overflow() {
        let _ = narayana_polynomials(67);
    }

    #[test]
    fn test_type_b_eulerian_values() {
        let polys = type_b_eulerian_polynomials(4);
        assert_eq!(polys[0], vec![1]); // B_0 = 1
        assert_eq!(polys[1], vec![1, 1]); // B_1 = 1 + t
        assert_eq!(polys[2], vec![1, 6, 1]); // B_2 = 1 + 6t + t^2
        assert_eq!(polys[3], vec![1, 23, 23, 1]); // B_3
                                                  // Row sums: (2n)! / 2^n = 1, 1, 3, 15, 105, ...
                                                  // Actually B_n(1) = number of elements in hyperoctahedral group with...
                                                  // Actually sum of type B Eulerian numbers for B_n is 2^n n! / ... let me just verify
                                                  // B_0(1)=1, B_1(1)=2, B_2(1)=8, B_3(1)=48
        assert_eq!(polys[0].iter().sum::<i64>(), 1);
        assert_eq!(polys[1].iter().sum::<i64>(), 2);
        assert_eq!(polys[2].iter().sum::<i64>(), 8);
        assert_eq!(polys[3].iter().sum::<i64>(), 48);
    }

    #[test]
    fn test_type_b_eulerian_properties() {
        let polys = type_b_eulerian_polynomials(8);
        for (i, p) in polys.iter().enumerate() {
            assert!(is_palindromic(p), "B_{} not palindromic", i);
            assert!(is_real_rooted(p), "B_{} not real-rooted", i);
            assert!(is_gamma_positive(p), "B_{} not gamma-positive", i);
        }
    }

    #[test]
    fn test_chebyshev_t_values() {
        let polys = chebyshev_polynomials_t(5);
        assert_eq!(polys[0], vec![1]); // T_0 = 1
        assert_eq!(polys[1], vec![0, 1]); // T_1 = t
        assert_eq!(polys[2], vec![-1, 0, 2]); // T_2 = 2t^2 - 1
        assert_eq!(polys[3], vec![0, -3, 0, 4]); // T_3 = 4t^3 - 3t
        assert_eq!(polys[4], vec![1, 0, -8, 0, 8]); // T_4 = 8t^4 - 8t^2 + 1
    }

    #[test]
    fn test_chebyshev_real_rooted() {
        let t_polys = chebyshev_polynomials_t(12);
        let u_polys = chebyshev_polynomials_u(12);
        for (i, p) in t_polys.iter().enumerate() {
            assert!(is_real_rooted(p), "T_{} not real-rooted", i);
        }
        for (i, p) in u_polys.iter().enumerate() {
            assert!(is_real_rooted(p), "U_{} not real-rooted", i);
        }
    }

    #[test]
    #[should_panic(expected = "sequence coefficient overflow")]
    fn test_chebyshev_t_rejects_i64_overflow() {
        let _ = chebyshev_polynomials_t(64);
    }

    #[test]
    fn test_hermite_values() {
        let polys = hermite_polynomials(5);
        assert_eq!(polys[0], vec![1]); // He_0 = 1
        assert_eq!(polys[1], vec![0, 1]); // He_1 = t
        assert_eq!(polys[2], vec![-1, 0, 1]); // He_2 = t^2 - 1
        assert_eq!(polys[3], vec![0, -3, 0, 1]); // He_3 = t^3 - 3t
        assert_eq!(polys[4], vec![3, 0, -6, 0, 1]); // He_4 = t^4 - 6t^2 + 3
    }

    #[test]
    fn test_hermite_real_rooted() {
        let polys = hermite_polynomials(12);
        for (i, p) in polys.iter().enumerate() {
            assert!(is_real_rooted(p), "He_{} not real-rooted", i);
        }
    }

    #[test]
    fn test_chebyshev_interlacing() {
        // T'_n = n * U_{n-1}, so T_n and U_{n-1} interlace.
        let t_polys = chebyshev_polynomials_t(10);
        let u_polys = chebyshev_polynomials_u(10);
        for n in 2..t_polys.len() {
            let result = check_interlacing(&u_polys[n - 1], &t_polys[n]);
            assert_eq!(
                result,
                Some(true),
                "T_{} and U_{} don't interlace",
                n,
                n - 1
            );
        }
    }

    #[test]
    fn test_hermite_interlacing() {
        // Consecutive Hermite polynomials interlace.
        let polys = hermite_polynomials(10);
        for i in 2..polys.len() {
            let result = check_interlacing(&polys[i - 1], &polys[i]);
            assert_eq!(
                result,
                Some(true),
                "He_{} and He_{} don't interlace",
                i,
                i - 1
            );
        }
    }

    // -- permutation statistics -------------------------------------------

    fn next_permutation(perm: &mut [usize]) -> bool {
        let n = perm.len();
        if n < 2 {
            return false;
        }
        let mut i = n - 1;
        while i > 0 && perm[i - 1] >= perm[i] {
            i -= 1;
        }
        if i == 0 {
            return false;
        }
        let mut j = n - 1;
        while perm[j] <= perm[i - 1] {
            j -= 1;
        }
        perm.swap(i - 1, j);
        perm[i..].reverse();
        true
    }

    /// Brute-force table `table[fix][des]` and `[exc]` counts over all permutations of `[n]`.
    fn brute_tables(n: usize) -> (Vec<Vec<BigInt>>, Vec<Vec<BigInt>>) {
        let mut fix_des = vec![vec![BigInt::from(0); n.max(1)]; n + 1];
        // derangement excedance counts, indexed by number of excedances
        let mut exc = vec![vec![BigInt::from(0); n + 1]];
        let mut perm: Vec<usize> = (1..=n).collect();
        loop {
            let fixed = perm.iter().enumerate().filter(|(i, &v)| v == i + 1).count();
            let des = perm.windows(2).filter(|w| w[0] > w[1]).count();
            fix_des[fixed][des] += 1;
            if fixed == 0 {
                let e = perm.iter().enumerate().filter(|(i, &v)| v > i + 1).count();
                exc[0][e] += 1;
            }
            if !next_permutation(&mut perm) {
                break;
            }
        }
        (fix_des, exc)
    }

    fn trimmed(rows: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
        rows.iter().map(|r| vec_poly::trim(r)).collect()
    }

    #[test]
    fn fixed_point_descent_matches_brute_force_up_to_n8() {
        let rows = fixed_point_descent_polynomials_bigint(8);
        assert_eq!(rows.len(), 9);
        assert_eq!(rows[0], vec![vec![BigInt::from(1)]]);
        for (n, row) in rows.iter().enumerate().skip(1) {
            let (brute, _) = brute_tables(n);
            assert_eq!(row.len(), n + 1);
            assert_eq!(*row, trimmed(&brute), "n = {n}");
        }
    }

    #[test]
    fn fixed_point_descent_specializations() {
        let rows = fixed_point_descent_polynomials_bigint(12);
        let euler = eulerian_polynomials_bigint(12);
        for n in 0..=12 {
            // a = 1 gives the Eulerian polynomial.
            let mut total = vec![BigInt::from(0)];
            for slice in &rows[n] {
                total = vec_poly::add(&total, slice);
            }
            assert_eq!(total, euler[n], "a = 1 at n = {n}");
            // Exactly n-1 fixed points is impossible.
            if n >= 1 {
                assert_eq!(rows[n][n - 1], vec![BigInt::from(0)]);
                assert_eq!(rows[n][n], vec![BigInt::from(1)]);
            }
        }
    }

    fn big_rows(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|r| r.iter().map(|&x| BigInt::from(x)).collect())
            .collect()
    }

    #[test]
    fn derangement_descent_values_and_brute_force() {
        let d = derangement_descent_polynomials_bigint(8);
        assert_eq!(
            d[..6],
            big_rows(&[&[1], &[0], &[0, 1], &[0, 2], &[0, 4, 4, 1], &[0, 8, 24, 12]])[..]
        );
        for (n, row) in d.iter().enumerate().skip(1) {
            let (brute, _) = brute_tables(n);
            assert_eq!(*row, vec_poly::trim(&brute[0]), "n = {n}");
        }
        // Row sums are the derangement numbers 1, 0, 1, 2, 9, 44, 265, 1854, 14833.
        let counts = [1i64, 0, 1, 2, 9, 44, 265, 1854, 14833];
        for (row, count) in d.iter().zip(counts) {
            assert_eq!(row.iter().sum::<BigInt>(), BigInt::from(count));
        }
    }

    #[test]
    fn derangement_descents_differ_from_excedances_a046739() {
        // OEIS A046739 counts derangements by *excedances*; rows n = 1.. are
        // 0; 1; 1,1; 1,7,1; 1,21,21,1 after dividing t^1 (every derangement
        // has at least one excedance).
        let a046739: [&[i64]; 5] = [&[0], &[1], &[1, 1], &[1, 7, 1], &[1, 21, 21, 1]];
        for (n, expected) in (1..=5).zip(a046739) {
            let (_, exc) = brute_tables(n);
            let mut row = vec_poly::trim(&exc[0]);
            if n >= 2 {
                assert_eq!(row[0], BigInt::from(0));
                row.remove(0);
            }
            let expected: Vec<BigInt> = expected.iter().map(|&x| BigInt::from(x)).collect();
            assert_eq!(row, expected, "excedance row n = {n}");
        }
        // The descent distribution is a different triangle.
        let d = derangement_descent_polynomials_bigint(4);
        assert_eq!(d[3], big_rows(&[&[0, 2]])[0]);
        assert_ne!(d[3], big_rows(&[&[0, 1, 1]])[0]);
    }

    #[test]
    fn derangement_descent_matches_oeis_a219836_and_its_recurrence() {
        // OEIS A219836: T(n,k) = derangements of [n] with k descents, n >= 0,
        // 0 <= k <= n, read by rows (data for rows 0..=10).
        const FLAT: [i64; 66] = [
            1, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 4, 4, 1, 0, 0, 8, 24, 12, 0, 0, 0, 16, 104, 120, 24,
            1, 0, 0, 32, 392, 896, 480, 54, 0, 0, 0, 64, 1368, 5544, 5984, 1764, 108, 1, 0, 0, 128,
            4552, 30384, 57640, 34520, 6048, 224, 0, 0, 0, 256, 14680, 153400, 470504, 495320,
            180416, 19936, 448, 1, 0,
        ];
        let rows = derangement_descent_polynomials_bigint(10);
        let mut offset = 0;
        for (n, row) in rows.iter().enumerate() {
            let expected = &FLAT[offset..offset + n + 1];
            offset += n + 1;
            for (k, &value) in expected.iter().enumerate() {
                let got = row.get(k).cloned().unwrap_or_else(|| BigInt::from(0));
                assert_eq!(got, BigInt::from(value), "T({n},{k})");
            }
        }
        assert_eq!(offset, FLAT.len());

        // T(n,k) = (k+1) T(n-1,k) + (n-k) T(n-1,k-1) + (-1)^n [k = n-1], n >= 2.
        let rows = derangement_descent_polynomials_bigint(25);
        let at = |n: usize, k: usize| rows[n].get(k).cloned().unwrap_or_else(|| BigInt::from(0));
        for n in 2..=25usize {
            for k in 0..n {
                let mut expected = BigInt::from(k + 1) * at(n - 1, k);
                if k >= 1 {
                    expected += BigInt::from(n - k) * at(n - 1, k - 1);
                }
                if k == n - 1 {
                    expected += if n % 2 == 0 { 1 } else { -1 };
                }
                assert_eq!(at(n, k), expected, "recurrence at ({n},{k})");
            }
        }
    }

    #[test]
    fn binomial_eulerian_values_and_properties() {
        let rows = binomial_eulerian_polynomials_bigint(7);
        assert_eq!(rows[0], big_rows(&[&[1]])[0]);
        assert_eq!(rows[1], big_rows(&[&[1, 1]])[0]);
        assert_eq!(rows[2], big_rows(&[&[1, 3, 1]])[0]);
        assert_eq!(rows[3], big_rows(&[&[1, 7, 7, 1]])[0]);
        assert_eq!(rows[4], big_rows(&[&[1, 15, 33, 15, 1]])[0]);
        // Totals are sum_m C(n,m) m! = 1, 2, 5, 16, 65, 326, 1957, 13700.
        let totals = [1i64, 2, 5, 16, 65, 326, 1957, 13700];
        for (n, row) in rows.iter().enumerate() {
            assert_eq!(row.iter().sum::<BigInt>(), BigInt::from(totals[n]));
            assert_eq!(row.len(), n + 1);
            let reversed: Vec<BigInt> = row.iter().rev().cloned().collect();
            assert_eq!(*row, reversed);
        }
    }

    #[test]
    fn refined_eulerian_matches_brute_force_and_sums_to_eulerian() {
        let rows = refined_eulerian_polynomials_bigint(8);
        let euler = eulerian_polynomials_bigint(8);
        assert!(rows[0].is_empty());
        for n in 1..=8usize {
            assert_eq!(rows[n].len(), n);
            let mut brute = vec![vec![BigInt::from(0); n]; n];
            let mut perm: Vec<usize> = (1..=n).collect();
            loop {
                let des = perm.windows(2).filter(|w| w[0] > w[1]).count();
                brute[perm[0] - 1][des] += 1;
                if !next_permutation(&mut perm) {
                    break;
                }
            }
            assert_eq!(rows[n], trimmed(&brute), "n = {n}");
            let total = rows[n]
                .iter()
                .fold(vec![BigInt::from(0)], |acc, p| vec_poly::add(&acc, p));
            assert_eq!(total, euler[n]);
        }
    }

    fn s_eulerian_brute(s: &[usize]) -> Vec<BigInt> {
        let mut counts = vec![BigInt::from(0); s.len() + 1];
        let mut e = vec![0usize; s.len()];
        loop {
            let mut asc = usize::from(e[0] > 0);
            for i in 0..s.len() - 1 {
                asc += usize::from(e[i] * s[i + 1] < e[i + 1] * s[i]);
            }
            counts[asc] += 1;
            let mut i = s.len();
            loop {
                if i == 0 {
                    return vec_poly::trim(&counts);
                }
                i -= 1;
                e[i] += 1;
                if e[i] < s[i] {
                    break;
                }
                e[i] = 0;
            }
        }
    }

    #[test]
    fn s_eulerian_matches_brute_force_and_eulerian() {
        let euler = eulerian_polynomials_bigint(8);
        for (m, expected) in euler.iter().enumerate().skip(1) {
            let s: Vec<usize> = (1..=m).collect();
            assert_eq!(s_eulerian_polynomial_bigint(&s), *expected, "s = 1..{m}");
        }
        for s in [
            vec![1usize],
            vec![3],
            vec![1, 1, 1],
            vec![2, 3],
            vec![3, 2, 4],
            vec![2, 2, 2, 2],
            vec![1, 3, 2, 5],
            vec![4, 1, 3, 2],
            vec![2, 4, 3, 5, 3],
        ] {
            assert_eq!(
                s_eulerian_polynomial_bigint(&s),
                s_eulerian_brute(&s),
                "s = {s:?}"
            );
        }
        assert_eq!(s_eulerian_polynomial_bigint(&[]), vec![BigInt::from(0)]);
        assert_eq!(s_eulerian_polynomial_bigint(&[2, 0]), vec![BigInt::from(0)]);
    }

    #[test]
    fn derangement_descent_boundaries() {
        assert_eq!(derangement_descent_polynomials_bigint(0), big_rows(&[&[1]]));
        assert_eq!(
            fixed_point_descent_polynomials_bigint(1),
            vec![big_rows(&[&[1]]), big_rows(&[&[0], &[1]])]
        );
    }

    mod shifted_legendre {
        use super::super::shifted_legendre_polynomials_bigint;
        use crate::numeric::binomial_bigint;
        use crate::vec_poly;
        use num_bigint::BigInt;
        use num_rational::BigRational;
        use num_traits::{One, Zero};

        fn q(n: i64, d: i64) -> BigRational {
            BigRational::new(BigInt::from(n), BigInt::from(d))
        }

        fn to_q(p: &[BigInt]) -> Vec<BigRational> {
            p.iter()
                .map(|c| BigRational::from_integer(c.clone()))
                .collect()
        }

        /// Exact `∫_0^1 p(x) dx = sum_k p_k / (k + 1)`.
        fn integrate_unit_interval(p: &[BigRational]) -> BigRational {
            p.iter()
                .enumerate()
                .map(|(k, c)| c / BigRational::from_integer(BigInt::from(k + 1)))
                .fold(BigRational::zero(), |acc, x| acc + x)
        }

        #[test]
        fn first_rows_and_sign_pattern() {
            let rows = shifted_legendre_polynomials_bigint(4);
            let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
            assert_eq!(rows.len(), 5);
            assert_eq!(rows[0], big(&[1]));
            assert_eq!(rows[1], big(&[-1, 2]));
            assert_eq!(rows[2], big(&[1, -6, 6]));
            assert_eq!(rows[3], big(&[-1, 12, -30, 20]));
            assert_eq!(rows[4], big(&[1, -20, 90, -140, 70]));
            assert_eq!(shifted_legendre_polynomials_bigint(0), vec![big(&[1])]);
        }

        #[test]
        fn closed_form_and_squared_binomial_form() {
            // The squared-binomial form sum_j C(n,j)^2 x^j (x - 1)^(n - j) is
            // the same polynomial, not merely equal up to sign.
            for (n, row) in shifted_legendre_polynomials_bigint(14).iter().enumerate() {
                let mut squared_form = vec![BigInt::zero()];
                for j in 0..=n {
                    let c = binomial_bigint(n, j).pow(2);
                    let term = vec_poly::mul_linear_pow(
                        &vec_poly::shift_by(&[c], j),
                        &-BigInt::one(),
                        &BigInt::one(),
                        n - j,
                    );
                    squared_form = vec_poly::add(&squared_form, &term);
                }
                assert_eq!(row, &squared_form, "n = {n}");
                for (k, c) in row.iter().enumerate() {
                    let magnitude = binomial_bigint(n, k) * binomial_bigint(n + k, k);
                    let expected = if (n + k) % 2 == 0 {
                        magnitude
                    } else {
                        -magnitude
                    };
                    assert_eq!(c, &expected);
                }
                // P~_n(x) = P_n(2x - 1) with the classical recurrence for P_n.
                assert_eq!(vec_poly::evaluate(row, &BigInt::one()), BigInt::one());
                let sign = if n % 2 == 0 { 1 } else { -1 };
                assert_eq!(vec_poly::evaluate(row, &BigInt::zero()), BigInt::from(sign));
                assert_eq!(row[n], binomial_bigint(2 * n, n));
            }
        }

        #[test]
        fn matches_classical_legendre_composed_with_2x_minus_1() {
            // Bonnet recurrence over Q for P_n on [-1, 1].
            let mut classical = vec![
                vec![BigRational::one()],
                to_q(&[BigInt::zero(), BigInt::one()]),
            ];
            for n in 1..10 {
                let a = q(2 * n + 1, n + 1);
                let b = q(n, n + 1);
                let next = vec_poly::sub(
                    &vec_poly::scale(&vec_poly::shift(&classical[n as usize]), &a),
                    &vec_poly::scale(&classical[n as usize - 1], &b),
                );
                classical.push(next);
            }
            let two_x_minus_one = [q(-1, 1), q(2, 1)];
            for (n, row) in shifted_legendre_polynomials_bigint(10).iter().enumerate() {
                assert_eq!(
                    to_q(row),
                    vec_poly::compose(&classical[n], &two_x_minus_one)
                );
            }
        }

        #[test]
        fn three_term_recurrence() {
            let rows = shifted_legendre_polynomials_bigint(20);
            for n in 1..20 {
                // (n+1) P~_{n+1} = (2n+1)(2x-1) P~_n - n P~_{n-1}
                let lhs = vec_poly::scale(&rows[n + 1], &BigInt::from(n + 1));
                let rhs = vec_poly::sub(
                    &vec_poly::mul_linear(
                        &rows[n],
                        &-BigInt::from(2 * n + 1),
                        &BigInt::from(2 * (2 * n + 1)),
                    ),
                    &vec_poly::scale(&rows[n - 1], &BigInt::from(n)),
                );
                assert_eq!(lhs, rhs, "n = {n}");
            }
        }

        #[test]
        fn exact_orthogonality_on_unit_interval() {
            let rows = shifted_legendre_polynomials_bigint(6);
            for m in 0..=6 {
                for n in 0..=6 {
                    let product = vec_poly::mul(&to_q(&rows[m]), &to_q(&rows[n]));
                    let expected = if m == n {
                        q(1, 2 * n as i64 + 1)
                    } else {
                        BigRational::zero()
                    };
                    assert_eq!(
                        integrate_unit_interval(&product),
                        expected,
                        "m = {m}, n = {n}"
                    );
                }
            }
        }

        #[test]
        fn simple_real_roots_in_open_unit_interval() {
            for (n, row) in shifted_legendre_polynomials_bigint(16).iter().enumerate() {
                assert!(crate::is_real_rooted_bigint_coeffs(row), "n = {n}");
                assert_eq!(crate::count_real_roots_uspensky_bigint_coeffs(row), n);
                // x = 1 / (1 + u) maps u in (0, ∞) onto x in (0, 1), and
                // (1 + u)^n P~_n(1 / (1 + u)) = rev_n(P~_n)(1 + u).
                let reversed = vec_poly::reverse(row, n).unwrap();
                let substituted = vec_poly::compose(&reversed, &[BigInt::one(), BigInt::one()]);
                assert_eq!(vec_poly::degree(&substituted), Some(n));
                assert_eq!(
                    crate::count_positive_roots_uspensky_bigint_coeffs(&substituted),
                    n,
                    "n = {n}"
                );
            }
            let rows = shifted_legendre_polynomials_bigint(16);
            for pair in rows.windows(2) {
                assert_eq!(
                    crate::check_interlacing_bigint_coeffs(&pair[0], &pair[1]),
                    Some(true)
                );
            }
        }
    }
}
