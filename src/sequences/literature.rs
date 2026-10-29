//! Exact polynomial families for reproducible literature-conjecture checks.
//!
//! Coefficients are in ascending degree order. No real-rootedness assertion
//! for a conjectural family is part of these generation APIs.

use super::PolynomialFamilyRowIterator;
use crate::polynomial::Polynomial;
use crate::recurrence::{BivarPoly, Recurrence, RecurrenceTerm};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

/// Delannoy triangle rows D_0(x), ..., D_max_n(x); `rows[n]` = D_n.
///
/// D(n,k) counts paths to (n-k,k) with steps (1,0), (0,1), (1,1).
/// D_0=1, D_1=1+x, D_n=(1+x)D_(n-1)+x D_(n-2).
pub fn delannoy_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut rows = vec![vec![BigInt::one()]];
    for n in 1..=max_n {
        let mut row = vec![BigInt::zero(); n + 1];
        add_linear_multiple(&mut row, &rows[n - 1], 1, 1);
        if n >= 2 {
            add_linear_multiple(&mut row, &rows[n - 2], 0, 1);
        }
        rows.push(row);
    }
    rows
}

/// Lazily generate the Delannoy polynomials `D_0, D_1, ...`; item `n` equals
/// `delannoy_polynomials_bigint(max_n)[n]`.  Uses
/// `D_n = (1 + x) D_{n-1} + x D_{n-2}` with seeds `D_0 = 1`, `D_1 = 1 + x`; see
/// [`PolynomialFamilyRowIterator`].
pub fn delannoy_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[(1, 0, &[(0, 0, 1), (0, 1, 1)]), (2, 0, &[(0, 1, 1)])],
        &[&[1], &[1, 1]],
    )
}

/// Lazily generate the rows `G_0, G_1, ...` of
/// [`delannoy_square_polynomials_bigint`]; item `n` equals its `rows[n]`.
/// Uses the four-term recurrence of that function with its four seed rows; see
/// [`PolynomialFamilyRowIterator`].
pub fn delannoy_square_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[
            (1, 0, &[(0, 0, 3), (0, 1, 1)]),
            (2, 0, &[(0, 0, -1), (0, 1, 1)]),
            (3, 0, &[(0, 0, -1), (0, 1, 1)]),
            (4, 0, &[(0, 1, 1)]),
        ],
        &[&[1], &[2, 1], &[5, 6, 1], &[12, 25, 10, 1]],
    )
}

/// Row-generating polynomials of the square of the Delannoy triangle, n=0..max_n.
///
/// Uses the four-term recurrence from Section 4 of Mao--Wang,
/// "The Narayana transformation", arXiv:2607.01572v1:
/// G_n=(x+3)G_(n-1)+(x-1)G_(n-2)+(x-1)G_(n-3)+xG_(n-4).
/// Initial rows: 1; 2+x; 5+6x+x^2; 12+25x+10x^2+x^3.  `rows[n]` = G_n, `max_n + 1` rows.
pub fn delannoy_square_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let seeds: [&[i64]; 4] = [&[1], &[2, 1], &[5, 6, 1], &[12, 25, 10, 1]];
    let mut rows: Vec<Vec<BigInt>> = seeds
        .iter()
        .take(max_n.saturating_add(1))
        .map(|r| r.iter().map(|&x| BigInt::from(x)).collect())
        .collect();
    for n in 4..=max_n {
        let mut row = vec![BigInt::zero(); n + 1];
        add_linear_multiple(&mut row, &rows[n - 1], 3, 1);
        add_linear_multiple(&mut row, &rows[n - 2], -1, 1);
        add_linear_multiple(&mut row, &rows[n - 3], -1, 1);
        add_linear_multiple(&mut row, &rows[n - 4], 0, 1);
        rows.push(row);
    }
    rows
}

/// Eulerian matrix-square rows, n=0..max_n, in Mao--Wang's normalization.
///
/// A(0,0)=1; for n>=1, A(n,k) counts permutations of n with k-1 descents.
/// Thus the nonzero-row basis is x times the usual Eulerian polynomial.
/// This is triangular matrix multiplication, NOT coefficientwise squaring.
/// The output retains its structural zero constant term for n>=1.
/// `rows[n]` is the paper's row n, with `max_n + 1` rows.
pub fn eulerian_square_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut triangle = vec![vec![BigInt::one()]];
    for row in super::eulerian_polynomials_bigint(max_n)
        .into_iter()
        .skip(1)
    {
        let mut shifted = vec![BigInt::zero()];
        shifted.extend(row);
        triangle.push(shifted);
    }
    square_triangle(&triangle)
}

/// Hoggatt polynomials `H_0^[m](q,t), ..., H_max_n^[m](q,t)`, exact integer q>=1.
///
/// Indexing: `rows[n]` = `H_n` for `0 <= n <= max_n` (`max_n + 1` rows), like
/// every generator in `sequences`.  The paper's normalization starts at n=1,
/// and `H_n` has degree at most `n - 1`, so `rows[0]` is the zero polynomial
/// `[0]` (an empty sum over `k`); use `rows[1..]` for the paper's rows.
///
/// **Indexing change in 0.6.0:** earlier versions returned `max_n` rows with
/// `rows[i]` = `H_{i+1}` (and no rows for `max_n = 0`); use
/// `hoggatt_polynomials_bigint(n, m, q)[1..]` to reproduce that.
///
/// The coefficient of t^k in row n (1<=n<=max_n, 0<=k<n) is
/// q^(m*k*(k+1)/2) times the volume generating function of plane partitions
/// in a k by (n-1-k) by m box, evaluated at q. This is the normalization of
/// Fang--Zhang--Zhao, arXiv:2505.05873v1, Section 1 and Conjecture 4.1.
///
/// MacMahon's product is evaluated as a quotient of Gaussian binomial
/// products, also valid at q=1. Every division is checked to be exact.
/// m=1 gives product_(i=1)^(n-1)(1+q^i*t); m=2,q=1 gives Narayana rows;
/// m=3,q=1 gives Baxter descent rows. No root property is assumed.
///
/// # Panics
///
/// Panics for m=0 or q=0, or index/exponent arithmetic exceeding usize.
pub fn hoggatt_polynomials_bigint(max_n: usize, m: usize, q: u32) -> Vec<Vec<BigInt>> {
    assert!(m >= 1, "Hoggatt box height m must be positive");
    assert!(q >= 1, "Hoggatt specialization q must be positive");
    if max_n == 0 {
        return vec![vec![BigInt::zero()]];
    }
    let max_top = max_n.checked_add(m).expect("Hoggatt index overflow") - 2;
    let q = BigInt::from(q);
    let mut powers = vec![BigInt::one()];
    for i in 1..=max_top {
        powers.push(&powers[i - 1] * &q);
    }
    let mut gaussian = vec![vec![BigInt::one()]];
    for n in 1..=max_top {
        let mut row = vec![BigInt::one(); n + 1];
        for k in 1..n {
            row[k] = &gaussian[n - 1][k - 1] + &powers[k] * &gaussian[n - 1][k];
        }
        gaussian.push(row);
    }
    std::iter::once(vec![BigInt::zero()])
        .chain((1..=max_n).map(|n| {
            let top = n + m - 2;
            let denominator: BigInt = (1..m).map(|j| &gaussian[top][j]).product();
            (0..n)
                .map(|k| {
                    let numerator: BigInt = (0..m).map(|i| &gaussian[top][k + i]).product();
                    let (count, remainder) = numerator.div_rem(&denominator);
                    assert!(remainder.is_zero(), "nonintegral MacMahon evaluation");
                    let twice_triangular = k.checked_mul(k + 1).expect("Hoggatt exponent overflow");
                    let exponent = (twice_triangular / 2)
                        .checked_mul(m)
                        .expect("Hoggatt exponent overflow");
                    count * pow_usize(&q, exponent)
                })
                .collect()
        }))
        .collect()
}

/// The `1/k`-Eulerian polynomials `A_0^(k), ..., A_max_n^(k)` for an exact
/// rational `k`; `rows[n]` = `A_n^(k)`, with `max_n + 1` rows.
///
/// For a positive integer `k`, `A_n^(k)(x)` is the ascent polynomial of the
/// `k`-inversion sequences of length `n` (Savage--Viswanathan, "The
/// 1/k-Eulerian polynomials", Electron. J. Combin. 19 (2012), #P9), with
/// exponential generating function `((1 - x) / (e^(k z (x - 1)) - x))^(1/k)`.
/// The rows are generated by the first-order differential recurrence
///
/// ```text
/// A_0 = 1,   A_(n+1) = (1 + k n x) A_n + k x (1 - x) D A_n,
/// ```
///
/// which is used for every exact rational `k` (so `A_1 = 1`,
/// `A_2 = 1 + k x`, `A_3 = 1 + (k^2 + 3k) x + k^2 x^2`, and
/// `A_n(1) = prod_(i<n) (1 + i k)`).  `k = 1` gives the Eulerian polynomials
/// [`super::eulerian_polynomials_bigint`].  For `n >= 1` the degree is at most
/// `n - 1`; no root property is asserted here.
pub fn one_over_k_eulerian_polynomials_rational(
    max_n: usize,
    k: &BigRational,
) -> Vec<Vec<BigRational>> {
    let mut rows = vec![vec![BigRational::one()]];
    for n in 0..max_n {
        let previous = &rows[n];
        let kn = k * BigRational::from_integer(BigInt::from(n));
        // (1 + k n x) A_n
        let mut next = vec_poly::mul(previous, &[BigRational::one(), kn]);
        // k x (1 - x) A_n'
        let derivative = vec_poly::derivative(previous);
        let correction = vec_poly::mul(&derivative, &[BigRational::zero(), k.clone(), -k.clone()]);
        next = vec_poly::add(&next, &correction);
        rows.push(vec_poly::trim(&next));
    }
    rows
}

/// The symmetric decompositions of the `1/k`-Eulerian polynomials:
/// `rows[n] = (a_n, b_n)` with `A_n^(k) = a_n + x b_n`, `a_n` symmetric with
/// respect to `d_n` (`x^(d_n) a_n(1/x) = a_n`) and `b_n` with respect to
/// `d_n - 1`, where `d_n = n - 1` for `n >= 1` and `d_0 = 0`; `max_n + 1`
/// rows.
///
/// The rows are the Stapledon (Brändén--Solus) decompositions of
/// [`one_over_k_eulerian_polynomials_rational`]; a zero component is `[0]`.
/// For `n >= 2`, `(a_n, b_n)` satisfy the coupled recurrence
/// `a_(n+1) = U_n a_n + x b_n`, `b_(n+1) = (k - 1) a_n + V_n b_n` with
/// `U_n = [1 + (1 + k (n - 1)) x] + k x (1 - x) D` and
/// `V_n = [k + (k + k (n - 2)) x] + k x (1 - x) D`, and the symmetric parts
/// `a_n` satisfy the eliminated second-order recurrence of Ma--Ma--Yeh--Yeh,
/// arXiv:2102.00899; for `k = 1`, `b_n = 0` and `a_n` is the Eulerian
/// polynomial.  So `a_2 = 1 + x`, `b_2 = k - 1`, `a_3 = 1 + (3k + 1) x + x^2`.
pub fn one_over_k_eulerian_symmetric_parts_rational(
    max_n: usize,
    k: &BigRational,
) -> Vec<(Vec<BigRational>, Vec<BigRational>)> {
    one_over_k_eulerian_polynomials_rational(max_n, k)
        .into_iter()
        .enumerate()
        .map(|(n, row)| {
            let degree_bound = n.saturating_sub(1);
            let (a, b) = Polynomial::new(row)
                .stapledon_decomposition(degree_bound)
                .expect("the 1/k-Eulerian row n has degree at most max(n - 1, 0)");
            (vec_poly::trim(a.coeffs()), vec_poly::trim(b.coeffs()))
        })
        .collect()
}

/// A recurrence with its identifying seed rows and the index of the first
/// seed, as returned by the `*_recurrence` functions of this module: the
/// seeds are the members `first_index, first_index + 1, ...` and every later
/// member follows from the rule ([`Recurrence::generate_rows_rational`]).
pub type RecurrenceWithSeeds = (Recurrence, Vec<Vec<BigRational>>, usize);

/// The peak polynomials `W_0, ..., W_max_n`: `rows[n]` is
/// `W_n(z) = sum z^(pk(sigma))` over the permutations `sigma` of `[n]`,
/// where `pk` counts the interior peaks `sigma(i - 1) < sigma(i) > sigma(i + 1)`,
/// `1 < i < n`; `max_n + 1` rows.
///
/// For `n >= 1` the coefficients are row `n` of OEIS A008303 (permutations of
/// `[n]` by number of peaks); `W_0 = W_1 = 1`.  Complementing the values
/// (`sigma(i) -> n + 1 - sigma(i)`) exchanges peaks and valleys, so these are
/// also the valley polynomials.
///
/// The rows are generated by the first-order recurrence of
/// [`peak_polynomials_recurrence`],
///
/// ```text
/// W_1 = 1,   W_n = (2 + (n - 2) z) W_(n-1) + 2 z (1 - z) W'_(n-1)   (n >= 2),
/// ```
///
/// which is the image of the Eulerian recurrence under the gamma expansion:
/// with `A_n(t) = sum_j gamma_(n,j) t^j (1 + t)^(n-1-2j)` (Foata--Schützenberger;
/// see Shareshian--Wachs, arXiv:1702.06666, equation (1.3)), the number of
/// permutations of `[n]` with `j` peaks is `2^(n-1-2j) gamma_(n,j)`
/// (the modified Foata--Strehl action; Stembridge, "Enriched P-partitions",
/// 1997; Petersen, *Eulerian Numbers*, 2015, Chapter 4), that is
/// `A_n(t) = ((1 + t)/2)^(n-1) W_n(4t/(1 + t)^2)` and
/// `W_n(z) = 2^(n-1) Gamma_n(z/4)` with `Gamma_n` the gamma polynomial of
/// `A_n`.  Substituting this into the gamma transport of the Eulerian rule
/// ([`crate::recurrence::transport_gamma`]) gives the recurrence above; the
/// `polytool` tests compare the rows with a brute-force count of peaks.
/// [`crate::recurrence::find_gamma_lifts`] certifies the identity with the
/// Eulerian polynomials for every `n` and proves every `W_n` real-rooted.
///
/// ```
/// use polytool::sequences::literature::peak_polynomials_bigint;
/// use num_bigint::BigInt;
///
/// let rows = peak_polynomials_bigint(5);
/// let row = |values: &[i64]| values.iter().map(|&v| BigInt::from(v)).collect::<Vec<_>>();
/// assert_eq!(rows[3], row(&[4, 2]));
/// assert_eq!(rows[5], row(&[16, 88, 16]));
/// ```
pub fn peak_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut rows = vec![vec![BigInt::one()]];
    if max_n >= 1 {
        rows.push(vec![BigInt::one()]);
    }
    for n in 2..=max_n {
        let previous = &rows[n - 1];
        let mut next = vec![BigInt::zero(); previous.len() + 1];
        let n_minus_two = BigInt::from(n - 2);
        for (j, value) in previous.iter().enumerate() {
            // (2 + (n - 2) z) W + 2 z W' - 2 z^2 W' on the monomial value z^j.
            let j_big = BigInt::from(j);
            next[j] += value * (BigInt::from(2) + BigInt::from(2) * &j_big);
            next[j + 1] += value * (&n_minus_two - BigInt::from(2) * &j_big);
        }
        rows.push(vec_poly::trim(&next));
    }
    rows
}

/// Lazily generate the peak polynomials `W_0, W_1, ...`; item `n` equals
/// `peak_polynomials_bigint(max_n)[n]`.  Uses the seeds `W_0 = W_1 = 1` and
/// the recurrence of [`peak_polynomials_recurrence`] from `n = 2` on; see
/// [`PolynomialFamilyRowIterator`].
pub fn peak_polynomials_iter() -> PolynomialFamilyRowIterator {
    PolynomialFamilyRowIterator::from_sparse_recurrence(
        &[
            (1, 0, &[(0, 0, 2), (1, 1, 1), (0, 1, -2)]),
            (1, 1, &[(0, 1, 2), (0, 2, -2)]),
        ],
        &[&[1], &[1]],
    )
}

/// The defining recurrence of [`peak_polynomials_bigint`] from `n = 1` on:
/// `W_n = (2 + (n - 2) z) W_(n-1) + 2 z (1 - z) W'_(n-1)` with the seed
/// `W_1 = 1` and first index 1 (the row `W_0 = 1` is not produced by the
/// rule).  The generator is defined by this recurrence, so it reproduces
/// every row `n >= 1`.
pub fn peak_polynomials_recurrence() -> RecurrenceWithSeeds {
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                rational_bivariate(&[
                    (0, 0, rational_integer(2)),
                    (1, 1, rational_integer(1)),
                    (0, 1, rational_integer(-2)),
                ]),
            ),
            RecurrenceTerm::new(
                1,
                1,
                rational_bivariate(&[(0, 1, rational_integer(2)), (0, 2, rational_integer(-2))]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    (rule, vec![vec![BigRational::one()]], 1)
}

/// The second-order recurrence of the symmetric parts `a_n` of
/// [`one_over_k_eulerian_symmetric_parts_rational`], from `n = 1` on.
///
/// The step `A_(n+1) = (1 + k n x) A_n + k x (1 - x) A'_n` of the
/// `1/k`-Eulerian polynomials splits into the coupled recurrence
///
/// ```text
/// a_(n+1) = U_n a_n + x b_n,   b_(n+1) = (k - 1) a_n + V_n b_n   (n >= 1),
/// U_n = 1 + (1 + k (n - 1)) x + k x (1 - x) D,
/// V_n = k + k (n - 1) x + k x (1 - x) D:
/// ```
///
/// the sum `a_(n+1) + x b_(n+1)` is `A_(n+1)` by expanding, and the parts
/// have the right symmetry because `1 + (1 + kappa d) x + kappa x (1 - x) D`
/// maps a polynomial symmetric at index `d` to one symmetric at `d + 1`
/// (`U_n` with `kappa = k`, `d = n - 1`; `V_n / k` with `kappa = 1`,
/// `d = n - 2`), so by uniqueness of the symmetric decomposition at
/// `d_(n+1) = n` they are the parts of `A_(n+1)`.  Eliminating `b` with
/// `x b_n = a_(n+1) - U_n a_n` and the operator identity
/// `x V_n b = W_n (x b)`, `W_n = k n x + k x (1 - x) D`, gives, for target
/// index `m = n + 2 >= 3`,
///
/// ```text
/// a_m = (U_(m-1) + W_(m-2)) a_(m-1) + ((k - 1) x - W_(m-2) U_(m-2)) a_(m-2),
/// ```
///
/// a rule with second derivatives whose coefficients are polynomials in
/// `(m, x)` and `k` (compare the eliminated recurrence of Ma--Ma--Yeh--Yeh,
/// arXiv:2102.00899; this one is derived and tested here).  Seeds `a_1 = 1`,
/// `a_2 = 1 + x`, first index 1; the leader is 1, so the rule and the seeds
/// define `a_n` for every `n >= 1`, and [`one_over_k_eulerian_symmetric_parts_rational`]
/// could be computed from it.
pub fn one_over_k_eulerian_symmetric_part_a_recurrence(k: &BigRational) -> RecurrenceWithSeeds {
    let one = BigRational::one();
    let k2 = k * k;
    let int = rational_integer;
    let first = [
        (0, 0, one.clone()),
        (0, 1, &one - int(4) * k),
        (1, 1, int(2) * k),
    ];
    let second = [
        (2, 2, -k2.clone()),
        (1, 2, int(6) * &k2 - k),
        (1, 1, -&k2 - k),
        (0, 2, int(-9) * &k2 + int(3) * k),
        (0, 1, int(3) * &k2 + int(2) * k - &one),
    ];
    let second_derivative_one = [
        (1, 3, int(2) * &k2),
        (1, 2, int(-2) * &k2),
        (0, 3, int(-7) * &k2 + k),
        (0, 2, int(8) * &k2),
        (0, 1, -&k2 - k),
    ];
    let rule = symmetric_parts_rule(k, &first, &second, &second_derivative_one);
    let seeds = vec![vec![one.clone()], vec![one.clone(), one]];
    (rule, seeds, 1)
}

/// The second-order recurrence of the parts `b_n` of
/// [`one_over_k_eulerian_symmetric_parts_rational`], from `n = 1` on.
///
/// From the coupled recurrence of
/// [`one_over_k_eulerian_symmetric_part_a_recurrence`],
/// `(k - 1) a_n = b_(n+1) - V_n b_n` and
/// `b_(n+2) = V_(n+1) b_(n+1) + U_n ((k - 1) a_n) + (k - 1) x b_n`, so for target
/// index `m = n + 2 >= 3`
///
/// ```text
/// b_m = (V_(m-1) + U_(m-2)) b_(m-1) + ((k - 1) x - U_(m-2) V_(m-2)) b_(m-2).
/// ```
///
/// Seeds `b_1 = 0`, `b_2 = k - 1`, first index 1; the leader is 1, so the
/// rule and the seeds define `b_n` for every `n >= 1` (all zero for
/// `k = 1`).
pub fn one_over_k_eulerian_symmetric_part_b_recurrence(k: &BigRational) -> RecurrenceWithSeeds {
    let one = BigRational::one();
    let k2 = k * k;
    let int = rational_integer;
    let first = [
        (0, 0, &one + k),
        (0, 1, &one - int(5) * k),
        (1, 1, int(2) * k),
    ];
    let second = [
        (2, 2, -k2.clone()),
        (1, 2, int(7) * &k2 - k),
        (1, 1, int(-2) * &k2 - k),
        (0, 2, int(-12) * &k2 + int(3) * k),
        (0, 1, int(6) * &k2 + int(3) * k - &one),
        (0, 0, -k.clone()),
    ];
    let second_derivative_one = [
        (1, 3, int(2) * &k2),
        (1, 2, int(-2) * &k2),
        (0, 3, int(-8) * &k2 + k),
        (0, 2, int(10) * &k2),
        (0, 1, int(-2) * &k2 - k),
    ];
    let rule = symmetric_parts_rule(k, &first, &second, &second_derivative_one);
    let seeds = vec![vec![BigRational::zero()], vec![k - &one]];
    (rule, seeds, 1)
}

/// `P_m = (c_0 + 2k x (1 - x) D) P_(m-1) + (e_0 + e_1 D - k^2 x^2 (1 - x)^2 D^2) P_(m-2)`
/// from the monomials of `c_0`, `e_0` and `e_1`.
fn symmetric_parts_rule(
    k: &BigRational,
    first: &[(usize, usize, BigRational)],
    second: &[(usize, usize, BigRational)],
    second_derivative_one: &[(usize, usize, BigRational)],
) -> Recurrence {
    let int = rational_integer;
    let k2 = k * k;
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, rational_bivariate(first)),
            RecurrenceTerm::new(
                1,
                1,
                rational_bivariate(&[(0, 1, int(2) * k), (0, 2, int(-2) * k)]),
            ),
            RecurrenceTerm::new(2, 0, rational_bivariate(second)),
            RecurrenceTerm::new(2, 1, rational_bivariate(second_derivative_one)),
            RecurrenceTerm::new(
                2,
                2,
                rational_bivariate(&[
                    (0, 2, -k2.clone()),
                    (0, 3, int(2) * &k2),
                    (0, 4, -k2.clone()),
                ]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

fn rational_integer(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

/// `sum c n^i x^j` from the monomials `(i, j, c)` (repeated monomials add).
fn rational_bivariate(monomials: &[(usize, usize, BigRational)]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let x_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![BigRational::zero(); x_size]; n_size];
    for (i, j, c) in monomials {
        coeffs[*i][*j] += c;
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

fn pow_usize(base: &BigInt, mut exponent: usize) -> BigInt {
    let mut result = BigInt::one();
    let mut factor = base.clone();
    while exponent != 0 {
        if exponent & 1 != 0 {
            result *= &factor;
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = &factor * &factor;
        }
    }
    result
}

fn add_linear_multiple(out: &mut [BigInt], row: &[BigInt], constant: i32, linear: i32) {
    for (k, value) in row.iter().enumerate() {
        out[k] += value * constant;
        out[k + 1] += value * linear;
    }
}

fn square_triangle(triangle: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    triangle
        .iter()
        .enumerate()
        .map(|(n, row)| {
            let mut out = vec![BigInt::zero(); n + 1];
            for (j, weight) in row.iter().enumerate() {
                for (k, value) in triangle[j].iter().enumerate() {
                    out[k] += weight * value;
                }
            }
            out
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn big(row: &[i64]) -> Vec<BigInt> {
        row.iter().map(|&x| BigInt::from(x)).collect()
    }

    fn rational(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(BigInt::from(numerator), BigInt::from(denominator))
    }

    #[test]
    fn one_over_k_eulerian_rows_match_eulerian_closed_values_and_parts() {
        let one = rational(1, 1);
        let eulerian = super::super::eulerian_polynomials_bigint(9);
        let rows = one_over_k_eulerian_polynomials_rational(9, &one);
        assert_eq!(rows.len(), 10);
        for (row, expected) in rows.iter().zip(&eulerian) {
            let expected: Vec<BigRational> = expected
                .iter()
                .map(|c| BigRational::from_integer(c.clone()))
                .collect();
            assert_eq!(*row, vec_poly::trim(&expected));
        }
        for k in [rational(2, 1), rational(3, 2), rational(1, 3)] {
            let rows = one_over_k_eulerian_polynomials_rational(8, &k);
            assert_eq!(rows[2], vec![rational(1, 1), k.clone()]);
            let k_squared = &k * &k;
            assert_eq!(
                rows[3],
                vec![rational(1, 1), &k_squared + &k * rational(3, 1), k_squared]
            );
            let mut value_at_one = rational(1, 1);
            for (n, row) in rows.iter().enumerate() {
                let sum: BigRational = row.iter().cloned().sum();
                assert_eq!(sum, value_at_one, "A_{n}(1) for k = {k}");
                value_at_one *= rational(1, 1) + &k * BigRational::from_integer(BigInt::from(n));
            }
            let parts = one_over_k_eulerian_symmetric_parts_rational(8, &k);
            assert_eq!(parts[2].0, vec![rational(1, 1), rational(1, 1)]);
            assert_eq!(parts[2].1, vec![&k - rational(1, 1)]);
            assert_eq!(
                parts[3].0,
                vec![
                    rational(1, 1),
                    &k * rational(3, 1) + rational(1, 1),
                    rational(1, 1)
                ]
            );
            for (n, (a, b)) in parts.iter().enumerate() {
                let collapsed = vec_poly::add(a, &vec_poly::shift(b));
                assert_eq!(vec_poly::trim(&collapsed), rows[n]);
                let d = n.saturating_sub(1);
                assert_eq!(
                    Some(a.clone()),
                    vec_poly::reverse(a, d).map(|r| vec_poly::trim(&r))
                );
                if d >= 1 {
                    assert_eq!(
                        Some(b.clone()),
                        vec_poly::reverse(b, d - 1).map(|r| vec_poly::trim(&r))
                    );
                }
            }
        }
        let parts = one_over_k_eulerian_symmetric_parts_rational(6, &one);
        assert!(parts.iter().all(|(_, b)| vec_poly::is_zero(b)));
    }

    #[test]
    fn delannoy_square_recurrence_matches_matrix_product() {
        let d = delannoy_polynomials_bigint(25);
        assert_eq!(d[3], big(&[1, 5, 5, 1]));
        assert_eq!(delannoy_square_polynomials_bigint(25), square_triangle(&d));
    }

    #[test]
    fn delannoy_square_three_lag_differential_identity() {
        let rows = delannoy_square_polynomials_bigint(120);
        for n in 1..=120 {
            let mut right = vec![BigInt::zero(); n + 1];
            let n_i = n as i32;
            for (j, constant, linear, derivative) in [
                (1, 3 - n_i, 1, [4_i64, 11, 1]),
                (2, 3 * n_i - 2, -2, [-4, -6, -2]),
                (3, 5 - n_i, -1, [0, 5, -1]),
            ] {
                if n < j {
                    continue;
                }
                let previous = &rows[n - j];
                add_linear_multiple(&mut right, previous, constant, linear);
                for k in 1..previous.len() {
                    for (ell, factor) in derivative.iter().enumerate() {
                        if *factor != 0 {
                            right[k - 1 + ell] += &previous[k] * BigInt::from(k) * factor;
                        }
                    }
                }
            }
            let left: Vec<BigInt> = rows[n].iter().map(|c| c * n_i).collect();
            assert_eq!(left, right, "n={n}");
        }
    }

    #[test]
    fn eulerian_square_retains_paper_indexing() {
        let a = eulerian_square_polynomials_bigint(4);
        assert_eq!(a[0], big(&[1]));
        assert_eq!(a[1], big(&[0, 1]));
        assert_eq!(a[2], big(&[0, 2, 1]));
        assert_eq!(a[3], big(&[0, 6, 8, 1]));
        assert_eq!(a[4], big(&[0, 24, 66, 22, 1]));
    }

    #[test]
    fn hoggatt_at_one_matches_narayana_and_baxter() {
        assert_eq!(
            hoggatt_polynomials_bigint(8, 2, 1)[1..],
            super::super::narayana_polynomials_bigint(8)[1..]
        );
        let h = hoggatt_polynomials_bigint(5, 3, 1);
        assert_eq!(h.len(), 6);
        assert_eq!(h[0], big(&[0]));
        assert_eq!(h[1], big(&[1]));
        assert_eq!(h[3], big(&[1, 4, 1]));
        assert_eq!(h[4], big(&[1, 10, 10, 1]));
        assert_eq!(h[5], big(&[1, 20, 50, 20, 1]));
    }

    #[test]
    fn hoggatt_height_one_matches_q_binomial_product() {
        for q in [1, 2, 3, 10] {
            let h = hoggatt_polynomials_bigint(9, 1, q);
            let mut product = vec![BigInt::one()];
            assert_eq!(h.len(), 10);
            assert_eq!(h[0], big(&[0]));
            for (n, row) in h.iter().enumerate().skip(1) {
                assert_eq!(*row, product);
                let factor = pow_usize(&BigInt::from(q), n);
                let mut next = vec![BigInt::zero(); product.len() + 1];
                for (i, c) in product.iter().enumerate() {
                    next[i] += c;
                    next[i + 1] += c * &factor;
                }
                product = next;
            }
        }
    }

    #[test]
    fn hoggatt_matches_small_plane_partition_enumeration() {
        // Exhaust all 2x2 arrays with weakly decreasing rows/columns, entries 0..m.
        for m in 1..=3 {
            for q in [1, 2, 3] {
                let mut weight = BigInt::zero();
                for a in 0..=m {
                    for b in 0..=a {
                        for c in 0..=a {
                            for d in 0..=b.min(c) {
                                weight += pow_usize(&BigInt::from(q), a + b + c + d);
                            }
                        }
                    }
                }
                assert_eq!(
                    hoggatt_polynomials_bigint(5, m, q)[5][2],
                    weight * pow_usize(&BigInt::from(q), 3 * m)
                );
            }
        }
    }

    #[test]
    fn empty_family_and_zero_degree_boundaries() {
        assert_eq!(hoggatt_polynomials_bigint(0, 2, 2), vec![big(&[0])]);
        assert_eq!(hoggatt_polynomials_bigint(0, 1, 1), vec![big(&[0])]);
        assert_eq!(
            hoggatt_polynomials_bigint(1, 3, 2),
            vec![big(&[0]), big(&[1])]
        );
        assert_eq!(delannoy_square_polynomials_bigint(0), vec![big(&[1])]);
        assert_eq!(eulerian_square_polynomials_bigint(0), vec![big(&[1])]);
    }

    #[test]
    fn hoggatt_rank_shift_identity_including_vanishing_top_coefficient() {
        // Product_j [n-k+j]_q * c_(n+1,k) = Product_j [n+j]_q * c_(n,k).
        // At k=n, the first product vanishes; the new leading term is separate.
        for m in 1..=4 {
            for q in [1, 2, 3, 11] {
                let rows = hoggatt_polynomials_bigint(25, m, q);
                let q_integer = |r: usize| -> BigInt {
                    if q == 1 {
                        BigInt::from(r)
                    } else {
                        (pow_usize(&BigInt::from(q), r) - 1) / (q - 1)
                    }
                };
                for n in 1..25 {
                    let right_factor: BigInt = (0..m).map(|j| q_integer(n + j)).product();
                    for k in 0..=n {
                        let left_factor: BigInt = (0..m).map(|j| q_integer(n - k + j)).product();
                        let previous = rows[n].get(k).cloned().unwrap_or_else(BigInt::zero);
                        assert_eq!(left_factor * &rows[n + 1][k], &right_factor * previous);
                    }
                    assert_eq!(
                        rows[n + 1][n],
                        pow_usize(&BigInt::from(q), m * n * (n + 1) / 2)
                    );
                }
            }
        }
    }

    #[test]
    #[should_panic(expected = "must be positive")]
    fn zero_q_is_rejected() {
        hoggatt_polynomials_bigint(2, 2, 0);
    }

    /// Every permutation of `0..n` (Heap's algorithm).
    fn permutations(n: usize) -> Vec<Vec<usize>> {
        let mut current: Vec<usize> = (0..n).collect();
        let mut out = vec![current.clone()];
        let mut counters = vec![0; n];
        let mut i = 1;
        while i < n {
            if counters[i] < i {
                let j = if i % 2 == 0 { 0 } else { counters[i] };
                current.swap(j, i);
                out.push(current.clone());
                counters[i] += 1;
                i = 1;
            } else {
                counters[i] = 0;
                i += 1;
            }
        }
        out
    }

    #[test]
    fn peak_and_valley_polynomials_match_brute_force_counts() {
        let rows = peak_polynomials_bigint(8);
        for (n, row) in rows.iter().enumerate() {
            let mut peaks = vec![BigInt::zero(); n / 2 + 1];
            let mut valleys = vec![BigInt::zero(); n / 2 + 1];
            for sigma in permutations(n) {
                let interior = |test: fn(usize, usize, usize) -> bool| {
                    (1..n.saturating_sub(1))
                        .filter(|&i| test(sigma[i - 1], sigma[i], sigma[i + 1]))
                        .count()
                };
                peaks[interior(|a, b, c| a < b && b > c)] += 1;
                valleys[interior(|a, b, c| a > b && b < c)] += 1;
            }
            assert_eq!(*row, vec_poly::trim(&peaks), "peaks at n = {n}");
            assert_eq!(*row, vec_poly::trim(&valleys), "valleys at n = {n}");
        }
        // A008303 rows 6 and 7.
        assert_eq!(rows[6], big(&[32, 416, 272]));
        assert_eq!(rows[7], big(&[64, 1824, 2880, 272]));
        let lazy: Vec<_> = peak_polynomials_iter().take(9).collect();
        assert_eq!(lazy, rows);
    }

    fn rational_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<BigRational>> {
        rows.iter()
            .map(|row| row.iter().cloned().map(BigRational::from_integer).collect())
            .collect()
    }

    #[test]
    fn peak_recurrence_reproduces_the_generator_and_the_eulerian_gamma_identity() {
        let (rule, seeds, first_index) = peak_polynomials_recurrence();
        let rows = rational_rows(&peak_polynomials_bigint(20));
        let generated = rule
            .generate_rows_rational(&seeds, first_index, 20)
            .unwrap();
        for (left, right) in generated.iter().zip(&rows[1..]) {
            assert_eq!(vec_poly::trim(left), vec_poly::trim(right));
        }
        // A_n(t) = ((1 + t)/2)^(n-1) W_n(4t/(1 + t)^2): with W_n = sum w_j z^j,
        // 2^(n-1) A_n = sum_j w_j 4^j t^j (1 + t)^(n-1-2j).
        let eulerian = super::super::eulerian_polynomials_bigint(20);
        for n in 1..=20 {
            let mut sum = vec![BigInt::zero()];
            for (j, w) in rows[n].iter().enumerate() {
                let w = w.numer() * pow_usize(&BigInt::from(4), j);
                let mut term = vec![w];
                for _ in 0..n - 1 - 2 * j {
                    term = vec_poly::mul(&term, &[BigInt::one(), BigInt::one()]);
                }
                let mut shifted = vec![BigInt::zero(); j];
                shifted.extend(term);
                sum = vec_poly::add(&sum, &shifted);
            }
            let scaled: Vec<BigInt> = eulerian[n]
                .iter()
                .map(|c| c * pow_usize(&BigInt::from(2), n - 1))
                .collect();
            assert_eq!(vec_poly::trim(&sum), vec_poly::trim(&scaled), "n = {n}");
        }
    }

    #[test]
    fn symmetric_part_recurrences_reproduce_the_decompositions() {
        for k in [
            rational(1, 1),
            rational(2, 1),
            rational(3, 1),
            rational(1, 2),
            rational(-5, 3),
        ] {
            let parts = one_over_k_eulerian_symmetric_parts_rational(16, &k);
            for (component, (rule, seeds, first_index)) in [
                one_over_k_eulerian_symmetric_part_a_recurrence(&k),
                one_over_k_eulerian_symmetric_part_b_recurrence(&k),
            ]
            .into_iter()
            .enumerate()
            {
                let generated = rule
                    .generate_rows_rational(&seeds, first_index, 16)
                    .unwrap();
                for (offset, row) in generated.iter().enumerate() {
                    let (a, b) = &parts[first_index + offset];
                    let expected = if component == 0 { a } else { b };
                    assert_eq!(
                        vec_poly::trim(row),
                        vec_poly::trim(expected),
                        "k = {k}, component {component}, n = {}",
                        first_index + offset
                    );
                }
            }
        }
    }
}
