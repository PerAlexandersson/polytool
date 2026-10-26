//! Exact integer combinatorial helpers: binomials, factorials, falling
//! factorials, multinomials, and their q-analogues ([`q_integer`],
//! [`q_factorial`], [`q_binomial`], [`q_catalan`], [`q_narayana`]).
//!
//! The arbitrary-precision functions return [`BigInt`] and never overflow.
//! [`binomial_u128`] is the only machine-integer function and reports overflow
//! as `None` instead of wrapping or panicking.
//!
//! The same names and semantics are used by the `combinatoric-core` crate's
//! `numeric` module, so code can move between the two crates unchanged.
//!
//! # Conventions
//!
//! * `binomial_bigint(n, k)` is `0` when `k > n` (no panic).
//! * `factorial_bigint(0) = 1` and `falling_factorial_bigint(x, 0) = 1`.
//! * The falling factorial is `x (x-1) ... (x-k+1)` for any integer `x`,
//!   including negative `x`; `x` is passed as `&BigInt`.
//! * `multinomial_bigint(&[])` is `1`, and zero parts are allowed.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

/// The binomial coefficient `C(n, k)` as an exact [`BigInt`].
///
/// Returns `0` when `k > n`. Uses the multiplicative formula with the smaller
/// of `k` and `n - k` factors, so it takes `O(min(k, n-k))` big-integer
/// operations.
///
/// ```
/// use polytool::numeric::binomial_bigint;
/// assert_eq!(binomial_bigint(5, 2).to_string(), "10");
/// assert_eq!(binomial_bigint(3, 5).to_string(), "0");
/// ```
pub fn binomial_bigint(n: usize, k: usize) -> BigInt {
    if k > n {
        return BigInt::from(0);
    }
    let k = k.min(n - k);
    let mut result = BigInt::one();
    for i in 1..=k {
        // After step i, `result == C(n - k + i, i)`, so the division is exact.
        result = result * BigInt::from(n - k + i) / BigInt::from(i);
    }
    result
}

/// The checked binomial coefficient `C(n, k)` as a `u128`.
///
/// Returns `Some(0)` when `k > n` and `None` if the value does not fit in a
/// `u128`. Intermediate values never exceed the result, and every product is
/// checked, so the function never panics.
///
/// ```
/// use polytool::numeric::binomial_u128;
/// assert_eq!(binomial_u128(10, 3), Some(120));
/// assert_eq!(binomial_u128(200, 100), None);
/// ```
pub fn binomial_u128(n: usize, k: usize) -> Option<u128> {
    if k > n {
        return Some(0);
    }
    let k = k.min(n - k);
    let mut result = 1u128;
    for i in 1..=k {
        // C(n-k+i, i) = C(n-k+i-1, i-1) * factor / i.  After cancelling
        // g = gcd(factor, i), the reduced divisor is coprime to the reduced
        // factor and therefore divides the previous value exactly.
        let mut factor = (n - k + i) as u128;
        let mut divisor = i as u128;
        let g = Integer::gcd(&factor, &divisor);
        factor /= g;
        divisor /= g;
        result = (result / divisor).checked_mul(factor)?;
    }
    Some(result)
}

/// The factorial `n!` as an exact [`BigInt`] (`0! = 1`).
///
/// ```
/// use polytool::numeric::factorial_bigint;
/// assert_eq!(factorial_bigint(20).to_string(), "2432902008176640000");
/// ```
pub fn factorial_bigint(n: usize) -> BigInt {
    let mut result = BigInt::one();
    for i in 2..=n {
        result *= BigInt::from(i);
    }
    result
}

/// The falling factorial `x (x-1) ... (x-k+1)` as an exact [`BigInt`].
///
/// `k = 0` gives `1`. The argument `x` may be negative or smaller than `k`
/// (then the product contains a zero factor exactly when `0 <= x < k`).
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::numeric::falling_factorial_bigint;
/// assert_eq!(falling_factorial_bigint(&BigInt::from(-2), 3), BigInt::from(-24));
/// assert_eq!(falling_factorial_bigint(&BigInt::from(5), 2), BigInt::from(20));
/// ```
pub fn falling_factorial_bigint(x: &BigInt, k: usize) -> BigInt {
    let mut result = BigInt::one();
    for offset in 0..k {
        result *= x - BigInt::from(offset);
    }
    result
}

/// The multinomial coefficient `(p_1 + ... + p_m)! / (p_1! ... p_m!)`.
///
/// The empty slice and all-zero parts give `1`. Computed as a product of
/// binomials `C(p_1 + ... + p_i, p_i)`, so no huge factorial is formed.
///
/// ```
/// use polytool::numeric::multinomial_bigint;
/// assert_eq!(multinomial_bigint(&[2, 1, 1]).to_string(), "12");
/// ```
pub fn multinomial_bigint(parts: &[usize]) -> BigInt {
    let mut total = 0usize;
    let mut result = BigInt::one();
    for &part in parts {
        total += part;
        result *= binomial_bigint(total, part);
    }
    result
}

// ---------------------------------------------------------------------------
// q-analogues
// ---------------------------------------------------------------------------

/// The zero polynomial in canonical form.
fn zero_poly() -> Vec<BigInt> {
    vec![BigInt::from(0)]
}

/// Exact division of `numerator` by a *monic* `denominator` over `Z`.
///
/// # Panics
/// Panics if the denominator is not monic or does not divide the numerator.
fn exact_div_monic(numerator: &[BigInt], denominator: &[BigInt]) -> Vec<BigInt> {
    let den = crate::vec_poly::trim(denominator);
    let dd = den.len() - 1;
    assert!(
        den[dd].is_one(),
        "q-analogue division needs a monic divisor"
    );
    let mut rem = crate::vec_poly::trim(numerator);
    if rem.len() < den.len() {
        assert!(rem.iter().all(Zero::is_zero), "inexact q-analogue division");
        return zero_poly();
    }
    let mut quotient = vec![BigInt::from(0); rem.len() - dd];
    for i in (0..quotient.len()).rev() {
        let c = rem[i + dd].clone();
        for (j, d) in den.iter().enumerate() {
            rem[i + j] -= &c * d;
        }
        quotient[i] = c;
    }
    assert!(rem.iter().all(Zero::is_zero), "inexact q-analogue division");
    crate::vec_poly::trim(&quotient)
}

/// The q-integer `[n]_q = 1 + q + ... + q^(n-1)`, as ascending coefficients.
///
/// `[0]_q = 0` is returned as the canonical zero polynomial `[0]`.
///
/// ```
/// use polytool::numeric::q_integer;
/// assert_eq!(q_integer(3).iter().map(|c| c.to_string()).collect::<Vec<_>>(), ["1", "1", "1"]);
/// ```
pub fn q_integer(n: usize) -> Vec<BigInt> {
    if n == 0 {
        zero_poly()
    } else {
        vec![BigInt::one(); n]
    }
}

/// The q-factorial `[n]_q! = [1]_q [2]_q ... [n]_q`, with `[0]_q! = 1`.
///
/// Ascending coefficients; the degree is `n (n - 1) / 2`.
///
/// ```
/// use polytool::numeric::q_factorial;
/// // [3]_q! = (1 + q)(1 + q + q^2) = 1 + 2q + 2q^2 + q^3.
/// let coefficients: Vec<String> = q_factorial(3).iter().map(|c| c.to_string()).collect();
/// assert_eq!(coefficients, ["1", "2", "2", "1"]);
/// ```
pub fn q_factorial(n: usize) -> Vec<BigInt> {
    let mut result = vec![BigInt::one()];
    for i in 2..=n {
        result = crate::vec_poly::mul(&result, &q_integer(i));
    }
    result
}

/// The Gaussian binomial coefficient `[n choose k]_q`, as ascending
/// coefficients of a polynomial of degree `k (n - k)`.
///
/// Returns the zero polynomial `[0]` when `k > n`.  It generates partitions
/// in a `k x (n - k)` box by size (equivalently, `k`-subsets by inversions or
/// words with `k` ones by inversions).  At `q = 1` it is the ordinary binomial
/// coefficient.
///
/// ```
/// use polytool::numeric::q_binomial;
/// // [4 choose 2]_q = 1 + q + 2q^2 + q^3 + q^4.
/// let coefficients: Vec<String> = q_binomial(4, 2).iter().map(|c| c.to_string()).collect();
/// assert_eq!(coefficients, ["1", "1", "2", "1", "1"]);
/// ```
pub fn q_binomial(n: usize, k: usize) -> Vec<BigInt> {
    if k > n {
        return zero_poly();
    }
    let k = k.min(n - k);
    // After step i the value is [n - k + i choose i]_q, so each division is exact.
    let mut result = vec![BigInt::one()];
    for i in 1..=k {
        result = exact_div_monic(
            &crate::vec_poly::mul(&result, &q_integer(n - k + i)),
            &q_integer(i),
        );
    }
    result
}

/// MacMahon's q-Catalan polynomial `C_n(q) = [2n choose n]_q / [n + 1]_q`.
///
/// This is the *maj/inv* q-Catalan (the hook-formula q-analogue for the
/// `2 x n` rectangle: up to a power-of-`q` normalization it is the
/// major-index generating function of Dyck paths, equivalently of standard
/// Young tableaux of shape `(n, n)`), and is **not** the Carlitz--Riordan
/// area q-Catalan.  Ascending coefficients, constant term `1`, degree
/// `n (n - 1)`, palindromic, and `C_n(1)` is the Catalan number.  Also
/// `C_n(q) = [2n choose n]_q - q [2n choose n + 1]_q`.
///
/// ```
/// use polytool::numeric::q_catalan;
/// // C_2(q) = 1 + q^2.
/// let coefficients: Vec<String> = q_catalan(2).iter().map(|c| c.to_string()).collect();
/// assert_eq!(coefficients, ["1", "0", "1"]);
/// ```
pub fn q_catalan(n: usize) -> Vec<BigInt> {
    exact_div_monic(&q_binomial(2 * n, n), &q_integer(n + 1))
}

/// The q-Narayana polynomial `q^(k(k-1)) [n choose k]_q [n choose k-1]_q / [n]_q`.
///
/// This is exactly the normalization `q^(k(k-1))` of the q-Narayana
/// polynomial in `qNarayana`; at `q = 1` it is the Narayana number
/// `N(n,k) = C(n,k) C(n,k-1) / n`.  Ascending coefficients.  It is the zero
/// polynomial `[0]` unless `1 <= k <= n` (for example when `k = 0`).  For
/// `1 <= k <= n` it is a genuine polynomial with lowest term `q^(k(k-1))`.
///
/// # Panics
/// Panics if `n == 0`, where `1 / [0]_q` is undefined.
///
/// ```
/// use polytool::numeric::q_narayana;
/// // q^2 [3 choose 2]_q [3 choose 1]_q / [3]_q = q^2 + q^3 + q^4.
/// let coefficients: Vec<String> = q_narayana(3, 2).iter().map(|c| c.to_string()).collect();
/// assert_eq!(coefficients, ["0", "0", "1", "1", "1"]);
/// ```
pub fn q_narayana(n: usize, k: usize) -> Vec<BigInt> {
    assert!(n >= 1, "q_narayana requires n >= 1");
    if k == 0 || k > n {
        return zero_poly();
    }
    let product = crate::vec_poly::mul(&q_binomial(n, k), &q_binomial(n, k - 1));
    let quotient = exact_div_monic(&product, &q_integer(n));
    crate::vec_poly::shift_by(&quotient, k * (k - 1))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pascal's triangle as an independent oracle.
    fn pascal(rows: usize) -> Vec<Vec<BigInt>> {
        let mut table = vec![vec![BigInt::one()]];
        for n in 1..rows {
            let prev = &table[n - 1];
            let mut row = vec![BigInt::one(); n + 1];
            for k in 1..n {
                row[k] = &prev[k - 1] + &prev[k];
            }
            table.push(row);
        }
        table
    }

    #[test]
    fn binomial_matches_pascal_triangle() {
        let table = pascal(70);
        for (n, row) in table.iter().enumerate() {
            for (k, value) in row.iter().enumerate() {
                assert_eq!(&binomial_bigint(n, k), value, "C({n}, {k})");
                assert_eq!(
                    binomial_u128(n, k).map(BigInt::from).as_ref(),
                    Some(value),
                    "u128 C({n}, {k})"
                );
            }
            assert_eq!(binomial_bigint(n, n + 1), BigInt::from(0));
            assert_eq!(binomial_u128(n, n + 3), Some(0));
        }
    }

    #[test]
    fn binomial_u128_overflow_boundary() {
        // C(130, 65) < 2^128 <= C(132, 66).
        let fits = binomial_bigint(130, 65);
        assert_eq!(binomial_u128(130, 65).map(BigInt::from), Some(fits));
        assert!(binomial_bigint(132, 66) > BigInt::from(u128::MAX));
        assert_eq!(binomial_u128(132, 66), None);
        // Values near the limit at the edges of the row still fit.
        assert_eq!(binomial_u128(1000, 1), Some(1000));
        assert_eq!(binomial_u128(1000, 999), Some(1000));
        // Every u128 result agrees with the BigInt value whenever it fits.
        for n in 120..=140 {
            for k in 0..=n {
                let exact = binomial_bigint(n, k);
                match binomial_u128(n, k) {
                    Some(v) => assert_eq!(BigInt::from(v), exact),
                    None => assert!(exact > BigInt::from(u128::MAX)),
                }
            }
        }
    }

    #[test]
    fn factorial_and_falling_factorial() {
        assert_eq!(factorial_bigint(0), BigInt::one());
        assert_eq!(factorial_bigint(1), BigInt::one());
        assert_eq!(factorial_bigint(5), BigInt::from(120));
        assert_eq!(
            factorial_bigint(25).to_string(),
            "15511210043330985984000000"
        );
        for n in 0..20usize {
            let x = BigInt::from(n);
            assert_eq!(falling_factorial_bigint(&x, n), factorial_bigint(n));
            for k in 0..=n {
                // (n)_k = k! C(n, k).
                assert_eq!(
                    falling_factorial_bigint(&x, k),
                    factorial_bigint(k) * binomial_bigint(n, k)
                );
            }
            assert_eq!(falling_factorial_bigint(&x, n + 1), BigInt::from(0));
        }
        // (-x)_k = (-1)^k x (x+1) ... (x+k-1).
        assert_eq!(
            falling_factorial_bigint(&BigInt::from(-1), 4),
            BigInt::from(24)
        );
        assert_eq!(
            falling_factorial_bigint(&BigInt::from(-3), 0),
            BigInt::one()
        );
    }

    #[test]
    fn multinomial_matches_factorial_formula() {
        assert_eq!(multinomial_bigint(&[]), BigInt::one());
        assert_eq!(multinomial_bigint(&[0, 0]), BigInt::one());
        assert_eq!(multinomial_bigint(&[3]), BigInt::one());
        for parts in [vec![2, 1, 1], vec![3, 0, 4], vec![5, 5, 5], vec![1; 8]] {
            let total: usize = parts.iter().sum();
            let denominator = parts
                .iter()
                .fold(BigInt::one(), |acc, &p| acc * factorial_bigint(p));
            assert_eq!(
                multinomial_bigint(&parts),
                factorial_bigint(total) / denominator,
                "{parts:?}"
            );
        }
    }

    // -- q-analogues -------------------------------------------------------

    fn qp(values: &[i64]) -> Vec<BigInt> {
        values.iter().map(|&v| BigInt::from(v)).collect()
    }

    fn eval_at_one(p: &[BigInt]) -> BigInt {
        p.iter().sum()
    }

    #[test]
    fn q_integers_and_factorials() {
        assert_eq!(q_integer(0), qp(&[0]));
        assert_eq!(q_integer(1), qp(&[1]));
        assert_eq!(q_integer(4), qp(&[1, 1, 1, 1]));
        assert_eq!(q_factorial(0), qp(&[1]));
        assert_eq!(q_factorial(1), qp(&[1]));
        assert_eq!(q_factorial(3), qp(&[1, 2, 2, 1]));
        assert_eq!(q_factorial(4), qp(&[1, 3, 5, 6, 5, 3, 1]));
        for n in 0..9 {
            assert_eq!(eval_at_one(&q_factorial(n)), factorial_bigint(n));
        }
    }

    #[test]
    fn q_binomials_match_pascal_recurrence_and_binomials() {
        // [n,k] = [n-1,k-1] + q^k [n-1,k].
        for n in 1..=10usize {
            for k in 1..n {
                let left = q_binomial(n - 1, k - 1);
                let right = crate::vec_poly::shift_by(&q_binomial(n - 1, k), k);
                assert_eq!(
                    q_binomial(n, k),
                    crate::vec_poly::add(&left, &right),
                    "n={n} k={k}"
                );
            }
        }
        for n in 0..=10usize {
            assert_eq!(q_binomial(n, 0), qp(&[1]));
            assert_eq!(q_binomial(n, n), qp(&[1]));
            for k in 0..=n {
                assert_eq!(eval_at_one(&q_binomial(n, k)), binomial_bigint(n, k));
                assert_eq!(q_binomial(n, k), q_binomial(n, n - k));
            }
            assert_eq!(q_binomial(n, n + 1), qp(&[0]));
        }
        assert_eq!(q_binomial(4, 2), qp(&[1, 1, 2, 1, 1]));
        // [n,k] [k]! [n-k]! = [n]!.
        let lhs = crate::vec_poly::mul(
            &q_binomial(7, 3),
            &crate::vec_poly::mul(&q_factorial(3), &q_factorial(4)),
        );
        assert_eq!(lhs, q_factorial(7));
    }

    #[test]
    fn q_catalan_values_and_identities() {
        assert_eq!(q_catalan(0), qp(&[1]));
        assert_eq!(q_catalan(1), qp(&[1]));
        assert_eq!(q_catalan(2), qp(&[1, 0, 1]));
        assert_eq!(q_catalan(3), qp(&[1, 0, 1, 1, 1, 0, 1]));
        let catalan = [1i64, 1, 2, 5, 14, 42, 132, 429, 1430];
        for (n, &c) in catalan.iter().enumerate() {
            let p = q_catalan(n);
            assert_eq!(eval_at_one(&p), BigInt::from(c));
            assert_eq!(p.len(), n * n.saturating_sub(1) + 1);
            let reversed: Vec<BigInt> = p.iter().rev().cloned().collect();
            assert_eq!(p, reversed, "palindromic n={n}");
            // C_n(q) = [2n,n] - q [2n,n+1].
            let difference = crate::vec_poly::sub(
                &q_binomial(2 * n, n),
                &crate::vec_poly::shift(&q_binomial(2 * n, n + 1)),
            );
            assert_eq!(p, difference, "difference identity n={n}");
        }
    }

    #[test]
    fn q_narayana_values_and_sums() {
        assert_eq!(q_narayana(1, 1), qp(&[1]));
        assert_eq!(q_narayana(3, 1), qp(&[1]));
        assert_eq!(q_narayana(3, 2), qp(&[0, 0, 1, 1, 1]));
        assert_eq!(q_narayana(3, 3), qp(&[0, 0, 0, 0, 0, 0, 1]));
        assert_eq!(q_narayana(4, 0), qp(&[0]));
        assert_eq!(q_narayana(4, 5), qp(&[0]));
        // N(n,k) at q = 1, and row sums are Catalan numbers.
        for n in 1..=8usize {
            let mut total = BigInt::from(0);
            for k in 1..=n {
                let p = q_narayana(n, k);
                let expected = binomial_bigint(n, k) * binomial_bigint(n, k - 1) / BigInt::from(n);
                assert_eq!(eval_at_one(&p), expected, "N({n},{k})");
                assert!(p[..k * (k - 1)].iter().all(|c| c.is_zero()));
                total += eval_at_one(&p);
            }
            assert_eq!(total, binomial_bigint(2 * n, n) / BigInt::from(n + 1));
        }
    }

    #[test]
    #[should_panic(expected = "n >= 1")]
    fn q_narayana_rejects_zero_n() {
        let _ = q_narayana(0, 0);
    }
}
