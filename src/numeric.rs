//! Exact integer combinatorial helpers: binomials, factorials, falling
//! factorials, and multinomials.
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
use num_traits::One;

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
}
