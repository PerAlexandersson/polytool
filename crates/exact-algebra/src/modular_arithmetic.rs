//! Overflow-safe word-size residue arithmetic.
//!
//! Arithmetic arguments are canonical residues in `0..modulus`, with nonzero
//! modulus. Prime-specific inversion uses a validated prime at public boundaries.

use crate::primes::ValidatedPrime;
use num_bigint::BigInt;
use num_traits::{Signed, ToPrimitive};

/// Multiply canonical residues modulo a nonzero modulus, widening when needed.
///
/// # Panics
///
/// May panic if `modulus` is zero or an argument is not a canonical residue.
#[inline]
pub fn multiply_modulo_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if modulus <= u32::MAX as u64 {
        return (a * b) % modulus;
    }
    ((a as u128 * b as u128) % modulus as u128) as u64
}

/// Subtract canonical residues modulo a nonzero modulus.
///
/// # Panics
///
/// May panic if `modulus` is zero or an argument is not a canonical residue.
#[inline]
pub fn subtract_modulo_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

/// Add canonical residues modulo a nonzero modulus without overflow.
///
/// # Panics
///
/// May panic if `modulus` is zero or an argument is not a canonical residue.
#[inline]
pub fn add_modulo_u64(a: u64, b: u64, modulus: u64) -> u64 {
    let (sum, overflow) = a.overflowing_add(b);
    if overflow {
        ((a as u128 + b as u128) % modulus as u128) as u64
    } else if sum >= modulus {
        sum - modulus
    } else {
        sum
    }
}

/// Raise a canonical residue to a nonnegative power modulo a nonzero modulus.
///
/// # Panics
///
/// Panics if `modulus` is zero.
#[inline]
pub fn power_modulo_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64 % modulus;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = multiply_modulo_u64(result, base, modulus);
        }
        base = multiply_modulo_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

/// Invert a canonical residue modulo an assumed prime, or return `None` for zero.
/// Prefer [`ValidatedPrime::inverse`] at checked public boundaries.
///
/// # Panics
///
/// May panic if `prime < 2` and `value` is nonzero.
#[inline]
pub fn inverse_modulo_assumed_prime(value: u64, prime: u64) -> Option<u64> {
    (value != 0).then(|| power_modulo_u64(value, prime - 2, prime))
}

/// Reduce a signed machine integer modulo a nonzero modulus.
///
/// # Panics
///
/// Panics if `modulus` is zero.
#[inline]
pub fn signed_residue_u64(value: i64, modulus: u64) -> u64 {
    (i128::from(value).rem_euclid(i128::from(modulus))) as u64
}

/// Reduce an arbitrary integer to a word-size canonical residue.
///
/// # Panics
///
/// Panics if `modulus` is zero.
#[inline]
pub fn bigint_residue_u64(value: &BigInt, modulus: u64) -> u64 {
    let modulus_big = BigInt::from(modulus);
    let mut reduced = value % &modulus_big;
    if reduced.is_negative() {
        reduced += &modulus_big;
    }
    reduced.to_u64().expect("reduced residue should fit in u64")
}

/// Invert a residue using a validated prime; inputs are reduced first.
#[inline]
pub fn inverse_modulo_prime(value: u64, prime: ValidatedPrime) -> Option<u64> {
    inverse_modulo_assumed_prime(value % prime.get(), prime.get())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn arithmetic_near_word_limit() {
        assert_eq!(power_modulo_u64(0, 0, 1), 0);
        let p = 18_446_744_073_709_551_557;
        assert_eq!(add_modulo_u64(p - 1, p - 1, p), p - 2);
        assert_eq!(multiply_modulo_u64(p - 1, p - 1, p), 1);
        assert_eq!(subtract_modulo_u64(0, p - 1, p), 1);
        assert_eq!(signed_residue_u64(i64::MIN, p), p - (1u64 << 63));
        assert_eq!(bigint_residue_u64(&BigInt::from(-1), p), p - 1);
        assert_eq!(
            inverse_modulo_prime(0, ValidatedPrime::new(p).unwrap()),
            None
        );
        assert_eq!(
            inverse_modulo_prime(p + 1, ValidatedPrime::new(p).unwrap()),
            Some(1)
        );
    }
    #[test]
    #[should_panic]
    fn zero_modulus_is_rejected() {
        multiply_modulo_u64(1, 1, 0);
    }
}
