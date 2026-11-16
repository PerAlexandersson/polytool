//! Deterministic primality checks and descending streams of validated primes.

use crate::modular_arithmetic::{
    multiply_modulo_u64 as mul_mod_u64, power_modulo_u64 as pow_mod_u64,
};

/// A rejected modulus that is not prime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InvalidPrime {
    /// The rejected modulus.
    pub modulus: u64,
}
impl std::fmt::Display for InvalidPrime {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "modulus {} is not prime", self.modulus)
    }
}
impl std::error::Error for InvalidPrime {}

/// A word-size prime checked by deterministic primality testing.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ValidatedPrime(u64);
impl ValidatedPrime {
    /// Check a modulus for primality.
    pub fn new(modulus: u64) -> Result<Self, InvalidPrime> {
        if is_prime_u64(modulus) {
            Ok(Self(modulus))
        } else {
            Err(InvalidPrime { modulus })
        }
    }
    /// Return the prime modulus.
    #[inline]
    pub const fn get(self) -> u64 {
        self.0
    }
    /// Invert a residue after reduction; zero has no inverse.
    #[inline]
    pub fn inverse(self, value: u64) -> Option<u64> {
        crate::modular_arithmetic::inverse_modulo_prime(value, self)
    }
}

/// A descending iterator over all primes at or below its starting bound.
#[derive(Clone, Debug)]
pub struct DescendingPrimes {
    next_bound: Option<u64>,
}
impl DescendingPrimes {
    /// Start at the greatest prime no larger than `upper_bound`.
    pub fn new(upper_bound: u64) -> Self {
        Self {
            next_bound: Some(upper_bound),
        }
    }
}
impl Iterator for DescendingPrimes {
    type Item = ValidatedPrime;
    fn next(&mut self) -> Option<Self::Item> {
        let prime = previous_prime_at_or_below(self.next_bound?)?;
        self.next_bound = prime.checked_sub(1);
        Some(ValidatedPrime(prime))
    }
}
impl std::iter::FusedIterator for DescendingPrimes {}

fn previous_prime_at_or_below(mut n: u64) -> Option<u64> {
    if n < 2 {
        return None;
    }
    if n == 2 {
        return Some(2);
    }
    if n & 1 == 0 {
        n -= 1;
    }
    while n >= 3 {
        if is_prime_u64(n) {
            return Some(n);
        }
        n = n.saturating_sub(2);
    }
    None
}

/// Deterministic primality test for every `u64`, using twelve Miller-Rabin bases.
pub fn is_prime_u64(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for p in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n == p {
            return true;
        }
        if n % p == 0 {
            return false;
        }
    }

    let mut d = n - 1;
    let mut s = 0;
    while d & 1 == 0 {
        d /= 2;
        s += 1;
    }

    for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if a >= n {
            continue;
        }
        let mut x = pow_mod_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        let mut probably_prime_for_base = false;
        for _ in 1..s {
            x = mul_mod_u64(x, x, n);
            if x == n - 1 {
                probably_prime_for_base = true;
                break;
            }
        }
        if !probably_prime_for_base {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn primes_and_pseudoprimes() {
        for n in [
            2,
            3,
            37,
            2_305_843_009_213_693_951,
            18_446_744_073_709_551_557,
        ] {
            assert!(is_prime_u64(n));
        }
        for n in [
            0,
            1,
            4,
            561,
            3_215_031_751,
            341_550_071_728_321,
            3_825_123_056_546_413_051,
            u64::MAX,
        ] {
            assert!(!is_prime_u64(n));
            assert!(ValidatedPrime::new(n).is_err());
        }
    }
    #[test]
    fn descending_stream_terminates() {
        assert_eq!(
            DescendingPrimes::new(12)
                .map(ValidatedPrime::get)
                .collect::<Vec<_>>(),
            vec![11, 7, 5, 3, 2]
        );
        let mut empty = DescendingPrimes::new(1);
        assert_eq!(empty.next(), None);
        assert_eq!(empty.next(), None);
    }
}
