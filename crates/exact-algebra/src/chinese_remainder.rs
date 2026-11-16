//! Arbitrary-size Chinese remaindering with explicit reconstruction bounds.
//!
//! Bounds are caller-supplied mathematical premises. A modulus exceeding the
//! uniqueness threshold certifies reconstruction only when the target satisfies
//! those premises; it does not certify correctness of the modular computations.

use crate::modular_arithmetic::{
    bigint_residue_u64 as bigint_mod_u64, multiply_modulo_u64 as mul_mod_u64,
    subtract_modulo_u64 as sub_mod_u64,
};
use crate::primes::ValidatedPrime;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Failure of a congruence update or a bounded reconstruction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChineseRemainderError {
    /// Moduli must be positive.
    InvalidModulus,
    /// A new modulus shares a factor with the accumulated modulus.
    NonCoprimeModuli,
    /// Residue or bound vector lengths disagree with the accumulator width.
    InvalidDimensions,
    /// Bounds must be nonnegative.
    InvalidBound,
    /// The accumulated modulus does not yet prove uniqueness.
    InsufficientModulus,
    /// The reconstructed representative does not satisfy its supplied bound.
    BoundViolation,
    /// A requested fixed-width conversion does not fit.
    Overflow,
}
impl std::fmt::Display for ChineseRemainderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidModulus => "modulus must be positive",
            Self::NonCoprimeModuli => "moduli are not coprime",
            Self::InvalidDimensions => "residue or bound vector has the wrong length",
            Self::InvalidBound => "reconstruction bound must be nonnegative",
            Self::InsufficientModulus => "accumulated modulus does not certify uniqueness",
            Self::BoundViolation => "reconstructed value exceeds its supplied bound",
            Self::Overflow => "reconstructed value does not fit the requested integer type",
        })
    }
}
impl std::error::Error for ChineseRemainderError {}

/// Batched congruences with a common arbitrary-size accumulated modulus.
/// Failed updates leave both residues and modulus unchanged.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChineseRemainderAccumulator {
    residues: Vec<BigInt>,
    modulus: BigInt,
}
impl ChineseRemainderAccumulator {
    /// Create `width` zero residues modulo one (no congruences yet).
    pub fn new(width: usize) -> Self {
        Self {
            residues: vec![BigInt::zero(); width],
            modulus: BigInt::one(),
        }
    }
    /// Canonical residues in `[0, modulus)`.
    pub fn residues(&self) -> &[BigInt] {
        &self.residues
    }
    /// Product of all moduli supplied so far; initially one.
    pub fn modulus(&self) -> &BigInt {
        &self.modulus
    }
    /// Combine word-size residues using a validated prime.
    /// Inputs are reduced modulo the prime, and a repeated prime is rejected.
    pub fn combine_prime(
        &mut self,
        values: &[u64],
        prime: ValidatedPrime,
    ) -> Result<(), ChineseRemainderError> {
        if values.len() != self.residues.len() {
            return Err(ChineseRemainderError::InvalidDimensions);
        }
        let p = prime.get();
        let inverse = prime
            .inverse(bigint_mod_u64(&self.modulus, p))
            .ok_or(ChineseRemainderError::NonCoprimeModuli)?;
        for (residue, &value) in self.residues.iter_mut().zip(values) {
            let correction = sub_mod_u64(value % p, bigint_mod_u64(residue, p), p);
            *residue += &self.modulus * BigInt::from(mul_mod_u64(correction, inverse, p));
        }
        self.modulus *= BigInt::from(p);
        Ok(())
    }
    /// Combine arbitrary integer residues modulo a new coprime positive modulus.
    pub fn combine(
        &mut self,
        values: &[BigInt],
        modulus: &BigInt,
    ) -> Result<(), ChineseRemainderError> {
        if values.len() != self.residues.len() {
            return Err(ChineseRemainderError::InvalidDimensions);
        }
        if modulus <= &BigInt::zero() {
            return Err(ChineseRemainderError::InvalidModulus);
        }
        let extended = self.modulus.extended_gcd(modulus);
        if !extended.gcd.is_one() {
            return Err(ChineseRemainderError::NonCoprimeModuli);
        }
        for (residue, value) in self.residues.iter_mut().zip(values) {
            let correction = ((value - &*residue) * &extended.x).mod_floor(modulus);
            *residue += &self.modulus * correction;
        }
        self.modulus *= modulus;
        Ok(())
    }
    fn validate_bounds(&self, bounds: &[BigInt]) -> Result<(), ChineseRemainderError> {
        if bounds.len() != self.residues.len() {
            return Err(ChineseRemainderError::InvalidDimensions);
        }
        if bounds.iter().any(Signed::is_negative) {
            return Err(ChineseRemainderError::InvalidBound);
        }
        Ok(())
    }
    /// Whether `modulus > 2 * bound` for every absolute-value bound.
    pub fn certifies_symmetric_bounds(
        &self,
        bounds: &[BigInt],
    ) -> Result<bool, ChineseRemainderError> {
        self.validate_bounds(bounds)?;
        Ok(bounds.iter().all(|bound| self.modulus > bound * 2))
    }
    /// Whether `modulus^2 > 4 * bound` for every squared absolute-value bound.
    /// This avoids rounding a square root when using Hadamard bounds.
    pub fn certifies_symmetric_squared_bounds(
        &self,
        bounds: &[BigInt],
    ) -> Result<bool, ChineseRemainderError> {
        self.validate_bounds(bounds)?;
        let squared = &self.modulus * &self.modulus;
        Ok(bounds.iter().all(|bound| squared > bound * 4))
    }
    /// Whether `modulus > bound` for every inclusive nonnegative upper bound.
    pub fn certifies_nonnegative_bounds(
        &self,
        bounds: &[BigInt],
    ) -> Result<bool, ChineseRemainderError> {
        self.validate_bounds(bounds)?;
        Ok(bounds.iter().all(|bound| &self.modulus > bound))
    }
    /// Return signed representatives in `(-modulus/2, modulus/2]` after certification.
    pub fn reconstruct_symmetric(
        &self,
        bounds: &[BigInt],
    ) -> Result<Vec<BigInt>, ChineseRemainderError> {
        if !self.certifies_symmetric_bounds(bounds)? {
            return Err(ChineseRemainderError::InsufficientModulus);
        }
        let values = self.symmetric_values();
        if values
            .iter()
            .zip(bounds)
            .any(|(value, bound)| value.abs() > *bound)
        {
            return Err(ChineseRemainderError::BoundViolation);
        }
        Ok(values)
    }
    /// Return signed representatives after certification from squared bounds.
    pub fn reconstruct_symmetric_with_squared_bounds(
        &self,
        bounds: &[BigInt],
    ) -> Result<Vec<BigInt>, ChineseRemainderError> {
        if !self.certifies_symmetric_squared_bounds(bounds)? {
            return Err(ChineseRemainderError::InsufficientModulus);
        }
        let values = self.symmetric_values();
        if values
            .iter()
            .zip(bounds)
            .any(|(value, bound)| value * value > *bound)
        {
            return Err(ChineseRemainderError::BoundViolation);
        }
        Ok(values)
    }
    fn symmetric_values(&self) -> Vec<BigInt> {
        self.residues
            .iter()
            .map(|value| {
                if value * 2 > self.modulus {
                    value - &self.modulus
                } else {
                    value.clone()
                }
            })
            .collect()
    }
    /// Return values in `[0, bound]` once the modulus exceeds every inclusive bound.
    pub fn reconstruct_nonnegative(
        &self,
        bounds: &[BigInt],
    ) -> Result<Vec<BigInt>, ChineseRemainderError> {
        if !self.certifies_nonnegative_bounds(bounds)? {
            return Err(ChineseRemainderError::InsufficientModulus);
        }
        if self
            .residues
            .iter()
            .zip(bounds)
            .any(|(value, bound)| value > bound)
        {
            return Err(ChineseRemainderError::BoundViolation);
        }
        Ok(self.residues.clone())
    }
    /// Reconstruct signed values and narrow explicitly, reporting overflow separately.
    pub fn reconstruct_symmetric_i128(
        &self,
        bounds: &[BigInt],
    ) -> Result<Vec<i128>, ChineseRemainderError> {
        self.reconstruct_symmetric(bounds)?
            .iter()
            .map(|value| value.to_i128().ok_or(ChineseRemainderError::Overflow))
            .collect()
    }
}

/// Normalize an arbitrary residue into `(-modulus/2, modulus/2]`; the even tie is positive.
pub fn symmetric_residue(
    residue: &BigInt,
    modulus: &BigInt,
) -> Result<BigInt, ChineseRemainderError> {
    if modulus <= &BigInt::zero() {
        return Err(ChineseRemainderError::InvalidModulus);
    }
    let value = residue.mod_floor(modulus);
    Ok(if &value * 2 > *modulus {
        value - modulus
    } else {
        value
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use ChineseRemainderError::*;
    #[test]
    fn updates_bounds_and_atomic_errors() {
        let mut crt = ChineseRemainderAccumulator::new(2);
        crt.combine(&[(-17).into(), 99.into()], &1.into()).unwrap();
        assert_eq!(crt.modulus(), &BigInt::one());
        assert_eq!(crt.residues(), &[BigInt::zero(), BigInt::zero()]);
        crt.combine(&[(-2).into(), 10.into()], &5.into()).unwrap();
        assert_eq!(
            crt.reconstruct_symmetric(&[2.into(), 10.into()]),
            Err(InsufficientModulus)
        );
        let saved = crt.clone();
        assert_eq!(
            crt.combine(&[0.into(), 0.into()], &10.into()),
            Err(NonCoprimeModuli)
        );
        assert_eq!(crt, saved);
        assert_eq!(crt.combine(&[], &7.into()), Err(InvalidDimensions));
        assert_eq!(
            crt.combine(&[0.into(), 0.into()], &0.into()),
            Err(InvalidModulus)
        );
        crt.combine_prime(&[5, 3], ValidatedPrime::new(7).unwrap())
            .unwrap();
        assert_eq!(
            crt.reconstruct_symmetric(&[2.into(), 10.into()]).unwrap(),
            vec![BigInt::from(-2), BigInt::from(10)]
        );
        assert_eq!(
            crt.reconstruct_symmetric_with_squared_bounds(&[4.into(), 100.into()])
                .unwrap(),
            vec![BigInt::from(-2), BigInt::from(10)]
        );
        assert_eq!(
            crt.reconstruct_nonnegative(&[33.into(), 10.into()])
                .unwrap(),
            vec![BigInt::from(33), BigInt::from(10)]
        );
        assert_eq!(
            crt.reconstruct_nonnegative(&[2.into(), 10.into()]),
            Err(BoundViolation)
        );
        assert_eq!(
            crt.reconstruct_symmetric(&[1.into(), 10.into()]),
            Err(BoundViolation)
        );
        assert_eq!(
            crt.reconstruct_symmetric_with_squared_bounds(&[1.into(), 100.into()]),
            Err(BoundViolation)
        );
        assert_eq!(
            crt.certifies_symmetric_bounds(&[(-1).into(), 0.into()]),
            Err(InvalidBound)
        );
        assert_eq!(
            crt.certifies_nonnegative_bounds(&[]),
            Err(InvalidDimensions)
        );
        assert_eq!(
            crt.certifies_symmetric_squared_bounds(&[(-1).into(), 0.into()]),
            Err(InvalidBound)
        );
        let saved = crt.clone();
        assert_eq!(
            crt.combine_prime(&[0, 0], ValidatedPrime::new(7).unwrap()),
            Err(NonCoprimeModuli)
        );
        assert_eq!(crt, saved);
    }
    #[test]
    fn strict_thresholds_and_even_tie() {
        let mut crt = ChineseRemainderAccumulator::new(1);
        crt.combine(&[3.into()], &6.into()).unwrap();
        assert!(!crt.certifies_symmetric_bounds(&[3.into()]).unwrap());
        assert!(!crt.certifies_symmetric_squared_bounds(&[9.into()]).unwrap());
        assert!(!crt.certifies_nonnegative_bounds(&[6.into()]).unwrap());
        assert_eq!(
            symmetric_residue(&(-3).into(), &6.into()).unwrap(),
            BigInt::from(3)
        );
        assert!(symmetric_residue(&0.into(), &0.into()).is_err());
    }
    #[test]
    fn prime_product_exceeds_i128_without_false_overflow() {
        let mut reconstruction = ChineseRemainderAccumulator::new(1);
        for modulus in [18_446_744_073_709_551_557, 18_446_744_073_709_551_533] {
            let prime = ValidatedPrime::new(modulus).unwrap();
            reconstruction.combine_prime(&[modulus - 5], prime).unwrap();
        }
        assert!(reconstruction.modulus() > &BigInt::from(i128::MAX));
        assert_eq!(
            reconstruction.reconstruct_symmetric(&[5.into()]).unwrap(),
            vec![BigInt::from(-5)]
        );
        assert_eq!(
            reconstruction
                .reconstruct_symmetric_i128(&[5.into()])
                .unwrap(),
            vec![-5i128]
        );
    }

    #[test]
    fn arbitrary_size_and_explicit_overflow() {
        let value: BigInt = BigInt::one() << 200;
        let bound = value.clone();
        let modulus = (&value << 2) + 1;
        let mut crt = ChineseRemainderAccumulator::new(1);
        crt.combine(std::slice::from_ref(&value), &modulus).unwrap();
        assert_eq!(
            crt.reconstruct_symmetric(std::slice::from_ref(&bound))
                .unwrap(),
            vec![value]
        );
        assert_eq!(crt.reconstruct_symmetric_i128(&[bound]), Err(Overflow));
    }
}
