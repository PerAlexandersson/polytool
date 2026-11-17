//! Rank-profile agreement and independent-prime verification for CRT outputs.
//!
//! Agreement across primes is a safeguard, not a characteristic-zero rank proof.
//! Callers supply the coefficient bound and a profile covering every elimination
//! that affects the residues. A disagreement is returned explicitly.

use crate::chinese_remainder::{ChineseRemainderAccumulator, ChineseRemainderError};
use crate::modular_arithmetic::bigint_residue_u64;
use crate::primes::ValidatedPrime;
use num_bigint::BigInt;
use num_traits::Zero;

/// Failure of profile-checked modular reconstruction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ProfiledReconstructionError {
    /// Two modular computations used different pivot or free-column profiles.
    RankProfileMismatch,
    /// The reconstructed coefficients disagree with the independent spare prime.
    SparePrimeMismatch,
    /// Congruence dimensions, bounds, or moduli are invalid.
    ChineseRemainder(ChineseRemainderError),
}

impl std::fmt::Display for ProfiledReconstructionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RankProfileMismatch => f.write_str("modular rank profiles disagree"),
            Self::SparePrimeMismatch => f.write_str("independent spare-prime verification failed"),
            Self::ChineseRemainder(error) => error.fmt(f),
        }
    }
}
impl std::error::Error for ProfiledReconstructionError {}
impl From<ChineseRemainderError> for ProfiledReconstructionError {
    fn from(error: ChineseRemainderError) -> Self {
        Self::ChineseRemainder(error)
    }
}

/// CRT reconstruction requiring two agreeing profiles and a distinct spare prime.
///
/// Profiles consist of caller-defined lists of zero-based pivot/free columns,
/// including empty lists and ambient dimensions when needed. Every list must
/// agree exactly. A mathematical absolute-value bound remains a caller premise.
#[derive(Clone, Debug)]
pub struct ProfiledChineseRemainder {
    accumulator: ChineseRemainderAccumulator,
    bounds: Vec<BigInt>,
    profile: Option<Vec<Vec<usize>>>,
    combined_primes: usize,
    candidate: Option<Vec<BigInt>>,
}

impl ProfiledChineseRemainder {
    /// Start a reconstruction with an absolute-value bound for each coefficient.
    pub fn new(bounds: Vec<BigInt>) -> Result<Self, ChineseRemainderError> {
        let accumulator = ChineseRemainderAccumulator::new(bounds.len());
        accumulator.certifies_symmetric_bounds(&bounds)?;
        Ok(Self {
            accumulator,
            bounds,
            profile: None,
            combined_primes: 0,
            candidate: None,
        })
    }

    /// Observe one prime, returning coefficients only after spare-prime agreement.
    ///
    /// The spare prime is never combined into the candidate being verified.
    /// On any error discard the reconstruction; no rank is inferred from one prime.
    pub fn push_prime(
        &mut self,
        values: &[u64],
        profile: &[Vec<usize>],
        prime: ValidatedPrime,
    ) -> Result<Option<Vec<BigInt>>, ProfiledReconstructionError> {
        if values.len() != self.bounds.len() {
            return Err(ChineseRemainderError::InvalidDimensions.into());
        }
        if self
            .profile
            .as_deref()
            .is_some_and(|expected| expected != profile)
        {
            return Err(ProfiledReconstructionError::RankProfileMismatch);
        }
        if bigint_residue_u64(self.accumulator.modulus(), prime.get()).is_zero() {
            return Err(ChineseRemainderError::NonCoprimeModuli.into());
        }
        if let Some(candidate) = &self.candidate {
            if candidate.iter().zip(values).all(|(value, residue)| {
                bigint_residue_u64(value, prime.get()) == residue % prime.get()
            }) {
                return Ok(Some(candidate.clone()));
            }
            return Err(ProfiledReconstructionError::SparePrimeMismatch);
        }
        self.accumulator.combine_prime(values, prime)?;
        self.profile = Some(profile.to_vec());
        self.combined_primes += 1;
        if self.combined_primes >= 2 && self.accumulator.certifies_symmetric_bounds(&self.bounds)? {
            self.candidate = Some(self.accumulator.reconstruct_symmetric(&self.bounds)?);
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn needs_two_profiles_and_an_independent_spare() {
        let mut lift = ProfiledChineseRemainder::new(vec![10.into()]).unwrap();
        let profile = vec![vec![0]];
        for p in [101, 103] {
            assert_eq!(
                lift.push_prime(&[p - 3], &profile, ValidatedPrime::new(p).unwrap()),
                Ok(None)
            );
        }
        assert_eq!(
            lift.push_prime(&[104], &profile, ValidatedPrime::new(107).unwrap()),
            Ok(Some(vec![(-3).into()]))
        );
    }
    #[test]
    fn unlucky_rank_and_spare_are_rejected() {
        let mut lift = ProfiledChineseRemainder::new(vec![10.into()]).unwrap();
        lift.push_prime(&[0], &[vec![]], ValidatedPrime::new(7).unwrap())
            .unwrap();
        assert_eq!(
            lift.push_prime(&[0], &[vec![0]], ValidatedPrime::new(11).unwrap()),
            Err(ProfiledReconstructionError::RankProfileMismatch)
        );
        let mut lift = ProfiledChineseRemainder::new(vec![10.into()]).unwrap();
        for p in [101, 103] {
            lift.push_prime(&[1], &[], ValidatedPrime::new(p).unwrap())
                .unwrap();
        }
        assert_eq!(
            lift.push_prime(&[2], &[], ValidatedPrime::new(107).unwrap()),
            Err(ProfiledReconstructionError::SparePrimeMismatch)
        );
        assert!(ProfiledChineseRemainder::new(vec![(-1).into()]).is_err());
    }

    #[test]
    fn dimensions_and_reused_spare_primes_are_rejected() {
        let mut lift = ProfiledChineseRemainder::new(vec![10.into()]).unwrap();
        assert_eq!(
            lift.push_prime(&[], &[], ValidatedPrime::new(101).unwrap()),
            Err(ChineseRemainderError::InvalidDimensions.into())
        );
        for p in [101, 103] {
            assert_eq!(
                lift.push_prime(&[1], &[], ValidatedPrime::new(p).unwrap()),
                Ok(None)
            );
        }
        assert_eq!(
            lift.push_prime(&[1], &[], ValidatedPrime::new(101).unwrap()),
            Err(ChineseRemainderError::NonCoprimeModuli.into())
        );
    }
}
