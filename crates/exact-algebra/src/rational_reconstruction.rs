//! Bounded rational reconstruction over arbitrary-size integers.
//!
//! A successful result has coprime numerator and positive denominator, satisfies
//! both bounds and the modular congruence, and has a denominator coprime to the modulus.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

/// Failure to certify a unique bounded rational reconstruction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RationalReconstructionError {
    /// The modulus must be greater than one.
    InvalidModulus,
    /// Numerator bound must be nonnegative and denominator bound positive.
    InvalidBounds,
    /// The uniqueness inequality `2 * numerator_bound * denominator_bound < modulus` fails.
    InsufficientModulus,
    /// No fraction satisfies the bounds, coprimality and final congruence checks.
    NoReconstruction,
}
impl std::fmt::Display for RationalReconstructionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::InvalidModulus => "modulus must be greater than one",
            Self::InvalidBounds => "invalid numerator or denominator bound",
            Self::InsufficientModulus => "modulus does not certify uniqueness",
            Self::NoReconstruction => "no bounded rational reconstruction exists",
        })
    }
}
impl std::error::Error for RationalReconstructionError {}

/// Reconstruct `a/b` from `residue = a/b (mod modulus)` with `|a| <= numerator_bound`
/// and `1 <= b <= denominator_bound`, requiring the strict uniqueness inequality.
pub fn rational_reconstruction_with_bounds(
    residue: &BigInt,
    modulus: &BigInt,
    numerator_bound: &BigInt,
    denominator_bound: &BigInt,
) -> Result<BigRational, RationalReconstructionError> {
    use RationalReconstructionError::*;
    if modulus <= &BigInt::one() {
        return Err(InvalidModulus);
    }
    if numerator_bound.is_negative() || denominator_bound <= &BigInt::zero() {
        return Err(InvalidBounds);
    }
    if numerator_bound * denominator_bound * 2 >= *modulus {
        return Err(InsufficientModulus);
    }
    let target = residue.mod_floor(modulus);
    let (mut r0, mut r1) = (modulus.clone(), target);
    let (mut s0, mut s1) = (BigInt::zero(), BigInt::one());
    while r1.abs() > *numerator_bound {
        let quotient = &r0 / &r1;
        (r0, r1) = (r1.clone(), &r0 - &quotient * &r1);
        (s0, s1) = (s1.clone(), &s0 - &quotient * &s1);
    }
    let (mut numerator, mut denominator) = (r1, s1);
    if denominator.is_negative() {
        numerator = -numerator;
        denominator = -denominator;
    }
    if denominator.is_zero()
        || denominator > *denominator_bound
        || numerator.abs() > *numerator_bound
        || !numerator.gcd(&denominator).is_one()
        || !denominator.gcd(modulus).is_one()
        || !(&denominator * residue - &numerator)
            .mod_floor(modulus)
            .is_zero()
    {
        return Err(NoReconstruction);
    }
    Ok(BigRational::new(numerator, denominator))
}

/// Reconstruct using the symmetric bound `floor(sqrt((modulus-1)/2))` for both parts.
pub fn rational_reconstruction(
    residue: &BigInt,
    modulus: &BigInt,
) -> Result<BigRational, RationalReconstructionError> {
    if modulus <= &BigInt::one() {
        return Err(RationalReconstructionError::InvalidModulus);
    }
    let bound = ((modulus - 1u8) / 2u8).sqrt();
    // Modulus two admits only zero under the symmetric bound; use denominator one.
    let denominator_bound = bound.clone().max(BigInt::one());
    rational_reconstruction_with_bounds(residue, modulus, &bound, &denominator_bound)
}

#[cfg(test)]
mod tests {
    use super::*;
    use RationalReconstructionError::*;
    #[test]
    fn bounded_positive_negative_zero_and_failures() {
        let modulus = BigInt::from(101);
        let inverse = BigInt::from(3).extended_gcd(&modulus).x;
        for n in [-4, 0, 4] {
            let residue = BigInt::from(n) * &inverse;
            assert_eq!(
                rational_reconstruction(&residue, &modulus).unwrap(),
                BigRational::new(n.into(), 3.into())
            );
        }
        assert_eq!(
            rational_reconstruction(&0.into(), &1.into()),
            Err(InvalidModulus)
        );
        assert_eq!(
            rational_reconstruction_with_bounds(&1.into(), &101.into(), &(-1).into(), &1.into()),
            Err(InvalidBounds)
        );
        assert_eq!(
            rational_reconstruction_with_bounds(&1.into(), &101.into(), &1.into(), &0.into()),
            Err(InvalidBounds)
        );
        assert_eq!(
            rational_reconstruction_with_bounds(&1.into(), &100.into(), &5.into(), &10.into()),
            Err(InsufficientModulus)
        );
        assert_eq!(
            rational_reconstruction_with_bounds(&50.into(), &101.into(), &1.into(), &1.into()),
            Err(NoReconstruction)
        );
        assert_eq!(
            rational_reconstruction_with_bounds(&3.into(), &9.into(), &1.into(), &3.into()),
            Err(NoReconstruction)
        );
    }
    #[test]
    fn reconstructs_beyond_machine_integer_bounds() {
        let numerator: BigInt = -(BigInt::one() << 180usize) + 1;
        let denominator: BigInt = (BigInt::one() << 127usize) + 3;
        let modulus: BigInt = BigInt::one() << 400usize;
        let inverse = denominator.extended_gcd(&modulus).x;
        let residue = &numerator * inverse;
        assert_eq!(
            rational_reconstruction_with_bounds(&residue, &modulus, &numerator.abs(), &denominator)
                .unwrap(),
            BigRational::new(numerator, denominator)
        );
    }
    #[test]
    fn asymmetric_bounds_match_exhaustive_oracle() {
        for m in 2..70i64 {
            for nb in 0..=4i64 {
                for db in 1..=5i64 {
                    if 2 * nb * db >= m {
                        continue;
                    }
                    for residue in 0..m {
                        let expected = (-nb..=nb)
                            .flat_map(|n| (1..=db).map(move |d| (n, d)))
                            .find(|&(n, d)| {
                                n.gcd(&d) == 1
                                    && d.gcd(&m) == 1
                                    && (d * residue - n).rem_euclid(m) == 0
                            })
                            .map(|(n, d)| BigRational::new(n.into(), d.into()));
                        assert_eq!(
                            rational_reconstruction_with_bounds(
                                &residue.into(),
                                &m.into(),
                                &nb.into(),
                                &db.into()
                            )
                            .ok(),
                            expected,
                            "{residue} mod {m}, bounds {nb}, {db}"
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn exhaustive_small_bounded_fractions() {
        for modulus in 3..160i64 {
            let m = BigInt::from(modulus);
            let b = ((&m - 1u8) / 2u8).sqrt();
            for den in 1..=8i64 {
                for num in -8..=8i64 {
                    if BigInt::from(num).abs() > b
                        || BigInt::from(den) > b
                        || num.gcd(&den) != 1
                        || den.gcd(&modulus) != 1
                    {
                        continue;
                    }
                    let inv = BigInt::from(den).extended_gcd(&m).x;
                    let residue = BigInt::from(num) * inv;
                    assert_eq!(
                        rational_reconstruction(&residue, &m).unwrap(),
                        BigRational::new(num.into(), den.into())
                    );
                }
            }
        }
    }
}
