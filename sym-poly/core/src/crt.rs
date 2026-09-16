use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::Ratio;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Combine two congruences with coprime moduli.
///
/// Returns the unique residue modulo `m1 * m2`, normalized to
/// `0 <= residue < modulus`, or `None` if the combined modulus does not fit
/// in `i128`.
pub fn chinese_remainder_pair(
    residue1: i128,
    modulus1: i128,
    residue2: i128,
    modulus2: i128,
) -> Option<(i128, i128)> {
    assert!(modulus1 > 0, "first modulus must be positive");
    assert!(modulus2 > 0, "second modulus must be positive");
    if modulus1.gcd(&modulus2) != 1 {
        return None;
    }

    let modulus1 = BigInt::from(modulus1);
    let modulus2 = BigInt::from(modulus2);
    let r1 = residue_mod_big(&BigInt::from(residue1), &modulus1);
    let r2 = residue_mod_big(&BigInt::from(residue2), &modulus2);
    let inverse = modular_inverse_big(&modulus1, &modulus2)?;
    let adjustment = residue_mod_big(&((r2 - &r1) * inverse), &modulus2);
    let modulus = &modulus1 * &modulus2;
    let residue = residue_mod_big(&(r1 + &modulus1 * adjustment), &modulus);

    Some((residue.to_i128()?, modulus.to_i128()?))
}

/// Combine a nonempty list of pairwise coprime congruences.
pub fn chinese_remainder(congruences: &[(i128, i128)]) -> Option<(i128, i128)> {
    let (&(first_residue, first_modulus), rest) = congruences.split_first()?;
    assert!(first_modulus > 0, "first modulus must be positive");
    let mut combined = (residue_mod(first_residue, first_modulus), first_modulus);
    for &(residue, modulus) in rest {
        combined = chinese_remainder_pair(combined.0, combined.1, residue, modulus)?;
    }
    Some(combined)
}

/// Return the representative in `[-modulus/2, modulus/2]`.
pub fn symmetric_residue(residue: i128, modulus: i128) -> i128 {
    assert!(modulus > 0, "modulus must be positive");
    let normalized = residue_mod(residue, modulus);
    if normalized > modulus - normalized {
        normalized - modulus
    } else {
        normalized
    }
}

/// Rational reconstruction from a residue modulo `modulus`.
///
/// If `residue = a / b (mod modulus)` and `|a|, b <= sqrt(modulus / 2)`,
/// this returns `a / b`. Otherwise it returns `None`.
pub fn rational_reconstruction(residue: i128, modulus: i128) -> Option<Ratio<i128>> {
    assert!(modulus > 1, "modulus must be greater than 1");
    let modulus = BigInt::from(modulus);
    let bound = ((&modulus - BigInt::one()) / BigInt::from(2)).sqrt();
    let target = symmetric_residue_big(&BigInt::from(residue), &modulus);

    let (mut r0, mut r1) = (modulus.clone(), target);
    let (mut s0, mut s1) = (BigInt::zero(), BigInt::one());

    while r1.abs() > bound {
        if r1.is_zero() {
            return None;
        }
        let q = &r0 / &r1;
        let next_r = &r0 - &q * &r1;
        let next_s = &s0 - &q * &s1;
        (r0, r1) = (r1, next_r);
        (s0, s1) = (s1, next_s);
    }

    let mut numerator = r1;
    let mut denominator = s1;
    if denominator.is_negative() {
        numerator = -numerator;
        denominator = -denominator;
    }

    if denominator.is_zero() || numerator.abs() > bound || denominator > bound {
        return None;
    }
    if numerator.gcd(&denominator) != BigInt::one() {
        return None;
    }
    if !residue_mod_big(&(denominator.clone() * residue - &numerator), &modulus).is_zero() {
        return None;
    }

    Some(Ratio::new(numerator.to_i128()?, denominator.to_i128()?))
}

fn modular_inverse_big(value: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    let (gcd, inverse, _) = extended_gcd_big(value.clone(), modulus.clone());
    (gcd == BigInt::one()).then(|| residue_mod_big(&inverse, modulus))
}

#[cfg(test)]
fn modular_inverse(value: i128, modulus: i128) -> Option<i128> {
    modular_inverse_big(&BigInt::from(value), &BigInt::from(modulus))?.to_i128()
}

fn residue_mod(value: i128, modulus: i128) -> i128 {
    residue_mod_big(&BigInt::from(value), &BigInt::from(modulus))
        .to_i128()
        .expect("an i128 residue modulo a positive i128 modulus fits in i128")
}

fn residue_mod_big(value: &BigInt, modulus: &BigInt) -> BigInt {
    value.mod_floor(modulus)
}

fn symmetric_residue_big(residue: &BigInt, modulus: &BigInt) -> BigInt {
    let normalized = residue_mod_big(residue, modulus);
    if normalized > modulus - &normalized {
        normalized - modulus
    } else {
        normalized
    }
}

fn extended_gcd_big(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
    let (mut old_r, mut r) = (a, b);
    let (mut old_s, mut s) = (BigInt::one(), BigInt::zero());
    let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());

    while !r.is_zero() {
        let quotient = &old_r / &r;
        let next_r = &old_r - &quotient * &r;
        let next_s = &old_s - &quotient * &s;
        let next_t = &old_t - &quotient * &t;
        (old_r, r) = (r, next_r);
        (old_s, s) = (s, next_s);
        (old_t, t) = (t, next_t);
    }

    if old_r.is_negative() {
        (-old_r, -old_s, -old_t)
    } else {
        (old_r, old_s, old_t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chinese_remainder_pair() {
        let (residue, modulus) = chinese_remainder_pair(2, 5, 3, 7).unwrap();

        assert_eq!(modulus, 35);
        assert_eq!(residue % 5, 2);
        assert_eq!(residue % 7, 3);
        assert_eq!(residue, 17);
    }

    #[test]
    fn test_chinese_remainder_sequence() {
        let (residue, modulus) = chinese_remainder(&[(1, 5), (2, 7), (3, 11)]).unwrap();

        assert_eq!(modulus, 385);
        assert_eq!(residue % 5, 1);
        assert_eq!(residue % 7, 2);
        assert_eq!(residue % 11, 3);
    }

    #[test]
    fn test_symmetric_residue() {
        assert_eq!(symmetric_residue(98, 101), -3);
        assert_eq!(symmetric_residue(49, 101), 49);
    }

    #[test]
    fn test_symmetric_residue_near_i128_limit() {
        assert_eq!(symmetric_residue(i128::MAX - 1, i128::MAX), -1);
    }

    #[test]
    fn test_chinese_remainder_rejects_unrepresentable_modulus() {
        assert_eq!(chinese_remainder_pair(0, i128::MAX, 1, 2), None);
    }

    #[test]
    fn test_rational_reconstruction() {
        let modulus = 1009;
        let inverse_7 = modular_inverse(7, modulus).unwrap();
        let residue = (5 * inverse_7) % modulus;

        assert_eq!(
            rational_reconstruction(residue, modulus),
            Some(Ratio::new(5, 7))
        );
    }

    #[test]
    fn test_rational_reconstruction_near_i128_limit() {
        let modulus = i128::MAX;
        // This is 5 / 2 modulo 2^127 - 1, without forming an overflowing
        // product with the modular inverse of 2.
        let residue = modulus / 2 + 3;

        assert_eq!(
            rational_reconstruction(residue, modulus),
            Some(Ratio::new(5, 2))
        );
    }
}
