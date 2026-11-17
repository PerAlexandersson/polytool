//! Differential checks for shared signed residue arithmetic.
use super::*;
use num_traits::ToPrimitive;
fn original_mod_norm(value: i64, modulus: i64) -> i64 {
    value.rem_euclid(modulus)
}

fn original_mod_mul(lhs: i64, rhs: i64, modulus: i64) -> i64 {
    if (0..modulus).contains(&lhs) && (0..modulus).contains(&rhs) && modulus <= 3_000_000_000 {
        return (lhs * rhs) % modulus;
    }
    ((lhs as i128 * rhs as i128).rem_euclid(modulus as i128)) as i64
}

fn original_mod_pow(mut base: i64, mut exp: i64, modulus: i64) -> i64 {
    let mut acc = 1;
    base = original_mod_norm(base, modulus);
    while exp > 0 {
        if exp & 1 == 1 {
            acc = original_mod_mul(acc, base, modulus);
        }
        base = original_mod_mul(base, base, modulus);
        exp >>= 1;
    }
    acc
}

fn original_mod_inv(value: i64, modulus: i64) -> Option<i64> {
    let value = original_mod_norm(value, modulus);
    (value != 0).then(|| original_mod_pow(value, modulus - 2, modulus))
}

fn original_bigint_mod_i64(value: &BigInt, modulus: i64) -> i64 {
    let modulus_big = BigInt::from(modulus);
    let residue = ((value % &modulus_big) + &modulus_big) % &modulus_big;
    residue
        .to_i64()
        .expect("residue modulo i64 prime fits in i64")
}
#[test]
fn shared_residues_match_former_recurrence_arithmetic() {
    let mut seed = 0xa198_0f32_76c4_534du64;
    for prime in MODULAR_PREFILTER_PRIMES {
        assert!(exact_algebra::primes::ValidatedPrime::new(prime as u64).is_ok());
        for _ in 0..1000 {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let lhs = seed as i64;
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            let rhs = seed as i64;
            assert_eq!(mod_mul(lhs, rhs, prime), original_mod_mul(lhs, rhs, prime));
            assert_eq!(mod_inv(lhs, prime), original_mod_inv(lhs, prime));
            let big = (BigInt::from(lhs) << 180usize) + BigInt::from(rhs);
            assert_eq!(
                bigint_mod_i64(&big, prime),
                original_bigint_mod_i64(&big, prime)
            );
        }
        for value in [i64::MIN, -prime, -1, 0, 1, prime, i64::MAX] {
            assert_eq!(mod_inv(value, prime), original_mod_inv(value, prime));
        }
    }
    let mut primes = MODULAR_PREFILTER_PRIMES.to_vec();
    primes.sort_unstable();
    primes.dedup();
    assert_eq!(primes.len(), MODULAR_PREFILTER_PRIMES.len());
}
