//! Modular arithmetic on `u64` residues and CRT reconstruction of leading
//! principal minors.
//!
//! [`modular_leading_principal_minors_bigint`]
//! combines determinants modulo word-sized primes until the Hadamard bounds
//! certify the reconstruction.  The private residue helpers are shared with the
//! sparse modular solver.

use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// Compute leading principal determinants by modular determinant computation
/// and Chinese remaindering.
///
/// For a square integer matrix, this returns the determinant of each leading
/// principal submatrix.  Each determinant is computed modulo a sequence of
/// large primes and reconstructed by CRT once the accumulated modulus exceeds
/// twice a Hadamard bound for every leading minor.
///
/// Unlike [`bareiss_leading_principal_minors_bigint`], this method can handle
/// zero Bareiss pivots because each leading determinant is computed separately
/// with modular pivoting.  It is intended as a comparison/large-entry path, not
/// as a replacement for the simpler Bareiss implementation.
///
/// [`bareiss_leading_principal_minors_bigint`]: super::bareiss_leading_principal_minors_bigint
pub fn modular_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    const MAX_MODULAR_PRIMES: usize = 256;

    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let hadamard_squared_bounds = leading_principal_hadamard_squared_bounds(mat)?;
    let mut residues = vec![BigInt::zero(); n];
    let mut modulus = BigInt::one();
    let mut prime_search_start = (1u64 << 61) - 1;

    for _ in 0..MAX_MODULAR_PRIMES {
        let prime = previous_prime_at_or_below(prime_search_start)?;
        prime_search_start = prime.saturating_sub(2);
        let modular_minors = leading_principal_minors_mod_prime(mat, prime);
        crt_update_residues(&mut residues, &mut modulus, &modular_minors, prime)?;

        if crt_modulus_certifies_bounds(&modulus, &hadamard_squared_bounds) {
            return Some(
                residues
                    .into_iter()
                    .map(|r| symmetric_residue(r, &modulus))
                    .collect(),
            );
        }
    }

    None
}

fn leading_principal_hadamard_squared_bounds(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    let mut bounds = Vec::with_capacity(n);
    for k in 1..=n {
        let mut product = BigInt::one();
        for row in &mat[..k] {
            let row_square_sum: BigInt = row[..k].iter().map(|entry| entry * entry).sum();
            product *= row_square_sum;
        }
        bounds.push(product);
    }
    Some(bounds)
}

fn crt_modulus_certifies_bounds(modulus: &BigInt, squared_bounds: &[BigInt]) -> bool {
    let modulus_squared = modulus * modulus;
    squared_bounds.iter().all(|bound| {
        let threshold = BigInt::from(4) * bound;
        modulus_squared > threshold
    })
}

fn symmetric_residue(residue: BigInt, modulus: &BigInt) -> BigInt {
    if &residue * BigInt::from(2) > *modulus {
        residue - modulus
    } else {
        residue
    }
}

fn crt_update_residues(
    residues: &mut [BigInt],
    modulus: &mut BigInt,
    modular_values: &[u64],
    prime: u64,
) -> Option<()> {
    let modulus_mod_prime = bigint_mod_u64(modulus, prime);
    let inverse = mod_inverse_prime(modulus_mod_prime, prime)?;
    let old_modulus = modulus.clone();

    for (residue, &value_mod_prime) in residues.iter_mut().zip(modular_values.iter()) {
        let residue_mod_prime = bigint_mod_u64(residue, prime);
        let correction = if value_mod_prime >= residue_mod_prime {
            value_mod_prime - residue_mod_prime
        } else {
            prime - (residue_mod_prime - value_mod_prime)
        };
        let t = mul_mod_u64(correction, inverse, prime);
        *residue += &old_modulus * BigInt::from(t);
    }

    *modulus = old_modulus * BigInt::from(prime);
    Some(())
}

fn bigint_mod_u64(value: &BigInt, modulus: u64) -> u64 {
    let modulus_big = BigInt::from(modulus);
    let mut reduced = value % &modulus_big;
    if reduced.is_negative() {
        reduced += &modulus_big;
    }
    reduced.to_u64().expect("reduced residue should fit in u64")
}

fn leading_principal_minors_mod_prime(mat: &[Vec<BigInt>], prime: u64) -> Vec<u64> {
    let reduced: Vec<Vec<u64>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| bigint_mod_u64(entry, prime))
                .collect()
        })
        .collect();
    if let Some(minors) = leading_principal_minors_mod_prime_bareiss(&reduced, prime) {
        return minors;
    }

    let n = mat.len();
    (1..=n)
        .map(|prefix_size| determinant_mod_prime_prefix(&reduced, prefix_size, prime))
        .collect()
}

fn leading_principal_minors_mod_prime_bareiss(mat: &[Vec<u64>], prime: u64) -> Option<Vec<u64>> {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut denom = 1u64;
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k];
        if pivot == 0 {
            if k == n - 1 {
                minors.push(0);
                break;
            }
            return None;
        }
        minors.push(pivot);
        if k == n - 1 {
            break;
        }

        let denom_inverse = mod_inverse_prime(denom, prime)?;
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let left = mul_mod_u64(a[i][j], pivot, prime);
                let right = mul_mod_u64(a[i][k], a[k][j], prime);
                let numerator = sub_mod_u64(left, right, prime);
                a[i][j] = mul_mod_u64(numerator, denom_inverse, prime);
            }
        }
        denom = pivot;
    }

    Some(minors)
}

fn determinant_mod_prime_prefix(mat: &[Vec<u64>], size: usize, prime: u64) -> u64 {
    let mut a: Vec<Vec<u64>> = mat
        .iter()
        .take(size)
        .map(|row| row.iter().take(size).copied().collect())
        .collect();
    let mut det = 1u64;
    let mut sign_is_negative = false;

    for k in 0..size {
        let Some(pivot_row) = (k..size).find(|&i| a[i][k] != 0) else {
            return 0;
        };
        if pivot_row != k {
            a.swap(k, pivot_row);
            sign_is_negative = !sign_is_negative;
        }
        let pivot = a[k][k];
        det = mul_mod_u64(det, pivot, prime);
        let inverse = mod_inverse_prime(pivot, prime).expect("nonzero pivot modulo prime");
        let (pivot_rows, below) = a.split_at_mut(k + 1);
        let pivot_row = &pivot_rows[k];
        for row in &mut below[..size - k - 1] {
            if row[k] == 0 {
                continue;
            }
            let factor = mul_mod_u64(row[k], inverse, prime);
            for (target, &source) in row[(k + 1)..size].iter_mut().zip(&pivot_row[(k + 1)..size]) {
                let sub = mul_mod_u64(factor, source, prime);
                *target = sub_mod_u64(*target, sub, prime);
            }
        }
    }

    if sign_is_negative && det != 0 {
        prime - det
    } else {
        det
    }
}

pub(super) fn mul_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if modulus <= u32::MAX as u64 {
        return (a * b) % modulus;
    }
    ((a as u128 * b as u128) % modulus as u128) as u64
}

pub(super) fn sub_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    if a >= b {
        a - b
    } else {
        modulus - (b - a)
    }
}

pub(super) fn add_mod_u64(a: u64, b: u64, modulus: u64) -> u64 {
    let (sum, overflow) = a.overflowing_add(b);
    if overflow {
        ((a as u128 + b as u128) % modulus as u128) as u64
    } else if sum >= modulus {
        sum - modulus
    } else {
        sum
    }
}

fn pow_mod_u64(mut base: u64, mut exponent: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = mul_mod_u64(result, base, modulus);
        }
        base = mul_mod_u64(base, base, modulus);
        exponent >>= 1;
    }
    result
}

pub(super) fn mod_inverse_prime(value: u64, prime: u64) -> Option<u64> {
    (value != 0).then(|| pow_mod_u64(value, prime - 2, prime))
}

pub(super) fn signed_mod_u64(value: i64, prime: u64) -> u64 {
    (i128::from(value).rem_euclid(i128::from(prime))) as u64
}

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

pub(super) fn is_prime_u64(n: u64) -> bool {
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

    use crate::linalg::{
        determinants::bareiss_leading_principal_minors_bigint,
        test_support::{bi, bi_mat},
    };

    #[test]
    fn test_modular_leading_principal_minors_match_bareiss() {
        let m = bi_mat(&[&[4, 2, 1], &[2, 5, 3], &[1, 3, 6]]);
        let bareiss = bareiss_leading_principal_minors_bigint(&m).unwrap();
        let modular = modular_leading_principal_minors_bigint(&m).unwrap();
        assert_eq!(modular, bareiss);
    }

    #[test]
    fn test_modular_leading_principal_minors_handles_zero_bareiss_pivot() {
        let m = bi_mat(&[&[0, 1], &[1, 0]]);
        assert_eq!(bareiss_leading_principal_minors_bigint(&m), None);
        assert_eq!(
            modular_leading_principal_minors_bigint(&m),
            Some(vec![bi(0), bi(-1)])
        );
    }
}
