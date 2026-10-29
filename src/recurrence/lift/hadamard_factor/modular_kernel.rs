//! Exact kernel vectors of integer matrices through modular elimination:
//! the reduced row-echelon form modulo several word-size primes, Chinese
//! remaindering and rational reconstruction, and an exact check of every
//! reconstructed vector against every row.  Used to guess row annihilators
//! (`annihilator_system`); a vector is returned only after the exact check,
//! and nothing in a certificate depends on how it was found.

use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};

type Q = BigRational;

/// Largest number of primes combined before giving up.
const MAX_PRIMES: usize = 16;

/// The primes below `2^31` in descending order, found by trial division.
fn primes() -> impl Iterator<Item = u64> {
    (1u64 << 30..1u64 << 31).rev().filter(|&candidate| {
        candidate % 2 == 1
            && (3..)
                .step_by(2)
                .take_while(|d| d * d <= candidate)
                .all(|d| candidate % d != 0)
    })
}

fn residue(value: &BigInt, prime: u64) -> u64 {
    value
        .mod_floor(&BigInt::from(prime))
        .to_u64()
        .unwrap_or_default()
}

fn inverse(value: u64, prime: u64) -> u64 {
    let (mut base, mut exponent, mut acc) = (value % prime, prime - 2, 1u64);
    while exponent > 0 {
        if exponent & 1 == 1 {
            acc = acc * base % prime;
        }
        base = base * base % prime;
        exponent >>= 1;
    }
    acc
}

/// The reduced row-echelon form modulo `prime`: pivot columns and the
/// reduced pivot rows.
fn reduced_echelon(
    rows: &[Vec<BigInt>],
    columns: usize,
    prime: u64,
) -> (Vec<usize>, Vec<Vec<u64>>) {
    let mut matrix: Vec<Vec<u64>> = rows
        .iter()
        .map(|row| row.iter().map(|value| residue(value, prime)).collect())
        .collect();
    let mut pivots = Vec::new();
    for column in 0..columns {
        let rank = pivots.len();
        let Some(found) = (rank..matrix.len()).find(|&r| matrix[r][column] != 0) else {
            continue;
        };
        matrix.swap(rank, found);
        let scale = inverse(matrix[rank][column], prime);
        for value in &mut matrix[rank][column..] {
            *value = *value * scale % prime;
        }
        let pivot_row = matrix[rank].clone();
        for (r, row) in matrix.iter_mut().enumerate() {
            if r == rank || row[column] == 0 {
                continue;
            }
            let factor = row[column];
            for (value, pivot) in row[column..].iter_mut().zip(&pivot_row[column..]) {
                *value = (*value + (prime - factor) * pivot) % prime;
            }
        }
        pivots.push(column);
    }
    matrix.truncate(pivots.len());
    (pivots, matrix)
}

/// `p / q` with `|p|, q <= sqrt(modulus / 2)` congruent to `value`, if any.
fn reconstruct(value: &BigInt, modulus: &BigInt) -> Option<Q> {
    let bound = (modulus / BigInt::from(2)).sqrt();
    let (mut r0, mut r1) = (modulus.clone(), value.mod_floor(modulus));
    let (mut t0, mut t1) = (BigInt::zero(), BigInt::one());
    while r1 > bound {
        let quotient = &r0 / &r1;
        let r2 = &r0 - &quotient * &r1;
        let t2 = &t0 - &quotient * &t1;
        (r0, r1, t0, t1) = (r1, r2, t1, t2);
    }
    if t1.is_zero() || t1.abs() > bound {
        return None;
    }
    Some(Q::new(r1, t1))
}

/// Every row of the integer matrix annihilates `vector` exactly.
fn kills_every_row(rows: &[Vec<BigInt>], vector: &[Q]) -> bool {
    let denominator = vector
        .iter()
        .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
    let integral: Vec<BigInt> = vector
        .iter()
        .map(|value| value.numer() * (&denominator / value.denom()))
        .collect();
    rows.iter().all(|row| {
        row.iter()
            .zip(&integral)
            .filter(|(a, b)| !a.is_zero() && !b.is_zero())
            .fold(BigInt::zero(), |acc, (a, b)| acc + a * b)
            .is_zero()
    })
}

/// For every column group (`group[column]`, groups `0..groups`), the kernel
/// vector of the integer matrix `rows` whose largest nonzero column is the
/// smallest non-pivot column of that group (modulo the first prime), with
/// one there and zero in every other non-pivot column; `None` for a group
/// without one or when no exact vector was reconstructed.
pub(super) fn kernel_vectors_by_group(
    rows: &[Vec<BigInt>],
    group: &[usize],
    groups: usize,
) -> Vec<Option<Vec<Q>>> {
    let columns = group.len();
    let mut out = vec![None; groups];
    let mut pivots: Option<Vec<usize>> = None;
    let mut targets: Vec<Option<usize>> = vec![None; groups];
    let mut modulus = BigInt::one();
    let mut combined = vec![vec![BigInt::zero(); columns]; groups];
    for prime in primes().take(MAX_PRIMES) {
        let (prime_pivots, reduced) = reduced_echelon(rows, columns, prime);
        match &pivots {
            None => {
                for column in (0..columns).rev() {
                    if !prime_pivots.contains(&column) {
                        targets[group[column]] = Some(column);
                    }
                }
                pivots = Some(prime_pivots);
            }
            // A prime with another rank profile is skipped.
            Some(expected) if *expected != prime_pivots => continue,
            Some(_) => {}
        }
        let prime_big = BigInt::from(prime);
        let lift = BigInt::from(inverse(residue(&modulus, prime), prime));
        let mut open = false;
        for (slot, target) in targets.iter().enumerate() {
            let Some(free) = *target else {
                continue;
            };
            if out[slot].is_some() {
                continue;
            }
            // v[free] = 1, v[pivot j] = -reduced[j][free], zero elsewhere.
            let mut residues = vec![0u64; columns];
            residues[free] = 1;
            for (row, &column) in reduced.iter().zip(pivots.iter().flatten()) {
                residues[column] = (prime - row[free]) % prime;
            }
            // value' = value + modulus ((residue - value) / modulus mod prime).
            for (value, residue) in combined[slot].iter_mut().zip(&residues) {
                let difference = (BigInt::from(*residue) - &*value).mod_floor(&prime_big);
                *value += &modulus * (difference * &lift % &prime_big);
            }
            let next_modulus = &modulus * &prime_big;
            let vector: Option<Vec<Q>> = combined[slot]
                .iter()
                .map(|value| reconstruct(value, &next_modulus))
                .collect();
            match vector {
                Some(vector) if kills_every_row(rows, &vector) => out[slot] = Some(vector),
                _ => open = true,
            }
        }
        modulus *= &prime_big;
        if !open {
            break;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn int_rows(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|v| BigInt::from(*v)).collect())
            .collect()
    }

    #[test]
    fn kernel_vectors_are_exact_and_echelon() {
        // Kernel of [[1, 2, 3, 4], [2, 4, 7, 9]]: columns grouped {0, 1} and {2, 3}.
        // Free columns modulo p: 1 and 3.
        let rows = int_rows(&[&[1, 2, 3, 4], &[2, 4, 7, 9]]);
        let vectors = kernel_vectors_by_group(&rows, &[0, 0, 1, 1], 2);
        let first = vectors[0].as_ref().expect("a vector led by column 1");
        assert_eq!(first[1], Q::one());
        assert!(first[2].is_zero() && first[3].is_zero());
        let second = vectors[1].as_ref().expect("a vector led by column 3");
        assert_eq!(second[3], Q::one());
        assert!(second[1].is_zero());
        assert!(kills_every_row(&rows, first) && kills_every_row(&rows, second));
        // A nonsingular matrix has no kernel vectors.
        let rows = int_rows(&[&[1, 0], &[0, 3]]);
        assert_eq!(kernel_vectors_by_group(&rows, &[0, 0], 1), vec![None]);
        // Large rational entries are reconstructed.
        let big = BigInt::from(10).pow(15u32) + BigInt::from(7);
        let rows = vec![vec![big.clone(), BigInt::from(-3)]];
        let vector = kernel_vectors_by_group(&rows, &[0, 0], 1)[0]
            .clone()
            .unwrap();
        assert_eq!(vector[0], Q::new(BigInt::from(3), big));
    }
}
