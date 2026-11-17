//! Differential checks against the former modular kernel implementation.
use super::*;

#[test]
fn overwide_rows_keep_the_former_prefix_kernel() {
    let rows = vec![vec![1.into(), 2.into(), 17.into(), 29.into()]];
    let expected = original::exact_kernel_vectors_by_group(&rows, &[0, 1], 2);
    assert!(expected[1].is_some());
    assert_eq!(exact_kernel_vectors_by_group(&rows, &[0, 1], 2), expected);
}

#[test]
#[should_panic(expected = "kernel rows have the requested width")]
fn short_rows_are_rejected_at_the_shared_shape_boundary() {
    exact_kernel_vectors_by_group(&[vec![1.into()]], &[0, 1], 2);
}
mod original {
    //! Exact kernel vectors of integer matrices through modular elimination:
    //! the reduced row-echelon form modulo several word-size primes, Chinese
    //! remaindering and rational reconstruction, and an exact check of every
    //! reconstructed vector against every row.  A vector is returned only after
    //! the exact check, so nothing that uses it depends on how it was found;
    //! a vector that is not found (an unlucky rank profile, too many primes) is
    //! reported as `None`, never guessed.

    #![allow(dead_code)]
    use crate::Rational;
    use num_bigint::BigInt;
    use num_integer::Integer;
    use num_traits::{One, Signed, Zero};
    use std::collections::BTreeMap;
    use std::sync::OnceLock;

    type Q = Rational;

    /// Largest number of primes combined before giving up.
    const MAX_PRIMES: usize = 16;

    /// The [`MAX_PRIMES`] largest primes below `2^31` in descending order,
    /// found once by trial division.
    pub(super) fn primes() -> impl Iterator<Item = u64> {
        static PRIMES: OnceLock<Vec<u64>> = OnceLock::new();
        PRIMES
            .get_or_init(|| {
                (1u64 << 30..1u64 << 31)
                    .rev()
                    .filter(|&candidate| {
                        candidate % 2 == 1
                            && (3..)
                                .step_by(2)
                                .take_while(|d| d * d <= candidate)
                                .all(|d| candidate % d != 0)
                    })
                    .take(MAX_PRIMES)
                    .collect()
            })
            .iter()
            .copied()
    }

    /// `value mod prime` in `[0, prime)` (the floor residue), from the 64-bit
    /// digits of the magnitude.
    fn residue(value: &BigInt, prime: u64) -> u64 {
        if value.is_zero() {
            return 0;
        }
        let modulus = u128::from(prime);
        let magnitude = value
            .magnitude()
            .iter_u64_digits()
            .rev()
            .fold(0u128, |acc, digit| {
                ((acc << 64) | u128::from(digit)) % modulus
            }) as u64;
        if value.is_negative() && magnitude != 0 {
            prime - magnitude
        } else {
            magnitude
        }
    }

    /// The integer row `L v` of a sparse rational equation `v` (`L` the least
    /// common multiple of the denominators), dense over `columns` columns.
    pub(crate) fn integer_equation_row(
        entries: &BTreeMap<usize, Q>,
        columns: usize,
    ) -> Vec<BigInt> {
        let denominator = entries.values().fold(BigInt::one(), |acc, value| {
            if value.denom().is_one() {
                acc
            } else {
                acc.lcm(value.denom())
            }
        });
        let mut row = vec![BigInt::zero(); columns];
        for (column, value) in entries {
            row[*column] = if value.denom().is_one() {
                if denominator.is_one() {
                    value.numer().clone()
                } else {
                    value.numer() * &denominator
                }
            } else {
                value.numer() * (&denominator / value.denom())
            };
        }
        row
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
            // Only the nonzero entries of the pivot row change other rows.
            let support: Vec<usize> = (column..pivot_row.len())
                .filter(|&c| pivot_row[c] != 0)
                .collect();
            for (r, row) in matrix.iter_mut().enumerate() {
                if r == rank || row[column] == 0 {
                    continue;
                }
                let factor = prime - row[column];
                for &c in &support {
                    row[c] = (row[c] + factor * pivot_row[c]) % prime;
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
    pub fn exact_kernel_vectors_by_group(
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
}

#[test]
fn cached_validated_primes_keep_every_former_value_and_order() {
    let actual: Vec<_> = validated_primes().map(ValidatedPrime::get).collect();
    assert_eq!(actual, original::primes().collect::<Vec<_>>());
    assert_eq!(actual.len(), MAX_PRIMES);
    assert!(actual
        .iter()
        .all(|p| ((1u64 << 30)..(1u64 << 31)).contains(p)));
    assert!(actual.windows(2).all(|pair| pair[0] > pair[1]));
}
#[test]
fn grouped_kernel_matches_former_random_and_unlucky_prime_cases() {
    let mut seed = 0x7429_813b_63d0_1a95u64;
    for columns in 0..=8 {
        for height in 0..=8 {
            let rows: Vec<Vec<BigInt>> = (0..height)
                .map(|_| {
                    (0..columns)
                        .map(|_| {
                            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                            BigInt::from((seed >> 32) as i64 % 11 - 5)
                        })
                        .collect()
                })
                .collect();
            let group: Vec<usize> = (0..columns).map(|c| c % 2).collect();
            assert_eq!(
                exact_kernel_vectors_by_group(&rows, &group, 2),
                original::exact_kernel_vectors_by_group(&rows, &group, 2)
            );
        }
    }
    let p = primes().next().unwrap();
    for rows in [
        vec![vec![BigInt::from(p), 1.into()]],
        vec![vec![BigInt::from(p), 0.into()], vec![0.into(), 1.into()]],
        vec![vec![0.into(), 0.into()]],
        vec![vec![1.into(), BigInt::one() << 200usize]],
    ] {
        assert_eq!(
            exact_kernel_vectors_by_group(&rows, &[0, 1], 2),
            original::exact_kernel_vectors_by_group(&rows, &[0, 1], 2)
        );
    }
}
