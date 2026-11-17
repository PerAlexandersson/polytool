//! Division-free integer matrix products with certified nonnegative reconstruction.
//!
//! The caller supplies an inclusive upper bound for each output column. Bounds
//! are mathematical premises, not estimates: finite modular agreement alone is
//! not an exact certificate. An unused prime also checks the reconstructed product.

use std::sync::{Mutex, OnceLock};

use num_bigint::BigInt;
use num_traits::{Signed, Zero};

use crate::chinese_remainder::{ChineseRemainderAccumulator, ChineseRemainderError};
use crate::modular_arithmetic::bigint_residue_u64;
use crate::primes::{DescendingPrimes, ValidatedPrime};

/// Failure of a bounded nonnegative matrix product.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NonnegativeMatrixProductError {
    /// Matrices are ragged, inner dimensions disagree, or column bounds have the wrong length.
    InvalidDimensions,
    /// An inclusive upper bound is negative.
    InvalidBound,
    /// The reconstructed product fails its bounds or congruences.
    Reconstruction(ChineseRemainderError),
    /// The reconstructed product disagrees with an independent unused prime.
    SparePrimeMismatch,
    /// All available primes have been exhausted before certification.
    InsufficientPrimes,
}

impl std::fmt::Display for NonnegativeMatrixProductError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDimensions => f.write_str("incompatible matrix or bound dimensions"),
            Self::InvalidBound => f.write_str("matrix product bounds must be nonnegative"),
            Self::Reconstruction(error) => error.fmt(f),
            Self::SparePrimeMismatch => f.write_str("matrix product spare-prime mismatch"),
            Self::InsufficientPrimes => f.write_str("insufficient matrix product primes"),
        }
    }
}
impl std::error::Error for NonnegativeMatrixProductError {}
impl From<ChineseRemainderError> for NonnegativeMatrixProductError {
    fn from(error: ChineseRemainderError) -> Self {
        Self::Reconstruction(error)
    }
}

/// Compute `left * right` by modular multiplication and nonnegative bounded CRT.
///
/// `column_bounds[j]` must be a proved inclusive bound for **every** exact output
/// entry in column `j`, and every output entry must be nonnegative. Input entries
/// may have either sign. Reconstruction in `[0, M)` occurs only when the product
/// of distinct primes exceeds all bounds; an unused prime then verifies it.
/// Incorrect bounds are rejected when detected but cannot be diagnosed in general.
/// No modular divisions or rank assumptions enter the matrix products.
///
/// An empty left matrix has zero columns; an empty right matrix uses the bound
/// vector to specify its width. The returned matrix has `left.len()` rows
/// and `column_bounds.len()` columns, including when the inner dimension is zero.
/// Validated small primes are cached across calls. Without a proved bound, use
/// an exact integer product instead.
pub fn multiply_integer_matrices_nonnegative_bounded(
    left: &[Vec<BigInt>],
    right: &[Vec<BigInt>],
    column_bounds: &[BigInt],
) -> Result<Vec<Vec<BigInt>>, NonnegativeMatrixProductError> {
    let inner = left.first().map_or(0, Vec::len);
    let columns = right.first().map_or(column_bounds.len(), Vec::len);
    if left.iter().any(|row| row.len() != inner)
        || right.iter().any(|row| row.len() != columns)
        || inner != right.len()
        || columns != column_bounds.len()
    {
        return Err(NonnegativeMatrixProductError::InvalidDimensions);
    }
    if column_bounds.iter().any(Signed::is_negative) {
        return Err(NonnegativeMatrixProductError::InvalidBound);
    }
    if left.is_empty() || columns == 0 || inner == 0 {
        return Ok(vec![vec![BigInt::zero(); columns]; left.len()]);
    }
    let mut accumulators: Vec<_> = (0..left.len())
        .map(|_| ChineseRemainderAccumulator::new(columns))
        .collect();
    let mut candidate: Option<Vec<Vec<BigInt>>> = None;
    for index in 0.. {
        let prime = cached_prime(index).ok_or(NonnegativeMatrixProductError::InsufficientPrimes)?;
        let p = prime.get();
        let a = reduce_matrix(left, p);
        let b = reduce_matrix(right, p);
        if let Some(result) = candidate {
            verify_product(&a, &b, &result, p)?;
            return Ok(result);
        }
        for (row, accumulator) in a.iter().zip(&mut accumulators) {
            accumulator.combine_prime(&product_row(row, &b, p), prime)?;
        }
        if accumulators[0].certifies_nonnegative_bounds(column_bounds)? {
            candidate = Some(
                accumulators
                    .iter()
                    .map(|accumulator| accumulator.reconstruct_nonnegative(column_bounds))
                    .collect::<Result<_, _>>()?,
            );
        }
    }
    unreachable!("unbounded prime index iterator")
}

fn reduce_matrix(matrix: &[Vec<BigInt>], prime: u64) -> Vec<Vec<u64>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|entry| bigint_residue_u64(entry, prime))
                .collect()
        })
        .collect()
}

fn verify_product(
    left: &[Vec<u64>],
    right: &[Vec<u64>],
    candidate: &[Vec<BigInt>],
    prime: u64,
) -> Result<(), NonnegativeMatrixProductError> {
    for (row, expected) in left.iter().zip(candidate) {
        if product_row(row, right, prime)
            .iter()
            .zip(expected)
            .any(|(&value, expected)| value != bigint_residue_u64(expected, prime))
        {
            return Err(NonnegativeMatrixProductError::SparePrimeMismatch);
        }
    }
    Ok(())
}

/// Delay reductions over a whole block, keeping every sum within u64.
fn product_row(left: &[u64], right: &[Vec<u64>], prime: u64) -> Vec<u64> {
    let mut output = vec![0u64; right[0].len()];
    let maximum = prime - 1;
    let block = usize::try_from((u64::MAX - maximum) / (maximum * maximum)).unwrap_or(usize::MAX);
    for (coefficients, rows) in left.chunks(block).zip(right.chunks(block)) {
        for (&coefficient, row) in coefficients.iter().zip(rows) {
            if coefficient == 0 {
                continue;
            }
            for (value, &entry) in output.iter_mut().zip(row) {
                // At block entry value < p. At most `block` products <= (p-1)^2
                // are added, so no addition or multiplication actually wraps.
                // Explicit wrapping operations allow vectorization even in checked builds.
                *value = value.wrapping_add(coefficient.wrapping_mul(entry));
            }
        }
        for value in &mut output {
            *value %= prime;
        }
    }
    output
}

fn cached_prime(index: usize) -> Option<ValidatedPrime> {
    type PrimeCache = Mutex<(Vec<ValidatedPrime>, DescendingPrimes)>;
    static PRIMES: OnceLock<PrimeCache> = OnceLock::new();
    let mut cache = PRIMES
        .get_or_init(|| Mutex::new((Vec::new(), DescendingPrimes::new((1 << 26) - 1))))
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    while cache.0.len() <= index {
        let prime = cache.1.next()?;
        cache.0.push(prime);
    }
    Some(cache.0[index])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn integer_matrix(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
        rows.iter()
            .map(|row| row.iter().map(|&v| BigInt::from(v)).collect())
            .collect()
    }

    #[test]
    fn bounded_products_match_exact_rectangular_signed_inputs() {
        for inner in 1..=9 {
            let left: Vec<Vec<BigInt>> = (0..4)
                .map(|i| (0..inner).map(|k| BigInt::from((i + k) % 7) - 3).collect())
                .collect();
            let right: Vec<_> = (0..inner)
                .map(|k| (0..4).map(|i| left[i][k].clone()).collect())
                .collect();
            // Shift by a positive Gram contribution so all products are nonnegative.
            let mut a = left.clone();
            for row in &mut a {
                row.push(100.into());
            }
            let mut b = right;
            b.push(vec![100.into(); 4]);
            let expected: Vec<Vec<BigInt>> = a
                .iter()
                .map(|row| {
                    (0..4)
                        .map(|j| row.iter().enumerate().map(|(k, v)| v * &b[k][j]).sum())
                        .collect()
                })
                .collect();
            assert_eq!(
                multiply_integer_matrices_nonnegative_bounded(
                    &a,
                    &b,
                    &[11000.into(), 11000.into(), 11000.into(), 11000.into()]
                )
                .unwrap(),
                expected
            );
        }
        let huge = BigInt::from(1) << 150;
        let a = vec![vec![huge]];
        let b = integer_matrix(&[&[7, 0]]);
        let bound: BigInt = &a[0][0] * 7;
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &b, &[bound.clone(), 0.into()])
                .unwrap(),
            vec![vec![bound, 0.into()]]
        );
    }

    #[test]
    fn malformed_bounds_shapes_and_spare_mismatch_are_errors() {
        use NonnegativeMatrixProductError::*;
        let a = integer_matrix(&[&[1]]);
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &a, &[]),
            Err(InvalidDimensions)
        );
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(
                &integer_matrix(&[&[1], &[]]),
                &a,
                &[1.into()]
            ),
            Err(InvalidDimensions)
        );
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &[], &[1.into()]),
            Err(InvalidDimensions)
        );
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &a, &[(-1).into()]),
            Err(InvalidBound)
        );
        assert!(matches!(
            multiply_integer_matrices_nonnegative_bounded(&a, &a, &[0.into()]),
            Err(Reconstruction(ChineseRemainderError::BoundViolation))
        ));
        let hidden = vec![vec![BigInt::from(cached_prime(0).unwrap().get()) + 1]];
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &hidden, &[1.into()]),
            Err(SparePrimeMismatch)
        );
        let boundary = vec![vec![BigInt::from(cached_prime(0).unwrap().get())]];
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&a, &boundary, &[boundary[0][0].clone()])
                .unwrap(),
            boundary
        );
        assert_eq!(
            verify_product(&[vec![1]], &[vec![1]], &[vec![2.into()]], 7),
            Err(SparePrimeMismatch)
        );
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&[vec![]], &[], &[0.into()]).unwrap(),
            vec![vec![BigInt::zero()]]
        );
        assert_eq!(
            multiply_integer_matrices_nonnegative_bounded(&[], &[], &[]).unwrap(),
            Vec::<Vec<BigInt>>::new()
        );
    }

    #[test]
    fn delayed_reduction_handles_block_boundaries_and_cached_primes() {
        let p = cached_prime(0).unwrap().get();
        let length = 10000;
        let a = vec![p - 1; length];
        let b = vec![vec![p - 1, 0]; length];
        assert_eq!(product_row(&a, &b, p), vec![length as u64 % p, 0]);
        assert_eq!(cached_prime(0), cached_prime(0));
        assert_ne!(cached_prime(0), cached_prime(1));
    }
}
