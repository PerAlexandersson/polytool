//! Dense finite-field row reduction shared by modular kernel consumers.

use crate::modular_arithmetic::{add_modulo_u64, multiply_modulo_u64};
use crate::primes::ValidatedPrime;

/// A rectangular matrix has a row of the wrong length.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ModularMatrixDimensionsError;

impl std::fmt::Display for ModularMatrixDimensionsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("modular matrix rows have inconsistent dimensions")
    }
}

impl std::error::Error for ModularMatrixDimensionsError {}

/// Reduced pivot rows and their zero-based pivot columns, in increasing order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ModularRowEchelon {
    /// The nonzero rows in reduced row-echelon form, with canonical residues.
    pub rows: Vec<Vec<u64>>,
    /// The column of the leading one in each returned row.
    pub pivot_columns: Vec<usize>,
}

/// Reduce a rectangular matrix over a validated prime field.
///
/// Entries are reduced first. `columns` also specifies the width of an empty
/// matrix. The first available pivot row is selected in input order.
pub fn reduced_row_echelon_modulo_prime(
    mut matrix: Vec<Vec<u64>>,
    columns: usize,
    prime: ValidatedPrime,
) -> Result<ModularRowEchelon, ModularMatrixDimensionsError> {
    if matrix.iter().any(|row| row.len() != columns) {
        return Err(ModularMatrixDimensionsError);
    }
    let p = prime.get();
    for value in matrix.iter_mut().flatten() {
        *value %= p;
    }
    let mut pivots = Vec::new();
    for column in 0..columns {
        let rank = pivots.len();
        let Some(found) = (rank..matrix.len()).find(|&r| matrix[r][column] != 0) else {
            continue;
        };
        matrix.swap(rank, found);
        let scale = prime
            .inverse(matrix[rank][column])
            .expect("nonzero field pivot");
        for value in &mut matrix[rank][column..] {
            *value = multiply_modulo_u64(*value, scale, p);
        }
        let pivot_row = matrix[rank].clone();
        let support: Vec<usize> = (column..columns).filter(|&c| pivot_row[c] != 0).collect();
        for (r, row) in matrix.iter_mut().enumerate() {
            if r == rank || row[column] == 0 {
                continue;
            }
            let factor = p - row[column];
            for &c in &support {
                row[c] = add_modulo_u64(row[c], multiply_modulo_u64(factor, pivot_row[c], p), p);
            }
        }
        pivots.push(column);
    }
    matrix.truncate(pivots.len());
    Ok(ModularRowEchelon {
        rows: matrix,
        pivot_columns: pivots,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn shape_validation_and_unlucky_rank() {
        let p = ValidatedPrime::new(7).unwrap();
        assert_eq!(
            reduced_row_echelon_modulo_prime(vec![vec![1]], 2, p),
            Err(ModularMatrixDimensionsError)
        );
        assert!(reduced_row_echelon_modulo_prime(vec![], 3, p)
            .unwrap()
            .rows
            .is_empty());
        let rows = vec![vec![7, 1], vec![0, 1]];
        assert_eq!(
            reduced_row_echelon_modulo_prime(rows.clone(), 2, p)
                .unwrap()
                .pivot_columns,
            vec![1]
        );
        assert_eq!(
            reduced_row_echelon_modulo_prime(rows, 2, ValidatedPrime::new(11).unwrap())
                .unwrap()
                .pivot_columns,
            vec![0, 1]
        );
    }
}
