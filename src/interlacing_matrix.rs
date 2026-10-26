//! Athanasiadis--Wagner interlacing matrices.
//!
//! For a `p x q` matrix `A` of formal power series, Athanasiadis and Wagner
//! define an infinite block Toeplitz matrix `Lace(A)`.  If
//! `A_ij(x) = sum_n a_ij(n) x^n`, then the entry with row `u = p u' + i`
//! and column `v = q v' + j` is
//!
//! ```text
//! [Lace(A)]_{u,v} = a_ij(v' - u').
//! ```
//!
//! This module implements finite polynomial truncations of this infinite
//! matrix.  A successful finite total-nonnegativity check is useful
//! computational evidence; it is not, by itself, a proof that the infinite
//! matrix is totally nonnegative unless a separate finite criterion applies.
//!
//! Entry points: [`lace_matrix`] (generic) with the `_i64`/`_bigint` wrappers,
//! [`lace_matrix_sequence`] for a column vector of polynomials, and the TNN
//! checks [`check_lace_tnn_neville_bigint`] / [`check_lace_tnn_neville_i64`]
//! (full total nonnegativity of the finite matrix by Neville elimination) and
//! [`check_lace_total_nonnegative_i64`] (brute force over minors up to a given
//! size).
//!
//! Conventions: a polynomial is a `Vec<C>` of ascending coefficients (an empty
//! vector is the zero polynomial); `block_rows` and `block_cols` count how
//! many `p`-row and `q`-column blocks of the infinite matrix are kept.  Matrix
//! rows and columns are zero-based.  The checks return `Result<(), String>`
//! with a message describing the failing minor or malformed input.  The TNN
//! machinery itself is in [`crate::linalg`].

use crate::linalg::{check_tnn_neville, check_tnn_neville_bigint, check_total_positivity};
use num_bigint::BigInt;
use num_traits::Zero;
use std::error::Error;
use std::fmt;

/// Error returned when an interlacing matrix input is not rectangular.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InterlacingMatrixError {
    /// The matrix has no rows.
    EmptyMatrix,
    /// Row `row` has no columns.
    EmptyRow {
        /// Zero-based index of the empty row.
        row: usize,
    },
    /// Row `row` has a different number of columns than the first row.
    RaggedRows {
        /// Zero-based index of the offending row.
        row: usize,
        /// Number of columns of the first row.
        expected: usize,
        /// Number of columns in row `row`.
        found: usize,
    },
}

impl fmt::Display for InterlacingMatrixError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMatrix => write!(f, "polynomial matrix must have at least one row"),
            Self::EmptyRow { row } => {
                write!(
                    f,
                    "polynomial matrix row {row} must have at least one column"
                )
            }
            Self::RaggedRows {
                row,
                expected,
                found,
            } => write!(
                f,
                "polynomial matrix row {row} has {found} columns, expected {expected}"
            ),
        }
    }
}

impl Error for InterlacingMatrixError {}

/// Build a finite truncation of the Athanasiadis--Wagner `Lace(A)` matrix.
///
/// The input is a rectangular matrix of polynomials in ascending coefficient
/// order.  If `matrix` is `p x q`, the output has
/// `p * block_rows` rows and `q * block_cols` columns.  Its entry in row
/// `p * rb + i` and column `q * cb + j` is the coefficient of
/// `x^(cb-rb)` in `matrix[i][j]`, with negative and out-of-range coefficient
/// indices interpreted as zero.
///
/// Empty coefficient vectors are allowed and represent the zero polynomial.
pub fn lace_matrix<C>(
    matrix: &[Vec<Vec<C>>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<C>>, InterlacingMatrixError>
where
    C: Clone + Zero,
{
    let (p, q) = validate_rectangular_matrix(matrix)?;
    let mut result = vec![vec![C::zero(); q * block_cols]; p * block_rows];

    for rb in 0..block_rows {
        for cb in rb..block_cols {
            let coeff_index = cb - rb;
            for i in 0..p {
                for j in 0..q {
                    if let Some(value) = matrix[i][j].get(coeff_index) {
                        result[p * rb + i][q * cb + j] = value.clone();
                    }
                }
            }
        }
    }

    Ok(result)
}

/// Build a finite `Lace(A)` truncation for integer coefficient polynomials.
pub fn lace_matrix_i64(
    matrix: &[Vec<Vec<i64>>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<i64>>, InterlacingMatrixError> {
    lace_matrix(matrix, block_rows, block_cols)
}

/// Build a finite `Lace(A)` truncation for `BigInt` coefficient polynomials.
pub fn lace_matrix_bigint(
    matrix: &[Vec<Vec<BigInt>>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<BigInt>>, InterlacingMatrixError> {
    lace_matrix(matrix, block_rows, block_cols)
}

/// Build a finite `Lace(A)` truncation for a column vector of polynomials.
///
/// Athanasiadis--Wagner call a column vector fully interlacing when the
/// corresponding infinite `Lace(A)` is totally nonnegative.  This function
/// returns the finite truncation used for experiments with that condition.
pub fn lace_matrix_sequence<C>(
    polynomials: &[Vec<C>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<C>>, InterlacingMatrixError>
where
    C: Clone + Zero,
{
    if polynomials.is_empty() {
        return Err(InterlacingMatrixError::EmptyMatrix);
    }
    let matrix: Vec<Vec<Vec<C>>> = polynomials.iter().map(|p| vec![p.clone()]).collect();
    lace_matrix(&matrix, block_rows, block_cols)
}

/// Build a finite `Lace(A)` truncation for an integer polynomial sequence.
pub fn lace_matrix_sequence_i64(
    polynomials: &[Vec<i64>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<i64>>, InterlacingMatrixError> {
    lace_matrix_sequence(polynomials, block_rows, block_cols)
}

/// Build a finite `Lace(A)` truncation for a `BigInt` polynomial sequence.
pub fn lace_matrix_sequence_bigint(
    polynomials: &[Vec<BigInt>],
    block_rows: usize,
    block_cols: usize,
) -> Result<Vec<Vec<BigInt>>, InterlacingMatrixError> {
    lace_matrix_sequence(polynomials, block_rows, block_cols)
}

/// Brute-force finite total-nonnegativity check for a `Lace(A)` truncation.
///
/// This checks all minors up to `max_minor_size` in the finite truncation.
/// The cost grows combinatorially with the matrix size and `max_minor_size`;
/// use [`check_lace_tnn_neville_i64`] for a full test of the finite matrix.
/// Returns `Err` with a message for malformed input or a negative minor.
pub fn check_lace_total_nonnegative_i64(
    matrix: &[Vec<Vec<i64>>],
    block_rows: usize,
    block_cols: usize,
    max_minor_size: usize,
) -> Result<(), String> {
    let lace = lace_matrix_i64(matrix, block_rows, block_cols).map_err(|err| err.to_string())?;
    check_total_positivity(&lace, max_minor_size, false)
        .map_err(|err| format!("finite Lace truncation is not TNN: {err}"))
}

/// Brute-force finite total-nonnegativity check for a polynomial sequence.
pub fn check_lace_sequence_total_nonnegative_i64(
    polynomials: &[Vec<i64>],
    block_rows: usize,
    block_cols: usize,
    max_minor_size: usize,
) -> Result<(), String> {
    let lace = lace_matrix_sequence_i64(polynomials, block_rows, block_cols)
        .map_err(|err| err.to_string())?;
    check_total_positivity(&lace, max_minor_size, false)
        .map_err(|err| format!("finite Lace truncation is not TNN: {err}"))
}

/// Return whether a finite `Lace(A)` truncation is TNN up to `max_minor_size`.
pub fn is_lace_totally_nonnegative_i64(
    matrix: &[Vec<Vec<i64>>],
    block_rows: usize,
    block_cols: usize,
    max_minor_size: usize,
) -> bool {
    check_lace_total_nonnegative_i64(matrix, block_rows, block_cols, max_minor_size).is_ok()
}

/// Return whether a sequence's finite `Lace(A)` truncation is TNN up to
/// `max_minor_size`.
pub fn is_lace_sequence_totally_nonnegative_i64(
    polynomials: &[Vec<i64>],
    block_rows: usize,
    block_cols: usize,
    max_minor_size: usize,
) -> bool {
    check_lace_sequence_total_nonnegative_i64(polynomials, block_rows, block_cols, max_minor_size)
        .is_ok()
}

/// Neville-elimination TNN check for a finite integer `Lace(A)` truncation.
pub fn check_lace_tnn_neville_i64(
    matrix: &[Vec<Vec<i64>>],
    block_rows: usize,
    block_cols: usize,
) -> Result<(), String> {
    let lace = lace_matrix_i64(matrix, block_rows, block_cols).map_err(|err| err.to_string())?;
    check_tnn_neville(&lace).map_err(|err| format!("finite Lace truncation is not TNN: {err}"))
}

/// Neville-elimination TNN check for a finite `BigInt` `Lace(A)` truncation.
pub fn check_lace_tnn_neville_bigint(
    matrix: &[Vec<Vec<BigInt>>],
    block_rows: usize,
    block_cols: usize,
) -> Result<(), String> {
    let lace = lace_matrix_bigint(matrix, block_rows, block_cols).map_err(|err| err.to_string())?;
    check_tnn_neville_bigint(&lace)
        .map_err(|err| format!("finite Lace truncation is not TNN: {err}"))
}

fn validate_rectangular_matrix<C>(
    matrix: &[Vec<Vec<C>>],
) -> Result<(usize, usize), InterlacingMatrixError> {
    if matrix.is_empty() {
        return Err(InterlacingMatrixError::EmptyMatrix);
    }
    let q = matrix[0].len();
    if q == 0 {
        return Err(InterlacingMatrixError::EmptyRow { row: 0 });
    }
    for (row, entries) in matrix.iter().enumerate().skip(1) {
        if entries.is_empty() {
            return Err(InterlacingMatrixError::EmptyRow { row });
        }
        if entries.len() != q {
            return Err(InterlacingMatrixError::RaggedRows {
                row,
                expected: q,
                found: entries.len(),
            });
        }
    }
    Ok((matrix.len(), q))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::check_total_positivity;

    #[test]
    fn sequence_truncation_uses_athanasiadis_wagner_indexing() {
        let polynomials = vec![vec![1, 1], vec![1, 2, 1]];
        let lace = lace_matrix_sequence_i64(&polynomials, 3, 4).unwrap();

        assert_eq!(
            lace,
            vec![
                vec![1, 1, 0, 0],
                vec![1, 2, 1, 0],
                vec![0, 1, 1, 0],
                vec![0, 1, 2, 1],
                vec![0, 0, 1, 1],
                vec![0, 0, 1, 2],
            ]
        );
    }

    #[test]
    fn polynomial_matrix_truncation_uses_block_toeplitz_form() {
        let matrix = vec![vec![vec![1], vec![2, 3]], vec![vec![4, 5], vec![6, 7, 8]]];
        let lace = lace_matrix_i64(&matrix, 2, 3).unwrap();

        assert_eq!(
            lace,
            vec![
                vec![1, 2, 0, 3, 0, 0],
                vec![4, 6, 5, 7, 0, 8],
                vec![0, 0, 1, 2, 0, 3],
                vec![0, 0, 4, 6, 5, 7],
            ]
        );
    }

    #[test]
    fn finite_lace_detects_athanasiadis_wagner_pairwise_not_full_example() {
        // Athanasiadis--Wagner, Example 3.4:
        // P=t+x, Q=(b+x)(d+x), R=(a+x)(c+x) with
        // a <= b <= t <= c <= d can be pairwise interlacing without being
        // fully interlacing.  Here a=1, b=2, t=2, c=3, d=4 gives determinant -1.
        let polynomials = vec![vec![2, 1], vec![8, 6, 1], vec![3, 4, 1]];
        let lace = lace_matrix_sequence_i64(&polynomials, 1, 3).unwrap();

        assert_eq!(lace, vec![vec![2, 1, 0], vec![8, 6, 1], vec![3, 4, 1]]);
        assert!(check_total_positivity(&lace, 3, false).is_err());
        let err = check_lace_sequence_total_nonnegative_i64(&polynomials, 1, 3, 3).unwrap_err();
        assert!(err.contains("det = -1"));
    }

    #[test]
    fn finite_lace_tnn_check_accepts_pf_toeplitz_example() {
        let matrix = vec![vec![vec![1, 1]]];
        let lace = lace_matrix_i64(&matrix, 4, 5).unwrap();

        assert_eq!(
            lace,
            vec![
                vec![1, 1, 0, 0, 0],
                vec![0, 1, 1, 0, 0],
                vec![0, 0, 1, 1, 0],
                vec![0, 0, 0, 1, 1],
            ]
        );
        assert!(check_lace_total_nonnegative_i64(&matrix, 4, 5, 4).is_ok());
        assert!(check_lace_tnn_neville_i64(&matrix, 4, 5).is_ok());
    }

    #[test]
    fn validates_rectangular_input() {
        let empty: Vec<Vec<Vec<i64>>> = Vec::new();
        assert_eq!(
            lace_matrix_i64(&empty, 1, 1).unwrap_err(),
            InterlacingMatrixError::EmptyMatrix
        );

        let ragged = vec![vec![vec![1]], vec![vec![1], vec![2]]];
        assert_eq!(
            lace_matrix_i64(&ragged, 1, 1).unwrap_err(),
            InterlacingMatrixError::RaggedRows {
                row: 1,
                expected: 1,
                found: 2,
            }
        );
    }
}
