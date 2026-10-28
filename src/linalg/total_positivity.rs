//! Total non-negativity: brute-force minor enumeration and Neville elimination.
//!
//! [`check_total_positivity`] enumerates minors up
//! to a given size (exponential); [`check_totally_nonnegative_neville`]
//! and its `BigInt` variant decide total non-negativity of all minors in
//! `O(n^2 m)` exact rational operations.  All checks report failures as a
//! [`TotalPositivityError`]; [`is_totally_nonnegative_neville`] and
//! [`is_totally_nonnegative_by_minors`] are the `bool` forms.

use super::determinants::determinant;
use super::Q;
use num_bigint::BigInt;
use num_traits::Zero;

/// Exact failure reported by the total-positivity checks: brute-force minor
/// enumeration ([`check_total_positivity`]) and Neville elimination
/// ([`check_totally_nonnegative_neville`], [`check_totally_nonnegative_neville_bigint`]).
///
/// The first three variants come from minor enumeration, the remaining ones
/// from Neville elimination (a failed elimination step is reported instead of
/// an explicit minor).
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TotalPositivityError {
    /// The rows do not all have the same length.
    RaggedMatrix {
        /// Zero-based index of the first row whose length differs from that of row 0.
        row: usize,
        /// Length of row 0.
        expected_columns: usize,
        /// Length of the offending row.
        found_columns: usize,
    },
    /// Some minor is negative.
    NegativeMinor {
        /// Order of the minor.
        size: usize,
        /// Zero-based indices of the rows of the minor, in increasing order.
        rows: Vec<usize>,
        /// Zero-based indices of the columns of the minor, in increasing order.
        columns: Vec<usize>,
        /// Value of the minor.
        determinant: BigInt,
    },
    /// In strict mode, a minor that is zero although its submatrix is not entirely zero.
    ZeroMinor {
        /// Order of the minor.
        size: usize,
        /// Zero-based indices of the rows of the minor, in increasing order.
        rows: Vec<usize>,
        /// Zero-based indices of the columns of the minor, in increasing order.
        columns: Vec<usize>,
    },
    /// Neville elimination: the rows do not all have the same length.
    InconsistentRowLengths,
    /// Neville elimination: an entry of the input matrix is negative.
    NegativeEntry {
        /// Zero-based row index of the entry.
        row: usize,
        /// Zero-based column index of the entry.
        column: usize,
        /// The entry.  The `i64` entry point ([`check_totally_nonnegative_neville`]) reports it
        /// and the `BigInt` entry point leaves it `None`, matching the
        /// original error messages.
        value: Option<BigInt>,
    },
    /// Neville elimination: a positive entry lies directly below a zero pivot
    /// whose row has a later nonzero entry, which forces a negative `2 x 2`
    /// minor.
    PositiveEntryBelowZeroPivot {
        /// Zero-based row index of the positive entry.
        row: usize,
        /// Zero-based column index of the positive entry (the pivot column).
        column: usize,
    },
    /// Neville elimination: subtracting a multiple of row `row - 1` from row
    /// `row` made an entry negative.
    NegativeEliminationEntry {
        /// Zero-based index of the lower of the two adjacent rows.
        row: usize,
        /// Zero-based column index of the entry that became negative.
        column: usize,
        /// Zero-based column of the pivot used in the elimination step.
        pivot_column: usize,
    },
}

impl std::fmt::Display for TotalPositivityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::RaggedMatrix {
                row,
                expected_columns,
                found_columns,
            } => write!(
                f,
                "matrix row {row} has {found_columns} columns; expected {expected_columns}"
            ),
            Self::NegativeMinor {
                size,
                rows,
                columns,
                determinant,
            } => write!(
                f,
                "{size}x{size} minor rows {rows:?} cols {columns:?} has det = {determinant} < 0"
            ),
            Self::ZeroMinor {
                size,
                rows,
                columns,
            } => write!(
                f,
                "{size}x{size} minor rows {rows:?} cols {columns:?} has det = 0 (strict mode)"
            ),
            Self::InconsistentRowLengths => write!(f, "matrix rows have inconsistent lengths"),
            Self::NegativeEntry {
                row,
                column,
                value: Some(value),
            } => write!(f, "Entry [{row},{column}] = {value} < 0"),
            Self::NegativeEntry {
                row,
                column,
                value: None,
            } => write!(f, "Entry [{row},{column}] < 0"),
            Self::PositiveEntryBelowZeroPivot { row, column } => write!(
                f,
                "Neville elimination: positive entry [{row},{column}] below zero pivot \
                 with nonzero pivot-row suffix"
            ),
            Self::NegativeEliminationEntry {
                row,
                column,
                pivot_column,
            } => write!(
                f,
                "Neville elimination: entry [{row},{column}] became negative \
                 (multiplier at column {pivot_column}, rows [{},{row}])",
                row - 1
            ),
        }
    }
}

impl std::error::Error for TotalPositivityError {}

/// Check if an integer matrix is totally positive (all minors strictly positive,
/// except those that are trivially zero due to zero rows/columns).
///
/// More precisely, checks that all minors are **non-negative** (totally non-negative, TNN).
/// Set `strict` to require all non-trivially-zero minors to be strictly positive.
///
/// `max_minor_size` limits the largest minor checked. For a full check, pass
/// `min(nrows, ncols)`. Returns `Ok(())` if all minors pass, or a structured
/// [`TotalPositivityError`] describing ragged input or the first violation.
///
/// # Example
/// ```
/// use polytool::linalg::check_total_positivity;
/// // Pascal's triangle is totally positive
/// let mat = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_total_positivity(&mat, 3, false).is_ok());
/// ```
pub fn check_total_positivity(
    mat: &[Vec<i64>],
    max_minor_size: usize,
    strict: bool,
) -> Result<(), TotalPositivityError> {
    check_total_positivity_entries(mat, max_minor_size, strict)
}

/// `BigInt` version of [`check_total_positivity`], for matrices whose entries
/// may not fit in `i64`.
///
/// The contract, minor order and reported errors are those of the `i64`
/// function.
///
/// # Example
/// ```
/// use num_bigint::BigInt;
/// use polytool::linalg::check_total_positivity_bigint;
/// // Rows of (1 + 10^20 x)^n are far beyond i64.
/// let big = BigInt::from(10).pow(20);
/// let mat = vec![
///     vec![BigInt::from(1), BigInt::from(0)],
///     vec![BigInt::from(1), big.clone()],
///     vec![BigInt::from(1), big * 2],
/// ];
/// assert!(check_total_positivity_bigint(&mat, 2, false).is_ok());
/// ```
pub fn check_total_positivity_bigint(
    mat: &[Vec<BigInt>],
    max_minor_size: usize,
    strict: bool,
) -> Result<(), TotalPositivityError> {
    check_total_positivity_entries(mat, max_minor_size, strict)
}

// Shared kernel of the `i64` and `BigInt` entry points; the minors are exact
// `BigInt` determinants in both cases.
fn check_total_positivity_entries<T>(
    mat: &[Vec<T>],
    max_minor_size: usize,
    strict: bool,
) -> Result<(), TotalPositivityError>
where
    T: Clone + Into<BigInt>,
{
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if let Some((row, found_columns)) = mat
        .iter()
        .enumerate()
        .find_map(|(row, values)| (values.len() != ncols).then_some((row, values.len())))
    {
        return Err(TotalPositivityError::RaggedMatrix {
            row,
            expected_columns: ncols,
            found_columns,
        });
    }
    let max_k = max_minor_size.min(nrows).min(ncols);

    for k in 1..=max_k {
        for rows in CombinationIterator::new(nrows, k) {
            for cols in CombinationIterator::new(ncols, k) {
                let sub = extract_submatrix(mat, &rows, &cols);
                let det = determinant(&sub);
                if det < BigInt::zero() {
                    return Err(TotalPositivityError::NegativeMinor {
                        size: k,
                        rows,
                        columns: cols,
                        determinant: det,
                    });
                }
                if strict && det.is_zero() {
                    // Check if this is a non-trivial zero
                    let all_zero = rows
                        .iter()
                        .all(|&r| cols.iter().all(|&c| mat[r][c].clone().into().is_zero()));
                    if !all_zero {
                        return Err(TotalPositivityError::ZeroMinor {
                            size: k,
                            rows,
                            columns: cols,
                        });
                    }
                }
            }
        }
    }
    Ok(())
}

/// Check total non-negativity by brute force: all minors up to size
/// `max_minor_size` are >= 0.
/// Convenience wrapper around [`check_total_positivity`]; see
/// [`is_totally_nonnegative_neville`] for the polynomial-time test of all minors.
pub fn is_totally_nonnegative_by_minors(mat: &[Vec<i64>], max_minor_size: usize) -> bool {
    check_total_positivity(mat, max_minor_size, false).is_ok()
}

/// Deprecated name of [`is_totally_nonnegative_by_minors`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `is_totally_nonnegative_by_minors`"
)]
pub fn is_totally_nonnegative(mat: &[Vec<i64>], max_minor_size: usize) -> bool {
    is_totally_nonnegative_by_minors(mat, max_minor_size)
}

// Neville elimination on the rows of an exact rational matrix: the shared
// kernel of the `i64` and `BigInt` entry points, which also apply it to the transpose.
fn check_totally_nonnegative_neville_rational_rows(
    a: &mut [Vec<Q>],
) -> Result<(), TotalPositivityError> {
    let nrows = a.len();
    let ncols = a.first().map_or(0, Vec::len);
    let min_dim = nrows.min(ncols);
    for k in 0..min_dim {
        for i in (k + 1..nrows).rev() {
            if a[i][k].is_zero() {
                continue;
            }
            if a[i - 1][k].is_zero() {
                // A positive entry below a zero pivot gives a negative 2x2
                // minor if the pivot row has any later positive entry. A
                // zero suffix row can be moved down without changing a minor.
                if ((k + 1)..ncols).any(|j| !a[i - 1][j].is_zero()) {
                    return Err(TotalPositivityError::PositiveEntryBelowZeroPivot {
                        row: i,
                        column: k,
                    });
                }
                a.swap(i - 1, i);
                continue;
            }
            let multiplier = a[i][k].clone() / a[i - 1][k].clone();
            debug_assert!(multiplier >= Q::zero());

            let (above, below) = a.split_at_mut(i);
            for (target, source) in below[0][k..ncols].iter_mut().zip(&above[i - 1][k..ncols]) {
                let sub = multiplier.clone() * source.clone();
                *target -= sub;
            }
            a[i][k] = Q::zero();

            for (j, entry) in a[i].iter().enumerate().take(ncols).skip(k + 1) {
                if *entry < Q::zero() {
                    return Err(TotalPositivityError::NegativeEliminationEntry {
                        row: i,
                        column: j,
                        pivot_column: k,
                    });
                }
            }
        }
    }
    Ok(())
}

/// Check total non-negativity of an integer matrix via Neville elimination
/// (adjacent-row subtraction), in exact rational arithmetic.
///
/// This is much faster than [`check_total_positivity`]: O(n²m) instead of
/// exponential in `min(n,m)`. It checks all minors, not just up to a given size.
/// The matrix is TNN iff the elimination, applied to the matrix and to its
/// transpose, completes with all multipliers ≥ 0 and no negative entry arises.
///
/// Returns `Ok(())` if the matrix is TNN, or a [`TotalPositivityError`]
/// describing the first failure: ragged rows
/// ([`TotalPositivityError::InconsistentRowLengths`]), a negative entry, or a
/// failed elimination step. An empty
/// matrix and a matrix with no columns are TNN. Use [`check_totally_nonnegative_neville_bigint`]
/// when entries may not fit in `i64`.
///
/// # Example
/// ```
/// use polytool::linalg::check_totally_nonnegative_neville;
/// let pascal = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_totally_nonnegative_neville(&pascal).is_ok());
/// ```
pub fn check_totally_nonnegative_neville(mat: &[Vec<i64>]) -> Result<(), TotalPositivityError> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err(TotalPositivityError::InconsistentRowLengths);
    }
    if ncols == 0 {
        return Ok(());
    }

    // Convert to rationals for exact arithmetic
    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| {
            row.iter()
                .map(|&v| Q::from_integer(BigInt::from(v)))
                .collect()
        })
        .collect();

    // Check all entries are non-negative
    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(TotalPositivityError::NegativeEntry {
                    row: i,
                    column: j,
                    value: Some(BigInt::from(mat[i][j])),
                });
            }
        }
    }

    check_totally_nonnegative_neville_rational_rows(&mut a)?;

    // Row elimination detects the initial minors involving the first column;
    // applying the same test to the transpose checks the complementary family
    // of initial minors. Both families are required in the Neville criterion.
    let mut transpose = (0..ncols)
        .map(|j| {
            (0..nrows)
                .map(|i| Q::from_integer(BigInt::from(mat[i][j])))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    check_totally_nonnegative_neville_rational_rows(&mut transpose)
}

/// Deprecated name of [`check_totally_nonnegative_neville`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `check_totally_nonnegative_neville`"
)]
pub fn check_tnn_neville(mat: &[Vec<i64>]) -> Result<(), TotalPositivityError> {
    check_totally_nonnegative_neville(mat)
}

/// Fast TNN check via Neville elimination. Returns bool.
/// Equivalent to `check_totally_nonnegative_neville(mat).is_ok()`.
pub fn is_totally_nonnegative_neville(mat: &[Vec<i64>]) -> bool {
    check_totally_nonnegative_neville(mat).is_ok()
}

/// Deprecated name of [`is_totally_nonnegative_neville`].
#[deprecated(since = "0.6.0", note = "renamed to `is_totally_nonnegative_neville`")]
pub fn is_tnn(mat: &[Vec<i64>]) -> bool {
    is_totally_nonnegative_neville(mat)
}

/// Check TNN via Neville elimination on a BigInt matrix.
/// Use this when entries may overflow i64.
pub fn check_totally_nonnegative_neville_bigint(
    mat: &[Vec<BigInt>],
) -> Result<(), TotalPositivityError> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err(TotalPositivityError::InconsistentRowLengths);
    }
    if ncols == 0 {
        return Ok(());
    }

    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect();

    for (i, row) in a.iter().enumerate() {
        for (j, entry) in row.iter().enumerate() {
            if *entry < Q::zero() {
                return Err(TotalPositivityError::NegativeEntry {
                    row: i,
                    column: j,
                    value: None,
                });
            }
        }
    }

    check_totally_nonnegative_neville_rational_rows(&mut a)?;
    let mut transpose = (0..ncols)
        .map(|j| {
            (0..nrows)
                .map(|i| Q::from_integer(mat[i][j].clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    check_totally_nonnegative_neville_rational_rows(&mut transpose)
}

/// Deprecated name of [`check_totally_nonnegative_neville_bigint`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `check_totally_nonnegative_neville_bigint`"
)]
pub fn check_tnn_neville_bigint(mat: &[Vec<BigInt>]) -> Result<(), TotalPositivityError> {
    check_totally_nonnegative_neville_bigint(mat)
}

fn extract_submatrix<T: Clone + Into<BigInt>>(
    mat: &[Vec<T>],
    rows: &[usize],
    cols: &[usize],
) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|&r| cols.iter().map(|&c| mat[r][c].clone().into()).collect())
        .collect()
}

struct CombinationIterator {
    n: usize,
    current: Option<Vec<usize>>,
}

impl CombinationIterator {
    fn new(n: usize, k: usize) -> Self {
        let current = (k <= n).then(|| (0..k).collect());
        Self { n, current }
    }
}

impl Iterator for CombinationIterator {
    type Item = Vec<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current.clone()?;
        let k = result.len();
        if k == 0 {
            self.current = None;
            return Some(result);
        }

        let current = self.current.as_mut().expect("combination exists");
        let Some(index) = (0..k).rev().find(|&i| current[i] < self.n - k + i) else {
            self.current = None;
            return Some(result);
        };
        current[index] += 1;
        for i in index + 1..k {
            current[i] = current[i - 1] + 1;
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::linalg::test_support::bi;

    // -----------------------------------------------------------------------
    // Total non-negativity
    // -----------------------------------------------------------------------

    #[test]
    fn bigint_minor_check_agrees_with_the_i64_check() {
        let matrices: Vec<Vec<Vec<i64>>> = vec![
            vec![],
            vec![vec![1, 0, 0], vec![1, 1, 0], vec![1, 2, 1], vec![1, 3, 3]],
            vec![vec![1, 2], vec![3, 4]],
            vec![vec![0, 0], vec![0, 1]],
            vec![vec![1, 1], vec![1]],
        ];
        for matrix in matrices {
            let big: Vec<Vec<BigInt>> = matrix
                .iter()
                .map(|row| row.iter().map(|&v| BigInt::from(v)).collect())
                .collect();
            for size in 0..4 {
                for strict in [false, true] {
                    assert_eq!(
                        check_total_positivity_bigint(&big, size, strict),
                        check_total_positivity(&matrix, size, strict),
                        "{matrix:?} size {size} strict {strict}"
                    );
                }
            }
        }
    }

    #[test]
    fn bigint_minor_check_handles_entries_beyond_i64() {
        // Binomial rows of (1 + 10^10 x)^n for n = 0..=3 have entries up to 10^30.
        let q = BigInt::from(10).pow(10);
        let rows: Vec<Vec<BigInt>> = (0..4usize)
            .map(|n| {
                (0..4usize)
                    .map(|k| {
                        if k > n {
                            return BigInt::from(0);
                        }
                        let mut binomial = BigInt::from(1);
                        for i in 0..k {
                            binomial = binomial * BigInt::from(n - i) / BigInt::from(i + 1);
                        }
                        binomial * q.pow(k as u32)
                    })
                    .collect()
            })
            .collect();
        assert!(rows[3][3] > BigInt::from(i64::MAX));
        assert!(check_total_positivity_bigint(&rows, 3, false).is_ok());

        // Negating one huge entry creates a negative 1 x 1 minor.
        let mut bad = rows.clone();
        bad[3][3] = -bad[3][3].clone();
        assert!(matches!(
            check_total_positivity_bigint(&bad, 1, false),
            Err(TotalPositivityError::NegativeMinor { size: 1, .. })
        ));
    }

    #[test]
    fn test_tnn_pascal() {
        let pascal = vec![
            vec![1, 0, 0],
            vec![1, 1, 0],
            vec![1, 2, 1],
            vec![1, 3, 3],
            vec![1, 4, 6],
        ];
        assert!(is_totally_nonnegative_neville(&pascal));
        assert!(check_totally_nonnegative_neville(&pascal).is_ok());
        assert!(check_total_positivity(&pascal, 3, false).is_ok());
    }

    #[test]
    fn test_tnn_identity() {
        let id = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        assert!(is_totally_nonnegative_neville(&id));
    }

    #[test]
    fn test_tnn_negative_entry() {
        let m = vec![vec![1, -1], vec![0, 1]];
        assert!(!is_totally_nonnegative_neville(&m));
    }

    #[test]
    fn test_tnn_rejects_row_swap_false_positive() {
        let permutation = vec![vec![0, 1], vec![1, 0]];
        assert!(check_totally_nonnegative_neville(&permutation).is_err());

        let permutation_bigint = vec![vec![bi(0), bi(1)], vec![bi(1), bi(0)]];
        assert!(check_totally_nonnegative_neville_bigint(&permutation_bigint).is_err());

        // Moving an all-zero suffix row remains safe and avoids rejecting a
        // genuinely TNN rank-deficient matrix.
        assert!(check_totally_nonnegative_neville(&[vec![0, 0], vec![1, 1]]).is_ok());
    }

    #[test]
    fn test_tnn_rejects_ragged_matrices() {
        assert_eq!(
            check_total_positivity(&[vec![1, 0], vec![1]], 2, false),
            Err(TotalPositivityError::RaggedMatrix {
                row: 1,
                expected_columns: 2,
                found_columns: 1,
            })
        );
        assert!(check_totally_nonnegative_neville(&[vec![1, 0], vec![1]]).is_err());
        assert!(
            check_totally_nonnegative_neville_bigint(&[vec![bi(1), bi(0)], vec![bi(1)]]).is_err()
        );
    }

    #[test]
    fn lazy_combinations_match_eager_reference() {
        fn eager(n: usize, k: usize) -> Vec<Vec<usize>> {
            fn generate(
                start: usize,
                remaining: usize,
                current: &mut Vec<usize>,
                output: &mut Vec<Vec<usize>>,
                n: usize,
            ) {
                if remaining == 0 {
                    output.push(current.clone());
                    return;
                }
                for value in start..=n - remaining {
                    current.push(value);
                    generate(value + 1, remaining - 1, current, output, n);
                    current.pop();
                }
            }

            if k > n {
                return Vec::new();
            }
            let mut output = Vec::new();
            generate(0, k, &mut Vec::new(), &mut output, n);
            output
        }

        for n in 0..=9 {
            for k in 0..=n + 1 {
                assert_eq!(
                    CombinationIterator::new(n, k).collect::<Vec<_>>(),
                    eager(n, k)
                );
            }
        }
    }

    #[test]
    fn test_tnn_not_tnn_positive_entries() {
        // All entries ≥ 0 but 2x2 minor det = 1*1 - 2*2 = -3 < 0
        let m = vec![vec![1, 2], vec![2, 1]];
        assert!(!is_totally_nonnegative_neville(&m));
    }

    #[test]
    fn test_tnn_path_matrix() {
        // Path matrix for P3: entry (i,j) = 1 if i <= j, 0 otherwise
        // These are always TNN
        let m = vec![vec![1, 1, 1], vec![0, 1, 1], vec![0, 0, 1]];
        assert!(is_totally_nonnegative_neville(&m));
    }

    #[test]
    fn test_tnn_neville_vs_brute() {
        // Verify Neville and brute-force agree on a 3x3 matrix
        let m = vec![vec![1, 1, 0], vec![0, 1, 1], vec![0, 0, 1]];
        let neville = check_totally_nonnegative_neville(&m).is_ok();
        let brute = check_total_positivity(&m, 3, false).is_ok();
        assert_eq!(neville, brute);
    }

    #[test]
    fn test_tnn_neville_matches_brute_force_on_binary_3x3_matrices() {
        for bits in 0_u16..(1 << 9) {
            let mut matrix = vec![vec![0; 3]; 3];
            for entry in 0..9 {
                matrix[entry / 3][entry % 3] = i64::from((bits >> entry) & 1);
            }
            assert_eq!(
                check_totally_nonnegative_neville(&matrix).is_ok(),
                check_total_positivity(&matrix, 3, false).is_ok(),
                "Neville and brute-force checks disagree for {matrix:?}"
            );
        }
    }

    #[test]
    fn test_tnn_bigint() {
        let m = vec![vec![bi(1), bi(1)], vec![bi(0), bi(1)]];
        assert!(check_totally_nonnegative_neville_bigint(&m).is_ok());
    }

    #[test]
    fn neville_error_messages_are_stable() {
        assert_eq!(
            check_totally_nonnegative_neville(&[vec![1, 2], vec![3]]).unwrap_err(),
            TotalPositivityError::InconsistentRowLengths
        );
        assert_eq!(
            check_totally_nonnegative_neville(&[vec![1, 2], vec![3]])
                .unwrap_err()
                .to_string(),
            "matrix rows have inconsistent lengths"
        );
        assert_eq!(
            check_totally_nonnegative_neville(&[vec![1, -2], vec![3, 4]])
                .unwrap_err()
                .to_string(),
            "Entry [0,1] = -2 < 0"
        );
        assert_eq!(
            check_totally_nonnegative_neville_bigint(&[vec![bi(1), bi(-2)], vec![bi(3), bi(4)]])
                .unwrap_err()
                .to_string(),
            "Entry [0,1] < 0"
        );
        assert_eq!(
            check_totally_nonnegative_neville(&[vec![0, 1], vec![1, 0]])
                .unwrap_err()
                .to_string(),
            "Neville elimination: positive entry [1,0] below zero pivot with nonzero \
             pivot-row suffix"
        );
        let error = check_totally_nonnegative_neville(&[vec![1, 2], vec![3, 4]]).unwrap_err();
        assert_eq!(
            error,
            TotalPositivityError::NegativeEliminationEntry {
                row: 1,
                column: 1,
                pivot_column: 0
            }
        );
        assert_eq!(
            error.to_string(),
            "Neville elimination: entry [1,1] became negative (multiplier at column 0, \
             rows [0,1])"
        );
    }
}
