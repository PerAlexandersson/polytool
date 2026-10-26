//! Total non-negativity: brute-force minor enumeration and Neville elimination.
//!
//! [`check_total_positivity`] enumerates minors up
//! to a given size (exponential); [`check_tnn_neville`]
//! and its `BigInt` variant decide total non-negativity of all minors in
//! `O(n^2 m)` exact rational operations.

use super::determinants::determinant;
use super::Q;
use num_bigint::BigInt;
use num_traits::Zero;

/// Exact failure reported by brute-force total-positivity checking.
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
                let sub = extract_submatrix_i64(mat, &rows, &cols);
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
                    let all_zero = rows.iter().all(|&r| cols.iter().all(|&c| mat[r][c] == 0));
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

/// Check total non-negativity: all minors >= 0.
/// Convenience wrapper around [`check_total_positivity`].
pub fn is_totally_nonnegative(mat: &[Vec<i64>], max_minor_size: usize) -> bool {
    check_total_positivity(mat, max_minor_size, false).is_ok()
}

// Neville elimination on the rows of an exact rational matrix: the shared
// kernel of the `i64` and `BigInt` entry points, which also apply it to the transpose.
fn check_tnn_neville_rational_rows(a: &mut [Vec<Q>]) -> Result<(), String> {
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
                    return Err(format!(
                        "Neville elimination: positive entry [{},{}] below zero pivot \
                         with nonzero pivot-row suffix",
                        i, k
                    ));
                }
                a.swap(i - 1, i);
                continue;
            }
            let multiplier = a[i][k].clone() / a[i - 1][k].clone();
            debug_assert!(multiplier >= Q::zero());

            for j in k..ncols {
                let sub = multiplier.clone() * a[i - 1][j].clone();
                a[i][j] -= sub;
            }
            a[i][k] = Q::zero();

            for j in (k + 1)..ncols {
                if a[i][j] < Q::zero() {
                    return Err(format!(
                        "Neville elimination: entry [{},{}] became negative \
                         (multiplier at column {}, rows [{},{}])",
                        i,
                        j,
                        k,
                        i - 1,
                        i
                    ));
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
/// Returns `Ok(())` if the matrix is TNN, or `Err(msg)` describing the first
/// failure: ragged rows, a negative entry, or a failed elimination step. An empty
/// matrix and a matrix with no columns are TNN. Use [`check_tnn_neville_bigint`]
/// when entries may not fit in `i64`.
///
/// # Example
/// ```
/// use polytool::linalg::check_tnn_neville;
/// let pascal = vec![
///     vec![1, 0, 0],
///     vec![1, 1, 0],
///     vec![1, 2, 1],
///     vec![1, 3, 3],
///     vec![1, 4, 6],
/// ];
/// assert!(check_tnn_neville(&pascal).is_ok());
/// ```
pub fn check_tnn_neville(mat: &[Vec<i64>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err("matrix rows have inconsistent lengths".to_string());
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
                return Err(format!("Entry [{},{}] = {} < 0", i, j, mat[i][j]));
            }
        }
    }

    check_tnn_neville_rational_rows(&mut a)?;

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
    check_tnn_neville_rational_rows(&mut transpose)
}

/// Fast TNN check via Neville elimination. Returns bool.
/// Equivalent to `check_tnn_neville(mat).is_ok()`.
pub fn is_tnn(mat: &[Vec<i64>]) -> bool {
    check_tnn_neville(mat).is_ok()
}

/// Check TNN via Neville elimination on a BigInt matrix.
/// Use this when entries may overflow i64.
pub fn check_tnn_neville_bigint(mat: &[Vec<BigInt>]) -> Result<(), String> {
    let nrows = mat.len();
    if nrows == 0 {
        return Ok(());
    }
    let ncols = mat[0].len();
    if mat.iter().any(|row| row.len() != ncols) {
        return Err("matrix rows have inconsistent lengths".to_string());
    }
    if ncols == 0 {
        return Ok(());
    }

    let mut a: Vec<Vec<Q>> = mat
        .iter()
        .map(|row| row.iter().map(|v| Q::from_integer(v.clone())).collect())
        .collect();

    for i in 0..nrows {
        for j in 0..ncols {
            if a[i][j] < Q::zero() {
                return Err(format!("Entry [{},{}] < 0", i, j));
            }
        }
    }

    check_tnn_neville_rational_rows(&mut a)?;
    let mut transpose = (0..ncols)
        .map(|j| {
            (0..nrows)
                .map(|i| Q::from_integer(mat[i][j].clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    check_tnn_neville_rational_rows(&mut transpose)
}

fn extract_submatrix_i64(mat: &[Vec<i64>], rows: &[usize], cols: &[usize]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|&r| cols.iter().map(|&c| BigInt::from(mat[r][c])).collect())
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
    fn test_tnn_pascal() {
        let pascal = vec![
            vec![1, 0, 0],
            vec![1, 1, 0],
            vec![1, 2, 1],
            vec![1, 3, 3],
            vec![1, 4, 6],
        ];
        assert!(is_tnn(&pascal));
        assert!(check_tnn_neville(&pascal).is_ok());
        assert!(check_total_positivity(&pascal, 3, false).is_ok());
    }

    #[test]
    fn test_tnn_identity() {
        let id = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
        assert!(is_tnn(&id));
    }

    #[test]
    fn test_tnn_negative_entry() {
        let m = vec![vec![1, -1], vec![0, 1]];
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_rejects_row_swap_false_positive() {
        let permutation = vec![vec![0, 1], vec![1, 0]];
        assert!(check_tnn_neville(&permutation).is_err());

        let permutation_bigint = vec![vec![bi(0), bi(1)], vec![bi(1), bi(0)]];
        assert!(check_tnn_neville_bigint(&permutation_bigint).is_err());

        // Moving an all-zero suffix row remains safe and avoids rejecting a
        // genuinely TNN rank-deficient matrix.
        assert!(check_tnn_neville(&[vec![0, 0], vec![1, 1]]).is_ok());
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
        assert!(check_tnn_neville(&[vec![1, 0], vec![1]]).is_err());
        assert!(check_tnn_neville_bigint(&[vec![bi(1), bi(0)], vec![bi(1)]]).is_err());
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
        assert!(!is_tnn(&m));
    }

    #[test]
    fn test_tnn_path_matrix() {
        // Path matrix for P3: entry (i,j) = 1 if i <= j, 0 otherwise
        // These are always TNN
        let m = vec![vec![1, 1, 1], vec![0, 1, 1], vec![0, 0, 1]];
        assert!(is_tnn(&m));
    }

    #[test]
    fn test_tnn_neville_vs_brute() {
        // Verify Neville and brute-force agree on a 3x3 matrix
        let m = vec![vec![1, 1, 0], vec![0, 1, 1], vec![0, 0, 1]];
        let neville = check_tnn_neville(&m).is_ok();
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
                check_tnn_neville(&matrix).is_ok(),
                check_total_positivity(&matrix, 3, false).is_ok(),
                "Neville and brute-force checks disagree for {matrix:?}"
            );
        }
    }

    #[test]
    fn test_tnn_bigint() {
        let m = vec![vec![bi(1), bi(1)], vec![bi(0), bi(1)]];
        assert!(check_tnn_neville_bigint(&m).is_ok());
    }
}
