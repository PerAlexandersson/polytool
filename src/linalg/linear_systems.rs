//! Exact linear systems and kernels over the rationals.
//!
//! Solving, rank/nullity reports and complete kernels with a rank certificate;
//! no modular or floating-point rank decisions are used.

use super::Q;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// A particular solution of a consistent exact linear system, together with
/// the rank information needed to decide whether it is identifiable.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearSystemSolution {
    /// One solution, with every free variable set to zero.
    pub particular: Vec<BigRational>,
    /// Rank of the coefficient matrix.
    pub rank: usize,
    /// Dimension of the solution space of the associated homogeneous system.
    pub nullity: usize,
}

/// Resource limits for a complete rational kernel computation.
#[derive(Debug, Clone, Copy)]
pub struct RationalKernelLimits {
    /// Maximum number of matrix entries (rows times columns); the kernel computation also
    /// requires `columns^2 <= max_entries / 2`.
    pub max_entries: usize,
    /// Maximum number of exact arithmetic operations during elimination.
    pub max_operations: usize,
    /// Maximum bit length of a numerator or denominator of any rational that occurs.
    pub max_coefficient_bits: u64,
}

/// Failure to compute a kernel, distinguished from a proved zero kernel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RationalKernelError {
    /// The rows do not all have the stated number of columns, or the rank certificate is
    /// inconsistent.
    InvalidDimensions,
    /// A resource limit was reached; the payload names it (`matrix entries`,
    /// `coefficient bits` or `arithmetic operations`). This is not a statement about the kernel.
    BudgetExceeded(&'static str),
}

impl std::fmt::Display for RationalKernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidDimensions => write!(
                f,
                "matrix rows do not all have the stated number of columns, or the rank \
                 certificate is inconsistent"
            ),
            Self::BudgetExceeded(limit) => {
                write!(f, "kernel computation exceeded its limit on {limit}")
            }
        }
    }
}

impl std::error::Error for RationalKernelError {}

/// A complete nullspace, with a nonsingular minor proving its dimension.
/// Basis vectors are stored separately (as columns mathematically).
#[derive(Debug, Clone)]
pub struct RationalKernel {
    /// Kernel basis; `basis[i]` has one entry per column of the matrix. Entry
    /// `free_columns[i]` of `basis[i]` is one and the other free-column entries are zero.
    pub basis: Vec<Vec<BigRational>>,
    /// Free columns in ascending order, one per basis vector.
    pub free_columns: Vec<usize>,
    /// Pivot columns in ascending order; their number is the rank.
    pub pivot_columns: Vec<usize>,
    /// Row indices of the input matrix forming an independent set, one per pivot column in
    /// the same order.
    pub independent_rows: Vec<usize>,
    /// Inverse of `A[independent_rows, pivot_columns]`.
    pub inverse_minor: Vec<Vec<BigRational>>,
}

struct EliminationBudget {
    remaining: usize,
    bits: u64,
}

impl EliminationBudget {
    fn check(&mut self, value: &Q) -> Result<(), RationalKernelError> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or(RationalKernelError::BudgetExceeded("arithmetic operations"))?;
        if value.numer().bits().max(value.denom().bits()) > self.bits {
            return Err(RationalKernelError::BudgetExceeded("coefficient bits"));
        }
        Ok(())
    }
}

// Shared elimination for the legacy solver and complete-kernel API. Columns
// after `columns` are carried through the row operations, never pivoted on.
fn rational_rref(
    aug: &mut [Vec<Q>],
    columns: usize,
    mut budget: Option<&mut EliminationBudget>,
) -> Result<(Vec<usize>, Vec<usize>), RationalKernelError> {
    let mut origins: Vec<_> = (0..aug.len()).collect();
    let mut pivots = Vec::new();
    for col in 0..columns {
        let row = pivots.len();
        let Some(pr) = (row..aug.len()).find(|&r| !aug[r][col].is_zero()) else {
            continue;
        };
        aug.swap(row, pr);
        origins.swap(row, pr);
        let value = aug[row][col].clone();
        let snapshot = aug[row][col..].to_vec();
        for (r, target_row) in aug.iter_mut().enumerate() {
            if r == row || target_row[col].is_zero() {
                continue;
            }
            let factor = &target_row[col] / &value;
            if let Some(b) = budget.as_deref_mut() {
                b.check(&factor)?;
            }
            for (entry, pivot) in target_row[col..].iter_mut().zip(&snapshot) {
                *entry -= pivot * &factor;
                if let Some(b) = budget.as_deref_mut() {
                    b.check(entry)?;
                }
            }
        }
        pivots.push(col);
    }
    origins.truncate(pivots.len());
    Ok((pivots, origins))
}

/// Compute every rational null vector, including a rank certificate.
/// `columns` is explicit so a zero-equation matrix has kernel Q^columns.
/// No modular or floating-point rank decisions are used.
pub fn rational_kernel_with_certificate(
    a: &[Vec<BigRational>],
    columns: usize,
    limits: RationalKernelLimits,
) -> Result<RationalKernel, RationalKernelError> {
    use RationalKernelError::{BudgetExceeded, InvalidDimensions};
    if a.iter().any(|r| r.len() != columns) {
        return Err(InvalidDimensions);
    }
    if a.len()
        .checked_mul(columns)
        .is_none_or(|n| n > limits.max_entries)
        || columns
            .checked_mul(columns)
            .is_none_or(|n| n > limits.max_entries / 2)
    {
        return Err(BudgetExceeded("matrix entries"));
    }
    if a.iter()
        .flatten()
        .any(|x| x.numer().bits().max(x.denom().bits()) > limits.max_coefficient_bits)
    {
        return Err(BudgetExceeded("coefficient bits"));
    }
    let mut budget = EliminationBudget {
        remaining: limits.max_operations,
        bits: limits.max_coefficient_bits,
    };
    let mut reduced = a.to_vec();
    let (pivots, rows) = rational_rref(&mut reduced, columns, Some(&mut budget))?;
    let rank = pivots.len();
    let free: Vec<_> = (0..columns).filter(|c| !pivots.contains(c)).collect();
    let mut basis = Vec::new();
    for &f in &free {
        let mut v = vec![Q::zero(); columns];
        v[f] = Q::one();
        for (i, &p) in pivots.iter().enumerate() {
            v[p] = -&reduced[i][f] / &reduced[i][p];
            budget.check(&v[p])?;
        }
        basis.push(v);
    }
    let mut minor: Vec<Vec<Q>> = rows
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            let mut row: Vec<_> = pivots.iter().map(|&c| a[r][c].clone()).collect();
            row.extend((0..rank).map(|j| if i == j { Q::one() } else { Q::zero() }));
            row
        })
        .collect();
    let (minor_pivots, _) = rational_rref(&mut minor, rank, Some(&mut budget))?;
    if minor_pivots.len() != rank {
        return Err(InvalidDimensions);
    }
    let mut inverse = Vec::new();
    for (i, row) in minor.iter().enumerate() {
        let mut out = Vec::new();
        for x in &row[rank..] {
            let v = x / &row[i];
            budget.check(&v)?;
            out.push(v);
        }
        inverse.push(out);
    }
    Ok(RationalKernel {
        basis,
        free_columns: free,
        pivot_columns: pivots,
        independent_rows: rows,
        inverse_minor: inverse,
    })
}

/// Solve `Ax = b` over ℚ and report the rank and nullity of `A`.
///
/// Returns `None` if the dimensions are invalid or the system is
/// inconsistent. For an underdetermined consistent system, free variables in
/// [`LinearSystemSolution::particular`] are set to zero.
pub fn solve_linear_system_with_rank(
    a: &[Vec<BigRational>],
    b: &[BigRational],
) -> Option<LinearSystemSolution> {
    let num_rows = a.len();
    if b.len() != num_rows {
        return None;
    }
    if num_rows == 0 {
        return Some(LinearSystemSolution {
            particular: vec![],
            rank: 0,
            nullity: 0,
        });
    }
    let num_cols = a[0].len();
    if a.iter().any(|row| row.len() != num_cols) {
        return None;
    }

    // Build augmented matrix [A | b].
    let mut aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    let (pivot_columns, _) = rational_rref(&mut aug, num_cols, None).ok()?;
    let pivot_row = pivot_columns.len();

    // Consistency check: any row [0 ... 0 | nonzero] means no solution.
    for aug_row in &aug[pivot_row..num_rows] {
        if !aug_row[num_cols].is_zero() {
            return None;
        }
    }

    // Extract solution (free variables = 0).
    let mut x = vec![Q::zero(); num_cols];
    for (pr, &pc) in pivot_columns.iter().enumerate() {
        x[pc] = aug[pr][num_cols].clone() / aug[pr][pc].clone();
    }
    let rank = pivot_columns.len();
    Some(LinearSystemSolution {
        particular: x,
        rank,
        nullity: num_cols - rank,
    })
}

/// Solve Ax = b via Gaussian elimination over ℚ.
///
/// Returns `None` if the system is inconsistent. For underdetermined systems,
/// free variables are set to zero. Use [`solve_linear_system_with_rank`] when
/// the caller needs to distinguish a unique solution from an arbitrary
/// particular solution.
pub fn solve_linear_system(a: &[Vec<BigRational>], b: &[BigRational]) -> Option<Vec<BigRational>> {
    solve_linear_system_with_rank(a, b).map(|solution| solution.particular)
}

/// Solve a nonsingular square system `Ax = b` over `Q`.
///
/// This uses ordinary forward elimination followed by back-substitution.  It is
/// faster than [`solve_linear_system`] when the caller already knows that the
/// chosen rows form a full-rank square subsystem, because it does not compute a
/// full reduced row-echelon form.  Returns `None` if the input is not square or
/// if a zero pivot proves that the matrix is singular.
pub fn solve_full_rank_square_linear_system(
    a: &[Vec<BigRational>],
    b: &[BigRational],
) -> Option<Vec<BigRational>> {
    let n = a.len();
    if b.len() != n || a.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(vec![]);
    }

    if let Some(integer_aug) = square_rational_system_as_bigint_augmented(a, b) {
        if let Some(solution) = solve_full_rank_square_integer_augmented(integer_aug) {
            return Some(solution);
        }
    }

    let mut aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    for col in 0..n {
        let pivot_row = (col..n).find(|&row| !aug[row][col].is_zero())?;
        aug.swap(col, pivot_row);

        let pivot_val = aug[col][col].clone();
        let (pivot_block, lower_rows) = aug.split_at_mut(col + 1);
        let pivot_row = &pivot_block[col];
        for row in lower_rows {
            if row[col].is_zero() {
                continue;
            }
            let factor = row[col].clone() / &pivot_val;
            row[col] = Q::zero();
            for (aug_j, pivot_j) in row[col + 1..=n]
                .iter_mut()
                .zip(pivot_row[col + 1..=n].iter())
            {
                let sub = pivot_j.clone() * &factor;
                *aug_j -= sub;
            }
        }
    }

    let mut x = vec![Q::zero(); n];
    for row in (0..n).rev() {
        let mut value = aug[row][n].clone();
        for (col, coeff) in aug[row].iter().enumerate().take(n).skip(row + 1) {
            if !coeff.is_zero() {
                value -= coeff.clone() * &x[col];
            }
        }
        x[row] = value / &aug[row][row];
    }

    Some(x)
}

fn square_rational_system_as_bigint_augmented(a: &[Vec<Q>], b: &[Q]) -> Option<Vec<Vec<BigInt>>> {
    a.iter()
        .zip(b.iter())
        .map(|row| {
            let (row, rhs) = row;
            let mut augmented = row
                .iter()
                .map(|entry| entry.denom().is_one().then(|| entry.numer().clone()))
                .collect::<Option<Vec<_>>>()?;
            augmented.push(rhs.denom().is_one().then(|| rhs.numer().clone())?);
            Some(augmented)
        })
        .collect()
}

fn solve_full_rank_square_integer_augmented(mut aug: Vec<Vec<BigInt>>) -> Option<Vec<Q>> {
    let n = aug.len();
    if aug.iter().any(|row| row.len() != n + 1) {
        return None;
    }
    if n == 0 {
        return Some(vec![]);
    }

    let mut denom = BigInt::one();

    for col in 0..n {
        let pivot_row = (col..n).find(|&row| !aug[row][col].is_zero())?;
        aug.swap(col, pivot_row);
        let pivot = aug[col][col].clone();
        if pivot.is_zero() {
            return None;
        }
        if col == n - 1 {
            break;
        }

        let pivot_tail = aug[col][col + 1..=n].to_vec();
        for target_row in &mut aug[(col + 1)..n] {
            let row_pivot = target_row[col].clone();
            target_row[col] = BigInt::zero();
            for (cell, tail) in target_row[(col + 1)..=n].iter_mut().zip(&pivot_tail) {
                let numerator = if row_pivot.is_zero() {
                    &*cell * &pivot
                } else {
                    &*cell * &pivot - &row_pivot * tail
                };
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                *cell = numerator / &denom;
            }
        }
        denom = pivot;
    }

    let mut x = vec![Q::zero(); n];
    for row in (0..n).rev() {
        let diagonal = aug[row][row].clone();
        if diagonal.is_zero() {
            return None;
        }
        let mut value = Q::from_integer(aug[row][n].clone());
        for (col, coeff) in aug[row].iter().enumerate().take(n).skip(row + 1) {
            if !coeff.is_zero() {
                value -= Q::from_integer(coeff.clone()) * &x[col];
            }
        }
        x[row] = value / Q::from_integer(diagonal);
    }

    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::linalg::test_support::bi;

    // -----------------------------------------------------------------------
    // Solve linear system
    // -----------------------------------------------------------------------

    #[test]
    fn linear_system_reports_rank_and_nullity() {
        let a = vec![
            vec![Q::one(), Q::one(), Q::zero()],
            vec![
                Q::from_integer(BigInt::from(2)),
                Q::from_integer(BigInt::from(2)),
                Q::zero(),
            ],
        ];
        let b = vec![Q::one(), Q::from_integer(BigInt::from(2))];
        let solution = solve_linear_system_with_rank(&a, &b).unwrap();

        assert_eq!(solution.rank, 1);
        assert_eq!(solution.nullity, 2);
        assert_eq!(solution.particular, vec![Q::one(), Q::zero(), Q::zero()]);
    }

    #[test]
    fn test_solve_identity() {
        let a = vec![vec![Q::one(), Q::zero()], vec![Q::zero(), Q::one()]];
        let b = vec![Q::from_integer(bi(3)), Q::from_integer(bi(5))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(3)));
        assert_eq!(x[1], Q::from_integer(bi(5)));
    }

    #[test]
    fn test_solve_2x2() {
        // x + 2y = 5, 3x + 4y = 11 -> x = 1, y = 2
        let a = vec![
            vec![Q::from_integer(bi(1)), Q::from_integer(bi(2))],
            vec![Q::from_integer(bi(3)), Q::from_integer(bi(4))],
        ];
        let b = vec![Q::from_integer(bi(5)), Q::from_integer(bi(11))];
        let x = solve_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(1)));
        assert_eq!(x[1], Q::from_integer(bi(2)));
    }

    #[test]
    fn test_solve_full_rank_square_2x2() {
        // x + 2y = 5, 3x + 4y = 11 -> x = 1, y = 2
        let a = vec![
            vec![Q::from_integer(bi(1)), Q::from_integer(bi(2))],
            vec![Q::from_integer(bi(3)), Q::from_integer(bi(4))],
        ];
        let b = vec![Q::from_integer(bi(5)), Q::from_integer(bi(11))];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(x[0], Q::from_integer(bi(1)));
        assert_eq!(x[1], Q::from_integer(bi(2)));
    }

    #[test]
    fn test_solve_full_rank_square_integer_row_swap() {
        // Solution is (3, -1, 2); the first pivot requires a row swap.
        let a = vec![
            vec![
                Q::from_integer(bi(0)),
                Q::from_integer(bi(2)),
                Q::from_integer(bi(1)),
            ],
            vec![
                Q::from_integer(bi(1)),
                Q::from_integer(bi(1)),
                Q::from_integer(bi(0)),
            ],
            vec![
                Q::from_integer(bi(2)),
                Q::from_integer(bi(3)),
                Q::from_integer(bi(1)),
            ],
        ];
        let b = vec![
            Q::from_integer(bi(0)),
            Q::from_integer(bi(2)),
            Q::from_integer(bi(5)),
        ];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(
            x,
            vec![
                Q::from_integer(bi(3)),
                Q::from_integer(bi(-1)),
                Q::from_integer(bi(2)),
            ]
        );
    }

    #[test]
    fn test_solve_full_rank_square_rational_fallback() {
        // x - y = 0, x/2 + y = 2 -> x = y = 4/3.
        let half = Q::new(bi(1), bi(2));
        let a = vec![
            vec![half, Q::one()],
            vec![Q::one(), Q::from_integer(bi(-1))],
        ];
        let b = vec![Q::from_integer(bi(2)), Q::zero()];
        let x = solve_full_rank_square_linear_system(&a, &b).unwrap();
        assert_eq!(x, vec![Q::new(bi(4), bi(3)), Q::new(bi(4), bi(3))]);
    }

    #[test]
    fn test_solve_full_rank_square_rejects_singular() {
        let a = vec![vec![Q::one(), Q::one()], vec![Q::one(), Q::one()]];
        let b = vec![Q::one(), Q::one()];
        assert!(solve_full_rank_square_linear_system(&a, &b).is_none());
    }

    #[test]
    fn test_solve_inconsistent() {
        // x + y = 1, x + y = 2 -> no solution
        let a = vec![vec![Q::one(), Q::one()], vec![Q::one(), Q::one()]];
        let b = vec![Q::one(), Q::from_integer(bi(2))];
        assert!(solve_linear_system(&a, &b).is_none());
    }

    #[test]
    fn test_solve_rejects_mismatched_shapes() {
        assert!(solve_linear_system(&[vec![Q::one()]], &[]).is_none());
        assert!(solve_linear_system(&[], &[Q::one()]).is_none());
        assert!(solve_linear_system(
            &[vec![Q::one(), Q::zero()], vec![Q::one()]],
            &[Q::one(), Q::one()]
        )
        .is_none());
    }

    #[test]
    fn kernel_errors_display_and_implement_error() {
        let error: Box<dyn std::error::Error> = Box::new(RationalKernelError::InvalidDimensions);
        assert!(error.to_string().contains("number of columns"));
        let limits = RationalKernelLimits {
            max_entries: 1,
            max_operations: 10,
            max_coefficient_bits: 64,
        };
        let error =
            rational_kernel_with_certificate(&[vec![Q::one(), Q::one()]], 2, limits).unwrap_err();
        assert_eq!(error, RationalKernelError::BudgetExceeded("matrix entries"));
        assert_eq!(
            error.to_string(),
            "kernel computation exceeded its limit on matrix entries"
        );
    }
}
