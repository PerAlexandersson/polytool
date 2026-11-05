//! Exact linear systems and kernels over the rationals.
//!
//! Solving, rank/nullity reports and complete kernels with a rank certificate;
//! no modular or floating-point rank decisions are used.

use super::integer_elimination::{
    back_substitute, integer_echelon_for_solve, integer_row, integer_rref, Elimination,
    EliminationBudget, IntegerEliminationResult, IntegerRref,
};
use super::Q;
use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
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

/// Gauss-Jordan elimination of `aug` (columns after `columns` are carried through the
/// row operations, never pivoted on), on integer rows (see `integer_elimination`).
///
/// Every returned row is a nonzero rational multiple of the corresponding row of the
/// unnormalised reduced row-echelon form, so ratios of entries of one row are those
/// of the rational elimination; the pivot columns and rows and the budget accounting
/// are the same as those of [`rational_rref`], which is the fallback should an exact
/// division of the integer elimination ever fail.
fn gauss_jordan(
    aug: &[Vec<Q>],
    columns: usize,
    mut budget: Option<&mut EliminationBudget>,
) -> Result<IntegerRref, RationalKernelError> {
    let remaining = budget.as_deref().map(|b| b.remaining);
    if let Some(reduced) = integer_rref(
        aug,
        columns,
        Elimination::GaussJordan,
        budget.as_deref_mut(),
    )? {
        return Ok(reduced);
    }
    if let (Some(b), Some(remaining)) = (budget.as_deref_mut(), remaining) {
        b.remaining = remaining;
    }
    let mut rows = aug.to_vec();
    let (pivots, origins) = rational_rref(&mut rows, columns, budget)?;
    Ok(IntegerRref {
        rows: rows.iter().map(|row| integer_row(row).0).collect(),
        pivots,
        origins,
    })
}

// The elimination on rational entries, with a gcd reduction after every operation.
// It is the reference of `gauss_jordan` and its fallback.
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
    let reduced = gauss_jordan(a, columns, Some(&mut budget))?;
    let IntegerRref {
        rows: reduced,
        pivots,
        origins: rows,
    } = reduced;
    let rank = pivots.len();
    let free: Vec<_> = (0..columns).filter(|c| !pivots.contains(c)).collect();
    let mut basis = Vec::new();
    for &f in &free {
        let mut v = vec![Q::zero(); columns];
        v[f] = Q::one();
        for (i, &p) in pivots.iter().enumerate() {
            v[p] = Q::new(-&reduced[i][f], reduced[i][p].clone());
            budget.check(&v[p])?;
        }
        basis.push(v);
    }
    let minor: Vec<Vec<Q>> = rows
        .iter()
        .enumerate()
        .map(|(i, &r)| {
            let mut row: Vec<_> = pivots.iter().map(|&c| a[r][c].clone()).collect();
            row.extend((0..rank).map(|j| if i == j { Q::one() } else { Q::zero() }));
            row
        })
        .collect();
    let minor = gauss_jordan(&minor, rank, Some(&mut budget))?;
    if minor.pivots.len() != rank {
        return Err(InvalidDimensions);
    }
    let mut inverse = Vec::new();
    for (i, row) in minor.rows.iter().enumerate() {
        let mut out = Vec::new();
        for x in &row[rank..] {
            let v = Q::new(x.clone(), row[i].clone());
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
    let aug: Vec<Vec<Q>> = a
        .iter()
        .zip(b.iter())
        .map(|(row, bi)| {
            let mut r = row.clone();
            r.push(bi.clone());
            r
        })
        .collect();

    // The solution, the rank and the consistency do not depend on the pivot rows, so
    // the integer forward elimination with back substitution decides them.
    if let Some(outcome) = integer_echelon_for_solve(&aug, num_cols) {
        let echelon = match outcome {
            IntegerEliminationResult::Inconsistent => return None,
            IntegerEliminationResult::Echelon(echelon) => echelon,
        };
        let rank = echelon.pivots.len();
        if echelon.rows[rank..]
            .iter()
            .any(|row| !row[num_cols].is_zero())
        {
            return None;
        }
        if let Some(particular) = back_substitute(&echelon, num_cols) {
            return Some(LinearSystemSolution {
                particular,
                rank,
                nullity: num_cols - rank,
            });
        }
    }

    let IntegerRref {
        rows: reduced,
        pivots: pivot_columns,
        ..
    } = gauss_jordan(&aug, num_cols, None).ok()?;
    let pivot_row = pivot_columns.len();

    // Consistency check: any row [0 ... 0 | nonzero] means no solution.
    for reduced_row in &reduced[pivot_row..num_rows] {
        if !reduced_row[num_cols].is_zero() {
            return None;
        }
    }

    // Extract solution (free variables = 0).
    let mut x = vec![Q::zero(); num_cols];
    for (pr, &pc) in pivot_columns.iter().enumerate() {
        x[pc] = Q::new(reduced[pr][num_cols].clone(), reduced[pr][pc].clone());
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
/// Each equation's denominators are cleared once, including its right-hand
/// side. Fraction-free forward elimination and integer back substitution then
/// reduce rationals only for the final solution coordinates. Rational
/// elimination remains the fallback if the integer solver cannot finish.
/// Returns `None` if the input is not square or if a zero pivot proves that the
/// matrix is singular.
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

    let integer_aug = square_rational_system_as_bigint_augmented(a, b);
    if let Some(solution) = solve_full_rank_square_integer_augmented(integer_aug) {
        return Some(solution);
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

fn square_rational_system_as_bigint_augmented(a: &[Vec<Q>], b: &[Q]) -> Vec<Vec<BigInt>> {
    a.iter()
        .zip(b.iter())
        .map(|(row, rhs)| {
            // Keep the integer-input path free of temporary rational copies.
            let integer = row
                .iter()
                .chain(std::iter::once(rhs))
                .map(|entry| entry.denom().is_one().then(|| entry.numer().clone()))
                .collect::<Option<Vec<_>>>();
            integer.unwrap_or_else(|| {
                let augmented: Vec<_> = row.iter().cloned().chain([rhs.clone()]).collect();
                integer_row(&augmented).0
            })
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
                if denom.is_one() {
                    *cell = numerator;
                } else {
                    let (quotient, remainder) = numerator.div_rem(&denom);
                    if !remainder.is_zero() {
                        return None;
                    }
                    *cell = quotient;
                }
            }
        }
        denom = pivot;
    }

    back_substitute(
        &IntegerRref {
            rows: aug,
            pivots: (0..n).collect(),
            origins: Vec::new(),
        },
        n,
    )
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
    fn test_solve_full_rank_square_rational_inputs() {
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

    /// A deterministic pseudo-random matrix with small entries (some zero), some
    /// rows repeated or combined so that the rank is deficient, and optional
    /// denominators.
    fn pseudo_random_matrix(
        seed: u64,
        rows: usize,
        columns: usize,
        fractions: bool,
    ) -> Vec<Vec<Q>> {
        let mut state = seed;
        let mut next = move |modulus: u64| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (state >> 33) % modulus
        };
        let mut matrix: Vec<Vec<Q>> = Vec::new();
        for r in 0..rows {
            let row = if r >= 2 && next(4) == 0 {
                // A combination of two earlier rows.
                let (a, b) = (next(r as u64) as usize, next(r as u64) as usize);
                let (x, y) = (next(5) as i64 - 2, next(5) as i64 - 2);
                (0..columns)
                    .map(|c| {
                        &matrix[a][c] * Q::from_integer(bi(x))
                            + &matrix[b][c] * Q::from_integer(bi(y))
                    })
                    .collect()
            } else {
                (0..columns)
                    .map(|_| {
                        if next(3) == 0 {
                            Q::zero()
                        } else if fractions && next(3) == 0 {
                            Q::new(bi(next(9) as i64 - 4), bi(1 + next(4) as i64))
                        } else {
                            Q::from_integer(bi(next(11) as i64 - 5))
                        }
                    })
                    .collect()
            };
            matrix.push(row);
        }
        matrix
    }

    /// Rows, pivot columns and pivot-row origins of the rational reference elimination.
    type ReferenceReduction = (Vec<Vec<Q>>, Vec<usize>, Vec<usize>);

    /// The same value as the integer elimination, from the rational reference.
    fn reference_reduction(
        matrix: &[Vec<Q>],
        columns: usize,
        budget: Option<&mut EliminationBudget>,
    ) -> Result<ReferenceReduction, RationalKernelError> {
        let mut rows = matrix.to_vec();
        let (pivots, origins) = rational_rref(&mut rows, columns, budget)?;
        Ok((rows, pivots, origins))
    }

    #[test]
    fn integer_elimination_matches_the_rational_reduction() {
        let mut fallbacks = 0;
        for seed in 0..300u64 {
            let rows = 1 + (seed % 9) as usize;
            let columns = 1 + (seed / 3 % 8) as usize;
            let matrix = pseudo_random_matrix(seed, rows, columns + 1, seed % 2 == 0);
            let (reference, pivots, origins) = reference_reduction(&matrix, columns, None).unwrap();
            let Some(fast) =
                integer_rref(&matrix, columns, Elimination::GaussJordan, None).unwrap()
            else {
                fallbacks += 1;
                continue;
            };
            assert_eq!(fast.pivots, pivots, "seed {seed}");
            assert_eq!(fast.origins, origins, "seed {seed}");
            // Each integer row is a nonzero multiple of the reference row.
            for (integer, exact) in fast.rows.iter().zip(&reference) {
                let ratio = integer.iter().zip(exact).find(|(_, e)| !e.is_zero());
                match ratio {
                    None => assert!(integer.iter().all(Zero::is_zero), "seed {seed}"),
                    Some((i, e)) => {
                        let scale = Q::from_integer(i.clone()) / e;
                        assert!(!scale.is_zero(), "seed {seed}");
                        for (i, e) in integer.iter().zip(exact) {
                            assert_eq!(Q::from_integer(i.clone()), e * &scale, "seed {seed}");
                        }
                    }
                }
            }
        }
        assert_eq!(fallbacks, 0, "the exact divisions never fail");
    }

    #[test]
    fn integer_elimination_charges_the_same_budget_as_the_rational_reduction() {
        for seed in 0..60u64 {
            let rows = 2 + (seed % 6) as usize;
            let columns = 2 + (seed / 2 % 5) as usize;
            let matrix = pseudo_random_matrix(1000 + seed, rows, columns, seed % 3 != 0);
            let unlimited = |bits| EliminationBudget {
                remaining: usize::MAX,
                bits,
            };
            // Smallest sufficient operation count and bit limit, for both.
            let minimal = |reference: bool, by_operations: bool| {
                let succeeds = |value: u64| {
                    let mut budget = if by_operations {
                        EliminationBudget {
                            remaining: value as usize,
                            bits: u64::MAX,
                        }
                    } else {
                        unlimited(value)
                    };
                    if reference {
                        rational_rref(&mut matrix.clone(), columns, Some(&mut budget)).is_ok()
                    } else {
                        integer_rref(
                            &matrix,
                            columns,
                            Elimination::GaussJordan,
                            Some(&mut budget),
                        )
                        .is_ok()
                    }
                };
                let (mut low, mut high) = (0u64, 1_000_000u64);
                while low < high {
                    let mid = low + (high - low) / 2;
                    if succeeds(mid) {
                        high = mid;
                    } else {
                        low = mid + 1;
                    }
                }
                low
            };
            for by_operations in [true, false] {
                assert_eq!(
                    minimal(false, by_operations),
                    minimal(true, by_operations),
                    "seed {seed}, by operations: {by_operations}"
                );
            }
            // The error reported below the limit is the same as well.
            let run = |reference: bool, remaining: usize, bits: u64| {
                let mut budget = EliminationBudget { remaining, bits };
                if reference {
                    rational_rref(&mut matrix.clone(), columns, Some(&mut budget)).err()
                } else {
                    integer_rref(
                        &matrix,
                        columns,
                        Elimination::GaussJordan,
                        Some(&mut budget),
                    )
                    .err()
                }
            };
            for remaining in [0usize, 3, 10, 30, 100] {
                for bits in [1u64, 3, 5, 8, 20, 64] {
                    assert_eq!(
                        run(false, remaining, bits),
                        run(true, remaining, bits),
                        "seed {seed}, remaining {remaining}, bits {bits}"
                    );
                }
            }
        }
    }

    fn q(value: i64) -> Q {
        Q::from_integer(bi(value))
    }

    /// The solver as it was before the integer elimination: rational Gauss-Jordan
    /// elimination, the solution read off the unnormalised pivot rows.
    fn reference_solution(a: &[Vec<Q>], b: &[Q]) -> Option<LinearSystemSolution> {
        let columns = a[0].len();
        let mut aug: Vec<Vec<Q>> = a
            .iter()
            .zip(b)
            .map(|(row, bi)| row.iter().cloned().chain([bi.clone()]).collect())
            .collect();
        let (pivots, _) = rational_rref(&mut aug, columns, None).ok()?;
        if aug[pivots.len()..]
            .iter()
            .any(|row| !row[columns].is_zero())
        {
            return None;
        }
        let mut x = vec![Q::zero(); columns];
        for (row, &column) in pivots.iter().enumerate() {
            x[column] = &aug[row][columns] / &aug[row][column];
        }
        Some(LinearSystemSolution {
            particular: x,
            rank: pivots.len(),
            nullity: columns - pivots.len(),
        })
    }

    #[test]
    fn exact_contradictions_stop_solve_but_not_full_reduction() {
        for (a, b) in [
            (vec![vec![q(0), q(0)]], vec![q(1)]),
            (vec![vec![q(1), q(1)], vec![q(2), q(2)]], vec![q(3), q(7)]),
            (
                vec![vec![q(0), q(1) / q(2)], vec![q(0), q(1) / q(3)]],
                vec![q(1), q(1)],
            ),
        ] {
            let columns = a[0].len();
            let augmented = a
                .iter()
                .zip(&b)
                .map(|(row, rhs)| row.iter().cloned().chain([rhs.clone()]).collect::<Vec<_>>())
                .collect::<Vec<_>>();
            assert!(matches!(
                integer_echelon_for_solve(&augmented, columns),
                Some(IntegerEliminationResult::Inconsistent)
            ));
            assert_eq!(
                solve_linear_system_with_rank(&a, &b),
                reference_solution(&a, &b)
            );
            // Kernel/certificate consumers still receive complete reduction,
            // including when the trailing column is not a system RHS.
            let full = integer_rref(&augmented, columns, Elimination::GaussJordan, None)
                .unwrap()
                .unwrap();
            let (expected, pivots, _) = reference_reduction(&augmented, columns, None).unwrap();
            assert_eq!(full.pivots, pivots);
            assert_eq!(full.rows.len(), expected.len());
        }
        // Zero rows and rank-deficient consistent systems must not be rejected.
        let a = vec![vec![q(0), q(0)], vec![q(1), q(2)], vec![q(2), q(4)]];
        let b = vec![q(0), q(3), q(6)];
        assert_eq!(
            solve_linear_system_with_rank(&a, &b),
            reference_solution(&a, &b)
        );
    }

    #[test]
    fn sparse_solves_with_delayed_row_updates_stay_fraction_free() {
        // Skipped pivot columns leave rows at different integer levels. Include
        // negative pivots, rational row scalings, free columns and dependent rows.
        let base = [
            [2, 0, 0, 3, 0, 1],
            [0, 0, -3, 0, 2, 1],
            [0, 0, 0, 5, 0, 2],
            [0, 0, 0, 0, 0, 7],
            [0, 0, 0, 10, 0, 4],
            [0, 0, 0, 0, 0, 0],
        ];
        for (rotation, shift) in
            (0..base.len()).flat_map(|rotation| [0, 96].map(move |shift| (rotation, shift)))
        {
            let scale = Q::from_integer(BigInt::one() << shift);
            let mut a: Vec<Vec<Q>> = base
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    row.iter()
                        .map(|&v| q(v) * &scale / q(i as i64 + 2))
                        .collect()
                })
                .collect();
            a.rotate_left(rotation);
            let b: Vec<Q> = a
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .map(|(j, x)| x / q(j as i64 + 2))
                        .sum()
                })
                .collect();
            let augmented: Vec<Vec<Q>> = a
                .iter()
                .zip(&b)
                .map(|(row, rhs)| row.iter().cloned().chain([rhs.clone()]).collect())
                .collect();
            let Some(IntegerEliminationResult::Echelon(echelon)) =
                integer_echelon_for_solve(&augmented, 6)
            else {
                panic!("consistent delayed updates must not need the rational fallback");
            };
            let expected = reference_solution(&a, &b).unwrap();
            assert_eq!(echelon.pivots, [0, 2, 3, 5]);
            assert_eq!(back_substitute(&echelon, 6), Some(expected.particular));
        }
    }

    #[test]
    fn square_integer_back_substitution_matches_rational_reference() {
        for seed in 0..100u64 {
            let n = 1 + (seed % 10) as usize;
            let a = pseudo_random_matrix(7000 + seed, n, n, false);
            let b: Vec<_> = (0..n).map(|i| q(i as i64 - 3)).collect();
            let expected = reference_solution(&a, &b)
                .filter(|solution| solution.rank == n)
                .map(|solution| solution.particular);
            let integer_aug = square_rational_system_as_bigint_augmented(&a, &b);
            assert_eq!(
                solve_full_rank_square_integer_augmented(integer_aug),
                expected,
                "seed {seed}"
            );
        }
    }

    #[test]
    fn square_denominator_clearing_includes_rhs_and_preserves_row_swaps() {
        let a = vec![
            vec![q(0), q(1) / q(2), -q(1) / q(3)],
            vec![q(2) / q(5), q(0), q(3) / q(7)],
            vec![q(1), q(0), q(1) / q(11)],
        ];
        let b = vec![q(1) / q(13), -q(2) / q(17), q(3) / q(19)];
        let augmented = square_rational_system_as_bigint_augmented(&a, &b);
        assert_eq!(augmented[0], [bi(0), bi(39), bi(-26), bi(6)]);
        let expected = reference_solution(&a, &b).unwrap().particular;
        assert_eq!(
            solve_full_rank_square_integer_augmented(augmented),
            Some(expected.clone())
        );
        assert_eq!(solve_full_rank_square_linear_system(&a, &b), Some(expected));

        // Only the RHS needs denominator clearing in this one-equation system.
        assert_eq!(
            solve_full_rank_square_linear_system(&[vec![q(2)]], &[q(1) / q(3)]),
            Some(vec![q(1) / q(6)])
        );
    }

    #[test]
    fn square_rational_inputs_match_reference_without_rational_fallback() {
        let (mut nonsingular, mut singular) = (0, 0);
        for seed in 0..200u64 {
            let n = 1 + (seed % 10) as usize;
            let a = pseudo_random_matrix(8000 + seed, n, n, true);
            let b: Vec<_> = (0..n).map(|i| q(i as i64 - 3) / q(i as i64 + 2)).collect();
            let expected = reference_solution(&a, &b)
                .filter(|solution| solution.rank == n)
                .map(|solution| solution.particular);
            if expected.is_some() {
                nonsingular += 1;
            } else {
                singular += 1;
            }
            assert_eq!(
                solve_full_rank_square_integer_augmented(
                    square_rational_system_as_bigint_augmented(&a, &b)
                ),
                expected,
                "seed {seed}"
            );
            assert_eq!(solve_full_rank_square_linear_system(&a, &b), expected);
        }
        assert!(nonsingular > 0 && singular > 0);
    }

    #[test]
    fn solving_matches_the_rational_reference() {
        let (mut consistent, mut inconsistent) = (0, 0);
        for seed in 0..400u64 {
            let rows = 1 + (seed % 10) as usize;
            let columns = 1 + (seed / 4 % 8) as usize;
            let a = pseudo_random_matrix(5000 + seed, rows, columns, seed % 2 == 1);
            // Right-hand sides: one in the column space, one arbitrary.
            let witness: Vec<Q> = (0..columns)
                .map(|c| q((c as i64 * 7 + seed as i64) % 5 - 2))
                .collect();
            let in_span: Vec<Q> = a
                .iter()
                .map(|row| row.iter().zip(&witness).map(|(x, y)| x * y).sum())
                .collect();
            let arbitrary: Vec<Q> = (0..rows)
                .map(|r| q((r as i64 * 3 + seed as i64) % 7 - 3))
                .collect();
            for b in [in_span, arbitrary] {
                let expected = reference_solution(&a, &b);
                assert_eq!(
                    solve_linear_system_with_rank(&a, &b),
                    expected,
                    "seed {seed}"
                );
                if expected.is_some() {
                    consistent += 1;
                } else {
                    inconsistent += 1;
                }
            }
        }
        assert!(
            consistent > 200 && inconsistent > 50,
            "{consistent} {inconsistent}"
        );
    }

    /// The kernel computation as it was before the integer elimination.
    fn reference_kernel(
        a: &[Vec<Q>],
        columns: usize,
        limits: RationalKernelLimits,
    ) -> Result<String, RationalKernelError> {
        if a.iter()
            .flatten()
            .any(|x| x.numer().bits().max(x.denom().bits()) > limits.max_coefficient_bits)
        {
            return Err(RationalKernelError::BudgetExceeded("coefficient bits"));
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
        rational_rref(&mut minor, rank, Some(&mut budget))?;
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
        Ok(format!(
            "{basis:?} {free:?} {pivots:?} {rows:?} {inverse:?}"
        ))
    }

    #[test]
    fn kernels_match_the_rational_reference_including_budgets() {
        for seed in 0..150u64 {
            let rows = 1 + (seed % 8) as usize;
            let columns = 1 + (seed / 3 % 7) as usize;
            let matrix = pseudo_random_matrix(9000 + seed, rows, columns, seed % 2 == 0);
            for (max_operations, max_coefficient_bits) in [
                (usize::MAX, u64::MAX),
                (usize::MAX, 6),
                (150, u64::MAX),
                (40, 12),
                (usize::MAX, 3),
            ] {
                let limits = RationalKernelLimits {
                    max_entries: usize::MAX / 4,
                    max_operations,
                    max_coefficient_bits,
                };
                let actual = rational_kernel_with_certificate(&matrix, columns, limits).map(|k| {
                    format!(
                        "{:?} {:?} {:?} {:?} {:?}",
                        k.basis,
                        k.free_columns,
                        k.pivot_columns,
                        k.independent_rows,
                        k.inverse_minor
                    )
                });
                assert_eq!(
                    actual,
                    reference_kernel(&matrix, columns, limits),
                    "seed {seed}, limits {max_operations} {max_coefficient_bits}"
                );
            }
        }
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
