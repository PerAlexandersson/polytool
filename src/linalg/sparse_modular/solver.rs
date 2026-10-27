use super::super::modular::{is_prime_u64, mul_mod_u64, sub_mod_u64};
use super::SparseModRow;

/// Row ordering strategy for sparse modular elimination.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum SparseModRowOrder {
    /// Process rows in input order.  This streams rows without collecting them.
    #[default]
    Input,
    /// Process rows by increasing number of nonzero coefficients.
    ///
    /// This is a cheap fill-in heuristic.  The result still reports original
    /// zero-based input row indices in `pivot_rows` and `inconsistent_row`.
    IncreasingNonzeros,
    /// Process rows by an initial Markowitz-style fill-in score.
    ///
    /// The score is `(row_nnz - 1) * (leading_column_nnz - 1)`, where
    /// `leading_column_nnz` is computed from the unreduced input system.  This
    /// is cheaper than dynamic Markowitz pivoting but tends to prefer sparse
    /// rows whose leading pivot column is also sparse.
    IncreasingMarkowitzCost,
    /// Repeatedly choose a pivot entry with low current Markowitz fill score.
    ///
    /// This collects the matrix and updates all active rows after each pivot,
    /// so it has higher overhead than the streaming backends. It can preserve
    /// sparsity much better on systems where the leading-column order produces
    /// unnecessary fill-in.
    DynamicMarkowitz,
}

/// Options for sparse modular elimination.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct SparseModEliminationOptions {
    /// Row ordering strategy.
    pub row_order: SparseModRowOrder,
    /// Whether to compute one solution when the system is consistent.
    pub compute_solution: bool,
}

impl Default for SparseModEliminationOptions {
    fn default() -> Self {
        Self {
            row_order: SparseModRowOrder::Input,
            compute_solution: true,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SparseModSolutionMode {
    None,
    AnyConsistent,
    FullRank,
}

impl SparseModSolutionMode {
    fn from_compute_solution(compute_solution: bool) -> Self {
        if compute_solution {
            Self::AnyConsistent
        } else {
            Self::None
        }
    }
}

/// Lightweight statistics from sparse modular elimination.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SparseModEliminationStats {
    /// Number of input rows seen by the solver.
    pub input_rows: usize,
    /// Number of nonzero coefficients in the input rows.
    pub input_nonzeros: usize,
    /// Number of tautological rows reduced to `0 = 0`.
    pub zero_rows: usize,
    /// Number of sparse row reductions against existing pivot rows.
    pub row_reductions: usize,
    /// Number of rows checked by evaluating the unique solution after full
    /// column rank was reached.
    pub full_rank_checks: usize,
    /// Sum of nonzero coefficients in stored pivot rows.
    pub pivot_nonzeros: usize,
    /// Maximum number of nonzero coefficients in any active row during
    /// elimination.
    pub max_active_nonzeros: usize,
}

/// Summary of sparse Gaussian elimination over a prime field.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseModEliminationResult {
    /// Whether the input system is consistent.
    pub consistent: bool,
    /// Number of pivot rows found before termination.
    pub rank: usize,
    /// Zero-based input row indices that introduced new pivots.
    pub pivot_rows: Vec<usize>,
    /// Pivot columns in the order they were inserted.
    pub pivot_columns: Vec<usize>,
    /// Zero-based input row index witnessing inconsistency, if any.
    pub inconsistent_row: Option<usize>,
    /// One solution with all free variables set to zero, if the system is
    /// consistent and solution computation was requested.
    pub solution: Option<Vec<u64>>,
    /// Sparse elimination statistics useful for tuning row-order and fill-in
    /// heuristics.
    pub stats: SparseModEliminationStats,
}

/// Errors for sparse modular linear-system elimination.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SparseModEliminationError {
    /// The modulus is not prime.  The elimination routine divides by nonzero
    /// pivots, so it must run over a field.
    ModulusNotPrime {
        /// The offending modulus.
        modulus: u64,
    },
    /// A row contains a column outside `0..num_vars`.
    ColumnOutOfRange {
        /// Caller-supplied index of the offending row.
        row: usize,
        /// The out-of-range column.
        column: usize,
        /// Number of variables; valid columns are `0..num_vars`.
        num_vars: usize,
    },
}

impl std::fmt::Display for SparseModEliminationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ModulusNotPrime { modulus } => {
                write!(f, "modulus {modulus} is not prime")
            }
            Self::ColumnOutOfRange {
                row,
                column,
                num_vars,
            } => write!(
                f,
                "row {row} contains column {column}, but the system has {num_vars} variables"
            ),
        }
    }
}

impl std::error::Error for SparseModEliminationError {}

fn sparse_modular_solution_from_pivots(
    pivots: &[Option<SparseModRow>],
    num_vars: usize,
    prime: u64,
) -> Vec<u64> {
    let mut solution = vec![0; num_vars];
    for pivot_col in (0..num_vars).rev() {
        let Some(row) = &pivots[pivot_col] else {
            continue;
        };
        let mut value = row.rhs;
        for &(col, coeff) in row.entries.iter().skip(1) {
            value = sub_mod_u64(value, mul_mod_u64(coeff, solution[col], prime), prime);
        }
        solution[pivot_col] = value;
    }
    solution
}

fn sparse_modular_solution_from_ordered_pivots(
    pivot_rows: &[SparseModRow],
    pivot_columns: &[usize],
    num_vars: usize,
    prime: u64,
) -> Vec<u64> {
    let mut solution = vec![0; num_vars];
    for (row, &pivot_col) in pivot_rows.iter().zip(pivot_columns.iter()).rev() {
        let mut value = row.rhs;
        for &(col, coeff) in &row.entries {
            if col != pivot_col {
                value = sub_mod_u64(value, mul_mod_u64(coeff, solution[col], prime), prime);
            }
        }
        solution[pivot_col] = value;
    }
    solution
}

fn sparse_modular_row_matches_solution(row: &SparseModRow, solution: &[u64], prime: u64) -> bool {
    row.is_satisfied_by_validated_solution(solution, prime)
}

fn sparse_modular_column_counts(
    indexed_rows: &[(usize, SparseModRow)],
    num_vars: usize,
) -> Vec<usize> {
    let mut column_counts = vec![0; num_vars];
    for (_, row) in indexed_rows {
        for &(col, _) in &row.entries {
            if col < num_vars {
                column_counts[col] += 1;
            }
        }
    }
    column_counts
}

fn sparse_modular_initial_markowitz_cost(
    row: &SparseModRow,
    column_counts: &[usize],
) -> (usize, usize, usize) {
    let Some(&(col, _)) = row.entries.first() else {
        return (0, 0, 0);
    };
    let row_fill = row.entries.len().saturating_sub(1);
    let column_fill = column_counts
        .get(col)
        .copied()
        .unwrap_or(usize::MAX)
        .saturating_sub(1);
    (
        row_fill.saturating_mul(column_fill),
        row.entries.len(),
        column_fill,
    )
}

fn sparse_modular_dynamic_markowitz_pivot(
    active_rows: &[(usize, SparseModRow)],
    num_vars: usize,
) -> Option<(usize, usize)> {
    let column_counts = sparse_modular_column_counts(active_rows, num_vars);
    let mut best: Option<(usize, usize, usize, usize, usize, usize)> = None;

    for (row_pos, (row_index, row)) in active_rows.iter().enumerate() {
        let row_fill = row.entries.len().saturating_sub(1);
        for &(col, _) in &row.entries {
            let column_fill = column_counts[col].saturating_sub(1);
            let score = row_fill.saturating_mul(column_fill);
            let candidate = (
                score,
                row.entries.len(),
                column_counts[col],
                *row_index,
                col,
                row_pos,
            );
            if best.is_none_or(|best| candidate < best) {
                best = Some(candidate);
            }
        }
    }

    best.map(|(_, _, _, _, col, row_pos)| (row_pos, col))
}

fn sparse_modular_linear_system_consistency_indexed<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = (usize, SparseModRow)>,
{
    let mut pivots = vec![None; num_vars];
    let mut pivot_rows = Vec::new();
    let mut pivot_columns = Vec::new();
    let mut stats = SparseModEliminationStats::default();
    let mut full_rank_solution: Option<Vec<u64>> = None;

    for (row_index, mut row) in rows {
        stats.input_rows += 1;
        stats.input_nonzeros += row.entries.len();
        stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        for &(column, _) in &row.entries {
            if column >= num_vars {
                return Err(SparseModEliminationError::ColumnOutOfRange {
                    row: row_index,
                    column,
                    num_vars,
                });
            }
        }

        if let Some(solution) = &full_rank_solution {
            stats.full_rank_checks += 1;
            if !sparse_modular_row_matches_solution(&row, solution, prime) {
                return Ok(SparseModEliminationResult {
                    consistent: false,
                    rank: pivot_columns.len(),
                    pivot_rows,
                    pivot_columns,
                    inconsistent_row: Some(row_index),
                    solution: None,
                    stats,
                });
            }
            continue;
        }

        loop {
            let Some((pivot_col, pivot_coeff)) = row.leading_entry() else {
                if row.rhs != 0 {
                    return Ok(SparseModEliminationResult {
                        consistent: false,
                        rank: pivot_columns.len(),
                        pivot_rows,
                        pivot_columns,
                        inconsistent_row: Some(row_index),
                        solution: None,
                        stats,
                    });
                }
                stats.zero_rows += 1;
                break;
            };

            let Some(pivot) = &pivots[pivot_col] else {
                row.normalize_pivot(prime);
                stats.pivot_nonzeros += row.entries.len();
                pivots[pivot_col] = Some(row);
                pivot_rows.push(row_index);
                pivot_columns.push(pivot_col);
                if pivot_columns.len() == num_vars {
                    full_rank_solution = Some(sparse_modular_solution_from_pivots(
                        &pivots, num_vars, prime,
                    ));
                }
                break;
            };
            row.subtract_scaled(pivot, pivot_coeff, prime);
            stats.row_reductions += 1;
            stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        }
    }
    let solution = match solution_mode {
        SparseModSolutionMode::None => None,
        SparseModSolutionMode::AnyConsistent => full_rank_solution.or_else(|| {
            Some(sparse_modular_solution_from_pivots(
                &pivots, num_vars, prime,
            ))
        }),
        SparseModSolutionMode::FullRank => full_rank_solution,
    };

    Ok(SparseModEliminationResult {
        consistent: true,
        rank: pivot_columns.len(),
        pivot_rows,
        pivot_columns,
        inconsistent_row: None,
        solution,
        stats,
    })
}

fn sparse_modular_linear_system_consistency_dynamic_markowitz(
    rows: Vec<(usize, SparseModRow)>,
    num_vars: usize,
    prime: u64,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError> {
    let mut active_rows = Vec::with_capacity(rows.len());
    let mut pivot_rows = Vec::new();
    let mut pivot_columns = Vec::new();
    let mut ordered_pivots = Vec::new();
    let mut stats = SparseModEliminationStats::default();

    for (row_index, row) in rows {
        stats.input_rows += 1;
        stats.input_nonzeros += row.entries.len();
        stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        for &(column, _) in &row.entries {
            if column >= num_vars {
                return Err(SparseModEliminationError::ColumnOutOfRange {
                    row: row_index,
                    column,
                    num_vars,
                });
            }
        }
        if row.is_contradiction() {
            return Ok(SparseModEliminationResult {
                consistent: false,
                rank: 0,
                pivot_rows,
                pivot_columns,
                inconsistent_row: Some(row_index),
                solution: None,
                stats,
            });
        }
        if row.is_tautology() {
            stats.zero_rows += 1;
        } else {
            active_rows.push((row_index, row));
        }
    }

    loop {
        let mut row_pos = 0;
        while row_pos < active_rows.len() {
            if active_rows[row_pos].1.is_contradiction() {
                return Ok(SparseModEliminationResult {
                    consistent: false,
                    rank: pivot_columns.len(),
                    pivot_rows,
                    pivot_columns,
                    inconsistent_row: Some(active_rows[row_pos].0),
                    solution: None,
                    stats,
                });
            }
            if active_rows[row_pos].1.is_tautology() {
                stats.zero_rows += 1;
                active_rows.swap_remove(row_pos);
            } else {
                row_pos += 1;
            }
        }

        let Some((pivot_pos, pivot_col)) =
            sparse_modular_dynamic_markowitz_pivot(&active_rows, num_vars)
        else {
            break;
        };
        let (pivot_row_index, mut pivot) = active_rows.swap_remove(pivot_pos);
        pivot.normalize_pivot_at(pivot_col, prime);
        stats.pivot_nonzeros += pivot.entries.len();
        pivot_rows.push(pivot_row_index);
        pivot_columns.push(pivot_col);
        ordered_pivots.push(pivot);

        if pivot_columns.len() == num_vars {
            let solution = sparse_modular_solution_from_ordered_pivots(
                &ordered_pivots,
                &pivot_columns,
                num_vars,
                prime,
            );
            stats.full_rank_checks += active_rows.len();
            for (row_index, row) in &active_rows {
                if !sparse_modular_row_matches_solution(row, &solution, prime) {
                    return Ok(SparseModEliminationResult {
                        consistent: false,
                        rank: pivot_columns.len(),
                        pivot_rows,
                        pivot_columns,
                        inconsistent_row: Some(*row_index),
                        solution: None,
                        stats,
                    });
                }
            }
            return Ok(SparseModEliminationResult {
                consistent: true,
                rank: pivot_columns.len(),
                pivot_rows,
                pivot_columns,
                inconsistent_row: None,
                solution: match solution_mode {
                    SparseModSolutionMode::None => None,
                    SparseModSolutionMode::AnyConsistent | SparseModSolutionMode::FullRank => {
                        Some(solution)
                    }
                },
                stats,
            });
        }

        let pivot = ordered_pivots
            .last()
            .expect("pivot was just pushed before active-row reduction");
        for (_, row) in &mut active_rows {
            let Some(factor) = row.coefficient_at(pivot_col) else {
                continue;
            };
            row.subtract_scaled_general(pivot, factor, prime);
            stats.row_reductions += 1;
            stats.max_active_nonzeros = stats.max_active_nonzeros.max(row.entries.len());
        }
    }

    let solution = match solution_mode {
        SparseModSolutionMode::None | SparseModSolutionMode::FullRank => None,
        SparseModSolutionMode::AnyConsistent => Some(sparse_modular_solution_from_ordered_pivots(
            &ordered_pivots,
            &pivot_columns,
            num_vars,
            prime,
        )),
    };

    Ok(SparseModEliminationResult {
        consistent: true,
        rank: pivot_columns.len(),
        pivot_rows,
        pivot_columns,
        inconsistent_row: None,
        solution,
        stats,
    })
}

/// Check consistency of a sparse linear system over the prime field `F_p`.
///
/// The rows are augmented rows `A_i x = b_i`.  The implementation performs
/// incremental sparse Gaussian elimination, normalizing each stored pivot row.
/// It is intended for sparse systems where one mostly needs rank/consistency
/// information or a single solution over `F_p` before a more expensive exact
/// solve.
pub fn sparse_modular_linear_system_consistency<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_with_options(
        rows,
        num_vars,
        prime,
        SparseModEliminationOptions::default(),
    )
}

/// Check consistency of a sparse linear system over `F_p` with explicit
/// elimination options.
pub fn sparse_modular_linear_system_consistency_with_options<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    options: SparseModEliminationOptions,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    if !is_prime_u64(prime) {
        return Err(SparseModEliminationError::ModulusNotPrime { modulus: prime });
    }

    sparse_modular_linear_system_consistency_prime_unchecked_with_options(
        rows, num_vars, prime, options,
    )
}

pub(crate) fn sparse_modular_linear_system_consistency_prime_unchecked_with_options<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    options: SparseModEliminationOptions,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
        rows,
        num_vars,
        prime,
        options.row_order,
        SparseModSolutionMode::from_compute_solution(options.compute_solution),
    )
}

pub(crate) fn sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
    row_order: SparseModRowOrder,
    solution_mode: SparseModSolutionMode,
) -> Result<SparseModEliminationResult, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    match row_order {
        SparseModRowOrder::Input => sparse_modular_linear_system_consistency_indexed(
            rows.into_iter().enumerate(),
            num_vars,
            prime,
            solution_mode,
        ),
        SparseModRowOrder::IncreasingNonzeros => {
            let mut indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            indexed_rows
                .sort_unstable_by_key(|(index, row)| (row.entries.len(), row.rhs == 0, *index));
            sparse_modular_linear_system_consistency_indexed(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
        SparseModRowOrder::IncreasingMarkowitzCost => {
            let mut indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            let column_counts = sparse_modular_column_counts(&indexed_rows, num_vars);
            indexed_rows.sort_unstable_by_key(|(index, row)| {
                (
                    sparse_modular_initial_markowitz_cost(row, &column_counts),
                    row.rhs == 0,
                    *index,
                )
            });
            sparse_modular_linear_system_consistency_indexed(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
        SparseModRowOrder::DynamicMarkowitz => {
            let indexed_rows = rows.into_iter().enumerate().collect::<Vec<_>>();
            sparse_modular_linear_system_consistency_dynamic_markowitz(
                indexed_rows,
                num_vars,
                prime,
                solution_mode,
            )
        }
    }
}

/// Convenience boolean wrapper around
/// [`sparse_modular_linear_system_consistency`].
pub fn is_sparse_modular_linear_system_consistent<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<bool, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency_with_options(
        rows,
        num_vars,
        prime,
        SparseModEliminationOptions {
            compute_solution: false,
            ..SparseModEliminationOptions::default()
        },
    )
    .map(|result| result.consistent)
}

/// Deprecated name of [`is_sparse_modular_linear_system_consistent`].
#[deprecated(
    since = "0.6.0",
    note = "renamed to `is_sparse_modular_linear_system_consistent`"
)]
pub fn sparse_modular_linear_system_consistent<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<bool, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    is_sparse_modular_linear_system_consistent(rows, num_vars, prime)
}

/// Return one solution of a sparse linear system over `F_p`, if the system is
/// consistent.
///
/// Free variables are set to zero.  Use
/// [`sparse_modular_linear_system_consistency`] when rank, pivot row indices,
/// pivot columns, or an inconsistency witness are also needed.
pub fn sparse_modular_linear_system_solution<I>(
    rows: I,
    num_vars: usize,
    prime: u64,
) -> Result<Option<Vec<u64>>, SparseModEliminationError>
where
    I: IntoIterator<Item = SparseModRow>,
{
    sparse_modular_linear_system_consistency(rows, num_vars, prime).map(|result| result.solution)
}
