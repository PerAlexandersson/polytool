//! Sparse Gaussian elimination over a prime field `F_p`.
//!
//! [`SparseModRow`] stores augmented rows sorted by column;
//! the `sparse_modular_linear_system_*` entry points decide consistency, report
//! rank statistics and optionally return one solution.  Row order and pivoting
//! are chosen with [`SparseModEliminationOptions`].

use super::modular::{
    add_mod_u64, is_prime_u64, mod_inverse_prime, mul_mod_u64, signed_mod_u64, sub_mod_u64,
};

/// Errors raised while constructing or updating a sparse modular row.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SparseModRowError {
    /// The modulus is zero.
    ZeroModulus,
}

impl std::fmt::Display for SparseModRowError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::ZeroModulus => write!(f, "sparse modular row modulus must be nonzero"),
        }
    }
}

impl std::error::Error for SparseModRowError {}

/// Sparse augmented row over the prime field `F_p`.
///
/// The row represents
///
/// ```text
///   sum_i entries[i].1 * x_{entries[i].0} = rhs   (mod p).
/// ```
///
/// Entries are stored sorted by column and normalized into `0..p`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SparseModRow {
    entries: Vec<(usize, u64)>,
    rhs: u64,
}

impl SparseModRow {
    /// Create an empty sparse row with right-hand side `rhs mod prime`.
    ///
    /// The `prime` argument is not stored; it is used only to normalize input.
    pub fn new(rhs: u64, prime: u64) -> Result<Self, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        Ok(Self {
            entries: Vec::new(),
            rhs: rhs % prime,
        })
    }

    pub(crate) fn with_capacity(rhs: u64, prime: u64, capacity: usize) -> Self {
        debug_assert_ne!(prime, 0);
        Self {
            entries: Vec::with_capacity(capacity),
            rhs: rhs % prime,
        }
    }

    /// Create a sparse row from unsigned `(column, value)` entries.
    ///
    /// Duplicate columns are combined modulo `prime`, and zero coefficients are
    /// removed.
    pub fn from_entries<I>(entries: I, rhs: u64, prime: u64) -> Result<Self, SparseModRowError>
    where
        I: IntoIterator<Item = (usize, u64)>,
    {
        let mut row = Self::new(rhs, prime)?;
        for (col, value) in entries {
            row.add_entry_unchecked(col, value, prime);
        }
        Ok(row)
    }

    /// Create an empty sparse row with a signed right-hand side.
    pub fn from_signed_rhs(rhs: i64, prime: u64) -> Result<Self, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        Ok(Self {
            entries: Vec::new(),
            rhs: signed_mod_u64(rhs, prime),
        })
    }

    /// Create a sparse row from signed `(column, value)` entries.
    ///
    /// Duplicate columns are combined modulo `prime`, and zero coefficients are
    /// removed.
    pub fn from_signed_entries<I>(
        entries: I,
        rhs: i64,
        prime: u64,
    ) -> Result<Self, SparseModRowError>
    where
        I: IntoIterator<Item = (usize, i64)>,
    {
        let mut row = Self::from_signed_rhs(rhs, prime)?;
        for (col, value) in entries {
            row.add_entry_unchecked(col, signed_mod_u64(value, prime), prime);
        }
        Ok(row)
    }

    /// Return the nonzero column/value entries, sorted by column.
    pub fn entries(&self) -> &[(usize, u64)] {
        &self.entries
    }

    /// Return the right-hand side.
    pub fn rhs(&self) -> u64 {
        self.rhs
    }

    /// Return true if `solution` satisfies this augmented row over `F_p`.
    pub fn is_satisfied_by(&self, solution: &[u64], prime: u64) -> Result<bool, SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        let mut lhs = 0;
        for &(col, coeff) in &self.entries {
            let Some(&value) = solution.get(col) else {
                return Ok(false);
            };
            lhs = add_mod_u64(lhs, mul_mod_u64(coeff, value, prime), prime);
        }
        Ok(lhs == self.rhs)
    }

    fn is_satisfied_by_validated_solution(&self, solution: &[u64], prime: u64) -> bool {
        let mut lhs = 0;
        for &(col, coeff) in &self.entries {
            debug_assert!(col < solution.len());
            lhs = add_mod_u64(lhs, mul_mod_u64(coeff, solution[col], prime), prime);
        }
        lhs == self.rhs
    }

    /// Add `value` to the coefficient of `col`, reducing modulo `prime`.
    pub fn add_entry(
        &mut self,
        col: usize,
        value: u64,
        prime: u64,
    ) -> Result<(), SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        self.add_entry_unchecked(col, value, prime);
        Ok(())
    }

    fn add_entry_unchecked(&mut self, col: usize, value: u64, prime: u64) {
        debug_assert_ne!(prime, 0);
        let value = value % prime;
        if value == 0 {
            return;
        }
        match self.entries.binary_search_by_key(&col, |&(col, _)| col) {
            Ok(index) => {
                let new_value = add_mod_u64(self.entries[index].1, value, prime);
                if new_value == 0 {
                    self.entries.remove(index);
                } else {
                    self.entries[index].1 = new_value;
                }
            }
            Err(index) => self.entries.insert(index, (col, value)),
        }
    }

    /// Add a signed value to the coefficient of `col`.
    pub fn add_signed_entry(
        &mut self,
        col: usize,
        value: i64,
        prime: u64,
    ) -> Result<(), SparseModRowError> {
        if prime == 0 {
            return Err(SparseModRowError::ZeroModulus);
        }
        self.add_entry_unchecked(col, signed_mod_u64(value, prime), prime);
        Ok(())
    }

    pub(crate) fn add_reduced_entry_sorted(&mut self, col: usize, value: u64, prime: u64) {
        debug_assert!(value < prime);
        if value == 0 {
            return;
        }
        let Some(&(last_col, _)) = self.entries.last() else {
            self.entries.push((col, value));
            return;
        };
        if col > last_col {
            self.entries.push((col, value));
        } else {
            self.add_entry_unchecked(col, value, prime);
        }
    }

    /// Return true if the row is the tautology `0 = 0`.
    pub fn is_tautology(&self) -> bool {
        self.entries.is_empty() && self.rhs == 0
    }

    /// Return true if the row is the contradiction `0 = c` with `c != 0`.
    pub fn is_contradiction(&self) -> bool {
        self.entries.is_empty() && self.rhs != 0
    }

    fn leading_entry(&self) -> Option<(usize, u64)> {
        self.entries.first().copied()
    }

    fn coefficient_at(&self, col: usize) -> Option<u64> {
        self.entries
            .binary_search_by_key(&col, |&(entry_col, _)| entry_col)
            .ok()
            .map(|index| self.entries[index].1)
    }

    fn normalize_pivot(&mut self, prime: u64) {
        let Some((_, pivot)) = self.leading_entry() else {
            return;
        };
        let pivot_inv = mod_inverse_prime(pivot, prime)
            .expect("nonzero element modulo prime should be invertible");
        for (_, value) in &mut self.entries {
            *value = mul_mod_u64(*value, pivot_inv, prime);
        }
        self.rhs = mul_mod_u64(self.rhs, pivot_inv, prime);
    }

    fn normalize_pivot_at(&mut self, pivot_col: usize, prime: u64) {
        let Some(pivot) = self.coefficient_at(pivot_col) else {
            return;
        };
        let pivot_inv = mod_inverse_prime(pivot, prime)
            .expect("nonzero element modulo prime should be invertible");
        for (_, value) in &mut self.entries {
            *value = mul_mod_u64(*value, pivot_inv, prime);
        }
        self.rhs = mul_mod_u64(self.rhs, pivot_inv, prime);
    }

    fn subtract_entry(&mut self, col: usize, value: u64, prime: u64) {
        let value = value % prime;
        if value == 0 {
            return;
        }
        match self.entries.binary_search_by_key(&col, |&(col, _)| col) {
            Ok(index) => {
                let new_value = sub_mod_u64(self.entries[index].1, value, prime);
                if new_value == 0 {
                    self.entries.remove(index);
                } else {
                    self.entries[index].1 = new_value;
                }
            }
            Err(index) => {
                let value = sub_mod_u64(0, value, prime);
                self.entries.insert(index, (col, value));
            }
        }
    }

    fn subtract_scaled_small_pivot(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        if self
            .entries
            .first()
            .is_some_and(|&(col, value)| col == pivot.entries[0].0 && value == factor)
        {
            self.entries.remove(0);
            for &(col, value) in pivot.entries.iter().skip(1) {
                self.subtract_entry(col, mul_mod_u64(factor, value, prime), prime);
            }
        } else {
            for &(col, value) in &pivot.entries {
                self.subtract_entry(col, mul_mod_u64(factor, value, prime), prime);
            }
        }
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled_leading_pivot(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let mut merged = Vec::with_capacity(
            self.entries
                .len()
                .saturating_add(pivot.entries.len().saturating_sub(2)),
        );
        let mut lhs = 1;
        let mut rhs = 1;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let factor = factor % prime;
        if factor == 0 {
            return;
        }
        if self
            .entries
            .first()
            .is_some_and(|&(col, value)| col == pivot.entries[0].0 && value == factor)
        {
            self.subtract_scaled_leading_pivot(pivot, factor, prime);
            return;
        }
        if pivot.entries.len() <= 4 {
            self.subtract_scaled_small_pivot(pivot, factor, prime);
            return;
        }

        let mut merged = Vec::with_capacity(self.entries.len() + pivot.entries.len());
        let mut lhs = 0;
        let mut rhs = 0;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }

    fn subtract_scaled_general(&mut self, pivot: &SparseModRow, factor: u64, prime: u64) {
        let factor = factor % prime;
        if factor == 0 {
            return;
        }

        let mut merged = Vec::with_capacity(self.entries.len() + pivot.entries.len());
        let mut lhs = 0;
        let mut rhs = 0;
        while lhs < self.entries.len() || rhs < pivot.entries.len() {
            let next_col = match (self.entries.get(lhs), pivot.entries.get(rhs)) {
                (Some(&(lhs_col, _)), Some(&(rhs_col, _))) => lhs_col.min(rhs_col),
                (Some(&(lhs_col, _)), None) => lhs_col,
                (None, Some(&(rhs_col, _))) => rhs_col,
                (None, None) => break,
            };

            let mut value = 0;
            if lhs < self.entries.len() && self.entries[lhs].0 == next_col {
                value = self.entries[lhs].1;
                lhs += 1;
            }
            if rhs < pivot.entries.len() && pivot.entries[rhs].0 == next_col {
                value = sub_mod_u64(
                    value,
                    mul_mod_u64(factor, pivot.entries[rhs].1, prime),
                    prime,
                );
                rhs += 1;
            }
            if value != 0 {
                merged.push((next_col, value));
            }
        }

        self.entries = merged;
        self.rhs = sub_mod_u64(self.rhs, mul_mod_u64(factor, pivot.rhs, prime), prime);
    }
}

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
pub fn sparse_modular_linear_system_consistent<I>(
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

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------------
    // Sparse modular linear systems
    // -----------------------------------------------------------------------

    fn dense_modular_linear_system_consistent(
        mut aug: Vec<Vec<u64>>,
        num_vars: usize,
        prime: u64,
    ) -> bool {
        let num_rows = aug.len();
        let mut pivot_row = 0;

        for col in 0..num_vars {
            let Some(pr) = (pivot_row..num_rows).find(|&row| aug[row][col] != 0) else {
                continue;
            };
            aug.swap(pivot_row, pr);
            let pivot_inv = mod_inverse_prime(aug[pivot_row][col], prime)
                .expect("nonzero element modulo prime is invertible");
            let pivot_snapshot = aug[pivot_row][col..=num_vars].to_vec();

            for aug_row in aug.iter_mut().take(num_rows).skip(pivot_row + 1) {
                if aug_row[col] == 0 {
                    continue;
                }
                let factor = mul_mod_u64(aug_row[col], pivot_inv, prime);
                for (entry, pivot_entry) in aug_row[col..=num_vars]
                    .iter_mut()
                    .zip(pivot_snapshot.iter())
                {
                    *entry = sub_mod_u64(*entry, mul_mod_u64(factor, *pivot_entry, prime), prime);
                }
            }

            pivot_row += 1;
            if pivot_row == num_rows {
                return true;
            }
        }

        aug[pivot_row..]
            .iter()
            .all(|row| row[..num_vars].iter().any(|&entry| entry != 0) || row[num_vars] == 0)
    }

    fn sparse_rows_from_dense_aug(
        aug: &[Vec<u64>],
        num_vars: usize,
        prime: u64,
    ) -> Vec<SparseModRow> {
        aug.iter()
            .map(|dense_row| {
                let mut row = SparseModRow::new(dense_row[num_vars], prime).unwrap();
                for (col, &value) in dense_row.iter().take(num_vars).enumerate() {
                    row.add_entry(col, value, prime).unwrap();
                }
                row
            })
            .filter(|row| !row.is_tautology())
            .collect()
    }

    fn solution_satisfies_dense_aug(
        aug: &[Vec<u64>],
        num_vars: usize,
        prime: u64,
        solution: &[u64],
    ) -> bool {
        aug.iter().all(|row| {
            let lhs = row
                .iter()
                .take(num_vars)
                .zip(solution.iter())
                .fold(0, |acc, (&coeff, &value)| {
                    add_mod_u64(acc, mul_mod_u64(coeff, value, prime), prime)
                });
            lhs == row[num_vars] % prime
        })
    }

    #[test]
    fn test_sparse_mod_row_normalizes_and_cancels_entries() {
        let prime = 101;
        let mut row = SparseModRow::from_signed_rhs(-3, prime).unwrap();
        row.add_signed_entry(2, -2, prime).unwrap();
        row.add_entry(2, 103, prime).unwrap();
        row.add_signed_entry(1, -1, prime).unwrap();

        assert_eq!(row.rhs(), 98);
        assert_eq!(row.entries(), &[(1, 100)]);

        let row = SparseModRow::from_entries([(2, 5), (1, 3), (2, 96)], 204, prime).unwrap();
        assert_eq!(row.rhs(), 2);
        assert_eq!(row.entries(), &[(1, 3)]);
        assert!(!row.is_tautology());
        assert!(!row.is_contradiction());

        let row = SparseModRow::new(7, prime).unwrap();
        assert!(row.is_contradiction());
    }

    #[test]
    fn test_sparse_mod_sorted_add_fast_path_and_fallback() {
        let prime = 101;
        let mut row = SparseModRow::new(0, prime).unwrap();

        row.add_reduced_entry_sorted(1, 4, prime);
        row.add_reduced_entry_sorted(3, 8, prime);
        row.add_reduced_entry_sorted(2, 5, prime);
        row.add_reduced_entry_sorted(3, 93, prime);

        assert_eq!(row.entries(), &[(1, 4), (2, 5)]);
    }

    #[test]
    fn test_sparse_mod_small_pivot_reduction_cancels_in_place() {
        let prime = 101;
        let pivot = SparseModRow::from_entries([(0, 1), (2, 4)], 7, prime).unwrap();
        let mut row = SparseModRow::from_entries([(0, 5), (1, 3), (2, 20)], 9, prime).unwrap();

        row.subtract_scaled(&pivot, 5, prime);

        assert_eq!(row.rhs(), 75);
        assert_eq!(row.entries(), &[(1, 3)]);
    }

    #[test]
    fn test_sparse_modular_solver_matches_dense_consistency() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![
            vec![1, 2, 0, 5],
            vec![0, 3, 4, 7],
            vec![2, 7, 4, 17],
            vec![0, 0, 0, 0],
        ];
        let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
        let result =
            sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

        assert_eq!(
            result.consistent,
            dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
        );
        assert_eq!(result.rank, 2);
        assert_eq!(result.pivot_rows, vec![0, 1]);
        assert_eq!(result.pivot_columns, vec![0, 1]);
        assert_eq!(result.inconsistent_row, None);
        let solution = result.solution.as_ref().expect("consistent system");
        assert_eq!(solution.len(), num_vars);
        assert!(solution_satisfies_dense_aug(
            &aug, num_vars, prime, solution
        ));

        let direct_solution = sparse_modular_linear_system_solution(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
        )
        .unwrap()
        .expect("consistent system");
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            &direct_solution
        ));
    }

    #[test]
    fn test_sparse_modular_solver_can_process_sparsest_rows_first() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![vec![1, 1, 1, 1], vec![0, 1, 0, 2], vec![0, 0, 1, 3]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingNonzeros,
                compute_solution: true,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.rank, 3);
        assert_eq!(result.pivot_rows, vec![1, 2, 0]);
        assert_eq!(result.pivot_columns, vec![1, 2, 0]);
        assert_eq!(result.stats.input_rows, 3);
        assert_eq!(result.stats.input_nonzeros, 5);
        assert_eq!(result.stats.row_reductions, 0);
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            result.solution.as_ref().expect("consistent system")
        ));

        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingNonzeros,
                compute_solution: false,
            },
        )
        .unwrap();
        assert!(result.consistent);
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_checks_tail_rows_after_full_rank() {
        let prime = 101;
        let num_vars = 2;
        let aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 5], vec![2, 3, 13]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::Input,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.pivot_rows, vec![0, 1]);
        assert_eq!(result.stats.row_reductions, 0);
        assert_eq!(result.stats.full_rank_checks, 2);
        assert_eq!(result.solution, None);

        let inconsistent_aug = vec![vec![1, 0, 2], vec![0, 1, 3], vec![1, 1, 6]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&inconsistent_aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::Input,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.inconsistent_row, Some(2));
        assert_eq!(result.stats.full_rank_checks, 1);
    }

    #[test]
    fn test_sparse_modular_solver_full_rank_solution_mode() {
        let prime = 101;
        let num_vars = 3;
        let rank_deficient_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3]];
        let full_rank_only =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::FullRank,
            )
            .unwrap();
        assert!(full_rank_only.consistent);
        assert_eq!(full_rank_only.rank, 2);
        assert_eq!(full_rank_only.solution, None);

        let any_consistent =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&rank_deficient_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::AnyConsistent,
            )
            .unwrap();
        assert_eq!(any_consistent.rank, 2);
        assert!(solution_satisfies_dense_aug(
            &rank_deficient_aug,
            num_vars,
            prime,
            any_consistent
                .solution
                .as_ref()
                .expect("rank-deficient solution")
        ));

        let full_rank_aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![0, 0, 1, 4]];
        let full_rank =
            sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode(
                sparse_rows_from_dense_aug(&full_rank_aug, num_vars, prime),
                num_vars,
                prime,
                SparseModRowOrder::IncreasingNonzeros,
                SparseModSolutionMode::FullRank,
            )
            .unwrap();
        assert_eq!(full_rank.rank, 3);
        assert!(solution_satisfies_dense_aug(
            &full_rank_aug,
            num_vars,
            prime,
            full_rank.solution.as_ref().expect("full-rank solution")
        ));
    }

    #[test]
    fn test_sparse_modular_solver_can_use_initial_markowitz_cost() {
        let prime = 101;
        let num_vars = 5;
        let rows = vec![
            SparseModRow::from_entries([(0, 1), (3, 1)], 1, prime).unwrap(),
            SparseModRow::from_entries([(0, 1), (4, 1)], 2, prime).unwrap(),
            SparseModRow::from_entries([(1, 1), (3, 1)], 3, prime).unwrap(),
            SparseModRow::from_entries([(1, 1), (4, 1)], 4, prime).unwrap(),
            SparseModRow::from_entries([(2, 1), (3, 1)], 5, prime).unwrap(),
        ];
        let result = sparse_modular_linear_system_consistency_with_options(
            rows,
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::IncreasingMarkowitzCost,
                compute_solution: false,
            },
        )
        .unwrap();

        assert!(result.consistent);
        assert_eq!(result.pivot_rows[0], 4);
        assert_eq!(result.pivot_columns[0], 2);
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_dynamic_markowitz_matches_dense() {
        let prime = 101;
        let num_vars = 4;
        let aug = vec![
            vec![1, 1, 1, 1, 10],
            vec![0, 0, 1, 1, 7],
            vec![0, 1, 0, 1, 8],
            vec![1, 0, 0, 0, 2],
            vec![2, 1, 1, 2, 19],
            vec![0, 0, 0, 0, 0],
        ];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::DynamicMarkowitz,
                compute_solution: true,
            },
        )
        .unwrap();

        assert_eq!(
            result.consistent,
            dense_modular_linear_system_consistent(aug.clone(), num_vars, prime)
        );
        assert_eq!(result.rank, 4);
        assert_eq!(result.solution.as_ref().map(Vec::len), Some(num_vars));
        assert!(solution_satisfies_dense_aug(
            &aug,
            num_vars,
            prime,
            result.solution.as_ref().expect("consistent system")
        ));
        assert_eq!(result.stats.zero_rows, 0);
        assert!(result.stats.row_reductions > 0);
    }

    #[test]
    fn test_sparse_modular_solver_dynamic_markowitz_detects_inconsistency() {
        let prime = 101;
        let num_vars = 3;
        let aug = vec![vec![1, 0, 0, 2], vec![0, 1, 0, 3], vec![1, 1, 0, 6]];
        let result = sparse_modular_linear_system_consistency_with_options(
            sparse_rows_from_dense_aug(&aug, num_vars, prime),
            num_vars,
            prime,
            SparseModEliminationOptions {
                row_order: SparseModRowOrder::DynamicMarkowitz,
                compute_solution: true,
            },
        )
        .unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 2);
        assert_eq!(result.inconsistent_row, Some(2));
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_detects_inconsistency() {
        let prime = 101;
        let num_vars = 2;
        let aug = vec![vec![1, 1, 1], vec![2, 2, 3]];
        let sparse_rows = sparse_rows_from_dense_aug(&aug, num_vars, prime);
        let result =
            sparse_modular_linear_system_consistency(sparse_rows, num_vars, prime).unwrap();

        assert!(!result.consistent);
        assert_eq!(result.rank, 1);
        assert_eq!(result.pivot_rows, vec![0]);
        assert_eq!(result.pivot_columns, vec![0]);
        assert_eq!(result.inconsistent_row, Some(1));
        assert_eq!(result.solution, None);
    }

    #[test]
    fn test_sparse_modular_solver_rejects_bad_inputs() {
        let row = {
            let mut row = SparseModRow::new(0, 101).unwrap();
            row.add_entry(3, 1, 101).unwrap();
            row
        };
        assert_eq!(
            sparse_modular_linear_system_consistency(vec![row], 3, 101),
            Err(SparseModEliminationError::ColumnOutOfRange {
                row: 0,
                column: 3,
                num_vars: 3
            })
        );
        assert_eq!(
            sparse_modular_linear_system_consistency(Vec::new(), 3, 100),
            Err(SparseModEliminationError::ModulusNotPrime { modulus: 100 })
        );
        assert_eq!(SparseModRow::new(1, 0), Err(SparseModRowError::ZeroModulus));
        assert_eq!(
            SparseModRow::from_signed_rhs(-1, 0),
            Err(SparseModRowError::ZeroModulus)
        );
        let mut row = SparseModRow::new(0, 101).unwrap();
        assert_eq!(row.add_entry(0, 1, 0), Err(SparseModRowError::ZeroModulus));
    }
}
