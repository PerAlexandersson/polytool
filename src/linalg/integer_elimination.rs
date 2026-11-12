//! Fraction-free Gauss-Jordan elimination of a rational matrix on integer rows.
//!
//! The reduced row-echelon computation of [`super::linear_systems`] used to run on
//! `Ratio<BigInt>` entries, where every product and difference ends in a gcd
//! reduction even though the entries of a recurrence-fit system are integers
//! (or become integers after one scaling of each row).  This module keeps every row
//! as an integer vector that is a nonzero rational multiple of the true row, scales
//! it once at the start (least common multiple of the denominators), and eliminates
//! with the Bareiss/Montante update
//!
//! ```text
//! row_r  <-  (v * row_r - t * row_pivot) / prev
//! ```
//!
//! where `v` is the pivot, `t = row_r[col]` and `prev` the previous pivot; the
//! division is exact (the entries are minors of the input).  The pivot columns, the
//! pivot-row order and the zero pattern of every entry are those of the rational
//! elimination, because each integer row stays a nonzero multiple of the rational
//! row it stands for; only the scale of a row differs.  Every output of the callers
//! is a ratio of two entries of one row (a solution coordinate, a kernel
//! coordinate, an inverse-minor entry), so the scale cancels.
//!
//! Rows that are not touched by a step (their entry in the pivot column is zero) are
//! not rescaled; each row records the pivot (its *level*) at which its integers are
//! valid. A pivot row is brought to the current level when it is next used.
//! For a target row at level `L`, we combine that rescaling with elimination:
//! `(v * (prev/L) * row - (prev/L) * t * pivot_row) / prev` simplifies to
//! `(v * row - t * pivot_row) / L`. This avoids separately multiplying and
//! dividing every entry of a stale target row, with the same integer result.
//!
//! With a [`EliminationBudget`], the true rational entries are reconstructed lazily:
//! every row carries the exact rational scalar `c` with `true row = c * integer
//! row`, and a value is reduced to lowest terms only when the cheap size bound
//! (bits of the unreduced numerator and denominator) exceeds the limit, so the
//! accounting is the same as that of the rational elimination.

use super::linear_systems::RationalKernelError;
use super::Q;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

// Test hook: while a `ForcedFallback` guard lives on this thread, every integer
// elimination reports an unfinished (inexact) division, so the rational fallbacks
// run on inputs where they would otherwise never be reached.
#[cfg(test)]
thread_local! {
    static FORCE_FALLBACK: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

#[cfg(test)]
pub(super) struct ForcedFallback;

#[cfg(test)]
impl ForcedFallback {
    pub(super) fn engage() -> Self {
        FORCE_FALLBACK.with(|flag| flag.set(true));
        Self
    }
}

#[cfg(test)]
impl Drop for ForcedFallback {
    fn drop(&mut self) {
        FORCE_FALLBACK.with(|flag| flag.set(false));
    }
}

/// Whether the integer elimination must pretend that an exact division failed.
pub(super) fn fallback_forced() -> bool {
    #[cfg(test)]
    {
        FORCE_FALLBACK.with(std::cell::Cell::get)
    }
    #[cfg(not(test))]
    {
        false
    }
}

/// Operation and size limits of one elimination.
pub(super) struct EliminationBudget {
    pub(super) remaining: usize,
    pub(super) bits: u64,
}

impl EliminationBudget {
    /// Count one arithmetic operation.
    fn spend(&mut self) -> Result<(), RationalKernelError> {
        self.remaining = self
            .remaining
            .checked_sub(1)
            .ok_or(RationalKernelError::BudgetExceeded("arithmetic operations"))?;
        Ok(())
    }

    /// Count one operation and check the size of a computed rational.
    pub(super) fn check(&mut self, value: &Q) -> Result<(), RationalKernelError> {
        self.spend()?;
        if value.numer().bits().max(value.denom().bits()) > self.bits {
            return Err(RationalKernelError::BudgetExceeded("coefficient bits"));
        }
        Ok(())
    }

    /// Count one operation and check the size of `entry * scale` in lowest terms
    /// without reducing it unless the unreduced bound exceeds the limit.
    fn check_scaled_entry(&mut self, entry: &BigInt, scale: &Q) -> Result<(), RationalKernelError> {
        self.spend()?;
        if entry.is_zero() {
            return Ok(());
        }
        let numerator_bound = entry.bits() + scale.numer().bits();
        if numerator_bound.max(scale.denom().bits()) <= self.bits {
            return Ok(());
        }
        let value = Q::new(entry * scale.numer(), scale.denom().clone());
        if value.numer().bits().max(value.denom().bits()) > self.bits {
            return Err(RationalKernelError::BudgetExceeded("coefficient bits"));
        }
        Ok(())
    }

    /// Count one operation and check the elimination factor
    /// `(t * c_target) / (v * c_pivot)`.
    fn check_factor(
        &mut self,
        t: &BigInt,
        target_scale: &Q,
        v: &BigInt,
        pivot_scale: &Q,
    ) -> Result<(), RationalKernelError> {
        self.spend()?;
        let numerator_bound = t.bits() + target_scale.numer().bits() + pivot_scale.denom().bits();
        let denominator_bound = target_scale.denom().bits() + v.bits() + pivot_scale.numer().bits();
        if numerator_bound.max(denominator_bound) <= self.bits {
            return Ok(());
        }
        let value = Q::new(
            t * target_scale.numer() * pivot_scale.denom(),
            target_scale.denom() * v * pivot_scale.numer(),
        );
        if value.numer().bits().max(value.denom().bits()) > self.bits {
            return Err(RationalKernelError::BudgetExceeded("coefficient bits"));
        }
        Ok(())
    }
}

/// The Gauss-Jordan form of a matrix: integer rows, each a nonzero rational
/// multiple of the corresponding reduced row-echelon row of the rational
/// elimination (without normalising the pivots).
pub(super) struct IntegerRref {
    /// All rows, the first `pivots.len()` being the pivot rows.
    pub(super) rows: Vec<Vec<BigInt>>,
    /// Pivot columns in ascending order.
    pub(super) pivots: Vec<usize>,
    /// Input row index of each pivot row.
    pub(super) origins: Vec<usize>,
}

/// Replace `row` by `row * numerator / denominator` (an exact integer division by
/// the theory above); `false` if a division was not exact.
pub(super) fn rescale_row(row: &mut [BigInt], numerator: &BigInt, denominator: &BigInt) -> bool {
    for entry in row.iter_mut().filter(|entry| !entry.is_zero()) {
        *entry *= numerator;
        let (quotient, remainder) = entry.div_rem(denominator);
        if !remainder.is_zero() {
            return false;
        }
        *entry = quotient;
    }
    true
}

/// An integer row that is `scale` times the rational row, and `scale`.
pub(super) fn integer_row(row: &[Q]) -> (Vec<BigInt>, BigInt) {
    let mut scale = BigInt::one();
    for entry in row {
        if !entry.denom().is_one() {
            scale = scale.lcm(entry.denom());
        }
    }
    let entries = if scale.is_one() {
        row.iter().map(|entry| entry.numer().clone()).collect()
    } else {
        row.iter()
            .map(|entry| entry.numer() * (&scale / entry.denom()))
            .collect()
    };
    (entries, scale)
}

/// How [`integer_rref`] chooses pivot rows and which rows it updates.
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum Elimination {
    /// Gauss-Jordan elimination of every other row, the pivot row of a column being the
    /// first row with a nonzero entry (the pivot rows and their order are those of the
    /// rational reference elimination).
    GaussJordan,
    /// Fraction-free forward elimination (Bareiss echelon form: only the rows below the
    /// pivot row are updated), the pivot row of a column being the row with a nonzero
    /// entry there whose entries from that column on have the smallest total bit
    /// length, which keeps the integers small.  The pivot columns, the rank and the
    /// solutions do not depend on the choice of pivot rows.
    ForwardSmallestRow,
}

/// Elimination of `matrix` over its first `columns` columns (later columns are
/// carried along), as chosen by `mode`.  Returns `Ok(None)` if an exact division
/// fails (this does not happen by the Bareiss identities; the caller then uses the
/// rational elimination).
pub(super) fn integer_rref(
    matrix: &[Vec<Q>],
    columns: usize,
    mode: Elimination,
    budget: Option<&mut EliminationBudget>,
) -> Result<Option<IntegerRref>, RationalKernelError> {
    match integer_eliminate(matrix, columns, mode, budget, false)? {
        Some(IntegerEliminationResult::Echelon(reduced)) => Ok(Some(reduced)),
        Some(IntegerEliminationResult::Inconsistent) => {
            unreachable!("inconsistency detection disabled")
        }
        None => Ok(None),
    }
}

/// A solve may terminate before echelon form when an exact contradictory row appears.
pub(super) enum IntegerEliminationResult {
    Echelon(IntegerRref),
    Inconsistent,
}

/// Forward elimination of a validated augmented system, stopping as soon as an
/// integer row proves `0 = c` with `c != 0`. No rank is claimed on that path.
/// `None` requests the existing exact fallback after a non-exact division.
pub(super) fn integer_echelon_for_solve(
    matrix: &[Vec<Q>],
    columns: usize,
) -> Option<IntegerEliminationResult> {
    integer_eliminate(matrix, columns, Elimination::ForwardSmallestRow, None, true)
        .ok()
        .flatten()
}

fn integer_eliminate(
    matrix: &[Vec<Q>],
    columns: usize,
    mode: Elimination,
    budget: Option<&mut EliminationBudget>,
    stop_on_inconsistency: bool,
) -> Result<Option<IntegerEliminationResult>, RationalKernelError> {
    let row_count = matrix.len();
    let mut rows = Vec::with_capacity(row_count);
    // `scales[r]`: the rational `c` with true row = `c * rows[r]` (budgeted runs only).
    let mut scales: Vec<Q> = Vec::new();
    for row in matrix {
        let (entries, scale) = integer_row(row);
        if stop_on_inconsistency
            && !entries[columns].is_zero()
            && entries[..columns].iter().all(Zero::is_zero)
        {
            return Ok(Some(IntegerEliminationResult::Inconsistent));
        }
        rows.push(entries);
        if budget.is_some() {
            scales.push(Q::new_raw(BigInt::one(), scale));
        }
    }
    integer_eliminate_rows(rows, scales, columns, mode, budget, stop_on_inconsistency)
}

/// Forward elimination of an already integral augmented system.
pub(super) fn integer_echelon_from_integer_rows(
    rows: Vec<Vec<BigInt>>,
    columns: usize,
) -> Option<IntegerEliminationResult> {
    if rows
        .iter()
        .any(|row| !row[columns].is_zero() && row[..columns].iter().all(Zero::is_zero))
    {
        return Some(IntegerEliminationResult::Inconsistent);
    }
    integer_eliminate_rows(
        rows,
        Vec::new(),
        columns,
        Elimination::ForwardSmallestRow,
        None,
        true,
    )
    .ok()
    .flatten()
}

fn integer_eliminate_rows(
    mut rows: Vec<Vec<BigInt>>,
    mut scales: Vec<Q>,
    columns: usize,
    mode: Elimination,
    mut budget: Option<&mut EliminationBudget>,
    stop_on_inconsistency: bool,
) -> Result<Option<IntegerEliminationResult>, RationalKernelError> {
    if fallback_forced() {
        return Ok(None);
    }
    let row_count = rows.len();
    let width = rows.first().map_or(0, Vec::len);
    // `levels[r]`: the pivot at which `rows[r]` is valid (1 initially).
    let mut levels: Vec<BigInt> = vec![BigInt::one(); row_count];
    let mut origins: Vec<usize> = (0..row_count).collect();
    let mut pivots: Vec<usize> = Vec::new();
    let mut previous = BigInt::one();

    for col in 0..columns {
        let row = pivots.len();
        let pivot_index = match mode {
            Elimination::GaussJordan => (row..row_count).find(|&r| !rows[r][col].is_zero()),
            Elimination::ForwardSmallestRow => (row..row_count)
                .filter(|&r| !rows[r][col].is_zero())
                .min_by_key(|&r| rows[r][col..].iter().map(|x| x.bits() + 1).sum::<u64>()),
        };
        let Some(pivot_index) = pivot_index else {
            continue;
        };
        rows.swap(row, pivot_index);
        origins.swap(row, pivot_index);
        levels.swap(row, pivot_index);
        if budget.is_some() {
            scales.swap(row, pivot_index);
        }
        // Bring the pivot row to the current level.
        if levels[row] != previous {
            if !rescale_row(&mut rows[row], &previous, &levels[row]) {
                return Ok(None);
            }
            if budget.is_some() {
                scales[row] *= Q::new(levels[row].clone(), previous.clone());
            }
            levels[row] = previous.clone();
        }
        let pivot_row = std::mem::take(&mut rows[row]);
        let value = pivot_row[col].clone();
        let pivot_scale = if budget.is_some() {
            scales[row].clone()
        } else {
            Q::zero()
        };

        let first_target = if mode == Elimination::GaussJordan {
            0
        } else {
            row + 1
        };
        for r in first_target..row_count {
            if r == row || rows[r][col].is_zero() {
                continue;
            }
            // Fuse rescaling from the target's level with elimination. The
            // common `previous` factor cancels, so we divide by its own level.
            let divisor = &levels[r];
            let target = &mut rows[r];
            let t = target[col].clone();
            if let Some(b) = budget.as_deref_mut() {
                b.check_factor(&t, &scales[r], &value, &pivot_scale)?;
                // true row' = c * target_level / value * row'.
                scales[r] *= Q::new(divisor.clone(), value.clone());
            }
            let start = if r < row { 0 } else { col };
            for j in start..width {
                let pivot_entry = &pivot_row[j];
                let entry = &mut target[j];
                if j == col {
                    // v * t - t * v is zero without either multiplication.
                    *entry = BigInt::zero();
                } else if !(entry.is_zero() && pivot_entry.is_zero()) {
                    let mut updated = if entry.is_zero() {
                        BigInt::zero()
                    } else {
                        &*entry * &value
                    };
                    if !pivot_entry.is_zero() {
                        updated -= &t * pivot_entry;
                    }
                    if !divisor.is_one() {
                        let (quotient, remainder) = updated.div_rem(divisor);
                        if !remainder.is_zero() {
                            return Ok(None);
                        }
                        updated = quotient;
                    }
                    *entry = updated;
                }
                if j >= col {
                    if let Some(b) = budget.as_deref_mut() {
                        b.check_scaled_entry(entry, &scales[r])?;
                    }
                }
            }
            if stop_on_inconsistency
                && !target[columns].is_zero()
                && target[..columns].iter().all(Zero::is_zero)
            {
                // Every operation so far is exact and preserves the solution
                // set. A contradictory row suffices; finishing the rank does not.
                return Ok(Some(IntegerEliminationResult::Inconsistent));
            }
            levels[r] = value.clone();
        }
        // The pivot row is unchanged; as an upper row from now on its integers
        // are those of the new level (its diagonal entry is the new pivot).
        rows[row] = pivot_row;
        levels[row] = value.clone();
        previous = value;
        pivots.push(col);
    }
    origins.truncate(pivots.len());
    Ok(Some(IntegerEliminationResult::Echelon(IntegerRref {
        rows,
        pivots,
        origins,
    })))
}

/// The solution of the consistent system `[A | b]` (`columns` columns of `A`) from the
/// forward echelon form of `integer_rref`, with every free variable zero: back
/// substitution on integers.  With `D` the last pivot (the determinant of the pivot
/// minor), the entries `N_j = D * x_j` are integers by Cramer's rule, and
/// `N_i = (D * b_i - sum_(j > i) a_ij N_j) / a_ii` is an exact division.  Only the pivot
/// rows are read: the rows below them are not inspected, so the caller must already
/// have excluded a nonzero right-hand side there (an inconsistent system).  `None` only
/// if a division is not exact, which Cramer's rule rules out for a consistent system.
pub(super) fn back_substitute(echelon: &IntegerRref, columns: usize) -> Option<Vec<Q>> {
    let rank = echelon.pivots.len();
    let mut solution = vec![Q::zero(); columns];
    let Some(determinant) = rank
        .checked_sub(1)
        .map(|last| &echelon.rows[last][echelon.pivots[last]])
    else {
        return Some(solution);
    };
    let mut numerators: Vec<BigInt> = vec![BigInt::zero(); rank];
    for i in (0..rank).rev() {
        let row = &echelon.rows[i];
        let mut total = determinant * &row[columns];
        for (j, numerator) in numerators.iter().enumerate().skip(i + 1) {
            let entry = &row[echelon.pivots[j]];
            if !entry.is_zero() && !numerator.is_zero() {
                total -= entry * numerator;
            }
        }
        let (quotient, remainder) = total.div_rem(&row[echelon.pivots[i]]);
        if !remainder.is_zero() {
            return None;
        }
        numerators[i] = quotient;
    }
    for (i, numerator) in numerators.into_iter().enumerate() {
        solution[echelon.pivots[i]] = Q::new(numerator, determinant.clone());
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_row_keeps_its_raw_negative_denominator_fast_path() {
        let row = [Q::new_raw(BigInt::one(), -BigInt::one())];
        assert_eq!(integer_row(&row), (vec![BigInt::one()], BigInt::one()));
        assert_eq!(
            crate::vec_poly::clear_denominators(&row),
            vec![-BigInt::one()]
        );
    }
}
