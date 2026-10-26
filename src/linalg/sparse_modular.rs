//! Sparse Gaussian elimination over a prime field `F_p`.
//!
//! [`SparseModRow`] stores augmented rows sorted by column;
//! the `sparse_modular_linear_system_*` entry points decide consistency, report
//! rank statistics and optionally return one solution.  Row order and pivoting
//! are chosen with [`SparseModEliminationOptions`].
//!
//! This file holds the row type; the elimination driver, its options, errors
//! and statistics live in the private `solver` submodule (a child module, so it
//! reads the row's private fields directly) and are re-exported here.

use super::modular::{add_mod_u64, mod_inverse_prime, mul_mod_u64, signed_mod_u64, sub_mod_u64};

mod solver;
#[cfg(test)]
mod tests;

pub use solver::{
    sparse_modular_linear_system_consistency,
    sparse_modular_linear_system_consistency_with_options, sparse_modular_linear_system_consistent,
    sparse_modular_linear_system_solution, SparseModEliminationError, SparseModEliminationOptions,
    SparseModEliminationResult, SparseModEliminationStats, SparseModRowOrder,
};
pub(crate) use solver::{
    sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode,
    SparseModSolutionMode,
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
