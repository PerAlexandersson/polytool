//! Deterministic work accounting for the exact solves of the recurrence search.
//!
//! A [`WorkMeter`] counts abstract work units, never time, so a search that
//! stops when its limit is reached stops at the same step on every machine.
//! The units are a fixed model of the arithmetic of each stage, computed from
//! the dimensions of the system and not from the data:
//!
//! * a linear solve of an `r x c` system (modular or exact) costs
//!   [`elimination_work`]`(r, c) = r * c * (min(r, c) + 1)`, the field
//!   operations of dense Gaussian elimination plus the assembly of the matrix;
//! * checking a candidate against `r` equations in `c` unknowns costs `r * c`;
//! * precomputing the derivative rows and the modular images costs one unit
//!   per coefficient and derivative order (and prime);
//! * every candidate costs one unit of bookkeeping.
//!
//! A stage is charged *before* it runs and runs only when its whole cost
//! fits in the remaining limit, so the units used never exceed the limit.
//! The first stage that does not fit marks the meter exhausted; the callers
//! stop as soon as they see it.

use std::cell::Cell;

/// A shared counter of deterministic work units with an optional limit.
#[derive(Debug)]
pub(in crate::recurrence) struct WorkMeter {
    used: Cell<u64>,
    limit: Option<u64>,
    exhausted: Cell<bool>,
}

impl WorkMeter {
    /// A meter that stops charging once `limit` units would be exceeded
    /// (`None`: only count).
    pub(in crate::recurrence) fn new(limit: Option<u64>) -> Self {
        Self {
            used: Cell::new(0),
            limit,
            exhausted: Cell::new(false),
        }
    }

    /// Charge `units`.  Returns `false`, without charging, when they do not
    /// fit in the remaining limit (and from then on for every charge).
    pub(in crate::recurrence) fn try_charge(&self, units: u64) -> bool {
        if self.exhausted.get() {
            return false;
        }
        let total = self.used.get().saturating_add(units);
        if self.limit.is_some_and(|limit| total > limit) {
            self.exhausted.set(true);
            return false;
        }
        self.used.set(total);
        true
    }

    /// Units charged so far.
    pub(in crate::recurrence) fn used(&self) -> u64 {
        self.used.get()
    }

    /// True once a charge did not fit.
    pub(in crate::recurrence) fn is_exhausted(&self) -> bool {
        self.exhausted.get()
    }
}

/// Work units of one dense linear solve of `rows x columns` (see the module
/// documentation).
pub(in crate::recurrence) fn elimination_work(rows: usize, columns: usize) -> u64 {
    let rows = rows as u64;
    let columns = columns as u64;
    rows.saturating_mul(columns)
        .saturating_mul(rows.min(columns).saturating_add(1))
}

/// Work units of evaluating `rows` equations in `columns` unknowns.
pub(in crate::recurrence) fn evaluation_work(rows: usize, columns: usize) -> u64 {
    (rows as u64).saturating_mul(columns as u64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_charge_that_does_not_fit_is_refused_and_sticky() {
        let meter = WorkMeter::new(Some(10));
        assert!(meter.try_charge(6));
        assert!(!meter.try_charge(5));
        assert!(meter.is_exhausted());
        assert!(!meter.try_charge(1));
        assert_eq!(meter.used(), 6);
        let counting = WorkMeter::new(None);
        assert!(counting.try_charge(u64::MAX));
        assert!(counting.try_charge(1));
        assert!(!counting.is_exhausted());
    }

    #[test]
    fn solve_cost_is_dimension_based() {
        assert_eq!(elimination_work(4, 3), 4 * 3 * 4);
        assert_eq!(evaluation_work(4, 3), 12);
    }
}
