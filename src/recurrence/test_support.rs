//! Shared helpers for the recurrence unit tests.

use super::*;

/// Helper: assert a recurrence is found and matches expected display.
pub(super) fn assert_recurrence(polys: &[Vec<i64>], opts: &RecurrenceOptions, expected: &str) {
    let rec = find_polynomial_recurrence(polys, opts)
        .unwrap_or_else(|| panic!("Expected recurrence, got None"));
    let display = format!("{rec}");
    assert_eq!(display, expected, "\npolys: {:?}\nopts: {:?}", polys, opts);
}

pub(super) fn br(n: i64, d: i64) -> BigRational {
    BigRational::new(BigInt::from(n), BigInt::from(d))
}
