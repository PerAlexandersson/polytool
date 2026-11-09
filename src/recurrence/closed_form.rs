//! Recurrences derived from closed forms of the coefficients.
//!
//! - `hypergeometric_triangle`: certified annihilating operators and
//!   recurrence records for the row polynomials of finite sums of
//!   hypergeometric terms `T(n, k)` given by their shift ratios
//!   ([`derive_hypergeometric_row_recurrence`],
//!   [`certify_hypergeometric_row_operator`]).

use super::*;

mod hypergeometric_triangle;
pub use hypergeometric_triangle::*;
