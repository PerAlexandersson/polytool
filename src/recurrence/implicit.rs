//! Implicit differential relations.
//!
//! Relations with derivatives on any row, the bivariate algebra behind them,
//! diagonal maps and derivative antecedents with their transport, reports with
//! reconstruction certificates, the direct search for implicit relations, the
//! lossy derivative map, and exact image recurrences of operators
//! ([`image_recurrence_of_operator`]).

use super::*;

pub(super) mod algebra;
mod diagonal;
pub use diagonal::*;
mod image_rule;
pub use image_rule::*;
mod lossy;
pub use lossy::*;
mod relation;
pub use relation::*;
mod report;
pub use report::*;
mod search;
pub use search::*;
