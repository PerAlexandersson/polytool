//! Implicit differential relations.
//!
//! Relations with derivatives on any row, the bivariate algebra behind them,
//! diagonal maps and derivative antecedents with their transport, reports with
//! reconstruction certificates, the direct search for implicit relations, and
//! the lossy derivative map.

use super::*;

pub(super) mod algebra;
mod diagonal;
pub use diagonal::*;
mod lossy;
pub use lossy::*;
mod relation;
pub use relation::*;
mod report;
pub use report::*;
mod search;
pub use search::*;
