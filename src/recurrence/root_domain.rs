//! Optional root-domain information for transformed recurrence results.
//!
//! Theorem adapters (including Rolle-type arguments and implicit maps) and
//! exact root enclosure reports. None of this is needed for algebraic
//! discovery; it is attached only when a request asks for it.

use super::*;

mod analysis;
pub use analysis::*;
mod enclosure_report;
pub use enclosure_report::*;
mod implicit_maps;
pub use implicit_maps::*;
mod rolle;
use rolle::*;
