//! Bounded search over equivalent representations of one sequence.
//!
//! Request and report types, the map chains used to move between
//! representations, and the queue search that composes the fitting,
//! transform and implicit layers.

use super::*;

mod map;
pub use map::RepresentationMapRecord;
mod report;
pub use report::*;
mod search;
pub use search::*;
