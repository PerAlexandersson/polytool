//! Core types and kernels shared by every recurrence search.
//!
//! Recurrence and vector-system types, polynomial helpers, exact
//! verification, modular prefiltering, evaluation, display, export, JSON
//! encoding and parameter counting.

use super::*;

mod counting;
pub use counting::*;
mod display;
pub use display::*;
mod evaluate;
pub use evaluate::*;
mod export;
mod json;
pub use json::*;
mod modular;
pub(super) use modular::*;
mod poly_ops;
pub(super) use poly_ops::*;
mod types;
pub use types::*;
mod vector;
pub use vector::*;
mod verify;
pub use verify::*;
