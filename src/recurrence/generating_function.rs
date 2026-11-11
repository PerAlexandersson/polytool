//! Exact rational row generating functions for constant-coefficient recurrences.
//!
//! [`rational_generating_function`] includes **all** prescribed rows, including
//! prefix rows before the rule starts. The coefficient of `z^m` is
//! `P_(first_index + m)(x)`. [`riordan_decomposition`] recognizes
//! `d(z)/(1 - x h(z))`, without assuming `h'(0) != 0`.
//!
//! JSON uses ascending coefficient arrays of rational strings; bivariate
//! arrays are indexed `[z_power][x_power]`. This module is independent of the
//! proof routes and is suitable for catalog indexing and output adapters.

mod rational;
pub(in crate::recurrence) use rational::{decode_bivariate, encode_bivariate};
pub use rational::{
    rational_generating_function, riordan_decomposition, RationalGeneratingFunction,
    RationalGeneratingSeries, RiordanDecomposition, RATIONAL_GENERATING_FUNCTION_SCHEMA,
};
