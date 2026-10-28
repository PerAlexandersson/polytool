//! Veronese-section lifts: recognize a recurrence-defined sequence as a
//! coefficient section of a parent sequence whose real-rootedness transfers.
//!
//! The entry point is [`find_veronese_section_lifts`]; its documentation
//! states the transfer theorem, the supported input class, the algorithm and
//! the argument that the finitely many identifying rows prove the identity
//! for all indices.  The group is split into
//!
//! - `veronese_algebra`: coefficient sections, bivariate polynomials in the
//!   shift variable `L`, the Veronese norm (a determinant) and monic division;
//! - `veronese_search`: input validation, factorization of `chi(L, x^r)`
//!   ([`crate::factorization`]), divisor enumeration and the per-candidate
//!   work behind a crate-internal entry point that takes `G` directly;
//! - `veronese_report`: the report types and their
//!   `polytool.recurrence-veronese-lifts.v1` JSON form.
//!
//! Related: [`RecurrenceJson`] (a lift's parent is exported in that format),
//! [`crate::is_real_rooted_bigint_coeffs`] (the exact parent-row checks).

use super::*;

mod veronese_algebra;
mod veronese_report;
pub use veronese_report::*;
mod veronese_search;
pub use veronese_search::{
    find_veronese_section_lifts, VeroneseSectionLiftOptions, MAX_VERONESE_SECTION_MODULUS,
};

#[cfg(test)]
mod veronese_tests;
