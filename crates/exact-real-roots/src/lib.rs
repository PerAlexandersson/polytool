//! Exact real-root counting over the integers.
//!
//! This module provides the default exact backend for boolean
//! real-rootedness. Degrees at most four use exact discriminant and quartic
//! invariant formulas. Degree five uses an unreduced signed pseudo-remainder
//! sequence; from degree six, a Collins--Brown subresultant PRS removes known
//! exact factors without computing coefficient-content GCDs. One-signed
//! inputs first try cheap coefficient filters. These paths avoid the large
//! exact PSD matrices used by the Bézout criterion.
//!
//! The positive-coefficient path is intended for combinatorial polynomials.  If
//! all coefficients have one sign, roots can only be non-positive; after
//! removing powers of `t`, real-rootedness of `f(t)` is equivalent to all roots
//! of `f(-t)` being positive.
//!
//! Public functions in this module count distinct roots of the square-free
//! part.  Real-rootedness tests compare that count with the square-free degree,
//! so repeated real roots are handled without needing multiplicities from the
//! Sturm sequence.
//!
//! An independent Uspensky/Descartes implementation is also available for
//! exact comparison.  It bounds the roots, applies homographic subdivision,
//! and uses Descartes' rule of signs to discard or certify open intervals.
//!
//! # Entry points
//!
//! - Decide real-rootedness: [`is_real_rooted_fast_bigint_coeffs`] (adaptive
//!   default, also [`is_real_rooted_fast_i64`]);
//!   [`is_real_rooted_prs_bigint_coeffs`] and
//!   [`is_real_rooted_uspensky_bigint_coeffs`] force one backend.
//! - Count distinct real roots: [`count_real_roots_prs_bigint_coeffs`],
//!   [`count_positive_roots_prs_bigint_coeffs`] and the Uspensky counterparts;
//!   [`squarefree_degree_bigint_coeffs`] gives the square-free degree to compare
//!   against.
//! - Isolate roots: [`isolate_real_roots_bigint_coeffs`] returns sorted
//!   [`IsolatedRealRoot`] values with exact dyadic endpoints.
//! - Signs at roots: [`count_root_signs_prs_bigint_coeffs`] and
//!   [`tarski_query_prs_bigint_coeffs`].
//! - Coefficient-only filters: [`satisfies_newton_inequalities_bigint_coeffs`]
//!   (necessary) and [`satisfies_kurtz_condition_bigint_coeffs`] (sufficient).
//!
//! # Conventions
//!
//! Polynomials are `&[BigInt]` (or `&[i64]` in the convenience wrappers) in
//! ascending degree order: `coeffs[i]` is the coefficient of `t^i`.  Trailing
//! zeros are ignored, and the zero polynomial and nonzero constants have no
//! roots and count as real-rooted.  Repeated roots are counted once, so counts
//! are numbers of *distinct* roots.  All results are exact.

mod arithmetic;
mod coefficient_tests;
mod isolation;
mod prs;
mod uspensky;

pub use arithmetic::bigint_coeffs_to_i64;
pub use coefficient_tests::{
    is_real_rooted_fast_bigint_coeffs, is_real_rooted_fast_i64,
    is_real_rooted_one_signed_bigint_coeffs, satisfies_kurtz_condition_bigint_coeffs,
    satisfies_newton_inequalities_bigint_coeffs,
};
pub use isolation::{isolate_real_roots_bigint_coeffs, DyadicRational, IsolatedRealRoot};
pub use prs::{
    count_positive_roots_prs_bigint_coeffs, count_real_roots_prs_bigint_coeffs,
    count_root_signs_prs_bigint_coeffs, is_real_rooted_prs_bigint_coeffs,
    primitive_sturm_max_coefficient_bits, squarefree_degree_bigint_coeffs,
    tarski_query_prs_bigint_coeffs, RootSignCounts,
};
pub use uspensky::{
    count_positive_roots_uspensky_bigint_coeffs, count_real_roots_uspensky_bigint_coeffs,
    is_real_rooted_uspensky_bigint_coeffs, is_real_rooted_uspensky_i64,
};

#[doc(hidden)]
pub use arithmetic::squarefree_part_primitive_bigint_coeffs;
#[allow(deprecated)]
pub use coefficient_tests::{
    satisfies_kurtz_condition_bigint, satisfies_newton_inequalities_bigint,
};

#[cfg(test)]
mod tests;
