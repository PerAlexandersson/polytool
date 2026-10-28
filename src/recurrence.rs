//! Find scalar and coupled linear recurrences among polynomial sequences.
//!
//! Given polynomials P_1(t), P_2(t), ..., P_m(t), searches for a recurrence:
//!
//!   f(n,t) P_n(t) = sum_{r,d} c_{r,d}(n,t) D^d P_{n-r}(t)  [+ g(n,t)]
//!
//! where c_{r,d}(n,t) are polynomial coefficients in n and t, D^d is the
//! d-th derivative in t, and f(n,t) is an optional LHS denominator.  When
//! `alternating_sign` is enabled, the search also allows terms of the form
//! (-1)^n c_{r,d}(n,t) D^d P_{n-r}(t).
//!
//! Coupled states `F_n = (P_n, Q_n, ...)^T` can also be fitted to a first-order
//! affine system `F_{n+1} = M(n,t,D) F_n + G(n,t)`. Matrix entries use the
//! normal-ordered Weyl form `sum_d c_d(n,t) D^d`; exact row ranks and nullities
//! are reported, and complete final transitions can be held out for exact
//! verification.
//!
//! Both searches reduce to solving linear systems over the rationals.
//!
//! The implementation lives in private submodule groups whose public items
//! are re-exported here:
//!
//! - `base`: core types, vector systems, polynomial helpers, exact
//!   verification, modular prefiltering, evaluation (including the lazy
//!   [`RecurrenceRowIterator`]), display/export, JSON and parameter counting;
//! - `fit`: fixed-shape and vector search, candidate enumeration, adaptive
//!   search, representation menus, and complete bounded spaces with
//!   certificates (also after factorial coefficient rescaling);
//! - `transform`: exact coordinate transformations with transport of
//!   recurrences, gamma coordinates, vector transport and reports;
//! - `implicit`: implicit differential relations, diagonal transport,
//!   derivative antecedents, reconstruction certificates, the direct search
//!   [`search_implicit_relations`] and the lossy derivative map
//!   ([`lossy_reconstruction`]);
//! - `root_domain`: optional root-domain information (theorem adapters, exact
//!   root enclosures) for transformed results;
//! - `representation`: the bounded queue search over equivalent
//!   representations (scalar and vector) that composes all of these
//!   ([`search_representations`]);
//! - `lift`: parent sequences whose real-rootedness transfers to the given
//!   sequence: Veronese-section lifts ([`find_veronese_section_lifts`]) and
//!   symmetric-companion lifts ([`find_symmetric_companion_lifts`]).
//!
//! Report, search and replay functions that can reject their input return
//! [`RecurrenceError`] (and [`parse_rational_coeff`] a [`ParseRationalError`]);
//! the core types keep their own finer-grained errors
//! ([`RecurrenceEvaluationError`], [`VectorRecurrenceFitError`],
//! [`TransformError`], ...).  Each `Display` text is the full message.

use num_bigint::BigInt;

/// Exact rational number used for every coefficient (`num_rational::Ratio<BigInt>`).
pub use crate::BigRational;

mod base;
pub use base::*;
mod error;
pub use error::{ParseRationalError, RecurrenceError};
mod fit;
pub use fit::*;
mod implicit;
pub use implicit::*;
mod lift;
pub use lift::*;
mod representation;
pub use representation::*;
mod root_domain;
pub use root_domain::*;
mod transform;
pub use transform::*;

#[cfg(test)]
mod test_support;
