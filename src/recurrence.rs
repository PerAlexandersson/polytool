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
//! The implementation is split by concern into private submodules whose
//! public items are re-exported here: core types, vector systems, exact
//! verification, modular prefiltering, fixed-shape and vector search,
//! display/export, row generation, JSON, parameter counting, candidate
//! enumeration, adaptive search, representation menus, complete spaces
//! (also after factorial coefficient rescaling), exact coordinate
//! transformations with recurrence transport, implicit differential
//! relations (derivatives on any row, row equations) with diagonal transport,
//! derivative antecedents, reconstruction data and certificates, and optional
//! root-domain information (theorem adapters, exact root enclosures) for
//! transformed results.

use num_bigint::BigInt;
use num_rational::Ratio;

/// Exact rational number used for every coefficient (`num_rational::Ratio<BigInt>`).
pub type BigRational = Ratio<BigInt>;

mod types;
pub use types::*;
mod vector;
pub use vector::*;
mod poly_ops;
use poly_ops::*;
mod verify;
pub use verify::*;
mod modular;
use modular::*;
mod search;
pub use search::*;
mod vector_search;
pub use vector_search::*;
mod display;
pub use display::*;
mod evaluate;
mod export;
pub use evaluate::*;
mod json;
pub use json::*;
mod counting;
pub use counting::*;
mod candidates;
use candidates::*;
mod adaptive;
pub use adaptive::*;
mod representations;
pub use representations::*;
mod space;
pub use space::*;
mod rescaling;
pub use rescaling::*;
mod transform;
pub use transform::*;
mod root_domain;
pub use root_domain::*;
mod root_enclosure_report;
pub use root_enclosure_report::*;
mod implicit;
mod implicit_algebra;
pub use implicit::*;
mod diagonal;
pub use diagonal::*;
mod implicit_report;
pub use implicit_report::*;

#[cfg(test)]
mod test_support;
