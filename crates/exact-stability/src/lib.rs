//! Exact root counts for Gaussian-rational univariate polynomials.
//!
//! [`check_univariate_stability`] counts roots with multiplicity in an open
//! circular region, on its boundary, and outside. Coefficients are ascending;
//! trailing zeros are ignored. All decisions use rational arithmetic.
//!
//! Upper/lower half-plane stability means avoidance; left/right half-plane,
//! disk and exterior stability mean that every root lies strictly inside.
//! [`StabilityConvention`] allows either reading of the same exact counts.
//! The zero polynomial is stable by the avoidance convention but has no
//! finite root counts (`counts: None`) and fails the all-inside convention.

pub mod gaussian;
pub mod parse;
pub mod univariate;

pub use exact_algebra::BigRational;
pub use gaussian::GaussianRational;
pub use univariate::mobius::{MobiusTransform, TransformedPolynomial};
pub use univariate::region::StabilityRegion;
pub use univariate::{
    check_univariate_stability, check_univariate_stability_with_checkpoint,
    check_univariate_stability_with_convention, RootCounts, StabilityCertificate,
    StabilityConvention, StabilityError, StabilityMethod, StabilityReport,
};
