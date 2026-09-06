//! Shared foundation for symmetric polynomial algebras.
//!
//! This crate provides the common building blocks used by `sym-poly-sym` (symmetric functions),
//! `sym-poly-qsym` (quasisymmetric functions), and `sym-poly-multipoly` (multivariate polynomials):
//!
//! - [`Ring`] trait for generic coefficients (i64, BigInt, rationals)
//! - [`Partition`] and [`Composition`] types for basis indexing
//! - [`UnivariatePolynomial`] for polynomial-valued coefficients
//! - [`Tableau`], [`SkewTableau`], and their lazy standard-tableau iterators
//! - traced inverse `P`-RS column insertion for abstract ordered alphabets
//! - [`BasisIndex`] trait abstracting over index types
//! - [`PrimeField`] for modular computations over finite prime fields
//! - [`FormalSum`] for generic linear combinations of basis elements
//! - [`matrix`] utilities for transition matrix computation
//! - [`hamel_goulden`] utilities for cutting-strip determinant examples
//! - [`linear_algebra`] utilities for exact row reduction and quotient spaces
//! - [`sparse_linear_algebra`] utilities for sparse modular row reduction
//! - [`packed_sparse_linear_algebra`] experiments for byte-sized small-prime rows
//! - [`sn_action`] utilities for small exact `S_n`-module computations
//! - [`RationalFunction`] and [`QtRationalFunction`] for `Q(q,t)` coefficient fields
//! - [`TransitionCache`] for per-algebra cached transition matrices

pub mod composition {
    pub use combinatoric_core::composition::*;
}
pub mod crt;
pub mod field;
pub mod finite_field;
pub mod finite_sn_module;
pub mod formal_sum;
pub mod hamel_goulden;
pub mod index;
pub mod linear_algebra;
pub mod matrix;
pub mod p_rs;
pub mod packed_sparse_linear_algebra;
pub mod partition {
    pub use combinatoric_core::partition::*;
}
pub mod polynomial;
pub mod rational_function;
pub mod ring {
    pub use combinatoric_core::ring::*;
}
pub mod sn_action;
pub mod sparse_linear_algebra;
pub mod ssaf;
pub mod tableau;
pub mod transition_cache;

pub use combinatoric_core::{Composition, Partition, Ring, WeakComposition};
pub use crt::{
    chinese_remainder, chinese_remainder_pair, rational_reconstruction, symmetric_residue,
};
pub use field::Field;
pub use finite_field::PrimeField;
pub use finite_sn_module::{
    left_permutation_basis_action, permutation_basis_module, right_permutation_basis_action,
    symmetric_group_permutation_basis, ungraded_permutation_basis_module, FiniteSnModule,
    PermutationBasis,
};
pub use formal_sum::FormalSum;
pub use hamel_goulden::{ContentInterval, CuttingStripSegment, OutsideDecomposition};
pub use index::BasisIndex;
pub use p_rs::{
    inverse_column_step, inverse_column_step_with_trace, reverse_complement, BoundaryCertificate,
    InverseColumnEvent, InverseColumnEventKind, InverseColumnStep, PInsertionError,
    PInsertionOrder,
};
pub use polynomial::UnivariatePolynomial;
pub use rational_function::{
    qt_coefficient, qt_constant, qt_monomial, qt_polynomial_constant, qt_rational_monomial,
    qt_unit_monomial, QtPolynomial, QtRationalFunction, RationalFunction,
};
pub use ssaf::Ssaf;
pub use tableau::{SkewTableau, StandardSkewTableauxIter, StandardTableauxIter, Tableau};
pub use transition_cache::TransitionCache;
