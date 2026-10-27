//! Linear algebra over exact rationals (BigInt/BigRational).
//!
//! Provides fraction-free Bareiss elimination, Gaussian elimination, positive
//! definiteness checks, determinants, exact characteristic polynomials
//! ([`characteristic_polynomial_bigint`]) and closed-walk counts
//! ([`closed_walk_counts_bigint`]), linear system solving, and total
//! non-negativity (TNN) checking.
//!
//! # Total non-negativity
//!
//! A matrix is *totally non-negative* (TNN) if every minor (determinant of every
//! square submatrix) is non-negative. The brute-force check
//! ([`check_total_positivity`]) enumerates all submatrices up to a given size,
//! which is exponential. The Neville elimination method ([`check_tnn_neville`])
//! checks TNN for *all* minors in O(n²m) time by performing adjacent-row
//! subtraction and verifying that no entry becomes negative.

// The implementation lives in private submodules; everything public is
// re-exported below at its historical `linalg::*` path.
mod definiteness;
mod determinants;
mod linear_systems;
mod modular;
mod sparse_modular;
#[cfg(test)]
mod test_support;
mod total_positivity;

pub use definiteness::{
    is_positive_definite, is_positive_definite_bareiss, is_positive_definite_modular,
    is_positive_semidefinite, MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD,
};
pub use determinants::{
    bareiss_determinant_bigint, bareiss_determinant_polynomial_bigint,
    bareiss_leading_principal_minors_bigint, bareiss_leading_principal_minors_polynomial_bigint,
    characteristic_polynomial_bigint, closed_walk_counts_bigint, determinant,
};
pub use linear_systems::{
    rational_kernel_with_certificate, solve_full_rank_square_linear_system, solve_linear_system,
    solve_linear_system_with_rank, LinearSystemSolution, RationalKernel, RationalKernelError,
    RationalKernelLimits,
};
pub use modular::modular_leading_principal_minors_bigint;
#[allow(deprecated)]
pub use sparse_modular::sparse_modular_linear_system_consistent;
pub use sparse_modular::{
    is_sparse_modular_linear_system_consistent, sparse_modular_linear_system_consistency,
    sparse_modular_linear_system_consistency_with_options, sparse_modular_linear_system_solution,
    SparseModEliminationError, SparseModEliminationOptions, SparseModEliminationResult,
    SparseModEliminationStats, SparseModRow, SparseModRowError, SparseModRowOrder,
};
pub(crate) use sparse_modular::{
    sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode,
    SparseModSolutionMode,
};
pub use total_positivity::{
    check_tnn_neville, check_tnn_neville_bigint, check_total_positivity,
    is_totally_nonnegative_by_minors, is_totally_nonnegative_neville, TotalPositivityError,
};
#[allow(deprecated)]
pub use total_positivity::{is_tnn, is_totally_nonnegative};

type Q = crate::BigRational;
