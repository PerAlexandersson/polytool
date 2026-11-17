//! Compatibility re-exports of the shared sparse modular engine.
#[allow(deprecated)]
pub use exact_algebra::sparse_modular::sparse_modular_linear_system_consistent;
pub use exact_algebra::sparse_modular::{
    is_sparse_modular_linear_system_consistent, sparse_modular_linear_system_consistency,
    sparse_modular_linear_system_consistency_with_options, sparse_modular_linear_system_solution,
    SparseModEliminationError, SparseModEliminationOptions, SparseModEliminationResult,
    SparseModEliminationStats, SparseModRow, SparseModRowError, SparseModRowOrder,
};
pub(crate) use exact_algebra::sparse_modular::{
    sparse_modular_linear_system_consistency_prime_unchecked_with_solution_mode,
    SparseModSolutionMode,
};
