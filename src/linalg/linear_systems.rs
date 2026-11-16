//! Exact rational solving and certified kernels, provided by the shared arithmetic crate.
pub(crate) use exact_algebra::linear_systems::solve_integer_fit_system;
pub use exact_algebra::linear_systems::{
    rational_kernel_with_certificate, solve_full_rank_square_linear_system, solve_linear_system,
    solve_linear_system_with_rank, LinearSystemSolution, RationalKernel, RationalKernelError,
    RationalKernelLimits,
};
