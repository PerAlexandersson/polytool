//! Exact arithmetic kernels independent of consumer polynomial and combinatorial types.
//!
//! Matrices use rows of `BigInt` or `BigRational` entries. Dense solvers retain
//! exact fraction-free elimination and rational fallbacks; reconstruction APIs
//! require explicit mathematical bounds and return typed errors. Sparse Gröbner
//! kernels accept a borrowed consumer representation or the owned rational type.

#![doc = include_str!("../README.md")]

pub mod chinese_remainder;
pub mod denominators;
pub mod dense_modular;
pub mod groebner;
pub mod integer_elimination;
pub mod linear_systems;
pub mod modular_arithmetic;
pub mod multimodular;
pub mod nonnegative_matrix_product;
pub mod primes;
pub mod rational_reconstruction;
pub mod sparse_modular;

pub use dense_modular::modular_leading_principal_minors_bigint;
pub use linear_systems::{
    rational_kernel_with_certificate, solve_full_rank_square_linear_system, solve_linear_system,
    solve_linear_system_with_rank, LinearSystemSolution, RationalKernel, RationalKernelError,
    RationalKernelLimits,
};
pub use num_rational::BigRational;
type Q = BigRational;

#[cfg(test)]
mod differential_tests;
#[cfg(test)]
mod test_support;
