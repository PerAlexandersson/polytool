//! Shared dense modular kernel and residue helpers for the sparse adapter.
pub use exact_algebra::dense_modular::modular_leading_principal_minors_bigint;
pub(super) use exact_algebra::modular_arithmetic::{
    add_modulo_u64 as add_mod_u64, inverse_modulo_assumed_prime as mod_inverse_prime,
    multiply_modulo_u64 as mul_mod_u64, signed_residue_u64 as signed_mod_u64,
    subtract_modulo_u64 as sub_mod_u64,
};
pub(super) use exact_algebra::primes::is_prime_u64;
