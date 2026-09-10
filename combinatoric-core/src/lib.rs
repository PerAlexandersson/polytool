//! Shared library for foundational combinatorics: partitions, set partitions,
//! compositions, permutations, graphs, posets, and related utilities.

pub mod chain_complex;
pub mod composition;
pub mod graph;
pub mod integer_linear_algebra;
pub mod key_polynomial;
pub mod meander;
pub mod partition;
pub mod permutation;
pub mod poset;
pub mod ring;
pub mod set_partition;
pub mod sparse_matrix;

// Top-level re-exports for convenience
pub use chain_complex::{
    cancel_units, field_betti_number, integral_homology, is_prime_u64, replay_unit_cancellation,
    replay_unit_cancellation_and_verify, universal_coefficient_dimension, AbelianGroup,
    ChainComplexError, FiniteChainComplex, UnitCancellationCertificate, UnitPivot,
    UnitReductionOptions, UnitReductionResult, UnitReductionStats,
};
pub use composition::{Composition, WeakComposition};
pub use graph::Graph;
pub use integer_linear_algebra::{
    replay_smith_operations, replay_smith_operations_bounded, smith_normal_form,
    verify_smith_certificate, SmithError, SmithLimit, SmithNormalForm, SmithOperation,
    SmithOptions, SmithReplayOptions,
};
pub use meander::{
    is_connected_arch_pair, noncrossing_perfect_matchings, rooted_meandric_permutation_count,
    rooted_meandric_permutation_from_arch_pair, rooted_meandric_permutations,
    CLOSED_MEANDRIC_NUMBERS_INITIAL,
};
pub use partition::Partition;
pub use permutation::{
    all_permutations_one_indexed, all_permutations_zero_indexed, assert_one_indexed_permutation,
    compose_permutations, inverse_permutation, is_one_indexed_permutation, longest_permutation,
    next_permutation, optimist_sort_derangement_step_distribution, optimist_sort_step,
    optimist_sort_step_distribution, optimist_sort_step_distribution_via_derangements,
    optimist_sort_steps, optimist_sort_steps_word, permutation_from_simple_transpositions,
    reduced_word, stable_standardization, unfixed_standardization,
};
pub use ring::Ring;
pub use set_partition::{
    ordered_set_partitions, set_partitions, OrderedSetPartition, SetPartition,
};
pub use sparse_matrix::{
    MutableSparseMatrix, SparseMatrix, SparseMatrixBuilder, SparseMatrixError, SparseMatrixLimits,
};
