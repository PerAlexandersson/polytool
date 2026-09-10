//! Shared library for foundational combinatorics: partitions, set partitions,
//! compositions, permutations, graphs, posets, and related utilities.

pub mod composition;
pub mod graph;
pub mod key_polynomial;
pub mod meander;
pub mod partition;
pub mod permutation;
pub mod poset;
pub mod ring;
pub mod set_partition;
pub mod sparse_matrix;

// Top-level re-exports for convenience
pub use composition::{Composition, WeakComposition};
pub use graph::Graph;
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
pub use sparse_matrix::{MutableSparseMatrix, SparseMatrix, SparseMatrixBuilder, SparseMatrixError};
pub use set_partition::{
    ordered_set_partitions, set_partitions, OrderedSetPartition, SetPartition,
};
