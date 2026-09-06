//! Multivariate polynomials with divided difference operators and nonsymmetric bases.
//!
//! This crate provides:
//! - [`MultiPoly<C>`] — sparse multivariate polynomials
//! - [`MultiPolyFunction<C>`] — polynomials in nonsymmetric bases (key, atom, slide)
//! - [`MultiPolyBasis`] — basis enum (Monomial, Key, Atom, MonSlide, FundSlide)
//! - [`monomial_order`] — monomial orders and leading terms
//! - [`division`] — multivariate division and normal forms
//! - [`groebner`] — basic Buchberger algorithm for small exact quotients
//! - [`quotient`] — finite standard-monomial quotient bases
//! - [`operators`] — simple, Demazure, and t-deformed operators (∂_i, π_i, θ_i)
//! - [`key_polynomial`] — Demazure characters via π operators
//! - [`atom_polynomial`] — Demazure atoms via θ operators
//! - [`borodin_wheeler`] — local weights for the Borodin--Wheeler vertex model
//! - [`flagged_schur`] — flagged Schur and flagged skew Schur polynomials
//! - [`kohnert`] — Kohnert diagrams and Assaf Yamanouchi tests
//! - [`lock_polynomial`] — finite lock polynomials from lock fillings
//! - [`multiline_queue`] — multiline queues and Ferrari--Martin labelings
//! - [`nonsymmetric_macdonald`] — nonsymmetric Macdonald filling formulas
//!   and the operator-side `q = 0` Hall-Littlewood specialization
//! - [`schubert_polynomial`] — Schubert polynomials via divided differences
//! - [`beta_grothendieck_polynomial`] — connective-K Grothendieck polynomials
//! - [`grothendieck_to_lascoux`] — Grothendieck-to-Lascoux expansion
//! - [`slide_polynomial`] — monomial slide, fundamental slide, and glide polynomials

pub mod atom_polynomial;
pub mod basis;
pub mod borodin_wheeler;
pub mod division;
pub mod flagged_schur;
pub mod graded_quotient;
pub mod groebner;
pub mod grothendieck_polynomial;
pub mod indexed_variables;
pub mod key_polynomial;
pub mod kohnert;
pub mod lock_polynomial;
pub mod lorentzian;
pub mod modular_groebner;
pub mod monomial_order;
pub mod multiline_queue;
pub mod multipoly;
pub mod multipoly_function;
pub mod nonsymmetric_macdonald;
pub mod operators;
pub mod quotient;
pub mod quotient_module;
pub mod schubert_polynomial;
pub mod slide_polynomial;
pub mod symmetric_polynomials;
pub mod transition;

pub use atom_polynomial::{atom_polynomial, t_atom_polynomial};
pub use basis::MultiPolyBasis;
pub use borodin_wheeler::{borodin_wheeler_l_weight, BorodinWheelerFaceWeight};
pub use division::{
    divide_by_polynomials, matrix_normal_forms_with_leading_terms, multiply_by_monomial,
    normal_form, normal_form_with_leading_terms, DivisionResult,
};
pub use flagged_schur::{
    flagged_schur, flagged_skew_schur, flagged_skew_tableaux, flagged_tableaux,
    row_interval_flagged_schur, row_interval_flagged_skew_schur,
    row_interval_flagged_skew_tableaux, row_interval_flagged_tableaux,
};
pub use graded_quotient::{
    graded_quotient_component, monomials_with_multidegree, polynomial_multidegree,
    GradedQuotientComponent,
};
pub use groebner::{
    buchberger_basis, buchberger_basis_with_options, buchberger_basis_with_stats,
    is_groebner_basis, make_monic, reduced_groebner_basis, reduced_groebner_basis_with_options,
    reduced_groebner_basis_with_stats, s_polynomial, BuchbergerStats, GroebnerBasis,
    GroebnerComputation, GroebnerOptions,
};
pub use grothendieck_polynomial::{
    beta_grothendieck_polynomial, beta_grothendieck_to_lascoux, grothendieck_polynomial,
    grothendieck_to_lascoux, lascoux_polynomial_by_operators, polynomial_to_lascoux,
    LascouxExpansion,
};
pub use indexed_variables::{
    ideal_generators_are_invariant_under_index_permutation, ideal_generators_are_sn_invariant,
    is_multidegree_preserving_action_matrix,
    quotient_action_matrices_by_index_permutation_and_multidegree,
    quotient_action_matrices_by_multidegree_and_cycle_type,
    quotient_action_matrix_multidegree_blocks, quotient_basis_multidegrees, IndexedVariables,
};
pub use key_polynomial::{key_polynomial, t_key_polynomial};
pub use kohnert::{
    canonical_labeling, cells_in_col, column_pairing, diagram_from_labeling, diagram_weight,
    format_diagram, ghost_diagram_weight, is_yamanouchi, k_kohnert_diagrams,
    k_kohnert_diagrams_for_composition, k_kohnert_moves, k_kohnert_weight_counts, key_diagram,
    kohnert_diagrams, kohnert_diagrams_for_composition, kohnert_moves, kohnert_polynomial,
    kohnert_polynomial_for_composition, kohnert_weight_counts, label_pairing, lascoux_polynomial,
    max_col, rectify_labeled, rectify_labeled_column_star, rothe_diagram, sorted_rows_in_col,
    yamanouchi_diagrams, Cell, Diagram, GhostDiagram, Labeling,
};
pub use lock_polynomial::lock_polynomial;
pub use lorentzian::{
    is_lorentzian, is_lorentzian_bool, is_m_convex, is_normalized_lorentzian,
    is_normalized_lorentzian_bool, is_strictly_lorentzian, is_strictly_normalized_lorentzian,
    support_is_m_convex, LorentzianResult,
};
pub use modular_groebner::{
    crt_lift_prime_field_basis_pair, crt_lift_prime_field_polynomial_pair,
    groebner_leading_monomials, modular_groebner_basis_i64_mod_prime,
    modular_groebner_basis_rational_i64_mod_prime, modular_leading_monomials_match_i64,
    modular_leading_monomials_match_rational_i64, rational_reconstruct_prime_field_basis_pair,
    rational_reconstruct_prime_field_polynomial_pair, reduce_i64_polynomial_mod_prime,
    reduce_i64_polynomials_mod_prime, reduce_rational_i64_mod_prime,
    reduce_rational_i64_polynomial_mod_prime, reduce_rational_i64_polynomials_mod_prime,
    ModularGroebnerError,
};
pub use monomial_order::{
    leading_term, monomial_divides, monomial_quotient, LeadingTerm, MonomialOrder,
    ParseMonomialOrderError,
};
pub use multiline_queue::{
    multiline_queue_weight_counts, multiline_queues_with_row_sizes, FerrariMartinLabeling,
    MultilineQueue, MultilineQueuePairing,
};
pub use multipoly::MultiPoly;
pub use multipoly_function::MultiPolyFunction;
pub use nonsymmetric_macdonald::{
    nonsymmetric_hall_littlewood, nonsymmetric_macdonald_filling_formula,
    nonsymmetric_macdonald_q0, permuted_basement_macdonald_filling_formula,
};
pub use operators::{
    demazure_lascoux_partial_i, demazure_lascoux_partial_word, demazure_lascoux_pi_i,
    demazure_lascoux_pi_word, partial_i, partial_word, pi_i, pi_word, theta_i, theta_word, tpi_i,
    tpi_word, ttheta_i, ttheta_word,
};
pub use quotient::{
    is_degree_preserving_action_matrix, normal_form_in_basis, permute_monomial, permute_variables,
    pure_power_bounds, quotient_action_matrices_by_permutation_and_degree,
    quotient_action_matrix_by_permutation, quotient_action_matrix_degree_blocks, quotient_basis,
    quotient_basis_degrees, quotient_basis_index, quotient_coordinates, restrict_matrix_to_indices,
    standard_monomials_from_leading_monomials, QuotientBasis,
};
pub use quotient_module::{PolynomialQuotientSnModule, PolynomialQuotientSnModuleError};
pub use schubert_polynomial::{
    schubert_polynomial, schubert_to_atom, schubert_to_fund_slide, schubert_to_key,
    schubert_to_monomial,
};
pub use slide_polynomial::{
    fundamental_slide_expansion, fundamental_slide_polynomial, glide_polynomial,
    monomial_slide_polynomial,
};
pub use symmetric_polynomials::{elementary_symmetric_generators, elementary_symmetric_polynomial};
