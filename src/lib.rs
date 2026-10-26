#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Dense univariate polynomial toolkit for combinatorial research.
//!
//! Provides [`Polynomial<C>`] with generic coefficients, plus specialized routines
//! for real-rootedness, interlacing, gamma-positivity, Stapledon decomposition, resultants, Ehrhart theory,
//! recurrence search, and standard polynomial sequences — all with exact arithmetic.
//!
//! # Quick start
//!
//! Most property-checking functions accept `&[i64]` coefficient vectors in ascending
//! degree order (`coeffs[i]` = coefficient of t^i):
//!
//! ```
//! use polytool::*;
//!
//! // Eulerian polynomial A_4(t) = 1 + 11t + 11t^2 + t^3
//! assert!(is_real_rooted(&[1, 11, 11, 1]));     // exact default path
//! assert!(is_palindromic(&[1, 11, 11, 1]));
//! assert!(is_gamma_positive(&[1, 11, 11, 1]));   // gamma = [1, 8]
//!
//! // Strict interlacing (degree diff = 1, Bézout-based)
//! assert_eq!(check_interlacing(&[8, -6, 1], &[-15, 23, -9, 1]), Some(true));
//! ```
//!
//! For polynomial arithmetic, use [`Polynomial<C>`]:
//!
//! ```
//! use polytool::Polynomial;
//!
//! let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]); // 1 + t
//! let q = p.clone() * p.clone(); // 1 + 2t + t^2
//! assert_eq!(format!("{}", q), "1 + 2t + t^2");
//! assert!(q.is_palindromic());
//! ```
//!
//! # API conventions
//!
//! Coefficient vectors are always in ascending degree order. Public `i64`
//! functions are ergonomic wrappers around exact integer implementations where
//! practical; use the `*_bigint_coeffs` variants when coefficients may exceed
//! `i64`.
//!
//! Interlacing checks return `Option<bool>` because the degree relation is part
//! of the contract. `Some(true)` and `Some(false)` mean the degree relation is
//! valid and the property was decided exactly; `None` means the input degrees
//! are not valid for that directed interlacing test. For the undirected
//! question "do these two polynomials interlace in some order?", use
//! [`interlaces_either_order`], which returns a plain `bool`.
//!
//! ```
//! use num_bigint::BigInt;
//! use polytool::check_interlacing_bigint_coeffs;
//!
//! let center = BigInt::from(10).pow(20);
//! let f = vec![-&center, BigInt::from(1)]; // t - center
//! let g = vec![
//!     (&center - 1u32) * (&center + 1u32),
//!     -BigInt::from(2) * &center,
//!     BigInt::from(1),
//! ]; // (t - center + 1)(t - center - 1)
//! assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
//! ```
//!
//! # Modules
//!
//! - [`polynomial`] — `Polynomial<C>` with `CoeffRing`/`FieldRing` traits, arithmetic,
//!   derivative, evaluate, shift, reverse, dilate, GCD, division, Lagrange interpolation
//! - [`vec_poly`] — Arithmetic on plain coefficient vectors `&[C]`, generic
//!   over [`CoeffRing`]: add/sub/mul, powers, composition, derivative,
//!   reversal, Wronskian, polynomial-matrix permanent
//! - [`numeric`] — Exact binomials, factorials, falling factorials,
//!   multinomials, and a checked `u128` binomial
//! - [`linalg`] — Exact linear algebra: fraction-free Bareiss elimination,
//!   modular/CRT determinant reconstruction, Gaussian elimination over ℚ,
//!   sparse finite-field consistency/solution checks, positive
//!   definiteness/semi-definiteness, determinants, linear system solving, total
//!   non-negativity via Neville elimination
//! - [`real_rootedness`] — exact real-rootedness wrappers, Bézout/Hermite
//!   comparison paths, strict/weak interlacing (including same-degree via
//!   Cauchy bound reduction), log-concavity, ultra-log-concavity, palindromic
//!   check, gamma-positivity, Stapledon decomposition, resultant,
//!   discriminant, Ehrhart ↔ h*-vector, display utilities
//! - [`real_rooted_enumeration`] — exact enumeration of all real-rooted
//!   integer polynomials of fixed degree with at least three prescribed
//!   leading coefficients (Greaves--Syatriadi `AllRealRooted`, certified
//!   endpoints), and the antiderivative-shift kernel behind it
//! - [`root_count`] — adaptive exact real-rootedness using primitive integer
//!   PRS and Uspensky/Descartes root counting
//! - [`sturm`] — Sturm chains for exact root counting and isolation (the slower
//!   alternative to the default real-rootedness path)
//! - [`sturm_cf`] — signed Euclidean/Sturm continued-fraction certificates
//! - [`recurrence`] — Adaptive recurrence search for polynomial sequences
//! - [`renewal`] — Exact recognition of renewal-type recurrences
//!   `P_n = sum_j r_j t^j P_(n-1-ell*j)` with interlacing consequences
//! - [`oeis`] — Bundled recurrence-backed OEIS polynomial families
//! - [`sequences`] — Standard sequences: Eulerian, Narayana, type B Eulerian,
//!   Chebyshev T/U, Hermite, shifted Legendre
//! - [`parse`] — Flexible polynomial parsing (comma/space-separated, bracketed,
//!   expanded polynomial notation)
//! - [`basis`] — Exact expansion of polynomials in prescribed bases, including
//!   the magic basis `{t^i (1+t)^{d-i}}`
//! - [`bkw`] — Beraha--Kahane--Weiss scout tools for characteristic symbols of
//!   polynomial recurrences
//! - [`decomposition`] — `I_d` / `R_d` symmetric decompositions, `f`-polynomials,
//!   alternatingly increasing checks, and Brandén--Solus-style magic-basis analysis
//! - [`brenti_sequence`] — Planar strip digraphs with height-dependent
//!   transitions; their path matrices are totally nonnegative (they do not by
//!   themselves certify real-rooted rows)
//! - [`brenti_recognition`] — Exact search for a planar strip digraph that
//!   reproduces a given prefix of polynomial rows
//! - [`interlacing_matrix`] — Finite Athanasiadis--Wagner `Lace(A)`
//!   interlacing-matrix truncations and exact finite TNN checks
//! - [`cyclic_sieving`] — Exact cyclic and bicyclic sieving checks via
//!   root-of-unity evaluation, orbit decompositions and fixed-point tables
//! - [`hstar_inequalities`] — Named inequalities for Ehrhart h*-vectors
//!   (Hibi, Stanley, Balletti--Higashitani, Stapledon)
//! - [`coefficient_tests`] — Coefficient-only real-rootedness criteria
//!   (Newton inequalities, Kurtz condition) with structured reports
//! - [`hurwitz`] — Exact Hurwitz matrices, Hurwitz determinants and the
//!   Routh--Hurwitz stability test (all roots in the open left half-plane)
//! - [`version`] — build-version formatting shared by the CLI
//! - [`tnn_network`] — Constructive planar-network certificates for
//!   lower-unitriangular totally nonnegative matrices and monic polynomial
//!   sequences, structured TNN diagnostics, and Brändén--Saud Maia Leite chain
//!   polynomials

pub mod basis;
pub mod bkw;
pub mod brenti_recognition;
pub mod brenti_sequence;
pub mod coefficient_tests;
pub mod cyclic_sieving;
pub mod decomposition;
pub mod hstar_inequalities;
pub mod hurwitz;
pub mod interlacing_matrix;
pub mod linalg;
pub mod numeric;
pub mod oeis;
pub mod polynomial;
pub mod sturm;
pub mod sturm_cf;
pub mod tnn_network;
pub mod vec_poly;
pub mod version;

pub use linalg::{
    bareiss_determinant_bigint, bareiss_determinant_polynomial_bigint,
    bareiss_leading_principal_minors_bigint, bareiss_leading_principal_minors_polynomial_bigint,
    check_tnn_neville, check_tnn_neville_bigint, check_total_positivity, determinant,
    is_positive_definite, is_positive_definite_bareiss, is_positive_definite_modular,
    is_positive_semidefinite, is_tnn, is_totally_nonnegative,
    modular_leading_principal_minors_bigint, sparse_modular_linear_system_consistency,
    sparse_modular_linear_system_consistency_with_options, sparse_modular_linear_system_consistent,
    sparse_modular_linear_system_solution, SparseModEliminationError, SparseModEliminationOptions,
    SparseModEliminationResult, SparseModEliminationStats, SparseModRow, SparseModRowError,
    SparseModRowOrder, TotalPositivityError, MODULAR_POSITIVE_DEFINITE_DIMENSION_THRESHOLD,
};

pub mod parse;
pub mod real_rooted_enumeration;
pub mod real_rootedness;
pub mod recurrence;
pub mod renewal;
pub mod root_count;
pub mod sequences;

pub use basis::{
    analyze_magic_basis_bigint, analyze_magic_basis_i64, bernstein_basis_bigint,
    bernstein_basis_coordinates, bernstein_basis_coordinates_bigint,
    bernstein_basis_coordinates_i64, coordinates_in_basis, coordinates_in_basis_bigint,
    coordinates_in_basis_i64, is_magic_positive_bigint, is_magic_positive_i64, magic_basis,
    magic_basis_coordinates_bigint, magic_basis_coordinates_i64, BasisError, MagicBasisAnalysis,
};
pub use bkw::{
    BkwError, BkwRootComputation, BkwRootInfo, BkwRootOptions, BkwScoutCandidate, BkwScoutOptions,
    BkwSymbol, Complex64,
};
pub use brenti_recognition::{
    recognize_brenti_strip_digraph, BrentiRecognition, BrentiRecognitionError,
    BrentiRecognitionOptions, BrentiWeightDomain,
};
pub use brenti_sequence::{
    build_brenti_sequence_certificate, BrentiEdge, BrentiError, BrentiSequenceCertificate,
    BrentiStripDigraph,
};
pub use coefficient_tests::{
    coefficient_test_report_bigint, kurtz_condition_report_bigint, newton_inequality_report_bigint,
    CoefficientCriterionReport, CoefficientInequalityCheck, CoefficientTestReport,
};
pub use cyclic_sieving::{
    bicyclic_sieving_report_bigint, cyclic_sieving_check_orbit_sizes,
    cyclic_sieving_orbit_decomposition_bigint, cyclic_sieving_report_bigint,
    cyclic_sieving_sequence_reports_bigint, cyclotomic_polynomial_bigint,
    fixed_points_from_orbit_counts, fixed_points_from_orbit_sizes, orbit_counts_from_fixed_points,
    orbit_counts_from_orbit_sizes, reduce_modulo_cyclic_bigint, root_of_unity_evaluation_bigint,
    BicyclicSievingReport, CyclicSievingPowerCheck, CyclicSievingReport, CyclicSievingSequenceItem,
    OrbitDecomposition, OrbitSizeEntry, OrbitSizesCheck, ResidueClass, RootOfUnityEvaluation,
    BICYCLIC_MAX_GROUP_SIZE,
};
pub use decomposition::{
    analyze_symmetric_decomposition_bigint, analyze_symmetric_decomposition_i64, f_polynomial,
    f_polynomial_i64, is_alternatingly_increasing, is_alternatingly_increasing_bigint,
    r_decomposition, r_decomposition_i64, r_transform, r_transform_i64,
    SymmetricDecompositionAnalysis, SymmetricDecompositionAnalysisBigInt,
};
pub use hstar_inequalities::{
    hstar_inequality_report_bigint, HStarInequalityCheck, HStarInequalityReport,
};
pub use hurwitz::{hurwitz_leading_minors_bigint, hurwitz_matrix_bigint, is_hurwitz_stable_bigint};
pub use interlacing_matrix::{
    check_lace_sequence_total_nonnegative_i64, check_lace_tnn_neville_bigint,
    check_lace_tnn_neville_i64, check_lace_total_nonnegative_i64,
    is_lace_sequence_totally_nonnegative_i64, is_lace_totally_nonnegative_i64, lace_matrix,
    lace_matrix_bigint, lace_matrix_i64, lace_matrix_sequence, lace_matrix_sequence_bigint,
    lace_matrix_sequence_i64, InterlacingMatrixError,
};
pub use parse::{
    parse_polynomial, parse_polynomial_bigint, parse_polynomials, parse_polynomials_bigint,
    MAX_POLYNOMIAL_COEFFICIENTS, MAX_POLYNOMIAL_INPUT_BYTES,
};
pub use polynomial::{CoeffRing, FieldRing, Polynomial};
pub use real_rooted_enumeration::{
    enumerate_real_rooted_integer_polynomials_bigint,
    real_rooted_antiderivative_constant_multiples_bigint,
    real_rooted_antiderivative_integer_constants_bigint, AntiderivativeIntervalError,
    InclusiveIntegerRange, RealRootedEnumerationError, RealRootedEnumerationLimit,
    RealRootedEnumerationOptions, RealRootedEnumerationResult, RealRootedEnumerationStats,
    RealRootedEnumerationStatus,
};
pub use real_rootedness::{
    // Bézout matrix directly
    bezout_matrix,
    bezout_matrix_bigint_coeffs,
    bezoutian_bigint_coeffs,
    check_interlacing,
    check_interlacing_bigint_coeffs,
    check_interlacing_sturm,
    check_weak_interlacing,
    check_weak_interlacing_bigint_coeffs,
    discriminant,
    // Ehrhart polynomial <-> h*-vector
    discriminant_bigint_coeffs,
    ehrhart_to_hstar,
    ehrhart_to_hstar_bigint,
    ehrhart_to_hstar_with_denom,
    // Display
    format_poly,
    format_poly_bigint_coeffs,
    format_poly_bigint_var,
    format_poly_var,
    gamma_coefficients,
    gamma_coefficients_bigint_coeffs,
    gamma_coefficients_ignoring_initial_zeros,
    gamma_coefficients_ignoring_initial_zeros_bigint_coeffs,
    // Concavity and symmetry
    has_simple_roots,
    has_simple_roots_bigint_coeffs,
    hermite_biehler_parts,
    hermite_matrix_bigint_coeffs,
    hstar_to_ehrhart,
    hstar_to_ehrhart_bigint_coeffs,
    // Undirected interlacing
    interlaces_either_order,
    interlaces_either_order_bigint_coeffs,
    is_gamma_positive,
    is_gamma_positive_bigint_coeffs,
    is_gamma_positive_ignoring_initial_zeros,
    is_gamma_positive_ignoring_initial_zeros_bigint_coeffs,
    is_log_concave,
    is_log_concave_bigint_coeffs,
    is_palindromic,
    is_palindromic_bigint_coeffs,
    is_palindromic_ignoring_initial_zeros,
    is_palindromic_ignoring_initial_zeros_bigint_coeffs,
    // Default methods and explicit Bézout path
    is_real_rooted,
    is_real_rooted_bezout_bigint_coeffs,
    is_real_rooted_bezout_squarefree,
    is_real_rooted_bezout_squarefree_bigint_coeffs,
    is_real_rooted_bigint_coeffs,
    is_real_rooted_hermite,
    is_real_rooted_hermite_bigint_coeffs,
    // Sturm-chain methods (slower, but can isolate roots)
    is_real_rooted_sturm,
    is_real_rooted_sturm_bigint_coeffs,
    is_strictly_real_rooted_bezout,
    is_strictly_real_rooted_bezout_bigint_coeffs,
    is_strictly_real_rooted_hermite,
    is_strictly_real_rooted_hermite_bigint_coeffs,
    is_ultra_log_concave,
    is_ultra_log_concave_bigint_coeffs,
    is_ultra_log_concave_infinite,
    is_ultra_log_concave_infinite_bigint_coeffs,
    is_ultra_log_concave_of_degree,
    is_ultra_log_concave_of_degree_bigint_coeffs,
    is_unimodal,
    is_unimodal_bigint_coeffs,
    real_roots,
    // Resultant and discriminant
    resultant,
    resultant_bigint_coeffs,
    squarefree_part_bigint_coeffs,
    stapledon_decomposition,
    stapledon_decomposition_bigint_coeffs,
    strip_initial_zeros,
    strip_initial_zeros_bigint_coeffs,
    sylvester_matrix,
};
pub use real_rootedness::{
    check_interlacing_bezout, check_weak_interlacing_bezout, is_real_rooted_bezout,
};
pub use root_count::{
    bigint_coeffs_to_i64, count_positive_roots_prs_bigint_coeffs,
    count_positive_roots_uspensky_bigint_coeffs, count_real_roots_prs_bigint_coeffs,
    count_real_roots_uspensky_bigint_coeffs, count_root_signs_prs_bigint_coeffs,
    is_real_rooted_fast_bigint_coeffs, is_real_rooted_fast_i64,
    is_real_rooted_one_signed_bigint_coeffs, is_real_rooted_prs_bigint_coeffs,
    is_real_rooted_uspensky_bigint_coeffs, is_real_rooted_uspensky_i64,
    primitive_sturm_max_coefficient_bits, satisfies_kurtz_condition_bigint,
    satisfies_newton_inequalities_bigint, squarefree_degree_bigint_coeffs,
    tarski_query_prs_bigint_coeffs, RootSignCounts,
};
pub use sturm_cf::{
    check_sturm_continued_fraction_certificate,
    check_sturm_continued_fraction_certificate_with_options, SturmContinuedFractionError,
    SturmContinuedFractionOptions, SturmContinuedFractionSummary,
};
pub use tnn_network::{
    build_tnn_certificate_from_monic_polynomials, chain_polynomials_bigint,
    chain_polynomials_rational, check_chain_polynomials_bigint, check_lower_unitriangular_tnn,
    coefficient_matrix_from_monic_polynomials, evaluate_path_matrix,
    lower_triangular_rows_from_square, reconstruct_canonical_tnn_network,
    verify_path_matrix_certificate, BigRational, CanonicalPlanarNetwork, CanonicalTnnProof,
    ChainPolynomialReport, NegativeMinor, NetworkError,
};
