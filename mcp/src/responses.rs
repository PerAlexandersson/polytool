//! Tool response types.  Exact integers and rationals are serialized as strings.

use crate::requests::{EhrhartHstarMode, SequenceKind};
use schemars::JsonSchema;
use serde::Serialize;

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct NormalizedPolynomial {
    pub polynomial: String,
    pub coefficients: Vec<i64>,
    pub degree: usize,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct BigIntNormalizedPolynomial {
    pub polynomial: String,
    pub coefficients: Vec<String>,
    pub degree: usize,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct OeisSequenceSummary {
    pub id: String,
    pub name: String,
    pub status: String,
    pub layout: String,
    pub first_row: i64,
    pub flattened_offset: i64,
    pub bfile_available: bool,
    pub source_rows: usize,
    pub verification_rows: usize,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ListOeisSequencesResponse {
    pub count: usize,
    pub sequences: Vec<OeisSequenceSummary>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct GetOeisSequenceResponse {
    pub sequence: OeisSequenceSummary,
    pub recurrence: String,
    pub latex: String,
    pub mathematica: String,
    pub sage: String,
    pub python: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct OeisPolynomialRow {
    pub n: i64,
    pub polynomial: String,
    pub coefficients: Vec<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct GenerateOeisRowsResponse {
    pub id: String,
    pub first_row: i64,
    pub row_count: usize,
    pub rows: Vec<OeisPolynomialRow>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ParsePolynomialItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coefficients: Option<Vec<i64>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct BigIntParsePolynomialItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<BigIntNormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ParsePolynomialsResponse {
    pub items: Vec<ParsePolynomialItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PolynomialPropertiesItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<BigIntNormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub simple_roots: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub palindromic: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_positive: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_coefficients: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unimodal: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub log_concave: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ultra_log_concave: Option<bool>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PolynomialPropertiesResponse {
    pub items: Vec<PolynomialPropertiesItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingPairResult {
    pub p: BigIntNormalizedPolynomial,
    pub q: BigIntNormalizedPolynomial,
    pub strict: Option<bool>,
    pub weak: Option<bool>,
    pub status: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingPairItem {
    pub pair_index: usize,
    pub left_index: usize,
    pub right_index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub result: Option<InterlacingPairResult>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingSequenceResponse {
    pub items: Vec<BigIntParsePolynomialItem>,
    pub pairs: Vec<InterlacingPairItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingProfileItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<BigIntNormalizedPolynomial>,
    pub previous_count: usize,
    pub checked_previous_count: usize,
    pub interlacing_previous_count: usize,
    pub strict_previous_count: usize,
    pub weak_previous_count: usize,
    pub previous_interlacing_indices: Vec<usize>,
    pub previous: Vec<InterlacingPairItem>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct InterlacingProfileResponse {
    pub items: Vec<BigIntParsePolynomialItem>,
    pub profile: Vec<InterlacingProfileItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RealRootsItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<NormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub roots: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RealRootsResponse {
    pub items: Vec<RealRootsItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RecurrenceSearchStatus {
    Found,
    NotFound,
    BudgetExhausted,
    InvalidInput,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct FindRecurrenceResponse {
    pub status: RecurrenceSearchStatus,
    pub found: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recurrence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub latex: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub mathematica: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub sage: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub python: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recurrence_json: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub weighted_unknowns: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub equations: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fit_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verification_polynomials: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub candidates_tried: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub candidates_considered: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub parse_errors: Vec<ParsePolynomialItem>,
    /// Conditional theorem matches of a found recurrence (objects of
    /// `polytool.recurrence-theorems.v1`); omitted when none match.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub theorem_matches: Option<Vec<serde_json::Value>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct GenerateRecurrenceRowsResponse {
    pub first_index: usize,
    pub row_count: usize,
    pub polynomials: Vec<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ResultantResponse {
    pub p: BigIntNormalizedPolynomial,
    pub q: BigIntNormalizedPolynomial,
    pub resultant: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DiscriminantItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<BigIntNormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub discriminant: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DiscriminantResponse {
    pub items: Vec<DiscriminantItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct EhrhartHstarResponse {
    pub mode: EhrhartHstarMode,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hstar: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ehrhart_coefficients: Option<Vec<String>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct HStarInequalityCheckItem {
    pub family: String,
    pub name: String,
    pub formula: String,
    pub reference: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index: Option<usize>,
    pub applicable: bool,
    pub holds: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lhs: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub rhs: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub value: Option<String>,
    pub details: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct HStarInequalitiesItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub hstar: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub dimension: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub codegree: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub all_applicable_hold: Option<bool>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub checks: Vec<HStarInequalityCheckItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct HStarInequalitiesResponse {
    pub items: Vec<HStarInequalitiesItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CoefficientInequalityCheckItem {
    pub index: usize,
    pub lhs: String,
    pub rhs: String,
    pub comparison: String,
    pub holds: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CoefficientCriterionItem {
    pub name: String,
    pub reference: String,
    pub applicable: bool,
    pub holds: bool,
    pub implies_real_rooted: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    pub checks: Vec<CoefficientInequalityCheckItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CoefficientTestsItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub polynomial: Option<BigIntNormalizedPolynomial>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub newton: Option<CoefficientCriterionItem>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kurtz: Option<CoefficientCriterionItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CoefficientTestsResponse {
    pub items: Vec<CoefficientTestsItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RootOfUnityEvaluationItem {
    pub group_order: usize,
    pub power: usize,
    pub root_order: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub integer_value: Option<String>,
    pub remainder: Vec<String>,
    pub remainder_polynomial: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CyclicSievingPowerCheckItem {
    pub power: usize,
    pub root_order: usize,
    pub expected_fixed_points: String,
    pub evaluation: RootOfUnityEvaluationItem,
    pub holds: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CyclicSievingReportItem {
    pub order: usize,
    pub polynomial: BigIntNormalizedPolynomial,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_counts: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub holds: Option<bool>,
    pub evaluations: Vec<RootOfUnityEvaluationItem>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub checks: Vec<CyclicSievingPowerCheckItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ResidueClassItem {
    pub gcd_class: usize,
    pub indices: Vec<usize>,
    pub residues: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub common_value: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct OrbitSizeEntryItem {
    pub orbit_size: usize,
    /// `r = n / orbit_size`, the index of the orbit-basis polynomial `h_r`.
    pub basis_index: usize,
    /// Number of orbits of this size (negative: not realizable).
    pub orbit_count: String,
    /// Number of elements in orbits of this size, `orbit_size * orbit_count`.
    pub element_count: String,
}

/// Orbit decomposition (Alexandersson-Amini, arXiv:1804.01447).
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct OrbitDecompositionItem {
    pub order: usize,
    /// Coefficients of `f mod q^n - 1`.
    pub reduced_residues: Vec<String>,
    pub classes: Vec<ResidueClassItem>,
    /// `f(zeta^d)` is an integer for every `d` (iff `c_j = c_gcd(j,n)`).
    pub integral: bool,
    pub violating_classes: Vec<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orbits: Option<Vec<OrbitSizeEntryItem>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_points: Option<Vec<String>>,
    /// Hypothesis of the theorem: all input coefficients are nonnegative.
    pub natural_coefficients: bool,
    /// All orbit counts are nonnegative (all `S_k >= 0`); present iff integral.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cone_condition: Option<bool>,
    /// Some cyclic action has the fixed-point counts `f(zeta^d)`.
    pub realizable: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct OrbitSizesCheckItem {
    pub orbit_sizes: Vec<usize>,
    pub num_elements: String,
    pub actual_orbits: Vec<OrbitSizeEntryItem>,
    pub actual_fixed_points: Vec<String>,
    pub mismatched_orbit_sizes: Vec<usize>,
    pub holds: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CyclicSievingResponse {
    pub report: CyclicSievingReportItem,
    pub orbit_decomposition: OrbitDecompositionItem,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub orbit_sizes_check: Option<OrbitSizesCheckItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct BicyclicSievingResponse {
    pub order1: usize,
    pub order2: usize,
    pub reduced_matrix: Vec<Vec<String>>,
    /// `f(zeta1^d1, zeta2^d2)`; `null` where the value is not an integer.
    pub evaluations: Vec<Vec<Option<String>>>,
    pub integral: bool,
    pub natural_coefficients: bool,
    pub divisors1: Vec<usize>,
    pub divisors2: Vec<usize>,
    /// `c_{j1,j2} = c_{gcd(j1,n1),gcd(j2,n2)}` for all `j1, j2`.
    pub product_basis_form: bool,
    /// Coordinates `a_{r1,r2}` in the basis `h_{r1}(q) h_{r2}(t)` (rows `r1` in
    /// `divisors1`, columns `r2` in `divisors2`).  They are NOT orbit counts and
    /// may be negative (e.g. `1 + q t` for `C_2 x C_2` gives 1, -1, -1, 2).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub product_coordinates: Option<Vec<Vec<String>>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub product_formula_evaluations: Option<Vec<Vec<String>>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fixed_table: Option<Vec<Vec<String>>>,
    pub mismatches: Vec<[usize; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub holds: Option<bool>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CyclicSievingSequenceItemResponse {
    pub row: usize,
    pub index: isize,
    pub polynomial: BigIntNormalizedPolynomial,
    pub candidate_orders: Vec<CyclicSievingReportItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CyclicSievingSequenceResponse {
    pub first_index: isize,
    pub offsets: Vec<isize>,
    pub items: Vec<CyclicSievingSequenceItemResponse>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DisplayedPolynomial {
    pub polynomial: String,
    pub coefficients: Vec<i64>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct MagicBasisReport {
    pub coordinates: Vec<String>,
    pub left_partial_sums: Vec<String>,
    pub right_partial_sums: Vec<String>,
    pub partial_sum_checks: Vec<bool>,
    pub all_nonnegative: bool,
    pub left_leq_right: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionReport {
    #[serde(flatten)]
    pub polynomial: NormalizedPolynomial,
    pub reciprocal: DisplayedPolynomial,
    pub a: DisplayedPolynomial,
    pub b: DisplayedPolynomial,
    pub a_real_rooted: bool,
    pub b_real_rooted: bool,
    pub b_interlaces_a: Option<bool>,
    pub reciprocal_interlaces_input: Option<bool>,
    pub alternatingly_increasing: bool,
    pub f_polynomial: DisplayedPolynomial,
    pub r_transform_of_f: DisplayedPolynomial,
    pub r_a: DisplayedPolynomial,
    pub r_b: DisplayedPolynomial,
    pub r_interlaces_f: Option<bool>,
    pub magic: MagicBasisReport,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionItem {
    pub index: usize,
    pub ok: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub report: Option<DecompositionReport>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct DecompositionResponse {
    pub items: Vec<DecompositionItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct GenerateSequenceResponse {
    pub sequence: SequenceKind,
    pub max_n: usize,
    pub polynomials: Vec<BigIntNormalizedPolynomial>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct FamilyPolynomialReport {
    pub index: usize,
    #[serde(flatten)]
    pub polynomial: NormalizedPolynomial,
    pub real_rooted: bool,
    pub simple_roots: bool,
    pub palindromic: bool,
    pub gamma_positive: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gamma_coefficients: Option<Vec<i64>>,
    pub unimodal: bool,
    pub log_concave: bool,
    pub ultra_log_concave: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct LaceCheckResponse {
    pub block_rows: usize,
    pub block_cols: usize,
    pub max_minor_size: usize,
    pub rows: usize,
    pub columns: usize,
    pub tnn: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub matrix: Option<Vec<Vec<i64>>>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct CheckPolynomialFamilyResponse {
    pub source: String,
    pub item_count: usize,
    pub all_required_checks_passed: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_failure: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub parse_errors: Vec<ParsePolynomialItem>,
    pub items: Vec<FamilyPolynomialReport>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub consecutive_pairs: Vec<InterlacingPairItem>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lace: Option<LaceCheckResponse>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub recurrence: Option<FindRecurrenceResponse>,
    pub markdown: String,
}

/// An explicit negative minor (zero-based, strictly increasing indices into
/// the zero-padded square matrix).
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct NegativeMinorItem {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    /// Exact (negative) value of the minor.
    pub value: String,
}

/// Why a matrix is not lower-unitriangular TNN.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct TnnFailureItem {
    /// Snake-case `NetworkError` variant, e.g. `not_tnn`, `negative_entry`,
    /// `non_monic_polynomial`, `nonzero_above_diagonal`.
    pub kind: String,
    pub message: String,
    /// Zero-based row (or sequence position) of the offending entry.
    pub row: Option<usize>,
    /// Zero-based column of the offending entry.
    pub col: Option<usize>,
    /// The offending value (exact integer or rational string).
    pub value: Option<String>,
    /// For `not_tnn`: an explicit negative minor.
    pub witness: Option<NegativeMinorItem>,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ChainPolynomialItem {
    pub n: usize,
    pub polynomial: String,
    pub coefficients: Vec<String>,
    pub real_rooted: bool,
    /// Real-rooted with every zero in `[-1, 0]`.
    pub roots_in_minus_one_zero: bool,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ChainInterlacingItem {
    pub left: usize,
    pub right: usize,
    /// `p_left` interlaces `p_right` (zeros of `p_left` weakly to the left).
    pub interlaces: bool,
}

/// Chain polynomials and the checks of Brändén--Saud Maia Leite,
/// arXiv:2412.06595, Theorem 3.7.  These concern the chain polynomials, not
/// the row polynomials of the input.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ChainPolynomialsItem {
    pub theorem: String,
    /// Set when the matrix shape does not define chain polynomials.
    pub error: Option<String>,
    pub polynomials: Vec<ChainPolynomialItem>,
    pub interlacing: Vec<ChainInterlacingItem>,
    /// Every `p_n` is real-rooted in `[-1, 0]` and consecutive ones interlace.
    pub conclusion_holds: Option<bool>,
    /// The matrix is TNN but the conclusion fails (would contradict the
    /// theorem).
    pub contradicts_theorem: Option<bool>,
}

/// Result of `check_tnn_network`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct TnnNetworkResponse {
    /// `polynomials` or `matrix`.
    pub input_mode: String,
    pub input_rows: usize,
    /// `P_0 = 1` was prepended to the polynomial rows.
    pub prepended_p0: bool,
    /// Triangular rows `r_{n,k}` (`k <= n`) as exact integer strings.
    pub coefficient_matrix: Option<Vec<Vec<String>>>,
    /// The coefficient matrix is lower-unitriangular and totally nonnegative.
    pub tnn: bool,
    pub failure: Option<TnnFailureItem>,
    /// Canonical network weights: `lambda[n][k]` is the weight of the
    /// vertical edge `(n+1,k) -> (n,k)` of `Gamma_N(lambda)` (exact rationals).
    pub lambda: Option<Vec<Vec<String>>>,
    /// The network's path matrix re-evaluates exactly to the input matrix.
    pub certificate_verified: Option<bool>,
    pub chain: Option<ChainPolynomialsItem>,
    /// Fixed guardrail: TNN of the matrix says nothing about the rows.
    pub note: String,
}

#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct StripDigraphEdgeItem {
    /// Target height in the next column.
    pub to: usize,
    /// Exact nonnegative weight.
    pub weight: String,
}

/// A planar strip digraph: `outgoing[y]` lists the edges leaving height `y`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct StripDigraphItem {
    pub num_heights: usize,
    pub outgoing: Vec<Vec<StripDigraphEdgeItem>>,
}

/// The options actually used by `recognize_brenti_strip_digraph`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct BrentiRecognitionOptionsItem {
    pub first_row_index: usize,
    pub max_target_span: Option<usize>,
    pub weights: String,
    pub max_search_nodes: usize,
    pub max_integer_search_nodes: usize,
}

/// Result of `recognize_brenti_strip_digraph`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct BrentiRecognitionResponse {
    pub options: BrentiRecognitionOptionsItem,
    pub rows: usize,
    /// `certified`, `impossible_within_bounds` or `inconclusive`.
    pub verdict: String,
    pub certified: bool,
    pub digraph: Option<StripDigraphItem>,
    /// Digraph row of each input polynomial.
    pub row_indices: Option<Vec<usize>>,
    /// The certificate was re-checked against every input row.
    pub certificate_verified: Option<bool>,
    pub supports_checked: Option<usize>,
    pub obstruction: Option<String>,
    pub reason: Option<String>,
    pub message: String,
    /// Fixed guardrail: a certificate proves TNN of the coefficient matrix,
    /// not real-rootedness of the rows.
    pub note: String,
}

/// One step of a height-invariant network: for every height `y`, an edge
/// `(n - length, y) -> (n, y + rise)` with this weight.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PfNetworkStepItem {
    pub length: usize,
    pub rise: usize,
    /// Exact nonnegative weight.
    pub weight: String,
}

/// A height-invariant planar network: `steps[n - 1]` lists the steps into
/// column `n`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PfNetworkItem {
    pub last_column: usize,
    pub steps: Vec<Vec<PfNetworkStepItem>>,
}

/// The options actually used by `recognize_pf_network`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PfNetworkOptionsItem {
    pub first_row_index: usize,
    pub max_step_length: usize,
    pub max_rise: usize,
    pub weights: String,
    pub max_search_nodes: usize,
    pub max_integer_search_nodes: usize,
}

/// Result of `recognize_pf_network`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct PfNetworkResponse {
    pub options: PfNetworkOptionsItem,
    pub rows: usize,
    /// `certified`, `impossible_within_bounds` or `inconclusive`.
    pub verdict: String,
    /// Every supplied row is certified PF (real-rooted, nonpositive zeros).
    pub certified: bool,
    pub network: Option<PfNetworkItem>,
    /// Network row of each input polynomial.
    pub row_indices: Option<Vec<usize>>,
    /// The certificate was re-checked (planarity, rows, Toeplitz windows).
    pub certificate_verified: Option<bool>,
    pub supports_checked: Option<usize>,
    pub obstruction: Option<String>,
    pub reason: Option<String>,
    pub message: String,
    /// Fixed guardrail: a certificate proves PF rows, not TNN of the
    /// coefficient matrix.
    pub note: String,
}

/// A negative minor of a rational matrix (zero-based rows and columns).
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct RationalNegativeMinorItem {
    pub rows: Vec<usize>,
    pub cols: Vec<usize>,
    /// Exact negative value.
    pub value: String,
}

/// A planar network for a TNN production matrix: `Q_N = L D`, with
/// `lambda` the canonical network of `L` and `D = diag(column_scales)`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct ProductionNetworkItem {
    /// `lambda[n][k]`: weight of the vertical edge `(n+1,k) -> (n,k)`.
    pub lambda: Vec<Vec<String>>,
    /// The diagonal of `Q_N`.
    pub column_scales: Vec<String>,
}

/// Result of `check_left_production_matrix`.
#[derive(Debug, Clone, Serialize, JsonSchema, PartialEq, Eq)]
pub struct LeftProductionResponse {
    pub first_row_index: usize,
    pub rows: usize,
    /// The triangle `A_N` (triangular rows, exact integer strings).
    pub triangle: Vec<Vec<String>>,
    /// The left production matrix `Q_N` (triangular rows, exact rationals).
    pub production_matrix: Vec<Vec<String>>,
    /// `Q_N` is totally nonnegative.
    pub tnn: bool,
    pub failure_reason: Option<String>,
    pub negative_minor: Option<RationalNegativeMinorItem>,
    /// When TNN: the certificate re-verified (regeneration, TNN, Toeplitz
    /// windows of the unified network).
    pub certificate_verified: Option<bool>,
    /// Rows `P_0..P_N` are certified PF and the triangle and its reversal TNN.
    pub rows_certified_pf: bool,
    pub production_network: Option<ProductionNetworkItem>,
    pub message: String,
    /// Fixed guardrail: the criterion is sufficient, not necessary.
    pub note: String,
}
