use num_bigint::BigInt;
use polytool::recurrence::{
    candidate_complexity, find_recurrence_adaptive_rational_with_budget,
    find_recurrence_adaptive_with_budget, format_rational_coeff, parse_rational_coeff,
    AdaptiveSearchBudget, AdaptiveSearchOptions, AdaptiveSearchOutcome, AdaptiveSearchResult,
    AdaptiveSearchSummary, BigRational, BivarPoly, Recurrence, RecurrenceJson,
    RecurrenceJsonSearch, RecurrenceOptions, RecurrenceOptionsJson,
};
use polytool::sequences::{
    chebyshev_polynomials_t_bigint, chebyshev_polynomials_u_bigint, eulerian_polynomials_bigint,
    hermite_polynomials_bigint, narayana_polynomials_bigint, type_b_eulerian_polynomials_bigint,
};
use polytool::*;
use rmcp::{
    handler::server::{router::tool::ToolRouter, wrapper::Parameters},
    model::{Implementation, ServerCapabilities, ServerInfo},
    tool, tool_handler, tool_router, ErrorData as McpError, Json, ServerHandler,
};
use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::borrow::Cow;
use std::collections::BTreeMap;

const MCP_MAX_TEXT_BYTES: usize = 1_048_576;
const MCP_MAX_BATCH_POLYNOMIALS: usize = 512;
const MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL: usize = 4_096;
const MCP_MAX_TOTAL_COEFFICIENTS: usize = 65_536;
const MCP_MAX_COEFFICIENT_TEXT_BYTES: usize = 16_384;
const MCP_MAX_SEQUENCE_N: usize = 200;
const MCP_MAX_RECURRENCE_ROWS: usize = 200;
const MCP_MAX_RECURRENCE_INPUT_ROWS: usize = 256;
const MCP_MAX_RECURRENCE_DEGREE: usize = 16;
const MCP_MAX_RECURRENCE_LENGTH: usize = 32;
const MCP_MAX_RECURRENCE_CANDIDATES: usize = 50_000;
const MCP_MAX_RECURRENCE_UNKNOWNS: usize = 20_000;
const MCP_MAX_LACE_MATRIX_CELLS: usize = 16_384;
const MCP_MAX_LACE_MINORS: usize = 250_000;

#[derive(Debug, Clone)]
pub struct PolynomialToolsServer {
    tool_router: ToolRouter<Self>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialInput {
    pub coefficients: Option<Vec<i64>>,
    pub expression: Option<String>,
}

impl JsonSchema for PolynomialInput {
    fn schema_name() -> Cow<'static, str> {
        "PolynomialInput".into()
    }

    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "coefficients": <Vec<i64>>::json_schema(generator),
                "expression": String::json_schema(generator)
            },
            "oneOf": [
                { "required": ["coefficients"] },
                { "required": ["expression"] }
            ]
        }))
        .expect("valid PolynomialInput schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(untagged)]
pub enum BigIntCoefficientInput {
    Integer(i64),
    Text(String),
}

impl From<i64> for BigIntCoefficientInput {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}

impl JsonSchema for BigIntCoefficientInput {
    fn schema_name() -> Cow<'static, str> {
        "BigIntCoefficientInput".into()
    }

    fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "anyOf": [
                { "type": "integer" },
                {
                    "type": "string",
                    "description": "Exact integer coefficient, e.g. \"42\" or \"1267650600228229401496703205376\"."
                }
            ]
        }))
        .expect("valid BigIntCoefficientInput schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct BigIntPolynomialInput {
    pub coefficients: Option<Vec<BigIntCoefficientInput>>,
    pub expression: Option<String>,
}

impl JsonSchema for BigIntPolynomialInput {
    fn schema_name() -> Cow<'static, str> {
        "BigIntPolynomialInput".into()
    }

    fn json_schema(generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "type": "object",
            "additionalProperties": false,
            "properties": {
                "coefficients": <Vec<BigIntCoefficientInput>>::json_schema(generator),
                "expression": String::json_schema(generator)
            },
            "oneOf": [
                { "required": ["coefficients"] },
                { "required": ["expression"] }
            ]
        }))
        .expect("valid BigIntPolynomialInput schema")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BigIntPolynomialBatchInput {
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    pub text: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PolynomialBatchInput {
    pub polynomials: Option<Vec<PolynomialInput>>,
    pub text: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FamilyCheckOptions {
    pub require_real_rooted: Option<bool>,
    pub require_simple_roots: Option<bool>,
    pub require_palindromic: Option<bool>,
    pub require_gamma_positive: Option<bool>,
    pub require_unimodal: Option<bool>,
    pub require_log_concave: Option<bool>,
    pub require_ultra_log_concave: Option<bool>,
    pub check_consecutive_interlacing: Option<bool>,
    pub require_consecutive_weak_interlacing: Option<bool>,
    pub find_recurrence: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LaceCheckRequest {
    pub block_rows: Option<usize>,
    pub block_cols: Option<usize>,
    pub max_minor_size: Option<usize>,
    pub include_matrix: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CheckPolynomialFamilyRequest {
    pub polynomials: Option<Vec<PolynomialInput>>,
    pub text: Option<String>,
    pub sequence: Option<SequenceKind>,
    pub max_n: Option<usize>,
    pub options: Option<FamilyCheckOptions>,
    pub recurrence_options: Option<RecurrenceSearchOptionsInput>,
    pub lace: Option<LaceCheckRequest>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BigIntPolynomialPairRequest {
    pub p: BigIntPolynomialInput,
    pub q: BigIntPolynomialInput,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSearchOptionsInput {
    /// Ignore this many leading polynomial rows before searching.
    pub skip_prefix: Option<usize>,
    /// Minimum recurrence length (largest lag).
    pub min_rec_len: Option<usize>,
    /// Maximum recurrence length (largest lag).
    pub max_rec_len: Option<usize>,
    /// Minimum degree in the polynomial variable.
    pub min_var_deg: Option<usize>,
    /// Maximum degree in the polynomial variable.
    pub max_var_deg: Option<usize>,
    /// Minimum degree in the recurrence index.
    pub min_idx_deg: Option<usize>,
    /// Maximum degree in the recurrence index.
    pub max_idx_deg: Option<usize>,
    /// Minimum derivative order allowed in recurrence terms.
    pub min_diff_deg: Option<usize>,
    /// Maximum derivative order allowed in recurrence terms.
    pub max_diff_deg: Option<usize>,
    /// Search for an additive inhomogeneous term.
    pub try_inhomogeneous: Option<bool>,
    /// Minimum polynomial-variable degree of the inhomogeneous term.
    pub min_inhomo_var_deg: Option<usize>,
    /// Maximum polynomial-variable degree of the inhomogeneous term.
    pub max_inhomo_var_deg: Option<usize>,
    /// Minimum recurrence-index degree of the inhomogeneous term.
    pub min_inhomo_idx_deg: Option<usize>,
    /// Maximum recurrence-index degree of the inhomogeneous term.
    pub max_inhomo_idx_deg: Option<usize>,
    /// Search for a nonconstant denominator on the left-hand side.
    pub try_denominator: Option<bool>,
    /// Search recurrence terms multiplied by `(-1)^n`.
    pub try_alternating_sign: Option<bool>,
    /// Maximum polynomial-variable degree of the denominator.
    pub max_denom_var_deg: Option<usize>,
    /// Maximum recurrence-index degree of the denominator.
    pub max_denom_idx_deg: Option<usize>,
    /// Minimum excess of equations over unknowns.
    pub min_margin: Option<usize>,
    /// Disable holdout verification of fitted recurrences.
    pub no_verify: Option<bool>,
    /// Add this many rows to the fit before holdout verification.
    pub fit_extra_rows: Option<usize>,
    /// Use the modular inconsistency prefilter (default `true`).
    pub modular_prefilter: Option<bool>,
    /// Inspect at most this many adaptive candidate configurations. Zero
    /// returns budget exhaustion before the first candidate.
    pub max_candidates: Option<usize>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(untagged)]
pub enum RationalCoefficientInput {
    Integer(i64),
    Text(String),
}

impl From<i64> for RationalCoefficientInput {
    fn from(value: i64) -> Self {
        Self::Integer(value)
    }
}

impl JsonSchema for RationalCoefficientInput {
    fn schema_name() -> Cow<'static, str> {
        "RationalCoefficientInput".into()
    }

    fn json_schema(_generator: &mut SchemaGenerator) -> Schema {
        Schema::try_from(json!({
            "anyOf": [
                { "type": "integer" },
                {
                    "type": "string",
                    "description": "Exact rational coefficient, e.g. \"42\", \"-17\", or \"3/7\"."
                }
            ]
        }))
        .expect("valid RationalCoefficientInput schema")
    }
}

#[derive(Debug, Clone, Default, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FindRecurrenceRequest {
    /// Polynomial objects; provide exactly one polynomial input form.
    pub polynomials: Option<Vec<PolynomialInput>>,
    /// Dense exact rational coefficient rows; integers or rational strings.
    pub coefficients: Option<Vec<Vec<RationalCoefficientInput>>>,
    /// Polynomial expressions, one expression per sequence row.
    pub expressions: Option<Vec<String>>,
    /// Newline-separated polynomial expressions or coefficient rows.
    pub text: Option<String>,
    /// Legacy nested search controls. Top-level controls take precedence.
    pub options: Option<RecurrenceSearchOptionsInput>,
    /// Omit generated Mathematica, Python, Sage, and recurrence JSON when false.
    /// Defaults to true for backward compatibility.
    pub include_code: Option<bool>,
    /// Ignore this many leading polynomial rows before searching.
    pub skip_prefix: Option<usize>,
    /// Minimum recurrence length (largest lag).
    pub min_rec_len: Option<usize>,
    /// Maximum recurrence length (largest lag).
    pub max_rec_len: Option<usize>,
    /// Minimum degree in the polynomial variable.
    pub min_var_deg: Option<usize>,
    /// Maximum degree in the polynomial variable.
    pub max_var_deg: Option<usize>,
    /// Minimum degree in the recurrence index.
    pub min_idx_deg: Option<usize>,
    /// Maximum degree in the recurrence index.
    pub max_idx_deg: Option<usize>,
    /// Minimum derivative order allowed in recurrence terms.
    pub min_diff_deg: Option<usize>,
    /// Maximum derivative order allowed in recurrence terms.
    pub max_diff_deg: Option<usize>,
    /// Search for an additive inhomogeneous term.
    pub try_inhomogeneous: Option<bool>,
    /// Minimum polynomial-variable degree of the inhomogeneous term.
    pub min_inhomo_var_deg: Option<usize>,
    /// Maximum polynomial-variable degree of the inhomogeneous term.
    pub max_inhomo_var_deg: Option<usize>,
    /// Minimum recurrence-index degree of the inhomogeneous term.
    pub min_inhomo_idx_deg: Option<usize>,
    /// Maximum recurrence-index degree of the inhomogeneous term.
    pub max_inhomo_idx_deg: Option<usize>,
    /// Search for a nonconstant denominator on the left-hand side.
    pub try_denominator: Option<bool>,
    /// Search recurrence terms multiplied by `(-1)^n`.
    pub try_alternating_sign: Option<bool>,
    /// Maximum polynomial-variable degree of the denominator.
    pub max_denom_var_deg: Option<usize>,
    /// Maximum recurrence-index degree of the denominator.
    pub max_denom_idx_deg: Option<usize>,
    /// Minimum excess of equations over unknowns.
    pub min_margin: Option<usize>,
    /// Disable holdout verification of fitted recurrences.
    pub no_verify: Option<bool>,
    /// Add this many rows to the fit before holdout verification.
    pub fit_extra_rows: Option<usize>,
    /// Use the modular inconsistency prefilter (default `true`).
    pub modular_prefilter: Option<bool>,
    /// Inspect at most this many adaptive candidate configurations.
    pub max_candidates: Option<usize>,
}

impl FindRecurrenceRequest {
    fn effective_options(&self) -> Option<RecurrenceSearchOptionsInput> {
        let mut options = self.options.clone().unwrap_or_default();
        let mut present = self.options.is_some();
        macro_rules! prefer_top_level {
            ($field:ident) => {
                if let Some(value) = self.$field {
                    options.$field = Some(value);
                    present = true;
                }
            };
        }
        prefer_top_level!(skip_prefix);
        prefer_top_level!(min_rec_len);
        prefer_top_level!(max_rec_len);
        prefer_top_level!(min_var_deg);
        prefer_top_level!(max_var_deg);
        prefer_top_level!(min_idx_deg);
        prefer_top_level!(max_idx_deg);
        prefer_top_level!(min_diff_deg);
        prefer_top_level!(max_diff_deg);
        prefer_top_level!(try_inhomogeneous);
        prefer_top_level!(min_inhomo_var_deg);
        prefer_top_level!(max_inhomo_var_deg);
        prefer_top_level!(min_inhomo_idx_deg);
        prefer_top_level!(max_inhomo_idx_deg);
        prefer_top_level!(try_denominator);
        prefer_top_level!(try_alternating_sign);
        prefer_top_level!(max_denom_var_deg);
        prefer_top_level!(max_denom_idx_deg);
        prefer_top_level!(min_margin);
        prefer_top_level!(no_verify);
        prefer_top_level!(fit_extra_rows);
        prefer_top_level!(modular_prefilter);
        prefer_top_level!(max_candidates);
        present.then_some(options)
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EhrhartHstarRequest {
    pub mode: EhrhartHstarMode,
    pub hstar: Option<Vec<BigIntCoefficientInput>>,
    pub ehrhart_coefficients: Option<Vec<String>>,
    pub numerator_coefficients: Option<Vec<BigIntCoefficientInput>>,
    pub denominator: Option<BigIntCoefficientInput>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct HStarInequalitiesRequest {
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    pub text: Option<String>,
    pub dimension: Option<usize>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CyclicSievingRequest {
    pub polynomial: BigIntPolynomialInput,
    pub order: usize,
    pub fixed_counts: Option<Vec<BigIntCoefficientInput>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CyclicSievingFixedCountsInput {
    pub index: isize,
    pub order: usize,
    pub counts: Vec<BigIntCoefficientInput>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CyclicSievingSequenceRequest {
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    pub text: Option<String>,
    pub first_index: Option<isize>,
    pub offsets: Option<Vec<isize>>,
    pub fixed_counts: Option<Vec<CyclicSievingFixedCountsInput>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum EhrhartHstarMode {
    HstarToEhrhart,
    EhrhartToHstar,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateSequenceRequest {
    pub sequence: SequenceKind,
    pub max_n: usize,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, Default)]
#[serde(deny_unknown_fields)]
pub struct ListOeisSequencesRequest {
    pub include_experimental: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GetOeisSequenceRequest {
    pub id: String,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateOeisRowsRequest {
    pub id: String,
    pub rows: usize,
    pub first_row: Option<i64>,
    pub include_experimental: Option<bool>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SequenceKind {
    Eulerian,
    Narayana,
    TypeBEulerian,
    ChebyshevT,
    ChebyshevU,
    Hermite,
}

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

#[derive(Debug, Clone, PartialEq, Eq)]
struct NormalizedBigIntPolynomial {
    display: BigIntNormalizedPolynomial,
    coefficients: Vec<BigInt>,
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
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateRecurrenceRowsRequest {
    pub recurrence_json: String,
    pub rows: Option<usize>,
    pub additional: Option<usize>,
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
pub struct CyclicSievingResponse {
    pub report: CyclicSievingReportItem,
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

type ParsedBatch = Vec<Result<NormalizedPolynomial, String>>;
type BigIntParsedBatch = Vec<Result<NormalizedBigIntPolynomial, String>>;
type RationalParsedBatch = Vec<Result<Vec<BigRational>, String>>;
type FamilySourceParse = (
    String,
    Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>>,
);

fn invalid_params(message: impl Into<String>) -> McpError {
    McpError::invalid_params(message.into(), None)
}

fn validate_mcp_text(text: &str, description: &str) -> Result<(), String> {
    if text.len() > MCP_MAX_TEXT_BYTES {
        return Err(format!(
            "{description} is {} bytes; the MCP limit is {MCP_MAX_TEXT_BYTES}",
            text.len()
        ));
    }
    Ok(())
}

fn validate_mcp_polynomial_length<T>(coefficients: Vec<T>) -> Result<Vec<T>, String> {
    if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            coefficients.len()
        ));
    }
    Ok(coefficients)
}

fn validate_mcp_batch_size(item_count: usize) -> Result<(), McpError> {
    if item_count > MCP_MAX_BATCH_POLYNOMIALS {
        return Err(invalid_params(format!(
            "batch has {item_count} polynomials; the MCP limit is {MCP_MAX_BATCH_POLYNOMIALS}"
        )));
    }
    Ok(())
}

fn parse_polynomial_input(input: &PolynomialInput) -> Result<Vec<i64>, String> {
    match (&input.coefficients, &input.expression) {
        (Some(coefficients), None) => validate_mcp_polynomial_length(coefficients.clone()),
        (None, Some(expression)) => {
            validate_mcp_text(expression, "polynomial expression")?;
            validate_mcp_polynomial_length(parse_polynomial(expression)?)
        }
        (Some(_), Some(_)) => {
            Err("expected exactly one of `coefficients` or `expression`, got both".to_string())
        }
        (None, None) => Err("expected exactly one of `coefficients` or `expression`".to_string()),
    }
}

fn parse_bigint_coefficient(input: &BigIntCoefficientInput) -> Result<BigInt, String> {
    match input {
        BigIntCoefficientInput::Integer(value) => Ok(BigInt::from(*value)),
        BigIntCoefficientInput::Text(value) => {
            if value.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                return Err(format!(
                    "integer coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                    value.len()
                ));
            }
            value
                .parse::<BigInt>()
                .map_err(|e| format!("invalid integer '{}': {}", value, e))
        }
    }
}

fn parse_bigint_coefficients(input: &[BigIntCoefficientInput]) -> Result<Vec<BigInt>, McpError> {
    if input.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(invalid_params(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            input.len()
        )));
    }
    input
        .iter()
        .map(parse_bigint_coefficient)
        .collect::<Result<Vec<_>, _>>()
        .map_err(invalid_params)
}

fn parse_bigint_polynomial_input(input: &BigIntPolynomialInput) -> Result<Vec<BigInt>, String> {
    match (&input.coefficients, &input.expression) {
        (Some(coefficients), None) => validate_mcp_polynomial_length(
            coefficients
                .iter()
                .enumerate()
                .map(|(index, coefficient)| {
                    parse_bigint_coefficient(coefficient)
                        .map_err(|error| format!("coefficient {index}: {error}"))
                })
                .collect::<Result<Vec<_>, _>>()?,
        ),
        (None, Some(expression)) => {
            validate_mcp_text(expression, "polynomial expression")?;
            validate_mcp_polynomial_length(parse_polynomial_bigint(expression)?)
        }
        (Some(_), Some(_)) => {
            Err("expected exactly one of `coefficients` or `expression`, got both".to_string())
        }
        (None, None) => Err("expected exactly one of `coefficients` or `expression`".to_string()),
    }
}

fn parse_required_bigint_polynomial(
    input: &BigIntPolynomialInput,
    name: &str,
) -> Result<Vec<BigInt>, McpError> {
    parse_bigint_polynomial_input(input).map_err(|e| invalid_params(format!("{name}: {e}")))
}

fn normalize_coefficients(mut coefficients: Vec<i64>) -> Vec<i64> {
    while coefficients.len() > 1 && coefficients.last() == Some(&0) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(0);
    }
    coefficients
}

fn normalize_rational_coefficients(mut coefficients: Vec<BigRational>) -> Vec<BigRational> {
    let zero = parse_rational("0").expect("0 parses as BigRational");
    while coefficients.len() > 1 && coefficients.last() == Some(&zero) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(zero);
    }
    coefficients
}

fn normalize_bigint_coefficients(mut coefficients: Vec<BigInt>) -> Vec<BigInt> {
    while coefficients.len() > 1 && coefficients.last() == Some(&BigInt::from(0)) {
        coefficients.pop();
    }
    if coefficients.is_empty() {
        coefficients.push(BigInt::from(0));
    }
    coefficients
}

fn degree(coefficients: &[i64]) -> usize {
    coefficients.iter().rposition(|&c| c != 0).unwrap_or(0)
}

fn bigint_degree(coefficients: &[BigInt]) -> usize {
    coefficients
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .unwrap_or(0)
}

fn bigint_strings(coefficients: &[BigInt]) -> Vec<String> {
    coefficients.iter().map(ToString::to_string).collect()
}

fn normalize_polynomial(coefficients: Vec<i64>) -> NormalizedPolynomial {
    let coefficients = normalize_coefficients(coefficients);
    NormalizedPolynomial {
        polynomial: format_poly(&coefficients),
        degree: degree(&coefficients),
        coefficients,
    }
}

fn normalize_bigint_polynomial(coefficients: Vec<BigInt>) -> NormalizedBigIntPolynomial {
    let coefficients = normalize_bigint_coefficients(coefficients);
    NormalizedBigIntPolynomial {
        display: BigIntNormalizedPolynomial {
            polynomial: format_poly_bigint_coeffs(&coefficients),
            coefficients: bigint_strings(&coefficients),
            degree: bigint_degree(&coefficients),
        },
        coefficients,
    }
}

fn parse_batch(input: &PolynomialBatchInput) -> Result<ParsedBatch, McpError> {
    let parsed: ParsedBatch = match (&input.polynomials, &input.text) {
        (Some(polynomials), None) => {
            validate_mcp_batch_size(polynomials.len())?;
            polynomials
                .iter()
                .map(|p| parse_polynomial_input(p).map(normalize_polynomial))
                .collect()
        }
        (None, Some(text)) => {
            validate_mcp_text(text, "polynomial batch").map_err(invalid_params)?;
            validate_mcp_batch_size(
                text.lines()
                    .filter(|line| {
                        let line = line.trim();
                        !line.is_empty() && !line.starts_with('#')
                    })
                    .count(),
            )?;
            parse_polynomials(text)
                .into_iter()
                .map(|result| {
                    result
                        .and_then(validate_mcp_polynomial_length)
                        .map(normalize_polynomial)
                })
                .collect()
        }
        (Some(_), Some(_)) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`, got both",
        ))?,
        (None, None) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`",
        ))?,
    };
    let total = parsed
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.coefficients.len())
        });
    let total = total.ok_or_else(|| invalid_params("total coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "batch has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(parsed)
}

fn parse_bigint_batch(input: &BigIntPolynomialBatchInput) -> Result<BigIntParsedBatch, McpError> {
    let parsed: BigIntParsedBatch = match (&input.polynomials, &input.text) {
        (Some(polynomials), None) => {
            validate_mcp_batch_size(polynomials.len())?;
            polynomials
                .iter()
                .map(|p| parse_bigint_polynomial_input(p).map(normalize_bigint_polynomial))
                .collect()
        }
        (None, Some(text)) => {
            validate_mcp_text(text, "polynomial batch").map_err(invalid_params)?;
            validate_mcp_batch_size(
                text.lines()
                    .filter(|line| {
                        let line = line.trim();
                        !line.is_empty() && !line.starts_with('#')
                    })
                    .count(),
            )?;
            parse_polynomials_bigint(text)
                .into_iter()
                .map(|result| {
                    result
                        .and_then(validate_mcp_polynomial_length)
                        .map(normalize_bigint_polynomial)
                })
                .collect()
        }
        (Some(_), Some(_)) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`, got both",
        ))?,
        (None, None) => Err(invalid_params(
            "expected exactly one of `polynomials` or `text`",
        ))?,
    };
    let total = parsed
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.coefficients.len())
        });
    let total = total.ok_or_else(|| invalid_params("total coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "batch has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(parsed)
}

fn i64_coefficients_to_rational(coefficients: &[i64]) -> Vec<BigRational> {
    coefficients
        .iter()
        .map(|coeff| parse_rational(&coeff.to_string()).expect("i64 parses as BigRational"))
        .collect()
}

fn parse_rational_coefficient(input: &RationalCoefficientInput) -> Result<BigRational, String> {
    match input {
        RationalCoefficientInput::Integer(value) => parse_rational(&value.to_string()),
        RationalCoefficientInput::Text(value) => {
            if value.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                return Err(format!(
                    "rational coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                    value.len()
                ));
            }
            parse_rational(value)
        }
    }
}

fn parse_rational_coefficients(
    coefficients: &[RationalCoefficientInput],
) -> Result<Vec<BigRational>, String> {
    if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
        return Err(format!(
            "polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
            coefficients.len()
        ));
    }
    coefficients
        .iter()
        .enumerate()
        .map(|(index, coefficient)| {
            parse_rational_coefficient(coefficient)
                .map_err(|error| format!("coefficient {index}: {error}"))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(normalize_rational_coefficients)
}

fn validate_recurrence_batch(batch: RationalParsedBatch) -> Result<RationalParsedBatch, McpError> {
    if batch.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
        return Err(invalid_params(format!(
            "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
            batch.len()
        )));
    }
    let total = batch
        .iter()
        .filter_map(|item| item.as_ref().ok())
        .try_fold(0usize, |total, polynomial| {
            total.checked_add(polynomial.len())
        })
        .ok_or_else(|| invalid_params("total recurrence coefficient count overflow"))?;
    if total > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "recurrence input has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(batch)
}

fn parse_recurrence_batch_rational(
    input: &FindRecurrenceRequest,
) -> Result<RationalParsedBatch, McpError> {
    let source_count = [
        input.polynomials.is_some(),
        input.coefficients.is_some(),
        input.expressions.is_some(),
        input.text.is_some(),
    ]
    .into_iter()
    .filter(|present| *present)
    .count();

    if source_count != 1 {
        return Err(invalid_params(
            "expected exactly one of `polynomials`, `coefficients`, `expressions`, or `text`",
        ));
    }

    if let Some(polynomials) = &input.polynomials {
        if polynomials.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                polynomials.len()
            )));
        }
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials.clone()),
            text: None,
        })?;
        return validate_recurrence_batch(
            batch
                .into_iter()
                .map(|item| {
                    item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
                })
                .collect(),
        );
    }

    if let Some(coefficients) = &input.coefficients {
        if coefficients.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                coefficients.len()
            )));
        }
        let total = coefficients
            .iter()
            .try_fold(0usize, |total, row| total.checked_add(row.len()))
            .ok_or_else(|| invalid_params("total recurrence coefficient count overflow"))?;
        if total > MCP_MAX_TOTAL_COEFFICIENTS {
            return Err(invalid_params(format!(
                "recurrence input has {total} coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
            )));
        }
        return validate_recurrence_batch(
            coefficients
                .iter()
                .map(|coefficients| parse_rational_coefficients(coefficients))
                .collect(),
        );
    }

    if let Some(expressions) = &input.expressions {
        if expressions.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
                expressions.len()
            )));
        }
        let polynomials = expressions
            .iter()
            .cloned()
            .map(|expression| PolynomialInput {
                coefficients: None,
                expression: Some(expression),
            })
            .collect();
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: Some(polynomials),
            text: None,
        })?;
        return validate_recurrence_batch(
            batch
                .into_iter()
                .map(|item| {
                    item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
                })
                .collect(),
        );
    }

    if let Some(text) = &input.text {
        validate_mcp_text(text, "recurrence polynomial batch").map_err(invalid_params)?;
        let row_count = text
            .lines()
            .filter(|line| {
                let line = line.trim();
                !line.is_empty() && !line.starts_with('#')
            })
            .count();
        if row_count > MCP_MAX_RECURRENCE_INPUT_ROWS {
            return Err(invalid_params(format!(
                "recurrence input has {row_count} rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}"
            )));
        }
    }
    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: None,
        text: input.text.clone(),
    })?;
    validate_recurrence_batch(
        batch
            .into_iter()
            .map(|item| {
                item.map(|polynomial| i64_coefficients_to_rational(&polynomial.coefficients))
            })
            .collect(),
    )
}

fn parse_items(batch: &ParsedBatch) -> Vec<ParsePolynomialItem> {
    batch
        .iter()
        .enumerate()
        .map(|(index, item)| match item {
            Ok(polynomial) => ParsePolynomialItem {
                index,
                ok: true,
                polynomial: Some(polynomial.polynomial.clone()),
                coefficients: Some(polynomial.coefficients.clone()),
                degree: Some(polynomial.degree),
                error: None,
            },
            Err(error) => ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            },
        })
        .collect()
}

fn parse_bigint_items(batch: &BigIntParsedBatch) -> Vec<BigIntParsePolynomialItem> {
    batch
        .iter()
        .enumerate()
        .map(|(index, item)| match item {
            Ok(polynomial) => BigIntParsePolynomialItem {
                index,
                ok: true,
                polynomial: Some(polynomial.display.clone()),
                error: None,
            },
            Err(error) => BigIntParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                error: Some(error.clone()),
            },
        })
        .collect()
}

fn hstar_check_item(check: &HStarInequalityCheck) -> HStarInequalityCheckItem {
    HStarInequalityCheckItem {
        family: check.family.clone(),
        name: check.name.clone(),
        formula: check.formula.clone(),
        reference: check.reference.clone(),
        url: check.url.clone(),
        index: check.index,
        applicable: check.applicable,
        holds: check.holds,
        lhs: check.lhs.as_ref().map(ToString::to_string),
        rhs: check.rhs.as_ref().map(ToString::to_string),
        value: check.value.clone(),
        details: check.details.clone(),
    }
}

fn hstar_inequalities_item(
    index: usize,
    item: Result<NormalizedBigIntPolynomial, String>,
    dimension: Option<usize>,
) -> HStarInequalitiesItem {
    match item {
        Ok(polynomial) => {
            let report = hstar_inequality_report_bigint(&polynomial.coefficients, dimension);
            HStarInequalitiesItem {
                index,
                ok: true,
                error: None,
                hstar: Some(bigint_strings(&report.hstar)),
                dimension: Some(report.dimension),
                degree: Some(report.degree),
                codegree: report.codegree,
                all_applicable_hold: Some(report.all_applicable_hold),
                checks: report.checks.iter().map(hstar_check_item).collect(),
            }
        }
        Err(error) => HStarInequalitiesItem {
            index,
            ok: false,
            error: Some(error),
            hstar: None,
            dimension: None,
            degree: None,
            codegree: None,
            all_applicable_hold: None,
            checks: Vec::new(),
        },
    }
}

fn coefficient_check_item(check: &CoefficientInequalityCheck) -> CoefficientInequalityCheckItem {
    CoefficientInequalityCheckItem {
        index: check.index,
        lhs: check.lhs.to_string(),
        rhs: check.rhs.to_string(),
        comparison: check.comparison.to_string(),
        holds: check.holds,
    }
}

fn coefficient_criterion_item(report: &CoefficientCriterionReport) -> CoefficientCriterionItem {
    CoefficientCriterionItem {
        name: report.name.to_string(),
        reference: report.reference.to_string(),
        applicable: report.applicable,
        holds: report.holds,
        implies_real_rooted: report.implies_real_rooted,
        reason: report.reason.clone(),
        checks: report.checks.iter().map(coefficient_check_item).collect(),
    }
}

fn coefficient_tests_item(
    index: usize,
    item: Result<NormalizedBigIntPolynomial, String>,
) -> CoefficientTestsItem {
    match item {
        Ok(polynomial) => {
            let report = coefficient_test_report_bigint(&polynomial.coefficients);
            CoefficientTestsItem {
                index,
                ok: true,
                error: None,
                polynomial: Some(polynomial.display),
                newton: Some(coefficient_criterion_item(&report.newton)),
                kurtz: Some(coefficient_criterion_item(&report.kurtz)),
            }
        }
        Err(error) => CoefficientTestsItem {
            index,
            ok: false,
            error: Some(error),
            polynomial: None,
            newton: None,
            kurtz: None,
        },
    }
}

fn root_evaluation_item(evaluation: &RootOfUnityEvaluation) -> RootOfUnityEvaluationItem {
    RootOfUnityEvaluationItem {
        group_order: evaluation.group_order,
        power: evaluation.power,
        root_order: evaluation.root_order,
        integer_value: evaluation.integer_value.as_ref().map(ToString::to_string),
        remainder: bigint_strings(&evaluation.remainder),
        remainder_polynomial: format_poly_bigint_coeffs(&evaluation.remainder),
    }
}

fn cyclic_power_check_item(check: &CyclicSievingPowerCheck) -> CyclicSievingPowerCheckItem {
    CyclicSievingPowerCheckItem {
        power: check.power,
        root_order: check.root_order,
        expected_fixed_points: check.expected_fixed_points.to_string(),
        evaluation: root_evaluation_item(&check.evaluation),
        holds: check.holds,
    }
}

fn cyclic_sieving_report_item(report: &CyclicSievingReport) -> CyclicSievingReportItem {
    CyclicSievingReportItem {
        order: report.order,
        polynomial: normalize_bigint_polynomial(report.coefficients.clone()).display,
        fixed_counts: report
            .fixed_counts
            .as_ref()
            .map(|counts| bigint_strings(counts)),
        holds: report.holds,
        evaluations: report
            .evaluations
            .iter()
            .map(root_evaluation_item)
            .collect(),
        checks: report.checks.iter().map(cyclic_power_check_item).collect(),
    }
}

fn cyclic_sequence_item_response(
    item: &CyclicSievingSequenceItem,
) -> CyclicSievingSequenceItemResponse {
    CyclicSievingSequenceItemResponse {
        row: item.row,
        index: item.index,
        polynomial: normalize_bigint_polynomial(item.coefficients.clone()).display,
        candidate_orders: item
            .candidate_orders
            .iter()
            .map(cyclic_sieving_report_item)
            .collect(),
    }
}

fn collect_polynomials_or_errors(
    batch: ParsedBatch,
) -> Result<Vec<NormalizedPolynomial>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

fn collect_rational_polynomials_or_errors(
    batch: RationalParsedBatch,
) -> Result<Vec<Vec<BigRational>>, Vec<ParsePolynomialItem>> {
    let errors: Vec<ParsePolynomialItem> = batch
        .iter()
        .enumerate()
        .filter_map(|(index, item)| match item {
            Ok(_) => None,
            Err(error) => Some(ParsePolynomialItem {
                index,
                ok: false,
                polynomial: None,
                coefficients: None,
                degree: None,
                error: Some(error.clone()),
            }),
        })
        .collect();
    if errors.is_empty() {
        Ok(batch.into_iter().map(Result::unwrap).collect())
    } else {
        Err(errors)
    }
}

fn interlacing_result(
    p: NormalizedBigIntPolynomial,
    q: NormalizedBigIntPolynomial,
) -> InterlacingPairResult {
    let strict = check_interlacing_bigint_coeffs(&p.coefficients, &q.coefficients);
    let weak = check_weak_interlacing_bigint_coeffs(&p.coefficients, &q.coefficients);
    let status = match (strict, weak) {
        (Some(true), _) => "strictly_interlace",
        (_, Some(true)) => "weakly_interlace",
        (Some(false), Some(false)) => "do_not_interlace",
        (Some(false), None) => "not_real_rooted_or_incompatible",
        (None, Some(false)) => "not_real_rooted_or_incompatible",
        (None, None) => "not_real_rooted_or_incompatible",
    }
    .to_string();
    InterlacingPairResult {
        p: p.display,
        q: q.display,
        strict,
        weak,
        status,
    }
}

fn interlacing_result_i64(
    p: NormalizedPolynomial,
    q: NormalizedPolynomial,
) -> InterlacingPairResult {
    let p = normalize_bigint_polynomial(p.coefficients.into_iter().map(BigInt::from).collect());
    let q = normalize_bigint_polynomial(q.coefficients.into_iter().map(BigInt::from).collect());
    interlacing_result(p, q)
}

fn interlacing_pair_has_interlacing(item: &InterlacingPairItem) -> bool {
    item.result.as_ref().is_some_and(|result| {
        matches!(
            result.status.as_str(),
            "strictly_interlace" | "weakly_interlace"
        )
    })
}

fn interlacing_pair_item(
    pair_index: usize,
    left_index: usize,
    right_index: usize,
    left: &Result<NormalizedBigIntPolynomial, String>,
    right: &Result<NormalizedBigIntPolynomial, String>,
) -> InterlacingPairItem {
    match (left, right) {
        (Ok(p), Ok(q)) => InterlacingPairItem {
            pair_index,
            left_index,
            right_index,
            ok: true,
            result: Some(interlacing_result(p.clone(), q.clone())),
            error: None,
        },
        _ => {
            let mut messages = Vec::new();
            if let Err(error) = left {
                messages.push(format!("polynomial {left_index}: {error}"));
            }
            if let Err(error) = right {
                messages.push(format!("polynomial {right_index}: {error}"));
            }
            InterlacingPairItem {
                pair_index,
                left_index,
                right_index,
                ok: false,
                result: None,
                error: Some(messages.join("; ")),
            }
        }
    }
}

fn integer_polys_to_rational(coefficients: &[Vec<i64>]) -> Vec<Vec<BigRational>> {
    coefficients
        .iter()
        .map(|row| {
            row.iter()
                .map(|coeff| {
                    parse_rational_coeff(&coeff.to_string())
                        .expect("integer coefficients parse as rationals")
                })
                .collect()
        })
        .collect()
}

fn recurrence_json_string(
    result: &AdaptiveSearchResult,
    polynomials: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    source_rows: usize,
) -> String {
    let searched_polys = polynomials.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result
        .recurrence
        .generation_initial_count(1, searched_polys.len());
    let initial_polys = &searched_polys[..initial_count];
    let recurrence_json = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        initial_polys,
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows,
            skip_prefix: search.skip_prefix,
            unknowns: result.num_unknowns,
            weighted_unknowns: result.weighted_unknowns,
            equations: result.num_equations,
            fit_polynomials: result.fit_polynomials,
            verification_polynomials: result.verification_polynomials,
            candidates_tried: result.candidates_tried,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    );
    serde_json::to_string_pretty(&recurrence_json).expect("serialize recurrence JSON")
}

fn family_report(index: usize, polynomial: NormalizedPolynomial) -> FamilyPolynomialReport {
    let coefficients = &polynomial.coefficients;
    FamilyPolynomialReport {
        index,
        real_rooted: is_real_rooted(coefficients),
        simple_roots: has_simple_roots(coefficients),
        palindromic: is_palindromic_ignoring_initial_zeros(coefficients),
        gamma_positive: is_gamma_positive_ignoring_initial_zeros(coefficients),
        gamma_coefficients: gamma_coefficients_ignoring_initial_zeros(coefficients),
        unimodal: is_unimodal(coefficients),
        log_concave: is_log_concave(coefficients),
        ultra_log_concave: is_ultra_log_concave(coefficients),
        polynomial,
    }
}

fn polynomial_properties_item_i64(
    index: usize,
    polynomial: BigIntNormalizedPolynomial,
    coefficients: &[i64],
) -> PolynomialPropertiesItem {
    let gamma = gamma_coefficients_ignoring_initial_zeros(coefficients);
    let gamma_positive = gamma
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|&value| value >= 0));
    let gamma_coefficients =
        gamma.map(|gamma| gamma.into_iter().map(|value| value.to_string()).collect());
    PolynomialPropertiesItem {
        index,
        ok: true,
        error: None,
        polynomial: Some(polynomial),
        real_rooted: Some(is_real_rooted(coefficients)),
        simple_roots: Some(has_simple_roots(coefficients)),
        palindromic: Some(is_palindromic_ignoring_initial_zeros(coefficients)),
        gamma_positive: Some(gamma_positive),
        gamma_coefficients,
        unimodal: Some(is_unimodal(coefficients)),
        log_concave: Some(is_log_concave(coefficients)),
        ultra_log_concave: Some(is_ultra_log_concave(coefficients)),
    }
}

fn polynomial_properties_item_bigint(
    index: usize,
    polynomial: NormalizedBigIntPolynomial,
) -> PolynomialPropertiesItem {
    if let Some(coefficients) = bigint_coeffs_to_i64(&polynomial.coefficients) {
        return polynomial_properties_item_i64(index, polynomial.display, &coefficients);
    }

    let zero = BigInt::from(0);
    let gamma_coefficients =
        gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(&polynomial.coefficients);
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|value| value >= &zero));
    PolynomialPropertiesItem {
        index,
        ok: true,
        error: None,
        polynomial: Some(polynomial.display),
        real_rooted: Some(is_real_rooted_bigint_coeffs(&polynomial.coefficients)),
        simple_roots: Some(has_simple_roots_bigint_coeffs(&polynomial.coefficients)),
        palindromic: Some(is_palindromic_ignoring_initial_zeros_bigint_coeffs(
            &polynomial.coefficients,
        )),
        gamma_positive: Some(gamma_positive),
        gamma_coefficients: gamma_coefficients.map(|gamma| bigint_strings(&gamma)),
        unimodal: Some(is_unimodal_bigint_coeffs(&polynomial.coefficients)),
        log_concave: Some(is_log_concave_bigint_coeffs(&polynomial.coefficients)),
        ultra_log_concave: Some(is_ultra_log_concave_bigint_coeffs(&polynomial.coefficients)),
    }
}

fn family_source_and_polynomials(
    input: &CheckPolynomialFamilyRequest,
) -> Result<FamilySourceParse, McpError> {
    let explicit_sources = usize::from(input.polynomials.is_some())
        + usize::from(input.text.is_some())
        + usize::from(input.sequence.is_some());
    if explicit_sources != 1 {
        return Err(invalid_params(
            "expected exactly one source: `polynomials`, `text`, or `sequence`",
        ));
    }

    if let Some(sequence) = input.sequence.clone() {
        let max_n = input
            .max_n
            .ok_or_else(|| invalid_params("`max_n` is required with `sequence`"))?;
        if max_n > MCP_MAX_SEQUENCE_N {
            return Err(invalid_params(format!(
                "max_n must be at most {MCP_MAX_SEQUENCE_N} for MCP sequence generation"
            )));
        }
        let bigint_polynomials = generated_sequence_polynomials_bigint(&sequence, max_n);
        let polynomials = bigint_polynomials
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|coefficient| {
                        i64::try_from(&coefficient).map_err(|_| {
                            invalid_params(format!(
                                "sequence {sequence:?} at max_n={max_n} exceeds the i64 coefficient range supported by check_polynomial_family"
                            ))
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()
                    .map(normalize_polynomial)
            })
            .collect::<Result<Vec<_>, _>>()?;
        return Ok((
            format!("sequence:{sequence:?}, max_n={max_n}"),
            Ok(polynomials),
        ));
    }

    if input.max_n.is_some() {
        return Err(invalid_params(
            "`max_n` is only valid when `sequence` is provided",
        ));
    }

    let batch = parse_batch(&PolynomialBatchInput {
        polynomials: input.polynomials.clone(),
        text: input.text.clone(),
    })?;
    let source = if input.text.is_some() {
        "text".to_string()
    } else {
        "polynomials".to_string()
    };
    Ok((source, collect_polynomials_or_errors(batch)))
}

fn generated_sequence_polynomials_bigint(
    sequence: &SequenceKind,
    max_n: usize,
) -> Vec<Vec<BigInt>> {
    match sequence {
        SequenceKind::Eulerian => eulerian_polynomials_bigint(max_n),
        SequenceKind::Narayana => narayana_polynomials_bigint(max_n),
        SequenceKind::TypeBEulerian => type_b_eulerian_polynomials_bigint(max_n),
        SequenceKind::ChebyshevT => chebyshev_polynomials_t_bigint(max_n),
        SequenceKind::ChebyshevU => chebyshev_polynomials_u_bigint(max_n),
        SequenceKind::Hermite => hermite_polynomials_bigint(max_n),
    }
}

fn oeis_sequence_summary(entry: &polytool::oeis::OeisSequenceDefinition) -> OeisSequenceSummary {
    OeisSequenceSummary {
        id: entry.id.to_string(),
        name: entry.name.to_string(),
        status: entry.status.as_str().to_string(),
        layout: entry.layout.as_str().to_string(),
        first_row: entry.first_row,
        flattened_offset: entry.flattened_offset,
        bfile_available: entry.bfile_prefix_verified,
        source_rows: entry.source_rows,
        verification_rows: entry.verification_rows,
    }
}

fn check_family_lace(
    polynomials: &[NormalizedPolynomial],
    request: &LaceCheckRequest,
) -> LaceCheckResponse {
    let coefficients: Vec<Vec<i64>> = polynomials.iter().map(|p| p.coefficients.clone()).collect();
    let max_degree = polynomials.iter().map(|p| p.degree).max().unwrap_or(0);
    let block_rows = request
        .block_rows
        .unwrap_or_else(|| polynomials.len().max(1));
    let block_cols = request
        .block_cols
        .unwrap_or_else(|| block_rows + max_degree + 1);
    let rows = match polynomials.len().checked_mul(block_rows) {
        Some(rows) => rows,
        None => {
            return LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size: request.max_minor_size.unwrap_or(0),
                rows: 0,
                columns: block_cols,
                tnn: false,
                error: Some("finite Lace matrix row count overflow".to_string()),
                matrix: None,
            };
        }
    };
    let max_minor_size = request
        .max_minor_size
        .unwrap_or_else(|| rows.min(block_cols).min(4));
    let cells = rows.checked_mul(block_cols);
    let minor_count =
        (1..=max_minor_size.min(rows).min(block_cols)).try_fold(0usize, |total, size| {
            binomial_capped(rows, size, MCP_MAX_LACE_MINORS)
                .checked_mul(binomial_capped(block_cols, size, MCP_MAX_LACE_MINORS))
                .and_then(|count| total.checked_add(count))
        });
    if cells.is_none_or(|cells| cells > MCP_MAX_LACE_MATRIX_CELLS)
        || minor_count.is_none_or(|count| count > MCP_MAX_LACE_MINORS)
    {
        return LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size,
            rows,
            columns: block_cols,
            tnn: false,
            error: Some(format!(
                "finite Lace request exceeds MCP limits ({MCP_MAX_LACE_MATRIX_CELLS} matrix cells and {MCP_MAX_LACE_MINORS} checked minors)"
            )),
            matrix: None,
        };
    }

    match lace_matrix_sequence_i64(&coefficients, block_rows, block_cols) {
        Ok(matrix) => {
            let rows = matrix.len();
            let columns = matrix.first().map(Vec::len).unwrap_or(0);
            let check = check_lace_sequence_total_nonnegative_i64(
                &coefficients,
                block_rows,
                block_cols,
                max_minor_size,
            );
            LaceCheckResponse {
                block_rows,
                block_cols,
                max_minor_size,
                rows,
                columns,
                tnn: check.is_ok(),
                error: check.err(),
                matrix: request.include_matrix.unwrap_or(false).then_some(matrix),
            }
        }
        Err(error) => LaceCheckResponse {
            block_rows,
            block_cols,
            max_minor_size: request.max_minor_size.unwrap_or(0),
            rows: 0,
            columns: 0,
            tnn: false,
            error: Some(error.to_string()),
            matrix: None,
        },
    }
}

fn binomial_capped(n: usize, k: usize, cap: usize) -> usize {
    if k > n {
        return 0;
    }
    let k = k.min(n - k);
    let mut result = 1usize;
    for divisor in 1..=k {
        result = result
            .checked_mul(n - k + divisor)
            .and_then(|value| value.checked_div(divisor))
            .unwrap_or(cap.saturating_add(1));
        if result > cap {
            return cap.saturating_add(1);
        }
    }
    result
}

fn recurrence_failure_response(
    status: RecurrenceSearchStatus,
    error: impl Into<String>,
    summary: Option<AdaptiveSearchSummary>,
    parse_errors: Vec<ParsePolynomialItem>,
) -> FindRecurrenceResponse {
    let (candidates_tried, candidates_considered) = summary.map_or((None, None), |summary| {
        (
            Some(summary.candidates_tried),
            Some(summary.diagnostics.considered_candidates),
        )
    });
    FindRecurrenceResponse {
        status,
        found: false,
        recurrence: None,
        latex: None,
        mathematica: None,
        sage: None,
        python: None,
        recurrence_json: None,
        unknowns: None,
        weighted_unknowns: None,
        equations: None,
        fit_polynomials: None,
        verification_polynomials: None,
        candidates_tried,
        candidates_considered,
        error: Some(error.into()),
        parse_errors,
    }
}

fn family_recurrence_response(
    coefficients: &[Vec<i64>],
    search: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> FindRecurrenceResponse {
    if coefficients.len() < 3 {
        return recurrence_failure_response(
            RecurrenceSearchStatus::InvalidInput,
            "need at least 3 polynomials",
            None,
            Vec::new(),
        );
    }

    match find_recurrence_adaptive_with_budget(coefficients, search, budget) {
        AdaptiveSearchOutcome::Found(result) => {
            let rational_polys = integer_polys_to_rational(coefficients);
            FindRecurrenceResponse {
                status: RecurrenceSearchStatus::Found,
                found: true,
                recurrence: Some(format!("{}", result.recurrence)),
                latex: Some(result.recurrence.to_latex()),
                mathematica: Some(result.recurrence.to_mathematica_definition(coefficients)),
                sage: Some(result.recurrence.to_sage_definition(coefficients)),
                python: Some(result.recurrence.to_python_definition(coefficients)),
                recurrence_json: Some(recurrence_json_string(
                    &result,
                    &rational_polys,
                    search,
                    coefficients.len(),
                )),
                unknowns: Some(result.num_unknowns),
                weighted_unknowns: Some(result.weighted_unknowns),
                equations: Some(result.num_equations),
                fit_polynomials: Some(result.fit_polynomials),
                verification_polynomials: Some(result.verification_polynomials),
                candidates_tried: Some(result.candidates_tried),
                candidates_considered: Some(result.diagnostics.considered_candidates),
                error: None,
                parse_errors: Vec::new(),
            }
        }
        AdaptiveSearchOutcome::NoRecurrence(summary) => recurrence_failure_response(
            RecurrenceSearchStatus::NotFound,
            "no recurrence found within the search bounds",
            Some(summary),
            Vec::new(),
        ),
        AdaptiveSearchOutcome::BudgetExhausted(summary) => recurrence_failure_response(
            RecurrenceSearchStatus::BudgetExhausted,
            format!(
                "recurrence candidate budget exhausted after {} candidates; unsearched candidates remain",
                summary.diagnostics.considered_candidates
            ),
            Some(summary),
            Vec::new(),
        ),
    }
}

fn first_family_failure(
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
    options: Option<&FamilyCheckOptions>,
) -> Option<String> {
    let require_real_rooted = options.and_then(|o| o.require_real_rooted).unwrap_or(true);
    let require_simple_roots = options
        .and_then(|o| o.require_simple_roots)
        .unwrap_or(false);
    let require_palindromic = options.and_then(|o| o.require_palindromic).unwrap_or(false);
    let require_gamma_positive = options
        .and_then(|o| o.require_gamma_positive)
        .unwrap_or(false);
    let require_unimodal = options.and_then(|o| o.require_unimodal).unwrap_or(false);
    let require_log_concave = options.and_then(|o| o.require_log_concave).unwrap_or(false);
    let require_ultra_log_concave = options
        .and_then(|o| o.require_ultra_log_concave)
        .unwrap_or(false);
    let require_interlacing = options
        .and_then(|o| o.require_consecutive_weak_interlacing)
        .unwrap_or(false);

    for item in items {
        if require_real_rooted && !item.real_rooted {
            return Some(format!("polynomial {} is not real-rooted", item.index));
        }
        if require_simple_roots && !item.simple_roots {
            return Some(format!(
                "polynomial {} does not have simple roots",
                item.index
            ));
        }
        if require_palindromic && !item.palindromic {
            return Some(format!("polynomial {} is not palindromic", item.index));
        }
        if require_gamma_positive && !item.gamma_positive {
            return Some(format!("polynomial {} is not gamma-positive", item.index));
        }
        if require_unimodal && !item.unimodal {
            return Some(format!("polynomial {} is not unimodal", item.index));
        }
        if require_log_concave && !item.log_concave {
            return Some(format!("polynomial {} is not log-concave", item.index));
        }
        if require_ultra_log_concave && !item.ultra_log_concave {
            return Some(format!(
                "polynomial {} is not ultra-log-concave",
                item.index
            ));
        }
    }

    if require_interlacing {
        for pair in pairs {
            if pair.result.as_ref().and_then(|r| r.weak) != Some(true) {
                return Some(format!(
                    "consecutive pair {}-{} is not weakly interlacing",
                    pair.left_index, pair.right_index
                ));
            }
        }
    }

    if let Some(lace) = lace {
        if !lace.tnn {
            return Some(format!(
                "finite Lace truncation is not TNN up to minors of size {}",
                lace.max_minor_size
            ));
        }
    }

    if let Some(recurrence) = recurrence {
        if !recurrence.found {
            return Some("recurrence search did not find a recurrence".to_string());
        }
    }

    None
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no"
    }
}

fn opt_yes_no(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "yes",
        Some(false) => "no",
        None => "n/a",
    }
}

fn family_markdown(
    source: &str,
    all_required_checks_passed: bool,
    first_failure: Option<&str>,
    items: &[FamilyPolynomialReport],
    pairs: &[InterlacingPairItem],
    lace: Option<&LaceCheckResponse>,
    recurrence: Option<&FindRecurrenceResponse>,
) -> String {
    let mut out = String::new();
    out.push_str("# Polynomial family check\n\n");
    out.push_str(&format!("- Source: `{source}`\n"));
    out.push_str(&format!("- Polynomials: {}\n", items.len()));
    out.push_str(&format!(
        "- Required checks: {}\n",
        if all_required_checks_passed {
            "passed"
        } else {
            "failed"
        }
    ));
    if let Some(first_failure) = first_failure {
        out.push_str(&format!("- First failure: {first_failure}\n"));
    }

    out.push_str("\n## Polynomial properties\n\n");
    out.push_str("| i | polynomial | degree | RR | gamma+ | unimodal | LC | ULC | pal |\n");
    out.push_str("|---:|---|---:|:---:|:---:|:---:|:---:|:---:|:---:|\n");
    for item in items {
        out.push_str(&format!(
            "| {} | `{}` | {} | {} | {} | {} | {} | {} | {} |\n",
            item.index,
            item.polynomial.polynomial,
            item.polynomial.degree,
            yes_no(item.real_rooted),
            yes_no(item.gamma_positive),
            yes_no(item.unimodal),
            yes_no(item.log_concave),
            yes_no(item.ultra_log_concave),
            yes_no(item.palindromic)
        ));
    }

    if !pairs.is_empty() {
        out.push_str("\n## Consecutive interlacing\n\n");
        out.push_str("| pair | weak | strict | status |\n");
        out.push_str("|---:|:---:|:---:|---|\n");
        for pair in pairs {
            if let Some(result) = &pair.result {
                out.push_str(&format!(
                    "| {}-{} | {} | {} | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    opt_yes_no(result.weak),
                    opt_yes_no(result.strict),
                    result.status
                ));
            } else {
                out.push_str(&format!(
                    "| {}-{} | n/a | n/a | `{}` |\n",
                    pair.left_index,
                    pair.right_index,
                    pair.error.as_deref().unwrap_or("parse error")
                ));
            }
        }
    }

    if let Some(lace) = lace {
        out.push_str("\n## Finite Lace check\n\n");
        out.push_str(&format!(
            "- Truncation: block rows `{}`, block cols `{}`; matrix size `{} x {}`.\n",
            lace.block_rows, lace.block_cols, lace.rows, lace.columns
        ));
        out.push_str(&format!(
            "- Checked minors up to size `{}`: {}.\n",
            lace.max_minor_size,
            if lace.tnn { "TNN" } else { "not TNN" }
        ));
        if let Some(error) = &lace.error {
            out.push_str(&format!("- Failure: `{error}`\n"));
        }
    }

    if let Some(recurrence) = recurrence {
        out.push_str("\n## Recurrence search\n\n");
        if let Some(found) = &recurrence.recurrence {
            out.push_str(&format!("- Found: `{found}`\n"));
        } else if let Some(error) = &recurrence.error {
            out.push_str(&format!("- Not found: {error}\n"));
        } else {
            out.push_str("- Not found.\n");
        }
    }

    out
}

fn displayed(coefficients: &[i64]) -> DisplayedPolynomial {
    DisplayedPolynomial {
        polynomial: format_poly(coefficients),
        coefficients: coefficients.to_vec(),
    }
}

fn format_rational<T: ToString>(rational: T) -> String {
    rational.to_string()
}

fn parse_rational(input: &str) -> Result<BigRational, String> {
    input
        .parse::<BigRational>()
        .map_err(|e| format!("invalid rational `{input}`: {e}"))
}

fn inclusive_option_count(minimum: usize, maximum: usize) -> Result<usize, McpError> {
    maximum
        .checked_sub(minimum)
        .and_then(|difference| difference.checked_add(1))
        .ok_or_else(|| invalid_params(format!("invalid recurrence range {minimum}..={maximum}")))
}

fn validate_recurrence_options(
    options: &AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
) -> Result<(), McpError> {
    let ranges = [
        (
            "recurrence length",
            options.min_rec_len,
            options.max_rec_len,
        ),
        ("variable degree", options.min_var_deg, options.max_var_deg),
        ("index degree", options.min_idx_deg, options.max_idx_deg),
        (
            "derivative degree",
            options.min_diff_deg,
            options.max_diff_deg,
        ),
        (
            "inhomogeneous variable degree",
            options.min_inhomo_var_deg,
            options.max_inhomo_var_deg,
        ),
        (
            "inhomogeneous index degree",
            options.min_inhomo_idx_deg,
            options.max_inhomo_idx_deg,
        ),
    ];
    for (name, minimum, maximum) in ranges {
        if minimum > maximum {
            return Err(invalid_params(format!(
                "minimum {name} {minimum} exceeds maximum {maximum}"
            )));
        }
    }
    if options.max_rec_len > MCP_MAX_RECURRENCE_LENGTH {
        return Err(invalid_params(format!(
            "maximum recurrence length {} exceeds the MCP limit {MCP_MAX_RECURRENCE_LENGTH}",
            options.max_rec_len
        )));
    }
    let degree_bounds = [
        ("variable", options.max_var_deg),
        ("index", options.max_idx_deg),
        ("derivative", options.max_diff_deg),
        ("inhomogeneous variable", options.max_inhomo_var_deg),
        ("inhomogeneous index", options.max_inhomo_idx_deg),
        ("denominator variable", options.max_denom_var_deg),
        ("denominator index", options.max_denom_idx_deg),
    ];
    for (name, degree) in degree_bounds {
        if degree > MCP_MAX_RECURRENCE_DEGREE {
            return Err(invalid_params(format!(
                "maximum {name} degree {degree} exceeds the MCP limit {MCP_MAX_RECURRENCE_DEGREE}"
            )));
        }
    }
    if options.skip_prefix > MCP_MAX_RECURRENCE_INPUT_ROWS
        || options.fit_extra_rows > MCP_MAX_RECURRENCE_INPUT_ROWS
    {
        return Err(invalid_params(format!(
            "skip_prefix and fit_extra_rows must be at most {MCP_MAX_RECURRENCE_INPUT_ROWS}"
        )));
    }
    if options.min_margin > MCP_MAX_RECURRENCE_UNKNOWNS {
        return Err(invalid_params(format!(
            "min_margin {} exceeds the MCP limit {MCP_MAX_RECURRENCE_UNKNOWNS}",
            options.min_margin
        )));
    }

    let base_count = inclusive_option_count(options.min_rec_len.max(1), options.max_rec_len)?
        .checked_mul(inclusive_option_count(
            options.min_diff_deg,
            options.max_diff_deg,
        )?)
        .and_then(|count| {
            count
                .checked_mul(inclusive_option_count(options.min_idx_deg, options.max_idx_deg).ok()?)
        })
        .and_then(|count| {
            count
                .checked_mul(inclusive_option_count(options.min_var_deg, options.max_var_deg).ok()?)
        })
        .and_then(|count| count.checked_mul(if options.try_alternating_sign { 2 } else { 1 }))
        .ok_or_else(|| invalid_params("recurrence candidate count overflow"))?;
    let inhomogeneous_count = if options.try_inhomogeneous {
        inclusive_option_count(options.min_inhomo_idx_deg, options.max_inhomo_idx_deg)?
            .checked_mul(inclusive_option_count(
                options.min_inhomo_var_deg,
                options.max_inhomo_var_deg,
            )?)
            .ok_or_else(|| invalid_params("inhomogeneous candidate count overflow"))?
    } else {
        0
    };
    let denominator_count = if options.try_denominator {
        options
            .max_denom_idx_deg
            .checked_add(1)
            .and_then(|count| count.checked_mul(options.max_denom_var_deg.checked_add(1)?))
            .and_then(|count| count.checked_sub(1))
            .ok_or_else(|| invalid_params("denominator candidate count overflow"))?
    } else {
        0
    };
    let candidate_count = 1usize
        .checked_add(inhomogeneous_count)
        .and_then(|count| count.checked_add(denominator_count))
        .and_then(|count| count.checked_mul(base_count))
        .ok_or_else(|| invalid_params("recurrence candidate count overflow"))?;
    let effective_candidate_count = budget
        .max_candidates
        .map_or(candidate_count, |limit| candidate_count.min(limit));
    if effective_candidate_count > MCP_MAX_RECURRENCE_CANDIDATES {
        return Err(invalid_params(format!(
            "recurrence search may inspect {effective_candidate_count} candidates; the MCP limit is {MCP_MAX_RECURRENCE_CANDIDATES}"
        )));
    }

    let base_options = RecurrenceOptions {
        rec_len: options.max_rec_len,
        var_deg: options.max_var_deg,
        idx_deg: options.max_idx_deg,
        diff_deg: options.max_diff_deg,
        homogeneous: true,
        inhomo_var_deg: 0,
        inhomo_idx_deg: 0,
        denom_var_deg: 0,
        denom_idx_deg: 0,
        alternating_sign: options.try_alternating_sign,
        modular_prefilter: options.modular_prefilter,
    };
    let mut maximum_unknowns = candidate_complexity(&base_options).raw_unknowns;
    if options.try_inhomogeneous {
        maximum_unknowns = maximum_unknowns.max(
            candidate_complexity(&RecurrenceOptions {
                homogeneous: false,
                inhomo_var_deg: options.max_inhomo_var_deg,
                inhomo_idx_deg: options.max_inhomo_idx_deg,
                ..base_options.clone()
            })
            .raw_unknowns,
        );
    }
    if options.try_denominator {
        maximum_unknowns = maximum_unknowns.max(
            candidate_complexity(&RecurrenceOptions {
                denom_var_deg: options.max_denom_var_deg,
                denom_idx_deg: options.max_denom_idx_deg,
                ..base_options
            })
            .raw_unknowns,
        );
    }
    if maximum_unknowns > MCP_MAX_RECURRENCE_UNKNOWNS {
        return Err(invalid_params(format!(
            "recurrence search may create {maximum_unknowns} unknowns; the MCP limit is {MCP_MAX_RECURRENCE_UNKNOWNS}"
        )));
    }
    Ok(())
}

fn validate_recurrence_for_generation(
    recurrence: &Recurrence,
    initial_polys: &[Vec<BigRational>],
) -> Result<(), McpError> {
    if recurrence.terms.len() > MCP_MAX_RECURRENCE_LENGTH {
        return Err(invalid_params(format!(
            "recurrence has {} terms; the MCP limit is {MCP_MAX_RECURRENCE_LENGTH}",
            recurrence.terms.len()
        )));
    }
    let mut total_cells = initial_polys
        .iter()
        .try_fold(0usize, |total, row| total.checked_add(row.len()))
        .ok_or_else(|| invalid_params("recurrence coefficient count overflow"))?;
    if initial_polys.len() > MCP_MAX_RECURRENCE_INPUT_ROWS {
        return Err(invalid_params(format!(
            "recurrence JSON has {} initial rows; the MCP limit is {MCP_MAX_RECURRENCE_INPUT_ROWS}",
            initial_polys.len()
        )));
    }
    if initial_polys
        .iter()
        .any(|row| row.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL)
    {
        return Err(invalid_params(format!(
            "an initial polynomial exceeds the MCP limit of {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL} coefficients"
        )));
    }

    let mut validate_bivar = |name: &str, polynomial: &BivarPoly| -> Result<(), McpError> {
        if polynomial.coeffs.len() > MCP_MAX_RECURRENCE_DEGREE + 1
            || polynomial
                .coeffs
                .iter()
                .any(|row| row.len() > MCP_MAX_RECURRENCE_DEGREE + 1)
        {
            return Err(invalid_params(format!(
                "{name} exceeds the MCP bivariate degree limit {MCP_MAX_RECURRENCE_DEGREE}"
            )));
        }
        for row in &polynomial.coeffs {
            total_cells = total_cells
                .checked_add(row.len())
                .ok_or_else(|| invalid_params("recurrence coefficient count overflow"))?;
        }
        Ok(())
    };
    for (index, term) in recurrence.terms.iter().enumerate() {
        if term.offset > MCP_MAX_RECURRENCE_LENGTH {
            return Err(invalid_params(format!(
                "recurrence term {index} has offset {}; the MCP limit is {MCP_MAX_RECURRENCE_LENGTH}",
                term.offset
            )));
        }
        if term.deriv_order > MCP_MAX_RECURRENCE_DEGREE {
            return Err(invalid_params(format!(
                "recurrence term {index} has derivative order {}; the MCP limit is {MCP_MAX_RECURRENCE_DEGREE}",
                term.deriv_order
            )));
        }
        validate_bivar(&format!("recurrence term {index}"), &term.coeff)?;
    }
    if let Some(denominator) = &recurrence.denominator {
        validate_bivar("recurrence denominator", denominator)?;
    }
    if let Some(inhomogeneous) = &recurrence.inhomogeneous {
        validate_bivar("recurrence inhomogeneous term", inhomogeneous)?;
    }
    if total_cells > MCP_MAX_TOTAL_COEFFICIENTS {
        return Err(invalid_params(format!(
            "recurrence JSON has {total_cells} stored coefficients; the MCP limit is {MCP_MAX_TOTAL_COEFFICIENTS}"
        )));
    }
    Ok(())
}

fn apply_recurrence_options(
    input: Option<RecurrenceSearchOptionsInput>,
) -> Result<(AdaptiveSearchOptions, AdaptiveSearchBudget), McpError> {
    let mut options = AdaptiveSearchOptions::default();
    let Some(input) = input else {
        options.verbose = false;
        let budget = AdaptiveSearchBudget::unbounded();
        validate_recurrence_options(&options, budget)?;
        return Ok((options, budget));
    };
    let budget = AdaptiveSearchBudget {
        max_candidates: input.max_candidates,
    };

    if let Some(value) = input.skip_prefix {
        options.skip_prefix = value;
    }
    if let Some(value) = input.min_rec_len {
        options.min_rec_len = value;
    }
    if let Some(value) = input.max_rec_len {
        options.max_rec_len = value;
    }
    if let Some(value) = input.min_var_deg {
        options.min_var_deg = value;
    }
    if let Some(value) = input.max_var_deg {
        options.max_var_deg = value;
    }
    if let Some(value) = input.min_idx_deg {
        options.min_idx_deg = value;
    }
    if let Some(value) = input.max_idx_deg {
        options.max_idx_deg = value;
    }
    if let Some(value) = input.min_diff_deg {
        options.min_diff_deg = value;
    }
    if let Some(value) = input.max_diff_deg {
        options.max_diff_deg = value;
    }
    if let Some(value) = input.try_inhomogeneous {
        options.try_inhomogeneous = value;
    }
    if let Some(value) = input.min_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_var_deg = value;
    }
    if let Some(value) = input.max_inhomo_var_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_var_deg = value;
    }
    if let Some(value) = input.min_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.min_inhomo_idx_deg = value;
    }
    if let Some(value) = input.max_inhomo_idx_deg {
        options.try_inhomogeneous = true;
        options.max_inhomo_idx_deg = value;
    }
    if let Some(value) = input.try_denominator {
        options.try_denominator = value;
    }
    if let Some(value) = input.try_alternating_sign {
        options.try_alternating_sign = value;
    }
    if let Some(value) = input.max_denom_var_deg {
        options.try_denominator = true;
        options.max_denom_var_deg = value;
    }
    if let Some(value) = input.max_denom_idx_deg {
        options.try_denominator = true;
        options.max_denom_idx_deg = value;
    }
    if let Some(value) = input.min_margin {
        options.min_margin = value;
    }
    if let Some(value) = input.no_verify {
        options.no_verify = value;
    }
    if let Some(value) = input.fit_extra_rows {
        options.fit_extra_rows = value;
    }
    if let Some(value) = input.modular_prefilter {
        options.modular_prefilter = value;
    }
    options.verbose = false;
    validate_recurrence_options(&options, budget)?;
    Ok((options, budget))
}

#[tool_router(router = tool_router)]
impl PolynomialToolsServer {
    pub fn new() -> Self {
        Self {
            tool_router: Self::tool_router(),
        }
    }

    #[tool(description = "Parse and format one or more dense univariate polynomials.")]
    pub fn parse_polynomials(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<ParsePolynomialsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        Ok(Json(ParsePolynomialsResponse {
            items: parse_items(&batch),
        }))
    }

    #[tool(
        description = "Check a polynomial family in one exact batch: properties, consecutive interlacing, optional finite Lace TNN, and optional recurrence search."
    )]
    pub fn check_polynomial_family(
        &self,
        Parameters(input): Parameters<CheckPolynomialFamilyRequest>,
    ) -> Result<Json<CheckPolynomialFamilyResponse>, McpError> {
        let (source, polynomials) = family_source_and_polynomials(&input)?;
        let polynomials = match polynomials {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                let first_failure = Some("one or more polynomials failed to parse".to_string());
                let mut markdown = String::new();
                markdown.push_str("# Polynomial family check\n\n");
                markdown.push_str(&format!("- Source: `{source}`\n"));
                markdown.push_str("- Required checks: failed\n");
                markdown.push_str("- First failure: one or more polynomials failed to parse\n\n");
                markdown.push_str("## Parse errors\n\n");
                for item in &parse_errors {
                    markdown.push_str(&format!(
                        "- polynomial {}: {}\n",
                        item.index,
                        item.error.as_deref().unwrap_or("parse error")
                    ));
                }
                return Ok(Json(CheckPolynomialFamilyResponse {
                    source,
                    item_count: 0,
                    all_required_checks_passed: false,
                    first_failure,
                    parse_errors,
                    items: Vec::new(),
                    consecutive_pairs: Vec::new(),
                    lace: None,
                    recurrence: None,
                    markdown,
                }));
            }
        };

        let options = input.options.as_ref();
        let find_recurrence = options.and_then(|o| o.find_recurrence).unwrap_or(false);
        let recurrence_search = if find_recurrence {
            Some(apply_recurrence_options(input.recurrence_options.clone())?)
        } else {
            None
        };
        let check_consecutive_interlacing = options
            .and_then(|o| o.check_consecutive_interlacing)
            .unwrap_or(true);

        let items: Vec<_> = polynomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, polynomial)| family_report(index, polynomial))
            .collect();

        let consecutive_pairs = if check_consecutive_interlacing {
            polynomials
                .windows(2)
                .enumerate()
                .map(|(pair_index, pair)| InterlacingPairItem {
                    pair_index,
                    left_index: pair_index,
                    right_index: pair_index + 1,
                    ok: true,
                    result: Some(interlacing_result_i64(pair[0].clone(), pair[1].clone())),
                    error: None,
                })
                .collect()
        } else {
            Vec::new()
        };

        let lace = input
            .lace
            .as_ref()
            .map(|request| check_family_lace(&polynomials, request));

        let coefficients: Vec<Vec<i64>> =
            polynomials.iter().map(|p| p.coefficients.clone()).collect();
        let recurrence = recurrence_search
            .as_ref()
            .map(|(search, budget)| family_recurrence_response(&coefficients, search, *budget));

        let first_failure = first_family_failure(
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
            options,
        );
        let all_required_checks_passed = first_failure.is_none();
        let markdown = family_markdown(
            &source,
            all_required_checks_passed,
            first_failure.as_deref(),
            &items,
            &consecutive_pairs,
            lace.as_ref(),
            recurrence.as_ref(),
        );

        Ok(Json(CheckPolynomialFamilyResponse {
            source,
            item_count: items.len(),
            all_required_checks_passed,
            first_failure,
            parse_errors: Vec::new(),
            items,
            consecutive_pairs,
            lace,
            recurrence,
            markdown,
        }))
    }

    #[tool(
        description = "Check real-rootedness, gamma-positivity, unimodality, and concavity properties."
    )]
    pub fn polynomial_properties(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<PolynomialPropertiesResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => polynomial_properties_item_bigint(index, polynomial),
                Err(error) => PolynomialPropertiesItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    simple_roots: None,
                    palindromic: None,
                    gamma_positive: None,
                    gamma_coefficients: None,
                    unimodal: None,
                    log_concave: None,
                    ultra_log_concave: None,
                },
            })
            .collect();
        Ok(Json(PolynomialPropertiesResponse { items }))
    }

    #[tool(
        description = "Check Newton inequalities and Kurtz's sufficient real-rootedness criterion for each polynomial."
    )]
    pub fn coefficient_tests(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<CoefficientTestsResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| coefficient_tests_item(index, item))
            .collect();
        Ok(Json(CoefficientTestsResponse { items }))
    }

    #[tool(
        description = "Check named Ehrhart h*-vector inequalities, including Stanley, Hibi, Balletti-Higashitani, and Stapledon-derived checks."
    )]
    pub fn hstar_inequalities(
        &self,
        Parameters(input): Parameters<HStarInequalitiesRequest>,
    ) -> Result<Json<HStarInequalitiesResponse>, McpError> {
        let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: input.polynomials,
            text: input.text,
        })?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| hstar_inequalities_item(index, item, input.dimension))
            .collect();
        Ok(Json(HStarInequalitiesResponse { items }))
    }

    #[tool(
        description = "Check strict and weak interlacing for a pair of polynomials with arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_pair(
        &self,
        Parameters(input): Parameters<BigIntPolynomialPairRequest>,
    ) -> Result<Json<InterlacingPairResult>, McpError> {
        let p = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.p, "p")?);
        let q = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.q, "q")?);
        Ok(Json(interlacing_result(p, q)))
    }

    #[tool(
        description = "Check strict and weak interlacing for consecutive polynomial pairs with arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_sequence(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<InterlacingSequenceResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = parse_bigint_items(&batch);
        let pairs = batch
            .windows(2)
            .enumerate()
            .map(|(pair_index, pair)| {
                interlacing_pair_item(pair_index, pair_index, pair_index + 1, &pair[0], &pair[1])
            })
            .collect();
        Ok(Json(InterlacingSequenceResponse { items, pairs }))
    }

    #[tool(
        description = "For each polynomial, count backward consecutive previous interlacings until the first failure. Accepts arbitrary-size integer coefficients."
    )]
    pub fn check_interlacing_profile(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<InterlacingProfileResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = parse_bigint_items(&batch);
        let mut pair_index = 0;
        let profile = batch
            .iter()
            .enumerate()
            .map(|(index, item)| match item {
                Err(error) => InterlacingProfileItem {
                    index,
                    ok: false,
                    polynomial: None,
                    previous_count: index,
                    checked_previous_count: 0,
                    interlacing_previous_count: 0,
                    strict_previous_count: 0,
                    weak_previous_count: 0,
                    previous_interlacing_indices: Vec::new(),
                    previous: Vec::new(),
                    error: Some(error.clone()),
                },
                Ok(polynomial) => {
                    let mut previous = Vec::new();
                    for left_index in (0..index).rev() {
                        let pair = interlacing_pair_item(
                            pair_index,
                            left_index,
                            index,
                            &batch[left_index],
                            item,
                        );
                        pair_index += 1;
                        let interlaces = interlacing_pair_has_interlacing(&pair);
                        previous.push(pair);
                        if !interlaces {
                            break;
                        }
                    }
                    let previous_interlacing_indices = previous
                        .iter()
                        .filter(|pair| interlacing_pair_has_interlacing(pair))
                        .map(|pair| pair.left_index)
                        .collect::<Vec<_>>();
                    InterlacingProfileItem {
                        index,
                        ok: true,
                        polynomial: Some(polynomial.display.clone()),
                        previous_count: index,
                        checked_previous_count: previous.len(),
                        interlacing_previous_count: previous_interlacing_indices.len(),
                        strict_previous_count: previous
                            .iter()
                            .filter(|pair| {
                                pair.result.as_ref().and_then(|result| result.strict) == Some(true)
                            })
                            .count(),
                        weak_previous_count: previous
                            .iter()
                            .filter(|pair| {
                                pair.result.as_ref().and_then(|result| result.weak) == Some(true)
                            })
                            .count(),
                        previous_interlacing_indices,
                        previous,
                        error: None,
                    }
                }
            })
            .collect();
        Ok(Json(InterlacingProfileResponse { items, profile }))
    }

    #[tool(
        description = "Return rational midpoint representatives for isolated real roots, or real_rooted=false."
    )]
    pub fn real_roots(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<RealRootsResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => match polytool::real_roots(&polynomial.coefficients) {
                    Some(roots) => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(true),
                        roots: Some(roots.into_iter().map(format_rational).collect()),
                    },
                    None => RealRootsItem {
                        index,
                        ok: true,
                        error: None,
                        polynomial: Some(polynomial),
                        real_rooted: Some(false),
                        roots: None,
                    },
                },
                Err(error) => RealRootsItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    real_rooted: None,
                    roots: None,
                },
            })
            .collect();
        Ok(Json(RealRootsResponse { items }))
    }

    #[tool(
        description = "Search adaptively for a polynomial recurrence. Provide exactly one of `polynomials`, `coefficients`, `expressions`, or `text`. Pass search controls such as `min_rec_len`, `max_rec_len`, degree bounds, `try_denominator`, `fit_extra_rows`, `modular_prefilter`, and `max_candidates` as top-level arguments. Legacy nested `options` is accepted; a top-level control wins when both are supplied. Set `include_code` to false for a compact result without Mathematica, Python, Sage, or recurrence JSON."
    )]
    pub fn find_recurrence(
        &self,
        Parameters(input): Parameters<FindRecurrenceRequest>,
    ) -> Result<Json<FindRecurrenceResponse>, McpError> {
        let recurrence_options = input.effective_options();
        let include_code = input.include_code.unwrap_or(true);
        let batch = parse_recurrence_batch_rational(&input)?;
        let polynomials = match collect_rational_polynomials_or_errors(batch) {
            Ok(polynomials) => polynomials,
            Err(parse_errors) => {
                return Ok(Json(recurrence_failure_response(
                    RecurrenceSearchStatus::InvalidInput,
                    "one or more polynomials failed to parse",
                    None,
                    parse_errors,
                )));
            }
        };
        if polynomials.len() < 3 {
            return Ok(Json(recurrence_failure_response(
                RecurrenceSearchStatus::InvalidInput,
                "need at least 3 polynomials",
                None,
                Vec::new(),
            )));
        }
        let (search, budget) = apply_recurrence_options(recurrence_options)?;
        match find_recurrence_adaptive_rational_with_budget(&polynomials, &search, budget) {
            AdaptiveSearchOutcome::Found(result) => Ok(Json(FindRecurrenceResponse {
                status: RecurrenceSearchStatus::Found,
                found: true,
                recurrence: Some(format!("{}", result.recurrence)),
                latex: Some(result.recurrence.to_latex()),
                mathematica: include_code.then(|| {
                    result
                        .recurrence
                        .to_mathematica_definition_rational(&polynomials)
                }),
                sage: include_code
                    .then(|| result.recurrence.to_sage_definition_rational(&polynomials)),
                python: include_code.then(|| {
                    result
                        .recurrence
                        .to_python_definition_rational(&polynomials)
                }),
                recurrence_json: include_code.then(|| {
                    recurrence_json_string(&result, &polynomials, &search, polynomials.len())
                }),
                unknowns: Some(result.num_unknowns),
                weighted_unknowns: Some(result.weighted_unknowns),
                equations: Some(result.num_equations),
                fit_polynomials: Some(result.fit_polynomials),
                verification_polynomials: Some(result.verification_polynomials),
                candidates_tried: Some(result.candidates_tried),
                candidates_considered: Some(result.diagnostics.considered_candidates),
                error: None,
                parse_errors: Vec::new(),
            })),
            AdaptiveSearchOutcome::NoRecurrence(summary) => Ok(Json(recurrence_failure_response(
                RecurrenceSearchStatus::NotFound,
                "no recurrence found within the search bounds",
                Some(summary),
                Vec::new(),
            ))),
            AdaptiveSearchOutcome::BudgetExhausted(summary) => {
                let error = format!(
                    "recurrence candidate budget exhausted after {} candidates; unsearched candidates remain",
                    summary.diagnostics.considered_candidates
                );
                Ok(Json(recurrence_failure_response(
                    RecurrenceSearchStatus::BudgetExhausted,
                    error,
                    Some(summary),
                    Vec::new(),
                )))
            }
        }
    }

    #[tool(
        description = "Generate or extend coefficient rows from recurrence JSON returned by find_recurrence."
    )]
    pub fn generate_recurrence_rows(
        &self,
        Parameters(input): Parameters<GenerateRecurrenceRowsRequest>,
    ) -> Result<Json<GenerateRecurrenceRowsResponse>, McpError> {
        if input.rows.is_some() == input.additional.is_some() {
            return Err(invalid_params("provide exactly one of rows or additional"));
        }
        validate_mcp_text(&input.recurrence_json, "recurrence JSON").map_err(invalid_params)?;
        let recurrence_json: RecurrenceJson = serde_json::from_str(&input.recurrence_json)
            .map_err(|error| invalid_params(format!("failed to parse recurrence JSON: {error}")))?;
        let (recurrence, first_index, initial_polys) = recurrence_json
            .to_recurrence_parts()
            .map_err(|error| invalid_params(format!("invalid recurrence JSON: {error}")))?;
        let row_count = match (input.rows, input.additional) {
            (Some(rows), None) => rows,
            (None, Some(additional)) => initial_polys
                .len()
                .checked_add(additional)
                .ok_or_else(|| invalid_params("recurrence row count overflow"))?,
            _ => unreachable!("validated exactly one row-count option"),
        };
        if row_count > MCP_MAX_RECURRENCE_ROWS {
            return Err(invalid_params(format!(
                "requested {row_count} recurrence rows; the MCP limit is {MCP_MAX_RECURRENCE_ROWS}"
            )));
        }
        first_index
            .checked_add(row_count)
            .ok_or_else(|| invalid_params("recurrence index range overflow"))?;
        validate_recurrence_for_generation(&recurrence, &initial_polys)?;
        let generated = recurrence
            .generate_rows_rational(&initial_polys, first_index, row_count)
            .map_err(|error| invalid_params(format!("failed to generate rows: {error}")))?;
        let polynomials = generated
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect::<Vec<Vec<_>>>();
        Ok(Json(GenerateRecurrenceRowsResponse {
            first_index,
            row_count: generated.len(),
            polynomials,
        }))
    }

    #[tool(description = "Compute the exact resultant of two polynomials.")]
    pub fn resultant(
        &self,
        Parameters(input): Parameters<BigIntPolynomialPairRequest>,
    ) -> Result<Json<ResultantResponse>, McpError> {
        let p = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.p, "p")?);
        let q = normalize_bigint_polynomial(parse_required_bigint_polynomial(&input.q, "q")?);
        let resultant =
            polytool::resultant_bigint_coeffs(&p.coefficients, &q.coefficients).to_string();
        Ok(Json(ResultantResponse {
            p: p.display,
            q: q.display,
            resultant,
        }))
    }

    #[tool(description = "Compute exact discriminants for one or more polynomials.")]
    pub fn discriminant(
        &self,
        Parameters(input): Parameters<BigIntPolynomialBatchInput>,
    ) -> Result<Json<DiscriminantResponse>, McpError> {
        let batch = parse_bigint_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => DiscriminantItem {
                    index,
                    ok: true,
                    error: None,
                    discriminant: Some(
                        polytool::discriminant_bigint_coeffs(&polynomial.coefficients).to_string(),
                    ),
                    polynomial: Some(polynomial.display),
                },
                Err(error) => DiscriminantItem {
                    index,
                    ok: false,
                    error: Some(error),
                    polynomial: None,
                    discriminant: None,
                },
            })
            .collect();
        Ok(Json(DiscriminantResponse { items }))
    }

    #[tool(
        description = "Convert h*-vectors to Ehrhart polynomials or Ehrhart polynomials to h*-vectors."
    )]
    pub fn ehrhart_hstar(
        &self,
        Parameters(input): Parameters<EhrhartHstarRequest>,
    ) -> Result<Json<EhrhartHstarResponse>, McpError> {
        match input.mode {
            EhrhartHstarMode::HstarToEhrhart => {
                let hstar_input = input
                    .hstar
                    .ok_or_else(|| invalid_params("`hstar` is required for hstar_to_ehrhart"))?;
                let hstar = parse_bigint_coefficients(&hstar_input)?;
                let ehrhart = hstar_to_ehrhart_bigint_coeffs(&hstar)
                    .into_iter()
                    .map(format_rational)
                    .collect();
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::HstarToEhrhart,
                    hstar: Some(bigint_strings(&hstar)),
                    ehrhart_coefficients: Some(ehrhart),
                }))
            }
            EhrhartHstarMode::EhrhartToHstar => {
                let hstar = match (
                    input.ehrhart_coefficients,
                    input.numerator_coefficients,
                    input.denominator,
                ) {
                    (Some(coefficients), None, None) => {
                        if coefficients.len() > MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL {
                            return Err(invalid_params(format!(
                                "Ehrhart polynomial has {} coefficients; the MCP limit is {MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}",
                                coefficients.len()
                            )));
                        }
                        for coefficient in &coefficients {
                            if coefficient.len() > MCP_MAX_COEFFICIENT_TEXT_BYTES {
                                return Err(invalid_params(format!(
                                    "Ehrhart coefficient is {} bytes; the MCP limit is {MCP_MAX_COEFFICIENT_TEXT_BYTES}",
                                    coefficient.len()
                                )));
                            }
                        }
                        let coefficients: Result<Vec<_>, _> =
                            coefficients.iter().map(|c| parse_rational(c)).collect();
                        ehrhart_to_hstar_bigint(&coefficients.map_err(invalid_params)?)
                            .map_err(|error| invalid_params(error.to_string()))?
                    }
                    (None, Some(numerator_coefficients), Some(denominator)) => {
                        let numerator_coefficients =
                            parse_bigint_coefficients(&numerator_coefficients)?;
                        let denominator =
                            parse_bigint_coefficient(&denominator).map_err(invalid_params)?;
                        if denominator == BigInt::from(0) {
                            return Err(invalid_params("`denominator` must be nonzero"));
                        }
                        let coefficients = numerator_coefficients
                            .into_iter()
                            .map(|numerator| BigRational::new(numerator, denominator.clone()))
                            .collect::<Vec<_>>();
                        ehrhart_to_hstar_bigint(&coefficients)
                            .map_err(|error| invalid_params(error.to_string()))?
                    }
                    _ => {
                        return Err(invalid_params(
                            "for ehrhart_to_hstar, provide either `ehrhart_coefficients` or both `numerator_coefficients` and `denominator`",
                        ));
                    }
                };
                Ok(Json(EhrhartHstarResponse {
                    mode: EhrhartHstarMode::EhrhartToHstar,
                    hstar: Some(bigint_strings(&hstar)),
                    ehrhart_coefficients: None,
                }))
            }
        }
    }

    #[tool(
        description = "Profile or check a cyclic sieving polynomial for one cyclic group order using exact cyclotomic arithmetic."
    )]
    pub fn cyclic_sieving(
        &self,
        Parameters(input): Parameters<CyclicSievingRequest>,
    ) -> Result<Json<CyclicSievingResponse>, McpError> {
        if input.order == 0 {
            return Err(invalid_params("`order` must be positive"));
        }
        let polynomial = normalize_bigint_polynomial(parse_required_bigint_polynomial(
            &input.polynomial,
            "polynomial",
        )?);
        let fixed_counts = input
            .fixed_counts
            .as_ref()
            .map(|counts| parse_bigint_coefficients(counts))
            .transpose()?;
        let report = cyclic_sieving_report_bigint(
            &polynomial.coefficients,
            input.order,
            fixed_counts.as_deref(),
        )
        .map_err(invalid_params)?;
        Ok(Json(CyclicSievingResponse {
            report: cyclic_sieving_report_item(&report),
        }))
    }

    #[tool(
        description = "For a polynomial sequence P_n(q), profile candidate cyclic sieving group orders n-2,n-1,n,n+1,n+2,n+3 by default."
    )]
    pub fn cyclic_sieving_sequence(
        &self,
        Parameters(input): Parameters<CyclicSievingSequenceRequest>,
    ) -> Result<Json<CyclicSievingSequenceResponse>, McpError> {
        let first_index = input.first_index.unwrap_or(0);
        let offsets = input.offsets.unwrap_or_else(|| vec![-2, -1, 0, 1, 2, 3]);
        if offsets.is_empty() {
            return Err(invalid_params("`offsets` may not be empty"));
        }
        let batch = parse_bigint_batch(&BigIntPolynomialBatchInput {
            polynomials: input.polynomials,
            text: input.text,
        })?;
        let mut polynomials = Vec::new();
        for (index, item) in batch.into_iter().enumerate() {
            match item {
                Ok(polynomial) => polynomials.push(polynomial.coefficients),
                Err(error) => {
                    return Err(invalid_params(format!("polynomial {index}: {error}")));
                }
            }
        }
        let mut fixed_counts = BTreeMap::new();
        for item in input.fixed_counts.unwrap_or_default() {
            if item.order == 0 {
                return Err(invalid_params("fixed-count order must be positive"));
            }
            let counts = parse_bigint_coefficients(&item.counts)?;
            if counts.len() != item.order {
                return Err(invalid_params(format!(
                    "fixed counts for index {} order {}: expected {} counts, got {}",
                    item.index,
                    item.order,
                    item.order,
                    counts.len()
                )));
            }
            fixed_counts.insert((item.index, item.order), counts);
        }
        let reports = cyclic_sieving_sequence_reports_bigint(
            &polynomials,
            first_index,
            &offsets,
            &fixed_counts,
        );
        Ok(Json(CyclicSievingSequenceResponse {
            first_index,
            offsets,
            items: reports.iter().map(cyclic_sequence_item_response).collect(),
        }))
    }

    #[tool(
        description = "Analyze symmetric decomposition, R-transform, interlacing checks, and magic-basis coordinates."
    )]
    pub fn analyze_decomposition(
        &self,
        Parameters(input): Parameters<PolynomialBatchInput>,
    ) -> Result<Json<DecompositionResponse>, McpError> {
        let batch = parse_batch(&input)?;
        let items = batch
            .into_iter()
            .enumerate()
            .map(|(index, item)| match item {
                Ok(polynomial) => {
                    match analyze_symmetric_decomposition_i64(&polynomial.coefficients) {
                        Ok(analysis) => {
                            let partial_sum_checks = analysis
                                .magic
                                .left_partial_sums
                                .iter()
                                .zip(analysis.magic.right_partial_sums.iter())
                                .map(|(left, right)| left <= right)
                                .collect();
                            DecompositionItem {
                                index,
                                ok: true,
                                report: Some(DecompositionReport {
                                    polynomial,
                                    reciprocal: displayed(&analysis.reciprocal),
                                    a: displayed(&analysis.a),
                                    b: displayed(&analysis.b),
                                    a_real_rooted: analysis.a_real_rooted,
                                    b_real_rooted: analysis.b_real_rooted,
                                    b_interlaces_a: analysis.b_interlaces_a,
                                    reciprocal_interlaces_input: analysis
                                        .reciprocal_interlaces_input,
                                    alternatingly_increasing: analysis.alternatingly_increasing,
                                    f_polynomial: displayed(&analysis.f_polynomial),
                                    r_transform_of_f: displayed(&analysis.r_transform_of_f),
                                    r_a: displayed(&analysis.r_a),
                                    r_b: displayed(&analysis.r_b),
                                    r_interlaces_f: analysis.r_interlaces_f,
                                    magic: MagicBasisReport {
                                        coordinates: analysis
                                            .magic
                                            .coordinates
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        left_partial_sums: analysis
                                            .magic
                                            .left_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        right_partial_sums: analysis
                                            .magic
                                            .right_partial_sums
                                            .iter()
                                            .map(ToString::to_string)
                                            .collect(),
                                        partial_sum_checks,
                                        all_nonnegative: analysis.magic.all_nonnegative,
                                        left_leq_right: analysis.magic.left_leq_right,
                                    },
                                }),
                                error: None,
                            }
                        }
                        Err(error) => DecompositionItem {
                            index,
                            ok: false,
                            report: None,
                            error: Some(error.to_string()),
                        },
                    }
                }
                Err(error) => DecompositionItem {
                    index,
                    ok: false,
                    report: None,
                    error: Some(error),
                },
            })
            .collect();
        Ok(Json(DecompositionResponse { items }))
    }

    #[tool(description = "List bundled recurrence-backed OEIS polynomial families.")]
    pub fn list_oeis_sequences(
        &self,
        Parameters(input): Parameters<ListOeisSequencesRequest>,
    ) -> Result<Json<ListOeisSequencesResponse>, McpError> {
        let include_experimental = input.include_experimental.unwrap_or(false);
        let sequences = polytool::oeis::catalog()
            .iter()
            .filter(|entry| {
                include_experimental
                    || entry.status != polytool::oeis::OeisSequenceStatus::Experimental
            })
            .map(oeis_sequence_summary)
            .collect::<Vec<_>>();
        Ok(Json(ListOeisSequencesResponse {
            count: sequences.len(),
            sequences,
        }))
    }

    #[tool(
        description = "Get metadata and exact recurrence exports for one bundled OEIS polynomial family."
    )]
    pub fn get_oeis_sequence(
        &self,
        Parameters(input): Parameters<GetOeisSequenceRequest>,
    ) -> Result<Json<GetOeisSequenceResponse>, McpError> {
        let entry = polytool::oeis::by_id(&input.id).ok_or_else(|| {
            invalid_params(format!("unknown bundled OEIS sequence: {}", input.id))
        })?;
        let (recurrence, initial_rows) = entry
            .recurrence_parts()
            .map_err(|error| invalid_params(error.to_string()))?;
        Ok(Json(GetOeisSequenceResponse {
            sequence: oeis_sequence_summary(entry),
            recurrence: recurrence.to_string(),
            latex: recurrence.to_latex(),
            mathematica: recurrence.to_mathematica_definition_rational(&initial_rows),
            sage: recurrence.to_sage_definition_rational(&initial_rows),
            python: recurrence.to_python_definition_rational(&initial_rows),
        }))
    }

    #[tool(
        description = "Generate exact polynomial coefficient rows for one bundled recurrence-backed OEIS family."
    )]
    pub fn generate_oeis_rows(
        &self,
        Parameters(input): Parameters<GenerateOeisRowsRequest>,
    ) -> Result<Json<GenerateOeisRowsResponse>, McpError> {
        const MAX_ROWS: usize = 200;
        if input.rows > MAX_ROWS {
            return Err(invalid_params(format!(
                "rows must be at most {MAX_ROWS} for MCP generation"
            )));
        }
        let entry = polytool::oeis::by_id(&input.id).ok_or_else(|| {
            invalid_params(format!("unknown bundled OEIS sequence: {}", input.id))
        })?;
        if entry.status == polytool::oeis::OeisSequenceStatus::Experimental
            && !input.include_experimental.unwrap_or(false)
        {
            return Err(invalid_params(format!(
                "{} is experimental; set include_experimental to true",
                entry.id
            )));
        }
        let first_row = input.first_row.unwrap_or(entry.first_row);
        let generated = entry
            .generate_rows_from(first_row, input.rows)
            .map_err(|error| invalid_params(error.to_string()))?;
        let rows = generated
            .into_iter()
            .enumerate()
            .map(|(offset, coefficients)| {
                let offset = i64::try_from(offset)
                    .map_err(|_| invalid_params("OEIS row offset does not fit in i64"))?;
                let n = first_row
                    .checked_add(offset)
                    .ok_or_else(|| invalid_params("OEIS row index overflow"))?;
                Ok(OeisPolynomialRow {
                    n,
                    polynomial: format_poly_bigint_coeffs(&coefficients),
                    coefficients: coefficients.iter().map(ToString::to_string).collect(),
                })
            })
            .collect::<Result<Vec<_>, McpError>>()?;
        Ok(Json(GenerateOeisRowsResponse {
            id: entry.id.to_string(),
            first_row,
            row_count: rows.len(),
            rows,
        }))
    }

    #[tool(description = "Generate standard polynomial sequences.")]
    pub fn generate_sequence(
        &self,
        Parameters(input): Parameters<GenerateSequenceRequest>,
    ) -> Result<Json<GenerateSequenceResponse>, McpError> {
        if input.max_n > MCP_MAX_SEQUENCE_N {
            return Err(invalid_params(format!(
                "max_n must be at most {MCP_MAX_SEQUENCE_N} for MCP sequence generation"
            )));
        }
        let polynomials = generated_sequence_polynomials_bigint(&input.sequence, input.max_n)
            .into_iter()
            .map(|coefficients| normalize_bigint_polynomial(coefficients).display)
            .collect();

        Ok(Json(GenerateSequenceResponse {
            sequence: input.sequence,
            max_n: input.max_n,
            polynomials,
        }))
    }
}

#[tool_handler(router = self.tool_router)]
impl ServerHandler for PolynomialToolsServer {
    fn get_info(&self) -> ServerInfo {
        ServerInfo {
            protocol_version: Default::default(),
            capabilities: ServerCapabilities::builder().enable_tools().build(),
            server_info: Implementation {
                name: "polytool".to_string(),
                title: Some("polytool".to_string()),
                version: env!("CARGO_PKG_VERSION").to_string(),
                description: Some(
                    "Exact univariate polytool for combinatorial research.".to_string(),
                ),
                icons: None,
                website_url: None,
            },
            instructions: Some("Exact univariate polytool for combinatorial research.".to_string()),
        }
    }
}

impl Default for PolynomialToolsServer {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn coeffs(values: &[i64]) -> PolynomialInput {
        PolynomialInput {
            coefficients: Some(values.to_vec()),
            expression: None,
        }
    }

    fn expr(value: &str) -> PolynomialInput {
        PolynomialInput {
            coefficients: None,
            expression: Some(value.to_string()),
        }
    }

    fn bigint_coeffs(values: &[&str]) -> BigIntPolynomialInput {
        BigIntPolynomialInput {
            coefficients: Some(
                values
                    .iter()
                    .map(|value| BigIntCoefficientInput::Text((*value).to_string()))
                    .collect(),
            ),
            expression: None,
        }
    }

    #[test]
    fn parses_coefficients_and_expression() {
        assert_eq!(
            normalize_polynomial(parse_polynomial_input(&coeffs(&[1, 2, 3, 0])).unwrap()),
            NormalizedPolynomial {
                polynomial: "1 + 2t + 3t^2".to_string(),
                coefficients: vec![1, 2, 3],
                degree: 2,
            }
        );
        assert_eq!(
            normalize_polynomial(parse_polynomial_input(&expr("3t^2 + 2t + 1")).unwrap())
                .coefficients,
            vec![1, 2, 3]
        );
    }

    #[test]
    fn preserves_batch_parse_errors() {
        let batch = parse_batch(&PolynomialBatchInput {
            polynomials: None,
            text: Some("1,2,1\nbad + x\n".to_string()),
        })
        .unwrap();
        let items = parse_items(&batch);
        assert_eq!(items.len(), 2);
        assert!(items[0].ok);
        assert!(!items[1].ok);
    }

    #[test]
    fn checks_eulerian_properties() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![bigint_coeffs(&["1", "11", "11", "1"])]),
                text: None,
            }))
            .unwrap();
        let item = &response.items[0];
        assert_eq!(item.real_rooted, Some(true));
        assert_eq!(item.gamma_positive, Some(true));
        assert_eq!(
            item.gamma_coefficients,
            Some(vec!["1".to_string(), "8".to_string()])
        );
        assert_eq!(item.unimodal, Some(true));
    }

    #[test]
    fn polynomial_properties_accepts_bigint_coefficients() {
        let server = PolynomialToolsServer::new();
        let large = "1000000000000000000000000000000000000000000";
        let Json(response) = server
            .polynomial_properties(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![bigint_coeffs(&[large, "2", "1"])]),
                text: None,
            }))
            .unwrap();
        let item = &response.items[0];
        assert!(item.ok);
        assert_eq!(
            item.polynomial
                .as_ref()
                .and_then(|polynomial| polynomial.coefficients.first()),
            Some(&large.to_string())
        );
        let value = serde_json::to_value(response).unwrap();
        assert_eq!(value["items"][0]["polynomial"]["coefficients"][0], large);
    }

    #[test]
    fn coefficient_tests_report_kurtz() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .coefficient_tests(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![bigint_coeffs(&["1000", "1110", "111", "1"])]),
                text: None,
            }))
            .unwrap();
        let item = &response.items[0];
        assert!(item.ok);
        assert_eq!(item.kurtz.as_ref().unwrap().holds, true);
        assert_eq!(item.kurtz.as_ref().unwrap().implies_real_rooted, true);
    }

    #[test]
    fn hstar_inequalities_report_named_failures() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .hstar_inequalities(Parameters(HStarInequalitiesRequest {
                polynomials: Some(vec![bigint_coeffs(&["1", "20", "1", "0"])]),
                text: None,
                dimension: Some(3),
            }))
            .unwrap();
        let item = &response.items[0];
        assert!(item.ok);
        assert_eq!(item.all_applicable_hold, Some(false));
        assert!(item.checks.iter().any(|check| {
            check.family == "Balletti-Higashitani" && check.applicable && !check.holds
        }));
    }

    #[test]
    fn cyclic_sieving_checks_order_two() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .cyclic_sieving(Parameters(CyclicSievingRequest {
                polynomial: bigint_coeffs(&["1", "1"]),
                order: 2,
                fixed_counts: Some(vec![1.into(), 0.into()]),
            }))
            .unwrap();
        assert_eq!(response.report.holds, Some(false));

        let Json(response) = server
            .cyclic_sieving(Parameters(CyclicSievingRequest {
                polynomial: bigint_coeffs(&["1", "1"]),
                order: 2,
                fixed_counts: Some(vec![2.into(), 0.into()]),
            }))
            .unwrap();
        assert_eq!(response.report.holds, Some(true));
    }

    #[test]
    fn cyclic_sieving_sequence_profiles_default_offsets() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .cyclic_sieving_sequence(Parameters(CyclicSievingSequenceRequest {
                polynomials: Some(vec![
                    bigint_coeffs(&["1", "1"]),
                    bigint_coeffs(&["1", "1", "1"]),
                ]),
                text: None,
                first_index: Some(2),
                offsets: None,
                fixed_counts: None,
            }))
            .unwrap();
        assert_eq!(response.items[0].index, 2);
        assert!(response.items[0]
            .candidate_orders
            .iter()
            .any(|report| report.order == 2));
        assert!(response.items[0]
            .candidate_orders
            .iter()
            .any(|report| report.order == 5));
    }

    #[test]
    fn checks_strict_and_weak_interlacing() {
        let server = PolynomialToolsServer::new();
        let Json(strict) = server
            .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
                p: bigint_coeffs(&["-2", "1"]),
                q: bigint_coeffs(&["3", "-4", "1"]),
            }))
            .unwrap();
        assert_eq!(strict.strict, Some(true));

        let Json(weak) = server
            .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
                p: bigint_coeffs(&["-1", "1"]),
                q: bigint_coeffs(&["2", "-3", "1"]),
            }))
            .unwrap();
        assert_eq!(weak.weak, Some(true));
    }

    #[test]
    fn interlacing_accepts_bigint_coefficients() {
        let server = PolynomialToolsServer::new();
        let center = "100000000000000000000";
        let Json(result) = server
            .check_interlacing_pair(Parameters(BigIntPolynomialPairRequest {
                p: bigint_coeffs(&["-100000000000000000000", "1"]),
                q: bigint_coeffs(&[
                    "9999999999999999999999999999999999999999",
                    "-200000000000000000000",
                    "1",
                ]),
            }))
            .unwrap();

        assert_eq!(result.strict, Some(true));
        assert_eq!(result.p.coefficients[0], format!("-{center}"));
        assert_eq!(result.q.degree, 2);
    }

    #[test]
    fn interlacing_profile_stops_at_first_failure() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .check_interlacing_profile(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![
                    bigint_coeffs(&["1"]),
                    bigint_coeffs(&["1", "1"]),
                    bigint_coeffs(&["1", "2", "1"]),
                    bigint_coeffs(&["1", "0", "1"]),
                    bigint_coeffs(&["1", "3", "3", "1"]),
                ]),
                text: None,
            }))
            .unwrap();

        assert_eq!(response.profile[2].previous_count, 2);
        assert_eq!(response.profile[2].checked_previous_count, 2);
        assert_eq!(response.profile[2].interlacing_previous_count, 1);
        assert_eq!(response.profile[2].previous_interlacing_indices, vec![1]);

        assert_eq!(response.profile[3].previous_count, 3);
        assert_eq!(response.profile[3].checked_previous_count, 1);
        assert_eq!(response.profile[3].interlacing_previous_count, 0);
    }

    #[test]
    fn generate_sequence_returns_bigint_coefficients() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .generate_sequence(Parameters(GenerateSequenceRequest {
                sequence: SequenceKind::ChebyshevT,
                max_n: 64,
            }))
            .unwrap();

        let last = response.polynomials.last().expect("T_64");
        assert_eq!(last.degree, 64);
        assert_eq!(
            last.coefficients.last().map(String::as_str),
            Some("9223372036854775808")
        );
    }

    #[test]
    fn lists_and_generates_oeis_catalog_rows() {
        let server = PolynomialToolsServer::new();
        let Json(list) = server
            .list_oeis_sequences(Parameters(ListOeisSequencesRequest::default()))
            .unwrap();
        assert_eq!(list.count, 755);
        assert!(list.sequences.iter().any(|entry| entry.id == "A008292"));
        assert!(list
            .sequences
            .iter()
            .any(|entry| entry.id == "A390883" && entry.status == "validated"));

        let Json(response) = server
            .generate_oeis_rows(Parameters(GenerateOeisRowsRequest {
                id: "a008292".to_string(),
                rows: 3,
                first_row: None,
                include_experimental: None,
            }))
            .unwrap();
        assert_eq!(response.first_row, 1);
        assert_eq!(response.rows[2].n, 3);
        assert_eq!(response.rows[2].coefficients, vec!["1", "4", "1"]);
    }

    #[test]
    fn computes_resultant_and_discriminant() {
        let server = PolynomialToolsServer::new();
        let Json(resultant) = server
            .resultant(Parameters(BigIntPolynomialPairRequest {
                p: bigint_coeffs(&["2", "-3", "1"]),
                q: bigint_coeffs(&["-3", "1"]),
            }))
            .unwrap();
        assert_eq!(resultant.resultant, "2");

        let Json(discriminant) = server
            .discriminant(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![bigint_coeffs(&["-1", "0", "1"])]),
                text: None,
            }))
            .unwrap();
        assert_eq!(discriminant.items[0].discriminant.as_deref(), Some("4"));

        let large = "1000000000000000000000000000000";
        let Json(large_discriminant) = server
            .discriminant(Parameters(BigIntPolynomialBatchInput {
                polynomials: Some(vec![bigint_coeffs(&["1", "0", large])]),
                text: None,
            }))
            .unwrap();
        assert_eq!(
            large_discriminant.items[0].discriminant.as_deref(),
            Some("-4000000000000000000000000000000")
        );
    }

    #[test]
    fn roundtrips_hstar_and_ehrhart() {
        let server = PolynomialToolsServer::new();
        let Json(to_ehrhart) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::HstarToEhrhart,
                hstar: Some(vec![1.into(), 1.into(), 0.into()]),
                ehrhart_coefficients: None,
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        assert_eq!(
            to_ehrhart.ehrhart_coefficients.as_ref().unwrap(),
            &vec!["1".to_string(), "2".to_string(), "1".to_string()]
        );

        let Json(to_hstar) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::EhrhartToHstar,
                hstar: None,
                ehrhart_coefficients: Some(vec!["1".to_string(), "2".to_string(), "1".to_string()]),
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        assert_eq!(
            to_hstar.hstar,
            Some(vec!["1".to_string(), "1".to_string(), "0".to_string()])
        );

        let Json(common_denominator) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::EhrhartToHstar,
                hstar: None,
                ehrhart_coefficients: None,
                numerator_coefficients: Some(vec![2.into(), 3.into(), 1.into()]),
                denominator: Some(2.into()),
            }))
            .unwrap();
        assert_eq!(
            common_denominator.hstar,
            Some(vec!["1".to_string(), "0".to_string(), "0".to_string()])
        );
    }

    #[test]
    fn find_recurrence_schema_is_flat_and_documents_every_control() {
        let schema = serde_json::to_value(schemars::schema_for!(FindRecurrenceRequest)).unwrap();
        assert_eq!(schema["type"], "object");
        assert!(schema.get("oneOf").is_none());
        assert_eq!(schema["additionalProperties"], false);
        let properties = schema["properties"].as_object().expect("schema properties");
        for name in [
            "polynomials",
            "coefficients",
            "expressions",
            "text",
            "options",
            "include_code",
            "skip_prefix",
            "min_rec_len",
            "max_rec_len",
            "min_var_deg",
            "max_var_deg",
            "min_idx_deg",
            "max_idx_deg",
            "min_diff_deg",
            "max_diff_deg",
            "try_inhomogeneous",
            "min_inhomo_var_deg",
            "max_inhomo_var_deg",
            "min_inhomo_idx_deg",
            "max_inhomo_idx_deg",
            "try_denominator",
            "try_alternating_sign",
            "max_denom_var_deg",
            "max_denom_idx_deg",
            "min_margin",
            "no_verify",
            "fit_extra_rows",
            "modular_prefilter",
            "max_candidates",
        ] {
            let property = properties
                .get(name)
                .unwrap_or_else(|| panic!("missing {name}"));
            assert!(
                property["description"]
                    .as_str()
                    .is_some_and(|text| !text.is_empty()),
                "{name} lacks a description: {property}"
            );
        }
    }

    #[test]
    fn flat_and_legacy_recurrence_controls_work_with_top_level_precedence() {
        let server = PolynomialToolsServer::new();
        let flat: FindRecurrenceRequest = serde_json::from_value(json!({
            "coefficients": [[1], [2], [4], [8], [16]],
            "min_rec_len": 1,
            "max_rec_len": 1,
            "max_var_deg": 0,
            "max_idx_deg": 0,
            "max_diff_deg": 0,
            "max_candidates": 1
        }))
        .unwrap();
        let Json(flat_result) = server.find_recurrence(Parameters(flat)).unwrap();
        assert_eq!(flat_result.status, RecurrenceSearchStatus::Found);
        assert_eq!(flat_result.candidates_considered, Some(1));

        let legacy: FindRecurrenceRequest = serde_json::from_value(json!({
            "coefficients": [[1], [2], [4], [8], [16]],
            "options": {
                "min_rec_len": 1,
                "max_rec_len": 1,
                "max_var_deg": 0,
                "max_idx_deg": 0,
                "max_diff_deg": 0,
                "max_candidates": 1
            }
        }))
        .unwrap();
        let Json(legacy_result) = server.find_recurrence(Parameters(legacy)).unwrap();
        assert_eq!(legacy_result.status, RecurrenceSearchStatus::Found);

        let conflict: FindRecurrenceRequest = serde_json::from_value(json!({
            "coefficients": [[1], [2], [4], [8], [16]],
            "max_candidates": 1,
            "max_rec_len": 1,
            "options": { "max_candidates": 0, "max_rec_len": 7 }
        }))
        .unwrap();
        let effective = conflict.effective_options().unwrap();
        assert_eq!(effective.max_candidates, Some(1));
        assert_eq!(effective.max_rec_len, Some(1));
        let Json(conflict_result) = server.find_recurrence(Parameters(conflict)).unwrap();
        assert_eq!(conflict_result.status, RecurrenceSearchStatus::Found);
        assert_eq!(conflict_result.candidates_considered, Some(1));
    }

    #[test]
    fn compact_recurrence_result_omits_code_and_default_result_remains_full() {
        let server = PolynomialToolsServer::new();
        let request = || {
            json!({
                "coefficients": [[1], [2], [4], [8], [16]],
                "max_rec_len": 1,
                "max_var_deg": 0,
                "max_idx_deg": 0,
                "max_diff_deg": 0
            })
        };

        let full: FindRecurrenceRequest = serde_json::from_value(request()).unwrap();
        let Json(full) = server.find_recurrence(Parameters(full)).unwrap();
        assert!(full.recurrence.is_some());
        assert!(full.latex.is_some());
        assert!(full.mathematica.is_some());
        assert!(full.python.is_some());
        assert!(full.sage.is_some());
        assert!(full.recurrence_json.is_some());

        let mut compact_request = request();
        compact_request["include_code"] = json!(false);
        let compact: FindRecurrenceRequest = serde_json::from_value(compact_request).unwrap();
        let Json(compact) = server.find_recurrence(Parameters(compact)).unwrap();
        assert!(compact.recurrence.is_some());
        assert!(compact.latex.is_some());
        assert!(compact.unknowns.is_some());
        assert!(compact.candidates_considered.is_some());
        assert!(compact.mathematica.is_none());
        assert!(compact.python.is_none());
        assert!(compact.sage.is_none());
        assert!(compact.recurrence_json.is_none());
        let serialized = serde_json::to_value(compact).unwrap();
        for absent in ["mathematica", "python", "sage", "recurrence_json"] {
            assert!(serialized.get(absent).is_none(), "unexpected {absent}");
        }
    }

    #[test]
    fn recurrence_input_exclusivity_remains_runtime_validated() {
        let empty = FindRecurrenceRequest::default();
        assert!(
            format!("{:?}", parse_recurrence_batch_rational(&empty).unwrap_err())
                .contains("expected exactly one")
        );

        let conflict: FindRecurrenceRequest = serde_json::from_value(json!({
            "coefficients": [[1], [2], [4]],
            "expressions": ["1", "2", "4"]
        }))
        .unwrap();
        assert!(format!(
            "{:?}",
            parse_recurrence_batch_rational(&conflict).unwrap_err()
        )
        .contains("expected exactly one"));
    }

    #[test]
    fn recurrence_search_examples() {
        let server = PolynomialToolsServer::new();
        let Json(geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: Some(vec![
                    coeffs(&[1]),
                    coeffs(&[2]),
                    coeffs(&[4]),
                    coeffs(&[8]),
                    coeffs(&[16]),
                ]),
                coefficients: None,
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(geometric.recurrence.as_deref(), Some("P(n) = 2 P(n-1)"));
        assert!(geometric
            .python
            .as_deref()
            .is_some_and(|python| python.contains("def P(n):")));
        let recurrence_json = geometric
            .recurrence_json
            .clone()
            .expect("recurrence JSON output");
        let Json(generated) = server
            .generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
                recurrence_json,
                rows: Some(6),
                additional: None,
            }))
            .unwrap();
        assert_eq!(
            generated.polynomials,
            vec![
                vec!["1".to_string()],
                vec!["2".to_string()],
                vec!["4".to_string()],
                vec!["8".to_string()],
                vec!["16".to_string()],
                vec!["32".to_string()],
            ]
        );

        let Json(fibonacci) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![1.into()],
                    vec![1.into()],
                    vec![2.into()],
                    vec![3.into()],
                    vec![5.into()],
                    vec![8.into()],
                ]),
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(
            fibonacci.recurrence.as_deref(),
            Some("P(n) = P(n-1) + P(n-2)")
        );

        let Json(expression_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: None,
                expressions: Some(vec![
                    "1".to_string(),
                    "2".to_string(),
                    "4".to_string(),
                    "8".to_string(),
                    "16".to_string(),
                ]),
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(
            expression_geometric.recurrence.as_deref(),
            Some("P(n) = 2 P(n-1)")
        );

        let Json(rational_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![RationalCoefficientInput::Text("1/2".to_string())],
                    vec![RationalCoefficientInput::Text("1/4".to_string())],
                    vec![RationalCoefficientInput::Text("1/8".to_string())],
                    vec![RationalCoefficientInput::Text("1/16".to_string())],
                    vec![RationalCoefficientInput::Text("1/32".to_string())],
                ]),
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(
            rational_geometric.recurrence.as_deref(),
            Some("P(n) = 1/2 P(n-1)")
        );

        let large = "1267650600228229401496703205376";
        let double_large = "2535301200456458802993406410752";
        let Json(large_geometric) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![RationalCoefficientInput::Text(large.to_string())],
                    vec![RationalCoefficientInput::Text(double_large.to_string())],
                    vec![RationalCoefficientInput::Text(
                        "5070602400912917605986812821504".to_string(),
                    )],
                    vec![RationalCoefficientInput::Text(
                        "10141204801825835211973625643008".to_string(),
                    )],
                    vec![RationalCoefficientInput::Text(
                        "20282409603651670423947251286016".to_string(),
                    )],
                ]),
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(
            large_geometric.recurrence.as_deref(),
            Some("P(n) = 2 P(n-1)")
        );
        assert!(large_geometric
            .sage
            .as_deref()
            .is_some_and(|sage| sage.contains(large)));

        let Json(alternating) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![1.into()],
                    vec![1.into()],
                    vec![(-1).into()],
                    vec![(-1).into()],
                    vec![1.into()],
                    vec![1.into()],
                    vec![(-1).into()],
                    vec![(-1).into()],
                ]),
                expressions: None,
                text: None,
                options: Some(
                    serde_json::from_value(json!({
                        "try_alternating_sign": true,
                        "max_rec_len": 1,
                        "max_var_deg": 0,
                        "max_idx_deg": 0,
                        "max_diff_deg": 0
                    }))
                    .unwrap(),
                ),
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(
            alternating.recurrence.as_deref(),
            Some("P(n) = (-1)^n P(n-1)")
        );
        assert!(alternating
            .sage
            .as_deref()
            .is_some_and(|sage| sage.contains("(-1)**n*P(n - 1)")));

        let eulerian = vec![
            coeffs(&[1]),
            coeffs(&[1]),
            coeffs(&[1, 1]),
            coeffs(&[1, 4, 1]),
            coeffs(&[1, 11, 11, 1]),
            coeffs(&[1, 26, 66, 26, 1]),
        ];
        let Json(eulerian_result) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: Some(eulerian),
                coefficients: None,
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        assert!(eulerian_result.found);

        let conflict = parse_recurrence_batch_rational(&FindRecurrenceRequest {
            polynomials: Some(vec![coeffs(&[1]), coeffs(&[2]), coeffs(&[4])]),
            coefficients: Some(vec![vec![1.into()], vec![2.into()], vec![4.into()]]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        })
        .unwrap_err();
        assert!(format!("{conflict:?}").contains("expected exactly one"));
    }

    #[test]
    fn recurrence_candidate_budget_has_distinct_mcp_outcomes() {
        let server = PolynomialToolsServer::new();
        let fibonacci = || FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![1.into()],
                vec![1.into()],
                vec![2.into()],
                vec![3.into()],
                vec![5.into()],
            ]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        };

        let mut zero_budget = fibonacci();
        zero_budget.options = Some(
            serde_json::from_value(json!({ "max_candidates": 0 }))
                .expect("valid recurrence options"),
        );
        let Json(exhausted) = server.find_recurrence(Parameters(zero_budget)).unwrap();
        assert_eq!(exhausted.status, RecurrenceSearchStatus::BudgetExhausted);
        assert!(!exhausted.found);
        assert_eq!(exhausted.candidates_considered, Some(0));
        assert_eq!(exhausted.candidates_tried, Some(0));

        let mut exact_space = fibonacci();
        exact_space.options = Some(
            serde_json::from_value(json!({
                "max_candidates": 1,
                "min_rec_len": 1,
                "max_rec_len": 1,
                "max_var_deg": 0,
                "max_idx_deg": 0,
                "max_diff_deg": 0
            }))
            .expect("valid recurrence options"),
        );
        let Json(not_found) = server.find_recurrence(Parameters(exact_space)).unwrap();
        assert_eq!(not_found.status, RecurrenceSearchStatus::NotFound);
        assert_eq!(not_found.candidates_considered, Some(1));

        let Json(found) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: None,
                coefficients: Some(vec![
                    vec![1.into()],
                    vec![2.into()],
                    vec![4.into()],
                    vec![8.into()],
                    vec![16.into()],
                ]),
                expressions: None,
                text: None,
                options: Some(
                    serde_json::from_value(json!({ "max_candidates": 1 }))
                        .expect("valid recurrence options"),
                ),
                ..Default::default()
            }))
            .unwrap();
        assert_eq!(found.status, RecurrenceSearchStatus::Found);
        assert_eq!(found.candidates_considered, Some(1));
    }

    #[test]
    fn family_check_reports_properties_and_recurrence() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
                polynomials: Some(vec![
                    coeffs(&[1]),
                    coeffs(&[2]),
                    coeffs(&[4]),
                    coeffs(&[8]),
                    coeffs(&[16]),
                ]),
                text: None,
                sequence: None,
                max_n: None,
                options: Some(FamilyCheckOptions {
                    require_real_rooted: Some(true),
                    require_simple_roots: None,
                    require_palindromic: None,
                    require_gamma_positive: None,
                    require_unimodal: None,
                    require_log_concave: None,
                    require_ultra_log_concave: None,
                    check_consecutive_interlacing: Some(false),
                    require_consecutive_weak_interlacing: None,
                    find_recurrence: Some(true),
                }),
                recurrence_options: None,
                lace: None,
            }))
            .unwrap();

        assert!(response.all_required_checks_passed);
        assert_eq!(response.item_count, 5);
        assert_eq!(
            response
                .recurrence
                .as_ref()
                .and_then(|r| r.recurrence.as_deref()),
            Some("P(n) = 2 P(n-1)")
        );
        assert!(response.markdown.contains("Polynomial properties"));
    }

    #[test]
    fn family_check_can_require_finite_lace_tnn() {
        let server = PolynomialToolsServer::new();
        let Json(response) = server
            .check_polynomial_family(Parameters(CheckPolynomialFamilyRequest {
                polynomials: Some(vec![
                    coeffs(&[2, 1]),
                    coeffs(&[8, 6, 1]),
                    coeffs(&[3, 4, 1]),
                ]),
                text: None,
                sequence: None,
                max_n: None,
                options: None,
                recurrence_options: None,
                lace: Some(LaceCheckRequest {
                    block_rows: Some(1),
                    block_cols: Some(3),
                    max_minor_size: Some(3),
                    include_matrix: Some(true),
                }),
            }))
            .unwrap();

        assert!(!response.all_required_checks_passed);
        assert_eq!(
            response.first_failure.as_deref(),
            Some("finite Lace truncation is not TNN up to minors of size 3")
        );
        assert_eq!(response.lace.as_ref().map(|l| l.tnn), Some(false));
        assert!(response
            .lace
            .as_ref()
            .and_then(|l| l.error.as_deref())
            .is_some_and(|error| error.contains("det = -1")));
    }

    #[test]
    fn serializes_bigint_and_bigrational_as_strings() {
        let server = PolynomialToolsServer::new();
        let Json(resultant) = server
            .resultant(Parameters(BigIntPolynomialPairRequest {
                p: bigint_coeffs(&["1", "0", "1"]),
                q: bigint_coeffs(&["-1", "1"]),
            }))
            .unwrap();
        let value = serde_json::to_value(resultant).unwrap();
        assert_eq!(value["resultant"], "2");
        assert_eq!(value["p"]["coefficients"][0], "1");

        let Json(ehrhart) = server
            .ehrhart_hstar(Parameters(EhrhartHstarRequest {
                mode: EhrhartHstarMode::HstarToEhrhart,
                hstar: Some(vec![1.into(), 0.into(), 0.into()]),
                ehrhart_coefficients: None,
                numerator_coefficients: None,
                denominator: None,
            }))
            .unwrap();
        let value = serde_json::to_value(ehrhart).unwrap();
        assert_eq!(value["ehrhart_coefficients"][1], "3/2");
    }

    #[test]
    fn mcp_rejects_oversized_polynomial_and_recurrence_inputs() {
        let oversized_expression = PolynomialInput {
            coefficients: None,
            expression: Some(format!("x^{MCP_MAX_COEFFICIENTS_PER_POLYNOMIAL}")),
        };
        assert!(parse_polynomial_input(&oversized_expression)
            .unwrap_err()
            .contains("MCP limit"));

        let oversized_batch = PolynomialBatchInput {
            polynomials: Some(vec![coeffs(&[1]); MCP_MAX_BATCH_POLYNOMIALS + 1]),
            text: None,
        };
        assert!(format!("{:?}", parse_batch(&oversized_batch).unwrap_err()).contains("MCP limit"));

        let oversized_recurrence_batch = FindRecurrenceRequest {
            polynomials: None,
            coefficients: Some(vec![
                vec![RationalCoefficientInput::from(1)];
                MCP_MAX_RECURRENCE_INPUT_ROWS + 1
            ]),
            expressions: None,
            text: None,
            options: None,
            ..Default::default()
        };
        assert!(format!(
            "{:?}",
            parse_recurrence_batch_rational(&oversized_recurrence_batch).unwrap_err()
        )
        .contains("MCP limit"));

        let excessive_length: RecurrenceSearchOptionsInput =
            serde_json::from_value(json!({ "max_rec_len": MCP_MAX_RECURRENCE_LENGTH + 1 }))
                .unwrap();
        assert!(format!(
            "{:?}",
            apply_recurrence_options(Some(excessive_length)).unwrap_err()
        )
        .contains("maximum recurrence length"));

        let excessive_grid: RecurrenceSearchOptionsInput = serde_json::from_value(json!({
            "max_rec_len": MCP_MAX_RECURRENCE_LENGTH,
            "max_var_deg": MCP_MAX_RECURRENCE_DEGREE,
            "max_idx_deg": MCP_MAX_RECURRENCE_DEGREE,
            "max_diff_deg": MCP_MAX_RECURRENCE_DEGREE
        }))
        .unwrap();
        assert!(format!(
            "{:?}",
            apply_recurrence_options(Some(excessive_grid)).unwrap_err()
        )
        .contains("candidates"));
    }

    #[test]
    fn mcp_bounds_generation_and_reports_exact_conversion_errors() {
        let server = PolynomialToolsServer::new();
        assert!(server
            .generate_sequence(Parameters(GenerateSequenceRequest {
                sequence: SequenceKind::Eulerian,
                max_n: MCP_MAX_SEQUENCE_N + 1,
            }))
            .is_err());

        let nonintegral = server.ehrhart_hstar(Parameters(EhrhartHstarRequest {
            mode: EhrhartHstarMode::EhrhartToHstar,
            hstar: None,
            ehrhart_coefficients: Some(vec!["1/2".to_string()]),
            numerator_coefficients: None,
            denominator: None,
        }));
        assert!(
            format!("{:?}", nonintegral.err().expect("nonintegral error"))
                .contains("not an integer")
        );

        let Json(found) = server
            .find_recurrence(Parameters(FindRecurrenceRequest {
                polynomials: Some(vec![coeffs(&[1]), coeffs(&[2]), coeffs(&[4]), coeffs(&[8])]),
                coefficients: None,
                expressions: None,
                text: None,
                options: None,
                ..Default::default()
            }))
            .unwrap();
        let recurrence_json = found.recurrence_json.expect("small recurrence found");
        let too_many_rows =
            server.generate_recurrence_rows(Parameters(GenerateRecurrenceRowsRequest {
                recurrence_json,
                rows: Some(MCP_MAX_RECURRENCE_ROWS + 1),
                additional: None,
            }));
        assert!(
            format!("{:?}", too_many_rows.err().expect("row limit error")).contains("MCP limit")
        );

        let oversized_lace = check_family_lace(
            &vec![normalize_polynomial(vec![1]); 4],
            &LaceCheckRequest {
                block_rows: Some(MCP_MAX_LACE_MATRIX_CELLS),
                block_cols: Some(2),
                max_minor_size: Some(1),
                include_matrix: Some(false),
            },
        );
        assert!(!oversized_lace.tnn);
        assert!(oversized_lace
            .error
            .is_some_and(|error| error.contains("exceeds MCP limits")));
    }
}
