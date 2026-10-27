//! Tool request types and their JSON schemas.
//!
//! Coefficient inputs accept JSON integers or decimal strings; hand-written
//! `JsonSchema` impls document the accepted forms.

use schemars::{JsonSchema, Schema, SchemaGenerator};
use serde::{de, Deserialize, Deserializer, Serialize};
use serde_json::json;
use std::borrow::Cow;

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

#[derive(Debug, Clone, Serialize)]
#[serde(untagged)]
pub enum BigIntCoefficientInput {
    Integer(i64),
    Text(String),
}

impl<'de> Deserialize<'de> for BigIntCoefficientInput {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct BigIntCoefficientVisitor;

        impl de::Visitor<'_> for BigIntCoefficientVisitor {
            type Value = BigIntCoefficientInput;

            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("an integer or an exact integer encoded as a decimal string")
            }

            fn visit_i64<E>(self, value: i64) -> Result<Self::Value, E> {
                Ok(BigIntCoefficientInput::Integer(value))
            }

            fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E> {
                Ok(i64::try_from(value).map_or_else(
                    |_| BigIntCoefficientInput::Text(value.to_string()),
                    BigIntCoefficientInput::Integer,
                ))
            }

            fn visit_f64<E>(self, value: f64) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Err(E::custom(format!(
                    "integer coefficient was received as the floating-point JSON number {value}; \
                     send exact coefficients as quoted decimal strings, especially outside the \
                     JavaScript safe-integer range [-9007199254740991, 9007199254740991]"
                )))
            }

            fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                Ok(BigIntCoefficientInput::Text(value.to_string()))
            }

            fn visit_string<E>(self, value: String) -> Result<Self::Value, E> {
                Ok(BigIntCoefficientInput::Text(value))
            }
        }

        deserializer.deserialize_any(BigIntCoefficientVisitor)
    }
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
                    "description": "Exact integer coefficient. Use decimal strings for values outside a client's safe integer range, e.g. \"42\" or \"1267650600228229401496703205376\"."
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
    pub(crate) fn effective_options(&self) -> Option<RecurrenceSearchOptionsInput> {
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
    /// Fixed-point counts `fix(g^d)` for `d = 0..order-1`.
    pub fixed_counts: Option<Vec<BigIntCoefficientInput>>,
    /// Alternative to `fixed_counts`: the multiset of orbit sizes (each a
    /// divisor of `order`) of a concrete cyclic action.
    pub orbit_sizes: Option<Vec<usize>>,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BicyclicSievingRequest {
    /// Coefficient matrix: `matrix[i][j]` is the coefficient of `q^i t^j`.
    pub matrix: Vec<Vec<BigIntCoefficientInput>>,
    /// Order `n1` of the cyclic group acting through `q`.
    pub order1: usize,
    /// Order `n2` of the cyclic group acting through `t`.
    pub order2: usize,
    /// Optional `n1 x n2` table of fixed-point counts `fix(g1^d1 g2^d2)`.
    pub fixed_table: Option<Vec<Vec<BigIntCoefficientInput>>>,
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

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GenerateRecurrenceRowsRequest {
    pub recurrence_json: String,
    pub rows: Option<usize>,
    pub additional: Option<usize>,
}

/// A supplied `polytool.recurrence.v1` record to check against the
/// interlacing theorems.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTheoremsRequest {
    pub recurrence_json: String,
}

/// Structured recurrence-space request; the library validates the strict schema.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceSpacesToolRequest {
    pub request: serde_json::Value,
}

#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VerifyRecurrenceSpacesToolRequest {
    pub report: serde_json::Value,
    /// Optional independently supplied assumptions. Otherwise replay checks
    /// consistency with the assumptions embedded in the report.
    pub expected_request: Option<serde_json::Value>,
}

/// Structured coordinate-transformation request
/// (`polytool::recurrence::RecurrenceTransformRequest`); the library validates
/// the strict schema.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTransformToolRequest {
    pub request: serde_json::Value,
}

/// A `polytool.recurrence-transform.v1` report to replay exactly.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VerifyRecurrenceTransformToolRequest {
    pub report: serde_json::Value,
}

/// Exact root-enclosure request
/// (`polytool::recurrence::RootEnclosureRequest`); the library validates the
/// strict schema.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RootEnclosuresToolRequest {
    pub request: serde_json::Value,
}

/// Structured implicit-transport request
/// (`polytool::recurrence::ImplicitTransformRequest`); the library validates
/// the strict schema.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceImplicitToolRequest {
    pub request: serde_json::Value,
}

/// A `polytool.recurrence-implicit.v1` report to replay exactly.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VerifyRecurrenceImplicitToolRequest {
    pub report: serde_json::Value,
}

/// A generating description (`polytool::recurrence::GeneratingDescription`)
/// to assess.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct AssessGeneratingDescriptionToolRequest {
    pub description: serde_json::Value,
}

/// Input of `check_tnn_network`.  Give exactly one of `polynomials`, `text`
/// or `matrix`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct TnnNetworkRequest {
    /// Monic polynomial rows `P_0, P_1, ...` with `deg P_n = n` (ascending
    /// coefficients); `P_0 = 1` is prepended when the first row has degree 1.
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    /// The same rows as newline-separated coefficient lists or expressions.
    pub text: Option<String>,
    /// Matrix rows instead of polynomials: a square lower-triangular matrix
    /// (every row as long as the number of rows) or triangular rows (row `n`
    /// has `n + 1` entries).  Shape errors then name the offending entry.
    pub matrix: Option<Vec<Vec<BigIntCoefficientInput>>>,
    /// Also compute and check the chain polynomials `p_0 = 1`,
    /// `p_n = t * sum_{k<n} r_{n,k} p_k` (default `false`).
    pub chain: Option<bool>,
}

/// Edge-weight domain of `recognize_brenti_strip_digraph`.
#[derive(Debug, Clone, Copy, Default, Deserialize, Serialize, JsonSchema, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum BrentiWeightsInput {
    /// Nonnegative rational weights (default).
    #[default]
    Rational,
    /// Nonnegative integer weights.
    Integer,
}

/// Input of `recognize_brenti_strip_digraph`.  Give exactly one of
/// `polynomials` or `text`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct BrentiRecognitionRequest {
    /// Rows `P_1, P_2, ...` (or `P_0, P_1, ...` with `first_row_index = 0`),
    /// ascending coefficients.
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    /// The same rows as newline-separated coefficient lists or expressions.
    pub text: Option<String>,
    /// Digraph row of the first polynomial: 1 (default, `P_0 = 1` implicit)
    /// or 0 (the first row must be `1`).
    pub first_row_index: Option<usize>,
    /// At most this many consecutive target heights per height (bounds the
    /// out-degree; default: no bound).
    pub max_target_span: Option<usize>,
    /// Rational (default) or integer edge weights.
    pub weights: Option<BrentiWeightsInput>,
    /// Budget on the partial supports visited (default and maximum
    /// 1000000).
    pub max_search_nodes: Option<usize>,
    /// Budget on the integer search for one column (default 10000, maximum
    /// 100000).
    pub max_integer_search_nodes: Option<usize>,
}

/// Input of `recognize_pf_network`.  Give exactly one of `polynomials` or
/// `text`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PfNetworkRequest {
    /// Rows `P_1, P_2, ...` (or `P_0, P_1, ...` with `first_row_index = 0`),
    /// ascending coefficients.
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    /// The same rows as newline-separated coefficient lists or expressions.
    pub text: Option<String>,
    /// Network row of the first polynomial: 1 (default, `P_0 = 1` implicit)
    /// or 0 (the first row must be `1`).
    pub first_row_index: Option<usize>,
    /// Largest step length (default 3, maximum 8).
    pub max_step_length: Option<usize>,
    /// Largest rise (default 3, maximum 8).
    pub max_rise: Option<usize>,
    /// Rational (default) or integer edge weights.
    pub weights: Option<BrentiWeightsInput>,
    /// Budget on the column supports visited (default and maximum 1000000).
    pub max_search_nodes: Option<usize>,
    /// Budget on the integer search for one row (default 10000, maximum
    /// 100000).
    pub max_integer_search_nodes: Option<usize>,
}

/// Input of `check_left_production_matrix`.  Give exactly one of
/// `polynomials` or `text`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LeftProductionRequest {
    /// Rows `P_1, P_2, ...` with `deg P_n = n` (or `P_0, P_1, ...` with
    /// `first_row_index = 0`), ascending coefficients.
    pub polynomials: Option<Vec<BigIntPolynomialInput>>,
    /// The same rows as newline-separated coefficient lists or expressions.
    pub text: Option<String>,
    /// Index of the first row: 1 (default, `P_0 = 1` is prepended) or 0.
    pub first_row_index: Option<usize>,
}

/// Structured representation-search request
/// (`polytool::recurrence::RepresentationSearchRequest`); the library
/// validates the strict schema.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceRepresentationSearchToolRequest {
    pub request: serde_json::Value,
}

/// A `polytool.recurrence-representation-search.v1` report to replay exactly.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct VerifyRecurrenceRepresentationSearchToolRequest {
    pub report: serde_json::Value,
}

/// Input of `enumerate_real_rooted`: all real-rooted integer polynomials of
/// one degree with prescribed leading coefficients.  The fields mirror the
/// options of `polytool enumerate-real-rooted`.
#[derive(Debug, Clone, Deserialize, Serialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EnumerateRealRootedRequest {
    /// Degree `n` of the polynomials, from 2 to 32.
    pub degree: usize,
    /// The first `t` leading coefficients `a_0, a_1, ..., a_(t-1)` in
    /// **descending** order (`3 <= t <= n + 1`, `a_0 != 0`; JSON integers or
    /// exact decimal strings of at most 256 bits).
    #[serde(alias = "leading_desc")]
    pub leading_coefficients_descending: Vec<BigIntCoefficientInput>,
    /// Stop after this many polynomials (default 500, maximum 2000).
    pub max_results: Option<usize>,
    /// Stop after this many search nodes (default 5000, maximum 25000).
    pub max_nodes: Option<u64>,
    /// Stop when a derivative-tower coefficient needs more than this many
    /// bits (default 1024, maximum 2048).
    pub max_coefficient_bits: Option<u64>,
}
