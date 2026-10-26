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
