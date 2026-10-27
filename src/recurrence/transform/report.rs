//! Requests, reports and JSON forms of coordinate transformations: rows from
//! a defining rule, an OEIS family or explicit data; exact reconstruction;
//! condition statuses; transported rules and independent searches; replay.

use super::*;
use crate::vec_poly;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// JSON forms
// ---------------------------------------------------------------------------

fn bivar_json(p: &BivarPoly) -> BivarPolyJson {
    BivarPolyJson {
        coeffs: p
            .normalized_coefficients()
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect(),
    }
}

fn bivar_from_json(p: &BivarPolyJson) -> Result<BivarPoly, String> {
    let coeffs = p
        .coeffs
        .iter()
        .map(|row| row.iter().map(|c| parse_rational_coeff(c)).collect())
        .collect::<Result<Vec<Vec<Q>>, String>>()?;
    Ok(BivarPoly::from_coefficients_normalized(coeffs))
}

fn encode_row(row: &[Q]) -> Vec<String> {
    vec_poly::trim(row)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

fn decode_row(row: &[String]) -> Result<Vec<Q>, String> {
    row.iter()
        .map(|c| parse_rational_coeff(c))
        .collect::<Result<Vec<_>, _>>()
        .map(|r| vec_poly::trim(&r))
}

/// JSON form of a first-order vector system
/// `q(n,x) F_(n+1) = sum_j sum_d matrix[i][j][d](n,x) D^d F_(n,j) + forcing[i](n,x)`
/// in the source index `n`, with its initial state.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TransformVectorRecurrenceJson {
    /// Source index `n` of `initial_state`.
    pub first_index: usize,
    /// `F_(first_index)`: one ascending rational-string row per component.
    pub initial_state: Vec<Vec<String>>,
    /// Common left factor `q(n,x)`.
    pub denominator: BivarPolyJson,
    /// `matrix[i][j][d]` is the coefficient of `D^d` in entry `(i, j)`.
    pub matrix: Vec<Vec<Vec<BivarPolyJson>>>,
    /// Additive vector `G(n,x)`, one entry per component; omitted when homogeneous.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub forcing: Option<Vec<BivarPolyJson>>,
}

impl TransformVectorRecurrenceJson {
    /// Encode a vector system and its initial state.
    pub fn from_recurrence(
        recurrence: &SharedDenominatorVectorRecurrence,
        first_index: usize,
        initial_state: &[Vec<BigRational>],
    ) -> Self {
        Self {
            first_index,
            initial_state: initial_state.iter().map(|row| encode_row(row)).collect(),
            denominator: bivar_json(&recurrence.denominator),
            matrix: recurrence
                .matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|op| op.coefficients.iter().map(bivar_json).collect())
                        .collect()
                })
                .collect(),
            forcing: recurrence
                .forcing
                .as_ref()
                .map(|f| f.iter().map(bivar_json).collect()),
        }
    }

    /// Human-readable form `(q) F_i(n+1) = ...`, one equation per component,
    /// with the given component names.
    pub fn text(&self, names: &[String]) -> Result<String, String> {
        let (system, _, _) = self.to_recurrence()?;
        if names.len() != system.dimension() {
            return Err(format!(
                "{} component names for a {}-component system",
                names.len(),
                system.dimension()
            ));
        }
        Ok(vector_text(&system, names))
    }

    /// Decode the vector system, its first index and initial state.
    pub fn to_recurrence(
        &self,
    ) -> Result<
        (
            SharedDenominatorVectorRecurrence,
            usize,
            Vec<Vec<BigRational>>,
        ),
        String,
    > {
        let matrix = self
            .matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|op| {
                        op.iter()
                            .map(bivar_from_json)
                            .collect::<Result<Vec<_>, _>>()
                            .map(WeylOperator::new)
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let forcing = self
            .forcing
            .as_ref()
            .map(|f| f.iter().map(bivar_from_json).collect::<Result<Vec<_>, _>>())
            .transpose()?;
        let state = self
            .initial_state
            .iter()
            .map(|row| decode_row(row))
            .collect::<Result<Vec<_>, _>>()?;
        Ok((
            SharedDenominatorVectorRecurrence {
                denominator: bivar_from_json(&self.denominator)?,
                matrix,
                forcing,
            },
            self.first_index,
            state,
        ))
    }
}

fn vector_text(recurrence: &SharedDenominatorVectorRecurrence, names: &[String]) -> String {
    let mut lines = Vec::new();
    for (i, row) in recurrence.matrix.iter().enumerate() {
        let mut parts = Vec::new();
        for (j, op) in row.iter().enumerate() {
            for (d, c) in op.coefficients.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                let derivative = match d {
                    0 => String::new(),
                    1 => "'".into(),
                    _ => format!("^({d})"),
                };
                parts.push(format!("({c}) {}{derivative}(n)", names[j]));
            }
        }
        if let Some(forcing) = &recurrence.forcing {
            if !forcing[i].is_zero() {
                parts.push(format!("({})", forcing[i]));
            }
        }
        let rhs = if parts.is_empty() {
            "0".to_string()
        } else {
            parts.join(" + ")
        };
        lines.push(format!(
            "({}) {}(n+1) = {rhs}",
            recurrence.denominator, names[i]
        ));
    }
    lines.join("; ")
}

// ---------------------------------------------------------------------------
// Requests and reports
// ---------------------------------------------------------------------------

/// Bounds of the independent search in the new coordinates.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct TransformSearchBounds {
    /// Largest lag of a scalar rule.
    pub max_rec_len: usize,
    /// Largest `x`-degree of a coefficient.
    pub max_var_deg: usize,
    /// Largest index degree of a coefficient.
    pub max_idx_deg: usize,
    /// Largest derivative order.
    pub max_diff_deg: usize,
    /// Also try a left factor `f(n,x)` (scalar search).
    pub try_denominator: bool,
    /// Candidate budget of the adaptive scalar search.
    pub max_candidates: usize,
    /// Final transitions held out for exact verification (vector search).
    pub vector_held_out_transitions: usize,
    /// Largest `x`-degree of the common left factor (vector search).
    pub vector_denominator_var_deg: usize,
    /// Largest index degree of the common left factor (vector search).
    pub vector_denominator_idx_deg: usize,
}

impl Default for TransformSearchBounds {
    fn default() -> Self {
        Self {
            max_rec_len: 2,
            max_var_deg: 3,
            max_idx_deg: 2,
            max_diff_deg: 2,
            try_denominator: false,
            max_candidates: 5_000,
            vector_held_out_transitions: 2,
            vector_denominator_var_deg: 1,
            vector_denominator_idx_deg: 1,
        }
    }
}

/// Resource limits of [`transform_recurrence`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct TransformLimits {
    /// Largest number of source rows.
    pub max_rows: usize,
    /// Largest degree bound `d_n` or row degree.
    pub max_degree: usize,
    /// Largest derivative order in the source rule.
    pub max_rule_derivative_order: usize,
    /// Largest candidate budget of the independent search.
    pub max_search_candidates: usize,
    /// Largest degree or order bound of the independent search.
    pub max_search_degree: usize,
}

impl Default for TransformLimits {
    fn default() -> Self {
        Self {
            max_rows: 256,
            max_degree: 512,
            max_rule_derivative_order: 8,
            max_search_candidates: 20_000,
            max_search_degree: 8,
        }
    }
}

impl TransformLimits {
    /// Reject limits above a service ceiling.
    pub fn check_ceiling(&self, ceiling: &Self) -> Result<(), String> {
        let pairs = [
            ("max_rows", self.max_rows, ceiling.max_rows),
            ("max_degree", self.max_degree, ceiling.max_degree),
            (
                "max_rule_derivative_order",
                self.max_rule_derivative_order,
                ceiling.max_rule_derivative_order,
            ),
            (
                "max_search_candidates",
                self.max_search_candidates,
                ceiling.max_search_candidates,
            ),
            (
                "max_search_degree",
                self.max_search_degree,
                ceiling.max_search_degree,
            ),
        ];
        for (name, value, limit) in pairs {
            if value > limit {
                return Err(format!(
                    "{name} = {value} exceeds the service ceiling {limit}"
                ));
            }
        }
        Ok(())
    }
}

/// Request for [`transform_recurrence`].
///
/// The source is a defining `polytool.recurrence.v1` rule with seeds
/// (`source_recurrence`), a bundled OEIS identifier (`oeis_id`), or rows only.
/// With a source, `rows` (optional) must agree with the generated rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTransformRequest {
    /// The coordinate change.
    pub transformation: RecurrenceTransformation,
    /// Defining rule with seeds; it holds for targets `first_index + seeds..`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_recurrence: Option<RecurrenceJson>,
    /// Bundled OEIS family (displayed row indices), used as the defining rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oeis_id: Option<String>,
    /// Source rows `P_(row_origin + i)`, ascending rational strings.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rows: Option<Vec<Vec<String>>>,
    /// Original index of `rows[0]`; defaults to the source's first index, else 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_origin: Option<usize>,
    /// Number of rows to generate from the source when `rows` is absent (default 10).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_count: Option<usize>,
    /// Bounds of an independent search in the new coordinates; omitted means no search.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub search: Option<TransformSearchBounds>,
    /// Resource limits.
    #[serde(default)]
    pub limits: TransformLimits,
    /// Optional root-domain analysis (see [`RootDomainRequest`]); omitted means
    /// none, and algebraic discovery never depends on it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_domain: Option<RootDomainRequest>,
}

impl RecurrenceTransformRequest {
    /// A request transforming rows generated from a defining rule.
    pub fn from_source(
        transformation: RecurrenceTransformation,
        source: RecurrenceJson,
        row_count: usize,
    ) -> Self {
        Self {
            transformation,
            source_recurrence: Some(source),
            oeis_id: None,
            rows: None,
            row_origin: None,
            row_count: Some(row_count),
            search: None,
            limits: TransformLimits::default(),
            root_domain: None,
        }
    }
}

/// Status of one precondition or image condition.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TransformConditionStatus {
    /// The condition.
    pub condition: String,
    /// `certified_all_indices` or `checked_supplied_rows`.
    pub status: String,
    /// Why this status was assigned.
    pub detail: String,
}

/// A rule transported from the source rule.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TransportedRuleReport {
    /// Always `transported_identity`: an exact consequence of the source rule,
    /// conditional on that rule.
    pub evidence: String,
    /// First target index (scalar) or first target `n + 1` (vector) at which it holds.
    pub valid_from_index: usize,
    /// Last target index, when the domain is bounded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub valid_to_index: Option<usize>,
    /// Scalar rule with transformed seeds (`polytool.recurrence.v1`, original indices).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scalar: Option<RecurrenceJson>,
    /// Vector system with its initial state.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector: Option<TransformVectorRecurrenceJson>,
    /// The rank-two gamma case: the even part and its odd companion, two
    /// relations of `Q` that together carry the source rule on symmetric
    /// rows (no single scalar rule does).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub system: Option<GammaRelationSystemJson>,
    /// Human-readable form.
    pub text: String,
    /// Whether the rule holds exactly at every checked target.
    pub holds_on_transformed_rows: bool,
    /// Number of targets checked.
    pub checked_targets: usize,
}

/// The rank-two gamma case of a transported rule ([`GammaRelationSystem`]).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GammaRelationSystemJson {
    /// `[even part, odd companion]` as implicit relations of the gamma rows
    /// `Q[n+i]` (base index `n`), each with its first base and the evidence
    /// `transported_identity`: an exact consequence of the source rule at
    /// every target where the rows involved are symmetric at their indices
    /// (the map's image condition, reported in `conditions`).
    pub relations: Vec<DescribedRelation>,
    /// The split index `N` of the cleared identity.
    pub split_index: usize,
    /// The power of `1 + x` multiplied in after dividing by `(1 + x)^(d_m)`.
    pub cleared_power: usize,
}

/// An independent search in the new coordinates.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IndependentSearchReport {
    /// Always `finite_fit`: fitted to finitely many transformed rows.
    pub evidence: String,
    /// `found`, `not_found`, `budget_exhausted` or `failed`.
    pub status: String,
    /// Scalar fit (`polytool.recurrence.v1`, original indices, with search metadata).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scalar: Option<RecurrenceJson>,
    /// Vector fit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector: Option<TransformVectorRecurrenceJson>,
    /// Human-readable form of the fit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub text: Option<String>,
    /// Explanation when nothing was found.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    /// Whether the fit holds exactly on every applicable transformed row.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub holds_on_transformed_rows: Option<bool>,
}

/// Versioned report of [`transform_recurrence`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceTransformReport {
    /// `polytool.recurrence-transform.v1`.
    pub schema: String,
    /// The request, retained for replay.
    pub request: RecurrenceTransformRequest,
    /// Exact map, inverse, conditions and conventions.
    pub transformation: TransformationRecord,
    /// Resolved defining rule with seeds, when a source was given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_rule: Option<RecurrenceJson>,
    /// Original index of the first row.
    pub row_origin: usize,
    /// Source rows `P_(row_origin + i)`.
    pub source_rows: Vec<Vec<String>>,
    /// Transformed states: `[row][component][coefficient]`.
    pub transformed_states: Vec<Vec<Vec<String>>>,
    /// Whether the inverse reproduced every source row exactly.
    pub reconstruction_verified: bool,
    /// Status of every precondition and image condition.
    pub conditions: Vec<TransformConditionStatus>,
    /// The transported rule, when a source rule was given and transport applies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transported: Option<TransportedRuleReport>,
    /// Why transport was not performed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transport_error: Option<String>,
    /// Independent search in the new coordinates, when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub independent_search: Option<IndependentSearchReport>,
    /// Root-domain implications, claims and observations, when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_domain: Option<RootDomainReport>,
}

struct ResolvedSource {
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
    json: RecurrenceJson,
}

fn resolve_source(request: &RecurrenceTransformRequest) -> Result<Option<ResolvedSource>, String> {
    if request.source_recurrence.is_some() && request.oeis_id.is_some() {
        return Err("give at most one of source_recurrence and oeis_id".into());
    }
    let resolved = if let Some(json) = &request.source_recurrence {
        let (rule, first_index, seeds) = json.to_recurrence_parts()?;
        Some((rule, first_index, seeds))
    } else if let Some(id) = &request.oeis_id {
        let definition =
            crate::oeis::by_id(id).ok_or_else(|| format!("unknown OEIS family `{id}`"))?;
        let (rule, seeds, first_row) = definition
            .recurrence_export_parts()
            .map_err(|e| e.to_string())?;
        let first_index = usize::try_from(first_row)
            .map_err(|_| format!("OEIS family `{id}` starts at a negative row index"))?;
        Some((rule, first_index, seeds))
    } else {
        None
    };
    let Some((rule, first_index, seeds)) = resolved else {
        return Ok(None);
    };
    if rule.terms.iter().any(|t| t.offset == 0) {
        return Err("source rule terms need positive offsets".into());
    }
    if seeds.is_empty() || rule.max_offset() > seeds.len() {
        return Err("source rule needs at least max-offset seeds".into());
    }
    if rule
        .terms
        .iter()
        .any(|t| t.deriv_order > request.limits.max_rule_derivative_order)
    {
        return Err("source rule derivative order exceeds max_rule_derivative_order".into());
    }
    let json = RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None);
    Ok(Some(ResolvedSource {
        rule,
        first_index,
        seeds,
        json,
    }))
}

fn resolve_rows(
    request: &RecurrenceTransformRequest,
    source: Option<&ResolvedSource>,
) -> Result<(usize, Vec<Vec<Q>>), String> {
    let limit = request.limits.max_rows;
    let explicit = request
        .rows
        .as_ref()
        .map(|rows| {
            rows.iter()
                .map(|r| decode_row(r))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let origin = match (request.row_origin, source) {
        (Some(origin), Some(s)) if origin != s.first_index => {
            return Err(format!(
                "row_origin {origin} differs from the source first index {}",
                s.first_index
            ))
        }
        (Some(origin), _) => origin,
        (None, Some(s)) => s.first_index,
        (None, None) => 0,
    };
    let rows = match (explicit, source) {
        (None, None) => return Err("give rows, source_recurrence or oeis_id".into()),
        (Some(rows), None) => rows,
        (explicit, Some(s)) => {
            let count = explicit
                .as_ref()
                .map_or(request.row_count.unwrap_or(10), Vec::len);
            if count > limit {
                return Err(format!("{count} rows exceed max_rows = {limit}"));
            }
            let count = count.max(s.seeds.len());
            let generated = s
                .rule
                .generate_rows_rational(&s.seeds, s.first_index, count)
                .map_err(|e| format!("source rule generation failed: {e}"))?
                .into_iter()
                .map(|r| vec_poly::trim(&r))
                .collect::<Vec<_>>();
            if let Some(rows) = explicit {
                if rows != generated[..rows.len()] {
                    return Err("rows disagree with the rows generated by the source rule".into());
                }
            }
            generated
        }
    };
    if rows.is_empty() {
        return Err("at least one row is required".into());
    }
    if rows.len() > limit {
        return Err(format!("{} rows exceed max_rows = {limit}", rows.len()));
    }
    if let Some(degree) = rows.iter().filter_map(|r| vec_poly::degree(r)).max() {
        if degree > request.limits.max_degree {
            return Err(format!("row degree {degree} exceeds max_degree"));
        }
    }
    Ok((origin, rows))
}

/// The scalar reversal transport of `rule` (for symmetry certification).
fn reversal_rule(rule: &Recurrence, schedule: &DegreeSchedule) -> Option<Recurrence> {
    match (RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    })
    .transport(rule)
    {
        Ok(TransportedRecurrence::Scalar(r)) => Some(r),
        _ => None,
    }
}

fn condition_statuses(
    transformation: &RecurrenceTransformation,
    source: Option<&ResolvedSource>,
    transported: Option<&TransportedRecurrence>,
) -> Vec<TransformConditionStatus> {
    let checked = |condition: &str, detail: &str| TransformConditionStatus {
        condition: condition.to_string(),
        status: CHECKED_SUPPLIED_ROWS.to_string(),
        detail: detail.to_string(),
    };
    let certified = |condition: &str, detail: String| TransformConditionStatus {
        condition: condition.to_string(),
        status: CERTIFIED_ALL_INDICES.to_string(),
        detail,
    };
    let finite = "exact on every supplied row; finite evidence only";
    let polynomial_rule = |rule: Option<&Recurrence>| rule.is_some_and(|r| r.denominator.is_none());
    let scalar = match transported {
        Some(TransportedRecurrence::Scalar(r)) => Some(r),
        _ => None,
    };
    let induction = |what: &str, s: &ResolvedSource| {
        format!(
            "{what} is denominator-free with polynomial coefficients and the seeds at indices \
             {}..{} satisfy the condition; induction over the domain",
            s.first_index,
            s.first_index + s.seeds.len()
        )
    };
    let degree_condition = "deg P_n <= d_n";
    match transformation {
        RecurrenceTransformation::MultiplyPower { .. }
        | RecurrenceTransformation::Affine { .. }
        | RecurrenceTransformation::Derivative { .. } => Vec::new(),
        RecurrenceTransformation::Diagonal { .. }
        | RecurrenceTransformation::DerivativeAntecedent { .. } => {
            implicit_map_condition_statuses(transformation, source.map_or(0, |s| s.first_index))
        }
        RecurrenceTransformation::DividePower { power } => {
            let condition = format!("x^{power} divides P_n");
            vec![match source {
                Some(s) if polynomial_rule(scalar) => certified(
                    &condition,
                    induction("the transported rule for Q_n = P_n/x^r", s),
                ),
                _ => checked(&condition, finite),
            }]
        }
        RecurrenceTransformation::Reversal { .. }
        | RecurrenceTransformation::FPolynomial { .. }
        | RecurrenceTransformation::Mobius { .. } => {
            let pole_free = transformation
                .mobius_parameters()
                .ok()
                .flatten()
                .is_some_and(|m| !m.c.is_zero());
            vec![match source {
                Some(s) if pole_free && polynomial_rule(scalar) => certified(
                    degree_condition,
                    induction(
                        "the transported rule (Q_n polynomial iff deg P_n <= d_n)",
                        s,
                    ),
                ),
                _ => checked(degree_condition, finite),
            }]
        }
        RecurrenceTransformation::SymmetricDecomposition { schedule } => {
            let reversed = source.and_then(|s| reversal_rule(&s.rule, schedule));
            vec![match source {
                Some(s) if polynomial_rule(reversed.as_ref()) => certified(
                    degree_condition,
                    induction("the reversal-transported rule", s),
                ),
                _ => checked(degree_condition, finite),
            }]
        }
        RecurrenceTransformation::Gamma { schedule } => {
            let condition = "x^(d_n) P_n(1/x) = P_n";
            let certificate = source.and_then(|s| {
                let reversed = reversal_rule(&s.rule, schedule)?;
                let valid_from = s.first_index + s.seeds.len();
                let leader = canonical_recurrence(&s.rule).denominator;
                let leader_nonzero = leader.as_ref().is_none_or(|f| {
                    crate::recurrence::root_domain::bivariate_leader_nonvanishing(f, valid_from)
                });
                (leader_nonzero && super::gamma::rules_proportional(&reversed, &s.rule)).then(
                    || {
                        format!(
                            "the reversal-transported rule defines the same row map as the source \
                         rule, whose left factor is a nonzero polynomial at every target \
                         m >= {valid_from}, and the seeds at indices {}..{} are symmetric; the \
                         reversed rows satisfy the same rule with the same seeds, so they \
                         coincide with the rows by induction over the domain",
                            s.first_index, valid_from
                        )
                    },
                )
            });
            vec![match certificate {
                Some(detail) => certified(condition, detail),
                None => checked(condition, finite),
            }]
        }
    }
}

fn rule_text(rule: &Recurrence) -> String {
    rule.to_string()
}

fn search_bounds_check(
    bounds: &TransformSearchBounds,
    limits: &TransformLimits,
) -> Result<(), String> {
    let degree = [
        bounds.max_rec_len,
        bounds.max_var_deg,
        bounds.max_idx_deg,
        bounds.max_diff_deg,
        bounds.vector_denominator_var_deg,
        bounds.vector_denominator_idx_deg,
    ]
    .into_iter()
    .max()
    .unwrap_or(0);
    if degree > limits.max_search_degree {
        return Err(format!(
            "search bound {degree} exceeds max_search_degree = {}",
            limits.max_search_degree
        ));
    }
    if bounds.max_candidates > limits.max_search_candidates {
        return Err(format!(
            "max_candidates {} exceeds max_search_candidates = {}",
            bounds.max_candidates, limits.max_search_candidates
        ));
    }
    Ok(())
}

fn independent_search(
    transformation: &RecurrenceTransformation,
    bounds: &TransformSearchBounds,
    origin: usize,
    states: &[Vec<Vec<Q>>],
) -> IndependentSearchReport {
    let mut report = IndependentSearchReport {
        evidence: FINITE_FIT_EVIDENCE.to_string(),
        status: "not_found".to_string(),
        scalar: None,
        vector: None,
        text: None,
        reason: None,
        holds_on_transformed_rows: None,
    };
    if transformation.component_count() == 1 {
        let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
        let options = AdaptiveSearchOptions {
            max_rec_len: bounds.max_rec_len,
            max_var_deg: bounds.max_var_deg,
            max_idx_deg: bounds.max_idx_deg,
            max_diff_deg: bounds.max_diff_deg,
            try_denominator: bounds.try_denominator,
            ..AdaptiveSearchOptions::default()
        };
        let outcome = find_recurrence_adaptive_rational_with_budget(
            &rows,
            &options,
            AdaptiveSearchBudget::limited(bounds.max_candidates),
        );
        match outcome {
            AdaptiveSearchOutcome::Found(found) => {
                // Search indices are 1-based positions; shift to original indices.
                let shift = 1 - origin as i64;
                let rule = found.recurrence.shift_index(shift);
                let initial = rule.generation_initial_count(origin, rows.len());
                let holds = (origin + rule.max_offset()..origin + rows.len()).all(|target| {
                    transformed_rule_holds_at(&rule, &rows, origin, target) == Some(true)
                });
                let search = RecurrenceJsonSearch {
                    recurrence_text: rule_text(&rule),
                    source_rows: rows.len(),
                    skip_prefix: 0,
                    unknowns: found.num_unknowns,
                    weighted_unknowns: found.weighted_unknowns,
                    equations: found.num_equations,
                    fit_polynomials: found.fit_polynomials,
                    verification_polynomials: found.verification_polynomials,
                    candidates_tried: found.candidates_tried,
                    options: RecurrenceOptionsJson::from(&found.opts),
                };
                report.status = "found".into();
                report.text = Some(rule_text(&rule));
                report.scalar = Some(RecurrenceJson::from_recurrence_rational(
                    &rule,
                    origin,
                    &rows[..initial],
                    Some(search),
                ));
                report.holds_on_transformed_rows = Some(holds);
            }
            AdaptiveSearchOutcome::NoRecurrence(_) => {
                report.reason = Some("no rule within the search bounds".into());
            }
            AdaptiveSearchOutcome::BudgetExhausted(_) => {
                report.status = "budget_exhausted".into();
                report.reason = Some("candidate budget exhausted".into());
            }
        }
    } else {
        let options = VectorRecurrenceDenominatorOptions {
            var_deg: bounds.max_var_deg,
            idx_deg: bounds.max_idx_deg,
            diff_deg: bounds.max_diff_deg,
            homogeneous: true,
            forcing_var_deg: 0,
            forcing_idx_deg: 0,
            held_out_transitions: bounds.vector_held_out_transitions,
            require_unique: false,
            denominator_var_deg: bounds.vector_denominator_var_deg,
            denominator_idx_deg: bounds.vector_denominator_idx_deg,
        };
        match find_vector_recurrence_with_denominator_rational(states, origin, &options) {
            Ok(fit) => {
                let holds = states.windows(2).enumerate().all(|(i, pair)| {
                    fit.recurrence
                        .holds_transition_rational(&pair[0], &pair[1], origin + i)
                });
                report.status = "found".into();
                report.text = Some(vector_text(
                    &fit.recurrence,
                    &component_names(transformation),
                ));
                report.vector = Some(TransformVectorRecurrenceJson::from_recurrence(
                    &fit.recurrence,
                    origin,
                    &states[0],
                ));
                report.holds_on_transformed_rows = Some(holds);
            }
            Err(error) => {
                report.status = "failed".into();
                report.reason = Some(error.to_string());
            }
        }
    }
    report
}

fn component_names(transformation: &RecurrenceTransformation) -> Vec<String> {
    transformation
        .record()
        .components
        .into_iter()
        .map(|c| c.name)
        .collect()
}

/// Apply a transformation to a sequence, transport its defining rule and
/// optionally search independently in the new coordinates.
///
/// Every row is mapped and reconstructed exactly; a row violating a
/// precondition (degree bound, divisibility, symmetry) or outside the index
/// domain rejects the request.
pub fn transform_recurrence(
    request: &RecurrenceTransformRequest,
) -> Result<RecurrenceTransformReport, String> {
    let transformation = &request.transformation;
    transformation.validate().map_err(|e| e.to_string())?;
    if let Some(bounds) = &request.search {
        search_bounds_check(bounds, &request.limits)?;
    }
    let source = resolve_source(request)?;
    let (origin, rows) = resolve_rows(request, source.as_ref())?;

    let mut states = Vec::with_capacity(rows.len());
    let mut reconstruction_verified = true;
    for (i, row) in rows.iter().enumerate() {
        let index = origin + i;
        if let Some(schedule) = transformation.schedule() {
            let bound = schedule.degree_at(index).map_err(|e| e.to_string())?;
            if bound > request.limits.max_degree {
                return Err(format!(
                    "degree bound {bound} at index {index} exceeds max_degree"
                ));
            }
        }
        let state = transformation
            .map_row(index, row)
            .map_err(|e| e.to_string())?;
        let back = transformation
            .inverse_row(index, &state)
            .map_err(|e| e.to_string())?;
        reconstruction_verified &= vec_poly::trim(&back) == vec_poly::trim(row);
        states.push(state);
    }

    let mut transported_rule = None;
    let mut transport_error = None;
    let mut transported_report = None;
    if let Some(s) = &source {
        match transformation.transport(&s.rule) {
            Ok(transported) => {
                transported_report = Some(transported_report_for(
                    transformation,
                    &transported,
                    s,
                    origin,
                    &states,
                ));
                transported_rule = Some(transported);
            }
            Err(error) => match gamma_system_report(transformation, s, origin, &states) {
                Some(report) => transported_report = Some(report),
                None => transport_error = Some(error.to_string()),
            },
        }
    }
    let conditions = condition_statuses(transformation, source.as_ref(), transported_rule.as_ref());
    let independent = request
        .search
        .as_ref()
        .map(|bounds| independent_search(transformation, bounds, origin, &states));
    let root_domain =
        request
            .root_domain
            .as_ref()
            .map(|root_request| {
                crate::recurrence::root_domain::analyze_root_domain(
                    &crate::recurrence::root_domain::RootDomainInputs {
                        request: root_request,
                        transformation,
                        origin,
                        rows: &rows,
                        states: &states,
                        conditions: &conditions,
                        transported: source.as_ref().zip(transported_rule.as_ref()).map(
                            |(s, rule)| crate::recurrence::root_domain::TransportedInput {
                                rule,
                                first_index: s.first_index,
                                valid_from: s.first_index + s.seeds.len(),
                            },
                        ),
                    },
                )
            })
            .transpose()?;

    Ok(RecurrenceTransformReport {
        schema: RECURRENCE_TRANSFORM_SCHEMA.to_string(),
        request: request.clone(),
        transformation: transformation.record(),
        source_rule: source.as_ref().map(|s| s.json.clone()),
        row_origin: origin,
        source_rows: rows.iter().map(|r| encode_row(r)).collect(),
        transformed_states: states
            .iter()
            .map(|state| state.iter().map(|c| encode_row(c)).collect())
            .collect(),
        reconstruction_verified,
        conditions,
        transported: transported_report,
        transport_error,
        independent_search: independent,
        root_domain,
    })
}

/// The rank-two gamma case as a two-relation system, when it can be stated
/// exactly ([`GammaTransport::relation_system`]).
fn gamma_system_report(
    transformation: &RecurrenceTransformation,
    source: &ResolvedSource,
    origin: usize,
    states: &[Vec<Vec<Q>>],
) -> Option<TransportedRuleReport> {
    let RecurrenceTransformation::Gamma { schedule } = transformation else {
        return None;
    };
    let gamma = transport_gamma(&source.rule, schedule).ok()?;
    let system = gamma.relation_system().ok()?;
    let valid_from = source.first_index + source.seeds.len();
    let valid_to = schedule.last_index;
    let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
    let end = origin + rows.len();
    let mut holds = true;
    let mut relations = Vec::new();
    for (relation, shift, note) in [
        (
            &system.even,
            system.even_shift,
            "even part of the gamma-transported identity",
        ),
        (
            &system.odd,
            system.odd_shift,
            "odd companion of the gamma-transported identity",
        ),
    ] {
        let first_base = valid_from.checked_sub(shift)?;
        for base in first_base..end.saturating_sub(relation.max_shift()) {
            if valid_to.is_some_and(|last| base + shift > last) {
                break;
            }
            holds &= relation.holds_at(&rows, origin, base) != Some(false);
        }
        relations.push(
            DescribedRelation::new(relation, first_base, TRANSPORTED_IDENTITY_EVIDENCE, "Q")
                .with_note(note),
        );
    }
    let checked_targets = (valid_from..end)
        .filter(|t| valid_to.is_none_or(|last| *t <= last))
        .count();
    Some(TransportedRuleReport {
        evidence: TRANSPORTED_IDENTITY_EVIDENCE.to_string(),
        valid_from_index: valid_from,
        valid_to_index: valid_to,
        scalar: None,
        vector: None,
        text: format!("even: {}; odd: {}", relations[0].text, relations[1].text),
        system: Some(GammaRelationSystemJson {
            relations,
            split_index: gamma.split_index,
            cleared_power: gamma.cleared_power,
        }),
        holds_on_transformed_rows: holds,
        checked_targets,
    })
}

fn transported_report_for(
    transformation: &RecurrenceTransformation,
    transported: &TransportedRecurrence,
    source: &ResolvedSource,
    origin: usize,
    states: &[Vec<Vec<Q>>],
) -> TransportedRuleReport {
    let valid_from = source.first_index + source.seeds.len();
    let valid_to = transformation.schedule().and_then(|s| s.last_index);
    let end = origin + states.len();
    let targets = (valid_from..end)
        .filter(|t| valid_to.is_none_or(|last| *t <= last))
        .collect::<Vec<_>>();
    match transported {
        TransportedRecurrence::Scalar(rule) => {
            let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
            let holds = targets
                .iter()
                .all(|&t| transformed_rule_holds_at(rule, &rows, origin, t) == Some(true));
            let seeds = &rows[source.first_index - origin..valid_from - origin];
            TransportedRuleReport {
                evidence: TRANSPORTED_IDENTITY_EVIDENCE.to_string(),
                valid_from_index: valid_from,
                valid_to_index: valid_to,
                scalar: Some(RecurrenceJson::from_recurrence_rational(
                    rule,
                    source.first_index,
                    seeds,
                    None,
                )),
                vector: None,
                system: None,
                text: rule_text(rule),
                holds_on_transformed_rows: holds,
                checked_targets: targets.len(),
            }
        }
        TransportedRecurrence::Vector(system) => {
            let holds = targets.iter().all(|&t| {
                system.holds_transition_rational(
                    &states[t - 1 - origin],
                    &states[t - origin],
                    t - 1,
                )
            });
            TransportedRuleReport {
                evidence: TRANSPORTED_IDENTITY_EVIDENCE.to_string(),
                valid_from_index: valid_from,
                valid_to_index: valid_to,
                scalar: None,
                vector: Some(TransformVectorRecurrenceJson::from_recurrence(
                    system,
                    valid_from - 1,
                    &states[valid_from - 1 - origin],
                )),
                system: None,
                text: vector_text(system, &component_names(transformation)),
                holds_on_transformed_rows: holds,
                checked_targets: targets.len(),
            }
        }
    }
}

/// Replay a report exactly from its embedded request.
///
/// Recomputes the transformation, reconstruction, transport and search and
/// requires the result to equal the report.  The transported rule is replayed
/// as an identity of the supplied source rule; it is not evidence that the
/// source rule itself is true.
pub fn verify_recurrence_transform(report: &RecurrenceTransformReport) -> Result<(), String> {
    if report.schema != RECURRENCE_TRANSFORM_SCHEMA {
        return Err(format!("unsupported schema `{}`", report.schema));
    }
    let replay = transform_recurrence(&report.request)?;
    if &replay != report {
        return Err("report does not match its exact replay".into());
    }
    if !report.reconstruction_verified {
        return Err("exact reconstruction failed".into());
    }
    if report
        .transported
        .as_ref()
        .is_some_and(|t| !t.holds_on_transformed_rows)
    {
        return Err("transported rule fails on the transformed rows".into());
    }
    Ok(())
}
