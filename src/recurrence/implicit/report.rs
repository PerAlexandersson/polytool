//! Requests, reports and replay of implicit transports: a defining rule is
//! mapped by a diagonal map or a derivative antecedent, transported to an
//! implicit relation with exceptional-index conditions, assembled into a
//! generating description with a uniqueness argument, and optionally
//! complemented by a fitted successor rule that a reconstruction certificate
//! proves for every index.
//!
//! Evidence is kept apart throughout: `transported_identity` (exact
//! consequence of the defining rule), `checked_supplied_rows` (supplied,
//! verified on rows), `finite_fit` (fitted) and `certified_by_reconstruction`
//! (a fit proved by `L o T = sum H_nu S^nu M`, source uniqueness, initial
//! agreement and injectivity of the map).

use super::diagonal::{observed_multiplier_zeros, weight_regularity};
use super::*;
use crate::vec_poly;
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Schema tag of [`ImplicitTransformReport`].
pub const RECURRENCE_IMPLICIT_SCHEMA: &str = "polytool.recurrence-implicit.v1";

/// A supplied relation with the first base index at which it holds.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SuppliedRelation {
    /// The relation.
    pub relation: ImplicitRelationJson,
    /// It holds for base indices `n >= first_base`.
    pub first_base: usize,
}

/// A supplied coefficient condition with the first base index at which it holds.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct SuppliedCondition {
    /// The condition.
    pub condition: CoefficientConditionJson,
    /// It holds for base indices `n >= first_base`.
    pub first_base: usize,
}

/// Request for [`implicit_transform_recurrence`].
///
/// The source is a defining `polytool.recurrence.v1` rule with seeds
/// (`source_recurrence`) or a bundled OEIS family (`oeis_id`); the
/// transformation must be `diagonal` or `derivative_antecedent`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitTransformRequest {
    /// The coordinate change.
    pub transformation: RecurrenceTransformation,
    /// Defining rule with seeds; it holds for targets `first_index + seeds..`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_recurrence: Option<RecurrenceJson>,
    /// Bundled OEIS family (displayed row indices), used as the defining rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oeis_id: Option<String>,
    /// Number of rows generated from the source (default 10).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_count: Option<usize>,
    /// Further known relations of the source rows `P`, transported as well;
    /// they are checked on the generated rows only.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub source_relations: Vec<SuppliedRelation>,
    /// Relations of the new rows `Q` to store; checked on the rows only.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub target_relations: Vec<SuppliedRelation>,
    /// Coefficient boundary conditions of `Q`; checked on the rows only.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub boundary_conditions: Vec<SuppliedCondition>,
    /// Bounds of a successor-rule search in the new coordinates; omitted means
    /// no search.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub search: Option<TransformSearchBounds>,
    /// Resource limits.
    #[serde(default)]
    pub limits: TransformLimits,
}

/// One transported relation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitTransportReport {
    /// `defining_rule` or `source_relations[i]`.
    pub source: String,
    /// The transported relation in the new coordinates.
    pub relation: DescribedRelation,
    /// Clearing data of a diagonal map.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clearing: Option<ClearingCertificate>,
    /// Conditions at the certified exceptional indices.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exceptional_index_conditions: Vec<DescribedCondition>,
    /// `L o T = sum H_nu S^nu M`, when `T` is a differential operator.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reconstruction: Option<ReconstructionCertificate>,
    /// Whether the relation and its conditions hold at every checked base.
    pub holds_on_transformed_rows: bool,
    /// Number of base indices checked.
    pub checked_bases: usize,
}

/// A successor-rule search in the new coordinates and its upgrade attempt.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitSearchReport {
    /// `finite_fit`, or `certified_by_reconstruction` after a successful upgrade.
    pub evidence: String,
    /// `found`, `not_found` or `budget_exhausted`.
    pub status: String,
    /// The fitted rule (`polytool.recurrence.v1`, original indices).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scalar: Option<RecurrenceJson>,
    /// The fitted rule as an implicit relation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relation: Option<DescribedRelation>,
    /// `L o T = sum H_nu S^nu M` for the fitted `M`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reconstruction: Option<ReconstructionCertificate>,
    /// The checks of the upgrade, in order.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub upgrade_steps: Vec<String>,
    /// The fitted rule with its seeds as a generating description, when upgraded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<GeneratingDescription>,
    /// Assessment of that description.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assessment: Option<DescriptionAssessment>,
    /// Why nothing was found or why the fit was not upgraded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

/// Versioned report of [`implicit_transform_recurrence`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitTransformReport {
    /// `polytool.recurrence-implicit.v1`.
    pub schema: String,
    /// The request, retained for replay.
    pub request: ImplicitTransformRequest,
    /// Exact map, inverse, conditions and conventions.
    pub transformation: TransformationRecord,
    /// Hypotheses of the map (weights defined and nonzero) and their status.
    pub map_conditions: Vec<TransformConditionStatus>,
    /// The resolved defining rule with seeds.
    pub source_rule: RecurrenceJson,
    /// Original index of the first row.
    pub row_origin: usize,
    /// Source rows `P_(row_origin + i)`.
    pub source_rows: Vec<Vec<String>>,
    /// Transformed rows `Q_(row_origin + i)`.
    pub transformed_rows: Vec<Vec<String>>,
    /// Whether the inverse reproduced every source row exactly.
    pub reconstruction_verified: bool,
    /// The defining rule and every supplied source relation, transported.
    pub transports: Vec<ImplicitTransportReport>,
    /// Seeds, transported defining relation, exceptional-index conditions and
    /// supplied boundary conditions.
    pub description: GeneratingDescription,
    /// Annihilation evidence and uniqueness argument of `description`.
    pub assessment: DescriptionAssessment,
    /// Further relations of `Q` (transported supplied relations, supplied
    /// target relations) with their evidence; not part of `description`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub additional_relations: Vec<DescribedRelation>,
    /// The successor-rule search, when requested.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub independent_search: Option<ImplicitSearchReport>,
    /// Root-domain analysis, when requested
    /// ([`implicit_transform_recurrence_with_root_domain`]); it carries its
    /// own request for replay.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_domain: Option<ImplicitRootDomainReport>,
}

struct Source {
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
    json: RecurrenceJson,
}

fn resolve(request: &ImplicitTransformRequest) -> Result<Source, String> {
    let (rule, first_index, seeds) = match (&request.source_recurrence, &request.oeis_id) {
        (Some(_), Some(_)) => return Err("give only one of source_recurrence and oeis_id".into()),
        (None, None) => {
            return Err(
                "implicit transport needs a defining rule: source_recurrence or oeis_id".into(),
            )
        }
        (Some(json), None) => json.to_recurrence_parts()?,
        (None, Some(id)) => {
            let definition =
                crate::oeis::by_id(id).ok_or_else(|| format!("unknown OEIS family `{id}`"))?;
            let (rule, seeds, first_row) = definition
                .recurrence_export_parts()
                .map_err(|e| e.to_string())?;
            let first = usize::try_from(first_row)
                .map_err(|_| format!("OEIS family `{id}` starts at a negative row index"))?;
            (rule, first, seeds)
        }
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
    Ok(Source {
        rule,
        first_index,
        seeds,
        json,
    })
}

fn check_search_bounds(
    bounds: &TransformSearchBounds,
    limits: &TransformLimits,
) -> Result<(), String> {
    let degree = [
        bounds.max_rec_len,
        bounds.max_var_deg,
        bounds.max_idx_deg,
        bounds.max_diff_deg,
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

fn encode_rows(rows: &[Vec<Q>]) -> Vec<Vec<String>> {
    rows.iter().map(|r| encode_poly(r)).collect()
}

/// Bases `n >= first_base` whose rows `n..=n+r` are available.
fn checked_bases(first_base: usize, r: usize, origin: usize, len: usize) -> std::ops::Range<usize> {
    let end = (origin + len).saturating_sub(r);
    first_base.max(origin)..end.max(first_base.max(origin))
}

/// Apply a diagonal map or derivative antecedent to a sequence with a defining
/// rule, transport the rule (and supplied source relations) into implicit
/// relations with exceptional-index conditions, assess the resulting
/// generating description, and optionally fit a successor rule in the new
/// coordinates and try to prove it by a reconstruction certificate.
///
/// Supplied relations and conditions that fail on the rows reject the request.
pub fn implicit_transform_recurrence(
    request: &ImplicitTransformRequest,
) -> Result<ImplicitTransformReport, String> {
    let transformation = &request.transformation;
    transformation.validate().map_err(|e| e.to_string())?;
    if !matches!(
        transformation,
        RecurrenceTransformation::Diagonal { .. }
            | RecurrenceTransformation::DerivativeAntecedent { .. }
    ) {
        return Err(format!(
            "implicit transport needs a diagonal or derivative_antecedent map, not {}; use \
             transform_recurrence",
            transformation.kind_name()
        ));
    }
    if let Some(bounds) = &request.search {
        check_search_bounds(bounds, &request.limits)?;
    }
    let source = resolve(request)?;
    let origin = source.first_index;
    let count = request.row_count.unwrap_or(10);
    if count > request.limits.max_rows {
        return Err(format!(
            "{count} rows exceed max_rows = {}",
            request.limits.max_rows
        ));
    }
    let rows = source
        .rule
        .generate_rows_rational(&source.seeds, origin, count.max(source.seeds.len()))
        .map_err(|e| format!("source rule generation failed: {e}"))?
        .into_iter()
        .map(|r| vec_poly::trim(&r))
        .collect::<Vec<_>>();
    let max_degree = rows
        .iter()
        .filter_map(|r| vec_poly::degree(r))
        .max()
        .unwrap_or(0);
    if max_degree > request.limits.max_degree {
        return Err(format!("row degree {max_degree} exceeds max_degree"));
    }

    let map_conditions = match transformation {
        RecurrenceTransformation::Diagonal { weight } => {
            weight_regularity(weight, origin).map_err(|e| e.to_string())?
        }
        _ => vec![TransformConditionStatus {
            condition: "T_g bijective".into(),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: "the monomial matrix of T_g is triangular with diagonal lc(g) binom(k+m,m) \
                     != 0"
                .into(),
        }],
    };
    let map_certified = map_conditions
        .iter()
        .all(|c| c.status == CERTIFIED_ALL_INDICES);

    let mut q_rows = Vec::with_capacity(rows.len());
    let mut reconstruction_verified = true;
    for (i, row) in rows.iter().enumerate() {
        let index = origin + i;
        let state = transformation
            .map_row(index, row)
            .map_err(|e| e.to_string())?;
        let back = transformation
            .inverse_row(index, &state)
            .map_err(|e| e.to_string())?;
        reconstruction_verified &= vec_poly::trim(&back) == vec_poly::trim(row);
        q_rows.push(vec_poly::trim(&state[0]));
    }

    let (defining, r) = ImplicitRelation::from_recurrence(&source.rule)?;
    let valid_from = origin + source.seeds.len();
    let first_base = valid_from - r;
    let mut transports = Vec::new();
    let mut additional = Vec::new();
    // Without a certified weight hypothesis, Q_n = T_n P_n is only known to
    // exist on the supplied rows, so the transported identities are finite
    // evidence there.
    let defining_evidence = if map_certified {
        TRANSPORTED_IDENTITY_EVIDENCE
    } else {
        CHECKED_SUPPLIED_ROWS
    };
    let defining_report = transport_report(
        transformation,
        "defining_rule",
        &defining,
        first_base,
        defining_evidence,
        origin,
        &q_rows,
        max_degree,
    )?;

    for (i, supplied) in request.source_relations.iter().enumerate() {
        let relation = ImplicitRelation::from_json(&supplied.relation)?;
        check_on_rows(
            &relation,
            supplied.first_base,
            origin,
            &rows,
            &format!("source_relations[{i}]"),
        )?;
        let report = transport_report(
            transformation,
            &format!("source_relations[{i}]"),
            &relation,
            supplied.first_base,
            TRANSPORTED_SUPPLIED_EVIDENCE,
            origin,
            &q_rows,
            max_degree,
        )?;
        additional.push(report.relation.clone().with_note(format!(
            "transported from source_relations[{i}], which is checked on the supplied rows only"
        )));
        transports.push(report);
    }
    for (i, supplied) in request.target_relations.iter().enumerate() {
        let relation = ImplicitRelation::from_json(&supplied.relation)?;
        check_on_rows(
            &relation,
            supplied.first_base,
            origin,
            &q_rows,
            &format!("target_relations[{i}]"),
        )?;
        additional.push(
            DescribedRelation::new(&relation, supplied.first_base, CHECKED_SUPPLIED_ROWS, "Q")
                .with_note(format!("supplied as target_relations[{i}]")),
        );
    }
    let mut boundary = Vec::new();
    for (i, supplied) in request.boundary_conditions.iter().enumerate() {
        let condition = CoefficientCondition::from_json(&supplied.condition)?;
        let r = condition.max_shift();
        for base in checked_bases(supplied.first_base, r, origin, q_rows.len()) {
            if condition.holds_at(&q_rows, origin, base) == Some(false) {
                return Err(format!(
                    "boundary_conditions[{i}] fails on the transformed rows at base n = {base}"
                ));
            }
        }
        boundary.push(
            DescribedCondition::new(
                &condition,
                supplied.first_base,
                CHECKED_SUPPLIED_ROWS,
                None,
                "Q",
            )
            .with_note(format!("supplied as boundary_conditions[{i}]")),
        );
    }

    let description = GeneratingDescription {
        first_index: origin,
        seeds: encode_rows(&q_rows[..source.seeds.len()]),
        relations: vec![defining_report.relation.clone()],
        exceptional_index_conditions: defining_report.exceptional_index_conditions.clone(),
        boundary_conditions: boundary,
    };
    let assessment = assess_generating_description(&description)?;
    transports.insert(0, defining_report);

    let independent_search = request.search.as_ref().map(|bounds| {
        search_report(
            transformation,
            bounds,
            &SearchContext {
                source: &source,
                defining: &defining,
                defining_first_base: first_base,
                map_certified,
                origin,
                q_rows: &q_rows,
            },
        )
    });

    Ok(ImplicitTransformReport {
        schema: RECURRENCE_IMPLICIT_SCHEMA.to_string(),
        request: request.clone(),
        transformation: transformation.record(),
        map_conditions,
        source_rule: source.json.clone(),
        row_origin: origin,
        source_rows: encode_rows(&rows),
        transformed_rows: encode_rows(&q_rows),
        reconstruction_verified,
        transports,
        description,
        assessment,
        additional_relations: additional,
        independent_search,
        root_domain: None,
    })
}

/// [`implicit_transform_recurrence`] followed, when `root_domain` is given, by
/// [`analyze_implicit_root_domain`]; the analysis (with its request) is stored
/// in the report's `root_domain` field and replayed by
/// [`verify_recurrence_implicit`].
pub fn implicit_transform_recurrence_with_root_domain(
    request: &ImplicitTransformRequest,
    root_domain: Option<&RootDomainRequest>,
) -> Result<ImplicitTransformReport, String> {
    let mut report = implicit_transform_recurrence(request)?;
    if let Some(root_request) = root_domain {
        report.root_domain = Some(analyze_implicit_root_domain(&report, root_request)?);
    }
    Ok(report)
}

/// Parse a JSON implicit-transport request: the fields of
/// [`ImplicitTransformRequest`] plus an optional `root_domain`
/// ([`RootDomainRequest`], `{}` for the defaults).
pub fn parse_implicit_transform_input(
    mut value: serde_json::Value,
) -> Result<(ImplicitTransformRequest, Option<RootDomainRequest>), String> {
    let root_domain = value
        .as_object_mut()
        .and_then(|fields| fields.remove("root_domain"));
    let request = serde_json::from_value(value).map_err(|e| e.to_string())?;
    let root_domain = root_domain
        .map(serde_json::from_value)
        .transpose()
        .map_err(|e| format!("root_domain: {e}"))?;
    Ok((request, root_domain))
}

fn check_on_rows(
    relation: &ImplicitRelation,
    first_base: usize,
    origin: usize,
    rows: &[Vec<Q>],
    label: &str,
) -> Result<usize, String> {
    if relation.is_zero() {
        return Err(format!("{label} is the zero relation"));
    }
    let bases = checked_bases(first_base, relation.max_shift(), origin, rows.len());
    let count = bases.len();
    for base in bases {
        if relation.holds_at(rows, origin, base) == Some(false) {
            return Err(format!("{label} fails on the rows at base n = {base}"));
        }
    }
    Ok(count)
}

#[allow(clippy::too_many_arguments)]
fn transport_report(
    transformation: &RecurrenceTransformation,
    label: &str,
    relation: &ImplicitRelation,
    first_base: usize,
    evidence: &str,
    origin: usize,
    q_rows: &[Vec<Q>],
    max_degree: usize,
) -> Result<ImplicitTransportReport, String> {
    let transport = transformation
        .transport_implicit(relation, first_base)
        .map_err(|e| format!("{label}: {e}"))?;
    let r = transport.relation.max_shift();
    let bases = checked_bases(first_base, r, origin, q_rows.len());
    let mut holds = true;
    for base in bases.clone() {
        holds &= transport.relation.holds_at(q_rows, origin, base) != Some(false);
        for (_, condition) in &transport.exceptional_conditions {
            holds &= condition.holds_at(q_rows, origin, base) != Some(false);
        }
    }
    let mut clearing = transport.clearing.clone();
    if let Some(c) = clearing.as_mut() {
        if c.exceptional.status != CERTIFIED_ALL_INDICES {
            c.exceptional.observed_zeros = observed_multiplier_zeros(c, bases.clone(), max_degree);
        }
    }
    let conditions = transport
        .exceptional_conditions
        .iter()
        .map(|(k, condition)| {
            DescribedCondition::new(condition, first_base, evidence, Some(*k), "Q").with_note(
                format!("[x^{k}] of the uncleared residual: H(n,{k}) = 0 at the exceptional index"),
            )
        })
        .collect();
    let note = match &clearing {
        Some(c) if c.exceptional.status != CERTIFIED_ALL_INDICES => {
            "transported and cleared; the exceptional set is not certified, so equivalence with \
             the source relation is not claimed"
        }
        Some(_) => "transported through the diagonal map and cleared by one common multiplier",
        None => "transported by substitution P = T_g Q",
    };
    Ok(ImplicitTransportReport {
        source: label.to_string(),
        relation: DescribedRelation::new(&transport.relation, first_base, evidence, "Q")
            .with_note(note),
        clearing,
        exceptional_index_conditions: conditions,
        reconstruction: transport.reconstruction,
        holds_on_transformed_rows: holds,
        checked_bases: bases.len(),
    })
}

struct SearchContext<'a> {
    source: &'a Source,
    defining: &'a ImplicitRelation,
    defining_first_base: usize,
    map_certified: bool,
    origin: usize,
    q_rows: &'a [Vec<Q>],
}

fn search_report(
    transformation: &RecurrenceTransformation,
    bounds: &TransformSearchBounds,
    context: &SearchContext<'_>,
) -> ImplicitSearchReport {
    let mut report = ImplicitSearchReport {
        evidence: FINITE_FIT_EVIDENCE.to_string(),
        status: "not_found".into(),
        scalar: None,
        relation: None,
        reconstruction: None,
        upgrade_steps: Vec::new(),
        description: None,
        assessment: None,
        reason: None,
    };
    let rows = context.q_rows;
    let origin = context.origin;
    let options = AdaptiveSearchOptions {
        max_rec_len: bounds.max_rec_len,
        max_var_deg: bounds.max_var_deg,
        max_idx_deg: bounds.max_idx_deg,
        max_diff_deg: bounds.max_diff_deg,
        try_denominator: bounds.try_denominator,
        ..AdaptiveSearchOptions::default()
    };
    let found = match find_recurrence_adaptive_rational_with_budget(
        rows,
        &options,
        AdaptiveSearchBudget::limited(bounds.max_candidates),
    ) {
        AdaptiveSearchOutcome::Found(found) => found,
        AdaptiveSearchOutcome::NoRecurrence(_) => {
            report.reason = Some("no rule within the search bounds".into());
            return report;
        }
        AdaptiveSearchOutcome::BudgetExhausted(_) => {
            report.status = "budget_exhausted".into();
            report.reason = Some("candidate budget exhausted".into());
            return report;
        }
    };
    report.status = "found".into();
    let rule = found.recurrence.shift_index(1 - origin as i64);
    let initial = rule.generation_initial_count(origin, rows.len());
    report.scalar = Some(RecurrenceJson::from_recurrence_rational(
        &rule,
        origin,
        &rows[..initial],
        Some(RecurrenceJsonSearch {
            recurrence_text: rule.to_string(),
            source_rows: rows.len(),
            skip_prefix: 0,
            unknowns: found.num_unknowns,
            weighted_unknowns: found.weighted_unknowns,
            equations: found.num_equations,
            fit_polynomials: found.fit_polynomials,
            verification_polynomials: found.verification_polynomials,
            candidates_tried: found.candidates_tried,
            options: RecurrenceOptionsJson::from(&found.opts),
        }),
    ));
    let (fitted, r) = match ImplicitRelation::from_recurrence(&rule) {
        Ok(parts) => parts,
        Err(error) => {
            report.reason = Some(format!("the fit has no implicit form: {error}"));
            return report;
        }
    };
    let fitted_first_base = origin + initial - r;
    report.relation = Some(
        DescribedRelation::new(&fitted, fitted_first_base, FINITE_FIT_EVIDENCE, "Q")
            .with_note("fitted to the transformed rows"),
    );
    let mut steps = Vec::new();
    match upgrade(
        transformation,
        context,
        &rule,
        &fitted,
        fitted_first_base,
        initial,
        &mut steps,
    ) {
        Ok(certificate) => {
            report.evidence = CERTIFIED_BY_RECONSTRUCTION_EVIDENCE.to_string();
            report.reconstruction = Some(certificate);
            let relation = DescribedRelation::new(
                &fitted,
                fitted_first_base,
                CERTIFIED_BY_RECONSTRUCTION_EVIDENCE,
                "Q",
            )
            .with_note("fitted, then proved for every index by the reconstruction certificate");
            report.relation = Some(relation.clone());
            let description = GeneratingDescription {
                first_index: origin,
                seeds: encode_rows(&rows[..initial]),
                relations: vec![relation],
                exceptional_index_conditions: Vec::new(),
                boundary_conditions: Vec::new(),
            };
            report.assessment = assess_generating_description(&description).ok();
            report.description = Some(description);
        }
        Err(reason) => report.reason = Some(format!("not upgraded: {reason}")),
    }
    report.upgrade_steps = steps;
    report
}

fn upgrade(
    transformation: &RecurrenceTransformation,
    context: &SearchContext<'_>,
    rule: &Recurrence,
    fitted: &ImplicitRelation,
    fitted_first_base: usize,
    initial: usize,
    steps: &mut Vec<String>,
) -> Result<ReconstructionCertificate, String> {
    let source = context.source;
    let origin = context.origin;
    if !context.map_certified {
        return Err("the map hypotheses are not certified for all indices".into());
    }
    let t = transformation
        .reconstruction_operator()
        .ok_or("the map Q -> P is not a differential operator with polynomial coefficients")?;
    let certificate = find_reconstruction_certificate(context.defining, &t, fitted).ok_or(
        "no reconstruction certificate: the fit is not proportional to L o T and its leader is \
         not a nonzero constant",
    )?;
    certificate.verify_identity()?;
    steps.push(format!(
        "identity: L o T = sum_nu H_nu o S^nu o M holds in the Ore algebra ({} multipliers)",
        certificate.multipliers.len()
    ));
    let valid_from = origin + source.seeds.len();
    if let Some(f) = &source.rule.denominator {
        let analysis =
            super::algebra::analyze_zeros(&super::algebra::Biv::from_bivar(f), valid_from, None);
        if !analysis.singular_indices.is_empty() {
            return Err(format!(
                "the source leader vanishes identically at targets {:?}",
                analysis.singular_indices
            ));
        }
    }
    steps.push(format!(
        "source uniqueness: the leader of the defining rule is a nonzero polynomial at every \
         target m >= {valid_from}, so P is determined by its seeds"
    ));
    let r_l = context.defining.max_shift();
    let total = valid_from.max(fitted_first_base + r_l) + 1 - origin;
    let generated = rule
        .generate_rows_rational(&context.q_rows[..initial], origin, total.max(initial))
        .map_err(|e| format!("generation from the fit failed: {e}"))?
        .into_iter()
        .map(|r| vec_poly::trim(&r))
        .collect::<Vec<_>>();
    for (j, seed) in source.seeds.iter().enumerate() {
        let index = origin + j;
        let image = t.apply_to_row(&generated[j], index);
        if vec_poly::trim(&image) != vec_poly::trim(seed) {
            return Err(format!(
                "initial agreement fails at the source seed index {index}"
            ));
        }
    }
    steps.push(format!(
        "initial agreement: T Q*_j = P_j at the source seeds j = {origin}..{}, where Q* is \
         generated by the fit from its seeds",
        valid_from - 1
    ));
    let mapped = generated
        .iter()
        .enumerate()
        .map(|(j, row)| t.apply_to_row(row, origin + j))
        .collect::<Vec<_>>();
    let direct = context.defining_first_base..fitted_first_base.max(context.defining_first_base);
    for base in direct.clone() {
        if context.defining.holds_at(&mapped, origin, base) != Some(true) {
            return Err(format!("L (T Q*) fails at the uncovered base n = {base}"));
        }
    }
    steps.push(format!(
        "domain: for bases n >= {} every identity (M Q*)_(n+nu) lies in the fitted domain; the \
         {} earlier bases were checked directly on Q*",
        fitted_first_base.max(context.defining_first_base),
        direct.len()
    ));
    steps.push(format!(
        "conclusion: T Q* satisfies the defining rule from base {} and agrees with its seeds, \
         so T Q* = P; T is injective, so Q* = Q and the fit holds for every base \
         n >= {fitted_first_base}",
        context.defining_first_base
    ));
    Ok(certificate)
}

/// Replay a report exactly from its embedded request, then re-check every
/// reconstruction certificate symbolically and by deriving its identities on
/// the transformed rows.
pub fn verify_recurrence_implicit(report: &ImplicitTransformReport) -> Result<(), String> {
    if report.schema != RECURRENCE_IMPLICIT_SCHEMA {
        return Err(format!("unsupported schema `{}`", report.schema));
    }
    let replay = implicit_transform_recurrence_with_root_domain(
        &report.request,
        report.root_domain.as_ref().map(|r| &r.request),
    )?;
    if &replay != report {
        return Err("report does not match its exact replay".into());
    }
    if !report.reconstruction_verified {
        return Err("exact reconstruction failed".into());
    }
    let rows = report
        .transformed_rows
        .iter()
        .map(|r| decode_poly(r))
        .collect::<Result<Vec<_>, _>>()?;
    let certificates = report
        .transports
        .iter()
        .filter_map(|t| t.reconstruction.as_ref())
        .chain(
            report
                .independent_search
                .as_ref()
                .and_then(|s| s.reconstruction.as_ref()),
        );
    for certificate in certificates {
        certificate.verify_identity()?;
        certificate.replay_on_rows(&rows, report.row_origin)?;
    }
    if report
        .transports
        .iter()
        .any(|t| !t.holds_on_transformed_rows)
    {
        return Err("a transported relation fails on the transformed rows".into());
    }
    let assessment = assess_generating_description(&report.description)?;
    if assessment != report.assessment {
        return Err("the description assessment does not replay".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests;
