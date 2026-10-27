//! Bounded queue search over equivalent polynomial representations of one
//! sequence: alternative recurrences, changes of coordinates, searches in the
//! new coordinates, and certified relations fed into later steps.
//!
//! # Queue
//!
//! [`search_representations`] starts from a source: a defining
//! `polytool.recurrence.v1` rule with seeds and origin, a bundled OEIS family,
//! or rows only (exploratory fits). The source rows are representation 0.
//! The queue is breadth first. For every representation below
//! `limits.max_depth`, each transformation of the request-specified set is
//! applied to its rows:
//!
//! 1. the composite map from the original source is reduced to a certified
//!    normal form (module `representation_map`); a map certified equal to an
//!    existing representation's map is a duplicate (`duplicate_map`) and is
//!    not searched again; rows are never compared for deduplication;
//! 2. the rows are mapped exactly and the source rows are reconstructed by
//!    inverting the whole chain, so the original source stays the reference;
//! 3. known relations are transported: every certified scalar description of
//!    the parent (at most `limits.max_sources_per_edge`, simplest first)
//!    becomes the defining rule of a replayable step,
//!    [`transform_recurrence`] (Möbius class, gamma, symmetric decomposition)
//!    or [`implicit_transform_recurrence`] (diagonal maps, derivative
//!    antecedents, with a reconstruction-certified successor fit); without a
//!    certified parent relation the step maps rows only;
//! 4. for a factorial diagonal map on a constant-index monic certified source,
//!    complete recurrence spaces are certified ([`find_recurrence_spaces`]
//!    with `coefficient_rescaling`) on a bounded menu of supports;
//! 5. the new coordinates are searched by the bounded fit menu
//!    ([`find_recurrence_representations_rational`]).
//!
//! Every certified scalar description of the new representation is then
//! available to its children. A representation whose equivalent generating
//! description coincides with one of another representation (same relations,
//! seeds, first index and conditions) is the same sequence, certified by the
//! uniqueness arguments of both descriptions; it is kept (`equal_sequence_as`)
//! but not expanded.
//!
//! # Labels and ranking
//!
//! Each [`RepresentationDescription`] carries one of three labels:
//! `finite_fit` (fitted to finitely many rows, or depending on a condition
//! checked on rows only), `certified_annihilating_relation` (annihilation
//! proved for every index of the domain, uniqueness not certified) and
//! `equivalent_generating_description` (certified annihilation together with
//! sufficient initial or boundary data and a certified uniqueness argument,
//! [`assess_generating_description`]). Descriptions are ranked by shift length,
//! derivative order, coefficient (`x`) degree, index degree, component count
//! and reconstruction cost (maps to invert plus coefficient conditions), with
//! the label as the final tie break.
//!
//! # Limits and replay
//!
//! Depth, representation count, backend steps, space supports, rows, degrees,
//! fit candidates and elimination work are bounded; every limit that cuts work
//! is reported as a [`RepresentationBudgetEvent`] and makes the status
//! `budget_exhausted`. [`verify_representation_search`] replays the whole
//! search from the embedded request, replays every embedded step and space
//! certificate with its own verifier, re-derives every representation from
//! the original source rows, and re-checks every description on its rows.
//!
//! Related modules: `representation_report` (request and report types),
//! `transform` (registry and transport), `implicit_report`
//! (implicit transport and reconstruction), `space` and `rescaling`
//! (certified complete spaces), `representations` (fit menus).

use super::implicit::{decode_poly, encode_poly};
use super::representation_map::{apply_chain, invert_chain, MapNormalForm};
use super::representation_report::{conventions, label_rank};
use super::*;
use crate::vec_poly;
use num_traits::Zero;
use std::collections::VecDeque;

type Q = BigRational;

// ---------------------------------------------------------------------------
// Source resolution
// ---------------------------------------------------------------------------

struct Source {
    kind: &'static str,
    rule: Option<(Recurrence, Vec<Vec<Q>>, RecurrenceJson)>,
    origin: usize,
    rows: Vec<Vec<Q>>,
}

fn validate(request: &RepresentationSearchRequest) -> Result<(), String> {
    let limits = &request.limits;
    if request.transformations.len() > MAX_REPRESENTATION_TRANSFORMATIONS {
        return Err(format!(
            "at most {MAX_REPRESENTATION_TRANSFORMATIONS} transformations are supported"
        ));
    }
    for (i, t) in request.transformations.iter().enumerate() {
        t.validate()
            .map_err(|e| format!("transformations[{i}]: {e}"))?;
    }
    let bounds = &request.search;
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
    if limits.max_depth > 16 || limits.max_representations > 4096 {
        return Err("max_depth is at most 16 and max_representations at most 4096".into());
    }
    if let Some(root) = &request.root_domain {
        root.validate()?;
    }
    if request.space_shapes.len() > limits.max_space_shapes {
        return Err("more space_shapes than max_space_shapes".into());
    }
    Ok(())
}

fn resolve(request: &RepresentationSearchRequest) -> Result<Source, String> {
    let limits = &request.limits;
    let defining = match (&request.source_recurrence, &request.oeis_id) {
        (Some(_), Some(_)) => {
            return Err("give at most one of source_recurrence and oeis_id".into())
        }
        (Some(json), None) => {
            let (rule, first, seeds) = json.to_recurrence_parts()?;
            Some(("defining_rule", rule, first, seeds))
        }
        (None, Some(id)) => {
            let definition =
                crate::oeis::by_id(id).ok_or_else(|| format!("unknown OEIS family `{id}`"))?;
            let (rule, seeds, first_row) = definition
                .recurrence_export_parts()
                .map_err(|e| e.to_string())?;
            let first = usize::try_from(first_row)
                .map_err(|_| format!("OEIS family `{id}` starts at a negative row index"))?;
            Some(("oeis", rule, first, seeds))
        }
        (None, None) => None,
    };
    let explicit = request
        .rows
        .as_ref()
        .map(|rows| {
            rows.iter()
                .map(|r| decode_poly(r))
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let (kind, rule, origin, rows) = match defining {
        Some((kind, rule, first, seeds)) => {
            if rule.terms.iter().any(|t| t.offset == 0) {
                return Err("source rule terms need positive offsets".into());
            }
            if seeds.is_empty() || rule.max_offset() > seeds.len() {
                return Err("source rule needs at least max-offset seeds".into());
            }
            if rule
                .terms
                .iter()
                .any(|t| t.deriv_order > TransformLimits::default().max_rule_derivative_order)
            {
                return Err("source rule derivative order is too large".into());
            }
            if request.row_origin.is_some_and(|o| o != first) {
                return Err(format!(
                    "row_origin differs from the source first index {first}"
                ));
            }
            let count = explicit.as_ref().map_or(
                request
                    .row_count
                    .unwrap_or(DEFAULT_REPRESENTATION_ROW_COUNT),
                Vec::len,
            );
            if count > limits.max_rows {
                return Err(format!(
                    "{count} rows exceed max_rows = {}",
                    limits.max_rows
                ));
            }
            let generated = rule
                .generate_rows_rational(&seeds, first, count.max(seeds.len()))
                .map_err(|e| format!("source rule generation failed: {e}"))?
                .into_iter()
                .map(|r| vec_poly::trim(&r))
                .collect::<Vec<_>>();
            if let Some(rows) = &explicit {
                if rows[..] != generated[..rows.len()] {
                    return Err("rows disagree with the rows generated by the source rule".into());
                }
            }
            let json = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
            (kind, Some((rule, seeds, json)), first, generated)
        }
        None => {
            let rows = explicit.ok_or("give source_recurrence, oeis_id or rows")?;
            ("rows_only", None, request.row_origin.unwrap_or(0), rows)
        }
    };
    if rows.is_empty() {
        return Err("at least one row is required".into());
    }
    if rows.len() > limits.max_rows {
        return Err(format!(
            "{} rows exceed max_rows = {}",
            rows.len(),
            limits.max_rows
        ));
    }
    if let Some(degree) = rows.iter().filter_map(|r| vec_poly::degree(r)).max() {
        if degree > limits.max_degree {
            return Err(format!("row degree {degree} exceeds max_degree"));
        }
    }
    Ok(Source {
        kind,
        rule,
        origin,
        rows,
    })
}

// ---------------------------------------------------------------------------
// Description helpers
// ---------------------------------------------------------------------------

fn nonzero(text: &str) -> bool {
    parse_rational_coeff(text).map_or(true, |q| !q.is_zero())
}

/// `(index degree, x degree)` of the nonzero entries of `coeffs[n][x]`.
fn json_degrees(coeffs: &[Vec<String>]) -> (usize, usize) {
    let mut degrees = (0, 0);
    for (i, row) in coeffs.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if nonzero(c) {
                degrees = (degrees.0.max(i), degrees.1.max(j));
            }
        }
    }
    degrees
}

fn bivar_degrees(p: &BivarPoly) -> (usize, usize) {
    let mut degrees = (0, 0);
    for (i, row) in p.coeffs.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if !c.is_zero() {
                degrees = (degrees.0.max(i), degrees.1.max(j));
            }
        }
    }
    degrees
}

fn rule_cost(rule: &Recurrence) -> DescriptionCost {
    let mut cost = DescriptionCost {
        shift_length: rule.max_offset(),
        derivative_order: rule.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0),
        coefficient_degree: 0,
        index_degree: 0,
        component_count: 1,
        reconstruction_cost: 0,
    };
    let polys = rule
        .terms
        .iter()
        .map(|t| &t.coeff)
        .chain(rule.denominator.iter())
        .chain(rule.inhomogeneous.iter());
    for p in polys {
        let (n, x) = bivar_degrees(p);
        cost.index_degree = cost.index_degree.max(n);
        cost.coefficient_degree = cost.coefficient_degree.max(x);
    }
    cost
}

fn relation_cost(relation: &ImplicitRelationJson) -> DescriptionCost {
    let mut cost = DescriptionCost {
        shift_length: 0,
        derivative_order: 0,
        coefficient_degree: 0,
        index_degree: 0,
        component_count: 1,
        reconstruction_cost: 0,
    };
    for term in &relation.terms {
        let (n, x) = json_degrees(&term.coeff.coeffs);
        if n == 0 && x == 0 && !term.coeff.coeffs.iter().flatten().any(|c| nonzero(c)) {
            continue;
        }
        cost.shift_length = cost.shift_length.max(term.shift);
        cost.derivative_order = cost.derivative_order.max(term.deriv_order);
        cost.index_degree = cost.index_degree.max(n);
        cost.coefficient_degree = cost.coefficient_degree.max(x);
    }
    cost
}

fn vector_cost(json: &TransformVectorRecurrenceJson) -> DescriptionCost {
    let mut cost = DescriptionCost {
        shift_length: 1,
        derivative_order: 0,
        coefficient_degree: 0,
        index_degree: 0,
        component_count: json.initial_state.len(),
        reconstruction_cost: 0,
    };
    let entries = json
        .matrix
        .iter()
        .flatten()
        .flat_map(|op| op.iter().enumerate())
        .map(|(d, p)| (Some(d), p))
        .chain(std::iter::once((None, &json.denominator)))
        .chain(json.forcing.iter().flatten().map(|p| (None, p)));
    for (derivative, p) in entries {
        if !p.coeffs.iter().flatten().any(|c| nonzero(c)) {
            continue;
        }
        let (n, x) = json_degrees(&p.coeffs);
        cost.index_degree = cost.index_degree.max(n);
        cost.coefficient_degree = cost.coefficient_degree.max(x);
        if let Some(d) = derivative {
            cost.derivative_order = cost.derivative_order.max(d);
        }
    }
    cost
}

struct ScalarInput<'a> {
    kind: &'a str,
    rule: &'a Recurrence,
    json: RecurrenceJson,
    evidence: &'a str,
    valid_from: usize,
    bounded_domain: bool,
    symbol: &'a str,
    provenance: String,
}

fn scalar_description(input: ScalarInput<'_>) -> RepresentationDescription {
    let certified = evidence_is_certified(input.evidence);
    let mut item = RepresentationDescription {
        label: FINITE_FIT_LABEL.to_string(),
        kind: input.kind.to_string(),
        evidence: input.evidence.to_string(),
        text: input.rule.to_string(),
        scalar: Some(RecurrenceJson {
            search: None,
            ..input.json
        }),
        vector: None,
        description: None,
        assessment: None,
        cost: rule_cost(input.rule),
        provenance: input.provenance,
        also_found_by: Vec::new(),
        note: None,
    };
    let json = item.scalar.as_ref().expect("scalar");
    let assessed = ImplicitRelation::from_recurrence(input.rule).and_then(|(relation, r)| {
        let first_base = input
            .valid_from
            .checked_sub(r)
            .ok_or_else(|| "the rule starts before its leader shift".to_string())?;
        let description = GeneratingDescription {
            first_index: json.first_index,
            seeds: json.initial_polynomials.clone(),
            relations: vec![DescribedRelation::new(
                &relation,
                first_base,
                input.evidence,
                input.symbol,
            )],
            exceptional_index_conditions: Vec::new(),
            boundary_conditions: Vec::new(),
        };
        let assessment = assess_generating_description(&description)?;
        Ok((description, assessment))
    });
    match assessed {
        Ok((description, assessment)) => {
            item.label = representation_label(&assessment.status).to_string();
            if input.bounded_domain && item.label == EQUIVALENT_DESCRIPTION_LABEL {
                item.label = CERTIFIED_RELATION_LABEL.to_string();
                item.note = Some("the map's index domain is bounded above".into());
            }
            item.description = Some(description);
            item.assessment = Some(assessment);
        }
        Err(reason) => {
            if certified {
                item.label = CERTIFIED_RELATION_LABEL.to_string();
            }
            item.note = Some(format!(
                "not assessed as a generating description: {reason}"
            ));
        }
    }
    item
}

fn single_relation(description: &RepresentationDescription) -> Option<ImplicitRelation> {
    let d = description.description.as_ref()?;
    match d.relations.as_slice() {
        [only] => ImplicitRelation::from_json(&only.relation).ok(),
        _ => None,
    }
}

fn same_relation(a: &RepresentationDescription, b: &RepresentationDescription) -> bool {
    if let (Some(ra), Some(rb)) = (single_relation(a), single_relation(b)) {
        return ra.same_relation(&rb);
    }
    if let (Some(sa), Some(sb)) = (&a.scalar, &b.scalar) {
        if let (Ok((ra, _, _)), Ok((rb, _, _))) =
            (sa.to_recurrence_parts(), sb.to_recurrence_parts())
        {
            return same_canonical_recurrence(&ra, &rb);
        }
    }
    if let (Some(va), Some(vb)) = (&a.vector, &b.vector) {
        return va.matrix == vb.matrix && va.denominator == vb.denominator;
    }
    false
}

fn rank_key(d: &RepresentationDescription) -> (DescriptionCost, usize) {
    (d.cost.clone(), label_rank(&d.label))
}

/// Merge descriptions of the same relation (keeping the strongest label) and rank.
fn merge_descriptions(items: Vec<RepresentationDescription>) -> Vec<RepresentationDescription> {
    let mut merged: Vec<RepresentationDescription> = Vec::new();
    for item in items {
        match merged.iter_mut().find(|m| same_relation(m, &item)) {
            Some(existing) => {
                let (mut kept, other) = if label_rank(&item.label) < label_rank(&existing.label) {
                    (item, existing.clone())
                } else {
                    (existing.clone(), item)
                };
                for provenance in std::iter::once(other.provenance).chain(other.also_found_by) {
                    if provenance != kept.provenance && !kept.also_found_by.contains(&provenance) {
                        kept.also_found_by.push(provenance);
                    }
                }
                *existing = kept;
            }
            None => merged.push(item),
        }
    }
    merged.sort_by_key(rank_key);
    merged
}

fn identical_descriptions(a: &GeneratingDescription, b: &GeneratingDescription) -> bool {
    let relations = a.relations.len() == b.relations.len()
        && a.relations.iter().zip(&b.relations).all(|(x, y)| {
            x.first_base == y.first_base
                && match (
                    ImplicitRelation::from_json(&x.relation),
                    ImplicitRelation::from_json(&y.relation),
                ) {
                    (Ok(rx), Ok(ry)) => rx.same_relation(&ry),
                    _ => false,
                }
        });
    let conditions = |p: &[DescribedCondition], q: &[DescribedCondition]| {
        p.len() == q.len()
            && p.iter()
                .zip(q)
                .all(|(x, y)| x.condition == y.condition && x.first_base == y.first_base)
    };
    relations
        && a.first_index == b.first_index
        && a.seeds == b.seeds
        && conditions(
            &a.exceptional_index_conditions,
            &b.exceptional_index_conditions,
        )
        && conditions(&a.boundary_conditions, &b.boundary_conditions)
}

fn find_equal_sequence(
    nodes: &[RepresentationNode],
    node: &RepresentationNode,
) -> Option<RepresentationEquality> {
    let equivalent = |d: &RepresentationDescription| {
        d.label == EQUIVALENT_DESCRIPTION_LABEL && d.description.is_some()
    };
    for other in nodes {
        if other.components.len() != node.components.len() || other.row_origin != node.row_origin {
            continue;
        }
        for (i, mine) in node.descriptions.iter().enumerate() {
            if !equivalent(mine) {
                continue;
            }
            for (j, theirs) in other.descriptions.iter().enumerate() {
                if equivalent(theirs)
                    && identical_descriptions(
                        mine.description.as_ref().expect("description"),
                        theirs.description.as_ref().expect("description"),
                    )
                {
                    return Some(RepresentationEquality {
                        representation: other.id,
                        other_description: j,
                        description: i,
                        certificate: "both descriptions are equivalent generating descriptions \
                                      with the same relations, seeds, first index and \
                                      conditions; each determines its sequence uniquely"
                            .into(),
                    });
                }
            }
        }
    }
    None
}

fn encode_states(states: &[Vec<Vec<Q>>]) -> Vec<Vec<Vec<String>>> {
    states
        .iter()
        .map(|s| s.iter().map(|c| encode_poly(c)).collect())
        .collect()
}

fn decode_states(states: &[Vec<Vec<String>>]) -> Result<Vec<Vec<Vec<Q>>>, String> {
    states
        .iter()
        .map(|s| s.iter().map(|c| decode_poly(c)).collect())
        .collect()
}

fn scalar_rows(node: &RepresentationNode) -> Result<Vec<Vec<Q>>, String> {
    node.states
        .iter()
        .map(|s| decode_poly(&s[0]))
        .collect::<Result<Vec<_>, _>>()
}

fn component_names(t: Option<&RecurrenceTransformation>, chain_len: usize) -> Vec<String> {
    match t {
        None => vec![if chain_len == 0 { "P" } else { "Q" }.to_string()],
        Some(t) => t.record().components.into_iter().map(|c| c.name).collect(),
    }
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

struct Search<'a> {
    request: &'a RepresentationSearchRequest,
    origin: usize,
    source_rows: Vec<Vec<Q>>,
    steps_used: usize,
    budget: Vec<RepresentationBudgetEvent>,
}

impl Search<'_> {
    fn event(&mut self, limit: &str, representation: Option<usize>, detail: String) {
        self.budget.push(RepresentationBudgetEvent {
            limit: limit.to_string(),
            representation,
            detail,
        });
    }

    fn take_step(&mut self, representation: usize, what: &str) -> bool {
        if self.steps_used >= self.request.limits.max_steps {
            self.event(
                "max_steps",
                Some(representation),
                format!("{what} not run: the step budget is spent"),
            );
            return false;
        }
        self.steps_used += 1;
        true
    }

    fn fit_options(&self) -> AdaptiveSearchOptions {
        let bounds = &self.request.search;
        AdaptiveSearchOptions {
            max_rec_len: bounds.max_rec_len.max(1),
            max_var_deg: bounds.max_var_deg,
            max_idx_deg: bounds.max_idx_deg,
            max_diff_deg: bounds.max_diff_deg,
            try_denominator: bounds.try_denominator,
            // Left factors are polynomials in n only, as in the space menu:
            // generation then never divides by a polynomial in x.
            max_denom_var_deg: 0,
            max_denom_idx_deg: bounds.max_idx_deg.max(1),
            ..AdaptiveSearchOptions::default()
        }
    }

    /// The fit menu on scalar rows.
    fn fits(
        &mut self,
        id: usize,
        rows: &[Vec<Q>],
        symbol: &str,
    ) -> (
        Option<RepresentationFitSummary>,
        Vec<RepresentationDescription>,
    ) {
        if !self.request.fits {
            return (None, Vec::new());
        }
        let summary = |status: &str, reason: Option<String>| RepresentationFitSummary {
            status: status.into(),
            configurations_considered: 0,
            fits_found: 0,
            reason,
        };
        if rows.len() < 3 {
            return (Some(summary("insufficient_rows", None)), Vec::new());
        }
        if !self.take_step(id, "fit menu") {
            return (Some(summary("skipped_budget", None)), Vec::new());
        }
        let menu = match find_recurrence_representations_rational(
            rows,
            &self.fit_options(),
            AdaptiveSearchBudget::limited(self.request.search.max_candidates),
        ) {
            Ok(menu) => menu,
            Err(error) => return (Some(summary("failed", Some(error))), Vec::new()),
        };
        let exhausted = menu.status == RepresentationSearchStatus::BudgetExhausted;
        if exhausted {
            self.event(
                "fit_candidates",
                Some(id),
                "the fit menu exhausted its candidate budget; unresolved cells are incomplete"
                    .into(),
            );
        }
        let origin = self.origin;
        let mut items = Vec::new();
        for rep in &menu.representations {
            let Ok((rule, _, _)) = rep.recurrence_json.to_recurrence_parts() else {
                continue;
            };
            let rule = rule.shift_index(1 - origin as i64);
            let count = rep
                .recurrence_json
                .initial_polynomials
                .len()
                .max(rule.generation_initial_count(origin, rows.len()))
                .min(rows.len());
            let holds = (origin + count..origin + rows.len())
                .all(|target| transformed_rule_holds_at(&rule, rows, origin, target) == Some(true));
            if !holds {
                continue;
            }
            let json =
                RecurrenceJson::from_recurrence_rational(&rule, origin, &rows[..count], None);
            items.push(scalar_description(ScalarInput {
                kind: FINITE_FIT_LABEL,
                rule: &rule,
                json,
                evidence: FINITE_FIT_EVIDENCE,
                valid_from: origin + count,
                bounded_domain: false,
                symbol,
                provenance: "fits".into(),
            }));
        }
        let summary = RepresentationFitSummary {
            status: if exhausted {
                "budget_exhausted"
            } else {
                "complete"
            }
            .into(),
            configurations_considered: menu.configurations_considered,
            fits_found: menu.representations.len(),
            reason: None,
        };
        (Some(summary), items)
    }

    fn root(&mut self, source: &Source) -> RepresentationNode {
        let mut descriptions = Vec::new();
        if let Some((rule, seeds, json)) = &source.rule {
            descriptions.push(scalar_description(ScalarInput {
                kind: "defining_rule",
                rule,
                json: json.clone(),
                evidence: DEFINING_EVIDENCE,
                valid_from: source.origin + seeds.len(),
                bounded_domain: false,
                symbol: "P",
                provenance: "source".into(),
            }));
        }
        let (fits, fitted) = self.fits(0, &source.rows, "P");
        descriptions.extend(fitted);
        let form = MapNormalForm::from_chain(&[], source.origin).expect("identity");
        RepresentationNode {
            id: 0,
            depth: 0,
            parent: None,
            transformation: None,
            chain: Vec::new(),
            map: form.record(),
            components: component_names(None, 0),
            row_origin: source.origin,
            states: encode_states(
                &source
                    .rows
                    .iter()
                    .map(|r| vec![r.clone()])
                    .collect::<Vec<_>>(),
            ),
            reconstruction_verified: true,
            map_conditions: Vec::new(),
            descriptions: merge_descriptions(descriptions),
            steps: Vec::new(),
            spaces: Vec::new(),
            fits,
            equal_sequence_as: None,
            expanded: false,
            notes: Vec::new(),
        }
    }

    /// Certified scalar descriptions of the parent usable as defining rules.
    fn source_descriptions(&self, parent: &RepresentationNode) -> Vec<usize> {
        parent
            .descriptions
            .iter()
            .enumerate()
            .filter(|(_, d)| {
                d.label != FINITE_FIT_LABEL
                    && d.scalar
                        .as_ref()
                        .is_some_and(|json| json.first_index == self.origin)
            })
            .map(|(i, _)| i)
            .take(self.request.limits.max_sources_per_edge)
            .collect()
    }

    #[allow(clippy::too_many_arguments)]
    fn run_step(
        &mut self,
        id: usize,
        t: &RecurrenceTransformation,
        parent: &RepresentationNode,
        source: Option<usize>,
        states: &[Vec<Vec<Q>>],
        items: &mut Vec<RepresentationDescription>,
        conditions: &mut Option<Vec<TransformConditionStatus>>,
    ) -> RepresentationStep {
        let mut step = RepresentationStep {
            source_description: source,
            transform: None,
            implicit: None,
            error: None,
        };
        let row_count = states.len();
        let limits = self.request.limits.transform_limits();
        let source_json = source.map(|i| RecurrenceJson {
            search: None,
            ..parent.descriptions[i]
                .scalar
                .clone()
                .expect("scalar source")
        });
        let implicit_map = matches!(
            t,
            RecurrenceTransformation::Diagonal { .. }
                | RecurrenceTransformation::DerivativeAntecedent { .. }
        );
        let expected = encode_states(states);
        if implicit_map && source_json.is_some() {
            let request = ImplicitTransformRequest {
                transformation: t.clone(),
                source_recurrence: source_json,
                oeis_id: None,
                row_count: Some(row_count),
                source_relations: Vec::new(),
                target_relations: Vec::new(),
                boundary_conditions: Vec::new(),
                search: Some(self.request.search.clone()),
                limits,
            };
            match implicit_transform_recurrence(&request) {
                Ok(report) => {
                    let rows = expected.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
                    let n = report.transformed_rows.len().min(rows.len());
                    if report.transformed_rows[..n] != rows[..n] {
                        step.error = Some("the step's rows differ from the mapped rows".into());
                    } else {
                        if conditions.is_none() {
                            *conditions = Some(report.map_conditions.clone());
                        }
                        self.implicit_items(id, &report, items);
                        step.implicit = Some(report);
                    }
                }
                Err(error) => step.error = Some(error),
            }
        } else {
            let request = RecurrenceTransformRequest {
                transformation: t.clone(),
                source_recurrence: source_json.clone(),
                oeis_id: None,
                rows: source_json.is_none().then(|| {
                    scalar_rows(parent)
                        .unwrap_or_default()
                        .iter()
                        .map(|r| encode_poly(r))
                        .collect()
                }),
                row_origin: source_json.is_none().then_some(self.origin),
                row_count: source_json.is_some().then_some(row_count),
                search: None,
                limits,
                root_domain: self.request.root_domain.clone(),
            };
            match transform_recurrence(&request) {
                Ok(report) => {
                    let n = report.transformed_states.len().min(expected.len());
                    if report.transformed_states[..n] != expected[..n] {
                        step.error = Some("the step's rows differ from the mapped rows".into());
                    } else {
                        if conditions.is_none() {
                            *conditions = Some(report.conditions.clone());
                        }
                        self.transform_items(&report, source.is_some(), items);
                        step.transform = Some(report);
                    }
                }
                Err(error) => step.error = Some(error),
            }
        }
        step
    }

    fn transform_items(
        &self,
        report: &RecurrenceTransformReport,
        certified_source: bool,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let Some(transported) = &report.transported else {
            return;
        };
        let certified = certified_source
            && report.reconstruction_verified
            && transported.holds_on_transformed_rows
            && report
                .conditions
                .iter()
                .all(|c| c.status == CERTIFIED_ALL_INDICES);
        let evidence = if certified {
            TRANSPORTED_IDENTITY_EVIDENCE
        } else {
            CHECKED_SUPPLIED_ROWS
        };
        let provenance = "steps[STEP].transform.transported".to_string();
        if let Some(json) = &transported.scalar {
            if let Ok((rule, _, _)) = json.to_recurrence_parts() {
                items.push(scalar_description(ScalarInput {
                    kind: "transported_rule",
                    rule: &rule,
                    json: json.clone(),
                    evidence,
                    valid_from: transported.valid_from_index,
                    bounded_domain: transported.valid_to_index.is_some(),
                    symbol: "Q",
                    provenance,
                }));
            }
        } else if let Some(vector) = &transported.vector {
            items.push(RepresentationDescription {
                label: if certified {
                    CERTIFIED_RELATION_LABEL
                } else {
                    FINITE_FIT_LABEL
                }
                .into(),
                kind: "transported_vector_system".into(),
                evidence: evidence.into(),
                text: transported.text.clone(),
                scalar: None,
                vector: Some(vector.clone()),
                description: None,
                assessment: None,
                cost: vector_cost(vector),
                provenance,
                also_found_by: Vec::new(),
                note: Some(
                    "the complete component system; uniqueness of vector systems is not \
                     assessed"
                        .into(),
                ),
            });
        }
    }

    fn implicit_items(
        &mut self,
        id: usize,
        report: &ImplicitTransformReport,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let description = &report.description;
        let relation = &description.relations[0];
        let cost = relation_cost(&relation.relation);
        items.push(RepresentationDescription {
            label: representation_label(&report.assessment.status).into(),
            kind: "implicit_transport".into(),
            evidence: relation.evidence.clone(),
            text: relation.text.clone(),
            scalar: None,
            vector: None,
            description: Some(description.clone()),
            assessment: Some(report.assessment.clone()),
            cost,
            provenance: "steps[STEP].implicit.description".into(),
            also_found_by: Vec::new(),
            note: relation.note.clone(),
        });
        let Some(search) = &report.independent_search else {
            return;
        };
        if search.status == "budget_exhausted" {
            self.event(
                "implicit_fit_candidates",
                Some(id),
                "the successor fit of an implicit step exhausted its candidate budget".into(),
            );
        }
        let Some(json) = &search.scalar else {
            return;
        };
        let Ok((rule, _, _)) = json.to_recurrence_parts() else {
            return;
        };
        if search.evidence == CERTIFIED_BY_RECONSTRUCTION_EVIDENCE {
            if let (Some(description), Some(assessment)) = (&search.description, &search.assessment)
            {
                items.push(RepresentationDescription {
                    label: representation_label(&assessment.status).into(),
                    kind: "reconstruction_certified_fit".into(),
                    evidence: CERTIFIED_BY_RECONSTRUCTION_EVIDENCE.into(),
                    text: rule.to_string(),
                    scalar: Some(RecurrenceJson {
                        search: None,
                        ..json.clone()
                    }),
                    vector: None,
                    description: Some(description.clone()),
                    assessment: Some(assessment.clone()),
                    cost: rule_cost(&rule),
                    provenance: "steps[STEP].implicit.independent_search".into(),
                    also_found_by: Vec::new(),
                    note: Some(
                        "fitted, then proved for every index by a reconstruction certificate"
                            .into(),
                    ),
                });
                return;
            }
        }
        let valid_from = json.first_index + json.initial_polynomials.len();
        items.push(scalar_description(ScalarInput {
            kind: FINITE_FIT_LABEL,
            rule: &rule,
            json: json.clone(),
            evidence: FINITE_FIT_EVIDENCE,
            valid_from,
            bounded_domain: false,
            symbol: "Q",
            provenance: "steps[STEP].implicit.independent_search".into(),
        }));
    }

    /// Complete spaces after a factorial diagonal map.
    fn spaces(
        &mut self,
        id: usize,
        power: usize,
        parent: &RepresentationNode,
        node: &mut RepresentationNode,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let source = parent.descriptions.iter().find_map(|d| {
            let json = d.scalar.as_ref()?;
            (d.label != FINITE_FIT_LABEL
                && json.first_index == self.origin
                && ConstantLagSource::from_recurrence_json(json).is_ok())
            .then(|| RecurrenceJson {
                search: None,
                ..json.clone()
            })
        });
        let Some(source) = source else {
            node.notes.push(
                "complete spaces need a certified constant-index monic pure-lag rule of the \
                 parent; none is available"
                    .into(),
            );
            return;
        };
        let mut tried = 0usize;
        let mut run = |search: &mut Self,
                       node: &mut RepresentationNode,
                       shape: RecurrenceSpaceShape|
         -> Option<usize> {
            if tried >= search.request.limits.max_space_shapes {
                if tried == search.request.limits.max_space_shapes {
                    search.event(
                        "max_space_shapes",
                        Some(id),
                        "further complete-space supports were not tried".into(),
                    );
                    tried += 1;
                }
                return None;
            }
            tried += 1;
            if !search.take_step(id, "complete-space support") {
                return None;
            }
            let request = RecurrenceSpaceRequest {
                source: None,
                source_recurrence: Some(source.clone()),
                indexed_source: None,
                coefficient_rescaling: Some(RecurrenceCoefficientRescaling::FactorialPower {
                    power,
                }),
                include_leader_coverage: false,
                rows: None,
                row_origin: search.origin,
                domain_start: search.origin,
                shapes: vec![shape],
                limits: search.request.limits.space_limits(),
            };
            match find_recurrence_spaces(&request) {
                Ok(report) => {
                    if report.results[0].status == "incomplete_budget" {
                        search.event(
                            "space_operations",
                            Some(id),
                            format!(
                                "a complete-space support is incomplete: {}",
                                report.results[0].reason.clone().unwrap_or_default()
                            ),
                        );
                    }
                    node.spaces.push(report);
                    Some(node.spaces.len() - 1)
                }
                Err(error) => {
                    node.notes
                        .push(format!("complete-space support failed: {error}"));
                    None
                }
            }
        };
        if !self.request.space_shapes.is_empty() {
            for shape in self.request.space_shapes.clone() {
                run(self, node, shape);
            }
        } else {
            let bounds = self.request.search.clone();
            let mut pairs = (1..=bounds.max_rec_len.max(1))
                .flat_map(|l| (0..=bounds.max_diff_deg).map(move |j| (l, j)))
                .collect::<Vec<_>>();
            pairs.sort_by_key(|&(l, j)| (l + j, l));
            let top = bounds.max_idx_deg.max(bounds.max_var_deg);
            let mut successes: Vec<(usize, usize)> = Vec::new();
            'pairs: for (l, j) in pairs {
                if successes.iter().any(|&(sl, sj)| sl <= l && sj <= j) {
                    continue;
                }
                for d in 0..=top {
                    let u = d.min(bounds.max_idx_deg);
                    let shape = RecurrenceSpaceShape::Box {
                        length: l,
                        derivative_order: j,
                        index_degree: u,
                        variable_degree: d.min(bounds.max_var_deg),
                        leading_index_degree: u,
                        leading_variable_degree: 0,
                    };
                    let Some(index) = run(self, node, shape) else {
                        break 'pairs;
                    };
                    match node.spaces[index].results[0].status.as_str() {
                        "successor_exists" => {
                            successes.push((l, j));
                            break;
                        }
                        "absent_in_shape" => {}
                        _ => break,
                    }
                }
            }
        }
        for (s, report) in node.spaces.iter().enumerate() {
            let result = &report.results[0];
            if result.status != "successor_exists" || result.evidence != "conditional_on_source" {
                continue;
            }
            let Some(bound) = &result.rescaling_bound else {
                continue;
            };
            for (r, rep) in result.representatives.iter().enumerate() {
                let Ok(rule) = rep.to_recurrence(result) else {
                    continue;
                };
                let json = RecurrenceJson {
                    schema: RECURRENCE_JSON_SCHEMA.to_string(),
                    first_index: self.origin,
                    initial_polynomials: bound.initial_rows.clone(),
                    recurrence: RecurrenceJson::from_recurrence_rational(&rule, 0, &[], None)
                        .recurrence,
                    search: None,
                };
                items.push(scalar_description(ScalarInput {
                    kind: "certified_space",
                    rule: &rule,
                    json,
                    evidence: CERTIFIED_COMPLETE_SPACE_EVIDENCE,
                    valid_from: self.origin + result.target_lag,
                    bounded_domain: false,
                    symbol: "Q",
                    provenance: format!("spaces[{s}].results[0].representatives[{r}]"),
                }));
            }
        }
    }

    fn child(
        &mut self,
        nodes: &[RepresentationNode],
        parent_id: usize,
        ti: usize,
        form: &MapNormalForm,
    ) -> Result<RepresentationNode, String> {
        let request = self.request;
        let t = &request.transformations[ti];
        let parent = &nodes[parent_id];
        let id = nodes.len();
        let origin = self.origin;
        let parent_rows = scalar_rows(parent)?;
        let states = parent_rows
            .iter()
            .enumerate()
            .map(|(i, row)| t.map_row(origin + i, row))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| format!("{}: {e}", t.kind_name()))?;
        if let Some(degree) = states
            .iter()
            .flatten()
            .filter_map(|c| vec_poly::degree(c))
            .max()
        {
            if degree > request.limits.max_degree {
                return Err(format!("mapped row degree {degree} exceeds max_degree"));
            }
        }
        let mut chain = parent.chain.clone();
        chain.push(ti);
        let chain_refs = chain
            .iter()
            .map(|&i| &request.transformations[i])
            .collect::<Vec<_>>();
        let reconstruction_verified =
            invert_chain(&chain_refs, origin, &states).is_ok_and(|back| back == self.source_rows);
        if !reconstruction_verified {
            return Err("inverting the chain does not reproduce the source rows".into());
        }
        let mut node = RepresentationNode {
            id,
            depth: parent.depth + 1,
            parent: Some(parent_id),
            transformation: Some(ti),
            chain,
            map: form.record(),
            components: component_names(Some(t), 1),
            row_origin: origin,
            states: encode_states(&states),
            reconstruction_verified,
            map_conditions: Vec::new(),
            descriptions: Vec::new(),
            steps: Vec::new(),
            spaces: Vec::new(),
            fits: None,
            equal_sequence_as: None,
            expanded: false,
            notes: Vec::new(),
        };
        let mut items = Vec::new();
        let mut conditions = None;
        let sources = self.source_descriptions(parent);
        let step_sources = if sources.is_empty() {
            vec![None]
        } else {
            sources.into_iter().map(Some).collect()
        };
        for source in step_sources {
            if !self.take_step(id, "transport step") {
                break;
            }
            let before = items.len();
            let step = self.run_step(id, t, parent, source, &states, &mut items, &mut conditions);
            let index = node.steps.len();
            for item in &mut items[before..] {
                item.provenance = item.provenance.replace("STEP", &index.to_string());
            }
            node.steps.push(step);
        }
        node.map_conditions = conditions.unwrap_or_default();
        if let RecurrenceTransformation::Diagonal {
            weight: DiagonalWeight::FactorialPower { power },
        } = t
        {
            if request.certified_spaces {
                self.spaces(id, *power, parent, &mut node, &mut items);
            }
        }
        if node.components.len() == 1 {
            let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
            let (fits, fitted) = self.fits(id, &rows, "Q");
            node.fits = fits;
            items.extend(fitted);
        } else {
            node.notes
                .push("vector representations are not fitted or transformed further".into());
        }
        let segments = form.segment_count();
        for item in &mut items {
            let conditions = item.description.as_ref().map_or(0, |d| {
                d.exceptional_index_conditions.len() + d.boundary_conditions.len()
            });
            item.cost.reconstruction_cost = segments + conditions;
        }
        node.descriptions = merge_descriptions(items);
        Ok(node)
    }
}

/// Run the bounded representation search.
///
/// Errors are returned for invalid requests (malformed maps or sources,
/// limits, contradictory rows); a map that fails on some row is recorded as a
/// rejected edge, and every limit that cuts work is a budget event.
pub fn search_representations(
    request: &RepresentationSearchRequest,
) -> Result<RepresentationSearchReport, String> {
    validate(request)?;
    let source = resolve(request)?;
    let mut search = Search {
        request,
        origin: source.origin,
        source_rows: source.rows.clone(),
        steps_used: 0,
        budget: Vec::new(),
    };
    let mut nodes = vec![search.root(&source)];
    let mut forms = vec![MapNormalForm::from_chain(&[], source.origin)?];
    let mut edges = Vec::new();
    let mut depth_frontier = Vec::new();
    let mut queue = VecDeque::from([0usize]);
    while let Some(id) = queue.pop_front() {
        if let Some(other) = nodes[id]
            .equal_sequence_as
            .as_ref()
            .map(|e| e.representation)
        {
            nodes[id].notes.push(format!(
                "not expanded: the sequence equals representation {other}"
            ));
            continue;
        }
        if nodes[id].components.len() != 1 {
            continue;
        }
        if nodes[id].depth >= request.limits.max_depth {
            if !request.transformations.is_empty() {
                depth_frontier.push(id);
            }
            continue;
        }
        nodes[id].expanded = true;
        for ti in 0..request.transformations.len() {
            let mut edge = RepresentationEdge {
                parent: id,
                transformation: ti,
                outcome: "rejected".into(),
                representation: None,
                detail: None,
            };
            let mut chain = nodes[id].chain.clone();
            chain.push(ti);
            let chain_refs = chain
                .iter()
                .map(|&i| &request.transformations[i])
                .collect::<Vec<_>>();
            let form = match MapNormalForm::from_chain(&chain_refs, source.origin) {
                Ok(form) => form,
                Err(error) => {
                    edge.detail = Some(error);
                    edges.push(edge);
                    continue;
                }
            };
            if let Some(duplicate) = forms.iter().position(|f| f.certified_equal(&form)) {
                edge.outcome = "duplicate_map".into();
                edge.representation = Some(duplicate);
                edge.detail = Some(format!(
                    "certified equal map normal forms on {}: {}",
                    form.record().index_domain,
                    if form.record().identity {
                        "the identity".to_string()
                    } else {
                        form.record().segments.join("; then ")
                    }
                ));
                edges.push(edge);
                continue;
            }
            if nodes.len() >= request.limits.max_representations {
                search.event(
                    "max_representations",
                    Some(id),
                    format!("transformations[{ti}] was not applied"),
                );
                edge.outcome = "skipped_budget".into();
                edges.push(edge);
                continue;
            }
            match search.child(&nodes, id, ti, &form) {
                Ok(mut node) => {
                    node.equal_sequence_as = find_equal_sequence(&nodes, &node);
                    let new_id = nodes.len();
                    edge.outcome = "new_representation".into();
                    edge.representation = Some(new_id);
                    nodes.push(node);
                    forms.push(form);
                    queue.push_back(new_id);
                }
                Err(error) => edge.detail = Some(error),
            }
            edges.push(edge);
        }
    }
    let mut ranked = nodes
        .iter()
        .flat_map(|node| {
            node.descriptions
                .iter()
                .enumerate()
                .map(move |(i, d)| RankedDescription {
                    representation: node.id,
                    description: i,
                    label: d.label.clone(),
                    text: d.text.clone(),
                    cost: d.cost.clone(),
                })
        })
        .collect::<Vec<_>>();
    ranked.sort_by(|a, b| {
        (
            &a.cost,
            label_rank(&a.label),
            a.representation,
            a.description,
        )
            .cmp(&(
                &b.cost,
                label_rank(&b.label),
                b.representation,
                b.description,
            ))
    });
    Ok(RepresentationSearchReport {
        schema: RECURRENCE_REPRESENTATION_SEARCH_SCHEMA.to_string(),
        request: request.clone(),
        status: if search.budget.is_empty() {
            "complete"
        } else {
            "budget_exhausted"
        }
        .into(),
        source: RepresentationSource {
            kind: source.kind.to_string(),
            rule: source.rule.as_ref().map(|(_, _, json)| json.clone()),
            row_origin: source.origin,
            row_count: source.rows.len(),
        },
        representations: nodes,
        edges,
        ranked,
        depth_frontier,
        budget: search.budget,
        steps_used: search.steps_used,
        conventions: conventions(),
    })
}

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

fn check_description_on_rows(
    item: &RepresentationDescription,
    rows: &[Vec<Q>],
    origin: usize,
    label: &str,
) -> Result<(), String> {
    if let Some(json) = &item.scalar {
        let (rule, first, seeds) = json.to_recurrence_parts()?;
        if first < origin {
            return Err(format!("{label}: the rule starts before the rows"));
        }
        for (j, seed) in seeds.iter().enumerate() {
            if rows.get(first + j - origin).map(|r| vec_poly::trim(r)) != Some(vec_poly::trim(seed))
            {
                return Err(format!("{label}: seed {j} differs from the rows"));
            }
        }
        for target in first + seeds.len()..origin + rows.len() {
            if transformed_rule_holds_at(&rule, rows, origin, target) == Some(false) {
                return Err(format!("{label}: the rule fails at target {target}"));
            }
        }
    }
    if let Some(description) = &item.description {
        let assessment = assess_generating_description(description)?;
        if item.assessment.as_ref() != Some(&assessment) {
            return Err(format!("{label}: the assessment does not replay"));
        }
        if representation_label(&assessment.status) != item.label
            && !(item.label == CERTIFIED_RELATION_LABEL
                && assessment.status == EQUIVALENT_DESCRIPTION)
        {
            return Err(format!("{label}: the label does not match the assessment"));
        }
        for relation in &description.relations {
            let relation_rows = ImplicitRelation::from_json(&relation.relation)?;
            let end = (origin + rows.len()).saturating_sub(relation_rows.max_shift());
            for base in relation.first_base.max(origin)..end {
                if relation_rows.holds_at(rows, origin, base) == Some(false) {
                    return Err(format!("{label}: a relation fails at base {base}"));
                }
            }
        }
        for condition in description
            .exceptional_index_conditions
            .iter()
            .chain(&description.boundary_conditions)
        {
            let parsed = CoefficientCondition::from_json(&condition.condition)?;
            let end = (origin + rows.len()).saturating_sub(parsed.max_shift());
            for base in condition.first_base.max(origin)..end {
                if parsed.holds_at(rows, origin, base) == Some(false) {
                    return Err(format!("{label}: a condition fails at base {base}"));
                }
            }
        }
    }
    Ok(())
}

/// Replay a report: rerun the search from its embedded request and require
/// equality, replay every embedded step and complete-space certificate with
/// its own verifier, check that every step's defining rule is a certified
/// description of the parent, re-derive every representation from the
/// original source rows (map and exact inverse), and re-check every
/// description (seeds, rule, relations, conditions, assessment) on its rows.
pub fn verify_representation_search(report: &RepresentationSearchReport) -> Result<(), String> {
    if report.schema != RECURRENCE_REPRESENTATION_SEARCH_SCHEMA {
        return Err(format!("unsupported schema `{}`", report.schema));
    }
    let replay = search_representations(&report.request)?;
    if &replay != report {
        return Err("report does not match its exact replay".into());
    }
    let root = report
        .representations
        .first()
        .ok_or("the report has no source representation")?;
    let origin = report.source.row_origin;
    let source_rows = scalar_rows(root)?;
    for node in &report.representations {
        let label = format!("representations[{}]", node.id);
        let chain = node
            .chain
            .iter()
            .map(|&i| report.request.transformations.get(i))
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| format!("{label}: unknown transformation index"))?;
        let states = decode_states(&node.states)?;
        if apply_chain(&chain, origin, &source_rows)? != states {
            return Err(format!(
                "{label}: rows are not the image of the source rows"
            ));
        }
        if invert_chain(&chain, origin, &states)? != source_rows || !node.reconstruction_verified {
            return Err(format!("{label}: exact reconstruction of the source fails"));
        }
        for (s, step) in node.steps.iter().enumerate() {
            let parent = node
                .parent
                .and_then(|p| report.representations.get(p))
                .ok_or_else(|| format!("{label}: a step without a parent"))?;
            let source = step
                .source_description
                .map(|i| {
                    parent
                        .descriptions
                        .get(i)
                        .filter(|d| d.label != FINITE_FIT_LABEL)
                        .and_then(|d| d.scalar.clone())
                        .ok_or_else(|| format!("{label}.steps[{s}]: uncertified source"))
                })
                .transpose()?;
            if let Some(transform) = &step.transform {
                verify_recurrence_transform(transform)?;
                if transform.request.source_recurrence != source {
                    return Err(format!("{label}.steps[{s}]: source mismatch"));
                }
            }
            if let Some(implicit) = &step.implicit {
                verify_recurrence_implicit(implicit)?;
                if implicit.request.source_recurrence != source {
                    return Err(format!("{label}.steps[{s}]: source mismatch"));
                }
            }
        }
        for space in &node.spaces {
            verify_recurrence_spaces(&space.request, space)?;
        }
        if node.components.len() == 1 {
            let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
            for (i, item) in node.descriptions.iter().enumerate() {
                check_description_on_rows(
                    item,
                    &rows,
                    origin,
                    &format!("{label}.descriptions[{i}]"),
                )?;
            }
        }
    }
    Ok(())
}
