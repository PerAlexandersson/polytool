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
//!    certified parent relation the step maps rows only. A parent with
//!    certified implicit relations but no certified scalar rule feeds them to
//!    diagonal maps, derivative antecedents and the lossy `derivative` map
//!    through [`RecurrenceTransformation::transport_implicit`] (an
//!    [`ImplicitRelationStep`]), carrying its coefficient conditions where
//!    the map allows; the `derivative` map always takes this route;
//! 4. for a factorial diagonal map on a constant-index monic certified source,
//!    complete recurrence spaces are certified ([`find_recurrence_spaces`]
//!    with `coefficient_rescaling`) on a bounded menu of supports; after a
//!    Möbius-class map, a polynomial diagonal weight or `derivative`, the
//!    direct implicit search certifies them ([`search_implicit_relations`],
//!    transported-source and weight theorems; `implicit_spaces`, one
//!    successor support per report);
//! 5. the new coordinates are searched by the bounded fit menu
//!    ([`find_recurrence_representations_rational`]).
//!
//! # Vector representations
//!
//! The symmetric decomposition gives two-component states `(a_n, b_n)`.
//! They are expanded like scalar representations: every scalar bijective map
//! acts componentwise ([`RecurrenceTransformation::map_state`]; a schedule
//! acts on `b` at `d_n - 1`), the chain is inverted componentwise, and normal
//! forms keep one componentwise suffix per component. A certified vector
//! system of the parent is transported exactly through Möbius-class maps
//! ([`transport_vector_componentwise`], a [`VectorTransportStep`]); its
//! evidence is `transported_identity` when the parent system and every degree
//! condition are certified (degree bounds of the components are inherited
//! from a certified decomposition: `deg a_n <= d_n`, `deg b_n <= d_n - 1`).
//! Gamma coordinates of the components of a certified decomposition at the
//! same schedule have certified symmetry by construction. States are fitted
//! by a first-order vector system (`vector_fit`, finite evidence). A vector
//! system with its initial state carries a [`VectorSystemAssessment`]: it
//! determines the states when its shared left factor is a nonzero polynomial
//! at every index, certified by one `x`-coefficient without integer roots in
//! the domain; it is then an `equivalent_generating_description` when its
//! annihilation is certified and it starts at the row origin. Two vector
//! representations with identical equivalent systems and initial states are
//! the same sequence. Lossy and vector-producing maps are rejected on vector
//! states; the component count enters the ranking.
//!
//! # The rank-two gamma case
//!
//! When the even and odd parts of a gamma-transported rule are independent,
//! the step reports the pair (`transported.system`) and the representation
//! gets a `transported_gamma_system` description with both relations. They
//! are exact consequences of the source rule where the rows are symmetric,
//! so the description is certified only when that symmetry is: by its own
//! step, or by another step of the same edge whose source rule is invariant
//! under reversal (the symmetry is a property of the rows). Otherwise its
//! evidence is `checked_supplied_rows` and its label `finite_fit`.
//!
//! Every certified scalar description of the new representation is then
//! available to its children. A representation whose equivalent generating
//! description coincides with one of another representation (same relations,
//! seeds, first index and conditions) is the same sequence, certified by the
//! uniqueness arguments of both descriptions; it is kept (`equal_sequence_as`)
//! but not expanded.
//!
//! # Lossy maps
//!
//! A `derivative` map `Q_n = D^m P_n` forgets `[x^0..x^(m-1)] P_n`. The node
//! records these coefficients (`lost_data`), and the chain is inverted with
//! them. Its implicit-relation steps also record a [`LossyReconstruction`]:
//! the description of the parent through `Q` (pulled-back relations and
//! conditions) is assessed with the conditions `[x^k] (L P)_n = 0` that pin
//! the lost coefficients, and again without every condition that refers to
//! them; the first can be equivalent, the second is a consequence only.
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

use super::diagonal::{implicit_map_condition_statuses, weight_regularity};
use super::implicit::{biv_from_json, decode_poly, encode_poly, valuation_exceptional_indices};
use super::implicit_algebra::{binomial, roots_in_domain};
use super::implicit_search::finite_theorem_applies;
use super::lossy::coefficient_conditions_of;
use super::representation_map::{
    apply_chain, chain_lost_data, chain_offsets, invert_chain_with_lost_data, MapNormalForm,
};
use super::representation_report::{conventions, label_rank};
use super::*;
use crate::polynomial::Polynomial;
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
        vector_assessment: None,
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
            if mine.label == EQUIVALENT_DESCRIPTION_LABEL && mine.vector.is_some() {
                if let Some(j) = other.descriptions.iter().position(|theirs| {
                    theirs.label == EQUIVALENT_DESCRIPTION_LABEL && theirs.vector == mine.vector
                }) {
                    return Some(RepresentationEquality {
                        representation: other.id,
                        other_description: j,
                        description: i,
                        certificate: "both descriptions are equivalent vector systems with the \
                                      same system and initial state; each determines its states \
                                      uniquely"
                            .into(),
                    });
                }
            }
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

/// `(lag, derivative order)` pairs of the default complete-space menu,
/// simplest first.
fn menu_pairs(bounds: &TransformSearchBounds) -> Vec<(usize, usize)> {
    let mut pairs = (1..=bounds.max_rec_len.max(1))
        .flat_map(|l| (0..=bounds.max_diff_deg).map(move |j| (l, j)))
        .collect::<Vec<_>>();
    pairs.sort_by_key(|&(l, j)| (l + j, l));
    pairs
}

/// Maps whose complete spaces come from the direct implicit search.
fn implicit_space_map(t: &RecurrenceTransformation) -> bool {
    matches!(
        t,
        RecurrenceTransformation::MultiplyPower { .. }
            | RecurrenceTransformation::DividePower { .. }
            | RecurrenceTransformation::Affine { .. }
            | RecurrenceTransformation::Reversal { .. }
            | RecurrenceTransformation::Mobius { .. }
            | RecurrenceTransformation::FPolynomial { .. }
            | RecurrenceTransformation::Diagonal { .. }
            | RecurrenceTransformation::Derivative { .. }
    )
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
            lost_data: Vec::new(),
            map_conditions: Vec::new(),
            descriptions: merge_descriptions(descriptions),
            steps: Vec::new(),
            spaces: Vec::new(),
            implicit_spaces: Vec::new(),
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

    /// Certified descriptions of the parent with one implicit relation
    /// (scalar rules included, in implicit form), usable by the implicit
    /// machinery.
    fn implicit_source_descriptions(&self, parent: &RepresentationNode) -> Vec<usize> {
        parent
            .descriptions
            .iter()
            .enumerate()
            .filter(|(_, d)| {
                d.label != FINITE_FIT_LABEL
                    && d.description.as_ref().is_some_and(|g| {
                        g.first_index == self.origin
                            && g.relations.len() == 1
                            && evidence_is_certified(&g.relations[0].evidence)
                            && ImplicitRelation::from_json(&g.relations[0].relation).is_ok()
                    })
            })
            .map(|(i, _)| i)
            .take(self.request.limits.max_sources_per_edge)
            .collect()
    }

    fn implicit_relation_step(
        &mut self,
        t: &RecurrenceTransformation,
        parent: &RepresentationNode,
        source: usize,
        rows: &[Vec<Q>],
        items: &mut Vec<RepresentationDescription>,
        conditions: &mut Option<Vec<TransformConditionStatus>>,
    ) -> RepresentationStep {
        let mut step = RepresentationStep {
            source_description: Some(source),
            transform: None,
            implicit: None,
            implicit_relation: None,
            vector_transport: None,
            error: None,
        };
        let built = scalar_rows(parent).and_then(|parent_rows| {
            build_implicit_relation_step(
                t,
                &parent.descriptions[source],
                &parent_rows,
                rows,
                self.origin,
            )
        });
        match built {
            Ok((record, map_conditions)) => {
                if conditions.is_none() {
                    *conditions = Some(map_conditions);
                }
                let relation = &record.relation;
                items.push(RepresentationDescription {
                    label: representation_label(&record.assessment.status).into(),
                    kind: "implicit_transport".into(),
                    evidence: relation.evidence.clone(),
                    text: relation.text.clone(),
                    scalar: None,
                    vector: None,
                    description: Some(record.description.clone()),
                    assessment: Some(record.assessment.clone()),
                    cost: relation_cost(&relation.relation),
                    provenance: "steps[STEP].implicit_relation.description".into(),
                    also_found_by: Vec::new(),
                    note: relation.note.clone(),
                    vector_assessment: None,
                });
                step.implicit_relation = Some(record);
            }
            Err(error) => step.error = Some(error),
        }
        step
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
            implicit_relation: None,
            vector_transport: None,
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
            let assessment = assess_vector_system(vector).ok();
            items.push(RepresentationDescription {
                label: vector_label(
                    certified,
                    assessment.as_ref(),
                    vector.first_index == self.origin,
                    transported.valid_to_index.is_some(),
                )
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
                    "the complete component system with its initial state; uniqueness is \
                     vector_assessment"
                        .into(),
                ),
                vector_assessment: assessment,
            });
        } else if let Some(system) = &transported.system {
            if let Some(item) = gamma_system_item(report, transported, system, certified) {
                items.push(RepresentationDescription { provenance, ..item });
            }
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
            vector_assessment: None,
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
                    vector_assessment: None,
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
            let top = bounds.max_idx_deg.max(bounds.max_var_deg);
            let mut successes: Vec<(usize, usize)> = Vec::new();
            'pairs: for (l, j) in menu_pairs(&bounds) {
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

    /// Complete spaces after a map covered by a finite theorem of the direct
    /// implicit search (Möbius-class maps, polynomial diagonal weights,
    /// differentiation), one successor support per embedded report.
    fn implicit_spaces(
        &mut self,
        id: usize,
        t: &RecurrenceTransformation,
        parent: &RepresentationNode,
        node: &mut RepresentationNode,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let source = parent.descriptions.iter().find_map(|d| {
            let json = d.scalar.as_ref()?;
            (d.label != FINITE_FIT_LABEL
                && json.first_index == self.origin
                && finite_theorem_applies(t, json))
            .then(|| RecurrenceJson {
                search: None,
                ..json.clone()
            })
        });
        let Some(source) = source else {
            node.notes.push(format!(
                "complete spaces after {} need a certified constant-index monic pure-lag rule \
                 of the parent covered by a finite theorem; none is available",
                t.kind_name()
            ));
            return;
        };
        let mut tried = 0usize;
        let mut run = |search: &mut Self,
                       node: &mut RepresentationNode,
                       support: ImplicitSearchSupport|
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
            let request = ImplicitRelationSearchRequest {
                source_recurrence: Some(source.clone()),
                oeis_id: None,
                rows: None,
                row_origin: None,
                row_count: None,
                transformation: Some(t.clone()),
                domain_start: None,
                supports: vec![support],
                conditions: true,
                condition_index_degree: 1,
                limits: search.request.limits.space_limits(),
            };
            match search_implicit_relations(&request) {
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
                    node.implicit_spaces.push(report);
                    Some(node.implicit_spaces.len() - 1)
                }
                Err(error) => {
                    node.notes
                        .push(format!("complete-space support failed: {error}"));
                    None
                }
            }
        };
        let support = |l: usize, j: usize, u: usize, v: usize| ImplicitSearchSupport {
            role: SUCCESSOR_SUPPORT.into(),
            max_shift: l,
            derivative_order: j,
            index_degree: u,
            variable_degree: v,
        };
        if !self.request.space_shapes.is_empty() {
            for shape in self.request.space_shapes.clone() {
                match shape {
                    RecurrenceSpaceShape::Box {
                        length,
                        derivative_order,
                        index_degree,
                        variable_degree,
                        leading_index_degree,
                        leading_variable_degree,
                    } => {
                        run(
                            self,
                            node,
                            support(
                                length,
                                derivative_order,
                                index_degree.max(leading_index_degree),
                                variable_degree.max(leading_variable_degree),
                            ),
                        );
                    }
                    other => node.notes.push(format!(
                        "space shape {other:?} has no successor-support form; skipped"
                    )),
                }
            }
        } else {
            let bounds = self.request.search.clone();
            let top = bounds.max_idx_deg.max(bounds.max_var_deg);
            let mut successes: Vec<(usize, usize)> = Vec::new();
            'pairs: for (l, j) in menu_pairs(&bounds) {
                if successes.iter().any(|&(sl, sj)| sl <= l && sj <= j) {
                    continue;
                }
                for d in 0..=top {
                    let shape = support(l, j, d.min(bounds.max_idx_deg), d.min(bounds.max_var_deg));
                    let Some(index) = run(self, node, shape) else {
                        break 'pairs;
                    };
                    match node.implicit_spaces[index].results[0].status.as_str() {
                        "relations_found" => {
                            successes.push((l, j));
                            break;
                        }
                        "absent_in_support" => {}
                        _ => break,
                    }
                }
            }
        }
        let rows = scalar_rows(node).unwrap_or_default();
        for (s, report) in node.implicit_spaces.iter().enumerate() {
            let result = &report.results[0];
            if result.status != "relations_found"
                || result.evidence != CERTIFIED_COMPLETE_SPACE_EVIDENCE
            {
                continue;
            }
            for (r, rep) in result.representatives.iter().enumerate() {
                let provenance = format!("implicit_spaces[{s}].results[0].representatives[{r}]");
                let Ok(relation) = ImplicitRelation::from_json(&rep.relation.relation) else {
                    continue;
                };
                let label = representation_label(&rep.assessment.status).to_string();
                let mut item = RepresentationDescription {
                    label: label.clone(),
                    kind: "certified_space".into(),
                    evidence: CERTIFIED_COMPLETE_SPACE_EVIDENCE.into(),
                    text: rep.relation.text.clone(),
                    scalar: None,
                    vector: None,
                    description: Some(rep.description.clone()),
                    assessment: Some(rep.assessment.clone()),
                    cost: relation_cost(&rep.relation.relation),
                    provenance: provenance.clone(),
                    also_found_by: Vec::new(),
                    note: rep.relation.note.clone(),
                    vector_assessment: None,
                };
                if let Some((rule, shift)) = relation.to_successor_rule() {
                    if rep.relation.first_base == self.origin && rows.len() >= shift {
                        let json = RecurrenceJson::from_recurrence_rational(
                            &rule,
                            self.origin,
                            &rows[..shift],
                            None,
                        );
                        let scalar = scalar_description(ScalarInput {
                            kind: "certified_space",
                            rule: &rule,
                            json,
                            evidence: CERTIFIED_COMPLETE_SPACE_EVIDENCE,
                            valid_from: self.origin + shift,
                            bounded_domain: false,
                            symbol: "Q",
                            provenance,
                        });
                        // Keep the representative's description when its
                        // conditions give a stronger assessment.
                        if label_rank(&scalar.label) <= label_rank(&label) {
                            item = scalar;
                        } else {
                            item.scalar = scalar.scalar;
                            item.text = scalar.text;
                            item.cost = scalar.cost;
                        }
                    }
                }
                items.push(item);
            }
        }
    }

    /// Map conditions and transported vector systems of an edge from a
    /// vector parent (a componentwise map).
    #[allow(clippy::too_many_arguments)]
    fn vector_steps(
        &mut self,
        nodes: &[RepresentationNode],
        t: &RecurrenceTransformation,
        parent: &RepresentationNode,
        offsets: &[i64],
        states: &[Vec<Vec<Q>>],
        node: &mut RepresentationNode,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let id = node.id;
        let transformations = &self.request.transformations;
        let bounds = component_bounds(transformations, nodes, parent.id);
        let symmetric = symmetric_decomposition_schedule(transformations, nodes, parent.id);
        node.map_conditions = vector_map_conditions(
            t,
            bounds.as_deref(),
            symmetric.as_ref(),
            offsets,
            &parent.components,
            self.origin,
        );
        if !t.is_mobius_class() {
            node.notes.push(format!(
                "vector systems are transported through Möbius-class maps only; {} maps the \
                 rows componentwise",
                t.kind_name()
            ));
            return;
        }
        let sources = parent
            .descriptions
            .iter()
            .enumerate()
            .filter(|(_, d)| d.label != FINITE_FIT_LABEL && d.vector.is_some())
            .map(|(i, _)| i)
            .take(self.request.limits.max_sources_per_edge)
            .collect::<Vec<_>>();
        if sources.is_empty() {
            node.notes.push(
                "no certified vector system of the parent to transport; rows are mapped \
                 componentwise"
                    .into(),
            );
            return;
        }
        let bounded = offsets.iter().any(|&o| {
            t.with_degree_offset(o)
                .ok()
                .and_then(|m| m.schedule().and_then(|s| s.last_index))
                .is_some()
        });
        for source in sources {
            if !self.take_step(id, "vector transport step") {
                break;
            }
            let mut step = RepresentationStep {
                source_description: Some(source),
                transform: None,
                implicit: None,
                implicit_relation: None,
                vector_transport: None,
                error: None,
            };
            let built = parent.descriptions[source]
                .vector
                .as_ref()
                .ok_or_else(|| "the source has no vector system".to_string())
                .and_then(|json| {
                    build_vector_transport(
                        t,
                        offsets,
                        json,
                        states,
                        self.origin,
                        &node.map_conditions,
                        &node.components,
                    )
                });
            match built {
                Ok(record) => {
                    let assessment = assess_vector_system(&record.transported).ok();
                    items.push(RepresentationDescription {
                        label: vector_label(
                            evidence_is_certified(&record.evidence),
                            assessment.as_ref(),
                            record.transported.first_index == self.origin,
                            bounded,
                        )
                        .into(),
                        kind: "transported_vector_system".into(),
                        evidence: record.evidence.clone(),
                        text: record.text.clone(),
                        scalar: None,
                        vector: Some(record.transported.clone()),
                        description: None,
                        assessment: None,
                        cost: vector_cost(&record.transported),
                        provenance: format!(
                            "steps[{}].vector_transport.transported",
                            node.steps.len()
                        ),
                        also_found_by: Vec::new(),
                        note: Some(
                            "the parent's vector system transported componentwise; uniqueness \
                             is vector_assessment"
                                .into(),
                        ),
                        vector_assessment: assessment,
                    });
                    step.vector_transport = Some(record);
                }
                Err(error) => step.error = Some(error),
            }
            node.steps.push(step);
        }
    }

    /// A first-order vector system fitted to the states (finite evidence).
    fn vector_fits(
        &mut self,
        id: usize,
        states: &[Vec<Vec<Q>>],
        names: &[String],
    ) -> (
        Option<RepresentationFitSummary>,
        Vec<RepresentationDescription>,
    ) {
        if !self.request.fits {
            return (None, Vec::new());
        }
        let summary =
            |status: &str, found: usize, reason: Option<String>| RepresentationFitSummary {
                status: status.into(),
                configurations_considered: usize::from(status == "complete"),
                fits_found: found,
                reason,
            };
        if states.len() < 3 {
            return (Some(summary("insufficient_rows", 0, None)), Vec::new());
        }
        let request = self.request;
        let bounds = &request.search;
        let components = states[0].len();
        let unknowns = [
            components * components,
            bounds.max_idx_deg + 1,
            bounds.max_var_deg + 1,
            bounds.max_diff_deg + 1,
        ]
        .iter()
        .try_fold(1usize, |acc, &f| acc.checked_mul(f))
        .and_then(|matrix| {
            matrix.checked_add(
                (bounds.vector_denominator_idx_deg + 1)
                    .checked_mul(bounds.vector_denominator_var_deg + 1)?,
            )
        });
        if unknowns.is_none_or(|u| u > MAX_VECTOR_FIT_UNKNOWNS) {
            self.event(
                "vector_fit_unknowns",
                Some(id),
                format!("the vector fit has more than {MAX_VECTOR_FIT_UNKNOWNS} unknowns; not run"),
            );
            return (Some(summary("skipped_budget", 0, None)), Vec::new());
        }
        if !self.take_step(id, "vector fit") {
            return (Some(summary("skipped_budget", 0, None)), Vec::new());
        }
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
        let fit =
            match find_vector_recurrence_with_denominator_rational(states, self.origin, &options) {
                Ok(fit) => fit,
                Err(error) => {
                    return (
                        Some(summary("complete", 0, Some(error.to_string()))),
                        Vec::new(),
                    )
                }
            };
        let holds = states.windows(2).enumerate().all(|(i, pair)| {
            fit.recurrence
                .holds_transition_rational(&pair[0], &pair[1], self.origin + i)
        });
        let json = TransformVectorRecurrenceJson::from_recurrence(
            &fit.recurrence,
            self.origin,
            &states[0],
        );
        let text = json.text(names);
        let (true, Ok(text)) = (holds, text) else {
            return (
                Some(summary(
                    "complete",
                    0,
                    Some("the fit fails on a transition".into()),
                )),
                Vec::new(),
            );
        };
        let assessment = assess_vector_system(&json).ok();
        let item = RepresentationDescription {
            label: FINITE_FIT_LABEL.into(),
            kind: "vector_fit".into(),
            evidence: FINITE_FIT_EVIDENCE.into(),
            text,
            scalar: None,
            cost: vector_cost(&json),
            vector: Some(json),
            description: None,
            assessment: None,
            provenance: "fits".into(),
            also_found_by: Vec::new(),
            note: Some(format!(
                "a first-order system fitted to {} transitions ({} held out); finite evidence",
                fit.fit_transitions, fit.held_out_transitions
            )),
            vector_assessment: assessment,
        };
        (Some(summary("complete", 1, None)), vec![item])
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
        let parent_chain = parent
            .chain
            .iter()
            .map(|&i| &request.transformations[i])
            .collect::<Vec<_>>();
        let offsets = chain_offsets(&parent_chain)
            .pop()
            .expect("at least one offset list");
        let vector_parent = parent.components.len() > 1;
        let states = decode_states(&parent.states)?
            .iter()
            .enumerate()
            .map(|(i, state)| t.map_state(origin + i, state, &offsets))
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
        let lost_data = chain_lost_data(&chain_refs, origin, &self.source_rows)?;
        let reconstruction_verified =
            invert_chain_with_lost_data(&chain_refs, origin, &states, &lost_data)
                .is_ok_and(|back| back == self.source_rows);
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
            components: if vector_parent {
                parent.components.clone()
            } else {
                component_names(Some(t), 1)
            },
            row_origin: origin,
            states: encode_states(&states),
            reconstruction_verified,
            lost_data,
            map_conditions: Vec::new(),
            descriptions: Vec::new(),
            steps: Vec::new(),
            spaces: Vec::new(),
            implicit_spaces: Vec::new(),
            fits: None,
            equal_sequence_as: None,
            expanded: false,
            notes: Vec::new(),
        };
        let mut items = Vec::new();
        if vector_parent {
            self.vector_steps(nodes, t, parent, &offsets, &states, &mut node, &mut items);
        } else {
            self.scalar_steps(t, parent, &states, &mut node, &mut items);
        }
        if node.components.len() == 1 {
            let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
            let (fits, fitted) = self.fits(id, &rows, "Q");
            node.fits = fits;
            items.extend(fitted);
        } else {
            let (fits, fitted) = self.vector_fits(id, &states, &node.components);
            node.fits = fits;
            items.extend(fitted);
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

    /// Transport steps, complete spaces and map conditions of an edge from a
    /// scalar parent.
    fn scalar_steps(
        &mut self,
        t: &RecurrenceTransformation,
        parent: &RepresentationNode,
        states: &[Vec<Vec<Q>>],
        node: &mut RepresentationNode,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let id = node.id;
        let request = self.request;
        let mut conditions = None;
        let sources = self.source_descriptions(parent);
        let implicit_sources = self.implicit_source_descriptions(parent);
        let implicit_map = matches!(
            t,
            RecurrenceTransformation::Diagonal { .. }
                | RecurrenceTransformation::DerivativeAntecedent { .. }
                | RecurrenceTransformation::Derivative { .. }
        );
        if t.is_lossy() || (implicit_map && sources.is_empty() && !implicit_sources.is_empty()) {
            if implicit_sources.is_empty() {
                node.notes.push(
                    "no certified relation of the parent to transport through the lossy map; \
                     rows are mapped only"
                        .into(),
                );
            }
            let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
            for source in implicit_sources {
                if !self.take_step(id, "implicit transport step") {
                    break;
                }
                let before = items.len();
                let step =
                    self.implicit_relation_step(t, parent, source, &rows, items, &mut conditions);
                let index = node.steps.len();
                for item in &mut items[before..] {
                    item.provenance = item.provenance.replace("STEP", &index.to_string());
                }
                node.steps.push(step);
            }
            if t.is_lossy() && !node.steps.is_empty() {
                node.notes.push(format!(
                    "lossy map: the source is rebuilt only with lost_data; steps[k].implicit_relation.lossy \
                     assesses the description of the parent with and without it ({})",
                    t.kind_name()
                ));
            }
        } else {
            if sources.is_empty() && !implicit_sources.is_empty() {
                node.notes.push(format!(
                    "the parent has certified implicit relations only; they are not transported \
                     through {}",
                    t.kind_name()
                ));
            }
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
                let step = self.run_step(id, t, parent, source, states, items, &mut conditions);
                let index = node.steps.len();
                for item in &mut items[before..] {
                    item.provenance = item.provenance.replace("STEP", &index.to_string());
                }
                node.steps.push(step);
            }
        }
        node.map_conditions = conditions.unwrap_or_default();
        if request.certified_spaces && node.components.len() == 1 {
            if let RecurrenceTransformation::Diagonal {
                weight: DiagonalWeight::FactorialPower { power },
            } = t
            {
                self.spaces(id, *power, parent, node, items);
            } else if implicit_space_map(t) {
                self.implicit_spaces(id, t, parent, node, items);
            }
        }
        if matches!(t, RecurrenceTransformation::Gamma { .. }) {
            upgrade_gamma_systems(&node.steps, items);
        }
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
            let t = &request.transformations[ti];
            if nodes[id].components.len() > 1 && (t.component_count() > 1 || t.is_lossy()) {
                edge.detail = Some(format!(
                    "{} does not act componentwise on a vector representation",
                    t.kind_name()
                ));
                edges.push(edge);
                continue;
            }
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
// Vector representations and the two-relation gamma system
// ---------------------------------------------------------------------------

/// Uniqueness of a first-order vector system with its initial state: the
/// states are determined when `q(n, .)` is a nonzero polynomial at every
/// index `n >= first_index`, certified by a coefficient `[x^j] q(n,x)` (a
/// polynomial in `n`) with no integer root `n >= first_index`.
fn assess_vector_system(
    json: &TransformVectorRecurrenceJson,
) -> Result<VectorSystemAssessment, String> {
    let (system, first, _) = json.to_recurrence()?;
    let width = system
        .denominator
        .coeffs
        .iter()
        .map(Vec::len)
        .max()
        .unwrap_or(0);
    for j in 0..width {
        let column = vec_poly::trim(
            &system
                .denominator
                .coeffs
                .iter()
                .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
                .collect::<Vec<_>>(),
        );
        if vec_poly::is_zero(&column) {
            continue;
        }
        if roots_in_domain(&column, first, None).is_empty() {
            return Ok(VectorSystemAssessment {
                status: "certified".into(),
                denominator_coefficient: Some(j),
                statement: format!(
                    "[x^{j}] q(n,x) has no integer root n >= {first}, so q(n, .) is a nonzero \
                     polynomial and the system determines every state F_(n+1) from F_n; with \
                     the initial state at index {first} every later state is determined"
                ),
            });
        }
    }
    Ok(VectorSystemAssessment {
        status: "not_certified".into(),
        denominator_coefficient: None,
        statement: format!(
            "every coefficient of q(n,x) in x vanishes at some integer index n >= {first} (or \
             q = 0): uniqueness is not certified"
        ),
    })
}

/// The search label of a vector system: equivalent when its annihilation is
/// certified, its uniqueness is certified, it starts at the row origin and
/// the domain is unbounded; certified when annihilation is; else a fit.
fn vector_label(
    certified: bool,
    assessment: Option<&VectorSystemAssessment>,
    from_origin: bool,
    bounded: bool,
) -> &'static str {
    if !certified {
        FINITE_FIT_LABEL
    } else if assessment.is_some_and(|a| a.status == "certified") && from_origin && !bounded {
        EQUIVALENT_DESCRIPTION_LABEL
    } else {
        CERTIFIED_RELATION_LABEL
    }
}

/// The description of a rank-two gamma transport: seeds, the even part and
/// the odd companion, with the step's evidence.
fn gamma_system_item(
    report: &RecurrenceTransformReport,
    transported: &TransportedRuleReport,
    system: &GammaRelationSystemJson,
    certified: bool,
) -> Option<RepresentationDescription> {
    let origin = report.row_origin;
    let first = report.source_rule.as_ref()?.first_index;
    let seeds = report
        .transformed_states
        .get(first.checked_sub(origin)?..transported.valid_from_index.checked_sub(origin)?)?
        .iter()
        .map(|state| state.first().cloned())
        .collect::<Option<Vec<_>>>()?;
    let evidence = if certified {
        TRANSPORTED_IDENTITY_EVIDENCE
    } else {
        CHECKED_SUPPLIED_ROWS
    };
    let description = GeneratingDescription {
        first_index: first,
        seeds,
        relations: system
            .relations
            .iter()
            .map(|r| DescribedRelation {
                evidence: evidence.into(),
                ..r.clone()
            })
            .collect(),
        exceptional_index_conditions: Vec::new(),
        boundary_conditions: Vec::new(),
    };
    let mut cost = DescriptionCost {
        shift_length: 0,
        derivative_order: 0,
        coefficient_degree: 0,
        index_degree: 0,
        component_count: system.relations.len(),
        reconstruction_cost: 0,
    };
    for relation in &system.relations {
        let c = relation_cost(&relation.relation);
        cost.shift_length = cost.shift_length.max(c.shift_length);
        cost.derivative_order = cost.derivative_order.max(c.derivative_order);
        cost.coefficient_degree = cost.coefficient_degree.max(c.coefficient_degree);
        cost.index_degree = cost.index_degree.max(c.index_degree);
    }
    let mut item = RepresentationDescription {
        label: FINITE_FIT_LABEL.into(),
        kind: "transported_gamma_system".into(),
        evidence: evidence.into(),
        text: transported.text.clone(),
        scalar: None,
        vector: None,
        description: Some(description),
        assessment: None,
        cost,
        provenance: String::new(),
        also_found_by: Vec::new(),
        note: Some(
            "the even part and the odd companion of a rank-two gamma transport: together they \
             carry the source rule where the rows are symmetric (no single scalar rule does)"
                .into(),
        ),
        vector_assessment: None,
    };
    assess_gamma_item(&mut item, transported.valid_to_index.is_some());
    Some(item)
}

fn assess_gamma_item(item: &mut RepresentationDescription, bounded: bool) {
    let Some(description) = &item.description else {
        return;
    };
    match assess_generating_description(description) {
        Ok(assessment) => {
            let mut label = representation_label(&assessment.status);
            if bounded && label == EQUIVALENT_DESCRIPTION_LABEL {
                label = CERTIFIED_RELATION_LABEL;
            }
            item.label = label.into();
            item.assessment = Some(assessment);
        }
        Err(_) => {
            item.label = if evidence_is_certified(&item.evidence) {
                CERTIFIED_RELATION_LABEL
            } else {
                FINITE_FIT_LABEL
            }
            .into();
            item.assessment = None;
        }
    }
}

/// The index `k` of a provenance `steps[k]...`.
fn step_index_of(provenance: &str) -> Option<usize> {
    provenance
        .strip_prefix("steps[")?
        .split(']')
        .next()?
        .parse()
        .ok()
}

/// The step of a gamma edge whose report certifies the symmetry condition
/// (its source rule is invariant under reversal), if any.
fn symmetry_certifying_step(steps: &[RepresentationStep]) -> Option<usize> {
    steps.iter().position(|s| {
        s.source_description.is_some()
            && s.transform.as_ref().is_some_and(|r| {
                !r.conditions.is_empty()
                    && r.conditions
                        .iter()
                        .all(|c| c.status == CERTIFIED_ALL_INDICES)
            })
    })
}

/// A two-relation gamma system rests on its certified source rule and on the
/// symmetry of the rows. When another step of the same edge certifies that
/// symmetry, the system's relations are exact consequences for every index.
fn upgrade_gamma_systems(steps: &[RepresentationStep], items: &mut [RepresentationDescription]) {
    let Some(k) = symmetry_certifying_step(steps) else {
        return;
    };
    for item in items
        .iter_mut()
        .filter(|i| i.kind == "transported_gamma_system" && i.evidence == CHECKED_SUPPLIED_ROWS)
    {
        let Some(step) = step_index_of(&item.provenance).and_then(|s| steps.get(s)) else {
            continue;
        };
        let Some(report) = &step.transform else {
            continue;
        };
        let Some(transported) = report
            .transported
            .as_ref()
            .filter(|t| t.holds_on_transformed_rows && t.system.is_some())
        else {
            continue;
        };
        if step.source_description.is_none() || !report.reconstruction_verified {
            continue;
        }
        let Some(description) = item.description.as_mut() else {
            continue;
        };
        for relation in &mut description.relations {
            relation.evidence = TRANSPORTED_IDENTITY_EVIDENCE.into();
        }
        item.evidence = TRANSPORTED_IDENTITY_EVIDENCE.into();
        item.note = Some(format!(
            "the even part and the odd companion of a rank-two gamma transport of a certified \
             rule; the symmetry of the rows is certified by steps[{k}], so both relations hold \
             for every index"
        ));
        assess_gamma_item(item, transported.valid_to_index.is_some());
    }
}

/// `first_degree + step (n - first_index)`.
fn schedule_value(schedule: &DegreeSchedule, n: usize) -> i128 {
    schedule.first_degree as i128
        + schedule.step as i128 * (n as i128 - schedule.first_index as i128)
}

/// Whether `bound(n) <= schedule(n)` at every index `n >= from` of both
/// domains, with `bound` known from `from` on.
fn bound_at_most(bound: &DegreeSchedule, schedule: &DegreeSchedule, from: usize) -> bool {
    if bound.first_index > from {
        return false;
    }
    let start = from.max(schedule.first_index);
    let end = match (bound.last_index, schedule.last_index) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (a, b) => a.or(b),
    };
    if end.is_some_and(|e| e < start) {
        return true;
    }
    let difference = |n: usize| schedule_value(schedule, n) - schedule_value(bound, n);
    difference(start) >= 0
        && match end {
            Some(e) => difference(e) >= 0,
            None => schedule.step >= bound.step,
        }
}

fn shifted_schedule(schedule: &DegreeSchedule, offset: i64) -> Option<DegreeSchedule> {
    RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    }
    .with_degree_offset(offset)
    .ok()?
    .schedule()
    .cloned()
}

/// Certified degree bounds `deg F_(n,c) <= B_c(n)` of the components of a
/// vector representation: from a symmetric decomposition whose degree
/// condition is certified (`a` at `d_n`, `b` at `d_n - 1`), carried through
/// componentwise maps whose image degree is controlled.
fn component_bounds(
    transformations: &[RecurrenceTransformation],
    nodes: &[RepresentationNode],
    id: usize,
) -> Option<Vec<DegreeSchedule>> {
    let node = nodes.get(id)?;
    let t = transformations.get(node.transformation?)?;
    let parent = nodes.get(node.parent?)?;
    let certified = !node.map_conditions.is_empty()
        && node
            .map_conditions
            .iter()
            .all(|c| c.status == CERTIFIED_ALL_INDICES);
    if parent.components.len() == 1 {
        let RecurrenceTransformation::SymmetricDecomposition { schedule } = t else {
            return None;
        };
        if !certified {
            return None;
        }
        return t
            .component_degree_offsets()
            .iter()
            .map(|&o| shifted_schedule(schedule, o))
            .collect();
    }
    let bounds = component_bounds(transformations, nodes, parent.id)?;
    match t {
        RecurrenceTransformation::Reversal { schedule }
        | RecurrenceTransformation::Mobius { schedule, .. }
        | RecurrenceTransformation::FPolynomial { schedule } => {
            if !certified {
                return None;
            }
            let chain = node
                .chain
                .iter()
                .map(|&i| transformations.get(i))
                .collect::<Option<Vec<_>>>()?;
            let offsets = chain_offsets(&chain).pop()?;
            offsets
                .iter()
                .map(|&o| shifted_schedule(schedule, o))
                .collect()
        }
        RecurrenceTransformation::MultiplyPower { power } => bounds
            .iter()
            .map(|b| {
                Some(DegreeSchedule {
                    first_degree: b.first_degree.checked_add(*power)?,
                    ..b.clone()
                })
            })
            .collect(),
        RecurrenceTransformation::Affine { .. }
        | RecurrenceTransformation::Diagonal { .. }
        | RecurrenceTransformation::DerivativeAntecedent { .. } => Some(bounds),
        _ => None,
    }
}

/// The schedule of a symmetric decomposition of a scalar representation
/// whose degree condition is certified: its components are then symmetric
/// at `d_n` and `d_n - 1` by construction.
fn symmetric_decomposition_schedule(
    transformations: &[RecurrenceTransformation],
    nodes: &[RepresentationNode],
    id: usize,
) -> Option<DegreeSchedule> {
    let node = nodes.get(id)?;
    let parent = nodes.get(node.parent?)?;
    let RecurrenceTransformation::SymmetricDecomposition { schedule } =
        transformations.get(node.transformation?)?
    else {
        return None;
    };
    (parent.components.len() == 1
        && !node.map_conditions.is_empty()
        && node
            .map_conditions
            .iter()
            .all(|c| c.status == CERTIFIED_ALL_INDICES))
    .then(|| schedule.clone())
}

/// Preconditions of a componentwise map on a vector representation, one per
/// component where they depend on it.
fn vector_map_conditions(
    t: &RecurrenceTransformation,
    bounds: Option<&[DegreeSchedule]>,
    symmetric: Option<&DegreeSchedule>,
    offsets: &[i64],
    names: &[String],
    origin: usize,
) -> Vec<TransformConditionStatus> {
    let finite = "exact on every supplied row; finite evidence only";
    let checked = |condition: String, detail: &str| TransformConditionStatus {
        condition,
        status: CHECKED_SUPPLIED_ROWS.into(),
        detail: detail.into(),
    };
    let offset_text = |o: i64| match o {
        0 => String::new(),
        o if o < 0 => format!(" - {}", -o),
        o => format!(" + {o}"),
    };
    match t {
        RecurrenceTransformation::Reversal { schedule }
        | RecurrenceTransformation::Mobius { schedule, .. }
        | RecurrenceTransformation::FPolynomial { schedule } => names
            .iter()
            .zip(offsets)
            .enumerate()
            .map(|(c, (name, &o))| {
                let condition = format!("deg {name}_n <= d_n{}", offset_text(o));
                let shifted = shifted_schedule(schedule, o);
                match (bounds.and_then(|b| b.get(c)), shifted) {
                    (Some(bound), Some(shifted)) if bound_at_most(bound, &shifted, origin) => {
                        TransformConditionStatus {
                            condition,
                            status: CERTIFIED_ALL_INDICES.into(),
                            detail: format!(
                                "the parent's certified bound deg {name}_n <= {} + {}*(n - {}) \
                                 is at most d_n{} for every index of the domain",
                                bound.first_degree,
                                bound.step,
                                bound.first_index,
                                offset_text(o)
                            ),
                        }
                    }
                    _ => checked(condition, finite),
                }
            })
            .collect(),
        RecurrenceTransformation::DividePower { power } => names
            .iter()
            .map(|name| checked(format!("x^{power} divides {name}_n"), finite))
            .collect(),
        RecurrenceTransformation::Gamma { schedule } => names
            .iter()
            .zip(offsets)
            .map(|(name, &o)| {
                let condition = format!(
                    "x^(d_n{0}) {name}_n(1/x) = {name}_n (symmetry at d_n{0})",
                    offset_text(o)
                );
                if symmetric == Some(schedule) {
                    TransformConditionStatus {
                        condition,
                        status: CERTIFIED_ALL_INDICES.into(),
                        detail: "the components of a symmetric decomposition at the same \
                                 schedule are symmetric at d_n and d_n - 1 by construction, \
                                 and its degree condition is certified"
                            .into(),
                    }
                } else {
                    checked(condition, finite)
                }
            })
            .collect(),
        RecurrenceTransformation::Diagonal { .. }
        | RecurrenceTransformation::DerivativeAntecedent { .. } => {
            implicit_map_condition_statuses(t, origin)
        }
        _ => Vec::new(),
    }
}

/// Transport the parent's vector system through the componentwise map and
/// check it on the new states.
fn build_vector_transport(
    t: &RecurrenceTransformation,
    offsets: &[i64],
    source: &TransformVectorRecurrenceJson,
    states: &[Vec<Vec<Q>>],
    origin: usize,
    conditions: &[TransformConditionStatus],
    names: &[String],
) -> Result<VectorTransportStep, String> {
    let (system, first, initial) = source.to_recurrence()?;
    let start = first
        .checked_sub(origin)
        .filter(|i| *i < states.len())
        .ok_or("the source system starts outside the rows")?;
    if t.map_state(first, &initial, offsets)
        .map_err(|e| e.to_string())?
        != states[start]
    {
        return Err("the mapped initial state differs from the rows".into());
    }
    let transported =
        transport_vector_componentwise(&system, t, offsets).map_err(|e| e.to_string())?;
    let mut holds = true;
    let mut checked = 0;
    for i in start..states.len().saturating_sub(1) {
        holds &= transported.holds_transition_rational(&states[i], &states[i + 1], origin + i);
        checked += 1;
    }
    let component_maps = offsets
        .iter()
        .map(|&o| t.with_degree_offset(o))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| e.to_string())?;
    let json = TransformVectorRecurrenceJson::from_recurrence(&transported, first, &states[start]);
    let certified = holds && conditions.iter().all(|c| c.status == CERTIFIED_ALL_INDICES);
    Ok(VectorTransportStep {
        component_maps,
        text: json.text(names)?,
        transported: json,
        evidence: if certified {
            TRANSPORTED_IDENTITY_EVIDENCE
        } else {
            CHECKED_SUPPLIED_ROWS
        }
        .into(),
        holds_on_rows: holds,
        checked_transitions: checked,
    })
}

/// Re-check a vector description on the states: initial state, every
/// transition from its first index, uniqueness and label.
fn check_vector_description(
    item: &RepresentationDescription,
    states: &[Vec<Vec<Q>>],
    origin: usize,
    label: &str,
) -> Result<(), String> {
    let Some(json) = &item.vector else {
        return Ok(());
    };
    let (system, first, initial) = json.to_recurrence()?;
    let start = first
        .checked_sub(origin)
        .filter(|i| *i < states.len())
        .ok_or_else(|| format!("{label}: the system starts outside the rows"))?;
    if states[start] != initial {
        return Err(format!("{label}: the initial state differs from the rows"));
    }
    for i in start..states.len().saturating_sub(1) {
        if !system.holds_transition_rational(&states[i], &states[i + 1], origin + i) {
            return Err(format!("{label}: the system fails at n = {}", origin + i));
        }
    }
    let assessment = assess_vector_system(json)?;
    if item.vector_assessment.as_ref() != Some(&assessment) {
        return Err(format!("{label}: the vector assessment does not replay"));
    }
    if item.label == EQUIVALENT_DESCRIPTION_LABEL
        && (assessment.status != "certified"
            || first != origin
            || !evidence_is_certified(&item.evidence))
    {
        return Err(format!("{label}: the label does not match the assessment"));
    }
    if item.label != FINITE_FIT_LABEL && !evidence_is_certified(&item.evidence) {
        return Err(format!(
            "{label}: a certified label with uncertified evidence"
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Implicit-relation steps
// ---------------------------------------------------------------------------

/// `c / lambda(n + shift, k)` for a diagonal weight, as `(numerator(n),
/// denominator(n))` with `lambda = K (k!)^p prod f_j(n,k)^(e_j)`.
fn inverse_weight(weight: &DiagonalWeight, shift: usize, k: usize) -> Option<(Vec<Q>, Vec<Q>)> {
    let kf = super::implicit_algebra::factorial(k);
    let k_value = Q::from_integer(BigInt::from(k));
    let constant = |c: Q| Some((vec![c], vec![Q::from_integer(1.into())]));
    match weight {
        DiagonalWeight::FactorialPower { power } => {
            let mut value = Q::from_integer(1.into());
            for _ in 0..*power {
                value *= &kf;
            }
            constant(value)
        }
        DiagonalWeight::Binomial { m } => constant(binomial(k + m, *m).recip()),
        DiagonalWeight::InverseBinomial { m } => constant(binomial(k + m, *m)),
        DiagonalWeight::Product {
            scale,
            factorial_exponent,
            factors,
        } => {
            let mut value = parse_rational_coeff(scale).ok()?.recip();
            for _ in 0..factorial_exponent.unsigned_abs() {
                if *factorial_exponent > 0 {
                    value /= &kf;
                } else {
                    value *= &kf;
                }
            }
            let mut numerator = vec![value];
            let mut denominator = vec![Q::from_integer(1.into())];
            for factor in factors {
                let f = biv_from_json(&factor.coeffs).ok()?;
                let at = f
                    .shift(&Q::from_integer(BigInt::from(shift)), &Q::zero())
                    .eval_y(&k_value);
                let power = vec_poly::pow(&at, factor.exponent.unsigned_abs() as usize);
                if factor.exponent > 0 {
                    denominator = vec_poly::mul(&denominator, &power);
                } else {
                    numerator = vec_poly::mul(&numerator, &power);
                }
            }
            Some((numerator, denominator))
        }
    }
}

/// A coefficient condition of the parent rows `P` rewritten for `Q = T(P)`,
/// or `None` when the map forgets a coefficient it refers to (or is not an
/// implicit map).
fn transport_condition(
    t: &RecurrenceTransformation,
    condition: &CoefficientCondition,
) -> Option<CoefficientCondition> {
    let condition = condition.canonical();
    let transported = match t {
        RecurrenceTransformation::Derivative { order } => {
            push_condition_through_derivative(&condition, *order)?
        }
        RecurrenceTransformation::DerivativeAntecedent { g } => {
            // [x^e] P = binom(e+m,m) sum_j g_j [x^(e+m-j)] Q.
            let (g, m) = derivative_antecedent_parameters(g).ok()?;
            let mut terms = Vec::new();
            for term in &condition.terms {
                let e = term.coefficient_index;
                for (j, gj) in g.iter().enumerate() {
                    if gj.is_zero() || e + m < j {
                        continue;
                    }
                    terms.push(ConditionTerm {
                        shift: term.shift,
                        coefficient_index: e + m - j,
                        coeff: vec_poly::scale(&term.coeff, &(binomial(e + m, m) * gj)),
                    });
                }
            }
            CoefficientCondition { terms }.canonical()
        }
        RecurrenceTransformation::Diagonal { weight } => {
            // [x^k] P_(n+i) = [x^k] Q_(n+i) / lambda(n+i, k); clear the
            // polynomial denominators by their lcm (nonzero on the domain by
            // the weight hypothesis).
            let mut raw = Vec::new();
            for term in &condition.terms {
                let (numerator, denominator) =
                    inverse_weight(weight, term.shift, term.coefficient_index)?;
                raw.push((term, vec_poly::mul(&term.coeff, &numerator), denominator));
            }
            let lcm = raw.iter().fold(
                Polynomial::new(vec![Q::from_integer(1.into())]),
                |acc, (_, _, d)| {
                    let d = Polynomial::new(d.clone());
                    let g = acc.gcd(&d);
                    let (quotient, _) = (acc.clone() * d).div_rem(&g);
                    quotient
                },
            );
            let terms = raw
                .into_iter()
                .map(|(term, numerator, denominator)| {
                    let (factor, _) = lcm.div_rem(&Polynomial::new(denominator));
                    ConditionTerm {
                        shift: term.shift,
                        coefficient_index: term.coefficient_index,
                        coeff: vec_poly::trim(&vec_poly::mul(&numerator, factor.coeffs())),
                    }
                })
                .collect();
            CoefficientCondition { terms }.canonical()
        }
        _ => return None,
    };
    // Scale the leading term's leading coefficient to one.
    let lead = transported
        .terms
        .last()
        .and_then(|t| vec_poly::degree(&t.coeff).map(|d| t.coeff[d].clone()));
    Some(match lead {
        Some(lead) if !lead.is_zero() => CoefficientCondition {
            terms: transported
                .terms
                .into_iter()
                .map(|t| ConditionTerm {
                    coeff: vec_poly::scale(&t.coeff, &lead.recip()),
                    ..t
                })
                .collect(),
        },
        _ => transported,
    })
}

fn push_unique(list: &mut Vec<DescribedCondition>, record: DescribedCondition) {
    if !list.iter().any(|c| c.condition == record.condition) {
        list.push(record);
    }
}

/// Transport the single relation of a certified parent description (and its
/// conditions) through an implicit map.
fn build_implicit_relation_step(
    t: &RecurrenceTransformation,
    parent: &RepresentationDescription,
    parent_rows: &[Vec<Q>],
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(ImplicitRelationStep, Vec<TransformConditionStatus>), String> {
    let description = parent
        .description
        .as_ref()
        .ok_or("the parent description has no generating description")?;
    let [source] = description.relations.as_slice() else {
        return Err("the parent description needs exactly one relation".into());
    };
    let relation = ImplicitRelation::from_json(&source.relation)?;
    let map_conditions = match t {
        RecurrenceTransformation::Diagonal { weight } => {
            weight_regularity(weight, origin).map_err(|e| e.to_string())?
        }
        RecurrenceTransformation::DerivativeAntecedent { .. } => vec![TransformConditionStatus {
            condition: "T_g bijective".into(),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: "the monomial matrix of T_g is triangular with diagonal lc(g) binom(k+m,m) \
                     != 0"
                .into(),
        }],
        RecurrenceTransformation::Derivative { order } => vec![TransformConditionStatus {
            condition: format!(
                "Q_n = D^{order} P_n with lost data {} P_n",
                super::lossy::lost_range(*order)
            ),
            status: CERTIFIED_ALL_INDICES.into(),
            detail: "defined for every row; the lost coefficients are recorded as boundary data"
                .into(),
        }],
        other => {
            return Err(format!(
                "implicit relations are not transported through {}",
                other.kind_name()
            ))
        }
    };
    let map_certified = map_conditions
        .iter()
        .all(|c| c.status == CERTIFIED_ALL_INDICES);
    let certified = evidence_is_certified(&source.evidence) && map_certified;
    let evidence = if certified {
        TRANSPORTED_IDENTITY_EVIDENCE
    } else {
        CHECKED_SUPPLIED_ROWS
    };
    let first_base = source.first_base;
    let transport = t
        .transport_implicit(&relation, first_base)
        .map_err(|e| e.to_string())?;
    let mut exceptional = Vec::new();
    for (k, condition) in &transport.exceptional_conditions {
        push_unique(
            &mut exceptional,
            DescribedCondition::new(condition, first_base, evidence, Some(*k), "Q").with_note(
                format!("[x^{k}] of the uncleared residual: H(n,{k}) = 0 at the exceptional index"),
            ),
        );
    }
    let mut boundary = Vec::new();
    let mut dropped = Vec::new();
    let carried_evidence = |record: &DescribedCondition| {
        if evidence_is_certified(&record.evidence) && map_certified {
            record.evidence.clone()
        } else {
            CHECKED_SUPPLIED_ROWS.to_string()
        }
    };
    let carried_index = |k: Option<usize>| match t {
        RecurrenceTransformation::Diagonal { .. } => k,
        RecurrenceTransformation::Derivative { order } => k.and_then(|k| k.checked_sub(*order)),
        _ => None,
    };
    for (list, records, target) in [
        (
            "exceptional_index_conditions",
            &description.exceptional_index_conditions,
            true,
        ),
        (
            "boundary_conditions",
            &description.boundary_conditions,
            false,
        ),
    ] {
        for (i, record) in records.iter().enumerate() {
            let condition = CoefficientCondition::from_json(&record.condition)?;
            match transport_condition(t, &condition) {
                Some(carried) if !carried.is_vacuous() => {
                    let described = DescribedCondition::new(
                        &carried,
                        record.first_base,
                        &carried_evidence(record),
                        carried_index(record.coefficient_index),
                        "Q",
                    )
                    .with_note(format!("the parent's {list}[{i}] carried through the map"));
                    push_unique(
                        if target {
                            &mut exceptional
                        } else {
                            &mut boundary
                        },
                        described,
                    );
                }
                Some(_) => {}
                None => dropped.push(format!("{list}[{i}]")),
            }
        }
    }
    let mut lost_conditions = Vec::new();
    if let RecurrenceTransformation::Derivative { order } = t {
        let pulled = pull_back_relation_through_derivative(&transport.relation, *order);
        let indices = valuation_exceptional_indices(&pulled, first_base)
            .unwrap_or_default()
            .into_iter()
            .filter(|k| *k <= MAX_CONDITION_COEFFICIENT_INDEX)
            .collect::<Vec<_>>();
        lost_conditions = coefficient_conditions_of(source, &indices, "P")?;
        for record in &lost_conditions {
            let condition = CoefficientCondition::from_json(&record.condition)?;
            if let Some(pushed) = push_condition_through_derivative(&condition, *order) {
                if !pushed.is_vacuous() {
                    push_unique(
                        &mut boundary,
                        DescribedCondition::new(
                            &pushed,
                            record.first_base,
                            &carried_evidence(record),
                            carried_index(record.coefficient_index),
                            "Q",
                        )
                        .with_note("[x^k] of the parent relation, free of lost coefficients"),
                    );
                }
            }
        }
        lost_conditions.extend(description.exceptional_index_conditions.iter().cloned());
        lost_conditions.extend(description.boundary_conditions.iter().cloned());
    }
    let r = transport.relation.max_shift();
    let seed_end = (first_base + r).max(origin + 1);
    if rows.len() < seed_end - origin {
        return Err("too few rows for the seeds of the transported relation".into());
    }
    let note = match t {
        RecurrenceTransformation::Derivative { .. } => {
            "transported through the lossy derivative map (a consequence of the parent relation)"
        }
        _ if transport.clearing.is_some() => {
            "transported through the diagonal map and cleared by one common multiplier"
        }
        _ => "transported by substitution P = T_g Q",
    };
    let described = DescribedRelation::new(&transport.relation, first_base, evidence, "Q")
        .with_note(format!("{note}; source: {}", source.text));
    let q_description = GeneratingDescription {
        first_index: origin,
        seeds: rows[..seed_end - origin]
            .iter()
            .map(|r| encode_poly(r))
            .collect(),
        relations: vec![described.clone()],
        exceptional_index_conditions: exceptional,
        boundary_conditions: boundary,
    };
    let assessment = assess_generating_description(&q_description)?;
    let lossy = match t {
        RecurrenceTransformation::Derivative { order } => Some(lossy_reconstruction(
            *order,
            &q_description,
            origin,
            parent_rows,
            &lost_conditions,
        )?),
        _ => None,
    };
    let (holds, checked) = description_holds(&q_description, rows, origin)?;
    Ok((
        ImplicitRelationStep {
            source_relation: source.clone(),
            relation: described,
            clearing: transport.clearing,
            certificate: transport.reconstruction,
            dropped_conditions: dropped,
            description: q_description,
            assessment,
            lossy,
            holds_on_rows: holds,
            checked_bases: checked,
        },
        map_conditions,
    ))
}

/// Whether every relation and condition of a description holds on the rows
/// at each base of its domain, and the number of relation bases checked.
fn description_holds(
    description: &GeneratingDescription,
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(bool, usize), String> {
    let mut holds = true;
    let mut checked = 0;
    for record in &description.relations {
        let relation = ImplicitRelation::from_json(&record.relation)?;
        let end = (origin + rows.len()).saturating_sub(relation.max_shift());
        for base in record.first_base.max(origin)..end {
            holds &= relation.holds_at(rows, origin, base) != Some(false);
            checked += 1;
        }
    }
    for record in description
        .exceptional_index_conditions
        .iter()
        .chain(&description.boundary_conditions)
    {
        let condition = CoefficientCondition::from_json(&record.condition)?;
        let end = (origin + rows.len()).saturating_sub(condition.max_shift());
        for base in record.first_base.max(origin)..end {
            holds &= condition.holds_at(rows, origin, base) != Some(false);
        }
    }
    Ok((holds, checked))
}

/// Independent checks of an implicit-relation step on the new rows.
fn verify_implicit_relation_step(
    t: &RecurrenceTransformation,
    record: &ImplicitRelationStep,
    rows: &[Vec<Q>],
    origin: usize,
) -> Result<(), String> {
    let source = ImplicitRelation::from_json(&record.source_relation.relation)?;
    let transport = t
        .transport_implicit(&source, record.source_relation.first_base)
        .map_err(|e| e.to_string())?;
    let relation = ImplicitRelation::from_json(&record.relation.relation)?;
    if !relation.same_relation(&transport.relation) {
        return Err("the transported relation does not replay".into());
    }
    if let Some(certificate) = &record.certificate {
        certificate.verify_identity()?;
    }
    if !record.holds_on_rows || !description_holds(&record.description, rows, origin)?.0 {
        return Err("the transported description fails on the rows".into());
    }
    if assess_generating_description(&record.description)? != record.assessment {
        return Err("the assessment does not replay".into());
    }
    if let Some(lossy) = &record.lossy {
        if assess_generating_description(&lossy.with_lost_data)? != lossy.with_lost_data_assessment
            || assess_generating_description(&lossy.without_lost_data)?
                != lossy.without_lost_data_assessment
        {
            return Err("the lossy reconstruction assessments do not replay".into());
        }
    }
    Ok(())
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
        if chain_lost_data(&chain, origin, &source_rows)? != node.lost_data {
            return Err(format!("{label}: the lost data does not replay"));
        }
        if invert_chain_with_lost_data(&chain, origin, &states, &node.lost_data)? != source_rows
            || !node.reconstruction_verified
        {
            return Err(format!("{label}: exact reconstruction of the source fails"));
        }
        for (s, step) in node.steps.iter().enumerate() {
            let parent = node
                .parent
                .and_then(|p| report.representations.get(p))
                .ok_or_else(|| format!("{label}: a step without a parent"))?;
            if let Some(record) = &step.vector_transport {
                let t = node
                    .transformation
                    .and_then(|i| report.request.transformations.get(i))
                    .ok_or_else(|| format!("{label}: unknown transformation"))?;
                let source = step
                    .source_description
                    .and_then(|i| parent.descriptions.get(i))
                    .filter(|d| d.label != FINITE_FIT_LABEL)
                    .and_then(|d| d.vector.as_ref())
                    .ok_or_else(|| format!("{label}.steps[{s}]: uncertified source"))?;
                let parent_chain = parent
                    .chain
                    .iter()
                    .map(|&i| report.request.transformations.get(i))
                    .collect::<Option<Vec<_>>>()
                    .ok_or_else(|| format!("{label}: unknown transformation index"))?;
                let offsets = chain_offsets(&parent_chain)
                    .pop()
                    .expect("at least one offset list");
                let rebuilt = build_vector_transport(
                    t,
                    &offsets,
                    source,
                    &states,
                    origin,
                    &node.map_conditions,
                    &node.components,
                )
                .map_err(|e| format!("{label}.steps[{s}]: {e}"))?;
                if &rebuilt != record || !record.holds_on_rows {
                    return Err(format!(
                        "{label}.steps[{s}]: the vector transport does not replay"
                    ));
                }
                continue;
            }
            if let Some(record) = &step.implicit_relation {
                let t = node
                    .transformation
                    .and_then(|i| report.request.transformations.get(i))
                    .ok_or_else(|| format!("{label}: unknown transformation"))?;
                let described = step
                    .source_description
                    .and_then(|i| parent.descriptions.get(i))
                    .filter(|d| d.label != FINITE_FIT_LABEL)
                    .and_then(|d| d.description.as_ref())
                    .ok_or_else(|| format!("{label}.steps[{s}]: uncertified source"))?;
                if described.relations.as_slice() != std::slice::from_ref(&record.source_relation) {
                    return Err(format!("{label}.steps[{s}]: source mismatch"));
                }
                let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
                verify_implicit_relation_step(t, record, &rows, origin)
                    .map_err(|e| format!("{label}.steps[{s}]: {e}"))?;
                continue;
            }
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
        for (s, space) in node.implicit_spaces.iter().enumerate() {
            let label = format!("{label}.implicit_spaces[{s}]");
            verify_implicit_relation_search(space).map_err(|e| format!("{label}: {e}"))?;
            let parent = node
                .parent
                .and_then(|p| report.representations.get(p))
                .ok_or_else(|| format!("{label}: a space without a parent"))?;
            let certified_source = parent.descriptions.iter().any(|d| {
                d.label != FINITE_FIT_LABEL
                    && d.scalar.as_ref().map(|json| RecurrenceJson {
                        search: None,
                        ..json.clone()
                    }) == space.request.source_recurrence
            });
            if !certified_source {
                return Err(format!(
                    "{label}: the source is not a certified description of the parent"
                ));
            }
            let map = node
                .transformation
                .and_then(|i| report.request.transformations.get(i));
            if space.request.transformation.as_ref() != map {
                return Err(format!("{label}: the map is not the edge map"));
            }
        }
        let certifying = symmetry_certifying_step(&node.steps);
        for (i, item) in node.descriptions.iter().enumerate() {
            let item_label = format!("{label}.descriptions[{i}]");
            if node.components.len() == 1 {
                let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
                check_description_on_rows(item, &rows, origin, &item_label)?;
            }
            check_vector_description(item, &states, origin, &item_label)?;
            if item.kind == "transported_gamma_system"
                && evidence_is_certified(&item.evidence)
                && step_index_of(&item.provenance)
                    .and_then(|k| node.steps.get(k))
                    .and_then(|s| s.transform.as_ref())
                    .is_none_or(|r| {
                        !r.conditions
                            .iter()
                            .all(|c| c.status == CERTIFIED_ALL_INDICES)
                            && certifying.is_none()
                    })
            {
                return Err(format!(
                    "{item_label}: a certified gamma system without a certified symmetry \
                     condition"
                ));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(BigInt::from(n))
    }

    fn constant_condition() -> CoefficientCondition {
        CoefficientCondition {
            terms: vec![
                ConditionTerm {
                    shift: 1,
                    coefficient_index: 1,
                    coeff: vec![q(1)],
                },
                ConditionTerm {
                    shift: 0,
                    coefficient_index: 0,
                    coeff: vec![q(0), q(-1)],
                },
            ],
        }
    }

    #[test]
    fn conditions_are_carried_through_implicit_maps() {
        // [x^1] P[n+1] - n [x^0] P[n] = 0.
        let condition = constant_condition();
        // Q = T P with lambda = binom(k+1, 1): [x^k] P = [x^k] Q / (k + 1).
        let binomial = RecurrenceTransformation::Diagonal {
            weight: DiagonalWeight::Binomial { m: 1 },
        };
        assert_eq!(
            transport_condition(&binomial, &condition)
                .unwrap()
                .text("Q"),
            "[x^1] Q[n+1] - 2n [x^0] Q[n] = 0"
        );
        // P = D(x Q): [x^e] P = (e + 1) [x^e] Q.
        let antecedent = RecurrenceTransformation::DerivativeAntecedent {
            g: vec!["0".into(), "1".into()],
        };
        assert_eq!(
            transport_condition(&antecedent, &condition)
                .unwrap()
                .text("Q"),
            "[x^1] Q[n+1] - 1/2n [x^0] Q[n] = 0"
        );
        // Q = D P forgets [x^0] P.
        let derivative = RecurrenceTransformation::Derivative { order: 1 };
        assert!(transport_condition(&derivative, &condition).is_none());
        let reversal = RecurrenceTransformation::Reversal {
            schedule: DegreeSchedule::new(0, 0, 1),
        };
        assert!(transport_condition(&reversal, &condition).is_none());
    }

    fn pascal_source(rule: &Recurrence) -> RecurrenceJson {
        RecurrenceJson::from_recurrence_rational(rule, 0, &[vec![q(1)]], None)
    }

    fn bivar(rows: &[&[i64]]) -> BivarPoly {
        BivarPoly::from_coefficients_normalized(
            rows.iter()
                .map(|r| r.iter().map(|&c| q(c)).collect())
                .collect(),
        )
    }

    fn gamma_step(rule: &Recurrence, source: usize) -> RepresentationStep {
        let map = RecurrenceTransformation::Gamma {
            schedule: DegreeSchedule::new(0, 0, 1),
        };
        let request = RecurrenceTransformRequest::from_source(map, pascal_source(rule), 8);
        RepresentationStep {
            source_description: Some(source),
            transform: Some(transform_recurrence(&request).unwrap()),
            implicit: None,
            implicit_relation: None,
            vector_transport: None,
            error: None,
        }
    }

    #[test]
    fn gamma_systems_are_certified_only_with_a_certified_symmetry() {
        // (1 + x)^n: an invariant rule and a rank-two rule for the same rows.
        let invariant = Recurrence {
            terms: vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 1]]))],
            denominator: None,
            inhomogeneous: None,
        };
        let rank_two = Recurrence {
            terms: vec![
                RecurrenceTerm::new(1, 0, bivar(&[&[1, 2], &[0, -1]])),
                RecurrenceTerm::new(1, 1, bivar(&[&[0, 1, 1]])),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let certifying = gamma_step(&invariant, 0);
        let coupled = gamma_step(&rank_two, 1);
        let report = coupled.transform.as_ref().unwrap();
        assert_eq!(report.conditions[0].status, CHECKED_SUPPLIED_ROWS);
        assert_eq!(
            certifying.transform.as_ref().unwrap().conditions[0].status,
            CERTIFIED_ALL_INDICES
        );
        let transported = report.transported.as_ref().unwrap();
        let item = || RepresentationDescription {
            provenance: "steps[1].transform.transported".into(),
            ..gamma_system_item(
                report,
                transported,
                transported.system.as_ref().unwrap(),
                false,
            )
            .unwrap()
        };
        let mut alone = vec![item()];
        upgrade_gamma_systems(std::slice::from_ref(&coupled), &mut alone);
        assert_eq!(alone[0].label, FINITE_FIT_LABEL);
        assert_eq!(alone[0].evidence, CHECKED_SUPPLIED_ROWS);
        // Another step of the edge certifies the symmetry of the rows.
        let steps = vec![certifying, coupled.clone()];
        assert_eq!(symmetry_certifying_step(&steps), Some(0));
        let mut items = vec![item()];
        upgrade_gamma_systems(&steps, &mut items);
        assert_eq!(items[0].evidence, TRANSPORTED_IDENTITY_EVIDENCE);
        assert_eq!(items[0].label, EQUIVALENT_DESCRIPTION_LABEL);
        assert!(items[0].note.as_ref().unwrap().contains("steps[0]"));
        let description = items[0].description.as_ref().unwrap();
        assert!(description
            .relations
            .iter()
            .all(|r| r.evidence == TRANSPORTED_IDENTITY_EVIDENCE));
        assert_eq!(
            items[0].assessment.as_ref().unwrap().status,
            EQUIVALENT_DESCRIPTION
        );
        // A step without a certified source certifies nothing.
        let mut uncertified = steps.clone();
        uncertified[0].source_description = None;
        assert_eq!(symmetry_certifying_step(&uncertified), None);
    }

    #[test]
    fn degree_bounds_compare_on_the_whole_domain() {
        let bound = DegreeSchedule::new(1, 1, 1);
        assert!(bound_at_most(&bound, &DegreeSchedule::new(1, 1, 1), 1));
        assert!(bound_at_most(&bound, &DegreeSchedule::new(0, 0, 2), 1));
        // Smaller at the start, or a smaller slope on an unbounded domain.
        assert!(!bound_at_most(&bound, &DegreeSchedule::new(1, 0, 1), 1));
        assert!(!bound_at_most(&bound, &DegreeSchedule::new(1, 5, 0), 1));
        // ... but fine on a bounded domain that ends in time.
        assert!(bound_at_most(
            &bound,
            &DegreeSchedule::new(1, 5, 0).with_last_index(5),
            1
        ));
        // A bound known only from a later index proves nothing earlier.
        assert!(!bound_at_most(
            &DegreeSchedule::new(2, 0, 1),
            &DegreeSchedule::new(0, 9, 1),
            1
        ));
    }
}
