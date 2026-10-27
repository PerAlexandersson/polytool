//! Scalar descriptions: description costs, ranking, merging and row helpers.

use super::*;

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

pub(super) fn rule_cost(rule: &Recurrence) -> DescriptionCost {
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

pub(super) fn relation_cost(relation: &ImplicitRelationJson) -> DescriptionCost {
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

pub(super) fn vector_cost(json: &TransformVectorRecurrenceJson) -> DescriptionCost {
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

pub(super) struct ScalarInput<'a> {
    pub(super) kind: &'a str,
    pub(super) rule: &'a Recurrence,
    pub(super) json: RecurrenceJson,
    pub(super) evidence: &'a str,
    pub(super) valid_from: usize,
    pub(super) bounded_domain: bool,
    pub(super) symbol: &'a str,
    pub(super) provenance: String,
}

pub(super) fn scalar_description(input: ScalarInput<'_>) -> RepresentationDescription {
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
        let first_base = input.valid_from.checked_sub(r).ok_or_else(|| {
            RecurrenceError::InvalidInput("the rule starts before its leader shift".to_string())
        })?;
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
pub(super) fn merge_descriptions(
    items: Vec<RepresentationDescription>,
) -> Vec<RepresentationDescription> {
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

pub(super) fn find_equal_sequence(
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

pub(super) fn encode_states(states: &[Vec<Vec<Q>>]) -> Vec<Vec<Vec<String>>> {
    states
        .iter()
        .map(|s| s.iter().map(|c| encode_poly(c)).collect())
        .collect()
}

pub(super) fn decode_states(
    states: &[Vec<Vec<String>>],
) -> Result<Vec<Vec<Vec<Q>>>, RecurrenceError> {
    states
        .iter()
        .map(|s| s.iter().map(|c| decode_poly(c)).collect())
        .collect()
}

pub(super) fn scalar_rows(node: &RepresentationNode) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    node.states
        .iter()
        .map(|s| decode_poly(&s[0]))
        .collect::<Result<Vec<_>, _>>()
}

pub(super) fn component_names(
    t: Option<&RecurrenceTransformation>,
    chain_len: usize,
) -> Vec<String> {
    match t {
        None => vec![if chain_len == 0 { "P" } else { "Q" }.to_string()],
        Some(t) => t.record().components.into_iter().map(|c| c.name).collect(),
    }
}

/// `(lag, derivative order)` pairs of the default complete-space menu,
/// simplest first.
pub(super) fn menu_pairs(bounds: &TransformSearchBounds) -> Vec<(usize, usize)> {
    let mut pairs = (1..=bounds.max_rec_len.max(1))
        .flat_map(|l| (0..=bounds.max_diff_deg).map(move |j| (l, j)))
        .collect::<Vec<_>>();
    pairs.sort_by_key(|&(l, j)| (l + j, l));
    pairs
}

/// Maps whose complete spaces come from the direct implicit search.
pub(super) fn implicit_space_map(t: &RecurrenceTransformation) -> bool {
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
