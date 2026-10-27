//! Vector representations: fits and transport of first-order vector systems and the
//! two-relation gamma system.

use super::*;

impl Search<'_> {
    /// Map conditions and transported vector systems of an edge from a
    /// vector parent (a componentwise map).
    pub(super) fn vector_steps(
        &mut self,
        nodes: &[RepresentationNode],
        edge: &Edge<'_>,
        offsets: &[i64],
        states: &[Vec<Vec<Q>>],
        node: &mut RepresentationNode,
        items: &mut Vec<RepresentationDescription>,
    ) {
        let Edge { t, parent } = *edge;
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
    pub(super) fn vector_fits(
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
}

/// Uniqueness of a first-order vector system with its initial state: the
/// states are determined when `q(n, .)` is a nonzero polynomial at every
/// index `n >= first_index`, certified by a coefficient `[x^j] q(n,x)` (a
/// polynomial in `n`) with no integer root `n >= first_index`.
pub(super) fn assess_vector_system(
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
pub(super) fn vector_label(
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
pub(super) fn gamma_system_item(
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
pub(super) fn step_index_of(provenance: &str) -> Option<usize> {
    provenance
        .strip_prefix("steps[")?
        .split(']')
        .next()?
        .parse()
        .ok()
}

/// The step of a gamma edge whose report certifies the symmetry condition
/// (its source rule is invariant under reversal), if any.
pub(super) fn symmetry_certifying_step(steps: &[RepresentationStep]) -> Option<usize> {
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
pub(super) fn upgrade_gamma_systems(
    steps: &[RepresentationStep],
    items: &mut [RepresentationDescription],
) {
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
pub(super) fn bound_at_most(
    bound: &DegreeSchedule,
    schedule: &DegreeSchedule,
    from: usize,
) -> bool {
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
pub(super) fn build_vector_transport(
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
pub(super) fn check_vector_description(
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
