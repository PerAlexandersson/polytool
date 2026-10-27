//! Replay: re-derivation and verification of a whole representation search.

use super::*;

fn check_description_on_rows(
    item: &RepresentationDescription,
    rows: &[Vec<Q>],
    origin: usize,
    label: &str,
) -> Result<(), RecurrenceError> {
    if let Some(json) = &item.scalar {
        let (rule, first, seeds) = json.to_recurrence_parts()?;
        if first < origin {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: the rule starts before the rows"
            )));
        }
        for (j, seed) in seeds.iter().enumerate() {
            if rows.get(first + j - origin).map(|r| vec_poly::trim(r)) != Some(vec_poly::trim(seed))
            {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: seed {j} differs from the rows"
                )));
            }
        }
        for target in first + seeds.len()..origin + rows.len() {
            if transformed_rule_holds_at(&rule, rows, origin, target) == Some(false) {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: the rule fails at target {target}"
                )));
            }
        }
    }
    if let Some(description) = &item.description {
        let assessment = assess_generating_description(description)?;
        if item.assessment.as_ref() != Some(&assessment) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: the assessment does not replay"
            )));
        }
        if representation_label(&assessment.status) != item.label
            && !(item.label == CERTIFIED_RELATION_LABEL
                && assessment.status == EQUIVALENT_DESCRIPTION)
        {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: the label does not match the assessment"
            )));
        }
        for relation in &description.relations {
            let relation_rows = ImplicitRelation::from_json(&relation.relation)?;
            let end = (origin + rows.len()).saturating_sub(relation_rows.max_shift());
            for base in relation.first_base.max(origin)..end {
                if relation_rows.holds_at(rows, origin, base) == Some(false) {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}: a relation fails at base {base}"
                    )));
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
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}: a condition fails at base {base}"
                    )));
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
pub fn verify_representation_search(
    report: &RepresentationSearchReport,
) -> Result<(), RecurrenceError> {
    if report.schema != RECURRENCE_REPRESENTATION_SEARCH_SCHEMA {
        return Err(RecurrenceError::InvalidInput(format!(
            "unsupported schema `{}`",
            report.schema
        )));
    }
    let replay = search_representations(&report.request)?;
    if &replay != report {
        return Err(RecurrenceError::VerificationFailed(
            "report does not match its exact replay".into(),
        ));
    }
    let root = report
        .representations
        .first()
        .ok_or(RecurrenceError::VerificationFailed(
            "the report has no source representation".into(),
        ))?;
    let origin = report.source.row_origin;
    let source_rows = scalar_rows(root)?;
    for node in &report.representations {
        let label = format!("representations[{}]", node.id);
        let chain = node
            .chain
            .iter()
            .map(|&i| report.request.transformations.get(i))
            .collect::<Option<Vec<_>>>()
            .ok_or_else(|| {
                RecurrenceError::VerificationFailed(format!(
                    "{label}: unknown transformation index"
                ))
            })?;
        let states = decode_states(&node.states)?;
        if apply_chain(&chain, origin, &source_rows)? != states {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: rows are not the image of the source rows"
            )));
        }
        if chain_lost_data(&chain, origin, &source_rows)? != node.lost_data {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: the lost data does not replay"
            )));
        }
        if invert_chain_with_lost_data(&chain, origin, &states, &node.lost_data)? != source_rows
            || !node.reconstruction_verified
        {
            return Err(RecurrenceError::VerificationFailed(format!(
                "{label}: exact reconstruction of the source fails"
            )));
        }
        for (s, step) in node.steps.iter().enumerate() {
            let parent = node
                .parent
                .and_then(|p| report.representations.get(p))
                .ok_or_else(|| {
                    RecurrenceError::VerificationFailed(format!("{label}: a step without a parent"))
                })?;
            if let Some(record) = &step.vector_transport {
                let t = node
                    .transformation
                    .and_then(|i| report.request.transformations.get(i))
                    .ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}: unknown transformation"
                        ))
                    })?;
                let source = step
                    .source_description
                    .and_then(|i| parent.descriptions.get(i))
                    .filter(|d| d.label != FINITE_FIT_LABEL)
                    .and_then(|d| d.vector.as_ref())
                    .ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}.steps[{s}]: uncertified source"
                        ))
                    })?;
                let parent_chain = parent
                    .chain
                    .iter()
                    .map(|&i| report.request.transformations.get(i))
                    .collect::<Option<Vec<_>>>()
                    .ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}: unknown transformation index"
                        ))
                    })?;
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
                .map_err(|e| e.with_context(format!("{label}.steps[{s}]")))?;
                if &rebuilt != record || !record.holds_on_rows {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}.steps[{s}]: the vector transport does not replay"
                    )));
                }
                continue;
            }
            if let Some(record) = &step.implicit_relation {
                let t = node
                    .transformation
                    .and_then(|i| report.request.transformations.get(i))
                    .ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}: unknown transformation"
                        ))
                    })?;
                let described = step
                    .source_description
                    .and_then(|i| parent.descriptions.get(i))
                    .filter(|d| d.label != FINITE_FIT_LABEL)
                    .and_then(|d| d.description.as_ref())
                    .ok_or_else(|| {
                        RecurrenceError::VerificationFailed(format!(
                            "{label}.steps[{s}]: uncertified source"
                        ))
                    })?;
                if described.relations.as_slice() != std::slice::from_ref(&record.source_relation) {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}.steps[{s}]: source mismatch"
                    )));
                }
                let rows = states.iter().map(|s| s[0].clone()).collect::<Vec<_>>();
                verify_implicit_relation_step(t, record, &rows, origin)
                    .map_err(|e| e.with_context(format!("{label}.steps[{s}]")))?;
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
                        .ok_or_else(|| {
                            RecurrenceError::VerificationFailed(format!(
                                "{label}.steps[{s}]: uncertified source"
                            ))
                        })
                })
                .transpose()?;
            if let Some(transform) = &step.transform {
                verify_recurrence_transform(transform)?;
                if transform.request.source_recurrence != source {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}.steps[{s}]: source mismatch"
                    )));
                }
            }
            if let Some(implicit) = &step.implicit {
                verify_recurrence_implicit(implicit)?;
                if implicit.request.source_recurrence != source {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "{label}.steps[{s}]: source mismatch"
                    )));
                }
            }
        }
        for space in &node.spaces {
            verify_recurrence_spaces(&space.request, space)?;
        }
        for (s, space) in node.implicit_spaces.iter().enumerate() {
            let label = format!("{label}.implicit_spaces[{s}]");
            verify_implicit_relation_search(space).map_err(|e| e.with_context(&label))?;
            let parent = node
                .parent
                .and_then(|p| report.representations.get(p))
                .ok_or_else(|| {
                    RecurrenceError::VerificationFailed(format!(
                        "{label}: a space without a parent"
                    ))
                })?;
            let certified_source = parent.descriptions.iter().any(|d| {
                d.label != FINITE_FIT_LABEL
                    && d.scalar.as_ref().map(|json| RecurrenceJson {
                        search: None,
                        ..json.clone()
                    }) == space.request.source_recurrence
            });
            if !certified_source {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: the source is not a certified description of the parent"
                )));
            }
            let map = node
                .transformation
                .and_then(|i| report.request.transformations.get(i));
            if space.request.transformation.as_ref() != map {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{label}: the map is not the edge map"
                )));
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
                return Err(RecurrenceError::VerificationFailed(format!(
                    "{item_label}: a certified gamma system without a certified symmetry \
                     condition"
                )));
            }
        }
    }
    Ok(())
}
