//! Scalar steps of the queue: coordinate changes, complete spaces and implicit
//! relations applied to a representation.

use super::*;

impl Search<'_> {
    /// Certified descriptions of the parent with one implicit relation
    /// (scalar rules included, in implicit form), usable by the implicit
    /// machinery.
    pub(super) fn implicit_source_descriptions(&self, parent: &RepresentationNode) -> Vec<usize> {
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

    pub(super) fn implicit_relation_step(
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

    pub(super) fn transform_items(
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

    pub(super) fn implicit_items(
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
    pub(super) fn spaces(
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
    pub(super) fn implicit_spaces(
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
}
