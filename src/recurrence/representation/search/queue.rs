//! The breadth-first queue: search state, step budget, representation nodes and the
//! driver [`search_representations`].

use super::*;

/// An edge of the search tree being expanded: a transformation applied to its parent node.
pub(super) struct Edge<'a> {
    pub(super) t: &'a RecurrenceTransformation,
    pub(super) parent: &'a RepresentationNode,
}

pub(super) struct Search<'a> {
    pub(super) request: &'a RepresentationSearchRequest,
    pub(super) origin: usize,
    source_rows: Vec<Vec<Q>>,
    steps_used: usize,
    budget: Vec<RepresentationBudgetEvent>,
}

impl Search<'_> {
    pub(super) fn event(&mut self, limit: &str, representation: Option<usize>, detail: String) {
        self.budget.push(RepresentationBudgetEvent {
            limit: limit.to_string(),
            representation,
            detail,
        });
    }

    pub(super) fn take_step(&mut self, representation: usize, what: &str) -> bool {
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

    fn run_step(
        &mut self,
        id: usize,
        edge: &Edge<'_>,
        source: Option<usize>,
        states: &[Vec<Vec<Q>>],
        items: &mut Vec<RepresentationDescription>,
        conditions: &mut Option<Vec<TransformConditionStatus>>,
    ) -> RepresentationStep {
        let Edge { t, parent } = *edge;
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
            let edge = Edge { t, parent };
            self.vector_steps(nodes, &edge, &offsets, &states, &mut node, &mut items);
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
                let edge = Edge { t, parent };
                let step = self.run_step(id, &edge, source, states, items, &mut conditions);
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
