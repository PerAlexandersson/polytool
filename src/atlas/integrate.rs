//! Integration of derivations into the store
//! ([`integrate_derivations`]): match premises and conclusions to nodes,
//! create nodes, merge definitions with equality certificates, record
//! links, write edges and update node properties.
//!
//! # Rules
//!
//! 1. The certificate is analyzed ([`analyze_certificate`]; a derivation
//!    may carry the analysis its job computed, which the audit replays
//!    later).
//! 2. Every premise is matched to a node with a definition that agrees with
//!    the premise's definition ([`same_sequence_from`]) from an index at
//!    most the required first index, and with a proved property at least
//!    as strong as the requirement.  A derivation with an unmatched premise
//!    is rejected.
//! 3. The conclusion is the node one of whose definitions agrees with the
//!    certificate's input from an index at most every claim's first index
//!    (a tail of a node is the node).  Otherwise the input's first
//!    [`ATLAS_STORED_ROWS`] rows must be integers (else the derivation is
//!    rejected) and a node with the same fingerprint is looked up: the same
//!    rows from the same first index ask the equality callback for an
//!    all-`n` certificate (a merged definition) and otherwise give a new
//!    node with a candidate-equal link; different rows give a new node with
//!    a fingerprint link.  Only a derivation with a real-rootedness claim
//!    creates a node.
//! 4. Every claim is checked exactly on the node's stored rows
//!    ([`rows_satisfy_claim`]); a failure is a contradiction and the edge is
//!    rejected.
//! 5. The edge is written once (ids are content hashes) and every claim
//!    becomes (or extends) a property; a real-rootedness claim with every
//!    root in `(-inf, 0]` (or `[0, inf)`) also gives the derived properties
//!    of one-signed coefficients and, for `(-inf, 0]`, the Hadamard kernel.

use super::analysis::{analyze_certificate, canonical_definition, definition_rows};
use super::analysis::{same_sequence_from, CertificateAnalysis, PremiseRequirement};
use super::error::AtlasError;
use super::rows::{
    atlas_hash, fingerprint_rows, integer_rows, nonpositive_ray, observe_rows, region_within,
    rows_from_strings, rows_satisfy_claim, rows_to_strings,
};
use super::store::AtlasStore;
use super::types::*;
use crate::real_rootedness::RealInterval;
use crate::recurrence::{RecurrenceJson, RootInterval};
use crate::BigRational;
use num_traits::Zero;
use std::collections::{BTreeMap, HashMap};

/// Options of [`integrate_derivations`].
#[derive(Debug, Clone, Default)]
pub struct IntegrationOptions {
    /// Do not create nodes (only add edges to existing ones).
    pub no_new_nodes: bool,
}

/// What [`integrate_derivations`] did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntegrationSummary {
    /// Derivations examined.
    pub derivations: usize,
    /// Edges written.
    pub edges_added: usize,
    /// Edges already present.
    pub edges_present: usize,
    /// Nodes created.
    pub nodes_added: usize,
    /// Definitions merged into existing nodes by an equality certificate.
    pub definitions_merged: usize,
    /// Links written.
    pub links_added: usize,
    /// Rejected derivations by reason.
    pub rejected: BTreeMap<String, usize>,
    /// A few examples per rejection reason.
    pub examples: BTreeMap<String, Vec<String>>,
    /// Claims that failed on stored rows (contradictions; must stay empty).
    pub contradictions: Vec<String>,
    /// The node id of every derivation (`None` when rejected), in order.
    pub outcomes: Vec<Option<String>>,
}

impl IntegrationSummary {
    fn example(&mut self, reason: &str, text: String) {
        let examples = self.examples.entry(reason.to_string()).or_default();
        if examples.len() < 8 {
            examples.push(text);
        }
    }

    fn reject(&mut self, reason: &str) {
        *self.rejected.entry(reason.to_string()).or_default() += 1;
        self.outcomes.push(None);
    }
}

/// The key of a canonical definition.
pub(crate) fn definition_key(definition: &RecurrenceJson) -> String {
    serde_json::to_string(definition).unwrap_or_default()
}

/// Whether the proved claim `have` implies the requirement `need`.
pub(crate) fn claim_implies(have: &AtlasClaim, need: &AtlasClaim) -> bool {
    match (have, need) {
        (
            AtlasClaim::RealRooted {
                from_index: a,
                region: r,
            },
            AtlasClaim::RealRooted {
                from_index: b,
                region: s,
            },
        ) => a <= b && region_within(r, s),
        (
            AtlasClaim::Interlacing {
                from_index: a,
                kind,
                orientation,
            },
            AtlasClaim::Interlacing {
                from_index: b,
                kind: need_kind,
                orientation: need_orientation,
            },
        ) => {
            let kind_ok =
                *need_kind == crate::recurrence::RowInterlacing::Weak || kind == need_kind;
            let orientation_ok = *need_orientation
                == crate::recurrence::InterlacingOrientation::Unordered
                || orientation == need_orientation;
            a <= b && kind_ok && orientation_ok
        }
        _ => have == need,
    }
}

/// Whether the proved claims `haves` together imply `need`: one of them
/// implies it, or, for a real-rootedness requirement, the intersection of
/// the regions of every real-rootedness claim covering the required first
/// index lies in the required region (every covered row has its roots in
/// each of those regions, so in their intersection).
pub(crate) fn claims_imply<'a>(
    haves: impl IntoIterator<Item = &'a AtlasClaim>,
    need: &AtlasClaim,
) -> bool {
    let haves: Vec<&AtlasClaim> = haves.into_iter().collect();
    if haves.iter().any(|have| claim_implies(have, need)) {
        return true;
    }
    let AtlasClaim::RealRooted {
        from_index: needed_from,
        region: needed,
    } = need
    else {
        return false;
    };
    let Ok(needed) = super::rows::region_interval(needed) else {
        return false;
    };
    let mut intersection: Option<RealInterval> = None;
    for have in haves {
        if let AtlasClaim::RealRooted { from_index, region } = have {
            if from_index > needed_from {
                continue;
            }
            let Ok(interval) = super::rows::region_interval(region) else {
                continue;
            };
            intersection = Some(match intersection {
                None => interval,
                Some(current) => super::rows::intersect_intervals(&current, &interval),
            });
        }
    }
    intersection.is_some_and(|i| i.is_empty() || i.is_subset_of(&needed))
}

/// The first index of a claim, when it has one.
pub(crate) fn claim_from(claim: &AtlasClaim) -> Option<usize> {
    match claim {
        AtlasClaim::RealRooted { from_index, .. }
        | AtlasClaim::Interlacing { from_index, .. }
        | AtlasClaim::OneSignedCoefficients { from_index, .. }
        | AtlasClaim::HadamardKernel { from_index } => Some(*from_index),
        AtlasClaim::BasisMapPreserver { .. } => None,
    }
}

/// The derived claims of a real-rootedness claim (module documentation).
pub(crate) fn derived_claims(claim: &AtlasClaim) -> Vec<AtlasClaim> {
    let AtlasClaim::RealRooted { from_index, region } = claim else {
        return Vec::new();
    };
    let nonnegative_ray = RootInterval::from_real(&RealInterval::at_least(BigRational::zero()));
    if region_within(region, &nonpositive_ray()) {
        vec![
            AtlasClaim::OneSignedCoefficients {
                from_index: *from_index,
                alternating: false,
            },
            AtlasClaim::HadamardKernel {
                from_index: *from_index,
            },
        ]
    } else if region_within(region, &nonnegative_ray) {
        vec![AtlasClaim::OneSignedCoefficients {
            from_index: *from_index,
            alternating: true,
        }]
    } else {
        Vec::new()
    }
}

/// The node id of a first definition.
pub(crate) fn node_id(definition: &RecurrenceJson) -> String {
    format!("n-{}", atlas_hash(&definition_key(definition)))
}

/// Indexes over the store used while integrating.
pub(crate) struct NodeIndex {
    by_definition: HashMap<String, String>,
    by_fingerprint: BTreeMap<String, Vec<String>>,
}

impl NodeIndex {
    pub(crate) fn new(store: &AtlasStore) -> Self {
        let mut by_definition = HashMap::new();
        for node in store.nodes() {
            for definition in &node.definitions {
                by_definition.insert(definition_key(&definition.recurrence), node.id.clone());
            }
        }
        Self {
            by_definition,
            by_fingerprint: store.fingerprint_index(),
        }
    }

    fn add_definition(&mut self, definition: &RecurrenceJson, node: &str) {
        self.by_definition
            .insert(definition_key(definition), node.to_string());
    }

    fn add_node(&mut self, node: &AtlasNode) {
        for definition in &node.definitions {
            self.add_definition(&definition.recurrence, &node.id);
        }
        self.by_fingerprint
            .entry(node.fingerprint.clone())
            .or_default()
            .push(node.id.clone());
    }
}

/// Rows and fingerprint of a canonical definition, when integral.
pub(crate) fn definition_integer_rows(
    definition: &RecurrenceJson,
) -> Result<Option<Vec<Vec<num_bigint::BigInt>>>, AtlasError> {
    let rows = definition_rows(definition, ATLAS_STORED_ROWS)?;
    Ok(integer_rows(&rows))
}

/// A node `id` whose definition agrees with `definition` from an index at
/// most `at_most` (the agreement index is returned).
pub(crate) fn find_agreeing_node(
    store: &AtlasStore,
    index: &NodeIndex,
    definition: &RecurrenceJson,
    at_most: Option<usize>,
) -> Result<Option<(String, usize)>, AtlasError> {
    let canonical = canonical_definition(definition)?;
    if let Some(id) = index.by_definition.get(&definition_key(&canonical)) {
        return Ok(Some((id.clone(), canonical.first_index)));
    }
    let Some(rows) = definition_integer_rows(&canonical)? else {
        return Ok(None);
    };
    let fingerprint = fingerprint_rows(&rows);
    let Some(candidates) = index.by_fingerprint.get(&fingerprint) else {
        return Ok(None);
    };
    for id in candidates {
        let Some(node) = store.node(id) else { continue };
        for existing in &node.definitions {
            if let Some(start) = same_sequence_from(&canonical, &existing.recurrence)? {
                if at_most.is_none_or(|bound| start <= bound) && start >= node.first_index {
                    return Ok(Some((id.clone(), start)));
                }
            }
        }
    }
    Ok(None)
}

fn premise_node(
    store: &AtlasStore,
    index: &NodeIndex,
    requirement: &PremiseRequirement,
) -> Result<Result<String, &'static str>, AtlasError> {
    let needed = claim_from(&requirement.uses);
    let Some((id, _)) = find_agreeing_node(store, index, &requirement.definition, needed)? else {
        return Ok(Err("premise_not_in_atlas"));
    };
    let node = store.node(&id).expect("indexed node");
    let satisfied = claims_imply(node.properties.iter().map(|p| &p.claim), &requirement.uses);
    Ok(if satisfied {
        Ok(id)
    } else {
        Err("premise_claim_not_proved")
    })
}

/// A new node for a canonical definition with integer rows.
fn new_node(
    spec: &AtlasSequenceSpec,
    definition: &RecurrenceJson,
    rows: &[Vec<num_bigint::BigInt>],
) -> AtlasNode {
    let mut node = AtlasNode {
        schema: ATLAS_NODE_SCHEMA.into(),
        id: node_id(definition),
        revision: 0,
        first_index: definition.first_index,
        rows: rows_to_strings(rows),
        fingerprint: fingerprint_rows(rows),
        definitions: vec![AtlasDefinition {
            recurrence: definition.clone(),
            source: spec.source.clone(),
            equality_edge: None,
        }],
        description: AtlasDescription {
            name: spec.name.clone(),
            polytool_family: spec.polytool_family.clone(),
            oeis_ids: Vec::new(),
            oeis_tail_ids: Vec::new(),
            interpretation: None,
        },
        properties: Vec::new(),
        observations: observe_rows(definition.first_index, rows),
    };
    add_description(&mut node, spec);
    node
}

/// The edge id of its content.
pub(crate) fn edge_id(edge: &AtlasEdge) -> String {
    let text = format!(
        "{}|{}|{}|{}",
        serde_json::to_string(&edge.operation).unwrap_or_default(),
        serde_json::to_string(&edge.premises).unwrap_or_default(),
        edge.conclusion,
        serde_json::to_string(&edge.certificate).unwrap_or_default()
    );
    format!("e-{}", atlas_hash(&text))
}

/// The link id of its kind and nodes.
pub(crate) fn link_id(nodes: &[String], kind: &AtlasLinkKind) -> String {
    let text = format!(
        "{}|{}",
        nodes.join(","),
        serde_json::to_string(kind).unwrap_or_default()
    );
    format!("l-{}", atlas_hash(&text))
}

/// Add a link (once).
pub(crate) fn add_link(
    store: &mut AtlasStore,
    nodes: Vec<String>,
    kind: AtlasLinkKind,
    summary: &mut IntegrationSummary,
) {
    let id = link_id(&nodes, &kind);
    if store.put_link(AtlasLink {
        schema: ATLAS_LINK_SCHEMA.into(),
        id,
        nodes,
        kind,
    }) {
        summary.links_added += 1;
    }
}

fn add_property(node: &mut AtlasNode, claim: AtlasClaim, provenance: AtlasProvenance, edge: &str) {
    if let Some(property) = node
        .properties
        .iter_mut()
        .find(|p| p.claim == claim && p.provenance == provenance)
    {
        if !property.edges.iter().any(|e| e == edge) {
            property.edges.push(edge.to_string());
        }
        return;
    }
    node.properties.push(AtlasProperty {
        claim,
        provenance,
        edges: vec![edge.to_string()],
    });
}

/// The callback that certifies `new = existing` for every index (an
/// [`AtlasCertificate::DefinitionEquality`]), or `None`.
pub type EqualityCallback<'a> =
    dyn FnMut(&RecurrenceJson, &RecurrenceJson) -> Option<AtlasCertificate> + 'a;

/// Integrate the derivations in order (module documentation) and return
/// what was done.  Nothing is written to disk; call
/// [`AtlasStore::commit`].
///
/// # Errors
///
/// [`AtlasError`] only for internal failures; a derivation that cannot be
/// integrated is counted in [`IntegrationSummary::rejected`].
pub fn integrate_derivations(
    store: &mut AtlasStore,
    derivations: &[(AtlasDerivation, Option<CertificateAnalysis>)],
    options: &IntegrationOptions,
    equality: &mut EqualityCallback<'_>,
) -> Result<IntegrationSummary, AtlasError> {
    let mut summary = IntegrationSummary::default();
    let mut index = NodeIndex::new(store);
    for (derivation, analysis) in derivations {
        summary.derivations += 1;
        let analysis = match analysis {
            Some(analysis) => analysis.clone(),
            None => match analyze_certificate(&derivation.certificate) {
                Ok(analysis) => analysis,
                Err(_) => {
                    summary.reject("certificate_does_not_replay");
                    continue;
                }
            },
        };
        integrate_one(
            store,
            &mut index,
            derivation,
            &analysis,
            options,
            equality,
            &mut summary,
        )?;
    }
    Ok(summary)
}

#[allow(clippy::too_many_arguments)]
fn integrate_one(
    store: &mut AtlasStore,
    index: &mut NodeIndex,
    derivation: &AtlasDerivation,
    analysis: &CertificateAnalysis,
    options: &IntegrationOptions,
    equality: &mut EqualityCallback<'_>,
    summary: &mut IntegrationSummary,
) -> Result<(), AtlasError> {
    // Premises.
    let mut premises = Vec::new();
    for requirement in &analysis.premises {
        match premise_node(store, index, requirement)? {
            Ok(node) => premises.push(AtlasPremise {
                node,
                role: requirement.role,
                uses: requirement.uses.clone(),
            }),
            Err(reason) => {
                summary.example(
                    reason,
                    format!(
                        "{} ({}): needs {}",
                        derivation.conclusion.source,
                        derivation.operation.kind(),
                        serde_json::to_string(&requirement.uses).unwrap_or_default()
                    ),
                );
                summary.reject(reason);
                return Ok(());
            }
        }
    }
    let Ok(input) = canonical_definition(&analysis.input) else {
        summary.reject("input_not_runnable");
        return Ok(());
    };
    // Definition equalities.
    if let Some(equal_to) = &analysis.equal_to {
        let Some((id, start)) = find_agreeing_node(store, index, equal_to, None)? else {
            summary.reject("equality_target_not_in_atlas");
            return Ok(());
        };
        let node = store.node(&id).expect("indexed node").clone();
        if start != node.first_index || input.first_index != node.first_index {
            summary.reject("equality_index_mismatch");
            return Ok(());
        }
        let edge = make_edge(derivation, analysis, Vec::new(), &id);
        let edge_id = edge.id.clone();
        if store.put_edge(edge) {
            summary.edges_added += 1;
        } else {
            summary.edges_present += 1;
        }
        let mut node = node;
        if !node.definitions.iter().any(|d| d.recurrence == input) {
            node.definitions.push(AtlasDefinition {
                recurrence: input.clone(),
                source: derivation.conclusion.source.clone(),
                equality_edge: Some(edge_id),
            });
            summary.definitions_merged += 1;
        }
        add_description(&mut node, &derivation.conclusion);
        index.add_definition(&input, &id);
        store.put_node(node);
        summary.outcomes.push(Some(id));
        return Ok(());
    }
    let claims_from = analysis.claims.iter().filter_map(claim_from).min();
    let has_real_rootedness = analysis
        .claims
        .iter()
        .any(|c| matches!(c, AtlasClaim::RealRooted { .. }));
    // Conclusion node.
    let found = find_agreeing_node(store, index, &input, claims_from)?.map(|(id, _)| id);
    let node_id_value = match found {
        Some(id) => id,
        None => {
            if !has_real_rootedness || options.no_new_nodes {
                summary.reject("no_node_for_conclusion");
                return Ok(());
            }
            let Some(rows) = definition_integer_rows(&input)? else {
                summary.reject("rows_not_integral");
                return Ok(());
            };
            let fingerprint = fingerprint_rows(&rows);
            let related: Vec<String> = index
                .by_fingerprint
                .get(&fingerprint)
                .cloned()
                .unwrap_or_default();
            let mut merged_into = None;
            let mut candidate_equal = Vec::new();
            let mut fingerprint_related = Vec::new();
            for id in &related {
                let node = store.node(id).expect("indexed node");
                let same_rows = node.first_index == input.first_index
                    && rows_from_strings(&node.rows).ok().as_deref() == Some(&rows[..]);
                if !same_rows {
                    fingerprint_related.push(id.clone());
                    continue;
                }
                let existing = node.definitions[0].recurrence.clone();
                if let Some(certificate) = equality(&existing, &input) {
                    if let Ok(equality_analysis) = analyze_certificate(&certificate) {
                        let equality_derivation = AtlasDerivation {
                            conclusion: derivation.conclusion.clone(),
                            operation: AtlasOperation::DefinitionEquality,
                            certificate,
                        };
                        let before = summary.outcomes.len();
                        integrate_one(
                            store,
                            index,
                            &equality_derivation,
                            &equality_analysis,
                            options,
                            equality,
                            summary,
                        )?;
                        // The equality is bookkeeping of this derivation.
                        summary.outcomes.truncate(before);
                        if index.by_definition.contains_key(&definition_key(&input)) {
                            merged_into = Some(id.clone());
                            break;
                        }
                    }
                }
                candidate_equal.push(id.clone());
            }
            match merged_into {
                Some(id) => id,
                None => {
                    let node = new_node(&derivation.conclusion, &input, &rows);
                    let id = node.id.clone();
                    if store.node(&id).is_some() {
                        summary.reject("node_id_collision");
                        return Ok(());
                    }
                    index.add_node(&node);
                    store.put_node(node);
                    summary.nodes_added += 1;
                    for other in candidate_equal {
                        add_link(
                            store,
                            vec![other, id.clone()],
                            AtlasLinkKind::CandidateEqual {
                                rows_compared: rows.len(),
                                reason: "no all-index equality certificate was found".into(),
                            },
                            summary,
                        );
                    }
                    for other in fingerprint_related {
                        add_link(
                            store,
                            vec![other, id.clone()],
                            AtlasLinkKind::FingerprintRelated {
                                fingerprint: fingerprint.clone(),
                            },
                            summary,
                        );
                    }
                    id
                }
            }
        }
    };
    let mut node = store.node(&node_id_value).expect("node").clone();
    // Claims on the stored rows.
    let rows = rows_from_strings(&node.rows).map_err(AtlasError::InvalidRecord)?;
    for claim in &analysis.claims {
        if let Err(reason) = rows_satisfy_claim(node.first_index, &rows, claim) {
            summary.contradictions.push(format!(
                "{} ({}): {reason}",
                node.id,
                derivation.operation.kind()
            ));
            summary.reject("claim_fails_on_stored_rows");
            return Ok(());
        }
    }
    let edge = make_edge(derivation, analysis, premises, &node.id);
    let edge_id = edge.id.clone();
    if store.put_edge(edge) {
        summary.edges_added += 1;
    } else {
        summary.edges_present += 1;
    }
    for claim in &analysis.claims {
        add_property(&mut node, claim.clone(), AtlasProvenance::Edges, &edge_id);
        for derived in derived_claims(claim) {
            add_property(
                &mut node,
                derived,
                AtlasProvenance::DerivedFromEdges,
                &edge_id,
            );
        }
    }
    if node.definitions.iter().any(|d| d.recurrence == input) {
        add_description(&mut node, &derivation.conclusion);
    }
    store.put_node(node);
    summary.outcomes.push(Some(node_id_value));
    Ok(())
}

fn add_description(node: &mut AtlasNode, spec: &AtlasSequenceSpec) {
    if let Some(id) = &spec.oeis_id {
        let ids = if spec.oeis_prefix_rows == 0 {
            &mut node.description.oeis_ids
        } else {
            &mut node.description.oeis_tail_ids
        };
        if !ids.contains(id) {
            ids.push(id.clone());
            ids.sort();
        }
    }
    if node.description.polytool_family.is_none() {
        node.description.polytool_family = spec.polytool_family.clone();
    }
}

fn make_edge(
    derivation: &AtlasDerivation,
    analysis: &CertificateAnalysis,
    premises: Vec<AtlasPremise>,
    conclusion: &str,
) -> AtlasEdge {
    let mut edge = AtlasEdge {
        schema: ATLAS_EDGE_SCHEMA.into(),
        id: String::new(),
        operation: derivation.operation.clone(),
        premises,
        conclusion: conclusion.to_string(),
        claims: analysis.claims.clone(),
        theorem: analysis.theorem.clone(),
        lean_route: analysis.lean_route.clone(),
        lean_declarations: analysis.lean_declarations.clone(),
        certificate: derivation.certificate.clone(),
        verification: AtlasVerification {
            polytool: "replayed".into(),
            lean: ATLAS_LEAN_NOT_ATTEMPTED.into(),
        },
    };
    edge.id = edge_id(&edge);
    edge
}
