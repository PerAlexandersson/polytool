//! The audit of an atlas directory ([`audit_atlas`]): every edge is
//! replayed in file order from the roots, every claim re-derived, every
//! premise checked against the claims already established, every node
//! property checked against its cited edges, and every node's rows
//! regenerated from its definition and checked exactly against its claims.

use super::analysis::{analyze_certificate, canonical_definition, same_sequence_from};
use super::error::AtlasError;
use super::integrate::{
    claim_from, claims_imply, definition_integer_rows, derived_claims, node_id,
};
use super::rows::{fingerprint_rows, observe_rows, rows_from_strings, rows_satisfy_claim};
use super::store::AtlasStore;
use super::types::*;
use crate::recurrence::RecurrenceJson;
use std::collections::{BTreeMap, HashMap, HashSet};

/// The result of [`audit_atlas`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AtlasAuditReport {
    /// Nodes (current revisions).
    pub nodes: usize,
    /// Edges.
    pub edges: usize,
    /// Links.
    pub links: usize,
    /// Edges whose certificate replayed and whose claims and premises
    /// check.
    pub edges_accepted: usize,
    /// Edges taken as recorded by a sampled audit (not replayed).
    pub edges_trusted: usize,
    /// Nodes whose incoming edges were replayed (all nodes in a full audit).
    pub nodes_replayed: usize,
    /// Node properties checked against their edges.
    pub properties_checked: usize,
    /// Claims checked exactly on stored rows.
    pub row_claims_checked: usize,
    /// Edges per operation kind.
    pub edges_by_operation: BTreeMap<String, usize>,
    /// Every problem found (an empty list is a passing audit).
    pub problems: Vec<String>,
}

impl AtlasAuditReport {
    /// Whether the audit found no problem.
    pub fn passed(&self) -> bool {
        self.problems.is_empty()
    }
}

/// The first index from which `definition` agrees with an accepted
/// definition of the node.
fn agreement(
    accepted: &[RecurrenceJson],
    definition: &RecurrenceJson,
) -> Result<Option<usize>, AtlasError> {
    let canonical = canonical_definition(definition)?;
    let mut best: Option<usize> = None;
    for existing in accepted {
        if *existing == canonical {
            return Ok(Some(canonical.first_index));
        }
        if let Some(start) = same_sequence_from(&canonical, existing)? {
            best = Some(best.map_or(start, |b: usize| b.min(start)));
        }
    }
    Ok(best)
}

fn check_node_rows(node: &AtlasNode) -> Result<Vec<Vec<num_bigint::BigInt>>, String> {
    if node.schema != ATLAS_NODE_SCHEMA {
        return Err(format!("unsupported schema `{}`", node.schema));
    }
    let first = node
        .definitions
        .first()
        .ok_or("the node has no definition")?;
    let canonical = canonical_definition(&first.recurrence).map_err(|e| e.to_string())?;
    if canonical != first.recurrence {
        return Err("the first definition is not canonical".into());
    }
    if node.id != node_id(&canonical) {
        return Err("the id is not the hash of the first definition".into());
    }
    if node.first_index != canonical.first_index {
        return Err("the first index differs from the definition's".into());
    }
    let rows = definition_integer_rows(&canonical)
        .map_err(|e| e.to_string())?
        .ok_or("the definition's rows are not integral")?;
    let stored = rows_from_strings(&node.rows)?;
    if stored != rows {
        return Err("the stored rows differ from the definition's rows".into());
    }
    if node.fingerprint != fingerprint_rows(&rows) {
        return Err("the fingerprint differs from its recomputation".into());
    }
    if node.observations != observe_rows(node.first_index, &rows) {
        return Err("the observations differ from their recomputation".into());
    }
    Ok(rows)
}

/// What [`audit_atlas_with`] replays.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AtlasAuditOptions {
    /// `None`: replay every edge (the full audit).  `Some(k)`: replay the
    /// incoming edges of a pseudo-random sample of `k` nodes (the smallest
    /// hashes of `seed:id`), of every node with at least two incoming edges
    /// and of every node in a candidate-equal link; the other edges are
    /// taken as recorded (their claims count as established without a
    /// replay).  Every node's rows, id, fingerprint, observations and every
    /// property's exact check on the stored rows are checked in both modes.
    pub sample_nodes: Option<usize>,
    /// The seed of the sample.
    pub sample_seed: u64,
}

/// The full audit of every record of the store (module documentation).
/// The report lists every problem; nothing is changed.
pub fn audit_atlas(store: &AtlasStore) -> AtlasAuditReport {
    audit_atlas_with(store, &AtlasAuditOptions::default())
}

/// The nodes whose incoming edges a sampled audit replays.
fn selected_nodes(store: &AtlasStore, options: &AtlasAuditOptions) -> Option<HashSet<String>> {
    let k = options.sample_nodes?;
    let mut hashed: Vec<(String, String)> = store
        .nodes()
        .iter()
        .map(|n| {
            (
                super::rows::atlas_hash(&format!("{}:{}", options.sample_seed, n.id)),
                n.id.clone(),
            )
        })
        .collect();
    hashed.sort();
    let mut selected: HashSet<String> = hashed.into_iter().take(k).map(|(_, id)| id).collect();
    let mut incoming: HashMap<&str, usize> = HashMap::new();
    for edge in store.edges() {
        *incoming.entry(edge.conclusion.as_str()).or_default() += 1;
    }
    selected.extend(
        incoming
            .into_iter()
            .filter(|(_, count)| *count >= 2)
            .map(|(id, _)| id.to_string()),
    );
    for link in store.links() {
        if matches!(link.kind, AtlasLinkKind::CandidateEqual { .. }) {
            selected.extend(link.nodes.iter().cloned());
        }
    }
    Some(selected)
}

/// The audit of [`audit_atlas`], with a sampled replay on request
/// ([`AtlasAuditOptions`]).
pub fn audit_atlas_with(store: &AtlasStore, options: &AtlasAuditOptions) -> AtlasAuditReport {
    let selected = selected_nodes(store, options);
    let nodes_replayed = selected.as_ref().map_or(store.nodes().len(), HashSet::len);
    let mut report = AtlasAuditReport {
        nodes: store.nodes().len(),
        edges: store.edges().len(),
        links: store.links().len(),
        nodes_replayed,
        ..AtlasAuditReport::default()
    };
    let mut rows_of: HashMap<String, Vec<Vec<num_bigint::BigInt>>> = HashMap::new();
    let mut accepted: HashMap<String, Vec<RecurrenceJson>> = HashMap::new();
    for node in store.nodes() {
        match check_node_rows(node) {
            Ok(rows) => {
                rows_of.insert(node.id.clone(), rows);
                accepted.insert(
                    node.id.clone(),
                    vec![node.definitions[0].recurrence.clone()],
                );
            }
            Err(reason) => report.problems.push(format!("node {}: {reason}", node.id)),
        }
    }
    let mut established: HashMap<String, Vec<(AtlasClaim, String)>> = HashMap::new();
    let mut accepted_edges: HashSet<String> = HashSet::new();
    for edge in store.edges() {
        *report
            .edges_by_operation
            .entry(edge.operation.kind().to_string())
            .or_default() += 1;
        let replay = selected
            .as_ref()
            .is_none_or(|selected| selected.contains(&edge.conclusion));
        if !replay {
            report.edges_trusted += 1;
            accepted_edges.insert(edge.id.clone());
            if let Some(node) = store.node(&edge.conclusion) {
                let defs = accepted.entry(edge.conclusion.clone()).or_default();
                for definition in &node.definitions {
                    if !defs.contains(&definition.recurrence) {
                        defs.push(definition.recurrence.clone());
                    }
                }
            }
            let entry = established.entry(edge.conclusion.clone()).or_default();
            for claim in &edge.claims {
                entry.push((claim.clone(), edge.id.clone()));
            }
            continue;
        }
        match audit_edge(edge, &rows_of, &mut accepted, &established) {
            Ok(checked) => {
                report.row_claims_checked += checked;
                report.edges_accepted += 1;
                accepted_edges.insert(edge.id.clone());
                let entry = established.entry(edge.conclusion.clone()).or_default();
                for claim in &edge.claims {
                    entry.push((claim.clone(), edge.id.clone()));
                }
            }
            Err(reason) => report.problems.push(format!("edge {}: {reason}", edge.id)),
        }
    }
    for node in store.nodes() {
        let node_accepted = accepted.get(&node.id).cloned().unwrap_or_default();
        for (i, definition) in node.definitions.iter().enumerate().skip(1) {
            let backed = definition.equality_edge.as_ref().is_some_and(|id| {
                accepted_edges.contains(id)
                    && store.edge(id).is_some_and(|e| {
                        e.conclusion == node.id && e.operation == AtlasOperation::DefinitionEquality
                    })
            });
            if !backed || !node_accepted.contains(&definition.recurrence) {
                report.problems.push(format!(
                    "node {}: definition {i} is not backed by an accepted equality edge",
                    node.id
                ));
            }
        }
        for property in &node.properties {
            report.properties_checked += 1;
            if property.edges.is_empty() {
                report
                    .problems
                    .push(format!("node {}: a property cites no edge", node.id));
            }
            for id in &property.edges {
                let Some(edge) = store.edge(id).filter(|_| accepted_edges.contains(id)) else {
                    report.problems.push(format!(
                        "node {}: property {} cites edge {id}, which is missing or rejected",
                        node.id,
                        property.claim.kind()
                    ));
                    continue;
                };
                let backed = edge.conclusion == node.id
                    && match property.provenance {
                        AtlasProvenance::Edges => edge.claims.contains(&property.claim),
                        AtlasProvenance::DerivedFromEdges => edge
                            .claims
                            .iter()
                            .any(|c| derived_claims(c).contains(&property.claim)),
                    };
                if !backed {
                    report.problems.push(format!(
                        "node {}: property {} is not a conclusion of edge {id}",
                        node.id,
                        property.claim.kind()
                    ));
                }
            }
            if let Some(rows) = rows_of.get(&node.id) {
                report.row_claims_checked += 1;
                if let Err(reason) = rows_satisfy_claim(node.first_index, rows, &property.claim) {
                    report
                        .problems
                        .push(format!("node {}: contradiction: {reason}", node.id));
                }
            }
        }
    }
    for link in store.links() {
        if link.nodes.iter().any(|id| store.node(id).is_none()) {
            report
                .problems
                .push(format!("link {}: a node is missing", link.id));
        }
    }
    report
}

/// Audit one edge; returns the number of claims checked on rows.
fn audit_edge(
    edge: &AtlasEdge,
    rows_of: &HashMap<String, Vec<Vec<num_bigint::BigInt>>>,
    accepted: &mut HashMap<String, Vec<RecurrenceJson>>,
    established: &HashMap<String, Vec<(AtlasClaim, String)>>,
) -> Result<usize, String> {
    if edge.schema != ATLAS_EDGE_SCHEMA {
        return Err(format!("unsupported schema `{}`", edge.schema));
    }
    if edge.id != super::integrate::edge_id(edge) {
        return Err("the id is not the hash of the edge".into());
    }
    let analysis = analyze_certificate(&edge.certificate).map_err(|e| e.to_string())?;
    if analysis.claims != edge.claims {
        return Err("the claims differ from the certificate's".into());
    }
    let conclusion_rows = rows_of
        .get(&edge.conclusion)
        .ok_or("the conclusion node is missing or invalid")?;
    let node_accepted = accepted
        .get(&edge.conclusion)
        .cloned()
        .ok_or("the conclusion node is missing")?;
    let input = canonical_definition(&analysis.input).map_err(|e| e.to_string())?;
    if let Some(equal_to) = &analysis.equal_to {
        if edge.operation != AtlasOperation::DefinitionEquality {
            return Err("an equality certificate on another operation".into());
        }
        let first = node_accepted
            .first()
            .map(|d| d.first_index)
            .unwrap_or_default();
        let start = agreement(&node_accepted, equal_to).map_err(|e| e.to_string())?;
        if start != Some(first) || input.first_index != first {
            return Err("the equality does not cover the node from its first index".into());
        }
        accepted
            .get_mut(&edge.conclusion)
            .expect("accepted definitions")
            .push(input);
        return Ok(0);
    }
    let start = agreement(&node_accepted, &input)
        .map_err(|e| e.to_string())?
        .ok_or("the certificate's sequence is not a definition of the conclusion node")?;
    if let Some(from) = analysis.claims.iter().filter_map(claim_from).min() {
        if start > from {
            return Err("the certificate's sequence agrees with the node too late".into());
        }
    }
    if edge.premises.len() != analysis.premises.len() {
        return Err("the premises differ from the certificate's".into());
    }
    for (premise, requirement) in edge.premises.iter().zip(&analysis.premises) {
        if premise.uses != requirement.uses || premise.role != requirement.role {
            return Err("a premise differs from the certificate's requirement".into());
        }
        let premise_accepted = accepted
            .get(&premise.node)
            .ok_or_else(|| format!("premise node {} is missing", premise.node))?;
        let needed = claim_from(&requirement.uses);
        let start = agreement(premise_accepted, &requirement.definition)
            .map_err(|e| e.to_string())?
            .ok_or_else(|| format!("premise node {} does not define the premise", premise.node))?;
        if needed.is_some_and(|n| start > n) {
            return Err(format!(
                "premise node {} agrees with the premise only from {start}",
                premise.node
            ));
        }
        let holds = established
            .get(&premise.node)
            .is_some_and(|claims| claims_imply(claims.iter().map(|(c, _)| c), &requirement.uses));
        if !holds {
            return Err(format!(
                "premise node {} has no established claim implying {}",
                premise.node,
                requirement.uses.kind()
            ));
        }
    }
    let first_index = accepted
        .get(&edge.conclusion)
        .and_then(|d| d.first().map(|d| d.first_index))
        .unwrap_or_default();
    for claim in &edge.claims {
        rows_satisfy_claim(first_index, conclusion_rows, claim)
            .map_err(|reason| format!("contradiction on stored rows: {reason}"))?;
    }
    Ok(edge.claims.len())
}
