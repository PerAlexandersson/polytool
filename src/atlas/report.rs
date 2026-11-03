//! Summary statistics of an atlas directory ([`atlas_report`]).

use super::store::AtlasStore;
use super::types::{AtlasClaim, AtlasLinkKind, AtlasOperation, OeisMatchTarget};
use serde::Serialize;
use std::collections::{BTreeMap, BTreeSet};

/// Counts of an atlas directory.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize)]
pub struct AtlasReport {
    /// Nodes (current revisions).
    pub nodes: usize,
    /// Nodes by the source of their first definition (`named_family`,
    /// `oeis_catalog`, `derived`).
    pub nodes_by_source: BTreeMap<String, usize>,
    /// Distinct OEIS catalog ids that are definitions of a node.
    pub catalog_ids_identified: usize,
    /// Catalog entries flagged `lean_proved` (benchmarks and match targets).
    pub lean_flagged_total: usize,
    /// Of those, the ones that are definitions of a node (reached
    /// automatically; the Lean proof is an independent cross-check).
    pub lean_flagged_reached: usize,
    /// Of those, the ones whose recurrence part (after the entry's prefix
    /// rows) is a definition of a node, while the whole entry is not.
    pub lean_flagged_reached_as_tail: usize,
    /// Catalog ids whose recurrence part (after prefix rows) is a definition
    /// of a node.
    pub catalog_tail_ids_identified: usize,
    /// Tail identifications of entries flagged `known_not_real_rooted` (the
    /// failing rows are then prefix rows outside the node; listed for
    /// review).
    pub tail_identified_known_not_real_rooted: Vec<String>,
    /// Catalog ids that are definitions of a node but are not flagged
    /// `lean_proved`.
    pub identified_not_lean_flagged: Vec<String>,
    /// Flagged entries flagged `known_not_real_rooted` that are definitions
    /// of a node (a contradiction; must stay empty).
    pub identified_known_not_real_rooted: Vec<String>,
    /// Edges.
    pub edges: usize,
    /// Edges by operation kind.
    pub edges_by_operation: BTreeMap<String, usize>,
    /// Links by kind.
    pub links_by_kind: BTreeMap<String, usize>,
    /// OEIS links by target (`triangle`, `reversed_triangle`, `row_sums`)
    /// and whether the entry is in the bundled catalog.
    pub oeis_links: BTreeMap<String, usize>,
    /// Nodes with a triangle or reversed-triangle OEIS link to an entry that
    /// is not in the bundled catalog (with those A-numbers).
    pub nodes_matched_outside_catalog: Vec<(String, String, Vec<String>)>,
    /// Nodes with at least two incoming proof edges (equalities and map
    /// theorems not counted), with the operation labels of the edges
    /// (`theorem_check:<theorem>`, `parent_provers:<route>`, or the
    /// operation kind with the premise node ids).
    pub nodes_with_several_proofs: Vec<(String, String, Vec<String>)>,
    /// Of those, the nodes whose edges use at least two different routes
    /// (distinct labels after dropping the root theorem's variant): several
    /// genuinely different proofs.
    pub nodes_with_several_routes: usize,
    /// Nodes by their tightest proved root region (as text).
    pub nodes_by_region: BTreeMap<String, usize>,
    /// Nodes with a proved interlacing property.
    pub nodes_with_interlacing: usize,
    /// Nodes with a proved basis-map property.
    pub nodes_defining_basis_maps: usize,
}

fn source_kind(source: &str) -> String {
    source.split(':').next().unwrap_or("").to_string()
}

/// The counts of `store`.
pub fn atlas_report(store: &AtlasStore) -> AtlasReport {
    let mut report = AtlasReport {
        nodes: store.nodes().len(),
        edges: store.edges().len(),
        ..AtlasReport::default()
    };
    let mut catalog_ids = BTreeSet::new();
    let mut tail_ids = BTreeSet::new();
    for node in store.nodes() {
        *report
            .nodes_by_source
            .entry(source_kind(&node.definitions[0].source))
            .or_default() += 1;
        catalog_ids.extend(node.description.oeis_ids.iter().cloned());
        tail_ids.extend(node.description.oeis_tail_ids.iter().cloned());
        let mut regions: Vec<String> = node
            .properties
            .iter()
            .filter_map(|p| match &p.claim {
                AtlasClaim::RealRooted { region, .. } => Some(region.to_string()),
                _ => None,
            })
            .collect();
        regions.sort_by_key(|r| (r == "(-inf, +inf)", r.contains("inf"), r.clone()));
        if let Some(region) = regions.first() {
            *report.nodes_by_region.entry(region.clone()).or_default() += 1;
        }
        if node
            .properties
            .iter()
            .any(|p| matches!(p.claim, AtlasClaim::Interlacing { .. }))
        {
            report.nodes_with_interlacing += 1;
        }
        if node
            .properties
            .iter()
            .any(|p| matches!(p.claim, AtlasClaim::BasisMapPreserver { .. }))
        {
            report.nodes_defining_basis_maps += 1;
        }
    }
    report.catalog_ids_identified = catalog_ids.len();
    report.catalog_tail_ids_identified = tail_ids.difference(&catalog_ids).count();
    for (entry, flag) in crate::oeis::entries_with_real_rootedness_flag() {
        match flag.status {
            crate::oeis::OeisRealRootednessStatus::LeanProved => {
                report.lean_flagged_total += 1;
                if catalog_ids.contains(entry.id) {
                    report.lean_flagged_reached += 1;
                } else if tail_ids.contains(entry.id) {
                    report.lean_flagged_reached_as_tail += 1;
                }
            }
            crate::oeis::OeisRealRootednessStatus::KnownNotRealRooted => {
                if catalog_ids.contains(entry.id) {
                    report
                        .identified_known_not_real_rooted
                        .push(entry.id.to_string());
                }
                if tail_ids.contains(entry.id) {
                    report
                        .tail_identified_known_not_real_rooted
                        .push(entry.id.to_string());
                }
            }
        }
    }
    let flagged: BTreeSet<&str> = crate::oeis::entries_with_real_rootedness_flag()
        .filter(|(_, flag)| flag.status == crate::oeis::OeisRealRootednessStatus::LeanProved)
        .map(|(entry, _)| entry.id)
        .collect();
    report.identified_not_lean_flagged = catalog_ids
        .iter()
        .filter(|id| !flagged.contains(id.as_str()))
        .cloned()
        .collect();
    let mut incoming: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for edge in store.edges() {
        *report
            .edges_by_operation
            .entry(edge.operation.kind().to_string())
            .or_default() += 1;
        if !matches!(
            edge.operation,
            AtlasOperation::DefinitionEquality | AtlasOperation::BasisMapTheorem
        ) {
            let label = match &edge.operation {
                AtlasOperation::TheoremCheck { theorem } => format!("theorem_check:{theorem}"),
                AtlasOperation::ParentProvers { route } => format!("parent_provers:{route}"),
                other => {
                    let premises: BTreeSet<&str> =
                        edge.premises.iter().map(|p| p.node.as_str()).collect();
                    format!(
                        "{}<-{}",
                        other.kind(),
                        premises.into_iter().collect::<Vec<_>>().join("+")
                    )
                }
            };
            incoming
                .entry(edge.conclusion.clone())
                .or_default()
                .push(label);
        }
    }
    for (node, operations) in incoming {
        if operations.len() >= 2 {
            let routes: BTreeSet<&str> = operations
                .iter()
                .map(|label| label.split(':').next().unwrap_or(label))
                .collect();
            let distinct: BTreeSet<&String> = operations.iter().collect();
            if routes.len() >= 2 || (distinct.len() >= 2 && !routes.contains("theorem_check")) {
                report.nodes_with_several_routes += 1;
            }
            let name = store
                .node(&node)
                .map(|n| n.description.name.clone())
                .unwrap_or_default();
            report
                .nodes_with_several_proofs
                .push((node, name, operations));
        }
    }
    let mut outside: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for link in store.links() {
        *report
            .links_by_kind
            .entry(link.kind.kind().to_string())
            .or_default() += 1;
        if let AtlasLinkKind::OeisMatch {
            oeis_id,
            in_catalog,
            alignment,
        } = &link.kind
        {
            let target = match alignment.target {
                OeisMatchTarget::Triangle => "triangle",
                OeisMatchTarget::ReversedTriangle => "reversed_triangle",
                OeisMatchTarget::RowSums => "row_sums",
            };
            let key = format!(
                "{target}:{}",
                if *in_catalog {
                    "catalog"
                } else {
                    "not_in_catalog"
                }
            );
            *report.oeis_links.entry(key).or_default() += 1;
            if !*in_catalog && alignment.target != OeisMatchTarget::RowSums {
                outside
                    .entry(link.nodes[0].clone())
                    .or_default()
                    .insert(oeis_id.clone());
            }
        }
    }
    for (node, ids) in outside {
        let name = store
            .node(&node)
            .map(|n| n.description.name.clone())
            .unwrap_or_default();
        report
            .nodes_matched_outside_catalog
            .push((node, name, ids.into_iter().collect()));
    }
    report
}
