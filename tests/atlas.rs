//! The real-rooted atlas engine: storage, integration, audit, growth and
//! OEIS matching on small examples.

use polytool::atlas::grow::{forward_derivations, ForwardTask};
use polytool::atlas::obreschkoff::obreschkoff_operator;
use polytool::atlas::oeis_match::{
    match_atlas, OeisMatchOptions, OeisTermsIndex, OEIS_TERMS_CACHE_VERSION,
};
use polytool::atlas::seed::{
    certify_definition_equality, named_family_sources, root_derivations, AnalyzedDerivation,
};
use polytool::atlas::*;
use polytool::recurrence::{DegreeSchedule, RecurrenceJson};
use std::path::PathBuf;

fn scratch(name: &str) -> PathBuf {
    let directory =
        std::env::temp_dir().join(format!("polytool-atlas-test-{name}-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&directory);
    directory
}

fn named(name: &str) -> AtlasSequenceSpec {
    named_family_sources()
        .into_iter()
        .find(|s| s.name == name)
        .expect("a named family")
}

fn integrate(store: &mut AtlasStore, derivations: Vec<AnalyzedDerivation>) -> IntegrationSummary {
    let list: Vec<_> = derivations.into_iter().map(|(d, a)| (d, Some(a))).collect();
    let mut equality = |existing: &RecurrenceJson, new: &RecurrenceJson| {
        certify_definition_equality(existing, new).ok()
    };
    integrate_derivations(store, &list, &IntegrationOptions::default(), &mut equality).unwrap()
}

fn node_of(store: &AtlasStore, family: &str) -> AtlasNode {
    store
        .nodes()
        .iter()
        .find(|n| n.description.polytool_family.as_deref() == Some(family))
        .expect("node")
        .clone()
}

#[test]
fn eulerian_roots_round_trip_and_audit() {
    let directory = scratch("roots");
    let mut store = AtlasStore::open(&directory).unwrap();
    let derivations = root_derivations(&named("eulerian")).unwrap();
    assert!(!derivations.is_empty());
    let summary = integrate(&mut store, derivations);
    assert_eq!(summary.nodes_added, 1, "{summary:?}");
    assert!(summary.contradictions.is_empty());
    store.commit().unwrap();
    let reopened = AtlasStore::open(&directory).unwrap();
    assert_eq!(reopened.nodes().len(), 1);
    let node = &reopened.nodes()[0];
    assert_eq!(node.first_index, 1);
    assert_eq!(node.rows[3], vec!["1", "11", "11", "1"]);
    assert!(node.observations.palindromic && node.observations.gamma_positive);
    assert!(node.properties.iter().any(|p| matches!(
        &p.claim,
        AtlasClaim::RealRooted { region, .. } if region.to_string() == "(-inf, 0)"
    )));
    assert!(node
        .properties
        .iter()
        .any(|p| matches!(p.claim, AtlasClaim::Interlacing { .. })));
    assert!(node
        .properties
        .iter()
        .any(|p| matches!(p.claim, AtlasClaim::HadamardKernel { .. })));
    let audit = audit_atlas(&reopened);
    assert!(audit.passed(), "{:?}", audit.problems);
    assert_eq!(audit.edges_accepted, reopened.edges().len());
    // A sampled audit with an empty sample replays only the nodes with
    // several incoming edges (here the one node) and still checks every row.
    let sampled = audit_atlas_with(
        &reopened,
        &AtlasAuditOptions {
            sample_nodes: Some(0),
            sample_seed: 7,
        },
    );
    assert!(sampled.passed(), "{:?}", sampled.problems);
    assert_eq!(
        sampled.edges_accepted + sampled.edges_trusted,
        reopened.edges().len()
    );

    // Integrating the same derivations again adds nothing.
    let mut again = AtlasStore::open(&directory).unwrap();
    let summary = integrate(&mut again, root_derivations(&named("eulerian")).unwrap());
    assert_eq!(summary.edges_added, 0);
    assert_eq!(summary.nodes_added, 0);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn audit_rejects_tampered_records() {
    let directory = scratch("tamper");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(&named("narayana")).unwrap());
    store.commit().unwrap();
    // A claim that the certificate does not give.
    let edges = std::fs::read_to_string(directory.join("edges.jsonl")).unwrap();
    let tampered = edges.replacen("\"lower\":", "\"lowered\":", 0).replacen(
        "\"from_index\":1",
        "\"from_index\":0",
        1,
    );
    assert_ne!(edges, tampered);
    std::fs::write(directory.join("edges.jsonl"), tampered).unwrap();
    let audit = audit_atlas(&AtlasStore::open(&directory).unwrap());
    assert!(!audit.passed());
    // Restore the edges, then change a stored row.
    std::fs::write(directory.join("edges.jsonl"), &edges).unwrap();
    assert!(audit_atlas(&AtlasStore::open(&directory).unwrap()).passed());
    let nodes = std::fs::read_to_string(directory.join("nodes.jsonl")).unwrap();
    std::fs::write(
        directory.join("nodes.jsonl"),
        nodes.replacen("[\"1\",\"3\",\"1\"]", "[\"1\",\"3\",\"2\"]", 1),
    )
    .unwrap();
    let audit = audit_atlas(&AtlasStore::open(&directory).unwrap());
    assert!(
        audit.problems.iter().any(|p| p.contains("stored rows")),
        "{:?}",
        audit.problems
    );
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn forward_gamma_and_obreschkoff_from_eulerian() {
    let directory = scratch("forward");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(&named("eulerian")).unwrap());
    let eulerian = node_of(&store, "eulerian");
    let spec = AtlasSequenceSpec {
        definition: eulerian.definitions[0].recurrence.clone(),
        source: eulerian.definitions[0].source.clone(),
        name: "eulerian".into(),
        polytool_family: None,
        oeis_id: None,
        oeis_prefix_rows: 0,
    };
    // Gamma polynomials of A_n (n >= 1): d_n = n - 1.
    let gamma = forward_derivations(&ForwardTask::GammaPolynomial {
        node: eulerian.id.clone(),
        parent: spec.clone(),
        schedule: DegreeSchedule::new(1, 0, 1),
    })
    .unwrap();
    let summary = integrate(&mut store, gamma);
    assert_eq!(summary.nodes_added, 1, "{summary:?}");
    // P_n + P_(n+1) of interlacing Eulerian rows.
    let (_, start) = obreschkoff_operator(&spec.definition).unwrap();
    assert_eq!(start, 1);
    let combination = forward_derivations(&ForwardTask::Obreschkoff {
        node: eulerian.id.clone(),
        parent: spec,
        first_index: 1,
    })
    .unwrap();
    let summary = integrate(&mut store, combination);
    assert_eq!(summary.nodes_added, 1, "{summary:?}");
    store.commit().unwrap();
    let reopened = AtlasStore::open(&directory).unwrap();
    let combination = reopened
        .nodes()
        .iter()
        .find(|n| n.description.name.starts_with("obreschkoff"))
        .unwrap();
    // A_1 + A_2 = 2 + x, A_2 + A_3 = 2 + 5x + x^2.
    assert_eq!(combination.rows[0], vec!["2", "1"]);
    assert_eq!(combination.rows[1], vec!["2", "5", "1"]);
    let audit = audit_atlas(&reopened);
    assert!(audit.passed(), "{:?}", audit.problems);
    assert_eq!(audit.edges_by_operation.get("gamma_polynomial"), Some(&1));
    assert_eq!(
        audit.edges_by_operation.get("obreschkoff_combination"),
        Some(&1)
    );
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn equal_definitions_are_merged_with_a_certificate() {
    let directory = scratch("equality");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(&named("binomial")).unwrap());
    // (1 + x)^n by the rule n P_n = n (1 + x) P_(n-1), which differs from the
    // named rule P_n = (1 + x) P_(n-1) but gives the same rows.
    let alternative: RecurrenceJson = serde_json::from_value(serde_json::json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {
            "terms": [
                {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "0"], ["1", "1"]]}}
            ],
            "denominator": {"coeffs": [["0"], ["1"]]}
        }
    }))
    .unwrap();
    let spec = AtlasSequenceSpec {
        definition: canonical_definition(&alternative).unwrap(),
        source: "test:index_factor_binomial".into(),
        name: "binomial (index factor)".into(),
        polytool_family: None,
        oeis_id: None,
        oeis_prefix_rows: 0,
    };
    let derivations = root_derivations(&spec).unwrap();
    assert!(!derivations.is_empty());
    let summary = integrate(&mut store, derivations);
    assert_eq!(summary.nodes_added, 0, "{summary:?}");
    assert_eq!(summary.definitions_merged, 1, "{summary:?}");
    store.commit().unwrap();
    let reopened = AtlasStore::open(&directory).unwrap();
    assert_eq!(reopened.nodes().len(), 1);
    assert_eq!(reopened.nodes()[0].definitions.len(), 2);
    let audit = audit_atlas(&reopened);
    assert!(audit.passed(), "{:?}", audit.problems);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn oeis_matching_records_links() {
    let directory = scratch("oeis");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(&named("eulerian")).unwrap());
    store.commit().unwrap();
    let cache = directory.join("terms.tsv");
    let eulerian: Vec<String> = [
        1, 1, 1, 1, 4, 1, 1, 11, 11, 1, 1, 26, 66, 26, 1, 1, 57, 302, 302, 57, 1,
    ]
    .iter()
    .map(ToString::to_string)
    .collect();
    let factorials: Vec<String> = [
        1u64, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800,
    ]
    .iter()
    .map(ToString::to_string)
    .collect();
    std::fs::write(
        &cache,
        format!(
            "{OEIS_TERMS_CACHE_VERSION}\nA008292\t1\tnonn,tabl\t{}\nA000142\t0\tnonn\t{}\n",
            eulerian.join(","),
            factorials.join(",")
        ),
    )
    .unwrap();
    let options = OeisMatchOptions::default();
    let index = OeisTermsIndex::load(&cache, &options).unwrap();
    let summary = match_atlas(&mut store, &index, &options);
    store.commit().unwrap();
    assert!(
        summary.hits_by_target.contains_key("Triangle"),
        "{summary:?}"
    );
    assert!(
        summary.hits_by_target.contains_key("RowSums"),
        "{summary:?}"
    );
    let report = atlas_report(&AtlasStore::open(&directory).unwrap());
    assert!(
        report.oeis_links.contains_key("triangle:catalog"),
        "{report:?}"
    );
    assert!(
        report.oeis_links.contains_key("row_sums:not_in_catalog"),
        "{report:?}"
    );
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn catalog_entries_with_prefix_rows_are_tails() {
    let (sources, _) = polytool::atlas::seed::catalog_sources();
    // A122431, the Riordan array ((1+x)^3, x): its prefix rows 1, 3 + x,
    // 3 + 3x + x^2 come before the rule; 3 + 3x + x^2 is not real-rooted.
    let tail = sources
        .iter()
        .find(|s| s.oeis_id.as_deref() == Some("A122431"))
        .unwrap();
    assert_eq!(tail.oeis_prefix_rows, 3);
    let directory = scratch("tail");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(tail).unwrap());
    let node = &store.nodes()[0];
    assert!(node.description.oeis_ids.is_empty());
    assert_eq!(node.description.oeis_tail_ids, vec!["A122431".to_string()]);
    let report = atlas_report(&store);
    assert!(report.identified_known_not_real_rooted.is_empty());
    assert_eq!(
        report.tail_identified_known_not_real_rooted,
        vec!["A122431".to_string()]
    );
    let _ = std::fs::remove_dir_all(&directory);
}

/// Forward Veronese sections of the Eulerian node: the section rules are
/// computed from the parent rule, certified, integrated and audited.
#[test]
fn forward_veronese_sections_from_eulerian() {
    let directory = scratch("forward-veronese");
    let mut store = AtlasStore::open(&directory).unwrap();
    integrate(&mut store, root_derivations(&named("eulerian")).unwrap());
    let eulerian = node_of(&store, "eulerian");
    let spec = AtlasSequenceSpec {
        definition: eulerian.definitions[0].recurrence.clone(),
        source: eulerian.definitions[0].source.clone(),
        name: "eulerian".into(),
        polytool_family: None,
        oeis_id: None,
        oeis_prefix_rows: 0,
    };
    for residue in 0..2 {
        let sections = forward_derivations(&ForwardTask::VeroneseSection {
            node: eulerian.id.clone(),
            parent: spec.clone(),
            modulus: 2,
            residue,
        })
        .unwrap();
        let summary = integrate(&mut store, sections);
        assert_eq!(summary.nodes_added, 1, "{summary:?}");
    }
    store.commit().unwrap();
    let reopened = AtlasStore::open(&directory).unwrap();
    let even = reopened
        .nodes()
        .iter()
        .find(|n| n.description.name.starts_with("veronese_section_2_0"))
        .unwrap();
    // V_(2,0)(A_4) = V_(2,0)(1 + 11x + 11x^2 + x^3) = 1 + 11y.
    assert_eq!(even.rows[3], vec!["1", "11"]);
    let audit = audit_atlas(&reopened);
    assert!(audit.passed(), "{:?}", audit.problems);
    assert_eq!(audit.edges_by_operation.get("veronese_section"), Some(&2));
    let _ = std::fs::remove_dir_all(&directory);
}
