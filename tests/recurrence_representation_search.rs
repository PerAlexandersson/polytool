//! Bounded representation search (`polytool.recurrence-representation-search.v1`):
//! the acceptance examples of the combined transformation search.

use polytool::recurrence::*;
use serde_json::json;

fn request(value: serde_json::Value) -> RepresentationSearchRequest {
    serde_json::from_value(value).expect("request")
}

fn run(value: serde_json::Value) -> RepresentationSearchReport {
    let report = search_representations(&request(value)).unwrap();
    verify_representation_search(&report).unwrap();
    report
}

fn q(text: &str) -> BigRational {
    parse_rational_coeff(text).unwrap()
}

fn bivar(rows: &[&[&str]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|row| row.iter().map(|c| q(c)).collect())
            .collect(),
    )
}

fn node_with_chain(report: &RepresentationSearchReport, chain: &[usize]) -> RepresentationNode {
    report
        .representations
        .iter()
        .find(|n| n.chain == chain)
        .unwrap_or_else(|| panic!("no representation with chain {chain:?}"))
        .clone()
}

/// The description whose scalar rule has the canonical form of `expected`.
fn find_rule<'a>(
    node: &'a RepresentationNode,
    expected: &Recurrence,
) -> Option<&'a RepresentationDescription> {
    node.descriptions.iter().find(|d| {
        d.scalar.as_ref().is_some_and(|json| {
            let (rule, _, _) = json.to_recurrence_parts().unwrap();
            same_canonical_recurrence(&rule, expected)
        })
    })
}

fn factorial_request() -> serde_json::Value {
    json!({
        "oeis_id": "A007318",
        "row_count": 10,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
        ],
        "limits": {"max_depth": 1}
    })
}

#[test]
fn a007318_divided_by_factorials_recovers_derivative_and_three_term_rules() {
    let report = run(factorial_request());
    assert_eq!(report.schema, RECURRENCE_REPRESENTATION_SEARCH_SCHEMA);
    assert_eq!(report.status, "complete");
    assert_eq!(report.source.kind, "oeis");
    assert_eq!(report.source.row_origin, 0);
    let node = node_with_chain(&report, &[0]);
    // Q_3 = 1 + 3x + 3x^2/2 + x^3/6.
    assert_eq!(node.states[3][0], vec!["1", "3", "3/2", "1/6"]);
    assert!(node.reconstruction_verified);

    // n Q_n = (n + x) Q_(n-1) + x Q'_(n-1), the first-order derivative rule.
    let derivative = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&["0", "1"], &["1"]])),
            RecurrenceTerm::new(1, 1, bivar(&[&["0", "1"]])),
        ],
        denominator: Some(bivar(&[&["0"], &["1"]])),
        inhomogeneous: None,
    };
    // n Q_n = (2n - 1 + x) Q_(n-1) - (n - 1) Q_(n-2), the three-term pure-shift rule.
    let three_term = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&["-1", "1"], &["2"]])),
            RecurrenceTerm::new(2, 0, bivar(&[&["1"], &["-1"]])),
        ],
        denominator: Some(bivar(&[&["0"], &["1"]])),
        inhomogeneous: None,
    };
    for (expected, seeds) in [(&derivative, 1), (&three_term, 2)] {
        let found = find_rule(&node, expected).expect("certified rule");
        assert_eq!(found.kind, "certified_space");
        assert_eq!(found.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
        assert_eq!(found.label, EQUIVALENT_DESCRIPTION_LABEL);
        assert_eq!(
            found.scalar.as_ref().unwrap().initial_polynomials.len(),
            seeds
        );
        assert!(found.provenance.starts_with("spaces["));
        // The fit menu found the same relation; the certified record is kept.
        assert!(found.also_found_by.iter().any(|p| p == "fits"), "{found:?}");
    }
    let cost = |d: &RepresentationDescription| d.cost.clone();
    assert!(
        cost(find_rule(&node, &derivative).unwrap()) < cost(find_rule(&node, &three_term).unwrap())
    );
    // The complete spaces are kept: an absent shape has a dual certificate.
    assert!(node
        .spaces
        .iter()
        .any(|s| s.results[0].status == "absent_in_shape"));
    // The transported implicit rule with its exceptional constant condition.
    let implicit = node
        .descriptions
        .iter()
        .find(|d| d.kind == "implicit_transport")
        .unwrap();
    assert_eq!(implicit.text, "D Q[n+1] - D Q[n] - Q[n] = 0");
    assert_eq!(implicit.label, EQUIVALENT_DESCRIPTION_LABEL);
    // The global ranking starts with the simplest description.
    assert_eq!(report.ranked[0].representation, 0);
    assert_eq!(report.ranked[0].text, "P(n) = (1 + t) P(n-1)");
}

fn eulerian_rows(count: usize) -> Vec<Vec<String>> {
    polytool::oeis::by_id("A008292")
        .unwrap()
        .generate_rows(count)
        .unwrap()
        .into_iter()
        .map(|row| row.iter().map(|c| c.to_string()).collect())
        .collect()
}

#[test]
fn a065826_divided_by_k_plus_one_then_the_f_polynomial_map() {
    let report = run(json!({
        "oeis_id": "A065826",
        "row_count": 9,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "inverse_binomial", "m": 1}},
            {"kind": "f_polynomial",
             "schedule": {"first_index": 1, "first_degree": 0, "step": 1}}
        ]
    }));
    let eulerian = node_with_chain(&report, &[0]);
    let rows: Vec<Vec<String>> = eulerian.states.iter().map(|s| s[0].clone()).collect();
    assert_eq!(rows, eulerian_rows(9));
    let fitted = eulerian
        .descriptions
        .iter()
        .position(|d| d.kind == "reconstruction_certified_fit")
        .expect("certified Eulerian rule");
    assert_eq!(
        eulerian.descriptions[fitted].label,
        EQUIVALENT_DESCRIPTION_LABEL
    );
    assert_eq!(
        eulerian.descriptions[fitted].text,
        "P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)"
    );

    // Q_(n+1) = (x (1 + x) Q_n)' = (1 + 2x) Q_n + (x + x^2) Q_n'.
    let f = node_with_chain(&report, &[0, 1]);
    let expected = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&["1", "2"]])),
            RecurrenceTerm::new(1, 1, bivar(&[&["0", "1", "1"]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let found = find_rule(&f, &expected).expect("transported f-polynomial rule");
    assert_eq!(found.kind, "transported_rule");
    assert_eq!(found.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert_eq!(found.label, EQUIVALENT_DESCRIPTION_LABEL);
    assert_eq!(
        found.scalar.as_ref().unwrap().initial_polynomials,
        vec![vec!["1".to_string()]]
    );
    // The step used the certified Eulerian rule of its parent as defining rule.
    let step = f.steps.iter().find(|s| s.transform.is_some()).unwrap();
    let source = step.source_description.unwrap();
    assert_eq!(
        eulerian.descriptions[source].kind,
        "reconstruction_certified_fit"
    );
    assert_eq!(f.map.segments.len(), 2);
    assert_eq!(
        f.map.segments[1],
        "Q_n(x) = (x + 1)^(n - 1) P_n((x)/(x + 1))"
    );
}

fn power_source(base: &str) -> serde_json::Value {
    // P_n = (1 + base x)^n, P_0 = 1.
    json!({"schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
            "coeff": {"coeffs": [["1", base]]}}]}})
}

fn reversal(first_degree: usize) -> serde_json::Value {
    json!({"kind": "reversal",
        "schedule": {"first_index": 0, "first_degree": first_degree, "step": 1}})
}

#[test]
fn reversal_twice_deduplicates_only_at_the_same_degree_bound() {
    let report = run(json!({
        "source_recurrence": power_source("2"),
        "row_count": 8,
        "transformations": [reversal(0)]
    }));
    let twice = report
        .edges
        .iter()
        .find(|e| e.parent == 1 && e.transformation == 0)
        .unwrap();
    assert_eq!(twice.outcome, "duplicate_map");
    assert_eq!(twice.representation, Some(0));
    assert!(twice.detail.as_ref().unwrap().contains("identity"));
    assert_eq!(report.representations.len(), 2);

    // Degree bounds n and n + 1: the composite is x P_n, a new representation.
    let report = run(json!({
        "source_recurrence": power_source("2"),
        "row_count": 8,
        "transformations": [reversal(0), reversal(1)]
    }));
    let shifted = report
        .edges
        .iter()
        .find(|e| e.parent == 1 && e.transformation == 1)
        .unwrap();
    assert_eq!(shifted.outcome, "new_representation");
    let node = &report.representations[shifted.representation.unwrap()];
    assert_eq!(node.map.segments, vec!["Q_n(x) = x P_n(x)"]);
    assert_eq!(node.states[2][0], vec!["0", "1", "4", "4"]);
    for (parent, transformation) in [(1, 0), (2, 1)] {
        let edge = report
            .edges
            .iter()
            .find(|e| e.parent == parent && e.transformation == transformation)
            .unwrap();
        assert_eq!(edge.outcome, "duplicate_map", "{edge:?}");
        assert_eq!(edge.representation, Some(0));
    }
    // Reversal at bound n of a row of degree n + 1 is rejected, not dropped.
    let rejected = report
        .edges
        .iter()
        .find(|e| e.parent == 2 && e.transformation == 0)
        .unwrap();
    assert_eq!(rejected.outcome, "rejected");
}

#[test]
fn equal_rows_alone_never_deduplicate_a_map() {
    // (1 + x)^n is palindromic: its reversal has the same rows, but the map is
    // not the identity. The sequences are equal only by certified relations.
    let report = run(json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "transformations": [reversal(0)]
    }));
    assert_eq!(report.edges[0].outcome, "new_representation");
    let node = &report.representations[1];
    assert_eq!(node.states, report.representations[0].states);
    assert!(!node.map.identity);
    let equal = node.equal_sequence_as.as_ref().expect("certified equality");
    assert_eq!(equal.representation, 0);
    assert!(!node.expanded);
}

#[test]
fn lost_constant_is_a_consequence_and_root_domain_is_optional() {
    let report = run(factorial_request());
    let node = node_with_chain(&report, &[0]);
    let implicit = node
        .descriptions
        .iter()
        .find(|d| d.kind == "implicit_transport")
        .unwrap();
    let mut description = implicit.description.clone().unwrap();
    assert_eq!(description.exceptional_index_conditions.len(), 1);
    description.exceptional_index_conditions.clear();
    let assessment = assess_generating_description(&description).unwrap();
    assert_eq!(
        representation_label(&assessment.status),
        CERTIFIED_RELATION_LABEL
    );

    let eulerian = |root_domain: Option<serde_json::Value>| {
        let mut value = json!({
            "oeis_id": "A008292",
            "row_count": 7,
            "transformations": [{"kind": "f_polynomial",
                "schedule": {"first_index": 1, "first_degree": 0, "step": 1}}],
            "limits": {"max_depth": 1}
        });
        if let Some(root_domain) = root_domain {
            value["root_domain"] = root_domain;
        }
        run(value)
    };
    let plain = eulerian(None);
    assert!(plain.representations[1].steps.iter().all(|s| s
        .transform
        .as_ref()
        .unwrap()
        .root_domain
        .is_none()));
    let interval = json!({"lower": {"value": "-1"}, "upper": {"value": "0"}});
    let with_interval = eulerian(Some(json!({
        "candidate_intervals": [interval],
        "assumptions": [{"subject": "transformed", "property": "real_rooted_in",
            "interval": interval, "note": "supplied"}]
    })));
    let step = with_interval.representations[1].steps[0]
        .transform
        .as_ref()
        .unwrap();
    let root = step.root_domain.as_ref().unwrap();
    let evidence = root.claims.iter().map(|c| c.evidence).collect::<Vec<_>>();
    assert!(evidence.contains(&RootEvidence::SuppliedAssumption));
    assert!(evidence.contains(&RootEvidence::CompletedRecurrenceInduction));
    // Algebraic discovery does not depend on the interval.
    assert_eq!(
        plain.representations[1].descriptions,
        with_interval.representations[1].descriptions
    );
}

#[test]
fn reports_round_trip_and_tampering_is_detected() {
    let report = run(factorial_request());
    let text = serde_json::to_string(&report).unwrap();
    let back: RepresentationSearchReport = serde_json::from_str(&text).unwrap();
    assert_eq!(back, report);
    verify_representation_search(&back).unwrap();

    let mut tampered = report.clone();
    tampered.representations[1].states[2][0][1] = "3".into();
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report.clone();
    tampered.representations[1].descriptions[0].label = FINITE_FIT_LABEL.into();
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report.clone();
    tampered.request.row_count = Some(11);
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report;
    tampered.schema = "polytool.recurrence-representation-search.v0".into();
    assert!(verify_representation_search(&tampered).is_err());
}

#[test]
fn rows_only_sources_give_finite_fits_only() {
    let rows: Vec<Vec<String>> = (0..9usize)
        .map(|n| {
            (0..=n)
                .map(|k| polytool::numeric::binomial_bigint(n, k).to_string())
                .collect()
        })
        .collect();
    let report = run(json!({
        "rows": rows,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
        ],
        "limits": {"max_depth": 1}
    }));
    assert_eq!(report.source.kind, "rows_only");
    for node in &report.representations {
        assert!(!node.descriptions.is_empty());
        for d in &node.descriptions {
            assert_eq!(d.label, FINITE_FIT_LABEL, "{d:?}");
        }
    }
    let node = &report.representations[1];
    assert!(node.steps[0].source_description.is_none());
    assert!(node
        .notes
        .iter()
        .any(|n| n.contains("complete spaces need")));
}

#[test]
fn limits_are_reported_as_budget_outcomes() {
    let report = run(json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}},
            reversal(0)
        ],
        "limits": {"max_depth": 1, "max_representations": 2}
    }));
    assert_eq!(report.status, "budget_exhausted");
    assert_eq!(report.representations.len(), 2);
    assert_eq!(report.edges[1].outcome, "skipped_budget");
    assert!(report
        .budget
        .iter()
        .any(|b| b.limit == "max_representations"));

    let report = run(json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
        ],
        "limits": {"max_depth": 1, "max_steps": 3}
    }));
    assert_eq!(report.status, "budget_exhausted");
    assert_eq!(report.steps_used, 3);
    assert!(report.budget.iter().any(|b| b.limit == "max_steps"));

    // Depth is a requested bound: the frontier is listed, not a budget event.
    let report = run(json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "fits": false,
        "transformations": [reversal(0)],
        "limits": {"max_depth": 0}
    }));
    assert_eq!(report.status, "complete");
    assert_eq!(report.depth_frontier, vec![0]);

    // Invalid requests are errors.
    assert!(search_representations(&request(json!({
        "oeis_id": "A007318",
        "transformations": [{"kind": "mobius", "a": "1", "b": "1", "c": "1", "e": "1",
            "schedule": {"first_index": 0, "first_degree": 0, "step": 1}}]
    })))
    .is_err());
    assert!(search_representations(&request(json!({"transformations": []}))).is_err());
}

#[test]
fn lossy_derivative_of_eulerian_rows_needs_the_lost_constants() {
    let report = run(json!({
        "oeis_id": "A008292",
        "row_count": 9,
        "transformations": [{"kind": "derivative", "order": 1}],
        "fits": false,
        "limits": {"max_depth": 1}
    }));
    let node = node_with_chain(&report, &[0]);
    assert!(node.reconstruction_verified);
    // The lost data: the constant terms A_n(0) = 1.
    assert_eq!(node.lost_data.len(), 1);
    assert_eq!(node.lost_data[0].order, 1);
    assert!(node.lost_data[0].rows.iter().all(|r| r == &vec!["1"]));
    let step = node.steps[0].implicit_relation.as_ref().unwrap();
    assert_eq!(step.source_relation.evidence, DEFINING_EVIDENCE);
    assert_eq!(step.relation.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert_eq!(
        step.relation.text,
        "D Q[n+1] + (-x + x^2) D^2 Q[n] + (-3 + 4 x - n x) D Q[n] + (2 - 2 n) Q[n] = 0"
    );
    step.certificate
        .as_ref()
        .unwrap()
        .verify_identity()
        .unwrap();
    // The derivative rows alone: a certified consequence.
    let description = node
        .descriptions
        .iter()
        .find(|d| d.kind == "implicit_transport")
        .unwrap();
    assert_eq!(description.label, CERTIFIED_RELATION_LABEL);
    // The parent through Q: equivalent with the lost constants, a consequence without.
    let lossy = step.lossy.as_ref().unwrap();
    assert_eq!(
        lossy.with_lost_data_assessment.status,
        EQUIVALENT_DESCRIPTION
    );
    assert_eq!(
        lossy.without_lost_data_assessment.status,
        CERTIFIED_CONSEQUENCE
    );
    assert!(lossy
        .with_lost_data
        .boundary_conditions
        .iter()
        .any(|c| c.text == "[x^0] P[n+1] - [x^0] P[n] = 0"));
    assert!(!lossy.lost_data_conditions.is_empty());
    assert!(lossy.without_lost_data.boundary_conditions.is_empty());

    // The lost data is part of the replay.
    let mut tampered = report.clone();
    tampered.representations[1].lost_data[0].rows[2] = vec!["2".into()];
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report;
    tampered.representations[1].steps[0]
        .implicit_relation
        .as_mut()
        .unwrap()
        .lossy
        .as_mut()
        .unwrap()
        .without_lost_data_assessment
        .status = EQUIVALENT_DESCRIPTION.into();
    assert!(verify_representation_search(&tampered).is_err());
}

#[test]
fn implicit_only_parent_feeds_the_next_diagonal_step() {
    let report = run(json!({
        "oeis_id": "A007318",
        "row_count": 10,
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
        ],
        "fits": false,
        "certified_spaces": false,
        "search": {"max_rec_len": 1, "max_var_deg": 0, "max_idx_deg": 0,
                   "max_diff_deg": 0, "max_candidates": 50},
        "limits": {"max_depth": 2}
    }));
    // Representation 1 has the implicit rule only.
    let child = node_with_chain(&report, &[0]);
    assert!(child.descriptions.iter().all(|d| d.scalar.is_none()));
    assert_eq!(child.descriptions[0].text, "D Q[n+1] - D Q[n] - Q[n] = 0");
    // Representation 2 transports it (with its constant condition).
    let grandchild = node_with_chain(&report, &[0, 0]);
    assert_eq!(grandchild.steps.len(), 1);
    let step = grandchild.steps[0].implicit_relation.as_ref().unwrap();
    assert_eq!(step.source_relation.text, "D Q[n+1] - D Q[n] - Q[n] = 0");
    assert_eq!(step.relation.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert_eq!(
        step.relation.text,
        "x D^2 Q[n+1] + D Q[n+1] - x D^2 Q[n] - D Q[n] - Q[n] = 0"
    );
    assert!(step
        .description
        .exceptional_index_conditions
        .iter()
        .any(|c| c.text == "[x^0] Q[n+1] - [x^0] Q[n] = 0"));
    assert_eq!(step.assessment.status, EQUIVALENT_DESCRIPTION);
    let best = &grandchild.descriptions[0];
    assert_eq!(best.label, EQUIVALENT_DESCRIPTION_LABEL);
    assert_eq!(best.provenance, "steps[0].implicit_relation.description");

    // A step whose recorded source differs from the parent's relation is rejected.
    let mut tampered = report.clone();
    let index = tampered
        .representations
        .iter()
        .position(|n| n.chain == [0, 0])
        .unwrap();
    tampered.representations[index].steps[0]
        .implicit_relation
        .as_mut()
        .unwrap()
        .source_relation
        .first_base = 1;
    assert!(verify_representation_search(&tampered).is_err());
}

// ---------------------------------------------------------------------------
// Vector representations, the rank-two gamma system and complete spaces
// after Möbius-class maps
// ---------------------------------------------------------------------------

fn eulerian_split_schedule() -> serde_json::Value {
    json!({"first_index": 1, "first_degree": 1, "step": 1})
}

/// Eulerian rows split as `E_n = a_n + x b_n` at `d_n = n`, then mapped
/// componentwise.
fn eulerian_vector_request() -> serde_json::Value {
    let schedule = eulerian_split_schedule();
    json!({
        "oeis_id": "A008292",
        "row_count": 10,
        "transformations": [
            {"kind": "symmetric_decomposition", "schedule": schedule},
            {"kind": "reversal", "schedule": schedule},
            {"kind": "f_polynomial", "schedule": schedule},
            {"kind": "gamma", "schedule": schedule}
        ],
        "limits": {"max_depth": 2, "max_representations": 32, "max_steps": 400}
    })
}

fn decode(states: &[Vec<Vec<String>>]) -> Vec<Vec<Vec<BigRational>>> {
    states
        .iter()
        .map(|s| s.iter().map(|c| c.iter().map(|x| q(x)).collect()).collect())
        .collect()
}

#[test]
fn vector_representations_are_mapped_transported_and_fitted() {
    let report = run(eulerian_vector_request());
    assert_eq!(report.status, "complete");
    let split = node_with_chain(&report, &[0]);
    assert_eq!(split.components, vec!["a", "b"]);
    assert!(split.expanded);
    let system = split
        .descriptions
        .iter()
        .find(|d| d.kind == "transported_vector_system")
        .unwrap();
    assert_eq!(system.label, EQUIVALENT_DESCRIPTION_LABEL);
    assert_eq!(system.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert_eq!(system.cost.component_count, 2);
    let uniqueness = system.vector_assessment.as_ref().unwrap();
    assert_eq!(uniqueness.status, "certified");
    assert_eq!(uniqueness.denominator_coefficient, Some(0));
    let fit = split
        .descriptions
        .iter()
        .find(|d| d.kind == "vector_fit")
        .unwrap();
    assert_eq!(fit.label, FINITE_FIT_LABEL);
    assert_eq!(split.fits.as_ref().unwrap().fits_found, 1);
    // Maps that do not act componentwise are rejected on vector states.
    let rejected = report
        .edges
        .iter()
        .find(|e| e.parent == split.id && e.transformation == 0)
        .unwrap();
    assert_eq!(rejected.outcome, "rejected");
    assert!(rejected.detail.as_ref().unwrap().contains("componentwise"));

    // f-polynomials of a_n (bound n) and b_n (bound n - 1): the certified
    // system is transported exactly, with certified degree conditions.
    let f = node_with_chain(&report, &[0, 2]);
    assert_eq!(f.components, vec!["a", "b"]);
    assert!(f.reconstruction_verified);
    assert_eq!(
        f.map.segments[1..],
        [
            "componentwise on a: Q_n(x) = (x + 1)^(n) P_n((x)/(x + 1))",
            "componentwise on b: Q_n(x) = (x + 1)^(n - 1) P_n((x)/(x + 1))"
        ]
    );
    assert!(f
        .map_conditions
        .iter()
        .all(|c| c.status == CERTIFIED_ALL_INDICES));
    let step = f.steps[0].vector_transport.as_ref().unwrap();
    assert_eq!(step.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert!(step.holds_on_rows);
    assert_eq!(step.component_maps.len(), 2);
    assert_eq!(
        step.text,
        "(1) a(n+1) = (1 + 2t) a(n) + (t + t^2) a'(n) + (2t + 2t^2) b(n); \
         (1) b(n+1) = (-1) a(n) + (t + t^2) b'(n)"
    );
    let transported = f
        .descriptions
        .iter()
        .find(|d| d.kind == "transported_vector_system")
        .unwrap();
    assert_eq!(transported.label, EQUIVALENT_DESCRIPTION_LABEL);
    // The transported system regenerates the states from the initial state.
    let (system, first, initial) = transported
        .vector
        .as_ref()
        .unwrap()
        .to_recurrence()
        .unwrap();
    let states = decode(&f.states);
    let mut state = initial;
    for (i, expected) in states.iter().enumerate().skip(1) {
        state = system
            .evaluate_next_rational(&state, first + i - 1)
            .unwrap();
        assert_eq!(&state, expected);
    }

    // Componentwise reversal at the symmetry indices fixes (a_n, b_n): equal
    // rows and equal certified systems, so the sequences are certified equal.
    let reversed = node_with_chain(&report, &[0, 1]);
    assert_eq!(reversed.states, split.states);
    assert!(!reversed.map.identity);
    assert_eq!(
        reversed.equal_sequence_as.as_ref().unwrap().representation,
        split.id
    );

    // Gamma vectors of the components: their symmetry is certified by the
    // decomposition; rows only and a vector fit.
    let gamma = node_with_chain(&report, &[0, 3]);
    assert!(gamma
        .map_conditions
        .iter()
        .all(|c| c.status == CERTIFIED_ALL_INDICES));
    assert!(gamma.steps.is_empty());
    assert!(gamma
        .descriptions
        .iter()
        .all(|d| d.label == FINITE_FIT_LABEL));
    assert!(gamma.descriptions.iter().any(|d| d.kind == "vector_fit"));
    assert!(report
        .ranked
        .iter()
        .any(|r| r.representation == f.id && r.cost.component_count == 2));
}

#[test]
fn vector_reports_replay_and_tampering_is_detected() {
    let report = run(eulerian_vector_request());
    let f = node_with_chain(&report, &[0, 2]).id;
    let back: RepresentationSearchReport =
        serde_json::from_str(&serde_json::to_string(&report).unwrap()).unwrap();
    assert_eq!(back, report);
    let system = |r: &mut RepresentationSearchReport| -> usize {
        r.representations[f]
            .descriptions
            .iter()
            .position(|d| d.kind == "transported_vector_system")
            .unwrap()
    };
    // A changed matrix entry of the transported system.
    let mut tampered = report.clone();
    let i = system(&mut tampered);
    tampered.representations[f].descriptions[i]
        .vector
        .as_mut()
        .unwrap()
        .matrix[1][0][0]
        .coeffs[0][0] = "2".into();
    assert!(verify_representation_search(&tampered).is_err());
    // The same change inside the replayable step.
    let mut tampered = report.clone();
    tampered.representations[f].steps[0]
        .vector_transport
        .as_mut()
        .unwrap()
        .transported
        .matrix[1][0][0]
        .coeffs[0][0] = "2".into();
    assert!(verify_representation_search(&tampered).is_err());
    // A forged uniqueness argument or evidence.
    let mut tampered = report.clone();
    let i = system(&mut tampered);
    tampered.representations[f].descriptions[i]
        .vector_assessment
        .as_mut()
        .unwrap()
        .status = "not_certified".into();
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report.clone();
    tampered.representations[f].steps[0]
        .vector_transport
        .as_mut()
        .unwrap()
        .evidence = CHECKED_SUPPLIED_ROWS.into();
    assert!(verify_representation_search(&tampered).is_err());
    // A changed component of a mapped state.
    let mut tampered = report.clone();
    tampered.representations[f].states[3][1][0] = "5".into();
    assert!(verify_representation_search(&tampered).is_err());
    // A finite fit promoted to an equivalent description.
    let mut tampered = report;
    let split = node_with_chain(&tampered, &[0]).id;
    let fit = tampered.representations[split]
        .descriptions
        .iter()
        .position(|d| d.kind == "vector_fit")
        .unwrap();
    tampered.representations[split].descriptions[fit].label = EQUIVALENT_DESCRIPTION_LABEL.into();
    assert!(verify_representation_search(&tampered).is_err());
}

#[test]
fn componentwise_maps_compose_and_deduplicate_by_normal_forms() {
    let schedule = eulerian_split_schedule();
    let shifted = json!({"first_index": 1, "first_degree": 2, "step": 1});
    let report = run(json!({
        "oeis_id": "A008292",
        "row_count": 8,
        "fits": false,
        "transformations": [
            {"kind": "symmetric_decomposition", "schedule": schedule},
            {"kind": "reversal", "schedule": shifted}
        ],
        "limits": {"max_depth": 3, "max_representations": 32, "max_steps": 400}
    }));
    // Reversal at d_n + 1 maps (a_n, b_n) to (x a_n, x b_n) componentwise.
    let once = node_with_chain(&report, &[0, 1]);
    assert_eq!(
        once.map.segments[1..],
        [
            "componentwise on a: Q_n(x) = x^(n + 1) P_n((1)/(x))",
            "componentwise on b: Q_n(x) = x^(n) P_n((1)/(x))"
        ]
    );
    let split = node_with_chain(&report, &[0]);
    for (mapped, original) in once.states.iter().zip(&split.states) {
        for (m, o) in mapped.iter().zip(original) {
            let o = o.iter().map(|c| q(c)).collect::<Vec<_>>();
            let shifted = std::iter::once(BigRational::from_integer(0.into()))
                .chain(o)
                .collect::<Vec<_>>();
            let m = m.iter().map(|c| q(c)).collect::<Vec<_>>();
            assert_eq!(
                polytool::vec_poly::trim(&m),
                polytool::vec_poly::trim(&shifted)
            );
        }
    }
    // The same reversal twice is the identity on each component: duplicate.
    let twice = report
        .edges
        .iter()
        .find(|e| e.parent == once.id && e.transformation == 1)
        .unwrap();
    assert_eq!(twice.outcome, "duplicate_map");
    assert_eq!(twice.representation, Some(split.id));
}

/// `(1 + x)^n` from `P_m = (1 + 2x - m x) P_(m-1) + (x + x^2) P'_(m-1)`, a
/// rule that is not invariant under reversal at `d_n = n`.
fn non_invariant_pascal_rule() -> serde_json::Value {
    json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none",
             "coeff": {"coeffs": [["1", "2"], ["0", "-1"]]}},
            {"offset": 1, "deriv_order": 1, "sign": "none",
             "coeff": {"coeffs": [["0", "1", "1"]]}}]}
    })
}

#[test]
fn rank_two_gamma_transport_is_a_two_relation_system_with_finite_evidence() {
    let report = run(json!({
        "source_recurrence": non_invariant_pascal_rule(),
        "transformations": [{"kind": "gamma",
                             "schedule": {"first_index": 0, "first_degree": 0, "step": 1}}],
        "limits": {"max_depth": 1}
    }));
    let node = node_with_chain(&report, &[0]);
    let transform = node.steps[0].transform.as_ref().unwrap();
    assert!(transform.transport_error.is_none());
    let transported = transform.transported.as_ref().unwrap();
    let system = transported.system.as_ref().unwrap();
    assert_eq!(system.relations.len(), 2);
    assert_eq!(system.relations[0].text, "Q[n+1] - Q[n] = 0");
    assert_eq!(system.relations[1].text, "-x D Q[n] = 0");
    assert_eq!(system.relations[1].role, ROW_EQUATION_ROLE);
    assert!(transported.holds_on_transformed_rows);
    // The symmetry of the rows is only checked: the system is a finite fit.
    assert_eq!(node.map_conditions[0].status, CHECKED_SUPPLIED_ROWS);
    let item = node
        .descriptions
        .iter()
        .find(|d| d.kind == "transported_gamma_system")
        .unwrap();
    assert_eq!(item.label, FINITE_FIT_LABEL);
    assert_eq!(item.evidence, CHECKED_SUPPLIED_ROWS);
    assert_eq!(item.cost.component_count, 2);
    let description = item.description.as_ref().unwrap();
    assert!(description
        .relations
        .iter()
        .all(|r| r.evidence == CHECKED_SUPPLIED_ROWS));
    assert_eq!(item.assessment.as_ref().unwrap().status, FINITE_EVIDENCE);
    // Claiming certified evidence without a certified symmetry is rejected.
    let mut tampered = report.clone();
    let i = tampered.representations[node.id]
        .descriptions
        .iter()
        .position(|d| d.kind == "transported_gamma_system")
        .unwrap();
    let forged = &mut tampered.representations[node.id].descriptions[i];
    forged.evidence = TRANSPORTED_IDENTITY_EVIDENCE.into();
    for relation in &mut forged.description.as_mut().unwrap().relations {
        relation.evidence = TRANSPORTED_IDENTITY_EVIDENCE.into();
    }
    forged.assessment =
        Some(assess_generating_description(forged.description.as_ref().unwrap()).unwrap());
    forged.label = EQUIVALENT_DESCRIPTION_LABEL.into();
    assert!(verify_representation_search(&tampered).is_err());
    // A changed relation of the embedded system does not replay.
    let mut tampered = report;
    tampered.representations[node.id].steps[0]
        .transform
        .as_mut()
        .unwrap()
        .transported
        .as_mut()
        .unwrap()
        .system
        .as_mut()
        .unwrap()
        .relations
        .swap(0, 1);
    assert!(verify_representation_search(&tampered).is_err());
}

/// Fibonacci polynomials `P_(n+2) = P_(n+1) + x P_n`, `P_0 = P_1 = 1`.
fn fibonacci_rule() -> serde_json::Value {
    json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"], ["1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1"]]}},
            {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"]]}}]}
    })
}

#[test]
fn complete_spaces_after_mobius_class_maps_and_weights() {
    let report = run(json!({
        "source_recurrence": fibonacci_rule(),
        "transformations": [
            {"kind": "reversal", "schedule": {"first_index": 0, "first_degree": 0, "step": 1}},
            {"kind": "diagonal", "weight": {"kind": "binomial", "m": 1}}
        ],
        "limits": {"max_depth": 1}
    }));
    for chain in [[0], [1]] {
        let node = node_with_chain(&report, &chain);
        assert!(!node.implicit_spaces.is_empty(), "{chain:?}");
        assert!(node.spaces.is_empty());
        let certified = node
            .descriptions
            .iter()
            .filter(|d| d.kind == "certified_space")
            .collect::<Vec<_>>();
        assert!(!certified.is_empty(), "{chain:?}: {:?}", node.notes);
        for d in certified {
            assert_eq!(d.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
            assert_ne!(d.label, FINITE_FIT_LABEL);
            assert!(d.provenance.starts_with("implicit_spaces["));
        }
    }
    // The reversed rows Q_(n+2) = x Q_(n+1) + x Q_n: certified and equivalent.
    let reversed = node_with_chain(&report, &[0]);
    let expected = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&["0", "1"]])),
            RecurrenceTerm::new(2, 0, bivar(&[&["0", "1"]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let found = find_rule(&reversed, &expected).unwrap();
    assert_eq!(found.label, EQUIVALENT_DESCRIPTION_LABEL);
    assert!(
        found.kind == "certified_space"
            || found
                .also_found_by
                .iter()
                .any(|p| p.starts_with("implicit_spaces")),
        "{found:?}"
    );
    let space = &reversed.implicit_spaces[0];
    assert!(space
        .source
        .completeness_basis
        .contains("transported-source"));
    // A space report with a lowered residual count, or from another source,
    // does not replay.
    let mut tampered = report.clone();
    tampered.representations[reversed.id].implicit_spaces[0].results[0].residual_count -= 1;
    assert!(verify_representation_search(&tampered).is_err());
    let mut tampered = report;
    tampered.representations[reversed.id].implicit_spaces[0]
        .request
        .source_recurrence = Some(serde_json::from_value(non_invariant_pascal_rule()).unwrap());
    assert!(verify_representation_search(&tampered).is_err());
}

#[test]
fn oversized_vector_fits_are_budget_events() {
    let report = run(json!({
        "oeis_id": "A008292",
        "row_count": 6,
        "transformations": [
            {"kind": "symmetric_decomposition", "schedule": eulerian_split_schedule()}
        ],
        "search": {"max_rec_len": 1, "max_var_deg": 6, "max_idx_deg": 6, "max_diff_deg": 2,
                   "max_candidates": 50},
        "limits": {"max_depth": 1}
    }));
    let split = node_with_chain(&report, &[0]);
    assert_eq!(split.fits.as_ref().unwrap().status, "skipped_budget");
    assert!(split.descriptions.iter().all(|d| d.kind != "vector_fit"));
    assert_eq!(report.status, "budget_exhausted");
    assert!(report
        .budget
        .iter()
        .any(|b| b.limit == "vector_fit_unknowns" && b.representation == Some(split.id)));
    // The certified transported system does not depend on the fit.
    assert!(split
        .descriptions
        .iter()
        .any(|d| d.kind == "transported_vector_system" && d.label == EQUIVALENT_DESCRIPTION_LABEL));
}

#[test]
fn replay_budget_exhaustion_keeps_its_variant() {
    // With 64 operations the search still certifies its implicit spaces, but
    // replaying their kernel certificates needs more arithmetic: the replay
    // reports the exhausted budget, not a failed verification.
    let report = search_representations(&request(json!({
        "source_recurrence": fibonacci_rule(),
        "transformations": [
            {"kind": "reversal", "schedule": {"first_index": 0, "first_degree": 0, "step": 1}}
        ],
        "limits": {"max_depth": 1, "max_operations": 64}
    })))
    .unwrap();
    let error = verify_representation_search(&report).unwrap_err();
    assert!(
        matches!(&error, RecurrenceError::BudgetExhausted(message)
            if message.contains("implicit_spaces[") && message.contains("budget")),
        "{error:?}"
    );
}

#[test]
fn too_many_space_shapes_exceed_a_limit() {
    let mut value = factorial_request();
    value["limits"]["max_space_shapes"] = json!(0);
    value["space_shapes"] = json!([{"kind": "box", "length": 1, "derivative_order": 0,
        "index_degree": 0, "variable_degree": 0, "leading_index_degree": 0}]);
    let error = search_representations(&request(value)).unwrap_err();
    assert_eq!(
        error,
        RecurrenceError::LimitExceeded("more space_shapes than max_space_shapes".into())
    );
}

#[test]
fn typed_errors_are_wrapped_with_context_not_flattened() {
    use std::error::Error;
    // A malformed transformation keeps its `TransformError` as the source.
    let mut value = factorial_request();
    value["transformations"] = json!([
        {"kind": "reversal", "schedule": {"first_index": 0, "first_degree": 0, "step": 0}},
        {"kind": "mobius", "a": "0", "b": "0", "c": "0", "e": "0",
         "schedule": {"first_index": 0, "first_degree": 0, "step": 1}}
    ]);
    let error = search_representations(&request(value)).unwrap_err();
    let rendered = error.to_string();
    assert!(rendered.starts_with("transformations["), "{rendered}");
    assert!(
        error
            .source()
            .is_some_and(|source| source.downcast_ref::<TransformError>().is_some()),
        "{error:?}"
    );
}

#[test]
fn an_unknown_recurrence_sign_is_invalid_input() {
    let mut rule = fibonacci_rule();
    rule["recurrence"]["terms"][0]["sign"] = json!("sideways");
    let json: RecurrenceJson = serde_json::from_value(rule).unwrap();
    assert_eq!(
        json.to_recurrence_parts().unwrap_err(),
        RecurrenceError::InvalidInput("unsupported recurrence sign `sideways`".into())
    );
}
