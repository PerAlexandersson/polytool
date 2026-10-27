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
