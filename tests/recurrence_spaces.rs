//! Spaces of recurrences, through the library API.

use polytool::recurrence::*;
use serde_json::{json, Value};

fn request(source: Value, shapes: Value) -> RecurrenceSpaceRequest {
    serde_json::from_value(json!({"source":source,"shapes":shapes})).unwrap()
}
fn box_shape(r: usize, s: usize, n: usize, x: usize) -> Value {
    json!({"kind":"box","length":r,"derivative_order":s,"index_degree":n,"variable_degree":x,"leading_index_degree":n})
}
fn run(req: &RecurrenceSpaceRequest) -> RecurrenceSpaces {
    let report = find_recurrence_spaces(req).unwrap();
    verify_recurrence_spaces(req, &report).unwrap();
    report
}

// Source snapshot: generated recurrence definitions of the OEIS families
// below (rows checked against the OEIS data).
// The canonical generator is ConstantLagSource::generate_rows; no floating
// root calculation or independently reimplemented family generator is used.
fn fixture(id: &str) -> Value {
    match id {
        "A049325" => json!({"coefficients":[["0","16"],["0","16"],["0","6"],["0","1"]],
            "seeds":[["16","12","1"],["16","68","18","1"],["0","224","156","24","1"],["0","448","840","280","30","1"]]}),
        "A171608" => {
            json!({"coefficients":[["0","0","-1"],["0"],["0","2"],["0"]],"seeds":[["1"],["2"],["0","2"],["0","3"]]})
        }
        "A111049" => {
            json!({"coefficients":[["0"],["-2","-4","-2"],["3","3"]],"seeds":[["1"],["1","1"],["1","3","2"]]})
        }
        _ => panic!("unknown fixture"),
    }
}

#[test]
fn complete_kernel_has_zero_leading_relations_and_successors() {
    let req = request(
        json!({"coefficients":[["0","1"]],"seeds":[["1"]]}),
        json!([box_shape(1, 1, 1, 1)]),
    );
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    assert_eq!(r.evidence, "conditional_on_source");
    assert!(r.certificate.as_ref().unwrap().basis.len() > r.successor_projection_rank);
    assert_eq!(r.h, 3);
}

#[test]
fn finite_rows_can_certify_absence_without_a_source() {
    let req: RecurrenceSpaceRequest =
        serde_json::from_value(json!({"rows":[["1"],["2"],["5"]],"shapes":[box_shape(1,0,0,0)]}))
            .unwrap();
    let report = run(&req);
    assert_eq!(report.results[0].status, "absent_in_shape");
    assert_eq!(report.results[0].evidence, "rows_only_obstruction");
    assert_eq!(
        report.results[0]
            .certificate
            .as_ref()
            .unwrap()
            .absence_duals
            .len(),
        1
    );
    let mut bad = report.clone();
    bad.results[0].certificate.as_mut().unwrap().absence_duals[0].weights[0] = "0".into();
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
}

#[test]
fn all_zero_matrix_preserves_its_column_count() {
    let req = request(
        json!({"coefficients":[["2"]],"seeds":[["0"]]}),
        json!([box_shape(1, 0, 0, 0)]),
    );
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.equations, 0);
    assert_eq!(r.certificate.as_ref().unwrap().basis.len(), 2);
    assert_eq!(r.successor_projection_rank, 1);
    assert!(r.all_checked_rows_zero);
}

#[test]
fn finite_support_and_source_with_zero_oldest_coefficient() {
    let req = request(
        json!({"coefficients":[["0"]],"seeds":[["1"]]}),
        json!([box_shape(1, 1, 1, 1)]),
    );
    let report = run(&req);
    assert_eq!(report.results[0].evidence, "conditional_on_source");
    assert_eq!(report.results[0].residual_count, 3);
}

#[test]
fn coefficient_indices_remain_original_on_a_tail() {
    let mut req = request(
        json!({"coefficients":[["1"]],"seeds":[["1"]]}),
        json!([{"kind":"explicit","target_lag":1,"monomials":[
            {"lag":1,"index_degree":2},{"lag":0}
        ]}]),
    );
    req.domain_start = 3;
    let report = run(&req);
    assert_eq!(report.results[0].status, "absent_in_shape");
    assert_eq!(report.results[0].last_required_source_row, Some(6));
    let mut changed = req.clone();
    changed.domain_start = 4;
    assert!(verify_recurrence_spaces(&changed, &report).is_err());
}

#[test]
fn leading_degree_contributes_to_the_finite_bound() {
    let req = request(
        json!({"coefficients":[["1"]],"seeds":[["1"]]}),
        json!([{"kind":"explicit","target_lag":1,"monomials":[
            {"lag":1,"index_degree":0},{"lag":1,"index_degree":1},{"lag":1,"index_degree":2}
        ]}]),
    );
    let report = run(&req);
    assert_eq!(report.results[0].h, 3);
    assert_eq!(report.results[0].status, "absent_in_shape");
    // n(n-1) fools the first two ranks, but fails at the required last rank.
    let mut finite = req.clone();
    finite.source = None;
    finite.rows = Some(vec![vec!["1".into()]; 3]);
    assert_eq!(run(&finite).results[0].evidence, "exact_fit");
    let mut bad = report.clone();
    bad.results[0].residual_count -= 1;
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
}

#[test]
fn singular_rank_outside_checked_prefix_is_reported() {
    // P[n]=n-51: n(n+1) shifts are generated by a constant order-two rule.
    // Its one-row leader n-51 must vanish well beyond the 4 checked residuals.
    let req = request(
        json!({"coefficients":[["-1"],["2"]],"seeds":[["-51"],["-50"]]}),
        json!([box_shape(1, 0, 1, 0)]),
    );
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists");
    assert!(r.representatives.iter().any(|p| p
        .singular_ranks
        .as_ref()
        .unwrap()
        .contains(&"51".into())));
}

#[test]
fn mutation_of_basis_minor_and_source_is_rejected() {
    let req = request(
        json!({"coefficients":[["2"]],"seeds":[["1"]]}),
        json!([box_shape(1, 0, 0, 0)]),
    );
    let report = run(&req);
    let mut bad = report.clone();
    bad.results[0].certificate.as_mut().unwrap().basis[0][0] = "999".into();
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
    let mut bad = report.clone();
    bad.results[0].certificate.as_mut().unwrap().inverse_minor[0][0] = "0".into();
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
    let mut bad = report.clone();
    bad.request.source.as_mut().unwrap().seeds[0][0] = "2".into();
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
}

#[test]
fn budgets_and_missing_rows_never_mean_nonexistence() {
    let mut req = request(
        json!({"coefficients":[["2"]],"seeds":[["1"]]}),
        json!([box_shape(1, 0, 0, 0)]),
    );
    req.limits.max_operations = 0;
    let report = run(&req);
    assert_eq!(report.results[0].status, "incomplete_budget");
    req.limits = Default::default();
    req.limits.max_columns = 1;
    assert_eq!(run(&req).results[0].status, "incomplete_budget");
    let finite: RecurrenceSpaceRequest =
        serde_json::from_value(json!({"rows":[["1"]],"shapes":[box_shape(1,0,0,0)]})).unwrap();
    assert_eq!(run(&finite).results[0].status, "insufficient_rows");
}

#[test]
fn representations_replay_in_the_independent_legacy_evaluator() {
    for (id, shape) in [
        ("A049325", box_shape(1, 1, 1, 2)),
        ("A171608", box_shape(1, 2, 1, 3)),
    ] {
        let req = request(fixture(id), json!([shape]));
        let report = run(&req);
        let rows = req
            .source
            .as_ref()
            .unwrap()
            .generate_rows(30, &req.limits)
            .unwrap();
        let rows: Vec<Vec<BigRational>> = rows
            .iter()
            .map(|r| r.iter().map(|s| parse_rational_coeff(s).unwrap()).collect())
            .collect();
        let result = &report.results[0];
        for rep in &result.representatives {
            let rec = rep.to_recurrence(result).unwrap();
            for target in result.target_lag..rows.len() {
                let n = target - result.target_lag;
                if rep
                    .singular_ranks
                    .as_ref()
                    .unwrap()
                    .contains(&n.to_string())
                {
                    continue;
                }
                let calculated = rec.evaluate_next_rational(&rows[..target], 0).unwrap();
                assert_eq!(
                    polytool::Polynomial::new(calculated),
                    polytool::Polynomial::new(rows[target].clone())
                );
            }
        }
    }
}

#[test]
fn pinned_oeis_menus_and_initial_boundary() {
    let req = request(
        fixture("A171608"),
        json!([
            box_shape(3, 0, 1, 3),
            box_shape(2, 1, 1, 3),
            box_shape(1, 2, 1, 3),
            box_shape(1, 1, 1, 3),
            box_shape(1, 1, 2, 3)
        ]),
    );
    let report = run(&req);
    for &i in &[0, 1, 2, 4] {
        assert_eq!(
            report.results[i].status, "successor_exists",
            "{:?}",
            report.results[i]
        );
    }
    assert_eq!(report.results[3].status, "absent_in_shape");
    assert_eq!(report.frontier, vec![0, 1, 2, 4]);
    let req = request(fixture("A049325"), json!([box_shape(1, 1, 1, 2)]));
    assert_eq!(run(&req).results[0].status, "successor_exists");
    let mut req = request(fixture("A111049"), json!([box_shape(1, 1, 1, 2)]));
    let report = run(&req);
    assert!(report.results[0].representatives.iter().all(|p| p
        .singular_ranks
        .as_ref()
        .unwrap()
        .contains(&"0".into())));
    req.domain_start = 1;
    let report = run(&req);
    assert!(report.results[0].representatives.iter().any(|p| !p
        .singular_ranks
        .as_ref()
        .unwrap()
        .contains(&"0".into())));
}

#[test]
fn source_adapter_preserves_origin_and_rejects_inconsistent_extra_seeds() {
    let rec = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![BigRational::from_integer(
                2.into(),
            )]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let json = RecurrenceJson::from_recurrence_rational(
        &rec,
        7,
        &[vec![BigRational::from_integer(3.into())]],
        None,
    );
    let req: RecurrenceSpaceRequest = serde_json::from_value(
        json!({"source_recurrence":json,"domain_start":7,"shapes":[box_shape(1,0,0,0)]}),
    )
    .unwrap();
    assert_eq!(run(&req).results[0].status, "successor_exists");
    let mut bad = req.clone();
    bad.source_recurrence
        .as_mut()
        .unwrap()
        .initial_polynomials
        .push(vec!["7".into()]);
    assert!(find_recurrence_spaces(&bad)
        .unwrap_err()
        .contains("disagree"));
    // Additional supplied rows must not mask contradictory initial data.
    bad.rows = Some(vec![vec!["3".into()]]);
    bad.row_origin = 7;
    assert!(find_recurrence_spaces(&bad)
        .unwrap_err()
        .contains("disagree"));
}

#[test]
fn nonminimal_source_has_the_same_complete_kernel() {
    let a = request(
        json!({"coefficients":[["2"]],"seeds":[["1"]]}),
        json!([box_shape(1, 1, 1, 1)]),
    );
    let b = request(
        json!({"coefficients":[["-2"],["3"]],"seeds":[["1"],["2"]]}),
        json!([box_shape(1, 1, 1, 1)]),
    );
    assert_eq!(
        run(&a).results[0].certificate.as_ref().unwrap().basis,
        run(&b).results[0].certificate.as_ref().unwrap().basis
    );
}

#[test]
fn canonical_source_generator_matches_pinned_and_bundled_rows() {
    let source: ConstantLagSource = serde_json::from_value(fixture("A171608")).unwrap();
    let rows = source.generate_rows(12, &Default::default()).unwrap();
    let bundled = polytool::oeis::by_id("A171608")
        .unwrap()
        .generate_rows(12)
        .unwrap();
    for (a, b) in rows.iter().zip(bundled) {
        let mut b: Vec<_> = b.into_iter().map(|x| x.to_string()).collect();
        while b.len() > 1 && b.last().unwrap() == "0" {
            b.pop();
        }
        assert_eq!(a, &b);
    }
}

#[test]
fn automatic_variable_degree_preserves_the_leader_space() {
    let mut req = request(
        json!({"coefficients":[["0","1"]],"seeds":[["1"]]}),
        json!([{"kind":"automatic_variable_degree","length":1,"derivative_order":1,"index_degree":1,"leading_index_degree":1}]),
    );
    let report = run(&req);
    let result = &report.results[0];
    assert_eq!(result.status, "successor_exists");
    assert!(
        result
            .automatic_degree_bound
            .as_ref()
            .unwrap()
            .variable_degree
            >= 1
    );
    let mut tampered = report.clone();
    tampered.results[0]
        .automatic_degree_bound
        .as_mut()
        .unwrap()
        .variable_degree += 1;
    assert!(verify_recurrence_spaces(&req, &tampered).is_err());
    req.source = None;
    req.rows = Some(vec![vec!["1".into()]; 8]);
    assert!(find_recurrence_spaces(&req)
        .unwrap_err()
        .contains("defining source"));
}

#[test]
fn leader_combination_removes_each_avoidable_singularity() {
    // A constant sequence has independent leaders 1 and n. The combined
    // leader 1+n^2 has no real zero and preserves the whole coverage.
    let mut req = request(
        json!({"coefficients":[["1"]],"seeds":[["1"]]}),
        json!([box_shape(1, 0, 1, 0)]),
    );
    req.include_leader_coverage = true;
    let report = run(&req);
    let r = &report.results[0];
    let coverage = r.leader_coverage.as_ref().unwrap();
    assert_eq!(coverage.common_singular_ranks, Some(vec![]));
    assert_eq!(coverage.basis_indices.len(), r.successor_projection_rank);
    assert!(coverage.monomials.iter().any(|c| c.index_degree == 2));
    let mut bad = report.clone();
    bad.results[0]
        .leader_coverage
        .as_mut()
        .unwrap()
        .coefficients[0] = "999".into();
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
    let rec = coverage.to_recurrence(r.target_lag).unwrap();
    let one = vec![BigRational::from_integer(1.into())];
    for n in 1..8 {
        let rows = vec![one.clone(); n];
        assert_eq!(rec.evaluate_next_rational(&rows, 0).unwrap(), one);
    }
}

#[test]
fn leader_combination_reports_an_unavoidable_rank() {
    let mut req = request(
        json!({"coefficients":[["-1"],["2"]],"seeds":[["-51"],["-50"]]}),
        json!([box_shape(1, 0, 1, 0)]),
    );
    req.include_leader_coverage = true;
    let report = run(&req);
    assert_eq!(
        report.results[0]
            .leader_coverage
            .as_ref()
            .unwrap()
            .common_singular_ranks,
        Some(vec!["51".into()])
    );
}

#[test]
fn polynomial_index_source_has_a_replayable_common_propagator() {
    // Canonical generator: IndexedLagSource::generate_rows for P[n]=n! X^n.
    let req:RecurrenceSpaceRequest=serde_json::from_value(json!({
        "indexed_source":{"coefficients":[[["0","1"],["0","1"]]],"seeds":[["1"]]},
        "shapes":[{"kind":"box","length":1,"derivative_order":1,"index_degree":1,"variable_degree":1,"leading_index_degree":1}]
    })).unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists", "{:?}", r.reason);
    assert_eq!(r.evidence, "conditional_on_source");
    assert!(r.propagation.is_some());
    let rows = req
        .indexed_source
        .as_ref()
        .unwrap()
        .generate_rows(6, &Default::default())
        .unwrap();
    assert_eq!(rows[5], vec!["0", "0", "0", "0", "0", "120"]);
    let mut bad = report.clone();
    bad.results[0].propagation.as_mut().unwrap().coefficients[0] = vec![vec!["123".into()]];
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
    let mut bad = report.clone();
    bad.results[0].propagation.as_mut().unwrap().regular_from += 1;
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
}

#[test]
fn indexed_source_derivative_lift_for_hermite() {
    // H[n+2]=2X H[n+1]-2(n+1) H[n], seeds H[0]=1,H[1]=2X.
    let req: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "indexed_source":{"coefficients":[[["-2"],["-2"]],[["0","2"]]],"seeds":[["1"],["0","2"]]},
        "shapes":[{"kind":"box","length":1,"derivative_order":1,"variable_degree":1}]
    }))
    .unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists", "{:?}", r.reason);
    assert_eq!(r.successor_projection_rank, 1);
    let generated = req
        .indexed_source
        .as_ref()
        .unwrap()
        .generate_rows(5, &Default::default())
        .unwrap();
    assert_eq!(generated[4], vec!["12", "0", "-48", "0", "16"]);
    let rec = r.representatives[0].to_recurrence(r).unwrap();
    let exact: Vec<Vec<BigRational>> = generated
        .iter()
        .map(|row| {
            row.iter()
                .map(|c| parse_rational_coeff(c).unwrap())
                .collect()
        })
        .collect();
    assert_eq!(
        rec.evaluate_next_rational(&exact[..4], 0).unwrap(),
        exact[4]
    );
}

#[test]
fn indexed_propagation_checks_beyond_a_late_exception() {
    let req:RecurrenceSpaceRequest=serde_json::from_value(json!({
        "indexed_source":{"coefficients":[[["0","-17"],["0","1"]]],"seeds":[["1"]]},
        "shapes":[{"kind":"box","length":1,"index_degree":1,"variable_degree":1,"leading_index_degree":1}]
    })).unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists", "{:?}", r.reason);
    assert!(r.residual_count > r.h);
    assert!(r.propagation.as_ref().unwrap().regular_from > 17);
    let mut bad = report.clone();
    bad.results[0].residual_count = r.h;
    assert!(verify_recurrence_spaces(&req, &bad).is_err());
}

#[test]
fn automatic_bound_also_applies_to_indexed_sources() {
    let req: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "indexed_source":{"coefficients":[[["0","1"],["0","1"]]],"seeds":[["1"]]},
        "shapes":[{"kind":"automatic_variable_degree","length":1,"index_degree":1}]
    }))
    .unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists", "{:?}", r.reason);
    assert!(r.automatic_degree_bound.is_some() && r.propagation.is_some());
}

#[test]
fn automatic_degree_can_prove_absence_for_all_rhs_variable_degrees() {
    // P[n]=1+X^n has two active bases, hence no constant-index one-row rule.
    let req = request(
        json!({"coefficients":[["0","-1"],["1","1"]],"seeds":[["2"],["1","1"]]}),
        json!([{"kind":"automatic_variable_degree","length":1}]),
    );
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "absent_in_shape");
    assert_eq!(
        r.automatic_degree_bound.as_ref().unwrap().variable_degree,
        2
    );
    assert!(!r.certificate.as_ref().unwrap().absence_duals.is_empty());
}

#[test]
fn automatic_degree_finds_a_rule_missed_by_zero_variable_degree() {
    let mut req = request(
        json!({"coefficients":[["0","1"]],"seeds":[["1"]]}),
        json!([box_shape(1, 0, 0, 0)]),
    );
    assert_eq!(run(&req).results[0].status, "absent_in_shape");
    req.shapes =
        serde_json::from_value(json!([{"kind":"automatic_variable_degree","length":1}])).unwrap();
    assert_eq!(run(&req).results[0].status, "successor_exists");
    req.limits.max_columns = 2;
    let partial = run(&req);
    assert_eq!(partial.results[0].status, "incomplete_budget");
    assert!(partial.results[0]
        .reason
        .as_ref()
        .unwrap()
        .contains("automatic RHS X-degree 1"));
}

#[test]
fn propagation_uses_a_proved_root_bound_when_exact_root_analysis_is_unavailable() {
    let mut coefficient = vec![vec!["0"]; 10];
    coefficient[9] = vec!["1"];
    let req: RecurrenceSpaceRequest = serde_json::from_value(json!({
        "indexed_source":{"coefficients":[coefficient],"seeds":[["1"]]},
        "shapes":[{"kind":"explicit","target_lag":1,"monomials":[{"lag":1}]}]
    }))
    .unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "successor_exists", "{:?}", r.reason);
    let cert = r.propagation.as_ref().unwrap();
    assert_eq!(cert.cutoff_method, "coefficient_root_bound");
    assert_eq!(cert.exceptional_ranks, None);
    assert_eq!(cert.regular_from, 2);
}

#[test]
fn indexed_zero_order_and_original_origin_are_checked() {
    let req:RecurrenceSpaceRequest=serde_json::from_value(json!({
        "indexed_source":{"coefficients":[[["1"],["1"]]],"seeds":[["2"]],"origin":7},"domain_start":7,
        "shapes":[{"kind":"explicit","target_lag":1,"monomials":[]}]
    })).unwrap();
    let report = run(&req);
    let r = &report.results[0];
    assert_eq!(r.status, "absent_in_shape");
    assert_eq!(r.residual_count, 0);
    assert_eq!(r.propagation.as_ref().unwrap().coefficients.len(), 1);
    assert_eq!(
        req.indexed_source
            .as_ref()
            .unwrap()
            .generate_rows(3, &Default::default())
            .unwrap(),
        vec![vec!["2"], vec!["16"], vec!["144"]]
    );
}
