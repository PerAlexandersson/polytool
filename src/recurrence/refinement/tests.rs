//! Refinement discovery and exact rejection regressions.
use super::*;

fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn biv(rows: &[&[i64]]) -> BivarPolyJson {
    BivarPolyJson {
        coeffs: rows
            .iter()
            .map(|r| r.iter().map(ToString::to_string).collect())
            .collect(),
    }
}
fn request(source: TransformVectorRecurrenceJson) -> InterlacingRefinementRequest {
    serde_json::from_value(serde_json::json!({"source":source,"derive_candidates":false})).unwrap()
}
fn supplied(source: TransformVectorRecurrenceJson) -> InterlacingRefinementRequest {
    let mut r = request(source.clone());
    r.candidates.push(InterlacingRefinementCandidate {
        name: "supplied".into(),
        coordinates: None,
        transition: Some(source.clone()),
        projection: (0..source.matrix.len())
            .map(|i| usize::from(i == source.matrix.len() - 1).to_string())
            .collect(),
    });
    r
}
fn system(m: OperatorMatrix, seed: Vec<Vec<Q>>, first: usize) -> TransformVectorRecurrenceJson {
    encode_system(&m, &Biv::one(), &seed, first)
}
fn scalar(a: i64, b: i64, d: i64, e: i64) -> TransformVectorRecurrenceJson {
    system(
        vec![
            vec![Op::default(), Op::identity(), Op::default()],
            vec![Op::default(), Op::default(), Op::identity()],
            vec![
                Op::single(0, 0, Biv::new(vec![vec![q(0), q(0), q(d)]])),
                Op::single(0, 0, Biv::new(vec![vec![q(0), q(e)]])),
                Op::single(0, 0, Biv::new(vec![vec![q(a), q(b)]])),
            ],
        ],
        vec![vec![q(1)], vec![q(1)], vec![q(1), q(2)]],
        0,
    )
}
#[test]
fn discovers_the_sharp_scaled_three_state_template() {
    let mut r = request(scalar(1, 1, 4, 4));
    r.derive_candidates = true;
    r.limits.maximum_candidates = 2;
    let report = discover_interlacing_refinements(&r).unwrap();
    let interval = report.scaled_interval.as_ref().unwrap();
    assert!(interval.feasible);
    assert_eq!(interval.rational_scale.as_deref(), Some("2"));
    assert_eq!(interval.lower, interval.upper);
    assert_eq!(report.status, "proved", "{report:?}");
    verify_interlacing_refinements(&report).unwrap();
    let mut bad = report.clone();
    bad.candidates[0].identity = "finite_fit".into();
    assert!(verify_interlacing_refinements(&bad).is_err());
    let good = report.candidates[0].candidate.clone();
    let mut bad = good.clone();
    let t = bad.coordinates.as_mut().unwrap();
    t[1][0][0].coeffs[0][1] = "1".into();
    r.derive_candidates = false;
    r.candidates = vec![bad];
    let rejected = discover_interlacing_refinements(&r).unwrap();
    assert_ne!(rejected.status, "proved");
}
#[test]
fn scale_solver_has_all_zero_and_irrational_branches() {
    for (a, b, d, e, feasible, scale) in [
        (0, 0, 0, 0, true, Some("1")),
        (0, 0, 1, 0, false, None),
        (0, 1, 2, 1, true, Some("3")),
        (1, 0, 0, 2, true, Some("1")),
        (1, 1, 4, 3, false, None),
        (0, 0, 0, -1, false, None),
        (-1, 1, 1, 4, false, None),
    ] {
        let p = [a, b, d, e].map(|n| n.to_string());
        let s = solve_scaled_refinement_interval(&p).unwrap();
        assert_eq!(s.feasible, feasible, "{p:?}");
        assert_eq!(s.rational_scale.as_deref(), scale);
    }
    let s = solve_scaled_refinement_interval(&["1", "1", "1", "3"].map(str::to_string)).unwrap();
    assert!(s.feasible);
    assert_eq!(s.lower.unwrap().radicand(), &q(5));
}
#[test]
fn rejects_exactly_the_two_smooth_descents_blocks_at_six() {
    let matrix = vec![
        vec![biv(&[&[1]]); 3],
        vec![biv(&[&[0, 1031]]), biv(&[&[1249]]), biv(&[&[989]])],
        vec![biv(&[&[0, 491]]), biv(&[&[185]]), biv(&[&[1011]])],
    ];
    let mut matrix = matrix;
    for (i, j, divisor) in [
        (1, 0, 331),
        (1, 1, 1324),
        (1, 2, 1324),
        (2, 0, 331),
        (2, 1, 3972),
        (2, 2, 1324),
    ] {
        for row in &mut matrix[i][j].coeffs {
            for c in row {
                *c = (parse_rational_coeff(c).unwrap() / q(divisor)).to_string();
            }
        }
    }
    let check = check_affine_interlacing_preserver(&matrix, 6).unwrap();
    assert_eq!(check.status, "rejected");
    let pairs: Vec<_> = check
        .rejections
        .iter()
        .map(|w| (w.rows, w.columns))
        .collect();
    assert_eq!(pairs, vec![([0, 1], [1, 2]), ([1, 2], [0, 2])]);
    assert!(check.rejections.iter().all(|w| w.index == 6));
}
#[test]
fn repeated_indices_and_scalar_real_rootedness_are_not_vacuous() {
    let check = check_affine_interlacing_preserver(&[vec![biv(&[&[1, 1, 1]])]], 0).unwrap();
    assert_eq!(check.status, "rejected");
    let check =
        check_affine_interlacing_preserver(&[vec![biv(&[&[0, 1]]), biv(&[&[1]])]], 0).unwrap();
    assert_eq!(check.status, "proved");
    assert_eq!(check.pencils.len(), 3);
}
#[test]
fn symbolic_index_pencils_and_derivative_baseline() {
    // [[1,1],[4*n*x,n]] preserves for every n>=1.
    let check = check_affine_interlacing_preserver(
        &[
            vec![biv(&[&[1]]), biv(&[&[1]])],
            vec![biv(&[&[0], &[0, 4]]), biv(&[&[0], &[1]])],
        ],
        1,
    )
    .unwrap();
    assert_eq!(check.status, "proved", "{check:?}");
    assert_eq!(check.pencils.len(), 9);
    let entry = crate::oeis::by_id("A124324").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let record =
        RecurrenceJson::from_recurrence_rational(&rule, entry.recurrence_first_index, &seeds, None);
    let mut r: InterlacingRefinementRequest =
        serde_json::from_value(serde_json::json!({"recurrence":record})).unwrap();
    r.limits.maximum_candidates = 1;
    let report = discover_interlacing_refinements(&r).unwrap();
    assert_eq!(report.status, "proved", "{report:?}");
    assert!(report.candidates[0]
        .preservation
        .as_ref()
        .unwrap()
        .reason
        .contains("Euler"));
    verify_interlacing_refinements(&report).unwrap();
}

#[test]
fn finite_refinement_fit_never_becomes_an_all_index_proof() {
    let mut state = vec![vec![q(1)], vec![q(1)]];
    let matrix = vec![
        vec![Op::identity(), Op::identity()],
        vec![
            Op::single(0, 0, Biv::new(vec![vec![q(0), q(1)]])),
            Op::identity(),
        ],
    ];
    let mut states = Vec::new();
    let mut targets = Vec::new();
    for i in 0..7 {
        states.push(state.iter().map(|p| strings(p)).collect::<Vec<_>>());
        targets.push(strings(&vp::sum(&state)));
        state = apply(&matrix, &state, i);
    }
    let mut r:InterlacingRefinementRequest=serde_json::from_value(serde_json::json!({"states":states,"target_rows":targets,"derive_candidates":false,"candidates":[{"name":"finite vectors"}]})).unwrap();
    r.limits.index_degree = 0;
    let report = discover_interlacing_refinements(&r).unwrap();
    assert_eq!(report.status, "finite_fit");
    assert_eq!(report.candidates[0].identity, "finite_fit");
    assert_eq!(report.candidates[0].target_recovery, "finite_fit");
    assert_eq!(
        report.candidates[0].preservation.as_ref().unwrap().status,
        "proved"
    );
    verify_interlacing_refinements(&report).unwrap();
}
#[test]
fn eulerian_last_letter_refinement_has_rectangular_preserving_transfers() {
    let rows = crate::sequences::refined_eulerian_polynomials_bigint(5);
    let targets = crate::sequences::eulerian_polynomials_bigint(5);
    for n in 1..5 {
        let matrix: Vec<Vec<_>> = (0..=n)
            .map(|j| {
                (0..n)
                    .map(|k| if k < j { biv(&[&[0, 1]]) } else { biv(&[&[1]]) })
                    .collect()
            })
            .collect();
        let certificate = check_affine_interlacing_preserver(&matrix, n).unwrap();
        assert_eq!(certificate.status, "proved", "{n}: {certificate:?}");
        let operators = matrix
            .iter()
            .map(|r| {
                r.iter()
                    .map(|p| Op::single(0, 0, decode_biv(p).unwrap()))
                    .collect()
            })
            .collect();
        let input: Vec<_> = rows[n]
            .iter()
            .map(|p| p.iter().cloned().map(Q::from_integer).collect())
            .collect();
        let output = apply(&operators, &input, n);
        let expected: Vec<_> = rows[n + 1]
            .iter()
            .map(|p| p.iter().cloned().map(Q::from_integer).collect::<Vec<_>>())
            .collect();
        assert_eq!(
            output.iter().map(|p| vp::trim(p)).collect::<Vec<_>>(),
            expected.iter().map(|p| vp::trim(p)).collect::<Vec<_>>()
        );
        assert_eq!(vp::clear_denominators(&vp::sum(&output)), targets[n + 1]);
    }
}
#[test]
fn known_positive_rows_are_fitted_without_unsupported_preservation_claims() {
    for rows in [
        crate::sequences::narayana_polynomials_bigint(8),
        crate::sequences::derangement_descent_polynomials_bigint(8),
    ] {
        let states: Vec<_> = rows
            .iter()
            .skip(2)
            .map(|p| vec![p.iter().map(ToString::to_string).collect::<Vec<_>>()])
            .collect();
        let mut r:InterlacingRefinementRequest=serde_json::from_value(serde_json::json!({"states":states,"first_index":2,"candidates":[{"name":"one coordinate"}],"derive_candidates":false})).unwrap();
        r.limits.variable_degree = 2;
        r.limits.index_degree = 2;
        let report = discover_interlacing_refinements(&r).unwrap();
        assert_ne!(report.status, "proved");
        verify_interlacing_refinements(&report).unwrap();
    }
}
#[test]
fn stirling_bell_derivative_transition_has_exact_identity() {
    // Bell recurrence B_(n+1)=x*B_n+x*B_n' has a general PF-preserver.
    let m = vec![vec![Op::single(0, 0, Biv::new(vec![vec![q(0), q(1)]]))
        .add(&Op::single(0, 1, Biv::new(vec![vec![q(0), q(1)]])))]];
    let r = supplied(system(m, vec![vec![q(1)]], 0));
    let report = discover_interlacing_refinements(&r).unwrap();
    assert_eq!(report.candidates[0].identity, "exact");
    assert_eq!(report.candidates[0].target_recovery, "exact");
    assert_eq!(report.status, "proved");
}
#[test]
fn weyl_intertwining_uses_the_product_rule() {
    let d = Op::single(0, 1, Biv::one());
    let x = Op::single(0, 0, Biv::new(vec![vec![q(0), q(1)]]));
    assert_eq!(d.compose(&x), x.compose(&d).add(&Op::identity()));
    // T=D, C=X; fitting M=D X D^-1 would require a nonpolynomial term.
    let fit = fit_symbolic(
        &vec![vec![d]],
        &vec![vec![x]],
        &RefinementSearchLimits::default(),
    )
    .unwrap();
    assert!(fit.is_none());
}
#[test]
fn rejects_vacuous_target_and_enforces_limits() {
    let mut r = supplied(system(identity(1), vec![vec![q(1)]], 0));
    r.candidates[0].projection = vec!["0".into()];
    assert!(discover_interlacing_refinements(&r).is_err());
    r.limits.maximum_candidates = 0;
    assert!(discover_interlacing_refinements(&r).is_err());
}

#[test]
fn big_block_normalized_and_derivative_state_rejections_keep_exact_witnesses() {
    let theta =
        Op::single(0, 0, Biv::one()).add(&Op::single(0, 1, Biv::new(vec![vec![q(0), q(1)]])));
    // At source n=6, the normalized companion produces
    // 232+1828*x+105*x^2 on the ordered input listed in the issue.
    let matrix = vec![
        vec![Op::default(), Op::identity(), Op::default()],
        vec![Op::default(), Op::default(), Op::identity()],
        vec![
            Op::single(
                0,
                0,
                Biv::new(vec![
                    vec![q(0)],
                    vec![q(0), Q::new((-1).into(), 2.into())],
                    vec![q(0), Q::new(1.into(), 2.into())],
                ]),
            ),
            Op::single(0, 0, Biv::new(vec![vec![q(0)], vec![q(1)]])),
            theta.clone(),
        ],
    ];
    let report = preservation::check_operator_preserver(&matrix, 6).unwrap();
    assert_eq!(report.status, "rejected");
    assert!(report
        .operator_rejections
        .iter()
        .any(|w| w.output[2] == ["232", "1828", "105"].map(str::to_string)));
    let derivative = vec![
        vec![Op::default(); 3],
        vec![Op::default(), Op::default(), Op::identity()],
        vec![
            Op::default(),
            Op::single(0, 0, Biv::new(vec![vec![q(0)], vec![q(1)]])),
            theta,
        ],
    ];
    let report = preservation::check_operator_preserver(&derivative, 1).unwrap();
    assert_eq!(report.status, "rejected");
    assert!(!report.operator_rejections.is_empty());
}
#[test]
fn ordered_degree_phase_guard_rejects_a_lag_sum_state() {
    let seed = vec![vec![q(1)], vec![q(1), q(1)], vec![q(1)]];
    let report = discover_interlacing_refinements(&supplied(system(identity(3), seed, 2))).unwrap();
    assert!(!report.candidates[0].degree_conditions);
    assert!(!report.candidates[0].initial_invariant);
    assert_eq!(report.candidates[0].status, "rejected");
}

#[test]
fn singular_source_denominator_cannot_be_cancelled_into_a_proof() {
    let nminus1 = Biv::new(vec![vec![q(-1)], vec![q(1)]]);
    let mut source = system(
        vec![vec![Op::single(0, 0, nminus1.clone())]],
        vec![vec![q(1)]],
        0,
    );
    source.denominator = encode_biv(&nminus1);
    let mut r = supplied(source);
    r.candidates[0].transition = Some(system(identity(1), vec![vec![q(1)]], 0));
    let report = discover_interlacing_refinements(&r).unwrap();
    assert_eq!(report.candidates[0].identity, "exact");
    assert_ne!(report.status, "proved");
}

#[test]
fn varying_steps_require_one_common_scale() {
    let step = |a: i64, b: i64, d: i64, e: i64| [a, b, d, e].map(|n| n.to_string());
    assert!(
        solve_common_scaled_refinement_interval(&[step(1, 1, 4, 4), step(1, 1, 1, 2)])
            .unwrap()
            .is_none()
    );
    let interval = solve_common_scaled_refinement_interval(&[step(1, 1, 4, 4), step(0, 1, 1, 1)])
        .unwrap()
        .unwrap();
    assert_eq!(interval.lower.unwrap().value.as_rational(), Some(&q(2)));
    assert_eq!(interval.upper.unwrap().value.as_rational(), Some(&q(2)));
}

#[test]
fn four_coordinate_nonnegative_transition_can_break_real_rootedness() {
    let input = vec![
        vec![q(8), q(4)],
        vec![q(4), q(4)],
        vec![q(0), q(4)],
        vec![q(0), q(1)],
    ];
    assert!((0..4).all(|i| (i + 1..4).all(|j| preservation::directed(&input[i], &input[j]))));
    let mut matrix = identity(4);
    matrix[3] = vec![
        Op::identity(),
        Op::single(0, 0, Biv::new(vec![vec![q(0), q(1)]])),
        Op::default(),
        Op::identity().scale(&q(2)),
    ];
    let output = apply(&matrix, &input, 0);
    assert_eq!(output[3], vec![q(8), q(10), q(4)]);
    assert!(!crate::is_real_rooted_bigint_coeffs(
        &vp::clear_denominators(&output[3])
    ));
    assert_eq!(
        preservation::check_operator_preserver(&matrix, 0)
            .unwrap()
            .status,
        "rejected"
    );
}

#[test]
fn stored_transition_replay_does_not_rerun_a_successful_fit() {
    let source = system(
        vec![vec![Op::single(0, 0, Biv::new(vec![vec![q(1), q(1)]]))]],
        vec![vec![q(1)]],
        0,
    );
    let mut input = request(source);
    input.candidates.push(InterlacingRefinementCandidate {
        name: "fitted scalar".into(),
        coordinates: None,
        transition: None,
        projection: vec![],
    });
    let mut report = discover_interlacing_refinements(&input).unwrap();
    assert_eq!(report.status, "proved");
    // This work cap cannot fit the transition, but checking the stored
    // intertwiner needs no fit cells and still establishes the same theorem.
    report.request.limits.maximum_fit_cells = 1;
    assert_ne!(
        discover_interlacing_refinements(&report.request)
            .unwrap()
            .status,
        "proved"
    );
    verify_interlacing_refinements(&report).unwrap();
    let mut bad = report.clone();
    bad.candidates[0].degree_conditions = false;
    assert!(verify_interlacing_refinements(&bad).is_err());
    let mut bad = report;
    bad.candidates[0].transition.as_mut().unwrap().matrix[0][0][0].coeffs[0][0] = "2".into();
    assert!(verify_interlacing_refinements(&bad).is_err());
}
