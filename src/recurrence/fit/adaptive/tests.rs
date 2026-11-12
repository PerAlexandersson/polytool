//! Tests for the enclosing module.

use super::*;

#[test]
fn heldout_slots_include_denominators_derivatives_and_forcing() {
    let record: crate::recurrence::RecurrenceJson = serde_json::from_value(serde_json::json!({
        "schema": "polytool.recurrence.v1", "first_index": 1,
        "initial_polynomials": [["1"]],
        "recurrence": {
            "denominator": {"coeffs": [["0", "0", "1"]]},
            "inhomogeneous": {"coeffs": [["0", "0", "0", "0", "0", "0", "1"]]},
            "terms": [
                {"offset": 1, "deriv_order": 1, "sign": "none",
                 "coeff": {"coeffs": [["0", "0", "0", "0", "1"]]}},
                {"offset": 1, "deriv_order": 0, "sign": "none",
                 "coeff": {"coeffs": [["0", "0", "0", "0", "0", "0", "0", "-3"],
                                       ["0", "0", "0", "0", "0", "0", "0", "1"]]}}
            ]
        }
    }))
    .unwrap();
    let (mut rule, _, _) = record.to_recurrence_parts().unwrap();
    let rows = i64_polys_to_rational(&[vec![1], vec![0, 0], vec![1, 0], vec![1]]);
    // At n=3 both right-hand products vanish; forcing has degree 6.
    // At n=4 the (n-3)t^7 term survives, giving eight slots.
    assert_eq!(heldout_equation_count(&rows, &rule, 2), 7 + 8);
    rule.inhomogeneous = None;
    rule.terms.truncate(1);
    let rows = i64_polys_to_rational(&[vec![0, 0], vec![1], vec![0, 0, 1], vec![0]]);
    // Degrees: t^2*1; t^2*t^2; t^4*D(t^2).
    assert_eq!(heldout_equation_count(&rows, &rule, 1), 3 + 5 + 6);
    let zeros = i64_polys_to_rational(&[vec![0], vec![], vec![0, 0]]);
    assert_eq!(heldout_equation_count(&zeros, &rule, 1), 2);
    assert_eq!(heldout_equation_count(&zeros, &rule, zeros.len()), 0);
}

#[test]
fn adaptive_geometric() {
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
    assert!(result.verification_polynomials >= 1);
    assert!(result.fit_polynomials < polys.len());
}

#[test]
fn homogeneous_degree_bound_rejects_impossible_rhs_degree() {
    let polys = i64_polys_to_rational(&[vec![1], vec![1, 1], vec![1, 2, 1]]);
    let derivs = rational_derivatives_up_to(&polys, 0);
    let too_low = RecurrenceOptions {
        rec_len: 1,
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        ..Default::default()
    };
    let enough_t_degree = RecurrenceOptions {
        var_deg: 1,
        ..too_low.clone()
    };

    assert!(homogeneous_degree_bound_rejects(
        &polys, &derivs, &too_low, 3
    ));
    assert!(!homogeneous_degree_bound_rejects(
        &polys,
        &derivs,
        &enough_t_degree,
        3
    ));
}

#[test]
fn structural_zero_row_rejects_impossible_support_gap() {
    let polys = i64_polys_to_rational(&[vec![1, 0, 1], vec![0, 1]]);
    let derivs = rational_derivatives_up_to(&polys, 0);
    let prefix_max_degrees = prefix_max_degrees_rational(&polys);
    let no_t_shift = RecurrenceOptions {
        rec_len: 1,
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        ..Default::default()
    };
    let with_t_shift = RecurrenceOptions {
        var_deg: 1,
        ..no_t_shift.clone()
    };

    assert!(!homogeneous_degree_bound_rejects(
        &polys,
        &derivs,
        &no_t_shift,
        2
    ));
    assert!(structural_zero_row_rejects(
        &polys,
        &derivs,
        &no_t_shift,
        2,
        &prefix_max_degrees
    ));
    assert!(!structural_zero_row_rejects(
        &polys,
        &derivs,
        &with_t_shift,
        2,
        &prefix_max_degrees
    ));
}

#[test]
fn adaptive_no_verify_uses_all_rows() {
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
    let search = AdaptiveSearchOptions {
        no_verify: true,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
    assert_eq!(result.fit_polynomials, polys.len());
    assert_eq!(result.verification_polynomials, 0);
}

#[test]
fn adaptive_rejects_false_prefix_fit_on_heldout_row() {
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![17]];
    let search = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_var_deg: 0,
        max_idx_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    assert!(find_recurrence_adaptive(&polys, &search).is_none());
}

#[test]
fn weighted_complexity_delays_alternating_terms() {
    let ordinary = RecurrenceOptions {
        rec_len: 2,
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        ..Default::default()
    };
    let alternating = RecurrenceOptions {
        rec_len: 1,
        var_deg: 0,
        idx_deg: 0,
        diff_deg: 0,
        homogeneous: true,
        alternating_sign: true,
        ..Default::default()
    };
    assert_eq!(count_unknowns(&ordinary), count_unknowns(&alternating));
    assert!(count_weighted_unknowns(&ordinary) < count_weighted_unknowns(&alternating));
}

#[test]
fn adaptive_respects_min_bounds() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![2],
        vec![3],
        vec![5],
        vec![8],
        vec![13],
    ];
    let search = AdaptiveSearchOptions {
        min_rec_len: 2,
        max_rec_len: 2,
        min_var_deg: 0,
        max_var_deg: 0,
        min_idx_deg: 0,
        max_idx_deg: 0,
        min_diff_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
    assert_eq!(result.opts.rec_len, 2);
}

#[test]
fn adaptive_can_skip_prefix() {
    let polys: Vec<Vec<i64>> = vec![
        vec![9],
        vec![1],
        vec![2],
        vec![4],
        vec![8],
        vec![16],
        vec![32],
    ];
    let search = AdaptiveSearchOptions {
        skip_prefix: 1,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
}

// --- Adaptive search tests ---

#[test]
fn adaptive_fibonacci() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![2],
        vec![3],
        vec![5],
        vec![8],
        vec![13],
    ];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
    assert_eq!(result.opts.rec_len, 2);
    assert_eq!(result.opts.diff_deg, 0);
}

#[test]
fn adaptive_modular_prefilter_finds_fibonacci() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![2],
        vec![3],
        vec![5],
        vec![8],
        vec![13],
    ];
    let search = AdaptiveSearchOptions {
        modular_prefilter: true,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = P(n-1) + P(n-2)");
    assert!(result.opts.modular_prefilter);
}

#[test]
fn adaptive_modular_lift_skips_exact_solve_for_integer_solution() {
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![4], vec![8], vec![16], vec![32]];
    let search = AdaptiveSearchOptions {
        min_rec_len: 1,
        max_rec_len: 1,
        min_var_deg: 0,
        max_var_deg: 0,
        min_idx_deg: 0,
        max_idx_deg: 0,
        min_diff_deg: 0,
        max_diff_deg: 0,
        modular_prefilter: true,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = 2 P(n-1)");
    assert_eq!(result.diagnostics.modular_lift_attempts, 1);
    assert_eq!(result.diagnostics.modular_lift_successes, 1);
    assert_eq!(result.diagnostics.exact_solve_attempts, 0);
}

#[test]
fn adaptive_rank_deficient_modular_lift_can_skip_exact_solve() {
    let mut polys = Vec::new();
    let mut value = 1_i64;
    for _ in 0..32 {
        polys.push(vec![value]);
        value *= 2;
    }
    let search = AdaptiveSearchOptions {
        min_rec_len: 5,
        max_rec_len: 5,
        min_var_deg: 3,
        max_var_deg: 3,
        min_idx_deg: 3,
        max_idx_deg: 3,
        min_diff_deg: 0,
        max_diff_deg: 0,
        modular_prefilter: true,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert!(recurrence_holds_rational(
        &i64_polys_to_rational(&polys),
        &result.recurrence,
        result.opts.rec_len
    ));
    assert_eq!(result.diagnostics.modular_lift_attempts, 1);
    assert_eq!(result.diagnostics.modular_lift_successes, 1);
    assert_eq!(result.diagnostics.exact_solve_attempts, 0);
}

#[test]
fn adaptive_factorial() {
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![2], vec![6], vec![24], vec![120]];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = n P(n-1)");
    assert_eq!(result.opts.rec_len, 1);
    assert_eq!(result.opts.idx_deg, 1);
}

#[test]
fn adaptive_binomial() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1, 1],
        vec![1, 2, 1],
        vec![1, 3, 3, 1],
        vec![1, 4, 6, 4, 1],
    ];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = (1 + t) P(n-1)");
    assert_eq!(result.opts.rec_len, 1);
}

#[test]
fn adaptive_chebyshev() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![0, 1],
        vec![-1, 0, 2],
        vec![0, -3, 0, 4],
        vec![1, 0, -8, 0, 8],
    ];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(
        format!("{}", result.recurrence),
        "P(n) = 2t P(n-1) - P(n-2)"
    );
}

#[test]
fn adaptive_eulerian() {
    // With held-out verification, prefix-only non-derivative fits are
    // rejected, and the standard derivative recurrence is found.
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![1, 1],
        vec![1, 4, 1],
        vec![1, 11, 11, 1],
        vec![1, 26, 66, 26, 1],
    ];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default()).unwrap();
    assert_eq!(
        format!("{}", result.recurrence),
        "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
    );
    assert_eq!(result.opts.diff_deg, 1);
    assert_eq!(result.opts.rec_len, 1);

    // With max_rec_len=1, the same derivative recurrence is still found.
    let search = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_diff_deg: 2,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(
        format!("{}", result.recurrence),
        "P(n) = (1 - 2t + nt) P(n-1) + (t - t^2) P'(n-1)"
    );
    assert_eq!(result.opts.diff_deg, 1);
}

#[test]
fn adaptive_with_alternating_sign_builder_toggles_the_search_space() {
    let polys: Vec<Vec<i64>> = [1, 1, -1, -1, 1, 1, -1, -1]
        .iter()
        .map(|&c| vec![c])
        .collect();
    let base = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_var_deg: 0,
        max_idx_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    assert!(!base.try_alternating_sign);
    assert!(find_recurrence_adaptive(&polys, &base).is_none());
    let enabled = base.clone().with_alternating_sign(true);
    assert!(enabled.try_alternating_sign);
    let result = find_recurrence_adaptive(&polys, &enabled).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = (-1)^n P(n-1)");
    assert!(!enabled.with_alternating_sign(false).try_alternating_sign);
}

#[test]
fn adaptive_alternating_sign() {
    let polys: Vec<Vec<i64>> = vec![
        vec![1],
        vec![1],
        vec![-1],
        vec![-1],
        vec![1],
        vec![1],
        vec![-1],
        vec![-1],
    ];
    let search = AdaptiveSearchOptions {
        try_alternating_sign: true,
        max_rec_len: 1,
        max_var_deg: 0,
        max_idx_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    let result = find_recurrence_adaptive(&polys, &search).unwrap();
    assert_eq!(format!("{}", result.recurrence), "P(n) = (-1)^n P(n-1)");
    assert!(result.opts.alternating_sign);
}

#[test]
fn adaptive_short_sequence() {
    // m=3 with constant polys: too few equations for any 2-term recurrence.
    let polys: Vec<Vec<i64>> = vec![vec![1], vec![1], vec![2]];
    let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
    assert!(result.is_none());
}

fn fixed_adaptive_search(max_rec_len: usize) -> AdaptiveSearchOptions {
    AdaptiveSearchOptions {
        min_rec_len: 1,
        max_rec_len,
        min_var_deg: 0,
        max_var_deg: 0,
        min_idx_deg: 0,
        max_idx_deg: 0,
        min_diff_deg: 0,
        max_diff_deg: 0,
        try_denominator: false,
        try_inhomogeneous: false,
        try_alternating_sign: false,
        // A fixed space of `max_rec_len` candidates: no left-hand variants.
        try_laurent: false,
        ..Default::default()
    }
}

#[test]
fn adaptive_budget_zero_stops_before_the_first_candidate() {
    let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
    let outcome = find_recurrence_adaptive_with_budget(
        &polys,
        &fixed_adaptive_search(2),
        AdaptiveSearchBudget::limited(0),
    );

    let AdaptiveSearchOutcome::BudgetExhausted(summary) = outcome else {
        panic!("zero budget must stop a nonempty search");
    };
    assert_eq!(summary.diagnostics.considered_candidates, 0);
    assert_eq!(summary.candidates_tried, 0);
}

#[test]
fn adaptive_budget_allows_success_on_the_exact_boundary() {
    let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
    let outcome = find_recurrence_adaptive_with_budget(
        &polys,
        &fixed_adaptive_search(2),
        AdaptiveSearchBudget::limited(1),
    );

    let AdaptiveSearchOutcome::Found(result) = outcome else {
        panic!("the first candidate should be allowed to succeed");
    };
    assert_eq!(result.recurrence.to_string(), "P(n) = 2 P(n-1)");
    assert_eq!(result.diagnostics.considered_candidates, 1);
    assert_eq!(result.candidates_tried, 1);
}

#[test]
fn adaptive_budget_distinguishes_exhaustion_from_complete_failure() {
    let polys = vec![vec![1], vec![1], vec![2], vec![3], vec![5]];

    let complete = find_recurrence_adaptive_with_budget(
        &polys,
        &fixed_adaptive_search(1),
        AdaptiveSearchBudget::limited(1),
    );
    let AdaptiveSearchOutcome::NoRecurrence(summary) = complete else {
        panic!("one allowed candidate should exhaust a one-candidate search space");
    };
    assert_eq!(summary.diagnostics.considered_candidates, 1);

    let truncated = find_recurrence_adaptive_with_budget(
        &polys,
        &fixed_adaptive_search(2),
        AdaptiveSearchBudget::limited(1),
    );
    let AdaptiveSearchOutcome::BudgetExhausted(summary) = truncated else {
        panic!("a second uninspected candidate must produce budget exhaustion");
    };
    assert_eq!(summary.diagnostics.considered_candidates, 1);
}

#[test]
fn unbounded_budget_preserves_the_existing_search_result() {
    let polys = vec![vec![1], vec![2], vec![4], vec![8], vec![16]];
    let search = fixed_adaptive_search(2);
    let legacy = find_recurrence_adaptive(&polys, &search).unwrap();
    let outcome =
        find_recurrence_adaptive_with_budget(&polys, &search, AdaptiveSearchBudget::unbounded());
    let AdaptiveSearchOutcome::Found(unbounded) = outcome else {
        panic!("unbounded outcome API must find the legacy recurrence");
    };

    assert_eq!(
        unbounded.recurrence.to_string(),
        legacy.recurrence.to_string()
    );
    assert_eq!(unbounded.candidates_tried, legacy.candidates_tried);
    assert_eq!(
        unbounded.diagnostics.considered_candidates,
        legacy.diagnostics.considered_candidates
    );
}
