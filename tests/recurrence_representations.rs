//! Recurrence representation search through the library API.

use num_traits::Zero;
use polytool::recurrence::*;

fn rows(id: &str, count: usize) -> Vec<Vec<BigRational>> {
    polytool::oeis::by_id(id)
        .unwrap()
        .generate_rows(count)
        .unwrap()
        .into_iter()
        .map(|mut p| {
            while p.len() > 1 && p.last().unwrap().is_zero() {
                p.pop();
            }
            p.into_iter().map(BigRational::from_integer).collect()
        })
        .collect()
}

fn bounds() -> AdaptiveSearchOptions {
    AdaptiveSearchOptions {
        max_rec_len: 3,
        max_diff_deg: 2,
        max_var_deg: 3,
        max_idx_deg: 1,
        try_denominator: true,
        max_denom_idx_deg: 1,
        max_denom_var_deg: 0,
        ..Default::default()
    }
}

#[test]
fn three_lag_derivative_tradeoffs_replay_beyond_the_fitted_prefix() {
    let input = rows("A171608", 60);
    let menu = find_recurrence_representations_rational(
        &input[..40],
        &bounds(),
        AdaptiveSearchBudget::unbounded(),
    )
    .unwrap();
    assert_eq!(menu.status, RepresentationSearchStatus::Complete);
    for pair in [(3, 0), (2, 1), (1, 2)] {
        assert!(
            menu.frontier
                .iter()
                .map(|&i| &menu.representations[i])
                .any(|r| (r.complexity.rec_len, r.complexity.diff_deg) == pair),
            "missing {pair:?}: {:?}",
            menu.representations
        );
    }
    assert_eq!(menu.coverage.len(), 9);
    assert!(
        menu.representations.iter().any(|r| r
            .recurrence_json
            .recurrence
            .denominator
            .as_ref()
            .is_some_and(|p| p.coeffs[0][0] == "0")),
        "must support the denominator chart with zero constant term"
    );
    for item in &menu.representations {
        assert_eq!(item.evidence, "exact_fit");
        assert!(item.verification_polynomials > 0);
        let (rec, first, initial) = item.recurrence_json.to_recurrence_parts().unwrap();
        assert_eq!(initial.len() + 1, item.first_checked_index);
        assert!(recurrence_holds_from_rational(
            &input,
            &rec,
            rec.max_offset(),
            item.first_checked_index
        ));
        assert_eq!(
            rec.generate_rows_rational(&initial, first, 60).unwrap(),
            input
        );
    }
    let json = serde_json::to_string(&menu).unwrap();
    let decoded: RecurrenceRepresentations = serde_json::from_str(&json).unwrap();
    assert_eq!(menu.frontier, decoded.frontier);
}

#[test]
fn altered_held_out_row_rejects_the_prefix_fit() {
    let input: Vec<Vec<i64>> = [1, 2, 4, 8, 16, 33].into_iter().map(|n| vec![n]).collect();
    let opts = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_var_deg: 0,
        max_idx_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    let menu =
        find_recurrence_representations(&input, &opts, AdaptiveSearchBudget::unbounded()).unwrap();
    assert!(menu.representations.is_empty());
}

#[test]
fn global_budget_reports_partial_coverage_and_zero_budget_does_no_work() {
    let input = rows("A171608", 20);
    for limit in [0, 1, 8] {
        let menu = find_recurrence_representations_rational(
            &input,
            &bounds(),
            AdaptiveSearchBudget::limited(limit),
        )
        .unwrap();
        assert_eq!(menu.status, RepresentationSearchStatus::BudgetExhausted);
        assert_eq!(menu.configurations_considered, limit);
        assert!(menu.coverage.iter().any(|c| c.status == "incomplete"));
    }
}

#[test]
fn tail_export_retains_extra_initial_rows_when_actual_lag_shrinks() {
    // All but the first row form a geometric progression. The rank-two shape
    // must not export just its one actually used initial row.
    let input: Vec<Vec<i64>> = [99, 2, 4, 8, 16, 32, 64, 128]
        .into_iter()
        .map(|n| vec![n])
        .collect();
    let opts = AdaptiveSearchOptions {
        min_rec_len: 2,
        max_rec_len: 2,
        max_var_deg: 0,
        max_idx_deg: 0,
        max_diff_deg: 0,
        ..Default::default()
    };
    let menu =
        find_recurrence_representations(&input, &opts, AdaptiveSearchBudget::unbounded()).unwrap();
    let item = &menu.representations[0];
    assert_eq!(item.complexity.rec_len, 1);
    assert_eq!(item.first_checked_index, 3);
    let (rec, first, initial) = item.recurrence_json.to_recurrence_parts().unwrap();
    assert_eq!(initial.len(), 2);
    assert_eq!(
        rec.generate_rows_rational(&initial, first, 8).unwrap()[0],
        initial[0]
    );
    let all: Vec<_> = input
        .iter()
        .map(|p| {
            p.iter()
                .map(|n| BigRational::from_integer((*n).into()))
                .collect::<Vec<_>>()
        })
        .collect();
    assert_eq!(rec.generate_rows_rational(&initial, first, 8).unwrap(), all);
}

#[test]
fn skip_prefix_and_invalid_bounds_are_explicit() {
    let input: Vec<Vec<i64>> = [99, 1, 2, 4, 8, 16, 32]
        .into_iter()
        .map(|n| vec![n])
        .collect();
    let mut opts = AdaptiveSearchOptions {
        skip_prefix: 1,
        max_rec_len: 1,
        max_diff_deg: 0,
        max_var_deg: 0,
        max_idx_deg: 0,
        ..Default::default()
    };
    let menu =
        find_recurrence_representations(&input, &opts, AdaptiveSearchBudget::unbounded()).unwrap();
    assert_eq!(menu.representations[0].first_checked_input_row, 3);
    assert_eq!(menu.representations[0].first_checked_index, 2);
    opts.min_diff_deg = 2;
    assert!(
        find_recurrence_representations(&input, &opts, AdaptiveSearchBudget::unbounded()).is_err()
    );
}

#[test]
fn denominator_search_is_opt_in_and_empty_search_is_not_a_nonexistence_proof() {
    let input = rows("A171608", 30);
    let opts = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_diff_deg: 0,
        max_var_deg: 0,
        max_idx_deg: 0,
        try_denominator: false,
        ..Default::default()
    };
    let menu =
        find_recurrence_representations_rational(&input, &opts, AdaptiveSearchBudget::unbounded())
            .unwrap();
    assert_eq!(menu.status, RepresentationSearchStatus::Complete);
    assert!(menu.representations.is_empty());
    assert_eq!(menu.coverage[0].status, "not_found");
}

#[test]
fn projective_denominator_fallback_preserves_alternating_index_convention() {
    // Exact fixture: n P(n) = (-1)^n P(n-1), starting at n=2.
    let input: Vec<_> = [
        "1", "1/2", "-1/6", "-1/24", "1/120", "1/720", "-1/5040", "-1/40320",
    ]
    .into_iter()
    .map(|s| vec![parse_rational_coeff(s).unwrap()])
    .collect();
    let opts = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_diff_deg: 0,
        max_var_deg: 0,
        max_idx_deg: 0,
        try_alternating_sign: true,
        try_denominator: true,
        max_denom_idx_deg: 1,
        max_denom_var_deg: 0,
        ..Default::default()
    };
    let report =
        find_recurrence_representations_rational(&input, &opts, AdaptiveSearchBudget::unbounded())
            .unwrap();
    let candidate = &report.representations[0];
    assert!(candidate.complexity.alternating_sign);
    let (rec, first, initial) = candidate.recurrence_json.to_recurrence_parts().unwrap();
    assert_eq!(
        rec.generate_rows_rational(&initial, first, input.len())
            .unwrap(),
        input
    );
}

#[test]
fn variable_denominator_fit_does_not_promise_polynomial_continuation() {
    // Fixture x^6, ..., x: x P(n) = P(n-1). Polynomial continuation ends at 1.
    let input: Vec<Vec<i64>> = (1..=6)
        .rev()
        .map(|d| {
            let mut row = vec![0; d + 1];
            row[d] = 1;
            row
        })
        .collect();
    let opts = AdaptiveSearchOptions {
        max_rec_len: 1,
        max_diff_deg: 0,
        max_var_deg: 0,
        max_idx_deg: 0,
        try_denominator: true,
        max_denom_idx_deg: 0,
        max_denom_var_deg: 1,
        ..Default::default()
    };
    let report =
        find_recurrence_representations(&input, &opts, AdaptiveSearchBudget::unbounded()).unwrap();
    let candidate = &report.representations[0];
    assert_eq!(candidate.complexity.denom_var_deg, 1);
    assert_eq!(candidate.evidence, "exact_fit");
    let (rec, first, initial) = candidate.recurrence_json.to_recurrence_parts().unwrap();
    assert!(rec.generate_rows_rational(&initial, first, 7).is_ok());
    assert!(rec.generate_rows_rational(&initial, first, 8).is_err());
}
