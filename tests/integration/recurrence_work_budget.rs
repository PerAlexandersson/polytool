//! The deterministic work budget of the adaptive recurrence search.

use polytool::recurrence::{
    find_recurrence_adaptive_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome,
};

/// Rows of the Fibonacci-like recurrence `P(n) = P(n-1) + t P(n-2)`.
fn rows() -> Vec<Vec<i64>> {
    let mut rows: Vec<Vec<i64>> = vec![vec![1], vec![1]];
    while rows.len() < 10 {
        let a = rows[rows.len() - 1].clone();
        let b = rows[rows.len() - 2].clone();
        let mut next = vec![0; (b.len() + 1).max(a.len())];
        for (i, v) in a.iter().enumerate() {
            next[i] += v;
        }
        for (i, v) in b.iter().enumerate() {
            next[i + 1] += v;
        }
        rows.push(next);
    }
    rows
}

fn options() -> AdaptiveSearchOptions {
    AdaptiveSearchOptions {
        max_rec_len: 2,
        max_var_deg: 1,
        max_idx_deg: 1,
        max_diff_deg: 0,
        ..AdaptiveSearchOptions::default()
    }
}

fn run(budget: AdaptiveSearchBudget) -> AdaptiveSearchOutcome {
    find_recurrence_adaptive_with_budget(&rows(), &options(), budget)
}

fn work_units(outcome: &AdaptiveSearchOutcome) -> u64 {
    match outcome {
        AdaptiveSearchOutcome::Found(result) => result.diagnostics.work_units,
        AdaptiveSearchOutcome::NoRecurrence(summary)
        | AdaptiveSearchOutcome::BudgetExhausted(summary) => summary.diagnostics.work_units,
    }
}

fn describe(outcome: &AdaptiveSearchOutcome) -> String {
    match outcome {
        AdaptiveSearchOutcome::Found(result) => {
            format!(
                "found {} after {} work",
                result.recurrence,
                work_units(outcome)
            )
        }
        AdaptiveSearchOutcome::NoRecurrence(_) => format!("none after {}", work_units(outcome)),
        AdaptiveSearchOutcome::BudgetExhausted(summary) => format!(
            "exhausted after {} work, {} candidates",
            work_units(outcome),
            summary.diagnostics.considered_candidates
        ),
    }
}

#[test]
fn unbudgeted_search_reports_its_work_and_a_large_limit_changes_nothing() {
    let unbounded = run(AdaptiveSearchBudget::unbounded());
    let AdaptiveSearchOutcome::Found(reference) = &unbounded else {
        panic!("expected a recurrence: {}", describe(&unbounded));
    };
    let total = reference.diagnostics.work_units;
    assert!(total > 0);
    for limit in [total, total + 1, total * 10, u64::MAX] {
        let limited = run(AdaptiveSearchBudget::unbounded().with_max_work(limit));
        let AdaptiveSearchOutcome::Found(result) = &limited else {
            panic!("limit {limit}: {}", describe(&limited));
        };
        assert_eq!(
            result.recurrence.to_string(),
            reference.recurrence.to_string()
        );
        assert_eq!(result.candidates_tried, reference.candidates_tried);
        assert_eq!(result.diagnostics.work_units, total);
    }
}

#[test]
fn a_limit_below_the_total_exhausts_and_never_overshoots() {
    let total = work_units(&run(AdaptiveSearchBudget::unbounded()));
    let at_total_minus_one = run(AdaptiveSearchBudget::unbounded().with_max_work(total - 1));
    assert!(
        matches!(
            at_total_minus_one,
            AdaptiveSearchOutcome::BudgetExhausted(_)
        ),
        "{}",
        describe(&at_total_minus_one)
    );
    for limit in [0, 1, 5, 40, 200, total / 2, total - 1] {
        let outcome = run(AdaptiveSearchBudget::unbounded().with_max_work(limit));
        assert!(
            matches!(outcome, AdaptiveSearchOutcome::BudgetExhausted(_)),
            "limit {limit}: {}",
            describe(&outcome)
        );
        assert!(work_units(&outcome) <= limit, "limit {limit}");
    }
}

#[test]
fn the_same_limit_gives_the_same_outcome_every_time() {
    let total = work_units(&run(AdaptiveSearchBudget::unbounded()));
    for limit in [3, 77, total / 3, total / 2, total - 1, total] {
        let budget = AdaptiveSearchBudget::unbounded().with_max_work(limit);
        let first = run(budget);
        for _ in 0..3 {
            assert_eq!(describe(&run(budget)), describe(&first), "limit {limit}");
        }
    }
}

#[test]
fn larger_limits_never_stop_earlier() {
    let total = work_units(&run(AdaptiveSearchBudget::unbounded()));
    let mut previous_candidates = 0;
    for limit in (0..=total).step_by((total / 25).max(1) as usize) {
        let outcome = run(AdaptiveSearchBudget::unbounded().with_max_work(limit));
        if let AdaptiveSearchOutcome::BudgetExhausted(summary) = &outcome {
            assert!(summary.diagnostics.considered_candidates >= previous_candidates);
            previous_candidates = summary.diagnostics.considered_candidates;
        }
    }
}

#[test]
fn candidate_and_work_limits_combine() {
    let outcome = run(AdaptiveSearchBudget::limited(1).with_max_work(u64::MAX));
    assert!(matches!(outcome, AdaptiveSearchOutcome::BudgetExhausted(_)));
    let zero = run(AdaptiveSearchBudget::unbounded().with_max_work(0));
    let AdaptiveSearchOutcome::BudgetExhausted(summary) = zero else {
        panic!("zero work must exhaust");
    };
    assert_eq!(summary.diagnostics.considered_candidates, 0);
    assert_eq!(summary.diagnostics.work_units, 0);
}
