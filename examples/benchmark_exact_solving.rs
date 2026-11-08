//! Benchmark bounded recurrence searches that exercise exact linear solving.
//! Run with `--release`; `--repeat N` and `--only SUBSTRING` select runs.
//! JSON Lines include exact outcome keys for comparisons between builds.

use polytool::recurrence::{
    find_recurrence_adaptive_rational_with_budget, parse_rational_coeff, AdaptiveSearchBudget,
    AdaptiveSearchOptions, AdaptiveSearchOutcome,
};
use std::hint::black_box;
use std::time::Instant;

fn outcome_key(outcome: &AdaptiveSearchOutcome) -> String {
    let (tag, tried, d, found) = match outcome {
        AdaptiveSearchOutcome::Found(r) => (
            "found",
            r.candidates_tried,
            &r.diagnostics,
            format!(
                "{:?}|{:?}|{}|{}|{}|{}",
                r.recurrence,
                r.opts,
                r.num_unknowns,
                r.weighted_unknowns,
                r.fit_polynomials,
                r.verification_polynomials
            ),
        ),
        AdaptiveSearchOutcome::NoRecurrence(r) => {
            ("none", r.candidates_tried, &r.diagnostics, String::new())
        }
        AdaptiveSearchOutcome::BudgetExhausted(r) => {
            ("budget", r.candidates_tried, &r.diagnostics, String::new())
        }
    };
    format!(
        "{tag}|{tried}|{found}|{:?}",
        [
            d.generated_candidates,
            d.considered_candidates,
            d.insufficient_fit_rows,
            d.equation_bound_rejections,
            d.degree_bound_rejections,
            d.structural_zero_rejections,
            d.modular_prefilter_rejections,
            d.exact_solve_attempts,
            d.modular_lift_attempts,
            d.modular_lift_successes,
            d.modular_lift_failures,
            d.failed_exact_solves,
            d.heldout_verification_failures,
            usize::from(d.denominator_escalation_entered)
        ]
    )
}

fn main() {
    let mut repeat = 5usize;
    let mut only = String::new();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--repeat" => repeat = args.next().unwrap().parse().unwrap(),
            "--only" => only = args.next().unwrap(),
            _ => panic!("use --repeat N --only SUBSTRING"),
        }
    }
    assert!(repeat > 0);
    for (name, slug, large, unsuccessful) in [
        ("easy", "03_binomial_powers", false, false),
        ("derivative", "12_eulerian_derivative", false, false),
        ("large", "21_fifth_order_mixed_coefficients", true, false),
        ("wide", "21_fifth_order_mixed_coefficients", true, false),
        (
            "denominator",
            "29_oeis_a101920_three_lag_denominator_derivative",
            false,
            false,
        ),
        (
            "unsuccessful",
            "21_fifth_order_mixed_coefficients",
            false,
            true,
        ),
    ] {
        if !name.contains(&only) {
            continue;
        }
        let text = std::fs::read_to_string(format!(
            "{}/fixtures/recurrence-benchmarks/rows/{slug}.txt",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let rows = text
            .lines()
            .filter(|s| !s.trim().is_empty())
            .map(|line| {
                line.split(',')
                    .map(|s| parse_rational_coeff(s.trim()).unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut search = AdaptiveSearchOptions {
            max_rec_len: if large { 5 } else { 2 },
            max_idx_deg: if unsuccessful { 0 } else { 2 },
            max_var_deg: 2,
            max_diff_deg: if unsuccessful { 0 } else { 2 },
            try_laurent: false,
            ..Default::default()
        };
        if name == "denominator" {
            search.min_rec_len = 3;
            search.max_rec_len = 3;
            search.max_idx_deg = 1;
            search.max_diff_deg = 1;
            search.try_denominator = true;
            search.max_denom_idx_deg = 1;
            search.max_denom_var_deg = 0;
        }
        let candidate_limit = if name == "wide" {
            search.max_rec_len = 6;
            search.max_var_deg = 3;
            1500
        } else {
            500
        };
        let mut reference = None;
        // Warm once; paired comparisons should alternate executable order.
        for run in 0..=repeat {
            let start = Instant::now();
            let outcome = find_recurrence_adaptive_rational_with_budget(
                black_box(&rows),
                &search,
                AdaptiveSearchBudget::limited(candidate_limit),
            );
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            let key = outcome_key(&outcome);
            if let Some(expected) = &reference {
                assert_eq!(&key, expected);
            } else {
                reference = Some(key.clone());
            }
            let (status, d) = match &outcome {
                AdaptiveSearchOutcome::Found(r) => ("found", &r.diagnostics),
                AdaptiveSearchOutcome::NoRecurrence(r) => ("none", &r.diagnostics),
                AdaptiveSearchOutcome::BudgetExhausted(r) => ("budget", &r.diagnostics),
            };
            if run > 0 {
                println!(
                    "{}",
                    serde_json::json!({
                        "case": name, "run": run, "elapsed_ms": elapsed_ms,
                        "status": status, "considered": d.considered_candidates,
                        "search": search, "candidate_limit": candidate_limit, "outcome": key,
                        "modular_prefilter_ms": d.modular_prefilter_ms,
                        "exact_solve_ms": d.exact_solve_ms,
                    })
                );
            }
        }
    }
}
