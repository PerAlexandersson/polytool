//! A fixed set of real workloads of the exact linear algebra (`linalg`) and
//! of the recurrence fit built on it, with recorded output fingerprints.
//!
//! The regression test runs every workload and compares a fingerprint of its
//! complete output with `tests/fixtures/linalg-workload-fingerprints.tsv`,
//! so an optimisation of the elimination cannot change a result (the rule
//! found, the kernel basis and its certificate, the budget at which a kernel
//! computation stops) unnoticed.  The slowest workloads are compared only by
//! the ignored full test.  Re-bless only for an intentional change of
//! results: `POLYTOOL_BLESS_LINALG_WORKLOADS=1 cargo test --release --test
//! linalg_workloads -- --ignored every_`.
//!
//! The ignored timing test prints the median wall time of each workload:
//! `cargo test --release --test linalg_workloads -- --ignored --nocapture time_`
//! (`LINALG_WORKLOAD_REPEAT` sets the number of runs, `LINALG_WORKLOAD_ONLY` a
//! name filter).

use num_bigint::BigInt;
use polytool::linalg::{
    rational_kernel_with_certificate, solve_full_rank_square_linear_system, solve_linear_system,
    solve_linear_system_with_rank, RationalKernelLimits,
};
use polytool::recurrence::{
    find_hadamard_factor_lifts, find_recurrence_adaptive_rational_with_budget,
    AdaptiveSearchBudget, AdaptiveSearchOptions, AdaptiveSearchOutcome, HadamardFactorLiftOptions,
    RecurrenceJson,
};
use polytool::sequences;
use polytool::BigRational;
use std::collections::BTreeMap;
use std::time::Instant;

type Q = BigRational;

const FINGERPRINTS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/linalg-workload-fingerprints.tsv"
);

mod common;
use common::q;

/// 64-bit FNV-1a of the text (stable across platforms and toolchains).
fn fingerprint(text: &str) -> String {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in text.bytes() {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0100_0000_01b3);
    }
    format!("{hash:016x}:{}", text.len())
}

// ---------------------------------------------------------------------------
// Adaptive recurrence searches
// ---------------------------------------------------------------------------

fn rational_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| row.iter().cloned().map(Q::from_integer).collect())
        .collect()
}

/// Every deterministic part of the outcome (no timings).
fn adaptive_workload(rows: &[Vec<Q>], options: &AdaptiveSearchOptions) -> String {
    let outcome = find_recurrence_adaptive_rational_with_budget(
        rows,
        options,
        AdaptiveSearchBudget::unbounded(),
    );
    let (label, summary_text) = match outcome {
        AdaptiveSearchOutcome::Found(result) => {
            let d = &result.diagnostics;
            (
                "found",
                format!(
                    "{:?}\n{:?}\n{} {} {} {} {} {}\n{} {} {} {} {} {} {} {} {} {} {} {} {} {}",
                    result.recurrence,
                    result.opts,
                    result.num_unknowns,
                    result.weighted_unknowns,
                    result.num_equations,
                    result.fit_polynomials,
                    result.verification_polynomials,
                    result.candidates_tried,
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
                    d.denominator_escalation_entered,
                ),
            )
        }
        AdaptiveSearchOutcome::NoRecurrence(summary) => ("none", summary_text(&summary)),
        AdaptiveSearchOutcome::BudgetExhausted(_) => ("budget", String::new()),
    };
    format!("{label}\n{summary_text}")
}

fn summary_text(summary: &polytool::recurrence::AdaptiveSearchSummary) -> String {
    let d = &summary.diagnostics;
    format!(
        "{} {} {} {} {} {} {} {} {} {} {} {} {} {}",
        summary.candidates_tried,
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
    )
}

fn search_options(rec: usize, var: usize, idx: usize, diff: usize) -> AdaptiveSearchOptions {
    AdaptiveSearchOptions {
        max_rec_len: rec,
        max_var_deg: var,
        max_idx_deg: idx,
        max_diff_deg: diff,
        ..AdaptiveSearchOptions::default()
    }
}

// ---------------------------------------------------------------------------
// Hadamard-factor lifts (kernels of the row-annihilator systems)
// ---------------------------------------------------------------------------

/// The descent polynomials of 321-avoiding permutations, `n >= 1`.
const AV321_DESCENTS: &str = r#"{"schema":"polytool.recurrence.v1","first_index":1,
  "initial_polynomials":[["1"],["1","1"],["1","4"]],
  "recurrence":{"terms":[
    {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-1","-1"],["5","-1"]]}},
    {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["10","-10"],["-8","8"]]}},
    {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-8","16","-8"],["4","-8","4"]]}}],
  "denominator":{"coeffs":[["1"],["1"]]}}}"#;

/// The peak polynomials of 123-avoiding permutations, `n >= 1`.
const AV123_PEAKS: &str = r#"{"schema":"polytool.recurrence.v1","first_index":1,
  "initial_polynomials":[["1"],["2"],["3","2"]],
  "recurrence":{"terms":[
    {"offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["0","-2"],["5","-1"]]}},
    {"offset":2,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["6","-6"],["-8","8"]]}},
    {"offset":3,"deriv_order":0,"sign":"none","coeff":{"coeffs":[["-4","8","-4"],["4","-8","4"]]}}],
  "denominator":{"coeffs":[["1"],["1"]]}}}"#;

fn hadamard_workload(json: &str) -> String {
    let record: RecurrenceJson = serde_json::from_str(json).unwrap();
    let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
    let report =
        find_hadamard_factor_lifts(&rule, &seeds, first, &HadamardFactorLiftOptions::default());
    format!("{report:?}")
}

// ---------------------------------------------------------------------------
// Direct systems and kernels
// ---------------------------------------------------------------------------

struct Lcg(u64);

impl Lcg {
    fn next(&mut self, modulus: i64) -> i64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        ((self.0 >> 33) % modulus as u64) as i64
    }
}

/// A `rows x columns` integer matrix of rank at most `rank`, built as a
/// product of a `rows x rank` and a `rank x columns` matrix with small
/// entries, with a fraction of the entries zeroed (sparse like the residual
/// matrices of the fits).  `denominators` divides some entries by 2, 3.
fn low_rank_matrix(
    seed: u64,
    rows: usize,
    columns: usize,
    rank: usize,
    sparsity: i64,
    denominators: bool,
) -> Vec<Vec<Q>> {
    let mut rng = Lcg(seed);
    let left: Vec<Vec<i64>> = (0..rows)
        .map(|_| {
            (0..rank)
                .map(|_| {
                    if rng.next(sparsity) == 0 {
                        rng.next(9) - 4
                    } else {
                        0
                    }
                })
                .collect()
        })
        .collect();
    let right: Vec<Vec<i64>> = (0..rank)
        .map(|_| (0..columns).map(|_| rng.next(11) - 5).collect())
        .collect();
    (0..rows)
        .map(|i| {
            (0..columns)
                .map(|j| {
                    let value: i64 = (0..rank).map(|k| left[i][k] * right[k][j]).sum();
                    if denominators && (i + 2 * j) % 5 == 0 {
                        Q::new(BigInt::from(value), BigInt::from(1 + (i + j) as i64 % 3))
                    } else {
                        q(value)
                    }
                })
                .collect()
        })
        .collect()
}

fn matrix_times_vector(matrix: &[Vec<Q>], vector: &[Q]) -> Vec<Q> {
    matrix
        .iter()
        .map(|row| row.iter().zip(vector).map(|(a, b)| a * b).sum())
        .collect()
}

fn solve_workload(
    seed: u64,
    rows: usize,
    columns: usize,
    rank: usize,
    sparsity: i64,
    denominators: bool,
) -> String {
    let matrix = low_rank_matrix(seed, rows, columns, rank, sparsity, denominators);
    let mut rng = Lcg(seed ^ 0xabcd);
    let witness: Vec<Q> = (0..columns).map(|_| q(rng.next(7) - 3)).collect();
    let rhs = matrix_times_vector(&matrix, &witness);
    let mut inconsistent = rhs.clone();
    inconsistent[rows - 1] += q(1);
    format!(
        "{:?}\n{:?}\n{:?}\n{:?}",
        solve_linear_system_with_rank(&matrix, &rhs),
        solve_linear_system(&matrix, &rhs),
        solve_linear_system_with_rank(&matrix, &inconsistent),
        solve_linear_system_with_rank(&matrix, &vec![q(0); rows]),
    )
}

fn square_solve_workload(seed: u64, size: usize) -> String {
    let matrix = low_rank_matrix(seed, size, size, size, 1, false);
    let mut rng = Lcg(seed ^ 0x1234);
    let witness: Vec<Q> = (0..size).map(|_| q(rng.next(7) - 3)).collect();
    let rhs = matrix_times_vector(&matrix, &witness);
    format!("{:?}", solve_full_rank_square_linear_system(&matrix, &rhs))
}

const UNLIMITED: RationalKernelLimits = RationalKernelLimits {
    max_entries: usize::MAX / 4,
    max_operations: usize::MAX,
    max_coefficient_bits: u64::MAX,
};

/// The complete kernel, and the exact budgets at which it stops: the
/// smallest `max_operations` and `max_coefficient_bits` that still succeed
/// (so the budget accounting is pinned, not just the kernel).
fn kernel_workload(
    seed: u64,
    rows: usize,
    columns: usize,
    rank: usize,
    sparsity: i64,
    denominators: bool,
) -> String {
    let matrix = low_rank_matrix(seed, rows, columns, rank, sparsity, denominators);
    let kernel = rational_kernel_with_certificate(&matrix, columns, UNLIMITED);
    let succeeds = |limits: RationalKernelLimits| {
        rational_kernel_with_certificate(&matrix, columns, limits).is_ok()
    };
    // Bisect the smallest sufficient number of operations, then of bits.
    let (mut low, mut high) = (0usize, 100_000_000usize);
    while low < high {
        let mid = low + (high - low) / 2;
        let limits = RationalKernelLimits {
            max_operations: mid,
            ..UNLIMITED
        };
        if succeeds(limits) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    let minimal_operations = low;
    let (mut low, mut high) = (0u64, 100_000u64);
    while low < high {
        let mid = low + (high - low) / 2;
        let limits = RationalKernelLimits {
            max_coefficient_bits: mid,
            ..UNLIMITED
        };
        if succeeds(limits) {
            high = mid;
        } else {
            low = mid + 1;
        }
    }
    let minimal_bits = low;
    let below_operations = rational_kernel_with_certificate(
        &matrix,
        columns,
        RationalKernelLimits {
            max_operations: minimal_operations.saturating_sub(1),
            ..UNLIMITED
        },
    )
    .map(|_| ());
    let below_bits = rational_kernel_with_certificate(
        &matrix,
        columns,
        RationalKernelLimits {
            max_coefficient_bits: minimal_bits.saturating_sub(1),
            ..UNLIMITED
        },
    )
    .map(|_| ());
    format!(
        "{kernel:?}\nops {minimal_operations} bits {minimal_bits}\n{below_operations:?}\n{below_bits:?}"
    )
}

/// The kernel only (no budget bisection), for the larger systems.
fn kernel_only_workload(
    seed: u64,
    rows: usize,
    columns: usize,
    rank: usize,
    sparsity: i64,
    denominators: bool,
) -> String {
    let matrix = low_rank_matrix(seed, rows, columns, rank, sparsity, denominators);
    format!(
        "{:?}",
        rational_kernel_with_certificate(&matrix, columns, UNLIMITED)
    )
}

/// The system of the recurrence-fit residuals of a family: the kernel of the
/// matrix of the linear recurrence relations of order `order` among the
/// leading coefficient columns of the rows (a Hankel-like integer matrix).
fn hankel_kernel_workload(rows: &[Vec<BigInt>], order: usize, width: usize) -> String {
    let columns = order * width;
    let matrix: Vec<Vec<Q>> = (order..rows.len())
        .map(|n| {
            (0..columns)
                .map(|c| {
                    let (lag, coefficient) = (c / width + 1, c % width);
                    rows[n - lag]
                        .get(coefficient)
                        .cloned()
                        .map(Q::from_integer)
                        .unwrap_or_else(|| q(0))
                })
                .collect()
        })
        .collect();
    format!(
        "{:?}",
        rational_kernel_with_certificate(&matrix, columns, UNLIMITED)
    )
}

// ---------------------------------------------------------------------------
// The workload table
// ---------------------------------------------------------------------------

type Workload = (String, Box<dyn Fn() -> String>);

fn workloads() -> Vec<Workload> {
    let mut out: Vec<Workload> = Vec::new();
    let add = |out: &mut Vec<Workload>, name: &str, run: Box<dyn Fn() -> String>| {
        out.push((name.to_string(), run));
    };

    add(
        &mut out,
        "adaptive/eulerian/default",
        Box::new(|| {
            let rows = rational_rows(&sequences::eulerian_polynomials_bigint(12)[1..]);
            adaptive_workload(&rows, &AdaptiveSearchOptions::default())
        }),
    );
    add(
        &mut out,
        "adaptive/eulerian/laurent",
        Box::new(|| {
            let rows = rational_rows(&sequences::eulerian_polynomials_bigint(14)[1..]);
            adaptive_workload(&rows, &search_options(2, 2, 2, 1).with_laurent(true))
        }),
    );
    add(
        &mut out,
        "adaptive/narayana/laurent",
        Box::new(|| {
            let rows = rational_rows(&sequences::narayana_polynomials_bigint(14)[1..]);
            adaptive_workload(&rows, &search_options(2, 2, 1, 0).with_laurent(true))
        }),
    );
    add(
        &mut out,
        "adaptive/narayana/default",
        Box::new(|| {
            let rows = rational_rows(&sequences::narayana_polynomials_bigint(14)[1..]);
            adaptive_workload(&rows, &AdaptiveSearchOptions::default())
        }),
    );
    add(
        &mut out,
        "adaptive/type-b-eulerian/derivative",
        Box::new(|| {
            let rows = rational_rows(&sequences::type_b_eulerian_polynomials_bigint(13));
            adaptive_workload(&rows, &search_options(2, 3, 2, 2))
        }),
    );
    add(
        &mut out,
        "adaptive/type-b-eulerian/default",
        Box::new(|| {
            let rows = rational_rows(&sequences::type_b_eulerian_polynomials_bigint(13));
            adaptive_workload(&rows, &AdaptiveSearchOptions::default())
        }),
    );
    add(
        &mut out,
        "adaptive/shifted-legendre/laurent",
        Box::new(|| {
            let rows = rational_rows(&sequences::shifted_legendre_polynomials_bigint(14)[1..]);
            adaptive_workload(&rows, &search_options(2, 2, 1, 0).with_laurent(true))
        }),
    );
    add(
        &mut out,
        "adaptive/shifted-legendre/laurent-wide",
        Box::new(|| {
            let rows = rational_rows(&sequences::shifted_legendre_polynomials_bigint(16)[1..]);
            adaptive_workload(&rows, &search_options(3, 3, 2, 1).with_laurent(true))
        }),
    );
    add(
        &mut out,
        "adaptive/hermite/denominator",
        Box::new(|| {
            let rows = rational_rows(&sequences::hermite_polynomials_bigint(12));
            let mut options = search_options(3, 3, 3, 1);
            options.try_denominator = true;
            adaptive_workload(&rows, &options)
        }),
    );
    add(
        &mut out,
        "adaptive/no-recurrence/squares",
        Box::new(|| {
            // Polynomials with no short recurrence of the searched shape.
            let rows: Vec<Vec<Q>> = (0..12i64)
                .map(|n| vec![q(n * n * n + 1), q((n * 7919) % 13), q(n % 3)])
                .collect();
            adaptive_workload(&rows, &search_options(2, 2, 2, 1))
        }),
    );

    add(
        &mut out,
        "hadamard-lift/av321",
        Box::new(|| hadamard_workload(AV321_DESCENTS)),
    );
    add(
        &mut out,
        "hadamard-lift/av123",
        Box::new(|| hadamard_workload(AV123_PEAKS)),
    );

    add(
        &mut out,
        "solve/integer 50x32 rank 18",
        Box::new(|| solve_workload(1, 50, 32, 18, 1, false)),
    );
    add(
        &mut out,
        "solve/sparse 80x48 rank 26",
        Box::new(|| solve_workload(2, 80, 48, 26, 3, false)),
    );
    add(
        &mut out,
        "solve/rational 30x24 rank 14",
        Box::new(|| solve_workload(3, 30, 24, 14, 1, true)),
    );
    add(
        &mut out,
        "solve/square 40",
        Box::new(|| square_solve_workload(4, 40)),
    );
    add(
        &mut out,
        "kernel/integer 24x18 rank 10 with budgets",
        Box::new(|| kernel_workload(5, 24, 18, 10, 1, false)),
    );
    add(
        &mut out,
        "kernel/sparse 30x22 rank 12 with budgets",
        Box::new(|| kernel_workload(6, 30, 22, 12, 3, false)),
    );
    add(
        &mut out,
        "kernel/rational 18x15 rank 8 with budgets",
        Box::new(|| kernel_workload(7, 18, 15, 8, 1, true)),
    );
    add(
        &mut out,
        "kernel/integer 60x40 rank 22",
        Box::new(|| kernel_only_workload(8, 60, 40, 22, 1, false)),
    );
    add(
        &mut out,
        "kernel/sparse 90x56 rank 30",
        Box::new(|| kernel_only_workload(9, 90, 56, 30, 4, false)),
    );
    add(
        &mut out,
        "kernel/zero rows 0x5",
        Box::new(|| kernel_only_workload(10, 0, 5, 0, 1, false)),
    );
    add(
        &mut out,
        "kernel/eulerian hankel order 3 width 9",
        Box::new(|| hankel_kernel_workload(&sequences::eulerian_polynomials_bigint(30), 3, 9)),
    );
    add(
        &mut out,
        "kernel/type-b-eulerian hankel order 3 width 10",
        Box::new(|| {
            hankel_kernel_workload(&sequences::type_b_eulerian_polynomials_bigint(30), 3, 10)
        }),
    );
    out
}

fn recorded() -> BTreeMap<String, String> {
    std::fs::read_to_string(FINGERPRINTS)
        .unwrap_or_default()
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .filter_map(|line| {
            let (name, value) = line.split_once('\t')?;
            Some((name.to_string(), value.to_string()))
        })
        .collect()
}

/// Workloads that take seconds in a release build; checked by the ignored
/// full comparison only.
fn is_slow(name: &str) -> bool {
    name.starts_with("hadamard-lift/") || name.contains("default") || name.contains("wide")
}

fn compare_fingerprints(include_slow: bool) {
    let computed: BTreeMap<String, String> = workloads()
        .into_iter()
        .filter(|(name, _)| include_slow || !is_slow(name))
        .map(|(name, run)| (name, fingerprint(&run())))
        .collect();
    if include_slow && std::env::var_os("POLYTOOL_BLESS_LINALG_WORKLOADS").is_some() {
        let mut text = String::from(
            "# Fingerprints (FNV-1a 64 of the Debug output, and its length) of\n\
             # tests/linalg_workloads.rs; re-bless only for intentional changes.\n",
        );
        for (name, value) in &computed {
            text.push_str(&format!("{name}\t{value}\n"));
        }
        std::fs::write(FINGERPRINTS, text).unwrap();
        return;
    }
    let expected: BTreeMap<String, String> = recorded()
        .into_iter()
        .filter(|(name, _)| include_slow || !is_slow(name))
        .collect();
    assert_eq!(computed, expected);
}

#[test]
fn linalg_workload_outputs_are_unchanged() {
    compare_fingerprints(false);
}

/// Also the slow workloads; blesses with `POLYTOOL_BLESS_LINALG_WORKLOADS=1`.
#[test]
#[ignore = "slow in debug builds: cargo test --release --test linalg_workloads -- --ignored"]
fn every_linalg_workload_output_is_unchanged() {
    compare_fingerprints(true);
}

#[test]
#[ignore = "timing run: cargo test --release --test linalg_workloads -- --ignored --nocapture"]
fn time_linalg_workloads() {
    let repeat: usize = std::env::var("LINALG_WORKLOAD_REPEAT")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(3)
        .max(1);
    let only = std::env::var("LINALG_WORKLOAD_ONLY").unwrap_or_default();
    let mut total = 0.0;
    println!("workload\tmedian_ms\tmin_ms\tfingerprint");
    for (name, run) in workloads() {
        if !name.contains(&only) {
            continue;
        }
        let mut times = Vec::with_capacity(repeat);
        let mut print = String::new();
        for _ in 0..repeat {
            let start = Instant::now();
            print = fingerprint(&run());
            times.push(start.elapsed().as_secs_f64() * 1000.0);
        }
        times.sort_by(f64::total_cmp);
        let median = times[times.len() / 2];
        total += median;
        println!("{name}\t{median:.1}\t{:.1}\t{print}", times[0]);
    }
    println!("total\t{total:.1}");
}
