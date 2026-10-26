//! Benchmark adaptive recurrence search on the fixture suite.
//!
//! Run from the repository root:
//!
//! ```text
//! cargo run --release -p polytool --example benchmark_recurrence_fixtures
//! ```

use polytool::recurrence::{
    find_recurrence_adaptive_rational, parse_rational_coeff, AdaptiveSearchOptions, BigRational,
};
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

const BASE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/fixtures/recurrence-benchmarks"
);

#[derive(Debug)]
struct FixtureManifestRow {
    slug: String,
    suggested_args: String,
    rows_file: String,
}

#[derive(Debug)]
struct BenchOptions {
    repeat: usize,
    modular_prefilter: bool,
    only: Option<String>,
}

impl Default for BenchOptions {
    fn default() -> Self {
        Self {
            repeat: 1,
            modular_prefilter: true,
            only: None,
        }
    }
}

fn parse_bench_options() -> BenchOptions {
    let mut opts = BenchOptions::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--repeat" => {
                let value = args
                    .next()
                    .expect("--repeat requires a positive integer argument");
                opts.repeat = value
                    .parse::<usize>()
                    .expect("--repeat requires a positive integer argument")
                    .max(1);
            }
            "--no-modular-prefilter" => {
                opts.modular_prefilter = false;
            }
            "--only" => {
                opts.only = Some(args.next().expect("--only requires a slug substring"));
            }
            "-h" | "--help" => {
                print_help_and_exit();
            }
            other => {
                panic!("unknown benchmark option `{other}`");
            }
        }
    }
    opts
}

fn print_help_and_exit() -> ! {
    println!("Usage:");
    println!(
        "  cargo run --release -p polytool --example benchmark_recurrence_fixtures -- [options]"
    );
    println!();
    println!("Options:");
    println!("  --repeat <n>              Repeat each fixture n times");
    println!("  --no-modular-prefilter    Disable default modular recurrence prefilter");
    println!("  --only <substring>        Run only fixture slugs containing substring");
    println!("  -h, --help                Print this help");
    std::process::exit(0);
}

fn parse_manifest() -> Vec<FixtureManifestRow> {
    let manifest = fs::read_to_string(Path::new(BASE).join("manifest.tsv"))
        .expect("read recurrence benchmark manifest");
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(cols.len(), 7, "manifest row should have 7 columns: {line}");
            FixtureManifestRow {
                slug: cols[0].to_string(),
                suggested_args: cols[3].to_string(),
                rows_file: cols[5].to_string(),
            }
        })
        .collect()
}

fn parse_row(line: &str) -> Vec<BigRational> {
    line.split(',')
        .map(|coeff| parse_rational_coeff(coeff.trim()).expect("parse rational coefficient"))
        .collect()
}

fn read_rows(path: PathBuf) -> Vec<Vec<BigRational>> {
    fs::read_to_string(path)
        .expect("read coefficient rows")
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(parse_row)
        .collect()
}

fn parse_search_options(suggested_args: &str, modular_prefilter: bool) -> AdaptiveSearchOptions {
    let mut search = AdaptiveSearchOptions {
        modular_prefilter,
        ..Default::default()
    };
    let mut args = suggested_args.split_whitespace();
    while let Some(arg) = args.next() {
        match arg {
            "--skip-prefix" => search.skip_prefix = parse_usize_arg(arg, args.next()),
            "--min-rec-len" => search.min_rec_len = parse_usize_arg(arg, args.next()),
            "--max-rec-len" => search.max_rec_len = parse_usize_arg(arg, args.next()),
            "--min-var-deg" => search.min_var_deg = parse_usize_arg(arg, args.next()),
            "--max-var-deg" => search.max_var_deg = parse_usize_arg(arg, args.next()),
            "--min-idx-deg" => search.min_idx_deg = parse_usize_arg(arg, args.next()),
            "--max-idx-deg" => search.max_idx_deg = parse_usize_arg(arg, args.next()),
            "--min-diff-deg" => search.min_diff_deg = parse_usize_arg(arg, args.next()),
            "--max-diff-deg" => search.max_diff_deg = parse_usize_arg(arg, args.next()),
            "--inhomogeneous" => search.try_inhomogeneous = true,
            "--min-inhomo-var-deg" => {
                search.min_inhomo_var_deg = parse_usize_arg(arg, args.next());
            }
            "--max-inhomo-var-deg" => {
                search.max_inhomo_var_deg = parse_usize_arg(arg, args.next());
            }
            "--min-inhomo-idx-deg" => {
                search.min_inhomo_idx_deg = parse_usize_arg(arg, args.next());
            }
            "--max-inhomo-idx-deg" => {
                search.max_inhomo_idx_deg = parse_usize_arg(arg, args.next());
            }
            "--denominator" | "--try-denominator" => search.try_denominator = true,
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg = parse_usize_arg(arg, args.next());
            }
            "--max-denom-idx-deg" => {
                search.try_denominator = true;
                search.max_denom_idx_deg = parse_usize_arg(arg, args.next());
            }
            "--alternating-sign" => search.try_alternating_sign = true,
            "--min-margin" => search.min_margin = parse_usize_arg(arg, args.next()),
            "--fit-extra-rows" => search.fit_extra_rows = parse_usize_arg(arg, args.next()),
            "--no-verify" => search.no_verify = true,
            other => panic!("unsupported fixture search option `{other}`"),
        }
    }
    search
}

fn parse_usize_arg(flag: &str, value: Option<&str>) -> usize {
    value
        .unwrap_or_else(|| panic!("{flag} requires a value"))
        .parse()
        .unwrap_or_else(|_| panic!("{flag} requires a nonnegative integer value"))
}

fn main() {
    let bench_opts = parse_bench_options();
    println!(
        "slug\trun\tfound\telapsed_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows\trecurrence"
    );

    for fixture in parse_manifest() {
        if bench_opts
            .only
            .as_ref()
            .is_some_and(|needle| !fixture.slug.contains(needle))
        {
            continue;
        }
        let rows = read_rows(Path::new(BASE).join(&fixture.rows_file));
        let search = parse_search_options(&fixture.suggested_args, bench_opts.modular_prefilter);
        for run in 1..=bench_opts.repeat {
            let start = Instant::now();
            let result = find_recurrence_adaptive_rational(black_box(&rows), black_box(&search));
            let elapsed = start.elapsed();
            if let Some(result) = result {
                println!(
                    "{}\t{}\ttrue\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
                    fixture.slug,
                    run,
                    elapsed.as_secs_f64() * 1000.0,
                    result.candidates_tried,
                    result.num_unknowns,
                    result.weighted_unknowns,
                    result.fit_polynomials,
                    result.verification_polynomials,
                    result.recurrence,
                );
            } else {
                println!(
                    "{}\t{}\tfalse\t{:.3}\t0\t0\t0\t0\t0\t",
                    fixture.slug,
                    run,
                    elapsed.as_secs_f64() * 1000.0,
                );
            }
        }
    }
}
