//! `polytool bench`: recurrence-fixture timings, benchmark comparison and
//! interlacing timings.

use crate::bench_real_rootedness::{cmd_bench_real_rootedness, print_real_rootedness_help};
use crate::input::{polynomial_degree_bigint, strip_trailing_zeros_bigint};
use crate::json_out::json_bool_option;
use crate::sequence::SequenceKind;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::recurrence::*;
use polytool::*;
use serde_json::{json, Value};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::time::Instant;

pub(crate) fn print_bench_help() {
    println!("Usage:");
    println!("  polytool bench recurrence-fixtures [options]");
    println!("  polytool bench interlacing [options]");
    println!("  polytool bench real-rootedness [options]");
    println!("  polytool bench --help");
    println!();
    println!("Run built-in timing suites and print tab-separated results.");
    println!();
    println!("Subcommands:");
    println!("  recurrence-fixtures   Time adaptive recurrence search on fixture rows");
    println!("  compare               Compare two recurrence-fixture JSON benchmark runs");
    println!("  interlacing           Time consecutive interlacing checks on a sequence");
    println!("  real-rootedness       Verdict-checked timing suite on degree 20-80 families");
    println!();
    println!("Recurrence fixture options:");
    println!("  --base <dir>              Fixture directory containing manifest.tsv");
    println!("  --only <substring>        Run only fixture slugs containing substring");
    println!("  --repeat <n>              Repeat each fixture n times (default: 1)");
    println!("  --no-modular-prefilter    Disable default modular recurrence prefilter");
    println!("  --summary                 Append fixture and category summaries to stdout");
    println!("  --report <path.md>        Write a Markdown benchmark report");
    println!("  --format <tsv|json>       Output format (default: tsv)");
    println!("  --json                    Alias for --format json");
    println!("                            JSON includes adaptive counters and timing buckets");
    println!();
    println!("Compare options:");
    println!("  polytool bench compare [options] <old.json> <new.json>");
    println!("  --top <n>                 Number of worst regressions to show (default: 10)");
    println!("  --format <tsv|json>       Output format (default: tsv)");
    println!();
    println!("Interlacing options:");
    println!("  --sequence <name>         Sequence name (default: eulerian)");
    println!("  --max-n <n>               Maximum n/index to generate (default: 20)");
    println!("  --repeat <n>              Repeat each pair n times (default: 3)");
    println!();
    print_real_rootedness_help();
}

#[derive(Debug)]
struct BenchRecurrenceOptions {
    base: PathBuf,
    repeat: usize,
    modular_prefilter: bool,
    only: Option<String>,
    summary: bool,
    report: Option<PathBuf>,
    format: BenchOutputFormat,
}

impl Default for BenchRecurrenceOptions {
    fn default() -> Self {
        Self {
            base: Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("fixtures")
                .join("recurrence-benchmarks"),
            repeat: 1,
            modular_prefilter: true,
            only: None,
            summary: false,
            report: None,
            format: BenchOutputFormat::Tsv,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum BenchOutputFormat {
    Tsv,
    Json,
}

impl BenchOutputFormat {
    pub(crate) fn parse(value: &str) -> Result<Self, String> {
        match value {
            "tsv" | "text" => Ok(Self::Tsv),
            "json" => Ok(Self::Json),
            other => Err(format!("unknown bench output format `{other}`")),
        }
    }
}

#[derive(Debug)]
struct RecurrenceFixtureManifestRow {
    slug: String,
    title: String,
    description: String,
    suggested_args: String,
    rows_file: String,
}

#[derive(Debug)]
struct BenchRecurrenceRun {
    slug: String,
    title: String,
    description: String,
    categories: Vec<String>,
    run: usize,
    found: bool,
    elapsed_ms: f64,
    candidates: usize,
    unknowns: usize,
    weighted_unknowns: usize,
    fit_rows: usize,
    verify_rows: usize,
    recurrence: String,
    diagnostics: Option<AdaptiveSearchDiagnostics>,
}

#[derive(Debug)]
struct BenchRecurrenceFixtureSummary {
    slug: String,
    title: String,
    description: String,
    categories: Vec<String>,
    runs: usize,
    found_runs: usize,
    min_ms: f64,
    median_ms: f64,
    mean_ms: f64,
    max_ms: f64,
    candidates: usize,
    unknowns: usize,
    weighted_unknowns: usize,
    fit_rows: usize,
    verify_rows: usize,
    recurrence: String,
}

#[derive(Debug)]
struct BenchRecurrenceCategorySummary {
    category: String,
    fixtures: usize,
    runs: usize,
    found_runs: usize,
    mean_ms: f64,
    max_ms: f64,
}

pub(crate) fn cmd_bench(args: &[String]) {
    let Some((subcommand, rest)) = args.split_first() else {
        print_bench_help();
        std::process::exit(1);
    };
    match subcommand.as_str() {
        "recurrence-fixtures" | "recurrence" => cmd_bench_recurrence_fixtures(rest),
        "compare" => cmd_bench_compare(rest),
        "interlacing" => cmd_bench_interlacing(rest),
        "real-rootedness" | "real-rooted" => cmd_bench_real_rootedness(rest),
        _ => {
            eprintln!("unknown bench subcommand: {subcommand}");
            print_bench_help();
            std::process::exit(1);
        }
    }
}

fn parse_bench_recurrence_options(args: &[String]) -> Result<BenchRecurrenceOptions, String> {
    let mut options = BenchRecurrenceOptions::default();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--base" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--base expects a directory".to_string())?;
                options.base = PathBuf::from(value);
            }
            "--only" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--only expects a substring".to_string())?;
                options.only = Some(value.clone());
            }
            "--repeat" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--repeat expects a positive integer".to_string())?;
                options.repeat = value
                    .parse::<usize>()
                    .map_err(|_| format!("--repeat expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            "--no-modular-prefilter" => options.modular_prefilter = false,
            "--summary" => options.summary = true,
            "--report" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--report expects a path".to_string())?;
                options.report = Some(PathBuf::from(value));
            }
            "--format" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--format expects tsv or json".to_string())?;
                options.format = BenchOutputFormat::parse(value)?;
            }
            "--json" => options.format = BenchOutputFormat::Json,
            other => return Err(format!("unknown recurrence benchmark option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn parse_recurrence_fixture_manifest(
    base: &Path,
) -> Result<Vec<RecurrenceFixtureManifestRow>, String> {
    let manifest_path = base.join("manifest.tsv");
    let manifest = fs::read_to_string(&manifest_path)
        .map_err(|e| format!("failed to read `{}`: {e}", manifest_path.display()))?;
    manifest
        .lines()
        .skip(1)
        .map(|line| {
            let cols = line.split('\t').collect::<Vec<_>>();
            if cols.len() != 7 {
                return Err(format!("manifest row should have 7 columns: {line}"));
            }
            Ok(RecurrenceFixtureManifestRow {
                slug: cols[0].to_string(),
                title: cols[1].to_string(),
                description: cols[2].to_string(),
                suggested_args: cols[3].to_string(),
                rows_file: cols[5].to_string(),
            })
        })
        .collect()
}

fn read_recurrence_fixture_rows(path: &Path) -> Result<Vec<Vec<RecurrenceBigRational>>, String> {
    let input = fs::read_to_string(path)
        .map_err(|e| format!("failed to read `{}`: {e}", path.display()))?;
    input
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            line.split(',')
                .map(|coeff| {
                    parse_rational_coeff(coeff.trim())
                        .map_err(|e| format!("failed to parse `{}`: {e}", line.trim()))
                })
                .collect()
        })
        .collect()
}

fn parse_fixture_search_options(
    suggested_args: &str,
    modular_prefilter: bool,
) -> Result<AdaptiveSearchOptions, String> {
    let mut search = AdaptiveSearchOptions {
        modular_prefilter,
        ..Default::default()
    };
    let mut args = suggested_args.split_whitespace();
    while let Some(arg) = args.next() {
        match arg {
            "--skip-prefix" => search.skip_prefix = parse_fixture_usize(arg, args.next())?,
            "--min-rec-len" => search.min_rec_len = parse_fixture_usize(arg, args.next())?,
            "--max-rec-len" => search.max_rec_len = parse_fixture_usize(arg, args.next())?,
            "--min-var-deg" => search.min_var_deg = parse_fixture_usize(arg, args.next())?,
            "--max-var-deg" => search.max_var_deg = parse_fixture_usize(arg, args.next())?,
            "--min-idx-deg" => search.min_idx_deg = parse_fixture_usize(arg, args.next())?,
            "--max-idx-deg" => search.max_idx_deg = parse_fixture_usize(arg, args.next())?,
            "--min-diff-deg" => search.min_diff_deg = parse_fixture_usize(arg, args.next())?,
            "--max-diff-deg" => search.max_diff_deg = parse_fixture_usize(arg, args.next())?,
            "--inhomogeneous" => search.try_inhomogeneous = true,
            "--min-inhomo-var-deg" => {
                search.min_inhomo_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-inhomo-var-deg" => {
                search.max_inhomo_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--min-inhomo-idx-deg" => {
                search.min_inhomo_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-inhomo-idx-deg" => {
                search.max_inhomo_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--denominator" | "--try-denominator" => search.try_denominator = true,
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--max-denom-idx-deg" => {
                search.try_denominator = true;
                search.max_denom_idx_deg = parse_fixture_usize(arg, args.next())?;
            }
            "--alternating-sign" => search.try_alternating_sign = true,
            "--min-margin" => search.min_margin = parse_fixture_usize(arg, args.next())?,
            "--fit-extra-rows" => search.fit_extra_rows = parse_fixture_usize(arg, args.next())?,
            "--no-verify" => search.no_verify = true,
            other => return Err(format!("unsupported fixture search option `{other}`")),
        }
    }
    Ok(search)
}

fn recurrence_fixture_categories(
    fixture: &RecurrenceFixtureManifestRow,
    search: &AdaptiveSearchOptions,
) -> Vec<String> {
    let mut categories = vec![if fixture.slug.contains("_oeis_") {
        "oeis".to_string()
    } else {
        "synthetic".to_string()
    }];
    if search.skip_prefix > 0 {
        categories.push("skip-prefix".to_string());
    }
    if search.max_rec_len >= 4 {
        categories.push("high-order".to_string());
    }
    if search.max_diff_deg > 0 {
        categories.push("derivative".to_string());
    }
    if search.max_diff_deg >= 2 {
        categories.push("higher-derivative".to_string());
    }
    if search.max_idx_deg > 0
        || (search.try_denominator && search.max_denom_idx_deg > 0)
        || (search.try_inhomogeneous && search.max_inhomo_idx_deg > 0)
    {
        categories.push("index-dependent".to_string());
    }
    if search.try_denominator {
        categories.push("denominator".to_string());
    }
    if search.try_inhomogeneous {
        categories.push("inhomogeneous".to_string());
    }
    if search.try_alternating_sign {
        categories.push("alternating-sign".to_string());
    }
    if fixture.slug.contains("closed_form") {
        categories.push("closed-form".to_string());
    }
    categories
}

fn bench_mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}

fn bench_median(values: &mut [f64]) -> f64 {
    values.sort_by(|a, b| a.total_cmp(b));
    let middle = values.len() / 2;
    if values.len() & 1 == 0 {
        (values[middle - 1] + values[middle]) / 2.0
    } else {
        values[middle]
    }
}

fn markdown_cell(text: &str) -> String {
    text.replace('\n', " ").replace('|', "\\|")
}

fn recurrence_bench_fixture_summaries(
    records: &[BenchRecurrenceRun],
) -> Vec<BenchRecurrenceFixtureSummary> {
    let mut by_fixture: BTreeMap<&str, Vec<&BenchRecurrenceRun>> = BTreeMap::new();
    for record in records {
        by_fixture.entry(&record.slug).or_default().push(record);
    }

    by_fixture
        .into_iter()
        .map(|(slug, runs)| {
            let found_runs = runs.iter().filter(|record| record.found).count();
            let mut elapsed = runs
                .iter()
                .map(|record| record.elapsed_ms)
                .collect::<Vec<_>>();
            let min_ms = elapsed.iter().copied().fold(f64::INFINITY, f64::min);
            let max_ms = elapsed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mean_ms = bench_mean(&elapsed);
            let median_ms = bench_median(&mut elapsed);
            let sample = runs[0];
            BenchRecurrenceFixtureSummary {
                slug: slug.to_string(),
                title: sample.title.clone(),
                description: sample.description.clone(),
                categories: sample.categories.clone(),
                runs: runs.len(),
                found_runs,
                min_ms,
                median_ms,
                mean_ms,
                max_ms,
                candidates: sample.candidates,
                unknowns: sample.unknowns,
                weighted_unknowns: sample.weighted_unknowns,
                fit_rows: sample.fit_rows,
                verify_rows: sample.verify_rows,
                recurrence: sample.recurrence.clone(),
            }
        })
        .collect()
}

fn recurrence_bench_category_summary_values(
    records: &[BenchRecurrenceRun],
) -> Vec<BenchRecurrenceCategorySummary> {
    recurrence_bench_category_summaries(records)
        .into_iter()
        .map(|(category, fixtures, runs)| {
            let found_runs = runs.iter().filter(|record| record.found).count();
            let elapsed = runs
                .iter()
                .map(|record| record.elapsed_ms)
                .collect::<Vec<_>>();
            let max_ms = elapsed.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            BenchRecurrenceCategorySummary {
                category,
                fixtures: fixtures.len(),
                runs: runs.len(),
                found_runs,
                mean_ms: bench_mean(&elapsed),
                max_ms,
            }
        })
        .collect()
}

fn print_recurrence_bench_summary(records: &[BenchRecurrenceRun]) {
    if records.is_empty() {
        println!("# fixture_summary");
        println!("slug\ttitle\tcategories\truns\tfound_runs\tmin_ms\tmedian_ms\tmean_ms\tmax_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows");
        println!("# category_summary");
        println!("category\tfixtures\truns\tfound_runs\tmean_ms\tmax_ms");
        return;
    }

    println!("# fixture_summary");
    println!("slug\ttitle\tcategories\truns\tfound_runs\tmin_ms\tmedian_ms\tmean_ms\tmax_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows");
    for summary in recurrence_bench_fixture_summaries(records) {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\t{}\t{}",
            summary.slug,
            summary.title,
            summary.categories.join(","),
            summary.runs,
            summary.found_runs,
            summary.min_ms,
            summary.median_ms,
            summary.mean_ms,
            summary.max_ms,
            summary.candidates,
            summary.unknowns,
            summary.weighted_unknowns,
            summary.fit_rows,
            summary.verify_rows,
        );
    }

    println!("# category_summary");
    println!("category\tfixtures\truns\tfound_runs\tmean_ms\tmax_ms");
    for summary in recurrence_bench_category_summary_values(records) {
        println!(
            "{}\t{}\t{}\t{}\t{:.3}\t{:.3}",
            summary.category,
            summary.fixtures,
            summary.runs,
            summary.found_runs,
            summary.mean_ms,
            summary.max_ms,
        );
    }
}

fn recurrence_bench_category_summaries(
    records: &[BenchRecurrenceRun],
) -> Vec<(String, BTreeSet<&str>, Vec<&BenchRecurrenceRun>)> {
    let mut categories: BTreeMap<String, (BTreeSet<&str>, Vec<&BenchRecurrenceRun>)> =
        BTreeMap::new();
    for record in records {
        for category in &record.categories {
            let entry = categories.entry(category.clone()).or_default();
            entry.0.insert(&record.slug);
            entry.1.push(record);
        }
    }
    categories
        .into_iter()
        .map(|(category, (fixtures, runs))| (category, fixtures, runs))
        .collect()
}

fn write_recurrence_bench_report(
    options: &BenchRecurrenceOptions,
    records: &[BenchRecurrenceRun],
) -> Result<(), String> {
    let Some(path) = &options.report else {
        return Ok(());
    };

    let mut report = String::new();
    report.push_str("# Recurrence Fixture Benchmark Report\n\n");
    report.push_str(&format!("- Base: `{}`\n", options.base.display()));
    report.push_str(&format!("- Repeat: `{}`\n", options.repeat));
    report.push_str(&format!(
        "- Modular prefilter: `{}`\n",
        if options.modular_prefilter {
            "enabled"
        } else {
            "disabled"
        }
    ));
    if let Some(only) = &options.only {
        report.push_str(&format!("- Filter: `{}`\n", only));
    }
    report.push_str("\n## Category Summary\n\n");
    report.push_str("| category | fixtures | runs | found | mean ms | max ms |\n");
    report.push_str("|---|---:|---:|---:|---:|---:|\n");
    for summary in recurrence_bench_category_summary_values(records) {
        report.push_str(&format!(
            "| {} | {} | {} | {} | {:.3} | {:.3} |\n",
            markdown_cell(&summary.category),
            summary.fixtures,
            summary.runs,
            summary.found_runs,
            summary.mean_ms,
            summary.max_ms,
        ));
    }

    report.push_str("\n## Fixture Summary\n\n");
    report.push_str("| slug | title | description | categories | runs | found | min ms | median ms | mean ms | max ms | unknowns | weighted | recurrence |\n");
    report.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---|\n");
    for summary in recurrence_bench_fixture_summaries(records) {
        report.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {} | {} | {} |\n",
            markdown_cell(&summary.slug),
            markdown_cell(&summary.title),
            markdown_cell(&summary.description),
            markdown_cell(&summary.categories.join(", ")),
            summary.runs,
            summary.found_runs,
            summary.min_ms,
            summary.median_ms,
            summary.mean_ms,
            summary.max_ms,
            summary.unknowns,
            summary.weighted_unknowns,
            markdown_cell(&summary.recurrence),
        ));
    }

    if records.is_empty() {
        report.push_str("\nNo fixtures matched this benchmark filter.\n");
    }

    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)
            .map_err(|e| format!("failed to create `{}`: {e}", parent.display()))?;
    }
    fs::write(path, report).map_err(|e| format!("failed to write `{}`: {e}", path.display()))
}

fn recurrence_diagnostics_json(diagnostics: &AdaptiveSearchDiagnostics) -> Value {
    json!({
        "generated_candidates": diagnostics.generated_candidates,
        "considered_candidates": diagnostics.considered_candidates,
        "insufficient_fit_rows": diagnostics.insufficient_fit_rows,
        "equation_bound_rejections": diagnostics.equation_bound_rejections,
        "degree_bound_rejections": diagnostics.degree_bound_rejections,
        "structural_zero_rejections": diagnostics.structural_zero_rejections,
        "modular_prefilter_rejections": diagnostics.modular_prefilter_rejections,
        "exact_solve_attempts": diagnostics.exact_solve_attempts,
        "modular_lift_attempts": diagnostics.modular_lift_attempts,
        "modular_lift_successes": diagnostics.modular_lift_successes,
        "modular_lift_failures": diagnostics.modular_lift_failures,
        "failed_exact_solves": diagnostics.failed_exact_solves,
        "heldout_verification_failures": diagnostics.heldout_verification_failures,
        "denominator_escalation_entered": diagnostics.denominator_escalation_entered,
        "derivative_precompute_ms": diagnostics.derivative_precompute_ms,
        "modular_cache_build_ms": diagnostics.modular_cache_build_ms,
        "candidate_generation_ms": diagnostics.candidate_generation_ms,
        "fit_selection_ms": diagnostics.fit_selection_ms,
        "equation_count_ms": diagnostics.equation_count_ms,
        "degree_bound_ms": diagnostics.degree_bound_ms,
        "structural_zero_ms": diagnostics.structural_zero_ms,
        "modular_prefilter_ms": diagnostics.modular_prefilter_ms,
        "modular_lift_ms": diagnostics.modular_lift_ms,
        "exact_solve_ms": diagnostics.exact_solve_ms,
        "heldout_verify_ms": diagnostics.heldout_verify_ms,
    })
}

fn recurrence_bench_run_json(record: &BenchRecurrenceRun) -> Value {
    json!({
        "slug": &record.slug,
        "title": &record.title,
        "description": &record.description,
        "categories": &record.categories,
        "run": record.run,
        "found": record.found,
        "elapsed_ms": record.elapsed_ms,
        "candidates": record.candidates,
        "unknowns": record.unknowns,
        "weighted_unknowns": record.weighted_unknowns,
        "fit_rows": record.fit_rows,
        "verify_rows": record.verify_rows,
        "recurrence": &record.recurrence,
        "diagnostics": record
            .diagnostics
            .as_ref()
            .map(recurrence_diagnostics_json)
            .unwrap_or(Value::Null),
    })
}

fn recurrence_bench_fixture_summary_json(summary: &BenchRecurrenceFixtureSummary) -> Value {
    json!({
        "slug": &summary.slug,
        "title": &summary.title,
        "description": &summary.description,
        "categories": &summary.categories,
        "runs": summary.runs,
        "found_runs": summary.found_runs,
        "min_ms": summary.min_ms,
        "median_ms": summary.median_ms,
        "mean_ms": summary.mean_ms,
        "max_ms": summary.max_ms,
        "candidates": summary.candidates,
        "unknowns": summary.unknowns,
        "weighted_unknowns": summary.weighted_unknowns,
        "fit_rows": summary.fit_rows,
        "verify_rows": summary.verify_rows,
        "recurrence": &summary.recurrence,
    })
}

fn recurrence_bench_category_summary_json(summary: &BenchRecurrenceCategorySummary) -> Value {
    json!({
        "category": &summary.category,
        "fixtures": summary.fixtures,
        "runs": summary.runs,
        "found_runs": summary.found_runs,
        "mean_ms": summary.mean_ms,
        "max_ms": summary.max_ms,
    })
}

fn recurrence_bench_json(
    options: &BenchRecurrenceOptions,
    records: &[BenchRecurrenceRun],
) -> Value {
    let fixture_summaries = recurrence_bench_fixture_summaries(records);
    let category_summaries = recurrence_bench_category_summary_values(records);
    json!({
        "schema": "polytool.bench.recurrence-fixtures.v1",
        "command": "recurrence-fixtures",
        "options": {
            "base": options.base.display().to_string(),
            "repeat": options.repeat,
            "modular_prefilter": options.modular_prefilter,
            "only": &options.only,
        },
        "runs": records.iter().map(recurrence_bench_run_json).collect::<Vec<_>>(),
        "fixture_summaries": fixture_summaries
            .iter()
            .map(recurrence_bench_fixture_summary_json)
            .collect::<Vec<_>>(),
        "category_summaries": category_summaries
            .iter()
            .map(recurrence_bench_category_summary_json)
            .collect::<Vec<_>>(),
    })
}

#[derive(Debug)]
struct BenchCompareOptions {
    old_path: PathBuf,
    new_path: PathBuf,
    top: usize,
    format: BenchOutputFormat,
}

#[derive(Debug, Clone)]
struct BenchComparePoint {
    mean_ms: f64,
    found_runs: usize,
    runs: usize,
}

#[derive(Debug)]
struct BenchCompareRow {
    kind: String,
    name: String,
    old_mean_ms: Option<f64>,
    new_mean_ms: Option<f64>,
    speedup: Option<f64>,
    delta_ms: Option<f64>,
    old_found_runs: Option<usize>,
    old_runs: Option<usize>,
    new_found_runs: Option<usize>,
    new_runs: Option<usize>,
    status: String,
}

fn parse_bench_compare_options(args: &[String]) -> Result<BenchCompareOptions, String> {
    let mut top = 10;
    let mut format = BenchOutputFormat::Tsv;
    let mut paths = Vec::new();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--top" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--top expects a positive integer".to_string())?;
                top = value
                    .parse::<usize>()
                    .map_err(|_| format!("--top expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            "--format" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--format expects tsv or json".to_string())?;
                format = BenchOutputFormat::parse(value)?;
            }
            "--json" => format = BenchOutputFormat::Json,
            other if other.starts_with('-') => {
                return Err(format!("unknown benchmark compare option: {other}"));
            }
            path => paths.push(PathBuf::from(path)),
        }
        i += 1;
    }

    if paths.len() != 2 {
        return Err("bench compare expects exactly two JSON paths".to_string());
    }
    Ok(BenchCompareOptions {
        old_path: paths.remove(0),
        new_path: paths.remove(0),
        top,
        format,
    })
}

fn read_benchmark_json(path: &Path) -> Result<Value, String> {
    let text = fs::read_to_string(path)
        .map_err(|e| format!("failed to read benchmark JSON `{}`: {e}", path.display()))?;
    let value: Value = serde_json::from_str(&text)
        .map_err(|e| format!("failed to parse benchmark JSON `{}`: {e}", path.display()))?;
    let schema = value
        .get("schema")
        .and_then(Value::as_str)
        .ok_or_else(|| format!("benchmark JSON `{}` has no schema", path.display()))?;
    if schema != "polytool.bench.recurrence-fixtures.v1" {
        return Err(format!(
            "unsupported benchmark JSON schema `{schema}` in `{}`",
            path.display()
        ));
    }
    Ok(value)
}

fn value_usize(value: &Value, key: &str) -> Result<usize, String> {
    value
        .get(key)
        .and_then(Value::as_u64)
        .and_then(|n| usize::try_from(n).ok())
        .ok_or_else(|| format!("benchmark summary entry is missing integer `{key}`"))
}

fn value_f64(value: &Value, key: &str) -> Result<f64, String> {
    value
        .get(key)
        .and_then(Value::as_f64)
        .ok_or_else(|| format!("benchmark summary entry is missing number `{key}`"))
}

fn value_string(value: &Value, key: &str) -> Result<String, String> {
    value
        .get(key)
        .and_then(Value::as_str)
        .map(ToString::to_string)
        .ok_or_else(|| format!("benchmark summary entry is missing string `{key}`"))
}

fn extract_benchmark_points(
    report: &Value,
    array_key: &str,
    name_key: &str,
) -> Result<BTreeMap<String, BenchComparePoint>, String> {
    let entries = report
        .get(array_key)
        .and_then(Value::as_array)
        .ok_or_else(|| format!("benchmark JSON missing array `{array_key}`"))?;
    let mut map = BTreeMap::new();
    for entry in entries {
        let name = value_string(entry, name_key)?;
        map.insert(
            name,
            BenchComparePoint {
                mean_ms: value_f64(entry, "mean_ms")?,
                found_runs: value_usize(entry, "found_runs")?,
                runs: value_usize(entry, "runs")?,
            },
        );
    }
    Ok(map)
}

fn compare_status(
    old: Option<&BenchComparePoint>,
    new: Option<&BenchComparePoint>,
    speedup: Option<f64>,
) -> String {
    match (old, new) {
        (None, Some(_)) => "new".to_string(),
        (Some(_), None) => "removed".to_string(),
        (None, None) => "missing".to_string(),
        (Some(old), Some(new)) if old.found_runs > 0 && new.found_runs == 0 => "lost".to_string(),
        (Some(old), Some(new)) if old.found_runs == 0 && new.found_runs > 0 => {
            "newly_found".to_string()
        }
        (Some(_), Some(_)) => {
            let speedup = speedup.unwrap_or(1.0);
            if speedup >= 1.10 {
                "faster".to_string()
            } else if speedup <= 0.90 {
                "slower".to_string()
            } else {
                "flat".to_string()
            }
        }
    }
}

fn benchmark_compare_rows(
    kind: &str,
    old: &BTreeMap<String, BenchComparePoint>,
    new: &BTreeMap<String, BenchComparePoint>,
) -> Vec<BenchCompareRow> {
    let names = old
        .keys()
        .chain(new.keys())
        .cloned()
        .collect::<BTreeSet<_>>();
    names
        .into_iter()
        .map(|name| {
            let old_point = old.get(&name);
            let new_point = new.get(&name);
            let speedup = match (old_point, new_point) {
                (Some(old), Some(new)) if new.mean_ms > 0.0 => Some(old.mean_ms / new.mean_ms),
                _ => None,
            };
            let delta_ms = match (old_point, new_point) {
                (Some(old), Some(new)) => Some(new.mean_ms - old.mean_ms),
                _ => None,
            };
            BenchCompareRow {
                kind: kind.to_string(),
                name,
                old_mean_ms: old_point.map(|point| point.mean_ms),
                new_mean_ms: new_point.map(|point| point.mean_ms),
                speedup,
                delta_ms,
                old_found_runs: old_point.map(|point| point.found_runs),
                old_runs: old_point.map(|point| point.runs),
                new_found_runs: new_point.map(|point| point.found_runs),
                new_runs: new_point.map(|point| point.runs),
                status: compare_status(old_point, new_point, speedup),
            }
        })
        .collect()
}

fn format_optional_f64(value: Option<f64>) -> String {
    value.map(|value| format!("{value:.3}")).unwrap_or_default()
}

fn format_optional_usize(value: Option<usize>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn print_compare_rows(section: &str, rows: &[&BenchCompareRow]) {
    println!("# {section}");
    println!(
        "kind\tname\told_mean_ms\tnew_mean_ms\tspeedup\tdelta_ms\told_found\told_runs\tnew_found\tnew_runs\tstatus"
    );
    for row in rows {
        println!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.kind,
            row.name,
            format_optional_f64(row.old_mean_ms),
            format_optional_f64(row.new_mean_ms),
            format_optional_f64(row.speedup),
            format_optional_f64(row.delta_ms),
            format_optional_usize(row.old_found_runs),
            format_optional_usize(row.old_runs),
            format_optional_usize(row.new_found_runs),
            format_optional_usize(row.new_runs),
            row.status,
        );
    }
}

fn compare_row_json(row: &BenchCompareRow) -> Value {
    json!({
        "kind": &row.kind,
        "name": &row.name,
        "old_mean_ms": row.old_mean_ms,
        "new_mean_ms": row.new_mean_ms,
        "speedup": row.speedup,
        "delta_ms": row.delta_ms,
        "old_found_runs": row.old_found_runs,
        "old_runs": row.old_runs,
        "new_found_runs": row.new_found_runs,
        "new_runs": row.new_runs,
        "status": &row.status,
    })
}

fn cmd_bench_compare(args: &[String]) {
    let options = match parse_bench_compare_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let old_report = match read_benchmark_json(&options.old_path) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let new_report = match read_benchmark_json(&options.new_path) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };

    let old_fixtures = match extract_benchmark_points(&old_report, "fixture_summaries", "slug") {
        Ok(points) => points,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let new_fixtures = match extract_benchmark_points(&new_report, "fixture_summaries", "slug") {
        Ok(points) => points,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let old_categories =
        match extract_benchmark_points(&old_report, "category_summaries", "category") {
            Ok(points) => points,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };
    let new_categories =
        match extract_benchmark_points(&new_report, "category_summaries", "category") {
            Ok(points) => points,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };

    let fixture_rows = benchmark_compare_rows("fixture", &old_fixtures, &new_fixtures);
    let category_rows = benchmark_compare_rows("category", &old_categories, &new_categories);
    let mut worst_regressions = fixture_rows
        .iter()
        .filter(|row| row.old_mean_ms.is_some() && row.new_mean_ms.is_some())
        .collect::<Vec<_>>();
    worst_regressions.sort_by(|a, b| {
        let a_ratio = a.new_mean_ms.unwrap_or(0.0) / a.old_mean_ms.unwrap_or(1.0);
        let b_ratio = b.new_mean_ms.unwrap_or(0.0) / b.old_mean_ms.unwrap_or(1.0);
        b_ratio.total_cmp(&a_ratio)
    });
    worst_regressions.truncate(options.top);

    if options.format == BenchOutputFormat::Json {
        let output = json!({
            "schema": "polytool.bench.compare.v1",
            "old": options.old_path.display().to_string(),
            "new": options.new_path.display().to_string(),
            "fixture_compare": fixture_rows.iter().map(compare_row_json).collect::<Vec<_>>(),
            "category_compare": category_rows.iter().map(compare_row_json).collect::<Vec<_>>(),
            "worst_regressions": worst_regressions
                .iter()
                .map(|row| compare_row_json(row))
                .collect::<Vec<_>>(),
        });
        println!("{}", serde_json::to_string_pretty(&output).unwrap());
    } else {
        print_compare_rows("fixture_compare", &fixture_rows.iter().collect::<Vec<_>>());
        print_compare_rows(
            "category_compare",
            &category_rows.iter().collect::<Vec<_>>(),
        );
        print_compare_rows("worst_regressions", &worst_regressions);
    }
}

fn parse_fixture_usize(flag: &str, value: Option<&str>) -> Result<usize, String> {
    let value = value.ok_or_else(|| format!("{flag} expects a value"))?;
    value
        .parse()
        .map_err(|_| format!("{flag} expects a nonnegative integer, got '{value}'"))
}

fn cmd_bench_recurrence_fixtures(args: &[String]) {
    let options = match parse_bench_recurrence_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let manifest = match parse_recurrence_fixture_manifest(&options.base) {
        Ok(manifest) => manifest,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    if options.format == BenchOutputFormat::Tsv {
        println!(
            "slug\trun\tfound\telapsed_ms\tcandidates\tunknowns\tweighted\tfit_rows\tverify_rows\trecurrence"
        );
    }
    let mut records = Vec::new();
    for fixture in manifest {
        if options
            .only
            .as_ref()
            .is_some_and(|needle| !fixture.slug.contains(needle))
        {
            continue;
        }
        let rows_path = options.base.join(&fixture.rows_file);
        let rows = match read_recurrence_fixture_rows(&rows_path) {
            Ok(rows) => rows,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(1);
            }
        };
        let search = match parse_fixture_search_options(
            &fixture.suggested_args,
            options.modular_prefilter,
        ) {
            Ok(search) => search,
            Err(error) => {
                eprintln!("{}: {error}", fixture.slug);
                std::process::exit(1);
            }
        };
        let categories = recurrence_fixture_categories(&fixture, &search);
        for run in 1..=options.repeat {
            let start = Instant::now();
            let result = find_recurrence_adaptive_rational(black_box(&rows), black_box(&search));
            let elapsed_ms = start.elapsed().as_secs_f64() * 1000.0;
            if let Some(result) = result {
                let recurrence = result.recurrence.to_string();
                if options.format == BenchOutputFormat::Tsv {
                    println!(
                        "{}\t{}\ttrue\t{:.3}\t{}\t{}\t{}\t{}\t{}\t{}",
                        fixture.slug,
                        run,
                        elapsed_ms,
                        result.candidates_tried,
                        result.num_unknowns,
                        result.weighted_unknowns,
                        result.fit_polynomials,
                        result.verification_polynomials,
                        recurrence,
                    );
                }
                records.push(BenchRecurrenceRun {
                    slug: fixture.slug.clone(),
                    title: fixture.title.clone(),
                    description: fixture.description.clone(),
                    categories: categories.clone(),
                    run,
                    found: true,
                    elapsed_ms,
                    candidates: result.candidates_tried,
                    unknowns: result.num_unknowns,
                    weighted_unknowns: result.weighted_unknowns,
                    fit_rows: result.fit_polynomials,
                    verify_rows: result.verification_polynomials,
                    recurrence,
                    diagnostics: Some(result.diagnostics),
                });
            } else {
                if options.format == BenchOutputFormat::Tsv {
                    println!(
                        "{}\t{}\tfalse\t{:.3}\t0\t0\t0\t0\t0\t",
                        fixture.slug, run, elapsed_ms
                    );
                }
                records.push(BenchRecurrenceRun {
                    slug: fixture.slug.clone(),
                    title: fixture.title.clone(),
                    description: fixture.description.clone(),
                    categories: categories.clone(),
                    run,
                    found: false,
                    elapsed_ms,
                    candidates: 0,
                    unknowns: 0,
                    weighted_unknowns: 0,
                    fit_rows: 0,
                    verify_rows: 0,
                    recurrence: String::new(),
                    diagnostics: None,
                });
            }
        }
    }
    if options.format == BenchOutputFormat::Tsv && options.summary {
        print_recurrence_bench_summary(&records);
    }
    if let Err(error) = write_recurrence_bench_report(&options, &records) {
        eprintln!("{error}");
        std::process::exit(1);
    }
    if options.format == BenchOutputFormat::Json {
        println!(
            "{}",
            serde_json::to_string_pretty(&recurrence_bench_json(&options, &records)).unwrap()
        );
    }
}

#[derive(Debug)]
struct BenchInterlacingOptions {
    sequence: SequenceKind,
    max_n: usize,
    repeat: usize,
}

impl Default for BenchInterlacingOptions {
    fn default() -> Self {
        Self {
            sequence: SequenceKind::Eulerian,
            max_n: 20,
            repeat: 3,
        }
    }
}

fn parse_bench_interlacing_options(args: &[String]) -> Result<BenchInterlacingOptions, String> {
    let mut options = BenchInterlacingOptions::default();
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--sequence" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--sequence expects a sequence name".to_string())?;
                options.sequence = SequenceKind::parse(value)
                    .ok_or_else(|| format!("unknown sequence: {value}"))?;
            }
            "--max-n" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--max-n expects a nonnegative integer".to_string())?;
                options.max_n = value
                    .parse()
                    .map_err(|_| format!("--max-n expects a nonnegative integer, got '{value}'"))?;
            }
            "--repeat" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--repeat expects a positive integer".to_string())?;
                options.repeat = value
                    .parse::<usize>()
                    .map_err(|_| format!("--repeat expects a positive integer, got '{value}'"))?
                    .max(1);
            }
            other => return Err(format!("unknown interlacing benchmark option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn cmd_bench_interlacing(args: &[String]) {
    let options = match parse_bench_interlacing_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_bench_help();
            std::process::exit(1);
        }
    };
    let polynomials = options.sequence.polynomials(options.max_n);
    println!("sequence\tleft_index\tright_index\tdegree\trepeat\tavg_us\tresult");
    for right in 1..polynomials.len() {
        let left = right - 1;
        let p = strip_trailing_zeros_bigint(&polynomials[left]);
        let q = strip_trailing_zeros_bigint(&polynomials[right]);
        let mut result = None;
        let start = Instant::now();
        for _ in 0..options.repeat {
            result = check_interlacing_bigint_coeffs(black_box(p), black_box(q));
        }
        let avg_us = start.elapsed().as_secs_f64() * 1_000_000.0 / options.repeat as f64;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}",
            options.sequence.name(),
            left,
            right,
            polynomial_degree_bigint(q),
            options.repeat,
            avg_us,
            json_bool_option(result),
        );
    }
}
