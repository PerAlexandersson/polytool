//! `polytool bench real-rootedness`: a deterministic regression suite for exact
//! real-rootedness, strict interlacing, distinct-real-root counting and h\*
//! conversion on positive-coefficient families of degree 20--80.
//!
//! Every check carries an expected verdict that is known independently of the
//! implementation (a classical theorem, an explicit factorization, or a
//! violated Newton inequality), so the suite doubles as a correctness
//! regression test: a wrong verdict makes the command exit with status 1.
//! Timings are reported next to the verdicts and can be compared against a
//! saved JSON run with `--baseline`.

use crate::bench::BenchOutputFormat;
use crate::input::strip_trailing_zeros_bigint;
use num_bigint::BigInt;
use num_rational::Ratio;
use polytool::numeric::binomial_bigint;
use polytool::root_count::isolate_real_roots_bigint_coeffs;
use polytool::sequences::{
    eulerian_polynomials_bigint, narayana_polynomials_bigint, type_b_eulerian_polynomials_bigint,
};
use polytool::vec_poly::mul;
use polytool::*;
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

/// Schema tag of the JSON output; checked when a baseline is read.
const SCHEMA: &str = "polytool.bench.real-rootedness.v1";

/// Baseline medians below this many microseconds are too noisy for the
/// `--max-slowdown` gate.
const SLOWDOWN_FLOOR_US: f64 = 1000.0;

pub(crate) fn print_real_rootedness_help() {
    println!("Real-rootedness regression suite options:");
    println!("  polytool bench real-rootedness [options]");
    println!("  --quick                   Small degrees (about 8-14); fast enough for CI");
    println!("  --degrees <d1,d2,...>     Override the degree list of the scaling families");
    println!("                            (4..=200; the full default takes a few minutes)");
    println!("                            (default: 20,40,60,80; quick: 8,14)");
    println!("  --only <substring>        Run only cases whose name contains substring");
    println!("  --repeat <n>              Time each check n times (default: 1)");
    println!("  --format <tsv|json>       Output format (default: tsv)");
    println!("  --json                    Alias for --format json");
    println!("  --baseline <run.json>     Compare verdicts and timings with a saved JSON run");
    println!("  --max-slowdown <x>        With --baseline, fail when a check whose baseline");
    println!("                            median is at least 1 ms became more than x times slower");
    println!();
    println!("Each check has an exact expected verdict; any mismatch (or a verdict that");
    println!("differs from the baseline) is reported as a failure and exits with status 1.");
    println!("Check kinds: real_rooted, interlacing (strict, f << g), root_count (distinct");
    println!("real roots), hstar_identity (Ehrhart <-> h* conversion round trips).");
}

#[derive(Debug)]
struct Options {
    quick: bool,
    degrees: Option<Vec<usize>>,
    only: Option<String>,
    repeat: Option<usize>,
    format: BenchOutputFormat,
    baseline: Option<PathBuf>,
    max_slowdown: Option<f64>,
}

fn parse_options(args: &[String]) -> Result<Options, String> {
    let mut options = Options {
        quick: false,
        degrees: None,
        only: None,
        repeat: None,
        format: BenchOutputFormat::Tsv,
        baseline: None,
        max_slowdown: None,
    };
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        let mut value = |name: &str| -> Result<&String, String> {
            i += 1;
            args.get(i).ok_or_else(|| format!("{name} expects a value"))
        };
        match flag {
            "--quick" => options.quick = true,
            "--degrees" => {
                let text = value("--degrees")?;
                let mut degrees = Vec::new();
                for part in text.split(',') {
                    let degree: usize = part.trim().parse().map_err(|_| {
                        format!("--degrees expects comma-separated integers, got '{text}'")
                    })?;
                    if !(4..=200).contains(&degree) {
                        return Err(format!(
                            "--degrees entries must be in 4..=200, got {degree}"
                        ));
                    }
                    degrees.push(degree);
                }
                if degrees.is_empty() {
                    return Err("--degrees expects at least one degree".to_string());
                }
                options.degrees = Some(degrees);
            }
            "--only" => options.only = Some(value("--only")?.clone()),
            "--repeat" => {
                let text = value("--repeat")?;
                options.repeat = Some(
                    text.parse::<usize>()
                        .map_err(|_| format!("--repeat expects a positive integer, got '{text}'"))?
                        .max(1),
                );
            }
            "--format" => options.format = BenchOutputFormat::parse(value("--format")?)?,
            "--json" => options.format = BenchOutputFormat::Json,
            "--baseline" => options.baseline = Some(PathBuf::from(value("--baseline")?)),
            "--max-slowdown" => {
                let text = value("--max-slowdown")?;
                let factor: f64 = text
                    .parse()
                    .map_err(|_| format!("--max-slowdown expects a number, got '{text}'"))?;
                if !(factor.is_finite() && factor > 1.0) {
                    return Err("--max-slowdown expects a factor greater than 1".to_string());
                }
                options.max_slowdown = Some(factor);
            }
            other => return Err(format!("unknown real-rootedness benchmark option: {other}")),
        }
        i += 1;
    }
    if options.max_slowdown.is_some() && options.baseline.is_none() {
        return Err("--max-slowdown requires --baseline".to_string());
    }
    Ok(options)
}

/// One benchmark check: a named computation with a known expected verdict.
struct Check {
    case: String,
    family: &'static str,
    kind: &'static str,
    degree: usize,
    coeff_bits: u64,
    expected: String,
    run: Box<dyn Fn() -> String>,
}

fn bits(polys: &[&[BigInt]]) -> u64 {
    polys
        .iter()
        .flat_map(|p| p.iter())
        .map(|c| c.bits())
        .max()
        .unwrap_or(0)
}

fn verdict_bool(value: bool) -> String {
    value.to_string()
}

fn verdict_interlacing(value: Option<bool>) -> String {
    match value {
        Some(flag) => flag.to_string(),
        None => "none".to_string(),
    }
}

fn real_rooted_check(case: &str, family: &'static str, poly: Vec<BigInt>, expected: bool) -> Check {
    Check {
        case: case.to_string(),
        family,
        kind: "real_rooted",
        degree: poly.len() - 1,
        coeff_bits: bits(&[&poly]),
        expected: verdict_bool(expected),
        run: Box::new(move || verdict_bool(is_real_rooted_bigint_coeffs(&poly))),
    }
}

fn root_count_check(case: &str, family: &'static str, poly: Vec<BigInt>, expected: usize) -> Check {
    Check {
        case: case.to_string(),
        family,
        kind: "root_count",
        degree: poly.len() - 1,
        coeff_bits: bits(&[&poly]),
        expected: expected.to_string(),
        run: Box::new(move || isolate_real_roots_bigint_coeffs(&poly).len().to_string()),
    }
}

fn interlacing_check(
    case: &str,
    family: &'static str,
    f: Vec<BigInt>,
    g: Vec<BigInt>,
    expected: Option<bool>,
) -> Check {
    Check {
        case: case.to_string(),
        family,
        kind: "interlacing",
        degree: f.len().max(g.len()) - 1,
        coeff_bits: bits(&[&f, &g]),
        expected: verdict_interlacing(expected),
        run: Box::new(move || verdict_interlacing(check_interlacing_bigint_coeffs(&f, &g))),
    }
}

fn identity_check(
    case: &str,
    family: &'static str,
    degree: usize,
    coeff_bits: u64,
    run: impl Fn() -> bool + 'static,
) -> Check {
    Check {
        case: case.to_string(),
        family,
        kind: "hstar_identity",
        degree,
        coeff_bits,
        expected: "true".to_string(),
        run: Box::new(move || verdict_bool(run())),
    }
}

fn ints(values: &[i64]) -> Vec<BigInt> {
    values.iter().copied().map(BigInt::from).collect()
}

/// `p(t)^k`-free helper: the product of two coefficient vectors, trimmed of
/// trailing zeros.
fn product(a: &[BigInt], b: &[BigInt]) -> Vec<BigInt> {
    strip_trailing_zeros_bigint(&mul(a, b)).to_vec()
}

/// `(1 + t)^m` as an ascending coefficient vector.
fn one_plus_t_power(m: usize) -> Vec<BigInt> {
    (0..=m).map(|k| binomial_bigint(m, k)).collect()
}

/// The member of a generated sequence with exactly the given degree.
fn pick(table: &[Vec<BigInt>], degree: usize) -> Vec<BigInt> {
    table
        .iter()
        .find(|poly| poly.len() - 1 == degree)
        .unwrap_or_else(|| panic!("sequence has no member of degree {degree}"))
        .clone()
}

/// Ehrhart coefficients of the unit cube `[0,1]^dim`, namely `(t+1)^dim`.
fn cube_ehrhart(dim: usize) -> Vec<Ratio<BigInt>> {
    one_plus_t_power(dim)
        .into_iter()
        .map(Ratio::from_integer)
        .collect()
}

/// `true` when converting `hstar` to Ehrhart coefficients and back returns it.
fn hstar_round_trips(hstar: &[BigInt]) -> bool {
    let ehrhart = hstar_to_ehrhart_bigint_coeffs(hstar);
    match ehrhart_to_hstar_bigint(&ehrhart) {
        Ok(back) => strip_trailing_zeros_bigint(&back) == hstar,
        Err(_) => false,
    }
}

fn build_checks(degrees: &[usize]) -> Vec<Check> {
    let mut checks = Vec::new();
    let max_degree = degrees.iter().copied().max().unwrap_or(4) + 3;
    let eulerian = eulerian_polynomials_bigint(max_degree);
    let narayana = narayana_polynomials_bigint(max_degree);
    let type_b = type_b_eulerian_polynomials_bigint(max_degree);

    // Small counterexample from the unit tests: positive coefficients but not
    // real-rooted.  It is the only fixed-degree case; the rest scale with
    // the requested degrees.
    checks.push(real_rooted_check(
        "positive_quintic_counterexample",
        "counterexample",
        ints(&[1, 43, 196, 168, 23, 1]),
        false,
    ));

    for &d in degrees {
        // Eulerian polynomials A_{d+1} (degree d): real-rooted with d simple
        // negative roots (Frobenius); consecutive ones interlace strictly.
        let a = pick(&eulerian, d);
        let a_next = pick(&eulerian, d + 1);
        checks.push(real_rooted_check(
            &format!("eulerian_{d}"),
            "eulerian",
            a.clone(),
            true,
        ));
        checks.push(root_count_check(
            &format!("eulerian_{d}"),
            "eulerian",
            a.clone(),
            d,
        ));
        checks.push(interlacing_check(
            &format!("eulerian_{d}_next"),
            "eulerian",
            a.clone(),
            a_next.clone(),
            Some(true),
        ));
        checks.push(interlacing_check(
            &format!("eulerian_{d}_degree_gap"),
            "eulerian",
            a.clone(),
            pick(&eulerian, d + 2),
            None,
        ));

        // Narayana polynomials N_{d+1} (degree d): real-rooted with simple
        // roots (they are Jacobi-type orthogonal polynomials).
        let n = pick(&narayana, d);
        checks.push(real_rooted_check(
            &format!("narayana_{d}"),
            "narayana",
            n.clone(),
            true,
        ));
        checks.push(root_count_check(
            &format!("narayana_{d}"),
            "narayana",
            n.clone(),
            d,
        ));
        checks.push(interlacing_check(
            &format!("narayana_{d}_next"),
            "narayana",
            n,
            pick(&narayana, d + 1),
            Some(true),
        ));

        // Type B Eulerian polynomials: real-rooted (Brenti).
        checks.push(real_rooted_check(
            &format!("type_b_eulerian_{d}"),
            "type_b_eulerian",
            pick(&type_b, d),
            true,
        ));

        // h*-vectors.
        // Unit cube of dimension d+1: h* is the Eulerian polynomial A_{d+1}.
        let cube_hstar = a.clone();
        checks.push(identity_check(
            &format!("hstar_cube_dim{}", d + 1),
            "hstar",
            d,
            bits(&[&cube_hstar]),
            move || match ehrhart_to_hstar_bigint(&cube_ehrhart(d + 1)) {
                Ok(h) => strip_trailing_zeros_bigint(&h) == cube_hstar.as_slice(),
                Err(_) => false,
            },
        ));
        // Product of simplices Delta_d x Delta_{d+7}: h*_i = C(d,i) C(d+7,i),
        // a Hadamard product of (1+t)^d and (1+t)^{d+7}, hence real-rooted
        // (Malo--Schur--Szego); it round-trips through the Ehrhart polynomial.
        let simplices: Vec<BigInt> = (0..=d)
            .map(|i| binomial_bigint(d, i) * binomial_bigint(d + 7, i))
            .collect();
        checks.push(real_rooted_check(
            &format!("hstar_simplex_product_{d}x{}", d + 7),
            "hstar",
            simplices.clone(),
            true,
        ));
        let simplices_bits = bits(&[&simplices]);
        checks.push(identity_check(
            &format!("hstar_simplex_product_{d}x{}", d + 7),
            "hstar",
            d,
            simplices_bits,
            move || hstar_round_trips(&simplices),
        ));
        // Cross-polytope of dimension d: h* = (1+t)^d, a single root of
        // multiplicity d.
        let cross = one_plus_t_power(d);
        checks.push(real_rooted_check(
            &format!("hstar_cross_polytope_{d}"),
            "hstar",
            cross.clone(),
            true,
        ));
        checks.push(root_count_check(
            &format!("hstar_cross_polytope_{d}"),
            "hstar",
            cross.clone(),
            1,
        ));
        let cross_bits = bits(&[&cross]);
        checks.push(identity_check(
            &format!("hstar_cross_polytope_{d}"),
            "hstar",
            d,
            cross_bits,
            move || hstar_round_trips(&cross),
        ));

        // Positive-coefficient products with a known verdict.
        let trinomial = ints(&[1, 1, 1]); // roots: primitive cube roots of unity
        let bad_quadratic = ints(&[1, 2, 2]); // discriminant -4
        let good_quadratic = ints(&[1, 3, 1]); // two real roots
        let a_prev = pick(&eulerian, d - 2);

        // A_{d-1} * (1 + t + t^2) has degree d, positive coefficients and two
        // non-real roots; the Eulerian factor of degree d-2 contributes d-2
        // simple real roots.
        let times_trinomial = product(&a_prev, &trinomial);
        let label = format!("eulerian_{}_times_trinomial", d - 2);
        checks.push(real_rooted_check(
            &label,
            "counterexample",
            times_trinomial.clone(),
            false,
        ));
        checks.push(root_count_check(
            &label,
            "counterexample",
            times_trinomial.clone(),
            d - 2,
        ));
        // The degree d-1 Eulerian polynomial against that non-real-rooted
        // polynomial of degree d: strict interlacing must fail.
        checks.push(interlacing_check(
            &format!("eulerian_{}_vs_eulerian_{}_times_trinomial", d - 1, d - 2),
            "counterexample",
            pick(&eulerian, d - 1),
            times_trinomial,
            Some(false),
        ));
        // (1+t)^{d-2} (1 + 2t + 2t^2): one real root of high multiplicity and
        // a non-real pair.
        let multiple_root = product(&one_plus_t_power(d - 2), &bad_quadratic);
        checks.push(real_rooted_check(
            &format!("one_plus_t_{}_times_bad_quadratic", d - 2),
            "counterexample",
            multiple_root.clone(),
            false,
        ));
        checks.push(root_count_check(
            &format!("one_plus_t_{}_times_bad_quadratic", d - 2),
            "counterexample",
            multiple_root,
            1,
        ));
        // 1 + t^{d+1} with d even has d+1 odd: only the root -1 is real.
        let mut binomial_sum = vec![BigInt::from(0); d + 2];
        binomial_sum[0] = BigInt::from(1);
        binomial_sum[d + 1] = BigInt::from(1);
        checks.push(real_rooted_check(
            &format!("one_plus_t_pow_{}", d + 1),
            "counterexample",
            binomial_sum.clone(),
            false,
        ));
        checks.push(root_count_check(
            &format!("one_plus_t_pow_{}", d + 1),
            "counterexample",
            binomial_sum,
            (d + 1) % 2,
        ));
        // Eulerian polynomial with a middle coefficient pushed below the
        // geometric mean of its neighbours: positive coefficients, but the
        // log-concavity (hence Newton) inequality fails, so it is not
        // real-rooted.
        let mut dented = a.clone();
        let mid = d / 2;
        let geometric = (&dented[mid - 1] * &dented[mid + 1]).sqrt();
        dented[mid] = geometric - 1;
        checks.push(real_rooted_check(
            &format!("eulerian_{d}_dented"),
            "counterexample",
            dented,
            false,
        ));
        // A_d * (1 + 3t + t^2) is a product of real-rooted polynomials.
        checks.push(real_rooted_check(
            &format!("eulerian_{d}_times_good_quadratic"),
            "eulerian",
            product(&a, &good_quadratic),
            true,
        ));
    }
    checks
}

struct Record {
    case: String,
    family: &'static str,
    kind: &'static str,
    degree: usize,
    coeff_bits: u64,
    expected: String,
    verdict: String,
    repeat: usize,
    min_us: f64,
    median_us: f64,
    failure: Option<String>,
}

impl Record {
    fn key(&self) -> String {
        format!("{}|{}", self.case, self.kind)
    }
}

fn run_check(check: &Check, repeat: usize) -> Record {
    let mut times = Vec::with_capacity(repeat);
    let mut verdict = String::new();
    let mut failure = None;
    for run in 0..repeat {
        let start = Instant::now();
        let current = black_box((check.run)());
        times.push(start.elapsed().as_secs_f64() * 1_000_000.0);
        if run == 0 {
            verdict = current;
        } else if current != verdict {
            failure = Some(format!(
                "nondeterministic verdict: {verdict} then {current}"
            ));
        }
    }
    if failure.is_none() && verdict != check.expected {
        failure = Some(format!(
            "expected {} but computed {}",
            check.expected, verdict
        ));
    }
    times.sort_by(|a, b| a.partial_cmp(b).expect("timings are finite"));
    let median = if times.len() % 2 == 1 {
        times[times.len() / 2]
    } else {
        (times[times.len() / 2 - 1] + times[times.len() / 2]) / 2.0
    };
    Record {
        case: check.case.clone(),
        family: check.family,
        kind: check.kind,
        degree: check.degree,
        coeff_bits: check.coeff_bits,
        expected: check.expected.clone(),
        verdict,
        repeat,
        min_us: times[0],
        median_us: median,
        failure,
    }
}

/// Verdict and median of one check in a saved run.
struct BaselinePoint {
    verdict: String,
    median_us: f64,
}

fn read_baseline(path: &PathBuf) -> Result<BTreeMap<String, BaselinePoint>, String> {
    let text = fs::read_to_string(path)
        .map_err(|e| format!("failed to read baseline `{}`: {e}", path.display()))?;
    let value: Value = serde_json::from_str(&text)
        .map_err(|e| format!("failed to parse baseline `{}`: {e}", path.display()))?;
    if value.get("schema").and_then(Value::as_str) != Some(SCHEMA) {
        return Err(format!(
            "baseline `{}` is not a `{SCHEMA}` run",
            path.display()
        ));
    }
    let records = value
        .get("records")
        .and_then(Value::as_array)
        .ok_or_else(|| format!("baseline `{}` has no records", path.display()))?;
    let mut map = BTreeMap::new();
    for record in records {
        let field = |name: &str| record.get(name).and_then(Value::as_str);
        let (Some(case), Some(kind), Some(verdict)) =
            (field("case"), field("kind"), field("verdict"))
        else {
            return Err(format!("malformed record in baseline `{}`", path.display()));
        };
        let median_us = record
            .get("median_us")
            .and_then(Value::as_f64)
            .ok_or_else(|| format!("record without median_us in `{}`", path.display()))?;
        map.insert(
            format!("{case}|{kind}"),
            BaselinePoint {
                verdict: verdict.to_string(),
                median_us,
            },
        );
    }
    Ok(map)
}

pub(crate) fn cmd_bench_real_rootedness(args: &[String]) {
    let options = match parse_options(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            print_real_rootedness_help();
            std::process::exit(1);
        }
    };
    let baseline = match options.baseline.as_ref().map(read_baseline).transpose() {
        Ok(baseline) => baseline,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(1);
        }
    };
    let degrees = options.degrees.clone().unwrap_or_else(|| {
        if options.quick {
            vec![8, 14]
        } else {
            vec![20, 40, 60, 80]
        }
    });
    let repeat = options.repeat.unwrap_or(1);
    let mut records: Vec<Record> = build_checks(&degrees)
        .iter()
        .filter(|check| {
            options
                .only
                .as_ref()
                .is_none_or(|only| check.case.contains(only.as_str()))
        })
        .map(|check| run_check(check, repeat))
        .collect();

    // Baseline comparison: verdict changes are failures, slowdowns only when
    // a factor was requested.
    let mut ratios: Vec<Option<f64>> = vec![None; records.len()];
    if let Some(baseline) = &baseline {
        for (record, ratio) in records.iter_mut().zip(ratios.iter_mut()) {
            let Some(point) = baseline.get(&record.key()) else {
                continue;
            };
            if point.median_us > 0.0 {
                *ratio = Some(record.median_us / point.median_us);
            }
            if record.failure.is_none() && point.verdict != record.verdict {
                record.failure = Some(format!(
                    "verdict {} differs from baseline verdict {}",
                    record.verdict, point.verdict
                ));
            }
            if let (Some(limit), Some(ratio)) = (options.max_slowdown, *ratio) {
                if record.failure.is_none() && point.median_us >= SLOWDOWN_FLOOR_US && ratio > limit
                {
                    record.failure = Some(format!(
                        "{ratio:.2}x slower than baseline (limit {limit:.2}x)"
                    ));
                }
            }
        }
    }
    let failures = records.iter().filter(|r| r.failure.is_some()).count();

    match options.format {
        BenchOutputFormat::Tsv => print_tsv(&records, &ratios, baseline.is_some()),
        BenchOutputFormat::Json => print_json(&records, &ratios, &options, repeat, &degrees),
    }
    eprintln!(
        "real-rootedness suite: {} checks, {} failures",
        records.len(),
        failures
    );
    for record in records.iter().filter(|r| r.failure.is_some()) {
        eprintln!(
            "FAIL {} [{}]: {}",
            record.case,
            record.kind,
            record.failure.as_deref().unwrap_or("")
        );
    }
    if failures > 0 {
        std::process::exit(1);
    }
}

fn status(record: &Record) -> &'static str {
    if record.failure.is_some() {
        "FAIL"
    } else {
        "ok"
    }
}

fn print_tsv(records: &[Record], ratios: &[Option<f64>], with_baseline: bool) {
    let mut header = String::from(
        "case\tfamily\tkind\tdegree\tcoeff_bits\texpected\tverdict\tstatus\trepeat\tmin_us\tmedian_us",
    );
    if with_baseline {
        header.push_str("\tratio_vs_baseline");
    }
    println!("{header}");
    for (record, ratio) in records.iter().zip(ratios) {
        let mut line = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}",
            record.case,
            record.family,
            record.kind,
            record.degree,
            record.coeff_bits,
            record.expected,
            record.verdict,
            status(record),
            record.repeat,
            record.min_us,
            record.median_us,
        );
        if with_baseline {
            match ratio {
                Some(ratio) => line.push_str(&format!("\t{ratio:.2}")),
                None => line.push_str("\t-"),
            }
        }
        println!("{line}");
    }
}

fn print_json(
    records: &[Record],
    ratios: &[Option<f64>],
    options: &Options,
    repeat: usize,
    degrees: &[usize],
) {
    let rows: Vec<Value> = records
        .iter()
        .zip(ratios)
        .map(|(record, ratio)| {
            json!({
                "case": record.case,
                "family": record.family,
                "kind": record.kind,
                "degree": record.degree,
                "coeff_bits": record.coeff_bits,
                "expected": record.expected,
                "verdict": record.verdict,
                "status": status(record),
                "failure": record.failure,
                "repeat": record.repeat,
                "min_us": record.min_us,
                "median_us": record.median_us,
                "ratio_vs_baseline": ratio,
            })
        })
        .collect();
    let report = json!({
        "schema": SCHEMA,
        "quick": options.quick,
        "degrees": degrees,
        "repeat": repeat,
        "checks": records.len(),
        "failures": records.iter().filter(|r| r.failure.is_some()).count(),
        "records": rows,
    });
    println!(
        "{}",
        serde_json::to_string_pretty(&report).expect("benchmark report serializes")
    );
}
