//! `polytool family-check`: properties, consecutive interlacing and optional
//! recurrence search in one pass.

use crate::help::print_coefficient_input_help;
use crate::input::{read_poly_parse_results_bigint, strip_trailing_zeros_bigint, OutputFormat};
use crate::interlacing::{interlacing_report, interlacing_report_json, InterlacingReport};
use crate::json_out::{json_bool_option, json_string};
use crate::properties::{property_labels, property_report, property_report_json, PropertyReport};
use num_bigint::BigInt;
use polytool::recurrence::*;

pub(crate) fn print_family_check_help() {
    println!("Usage:");
    println!("  polytool family-check [options]");
    println!("  polytool family-check --help");
    println!();
    println!("Check a polynomial family in one pass: properties, consecutive");
    println!("interlacing, and optionally recurrence search.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --recurrence                  Search for an adaptive recurrence");
    println!("  --require-real-rooted         Fail if any row is not real-rooted");
    println!("  --require-palindromic         Fail if any row is not palindromic");
    println!("  --require-gamma-positive      Fail if any row is not gamma-positive");
    println!("  --require-unimodal            Fail if any row is not unimodal");
    println!("  --require-log-concave         Fail if any row is not log-concave");
    println!("  --require-ultra-log-concave   Fail if any row is not ultra-log-concave");
    println!("  --require-weak-interlacing    Fail if consecutive rows do not weakly interlace");
    println!("  --json                        Emit machine-readable JSON");
    println!("  -h, --help                    Print this help text");
    println!();
    println!("All checks, including --recurrence, accept arbitrary-size integer coefficients.");
}

#[derive(Clone, Debug)]
struct FamilyCheckOptions {
    format: OutputFormat,
    recurrence: bool,
    require_real_rooted: bool,
    require_palindromic: bool,
    require_gamma_positive: bool,
    require_unimodal: bool,
    require_log_concave: bool,
    require_ultra_log_concave: bool,
    require_weak_interlacing: bool,
}

#[derive(Clone, Debug)]
struct ParseErrorReport {
    index: usize,
    error: String,
}

impl Default for FamilyCheckOptions {
    fn default() -> Self {
        Self {
            format: OutputFormat::Text,
            recurrence: false,
            require_real_rooted: false,
            require_palindromic: false,
            require_gamma_positive: false,
            require_unimodal: false,
            require_log_concave: false,
            require_ultra_log_concave: false,
            require_weak_interlacing: false,
        }
    }
}

#[derive(Clone, Debug)]
struct RecurrenceSummary {
    found: bool,
    recurrence: Option<String>,
    unknowns: Option<usize>,
    weighted_unknowns: Option<usize>,
    equations: Option<usize>,
    fit_polynomials: Option<usize>,
    verification_polynomials: Option<usize>,
    candidates_tried: Option<usize>,
    error: Option<String>,
}

fn parse_family_check_args(args: &[String]) -> Result<FamilyCheckOptions, String> {
    let mut options = FamilyCheckOptions::default();
    for arg in args {
        match arg.as_str() {
            "--json" => options.format = OutputFormat::Json,
            "--recurrence" => options.recurrence = true,
            "--require-real-rooted" => options.require_real_rooted = true,
            "--require-palindromic" => options.require_palindromic = true,
            "--require-gamma-positive" => options.require_gamma_positive = true,
            "--require-unimodal" => options.require_unimodal = true,
            "--require-log-concave" => options.require_log_concave = true,
            "--require-ultra-log-concave" => options.require_ultra_log_concave = true,
            "--require-weak-interlacing" => options.require_weak_interlacing = true,
            other => return Err(format!("unknown family-check option: {other}")),
        }
    }
    Ok(options)
}

fn first_family_failure(
    reports: &[PropertyReport],
    pairs: &[InterlacingReport],
    recurrence: Option<&RecurrenceSummary>,
    options: &FamilyCheckOptions,
) -> Option<String> {
    for report in reports {
        if options.require_real_rooted && !report.real_rooted {
            return Some(format!("polynomial {} is not real-rooted", report.index));
        }
        if options.require_palindromic && !report.palindromic {
            return Some(format!("polynomial {} is not palindromic", report.index));
        }
        if options.require_gamma_positive && !report.gamma_positive {
            return Some(format!("polynomial {} is not gamma-positive", report.index));
        }
        if options.require_unimodal && !report.unimodal {
            return Some(format!("polynomial {} is not unimodal", report.index));
        }
        if options.require_log_concave && !report.log_concave {
            return Some(format!("polynomial {} is not log-concave", report.index));
        }
        if options.require_ultra_log_concave && !report.ultra_log_concave {
            return Some(format!(
                "polynomial {} is not ultra-log-concave",
                report.index
            ));
        }
    }
    if options.require_weak_interlacing {
        for pair in pairs {
            if pair.weak != Some(true) {
                return Some(format!(
                    "polynomials {} and {} do not weakly interlace",
                    pair.left_index, pair.right_index
                ));
            }
        }
    }
    if options.recurrence {
        if let Some(summary) = recurrence {
            if !summary.found {
                return Some(
                    summary
                        .error
                        .clone()
                        .unwrap_or_else(|| "no recurrence found".to_string()),
                );
            }
        }
    }
    None
}

fn find_family_recurrence(polys: &[Vec<BigInt>]) -> RecurrenceSummary {
    if polys.len() < 3 {
        return RecurrenceSummary {
            found: false,
            recurrence: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("need at least 3 polynomials".to_string()),
        };
    }
    let search = AdaptiveSearchOptions::default();
    let rational_polys: Vec<Vec<BigRational>> = polys
        .iter()
        .map(|row| row.iter().cloned().map(BigRational::from_integer).collect())
        .collect();
    match find_recurrence_adaptive_rational(&rational_polys, &search) {
        Some(result) => RecurrenceSummary {
            found: true,
            recurrence: Some(result.recurrence.to_string()),
            unknowns: Some(result.num_unknowns),
            weighted_unknowns: Some(result.weighted_unknowns),
            equations: Some(result.num_equations),
            fit_polynomials: Some(result.fit_polynomials),
            verification_polynomials: Some(result.verification_polynomials),
            candidates_tried: Some(result.candidates_tried),
            error: None,
        },
        None => RecurrenceSummary {
            found: false,
            recurrence: None,
            unknowns: None,
            weighted_unknowns: None,
            equations: None,
            fit_polynomials: None,
            verification_polynomials: None,
            candidates_tried: None,
            error: Some("no recurrence found within the search bounds".to_string()),
        },
    }
}

fn family_recurrence_unavailable(error: impl Into<String>) -> RecurrenceSummary {
    RecurrenceSummary {
        found: false,
        recurrence: None,
        unknowns: None,
        weighted_unknowns: None,
        equations: None,
        fit_polynomials: None,
        verification_polynomials: None,
        candidates_tried: None,
        error: Some(error.into()),
    }
}

fn parse_error_report_json(report: &ParseErrorReport) -> String {
    format!(
        "{{\"index\":{},\"error\":{}}}",
        report.index,
        json_string(&report.error)
    )
}

fn recurrence_summary_json(summary: &RecurrenceSummary) -> String {
    format!(
        "{{\"found\":{},\"recurrence\":{},\"unknowns\":{},\"weighted_unknowns\":{},\
         \"equations\":{},\"fit_polynomials\":{},\"verification_polynomials\":{},\
         \"candidates_tried\":{},\"error\":{}}}",
        summary.found,
        summary
            .recurrence
            .as_ref()
            .map(|r| json_string(r))
            .unwrap_or_else(|| "null".to_string()),
        summary
            .unknowns
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .weighted_unknowns
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .equations
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .fit_polynomials
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .verification_polynomials
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .candidates_tried
            .map(|n| n.to_string())
            .unwrap_or_else(|| "null".to_string()),
        summary
            .error
            .as_ref()
            .map(|e| json_string(e))
            .unwrap_or_else(|| "null".to_string())
    )
}

pub(crate) fn cmd_family_check(args: &[String]) {
    let options = match parse_family_check_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut parse_errors = Vec::new();
    let mut polys = Vec::new();
    for (index, item) in read_poly_parse_results_bigint().into_iter().enumerate() {
        match item {
            Ok(coeffs) => {
                polys.push((index, strip_trailing_zeros_bigint(&coeffs).to_vec()));
            }
            Err(error) => parse_errors.push(ParseErrorReport { index, error }),
        }
    }
    let reports = polys
        .iter()
        .map(|(index, coeffs)| property_report(*index, coeffs))
        .collect::<Vec<_>>();
    let pairs = if parse_errors.is_empty() {
        (0..polys.len().saturating_sub(1))
            .map(|pair_index| {
                let left_index = polys[pair_index].0;
                let right_index = polys[pair_index + 1].0;
                interlacing_report(
                    pair_index,
                    left_index,
                    right_index,
                    &polys[pair_index].1,
                    &polys[pair_index + 1].1,
                )
            })
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    let recurrence = options.recurrence.then(|| {
        if !parse_errors.is_empty() {
            return family_recurrence_unavailable(
                "recurrence skipped because one or more polynomials failed to parse",
            );
        }
        let rows = polys
            .iter()
            .map(|(_, coeffs)| coeffs.clone())
            .collect::<Vec<_>>();
        find_family_recurrence(&rows)
    });
    let first_failure = if let Some(error) = parse_errors.first() {
        Some(format!(
            "polynomial {} failed to parse: {}",
            error.index, error.error
        ))
    } else {
        first_family_failure(&reports, &pairs, recurrence.as_ref(), &options)
    };
    let passed = first_failure.is_none();

    if options.format == OutputFormat::Json {
        let parse_error_items = parse_errors
            .iter()
            .map(parse_error_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let items = reports
            .iter()
            .map(property_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let pair_items = pairs
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let recurrence_json = recurrence
            .as_ref()
            .map(recurrence_summary_json)
            .unwrap_or_else(|| "null".to_string());
        let failure_json = first_failure
            .as_ref()
            .map(|failure| json_string(failure))
            .unwrap_or_else(|| "null".to_string());
        println!(
            "{{\"item_count\":{},\"all_required_checks_passed\":{},\
             \"first_failure\":{},\"parse_errors\":[{}],\"items\":[{}],\
             \"consecutive_pairs\":[{}],\"recurrence\":{}}}",
            reports.len(),
            passed,
            failure_json,
            parse_error_items,
            items,
            pair_items,
            recurrence_json
        );
        if !passed {
            std::process::exit(1);
        }
        return;
    }

    println!("Polynomial family check");
    println!("- Polynomials: {}", reports.len());
    println!(
        "- Required checks: {}",
        if passed { "passed" } else { "failed" }
    );
    if let Some(failure) = &first_failure {
        println!("- First failure: {failure}");
    }
    if !parse_errors.is_empty() {
        println!();
        println!("Parse errors:");
        for error in &parse_errors {
            println!("{}: {}", error.index, error.error);
        }
    }
    println!();
    println!("Properties:");
    for report in &reports {
        println!("{}: {}", report.index, property_labels(report).join(", "));
    }
    if !pairs.is_empty() {
        println!();
        println!("Consecutive interlacing:");
        for pair in &pairs {
            println!(
                "{}-{}: weak={}, strict={}, status={}",
                pair.left_index,
                pair.right_index,
                json_bool_option(pair.weak),
                json_bool_option(pair.strict),
                pair.status
            );
        }
    }
    if let Some(summary) = &recurrence {
        println!();
        println!("Recurrence:");
        if let Some(rec) = &summary.recurrence {
            println!("- Found: {rec}");
        } else if let Some(error) = &summary.error {
            println!("- Not found: {error}");
        } else {
            println!("- Not found");
        }
    }
    if !passed {
        std::process::exit(1);
    }
}
