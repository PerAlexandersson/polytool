//! `polytool interlacing` and `polytool interlacing-profile`.

use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint, OutputFormat};
use crate::json_out::{json_bool_option, json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;

#[derive(Clone, Debug)]
pub(crate) struct InterlacingReport {
    pair_index: usize,
    pub(crate) left_index: usize,
    pub(crate) right_index: usize,
    p: String,
    q: String,
    pub(crate) strict: Option<bool>,
    pub(crate) weak: Option<bool>,
    pub(crate) status: String,
}

fn interlacing_status(strict: Option<bool>, weak: Option<bool>) -> String {
    match (strict, weak) {
        (Some(true), _) => "strictly_interlace",
        (_, Some(true)) => "weakly_interlace",
        (Some(false), Some(false)) => "do_not_interlace",
        (Some(false), None) => "not_real_rooted_or_incompatible",
        (None, Some(false)) => "not_real_rooted_or_incompatible",
        (None, None) => "not_real_rooted_or_incompatible",
    }
    .to_string()
}

pub(crate) fn interlacing_report(
    pair_index: usize,
    left_index: usize,
    right_index: usize,
    p: &[BigInt],
    q: &[BigInt],
) -> InterlacingReport {
    let strict = check_interlacing_bigint_coeffs(p, q);
    let weak = check_weak_interlacing_bigint_coeffs(p, q);
    InterlacingReport {
        pair_index,
        left_index,
        right_index,
        p: format_poly_bigint_coeffs(p),
        q: format_poly_bigint_coeffs(q),
        strict,
        weak,
        status: interlacing_status(strict, weak),
    }
}

pub(crate) fn interlacing_report_json(report: &InterlacingReport) -> String {
    format!(
        "{{\"pair_index\":{},\"left_index\":{},\"right_index\":{},\
         \"p\":{},\"q\":{},\"strict\":{},\"weak\":{},\"status\":{}}}",
        report.pair_index,
        report.left_index,
        report.right_index,
        json_string(&report.p),
        json_string(&report.q),
        json_bool_option(report.strict),
        json_bool_option(report.weak),
        json_string(&report.status)
    )
}

fn print_interlacing_reports_json(reports: &[InterlacingReport]) {
    println!(
        "{{\"pairs\":[{}]}}",
        reports
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

#[derive(Clone, Debug)]
struct InterlacingProfileReport {
    index: usize,
    polynomial: String,
    previous_count: usize,
    checked_previous_count: usize,
    interlacing_previous_count: usize,
    strict_previous_count: usize,
    weak_previous_count: usize,
    previous_interlacing_indices: Vec<usize>,
    previous: Vec<InterlacingReport>,
}

fn interlacing_report_has_interlacing(report: &InterlacingReport) -> bool {
    matches!(
        report.status.as_str(),
        "strictly_interlace" | "weakly_interlace"
    )
}

fn interlacing_profile_report(
    index: usize,
    polynomial: &[BigInt],
    previous: Vec<InterlacingReport>,
) -> InterlacingProfileReport {
    let previous_interlacing_indices = previous
        .iter()
        .filter(|report| interlacing_report_has_interlacing(report))
        .map(|report| report.left_index)
        .collect::<Vec<_>>();
    InterlacingProfileReport {
        index,
        polynomial: format_poly_bigint_coeffs(polynomial),
        previous_count: index,
        checked_previous_count: previous.len(),
        interlacing_previous_count: previous_interlacing_indices.len(),
        strict_previous_count: previous
            .iter()
            .filter(|report| report.strict == Some(true))
            .count(),
        weak_previous_count: previous
            .iter()
            .filter(|report| report.weak == Some(true))
            .count(),
        previous_interlacing_indices,
        previous,
    }
}

fn interlacing_profile_report_json(report: &InterlacingProfileReport) -> String {
    format!(
        "{{\"index\":{},\"polynomial\":{},\"previous_count\":{},\
         \"checked_previous_count\":{},\"interlacing_previous_count\":{},\
         \"strict_previous_count\":{},\"weak_previous_count\":{},\
         \"previous_interlacing_indices\":{},\"previous\":[{}]}}",
        report.index,
        json_string(&report.polynomial),
        report.previous_count,
        report.checked_previous_count,
        report.interlacing_previous_count,
        report.strict_previous_count,
        report.weak_previous_count,
        json_usize_vec(&report.previous_interlacing_indices),
        report
            .previous
            .iter()
            .map(interlacing_report_json)
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn print_interlacing_profile_reports_json(reports: &[InterlacingProfileReport]) {
    println!(
        "{{\"items\":[{}]}}",
        reports
            .iter()
            .map(interlacing_profile_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

pub(crate) fn cmd_interlacing(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let polys = read_polys_bigint();
    if polys.len() < 2 {
        eprintln!("Need at least two polynomials for interlacing check.");
        return;
    }
    let reports = polys
        .windows(2)
        .enumerate()
        .map(|(i, pair)| {
            let p = strip_trailing_zeros_bigint(&pair[0]);
            let q = strip_trailing_zeros_bigint(&pair[1]);
            interlacing_report(i, i, i + 1, p, q)
        })
        .collect::<Vec<_>>();
    if format == OutputFormat::Json {
        print_interlacing_reports_json(&reports);
        return;
    }
    for report in reports {
        let status = match report.status.as_str() {
            "strictly_interlace" => "strictly interlace",
            "weakly_interlace" => "weakly interlace (shared roots)",
            "do_not_interlace" => "do NOT interlace",
            _ => "incompatible degrees",
        };
        println!("{} & {}: {}", report.p, report.q, status);
    }
}

pub(crate) fn cmd_interlacing_profile(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let polys = read_polys_bigint()
        .into_iter()
        .map(|coeffs| strip_trailing_zeros_bigint(&coeffs).to_vec())
        .collect::<Vec<_>>();
    if polys.is_empty() {
        eprintln!("Need at least one polynomial for interlacing profile.");
        return;
    }

    let mut pair_index = 0;
    let mut reports = Vec::with_capacity(polys.len());
    for current in 0..polys.len() {
        let mut previous = Vec::new();
        for left in (0..current).rev() {
            let report =
                interlacing_report(pair_index, left, current, &polys[left], &polys[current]);
            pair_index += 1;
            let interlaces = interlacing_report_has_interlacing(&report);
            previous.push(report);
            if !interlaces {
                break;
            }
        }
        reports.push(interlacing_profile_report(
            current,
            &polys[current],
            previous,
        ));
    }

    if format == OutputFormat::Json {
        print_interlacing_profile_reports_json(&reports);
        return;
    }

    for report in reports {
        println!(
            "row {}: {} previous rows interlace before first failure (checked {}/{}); \
             strict={}, weak={}; indices={:?}; {}",
            report.index,
            report.interlacing_previous_count,
            report.checked_previous_count,
            report.previous_count,
            report.strict_previous_count,
            report.weak_previous_count,
            report.previous_interlacing_indices,
            report.polynomial
        );
    }
}
