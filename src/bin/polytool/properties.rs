//! `polytool real-rooted` and `polytool properties`.

use crate::input::{
    polynomial_degree, polynomial_degree_bigint, read_polys_bigint, strip_trailing_zeros_bigint,
    OutputFormat,
};
use crate::json_out::{json_bigint_vec, json_string};
use num_bigint::BigInt;
use polytool::*;

#[derive(Clone, Debug)]
pub(crate) struct PropertyReport {
    pub(crate) index: usize,
    coefficients: Vec<BigInt>,
    polynomial: String,
    degree: usize,
    pub(crate) real_rooted: bool,
    simple_roots: bool,
    pub(crate) palindromic: bool,
    pub(crate) gamma_positive: bool,
    gamma_coefficients: Option<Vec<BigInt>>,
    pub(crate) unimodal: bool,
    pub(crate) log_concave: bool,
    pub(crate) ultra_log_concave: bool,
}

fn property_report_i64(index: usize, coeffs: &[i64]) -> PropertyReport {
    // The `i64` routine also returns `None` when a gamma coefficient of a
    // palindromic row overflows `i64`; the `BigInt` routine then gives the
    // exact answer (and `None` again only for a non-palindromic row).
    let gamma_coefficients = gamma_coefficients_ignoring_initial_zeros(coeffs)
        .map(|gamma| gamma.into_iter().map(BigInt::from).collect::<Vec<_>>())
        .or_else(|| {
            let big: Vec<BigInt> = coeffs.iter().map(|&c| BigInt::from(c)).collect();
            gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(&big)
        });
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|g| g >= &BigInt::from(0)));
    PropertyReport {
        index,
        coefficients: coeffs.iter().map(|&c| BigInt::from(c)).collect(),
        polynomial: format_poly(coeffs),
        degree: polynomial_degree(coeffs),
        real_rooted: is_real_rooted(coeffs),
        simple_roots: has_simple_roots(coeffs),
        palindromic: is_palindromic_ignoring_initial_zeros(coeffs),
        gamma_positive,
        gamma_coefficients,
        unimodal: is_unimodal(coeffs),
        log_concave: is_log_concave(coeffs),
        ultra_log_concave: is_ultra_log_concave(coeffs),
    }
}

pub(crate) fn property_report(index: usize, coeffs: &[BigInt]) -> PropertyReport {
    if let Some(coeffs_i64) = bigint_coeffs_to_i64(coeffs) {
        return property_report_i64(index, &coeffs_i64);
    }

    let gamma_coefficients = gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(coeffs);
    let gamma_positive = gamma_coefficients
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|g| g >= &BigInt::from(0)));
    PropertyReport {
        index,
        coefficients: coeffs.to_vec(),
        polynomial: format_poly_bigint_coeffs(coeffs),
        degree: polynomial_degree_bigint(coeffs),
        real_rooted: is_real_rooted_bigint_coeffs(coeffs),
        simple_roots: has_simple_roots_bigint_coeffs(coeffs),
        palindromic: is_palindromic_ignoring_initial_zeros_bigint_coeffs(coeffs),
        gamma_positive,
        gamma_coefficients,
        unimodal: is_unimodal_bigint_coeffs(coeffs),
        log_concave: is_log_concave_bigint_coeffs(coeffs),
        ultra_log_concave: is_ultra_log_concave_bigint_coeffs(coeffs),
    }
}

pub(crate) fn property_labels(report: &PropertyReport) -> Vec<String> {
    let mut props = Vec::new();
    if report.real_rooted {
        props.push("real-rooted".to_string());
    }
    if report.palindromic {
        props.push("palindromic".to_string());
    }
    if report.gamma_positive {
        if let Some(gamma) = &report.gamma_coefficients {
            props.push(format!("gamma-positive {:?}", gamma));
        }
    }
    if report.unimodal {
        props.push("unimodal".to_string());
    }
    if report.log_concave {
        props.push("log-concave".to_string());
    }
    if report.ultra_log_concave {
        props.push("ultra-log-concave".to_string());
    }
    if props.is_empty() {
        props.push("(none)".to_string());
    }
    props
}

pub(crate) fn property_report_json(report: &PropertyReport) -> String {
    let gamma = report
        .gamma_coefficients
        .as_ref()
        .map(|coeffs| json_bigint_vec(coeffs))
        .unwrap_or_else(|| "null".to_string());
    format!(
        "{{\"index\":{},\"polynomial\":{},\"coefficients\":{},\"degree\":{},\
         \"real_rooted\":{},\"simple_roots\":{},\"palindromic\":{},\
         \"gamma_positive\":{},\"gamma_coefficients\":{},\"unimodal\":{},\
         \"log_concave\":{},\"ultra_log_concave\":{}}}",
        report.index,
        json_string(&report.polynomial),
        json_bigint_vec(&report.coefficients),
        report.degree,
        report.real_rooted,
        report.simple_roots,
        report.palindromic,
        report.gamma_positive,
        gamma,
        report.unimodal,
        report.log_concave,
        report.ultra_log_concave
    )
}

fn print_property_reports_json(reports: &[PropertyReport]) {
    println!(
        "{{\"items\":[{}]}}",
        reports
            .iter()
            .map(property_report_json)
            .collect::<Vec<_>>()
            .join(",")
    );
}

pub(crate) fn cmd_real_rooted() {
    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let rr = is_real_rooted_bigint_coeffs(c);
        println!(
            "{}: {}",
            format_poly_bigint_coeffs(c),
            if rr { "real-rooted" } else { "NOT real-rooted" }
        );
    }
}

pub(crate) fn cmd_properties(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let reports = read_polys_bigint()
        .into_iter()
        .enumerate()
        .map(|(index, coeffs)| {
            let c = strip_trailing_zeros_bigint(&coeffs);
            property_report(index, c)
        })
        .collect::<Vec<_>>();
    if format == OutputFormat::Json {
        print_property_reports_json(&reports);
        return;
    }
    for report in reports {
        println!(
            "{}: {}",
            report.polynomial,
            property_labels(&report).join(", ")
        );
    }
}
