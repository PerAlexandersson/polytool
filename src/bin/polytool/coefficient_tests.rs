//! `polytool coefficient-tests`: Newton and Kurtz coefficient criteria.

use crate::help::print_coefficient_input_help;
use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint, OutputFormat};
use crate::json_out::{json_bigint_vec, json_string};
use polytool::*;

pub(crate) fn print_coefficient_tests_help() {
    println!("Usage:");
    println!("  polytool coefficient-tests [options]");
    println!("  polytool coefficient-tests --help");
    println!();
    println!("Check Newton inequalities and Kurtz's sufficient real-rootedness condition.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --json    Emit machine-readable JSON");
}

fn coefficient_check_json(check: &CoefficientInequalityCheck) -> String {
    format!(
        "{{\"index\":{},\"lhs\":{},\"rhs\":{},\"comparison\":{},\"holds\":{}}}",
        check.index,
        json_string(&check.lhs.to_string()),
        json_string(&check.rhs.to_string()),
        json_string(check.comparison),
        check.holds
    )
}

fn coefficient_criterion_json(report: &CoefficientCriterionReport) -> String {
    let checks = report
        .checks
        .iter()
        .map(coefficient_check_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"name\":{},\"reference\":{},\"applicable\":{},\"holds\":{},\
         \"implies_real_rooted\":{},\"reason\":{},\"checks\":[{}]}}",
        json_string(report.name),
        json_string(report.reference),
        report.applicable,
        report.holds,
        report.implies_real_rooted,
        report
            .reason
            .as_ref()
            .map(|reason| json_string(reason))
            .unwrap_or_else(|| "null".to_string()),
        checks
    )
}

fn coefficient_test_report_json(index: usize, report: &CoefficientTestReport) -> String {
    format!(
        "{{\"index\":{},\"coefficients\":{},\"polynomial\":{},\"degree\":{},\
         \"newton\":{},\"kurtz\":{}}}",
        index,
        json_bigint_vec(&report.coefficients),
        json_string(&format_poly_bigint_coeffs(&report.coefficients)),
        report
            .degree
            .map(|degree| degree.to_string())
            .unwrap_or_else(|| "null".to_string()),
        coefficient_criterion_json(&report.newton),
        coefficient_criterion_json(&report.kurtz)
    )
}

fn print_coefficient_test_report_text(index: usize, report: &CoefficientTestReport) {
    println!(
        "row {index}: {}",
        format_poly_bigint_coeffs(&report.coefficients)
    );
    for criterion in [&report.newton, &report.kurtz] {
        println!(
            "  {}: {}{}",
            criterion.name,
            if criterion.holds { "passed" } else { "failed" },
            if criterion.implies_real_rooted {
                " (implies real-rooted)"
            } else {
                ""
            }
        );
        if let Some(reason) = &criterion.reason {
            println!("    {reason}");
        }
        for check in criterion.checks.iter().filter(|check| !check.holds) {
            println!(
                "    fail i={}: {} {} {} is false",
                check.index, check.lhs, check.comparison, check.rhs
            );
        }
    }
}

pub(crate) fn cmd_coefficient_tests(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let reports = read_polys_bigint()
        .into_iter()
        .map(|coeffs| coefficient_test_report_bigint(strip_trailing_zeros_bigint(&coeffs)))
        .collect::<Vec<_>>();
    if format == OutputFormat::Json {
        let items = reports
            .iter()
            .enumerate()
            .map(|(index, report)| coefficient_test_report_json(index, report))
            .collect::<Vec<_>>()
            .join(",");
        println!("{{\"items\":[{}]}}", items);
        return;
    }
    for (index, report) in reports.iter().enumerate() {
        print_coefficient_test_report_text(index, report);
    }
}
