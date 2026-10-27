//! `polytool hstar-to-ehrhart`, `ehrhart-to-hstar` and `hstar-inequalities`.

use crate::help::print_coefficient_input_help;
use crate::input::{
    read_polys_bigint, read_polys_rational, strip_trailing_zeros_bigint, OutputFormat,
};
use crate::json_out::{json_bigint_vec, json_string, json_string_vec};
use num_bigint::BigInt;
use polytool::*;

pub(crate) fn print_hstar_inequalities_help() {
    println!("Usage:");
    println!("  polytool hstar-inequalities [options]");
    println!("  polytool hstar-inequalities --help");
    println!();
    println!("Check named Ehrhart h*-vector inequalities with exact integer arithmetic.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --dimension <d>  Use dimension d instead of len(h*)-1 (at most 4096)");
    println!("  --json           Emit machine-readable JSON");
    println!();
    println!("Example:");
    println!("  echo '1, 4, 1, 0' | polytool hstar-inequalities --dimension 3");
}

pub(crate) fn cmd_hstar_to_ehrhart(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_bigint().into_iter().enumerate() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        let ehrhart = hstar_to_ehrhart_bigint_coeffs(c);
        let display: Vec<String> = ehrhart.iter().map(|r| format!("{}", r)).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"hstar\":{},\"ehrhart_coefficients\":{}}}",
                index,
                json_bigint_vec(c),
                json_string_vec(&display)
            ));
        } else {
            println!(
                "h*={} => L(n) coeffs: [{}]",
                format_poly_bigint_coeffs(c),
                display.join(", ")
            );
        }
    }
    if format == OutputFormat::Json {
        println!("{{\"items\":[{}]}}", json_items.join(","));
    }
}

pub(crate) fn cmd_ehrhart_to_hstar(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_rational().into_iter().enumerate() {
        let hstar = match ehrhart_to_hstar_bigint(&coeffs) {
            Ok(hstar) => hstar,
            Err(error) => {
                eprintln!("Invalid Ehrhart polynomial at input {index}: {error}");
                std::process::exit(2);
            }
        };
        let display: Vec<String> = coeffs.iter().map(ToString::to_string).collect();
        if format == OutputFormat::Json {
            json_items.push(format!(
                "{{\"index\":{},\"ehrhart_coefficients\":{},\"hstar\":{}}}",
                index,
                json_string_vec(&display),
                json_bigint_vec(&hstar)
            ));
        } else {
            println!(
                "L(n) coeffs: [{}] => h*={}",
                display.join(", "),
                format_poly_bigint_coeffs(&hstar)
            );
        }
    }
    if format == OutputFormat::Json {
        println!("{{\"items\":[{}]}}", json_items.join(","));
    }
}

#[derive(Clone, Debug)]
struct HStarInequalityOptions {
    format: OutputFormat,
    dimension: Option<usize>,
}

fn parse_hstar_inequality_args(args: &[String]) -> Result<HStarInequalityOptions, String> {
    let mut options = HStarInequalityOptions {
        format: OutputFormat::Text,
        dimension: None,
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => options.format = OutputFormat::Json,
            "--dimension" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--dimension expects a nonnegative integer".to_string())?;
                let dimension: usize = value.parse().map_err(|_| {
                    format!("--dimension expects a nonnegative integer, got '{value}'")
                })?;
                if dimension > HSTAR_INEQUALITIES_MAX_DIMENSION {
                    return Err(format!(
                        "--dimension {dimension} exceeds the maximum \
                         {HSTAR_INEQUALITIES_MAX_DIMENSION}"
                    ));
                }
                options.dimension = Some(dimension);
            }
            other => return Err(format!("unknown hstar-inequalities option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn json_bigint_option(value: Option<&BigInt>) -> String {
    value
        .map(|value| json_string(&value.to_string()))
        .unwrap_or_else(|| "null".to_string())
}

fn hstar_check_json(check: &HStarInequalityCheck) -> String {
    format!(
        "{{\"family\":{},\"name\":{},\"formula\":{},\"reference\":{},\"url\":{},\
         \"index\":{},\"applicable\":{},\"holds\":{},\"lhs\":{},\"rhs\":{},\
         \"value\":{},\"details\":{}}}",
        json_string(&check.family),
        json_string(&check.name),
        json_string(&check.formula),
        json_string(&check.reference),
        check
            .url
            .as_ref()
            .map(|url| json_string(url))
            .unwrap_or_else(|| "null".to_string()),
        check
            .index
            .map(|index| index.to_string())
            .unwrap_or_else(|| "null".to_string()),
        check.applicable,
        check.holds,
        json_bigint_option(check.lhs.as_ref()),
        json_bigint_option(check.rhs.as_ref()),
        check
            .value
            .as_ref()
            .map(|value| json_string(value))
            .unwrap_or_else(|| "null".to_string()),
        json_string(&check.details)
    )
}

fn hstar_report_json(index: usize, report: &HStarInequalityReport) -> String {
    let checks = report
        .checks
        .iter()
        .map(hstar_check_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"index\":{},\"hstar\":{},\"dimension\":{},\"degree\":{},\"codegree\":{},\
         \"all_applicable_hold\":{},\"checks\":[{}]}}",
        index,
        json_bigint_vec(&report.hstar),
        report.dimension,
        report.degree,
        report
            .codegree
            .map(|codegree| codegree.to_string())
            .unwrap_or_else(|| "null".to_string()),
        report.all_applicable_hold,
        checks
    )
}

fn print_hstar_report_text(index: usize, report: &HStarInequalityReport) {
    println!(
        "row {index}: h*={}, d={}, s={}, codegree={}, {}",
        format_poly_bigint_coeffs(&report.hstar),
        report.dimension,
        report.degree,
        report
            .codegree
            .map(|c| c.to_string())
            .unwrap_or_else(|| "n/a".to_string()),
        if report.all_applicable_hold {
            "all applicable checks passed"
        } else {
            "FAILED"
        }
    );
    for check in report
        .checks
        .iter()
        .filter(|check| check.applicable && !check.holds)
    {
        let index = check
            .index
            .map(|index| format!(", i={index}"))
            .unwrap_or_default();
        let lhs_rhs = match (&check.lhs, &check.rhs) {
            (Some(lhs), Some(rhs)) => format!(" lhs={lhs}, rhs={rhs};"),
            _ => check
                .value
                .as_ref()
                .map(|value| format!(" value={value};"))
                .unwrap_or_default(),
        };
        println!(
            "  FAIL {} / {}{}: {} ;{} reference: {}",
            check.family, check.name, index, check.formula, lhs_rhs, check.reference
        );
        if !check.details.is_empty() {
            println!("    {}", check.details);
        }
    }
}

pub(crate) fn cmd_hstar_inequalities(args: &[String]) {
    let options = match parse_hstar_inequality_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let reports = read_polys_bigint()
        .into_iter()
        .map(|coeffs| {
            let c = strip_trailing_zeros_bigint(&coeffs);
            hstar_inequality_report_bigint(c, options.dimension)
        })
        .collect::<Vec<_>>();
    if options.format == OutputFormat::Json {
        let items = reports
            .iter()
            .enumerate()
            .map(|(index, report)| hstar_report_json(index, report))
            .collect::<Vec<_>>()
            .join(",");
        println!("{{\"items\":[{}]}}", items);
        return;
    }
    for (index, report) in reports.iter().enumerate() {
        print_hstar_report_text(index, report);
    }
}
