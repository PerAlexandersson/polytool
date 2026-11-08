//! `polytool bicyclic-sieving`: bicyclic sieving for C_n1 x C_n2.

use crate::input::{read_stdin_or_exit, OutputFormat};
use crate::json_out::{json_bigint_vec, json_bool_option, json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;
use serde_json::Value;

fn parse_bigint_json_value(value: &Value) -> Result<BigInt, String> {
    match value {
        Value::Number(number) => {
            if let Some(v) = number.as_i64() {
                Ok(BigInt::from(v))
            } else if let Some(v) = number.as_u64() {
                Ok(BigInt::from(v))
            } else {
                Err(format!(
                    "matrix entry {number} is not an integer; write large integers as strings"
                ))
            }
        }
        Value::String(text) => text
            .trim()
            .parse::<BigInt>()
            .map_err(|_| format!("matrix entry '{text}' is not an integer")),
        other => Err(format!("matrix entry {other} is not an integer")),
    }
}

fn parse_bigint_matrix_json(input: &str, what: &str) -> Result<Vec<Vec<BigInt>>, String> {
    let value: Value = serde_json::from_str(input)
        .map_err(|error| format!("{what} is not valid JSON: {error}"))?;
    let rows = value
        .as_array()
        .ok_or_else(|| format!("{what} must be a JSON array of rows"))?;
    rows.iter()
        .map(|row| {
            row.as_array()
                .ok_or_else(|| format!("{what} must be a JSON array of arrays"))?
                .iter()
                .map(parse_bigint_json_value)
                .collect::<Result<Vec<_>, _>>()
        })
        .collect()
}

fn json_bigint_matrix(matrix: &[Vec<BigInt>]) -> String {
    format!(
        "[{}]",
        matrix
            .iter()
            .map(|row| json_bigint_vec(row))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn json_optional_bigint_matrix(matrix: &[Vec<Option<BigInt>>]) -> String {
    let rows = matrix
        .iter()
        .map(|row| {
            let cells = row
                .iter()
                .map(|cell| {
                    cell.as_ref()
                        .map(|v| json_string(&v.to_string()))
                        .unwrap_or_else(|| "null".to_string())
                })
                .collect::<Vec<_>>()
                .join(",");
            format!("[{cells}]")
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{rows}]")
}

fn bicyclic_sieving_report_json(report: &BicyclicSievingReport) -> String {
    let mismatches = report
        .mismatches
        .iter()
        .map(|(d1, d2)| format!("[{d1},{d2}]"))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"order1\":{},\"order2\":{},\"coefficient_matrix\":{},\
         \"reduced_matrix\":{},\"evaluations\":{},\"integral\":{},\
         \"natural_coefficients\":{},\"divisors1\":{},\"divisors2\":{},\
         \"product_basis_form\":{},\"product_coordinates\":{},\
         \"product_formula_evaluations\":{},\"fixed_table\":{},\"mismatches\":[{}],\
         \"holds\":{}}}",
        report.order1,
        report.order2,
        json_bigint_matrix(&report.coefficient_matrix),
        json_bigint_matrix(&report.reduced_matrix),
        json_optional_bigint_matrix(&report.evaluations),
        report.integral,
        report.natural_coefficients,
        json_usize_vec(&report.divisors1),
        json_usize_vec(&report.divisors2),
        report.product_basis_form,
        report
            .product_coordinates
            .as_deref()
            .map(json_bigint_matrix)
            .unwrap_or_else(|| "null".to_string()),
        report
            .product_formula_evaluations
            .as_deref()
            .map(json_bigint_matrix)
            .unwrap_or_else(|| "null".to_string()),
        report
            .fixed_table
            .as_deref()
            .map(json_bigint_matrix)
            .unwrap_or_else(|| "null".to_string()),
        mismatches,
        json_bool_option(report.holds)
    )
}

fn print_bicyclic_sieving_report_text(report: &BicyclicSievingReport) {
    println!(
        "C_{} x C_{}: {}",
        report.order1,
        report.order2,
        match report.holds {
            Some(true) => "bicyclic sieving checks passed",
            Some(false) => "bicyclic sieving checks failed",
            None => "profile only",
        }
    );
    println!(
        "  integral at all roots of unity: {}; natural coefficients: {}; product basis form: {}",
        if report.integral { "yes" } else { "no" },
        if report.natural_coefficients {
            "yes"
        } else {
            "no"
        },
        if report.product_basis_form {
            "yes"
        } else {
            "no"
        }
    );
    for (d1, row) in report.evaluations.iter().enumerate() {
        let cells = row
            .iter()
            .map(|cell| {
                cell.as_ref()
                    .map(ToString::to_string)
                    .unwrap_or_else(|| "non-integer".to_string())
            })
            .collect::<Vec<_>>()
            .join(", ");
        println!("  d1={d1}: f(zeta1^d1, zeta2^d2) for d2=0.. : {cells}");
    }
    if let Some(coords) = &report.product_coordinates {
        println!(
            "  product coordinates a_(r1,r2) (rows r1 in {:?}, columns r2 in {:?}); not orbit counts:",
            report.divisors1, report.divisors2
        );
        for row in coords {
            let cells = row
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(", ");
            println!("    {cells}");
        }
    }
    for (d1, d2) in &report.mismatches {
        let expected = report
            .fixed_table
            .as_ref()
            .map(|t| t[*d1][*d2].to_string())
            .unwrap_or_default();
        let got = report.evaluations[*d1][*d2]
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_else(|| "non-integer".to_string());
        println!("  FAIL (d1,d2)=({d1},{d2}): expected {expected}, got {got}");
    }
}

pub(crate) fn cmd_bicyclic_sieving(args: &[String]) {
    let mut format = OutputFormat::Text;
    let mut order1: Option<usize> = None;
    let mut order2: Option<usize> = None;
    let mut fixed_table: Option<Vec<Vec<BigInt>>> = None;
    let mut i = 0;
    while i < args.len() {
        let parse_order = |value: Option<&String>, name: &str| -> Result<usize, String> {
            let value = value.ok_or_else(|| format!("{name} expects a positive integer"))?;
            match value.parse::<usize>() {
                Ok(0) | Err(_) => Err(format!("{name} expects a positive integer, got '{value}'")),
                Ok(n) => Ok(n),
            }
        };
        let result = match args[i].as_str() {
            "--json" => {
                format = OutputFormat::Json;
                Ok(())
            }
            "--order1" => {
                i += 1;
                parse_order(args.get(i), "--order1").map(|n| order1 = Some(n))
            }
            "--order2" => {
                i += 1;
                parse_order(args.get(i), "--order2").map(|n| order2 = Some(n))
            }
            "--fixed-table" => {
                i += 1;
                match args.get(i) {
                    Some(value) => parse_bigint_matrix_json(value, "--fixed-table")
                        .map(|table| fixed_table = Some(table)),
                    None => Err("--fixed-table expects a JSON matrix".to_string()),
                }
            }
            other => Err(format!("unknown bicyclic-sieving option: {other}")),
        };
        if let Err(error) = result {
            eprintln!("{error}");
            std::process::exit(2);
        }
        i += 1;
    }
    let (Some(order1), Some(order2)) = (order1, order2) else {
        eprintln!("bicyclic-sieving requires --order1 <n1> and --order2 <n2>");
        std::process::exit(2);
    };
    let matrix = match parse_bigint_matrix_json(&read_stdin_or_exit(), "stdin") {
        Ok(matrix) => matrix,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };
    match bicyclic_sieving_report_bigint(&matrix, order1, order2, fixed_table.as_deref()) {
        Ok(report) => {
            if format == OutputFormat::Json {
                println!("{{\"report\":{}}}", bicyclic_sieving_report_json(&report));
            } else {
                print_bicyclic_sieving_report_text(&report);
            }
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}
