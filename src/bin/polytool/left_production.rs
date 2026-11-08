//! `polytool left-production`: the Chen--Fu--Ruan left production matrix of a
//! triangle, whose total nonnegativity certifies the triangle, its reversal
//! and Pólya frequency rows.

use crate::input::{parse_usize_option, read_polys_bigint};
use crate::json_out::{json_string, json_string_vec, json_usize_vec};
use num_bigint::BigInt;
use polytool::left_production::{ProductionNetwork, RationalNegativeMinor};
use polytool::*;

/// Fixed guardrail printed with every report.
pub(crate) const LEFT_PRODUCTION_NOTE: &str = "If the left production matrix Q_N is TNN, \
     Chen--Fu--Ruan (arXiv:2512.08369, Theorems 3 and 4) give one planar network whose \
     path matrices are the triangle, its reversal and the Toeplitz matrix of each row, so \
     the rows P_0, ..., P_N are Pólya frequency sequences (real-rooted, nonpositive zeros) \
     and the triangle and its reversal are TNN, for this prefix only. If Q_N is not TNN, \
     nothing follows about the rows (the Eulerian rows are real-rooted, yet q_{4,2} = -5).";

struct LeftProductionArgs {
    json: bool,
    first_row_index: usize,
}

fn parse_left_production_args(args: &[String]) -> Result<LeftProductionArgs, String> {
    let mut parsed = LeftProductionArgs {
        json: false,
        first_row_index: 1,
    };
    let mut index = 0;
    while index < args.len() {
        let arg = args[index].as_str();
        match arg {
            "--json" => parsed.json = true,
            "--first-index" => parsed.first_row_index = parse_usize_option(args, &mut index, arg)?,
            other => return Err(format!("unknown left-production option: {other}")),
        }
        index += 1;
    }
    Ok(parsed)
}

fn strings<T: ToString>(row: &[T]) -> Vec<String> {
    row.iter().map(ToString::to_string).collect()
}

fn matrix_json<T: ToString>(rows: &[Vec<T>]) -> String {
    format!(
        "[{}]",
        rows.iter()
            .map(|row| json_string_vec(&strings(row)))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn minor_json(minor: &RationalNegativeMinor) -> String {
    format!(
        "{{\"rows\":{},\"cols\":{},\"value\":{}}}",
        json_usize_vec(&minor.rows),
        json_usize_vec(&minor.cols),
        json_string(&minor.value.to_string())
    )
}

fn network_json(network: &ProductionNetwork) -> String {
    format!(
        "{{\"lambda\":{},\"column_scales\":{}}}",
        matrix_json(&network.network.lambda),
        json_string_vec(&strings(&network.column_scales))
    )
}

fn json_option<T>(value: Option<T>, render: impl FnOnce(T) -> String) -> String {
    value.map_or("null".to_string(), render)
}

fn report_json(report: &LeftProductionReport, first_row_index: usize, rows: usize) -> String {
    format!(
        "{{\"schema\":\"polytool.left-production.v1\",\"first_row_index\":{},\"first_index\":{},\"rows\":{},\
         \"triangle\":{},\"production_matrix\":{},\"tnn\":{},\"failure_reason\":{},\
         \"negative_minor\":{},\"certificate_verified\":{},\"rows_certified_pf\":{},\
         \"production_network\":{},\"message\":{},\"note\":{}}}",
        first_row_index,
        first_row_index,
        rows,
        matrix_json(&report.triangle),
        matrix_json(&report.production_matrix),
        report.totally_nonnegative,
        json_option(report.failure_reason.as_deref(), json_string),
        json_option(report.negative_minor.as_ref(), minor_json),
        json_option(report.certificate_verified, |v| v.to_string()),
        report.rows_certified_pf(),
        json_option(report.production_network.as_ref(), network_json),
        json_string(&report.to_string()),
        json_string(LEFT_PRODUCTION_NOTE)
    )
}

fn print_report_text(report: &LeftProductionReport, first_row_index: usize, rows: usize) {
    println!("Left production matrix (Chen--Fu--Ruan)");
    println!(
        "- Input: {rows} row(s){}",
        if first_row_index == 1 {
            "; P_0 = 1 prepended"
        } else {
            ""
        }
    );
    println!("- Production matrix Q_N (triangular rows):");
    for (n, row) in report.production_matrix.iter().enumerate() {
        println!("  {n}: {}", strings(row).join(", "));
    }
    let last = report.triangle.len() - 1;
    if report.rows_certified_pf() {
        println!("- Verdict: Q_{last} is TNN; certificate re-verified");
        println!(
            "- Conclusion: rows P_0..P_{last} are PF (real-rooted, nonpositive zeros); the \
             triangle and its reversal are TNN"
        );
    } else if report.totally_nonnegative {
        println!("- Verdict: Q_{last} is TNN, but the certificate FAILED re-verification");
    } else {
        println!("- Verdict: Q_{last} is not TNN; nothing follows about the rows");
        if let Some(reason) = &report.failure_reason {
            println!("- Reason: {reason}");
        }
        if let Some(minor) = &report.negative_minor {
            println!(
                "- Witness: minor with rows {:?} and columns {:?} equals {}",
                minor.rows, minor.cols, minor.value
            );
        }
    }
    if let Some(network) = &report.production_network {
        println!("- Planar network for Q_N = L D (lambda[n][k] of Gamma_N(lambda) for L):");
        for (n, row) in network.network.lambda.iter().enumerate() {
            println!("  n={n}: {}", strings(row).join(", "));
        }
        println!(
            "  column scales D = diag({})",
            strings(&network.column_scales).join(", ")
        );
    }
    println!();
    println!("Note: {LEFT_PRODUCTION_NOTE}");
}

pub(crate) fn cmd_left_production(args: &[String]) {
    let parsed = match parse_left_production_args(args) {
        Ok(parsed) => parsed,
        Err(error) => {
            eprintln!("{error}");
            eprintln!("Run `polytool help left-production` for usage.");
            std::process::exit(2);
        }
    };
    let polys: Vec<Polynomial<BigInt>> = read_polys_bigint()
        .into_iter()
        .map(Polynomial::new)
        .collect();
    let report = match check_left_production_matrix(&polys, parsed.first_row_index) {
        Ok(report) => report,
        Err(error) => {
            eprintln!("left-production: {error}");
            std::process::exit(2);
        }
    };
    if parsed.json {
        println!(
            "{}",
            report_json(&report, parsed.first_row_index, polys.len())
        );
    } else {
        print_report_text(&report, parsed.first_row_index, polys.len());
    }
    if !report.rows_certified_pf() {
        std::process::exit(1);
    }
}
