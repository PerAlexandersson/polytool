//! `polytool cyclic-sieving-sequence`: candidate cyclic sieving orders for a
//! polynomial sequence.

use crate::cyclic_sieving::{
    cyclic_sieving_report_json, parse_bigint_values, print_cyclic_sieving_report_text,
};
use crate::help::print_coefficient_input_help;
use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint, OutputFormat};
use crate::json_out::{json_bigint_vec, json_string, json_string_vec};
use num_bigint::BigInt;
use polytool::*;
use std::collections::BTreeMap;
use std::fs;

pub(crate) fn print_cyclic_sieving_sequence_help() {
    println!("Usage:");
    println!("  polytool cyclic-sieving-sequence [options]");
    println!("  polytool cyclic-sieving-sequence --help");
    println!();
    println!("For a polynomial sequence P_n(q), profile candidate cyclic group orders");
    println!("n-2, n-1, n, n+1, n+2, n+3 by default.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --first-index <n>         Index of the first input row (default: 0)");
    println!("  --offsets <list>          Candidate order offsets (default: -2,-1,0,1,2,3)");
    println!("  --fixed-counts-file <f>   Lines: <index> <order>: c_0,...,c_{{order-1}}");
    println!("  --json                    Emit machine-readable JSON");
}

#[derive(Clone, Debug)]
struct CyclicSievingSequenceOptions {
    format: OutputFormat,
    first_index: isize,
    offsets: Vec<isize>,
    fixed_counts: BTreeMap<(isize, usize), Vec<BigInt>>,
}

fn parse_offsets(input: &str) -> Result<Vec<isize>, String> {
    let mut offsets = Vec::new();
    for raw in input.split(',') {
        let part = raw.trim();
        if part.is_empty() {
            continue;
        }
        offsets.push(
            part.parse()
                .map_err(|_| format!("invalid offset '{part}'"))?,
        );
    }
    if offsets.is_empty() {
        return Err("offset list may not be empty".to_string());
    }
    Ok(offsets)
}

fn parse_fixed_counts_file(path: &str) -> Result<BTreeMap<(isize, usize), Vec<BigInt>>, String> {
    let text = fs::read_to_string(path).map_err(|e| format!("cannot read {path}: {e}"))?;
    let mut out = BTreeMap::new();
    for (line_index, raw_line) in text.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let normalized = line.replace(':', " ");
        let mut parts = normalized.split_whitespace();
        let index_text = parts
            .next()
            .ok_or_else(|| format!("line {}: expected index", line_index + 1))?;
        let order_text = parts
            .next()
            .ok_or_else(|| format!("line {}: expected order", line_index + 1))?;
        let index = index_text
            .parse::<isize>()
            .map_err(|_| format!("line {}: invalid index '{index_text}'", line_index + 1))?;
        let order = order_text
            .parse::<usize>()
            .map_err(|_| format!("line {}: invalid order '{order_text}'", line_index + 1))?;
        let counts_text = parts.collect::<Vec<_>>().join(" ");
        let counts = parse_bigint_values(&counts_text)
            .map_err(|e| format!("line {}: {e}", line_index + 1))?;
        if counts.len() != order {
            return Err(format!(
                "line {}: order {order} expects {order} counts, got {}",
                line_index + 1,
                counts.len()
            ));
        }
        out.insert((index, order), counts);
    }
    Ok(out)
}

fn parse_cyclic_sieving_sequence_args(
    args: &[String],
) -> Result<CyclicSievingSequenceOptions, String> {
    let mut options = CyclicSievingSequenceOptions {
        format: OutputFormat::Text,
        first_index: 0,
        offsets: vec![-2, -1, 0, 1, 2, 3],
        fixed_counts: BTreeMap::new(),
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => options.format = OutputFormat::Json,
            "--first-index" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--first-index expects an integer".to_string())?;
                options.first_index = value
                    .parse()
                    .map_err(|_| format!("--first-index expects an integer, got '{value}'"))?;
            }
            "--offsets" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--offsets expects a comma-separated list".to_string())?;
                options.offsets = parse_offsets(value)?;
            }
            "--fixed-counts-file" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--fixed-counts-file expects a path".to_string())?;
                options.fixed_counts = parse_fixed_counts_file(value)?;
            }
            other => return Err(format!("unknown cyclic-sieving-sequence option: {other}")),
        }
        i += 1;
    }
    Ok(options)
}

fn cyclic_sieving_sequence_item_json(item: &CyclicSievingSequenceItem) -> String {
    let reports = item
        .candidate_orders
        .iter()
        .map(cyclic_sieving_report_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"row\":{},\"index\":{},\"coefficients\":{},\"polynomial\":{},\
         \"candidate_orders\":[{}]}}",
        item.row,
        item.index,
        json_bigint_vec(&item.coefficients),
        json_string(&format_poly_bigint_coeffs(&item.coefficients)),
        reports
    )
}

pub(crate) fn cmd_cyclic_sieving_sequence(args: &[String]) {
    let options = match parse_cyclic_sieving_sequence_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let polynomials = read_polys_bigint()
        .into_iter()
        .map(|coeffs| strip_trailing_zeros_bigint(&coeffs).to_vec())
        .collect::<Vec<_>>();
    let reports = cyclic_sieving_sequence_reports_bigint(
        &polynomials,
        options.first_index,
        &options.offsets,
        &options.fixed_counts,
    );
    if options.format == OutputFormat::Json {
        let items = reports
            .iter()
            .map(cyclic_sieving_sequence_item_json)
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{{\"first_index\":{},\"offsets\":{},\"items\":[{}]}}",
            options.first_index,
            json_string_vec(
                &options
                    .offsets
                    .iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
            ),
            items
        );
        return;
    }
    for item in &reports {
        println!(
            "row {}, index {}: {}",
            item.row,
            item.index,
            format_poly_bigint_coeffs(&item.coefficients)
        );
        for report in &item.candidate_orders {
            print_cyclic_sieving_report_text(item.row, report);
        }
    }
}
