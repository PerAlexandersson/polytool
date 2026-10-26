//! `polytool oeis`: list, describe and generate bundled OEIS families; add new
//! families to the catalog data file and re-check entries against OEIS data.

mod add;
mod catalog_file;
mod check;
mod sha256;
mod source;

use crate::input::parse_usize_option;
use num_bigint::BigInt;
use polytool::recurrence::*;
use polytool::*;
use serde_json::{json, Value};

pub(crate) fn print_oeis_help() {
    println!("Usage:");
    println!("  polytool oeis list [--json] [--include-experimental]");
    println!("  polytool oeis info <A-number> [--json]");
    println!("  polytool oeis generate <A-number> --rows <n> [options]");
    println!("  polytool oeis add <A-number>... [options]   (see `oeis add --help`)");
    println!("  polytool oeis check [A-number...] [options] (see `oeis check --help`)");
    println!("  polytool oeis --help");
    println!();
    println!("Generate OEIS polynomial rows from bundled exact recurrences.");
    println!();
    println!("Generation options:");
    println!("  --rows <n>             Number of complete rows to emit");
    println!("  --start-row <n>         First displayed OEIS row (default: catalog start)");
    println!("  --format <format>       rows, triangle, polynomial, json, jsonl, csv, bfile");
    println!("                          (default: polynomial)");
    println!("  --max-terms <n>         B-file cap; never splits a row");
    println!("  --include-experimental  Permit an experimental catalog entry");
    println!();
    println!("Catalog maintenance (source checkout):");
    println!("  oeis add     Search a recurrence for OEIS rows, check it on held-out terms,");
    println!("               and record it with provenance in data/oeis-catalog.jsonl");
    println!("  oeis check   Re-validate catalog entries against a local OEIS copy");
    println!();
    println!("Examples:");
    println!("  polytool oeis generate A008292 --rows 6 --format triangle");
    println!("  polytool oeis generate A008292 --rows 30 --format bfile");
    println!("  polytool oeis add A008292 --oeis-data ~/oeis --dry-run");
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum OeisOutputFormat {
    Rows,
    Triangle,
    Polynomial,
    Json,
    JsonLines,
    Csv,
    BFile,
}

impl OeisOutputFormat {
    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "rows" | "coefficients" | "text" => Ok(Self::Rows),
            "triangle" | "table" => Ok(Self::Triangle),
            "polynomial" | "polynomials" | "poly" => Ok(Self::Polynomial),
            "json" => Ok(Self::Json),
            "jsonl" | "json-lines" => Ok(Self::JsonLines),
            "csv" => Ok(Self::Csv),
            "bfile" | "b-file" => Ok(Self::BFile),
            other => Err(format!("unknown OEIS output format: {other}")),
        }
    }
}

fn oeis_entry_json(entry: &polytool::oeis::OeisSequenceDefinition) -> Value {
    json!({
        "id": entry.id,
        "name": entry.name,
        "status": entry.status.as_str(),
        "layout": entry.layout.as_str(),
        "first_row": entry.first_row,
        "flattened_offset": entry.flattened_offset,
        "bfile_available": entry.bfile_prefix_verified,
        "source_rows": entry.source_rows,
        "verification_rows": entry.verification_rows,
        "rows_sha256": entry.rows_sha256,
    })
}

fn cmd_oeis_list(args: &[String]) {
    let mut json_output = false;
    let mut include_experimental = false;
    for arg in args {
        match arg.as_str() {
            "--json" => json_output = true,
            "--include-experimental" => include_experimental = true,
            other => {
                eprintln!("unknown oeis list option: {other}");
                std::process::exit(2);
            }
        }
    }
    let entries = polytool::oeis::catalog()
        .iter()
        .filter(|entry| {
            include_experimental || entry.status != polytool::oeis::OeisSequenceStatus::Experimental
        })
        .collect::<Vec<_>>();
    if json_output {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "schema": "polytool.oeis-catalog.v1",
                "count": entries.len(),
                "sequences": entries.iter().map(|entry| oeis_entry_json(entry)).collect::<Vec<_>>(),
            }))
            .unwrap()
        );
    } else {
        for entry in entries {
            let bfile = if entry.bfile_prefix_verified {
                "bfile"
            } else {
                "no-bfile"
            };
            println!(
                "{}\t{}\t{}\t{}",
                entry.id,
                entry.status.as_str(),
                bfile,
                entry.name
            );
        }
    }
}

fn cmd_oeis_info(args: &[String]) {
    let mut json_output = false;
    let mut id: Option<&str> = None;
    for arg in args {
        match arg.as_str() {
            "--json" => json_output = true,
            other if id.is_none() => id = Some(other),
            other => {
                eprintln!("unknown oeis info option: {other}");
                std::process::exit(2);
            }
        }
    }
    let Some(id) = id else {
        eprintln!("oeis info needs an A-number");
        std::process::exit(2);
    };
    let Some(entry) = polytool::oeis::by_id(id) else {
        eprintln!("unknown bundled OEIS sequence: {id}");
        std::process::exit(2);
    };
    let (recurrence, initial_rows, first_index) =
        entry.recurrence_export_parts().unwrap_or_else(|error| {
            eprintln!("{error}");
            std::process::exit(2);
        });
    if json_output {
        let mut value = oeis_entry_json(entry);
        let object = value.as_object_mut().unwrap();
        object.insert("recurrence".to_string(), json!(recurrence.to_string()));
        let recurrence_first_index = usize::try_from(first_index).unwrap_or_else(|_| {
            eprintln!(
                "{} has a negative displayed recurrence index, which cannot be represented by {}",
                entry.id, RECURRENCE_JSON_SCHEMA
            );
            std::process::exit(2);
        });
        object.insert(
            "recurrence_data".to_string(),
            serde_json::to_value(RecurrenceJson::from_recurrence_rational(
                &recurrence,
                recurrence_first_index,
                &initial_rows,
                None,
            ))
            .unwrap(),
        );
        object.insert("latex".to_string(), json!(recurrence.to_latex()));
        object.insert(
            "mathematica".to_string(),
            json!(recurrence
                .to_mathematica_definition_rational_indexed(&initial_rows, first_index)),
        );
        object.insert(
            "sage".to_string(),
            json!(recurrence.to_sage_definition_rational_indexed(&initial_rows, first_index)),
        );
        object.insert(
            "python".to_string(),
            json!(recurrence.to_python_definition_rational_indexed(&initial_rows, first_index)),
        );
        println!("{}", serde_json::to_string_pretty(&value).unwrap());
    } else {
        println!("{}: {}", entry.id, entry.name);
        println!("status: {}", entry.status.as_str());
        println!("layout: {}", entry.layout.as_str());
        println!("first row: {}", entry.first_row);
        println!("flattened offset: {}", entry.flattened_offset);
        println!(
            "b-file: {}",
            if entry.bfile_prefix_verified {
                "available"
            } else {
                "disabled (prefix mapping not verified)"
            }
        );
        println!("held-out verification rows: {}", entry.verification_rows);
        println!("recurrence: {recurrence}");
    }
}

fn bigint_rows_json(rows: &[Vec<BigInt>]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect()
}

fn cmd_oeis_generate(args: &[String]) {
    let mut id: Option<&str> = None;
    let mut row_count: Option<usize> = None;
    let mut first_row: Option<i64> = None;
    let mut format = OeisOutputFormat::Polynomial;
    let mut max_terms: Option<usize> = None;
    let mut include_experimental = false;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--rows" => {
                row_count = Some(
                    parse_usize_option(args, &mut index, "--rows").unwrap_or_else(|error| {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }),
                );
            }
            "--start-row" => {
                index += 1;
                let Some(value) = args.get(index) else {
                    eprintln!("--start-row expects a value");
                    std::process::exit(2);
                };
                first_row = Some(value.parse().unwrap_or_else(|_| {
                    eprintln!("--start-row expects an integer, got '{value}'");
                    std::process::exit(2);
                }));
            }
            "--format" => {
                index += 1;
                let Some(value) = args.get(index) else {
                    eprintln!("--format expects a value");
                    std::process::exit(2);
                };
                format = OeisOutputFormat::parse(value).unwrap_or_else(|error| {
                    eprintln!("{error}");
                    std::process::exit(2);
                });
            }
            "--json" => format = OeisOutputFormat::Json,
            "--max-terms" => {
                max_terms = Some(
                    parse_usize_option(args, &mut index, "--max-terms").unwrap_or_else(|error| {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }),
                );
            }
            "--include-experimental" => include_experimental = true,
            other if id.is_none() => id = Some(other),
            other => {
                eprintln!("unknown oeis generate option: {other}");
                std::process::exit(2);
            }
        }
        index += 1;
    }
    let Some(id) = id else {
        eprintln!("oeis generate needs an A-number");
        std::process::exit(2);
    };
    let Some(row_count) = row_count else {
        eprintln!("oeis generate needs --rows <n>");
        std::process::exit(2);
    };
    let Some(entry) = polytool::oeis::by_id(id) else {
        eprintln!("unknown bundled OEIS sequence: {id}");
        std::process::exit(2);
    };
    if entry.status == polytool::oeis::OeisSequenceStatus::Experimental && !include_experimental {
        eprintln!(
            "{} is experimental; pass --include-experimental to generate it",
            entry.id
        );
        std::process::exit(2);
    }
    let first_row = first_row.unwrap_or(entry.first_row);
    if format == OeisOutputFormat::BFile && first_row != entry.first_row {
        eprintln!("strict b-file output must begin at row {}", entry.first_row);
        std::process::exit(2);
    }
    if format != OeisOutputFormat::BFile && max_terms.is_some() {
        eprintln!("--max-terms is only valid with --format bfile");
        std::process::exit(2);
    }
    let rows = entry
        .generate_rows_from(first_row, row_count)
        .unwrap_or_else(|error| {
            eprintln!("{error}");
            std::process::exit(2);
        });

    match format {
        OeisOutputFormat::Rows => {
            for row in &rows {
                println!(
                    "{}",
                    row.iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
        }
        OeisOutputFormat::Triangle => {
            for (offset, row) in rows.iter().enumerate() {
                println!(
                    "{}: {}",
                    first_row + offset as i64,
                    row.iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(" ")
                );
            }
        }
        OeisOutputFormat::Polynomial => {
            for (offset, row) in rows.iter().enumerate() {
                println!(
                    "P_{}(t) = {}",
                    first_row + offset as i64,
                    format_poly_bigint_coeffs(row)
                );
            }
        }
        OeisOutputFormat::Json => {
            println!(
                "{}",
                serde_json::to_string_pretty(&json!({
                    "schema": "polytool.oeis-rows.v1",
                    "id": entry.id,
                    "first_row": first_row,
                    "row_count": rows.len(),
                    "rows": bigint_rows_json(&rows),
                }))
                .unwrap()
            );
        }
        OeisOutputFormat::JsonLines => {
            for (offset, row) in rows.iter().enumerate() {
                println!(
                    "{}",
                    json!({
                        "id": entry.id,
                        "n": first_row + offset as i64,
                        "coefficients": row.iter().map(ToString::to_string).collect::<Vec<_>>(),
                    })
                );
            }
        }
        OeisOutputFormat::Csv => {
            println!("n,k,value");
            for (row_offset, row) in rows.iter().enumerate() {
                for (column, value) in row.iter().enumerate() {
                    println!("{},{},{}", first_row + row_offset as i64, column, value);
                }
            }
        }
        OeisOutputFormat::BFile => {
            if !entry.bfile_prefix_verified {
                eprintln!(
                    "strict b-file output for {} is disabled because its OEIS prefix mapping is not verified",
                    entry.id
                );
                std::process::exit(2);
            }
            let mut included_terms = 0usize;
            let mut flattened_index = entry.flattened_offset;
            let mut included_rows = 0usize;
            for row in &rows {
                if max_terms.is_some_and(|limit| included_terms + row.len() > limit) {
                    break;
                }
                for value in row {
                    println!("{flattened_index} {value}");
                    flattened_index += 1;
                }
                included_terms += row.len();
                included_rows += 1;
            }
            if row_count > 0 && included_rows == 0 {
                eprintln!("--max-terms is smaller than the first complete row");
                std::process::exit(2);
            }
        }
    }
}

pub(crate) fn cmd_oeis(args: &[String]) {
    let Some((subcommand, rest)) = args.split_first() else {
        print_oeis_help();
        std::process::exit(2);
    };
    match subcommand.as_str() {
        "list" => cmd_oeis_list(rest),
        "info" => cmd_oeis_info(rest),
        "generate" => cmd_oeis_generate(rest),
        "add" => add::cmd_oeis_add(rest),
        "check" => check::cmd_oeis_check(rest),
        "--help" | "-h" | "help" => print_oeis_help(),
        other => {
            eprintln!("unknown oeis subcommand: {other}");
            print_oeis_help();
            std::process::exit(2);
        }
    }
}
