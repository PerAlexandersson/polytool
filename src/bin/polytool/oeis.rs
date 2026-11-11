//! `polytool oeis`: list, describe and generate bundled OEIS families; add new
//! families to the catalog data file and re-check entries against OEIS data.

mod add;
mod catalog_file;
mod check;
mod evidence;
mod rows;
mod sha256;
mod source;
mod supplied;

use crate::input::parse_usize_option;
use num_bigint::BigInt;
use polytool::recurrence::*;
use polytool::*;
use serde_json::{json, Value};

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
        "first_index": entry.first_row,
        "flattened_offset": entry.flattened_offset,
        "bfile_available": entry.bfile_prefix_verified,
        "source_rows": entry.source_rows,
        "verification_rows": entry.verification_rows,
        "rows_sha256": entry.rows_sha256,
        "real_rootedness": entry.real_rootedness().map(|flag| json!({
            "status": flag.status.as_str(),
            "as_of": flag.as_of,
        })),
        "evidence": entry.evidence().map(|evidence| evidence::evidence_json(
            evidence.strongest,
            evidence.rows_checked,
            &evidence.first_failure,
            evidence.as_of,
        )),
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
        match entry.real_rootedness() {
            Some(flag) => println!(
                "real-rootedness: {}{} (recorded flag, not proved by polytool)",
                flag.status.as_str(),
                flag.as_of
                    .map_or_else(String::new, |date| format!(" as of {date}"))
            ),
            None => println!("real-rootedness: unknown"),
        }
        match entry.evidence() {
            Some(evidence) => println!(
                "{}",
                evidence::describe_evidence(
                    evidence.strongest,
                    evidence.rows_checked,
                    &evidence.first_failure,
                    evidence.as_of,
                )
            ),
            None => println!("evidence: none recorded"),
        }
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
            "--count" => {
                row_count = Some(
                    parse_usize_option(args, &mut index, "--count").unwrap_or_else(|error| {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }),
                );
            }
            "--from-index" => {
                index += 1;
                let Some(value) = args.get(index) else {
                    eprintln!("--from-index expects a value");
                    std::process::exit(2);
                };
                first_row = Some(value.parse().unwrap_or_else(|_| {
                    eprintln!("--from-index expects an integer, got '{value}'");
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
        eprintln!("oeis generate needs --count <n>");
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
                    "first_index": first_row,
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
            let options = polytool::oeis::BFileOptions {
                first_index: entry.flattened_offset,
                max_terms,
                ..polytool::oeis::BFileOptions::default()
            };
            let file =
                polytool::oeis::flatten_rows_to_bfile(&rows, &options).unwrap_or_else(|error| {
                    eprintln!("{error}");
                    std::process::exit(2);
                });
            let included_rows = file.rows_included;
            print!("{}", file.text);
            if row_count > 0 && included_rows == 0 {
                eprintln!("--max-terms is smaller than the first complete row");
                std::process::exit(2);
            }
        }
    }
}

pub(crate) fn cmd_oeis(args: &[String]) {
    let Some((subcommand, rest)) = args.split_first() else {
        crate::input::usage_error("oeis", "expected an OEIS subcommand");
    };
    match subcommand.as_str() {
        "rows" => rows::run(rest),
        "list" => cmd_oeis_list(rest),
        "info" => cmd_oeis_info(rest),
        "generate" => cmd_oeis_generate(rest),
        "add" => add::cmd_oeis_add(rest),
        "check" => check::cmd_oeis_check(rest),
        "evidence" => evidence::cmd_oeis_evidence(rest),
        other => {
            crate::input::usage_error("oeis", &format!("unknown oeis subcommand: {other}"));
        }
    }
}
