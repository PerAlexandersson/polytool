//! `polytool oeis check`: re-validate catalog entries against local OEIS data.

use super::add::{normalize_id, read_ids_file, OEIS_DATA_ENV};
use super::catalog_file::{generate_flattened, leak_definition, load_catalog};
use super::source::read_local_record;
use crate::input::next_option_value;
use polytool::oeis::OeisSequenceDefinition;
use serde_json::json;
use std::path::PathBuf;

/// Upper bound on generated rows per entry while covering the OEIS terms.
const MAX_CHECK_ROWS: usize = 4096;

fn fail<T>(error: String) -> T {
    eprintln!("{error}");
    std::process::exit(2);
}

struct CheckResult {
    id: &'static str,
    outcome: &'static str,
    checked_terms: usize,
    detail: String,
}

fn check_entry(root: &std::path::Path, entry: &'static OeisSequenceDefinition) -> CheckResult {
    let result = |outcome, checked_terms, detail: String| CheckResult {
        id: entry.id,
        outcome,
        checked_terms,
        detail,
    };
    if !entry.bfile_prefix_verified {
        return result(
            "unchecked",
            0,
            "flattened OEIS alignment is not verified for this entry".to_string(),
        );
    }
    let source = match read_local_record(root, entry.id) {
        Ok(source) => source,
        Err(error) => return result("unchecked", 0, error),
    };
    let Ok(start) = usize::try_from(entry.flattened_offset - source.offset) else {
        return result(
            "mismatch",
            0,
            format!(
                "catalog flattened offset {} precedes the OEIS offset {}",
                entry.flattened_offset, source.offset
            ),
        );
    };
    let Some(expected) = source.terms.get(start..).filter(|terms| !terms.is_empty()) else {
        return result(
            "unchecked",
            0,
            "no OEIS terms at the catalog offset".to_string(),
        );
    };
    let rows = match generate_flattened(entry, expected.len(), MAX_CHECK_ROWS) {
        Ok(rows) => rows,
        Err(error) => return result("mismatch", 0, error),
    };
    let generated = rows.into_iter().flatten().collect::<Vec<_>>();
    let compared = expected.len().min(generated.len());
    if let Some(position) = (0..compared).find(|&index| generated[index] != expected[index]) {
        return result(
            "mismatch",
            position,
            format!(
                "OEIS term a({}) = {} but the catalog generates {}",
                entry.flattened_offset + position as i64,
                expected[position],
                generated[position]
            ),
        );
    }
    let mut detail = format!("{compared} terms match ({})", source.terms_source);
    if compared < expected.len() {
        detail.push_str(&format!(
            "; only {compared} of {} terms generated within {MAX_CHECK_ROWS} rows",
            expected.len()
        ));
    }
    result("ok", compared, detail)
}

pub(crate) fn cmd_oeis_check(args: &[String]) {
    let mut oeis_data: Option<PathBuf> = None;
    let mut catalog: Option<String> = None;
    let mut ids = Vec::new();
    let mut json_output = false;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--oeis-data" => {
                oeis_data = Some(
                    next_option_value(args, &mut index, "--oeis-data")
                        .unwrap_or_else(fail)
                        .into(),
                )
            }
            "--catalog" => {
                catalog = Some(
                    next_option_value(args, &mut index, "--catalog")
                        .unwrap_or_else(fail)
                        .to_string(),
                )
            }
            "--ids-file" => {
                let path = next_option_value(args, &mut index, "--ids-file").unwrap_or_else(fail);
                ids.extend(read_ids_file(path).unwrap_or_else(fail));
            }
            "--json" => json_output = true,
            other if other.starts_with('-') => fail(format!("unknown oeis check option: {other}")),
            other => ids.push(normalize_id(other).unwrap_or_else(fail)),
        }
        index += 1;
    }
    let root = oeis_data
        .or_else(|| std::env::var_os(OEIS_DATA_ENV).map(PathBuf::from))
        .unwrap_or_else(|| {
            fail(format!(
                "oeis check needs OEIS data: pass --oeis-data <dir> or set {OEIS_DATA_ENV}"
            ))
        });
    let entries: Vec<&'static OeisSequenceDefinition> = match &catalog {
        Some(path) => load_catalog(std::path::Path::new(path))
            .unwrap_or_else(fail)
            .iter()
            .map(leak_definition)
            .collect(),
        None => polytool::oeis::catalog().iter().collect(),
    };
    let selected = entries
        .into_iter()
        .filter(|entry| ids.is_empty() || ids.iter().any(|id| id == entry.id))
        .collect::<Vec<_>>();
    for id in &ids {
        if !selected.iter().any(|entry| entry.id == id) {
            fail::<()>(format!("{id} is not in the catalog"));
        }
    }

    let results = selected
        .iter()
        .map(|entry| check_entry(&root, entry))
        .collect::<Vec<_>>();
    let count = |outcome: &str| results.iter().filter(|r| r.outcome == outcome).count();
    let (ok, mismatched, unchecked) = (count("ok"), count("mismatch"), count("unchecked"));
    if json_output {
        let report = json!({
            "schema": "polytool.oeis-check.v1",
            "checked": ok,
            "mismatched": mismatched,
            "unchecked": unchecked,
            "results": results.iter().map(|result| json!({
                "id": result.id,
                "outcome": result.outcome,
                "checked_terms": result.checked_terms,
                "detail": result.detail,
            })).collect::<Vec<_>>(),
        });
        println!("{}", serde_json::to_string_pretty(&report).unwrap());
    } else {
        for result in &results {
            println!("{}\t{}\t{}", result.id, result.outcome, result.detail);
        }
        println!(
            "summary: {ok} ok, {mismatched} mismatched, {unchecked} unchecked ({} entries)",
            results.len()
        );
    }
    if mismatched > 0 {
        std::process::exit(1);
    }
}
