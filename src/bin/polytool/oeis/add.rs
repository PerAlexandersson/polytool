//! `polytool oeis add`: find a recurrence for OEIS rows and record it, with
//! provenance, in the catalog data file.

use super::catalog_file::{
    catalog_recurrence, leak_definition, load_catalog, resolve_catalog_path, write_catalog,
    DEFAULT_CATALOG_PATH,
};
use super::sha256::sha256_hex;
use super::source::{parse_bfile, read_local_record, OeisSourceRecord, RowWidths};
use crate::input::{next_option_value, parse_usize_option};
use num_bigint::BigInt;
use polytool::oeis::{
    OeisCatalogProvenance, OeisCatalogRecord, OeisLayout, OeisProvenanceSource, OeisSequenceStatus,
    OeisValidationRow,
};
use polytool::recurrence::{
    find_recurrence_adaptive_rational_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome, BigRational,
};
use serde_json::{json, Value};
use std::collections::BTreeMap;
use std::path::PathBuf;

/// Environment variable naming a local OEIS internal-format directory.
pub(crate) const OEIS_DATA_ENV: &str = "OEIS_DATA_DIR";

/// Fewest complete rows `oeis add` accepts.
const MIN_ROWS: usize = 5;

pub(crate) fn print_oeis_add_help() {
    println!("Usage:");
    println!("  polytool oeis add <A-number>... [options]");
    println!("  polytool oeis add --ids-file <file> [options]");
    println!(
        "  polytool oeis add <A-number> --bfile <file> --name <title> --layout <layout> [options]"
    );
    println!();
    println!("Read OEIS terms, cut them into complete rows, run the adaptive recurrence");
    println!("search (the same search and options as `polytool recurrence`) on a prefix");
    println!("of the rows, and check the recurrence exactly against every remaining");
    println!("term through the catalog's own row generator.  A passing entry is written");
    println!("to the catalog data file with provenance `polytool_recurrence_search`;");
    println!("its status is `verified` when held-out rows were checked, else");
    println!("`validated`.  Rebuild polytool to bundle new entries.");
    println!();
    println!("Input:");
    println!(
        "  --oeis-data <dir>       Local OEIS internal-format copy (default: ${OEIS_DATA_ENV});"
    );
    println!("                          reads <dir>/seq/A123/A123456.seq (or <dir>/A123/...),");
    println!("                          and a real b-file <dir>/files/A123/b123456.txt when it");
    println!("                          extends the %S/%T/%U terms");
    println!("  --bfile <file>          Explicit b-file (one A-number only); needs --name and");
    println!("                          --layout; the offset defaults to the first b-file index");
    println!("  --name <title>          Sequence title (with --bfile)");
    println!("  --offset <n>            First flattened index and displayed row (with --bfile)");
    println!("  --ids-file <file>       A-numbers, one per line (# comments allowed)");
    println!();
    println!("Rows:");
    println!("  OEIS `tabl` sequences are read as regular triangles (rows of length");
    println!("  1, 2, 3, ...).  Other sequences need an explicit rule:");
    println!("  --layout triangle       Treat the terms as a regular triangle");
    println!("  --layout table          Irregular table; requires --row-widths");
    println!("  --row-widths <W+S[/R]>  Row i has W + S*floor(i/R) terms (R defaults to 1):");
    println!("                          1+2 gives 1, 3, 5, ...; 1+1/2 gives 1, 1, 2, 2, ...");
    println!("  --fit-rows <n>          Rows given to the search (default: all but");
    println!("                          max(2, rows/4), which are held out)");
    println!();
    println!("Catalog:");
    println!("  --catalog <file>        Data file (default: ./{DEFAULT_CATALOG_PATH} if present)");
    println!("  --dry-run               Do not write; print each new entry as a JSONL line");
    println!("  --replace               Search again and update entries already present");
    println!("  --json                  Emit a polytool.oeis-add.v1 report");
    println!();
    println!("Search options (as in `polytool recurrence`): --skip-prefix, --min/max-rec-len,");
    println!("  --min/max-var-deg, --min/max-idx-deg, --min/max-diff-deg, --inhomogeneous,");
    println!("  --min/max-inhomo-var-deg, --min/max-inhomo-idx-deg, --denominator,");
    println!("  --max-denom-var-deg, --max-denom-idx-deg, --alternating-sign, --min-margin,");
    println!("  --fit-extra-rows, --no-verify, --no-modular-prefilter, --max-candidates <n>,");
    println!("  --verbose.");
    println!();
    println!("Outcomes per A-number: added, updated, already_present, no_recurrence");
    println!("(nothing within the bounds or budget), layout_unsupported, insufficient_data,");
    println!("mismatch (a held-out term differs), unavailable (no record).  The exit");
    println!("status is 0 when every outcome is added, updated or already_present.");
    println!();
    println!("Examples:");
    println!("  polytool oeis add A008292 --oeis-data ~/oeis --dry-run");
    println!("  polytool oeis add --ids-file ids.txt --max-candidates 20000");
}

/// Parse one search option shared with `polytool recurrence`; returns
/// `Ok(false)` when `args[*index]` is not a search option.
pub(crate) fn parse_search_option(
    args: &[String],
    index: &mut usize,
    search: &mut AdaptiveSearchOptions,
    budget: &mut AdaptiveSearchBudget,
) -> Result<bool, String> {
    let option = args[*index].as_str();
    let value = |index: &mut usize| parse_usize_option(args, index, option);
    match option {
        "--skip-prefix" => search.skip_prefix = value(index)?,
        "--min-rec-len" => search.min_rec_len = value(index)?,
        "--max-rec-len" => search.max_rec_len = value(index)?,
        "--min-var-deg" => search.min_var_deg = value(index)?,
        "--max-var-deg" => search.max_var_deg = value(index)?,
        "--min-idx-deg" => search.min_idx_deg = value(index)?,
        "--max-idx-deg" => search.max_idx_deg = value(index)?,
        "--min-diff-deg" => search.min_diff_deg = value(index)?,
        "--max-diff-deg" => search.max_diff_deg = value(index)?,
        "--inhomogeneous" => search.try_inhomogeneous = true,
        "--min-inhomo-var-deg"
        | "--max-inhomo-var-deg"
        | "--min-inhomo-idx-deg"
        | "--max-inhomo-idx-deg" => {
            search.try_inhomogeneous = true;
            let parsed = value(index)?;
            match option {
                "--min-inhomo-var-deg" => search.min_inhomo_var_deg = parsed,
                "--max-inhomo-var-deg" => search.max_inhomo_var_deg = parsed,
                "--min-inhomo-idx-deg" => search.min_inhomo_idx_deg = parsed,
                _ => search.max_inhomo_idx_deg = parsed,
            }
        }
        "--denominator" | "--try-denominator" => search.try_denominator = true,
        "--max-denom-var-deg" => {
            search.try_denominator = true;
            search.max_denom_var_deg = value(index)?;
        }
        "--max-denom-idx-deg" => {
            search.try_denominator = true;
            search.max_denom_idx_deg = value(index)?;
        }
        "--alternating-sign" => search.try_alternating_sign = true,
        "--min-margin" => search.min_margin = value(index)?,
        "--fit-extra-rows" => search.fit_extra_rows = value(index)?,
        "--no-verify" => search.no_verify = true,
        "--no-modular-prefilter" => search.modular_prefilter = false,
        "--modular-prefilter" => search.modular_prefilter = true,
        "--max-candidates" => budget.max_candidates = Some(value(index)?),
        "--verbose" => search.verbose = true,
        _ => return Ok(false),
    }
    Ok(true)
}

/// Normalize and check an A-number.
pub(crate) fn normalize_id(text: &str) -> Result<String, String> {
    let id = text.trim().to_ascii_uppercase();
    let bytes = id.as_bytes();
    if bytes.len() == 7 && bytes[0] == b'A' && bytes[1..].iter().all(u8::is_ascii_digit) {
        Ok(id)
    } else {
        Err(format!(
            "`{text}` is not an A-number (expected A followed by six digits)"
        ))
    }
}

/// A-numbers listed in a file, one per line; text after the first
/// whitespace and `#` comments are ignored.
pub(crate) fn read_ids_file(path: &str) -> Result<Vec<String>, String> {
    let text =
        std::fs::read_to_string(path).map_err(|error| format!("failed to read {path}: {error}"))?;
    text.lines()
        .map(|line| line.split('#').next().unwrap_or("").trim())
        .filter(|line| !line.is_empty())
        .map(|line| normalize_id(line.split_whitespace().next().unwrap_or(line)))
        .collect()
}

#[derive(Default)]
struct AddOptions {
    ids: Vec<String>,
    oeis_data: Option<PathBuf>,
    bfile: Option<String>,
    name: Option<String>,
    offset: Option<i64>,
    layout: Option<OeisLayout>,
    row_widths: Option<RowWidths>,
    fit_rows: Option<usize>,
    catalog: Option<String>,
    dry_run: bool,
    replace: bool,
    json: bool,
    search: AdaptiveSearchOptions,
    budget: AdaptiveSearchBudget,
}

fn parse_add_options(args: &[String]) -> Result<AddOptions, String> {
    let mut options = AddOptions::default();
    let mut index = 0;
    while index < args.len() {
        if parse_search_option(args, &mut index, &mut options.search, &mut options.budget)? {
            index += 1;
            continue;
        }
        match args[index].as_str() {
            "--oeis-data" => {
                options.oeis_data = Some(next_option_value(args, &mut index, "--oeis-data")?.into())
            }
            "--bfile" => {
                options.bfile = Some(next_option_value(args, &mut index, "--bfile")?.to_string())
            }
            "--name" => {
                options.name = Some(next_option_value(args, &mut index, "--name")?.to_string())
            }
            "--offset" => {
                let value = next_option_value(args, &mut index, "--offset")?;
                options.offset = Some(
                    value
                        .parse()
                        .map_err(|_| format!("--offset expects an integer, got '{value}'"))?,
                );
            }
            "--layout" => {
                options.layout = Some(match next_option_value(args, &mut index, "--layout")? {
                    "triangle" | "regular_triangle" | "tabl" => OeisLayout::RegularTriangle,
                    "table" | "tabf" => OeisLayout::Table,
                    other => return Err(format!("unknown layout `{other}` (triangle or table)")),
                })
            }
            "--row-widths" => {
                options.row_widths = Some(RowWidths::parse(next_option_value(
                    args,
                    &mut index,
                    "--row-widths",
                )?)?)
            }
            "--fit-rows" => {
                options.fit_rows = Some(parse_usize_option(args, &mut index, "--fit-rows")?)
            }
            "--catalog" => {
                options.catalog =
                    Some(next_option_value(args, &mut index, "--catalog")?.to_string())
            }
            "--ids-file" => {
                let path = next_option_value(args, &mut index, "--ids-file")?;
                options.ids.extend(read_ids_file(path)?);
            }
            "--dry-run" => options.dry_run = true,
            "--replace" => options.replace = true,
            "--json" => options.json = true,
            other if other.starts_with('-') => {
                return Err(format!("unknown oeis add option: {other}"))
            }
            other => options.ids.push(normalize_id(other)?),
        }
        index += 1;
    }
    if options.ids.is_empty() {
        return Err("oeis add needs at least one A-number (or --ids-file)".to_string());
    }
    let mut seen = std::collections::BTreeSet::new();
    options.ids.retain(|id| seen.insert(id.clone()));
    if options.bfile.is_some() {
        if options.ids.len() != 1 {
            return Err("--bfile applies to exactly one A-number".to_string());
        }
        if options.name.is_none() || options.layout.is_none() {
            return Err("--bfile needs --name and --layout".to_string());
        }
    } else {
        if options.name.is_some() || options.offset.is_some() {
            return Err("--name and --offset are only used with --bfile".to_string());
        }
        if options.oeis_data.is_none() {
            options.oeis_data = std::env::var_os(OEIS_DATA_ENV).map(PathBuf::from);
        }
        if options.oeis_data.is_none() {
            return Err(format!(
                "oeis add needs OEIS data: pass --oeis-data <dir>, set {OEIS_DATA_ENV}, or use --bfile"
            ));
        }
    }
    if options.layout == Some(OeisLayout::Table) && options.row_widths.is_none() {
        return Err("--layout table needs --row-widths W+S[/R]".to_string());
    }
    if options.layout == Some(OeisLayout::RegularTriangle)
        && options
            .row_widths
            .is_some_and(|widths| widths != RowWidths::TRIANGLE)
    {
        return Err("--layout triangle means --row-widths 1+1".to_string());
    }
    Ok(options)
}

/// Result of processing one A-number.
struct AddResult {
    id: String,
    outcome: &'static str,
    detail: String,
    record: Option<OeisCatalogRecord>,
}

impl AddResult {
    fn new(id: &str, outcome: &'static str, detail: impl Into<String>) -> Self {
        Self {
            id: id.to_string(),
            outcome,
            detail: detail.into(),
            record: None,
        }
    }

    fn succeeded(&self) -> bool {
        matches!(self.outcome, "added" | "updated" | "already_present")
    }
}

fn source_record(options: &AddOptions, id: &str) -> Result<OeisSourceRecord, String> {
    if let Some(path) = &options.bfile {
        let text = std::fs::read_to_string(path)
            .map_err(|error| format!("failed to read {path}: {error}"))?;
        let (first, terms) = parse_bfile(&text)?;
        let offset = options.offset.unwrap_or(first);
        if offset != first {
            return Err(format!(
                "--offset {offset} differs from the first b-file index {first}"
            ));
        }
        return Ok(OeisSourceRecord {
            name: options.name.clone().unwrap_or_default(),
            offset,
            keywords: Vec::new(),
            terms,
            terms_source: "bfile",
        });
    }
    let root = options
        .oeis_data
        .as_ref()
        .expect("checked by option parsing");
    read_local_record(root, id)
}

fn search_options_json(search: &AdaptiveSearchOptions, budget: AdaptiveSearchBudget) -> Value {
    let value = serde_json::to_value(search).expect("search options serialize");
    let mut entries = value
        .as_object()
        .map(|object| {
            object
                .iter()
                .filter(|(key, _)| key.as_str() != "verbose")
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect::<BTreeMap<_, _>>()
        })
        .unwrap_or_default();
    entries.insert("max_candidates".to_string(), json!(budget.max_candidates));
    // Sorted keys keep data-file lines canonical whatever map order serde_json uses.
    Value::Object(entries.into_iter().collect())
}

fn string_rows(rows: &[&[BigInt]]) -> Vec<Vec<String>> {
    rows.iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect()
}

/// Search, validate against every term, and build the record for one id.
fn add_one(options: &AddOptions, id: &str, present: bool) -> AddResult {
    if present && !options.replace {
        return AddResult::new(id, "already_present", "pass --replace to search again");
    }
    let source = match source_record(options, id) {
        Ok(source) => source,
        Err(error) => return AddResult::new(id, "unavailable", error),
    };
    let widths = match options.row_widths {
        Some(widths) => widths,
        None if options.layout == Some(OeisLayout::RegularTriangle)
            || source.has_keyword("tabl") =>
        {
            RowWidths::TRIANGLE
        }
        None => {
            return AddResult::new(
                id,
                "layout_unsupported",
                format!(
                    "keywords `{}` lack `tabl`, so the row structure is unknown; pass \
                     --layout triangle or --layout table --row-widths W+S[/R]",
                    source.keywords.join(",")
                ),
            )
        }
    };
    let layout = options.layout.unwrap_or(if widths == RowWidths::TRIANGLE {
        OeisLayout::RegularTriangle
    } else {
        OeisLayout::Table
    });
    let (rows, leftover) = widths.split(&source.terms);
    let skip = options.search.skip_prefix;
    if rows.len() < MIN_ROWS.max(skip + 3) {
        return AddResult::new(
            id,
            "insufficient_data",
            format!(
                "only {} complete rows from {} terms",
                rows.len(),
                source.terms.len()
            ),
        );
    }
    let default_fit = rows.len() - (rows.len() / 4).max(2);
    let fit_rows = options.fit_rows.unwrap_or(default_fit).min(rows.len());
    if fit_rows < skip + 3 {
        return AddResult::new(
            id,
            "insufficient_data",
            format!("--fit-rows {fit_rows} leaves fewer than 3 rows after --skip-prefix"),
        );
    }
    let fit = rows[..fit_rows]
        .iter()
        .map(|row| {
            row.iter()
                .map(|value| BigRational::from_integer(value.clone()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let found = match find_recurrence_adaptive_rational_with_budget(
        &fit,
        &options.search,
        options.budget,
    ) {
        AdaptiveSearchOutcome::Found(found) => found,
        AdaptiveSearchOutcome::NoRecurrence(summary) => {
            return AddResult::new(
                id,
                "no_recurrence",
                format!(
                    "no recurrence within the search bounds ({} candidates considered)",
                    summary.diagnostics.considered_candidates
                ),
            )
        }
        AdaptiveSearchOutcome::BudgetExhausted(summary) => {
            return AddResult::new(
                id,
                "no_recurrence",
                format!(
                    "candidate budget exhausted after {} candidates",
                    summary.diagnostics.considered_candidates
                ),
            )
        }
    };
    let searched = &rows[skip..fit_rows];
    let initial_count = found.recurrence.generation_initial_count(1, searched.len());
    let verification_rows = rows.len() - fit_rows;
    let status = if verification_rows > 0 {
        OeisSequenceStatus::Verified
    } else {
        OeisSequenceStatus::Validated
    };
    let checked_terms = source.terms.len();
    let terms_text = source
        .terms
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let last_row = rows.len() - 1;
    let record = OeisCatalogRecord {
        id: id.to_string(),
        name: source.name.clone(),
        layout,
        first_row: source.offset,
        flattened_offset: source.offset,
        bfile_prefix_verified: true,
        prefix_rows: string_rows(&rows[..skip]),
        recurrence_first_index: 1,
        recurrence_first_width: rows[skip].len(),
        initial_rows: string_rows(&searched[..initial_count]),
        recurrence: catalog_recurrence(&found.recurrence),
        provenance: OeisCatalogProvenance {
            source: OeisProvenanceSource::PolytoolRecurrenceSearch,
            status,
            source_rows: rows.len(),
            verification_rows,
            rows_sha256: sha256_hex(terms_text.as_bytes()),
            fixture_slug: None,
            fixture_row_offset: None,
            validation_row: Some(OeisValidationRow {
                index: last_row,
                terms: rows[last_row].iter().map(ToString::to_string).collect(),
            }),
            fit_rows: Some(fit_rows),
            checked_terms: Some(checked_terms),
            terms_source: Some(source.terms_source.to_string()),
            search_options: Some(search_options_json(&options.search, options.budget)),
            polytool_version: Some(env!("CARGO_PKG_VERSION").to_string()),
        },
    };
    if let Err(error) = record.validate() {
        return AddResult::new(id, "mismatch", format!("invalid catalog record: {error}"));
    }

    // Replay through the catalog's own generator and compare every term.
    let definition = leak_definition(&record);
    let wanted_rows = rows.len() + usize::from(leftover > 0);
    let generated = match definition.generate_rows(wanted_rows) {
        Ok(generated) => generated,
        Err(error) => return AddResult::new(id, "mismatch", error.to_string()),
    };
    for (index, row) in rows.iter().enumerate() {
        if generated[index].as_slice() != *row {
            let kind = if index < fit_rows {
                "fitted"
            } else {
                "held-out"
            };
            return AddResult::new(
                id,
                "mismatch",
                format!(
                    "{kind} row {} (displayed row {}) differs from the OEIS terms",
                    index,
                    source.offset + index as i64
                ),
            );
        }
    }
    if leftover > 0 {
        let partial = &source.terms[source.terms.len() - leftover..];
        if !generated[rows.len()].starts_with(partial) {
            return AddResult::new(
                id,
                "mismatch",
                "the trailing partial row differs from the OEIS terms",
            );
        }
    }
    let outcome = if present { "updated" } else { "added" };
    let mut result = AddResult::new(
        id,
        outcome,
        format!(
            "{} [{} fit rows, {} held-out rows, {} terms checked, {}]",
            found.recurrence,
            fit_rows,
            verification_rows,
            checked_terms,
            status.as_str()
        ),
    );
    result.record = Some(record);
    result
}

pub(crate) fn cmd_oeis_add(args: &[String]) {
    if args.iter().any(|arg| arg == "--help" || arg == "-h") {
        print_oeis_add_help();
        return;
    }
    let options = parse_add_options(args).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2);
    });
    let catalog_path = resolve_catalog_path(options.catalog.as_deref());
    if catalog_path.is_none() && !options.dry_run {
        eprintln!(
            "no catalog data file: run from the polytool source directory (./{DEFAULT_CATALOG_PATH}), \
             pass --catalog <file>, or use --dry-run"
        );
        std::process::exit(2);
    }
    let mut records: BTreeMap<String, OeisCatalogRecord> = BTreeMap::new();
    if let Some(path) = &catalog_path {
        for record in load_catalog(path).unwrap_or_else(|error| {
            eprintln!("{error}");
            std::process::exit(2);
        }) {
            records.insert(record.id.clone(), record);
        }
    }
    let is_present = |id: &str| match &catalog_path {
        Some(_) => records.contains_key(id),
        None => polytool::oeis::by_id(id).is_some(),
    };

    let mut results = Vec::new();
    for id in &options.ids {
        let result = add_one(&options, id, is_present(id));
        if !options.json {
            println!("{}\t{}\t{}", result.id, result.outcome, result.detail);
            if options.dry_run {
                if let Some(record) = &result.record {
                    println!("{}", record.to_jsonl_line());
                }
            }
        }
        results.push(result);
    }

    let changed = results
        .iter()
        .filter(|result| result.record.is_some())
        .count();
    if !options.dry_run && changed > 0 {
        let path = catalog_path.as_ref().expect("checked above");
        for result in &results {
            if let Some(record) = &result.record {
                records.insert(record.id.clone(), record.clone());
            }
        }
        let mut all = records.into_values().collect::<Vec<_>>();
        if let Err(error) = write_catalog(path, &mut all) {
            eprintln!("{error}");
            std::process::exit(2);
        }
        if !options.json {
            eprintln!(
                "wrote {changed} entr{} to {}; rebuild polytool to bundle them",
                if changed == 1 { "y" } else { "ies" },
                path.display()
            );
        }
    }
    if options.json {
        let report = json!({
            "schema": "polytool.oeis-add.v1",
            "dry_run": options.dry_run,
            "catalog": catalog_path.as_ref().map(|path| path.display().to_string()),
            "results": results.iter().map(|result| json!({
                "id": result.id,
                "outcome": result.outcome,
                "detail": result.detail,
                "record": result.record,
            })).collect::<Vec<_>>(),
        });
        println!("{}", serde_json::to_string_pretty(&report).unwrap());
    }
    if !results.iter().all(AddResult::succeeded) {
        std::process::exit(1);
    }
}
