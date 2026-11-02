//! `polytool oeis evidence`: compute, show and record the finite evidence
//! (strongest coefficient property on the first rows) of catalog entries.

use super::add::{normalize_id, read_ids_file};
use super::catalog_file::{leak_definition, load_catalog, resolve_catalog_path, write_catalog};
use crate::input::{next_option_value, parse_usize_option};
use polytool::oeis::{
    OeisCatalogEvidence, OeisCatalogRecord, OeisEvidenceAccumulator, OeisEvidenceFirstFailures,
    OeisEvidenceProperty, OeisRealRootednessStatus, OeisSequenceDefinition,
};
use serde_json::{json, Value};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

/// Default per-entry time limit in seconds.
const DEFAULT_ENTRY_TIMEOUT_SECONDS: u64 = 60;

pub(crate) fn print_oeis_evidence_help() {
    println!("Usage:");
    println!("  polytool oeis evidence [A-number...|--all] --rows <n> [options]");
    println!();
    println!("Compute, for each entry, the strongest coefficient property that holds on");
    println!("every one of its first <n> rows (fewer if the entry cannot generate that");
    println!("many), exactly, from the catalog's own rows.  Properties, strongest first:");
    println!("real_rooted, ultra_log_concave, log_concave, unimodal, none.  real_rooted");
    println!("is the exact polytool check (any coefficient signs); the other three need a");
    println!("row with nonnegative coefficients and no zero between two nonzero ones, and");
    println!("then use polytool's own log-concavity, ultra-log-concavity and unimodality");
    println!("checks (ultra-log-concavity at the row's degree).  The record lists the");
    println!("first failing row (displayed OEIS row index) of each stronger property.");
    println!();
    println!("This is evidence about finitely many rows, not a proof, and it is separate");
    println!("from the real_rootedness flag.");
    println!();
    println!("Options:");
    println!("  --rows <n>             Number of rows to check (required, at least 1)");
    println!("  --all                  Every catalog entry (instead of listing A-numbers)");
    println!("  --ids-file <file>      Restrict to the listed A-numbers");
    println!("  --write                Store the records in the catalog data file; only the");
    println!("                         `evidence` key of each computed entry changes");
    println!("  --per-entry-timeout-seconds <s>");
    println!("                         Limit per entry, checked before each row (default");
    println!("                         {DEFAULT_ENTRY_TIMEOUT_SECONDS}; 0 = none).  An entry over the limit gets no record");
    println!("                         (an existing record is kept) and is listed");
    println!("  --as-of <YYYY-MM-DD>   Date to record (default: today, UTC)");
    println!("  --catalog <file>       Data file to read (and to write with --write);");
    println!("                         default ./data/oeis-catalog.jsonl if present, else");
    println!("                         the bundled catalog (then --write is not possible)");
    println!("  --json                 Emit a polytool.oeis-evidence.v1 report");
    println!();
    println!("Besides the per-entry lines the report gives the distribution of `strongest`");
    println!("and lists every lean_proved entry whose evidence is weaker than real_rooted");
    println!("(a contradiction to investigate).");
}

fn fail<T>(error: String) -> T {
    eprintln!("{error}");
    std::process::exit(2);
}

/// Today's date (UTC) as `YYYY-MM-DD`.
fn today() -> String {
    let seconds = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |elapsed| elapsed.as_secs());
    let days = (seconds / 86_400) as i64;
    // Civil date from days since 1970-01-01 (Hinnant's algorithm).
    let shifted = days + 719_468;
    let era = shifted.div_euclid(146_097);
    let day_of_era = shifted.rem_euclid(146_097);
    let year_of_era =
        (day_of_era - day_of_era / 1_460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_index = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_index + 2) / 5 + 1;
    let month = if month_index < 10 {
        month_index + 3
    } else {
        month_index - 9
    };
    let year = year_of_era + era * 400 + i64::from(month <= 2);
    format!("{year:04}-{month:02}-{day:02}")
}

enum Outcome {
    /// A record covering `rows_checked` rows; `note` explains a short record.
    Computed {
        evidence: OeisCatalogEvidence,
        note: Option<String>,
    },
    /// The per-entry limit was exceeded (after this many rows).
    Timeout { rows_done: usize },
    /// Not even the first row could be generated.
    Error(String),
}

fn compute_entry(
    entry: &'static OeisSequenceDefinition,
    rows: usize,
    limit: Option<Duration>,
    as_of: &str,
) -> Outcome {
    let start = Instant::now();
    let mut iterator = match entry.generate_rows_iter() {
        Ok(iterator) => iterator,
        Err(error) => return Outcome::Error(error.to_string()),
    };
    let mut accumulator = OeisEvidenceAccumulator::new(entry.first_row);
    let mut note = None;
    while accumulator.rows_checked() < rows {
        if limit.is_some_and(|limit| start.elapsed() > limit) {
            return Outcome::Timeout {
                rows_done: accumulator.rows_checked(),
            };
        }
        match iterator.next() {
            Some(Ok(row)) => accumulator.push_row(&row),
            Some(Err(error)) if accumulator.rows_checked() == 0 => {
                return Outcome::Error(error.to_string())
            }
            Some(Err(error)) => {
                note = Some(format!("row generation stopped: {error}"));
                break;
            }
            None => {
                note = Some("no further rows".to_string());
                break;
            }
        }
    }
    if limit.is_some_and(|limit| start.elapsed() > limit) {
        return Outcome::Timeout {
            rows_done: accumulator.rows_checked(),
        };
    }
    match accumulator.finish(as_of) {
        Some(evidence) => Outcome::Computed { evidence, note },
        None => Outcome::Error("no rows requested".to_string()),
    }
}

/// JSON object of an evidence record (`first_failure` lists only failing properties).
pub(super) fn evidence_json(
    strongest: OeisEvidenceProperty,
    rows_checked: usize,
    first_failure: &OeisEvidenceFirstFailures,
    as_of: &str,
) -> Value {
    let mut failures = serde_json::Map::new();
    for property in OeisEvidenceProperty::ALL {
        if let Some(row) = first_failure.get(property) {
            failures.insert(property.as_str().to_string(), json!(row));
        }
    }
    json!({
        "strongest": strongest.as_str(),
        "rows_checked": rows_checked,
        "first_failure": failures,
        "as_of": as_of,
    })
}

/// The `oeis info` line describing an evidence record.
pub(super) fn describe_evidence(
    strongest: OeisEvidenceProperty,
    rows_checked: usize,
    first_failure: &OeisEvidenceFirstFailures,
    as_of: &str,
) -> String {
    let failures = OeisEvidenceProperty::ALL
        .iter()
        .filter_map(|&property| {
            first_failure
                .get(property)
                .map(|row| format!("{} fails first at row {row}", property.as_str()))
        })
        .collect::<Vec<_>>();
    format!(
        "evidence: strongest property on the first {rows_checked} rows: {} (as of {as_of}; \
         finite check, not a proof{}{})",
        strongest.as_str(),
        if failures.is_empty() { "" } else { "; " },
        failures.join("; ")
    )
}

/// The canonical line of `record` with its evidence removed.
fn line_without_evidence(record: &OeisCatalogRecord) -> String {
    let mut record = record.clone();
    record.evidence = None;
    record.to_jsonl_line()
}

pub(crate) fn cmd_oeis_evidence(args: &[String]) {
    if args.iter().any(|arg| arg == "--help" || arg == "-h") {
        print_oeis_evidence_help();
        return;
    }
    let mut ids: Vec<String> = Vec::new();
    let mut all = false;
    let mut rows: Option<usize> = None;
    let mut write = false;
    let mut json_output = false;
    let mut timeout_seconds = DEFAULT_ENTRY_TIMEOUT_SECONDS;
    let mut as_of: Option<String> = None;
    let mut catalog: Option<String> = None;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--all" => all = true,
            "--write" => write = true,
            "--json" => json_output = true,
            "--rows" => {
                rows = Some(parse_usize_option(args, &mut index, "--rows").unwrap_or_else(fail))
            }
            "--per-entry-timeout-seconds" => {
                timeout_seconds =
                    parse_usize_option(args, &mut index, "--per-entry-timeout-seconds")
                        .unwrap_or_else(fail) as u64
            }
            "--as-of" => {
                as_of = Some(
                    next_option_value(args, &mut index, "--as-of")
                        .unwrap_or_else(fail)
                        .to_string(),
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
            other if other.starts_with('-') => {
                fail(format!("unknown oeis evidence option: {other}"))
            }
            other => ids.push(normalize_id(other).unwrap_or_else(fail)),
        }
        index += 1;
    }
    let Some(rows) = rows.filter(|&rows| rows > 0) else {
        fail::<()>("oeis evidence needs --rows <n> with n at least 1".to_string());
        return;
    };
    if all == !ids.is_empty() {
        fail::<()>("oeis evidence needs A-numbers (or --ids-file), or --all, not both".to_string());
    }
    let as_of = as_of.unwrap_or_else(today);
    let limit = (timeout_seconds > 0).then(|| Duration::from_secs(timeout_seconds));

    // The records to work on: the data file when there is one, else the bundled catalog.
    let path = resolve_catalog_path(catalog.as_deref());
    let mut records = path
        .as_deref()
        .map(|path| load_catalog(path).unwrap_or_else(fail));
    if write && records.is_none() {
        fail::<()>("--write needs the catalog data file (--catalog <file>)".to_string());
    }
    let original_lines: Vec<String> = records
        .iter()
        .flatten()
        .map(OeisCatalogRecord::to_jsonl_line)
        .collect();
    if let (Some(path), true) = (&path, write) {
        let text = std::fs::read_to_string(path).unwrap_or_else(|error| fail(error.to_string()));
        let on_disk = text.lines().filter(|line| !line.trim().is_empty());
        if on_disk.eq(original_lines.iter().map(String::as_str)) && text.ends_with('\n') {
        } else {
            fail::<()>(format!(
                "{} is not in canonical form; refusing to rewrite it",
                path.display()
            ));
        }
    }
    let definitions: Vec<(
        &'static OeisSequenceDefinition,
        Option<OeisRealRootednessStatus>,
    )> = match &records {
        Some(records) => records
            .iter()
            .map(|record| {
                (
                    leak_definition(record),
                    record.real_rootedness.as_ref().map(|flag| flag.status),
                )
            })
            .collect(),
        None => polytool::oeis::catalog()
            .iter()
            .map(|entry| (entry, entry.real_rootedness().map(|flag| flag.status)))
            .collect(),
    };
    let selected: Vec<_> = if all {
        definitions.iter().collect()
    } else {
        ids.iter()
            .map(|id| {
                definitions
                    .iter()
                    .find(|(entry, _)| entry.id == id)
                    .unwrap_or_else(|| fail(format!("{id} is not in the catalog")))
            })
            .collect()
    };

    let mut distribution = [0usize; 5];
    let mut report = Vec::new();
    let mut unrecorded = Vec::new();
    let mut contradictions = Vec::new();
    let mut new_evidence: Vec<(&str, OeisCatalogEvidence)> = Vec::new();
    for &(entry, flag) in &selected {
        match compute_entry(entry, rows, limit, &as_of) {
            Outcome::Computed { evidence, note } => {
                distribution[evidence.strongest.rank()] += 1;
                if *flag == Some(OeisRealRootednessStatus::LeanProved)
                    && evidence.strongest != OeisEvidenceProperty::RealRooted
                {
                    contradictions.push(entry.id);
                }
                let failures = OeisEvidenceProperty::ALL
                    .iter()
                    .filter_map(|&property| {
                        evidence
                            .first_failure
                            .get(property)
                            .map(|row| format!("{}@{row}", property.as_str()))
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                if !json_output {
                    println!(
                        "{}\t{}\trows={}\t{}\t{}",
                        entry.id,
                        evidence.strongest.as_str(),
                        evidence.rows_checked,
                        if failures.is_empty() { "-" } else { &failures },
                        note.as_deref().unwrap_or("")
                    );
                }
                report.push(json!({
                    "id": entry.id,
                    "outcome": "computed",
                    "evidence": evidence_json(
                        evidence.strongest,
                        evidence.rows_checked,
                        &evidence.first_failure,
                        &evidence.as_of,
                    ),
                    "note": note,
                }));
                new_evidence.push((entry.id, evidence));
            }
            Outcome::Timeout { rows_done } => {
                unrecorded.push(entry.id);
                if !json_output {
                    println!("{}\ttimeout\tafter {rows_done} rows", entry.id);
                }
                report.push(json!({"id": entry.id, "outcome": "timeout", "rows_done": rows_done}));
            }
            Outcome::Error(message) => {
                unrecorded.push(entry.id);
                if !json_output {
                    println!("{}\terror\t{message}", entry.id);
                }
                report.push(json!({"id": entry.id, "outcome": "error", "message": message}));
            }
        }
    }

    let mut written = false;
    if let (Some(records), Some(path), true) = (records.as_mut(), path.as_deref(), write) {
        for (id, evidence) in &new_evidence {
            let record = records
                .iter_mut()
                .find(|record| record.id == *id)
                .expect("selected entries come from the records");
            record.evidence = Some(evidence.clone());
        }
        // Only the evidence key may change: every line, with its evidence removed, is
        // the original line with its evidence removed; unchanged records keep their bytes.
        let new_lines: Vec<String> = records
            .iter()
            .map(OeisCatalogRecord::to_jsonl_line)
            .collect();
        let originals = load_catalog(path).unwrap_or_else(fail);
        assert_eq!(originals.len(), records.len());
        for ((original, record), (old_line, new_line)) in originals
            .iter()
            .zip(records.iter())
            .zip(original_lines.iter().zip(&new_lines))
        {
            assert_eq!(
                line_without_evidence(original),
                line_without_evidence(record),
                "{}",
                record.id
            );
            if original.evidence == record.evidence {
                assert_eq!(old_line, new_line, "{}", record.id);
            }
        }
        write_catalog(path, records).unwrap_or_else(fail);
        written = true;
    }

    let names = [
        "real_rooted",
        "ultra_log_concave",
        "log_concave",
        "unimodal",
        "none",
    ];
    if json_output {
        let distribution_json: serde_json::Map<String, Value> = names
            .iter()
            .zip(distribution)
            .map(|(name, count)| (name.to_string(), json!(count)))
            .collect();
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "schema": "polytool.oeis-evidence.v1",
                "rows_requested": rows,
                "as_of": as_of,
                "written": written,
                "entries": report,
                "distribution": distribution_json,
                "unrecorded": unrecorded,
                "lean_proved_with_weaker_evidence": contradictions,
            }))
            .unwrap()
        );
    } else {
        println!();
        for (name, count) in names.iter().zip(distribution) {
            println!("strongest {name}: {count}");
        }
        println!("without a record (timeout or error): {}", unrecorded.len());
        println!(
            "lean_proved with evidence weaker than real_rooted: {}",
            if contradictions.is_empty() {
                "none".to_string()
            } else {
                contradictions.join(", ")
            }
        );
        if written {
            println!("wrote {} record(s)", new_evidence.len());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evidence_is_described_and_serialized_with_its_failures() {
        let failures = OeisEvidenceFirstFailures {
            real_rooted: Some(2),
            ultra_log_concave: Some(2),
            ..Default::default()
        };
        assert_eq!(
            describe_evidence(
                OeisEvidenceProperty::LogConcave,
                50,
                &failures,
                "2026-11-02"
            ),
            "evidence: strongest property on the first 50 rows: log_concave (as of 2026-11-02; \
             finite check, not a proof; real_rooted fails first at row 2; \
             ultra_log_concave fails first at row 2)"
        );
        assert_eq!(
            describe_evidence(
                OeisEvidenceProperty::RealRooted,
                7,
                &OeisEvidenceFirstFailures::default(),
                "2026-11-02"
            ),
            "evidence: strongest property on the first 7 rows: real_rooted (as of 2026-11-02; \
             finite check, not a proof)"
        );
        assert_eq!(
            evidence_json(
                OeisEvidenceProperty::LogConcave,
                50,
                &failures,
                "2026-11-02"
            ),
            json!({
                "strongest": "log_concave",
                "rows_checked": 50,
                "first_failure": {"real_rooted": 2, "ultra_log_concave": 2},
                "as_of": "2026-11-02",
            })
        );
    }
}
