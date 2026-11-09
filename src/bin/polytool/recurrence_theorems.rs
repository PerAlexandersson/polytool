//! `polytool recurrence-theorems` and the theorem-match output of
//! `polytool recurrence --theorem-matches`.

use crate::input::error_text;
use crate::input::parse_usize_option;
use crate::recurrence_input::read_recurrence_json;
use polytool::recurrence::{
    AdaptiveSearchOptions, AdaptiveSearchResult, BigRational, RecurrenceJson, RecurrenceJsonSearch,
    RecurrenceOptionsJson,
};
use polytool::renewal::RecurrenceTheoremReport;

/// Theorem report for an adaptive-search result; `polys` are all input rows.
/// Its `recurrence_json` is the record that `recurrence --json` prints.
pub(crate) fn search_theorem_report(
    result: &AdaptiveSearchResult,
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    first_index: usize,
) -> RecurrenceTheoremReport {
    let searched = polys.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result.export_initial_count(searched.len());
    let recurrence_json = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        &searched[..initial_count],
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows: polys.len(),
            skip_prefix: search.skip_prefix,
            unknowns: result.num_unknowns,
            weighted_unknowns: result.weighted_unknowns,
            equations: result.num_equations,
            fit_polynomials: result.fit_polynomials,
            verification_polynomials: result.verification_polynomials,
            verification_equations: result.verification_equations,
            candidates_tried: result.candidates_tried,
            work_units: None,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    )
    .reindexed(first_index)
    .expect("CLI index shift was validated");
    let mut shifted_result = result.clone();
    shifted_result.recurrence = result.recurrence.shift_index(1 - first_index as i64);
    RecurrenceTheoremReport::for_search(&shifted_result, recurrence_json, searched)
}

/// Print one stderr note per conditional theorem match (no-op without matches).
pub(crate) fn note_theorem_matches(report: &RecurrenceTheoremReport) {
    for theorem in &report.theorem_matches {
        eprintln!("Note: {}", theorem.summary());
    }
}

/// Text rendering of a theorem report.
pub(crate) fn print_theorem_report_text(report: &RecurrenceTheoremReport) {
    println!(
        "Finite evidence ({}): {}",
        report.finite_evidence.source, report.finite_evidence.note
    );
    if report.theorem_matches.is_empty() {
        println!("No conditional theorem match.");
    }
    for theorem in &report.theorem_matches {
        println!("Conditional theorem match: {}", theorem.theorem);
        println!("  family: {}", theorem.family);
        for hypothesis in &theorem.hypotheses {
            println!(
                "  hypothesis {} [{}]: {}",
                hypothesis.name, hypothesis.check, hypothesis.statement
            );
        }
        println!("  conclusion: {}", theorem.conclusion);
        println!(
            "  window: {} ({})",
            theorem.window.inequality, theorem.window.statement
        );
        println!(
            "  boundary: {:?}, {} supplied rows from index {} ({} required); contract: {}",
            theorem.boundary.status,
            theorem.boundary.supplied_rows,
            theorem.boundary.supplied_first_index,
            theorem.boundary.required_rows,
            theorem.boundary.contract
        );
        for obligation in &theorem.remaining_obligations {
            println!("  remaining: {obligation}");
        }
        println!("  status: {} ({})", theorem.proof_status, theorem.reference);
    }
    for rejection in &report.rejected_theorems {
        println!("Not matched {}: {}", rejection.theorem, rejection.reason);
    }
    if let Some(trivial) = report.trivially_real_rooted {
        println!(
            "Trivially real-rooted (product seed class, product proof): {}",
            if trivial { "yes" } else { "no" }
        );
    }
}

pub(crate) fn cmd_recurrence_theorems(args: &[String]) {
    let mut path = None;
    let mut candidate = None;
    let mut json = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--candidate" => {
                candidate = Some(
                    parse_usize_option(args, &mut i, "--candidate").unwrap_or_else(|error| {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }),
                );
            }
            "--recurrence" => {
                i += 1;
                path = args.get(i).cloned();
            }
            "--json" => json = true,
            other => {
                eprintln!("unknown recurrence-theorems option: {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }
    // Without --recurrence the record is the command input (stdin or
    // --input FILE).
    let path = path.unwrap_or_else(|| "-".to_string());
    let report = read_recurrence_json(&path, candidate)
        .and_then(|record| RecurrenceTheoremReport::for_supplied(&record).map_err(error_text));
    let report = match report {
        Ok(report) => report,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };
    if json {
        println!("{}", serde_json::to_string_pretty(&report).unwrap());
    } else {
        print_theorem_report_text(&report);
    }
    if report.theorem_matches.is_empty() {
        std::process::exit(1);
    }
}
