//! `polytool recurrence-theorems` and the theorem-match output of
//! `polytool recurrence --theorem-matches`.

use crate::input::error_text;
use crate::input::{is_help_arg, read_limited_text};
use polytool::recurrence::{
    AdaptiveSearchOptions, AdaptiveSearchResult, BigRational, RecurrenceJson, RecurrenceJsonSearch,
    RecurrenceOptionsJson,
};
use polytool::renewal::RecurrenceTheoremReport;
use std::fs;
use std::io;

pub(crate) fn print_recurrence_theorems_help() {
    println!("Usage:");
    println!("  polytool recurrence-theorems --recurrence <file|-> [--json]");
    println!("  polytool recurrence-theorems --help");
    println!();
    println!("Check a supplied recurrence JSON record (polytool.recurrence.v1, e.g. from");
    println!("`polytool recurrence --json`) against the supported interlacing theorems and");
    println!("report explicit conditional theorem matches.");
    println!();
    println!("Theorems (paper proofs, not Lean checked;");
    println!("https://github.com/PerAlexandersson/RealRooted/issues/2037):");
    println!("  weighted_renewal_window          P_n = sum_j r_j t^j P_(n-1-ell*j), r_0 = 1,");
    println!("                                   r_j >= 0, R(z) = sum r_j z^j real-rooted;");
    println!("                                   P_n << P_(n+s) for d*s <= d*ell+1, d = deg R.");
    println!(
        "  fixed_delay_predecessor_window   P_n = P_(n-1) + w(n) t P_(n-L), w polynomial in n"
    );
    println!("                                   with w(n) >= 0 for all integers n >= 1 (exact");
    println!("                                   certificate, not sampling); Q_A << Q_B for");
    println!("                                   B - A <= L with Q_(max(0, n-L)) boundary.");
    println!();
    println!("Hypotheses are checked exactly. The record's initial rows are compared with the");
    println!("canonical boundary rows (finite check; `canonical_after_relabel` when the first");
    println!("supplied row plays the role of index 0). A match is conditional: it never shows");
    println!("that the recurrence holds for the intended sequence for all n.");
    println!();
    println!("Options:");
    println!("  --recurrence <file|->   Read recurrence JSON from a file, or stdin with -");
    println!("  --json                  Emit the polytool.recurrence-theorems.v1 report");
    println!("  -h, --help              Print this help text");
    println!();
    println!("Exit status: 0 with at least one match, 1 with none, 2 on invalid input.");
    println!();
    println!("Example:");
    println!(
        "  polytool recurrence --json < rows.txt | polytool recurrence-theorems --recurrence -"
    );
}

/// Theorem report for an adaptive-search result; `polys` are all input rows.
/// Its `recurrence_json` is the record that `recurrence --json` prints.
pub(crate) fn search_theorem_report(
    result: &AdaptiveSearchResult,
    polys: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
) -> RecurrenceTheoremReport {
    let searched = polys.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result
        .recurrence
        .generation_initial_count(1, searched.len());
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
            candidates_tried: result.candidates_tried,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    );
    RecurrenceTheoremReport::for_search(result, recurrence_json, searched)
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
}

pub(crate) fn cmd_recurrence_theorems(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_theorems_help();
        return;
    }
    let mut path = None;
    let mut json = false;
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
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
    let Some(path) = path else {
        eprintln!("recurrence-theorems needs --recurrence <file|->");
        std::process::exit(2);
    };
    let input = if path == "-" {
        read_limited_text(io::stdin(), "recurrence JSON from stdin")
    } else {
        fs::File::open(&path)
            .map_err(|e| format!("failed to open recurrence JSON `{path}`: {e}"))
            .and_then(|file| read_limited_text(file, &format!("recurrence JSON `{path}`")))
    };
    let report = input
        .and_then(|text| {
            serde_json::from_str::<RecurrenceJson>(&text)
                .map_err(|e| format!("failed to parse recurrence JSON: {e}"))
        })
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
