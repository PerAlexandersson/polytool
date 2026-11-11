//! `polytool recurrence` (search, representations, complete spaces) and
//! `polytool recurrence-generate`.

use crate::input::error_text;
use crate::input::{
    next_option_value, parse_u64_option, parse_usize_option, read_polys_rational,
    read_stdin_or_exit,
};
use crate::recurrence_input::read_recurrence_json;
use crate::recurrence_theorems::{
    note_theorem_matches, print_theorem_report_text, search_theorem_report,
};
use num_bigint::BigInt;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::recurrence::*;
use serde_json::{json, Value};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RecurrenceOutputFormat {
    Text,
    Json,
    Python,
    Latex,
    Lean,
}

fn print_representations_text(report: &RecurrenceRepresentations) {
    println!("Exact fits to supplied rows (not all-rank proofs).");
    println!(
        "{:?}; {} configurations inspected.",
        report.status, report.configurations_considered
    );
    if report.representations.is_empty() {
        println!("No representations found within the inspected bounds.");
        return;
    }

    let mut table = vec![[
        "id", "", "lag", "deriv", "deg(n,t)", "lhs(n,t)", "coeffs", "fit", "check", "eqs",
    ]
    .map(String::from)
    .to_vec()];
    for (index, item) in report.representations.iter().enumerate() {
        let c = &item.complexity;
        table.push(vec![
            format!("[{index}]"),
            if report.frontier.contains(&index) {
                "*"
            } else {
                ""
            }
            .into(),
            c.rec_len.to_string(),
            c.diff_deg.to_string(),
            format!("{},{}", c.idx_deg, c.var_deg),
            format!("{},{}", c.denom_idx_deg, c.denom_var_deg),
            c.nonzero_coefficients.to_string(),
            item.fit_polynomials.to_string(),
            item.verification_polynomials.to_string(),
            item.verification_equations.to_string(),
        ]);
    }
    let widths: Vec<_> = (0..table[0].len())
        .map(|column| table.iter().map(|row| row[column].len()).max().unwrap_or(0))
        .collect();
    println!();
    for row in table {
        println!(
            "{}",
            row.iter()
                .zip(&widths)
                .map(|(cell, width)| format!("{cell:>width$}"))
                .collect::<Vec<_>>()
                .join("  ")
        );
    }
    println!("* Nondominated among returned candidates; indices match the JSON menu.");
    println!("deg/lhs: degrees in (n,t) of RHS/forcing coefficients and the left factor.");
    println!("coeffs: nonzero scalar coefficients; fit/check: fitted/held-out rows.");
    println!("eqs: held-out coefficient equations, not independent constraints.");
    println!();
    for (index, item) in report.representations.iter().enumerate() {
        println!("[{index}] {}", item.recurrence);
        println!(
            "    Checked n={}..={}.",
            item.first_checked_index, item.last_checked_index
        );
        if let Some(presentation) = &item.recurrence_json.presentation {
            if !presentation.left_hand_zeros.is_empty() {
                println!(
                    "    Left-hand zeros: [{}].",
                    presentation.left_hand_zeros.join(", ")
                );
            }
            if let Some(start) = presentation.start_index {
                println!("    Defining range: n >= {start}.");
            } else if let Some(reason) = &presentation.reason {
                println!("    No defining range: {reason}.");
            }
        }
        if !item.singular_checked_indices.is_empty() {
            println!(
                "    Cannot generate at checked singular indices: {:?}.",
                item.singular_checked_indices
            );
        }
    }
    println!("Use --json to save the menu; select an entry with --candidate N in recurrence-generate or recurrence-theorems.");
}

pub(crate) fn cmd_recurrence(args: &[String]) {
    if args.iter().any(|x| x == "--space" || x == "--verify-space") {
        cmd_recurrence_space(args);
        return;
    }

    // The command varies the left-hand side by default (`--no-laurent`).
    let mut search = AdaptiveSearchOptions::default().with_laurent(true);
    let mut budget = AdaptiveSearchBudget::unbounded();
    let mut format = RecurrenceOutputFormat::Text;
    let mut representations = false;
    let mut theorem_matches = false;
    let mut first_index = None;
    let mut lean_name = None;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--lean-name" => {
                let name = next_option_value(args, &mut i, "--lean-name")
                    .unwrap_or_else(|error| usage_error(&error));
                validate_lean_recurrence_name(name)
                    .unwrap_or_else(|error| usage_error(&error.to_string()));
                lean_name = Some(name.to_string());
            }
            "--first-index" => {
                first_index = Some(
                    parse_usize_option(args, &mut i, "--first-index")
                        .unwrap_or_else(|error| usage_error(&error)),
                );
            }
            "--representations" => representations = true,
            "--json" => {
                format = RecurrenceOutputFormat::Json;
            }
            "--python" => {
                format = RecurrenceOutputFormat::Python;
            }
            "--latex" => format = RecurrenceOutputFormat::Latex,
            "--format" => {
                let value = match next_option_value(args, &mut i, "--format") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
                format = match value {
                    "text" => RecurrenceOutputFormat::Text,
                    "json" => RecurrenceOutputFormat::Json,
                    "python" => RecurrenceOutputFormat::Python,
                    "latex" => RecurrenceOutputFormat::Latex,
                    "lean" => RecurrenceOutputFormat::Lean,
                    other => {
                        eprintln!("unknown recurrence output format: {other}");
                        std::process::exit(2);
                    }
                };
            }
            "--skip-prefix" => {
                search.skip_prefix = match parse_usize_option(args, &mut i, "--skip-prefix") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-recurrence-length" => {
                search.min_rec_len =
                    match parse_usize_option(args, &mut i, "--min-recurrence-length") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-recurrence-length" => {
                search.max_rec_len =
                    match parse_usize_option(args, &mut i, "--max-recurrence-length") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--min-variable-degree" => {
                search.min_var_deg = match parse_usize_option(args, &mut i, "--min-variable-degree")
                {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-variable-degree" => {
                search.max_var_deg = match parse_usize_option(args, &mut i, "--max-variable-degree")
                {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-index-degree" => {
                search.min_idx_deg = match parse_usize_option(args, &mut i, "--min-index-degree") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-index-degree" => {
                search.max_idx_deg = match parse_usize_option(args, &mut i, "--max-index-degree") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-derivative-order" => {
                search.min_diff_deg =
                    match parse_usize_option(args, &mut i, "--min-derivative-order") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-derivative-order" => {
                search.max_diff_deg =
                    match parse_usize_option(args, &mut i, "--max-derivative-order") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--try-inhomogeneous" => {
                search.try_inhomogeneous = true;
            }
            "--min-inhomogeneous-variable-degree" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomogeneous-variable-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-inhomogeneous-variable-degree" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomogeneous-variable-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--min-inhomogeneous-index-degree" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomogeneous-index-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-inhomogeneous-index-degree" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomogeneous-index-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--try-denominator" => {
                search.try_denominator = true;
            }
            "--try-alternating-sign" => {
                search.try_alternating_sign = true;
            }
            "--no-laurent" => {
                search.try_laurent = false;
            }
            "--max-laurent-power" => {
                search.max_laurent_x_power =
                    match parse_usize_option(args, &mut i, "--max-laurent-power") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-denominator-variable-degree" => {
                search.try_denominator = true;
                search.max_denom_var_deg =
                    match parse_usize_option(args, &mut i, "--max-denominator-variable-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-denominator-index-degree" => {
                search.try_denominator = true;
                search.max_denom_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-denominator-index-degree") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--min-equation-margin" => {
                search.min_margin = match parse_usize_option(args, &mut i, "--min-equation-margin")
                {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--fit-extra-rows" => {
                search.fit_extra_rows = match parse_usize_option(args, &mut i, "--fit-extra-rows") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--no-held-out-rows" => {
                search.no_verify = true;
            }
            "--no-modular-prefilter" => {
                search.modular_prefilter = false;
            }
            "--modular-prefilter" => {
                search.modular_prefilter = true;
            }
            "--max-candidates" => {
                budget.max_candidates = match parse_usize_option(args, &mut i, "--max-candidates") {
                    Ok(value) => Some(value),
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-work" => {
                budget.max_work = match parse_u64_option(args, &mut i, "--max-work") {
                    Ok(value) => Some(value),
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--verbose" => {
                search.verbose = true;
            }
            "--theorem-matches" => theorem_matches = true,
            other => {
                eprintln!("unknown recurrence option: {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }

    if representations && budget.max_work.is_some() {
        usage_error(
            "--max-work is not supported with --representations; bound it with --max-candidates",
        );
    }
    if representations && theorem_matches {
        usage_error("--representations cannot be combined with --theorem-matches; save the menu with --json, then run `polytool recurrence-theorems --recurrence menu.json --candidate N` (0-based)");
    }
    if representations && matches!(format, RecurrenceOutputFormat::Python) {
        usage_error("--representations supports text, JSON or LaTeX; use a candidate's recurrence_json for row generation");
    }

    let polys = read_polys_rational();
    if polys.len() < 3 {
        eprintln!("Need at least 3 polynomials for recurrence search.");
        std::process::exit(2);
    }

    let output_first_index = first_index
        .map_or(Some(1), |n: usize| n.checked_add(search.skip_prefix))
        .unwrap_or_else(|| usage_error("first index plus skipped rows overflows"));
    if output_first_index > i64::MAX as usize
        || output_first_index.checked_add(polys.len()).is_none()
    {
        usage_error("first index or final row index is too large");
    }
    if representations {
        match find_recurrence_representations_rational(&polys, &search, budget) {
            Ok(report) => {
                let mut report = report
                    .reindexed(output_first_index)
                    .unwrap_or_else(|error| usage_error(&error.to_string()));
                for (index, item) in report.representations.iter_mut().enumerate() {
                    item.recurrence_json = item
                        .recurrence_json
                        .clone()
                        .with_checked_presentation(&polys[search.skip_prefix..])
                        .and_then(RecurrenceJson::with_latex)
                        .unwrap_or_else(|error| usage_error(&error.to_string()));
                    if lean_name.is_some() || matches!(format, RecurrenceOutputFormat::Lean) {
                        item.recurrence_json = item.recurrence_json.clone().with_lean(
                            &format!("{}_{index}", lean_name.as_deref().unwrap_or("P")),
                            &polys[search.skip_prefix..],
                        );
                    }
                    let (rule, _, _) = item
                        .recurrence_json
                        .to_recurrence_parts()
                        .unwrap_or_else(|error| usage_error(&error.to_string()));
                    item.recurrence = rule.to_string();
                    item.latex = rule
                        .normalized_presentation()
                        .unwrap_or_else(|error| usage_error(&error.to_string()))
                        .to_latex();
                }
                if matches!(format, RecurrenceOutputFormat::Json) {
                    println!("{}", serde_json::to_string_pretty(&report).unwrap());
                } else if matches!(format, RecurrenceOutputFormat::Latex) {
                    println!(
                        "% {:?}; {} configurations inspected.",
                        report.status, report.configurations_considered
                    );
                    for (i, item) in report.representations.iter().enumerate() {
                        println!("% Representation {i}");
                        println!(
                            "{}",
                            item.recurrence_json
                                .latex
                                .as_deref()
                                .expect("LaTeX populated")
                        );
                    }
                } else if matches!(format, RecurrenceOutputFormat::Lean) {
                    let mut rendered = 0;
                    for (index, item) in report.representations.iter().enumerate() {
                        if let Some(lean) = &item.recurrence_json.lean {
                            println!("{lean}\n");
                            rendered += 1;
                        } else {
                            println!(
                                "-- Representation {index}: {}",
                                item.recurrence_json
                                    .lean_reason
                                    .as_deref()
                                    .unwrap_or("Lean output unavailable")
                            );
                        }
                    }
                    if rendered == 0 {
                        eprintln!("no representation can be rendered as Lean");
                        std::process::exit(2);
                    }
                } else {
                    print_representations_text(&report);
                }
                if report.status == RepresentationSearchStatus::BudgetExhausted {
                    if format == RecurrenceOutputFormat::Lean {
                        println!(
                            "-- Representation search budget exhausted; this menu is incomplete."
                        );
                    } else {
                        std::process::exit(3);
                    }
                }
            }
            Err(error) => usage_error(&error.to_string()),
        }
        return;
    }
    eprintln!(
        "Searching for recurrence among {} polynomials...",
        polys.len()
    );
    match find_recurrence_adaptive_rational_with_budget(&polys, &search, budget) {
        AdaptiveSearchOutcome::Found(res) => {
            let searched_polys = polys.get(search.skip_prefix..).unwrap_or(&[]);
            let rule = res.recurrence.shift_index(1 - output_first_index as i64);
            let mut theorem_report =
                search_theorem_report(&res, &polys, &search, output_first_index);
            if budget.max_work.is_some() {
                if let Some(metadata) = theorem_report.recurrence_json.search.as_mut() {
                    metadata.work_units = Some(res.diagnostics.work_units);
                }
            }
            let mut recurrence_json = theorem_report
                .recurrence_json
                .clone()
                .with_checked_presentation(searched_polys)
                .and_then(RecurrenceJson::with_latex)
                .unwrap_or_else(|error| usage_error(&error.to_string()));
            if matches!(format, RecurrenceOutputFormat::Lean)
                || (matches!(format, RecurrenceOutputFormat::Json) && lean_name.is_some())
            {
                recurrence_json =
                    recurrence_json.with_lean(lean_name.as_deref().unwrap_or("P"), searched_polys);
            }
            theorem_report.recurrence_json = recurrence_json.clone();
            match format {
                RecurrenceOutputFormat::Json if theorem_matches => {
                    println!("{}", serde_json::to_string_pretty(&theorem_report).unwrap());
                }
                RecurrenceOutputFormat::Json => {
                    println!(
                        "{}",
                        serde_json::to_string_pretty(&recurrence_json).unwrap()
                    );
                }
                RecurrenceOutputFormat::Latex => {
                    println!(
                        "{}",
                        recurrence_json.latex.as_deref().expect("LaTeX populated")
                    );
                }
                RecurrenceOutputFormat::Lean => {
                    let lean = recurrence_json.lean.as_deref().unwrap_or_else(|| {
                        usage_error(
                            recurrence_json
                                .lean_reason
                                .as_deref()
                                .unwrap_or("Lean output unavailable"),
                        )
                    });
                    println!("{lean}");
                }
                RecurrenceOutputFormat::Python => {
                    println!(
                        "{}",
                        rule.to_python_definition_rational_indexed(
                            &searched_polys[..recurrence_json.initial_polynomials.len()],
                            output_first_index as i64
                        )
                    );
                }
                RecurrenceOutputFormat::Text => {
                    println!("{rule}");
                    let presentation = recurrence_json
                        .presentation
                        .as_ref()
                        .expect("presentation populated");
                    if !presentation.left_hand_zeros.is_empty() {
                        println!(
                            "Left-hand zeros (nonnegative integers): [{}]",
                            presentation.left_hand_zeros.join(", ")
                        );
                    }
                    if let Some(start) = presentation.start_index {
                        println!("Range: n >= {start}");
                    } else if let Some(reason) = &presentation.reason {
                        println!("No defining range: {reason}");
                    }
                    for (offset, row) in recurrence_json.initial_polynomials.iter().enumerate() {
                        println!("P({}) = [{}]", output_first_index + offset, row.join(", "));
                    }
                    if theorem_matches {
                        print_theorem_report_text(&theorem_report);
                    }
                }
            }
            if !theorem_matches
                || matches!(
                    format,
                    RecurrenceOutputFormat::Python
                        | RecurrenceOutputFormat::Latex
                        | RecurrenceOutputFormat::Lean
                )
            {
                note_theorem_matches(&theorem_report);
            }
            eprintln!(
                "Found with {} unknowns, weighted score {}, {} equations, {} fit rows, \
                 {} verification rows ({} candidates tried)",
                res.num_unknowns,
                res.weighted_unknowns,
                res.num_equations,
                res.fit_polynomials,
                res.verification_polynomials,
                res.candidates_tried
            );
        }
        AdaptiveSearchOutcome::NoRecurrence(summary) => {
            if format == RecurrenceOutputFormat::Json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "status": "not_found",
                        "candidates_considered": summary.diagnostics.considered_candidates,
                        "candidates_tried": summary.candidates_tried,
                        "max_candidates": budget.max_candidates,
                        "work_units": summary.diagnostics.work_units,
                        "max_work": budget.max_work,
                    }))
                    .unwrap()
                );
            }
            eprintln!("No recurrence found within the search bounds.");
        }
        AdaptiveSearchOutcome::BudgetExhausted(summary) => {
            if format == RecurrenceOutputFormat::Json {
                println!(
                    "{}",
                    serde_json::to_string_pretty(&json!({
                        "status": "budget_exhausted",
                        "candidates_considered": summary.diagnostics.considered_candidates,
                        "candidates_tried": summary.candidates_tried,
                        "max_candidates": budget.max_candidates,
                        "work_units": summary.diagnostics.work_units,
                        "max_work": budget.max_work,
                    }))
                    .unwrap()
                );
            }
            eprintln!(
                "Recurrence search exhausted its budget after {} candidates and {} work units; \
                 unsearched candidates remain.",
                summary.diagnostics.considered_candidates, summary.diagnostics.work_units
            );
            std::process::exit(3);
        }
    }
}

fn cmd_recurrence_space(args: &[String]) {
    let verify = args.iter().any(|x| x == "--verify-space");
    let search = args.iter().any(|x| x == "--space");
    if verify == search
        || args
            .iter()
            .any(|x| !["--space", "--verify-space", "--json"].contains(&x.as_str()))
    {
        eprintln!("choose --space or --verify-space, with optional --json; bounds belong in the JSON request");
        std::process::exit(2);
    }
    let input = read_stdin_or_exit();
    let outcome: Result<(Value, bool), String> = if verify {
        serde_json::from_str::<RecurrenceSpaces>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_recurrence_spaces(&report.request, &report).map_err(error_text)?;
                let incomplete = report
                    .results
                    .iter()
                    .filter(|r| r.certificate.is_none())
                    .count();
                Ok((
                    json!({"verified_completed_claims":true,"incomplete_shapes":incomplete}),
                    incomplete > 0,
                ))
            })
    } else {
        serde_json::from_str::<RecurrenceSpaceRequest>(&input)
            .map_err(error_text)
            .and_then(|request| {
                let report = find_recurrence_spaces(&request).map_err(error_text)?;
                let incomplete = report.results.iter().any(|r| r.certificate.is_none());
                Ok((
                    serde_json::to_value(report).map_err(error_text)?,
                    incomplete,
                ))
            })
    };
    match outcome {
        Ok((value, incomplete)) => {
            println!("{}", serde_json::to_string_pretty(&value).unwrap());
            if incomplete {
                std::process::exit(3);
            }
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

pub(crate) fn format_rational_row(coeffs: &[RecurrenceBigRational]) -> String {
    coeffs
        .iter()
        .map(format_rational_coeff)
        .collect::<Vec<_>>()
        .join(", ")
}

fn usage_error(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(2);
}

/// A rational coefficient row as integers, or `None` if some coefficient is
/// not an integer.
fn integer_row(row: &[RecurrenceBigRational]) -> Option<Vec<BigInt>> {
    row.iter()
        .map(|value| value.is_integer().then(|| value.to_integer()))
        .collect()
}

/// Options of `recurrence-generate --format bfile`.
struct BFileRequest {
    options: polytool::oeis::BFileOptions,
    given: bool,
}

pub(crate) fn cmd_recurrence_generate(args: &[String]) {
    fn usage_error(message: &str) -> ! {
        eprintln!("{message}");
        std::process::exit(2);
    }

    let mut recurrence_path: Option<String> = None;
    let mut candidate = None;
    let mut first_index_override = None;
    let mut rows: Option<usize> = None;
    let mut additional: Option<usize> = None;
    let mut format = "text";
    let mut bfile = BFileRequest {
        options: polytool::oeis::BFileOptions::default(),
        given: false,
    };
    // Row reading shared by the b-file and the OEIS text formats.
    let mut reading_given = false;
    let mut offset_given = false;
    let mut first_row: Option<usize> = None;
    let mut regular_triangle = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--candidate" => {
                candidate = Some(
                    parse_usize_option(args, &mut i, "--candidate")
                        .unwrap_or_else(|e| usage_error(&e)),
                );
            }
            "--recurrence" => {
                i += 1;
                recurrence_path = args.get(i).cloned();
            }
            "--first-index" => {
                first_index_override = Some(
                    parse_usize_option(args, &mut i, "--first-index")
                        .unwrap_or_else(|error| usage_error(&error)),
                );
            }
            "--count" => {
                rows = Some(
                    parse_usize_option(args, &mut i, "--count").unwrap_or_else(|e| usage_error(&e)),
                );
            }
            "--additional-count" => {
                additional = Some(
                    parse_usize_option(args, &mut i, "--additional-count")
                        .unwrap_or_else(|e| usage_error(&e)),
                );
            }
            "--json" => format = "json",
            "--format" => {
                let value =
                    next_option_value(args, &mut i, "--format").unwrap_or_else(|e| usage_error(&e));
                format = match value {
                    "text" | "rows" => "text",
                    "json" => "json",
                    "bfile" | "b-file" => "bfile",
                    "oeis-example" => "oeis-example",
                    "oeis-formula" => "oeis-formula",
                    "mathematica" => "mathematica",
                    "pari" => "pari",
                    "oeis-json" => "oeis-json",
                    other => usage_error(&format!(
                        "unknown recurrence-generate format '{other}' (expected text, json, bfile, \
                         oeis-example, oeis-formula, mathematica, pari or oeis-json)"
                    )),
                };
            }
            "--first-term-index" => {
                let value = next_option_value(args, &mut i, "--first-term-index")
                    .unwrap_or_else(|e| usage_error(&e));
                bfile.options.first_index = value.parse().unwrap_or_else(|_| {
                    usage_error(&format!(
                        "--first-term-index expects an integer, got '{value}'"
                    ))
                });
                reading_given = true;
                offset_given = true;
            }
            "--from-index" => {
                first_row = Some(
                    parse_usize_option(args, &mut i, "--from-index")
                        .unwrap_or_else(|e| usage_error(&e)),
                );
                reading_given = true;
            }
            "--regular-triangle" => {
                regular_triangle = true;
                reading_given = true;
            }
            "--max-terms" => {
                bfile.options.max_terms = Some(
                    parse_usize_option(args, &mut i, "--max-terms")
                        .unwrap_or_else(|e| usage_error(&e)),
                );
                bfile.given = true;
            }
            "--descending" => {
                bfile.options.orientation = polytool::oeis::BFileRowOrientation::Reversed;
                reading_given = true;
            }
            "--nearest-row" => {
                bfile.options.stopping_rule = polytool::oeis::BFileStoppingRule::NearestWholeRow;
                bfile.given = true;
            }
            "--trim-zeros" => {
                let value = next_option_value(args, &mut i, "--trim-zeros")
                    .unwrap_or_else(|e| usage_error(&e));
                bfile.options.zero_policy = match value {
                    "none" => polytool::oeis::BFileZeroPolicy::KeepAll,
                    "low" => polytool::oeis::BFileZeroPolicy::TrimLowOrder,
                    "high" => polytool::oeis::BFileZeroPolicy::TrimHighOrder,
                    "both" => polytool::oeis::BFileZeroPolicy::TrimBothEnds,
                    other => usage_error(&format!(
                        "--trim-zeros expects none, low, high or both, got '{other}'"
                    )),
                };
                reading_given = true;
            }
            other => {
                eprintln!("unknown recurrence-generate option: {other}");
                std::process::exit(2);
            }
        }
        i += 1;
    }

    // Without --recurrence the record is the command input (stdin or
    // --input FILE).
    let path = recurrence_path.unwrap_or_else(|| "-".to_string());
    if rows.is_some() && additional.is_some() {
        eprintln!("use only one of --count or --additional-count");
        std::process::exit(2);
    }
    let oeis_text_format = matches!(
        format,
        "oeis-example" | "oeis-formula" | "mathematica" | "pari" | "oeis-json"
    );
    if format != "bfile" && bfile.given {
        usage_error("--max-terms and --nearest-row need --format bfile");
    }
    if format != "bfile" && !oeis_text_format && reading_given {
        usage_error(
            "--first-term-index, --from-index, --descending, --trim-zeros and --regular-triangle need \
             --format bfile or an OEIS text format",
        );
    }

    let mut recurrence_json = match read_recurrence_json(&path, candidate) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };
    if let Some(first_index) = first_index_override {
        recurrence_json = recurrence_json
            .reindexed(first_index)
            .unwrap_or_else(|error| usage_error(&error.to_string()));
    }
    let (recurrence, first_index, initial_polys) = match recurrence_json.to_recurrence_parts() {
        Ok(parts) => parts,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };

    let total_rows = if let Some(rows) = rows {
        Some(rows)
    } else if let Some(additional) = additional {
        match initial_polys.len().checked_add(additional) {
            Some(total) => Some(total),
            None => {
                eprintln!("recurrence row count overflow");
                std::process::exit(2);
            }
        }
    } else {
        None
    };

    let reading = polytool::oeis::OeisTriangleReading {
        first_row_index: if offset_given {
            bfile.options.first_index
        } else if oeis_text_format {
            i64::try_from(first_row.unwrap_or(first_index))
                .unwrap_or_else(|_| usage_error("--from-index does not fit in i64"))
        } else {
            0
        },
        orientation: bfile.options.orientation,
        zero_policy: bfile.options.zero_policy,
        regular_triangle,
    };

    if oeis_text_format {
        let presentation = polytool::oeis::OeisRecurrencePresentation::from_recurrence_json(
            &recurrence_json,
            first_row,
            reading,
        )
        .unwrap_or_else(|error| usage_error(&error.to_string()));
        let text = presentation
            .oeis_text(total_rows.unwrap_or(10))
            .unwrap_or_else(|error| usage_error(&error.to_string()));
        match format {
            "oeis-example" => print!("{}", text.triangle_example),
            "oeis-formula" => println!("{}", text.recurrence_sentence),
            "mathematica" => println!("{}", text.mathematica),
            "pari" => println!("{}", text.pari),
            _ => println!("{}", serde_json::to_string_pretty(&text).unwrap()),
        }
        return;
    }

    if format == "bfile" && (first_row.is_some() || regular_triangle) {
        // Read the rows from `--from-index` on as OEIS rows, then flatten them.
        let presentation = polytool::oeis::OeisRecurrencePresentation::from_recurrence_json(
            &recurrence_json,
            first_row,
            reading,
        )
        .unwrap_or_else(|error| usage_error(&error.to_string()));
        let mut iterator = presentation.rows_iter();
        let rows = collect_bfile_rows(&mut iterator, total_rows, &bfile.options);
        let read = presentation.reading().read_rows(&rows);
        let options = polytool::oeis::BFileOptions {
            orientation: polytool::oeis::BFileRowOrientation::AsStored,
            zero_policy: polytool::oeis::BFileZeroPolicy::KeepAll,
            ..bfile.options
        };
        print_bfile(&read, &options);
        return;
    }

    if format == "bfile" {
        let Some(max_terms) = bfile.options.max_terms.filter(|_| total_rows.is_none()) else {
            let Some(total_rows) = total_rows else {
                usage_error(
                    "recurrence-generate --format bfile needs --count <n>, --additional-count <n> or --max-terms <n>",
                );
            };
            let generated =
                match recurrence.generate_rows_rational(&initial_polys, first_index, total_rows) {
                    Ok(rows) => rows,
                    Err(error) => usage_error(&error.to_string()),
                };
            let rows = generated
                .iter()
                .enumerate()
                .map(|(index, row)| integer_row(row).unwrap_or_else(|| non_integral_error(index)))
                .collect::<Vec<_>>();
            print_bfile(&rows, &bfile.options);
            return;
        };
        // Only a term cap: pull rows until the cap is reached.
        let mut generated: Vec<Vec<BigInt>> = Vec::new();
        let mut terms = 0usize;
        let mut barren_rows = 0usize;
        let mut iterator = recurrence.generate_rows_rational_iter(&initial_polys, first_index);
        while terms < max_terms {
            let Some(row) = iterator.next() else { break };
            let row = match row {
                Ok(row) => row,
                Err(error) => usage_error(&error.to_string()),
            };
            let row = integer_row(&row).unwrap_or_else(|| non_integral_error(generated.len()));
            let kept = polytool::oeis::bfile_row_terms(&row, bfile.options.zero_policy).len();
            if kept == 0 {
                barren_rows += 1;
                if barren_rows > 10_000 {
                    usage_error("10000 consecutive rows without b-file terms; use --count");
                }
            } else {
                barren_rows = 0;
            }
            terms += kept;
            generated.push(row);
        }
        print_bfile(&generated, &bfile.options);
        return;
    }

    let Some(total_rows) = total_rows else {
        eprintln!("recurrence-generate needs --count <n> or --additional-count <n>");
        std::process::exit(2);
    };
    let generated = match recurrence.generate_rows_rational(&initial_polys, first_index, total_rows)
    {
        Ok(rows) => rows,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };

    if format == "json" {
        let rows_json: Vec<Vec<String>> = generated
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect();
        println!(
            "{}",
            serde_json::to_string_pretty(&serde_json::json!({
                "first_index": first_index,
                "polynomials": rows_json,
            }))
            .unwrap()
        );
    } else {
        for row in &generated {
            println!("{}", format_rational_row(row));
        }
    }
}

/// Integer rows for a b-file from `iterator`: `total_rows` rows, or without it
/// rows until the term cap of `options` is reached.  Exits with status 2 on
/// failure.
fn collect_bfile_rows(
    iterator: &mut RecurrenceRowIterator,
    total_rows: Option<usize>,
    options: &polytool::oeis::BFileOptions,
) -> Vec<Vec<BigInt>> {
    let mut generated: Vec<Vec<BigInt>> = Vec::new();
    let max_terms = match (total_rows, options.max_terms) {
        (Some(_), _) => usize::MAX,
        (None, Some(max_terms)) => max_terms,
        (None, None) => {
            eprintln!(
                "recurrence-generate --format bfile needs --count <n>, --additional-count <n> or --max-terms <n>"
            );
            std::process::exit(2);
        }
    };
    let mut terms = 0usize;
    let mut barren_rows = 0usize;
    while terms < max_terms && total_rows.is_none_or(|total| generated.len() < total) {
        let Some(row) = iterator.next() else { break };
        let row = row.unwrap_or_else(|error| {
            eprintln!("{error}");
            std::process::exit(2);
        });
        let row = integer_row(&row).unwrap_or_else(|| non_integral_error(generated.len()));
        let kept = polytool::oeis::bfile_row_terms(&row, options.zero_policy).len();
        if kept == 0 {
            barren_rows += 1;
            if barren_rows > 10_000 {
                eprintln!("10000 consecutive rows without b-file terms; use --count");
                std::process::exit(2);
            }
        } else {
            barren_rows = 0;
        }
        terms += kept;
        generated.push(row);
    }
    generated
}

fn non_integral_error(index: usize) -> ! {
    eprintln!("row {index} has a non-integral coefficient; a b-file needs integers");
    std::process::exit(2);
}

/// Write the b-file of `rows`; exits with status 2 on failure.
fn print_bfile(rows: &[Vec<BigInt>], options: &polytool::oeis::BFileOptions) {
    let file = polytool::oeis::flatten_rows_to_bfile(rows, options).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2);
    });
    if options.max_terms.is_some()
        && file.rows_included == 0
        && rows.iter().any(|row| !row.is_empty())
    {
        eprintln!("--max-terms is smaller than the first complete row (try --nearest-row)");
        std::process::exit(2);
    }
    print!("{}", file.text);
}
