//! `polytool recurrence` (search, representations, complete spaces) and
//! `polytool recurrence-generate`.

use crate::help::print_rational_coefficient_input_help;
use crate::input::error_text;
use crate::input::{
    is_help_arg, next_option_value, parse_usize_option, read_limited_text, read_polys_rational,
    read_stdin_or_exit, OutputFormat,
};
use crate::recurrence_theorems::{
    note_theorem_matches, print_theorem_report_text, search_theorem_report,
};
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::recurrence::*;
use serde_json::{json, Value};
use std::fs;
use std::io::{self};

pub(crate) fn print_recurrence_help() {
    println!("Usage:");
    println!("  polytool recurrence [options]");
    println!("  polytool recurrence --help");
    println!();
    println!("Search adaptively for a linear differential-polynomial recurrence");
    println!("among the input polynomial sequence.");
    println!();
    print_rational_coefficient_input_help();
    println!();
    println!("Search options:");
    println!("  --skip-prefix <n>          Ignore the first n input polynomials");
    println!("  --min-rec-len <n>          Minimum number of previous rows to use");
    println!("  --max-rec-len <n>          Maximum number of previous rows to use");
    println!("  --min-var-deg <d>          Minimum t-degree for recurrence coefficients");
    println!("  --max-var-deg <d>          Maximum t-degree for recurrence coefficients");
    println!("  --min-idx-deg <d>          Minimum n-degree for recurrence coefficients");
    println!("  --max-idx-deg <d>          Maximum n-degree for recurrence coefficients");
    println!("  --min-diff-deg <d>         Minimum derivative order");
    println!("  --max-diff-deg <d>         Maximum derivative order");
    println!("  --inhomogeneous            Try inhomogeneous terms");
    println!("  --min-inhomo-var-deg <d>   Minimum t-degree for inhomogeneous terms");
    println!("  --max-inhomo-var-deg <d>   Maximum t-degree for inhomogeneous terms");
    println!("  --min-inhomo-idx-deg <d>   Minimum n-degree for inhomogeneous terms");
    println!("  --max-inhomo-idx-deg <d>   Maximum n-degree for inhomogeneous terms");
    println!("  --denominator              Try LHS denominators");
    println!("  --try-denominator          Alias for --denominator");
    println!("  --max-denom-var-deg <d>    Maximum t-degree for denominators");
    println!("  --max-denom-idx-deg <d>    Maximum n-degree for denominators");
    println!("  --alternating-sign         Also try alternating-sign terms");
    println!("  --min-margin <n>           Require equations >= unknowns + n");
    println!("  --fit-extra-rows <n>       Extra rows beyond the first solvable prefix");
    println!("  --no-verify                Use all input rows for fitting");
    println!("  --no-modular-prefilter     Disable default modular candidate rejection");
    println!("  --max-candidates <n>       Inspect at most n candidate configurations");
    println!("  --representations          Return a menu of length/derivative tradeoffs");
    println!("  --space                    Complete bounded spaces; read a recurrence-space JSON request from stdin");
    println!(
        "  --verify-space             Replay a recurrence-space JSON report from stdin (no search)"
    );
    println!("                             These two modes emit exact JSON and accept only optional --json.");
    println!("                             Source-backed positives are conditional on the supplied source/seeds.");
    println!("                             Rows-only positives are fits; dual negatives prove absence in shape.");
    println!("                             Opt-in \"coefficient_rescaling\": {{\"kind\":\"factorial_power\",\"power\":m}}");
    println!("                             searches Q[n] = T(P[n]), T(X^k) = X^k/(k!)^m, m >= 1, and certifies");
    println!("                             complete spaces from r(U+m(A+J)+1) residuals (constant-index source).");
    println!("                             With --json: structured candidates and coverage");
    println!("  --json                     Emit recurrence JSON with initial conditions");
    println!("  --format json              Alias for --json");
    println!("  --python                   Emit exact Python code for the recurrence");
    println!("  --format python            Alias for --python");
    println!("  --verbose                  Print candidate search details");
    println!("  --theorem-matches          Also report conditional theorem matches");
    println!("                             (weighted renewal, fixed-delay predecessor),");
    println!("                             separate from the finite fit evidence. With --json,");
    println!("                             emit polytool.recurrence-theorems.v1 (the usual");
    println!("                             record is its recurrence_json). Without the flag a");
    println!("                             match is only noted on stderr. Not used with");
    println!("                             --representations; see recurrence-theorems.");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("Example:");
    println!("  printf '1\\n1\\n2\\n3\\n5\\n8\\n' | polytool recurrence");
    println!("  polytool recurrence --json < rows.txt > recurrence.json");
}

pub(crate) fn print_recurrence_generate_help() {
    println!("Usage:");
    println!("  polytool recurrence-generate --recurrence <file|-> --rows <n>");
    println!("  polytool recurrence-generate --recurrence <file|-> --additional <n>");
    println!("  polytool recurrence-generate --help");
    println!();
    println!("Generate coefficient rows from recurrence JSON emitted by");
    println!("`polytool recurrence --json`.");
    println!();
    println!("Options:");
    println!("  --recurrence <file|->      Read recurrence JSON from file, or stdin with -");
    println!("  --rows <n>                 Output n total rows, including initial rows");
    println!("  --additional <n>           Output initial rows plus n generated rows");
    println!("  --json                     Emit generated rows as JSON");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("Example:");
    println!("  polytool recurrence-generate --recurrence recurrence.json --rows 50");
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RecurrenceOutputFormat {
    Text,
    Json,
    Python,
}

pub(crate) fn cmd_recurrence(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_help();
        return;
    }

    if args.iter().any(|x| x == "--space" || x == "--verify-space") {
        cmd_recurrence_space(args);
        return;
    }

    let mut search = AdaptiveSearchOptions::default();
    let mut budget = AdaptiveSearchBudget::unbounded();
    let mut format = RecurrenceOutputFormat::Text;
    let mut representations = false;
    let mut theorem_matches = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--representations" => representations = true,
            "--json" => {
                format = RecurrenceOutputFormat::Json;
            }
            "--python" => {
                format = RecurrenceOutputFormat::Python;
            }
            "--format" => {
                let value = match next_option_value(args, &mut i, "--format") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
                format = match value {
                    "text" => RecurrenceOutputFormat::Text,
                    "json" => RecurrenceOutputFormat::Json,
                    "python" => RecurrenceOutputFormat::Python,
                    other => {
                        eprintln!("unknown recurrence output format: {other}");
                        return;
                    }
                };
            }
            "--skip-prefix" => {
                search.skip_prefix = match parse_usize_option(args, &mut i, "--skip-prefix") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-rec-len" => {
                search.min_rec_len = match parse_usize_option(args, &mut i, "--min-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-rec-len" => {
                search.max_rec_len = match parse_usize_option(args, &mut i, "--max-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-var-deg" => {
                search.min_var_deg = match parse_usize_option(args, &mut i, "--min-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-var-deg" => {
                search.max_var_deg = match parse_usize_option(args, &mut i, "--max-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-idx-deg" => {
                search.min_idx_deg = match parse_usize_option(args, &mut i, "--min-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-idx-deg" => {
                search.max_idx_deg = match parse_usize_option(args, &mut i, "--max-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-diff-deg" => {
                search.min_diff_deg = match parse_usize_option(args, &mut i, "--min-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-diff-deg" => {
                search.max_diff_deg = match parse_usize_option(args, &mut i, "--max-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--inhomogeneous" => {
                search.try_inhomogeneous = true;
            }
            "--min-inhomo-var-deg" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomo-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-inhomo-var-deg" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_var_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomo-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--min-inhomo-idx-deg" => {
                search.try_inhomogeneous = true;
                search.min_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--min-inhomo-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-inhomo-idx-deg" => {
                search.try_inhomogeneous = true;
                search.max_inhomo_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-inhomo-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--denominator" | "--try-denominator" => {
                search.try_denominator = true;
            }
            "--alternating-sign" => {
                search.try_alternating_sign = true;
            }
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg =
                    match parse_usize_option(args, &mut i, "--max-denom-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--max-denom-idx-deg" => {
                search.try_denominator = true;
                search.max_denom_idx_deg =
                    match parse_usize_option(args, &mut i, "--max-denom-idx-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            return;
                        }
                    };
            }
            "--min-margin" => {
                search.min_margin = match parse_usize_option(args, &mut i, "--min-margin") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--fit-extra-rows" => {
                search.fit_extra_rows = match parse_usize_option(args, &mut i, "--fit-extra-rows") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--no-verify" => {
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
                        return;
                    }
                };
            }
            "--verbose" => {
                search.verbose = true;
            }
            "--theorem-matches" => theorem_matches = true,
            other => {
                eprintln!("unknown recurrence option: {other}");
                return;
            }
        }
        i += 1;
    }

    let polys = read_polys_rational();
    if polys.len() < 3 {
        eprintln!("Need at least 3 polynomials for recurrence search.");
        return;
    }

    if representations {
        if matches!(format, RecurrenceOutputFormat::Python) {
            eprintln!("--representations supports text or JSON; use a candidate's recurrence_json for row generation");
            return;
        }
        match find_recurrence_representations_rational(&polys, &search, budget) {
            Ok(report) => {
                if matches!(format, RecurrenceOutputFormat::Json) {
                    println!("{}", serde_json::to_string_pretty(&report).unwrap());
                } else {
                    println!("Exact fits to supplied rows (not all-rank proofs).");
                    println!(
                        "{:?}; {} configurations inspected.",
                        report.status, report.configurations_considered
                    );
                    for (i, item) in report.representations.iter().enumerate() {
                        println!(
                            "{}[{}] length {}, derivative order {}, checked n={}..={}: {}",
                            if report.frontier.contains(&i) {
                                "* "
                            } else {
                                "  "
                            },
                            i,
                            item.complexity.rec_len,
                            item.complexity.diff_deg,
                            item.first_checked_index,
                            item.last_checked_index,
                            item.recurrence
                        );
                    }
                    println!("* Nondominated among returned candidates; use --json for bounds and coverage.");
                }
                if report.status == RepresentationSearchStatus::BudgetExhausted {
                    std::process::exit(3);
                }
            }
            Err(error) => eprintln!("{error}"),
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
            let initial_count = res
                .recurrence
                .generation_initial_count(1, searched_polys.len());
            let initial_polys = &searched_polys[..initial_count];
            let theorem_report = search_theorem_report(&res, &polys, &search);
            match format {
                RecurrenceOutputFormat::Json if theorem_matches => {
                    println!("{}", serde_json::to_string_pretty(&theorem_report).unwrap());
                }
                RecurrenceOutputFormat::Json => {
                    let recurrence_json = RecurrenceJson::from_recurrence_rational(
                        &res.recurrence,
                        1,
                        initial_polys,
                        Some(RecurrenceJsonSearch {
                            recurrence_text: res.recurrence.to_string(),
                            source_rows: polys.len(),
                            skip_prefix: search.skip_prefix,
                            unknowns: res.num_unknowns,
                            weighted_unknowns: res.weighted_unknowns,
                            equations: res.num_equations,
                            fit_polynomials: res.fit_polynomials,
                            verification_polynomials: res.verification_polynomials,
                            candidates_tried: res.candidates_tried,
                            options: RecurrenceOptionsJson::from(&res.opts),
                        }),
                    );
                    println!(
                        "{}",
                        serde_json::to_string_pretty(&recurrence_json).unwrap()
                    );
                }
                RecurrenceOutputFormat::Python => {
                    println!(
                        "{}",
                        res.recurrence.to_python_definition_rational(initial_polys)
                    );
                }
                RecurrenceOutputFormat::Text => {
                    println!("{}", res.recurrence);
                    if theorem_matches {
                        print_theorem_report_text(&theorem_report);
                    }
                }
            }
            if !theorem_matches || format == RecurrenceOutputFormat::Python {
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
                    }))
                    .unwrap()
                );
            }
            eprintln!(
                "Recurrence search exhausted its candidate budget after {} candidates; \
                 unsearched candidates remain.",
                summary.diagnostics.considered_candidates
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

fn read_recurrence_json(path: &str) -> Result<RecurrenceJson, String> {
    let input = if path == "-" {
        read_limited_text(io::stdin(), "recurrence JSON from stdin")?
    } else {
        let file = fs::File::open(path)
            .map_err(|e| format!("failed to open recurrence JSON `{path}`: {e}"))?;
        read_limited_text(file, &format!("recurrence JSON `{path}`"))?
    };
    serde_json::from_str(&input).map_err(|e| format!("failed to parse recurrence JSON: {e}"))
}

pub(crate) fn cmd_recurrence_generate(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_generate_help();
        return;
    }

    let mut recurrence_path: Option<String> = None;
    let mut rows: Option<usize> = None;
    let mut additional: Option<usize> = None;
    let mut format = OutputFormat::Text;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--recurrence" => {
                i += 1;
                recurrence_path = args.get(i).cloned();
            }
            "--rows" => {
                i += 1;
                rows = args.get(i).and_then(|s| s.parse().ok());
            }
            "--additional" => {
                i += 1;
                additional = args.get(i).and_then(|s| s.parse().ok());
            }
            "--json" => {
                format = OutputFormat::Json;
            }
            other => {
                eprintln!("unknown recurrence-generate option: {other}");
                return;
            }
        }
        i += 1;
    }

    let Some(path) = recurrence_path else {
        eprintln!("recurrence-generate needs --recurrence <file|->");
        return;
    };
    if rows.is_some() && additional.is_some() {
        eprintln!("use only one of --rows or --additional");
        return;
    }

    let recurrence_json = match read_recurrence_json(&path) {
        Ok(value) => value,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let (recurrence, first_index, initial_polys) = match recurrence_json.to_recurrence_parts() {
        Ok(parts) => parts,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };

    let total_rows = if let Some(rows) = rows {
        rows
    } else if let Some(additional) = additional {
        match initial_polys.len().checked_add(additional) {
            Some(total) => total,
            None => {
                eprintln!("recurrence row count overflow");
                return;
            }
        }
    } else {
        eprintln!("recurrence-generate needs --rows <n> or --additional <n>");
        return;
    };

    let generated = match recurrence.generate_rows_rational(&initial_polys, first_index, total_rows)
    {
        Ok(rows) => rows,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };

    if format == OutputFormat::Json {
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
