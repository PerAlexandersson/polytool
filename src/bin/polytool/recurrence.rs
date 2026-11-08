//! `polytool recurrence` (search, representations, complete spaces) and
//! `polytool recurrence-generate`.

use crate::help::print_rational_coefficient_input_help;
use crate::input::error_text;
use crate::input::{
    is_help_arg, next_option_value, parse_usize_option, read_polys_rational, read_stdin_or_exit,
};
use crate::recurrence_input::read_recurrence_json;
use crate::recurrence_theorems::{
    note_theorem_matches, print_theorem_report_text, search_theorem_report,
};
use num_bigint::BigInt;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::recurrence::*;
use serde_json::{json, Value};

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
    println!("  --first-index <n>          Index of the first input row (nonnegative)");
    println!("                             With --skip-prefix, discarded rows still count.");
    println!("                             Default: renumber the retained rows from 1.");
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
    println!(
        "  --no-laurent               Do not vary the left-hand side (x^e, n-factors, charts)"
    );
    println!("  --max-laurent-x-power <e>  Largest left-hand power x^e (default 2)");
    println!("  --min-margin <n>           Require equations >= unknowns + n");
    println!("  --fit-extra-rows <n>       Extra rows beyond the first solvable prefix");
    println!("  --no-verify                Use all input rows for fitting");
    println!("  --no-modular-prefilter     Disable default modular candidate rejection");
    println!("  --max-candidates <n>       Inspect at most n candidate configurations");
    println!("  --representations          Return a menu of length/derivative tradeoffs");
    println!("                             Text: comparison table and numbered formulas; JSON: full menu.");
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
    println!("  --latex                    Emit LaTeX with subscripts, seeds and starting index");
    println!("  --format latex             Alias for --latex; coefficients remain expanded");
    println!("                             JSON also includes this presentation in latex.");
    println!("  --python                   Emit exact Python code for the recurrence");
    println!("  --format python            Alias for --python");
    println!("  --verbose                  Print candidate search details");
    println!("  --theorem-matches          Also report conditional theorem matches");
    println!("                             (weighted renewal, fixed-delay predecessor),");
    println!("                             separate from the finite fit evidence. With --json,");
    println!("                             emit polytool.recurrence-theorems.v1 (the usual");
    println!("                             record is its recurrence_json). Without the flag a");
    println!("                             match is only noted on stderr. Cannot combine with");
    println!("                             --representations; save its JSON menu, then use");
    println!(
        "                             recurrence-theorems --recurrence menu.json --candidate N."
    );
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
    println!("  polytool recurrence-generate --recurrence <file|-> --format bfile --max-terms <n>");
    println!("  polytool recurrence-generate --help");
    println!("  --candidate <index>  Select a representation menu entry (0-based)");
    println!("  Accepts recurrence JSON, theorem wrappers and representation menus.");
    println!("  --format rows is an alias for the default coefficient-list output.");
    println!();
    println!("Generate coefficient rows from recurrence JSON emitted by");
    println!("`polytool recurrence --json`.");
    println!();
    println!("Options:");
    println!("  --recurrence <file|->      Read recurrence JSON from file, or stdin with -");
    println!("  --first-index <n>          Relabel the first row and shift the rule accordingly");
    println!("                             Default: keep the JSON record's first_index.");
    println!("  --rows <n>                 Output n total rows, including initial rows");
    println!("  --additional <n>           Output initial rows plus n generated rows");
    println!("  --json                     Emit generated rows as JSON");
    println!("  --format <format>          text (default), json (--json), bfile, or an OEIS");
    println!("                             text format: oeis-example, oeis-formula,");
    println!("                             mathematica, pari, oeis-json");
    println!("  -h, --help                 Print this help text");
    println!();
    println!("OEIS text formats: the rows from --first-row on, read as an OEIS triangle");
    println!("(--offset, --reversed, --trim-zeros, --regular-triangle), with the recurrence");
    println!("index n shifted to the OEIS row number.  --rows <n> rows are shown (default 10).");
    println!("  oeis-example               %e block \"Triangle begins:\" with aligned columns");
    println!("  oeis-formula               %F sentence: the recurrence in plain ASCII with its");
    println!("                             initial rows, e.g. P(n,x) = (2+x)*P(n-1,x) + ...");
    println!("  mathematica                %t program: memoized P[n_, x_] and the flattened rows");
    println!("  pari                       %o PARI/GP program listing the same rows");
    println!("  oeis-json                  all four as JSON (schema polytool.oeis-text.v1)");
    println!("Without --offset the first row keeps its own index as its OEIS row number.");
    println!();
    println!("B-file output (--format bfile): the rows are flattened into OEIS b-file");
    println!("lines `index value`, coefficients in ascending powers within each row, and");
    println!("every coefficient must be an integer.  A term cap never splits a row.");
    println!("  --offset <k>               Index of the first term (the OEIS offset; default 0)");
    println!("  --max-terms <n>            Write at most n terms, whole rows only.  Without");
    println!("                             --rows/--additional, rows are generated until the");
    println!("                             cap is reached");
    println!("  --nearest-row              With --max-terms, stop at the whole row nearest");
    println!("                             to n terms (it may exceed n; ties go to fewer)");
    println!("                             instead of the longest prefix with at most n");
    println!("  --reversed                 Write each row from its highest coefficient down");
    println!("  --trim-zeros <none|low|high|both>");
    println!("                             Drop zeros at the low (constant) end, the high");
    println!("                             (leading) end or both ends of each row, before");
    println!("                             --reversed is applied (default none)");
    println!("  --first-row <k>            Start at the row with absolute index k (default");
    println!("                             the record's first index); also for OEIS text");
    println!("  --regular-triangle         Row i (from the first) has exactly L0 + i entries:");
    println!("                             pad shorter rows with high zeros, stop before a");
    println!("                             longer row; also for OEIS text");
    println!();
    println!("Examples:");
    println!("  polytool recurrence-generate --recurrence recurrence.json --rows 50");
    println!("  polytool recurrence-generate --recurrence eulerian.json --format bfile \\");
    println!("      --offset 1 --max-terms 1000 > b008292.txt");
    println!("  polytool recurrence-generate --recurrence rec.json --format mathematica \\");
    println!("      --first-row 1 --offset 0 --reversed --rows 11");
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum RecurrenceOutputFormat {
    Text,
    Json,
    Python,
    Latex,
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
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_help();
        return;
    }

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

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
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
            "--min-rec-len" => {
                search.min_rec_len = match parse_usize_option(args, &mut i, "--min-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-rec-len" => {
                search.max_rec_len = match parse_usize_option(args, &mut i, "--max-rec-len") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-var-deg" => {
                search.min_var_deg = match parse_usize_option(args, &mut i, "--min-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-var-deg" => {
                search.max_var_deg = match parse_usize_option(args, &mut i, "--max-var-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-idx-deg" => {
                search.min_idx_deg = match parse_usize_option(args, &mut i, "--min-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-idx-deg" => {
                search.max_idx_deg = match parse_usize_option(args, &mut i, "--max-idx-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--min-diff-deg" => {
                search.min_diff_deg = match parse_usize_option(args, &mut i, "--min-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }
                };
            }
            "--max-diff-deg" => {
                search.max_diff_deg = match parse_usize_option(args, &mut i, "--max-diff-deg") {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        std::process::exit(2);
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
                            std::process::exit(2);
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
                            std::process::exit(2);
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
                            std::process::exit(2);
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
                            std::process::exit(2);
                        }
                    };
            }
            "--denominator" | "--try-denominator" => {
                search.try_denominator = true;
            }
            "--alternating-sign" => {
                search.try_alternating_sign = true;
            }
            "--no-laurent" => {
                search.try_laurent = false;
            }
            "--max-laurent-x-power" => {
                search.max_laurent_x_power =
                    match parse_usize_option(args, &mut i, "--max-laurent-x-power") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
                        }
                    };
            }
            "--max-denom-var-deg" => {
                search.try_denominator = true;
                search.max_denom_var_deg =
                    match parse_usize_option(args, &mut i, "--max-denom-var-deg") {
                        Ok(value) => value,
                        Err(error) => {
                            eprintln!("{error}");
                            std::process::exit(2);
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
                            std::process::exit(2);
                        }
                    };
            }
            "--min-margin" => {
                search.min_margin = match parse_usize_option(args, &mut i, "--min-margin") {
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
                for item in &mut report.representations {
                    item.recurrence_json = item
                        .recurrence_json
                        .clone()
                        .with_latex()
                        .unwrap_or_else(|error| usage_error(&error.to_string()));
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
                } else {
                    print_representations_text(&report);
                }
                if report.status == RepresentationSearchStatus::BudgetExhausted {
                    std::process::exit(3);
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
            let initial_count = res.export_initial_count(searched_polys.len());
            let initial_polys = &searched_polys[..initial_count];
            let rule = res.recurrence.shift_index(1 - output_first_index as i64);
            let mut theorem_report =
                search_theorem_report(&res, &polys, &search, output_first_index);
            let recurrence_json = theorem_report
                .recurrence_json
                .clone()
                .with_latex()
                .unwrap_or_else(|error| usage_error(&error.to_string()));
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
                RecurrenceOutputFormat::Python => {
                    println!(
                        "{}",
                        rule.to_python_definition_rational_indexed(
                            initial_polys,
                            output_first_index as i64
                        )
                    );
                }
                RecurrenceOutputFormat::Text => {
                    println!("{rule}");
                    if theorem_matches {
                        print_theorem_report_text(&theorem_report);
                    }
                }
            }
            if !theorem_matches
                || matches!(
                    format,
                    RecurrenceOutputFormat::Python | RecurrenceOutputFormat::Latex
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
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_generate_help();
        return;
    }

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
            "--rows" => {
                rows = Some(
                    parse_usize_option(args, &mut i, "--rows").unwrap_or_else(|e| usage_error(&e)),
                );
            }
            "--additional" => {
                additional = Some(
                    parse_usize_option(args, &mut i, "--additional")
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
            "--offset" => {
                let value =
                    next_option_value(args, &mut i, "--offset").unwrap_or_else(|e| usage_error(&e));
                bfile.options.first_index = value.parse().unwrap_or_else(|_| {
                    usage_error(&format!("--offset expects an integer, got '{value}'"))
                });
                reading_given = true;
                offset_given = true;
            }
            "--first-row" => {
                first_row = Some(
                    parse_usize_option(args, &mut i, "--first-row")
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
            "--reversed" => {
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

    let Some(path) = recurrence_path else {
        eprintln!("recurrence-generate needs --recurrence <file|->");
        std::process::exit(2);
    };
    if rows.is_some() && additional.is_some() {
        eprintln!("use only one of --rows or --additional");
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
            "--offset, --first-row, --reversed, --trim-zeros and --regular-triangle need \
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
                .unwrap_or_else(|_| usage_error("--first-row does not fit in i64"))
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
        // Read the rows from `--first-row` on as OEIS rows, then flatten them.
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
                    "recurrence-generate --format bfile needs --rows <n>, --additional <n> or --max-terms <n>",
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
                    usage_error("10000 consecutive rows without b-file terms; use --rows");
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
        eprintln!("recurrence-generate needs --rows <n> or --additional <n>");
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
                "recurrence-generate --format bfile needs --rows <n>, --additional <n> or --max-terms <n>"
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
                eprintln!("10000 consecutive rows without b-file terms; use --rows");
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
