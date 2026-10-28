//! `polytool recurrence-veronese-lifts`: parent sequences whose coefficient
//! sections are a recurrence-defined sequence (thin adapter over
//! `polytool::recurrence::find_veronese_section_lifts`).

use crate::input::{error_text, is_help_arg, parse_usize_option, read_limited_text};
use polytool::recurrence::{
    find_veronese_section_lifts, format_characteristic_polynomial, BigRational, RecurrenceJson,
    VeroneseSectionLiftOptions, VeroneseSectionLiftReport, MAX_VERONESE_SECTION_MODULUS,
};
use polytool::Polynomial;
use std::fs;
use std::io;

pub(crate) fn print_recurrence_veronese_lifts_help() {
    let defaults = VeroneseSectionLiftOptions::default();
    println!("Usage:");
    println!("  polytool recurrence-veronese-lifts [--recurrence <file|->] [options]");
    println!("  polytool recurrence-veronese-lifts --help");
    println!();
    println!("Search for parent sequences P_n(x) whose Veronese coefficient section is the");
    println!("sequence Q_n(y) given by a polytool.recurrence.v1 record (for example from");
    println!("`polytool recurrence --json`):");
    println!();
    println!(
        "  V_(r,j)(P)(y) = sum_k [x^(r k + j)] P(x) y^k,   Q_n = V_(r,j)(P_n) for all n >= n0."
    );
    println!();
    println!("If P_n is real-rooted with coefficients of one sign, every section is zero or");
    println!("real-rooted with nonpositive roots (Polya frequency sequences, Aissen--");
    println!("Schoenberg--Whitney--Edrei), so a lift transfers real-rootedness to Q_n.");
    println!();
    println!("Supported input: derivative-free terms without (-1)^n factors and with");
    println!("coefficients independent of n, a constant (or no) denominator, no");
    println!("inhomogeneous term. Anything else exits with status 2.");
    println!();
    println!("Method, for each modulus r = 2..max-modulus (at most the order d):");
    println!("  chi(L, y) = L^d - sum_k c_k(y) L^(d-k) is the characteristic polynomial of");
    println!("  Q_m = sum_k c_k(y) Q_(m-k). chi(L, x^r) is factored over Q; each monic divisor");
    println!("  G(L, x) with r deg_L G <= d is a candidate parent recurrence. G is kept when");
    println!("  N = prod_(w^r=1) G(L, w x) divides chi; for each residue j the parent seeds");
    println!("  are solved exactly (polynomial seeds required) and the first d rows are");
    println!("  compared. N | chi, chi monic and d matching rows prove the identity for all n.");
    println!();
    println!("Status of a lift:");
    println!("  proved_product_parent        order-one parent P_n = c^(n-n0) P_(n0) with c and");
    println!("                               P_(n0) real-rooted of one sign: every Q_n (n >= n0)");
    println!("                               is zero or real-rooted with nonpositive roots");
    println!("  parent_real_rootedness_open  all checked parent rows are real-rooted of one");
    println!("                               sign (finite evidence, not a proof)");
    println!("  parent_row_zero, parent_row_mixed_signs, parent_row_not_real_rooted");
    println!("                               a checked parent row fails: no proof from this lift");
    println!("Lifts are listed proved first, then open, then failing. Rejected candidates");
    println!("give a reason: norm_does_not_divide, seeds_not_determined,");
    println!("seeds_not_polynomial or row_mismatch. The search is bounded: no lift found");
    println!("does not mean that none exists.");
    println!();
    println!("Options:");
    println!("  --recurrence <file|->   Read the recurrence JSON from a file, or stdin with -");
    println!("                          (default: stdin)");
    println!(
        "  --max-modulus <r>       Largest modulus, 2..={MAX_VERONESE_SECTION_MODULUS} (default {})",
        defaults.max_modulus
    );
    println!(
        "  --parent-rows <n>       Parent rows checked exactly (default {})",
        defaults.parent_rows_to_check
    );
    println!(
        "  --max-candidates <n>    Candidate parents G tried per modulus (default {})",
        defaults.max_parent_candidates
    );
    println!("  --json                  Emit the polytool.recurrence-veronese-lifts.v1 report");
    println!("                          (exact rationals as strings; each parent as a");
    println!("                          polytool.recurrence.v1 record in parent_recurrence)");
    println!("  -h, --help              Print this help text");
    println!();
    println!("Exit status: 0 with at least one lift, 1 with none, 2 on invalid or");
    println!("unsupported input.");
    println!();
    println!("Example (Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n, Q_0 = Q_1 = 1: even sections of");
    println!("(1 + x)^n):");
    println!("  echo '{{\"schema\":\"polytool.recurrence.v1\",\"first_index\":0,");
    println!("    \"initial_polynomials\":[[\"1\"],[\"1\"]],\"recurrence\":{{\"terms\":[");
    println!("    {{\"offset\":1,\"deriv_order\":0,\"sign\":\"none\",\"coeff\":{{\"coeffs\":[[\"2\"]]}}}},");
    println!("    {{\"offset\":2,\"deriv_order\":0,\"sign\":\"none\",\"coeff\":{{\"coeffs\":[[\"-1\",\"1\"]]}}}}]}}}}' \\");
    println!("    | polytool recurrence-veronese-lifts");
}

fn polynomial_text(row: &[BigRational]) -> String {
    Polynomial::new(row.to_vec()).to_string()
}

fn print_report_text(report: &VeroneseSectionLiftReport) {
    println!(
        "Input characteristic polynomial: chi(L, y) = {}",
        format_characteristic_polynomial(&report.input_characteristic_polynomial, "y")
    );
    if report.lifts.is_empty() {
        println!("No Veronese-section lift found.");
    }
    for (number, lift) in report.lifts.iter().enumerate() {
        println!(
            "Lift {}: modulus {}, residue {} [{}]",
            number + 1,
            lift.modulus,
            lift.residue,
            lift.status.label()
        );
        println!(
            "  parent G(L, x) = {}",
            format_characteristic_polynomial(&lift.parent_characteristic_polynomial, "x")
        );
        println!(
            "  section N(L, y) = {}; cofactor chi/N = {}",
            format_characteristic_polynomial(&lift.section_characteristic_polynomial, "y"),
            format_characteristic_polynomial(&lift.cofactor, "y")
        );
        println!("  parent recurrence (t = x): {}", lift.parent_recurrence);
        for (offset, row) in lift.parent_initial_rows.iter().enumerate() {
            println!(
                "  parent seed P_{} = {}",
                lift.parent_first_index + offset,
                polynomial_text(row)
            );
        }
        let check = &lift.parent_row_check;
        println!(
            "  parent rows checked: {} (first zero: {}, first mixed signs: {}, first not \
             real-rooted: {})",
            check.rows_checked,
            optional_index(check.first_zero_row_index),
            optional_index(check.first_mixed_sign_row_index),
            optional_index(check.first_not_real_rooted_row_index)
        );
        println!("  {}", lift.status.statement());
    }
    for rejection in &report.rejected_candidates {
        let residue = rejection.residue.map_or_else(
            || "all residues".to_string(),
            |residue| format!("residue {residue}"),
        );
        println!(
            "Rejected: modulus {}, {residue}, G = {}: {}",
            rejection.modulus,
            format_characteristic_polynomial(&rejection.parent_characteristic_polynomial, "x"),
            rejection.reason.label()
        );
    }
    for skipped in &report.skipped_moduli {
        println!("Skipped modulus {}: {}", skipped.modulus, skipped.reason);
    }
    if report.exhausted_budget {
        println!("Candidate budget exhausted: raise --max-candidates to search further.");
    }
}

fn optional_index(index: Option<usize>) -> String {
    index.map_or_else(|| "none".to_string(), |index| index.to_string())
}

fn parse_args(
    args: &[String],
) -> Result<(Option<String>, VeroneseSectionLiftOptions, bool), String> {
    let mut path = None;
    let mut options = VeroneseSectionLiftOptions::default();
    let mut json = false;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--recurrence" => {
                index += 1;
                path = Some(
                    args.get(index)
                        .cloned()
                        .ok_or("--recurrence expects a file or -")?,
                );
            }
            "--max-modulus" => {
                options.max_modulus = parse_usize_option(args, &mut index, "--max-modulus")?
            }
            "--parent-rows" => {
                options.parent_rows_to_check =
                    parse_usize_option(args, &mut index, "--parent-rows")?
            }
            "--max-candidates" => {
                options.max_parent_candidates =
                    parse_usize_option(args, &mut index, "--max-candidates")?
            }
            "--json" => json = true,
            other => return Err(format!("unknown recurrence-veronese-lifts option: {other}")),
        }
        index += 1;
    }
    Ok((path, options, json))
}

pub(crate) fn cmd_recurrence_veronese_lifts(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_veronese_lifts_help();
        return;
    }
    let outcome = parse_args(args).and_then(|(path, options, json)| {
        let text = match path.as_deref() {
            None | Some("-") => read_limited_text(io::stdin(), "recurrence JSON from stdin"),
            Some(path) => fs::File::open(path)
                .map_err(|e| format!("failed to open recurrence JSON `{path}`: {e}"))
                .and_then(|file| read_limited_text(file, &format!("recurrence JSON `{path}`"))),
        }?;
        let record: RecurrenceJson = serde_json::from_str(&text)
            .map_err(|e| format!("failed to parse recurrence JSON: {e}"))?;
        let (recurrence, first_index, rows) = record.to_recurrence_parts().map_err(error_text)?;
        let report = find_veronese_section_lifts(&recurrence, &rows, first_index, &options)
            .map_err(error_text)?;
        Ok((report, json))
    });
    let (report, json) = match outcome {
        Ok(result) => result,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    };
    if json {
        println!(
            "{}",
            serde_json::to_string_pretty(&report).expect("serialize Veronese lift report")
        );
    } else {
        print_report_text(&report);
    }
    if report.lifts.is_empty() {
        std::process::exit(1);
    }
}
