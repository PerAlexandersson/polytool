//! `polytool recurrence-veronese-lifts`: parent sequences whose coefficient
//! sections are a recurrence-defined sequence (thin adapter over
//! `polytool::recurrence::find_veronese_section_lifts`).

use crate::input::{error_text, parse_usize_option};
use crate::recurrence_input::read_recurrence_json;
use polytool::recurrence::{
    find_veronese_section_lifts, format_characteristic_polynomial, BigRational,
    VeroneseSectionLiftOptions, VeroneseSectionLiftReport,
};
use polytool::Polynomial;

fn polynomial_text(row: &[BigRational]) -> String {
    Polynomial::new(row.to_vec()).to_string()
}

fn print_report_text(report: &VeroneseSectionLiftReport) {
    println!(
        "Input characteristic polynomial: chi(L, y) = {}",
        format_characteristic_polynomial(&report.input_characteristic_polynomial, "y")
    );
    if let Some(index) = report.tail_first_index {
        println!(
            "Supplied rows before index {index} do not follow the recurrence; the search uses \
             the rows from index {index} on."
        );
    }
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
) -> Result<
    (
        Option<String>,
        Option<usize>,
        VeroneseSectionLiftOptions,
        bool,
    ),
    String,
> {
    let mut path = None;
    let mut candidate = None;
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
            "--candidate" => {
                candidate = Some(parse_usize_option(args, &mut index, "--candidate")?);
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
    Ok((path, candidate, options, json))
}

pub(crate) fn cmd_recurrence_veronese_lifts(args: &[String]) {
    let outcome = parse_args(args).and_then(|(path, candidate, options, json)| {
        let record = read_recurrence_json(path.as_deref().unwrap_or("-"), candidate)?;
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
