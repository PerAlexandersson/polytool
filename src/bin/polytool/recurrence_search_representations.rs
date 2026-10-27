//! `polytool recurrence-search-representations`: the bounded queue search
//! over equivalent representations (coordinate changes, transported and
//! certified relations, fits) with exact replay.

use crate::input::error_text;
use crate::input::{is_help_arg, read_stdin_or_exit};
use polytool::recurrence::{
    search_representations, verify_representation_search, RepresentationSearchReport,
    RepresentationSearchRequest, EQUIVALENT_DESCRIPTION_LABEL,
};
use serde_json::json;

pub(crate) fn print_recurrence_search_representations_help() {
    println!("Usage:");
    println!("  polytool recurrence-search-representations [--compact] < request.json");
    println!("  polytool recurrence-search-representations --verify < report.json");
    println!("  polytool recurrence-search-representations --help");
    println!();
    println!("Search a bounded queue of equivalent representations of one polynomial");
    println!("sequence: apply the allowed coordinate changes, transport the known");
    println!("certified relations, certify complete spaces after factorial maps, fit the");
    println!("new coordinates, and feed certified relations into later steps. Prints a");
    println!("polytool.recurrence-representation-search.v1 report.");
    println!();
    println!("Request fields:");
    println!("  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or");
    println!("  oeis_id            a bundled OEIS family (e.g. \"A007318\"), or");
    println!("  rows               rows only (exploratory: every description is a finite fit)");
    println!(
        "  row_origin         original index of rows[0] (default: source first index, else 0)"
    );
    println!("  row_count          rows generated from a source (default 12)");
    println!("  transformations    the allowed maps (recurrence-transform / recurrence-implicit");
    println!("                     registry: multiply_power, divide_power, reversal, affine,");
    println!("                     mobius, f_polynomial, gamma, symmetric_decomposition,");
    println!("                     diagonal, derivative_antecedent)");
    println!(
        "                     and the lossy derivative {{\"kind\":\"derivative\",\"order\":m}}:"
    );
    println!("                     its lost coefficients are recorded (lost_data) and the parent");
    println!("                     is assessed through Q with and without them");
    println!("  search             fit bounds {{max_rec_len, max_var_deg, max_idx_deg,");
    println!("                     max_diff_deg, try_denominator, max_candidates}}");
    println!("                     (default 2, 2, 1, 1, true, 2000)");
    println!("  fits               run the fit menu in every representation (default true)");
    println!("  certified_spaces   certify complete spaces after factorial maps (default true),");
    println!("                     and after Möbius-class maps, polynomial diagonal weights and");
    println!("                     derivative maps (implicit_spaces)");
    println!("  space_shapes       explicit complete-space supports (default: a bounded menu)");
    println!("  root_domain        optional root-domain analysis of every coordinate change");
    println!("  limits             {{max_depth (2), max_representations (16), max_steps (160),");
    println!("                     max_sources_per_edge (2), max_space_shapes (16), max_rows,");
    println!("                     max_degree, max_search_candidates, max_search_degree,");
    println!("                     max_operations, max_coefficient_bits}}");
    println!();
    println!("Report:");
    println!("  representations    the source (0) and every new representation: chain of maps,");
    println!("                     certified map normal form, rows, exact reconstruction of");
    println!("                     the source, embedded replayable steps and complete spaces,");
    println!("                     and ranked descriptions labelled");
    println!("                       finite_fit | certified_annihilating_relation |");
    println!("                       equivalent_generating_description");
    println!("  edges              every map application: new_representation | duplicate_map");
    println!("                     (certified equal normal forms; rows are never compared) |");
    println!("                     rejected | skipped_budget");
    println!("  ranked             all descriptions by shift length, derivative order,");
    println!("                     coefficient degree, index degree, components, reconstruction");
    println!("  budget, status     limits that cut work; status complete | budget_exhausted");
    println!("  steps[k].implicit_relation  a certified (also implicit-only) parent relation");
    println!("                     transported through diagonal, derivative_antecedent or");
    println!("                     derivative maps, with carried conditions and, for lossy maps,");
    println!("                     the parent's description with and without the lost data");
    println!("  steps[k].vector_transport  a certified vector system of a vector parent");
    println!("                     transported exactly through a componentwise Möbius-class map");
    println!("  implicit_spaces    polytool.recurrence-implicit-search.v1 complete spaces after");
    println!("                     Möbius-class maps, polynomial diagonal weights and derivatives");
    println!();
    println!("Vector representations: after symmetric_decomposition the pair (a_n, b_n) is");
    println!("mapped componentwise by every scalar bijective map (a schedule acts on b at");
    println!("d_n - 1, its symmetry index), fitted by first-order vector systems (vector_fit),");
    println!("and its certified systems are transported through Möbius-class maps. Vector");
    println!("descriptions carry vector_assessment (uniqueness: the shared left factor is a");
    println!("nonzero polynomial at every index). A rank-two gamma transport is reported as");
    println!("transported_gamma_system (even part and odd companion), certified only where the");
    println!("symmetry of the rows is certified.");
    println!();
    println!("Options:");
    println!("  --verify    Replay a report: rerun the search, replay every step and space");
    println!("              certificate, re-derive every representation from the source");
    println!("  --compact   Print compact JSON");
    println!("  -h, --help  Print this help text");
    println!();
    println!("Exit status: 0 on success, 2 on invalid input or a failed replay.");
    println!();
    println!("Example (A007318 divided coefficientwise by k!):");
    println!("  echo '{{\"oeis_id\":\"A007318\",\"limits\":{{\"max_depth\":1}},");
    println!("    \"transformations\":[{{\"kind\":\"diagonal\",");
    println!("    \"weight\":{{\"kind\":\"factorial_power\",\"power\":1}}}}]}}' \\");
    println!("    | polytool recurrence-search-representations");
}

pub(crate) fn cmd_recurrence_search_representations(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_search_representations_help();
        return;
    }
    let mut verify = false;
    let mut compact = false;
    for arg in args {
        match arg.as_str() {
            "--verify" => verify = true,
            "--compact" => compact = true,
            other => {
                eprintln!("unknown option: {other}; bounds belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    let input = read_stdin_or_exit();
    let outcome = if verify {
        serde_json::from_str::<RepresentationSearchReport>(&input)
            .map_err(error_text)
            .and_then(|report| {
                verify_representation_search(&report).map_err(error_text)?;
                let equivalent = report
                    .representations
                    .iter()
                    .flat_map(|n| &n.descriptions)
                    .filter(|d| d.label == EQUIVALENT_DESCRIPTION_LABEL)
                    .count();
                Ok(json!({
                    "verified": true,
                    "status": report.status,
                    "representations": report.representations.len(),
                    "equivalent_descriptions": equivalent,
                }))
            })
    } else {
        serde_json::from_str::<RepresentationSearchRequest>(&input)
            .map_err(error_text)
            .and_then(|request| search_representations(&request).map_err(error_text))
            .and_then(|report| serde_json::to_value(report).map_err(error_text))
    };
    match outcome {
        Ok(value) => {
            let text = if compact {
                serde_json::to_string(&value)
            } else {
                serde_json::to_string_pretty(&value)
            };
            println!("{}", text.expect("serialize JSON value"));
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}
