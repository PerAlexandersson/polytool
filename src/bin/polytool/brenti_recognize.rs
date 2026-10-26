//! `polytool brenti-recognize`: exact search for a planar strip digraph whose
//! rows from `(0,0)` are the input polynomials.

use crate::help::print_coefficient_input_help;
use crate::input::{next_option_value, parse_usize_option, read_polys_bigint};
use crate::json_out::{json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;

/// Fixed guardrail printed with every report.
pub(crate) const BRENTI_RECOGNIZE_NOTE: &str = "A certificate shows that the supplied rows \
     are rows of the path matrix of a planar strip digraph with height-dependent \
     transitions, so by Lindström--Gessel--Viennot the coefficient matrix of the rows \
     P_0, ..., P_N is totally nonnegative. It does not show that the rows are real-rooted \
     or Pólya frequency sequences ([1 + t, 1 + t + t^2] is recognized), and it covers \
     only the supplied prefix.";

pub(crate) fn print_brenti_recognize_help() {
    println!("Usage:");
    println!("  polytool brenti-recognize [options]");
    println!("  polytool brenti-recognize --help");
    println!();
    println!("Search exactly for a planar strip digraph D (heights 0..H, edges from column");
    println!("x to column x + 1 that depend only on the height, monotone targets) whose");
    println!("weighted path counts from (0,0) to (n,k) are the coefficients [t^k] P_n of the");
    println!("input rows.  The verdict is one of:");
    println!("  certified                 a digraph, re-verified against every input row");
    println!("  impossible_within_bounds  an exact proof that no digraph within the bounds");
    println!("                            reproduces the rows");
    println!("  inconclusive              a search budget ran out");
    println!();
    print_coefficient_input_help();
    println!();
    println!("  By default the input is P_1, P_2, ... with P_0 = 1 implicit.");
    println!();
    println!("Options:");
    println!("  --first-row <0|1>          Digraph row of the first input polynomial");
    println!("                             (default 1; with 0 the first row must be 1)");
    println!("  --max-span <k>             At most k consecutive target heights per height");
    println!("                             (bounds the out-degree; default: no bound)");
    println!("  --weights <rational|integer>");
    println!("                             Nonnegative rational (default) or integer weights");
    println!("  --max-search-nodes <n>     Budget on partial supports (default 1000000)");
    println!("  --max-integer-nodes <n>    Budget on the integer search per column");
    println!("                             (default 10000)");
    println!("  --json                     Emit machine-readable JSON");
    println!("                             (schema polytool.brenti-recognize.v1)");
    println!();
    println!("Guardrails:");
    println!("  A certificate proves that the coefficient matrix of P_0, ..., P_N is totally");
    println!("  nonnegative (Lindström--Gessel--Viennot).  It does NOT prove that the rows");
    println!("  are real-rooted or PF: [1 + t, 1 + t + t^2] is certified.  The transitions");
    println!("  here depend on the height, so the digraph is invariant under shifts in x,");
    println!("  not in y; Brenti's PF theorem needs height-shift invariance.  Only the");
    println!("  supplied prefix is covered.");
    println!();
    println!("Exit status:");
    println!("  0  certified");
    println!("  1  impossible within the bounds");
    println!("  2  invalid input or options");
    println!("  3  inconclusive (a budget was exhausted)");
    println!();
    println!("Examples:");
    println!("  printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool brenti-recognize");
    println!("  printf '0,1\\n0,1,1\\n0,1,3,1\\n0,1,7,6,1\\n' | polytool brenti-recognize --weights integer");
    println!("  printf '1,1\\n1,4,1\\n1,11,11,1\\n' | polytool brenti-recognize --json");
}

struct BrentiRecognizeArgs {
    json: bool,
    options: BrentiRecognitionOptions,
}

fn parse_brenti_recognize_args(args: &[String]) -> Result<BrentiRecognizeArgs, String> {
    let mut json = false;
    let mut options = BrentiRecognitionOptions::default();
    let mut index = 0;
    while index < args.len() {
        let arg = args[index].as_str();
        match arg {
            "--json" => json = true,
            "--first-row" => options.first_row_index = parse_usize_option(args, &mut index, arg)?,
            "--max-span" => {
                options.max_target_span = Some(parse_usize_option(args, &mut index, arg)?)
            }
            "--weights" => {
                options.weights = match next_option_value(args, &mut index, arg)? {
                    "rational" => BrentiWeightDomain::Rational,
                    "integer" => BrentiWeightDomain::Integer,
                    other => {
                        return Err(format!(
                            "--weights expects 'rational' or 'integer', got '{other}'"
                        ))
                    }
                }
            }
            "--max-search-nodes" => {
                options.max_search_nodes = parse_usize_option(args, &mut index, arg)?
            }
            "--max-integer-nodes" => {
                options.max_integer_search_nodes = parse_usize_option(args, &mut index, arg)?
            }
            other => return Err(format!("unknown brenti-recognize option: {other}")),
        }
        index += 1;
    }
    Ok(BrentiRecognizeArgs { json, options })
}

pub(crate) fn weight_domain_name(domain: BrentiWeightDomain) -> &'static str {
    match domain {
        BrentiWeightDomain::Rational => "rational",
        BrentiWeightDomain::Integer => "integer",
    }
}

pub(crate) fn verdict_name(result: &BrentiRecognition) -> &'static str {
    match result {
        BrentiRecognition::Certified { .. } => "certified",
        BrentiRecognition::ImpossibleWithinBounds { .. } => "impossible_within_bounds",
        BrentiRecognition::Inconclusive { .. } => "inconclusive",
    }
}

fn digraph_json(digraph: &BrentiStripDigraph) -> String {
    let outgoing = digraph
        .outgoing
        .iter()
        .map(|edges| {
            format!(
                "[{}]",
                edges
                    .iter()
                    .map(|edge| format!(
                        "{{\"to\":{},\"weight\":{}}}",
                        edge.to,
                        json_string(&edge.weight.to_string())
                    ))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"num_heights\":{},\"outgoing\":[{}]}}",
        digraph.num_heights(),
        outgoing
    )
}

fn json_option<T>(value: Option<T>, render: impl FnOnce(T) -> String) -> String {
    value.map_or("null".to_string(), render)
}

fn result_json(
    result: &BrentiRecognition,
    options: &BrentiRecognitionOptions,
    rows: usize,
    certificate_verified: Option<bool>,
) -> String {
    let (digraph, row_indices, supports_checked, obstruction, reason) = match result {
        BrentiRecognition::Certified {
            digraph,
            certificate,
        } => (
            Some(digraph),
            Some(&certificate.row_indices),
            None,
            None,
            None,
        ),
        BrentiRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => (
            None,
            None,
            Some(*supports_checked),
            Some(obstruction.as_str()),
            None,
        ),
        BrentiRecognition::Inconclusive { reason } => {
            (None, None, None, None, Some(reason.as_str()))
        }
    };
    format!(
        "{{\"schema\":\"polytool.brenti-recognize.v1\",\"options\":{{\"first_row_index\":{},\
         \"max_target_span\":{},\"weights\":{},\"max_search_nodes\":{},\
         \"max_integer_search_nodes\":{}}},\"rows\":{},\"verdict\":{},\"certified\":{},\
         \"digraph\":{},\"row_indices\":{},\"certificate_verified\":{},\
         \"supports_checked\":{},\"obstruction\":{},\"reason\":{},\"message\":{},\"note\":{}}}",
        options.first_row_index,
        json_option(options.max_target_span, |v| v.to_string()),
        json_string(weight_domain_name(options.weights)),
        options.max_search_nodes,
        options.max_integer_search_nodes,
        rows,
        json_string(verdict_name(result)),
        result.is_certified(),
        json_option(digraph, digraph_json),
        json_option(row_indices, |r| json_usize_vec(r)),
        json_option(certificate_verified, |v| v.to_string()),
        json_option(supports_checked, |v| v.to_string()),
        json_option(obstruction, json_string),
        json_option(reason, json_string),
        json_string(&result.to_string()),
        json_string(BRENTI_RECOGNIZE_NOTE)
    )
}

fn print_result_text(
    result: &BrentiRecognition,
    options: &BrentiRecognitionOptions,
    rows: usize,
    certificate_verified: Option<bool>,
) {
    println!("Strip-digraph recognition");
    let first = options.first_row_index;
    println!(
        "- Input: {rows} row(s), digraph rows {first}..{}{}",
        first + rows.saturating_sub(1),
        if first == 1 { "; P_0 = 1 implicit" } else { "" }
    );
    println!(
        "- Bounds: {} weights, target span {}, {} search nodes, {} integer nodes per column",
        weight_domain_name(options.weights),
        options
            .max_target_span
            .map_or("unbounded".to_string(), |s| format!("at most {s}")),
        options.max_search_nodes,
        options.max_integer_search_nodes
    );
    match result {
        BrentiRecognition::Certified { digraph, .. } => {
            println!("- Verdict: certified");
            println!("- Digraph (height -> target: weight):");
            for (height, edges) in digraph.outgoing.iter().enumerate() {
                if edges.is_empty() {
                    println!("  {height} -> (no edges)");
                } else {
                    let targets = edges
                        .iter()
                        .map(|edge| format!("{}: {}", edge.to, edge.weight))
                        .collect::<Vec<_>>()
                        .join(", ");
                    println!("  {height} -> {targets}");
                }
            }
            println!(
                "- Certificate re-verified against the input rows: {}",
                if certificate_verified == Some(true) {
                    "yes"
                } else {
                    "NO"
                }
            );
        }
        BrentiRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => {
            println!(
                "- Verdict: impossible within bounds ({supports_checked} partial supports checked)"
            );
            println!("- Obstruction: {obstruction}");
        }
        BrentiRecognition::Inconclusive { reason } => {
            println!("- Verdict: inconclusive");
            println!("- Reason: {reason}");
        }
    }
    println!();
    println!("Note: {BRENTI_RECOGNIZE_NOTE}");
}

pub(crate) fn cmd_brenti_recognize(args: &[String]) {
    let parsed = match parse_brenti_recognize_args(args) {
        Ok(parsed) => parsed,
        Err(error) => {
            eprintln!("{error}");
            eprintln!("Run `polytool help brenti-recognize` for usage.");
            std::process::exit(2);
        }
    };
    let polys: Vec<Polynomial<BigInt>> = read_polys_bigint()
        .into_iter()
        .map(Polynomial::new)
        .collect();
    let result = match recognize_brenti_strip_digraph(&polys, &parsed.options) {
        Ok(result) => result,
        Err(error) => {
            eprintln!("brenti-recognize: {error}");
            std::process::exit(2);
        }
    };
    let certificate_verified = match &result {
        BrentiRecognition::Certified { certificate, .. } => Some(
            certificate
                .verify_against_bigint_polynomials(&polys)
                .is_ok(),
        ),
        _ => None,
    };
    if parsed.json {
        println!(
            "{}",
            result_json(&result, &parsed.options, polys.len(), certificate_verified)
        );
    } else {
        print_result_text(&result, &parsed.options, polys.len(), certificate_verified);
    }
    let code = match &result {
        BrentiRecognition::Certified { .. } if certificate_verified == Some(true) => 0,
        BrentiRecognition::Certified { .. } | BrentiRecognition::ImpossibleWithinBounds { .. } => 1,
        BrentiRecognition::Inconclusive { .. } => 3,
    };
    if code != 0 {
        std::process::exit(code);
    }
}
