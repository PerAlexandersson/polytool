//! `polytool pf-network`: exact search for a height-shift-invariant planar
//! network whose rows are the input polynomials (a Pólya frequency
//! certificate for each row).

use crate::brenti_recognize::weight_domain_name;
use crate::help::print_coefficient_input_help;
use crate::input::{next_option_value, parse_usize_option, read_polys_bigint};
use crate::json_out::{json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;

/// Fixed guardrail printed with every report.
pub(crate) const PF_NETWORK_NOTE: &str = "A certificate shows that each supplied row is a \
     row of a planar network that is invariant under height shifts, so the Toeplitz matrix \
     of each row is the path matrix of fully compatible sources and sinks and is totally \
     nonnegative (Lindström--Gessel--Viennot; Brenti 1995): every supplied row is a Pólya \
     frequency sequence, i.e. real-rooted with only nonpositive zeros \
     (Aissen--Schoenberg--Whitney). It covers only the supplied rows, and it does not \
     certify that the coefficient matrix of the rows is TNN (tnn-network and \
     brenti-recognize concern that matrix, not the rows).";

pub(crate) fn print_pf_network_help() {
    println!("Usage:");
    println!("  polytool pf-network [options]");
    println!("  polytool pf-network --help");
    println!();
    println!("Search exactly for a planar network D on {{0..N}} x Z that is invariant under");
    println!("height shifts: for each column n it has steps (length, rise, weight) standing");
    println!("for the edges (n - length, y) -> (n, y + rise) at every height y, so");
    println!("  P_n = sum over steps into n of weight * t^rise * P_(n - length),  P_0 = 1.");
    println!("Edges are straight segments; the network is valid when every positive-weight");
    println!("step is primitive (gcd(length, rise) = 1) and no two edges meet away from a");
    println!("common endpoint.  Then the Toeplitz matrix of every row is a path matrix of");
    println!("fully compatible sources (0, i) and sinks (n, j), hence totally nonnegative,");
    println!("so every row is a Pólya frequency (PF) sequence: real-rooted with only");
    println!("nonpositive zeros.  The certificate is re-verified exactly: structure,");
    println!("planarity, every row, and the explicitly rebuilt (deg+1) x (deg+1) path");
    println!("matrix against the Toeplitz window of the row.");
    println!();
    println!("The search tries, column by column, every inclusion-minimal set of steps that");
    println!("is planar with the edges already chosen and admits nonnegative weights for");
    println!("that row (exact linear feasibility).  The verdict is one of:");
    println!("  certified                 a network, re-verified against every input row");
    println!("  impossible_within_bounds  an exact proof that no network within the bounds");
    println!("                            reproduces the rows (a row that is not real-rooted");
    println!("                            rules out every network)");
    println!("  inconclusive              a search budget ran out");
    println!();
    print_coefficient_input_help();
    println!();
    println!("  By default the input is P_1, P_2, ... with P_0 = 1 implicit.");
    println!();
    println!("Options:");
    println!("  --first-row <0|1>          Network row of the first input polynomial");
    println!("                             (default 1; with 0 the first row must be 1)");
    println!("  --max-length <L>           Largest step length (default 3)");
    println!("  --max-rise <R>             Largest rise (default 3)");
    println!("  --weights <rational|integer>");
    println!("                             Nonnegative rational (default) or integer weights");
    println!("  --max-search-nodes <n>     Budget on column supports (default 1000000)");
    println!("  --max-integer-nodes <n>    Budget on the integer search per row");
    println!("                             (default 10000)");
    println!("  --json                     Emit machine-readable JSON");
    println!("                             (schema polytool.pf-network.v1)");
    println!();
    println!("Guardrails:");
    println!("  A certificate proves that every supplied row is PF (real-rooted, nonpositive");
    println!("  zeros).  It does NOT prove that the coefficient matrix of the rows is TNN,");
    println!("  and it says nothing about rows beyond the input.  Conversely, tnn-network");
    println!("  and brenti-recognize certify TNN of the coefficient matrix, not PF rows.");
    println!("  impossible_within_bounds is relative to --max-length and --max-rise unless");
    println!("  the obstruction is a row that is not real-rooted.  See also");
    println!("  left-production for the Chen--Fu--Ruan certificate of a whole triangle.");
    println!();
    println!("Exit status:");
    println!("  0  certified");
    println!("  1  impossible within the bounds");
    println!("  2  invalid input or options");
    println!("  3  inconclusive (a budget was exhausted)");
    println!();
    println!("Examples:");
    println!("  printf '1,1\\n1,2,1\\n1,3,3,1\\n' | polytool pf-network");
    println!("  printf '1,1\\n1,4,1\\n1,11,11,1\\n' | polytool pf-network --json");
    println!("  printf '1,1,1\\n' | polytool pf-network");
}

struct PfNetworkArgs {
    json: bool,
    options: PfNetworkOptions,
}

fn parse_pf_network_args(args: &[String]) -> Result<PfNetworkArgs, String> {
    let mut json = false;
    let mut options = PfNetworkOptions::default();
    let mut index = 0;
    while index < args.len() {
        let arg = args[index].as_str();
        match arg {
            "--json" => json = true,
            "--first-row" => options.first_row_index = parse_usize_option(args, &mut index, arg)?,
            "--max-length" => options.max_step_length = parse_usize_option(args, &mut index, arg)?,
            "--max-rise" => options.max_rise = parse_usize_option(args, &mut index, arg)?,
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
            other => return Err(format!("unknown pf-network option: {other}")),
        }
        index += 1;
    }
    Ok(PfNetworkArgs { json, options })
}

fn verdict_name(result: &PfNetworkRecognition) -> &'static str {
    match result {
        PfNetworkRecognition::Certified { .. } => "certified",
        PfNetworkRecognition::ImpossibleWithinBounds { .. } => "impossible_within_bounds",
        PfNetworkRecognition::Inconclusive { .. } => "inconclusive",
    }
}

fn network_json(network: &HeightInvariantNetwork) -> String {
    let steps = network
        .steps
        .iter()
        .map(|column| {
            format!(
                "[{}]",
                column
                    .iter()
                    .map(|step| format!(
                        "{{\"length\":{},\"rise\":{},\"weight\":{}}}",
                        step.length,
                        step.rise,
                        json_string(&step.weight.to_string())
                    ))
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"last_column\":{},\"steps\":[{}]}}",
        network.last_column(),
        steps
    )
}

fn json_option<T>(value: Option<T>, render: impl FnOnce(T) -> String) -> String {
    value.map_or("null".to_string(), render)
}

fn result_json(
    result: &PfNetworkRecognition,
    options: &PfNetworkOptions,
    rows: usize,
    certificate_verified: Option<bool>,
) -> String {
    let (network, row_indices, supports_checked, obstruction, reason) = match result {
        PfNetworkRecognition::Certified { certificate } => (
            Some(&certificate.network),
            Some(&certificate.row_indices),
            None,
            None,
            None,
        ),
        PfNetworkRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => (
            None,
            None,
            Some(*supports_checked),
            Some(obstruction.as_str()),
            None,
        ),
        PfNetworkRecognition::Inconclusive { reason } => {
            (None, None, None, None, Some(reason.as_str()))
        }
    };
    format!(
        "{{\"schema\":\"polytool.pf-network.v1\",\"options\":{{\"first_row_index\":{},\
         \"max_step_length\":{},\"max_rise\":{},\"weights\":{},\"max_search_nodes\":{},\
         \"max_integer_search_nodes\":{}}},\"rows\":{},\"verdict\":{},\"certified\":{},\
         \"network\":{},\"row_indices\":{},\"certificate_verified\":{},\
         \"supports_checked\":{},\"obstruction\":{},\"reason\":{},\"message\":{},\"note\":{}}}",
        options.first_row_index,
        options.max_step_length,
        options.max_rise,
        json_string(weight_domain_name(options.weights)),
        options.max_search_nodes,
        options.max_integer_search_nodes,
        rows,
        json_string(verdict_name(result)),
        result.is_certified(),
        json_option(network, network_json),
        json_option(row_indices, |r| json_usize_vec(r)),
        json_option(certificate_verified, |v| v.to_string()),
        json_option(supports_checked, |v| v.to_string()),
        json_option(obstruction, json_string),
        json_option(reason, json_string),
        json_string(&result.to_string()),
        json_string(PF_NETWORK_NOTE)
    )
}

fn print_result_text(
    result: &PfNetworkRecognition,
    options: &PfNetworkOptions,
    rows: usize,
    certificate_verified: Option<bool>,
) {
    println!("Height-invariant planar network (PF rows)");
    let first = options.first_row_index;
    println!(
        "- Input: {rows} row(s), network rows {first}..{}{}",
        first + rows.saturating_sub(1),
        if first == 1 { "; P_0 = 1 implicit" } else { "" }
    );
    println!(
        "- Bounds: {} weights, step length <= {}, rise <= {}, {} search nodes, {} integer \
         nodes per row",
        weight_domain_name(options.weights),
        options.max_step_length,
        options.max_rise,
        options.max_search_nodes,
        options.max_integer_search_nodes
    );
    match result {
        PfNetworkRecognition::Certified { certificate } => {
            println!("- Verdict: certified (every supplied row is PF)");
            println!("- Steps into column n (length, rise): weight");
            for (index, column) in certificate.network.steps.iter().enumerate() {
                let steps = if column.is_empty() {
                    "(no steps)".to_string()
                } else {
                    column
                        .iter()
                        .map(|step| format!("({}, {}): {}", step.length, step.rise, step.weight))
                        .collect::<Vec<_>>()
                        .join(", ")
                };
                println!("  n={}: {steps}", index + 1);
            }
            println!(
                "- Certificate re-verified (planarity, rows, Toeplitz windows): {}",
                if certificate_verified == Some(true) {
                    "yes"
                } else {
                    "NO"
                }
            );
        }
        PfNetworkRecognition::ImpossibleWithinBounds {
            supports_checked,
            obstruction,
        } => {
            println!(
                "- Verdict: impossible within bounds ({supports_checked} column supports checked)"
            );
            println!("- Obstruction: {obstruction}");
        }
        PfNetworkRecognition::Inconclusive { reason } => {
            println!("- Verdict: inconclusive");
            println!("- Reason: {reason}");
        }
    }
    println!();
    println!("Note: {PF_NETWORK_NOTE}");
}

pub(crate) fn cmd_pf_network(args: &[String]) {
    let parsed = match parse_pf_network_args(args) {
        Ok(parsed) => parsed,
        Err(error) => {
            eprintln!("{error}");
            eprintln!("Run `polytool help pf-network` for usage.");
            std::process::exit(2);
        }
    };
    let polys: Vec<Polynomial<BigInt>> = read_polys_bigint()
        .into_iter()
        .map(Polynomial::new)
        .collect();
    let result = match recognize_pf_network(&polys, &parsed.options) {
        Ok(result) => result,
        Err(error) => {
            eprintln!("pf-network: {error}");
            std::process::exit(2);
        }
    };
    let certificate_verified = match &result {
        PfNetworkRecognition::Certified { certificate } => Some(
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
        PfNetworkRecognition::Certified { .. } if certificate_verified == Some(true) => 0,
        PfNetworkRecognition::Certified { .. }
        | PfNetworkRecognition::ImpossibleWithinBounds { .. } => 1,
        PfNetworkRecognition::Inconclusive { .. } => 3,
    };
    if code != 0 {
        std::process::exit(code);
    }
}
