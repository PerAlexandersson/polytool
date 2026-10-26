//! `polytool recurrence-transform`: exact coordinate changes of a polynomial
//! sequence, transport of its defining recurrence, and exact replay.

use crate::input::{is_help_arg, read_stdin_or_exit};
use polytool::recurrence::{
    transform_recurrence, verify_recurrence_transform, RecurrenceTransformReport,
    RecurrenceTransformRequest,
};
use serde_json::json;

pub(crate) fn print_recurrence_transform_help() {
    println!("Usage:");
    println!("  polytool recurrence-transform [--compact] < request.json");
    println!("  polytool recurrence-transform --verify < report.json");
    println!("  polytool recurrence-transform --help");
    println!();
    println!("Apply an exact change of polynomial coordinates to a sequence, transport");
    println!("its defining recurrence into the new coordinates, and optionally search");
    println!("independently there. Reads one JSON request from stdin and prints a");
    println!("polytool.recurrence-transform.v1 report.");
    println!();
    println!("Request fields:");
    println!("  transformation   {{\"kind\": ...}} with kind one of");
    println!("                     multiply_power {{power}}      Q_n = x^r P_n");
    println!(
        "                     divide_power {{power}}        Q_n = P_n / x^r (needs x^r | P_n)"
    );
    println!("                     reversal {{schedule}}         Q_n = x^(d_n) P_n(1/x)");
    println!("                     affine {{a, b}}               Q_n = P_n(a x + b), a != 0");
    println!("                     mobius {{a, b, c, e, schedule}}");
    println!("                                     Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))");
    println!(
        "                     f_polynomial {{schedule}}     Q_n = (1 + x)^(d_n) P_n(x/(1 + x))"
    );
    println!("                     gamma {{schedule}}            gamma coordinates at symmetry index d_n");
    println!("                     symmetric_decomposition {{schedule}}");
    println!("                                     (a_n, b_n) with P_n = a_n + x b_n");
    println!("                   a, b, c, e are exact rationals as strings (\"3/2\").");
    println!("  schedule         {{first_index, first_degree, step, last_index?}}:");
    println!("                   d_n = first_degree + step*(n - first_index) on");
    println!("                   first_index <= n (<= last_index); decreasing needs last_index");
    println!("  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or");
    println!("  oeis_id          a bundled OEIS family (e.g. \"A008292\"), or");
    println!("  rows             explicit rows (ascending rational strings) with row_origin");
    println!("  row_count        rows generated from the source (default 10)");
    println!("  search           {{}} or bounds {{max_rec_len, max_var_deg, max_idx_deg,");
    println!("                   max_diff_deg, try_denominator, max_candidates, ...}}:");
    println!("                   also search independently in the new coordinates");
    println!("  limits           {{max_rows, max_degree, ...}} resource limits");
    println!();
    println!("Report:");
    println!("  transformation      exact row map, inverse, image conditions, domain");
    println!("  transformed_states  [row][component][coefficient]");
    println!("  reconstruction_verified  inverse(map(P_n)) == P_n for every row");
    println!("  conditions          certified_all_indices (from the source rule) or");
    println!("                      checked_supplied_rows (finite evidence)");
    println!("  transported         evidence transported_identity: exact consequence of");
    println!("                      the source rule, with transformed seeds");
    println!("  independent_search  evidence finite_fit, kept separately");
    println!();
    println!("Options:");
    println!("  --verify    Replay a report exactly from its embedded request");
    println!("  --compact   Print compact JSON");
    println!("  -h, --help  Print this help text");
    println!();
    println!("Exit status: 0 on success, 2 on invalid input or a failed replay.");
    println!();
    println!("Example:");
    println!("  echo '{{\"oeis_id\":\"A008292\",\"row_count\":8,\"search\":{{}},");
    println!("    \"transformation\":{{\"kind\":\"f_polynomial\",");
    println!("    \"schedule\":{{\"first_index\":1,\"first_degree\":0,\"step\":1}}}}}}' \\");
    println!("    | polytool recurrence-transform");
}

pub(crate) fn cmd_recurrence_transform(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_transform_help();
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
        serde_json::from_str::<RecurrenceTransformReport>(&input)
            .map_err(|e| e.to_string())
            .and_then(|report| {
                verify_recurrence_transform(&report)?;
                Ok(json!({
                    "verified": true,
                    "transported": report.transported.is_some(),
                    "independent_search": report.independent_search.is_some(),
                }))
            })
    } else {
        serde_json::from_str::<RecurrenceTransformRequest>(&input)
            .map_err(|e| e.to_string())
            .and_then(|request| transform_recurrence(&request))
            .and_then(|report| serde_json::to_value(report).map_err(|e| e.to_string()))
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
