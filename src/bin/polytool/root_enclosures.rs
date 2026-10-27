//! `polytool root-enclosures`: exact rational root enclosures with
//! multiplicities and exact interval tests.

use crate::input::{is_help_arg, read_stdin_or_exit};
use polytool::recurrence::{root_enclosures, RootEnclosureRequest};

pub(crate) fn print_root_enclosures_help() {
    println!("Usage:");
    println!("  polytool root-enclosures [--compact] < request.json");
    println!("  polytool root-enclosures --help");
    println!();
    println!("Isolate the real roots of each polynomial exactly. Every distinct real root");
    println!("gets a rational enclosure [lower, upper] and its multiplicity: lower == upper");
    println!("is an exact rational root, otherwise exactly one root lies strictly inside");
    println!("and neither endpoint is a root. Midpoints are representatives only, never");
    println!("certified bounds. Reads one JSON request on stdin and prints a");
    println!("polytool.root-enclosures.v1 report.");
    println!();
    println!("Request fields:");
    println!("  polynomials   [[\"c0\", \"c1\", ...], ...] ascending rational strings");
    println!("  max_width     optional rational: refine enclosures to at most this width");
    println!("                (default refinement: width <= 1/1024)");
    println!("  interval      optional {{lower?, upper?}} with endpoints");
    println!("                {{value, per_index?, closed?}}; an omitted end is infinite and");
    println!("                value + per_index * n uses n = first_index + position");
    println!("  first_index   index of the first polynomial (default 0)");
    println!();
    println!("Each item reports degree, valuation (multiplicity of 0), real_rooted,");
    println!("real_root_count (with multiplicity), the enclosures and, with an interval,");
    println!("real_roots_in_interval and real_rooted_in_interval (decided exactly). The");
    println!("zero polynomial has no degree and satisfies interval tests vacuously.");
    println!();
    println!("Options:");
    println!("  --compact   Print compact JSON");
    println!("  -h, --help  Print this help text");
    println!();
    println!("Exit status: 0 on success, 2 on invalid input.");
    println!();
    println!("Example:");
    println!("  echo '{{\"polynomials\":[[\"0\",\"0\",\"1\",\"3\",\"3\",\"1\"]],");
    println!(
        "    \"interval\":{{\"lower\":{{\"value\":\"-1\"}},\"upper\":{{\"value\":\"0\"}}}}}}' \\"
    );
    println!("    | polytool root-enclosures");
}

pub(crate) fn cmd_root_enclosures(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_root_enclosures_help();
        return;
    }
    let mut compact = false;
    for arg in args {
        match arg.as_str() {
            "--compact" => compact = true,
            other => {
                eprintln!("unknown option: {other}; parameters belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    let input = read_stdin_or_exit();
    let outcome = serde_json::from_str::<RootEnclosureRequest>(&input)
        .map_err(|e| e.to_string())
        .and_then(|request| root_enclosures(&request).map_err(|e| e.to_string()));
    match outcome {
        Ok(report) => {
            let text = if compact {
                serde_json::to_string(&report)
            } else {
                serde_json::to_string_pretty(&report)
            };
            println!("{}", text.expect("serialize root enclosures"));
        }
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}
