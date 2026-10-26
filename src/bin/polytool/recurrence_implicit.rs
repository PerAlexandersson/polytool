//! `polytool recurrence-implicit`: diagonal maps and derivative antecedents,
//! implicit differential relations with reconstruction data, uniqueness
//! arguments, reconstruction certificates and exact replay.

use crate::input::{is_help_arg, read_stdin_or_exit};
use polytool::recurrence::{
    assess_generating_description, implicit_transform_recurrence, verify_recurrence_implicit,
    GeneratingDescription, ImplicitTransformReport, ImplicitTransformRequest,
};
use serde_json::json;

pub(crate) fn print_recurrence_implicit_help() {
    println!("Usage:");
    println!("  polytool recurrence-implicit [--compact] < request.json");
    println!("  polytool recurrence-implicit --verify < report.json");
    println!("  polytool recurrence-implicit --assess < description.json");
    println!("  polytool recurrence-implicit --help");
    println!();
    println!("Map a sequence with a defining rule by a diagonal coefficient map or a");
    println!("derivative antecedent, transport the rule into an implicit differential");
    println!("relation (derivatives on any row) with exceptional-index conditions, and");
    println!("assess the resulting generating description: annihilation evidence and a");
    println!("uniqueness argument are reported separately. Prints a");
    println!("polytool.recurrence-implicit.v1 report.");
    println!();
    println!("Request fields:");
    println!("  transformation   {{\"kind\": \"diagonal\", \"weight\": W}} with W one of");
    println!(
        "                     {{\"kind\":\"factorial_power\",\"power\":m}}  lambda = 1/(k!)^m"
    );
    println!("                     {{\"kind\":\"binomial\",\"m\":m}}             lambda = binom(k+m,m), Q = W_m P");
    println!("                     {{\"kind\":\"inverse_binomial\",\"m\":m}}     lambda = 1/binom(k+m,m), P = W_m Q");
    println!(
        "                     {{\"kind\":\"product\",\"scale\":\"c\",\"factorial_exponent\":p,"
    );
    println!(
        "                      \"factors\":[{{\"coeffs\":{{\"coeffs\":[[..]]}},\"exponent\":e}}]}}"
    );
    println!("                                     lambda = c (k!)^p prod f(n,k)^e");
    println!(
        "                   or {{\"kind\": \"derivative_antecedent\", \"g\": [\"g0\", ...]}}:"
    );
    println!(
        "                     Q_n with P_n = D^m(g Q_n)/m!, m = deg g (exact back substitution)"
    );
    println!("  source_recurrence  defining polytool.recurrence.v1 rule with seeds, or");
    println!("  oeis_id          a bundled OEIS family (e.g. \"A065826\")");
    println!("  row_count        rows generated from the source (default 10)");
    println!("  source_relations [{{relation, first_base}}] further relations of P, transported");
    println!(
        "  target_relations [{{relation, first_base}}] relations of Q to store (checked on rows)"
    );
    println!("  boundary_conditions [{{condition, first_base}}] coefficient conditions of Q");
    println!("  search           {{}} or bounds: fit a successor rule to Q and try to prove it");
    println!("                   by a reconstruction certificate L o T = sum H_nu S^nu M");
    println!("  limits           {{max_rows, max_degree, ...}} resource limits");
    println!();
    println!("A relation is {{\"terms\": [{{\"shift\": i, \"deriv_order\": d, \"coeff\": {{\"coeffs\": ...}}}}]}}");
    println!("meaning sum c(n,x) D^d Q[n+i] = 0 (coeffs[a][b] multiplies n^a x^b); a");
    println!(
        "condition is {{\"terms\": [{{\"shift\": i, \"coefficient_index\": k, \"coeff\": [..]}}]}}"
    );
    println!("meaning sum c(n) [x^k] Q[n+i] = 0 (coeff ascending in n).");
    println!();
    println!("Report:");
    println!("  map_conditions     weights defined and nonzero: certified_all_indices or");
    println!("                     checked_supplied_rows (uncertified, visible)");
    println!("  transports         transported relation, clearing certificate (common");
    println!("                     multiplier H, exceptional indices), reconstruction certificate");
    println!("  description        seeds, relations, exceptional_index_conditions,");
    println!("                     boundary_conditions");
    println!("  assessment         equivalent_description | certified_consequence |");
    println!(
        "                     finite_evidence | inconsistent_seeds, with the uniqueness argument"
    );
    println!("  additional_relations  stored relations outside the description (row equations)");
    println!("  independent_search finite_fit, or certified_by_reconstruction after an upgrade");
    println!();
    println!("Options:");
    println!("  --verify    Replay a report exactly and re-check every certificate");
    println!("  --assess    Assess a generating description (for example after removing a");
    println!("              condition) and print the assessment");
    println!("  --compact   Print compact JSON");
    println!("  -h, --help  Print this help text");
    println!();
    println!("Exit status: 0 on success, 2 on invalid input or a failed replay.");
    println!();
    println!("Example (canonical A065826 divided coefficientwise by k+1 gives Eulerian rows):");
    println!("  echo '{{\"oeis_id\":\"A065826\",\"row_count\":9,\"search\":{{}},");
    println!("    \"transformation\":{{\"kind\":\"diagonal\",");
    println!("    \"weight\":{{\"kind\":\"inverse_binomial\",\"m\":1}}}}}}' \\");
    println!("    | polytool recurrence-implicit");
}

pub(crate) fn cmd_recurrence_implicit(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_recurrence_implicit_help();
        return;
    }
    let mut verify = false;
    let mut assess = false;
    let mut compact = false;
    for arg in args {
        match arg.as_str() {
            "--verify" => verify = true,
            "--assess" => assess = true,
            "--compact" => compact = true,
            other => {
                eprintln!("unknown option: {other}; bounds belong in the JSON request");
                std::process::exit(2);
            }
        }
    }
    if verify && assess {
        eprintln!("give at most one of --verify and --assess");
        std::process::exit(2);
    }
    let input = read_stdin_or_exit();
    let outcome = if verify {
        serde_json::from_str::<ImplicitTransformReport>(&input)
            .map_err(|e| e.to_string())
            .and_then(|report| {
                verify_recurrence_implicit(&report)?;
                Ok(json!({
                    "verified": true,
                    "assessment": report.assessment.status,
                    "independent_search": report
                        .independent_search
                        .as_ref()
                        .map(|s| s.evidence.clone()),
                }))
            })
    } else if assess {
        serde_json::from_str::<GeneratingDescription>(&input)
            .map_err(|e| e.to_string())
            .and_then(|description| assess_generating_description(&description))
            .and_then(|assessment| serde_json::to_value(assessment).map_err(|e| e.to_string()))
    } else {
        serde_json::from_str::<ImplicitTransformRequest>(&input)
            .map_err(|e| e.to_string())
            .and_then(|request| implicit_transform_recurrence(&request))
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
