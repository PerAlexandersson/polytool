//! `polytool enumerate-real-rooted`: all real-rooted integer polynomials of a
//! fixed degree with prescribed leading coefficients.

use crate::input::{next_option_value, OutputFormat};
use crate::json_out::{json_bigint_vec, json_string};
use num_bigint::BigInt;
use polytool::real_rooted_enumeration::{
    enumerate_real_rooted_integer_polynomials_bigint, RealRootedEnumerationError,
    RealRootedEnumerationOptions, RealRootedEnumerationResult, RealRootedEnumerationStatus,
};

const DEFAULT_MAX_NODES: u64 = 1_000_000;
const DEFAULT_MAX_RESULTS: usize = 10_000;

pub(crate) fn print_enumerate_real_rooted_help() {
    println!("Usage:");
    println!("  polytool enumerate-real-rooted --degree <n> --leading-desc <a0,a1,...> [options]");
    println!("  polytool enumerate-real-rooted --help");
    println!();
    println!("Enumerate every real-rooted integer polynomial");
    println!("  p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n");
    println!("of degree n whose first t leading coefficients a_0, ..., a_(t-1) are given,");
    println!("with 3 <= t <= n + 1 and a_0 != 0 (either sign).  Reducible polynomials and");
    println!("polynomials with repeated roots are included.  The search is the");
    println!("Greaves--Syatriadi derivative-tower algorithm (arXiv:2504.09241) with exact,");
    println!("certified endpoints: no floating point is used.  Nothing is read from stdin.");
    println!();
    println!("Coefficient order:");
    println!("  --leading-desc is DESCENDING: a_0 (leading), a_1, a_2, ...");
    println!("  each output polynomial is ASCENDING, as in all other polytool commands:");
    println!("  c_0, c_1, ..., c_n  represents  c_0 + c_1 x + ... + c_n x^n");
    println!("  (so c_k = a_(n-k)).  Output lines can be piped into `polytool real-rooted`.");
    println!("  Results are sorted lexicographically by (a_t, ..., a_n).");
    println!();
    println!("Options:");
    println!("  --degree <n>                 Degree n >= 2 (required)");
    println!("  --leading-desc <list>        Comma- or space-separated a_0, ..., a_(t-1)");
    println!("                               (required; arbitrary-size integers)");
    println!("  --max-results <k|none>       Stop after k results (default {DEFAULT_MAX_RESULTS})");
    println!("  --max-nodes <k|none>         Stop after examining k search nodes");
    println!("                               (default {DEFAULT_MAX_NODES})");
    println!("  --max-coefficient-bits <b|none>");
    println!("                               Stop if a derivative-tower coefficient needs more");
    println!("                               than b bits (default none)");
    println!("  --format <text|json>         Output format (default text)");
    println!("  --json                       Same as --format json");
    println!();
    println!("Status:");
    println!("  complete        the list is the full set (it may be empty)");
    println!("  resource_limit  a budget was hit; the list is an initial segment of the");
    println!("                  full sorted set, and limit_reached names the budget");
    println!("  invalid_input   (exit code 2) the arguments were rejected");
    println!("  internal_error  (exit code 1) an internal certification check failed");
    println!("Errors are reported on stderr; with JSON output a status object is also");
    println!("printed on stdout.");
    println!("JSON coefficients are strings, so they stay exact beyond 2^53.");
    println!();
    println!("Example:");
    println!("  polytool enumerate-real-rooted --degree 6 --leading-desc 1,0,-15 --format json");
}

struct EnumerateOptions {
    format: OutputFormat,
    degree: usize,
    leading: Vec<BigInt>,
    limits: RealRootedEnumerationOptions,
}

fn parse_limit<T: std::str::FromStr>(value: &str, option: &str) -> Result<Option<T>, String> {
    if value == "none" {
        return Ok(None);
    }
    value
        .parse()
        .map(Some)
        .map_err(|_| format!("{option} expects a nonnegative integer or 'none', got '{value}'"))
}

fn parse_leading(value: &str) -> Result<Vec<BigInt>, String> {
    let mut out = Vec::new();
    for token in value
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|token| !token.is_empty())
    {
        out.push(
            token
                .parse::<BigInt>()
                .map_err(|_| format!("--leading-desc: invalid integer '{token}'"))?,
        );
    }
    if out.is_empty() {
        return Err("--leading-desc expects at least one integer".to_string());
    }
    Ok(out)
}

fn parse_args(args: &[String]) -> Result<EnumerateOptions, String> {
    let mut format = OutputFormat::Text;
    let mut degree = None;
    let mut leading = None;
    let mut limits = RealRootedEnumerationOptions {
        max_nodes: Some(DEFAULT_MAX_NODES),
        max_results: Some(DEFAULT_MAX_RESULTS),
        max_coefficient_bits: None,
    };
    let mut i = 0;
    while i < args.len() {
        let option = args[i].as_str();
        match option {
            "--json" => format = OutputFormat::Json,
            "--format" => {
                format = match next_option_value(args, &mut i, option)? {
                    "json" => OutputFormat::Json,
                    "text" => OutputFormat::Text,
                    other => return Err(format!("--format expects text or json, got '{other}'")),
                }
            }
            "--degree" => {
                let value = next_option_value(args, &mut i, option)?;
                degree = Some(value.parse::<usize>().map_err(|_| {
                    format!("--degree expects a nonnegative integer, got '{value}'")
                })?);
            }
            "--leading-desc" => {
                leading = Some(parse_leading(next_option_value(args, &mut i, option)?)?);
            }
            "--max-results" => {
                limits.max_results = parse_limit(next_option_value(args, &mut i, option)?, option)?;
            }
            "--max-nodes" => {
                limits.max_nodes = parse_limit(next_option_value(args, &mut i, option)?, option)?;
            }
            "--max-coefficient-bits" => {
                limits.max_coefficient_bits =
                    parse_limit(next_option_value(args, &mut i, option)?, option)?;
            }
            other => return Err(format!("unknown enumerate-real-rooted option: {other}")),
        }
        i += 1;
    }
    Ok(EnumerateOptions {
        format,
        degree: degree.ok_or("--degree is required")?,
        leading: leading.ok_or("--leading-desc is required")?,
        limits,
    })
}

fn json_option<T: ToString>(value: Option<T>) -> String {
    value.map_or_else(|| "null".to_string(), |v| v.to_string())
}

fn limits_json(limits: &RealRootedEnumerationOptions) -> String {
    format!(
        "{{\"max_nodes\":{},\"max_results\":{},\"max_coefficient_bits\":{}}}",
        json_option(limits.max_nodes),
        json_option(limits.max_results),
        json_option(limits.max_coefficient_bits)
    )
}

fn limit_reached(status: RealRootedEnumerationStatus) -> Option<&'static str> {
    match status {
        RealRootedEnumerationStatus::Complete => None,
        RealRootedEnumerationStatus::ResourceLimit(limit) => Some(limit.as_str()),
    }
}

fn result_json(options: &EnumerateOptions, result: &RealRootedEnumerationResult) -> String {
    let s = &result.stats;
    let by_level = s
        .polynomials_by_level
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let polynomials = result
        .polynomials
        .iter()
        .map(|p| json_bigint_vec(p))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        concat!(
            "{{\"command\":\"enumerate-real-rooted\",\"degree\":{},",
            "\"leading_coefficients_descending\":{},\"output_coefficient_order\":\"ascending\",",
            "\"status\":{},\"limit_reached\":{},\"limits\":{},\"count\":{},",
            "\"polynomials\":[{}],",
            "\"stats\":{{\"derivative_levels\":{},\"nodes_examined\":{},\"nodes_pruned\":{},",
            "\"exact_root_isolations\":{},\"refinement_steps\":{},\"endpoint_tests\":{},",
            "\"final_verifications\":{},\"results_emitted\":{},",
            "\"polynomials_by_level\":[{}]}}}}"
        ),
        options.degree,
        json_bigint_vec(&options.leading),
        json_string(result.status.as_str()),
        limit_reached(result.status).map_or_else(|| "null".to_string(), json_string),
        limits_json(&options.limits),
        result.polynomials.len(),
        polynomials,
        s.derivative_levels,
        s.nodes_examined,
        s.nodes_pruned,
        s.exact_root_isolations,
        s.refinement_steps,
        s.endpoint_tests,
        s.final_verifications,
        s.results_emitted,
        by_level,
    )
}

fn join_bigints(values: &[BigInt]) -> String {
    values
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ")
}

fn print_text(options: &EnumerateOptions, result: &RealRootedEnumerationResult) {
    let s = &result.stats;
    let show = |value: Option<String>| value.unwrap_or_else(|| "none".to_string());
    println!(
        "# degree {}; leading coefficients (descending a_0, a_1, ...): {}",
        options.degree,
        join_bigints(&options.leading)
    );
    println!("# output: ascending coefficients c_0, ..., c_n of c_0 + c_1 x + ... + c_n x^n");
    match limit_reached(result.status) {
        None => println!(
            "# status: complete ({} polynomials)",
            result.polynomials.len()
        ),
        Some(limit) => println!(
            "# status: resource_limit ({limit} reached; {} polynomials found so far)",
            result.polynomials.len()
        ),
    }
    println!(
        "# limits: max_nodes={} max_results={} max_coefficient_bits={}",
        show(options.limits.max_nodes.map(|v| v.to_string())),
        show(options.limits.max_results.map(|v| v.to_string())),
        show(options.limits.max_coefficient_bits.map(|v| v.to_string()))
    );
    println!(
        "# stats: nodes_examined={} nodes_pruned={} root_isolations={} refinement_steps={} \
         endpoint_tests={} final_verifications={} polynomials_by_level=[{}]",
        s.nodes_examined,
        s.nodes_pruned,
        s.exact_root_isolations,
        s.refinement_steps,
        s.endpoint_tests,
        s.final_verifications,
        s.polynomials_by_level
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ")
    );
    for p in &result.polynomials {
        println!("{}", join_bigints(p));
    }
}

fn fail(format: OutputFormat, status: &str, message: &str, code: i32) -> ! {
    if format == OutputFormat::Json {
        println!(
            "{{\"command\":\"enumerate-real-rooted\",\"status\":{},\"error\":{}}}",
            json_string(status),
            json_string(message)
        );
    }
    eprintln!("enumerate-real-rooted: {message}");
    std::process::exit(code);
}

pub(crate) fn cmd_enumerate_real_rooted(args: &[String]) {
    let wants_json = args.iter().any(|a| a == "--json")
        || args
            .windows(2)
            .any(|w| w[0] == "--format" && w[1] == "json");
    let fallback_format = if wants_json {
        OutputFormat::Json
    } else {
        OutputFormat::Text
    };
    let options = match parse_args(args) {
        Ok(options) => options,
        Err(message) => fail(fallback_format, "invalid_input", &message, 2),
    };
    let result = match enumerate_real_rooted_integer_polynomials_bigint(
        options.degree,
        &options.leading,
        &options.limits,
    ) {
        Ok(result) => result,
        Err(error @ RealRootedEnumerationError::CertificationFailure { .. }) => {
            fail(options.format, "internal_error", &error.to_string(), 1)
        }
        Err(error) => fail(options.format, "invalid_input", &error.to_string(), 2),
    };
    match options.format {
        OutputFormat::Json => println!("{}", result_json(&options, &result)),
        OutputFormat::Text => print_text(&options, &result),
    }
}
