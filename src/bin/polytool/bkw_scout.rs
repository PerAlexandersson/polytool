//! `polytool bkw-scout`: numerical BKW equal-modulus scouting.

use crate::input::read_limited_text;
use crate::json_out::json_string;
use polytool::*;
use std::io::{self};

pub(crate) fn bkw_scout_usage() {
    println!("Usage:");
    println!("  polytool bkw-scout [options]");
    println!("  polytool bkw-scout --help");
    println!();
    println!("Input symbol: z-coefficient polynomials in x, ascending z-degree.");
    println!("Example: F(x,z)=1-xz+z^2 is `1; -x; 1`.");
    println!();
    println!("Options:");
    println!("  --symbol <s>             Symbol coefficients, e.g. '1; -x; 1'");
    println!("                           If omitted, read the symbol from stdin.");
    println!("  --box <r0> <r1> <i0> <i1>");
    println!("                           Complex x rectangle (default: -3 3 -3 3)");
    println!("  --grid <n>               Use an n by n grid (default: 61)");
    println!("  --grid-re <n>            Number of real-axis grid samples");
    println!("  --grid-im <n>            Number of imaginary-axis grid samples");
    println!("  --top <n>                Number of candidates to print (default: 20)");
    println!("  --include-real-axis      Include samples with Im(x)=0");
    println!("  --min-imag <eps>         Minimum |Im(x)| unless real axis is included");
    println!("  --refine-steps <n>       Local coordinate-refinement steps (default: 8)");
    println!("  --no-refine              Disable local refinement");
    println!("  --tol <eps>              Positive finite Durand--Kerner root tolerance");
    println!("  --max-iter <n>           Positive Durand--Kerner max iterations");
    println!("  --format <text|json>     Output format (default: text)");
    println!("  --json                  Alias for --format json");
    println!("  --mathematica            Also print an exact Mathematica Reduce skeleton");
    println!("                          With JSON, adds a mathematica string field");
    println!("  All floating-point options must be finite; --min-imag must be nonnegative.");
}

pub(crate) fn cmd_bkw_scout(args: &[String]) {
    let mut options = BkwScoutOptions::default();
    let mut symbol_input: Option<String> = None;
    let mut format = "text".to_string();
    let mut print_mathematica = false;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "-h" | "--help" => {
                bkw_scout_usage();
                std::process::exit(2);
            }
            "--symbol" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--symbol requires an argument.");
                    std::process::exit(2);
                }
                symbol_input = Some(args[i].clone());
            }
            "--box" => {
                if i + 4 >= args.len() {
                    eprintln!("--box requires four numeric arguments.");
                    std::process::exit(2);
                }
                options.re_min = parse_f64_arg("--box re_min", &args[i + 1]);
                options.re_max = parse_f64_arg("--box re_max", &args[i + 2]);
                options.im_min = parse_f64_arg("--box im_min", &args[i + 3]);
                options.im_max = parse_f64_arg("--box im_max", &args[i + 4]);
                i += 4;
            }
            "--grid" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid requires a positive integer.");
                    std::process::exit(2);
                }
                let grid = parse_usize_arg("--grid", &args[i]);
                options.grid_re = grid;
                options.grid_im = grid;
            }
            "--grid-re" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid-re requires a positive integer.");
                    std::process::exit(2);
                }
                options.grid_re = parse_usize_arg("--grid-re", &args[i]);
            }
            "--grid-im" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--grid-im requires a positive integer.");
                    std::process::exit(2);
                }
                options.grid_im = parse_usize_arg("--grid-im", &args[i]);
            }
            "--top" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--top requires a positive integer.");
                    std::process::exit(2);
                }
                options.max_results = parse_usize_arg("--top", &args[i]);
            }
            "--include-real-axis" => {
                options.include_real_axis = true;
            }
            "--min-imag" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--min-imag requires a numeric argument.");
                    std::process::exit(2);
                }
                options.min_imaginary_abs = parse_f64_arg("--min-imag", &args[i]);
            }
            "--refine-steps" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--refine-steps requires a nonnegative integer.");
                    std::process::exit(2);
                }
                options.refine_steps = parse_usize_arg("--refine-steps", &args[i]);
            }
            "--no-refine" => {
                options.refine_steps = 0;
            }
            "--tol" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--tol requires a numeric argument.");
                    std::process::exit(2);
                }
                options.root_options.tolerance = parse_f64_arg("--tol", &args[i]);
            }
            "--max-iter" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--max-iter requires a positive integer.");
                    std::process::exit(2);
                }
                options.root_options.max_iterations = parse_usize_arg("--max-iter", &args[i]);
            }
            "--format" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("--format requires text or json.");
                    std::process::exit(2);
                }
                format = args[i].clone();
            }
            "--mathematica" => {
                print_mathematica = true;
            }
            "--json" => format = "json".to_string(),
            other => {
                eprintln!("Unknown bkw-scout option: {}", other);
                bkw_scout_usage();
                std::process::exit(2);
            }
        }
        i += 1;
    }

    if !matches!(format.as_str(), "text" | "json") {
        eprintln!("Unknown output format: {format}; expected text or json.");
        std::process::exit(2);
    }
    if options.re_min > options.re_max || options.im_min > options.im_max {
        eprintln!("Invalid box: lower bounds must be <= upper bounds.");
        std::process::exit(2);
    }
    if options.grid_re == 0 || options.grid_im == 0 || options.max_results == 0 {
        eprintln!("Grid sizes and --top must be positive.");
        std::process::exit(2);
    }
    if !(options.re_max - options.re_min).is_finite()
        || !(options.im_max - options.im_min).is_finite()
        || !(options.re_max + options.re_min).is_finite()
        || !(options.im_max + options.im_min).is_finite()
    {
        eprintln!("Box endpoints overflow grid arithmetic; choose a smaller box.");
        std::process::exit(2);
    }
    if options.min_imaginary_abs < 0.0 {
        eprintln!("--min-imag must be nonnegative.");
        std::process::exit(2);
    }
    if options.root_options.tolerance <= 0.0 || options.root_options.max_iterations == 0 {
        eprintln!("--tol and --max-iter must be positive.");
        std::process::exit(2);
    }

    let input = match symbol_input {
        Some(s) => s,
        None => match read_limited_text(io::stdin(), "BKW symbol from stdin") {
            Ok(input) => input,
            Err(error) => {
                eprintln!("{error}");
                std::process::exit(2);
            }
        },
    };
    let symbol = match BkwSymbol::parse_z_coefficient_symbol(&input) {
        Ok(symbol) => symbol,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(2);
        }
    };

    let candidates = symbol.scout_equal_modulus_locus(&options);
    match format.as_str() {
        "text" => print_bkw_text(&symbol, &options, &candidates),
        "json" => {
            if !candidates.iter().all(candidate_has_finite_json_fields) {
                eprintln!("Numerical scout produced non-finite results; cannot emit JSON. Try a smaller box or a different symbol.");
                std::process::exit(2);
            }
            print_bkw_json(&symbol, &options, &candidates, print_mathematica);
        }
        other => {
            eprintln!("Unknown output format: {}", other);
            std::process::exit(2);
        }
    }

    if print_mathematica && format == "text" {
        println!();
        println!("Mathematica exact equal-modulus skeleton:");
        println!("{}", symbol.mathematica_equal_modulus_query());
    }
}

fn parse_f64_arg(name: &str, value: &str) -> f64 {
    value
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite())
        .unwrap_or_else(|| {
            eprintln!(
                "{} expects a finite floating-point number, got '{}'.",
                name, value
            );
            std::process::exit(2);
        })
}

fn candidate_has_finite_json_fields(candidate: &BkwScoutCandidate) -> bool {
    [
        candidate.x.re,
        candidate.x.im,
        candidate.relative_modulus_gap,
        candidate.log_modulus_gap,
        candidate.root_residual,
    ]
    .iter()
    .all(|value| value.is_finite())
        && candidate
            .dominance_ratio
            .is_none_or(|value| value.is_finite() || value == f64::INFINITY)
        && candidate.roots_by_modulus.iter().all(|root| {
            [root.root.re, root.root.im, root.modulus]
                .iter()
                .all(|value| value.is_finite())
        })
}

fn parse_usize_arg(name: &str, value: &str) -> usize {
    value.parse::<usize>().unwrap_or_else(|_| {
        eprintln!("{} expects a nonnegative integer, got '{}'.", name, value);
        std::process::exit(2);
    })
}

fn print_bkw_text(symbol: &BkwSymbol, options: &BkwScoutOptions, candidates: &[BkwScoutCandidate]) {
    println!("BKW equal-modulus scout");
    println!("Symbol: F(x,z) = {}", symbol.format_symbol());
    println!(
        "Box: Re(x) in [{:.6}, {:.6}], Im(x) in [{:.6}, {:.6}], grid {} x {}",
        options.re_min,
        options.re_max,
        options.im_min,
        options.im_max,
        options.grid_re,
        options.grid_im
    );
    println!(
        "Candidates: {} (ranked by dominant-root relative modulus gap)",
        candidates.len()
    );
    if candidates.is_empty() {
        println!("No candidate points found.  Try a larger box or finer grid.");
        return;
    }
    for (rank, candidate) in candidates.iter().enumerate() {
        let dominance = candidate
            .dominance_ratio
            .map(|r| format!("{:.6e}", r))
            .unwrap_or_else(|| "n/a".to_string());
        println!(
            "{}. x = {}   gap = {:.6e}   log_gap = {:.6e}   dominance(second/third) = {}",
            rank + 1,
            candidate.x,
            candidate.relative_modulus_gap,
            candidate.log_modulus_gap,
            dominance
        );
        println!(
            "   degree_at_x = {}, converged = {}, iterations = {}, residual = {:.3e}",
            candidate.z_degree_at_x,
            candidate.converged,
            candidate.iterations,
            candidate.root_residual
        );
        if let Some((left, right)) = candidate.tied_roots() {
            println!(
                "   tied roots: z{} = {} (|z|={:.12}), z{} = {} (|z|={:.12})",
                left.index, left.root, left.modulus, right.index, right.root, right.modulus
            );
        }
        let root_line = candidate
            .roots_by_modulus
            .iter()
            .map(|root| format!("z{}:{} |.|={:.6}", root.index, root.root, root.modulus))
            .collect::<Vec<_>>()
            .join("; ");
        println!("   roots by modulus: {}", root_line);
    }
    println!();
    println!(
        "{}",
        concat!(
            "Scout warning: this ranks numerical equal-modulus candidates only; ",
            "it does not prove BKW dominance, amplitude nonvanishing, ",
            "or eventual non-real-rootedness."
        )
    );
}

fn print_bkw_json(
    symbol: &BkwSymbol,
    options: &BkwScoutOptions,
    candidates: &[BkwScoutCandidate],
    mathematica: bool,
) {
    println!("{{");
    println!("  \"symbol\": {:?},", symbol.format_symbol());
    println!(
        "  \"box\": {{\"re_min\": {}, \"re_max\": {}, \"im_min\": {}, \"im_max\": {}}},",
        options.re_min, options.re_max, options.im_min, options.im_max
    );
    println!(
        "  \"grid\": {{\"re\": {}, \"im\": {}}},",
        options.grid_re, options.grid_im
    );
    println!("  \"candidates\": [");
    for (i, candidate) in candidates.iter().enumerate() {
        let comma = if i + 1 == candidates.len() { "" } else { "," };
        println!("    {{");
        println!(
            "      \"x\": {{\"re\": {}, \"im\": {}}},",
            candidate.x.re, candidate.x.im
        );
        println!("      \"z_degree_at_x\": {},", candidate.z_degree_at_x);
        println!(
            "      \"relative_modulus_gap\": {},",
            candidate.relative_modulus_gap
        );
        println!("      \"log_modulus_gap\": {},", candidate.log_modulus_gap);
        match candidate.dominance_ratio {
            Some(ratio) if ratio.is_finite() => println!("      \"dominance_ratio\": {},", ratio),
            Some(_) => println!("      \"dominance_ratio\": \"Infinity\","),
            None => println!("      \"dominance_ratio\": null,"),
        }
        println!("      \"root_residual\": {},", candidate.root_residual);
        println!("      \"converged\": {},", candidate.converged);
        println!("      \"iterations\": {},", candidate.iterations);
        println!("      \"roots_by_modulus\": [");
        for (j, root) in candidate.roots_by_modulus.iter().enumerate() {
            let root_comma = if j + 1 == candidate.roots_by_modulus.len() {
                ""
            } else {
                ","
            };
            println!(
                "        {{\"index\": {}, \"re\": {}, \"im\": {}, \"modulus\": {}}}{}",
                root.index, root.root.re, root.root.im, root.modulus, root_comma
            );
        }
        println!("      ]");
        println!("    }}{}", comma);
    }
    println!("  ],");
    if mathematica {
        println!(
            "  \"mathematica\": {},",
            json_string(&symbol.mathematica_equal_modulus_query())
        );
    }
    println!(
        "{}",
        concat!(
            "  \"warning\": \"numerical scout only; dominance and BKW ",
            "amplitude conditions are not certified\""
        )
    );
    println!("}}");
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn non_finite_results_cannot_reach_json_numeric_fields() {
        let symbol = BkwSymbol::from_z_coefficients(vec![vec![1], vec![0], vec![1]])
            .expect("quadratic symbol");
        let options = BkwScoutOptions {
            grid_re: 1,
            grid_im: 1,
            include_real_axis: true,
            refine_steps: 0,
            ..BkwScoutOptions::default()
        };
        let candidate = symbol.scout_equal_modulus_locus(&options).remove(0);
        assert!(candidate_has_finite_json_fields(&candidate));
        let mut altered = candidate.clone();
        altered.root_residual = f64::NAN;
        assert!(!candidate_has_finite_json_fields(&altered));
        let mut altered = candidate.clone();
        altered.roots_by_modulus[0].root.re = f64::INFINITY;
        assert!(!candidate_has_finite_json_fields(&altered));
        let mut altered = candidate;
        // The existing JSON convention represents this case as "Infinity".
        altered.dominance_ratio = Some(f64::INFINITY);
        assert!(candidate_has_finite_json_fields(&altered));
        for value in [f64::NAN, f64::NEG_INFINITY] {
            altered.dominance_ratio = Some(value);
            assert!(!candidate_has_finite_json_fields(&altered));
        }
    }
}
