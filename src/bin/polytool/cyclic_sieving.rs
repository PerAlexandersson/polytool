//! `polytool cyclic-sieving`: one cyclic sieving order.

use crate::help::print_coefficient_input_help;
use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint, OutputFormat};
use crate::json_out::{json_bigint_vec, json_bool_option, json_string, json_usize_vec};
use num_bigint::BigInt;
use polytool::*;

pub(crate) fn print_cyclic_sieving_help() {
    println!("Usage:");
    println!("  polytool cyclic-sieving --order <n> [options]");
    println!("  polytool cyclic-sieving --help");
    println!();
    println!("Profile or check P(q) at powers of a primitive n-th root of unity.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Options:");
    println!("  --order <n>          Cyclic group order");
    println!("  --fixed-counts <xs>  Fixed-point counts fix(g^0),...,fix(g^(n-1))");
    println!("  --orbit-sizes <xs>   Alternative to --fixed-counts: the multiset of orbit sizes");
    println!("                       (each a divisor of n) of a concrete C_n-action");
    println!("  --json               Emit machine-readable JSON");
    println!();
    println!("Without fixed counts or orbit sizes, this profiles exact root-of-unity");
    println!("evaluations only.");
    println!();
    println!("Every run also reports the orbit decomposition of the polynomial (exact integer");
    println!("arithmetic, Alexandersson-Amini, 'The cone of cyclic sieving phenomena',");
    println!("Discrete Math. 342 (2019), arXiv:1804.01447):");
    println!("  - reduce f modulo q^n-1 to residues c_0..c_{{n-1}}; f(zeta^d) is an integer for");
    println!("    all d iff c_j = c_gcd(j,n) (Desarmenien); violating gcd classes are listed;");
    println!("  - if integral, f = sum_(r|n) a_r h_r mod q^n-1 with h_r = 1+q^r+...+q^(n-r),");
    println!("    a_r = sum_(e|r) mu(r/e) c_e = number of orbits of size n/r a CSP action needs;");
    println!("    S_k = k*a_(n/k) = sum_(j|k) mu(k/j) f(zeta^j) elements lie in orbits of size k;");
    println!("  - cone condition: all a_r >= 0 (equivalently all S_k >= 0).  For f in N[q] this");
    println!("    is necessary and sufficient for a C_n-action making (X, C_n, f) a CSP.");
    println!("The hypothesis (natural coefficients) and the cone condition are reported");
    println!("separately.  With --orbit-sizes the predicted orbit counts are compared with");
    println!("the actual ones.  JSON adds 'orbit_decomposition' and 'orbit_sizes_check'");
    println!("next to the unchanged 'report' of every item; exact integers are strings.");
}

#[derive(Clone, Debug)]
struct CyclicSievingOptions {
    format: OutputFormat,
    order: Option<usize>,
    fixed_counts: Option<Vec<BigInt>>,
    orbit_sizes: Option<Vec<usize>>,
}

pub(crate) fn parse_bigint_values(input: &str) -> Result<Vec<BigInt>, String> {
    parse_polynomial_bigint(input).map_err(|error| error.to_string())
}

fn parse_cyclic_sieving_args(args: &[String]) -> Result<CyclicSievingOptions, String> {
    let mut options = CyclicSievingOptions {
        format: OutputFormat::Text,
        order: None,
        fixed_counts: None,
        orbit_sizes: None,
    };
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => options.format = OutputFormat::Json,
            "--orbit-sizes" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--orbit-sizes expects a list of orbit sizes".to_string())?;
                let sizes = value
                    .split(|c: char| c == ',' || c.is_whitespace())
                    .filter(|part| !part.is_empty())
                    .map(|part| {
                        part.parse::<usize>().map_err(|_| {
                            format!("--orbit-sizes expects positive integers, got '{part}'")
                        })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                options.orbit_sizes = Some(sizes);
            }
            "--order" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--order expects a positive integer".to_string())?;
                let order = value
                    .parse()
                    .map_err(|_| format!("--order expects a positive integer, got '{value}'"))?;
                if order == 0 {
                    return Err("--order must be positive".to_string());
                }
                options.order = Some(order);
            }
            "--fixed-counts" => {
                i += 1;
                let value = args
                    .get(i)
                    .ok_or_else(|| "--fixed-counts expects a coefficient list".to_string())?;
                options.fixed_counts = Some(parse_bigint_values(value)?);
            }
            other => return Err(format!("unknown cyclic-sieving option: {other}")),
        }
        i += 1;
    }
    if options.fixed_counts.is_some() && options.orbit_sizes.is_some() {
        return Err("--fixed-counts and --orbit-sizes are alternatives; give only one".to_string());
    }
    Ok(options)
}

fn json_orbit_entries(entries: &[OrbitSizeEntry]) -> String {
    let items = entries
        .iter()
        .map(|entry| {
            format!(
                "{{\"orbit_size\":{},\"basis_index\":{},\"orbit_count\":{},\
                 \"element_count\":{}}}",
                entry.orbit_size,
                entry.basis_index,
                json_string(&entry.orbit_count.to_string()),
                json_string(&entry.element_count.to_string())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!("[{items}]")
}

fn orbit_decomposition_json(dec: &OrbitDecomposition) -> String {
    let classes = dec
        .classes
        .iter()
        .map(|class| {
            format!(
                "{{\"gcd_class\":{},\"indices\":{},\"residues\":{},\"common_value\":{}}}",
                class.gcd_class,
                json_usize_vec(&class.indices),
                json_bigint_vec(&class.residues),
                class
                    .common_value
                    .as_ref()
                    .map(|v| json_string(&v.to_string()))
                    .unwrap_or_else(|| "null".to_string())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"order\":{},\"reduced_residues\":{},\"classes\":[{}],\"integral\":{},\
         \"violating_classes\":{},\"orbits\":{},\"fixed_points\":{},\
         \"natural_coefficients\":{},\"cone_condition\":{},\"realizable\":{}}}",
        dec.order,
        json_bigint_vec(&dec.reduced_residues),
        classes,
        dec.integral,
        json_usize_vec(&dec.violating_classes),
        dec.orbits
            .as_deref()
            .map(json_orbit_entries)
            .unwrap_or_else(|| "null".to_string()),
        dec.fixed_points
            .as_deref()
            .map(json_bigint_vec)
            .unwrap_or_else(|| "null".to_string()),
        dec.natural_coefficients,
        json_bool_option(dec.cone_condition),
        dec.realizable
    )
}

fn orbit_sizes_check_json(check: &OrbitSizesCheck) -> String {
    format!(
        "{{\"orbit_sizes\":{},\"num_elements\":{},\"actual_orbits\":{},\
         \"actual_fixed_points\":{},\"mismatched_orbit_sizes\":{},\"holds\":{}}}",
        json_usize_vec(&check.orbit_sizes),
        json_string(&check.num_elements.to_string()),
        json_orbit_entries(&check.actual_orbits),
        json_bigint_vec(&check.actual_fixed_points),
        json_usize_vec(&check.mismatched_orbit_sizes),
        check.holds
    )
}

fn print_orbit_decomposition_text(dec: &OrbitDecomposition) {
    let residues = dec
        .reduced_residues
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    println!("  f mod q^{}-1: [{residues}]", dec.order);
    if !dec.integral {
        let classes = dec
            .violating_classes
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(", ");
        println!(
            "  not integral: f(zeta^d) is not an integer for some d (non-constant gcd \
             classes: {classes}); no cyclic action can satisfy CSP"
        );
        return;
    }
    println!("  integral: yes (c_j = c_gcd(j,n))");
    for entry in dec.orbits.iter().flatten() {
        println!(
            "  orbit size {}: {} orbit(s), {} element(s)  [a_{} = {}]",
            entry.orbit_size,
            entry.orbit_count,
            entry.element_count,
            entry.basis_index,
            entry.orbit_count
        );
    }
    println!(
        "  Alexandersson-Amini: natural coefficients: {}; cone condition (all orbit counts >= 0): {}",
        if dec.natural_coefficients { "yes" } else { "no" },
        if dec.cone_condition == Some(true) { "yes" } else { "no" }
    );
    println!(
        "  verdict: {}",
        if dec.realizable {
            "a C_n-action making (X, C_n, f) a CSP exists"
        } else {
            "no C_n-action makes (X, C_n, f) a CSP"
        }
    );
}

fn print_orbit_sizes_check_text(check: &OrbitSizesCheck) {
    println!(
        "  action with orbit sizes {:?} ({} elements): predicted orbit counts {}",
        check.orbit_sizes,
        check.num_elements,
        if check.holds {
            "match the action".to_string()
        } else if !check.decomposition.integral {
            "do not exist (not integral)".to_string()
        } else {
            format!("differ for orbit sizes {:?}", check.mismatched_orbit_sizes)
        }
    );
}

fn root_evaluation_json(evaluation: &RootOfUnityEvaluation) -> String {
    format!(
        "{{\"group_order\":{},\"power\":{},\"root_order\":{},\"integer_value\":{},\
         \"remainder\":{},\"remainder_polynomial\":{}}}",
        evaluation.group_order,
        evaluation.power,
        evaluation.root_order,
        evaluation
            .integer_value
            .as_ref()
            .map(|value| json_string(&value.to_string()))
            .unwrap_or_else(|| "null".to_string()),
        json_bigint_vec(&evaluation.remainder),
        json_string(&format_poly_bigint_coeffs(&evaluation.remainder))
    )
}

fn cyclic_power_check_json(check: &CyclicSievingPowerCheck) -> String {
    format!(
        "{{\"power\":{},\"root_order\":{},\"expected_fixed_points\":{},\
         \"evaluation\":{},\"holds\":{}}}",
        check.power,
        check.root_order,
        json_string(&check.expected_fixed_points.to_string()),
        root_evaluation_json(&check.evaluation),
        check.holds
    )
}

pub(crate) fn cyclic_sieving_report_json(report: &CyclicSievingReport) -> String {
    let evaluations = report
        .evaluations
        .iter()
        .map(root_evaluation_json)
        .collect::<Vec<_>>()
        .join(",");
    let checks = report
        .checks
        .iter()
        .map(cyclic_power_check_json)
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"order\":{},\"coefficients\":{},\"polynomial\":{},\"fixed_counts\":{},\
         \"holds\":{},\"evaluations\":[{}],\"checks\":[{}]}}",
        report.order,
        json_bigint_vec(&report.coefficients),
        json_string(&format_poly_bigint_coeffs(&report.coefficients)),
        report
            .fixed_counts
            .as_ref()
            .map(|counts| json_bigint_vec(counts))
            .unwrap_or_else(|| "null".to_string()),
        json_bool_option(report.holds),
        evaluations,
        checks
    )
}

pub(crate) fn print_cyclic_sieving_report_text(row: usize, report: &CyclicSievingReport) {
    println!(
        "row {row}, order {}: {}",
        report.order,
        report
            .holds
            .map(|holds| if holds {
                "CSP checks passed"
            } else {
                "CSP checks failed"
            })
            .unwrap_or("profile only")
    );
    for evaluation in &report.evaluations {
        let value = evaluation
            .integer_value
            .as_ref()
            .map(ToString::to_string)
            .unwrap_or_else(|| format_poly_bigint_coeffs(&evaluation.remainder));
        println!(
            "  k={}: primitive order {}, value {}",
            evaluation.power, evaluation.root_order, value
        );
    }
    for check in report.checks.iter().filter(|check| !check.holds) {
        println!(
            "  FAIL k={}: expected {}, got {}",
            check.power,
            check.expected_fixed_points,
            check
                .evaluation
                .integer_value
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| format_poly_bigint_coeffs(&check.evaluation.remainder))
        );
    }
}

pub(crate) fn cmd_cyclic_sieving(args: &[String]) {
    let options = match parse_cyclic_sieving_args(args) {
        Ok(options) => options,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let Some(order) = options.order else {
        eprintln!("cyclic-sieving requires --order <n>");
        return;
    };
    let mut reports = Vec::new();
    let mut decompositions = Vec::new();
    let mut orbit_checks = Vec::new();
    for coeffs in read_polys_bigint() {
        let coeffs = strip_trailing_zeros_bigint(&coeffs);
        let check = match &options.orbit_sizes {
            Some(sizes) => match check_cyclic_sieving_orbit_sizes(coeffs, order, sizes) {
                Ok(check) => Some(check),
                Err(error) => {
                    eprintln!("{error}");
                    return;
                }
            },
            None => None,
        };
        let fixed_counts = match &check {
            Some(check) => Some(check.actual_fixed_points.as_slice()),
            None => options.fixed_counts.as_deref(),
        };
        let report = cyclic_sieving_report_bigint(coeffs, order, fixed_counts).and_then(|report| {
            cyclic_sieving_orbit_decomposition_bigint(coeffs, order)
                .map(|decomposition| (report, decomposition))
        });
        match report {
            Ok((report, decomposition)) => {
                reports.push(report);
                decompositions.push(decomposition);
                orbit_checks.push(check);
            }
            Err(error) => {
                eprintln!("{error}");
                return;
            }
        }
    }
    if options.format == OutputFormat::Json {
        let items = reports
            .iter()
            .enumerate()
            .map(|(index, report)| {
                format!(
                    "{{\"index\":{},\"report\":{},\"orbit_decomposition\":{},\
                     \"orbit_sizes_check\":{}}}",
                    index,
                    cyclic_sieving_report_json(report),
                    orbit_decomposition_json(&decompositions[index]),
                    orbit_checks[index]
                        .as_ref()
                        .map(orbit_sizes_check_json)
                        .unwrap_or_else(|| "null".to_string())
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        println!("{{\"items\":[{}]}}", items);
        return;
    }
    for (row, report) in reports.iter().enumerate() {
        print_cyclic_sieving_report_text(row, report);
        print_orbit_decomposition_text(&decompositions[row]);
        if let Some(check) = &orbit_checks[row] {
            print_orbit_sizes_check_text(check);
        }
    }
}
