//! `polytool pf-pencil`: built-in PF-bidiagonal/Jensen endpoint pencils.

use crate::input::{next_option_value, parse_usize_option, OutputFormat};
use crate::json_out::{json_bigint_vec, json_string};
use num_bigint::BigInt;
use polytool::numeric::binomial_bigint;
use polytool::*;

pub(crate) fn print_pf_pencil_help() {
    println!("Usage:");
    println!("  polytool pf-pencil --case <OEIS-id> --degree <d> [options]");
    println!("  polytool pf-pencil --case <OEIS-id> --min-degree <a> --max-degree <b> [options]");
    println!("  polytool pf-pencil --all-family-h --max-degree <b> [options]");
    println!("  polytool pf-pencil --help");
    println!();
    println!("Check built-in Family H PF-bidiagonal/Jensen endpoint pencils.");
    println!();
    println!("Options:");
    println!("  --case <OEIS-id>       Built-in case such as A036969");
    println!("  --degree <d>           Check one old-row degree");
    println!("  --min-degree <d>       First degree for a range (default: 1)");
    println!("  --max-degree <d>       Last degree for a range");
    println!("  --all-family-h         Check all built-in Family H cases");
    println!("  --lambdas <list>       Comma-separated nonnegative integer lambdas");
    println!("                         (default: 0,1,10,100)");
    println!("  --json                 Emit machine-readable JSON");
    println!("  -h, --help             Print this help text");
    println!();
    println!("Supported cases:");
    println!("  A036969 A071951 A080248 A156289 A160562 A269945");
    println!("  A166960 A166961 A166962 A166972 A191935");
    println!("  A371081 A371259 A390433 A198204");
}

#[derive(Clone, Copy, Debug)]
struct PfPencilCase {
    id: &'static str,
    alpha_label: &'static str,
    beta_label: &'static str,
}

const PF_PENCIL_CASES: &[PfPencilCase] = &[
    PfPencilCase {
        id: "A036969",
        alpha_label: "(k+1)^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A071951",
        alpha_label: "(k+1)(k+2)",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A080248",
        alpha_label: "(k+1)(k+2)/2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A156289",
        alpha_label: "(k+1)^2",
        beta_label: "2k+3",
    },
    PfPencilCase {
        id: "A160562",
        alpha_label: "(2k+1)^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A269945",
        alpha_label: "k^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A166960",
        alpha_label: "(k+1)^2",
        beta_label: "d+3-k",
    },
    PfPencilCase {
        id: "A166961",
        alpha_label: "(k+1)(2k+1)",
        beta_label: "2d+5-2k",
    },
    PfPencilCase {
        id: "A166962",
        alpha_label: "(k+1)(3k+1)",
        beta_label: "3d+1-3k",
    },
    PfPencilCase {
        id: "A166972",
        alpha_label: "(k+1)(3k+1)",
        beta_label: "d+1-k",
    },
    PfPencilCase {
        id: "A191935",
        alpha_label: "1",
        beta_label: "(d+2-k)(d+3-k)",
    },
    PfPencilCase {
        id: "A371081",
        alpha_label: "(d+k+2)^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A371259",
        alpha_label: "(d+k+4)^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A390433",
        alpha_label: "(d+k)^2",
        beta_label: "1",
    },
    PfPencilCase {
        id: "A198204",
        alpha_label: "(k+1)(k+d)/d",
        beta_label: "(k+d+1)(k+d)/d",
    },
];

#[derive(Clone, Debug)]
struct LambdaCheck {
    lambda: BigInt,
    coefficients: Vec<BigInt>,
    real_rooted: bool,
}

#[derive(Clone, Debug)]
struct PfPencilReport {
    case: PfPencilCase,
    degree: usize,
    common_denominator: BigInt,
    j_alpha: Vec<BigInt>,
    t_j_beta: Vec<BigInt>,
    j_alpha_real_rooted: bool,
    t_j_beta_real_rooted: bool,
    lambda_checks: Vec<LambdaCheck>,
    supported: bool,
    unsupported_reason: Option<String>,
}

impl PfPencilReport {
    fn finite_evidence_ok(&self) -> bool {
        self.supported
            && self.j_alpha_real_rooted
            && self.t_j_beta_real_rooted
            && self.lambda_checks.iter().all(|check| check.real_rooted)
    }
}

fn pf_pencil_case(id: &str) -> Option<PfPencilCase> {
    PF_PENCIL_CASES
        .iter()
        .copied()
        .find(|case| case.id.eq_ignore_ascii_case(id))
}

fn lcm_usize(a: usize, b: usize) -> usize {
    use num_integer::Integer;
    a.lcm(&b)
}

fn pf_value(case_id: &str, which: &str, k: usize, d: usize) -> (BigInt, usize) {
    let k = k as i128;
    let d_i = d as i128;
    let value = match (case_id, which) {
        ("A036969", "alpha") => ((k + 1) * (k + 1), 1),
        ("A036969", "beta") => (1, 1),
        ("A071951", "alpha") => ((k + 1) * (k + 2), 1),
        ("A071951", "beta") => (1, 1),
        ("A080248", "alpha") => ((k + 1) * (k + 2), 2),
        ("A080248", "beta") => (1, 1),
        ("A156289", "alpha") => ((k + 1) * (k + 1), 1),
        ("A156289", "beta") => (2 * k + 3, 1),
        ("A160562", "alpha") => ((2 * k + 1) * (2 * k + 1), 1),
        ("A160562", "beta") => (1, 1),
        ("A269945", "alpha") => (k * k, 1),
        ("A269945", "beta") => (1, 1),
        ("A166960", "alpha") => ((k + 1) * (k + 1), 1),
        ("A166960", "beta") => (d_i + 3 - k, 1),
        ("A166961", "alpha") => ((k + 1) * (2 * k + 1), 1),
        ("A166961", "beta") => (2 * d_i + 5 - 2 * k, 1),
        ("A166962", "alpha") => ((k + 1) * (3 * k + 1), 1),
        ("A166962", "beta") => (3 * d_i + 1 - 3 * k, 1),
        ("A166972", "alpha") => ((k + 1) * (3 * k + 1), 1),
        ("A166972", "beta") => (d_i + 1 - k, 1),
        ("A191935", "alpha") => (1, 1),
        ("A191935", "beta") => ((d_i + 2 - k) * (d_i + 3 - k), 1),
        ("A371081", "alpha") => ((d_i + k + 2) * (d_i + k + 2), 1),
        ("A371081", "beta") => (1, 1),
        ("A371259", "alpha") => ((d_i + k + 4) * (d_i + k + 4), 1),
        ("A371259", "beta") => (1, 1),
        ("A390433", "alpha") => ((d_i + k) * (d_i + k), 1),
        ("A390433", "beta") => (1, 1),
        ("A198204", "alpha") => ((k + 1) * (k + d_i), d),
        ("A198204", "beta") => ((k + d_i + 1) * (k + d_i), d),
        _ => unreachable!("unknown PF pencil case/action"),
    };
    (BigInt::from(value.0), value.1)
}

fn pf_common_denominator(case_id: &str, degree: usize) -> usize {
    let mut denom = 1;
    for k in 0..=degree {
        denom = lcm_usize(denom, pf_value(case_id, "alpha", k, degree).1);
        denom = lcm_usize(denom, pf_value(case_id, "beta", k, degree).1);
    }
    denom
}

fn pf_kernel(case: PfPencilCase, degree: usize, which: &str, common_den: usize) -> Vec<BigInt> {
    (0..=degree)
        .map(|k| {
            let (num, denom) = pf_value(case.id, which, k, degree);
            binomial_bigint(degree, k) * num * BigInt::from(common_den / denom)
        })
        .collect()
}

fn shift_by_t(coeffs: &[BigInt]) -> Vec<BigInt> {
    let mut out = Vec::with_capacity(coeffs.len() + 1);
    out.push(BigInt::from(0));
    out.extend_from_slice(coeffs);
    out
}

fn add_lambda_pencil(j_alpha: &[BigInt], t_j_beta: &[BigInt], lambda: &BigInt) -> Vec<BigInt> {
    let len = j_alpha.len().max(t_j_beta.len());
    (0..len)
        .map(|i| {
            let mut coeff = j_alpha.get(i).cloned().unwrap_or_else(|| BigInt::from(0));
            if let Some(beta_coeff) = t_j_beta.get(i) {
                coeff += lambda * beta_coeff;
            }
            coeff
        })
        .collect()
}

fn pf_pencil_report(case: PfPencilCase, degree: usize, lambdas: &[BigInt]) -> PfPencilReport {
    if case.id == "A198204" && degree == 0 {
        return PfPencilReport {
            case,
            degree,
            common_denominator: BigInt::from(0),
            j_alpha: Vec::new(),
            t_j_beta: Vec::new(),
            j_alpha_real_rooted: false,
            t_j_beta_real_rooted: false,
            lambda_checks: Vec::new(),
            supported: false,
            unsupported_reason: Some(
                "A198204 has denominator d and is unsupported for d=0".to_string(),
            ),
        };
    }
    let common_den = pf_common_denominator(case.id, degree);
    let j_alpha = pf_kernel(case, degree, "alpha", common_den);
    let beta = pf_kernel(case, degree, "beta", common_den);
    let t_j_beta = shift_by_t(&beta);
    let lambda_checks = lambdas
        .iter()
        .map(|lambda| {
            let coefficients = add_lambda_pencil(&j_alpha, &t_j_beta, lambda);
            let real_rooted = is_real_rooted_bigint_coeffs(&coefficients);
            LambdaCheck {
                lambda: lambda.clone(),
                coefficients,
                real_rooted,
            }
        })
        .collect();
    PfPencilReport {
        case,
        degree,
        common_denominator: BigInt::from(common_den),
        j_alpha_real_rooted: is_real_rooted_bigint_coeffs(&j_alpha),
        t_j_beta_real_rooted: is_real_rooted_bigint_coeffs(&t_j_beta),
        j_alpha,
        t_j_beta,
        lambda_checks,
        supported: true,
        unsupported_reason: None,
    }
}

fn pf_pencil_report_json(report: &PfPencilReport) -> String {
    let lambda_checks = report
        .lambda_checks
        .iter()
        .map(|check| {
            format!(
                "{{\"lambda\":{},\"coefficients\":{},\"polynomial\":{},\"real_rooted\":{}}}",
                json_string(&check.lambda.to_string()),
                json_bigint_vec(&check.coefficients),
                json_string(&format_poly_bigint_coeffs(&check.coefficients)),
                check.real_rooted
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"case\":{},\"degree\":{},\"alpha_label\":{},\"beta_label\":{},\
         \"common_denominator\":{},\"supported\":{},\"unsupported_reason\":{},\
         \"j_alpha_coefficients\":{},\"t_j_beta_coefficients\":{},\
         \"j_alpha_polynomial\":{},\"t_j_beta_polynomial\":{},\
         \"j_alpha_real_rooted\":{},\"t_j_beta_real_rooted\":{},\
         \"lambda_checks\":[{}],\"finite_evidence_ok\":{}}}",
        json_string(report.case.id),
        report.degree,
        json_string(report.case.alpha_label),
        json_string(report.case.beta_label),
        json_string(&report.common_denominator.to_string()),
        report.supported,
        report
            .unsupported_reason
            .as_ref()
            .map(|reason| json_string(reason))
            .unwrap_or_else(|| "null".to_string()),
        json_bigint_vec(&report.j_alpha),
        json_bigint_vec(&report.t_j_beta),
        json_string(&format_poly_bigint_coeffs(&report.j_alpha)),
        json_string(&format_poly_bigint_coeffs(&report.t_j_beta)),
        report.j_alpha_real_rooted,
        report.t_j_beta_real_rooted,
        lambda_checks,
        report.finite_evidence_ok()
    )
}

fn print_pf_pencil_text(reports: &[PfPencilReport]) {
    for report in reports {
        println!("{} d={}", report.case.id, report.degree);
        println!("  alpha(k,d) = {}", report.case.alpha_label);
        println!("  beta(k,d)  = {}", report.case.beta_label);
        if !report.supported {
            println!(
                "  unsupported: {}",
                report
                    .unsupported_reason
                    .as_deref()
                    .unwrap_or("unsupported case")
            );
            continue;
        }
        println!(
            "  common denominator cleared: {}",
            report.common_denominator
        );
        println!(
            "  J_alpha,d(t) = {} [{}]",
            format_poly_bigint_coeffs(&report.j_alpha),
            if report.j_alpha_real_rooted {
                "real-rooted"
            } else {
                "NOT real-rooted"
            }
        );
        println!(
            "  t J_beta,d(t) = {} [{}]",
            format_poly_bigint_coeffs(&report.t_j_beta),
            if report.t_j_beta_real_rooted {
                "real-rooted"
            } else {
                "NOT real-rooted"
            }
        );
        for check in &report.lambda_checks {
            println!(
                "  lambda={}: {} [{}]",
                check.lambda,
                format_poly_bigint_coeffs(&check.coefficients),
                if check.real_rooted {
                    "real-rooted"
                } else {
                    "NOT real-rooted"
                }
            );
        }
        println!("  finite evidence ok: {}", report.finite_evidence_ok());
    }
}

fn parse_lambdas(input: &str) -> Result<Vec<BigInt>, String> {
    let mut lambdas = Vec::new();
    for raw in input.split(',') {
        let value = raw.trim();
        if value.is_empty() {
            return Err("empty lambda in --lambdas list".to_string());
        }
        if value.starts_with('-') {
            return Err(format!(
                "--lambdas expects nonnegative integers, got '{value}'"
            ));
        }
        let lambda = value
            .parse::<BigInt>()
            .map_err(|_| format!("--lambdas expects nonnegative integers, got '{value}'"))?;
        lambdas.push(lambda);
    }
    if lambdas.is_empty() {
        return Err("--lambdas expects at least one value".to_string());
    }
    Ok(lambdas)
}

pub(crate) fn cmd_pf_pencil(args: &[String]) {
    let mut case_id: Option<String> = None;
    let mut all_family_h = false;
    let mut degree: Option<usize> = None;
    let mut min_degree: Option<usize> = None;
    let mut max_degree: Option<usize> = None;
    let mut lambdas = parse_lambdas("0,1,10,100").expect("default lambdas parse");
    let mut format = OutputFormat::Text;

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--case" => {
                case_id = match next_option_value(args, &mut i, "--case") {
                    Ok(value) => Some(value.to_string()),
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--all-family-h" => all_family_h = true,
            "--degree" => {
                degree = match parse_usize_option(args, &mut i, "--degree") {
                    Ok(value) => Some(value),
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--min-degree" => {
                min_degree = match parse_usize_option(args, &mut i, "--min-degree") {
                    Ok(value) => Some(value),
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--max-degree" => {
                max_degree = match parse_usize_option(args, &mut i, "--max-degree") {
                    Ok(value) => Some(value),
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--lambdas" => {
                lambdas = match next_option_value(args, &mut i, "--lambdas").and_then(parse_lambdas)
                {
                    Ok(value) => value,
                    Err(error) => {
                        eprintln!("{error}");
                        return;
                    }
                };
            }
            "--json" => format = OutputFormat::Json,
            other => {
                eprintln!("unknown option: {other}");
                print_pf_pencil_help();
                std::process::exit(1);
            }
        }
        i += 1;
    }

    if all_family_h == case_id.is_some() {
        eprintln!("choose exactly one of --case or --all-family-h");
        print_pf_pencil_help();
        std::process::exit(1);
    }

    let degrees = if let Some(d) = degree {
        if min_degree.is_some() || max_degree.is_some() {
            eprintln!("--degree cannot be combined with --min-degree or --max-degree");
            std::process::exit(1);
        }
        vec![d]
    } else {
        let min = min_degree.unwrap_or(1);
        let Some(max) = max_degree else {
            eprintln!("expected --degree or --max-degree");
            print_pf_pencil_help();
            std::process::exit(1);
        };
        if min > max {
            eprintln!("--min-degree must be <= --max-degree");
            std::process::exit(1);
        }
        (min..=max).collect()
    };

    let cases = if all_family_h {
        PF_PENCIL_CASES.to_vec()
    } else {
        let id = case_id.expect("case_id checked");
        let Some(case) = pf_pencil_case(&id) else {
            eprintln!("unknown PF pencil case: {id}");
            print_pf_pencil_help();
            std::process::exit(1);
        };
        vec![case]
    };

    let reports = cases
        .into_iter()
        .flat_map(|case| {
            degrees.iter().copied().map({
                let lambdas = lambdas.clone();
                move |degree| pf_pencil_report(case, degree, &lambdas)
            })
        })
        .collect::<Vec<_>>();

    if format == OutputFormat::Json {
        let items = reports
            .iter()
            .map(pf_pencil_report_json)
            .collect::<Vec<_>>()
            .join(",");
        let overall = reports.iter().all(PfPencilReport::finite_evidence_ok);
        println!(
            "{{\"schema\":\"polytool.pf-pencil.v1\",\"items\":[{}],\
             \"overall_finite_evidence_ok\":{}}}",
            items, overall
        );
    } else {
        print_pf_pencil_text(&reports);
    }
}
