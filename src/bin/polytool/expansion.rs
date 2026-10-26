//! `polytool gamma-expansion` and `polytool bernstein-expansion`.

use crate::input::{
    parse_usize_option, polynomial_degree_bigint, read_polys_bigint, read_polys_rational,
    strip_trailing_zeros_bigint, OutputFormat,
};
use crate::json_out::{json_bigint_vec, json_string};
use crate::recurrence::format_rational_row;
use num_bigint::BigInt;
use polytool::*;
use serde_json::json;

pub(crate) fn cmd_gamma_expansion(args: &[String]) {
    let format = match OutputFormat::from_args(args) {
        Ok(format) => format,
        Err(error) => {
            eprintln!("{error}");
            return;
        }
    };
    let mut had_error = false;
    let mut json_items = Vec::new();
    for (index, coeffs) in read_polys_bigint().into_iter().enumerate() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        match gamma_coefficients_bigint_coeffs(c) {
            Some(gamma) => {
                let polynomial = format_poly_bigint_coeffs(c);
                let expansion = format_gamma_expansion(&gamma, c.len().saturating_sub(1));
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":true,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"gamma\":{},\"expansion\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_bigint_vec(c),
                        polynomial_degree_bigint(c),
                        json_bigint_vec(&gamma),
                        json_string(&expansion)
                    ));
                } else {
                    println!("{polynomial}: gamma {:?}; expansion: {expansion}", gamma);
                }
            }
            None => {
                let polynomial = format_poly_bigint_coeffs(c);
                let error = "not palindromic; no gamma expansion";
                if format == OutputFormat::Json {
                    json_items.push(format!(
                        "{{\"index\":{},\"ok\":false,\"polynomial\":{},\"coefficients\":{},\
                         \"degree\":{},\"error\":{}}}",
                        index,
                        json_string(&polynomial),
                        json_bigint_vec(c),
                        polynomial_degree_bigint(c),
                        json_string(error)
                    ));
                } else {
                    eprintln!("{polynomial}: {error}");
                }
                had_error = true;
            }
        }
    }
    if format == OutputFormat::Json {
        println!(
            "{{\"ok\":{},\"items\":[{}]}}",
            !had_error,
            json_items.join(",")
        );
    }
    if had_error {
        std::process::exit(1);
    }
}

pub(crate) fn cmd_bernstein_expansion(args: &[String]) {
    let mut format = OutputFormat::Text;
    let mut ambient_degree = None;
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--json" => format = OutputFormat::Json,
            "--degree" => {
                ambient_degree = Some(
                    parse_usize_option(args, &mut index, "--degree").unwrap_or_else(|error| {
                        eprintln!("{error}");
                        std::process::exit(2);
                    }),
                );
            }
            other => {
                eprintln!("unknown option: {other}");
                std::process::exit(2);
            }
        }
        index += 1;
    }

    let mut items = Vec::new();
    let mut had_error = false;
    for (row_index, coefficients) in read_polys_rational().into_iter().enumerate() {
        let polynomial = Polynomial::new(coefficients.clone());
        let degree = ambient_degree.unwrap_or_else(|| polynomial.degree().unwrap_or(0));
        match bernstein_basis_coordinates(&polynomial, degree) {
            Ok(coordinates) => {
                if format == OutputFormat::Json {
                    items.push(json!({
                        "index": row_index,
                        "ok": true,
                        "coefficients": coefficients
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                        "degree": degree,
                        "coordinates": coordinates
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                    }));
                } else {
                    println!("{}", format_rational_row(&coordinates));
                }
            }
            Err(error) => {
                had_error = true;
                if format == OutputFormat::Json {
                    items.push(json!({
                        "index": row_index,
                        "ok": false,
                        "coefficients": coefficients
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>(),
                        "degree": degree,
                        "error": error.to_string(),
                    }));
                } else {
                    eprintln!("row {row_index}: {error}");
                }
            }
        }
    }
    if format == OutputFormat::Json {
        println!(
            "{}",
            serde_json::to_string_pretty(&json!({
                "ok": !had_error,
                "basis": "binomial(n,j) t^j (1-t)^(n-j)",
                "items": items,
            }))
            .unwrap()
        );
    }
    if had_error {
        std::process::exit(1);
    }
}

fn format_gamma_expansion(gamma: &[BigInt], degree: usize) -> String {
    let mut out = String::new();
    let zero = BigInt::from(0);
    let one = BigInt::from(1);
    for (i, coeff) in gamma.iter().enumerate() {
        if coeff == &zero {
            continue;
        }
        let negative = coeff < &zero;
        let abs_coeff = if negative {
            -coeff.clone()
        } else {
            coeff.clone()
        };
        let factor = gamma_basis_factor(i, degree.saturating_sub(2 * i));
        let body = if factor == "1" {
            abs_coeff.to_string()
        } else if abs_coeff == one {
            factor
        } else {
            format!("{abs_coeff} {factor}")
        };

        if out.is_empty() {
            if negative {
                out.push('-');
            }
            out.push_str(&body);
        } else if negative {
            out.push_str(" - ");
            out.push_str(&body);
        } else {
            out.push_str(" + ");
            out.push_str(&body);
        }
    }
    if out.is_empty() {
        "0".to_string()
    } else {
        out
    }
}

fn gamma_basis_factor(t_power: usize, one_plus_t_power: usize) -> String {
    let mut factors = Vec::new();
    match t_power {
        0 => {}
        1 => factors.push("t".to_string()),
        n => factors.push(format!("t^{n}")),
    }
    match one_plus_t_power {
        0 => {}
        1 => factors.push("(1+t)".to_string()),
        n => factors.push(format!("(1+t)^{n}")),
    }
    if factors.is_empty() {
        "1".to_string()
    } else {
        factors.join(" ")
    }
}
