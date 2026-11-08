//! `polytool resultant` and `polytool discriminant`.

use crate::input::{read_polys_bigint, strip_trailing_zeros_bigint, OutputFormat};
use polytool::*;

pub(crate) fn cmd_resultant(args: &[String]) {
    let format = output_format(args);
    let polys = read_polys_bigint();
    if polys.len() != 2 {
        eprintln!("Need exactly two polynomials for resultant.");
        std::process::exit(2);
    }
    let p = strip_trailing_zeros_bigint(&polys[0]);
    let q = strip_trailing_zeros_bigint(&polys[1]);
    let r = resultant_bigint_coeffs(p, q);
    if format == OutputFormat::Json {
        println!("{}", serde_json::json!({"resultant": r.to_string()}));
        return;
    }
    println!(
        "Res({}, {}) = {}",
        format_poly_bigint_coeffs(p),
        format_poly_bigint_coeffs(q),
        r
    );
}

pub(crate) fn cmd_discriminant(args: &[String]) {
    let format = output_format(args);
    let mut items = Vec::new();
    for (index, coeffs) in read_polys_bigint().iter().enumerate() {
        let c = strip_trailing_zeros_bigint(coeffs);
        let d = discriminant_bigint_coeffs(c);
        if format == OutputFormat::Json {
            items.push(serde_json::json!({"index": index, "discriminant": d.to_string()}));
        } else {
            println!("disc({}) = {}", format_poly_bigint_coeffs(c), d);
        }
    }
    if format == OutputFormat::Json {
        println!("{}", serde_json::json!({"items": items}));
    }
}

fn output_format(args: &[String]) -> OutputFormat {
    OutputFormat::from_args(args).unwrap_or_else(|error| {
        eprintln!("{error}");
        std::process::exit(2);
    })
}
