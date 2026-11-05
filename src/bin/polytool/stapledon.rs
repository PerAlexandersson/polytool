//! `polytool stapledon`: Stapledon decompositions.

use crate::help::print_coefficient_input_help;
use crate::input::{is_help_arg, read_polys_bigint, strip_trailing_zeros_bigint};
use polytool::*;

pub(crate) fn print_stapledon_help() {
    println!("Usage:");
    println!("  polytool stapledon <n>");
    println!("  polytool stapledon --help");
    println!();
    println!("Compute the Stapledon decomposition with respect to degree bound n.");
    println!();
    print_coefficient_input_help();
    println!();
    println!("Note: accepts arbitrary-size integer coefficients.");
}

pub(crate) fn cmd_stapledon(args: &[String]) {
    if args.iter().any(|arg| is_help_arg(arg)) {
        print_stapledon_help();
        return;
    }

    if args.len() != 1 {
        eprintln!("stapledon expects one degree bound");
        std::process::exit(2);
    }

    let n: usize = match args[0].parse() {
        Ok(n) => n,
        Err(_) => {
            eprintln!("Expected a nonnegative integer degree bound.");
            std::process::exit(2);
        }
    };

    for coeffs in read_polys_bigint() {
        let c = strip_trailing_zeros_bigint(&coeffs);
        match stapledon_decomposition_bigint_coeffs(c, n) {
            Some((a, b)) => {
                println!(
                    "{} = {} + x ({})",
                    format_poly_bigint_coeffs(c),
                    format_poly_bigint_coeffs(&a),
                    format_poly_bigint_coeffs(&b),
                );
            }
            None => {
                eprintln!(
                    "{} has degree greater than the requested bound {}.",
                    format_poly_bigint_coeffs(c),
                    n,
                );
            }
        }
    }
}
