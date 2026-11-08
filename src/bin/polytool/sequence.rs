//! `polytool sequence`: standard polynomial sequences.

use crate::input::{next_option_value, polynomial_degree_bigint, strip_trailing_zeros_bigint};
use crate::json_out::{json_bigint_vec, json_string};
use num_bigint::BigInt;
use polytool::*;

pub(crate) fn print_sequence_help() {
    println!("Usage:");
    println!("  polytool sequence <name> <max-n> [--format text|rows|json]");
    println!("  polytool sequence --help");
    println!();
    println!("Generate standard polynomial sequences.");
    println!();
    println!("Names:");
    println!("  eulerian");
    println!("  narayana");
    println!("  type-b-eulerian  (aliases: type-b, type_b_eulerian)");
    println!("  chebyshev-t      (aliases: chebyshev_t, t)");
    println!("  chebyshev-u      (aliases: chebyshev_u, u)");
    println!("  hermite");
    println!();
    println!("Options:");
    println!("  --format <name>  text (default), rows (coefficient lists), or json");
    println!("  --json           Alias for --format json");
    println!("  -h, --help       Print this help text");
    println!();
    println!("Note: generated coefficients use arbitrary-size integers.");
    println!();
    println!("Example:");
    println!("  polytool sequence eulerian 5");
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum SequenceKind {
    Eulerian,
    Narayana,
    TypeBEulerian,
    ChebyshevT,
    ChebyshevU,
    Hermite,
}

impl SequenceKind {
    pub(crate) fn parse(input: &str) -> Option<Self> {
        match input {
            "eulerian" => Some(Self::Eulerian),
            "narayana" => Some(Self::Narayana),
            "type-b-eulerian" | "type-b" | "type_b_eulerian" | "type_b" => {
                Some(Self::TypeBEulerian)
            }
            "chebyshev-t" | "chebyshev_t" | "t" => Some(Self::ChebyshevT),
            "chebyshev-u" | "chebyshev_u" | "u" => Some(Self::ChebyshevU),
            "hermite" => Some(Self::Hermite),
            _ => None,
        }
    }

    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::Eulerian => "eulerian",
            Self::Narayana => "narayana",
            Self::TypeBEulerian => "type_b_eulerian",
            Self::ChebyshevT => "chebyshev_t",
            Self::ChebyshevU => "chebyshev_u",
            Self::Hermite => "hermite",
        }
    }

    pub(crate) fn polynomials(self, max_n: usize) -> Vec<Vec<BigInt>> {
        match self {
            // The library families start at A_0 = N_0 = 1 (equal to A_1, N_1);
            // the CLI keeps the classical labels A_1..A_max_n and N_1..N_max_n.
            Self::Eulerian => polytool::sequences::eulerian_polynomials_bigint(max_n)
                .into_iter()
                .skip(1)
                .collect(),
            Self::Narayana => polytool::sequences::narayana_polynomials_bigint(max_n)
                .into_iter()
                .skip(1)
                .collect(),
            Self::TypeBEulerian => polytool::sequences::type_b_eulerian_polynomials_bigint(max_n),
            Self::ChebyshevT => polytool::sequences::chebyshev_polynomials_t_bigint(max_n),
            Self::ChebyshevU => polytool::sequences::chebyshev_polynomials_u_bigint(max_n),
            Self::Hermite => polytool::sequences::hermite_polynomials_bigint(max_n),
        }
    }

    fn label(self, index: usize) -> String {
        match self {
            Self::Eulerian => format!("A_{}(t)", index + 1),
            Self::Narayana => format!("N_{}(t)", index + 1),
            Self::TypeBEulerian => format!("B_{index}(t)"),
            Self::ChebyshevT => format!("T_{index}(t)"),
            Self::ChebyshevU => format!("U_{index}(t)"),
            Self::Hermite => format!("He_{index}(t)"),
        }
    }
}

pub(crate) fn cmd_sequence(args: &[String]) {
    let mut positional = Vec::new();
    let mut format = "text";
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--json" => format = "json",
            "--format" => {
                format = next_option_value(args, &mut i, "--format")
                    .unwrap_or_else(|error| fail(&error));
                if !matches!(format, "text" | "rows" | "json") {
                    fail(&format!("unknown sequence output format: {format}"));
                }
            }
            other if other.starts_with('-') => fail(&format!("unknown sequence option: {other}")),
            other => positional.push(other),
        }
        i += 1;
    }
    if positional.len() != 2 {
        crate::input::usage_error("sequence", "expected a sequence name and max-n");
    }
    let Some(kind) = SequenceKind::parse(positional[0]) else {
        crate::input::usage_error("sequence", &format!("unknown sequence: {}", positional[0]));
    };
    let max_n = positional[1].parse::<usize>().unwrap_or_else(|_| {
        eprintln!(
            "expected a nonnegative integer max-n, got '{}'",
            positional[1]
        );
        std::process::exit(2);
    });
    let polynomials = kind.polynomials(max_n);
    if format == "json" {
        let items = polynomials
            .iter()
            .enumerate()
            .map(|(index, coeffs)| {
                let c = strip_trailing_zeros_bigint(coeffs);
                format!(
                    "{{\"index\":{},\"label\":{},\"polynomial\":{},\"coefficients\":{},\
                     \"degree\":{}}}",
                    index,
                    json_string(&kind.label(index)),
                    json_string(&format_poly_bigint_coeffs(c)),
                    json_bigint_vec(c),
                    polynomial_degree_bigint(c)
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        println!(
            "{{\"sequence\":{},\"max_n\":{},\"polynomials\":[{}]}}",
            json_string(kind.name()),
            max_n,
            items
        );
        return;
    }
    for (index, coeffs) in polynomials.iter().enumerate() {
        let c = strip_trailing_zeros_bigint(coeffs);
        if format == "rows" {
            if c.is_empty() {
                println!("0");
            } else {
                println!(
                    "{}",
                    c.iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                );
            }
        } else {
            println!("{} = {}", kind.label(index), format_poly_bigint_coeffs(c));
        }
    }
}

fn fail(message: &str) -> ! {
    eprintln!("{message}");
    std::process::exit(2);
}
