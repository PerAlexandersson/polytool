//! Shared input handling: bounded stdin reads, coefficient-list parsing,
//! option-value parsing and the common `--json` output-format flag.

use num_bigint::BigInt;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::*;
use std::io::{self, Read};

const MAX_CLI_INPUT_BYTES: usize = 16 * 1024 * 1024;

pub(crate) fn read_limited_text(reader: impl Read, description: &str) -> Result<String, String> {
    let mut input = String::new();
    let mut limited = reader.take((MAX_CLI_INPUT_BYTES as u64) + 1);
    limited
        .read_to_string(&mut input)
        .map_err(|error| format!("failed to read {description}: {error}"))?;
    if input.len() > MAX_CLI_INPUT_BYTES {
        return Err(format!(
            "{description} exceeds the {MAX_CLI_INPUT_BYTES}-byte CLI input limit"
        ));
    }
    Ok(input)
}

pub(crate) fn read_stdin_or_exit() -> String {
    match read_limited_text(io::stdin(), "stdin") {
        Ok(input) => input,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

pub(crate) fn is_help_arg(arg: &str) -> bool {
    matches!(arg, "-h" | "--help" | "help")
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum OutputFormat {
    Text,
    Json,
}

impl OutputFormat {
    pub(crate) fn from_args(args: &[String]) -> Result<Self, String> {
        let mut format = Self::Text;
        for arg in args {
            match arg.as_str() {
                "--json" => format = Self::Json,
                other => return Err(format!("unknown option: {other}")),
            }
        }
        Ok(format)
    }
}

pub(crate) fn next_option_value<'a>(
    args: &'a [String],
    index: &mut usize,
    option: &str,
) -> Result<&'a str, String> {
    *index += 1;
    args.get(*index)
        .map(String::as_str)
        .ok_or_else(|| format!("{option} expects a value"))
}

pub(crate) fn parse_usize_option(
    args: &[String],
    index: &mut usize,
    option: &str,
) -> Result<usize, String> {
    let value = next_option_value(args, index, option)?;
    value
        .parse()
        .map_err(|_| format!("{option} expects a nonnegative integer, got '{value}'"))
}

pub(crate) fn read_polys_bigint() -> Vec<Vec<BigInt>> {
    let input = read_stdin_or_exit();
    parse_polynomials_bigint(&input)
        .into_iter()
        .map(|result| match result {
            Ok(polynomial) => polynomial,
            Err(error) => {
                eprintln!("Invalid polynomial input: {error}");
                std::process::exit(2);
            }
        })
        .collect()
}

pub(crate) fn read_poly_parse_results_bigint() -> Vec<Result<Vec<BigInt>, String>> {
    let input = read_stdin_or_exit();
    parse_polynomials_bigint(&input)
        .into_iter()
        .map(|result| result.map_err(|error| error.to_string()))
        .collect()
}

pub(crate) fn read_polys_rational() -> Vec<Vec<RecurrenceBigRational>> {
    let input = read_stdin_or_exit();
    input
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .map(|line| match parse_coeff_list_rational(line) {
            Ok(polynomial) => polynomial,
            Err(error) => {
                eprintln!("Invalid polynomial input: {error}");
                std::process::exit(2);
            }
        })
        .collect()
}

fn parse_coeff_list_rational(input: &str) -> Result<Vec<RecurrenceBigRational>, String> {
    let s = input.trim();
    if s.is_empty() {
        return Err("empty input".to_string());
    }
    if s.chars().any(|c| c.is_ascii_alphabetic()) {
        return Err(format!(
            "recurrence input expects coefficient lists, got expanded polynomial '{}'",
            input
        ));
    }

    let s = strip_coeff_list_brackets(s);
    let parts: Vec<&str> = if s.contains(',') {
        s.split(',').collect()
    } else {
        s.split_whitespace().collect()
    };
    if parts.len() > MAX_POLYNOMIAL_COEFFICIENTS {
        return Err(format!(
            "polynomial has {} coefficients; the limit is {MAX_POLYNOMIAL_COEFFICIENTS}",
            parts.len()
        ));
    }

    parts
        .iter()
        .map(|part| parse_big_rational(part.trim()))
        .collect()
}

fn parse_big_rational(token: &str) -> Result<RecurrenceBigRational, String> {
    if token.is_empty() {
        return Err("empty coefficient".to_string());
    }
    if let Some((num, den)) = token.split_once('/') {
        let num = parse_bigint(num.trim())?;
        let den = parse_bigint(den.trim())?;
        if den == BigInt::from(0) {
            return Err(format!("zero denominator in coefficient '{}'", token));
        }
        Ok(RecurrenceBigRational::new(num, den))
    } else {
        Ok(RecurrenceBigRational::from_integer(parse_bigint(token)?))
    }
}

fn parse_bigint(token: &str) -> Result<BigInt, String> {
    token
        .parse::<BigInt>()
        .map_err(|e| format!("invalid integer '{}': {}", token, e))
}

fn strip_coeff_list_brackets(s: &str) -> &str {
    let s = s.trim();
    if (s.starts_with('[') && s.ends_with(']'))
        || (s.starts_with('{') && s.ends_with('}'))
        || (s.starts_with('(') && s.ends_with(')'))
    {
        &s[1..s.len() - 1]
    } else {
        s
    }
}

pub(crate) fn strip_trailing_zeros_bigint(coeffs: &[BigInt]) -> &[BigInt] {
    let end = coeffs
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .map_or(0, |i| i + 1);
    &coeffs[..end]
}

pub(crate) fn polynomial_degree(coeffs: &[i64]) -> usize {
    coeffs.iter().rposition(|&c| c != 0).unwrap_or(0)
}

pub(crate) fn polynomial_degree_bigint(coeffs: &[BigInt]) -> usize {
    coeffs
        .iter()
        .rposition(|c| c != &BigInt::from(0))
        .unwrap_or(0)
}
