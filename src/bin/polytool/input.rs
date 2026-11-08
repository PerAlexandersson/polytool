//! Shared input handling: bounded stdin reads (or the `--input <file|->`
//! replacement), coefficient-list parsing, option-value parsing and the
//! common `--json` output-format flag.

use num_bigint::BigInt;
use polytool::recurrence::BigRational as RecurrenceBigRational;
use polytool::*;
use std::io::{self, Read};
use std::sync::OnceLock;

const MAX_CLI_INPUT_BYTES: usize = 16 * 1024 * 1024;

/// The file named by `--input`, read instead of standard input.
static INPUT_PATH: OnceLock<String> = OnceLock::new();

/// Read the command input from `path` instead of stdin (`-` keeps stdin).
pub(crate) fn set_input_path(path: &str) {
    if path != "-" {
        let _ = INPUT_PATH.set(path.to_string());
    }
}

/// Read the whole command input (stdin, or the `--input` file) within the
/// CLI size limit; `what` names the content in error messages.
pub(crate) fn read_command_input(what: &str) -> Result<String, String> {
    match INPUT_PATH.get() {
        Some(path) => {
            let file = std::fs::File::open(path)
                .map_err(|error| format!("failed to open input file `{path}`: {error}"))?;
            read_limited_text(file, &format!("{what} `{path}`"))
        }
        None => read_limited_text(io::stdin(), what),
    }
}

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
    match read_command_input("stdin") {
        Ok(input) => input,
        Err(error) => {
            eprintln!("{error}");
            std::process::exit(2);
        }
    }
}

pub(crate) fn is_help_arg(arg: &str) -> bool {
    matches!(arg, "-h" | "--help")
}

/// Reject a command's arguments without writing help text into its data stream.
pub(crate) fn usage_error(command: &str, message: &str) -> ! {
    eprintln!("{message}");
    if command.is_empty() {
        eprintln!("Run `polytool --help` for usage.");
    } else {
        eprintln!("Run `polytool help {command}` for usage.");
    }
    std::process::exit(2);
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

pub(crate) fn parse_u64_option(
    args: &[String],
    index: &mut usize,
    option: &str,
) -> Result<u64, String> {
    let value = next_option_value(args, index, option)?;
    value
        .parse()
        .map_err(|_| format!("{option} expects a nonnegative integer, got '{value}'"))
}

pub(crate) fn read_polys_bigint() -> Vec<Vec<BigInt>> {
    read_poly_parse_results_bigint()
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

/// The `Display` text of an error, for the `Result<_, String>` plumbing of the
/// subcommands.
pub(crate) fn error_text<E: std::fmt::Display>(error: E) -> String {
    error.to_string()
}

pub(crate) fn read_poly_parse_results_bigint() -> Vec<Result<Vec<BigInt>, String>> {
    let input = read_stdin_or_exit();
    polynomial_input_lines(&input)
        .map(|(number, line)| {
            parse_polynomial_bigint(line).map_err(|error| format!("line {number}: {error}"))
        })
        .collect()
}

/// Preserve physical source lines while ignoring comments and blank lines.
fn polynomial_input_lines(input: &str) -> impl Iterator<Item = (usize, &str)> {
    input.lines().enumerate().filter_map(|(index, line)| {
        let trimmed = line.trim();
        (!trimmed.is_empty() && !trimmed.starts_with('#')).then_some((index + 1, line))
    })
}

pub(crate) fn read_polys_rational() -> Vec<Vec<RecurrenceBigRational>> {
    let input = read_stdin_or_exit();
    polynomial_input_lines(&input)
        .map(|(number, line)| match parse_coeff_list_rational(line) {
            Ok(polynomial) => polynomial,
            Err(error) => {
                eprintln!("Invalid polynomial input at line {number}: {error}");
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
            "expected a coefficient list, got expanded polynomial '{}'",
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
        .enumerate()
        .map(|(index, part)| {
            polytool::recurrence::parse_rational_coeff(part)
                .map_err(|error| format!("coefficient {}: {error}", index + 1))
        })
        .collect()
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_lists_keep_exact_values_and_existing_syntax() {
        let expected = vec![
            RecurrenceBigRational::new(1.into(), 2.into()),
            RecurrenceBigRational::new((-3).into(), 4.into()),
        ];
        for input in ["1/2 -3/4", "[ 2 / 4, 3 / -4 ]", "{1/2,-3/4}", "(1/2,-3/4)"] {
            assert_eq!(
                parse_coeff_list_rational(input).unwrap(),
                expected,
                "{input}"
            );
        }
        assert!(parse_coeff_list_rational("[]").unwrap().is_empty());
        let huge = "1234567890123456789012345678901234567890";
        assert_eq!(
            parse_coeff_list_rational(huge).unwrap()[0].to_string(),
            huge
        );
        for input in ["1,,2", "1,2/0", "1,2/3/4", "1+x", "[1,2)"] {
            assert!(parse_coeff_list_rational(input).is_err(), "{input}");
        }
    }
}
