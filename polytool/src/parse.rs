//! Flexible polynomial parsing from human-readable text.
//!
//! Accepts many common formats:
//!
//! - Comma-separated coefficients: `1, 2, 3` → 1 + 2t + 3t²
//! - Space-separated coefficients: `1 2 3`
//! - Bracketed: `[1, 2, 3]` or `{1, 2, 3}` or `(1, 2, 3)`
//! - Expanded polynomial: `1 + 2t + 3t^2` (any single letter as variable,
//!   `*` optional for multiplication, `^` for exponent, `**` also accepted)
//!
//! All formats produce coefficients in ascending degree order.

use num_bigint::BigInt;
use num_traits::ToPrimitive;

/// Largest textual polynomial accepted by the general parser.
pub const MAX_POLYNOMIAL_INPUT_BYTES: usize = 1_048_576;

/// Largest dense coefficient vector produced by the general parser.
///
/// Polytool stores polynomials densely, so accepting a larger exponent would
/// allocate every intervening zero coefficient as well.
pub const MAX_POLYNOMIAL_COEFFICIENTS: usize = 100_001;

/// Parse a single line into polynomial coefficients (ascending degree).
///
/// Tries expanded polynomial format first (if a letter is found),
/// then falls back to coefficient list parsing.
pub fn parse_polynomial(input: &str) -> Result<Vec<i64>, String> {
    parse_polynomial_bigint(input)?
        .into_iter()
        .map(|coeff| {
            coeff
                .to_i64()
                .ok_or_else(|| format!("integer '{}' does not fit in i64", coeff))
        })
        .collect()
}

/// Parse a single line into arbitrary-size integer polynomial coefficients
/// (ascending degree).
///
/// Tries expanded polynomial format first (if a letter is found),
/// then falls back to coefficient list parsing.
pub fn parse_polynomial_bigint(input: &str) -> Result<Vec<BigInt>, String> {
    if input.len() > MAX_POLYNOMIAL_INPUT_BYTES {
        return Err(format!(
            "polynomial input is {} bytes; the limit is {MAX_POLYNOMIAL_INPUT_BYTES}",
            input.len()
        ));
    }
    let s = input.trim();
    if s.is_empty() {
        return Err("empty input".to_string());
    }

    // Strip outer brackets/braces/parens
    let s = strip_brackets(s);

    // Detect if this looks like an expanded polynomial (contains a letter)
    if has_variable(s) {
        parse_expanded_bigint(s)
    } else {
        parse_coeff_list_bigint(s)
    }
}

/// Parse multiple lines, skipping blanks and comments.
pub fn parse_polynomials(input: &str) -> Vec<Result<Vec<i64>, String>> {
    parse_polynomials_bigint(input)
        .into_iter()
        .map(|result| {
            result.and_then(|coeffs| {
                coeffs
                    .into_iter()
                    .map(|coeff| {
                        coeff
                            .to_i64()
                            .ok_or_else(|| format!("integer '{}' does not fit in i64", coeff))
                    })
                    .collect()
            })
        })
        .collect()
}

/// Parse multiple lines into arbitrary-size integer coefficient vectors,
/// skipping blanks and comments.
pub fn parse_polynomials_bigint(input: &str) -> Vec<Result<Vec<BigInt>, String>> {
    input
        .lines()
        .filter(|l| {
            let t = l.trim();
            !t.is_empty() && !t.starts_with('#')
        })
        .map(parse_polynomial_bigint)
        .collect()
}

// ---------------------------------------------------------------------------
// Coefficient list: "1, 2, 3" or "1 2 3"
// ---------------------------------------------------------------------------

fn parse_coeff_list_bigint(s: &str) -> Result<Vec<BigInt>, String> {
    // Determine separator: comma if present, else whitespace
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

    let coeffs: Result<Vec<BigInt>, _> = parts
        .iter()
        .map(|p| {
            p.trim()
                .parse::<BigInt>()
                .map_err(|e| format!("invalid integer '{}': {}", p.trim(), e))
        })
        .collect();

    coeffs
}

// ---------------------------------------------------------------------------
// Expanded polynomial: "3t^2 - 5t + 1"
// ---------------------------------------------------------------------------

fn has_variable(s: &str) -> bool {
    s.chars().any(|c| c.is_ascii_alphabetic())
}

/// Detect which single variable is used (the only ASCII letter in the string).
fn detect_variable(s: &str) -> Result<char, String> {
    let vars: std::collections::HashSet<char> =
        s.chars().filter(|c| c.is_ascii_alphabetic()).collect();

    match vars.len() {
        0 => Err("no variable found".to_string()),
        1 => Ok(*vars.iter().next().unwrap()),
        _ => Err(format!(
            "multiple variables found: {:?} (expected a single variable)",
            vars
        )),
    }
}

/// Parse an expanded polynomial string like "3t^2 - 5t + 1" or "x**3 + 2*x - 7".
fn parse_expanded_bigint(s: &str) -> Result<Vec<BigInt>, String> {
    let var = detect_variable(s)?;

    // Normalize: remove spaces around operators, but keep spacing for tokenizing
    // Strategy: split into terms by + and - (keeping the sign)
    let terms = split_into_terms(s);

    let mut coeffs: Vec<BigInt> = Vec::new();

    for term in &terms {
        let term = term.trim();
        if term.is_empty() {
            continue;
        }
        let (coeff, deg) = parse_term_bigint(term, var)?;

        // Extend coefficient vector if needed
        let required_len = deg
            .checked_add(1)
            .ok_or_else(|| format!("polynomial exponent {deg} is too large"))?;
        if required_len > MAX_POLYNOMIAL_COEFFICIENTS {
            return Err(format!(
                "polynomial exponent {deg} exceeds the maximum supported degree {}",
                MAX_POLYNOMIAL_COEFFICIENTS - 1
            ));
        }
        if required_len > coeffs.len() {
            coeffs.resize(required_len, BigInt::from(0));
        }
        coeffs[deg] += coeff;
    }

    if coeffs.is_empty() {
        coeffs.push(BigInt::from(0));
    }

    Ok(coeffs)
}

/// Split a polynomial string into signed terms.
/// "3t^2 - 5t + 1" → ["+3t^2", "-5t", "+1"]
fn split_into_terms(s: &str) -> Vec<String> {
    let mut terms = Vec::new();
    let mut current = String::new();

    for c in s.chars() {
        if (c == '+' || c == '-') && !current.trim().is_empty() {
            // Check if previous char was ^ or * (exponent sign, not term separator)
            let last_significant = current.trim_end().chars().last();
            if last_significant == Some('^') || last_significant == Some('*') {
                current.push(c);
                continue;
            }
            terms.push(current.clone());
            current.clear();
        }
        current.push(c);
    }
    if !current.trim().is_empty() {
        terms.push(current);
    }

    terms
}

/// Parse a single term like "3t^2", "-t", "5", "t", "-3*t**2".
/// Returns (coefficient, degree).
fn parse_term_bigint(term: &str, var: char) -> Result<(BigInt, usize), String> {
    let s: String = term.chars().filter(|c| !c.is_whitespace()).collect();

    if s.is_empty() {
        return Ok((BigInt::from(0), 0));
    }

    // Find the variable position
    let var_pos = s.find(var);

    match var_pos {
        None => {
            // Pure constant: no variable
            let val = s
                .parse::<BigInt>()
                .map_err(|e| format!("cannot parse '{}' as integer: {}", s, e))?;
            Ok((val, 0))
        }
        Some(pos) => {
            // Extract coefficient part (before the variable)
            let coeff_str = &s[..pos];
            // Remove trailing * if present (e.g., "3*t")
            let coeff_str = coeff_str.trim_end_matches('*');

            let coeff = if coeff_str.is_empty() || coeff_str == "+" {
                BigInt::from(1)
            } else if coeff_str == "-" {
                BigInt::from(-1)
            } else {
                coeff_str
                    .parse::<BigInt>()
                    .map_err(|e| format!("cannot parse coefficient '{}': {}", coeff_str, e))?
            };

            // Extract degree part (after the variable)
            let after_var = &s[pos + var.len_utf8()..];

            let degree = if after_var.is_empty() {
                1
            } else {
                // Strip ^ or **
                let exp_str = after_var.trim_start_matches("**").trim_start_matches('^');

                if exp_str.is_empty() {
                    1
                } else {
                    exp_str
                        .parse::<usize>()
                        .map_err(|e| format!("cannot parse exponent '{}': {}", exp_str, e))?
                }
            };

            Ok((coeff, degree))
        }
    }
}

fn strip_brackets(s: &str) -> &str {
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

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_comma_separated() {
        assert_eq!(parse_polynomial("1, 2, 3").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_polynomial("  1,2,3  ").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn test_space_separated() {
        assert_eq!(parse_polynomial("1 2 3").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_polynomial("  1   2   3  ").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn test_bracketed() {
        assert_eq!(parse_polynomial("[1, 2, 3]").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_polynomial("{1, 2, 3}").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_polynomial("(1, 2, 3)").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn test_expanded_basic() {
        assert_eq!(parse_polynomial("1 + 2t + 3t^2").unwrap(), vec![1, 2, 3]);
        assert_eq!(parse_polynomial("3t^2 + 2t + 1").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn test_expanded_negative() {
        assert_eq!(parse_polynomial("t^3 - 3t + 2").unwrap(), vec![2, -3, 0, 1]);
    }

    #[test]
    fn test_expanded_implicit_coeff() {
        assert_eq!(parse_polynomial("t").unwrap(), vec![0, 1]);
        assert_eq!(parse_polynomial("-t").unwrap(), vec![0, -1]);
        assert_eq!(parse_polynomial("t^2 + t + 1").unwrap(), vec![1, 1, 1]);
    }

    #[test]
    fn test_expanded_star_multiply() {
        assert_eq!(parse_polynomial("3*t^2 + 2*t + 1").unwrap(), vec![1, 2, 3]);
    }

    #[test]
    fn test_expanded_double_star() {
        assert_eq!(parse_polynomial("t**3 - 1").unwrap(), vec![-1, 0, 0, 1]);
    }

    #[test]
    fn test_expanded_variable_x() {
        assert_eq!(parse_polynomial("x^2 + 2x + 1").unwrap(), vec![1, 2, 1]);
    }

    #[test]
    fn test_expanded_variable_n() {
        assert_eq!(parse_polynomial("n^3 + n").unwrap(), vec![0, 1, 0, 1]);
    }

    #[test]
    fn test_constant() {
        assert_eq!(parse_polynomial("42").unwrap(), vec![42]);
    }

    #[test]
    fn test_bigint_coefficients() {
        let big = BigInt::from(10u64).pow(40);
        assert_eq!(
            parse_polynomial_bigint(&format!("{big}, 2, 1")).unwrap(),
            vec![big.clone(), BigInt::from(2), BigInt::from(1)]
        );
        assert_eq!(
            parse_polynomial_bigint(&format!("{big}t^2 + 2t + 1")).unwrap(),
            vec![BigInt::from(1), BigInt::from(2), big]
        );
    }

    #[test]
    fn rejects_exponents_that_cannot_be_allocated_safely() {
        let error = parse_polynomial_bigint(&format!("t^{}", usize::MAX)).unwrap_err();
        assert!(error.contains("too large"));

        let error =
            parse_polynomial_bigint(&format!("t^{}", MAX_POLYNOMIAL_COEFFICIENTS)).unwrap_err();
        assert!(error.contains("maximum supported degree"));
    }

    #[test]
    fn rejects_oversized_text_and_coefficient_lists() {
        let oversized_text = "1".repeat(MAX_POLYNOMIAL_INPUT_BYTES + 1);
        assert!(parse_polynomial_bigint(&oversized_text)
            .unwrap_err()
            .contains("bytes"));

        let oversized_list = std::iter::repeat_n("0", MAX_POLYNOMIAL_COEFFICIENTS + 1)
            .collect::<Vec<_>>()
            .join(",");
        assert!(parse_polynomial_bigint(&oversized_list)
            .unwrap_err()
            .contains("coefficients"));
    }

    #[test]
    fn test_negative_coeffs() {
        assert_eq!(parse_polynomial("-1, 0, 1").unwrap(), vec![-1, 0, 1]);
    }

    #[test]
    fn test_parse_multiple() {
        let input = "1, 11, 11, 1\n# comment\n1 + 4t + t^2\n\n[1, 2, 1]";
        let results = parse_polynomials(input);
        assert_eq!(results.len(), 3);
        assert_eq!(results[0].as_ref().unwrap(), &vec![1, 11, 11, 1]);
        assert_eq!(results[1].as_ref().unwrap(), &vec![1, 4, 1]);
        assert_eq!(results[2].as_ref().unwrap(), &vec![1, 2, 1]);
    }

    #[test]
    fn test_eulerian() {
        // Eulerian A_4(t) in various formats
        let expected = vec![1, 11, 11, 1];
        assert_eq!(parse_polynomial("1, 11, 11, 1").unwrap(), expected);
        assert_eq!(parse_polynomial("1 11 11 1").unwrap(), expected);
        assert_eq!(parse_polynomial("[1, 11, 11, 1]").unwrap(), expected);
        assert_eq!(parse_polynomial("{1, 11, 11, 1}").unwrap(), expected);
        assert_eq!(parse_polynomial("1 + 11t + 11t^2 + t^3").unwrap(), expected);
        assert_eq!(parse_polynomial("t^3 + 11t^2 + 11t + 1").unwrap(), expected);
    }
}
