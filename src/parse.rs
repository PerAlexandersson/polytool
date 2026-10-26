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
//!
//! Entry points: [`parse_polynomial`] / [`parse_polynomial_bigint`] for one
//! line (`i64` with an explicit "does not fit in `i64`" error, or `BigInt`) and
//! [`parse_polynomials`] / [`parse_polynomials_bigint`] for multi-line input
//! (blank lines and lines starting with `#` are skipped; the result has one
//! `Result` per remaining line).  Errors are human-readable `String`s.
//!
//! Limits and caveats:
//!
//! - Input longer than [`MAX_POLYNOMIAL_INPUT_BYTES`] and polynomials with more
//!   than [`MAX_POLYNOMIAL_COEFFICIENTS`] coefficients are rejected.
//! - Coefficients are integers only (no fractions, decimals or scientific
//!   notation) and a single variable letter is required for the expanded form;
//!   products and powers of parenthesized sums such as `(1+x)^2` are **not**
//!   expanded.
//! - Coefficient lists are returned as written (`"1, 2, 0"` gives `[1, 2, 0]`,
//!   so a nominal degree survives); the expanded form is trimmed (`"x - x"`
//!   gives `[0]`).
//! - Ambiguous input is rejected rather than guessed: digits separated by
//!   whitespace inside one term (`"2 3 x"`), a variable followed directly by
//!   digits (`"x2"`; write `x^2`), and scientific notation (`"1e3"`).

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

    // Cancellation (`x - x`) leaves trailing zeros; keep the canonical zero `[0]`.
    while coeffs.len() > 1 && coeffs.last().is_some_and(|c| c == &BigInt::from(0)) {
        coeffs.pop();
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

/// Parse a single signed term like "3t^2", "-t", "5", "t", "-3*t**2" or
/// "2 t ^ 3".  Returns (coefficient, degree).
///
/// Whitespace may separate the coefficient, `*`, the variable, `^`/`**` and
/// the exponent, but never the digits of one number.  A variable followed
/// directly by digits (`t2`, or `1e3` read as `e`) is rejected as ambiguous.
fn parse_term_bigint(term: &str, var: char) -> Result<(BigInt, usize), String> {
    let trimmed = term.trim();
    let (negative, rest) = match trimmed.chars().next() {
        Some('+') => (false, trimmed[1..].trim_start()),
        Some('-') => (true, trimmed[1..].trim_start()),
        _ => (false, trimmed),
    };
    if rest.is_empty() {
        return Err(format!("missing term after sign in '{trimmed}'"));
    }

    let digits_end = rest
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(rest.len());
    let (digits, after) = rest.split_at(digits_end);
    let after = after.trim_start();
    if after.starts_with(|c: char| c.is_ascii_digit()) {
        return Err(format!(
            "digits separated by whitespace in '{trimmed}'; write one integer"
        ));
    }
    let magnitude = if digits.is_empty() {
        BigInt::from(1)
    } else {
        digits
            .parse::<BigInt>()
            .map_err(|e| format!("cannot parse coefficient '{digits}': {e}"))?
    };
    let coeff = if negative { -magnitude } else { magnitude };

    if after.is_empty() {
        return Ok((coeff, 0));
    }
    let after = match after.strip_prefix('*') {
        Some(_) if digits.is_empty() => {
            return Err(format!("missing coefficient before '*' in '{trimmed}'"))
        }
        Some(rest) if rest.starts_with('*') => {
            return Err(format!(
                "unexpected '**' after a coefficient in '{trimmed}'"
            ))
        }
        Some(rest) => rest.trim_start(),
        None => after,
    };
    let Some(after_var) = after.strip_prefix(var) else {
        return Err(format!("unexpected '{after}' in term '{trimmed}'"));
    };
    let after_var = after_var.trim_start();

    let degree = if after_var.is_empty() {
        1
    } else if let Some(exponent) = after_var
        .strip_prefix("**")
        .or_else(|| after_var.strip_prefix('^'))
    {
        let exponent = exponent.trim();
        if exponent.is_empty() || !exponent.chars().all(|c| c.is_ascii_digit()) {
            return Err(format!(
                "exponent must be a nonnegative integer in '{trimmed}'"
            ));
        }
        exponent
            .parse::<usize>()
            .map_err(|_| format!("polynomial exponent {exponent} is too large"))?
    } else if after_var.starts_with(|c: char| c.is_ascii_digit()) {
        if (var == 'e' || var == 'E') && !digits.is_empty() {
            return Err(format!(
                "scientific notation such as '{trimmed}' is not supported; write the integer"
            ));
        }
        return Err(format!(
            "ambiguous '{var}{after_var}' in '{trimmed}'; write '{var}^{after_var}' for a power"
        ));
    } else {
        return Err(format!("unexpected '{after_var}' in term '{trimmed}'"));
    };

    Ok((coeff, degree))
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

    #[test]
    fn rejects_ambiguous_expanded_input() {
        let error = parse_polynomial("1e3").unwrap_err();
        assert!(error.contains("scientific notation"), "{error}");
        let error = parse_polynomial("2 3 x").unwrap_err();
        assert!(error.contains("whitespace"), "{error}");
        let error = parse_polynomial("x2 + 1").unwrap_err();
        assert!(error.contains("x^2"), "{error}");
        for bad in ["x^", "x^-1", "*x", "2**x", "2x x", "x +", "2y + x"] {
            assert!(parse_polynomial(bad).is_err(), "{bad} should be rejected");
        }
    }

    #[test]
    fn accepts_spacing_around_operators() {
        assert_eq!(
            parse_polynomial("2 x ^ 3 - 3 * x + 1").unwrap(),
            vec![1, -3, 0, 2]
        );
        assert_eq!(parse_polynomial("- x ** 2").unwrap(), vec![0, 0, -1]);
    }

    #[test]
    fn expanded_form_is_trimmed_but_lists_are_not() {
        assert_eq!(parse_polynomial("x - x").unwrap(), vec![0]);
        assert_eq!(parse_polynomial("1 + x^2 - x^2").unwrap(), vec![1]);
        assert_eq!(parse_polynomial("1, 2, 0").unwrap(), vec![1, 2, 0]);
    }
}
