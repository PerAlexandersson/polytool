//! Human-readable display of coefficient vectors (`format_poly*`).
//!
//! Coefficients are in ascending degree order; the zero polynomial (empty or
//! all-zero vector) is displayed as `0`.

use num_bigint::BigInt;
use num_traits::Zero;

/// Format an `i64` coefficient vector as a human-readable polynomial string.
///
/// Uses the variable name `t` by default. Coefficients of 1 or -1 are displayed
/// without the leading "1". Examples:
///
/// - `[1, 2, 1]` → `"1 + 2t + t^2"`
/// - `[0, -1, 0, 3]` → `"-t + 3t^3"`
/// - `[]` or `[0]` → `"0"`
pub fn format_poly(coeffs: &[i64]) -> String {
    format_poly_var(coeffs, "t")
}

/// Format an `i64` coefficient vector with a custom variable name.
pub fn format_poly_var(coeffs: &[i64], var: &str) -> String {
    let mut terms = Vec::new();
    for (i, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        let term = match (c, i) {
            (_, 0) => format!("{}", c),
            (1, 1) => var.to_string(),
            (-1, 1) => format!("-{}", var),
            (_, 1) => format!("{}{}", c, var),
            (1, e) => format!("{}^{}", var, e),
            (-1, e) => format!("-{}^{}", var, e),
            (_, e) => format!("{}{}^{}", c, var, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

/// Format a `BigInt` coefficient vector as a human-readable polynomial string.
pub fn format_poly_bigint_coeffs(coeffs: &[BigInt]) -> String {
    format_poly_bigint_var(coeffs, "t")
}

/// Format a `BigInt` coefficient vector with a custom variable name.
pub fn format_poly_bigint_var(coeffs: &[BigInt], var: &str) -> String {
    let one = BigInt::from(1);
    let minus_one = BigInt::from(-1);
    let mut terms = Vec::new();
    for (i, c) in coeffs.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let term = match i {
            0 => c.to_string(),
            1 if c == &one => var.to_string(),
            1 if c == &minus_one => format!("-{}", var),
            1 => format!("{}{}", c, var),
            e if c == &one => format!("{}^{}", var, e),
            e if c == &minus_one => format!("-{}^{}", var, e),
            e => format!("{}{}^{}", c, var, e),
        };
        terms.push(term);
    }
    if terms.is_empty() {
        return "0".to_string();
    }
    let mut result = terms[0].clone();
    for term in &terms[1..] {
        if let Some(rest) = term.strip_prefix('-') {
            result.push_str(" - ");
            result.push_str(rest);
        } else {
            result.push_str(" + ");
            result.push_str(term);
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_poly() {
        assert_eq!(format_poly(&[1, 2, 1]), "1 + 2t + t^2");
        assert_eq!(format_poly(&[0, -1, 0, 3]), "-t + 3t^3");
        assert_eq!(format_poly(&[]), "0");
        assert_eq!(format_poly(&[0]), "0");
        assert_eq!(format_poly(&[5]), "5");
        assert_eq!(format_poly(&[0, 0, -1]), "-t^2");
        assert_eq!(format_poly_var(&[1, 1], "q"), "1 + q");
    }
}
