//! Plain-text and LaTeX formatting of bivariate coefficients and recurrences,
//! plus exact rational coefficient formatting and parsing.

use super::*;
use num_traits::One;
use num_traits::Zero;
use std::fmt;

fn fmt_rational(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("{}/{}", r.numer(), r.denom())
    }
}

/// Format an exact rational coefficient as an integer or `num/den` string.
pub fn format_rational_coeff(r: &BigRational) -> String {
    fmt_rational(r)
}

/// Parse an exact rational coefficient from an integer or `num/den` string.
pub fn parse_rational_coeff(input: &str) -> Result<BigRational, ParseRationalError> {
    let token = input.trim();
    if token.is_empty() {
        return Err(ParseRationalError::Empty);
    }
    if let Some((num, den)) = token.split_once('/') {
        let numerator =
            num.trim()
                .parse::<BigInt>()
                .map_err(|e| ParseRationalError::InvalidNumerator {
                    text: num.trim().to_string(),
                    message: e.to_string(),
                })?;
        let denominator =
            den.trim()
                .parse::<BigInt>()
                .map_err(|e| ParseRationalError::InvalidDenominator {
                    text: den.trim().to_string(),
                    message: e.to_string(),
                })?;
        if denominator.is_zero() {
            return Err(ParseRationalError::ZeroDenominator {
                token: token.to_string(),
            });
        }
        Ok(BigRational::new(numerator, denominator))
    } else {
        let integer = token
            .parse::<BigInt>()
            .map_err(|e| ParseRationalError::InvalidInteger {
                token: token.to_string(),
                message: e.to_string(),
            })?;
        Ok(BigRational::from_integer(integer))
    }
}

fn fmt_monomial(c: &BigRational, n_pow: usize, t_pow: usize) -> String {
    let var = match (n_pow, t_pow) {
        (0, 0) => return fmt_rational(c),
        (1, 0) => "n".into(),
        (i, 0) => format!("n^{i}"),
        (0, 1) => "t".into(),
        (0, j) => format!("t^{j}"),
        (1, 1) => "nt".into(),
        (i, 1) => format!("n^{i}t"),
        (1, j) => format!("nt^{j}"),
        (i, j) => format!("n^{i}t^{j}"),
    };

    if c.is_one() {
        var
    } else if *c == BigRational::from(BigInt::from(-1)) {
        format!("-{var}")
    } else {
        format!("{}{var}", fmt_rational(c))
    }
}

impl fmt::Display for BivarPoly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mut terms: Vec<String> = Vec::new();
        // Iterate j (t-power) first, then i (n-power), for natural ordering.
        let width = self.coefficient_width();
        if width == 0 {
            return write!(f, "0");
        }
        let max_j = width - 1;
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let Some(c) = self.coefficient(i, j) else {
                    continue;
                };
                if c.is_zero() {
                    continue;
                }
                terms.push(fmt_monomial(c, i, j));
            }
        }

        if terms.is_empty() {
            return write!(f, "0");
        }

        let mut result = terms[0].clone();
        for t in &terms[1..] {
            if let Some(rest) = t.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(rest);
            } else {
                result.push_str(" + ");
                result.push_str(t);
            }
        }
        write!(f, "{result}")
    }
}

fn fmt_rational_latex(r: &BigRational) -> String {
    if r.denom() == &BigInt::one() {
        format!("{}", r.numer())
    } else {
        format!("\\frac{{{}}}{{{}}}", r.numer(), r.denom())
    }
}

fn fmt_monomial_latex(c: &BigRational, n_pow: usize, t_pow: usize) -> String {
    let var = match (n_pow, t_pow) {
        (0, 0) => return fmt_rational_latex(c),
        (1, 0) => "n".into(),
        (i, 0) => format!("n^{{{i}}}"),
        (0, 1) => "t".into(),
        (0, j) => format!("t^{{{j}}}"),
        (1, 1) => "nt".into(),
        (i, 1) => format!("n^{{{i}}}t"),
        (1, j) => format!("nt^{{{j}}}"),
        (i, j) => format!("n^{{{i}}}t^{{{j}}}"),
    };

    if c.is_one() {
        var
    } else if *c == BigRational::from(BigInt::from(-1)) {
        format!("-{var}")
    } else {
        format!("{}{var}", fmt_rational_latex(c))
    }
}

impl BivarPoly {
    /// Format as a LaTeX expression.
    pub fn to_latex(&self) -> String {
        let mut terms: Vec<String> = Vec::new();
        let width = self.coefficient_width();
        if width == 0 {
            return "0".into();
        }
        let max_j = width - 1;
        let max_i = self.coeffs.len().saturating_sub(1);
        for j in 0..=max_j {
            for i in 0..=max_i {
                let Some(c) = self.coefficient(i, j) else {
                    continue;
                };
                if c.is_zero() {
                    continue;
                }
                terms.push(fmt_monomial_latex(c, i, j));
            }
        }
        if terms.is_empty() {
            return "0".into();
        }
        let mut result = terms[0].clone();
        for t in &terms[1..] {
            if let Some(rest) = t.strip_prefix('-') {
                result.push_str(" - ");
                result.push_str(rest);
            } else {
                result.push_str(" + ");
                result.push_str(t);
            }
        }
        result
    }
}

fn fmt_poly_ref_latex(offset: usize, deriv_order: usize) -> String {
    let sub = if offset == 0 {
        "n".into()
    } else {
        format!("n-{offset}")
    };
    match deriv_order {
        0 => format!("P({sub})"),
        1 => format!("P'({sub})"),
        d => format!("P^{{({d})}}({sub})"),
    }
}

pub(super) fn fmt_poly_ref_latex_with_sign(
    sign: RecurrenceSign,
    offset: usize,
    deriv_order: usize,
) -> String {
    let pref = fmt_poly_ref_latex(offset, deriv_order);
    match sign {
        RecurrenceSign::None => pref,
        RecurrenceSign::AlternatingN => format!("(-1)^n {pref}"),
    }
}

fn fmt_poly_ref(offset: usize, deriv_order: usize) -> String {
    let sub = if offset == 0 {
        "n".into()
    } else {
        format!("n-{offset}")
    };
    match deriv_order {
        0 => format!("P({sub})"),
        1 => format!("P'({sub})"),
        d => format!("P^({d})({sub})"),
    }
}

fn fmt_poly_ref_with_sign(sign: RecurrenceSign, offset: usize, deriv_order: usize) -> String {
    let pref = fmt_poly_ref(offset, deriv_order);
    match sign {
        RecurrenceSign::None => pref,
        RecurrenceSign::AlternatingN => format!("(-1)^n {pref}"),
    }
}

impl fmt::Display for Recurrence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // LHS
        if let Some(ref denom) = self.denominator {
            write!(f, "({denom}) ")?;
        }
        write!(f, "P(n) = ")?;

        let mut first = true;
        for term in &self.terms {
            let cs = format!("{}", term.coeff);
            let pref = fmt_poly_ref_with_sign(term.sign, term.offset, term.deriv_order);

            if first {
                first = false;
                if cs == "1" {
                    write!(f, "{pref}")?;
                } else if cs == "-1" {
                    write!(f, "-{pref}")?;
                } else if term.coeff.num_terms() > 1 {
                    write!(f, "({cs}) {pref}")?;
                } else {
                    write!(f, "{cs} {pref}")?;
                }
            } else {
                // Determine sign for nice formatting.
                if cs.starts_with('-') && term.coeff.is_monomial() {
                    // Single negative term: pull out the minus sign.
                    let pos = &cs[1..];
                    if pos == "1" {
                        write!(f, " - {pref}")?;
                    } else {
                        write!(f, " - {pos} {pref}")?;
                    }
                } else if term.coeff.num_terms() > 1 {
                    write!(f, " + ({cs}) {pref}")?;
                } else if cs == "1" {
                    write!(f, " + {pref}")?;
                } else {
                    write!(f, " + {cs} {pref}")?;
                }
            }
        }

        if let Some(ref inh) = self.inhomogeneous {
            if !first {
                write!(f, " + ")?;
            }
            write!(f, "{inh}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    #[test]
    fn ragged_bivariate_polynomials_format_without_loss_or_panics() {
        let wider_later_row = BivarPoly {
            coeffs: vec![vec![br(1, 1)], vec![br(0, 1), br(1, 1)]],
        };
        assert_eq!(wider_later_row.to_string(), "1 + nt");
        assert_eq!(wider_later_row.to_latex(), "1 + nt");
        assert_eq!(wider_later_row.to_mathematica_code(), "1 + n*t");

        let shorter_later_row = BivarPoly {
            coeffs: vec![vec![br(1, 1), br(2, 1)], vec![br(3, 1)]],
        };
        assert_eq!(shorter_later_row.to_string(), "1 + 3n + 2t");

        let empty = BivarPoly { coeffs: vec![] };
        assert!(empty.is_zero());
        assert!(!empty.is_one());
        assert_eq!(empty.to_string(), "0");
        assert_eq!(empty.to_mathematica_code(), "0");
    }
}
