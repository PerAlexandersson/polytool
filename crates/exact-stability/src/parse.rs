//! Bounded exact Gaussian-rational and univariate expression parsing.

use crate::{BigRational as Q, GaussianRational as G};
use num_bigint::BigInt;
use num_traits::Zero;
use std::fmt;

/// Resource bounds checked during parsing, including intermediate expressions.
#[derive(Clone, Copy, Debug)]
pub struct ParseLimits {
    /// Largest allowed degree.
    pub max_degree: usize,
    /// Largest numerator or denominator bit length of either coefficient part.
    pub max_coefficient_bits: u64,
    /// Largest byte length of an input coefficient literal or integer token.
    pub max_literal_bytes: usize,
}

/// An invalid expression or an exceeded parsing bound.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ParseError(pub String);

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for ParseError {}

fn error(message: &str) -> ParseError {
    ParseError(message.into())
}

impl ParseLimits {
    /// Check each exact part before expensive root counting.
    pub fn validate(&self, coefficients: &[G]) -> Result<(), ParseError> {
        if coefficients
            .iter()
            .rposition(|c| !c.is_zero())
            .is_some_and(|d| d > self.max_degree)
        {
            return Err(error("polynomial exceeds max-degree"));
        }
        for c in coefficients {
            for q in [&c.real, &c.imaginary] {
                if q.numer().bits() > self.max_coefficient_bits
                    || q.denom().bits() > self.max_coefficient_bits
                {
                    return Err(error("coefficient exceeds max-coefficient-bits"));
                }
            }
        }
        Ok(())
    }
}

/// Parse an exact rational string; floating point and zero denominators fail.
pub fn parse_rational(input: &str, limits: ParseLimits) -> Result<Q, ParseError> {
    let parts: Vec<_> = input.trim().split('/').collect();
    if parts.is_empty() || parts.len() > 2 {
        return Err(error("expected an exact rational"));
    }
    let integer = |s: &str| -> Result<BigInt, ParseError> {
        let s = s.trim();
        if s.len() > limits.max_literal_bytes
            || s.len() as u64 > limits.max_coefficient_bits.saturating_add(2)
        {
            return Err(error("integer text exceeds coefficient bit budget"));
        }
        s.parse()
            .map_err(|_| error("expected an exact integer or rational"))
    };
    let numerator = integer(parts[0])?;
    let denominator = if parts.len() == 2 {
        integer(parts[1])?
    } else {
        1.into()
    };
    if denominator.is_zero() {
        return Err(error("zero denominator"));
    }
    if numerator.bits() > limits.max_coefficient_bits
        || denominator.bits() > limits.max_coefficient_bits
    {
        return Err(error("coefficient exceeds max-coefficient-bits"));
    }
    Ok(Q::new(numerator, denominator))
}

/// Parse a constant Gaussian expression, with `I` denoting the imaginary unit.
pub fn parse_gaussian(input: &str, limits: ParseLimits) -> Result<G, ParseError> {
    if input.len() > limits.max_literal_bytes {
        return Err(error("coefficient literal exceeds text limit"));
    }
    let mut parser = Parser::new(input, limits);
    let p = parser.expression(0)?;
    parser.space();
    if parser.position != parser.input.len() || parser.variable.is_some() {
        return Err(error("expected an exact Gaussian constant"));
    }
    Ok(p.first().cloned().unwrap_or_default())
}

/// Parse ascending comma-separated coefficients or one-variable expressions.
/// Integer/rational arithmetic, implicit products, parentheses and natural
/// powers are accepted; denominators must be syntactically constant.
pub fn parse_univariate(input: &str, limits: ParseLimits) -> Result<Vec<G>, ParseError> {
    let input = input.trim();
    if input.is_empty() {
        return Err(error("empty polynomial"));
    }
    let list = input.strip_prefix('[').and_then(|s| s.strip_suffix(']'));
    let mut result = if list.is_some() || input.contains(',') {
        let text = list.unwrap_or(input);
        if text.trim().is_empty() {
            vec![]
        } else {
            let mut coefficients = Vec::new();
            for part in text.split(',') {
                if coefficients.len() > limits.max_degree {
                    return Err(error("coefficient list exceeds max-degree"));
                }
                coefficients.push(parse_gaussian(part, limits)?);
            }
            coefficients
        }
    } else {
        let mut parser = Parser::new(input, limits);
        let p = parser.expression(0)?;
        parser.space();
        if parser.position != parser.input.len() {
            return Err(error("unexpected token in polynomial"));
        }
        p
    };
    while result.last().is_some_and(G::is_zero) {
        result.pop();
    }
    limits.validate(&result)?;
    Ok(result)
}

struct Parser<'a> {
    input: &'a [u8],
    position: usize,
    variable: Option<String>,
    limits: ParseLimits,
}
impl<'a> Parser<'a> {
    fn new(input: &'a str, limits: ParseLimits) -> Self {
        Self {
            input: input.as_bytes(),
            position: 0,
            variable: None,
            limits,
        }
    }
    fn space(&mut self) {
        while self
            .input
            .get(self.position)
            .is_some_and(u8::is_ascii_whitespace)
        {
            self.position += 1;
        }
    }
    fn take(&mut self, c: u8) -> bool {
        self.space();
        if self.input.get(self.position) == Some(&c) {
            self.position += 1;
            true
        } else {
            false
        }
    }
    fn expression(&mut self, depth: usize) -> Result<Vec<G>, ParseError> {
        if depth > 64 {
            return Err(error("expression nesting exceeds 64"));
        }
        let mut p = self.product(depth + 1)?;
        loop {
            let sign = if self.take(b'+') {
                1
            } else if self.take(b'-') {
                -1
            } else {
                break;
            };
            let q = self.product(depth + 1)?;
            p.resize(p.len().max(q.len()), G::zero());
            for (a, b) in p.iter_mut().zip(q) {
                *a = if sign == 1 { a.add(&b) } else { a.subtract(&b) };
            }
            self.limits.validate(&p)?;
        }
        Ok(p)
    }
    fn multiply(&self, p: &[G], q: &[G]) -> Result<Vec<G>, ParseError> {
        if p.is_empty() || q.is_empty() {
            return Ok(vec![]);
        }
        let size = p
            .len()
            .checked_add(q.len())
            .and_then(|n| n.checked_sub(1))
            .ok_or_else(|| error("degree overflow"))?;
        if size.saturating_sub(1) > self.limits.max_degree {
            return Err(error("product exceeds max-degree"));
        }
        let mut result = vec![G::zero(); size];
        for (i, a) in p.iter().enumerate() {
            for (j, b) in q.iter().enumerate() {
                result[i + j] = result[i + j].add(&a.multiply(b));
            }
        }
        self.limits.validate(&result)?;
        Ok(result)
    }
    fn product(&mut self, depth: usize) -> Result<Vec<G>, ParseError> {
        let mut p = self.signed_power(depth)?;
        loop {
            let explicit = self.take(b'*');
            if !explicit && self.take(b'/') {
                let start = self.position;
                let q = self.signed_power(depth)?;
                if self.input[start..self.position]
                    .iter()
                    .any(u8::is_ascii_alphabetic)
                    && self.input[start..self.position]
                        .iter()
                        .any(|&c| c != b'I' && c.is_ascii_alphabetic())
                {
                    return Err(error("nonconstant denominator"));
                }
                if q.len() > 1 {
                    return Err(error("nonconstant denominator"));
                }
                let inverse = G::one()
                    .checked_divide(q.first().unwrap_or(&G::zero()))
                    .ok_or_else(|| error("zero denominator"))?;
                for a in &mut p {
                    *a = a.multiply(&inverse);
                }
                self.limits.validate(&p)?;
                continue;
            }
            self.space();
            let implicit = self
                .input
                .get(self.position)
                .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'(' || c.is_ascii_digit());
            if !explicit && !implicit {
                break;
            }
            let q = self.signed_power(depth)?;
            p = self.multiply(&p, &q)?;
        }
        Ok(p)
    }
    fn signed_power(&mut self, depth: usize) -> Result<Vec<G>, ParseError> {
        if depth > 64 {
            return Err(error("expression nesting exceeds 64"));
        }
        if self.take(b'+') {
            return self.signed_power(depth + 1);
        }
        if self.take(b'-') {
            return Ok(self
                .signed_power(depth + 1)?
                .iter()
                .map(|x| G::zero().subtract(x))
                .collect());
        }
        let p = self.atom(depth + 1)?;
        if !self.take(b'^') {
            return Ok(p);
        }
        self.space();
        let start = self.position;
        while self
            .input
            .get(self.position)
            .is_some_and(u8::is_ascii_digit)
        {
            self.position += 1;
        }
        let text = std::str::from_utf8(&self.input[start..self.position])
            .map_err(|_| error("invalid exponent"))?;
        let mut n: usize = text
            .parse()
            .map_err(|_| error("expected a natural exponent"))?;
        if p.len() > 1 && n > self.limits.max_degree {
            return Err(error("power exceeds max-degree"));
        }
        let mut base = p;
        let mut result = vec![G::one()];
        while n > 0 {
            if n % 2 == 1 {
                result = self.multiply(&result, &base)?;
            }
            n /= 2;
            if n > 0 {
                base = self.multiply(&base, &base)?;
            }
        }
        Ok(result)
    }
    fn atom(&mut self, depth: usize) -> Result<Vec<G>, ParseError> {
        self.space();
        if self.take(b'(') {
            let p = self.expression(depth)?;
            if !self.take(b')') {
                return Err(error("missing closing parenthesis"));
            }
            return Ok(p);
        }
        let start = self.position;
        if self.input.get(start).is_some_and(u8::is_ascii_digit) {
            while self
                .input
                .get(self.position)
                .is_some_and(u8::is_ascii_digit)
            {
                self.position += 1;
            }
            if matches!(self.input.get(self.position), Some(b'e' | b'E')) {
                return Err(error("scientific notation is not an exact literal"));
            }
            let text = std::str::from_utf8(&self.input[start..self.position])
                .map_err(|_| error("invalid integer"))?;
            return Ok(vec![G::from(parse_rational(text, self.limits)?)]);
        }
        if self.input.get(start).is_some_and(u8::is_ascii_alphabetic) {
            while self
                .input
                .get(self.position)
                .is_some_and(|c| c.is_ascii_alphanumeric() || *c == b'_')
            {
                self.position += 1;
            }
            let name = std::str::from_utf8(&self.input[start..self.position])
                .map_err(|_| error("invalid variable"))?;
            if name == "I" {
                return Ok(vec![G::imaginary_unit()]);
            }
            if self.variable.as_deref().is_some_and(|v| v != name) {
                return Err(error("multivariate expressions are not supported"));
            }
            self.variable = Some(name.into());
            if self.limits.max_degree == 0 {
                return Err(error("variable exceeds max-degree"));
            }
            return Ok(vec![G::zero(), G::one()]);
        }
        Err(error(
            "expected integer, I, variable or parenthesis; floating inputs are not accepted",
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const LIMITS: ParseLimits = ParseLimits {
        max_degree: 20,
        max_coefficient_bits: 128,
        max_literal_bytes: 16384,
    };
    #[test]
    fn exact_syntax_and_rejections() {
        let p = parse_univariate("t^2 + (1+2I) t - I", LIMITS).unwrap();
        assert_eq!(p, parse_univariate("[-I, 1+2I, 1]", LIMITS).unwrap());
        assert_eq!(
            parse_gaussian("(1+I)/(1-I)", LIMITS).unwrap(),
            G::imaginary_unit()
        );
        assert_eq!(
            parse_univariate("(t+1)^3", LIMITS).unwrap(),
            [1, 3, 3, 1].map(G::from)
        );
        for s in [
            "1.2",
            "1e3",
            "x+y",
            "1/(t-t+1)",
            "t^-1",
            "t^21",
            "1/0",
            "1,,2",
            "I^2+t/0",
            "1/(x-x+1)",
        ] {
            assert!(parse_univariate(s, LIMITS).is_err(), "{s}");
        }
        assert!(parse_gaussian("t-t", LIMITS).is_err());
        assert!(parse_univariate(&"-".repeat(100), LIMITS).is_err());
    }

    #[test]
    fn bit_degree_and_literal_caps_are_checked_before_counting() {
        let small = ParseLimits {
            max_degree: 1,
            max_coefficient_bits: 3,
            max_literal_bytes: 16,
        };
        assert!(parse_univariate("7/3 + (1/7)I*t", small).is_ok());
        for input in ["8", "1/8", "I/8", "4+4", "4*4", "(t+1)^2"] {
            assert!(parse_univariate(input, small).is_err(), "{input}");
        }
        assert!(parse_rational("8/8", small).is_err());
        assert!(parse_gaussian("1".repeat(17).as_str(), small).is_err());
        let constants = ParseLimits {
            max_degree: 0,
            ..small
        };
        assert_eq!(parse_gaussian("I^2", constants).unwrap(), G::from(-1));
        assert!(parse_univariate("t", constants).is_err());
        let coefficients = parse_univariate("0, 0, 0", LIMITS).unwrap();
        assert!(coefficients.is_empty());
    }
}
