//! Products of powers of linear forms in `(n, k)`: the shift ratios of a
//! hypergeometric term, their parser, composition and exact evaluation.

use super::super::super::{BigRational as Q, RecurrenceError};
use crate::recurrence::implicit::algebra::Biv;
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};
use std::collections::BTreeMap;

/// The linear form `index n + column k + constant`.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct Linear {
    pub(super) index: Q,
    pub(super) column: Q,
    pub(super) constant: Q,
}

impl Linear {
    pub(super) fn constant(c: Q) -> Self {
        Self {
            index: Q::zero(),
            column: Q::zero(),
            constant: c,
        }
    }

    pub(super) fn index_variable() -> Self {
        Self {
            index: Q::one(),
            column: Q::zero(),
            constant: Q::zero(),
        }
    }

    pub(super) fn column_variable() -> Self {
        Self {
            index: Q::zero(),
            column: Q::one(),
            constant: Q::zero(),
        }
    }

    pub(super) fn is_constant(&self) -> bool {
        self.index.is_zero() && self.column.is_zero()
    }

    pub(super) fn add(&self, other: &Self) -> Self {
        Self {
            index: &self.index + &other.index,
            column: &self.column + &other.column,
            constant: &self.constant + &other.constant,
        }
    }

    pub(super) fn scale(&self, c: &Q) -> Self {
        Self {
            index: &self.index * c,
            column: &self.column * c,
            constant: &self.constant * c,
        }
    }

    /// The value at `(n, k)`.
    pub(super) fn evaluate(&self, n: &Q, k: &Q) -> Q {
        &self.index * n + &self.column * k + &self.constant
    }

    /// `phi(n + shift, k')` with `k'` the linear form `column_map(n, k)`.
    pub(super) fn compose(&self, shift: i64, column_map: &Linear) -> Self {
        let shifted = Self {
            index: self.index.clone(),
            column: Q::zero(),
            constant: &self.constant + &self.index * Q::from_integer(BigInt::from(shift)),
        };
        shifted.add(&column_map.scale(&self.column))
    }

    /// `(c, phi / c)` with `phi / c` monic in its first nonzero variable
    /// coefficient (`index`, then `column`); `None` for a constant form.
    fn normalized(&self) -> Option<(Q, Self)> {
        let lead = if !self.index.is_zero() {
            self.index.clone()
        } else if !self.column.is_zero() {
            self.column.clone()
        } else {
            return None;
        };
        let inverse = Q::one() / &lead;
        Some((lead, self.scale(&inverse)))
    }

    /// The form as a bivariate polynomial in `(n, k)`.
    pub(super) fn to_biv(&self) -> Biv {
        Biv::new(vec![
            vec![self.constant.clone(), self.column.clone()],
            vec![self.index.clone()],
        ])
    }

    pub(super) fn text(&self) -> String {
        let mut parts = Vec::new();
        for (c, name) in [(&self.index, "n"), (&self.column, "k")] {
            if c.is_zero() {
                continue;
            }
            let coefficient = if c.is_one() {
                String::new()
            } else if *c == -Q::one() {
                "-".into()
            } else {
                format!("{c}*")
            };
            parts.push(format!("{coefficient}{name}"));
        }
        if !self.constant.is_zero() || parts.is_empty() {
            parts.push(self.constant.to_string());
        }
        let mut out = parts[0].clone();
        for part in &parts[1..] {
            if let Some(rest) = part.strip_prefix('-') {
                out.push_str(" - ");
                out.push_str(rest);
            } else {
                out.push_str(" + ");
                out.push_str(part);
            }
        }
        out
    }
}

/// `constant * prod form^exponent` with monic, pairwise distinct forms
/// (unique factorization: two values are equal exactly when the rational
/// functions are).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Factored {
    pub(super) constant: Q,
    pub(super) factors: BTreeMap<Linear, i64>,
}

impl Factored {
    pub(super) fn constant(c: Q) -> Self {
        Self {
            constant: c,
            factors: BTreeMap::new(),
        }
    }

    pub(super) fn one() -> Self {
        Self::constant(Q::one())
    }

    /// `phi^exponent`; a constant form becomes a constant.
    pub(super) fn power_of(phi: &Linear, exponent: i64) -> Result<Self, RecurrenceError> {
        match phi.normalized() {
            Some((lead, monic)) => {
                let mut out = Self::constant(power(&lead, exponent)?);
                if exponent != 0 {
                    out.factors.insert(monic, exponent);
                }
                Ok(out)
            }
            None => Ok(Self::constant(power(&phi.constant, exponent)?)),
        }
    }

    pub(super) fn multiply(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.constant *= &other.constant;
        for (phi, e) in &other.factors {
            let entry = out.factors.entry(phi.clone()).or_insert(0);
            *entry += e;
            if *entry == 0 {
                out.factors.remove(phi);
            }
        }
        out
    }

    pub(super) fn inverse(&self) -> Result<Self, RecurrenceError> {
        if self.constant.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "a ratio is identically zero".into(),
            ));
        }
        Ok(Self {
            constant: Q::one() / &self.constant,
            factors: self.factors.iter().map(|(p, e)| (p.clone(), -e)).collect(),
        })
    }

    pub(super) fn power(&self, exponent: i64) -> Result<Self, RecurrenceError> {
        let mut out = Self::constant(power(&self.constant, exponent)?);
        for (phi, e) in &self.factors {
            if *e * exponent != 0 {
                out.factors.insert(phi.clone(), e * exponent);
            }
        }
        Ok(out)
    }

    /// Substitute `n -> n + shift`, `k -> column_map(n, k)` in every factor.
    pub(super) fn compose(&self, shift: i64, column_map: &Linear) -> Result<Self, RecurrenceError> {
        let mut out = Self::constant(self.constant.clone());
        for (phi, e) in &self.factors {
            let image = phi.compose(shift, column_map);
            if image.is_constant() && image.constant.is_zero() {
                return Err(RecurrenceError::Unsupported(format!(
                    "the factor {} vanishes identically after a substitution",
                    phi.text()
                )));
            }
            out = out.multiply(&Self::power_of(&image, *e)?);
        }
        Ok(out)
    }

    /// The exact value at `(n, k)`; `None` when a denominator factor
    /// vanishes there, zero when a numerator factor does.
    pub(super) fn evaluate(&self, n: &Q, k: &Q) -> Option<Q> {
        let mut value = self.constant.clone();
        for (phi, e) in &self.factors {
            let v = phi.evaluate(n, k);
            if v.is_zero() {
                if *e < 0 {
                    return None;
                }
                return Some(Q::zero());
            }
            let p = power(&v, *e).ok()?;
            value *= p;
        }
        Some(value)
    }

    /// The denominator exponents: `-e` for every factor with `e < 0`.
    pub(super) fn denominator_exponents(&self) -> BTreeMap<Linear, i64> {
        self.factors
            .iter()
            .filter(|(_, e)| **e < 0)
            .map(|(p, e)| (p.clone(), -e))
            .collect()
    }

    /// `self * prod phi^(d_phi)` as a polynomial in `(n, k)`, where
    /// `denominator` covers every negative exponent of `self`.
    pub(super) fn numerator_over(&self, denominator: &BTreeMap<Linear, i64>) -> Biv {
        let mut exponents: BTreeMap<&Linear, i64> = BTreeMap::new();
        for (phi, e) in &self.factors {
            *exponents.entry(phi).or_insert(0) += e;
        }
        for (phi, d) in denominator {
            *exponents.entry(phi).or_insert(0) += d;
        }
        let mut out = Biv::constant(self.constant.clone());
        for (phi, e) in exponents {
            debug_assert!(e >= 0);
            if e > 0 {
                out = out.multiply(&phi.to_biv().power(e as usize));
            }
        }
        out
    }

    pub(super) fn text(&self) -> String {
        let mut numerator = Vec::new();
        let mut denominator = Vec::new();
        for (phi, e) in &self.factors {
            let base = if phi.text().contains(' ') {
                format!("({})", phi.text())
            } else {
                phi.text()
            };
            let piece = if e.abs() == 1 {
                base
            } else {
                format!("{base}^{}", e.abs())
            };
            if *e > 0 {
                numerator.push(piece);
            } else {
                denominator.push(piece);
            }
        }
        let mut top = if numerator.is_empty() {
            self.constant.to_string()
        } else if self.constant.is_one() {
            numerator.join("*")
        } else {
            format!("{}*{}", self.constant, numerator.join("*"))
        };
        if !denominator.is_empty() {
            if numerator.is_empty() && !self.constant.is_one() {
                top = self.constant.to_string();
            }
            top = format!("{top}/({})", denominator.join("*"));
        }
        top
    }
}

fn power(value: &Q, exponent: i64) -> Result<Q, RecurrenceError> {
    if exponent < 0 && value.is_zero() {
        return Err(RecurrenceError::InvalidInput(
            "division by zero in a ratio".into(),
        ));
    }
    let magnitude = u32::try_from(exponent.unsigned_abs())
        .map_err(|_| RecurrenceError::LimitExceeded("exponent too large".into()))?;
    let p = num_traits::pow::Pow::pow(value, magnitude);
    Ok(if exponent < 0 { Q::one() / p } else { p })
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// A parsed value: a linear form while it is one, a product otherwise.
enum Value {
    Linear(Linear),
    Product(Factored),
}

impl Value {
    fn into_product(self) -> Result<Factored, RecurrenceError> {
        match self {
            Value::Linear(phi) => Factored::power_of(&phi, 1),
            Value::Product(p) => Ok(p),
        }
    }

    fn into_linear(self) -> Option<Linear> {
        match self {
            Value::Linear(phi) => Some(phi),
            Value::Product(p) => {
                if p.factors.is_empty() {
                    return Some(Linear::constant(p.constant));
                }
                if p.factors.len() == 1 {
                    let (phi, e) = p.factors.iter().next()?;
                    if *e == 1 {
                        return Some(phi.scale(&p.constant));
                    }
                }
                None
            }
        }
    }
}

struct Parser<'a> {
    text: &'a str,
    bytes: &'a [u8],
    position: usize,
}

impl<'a> Parser<'a> {
    fn error(&self, message: &str) -> RecurrenceError {
        RecurrenceError::InvalidInput(format!(
            "ratio `{}`: {message} at position {}",
            self.text, self.position
        ))
    }

    fn skip_spaces(&mut self) {
        while self.position < self.bytes.len() && self.bytes[self.position].is_ascii_whitespace() {
            self.position += 1;
        }
    }

    fn peek(&mut self) -> Option<u8> {
        self.skip_spaces();
        self.bytes.get(self.position).copied()
    }

    fn sum(&mut self) -> Result<Value, RecurrenceError> {
        let mut value = self.product()?;
        while let Some(op @ (b'+' | b'-')) = self.peek() {
            self.position += 1;
            let right = self.product()?;
            let (Some(a), Some(b)) = (value.into_linear(), right.into_linear()) else {
                return Err(self.error("a sum must be linear in n and k"));
            };
            let b = if op == b'-' { b.scale(&-Q::one()) } else { b };
            value = Value::Linear(a.add(&b));
        }
        Ok(value)
    }

    fn product(&mut self) -> Result<Value, RecurrenceError> {
        let mut value = self.unary()?;
        loop {
            match self.peek() {
                Some(b'*') => {
                    self.position += 1;
                    let right = self.unary()?;
                    value = multiply(value, right)?;
                }
                Some(b'/') => {
                    self.position += 1;
                    let right = self.unary()?;
                    value = multiply(value, Value::Product(right.into_product()?.inverse()?))?;
                }
                // Implicit multiplication: `2n`, `2(n+1)`, `(n+1)(k+2)`.
                Some(b'n' | b'k' | b'(' | b'0'..=b'9') => {
                    let right = self.unary()?;
                    value = multiply(value, right)?;
                }
                _ => return Ok(value),
            }
        }
    }

    fn unary(&mut self) -> Result<Value, RecurrenceError> {
        if self.peek() == Some(b'-') {
            self.position += 1;
            let inner = self.unary()?;
            return multiply(Value::Linear(Linear::constant(-Q::one())), inner);
        }
        if self.peek() == Some(b'+') {
            self.position += 1;
            return self.unary();
        }
        let base = self.primary()?;
        if self.peek() == Some(b'^') {
            self.position += 1;
            let negative = if self.peek() == Some(b'-') {
                self.position += 1;
                true
            } else {
                false
            };
            let exponent = self.integer()?;
            let exponent = i64::try_from(exponent).map_err(|_| self.error("exponent too large"))?;
            let exponent = if negative { -exponent } else { exponent };
            return Ok(Value::Product(base.into_product()?.power(exponent)?));
        }
        Ok(base)
    }

    fn integer(&mut self) -> Result<BigInt, RecurrenceError> {
        self.skip_spaces();
        let start = self.position;
        while self.position < self.bytes.len() && self.bytes[self.position].is_ascii_digit() {
            self.position += 1;
        }
        if start == self.position {
            return Err(self.error("expected an integer"));
        }
        self.text[start..self.position]
            .parse::<BigInt>()
            .map_err(|_| self.error("invalid integer"))
    }

    fn primary(&mut self) -> Result<Value, RecurrenceError> {
        match self.peek() {
            Some(b'(') => {
                self.position += 1;
                let inner = self.sum()?;
                if self.peek() != Some(b')') {
                    return Err(self.error("expected `)`"));
                }
                self.position += 1;
                Ok(inner)
            }
            Some(b'n') => {
                self.position += 1;
                Ok(Value::Linear(Linear::index_variable()))
            }
            Some(b'k') => {
                self.position += 1;
                Ok(Value::Linear(Linear::column_variable()))
            }
            Some(b'0'..=b'9') => {
                let value = self.integer()?;
                Ok(Value::Linear(Linear::constant(Q::from_integer(value))))
            }
            _ => Err(self.error("expected a number, `n`, `k` or `(`")),
        }
    }
}

fn multiply(a: Value, b: Value) -> Result<Value, RecurrenceError> {
    // Keep a linear form linear under multiplication by a constant.
    if let (Value::Linear(x), Value::Linear(y)) = (&a, &b) {
        if x.is_constant() {
            return Ok(Value::Linear(y.scale(&x.constant)));
        }
        if y.is_constant() {
            return Ok(Value::Linear(x.scale(&y.constant)));
        }
    }
    Ok(Value::Product(
        a.into_product()?.multiply(&b.into_product()?),
    ))
}

/// Parse a product of integer powers of linear forms in `n` and `k` with
/// rational constants, for example `(n-k)^2/((k+1)*(2n-2k+1))`.
pub(super) fn parse_factored(text: &str) -> Result<Factored, RecurrenceError> {
    let mut parser = Parser {
        text,
        bytes: text.as_bytes(),
        position: 0,
    };
    let value = parser.sum()?;
    if parser.peek().is_some() {
        return Err(parser.error("unexpected character"));
    }
    let product = value.into_product()?;
    if product.constant.is_zero() {
        return Err(RecurrenceError::InvalidInput(format!(
            "ratio `{text}` is identically zero"
        )));
    }
    Ok(product)
}

/// Parse a linear form (sums of rational multiples of `n`, `k` and
/// constants), for example `(n+1)/2`.
pub(super) fn parse_linear(text: &str) -> Result<Linear, RecurrenceError> {
    let mut parser = Parser {
        text,
        bytes: text.as_bytes(),
        position: 0,
    };
    let value = parser.sum()?;
    if parser.peek().is_some() {
        return Err(parser.error("unexpected character"));
    }
    value
        .into_linear()
        .ok_or_else(|| RecurrenceError::InvalidInput(format!("`{text}` is not linear in n and k")))
}

// ---------------------------------------------------------------------------
// Polynomials in (n, k)
// ---------------------------------------------------------------------------

struct PolynomialParser<'a> {
    text: &'a str,
    bytes: &'a [u8],
    position: usize,
}

impl PolynomialParser<'_> {
    fn error(&self, message: &str) -> RecurrenceError {
        RecurrenceError::InvalidInput(format!(
            "polynomial `{}`: {message} at position {}",
            self.text, self.position
        ))
    }

    fn peek(&mut self) -> Option<u8> {
        while self.position < self.bytes.len() && self.bytes[self.position].is_ascii_whitespace() {
            self.position += 1;
        }
        self.bytes.get(self.position).copied()
    }

    fn integer(&mut self) -> Result<BigInt, RecurrenceError> {
        self.peek();
        let start = self.position;
        while self.position < self.bytes.len() && self.bytes[self.position].is_ascii_digit() {
            self.position += 1;
        }
        if start == self.position {
            return Err(self.error("expected an integer"));
        }
        self.text[start..self.position]
            .parse::<BigInt>()
            .map_err(|_| self.error("invalid integer"))
    }

    fn sum(&mut self) -> Result<Biv, RecurrenceError> {
        let mut value = self.product()?;
        while let Some(op @ (b'+' | b'-')) = self.peek() {
            self.position += 1;
            let right = self.product()?;
            value = if op == b'+' {
                value.add(&right)
            } else {
                value.subtract(&right)
            };
        }
        Ok(value)
    }

    fn product(&mut self) -> Result<Biv, RecurrenceError> {
        let mut value = self.unary()?;
        loop {
            match self.peek() {
                Some(b'*') => {
                    self.position += 1;
                    value = value.multiply(&self.unary()?);
                }
                Some(b'/') => {
                    self.position += 1;
                    let divisor = self.integer()?;
                    if divisor.is_zero() {
                        return Err(self.error("division by zero"));
                    }
                    value = value.scale(&Q::new(BigInt::one(), divisor));
                }
                Some(b'n' | b'k' | b'(' | b'0'..=b'9') => {
                    value = value.multiply(&self.unary()?);
                }
                _ => return Ok(value),
            }
        }
    }

    fn unary(&mut self) -> Result<Biv, RecurrenceError> {
        if self.peek() == Some(b'-') {
            self.position += 1;
            return Ok(self.unary()?.scale(&-Q::one()));
        }
        let base = match self.peek() {
            Some(b'(') => {
                self.position += 1;
                let inner = self.sum()?;
                if self.peek() != Some(b')') {
                    return Err(self.error("expected `)`"));
                }
                self.position += 1;
                inner
            }
            Some(b'n') => {
                self.position += 1;
                Biv::monomial(1, 0, Q::one())
            }
            Some(b'k') => {
                self.position += 1;
                Biv::monomial(0, 1, Q::one())
            }
            Some(b'0'..=b'9') => Biv::constant(Q::from_integer(self.integer()?)),
            _ => return Err(self.error("expected a number, `n`, `k` or `(`")),
        };
        if self.peek() == Some(b'^') {
            self.position += 1;
            let exponent = self
                .integer()?
                .to_usize()
                .filter(|e| *e <= 64)
                .ok_or_else(|| self.error("exponent out of range"))?;
            return Ok(base.power(exponent));
        }
        Ok(base)
    }
}

/// Parse a polynomial in `n` and `k` with rational coefficients (integer
/// divisors only), for example `2n*k + n - 2k^2`.
pub(super) fn parse_polynomial(text: &str) -> Result<Biv, RecurrenceError> {
    let mut parser = PolynomialParser {
        text,
        bytes: text.as_bytes(),
        position: 0,
    };
    let value = parser.sum()?;
    if parser.peek().is_some() {
        return Err(parser.error("unexpected character"));
    }
    Ok(value)
}

/// `M(n + shift, column_map(n, k))` for a polynomial `M` in `(n, k)`.
pub(super) fn compose_polynomial(m: &Biv, shift: i64, column_map: &Linear) -> Biv {
    let index = Linear::index_variable()
        .add(&Linear::constant(Q::from_integer(shift.into())))
        .to_biv();
    let column = column_map.to_biv();
    let mut out = Biv::zero();
    for (a, row) in m.rows().iter().enumerate() {
        let index_power = index.power(a);
        for (b, c) in row.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            out = out.add(&index_power.multiply(&column.power(b)).scale(c));
        }
    }
    out
}
