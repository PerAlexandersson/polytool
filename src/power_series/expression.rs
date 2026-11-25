//! Bounded recursive-descent grammar and precision-aware evaluation.
use super::{SeriesError, TruncatedSeries};
use crate::numeric::factorial_bigint;
use crate::recurrence::parse_rational_coeff;
use num_rational::BigRational as Q;
use num_traits::{One, Zero};

/// Bounds checked before evaluation and after each expression operation.
#[derive(Clone, Debug)]
pub struct SeriesExpressionLimits {
    /// Maximum input bytes.
    pub bytes: usize,
    /// Maximum expression tree nodes.
    pub nodes: usize,
    /// Maximum tree depth, including unary operators.
    pub depth: usize,
    /// Maximum integer exponent or explicit truncation degree.
    pub exponent: usize,
    /// Maximum requested degree, including derivative lookahead.
    pub degree: usize,
    /// Maximum numerator or denominator bits in literals and intermediate results.
    pub coefficient_bits: u64,
}
impl Default for SeriesExpressionLimits {
    fn default() -> Self {
        Self {
            bytes: 8192,
            nodes: 512,
            depth: 64,
            exponent: 1024,
            degree: 4096,
            coefficient_bits: 65536,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Expr {
    Polynomial(Vec<Q>),
    Variable,
    Add(Box<Expr>, Box<Expr>),
    Subtract(Box<Expr>, Box<Expr>),
    Multiply(Box<Expr>, Box<Expr>),
    Divide(Box<Expr>, Box<Expr>),
    Negate(Box<Expr>),
    Power(Box<Expr>, usize),
    Exp(Box<Expr>),
    Log(Box<Expr>),
    Inverse(Box<Expr>),
    Compose(Box<Expr>, Box<Expr>),
    Truncate(Box<Expr>, usize),
    Derivative(Box<Expr>),
}
/// Parsed formal series expression, independent of evaluation precision.
///
/// Variables `x` and `z` are aliases. Functions: `exp`, `log`, `inv`,
/// `compose(outer,inner)`, `trunc(series,degree)`, `derivative(series)`.
/// Inline `coeffs:[c0,c1,...]` and `coeffs:exponential:[w0,w1,...]`
/// are finite polynomials, with an explicit degree-zero entry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeriesExpression {
    pub(crate) root: Expr,
    bytes: usize,
}
impl SeriesExpression {
    /// Parse using [`SeriesExpressionLimits::default`]. Errors identify a byte offset.
    pub fn parse(text: &str) -> Result<Self, SeriesError> {
        Self::parse_with_limits(text, &SeriesExpressionLimits::default())
    }
    /// Parse and validate with caller-supplied bounds before any series expansion.
    pub fn parse_with_limits(
        text: &str,
        limits: &SeriesExpressionLimits,
    ) -> Result<Self, SeriesError> {
        if text.len() > limits.bytes {
            return Err(SeriesError(format!(
                "expression exceeds {} bytes",
                limits.bytes
            )));
        }
        let mut parser = Parser {
            text,
            at: 0,
            nodes: 0,
            limits,
        };
        let root = parser.sum(0)?;
        parser.space();
        if parser.at != text.len() {
            return parser.error("unexpected trailing input");
        }
        let expression = Self {
            root,
            bytes: text.len(),
        };
        expression.validate_limits(limits)?;
        Ok(expression)
    }
    /// Recheck a parsed tree against bounds (including literal sizes).
    pub fn validate_limits(&self, limits: &SeriesExpressionLimits) -> Result<(), SeriesError> {
        if self.bytes > limits.bytes {
            return Err(SeriesError("expression byte limit exceeded".into()));
        }
        let mut nodes = 0;
        self.root.validate(limits, 0, &mut nodes)
    }
    /// Evaluate through inclusive `degree`, expanding no unused higher terms.
    pub fn evaluate(&self, degree: usize) -> Result<TruncatedSeries, SeriesError> {
        self.evaluate_with_limits(degree, &SeriesExpressionLimits::default())
    }
    /// Evaluate exactly with bounded intermediate precision and coefficient sizes.
    pub fn evaluate_with_limits(
        &self,
        degree: usize,
        limits: &SeriesExpressionLimits,
    ) -> Result<TruncatedSeries, SeriesError> {
        self.validate_limits(limits)?;
        self.root.evaluate(degree, limits)
    }
}
impl Expr {
    fn children(&self) -> Vec<&Self> {
        match self {
            Self::Add(a, b)
            | Self::Subtract(a, b)
            | Self::Multiply(a, b)
            | Self::Divide(a, b)
            | Self::Compose(a, b) => vec![a, b],
            Self::Negate(a)
            | Self::Power(a, _)
            | Self::Exp(a)
            | Self::Log(a)
            | Self::Inverse(a)
            | Self::Truncate(a, _)
            | Self::Derivative(a) => vec![a],
            _ => vec![],
        }
    }
    fn validate(
        &self,
        limits: &SeriesExpressionLimits,
        depth: usize,
        nodes: &mut usize,
    ) -> Result<(), SeriesError> {
        *nodes += 1;
        if depth > limits.depth || *nodes > limits.nodes {
            return Err(SeriesError(
                "expression depth or node limit exceeded".into(),
            ));
        }
        if let Self::Power(_, n) | Self::Truncate(_, n) = self {
            if *n > limits.exponent {
                return Err(SeriesError(format!(
                    "exponent/truncation degree exceeds {}",
                    limits.exponent
                )));
            }
        }
        if let Self::Polynomial(p) = self {
            if p.len() > limits.exponent.saturating_add(1) {
                return Err(SeriesError("coefficient list length limit exceeded".into()));
            }
            check_bits(p, limits)?;
        }
        for child in self.children() {
            child.validate(limits, depth + 1, nodes)?;
        }
        Ok(())
    }
    pub(crate) fn evaluate(
        &self,
        degree: usize,
        limits: &SeriesExpressionLimits,
    ) -> Result<TruncatedSeries, SeriesError> {
        if degree > limits.degree || degree == usize::MAX {
            return Err(SeriesError(format!(
                "evaluation degree exceeds {}",
                limits.degree
            )));
        }
        let order = degree + 1;
        let eval = |a: &Self| a.evaluate(degree, limits);
        let out = match self {
            Self::Polynomial(p) => {
                let mut p = p.clone();
                p.resize(order, Q::zero());
                p.truncate(order);
                TruncatedSeries::new(p)
            }
            Self::Variable => {
                let mut p = vec![Q::zero(); order];
                if degree > 0 {
                    p[1] = Q::one();
                }
                TruncatedSeries::new(p)
            }
            Self::Add(a, b) => eval(a)?.add(&eval(b)?),
            Self::Subtract(a, b) => eval(a)?.subtract(&eval(b)?),
            Self::Multiply(a, b) => eval(a)?.multiply(&eval(b)?),
            Self::Divide(a, b) => eval(a)?.multiply(&eval(b)?.inverse()?),
            Self::Negate(a) => eval(a)?.scale(&-Q::one()),
            Self::Power(a, n) => {
                // Bound every squaring, rather than only the final power.
                let mut base = eval(a)?;
                let mut out = TruncatedSeries::constant(Q::one(), order);
                let mut n = *n;
                while n > 0 {
                    if n % 2 == 1 {
                        out = out.multiply(&base);
                        check_bits(out.coefficients(), limits)?;
                    }
                    n /= 2;
                    if n > 0 {
                        base = base.multiply(&base);
                        check_bits(base.coefficients(), limits)?;
                    }
                }
                out
            }
            Self::Exp(a) => eval(a)?.exponential()?,
            Self::Log(a) => eval(a)?.logarithm()?,
            Self::Inverse(a) => eval(a)?.inverse()?,
            Self::Compose(a, b) => eval(a)?.compose(&eval(b)?)?,
            Self::Truncate(a, n) => {
                let mut p = a.evaluate(degree.min(*n), limits)?.coefficients;
                p.resize(order, Q::zero());
                TruncatedSeries::new(p)
            }
            Self::Derivative(a) => a.evaluate(degree + 1, limits)?.derivative(),
        };
        check_bits(out.coefficients(), limits)?;
        Ok(out)
    }
}
fn check_bits(coefficients: &[Q], limits: &SeriesExpressionLimits) -> Result<(), SeriesError> {
    if coefficients.iter().any(|q| {
        q.numer().bits() > limits.coefficient_bits || q.denom().bits() > limits.coefficient_bits
    }) {
        Err(SeriesError(format!(
            "series coefficient exceeds {} bits",
            limits.coefficient_bits
        )))
    } else {
        Ok(())
    }
}
struct Parser<'a> {
    text: &'a str,
    at: usize,
    nodes: usize,
    limits: &'a SeriesExpressionLimits,
}
impl Parser<'_> {
    fn error<T>(&self, message: &str) -> Result<T, SeriesError> {
        Err(SeriesError(format!("at byte {}: {message}", self.at)))
    }
    fn space(&mut self) {
        while self
            .text
            .as_bytes()
            .get(self.at)
            .is_some_and(u8::is_ascii_whitespace)
        {
            self.at += 1;
        }
    }
    fn eat(&mut self, token: &str) -> bool {
        self.space();
        if self.text[self.at..].starts_with(token) {
            self.at += token.len();
            true
        } else {
            false
        }
    }
    fn require(&mut self, token: &str) -> Result<(), SeriesError> {
        if self.eat(token) {
            Ok(())
        } else {
            self.error(&format!("expected '{token}'"))
        }
    }
    fn node(&mut self, expr: Expr) -> Result<Expr, SeriesError> {
        self.nodes += 1;
        if self.nodes > self.limits.nodes {
            self.error("expression node limit exceeded")
        } else {
            Ok(expr)
        }
    }
    fn guard(&self, depth: usize) -> Result<(), SeriesError> {
        if depth > self.limits.depth {
            self.error("expression nesting limit exceeded")
        } else {
            Ok(())
        }
    }
    fn sum(&mut self, depth: usize) -> Result<Expr, SeriesError> {
        self.guard(depth)?;
        let mut a = self.product(depth)?;
        loop {
            if self.eat("+") {
                let b = self.product(depth)?;
                a = self.node(Expr::Add(Box::new(a), Box::new(b)))?;
            } else if self.eat("-") {
                let b = self.product(depth)?;
                a = self.node(Expr::Subtract(Box::new(a), Box::new(b)))?;
            } else {
                return Ok(a);
            }
        }
    }
    fn product(&mut self, depth: usize) -> Result<Expr, SeriesError> {
        let mut a = self.unary(depth)?;
        loop {
            if self.eat("*") {
                let b = self.unary(depth)?;
                a = self.node(Expr::Multiply(Box::new(a), Box::new(b)))?;
            } else if self.eat("/") {
                let b = self.unary(depth)?;
                a = self.node(Expr::Divide(Box::new(a), Box::new(b)))?;
            } else {
                return Ok(a);
            }
        }
    }
    fn unary(&mut self, depth: usize) -> Result<Expr, SeriesError> {
        self.guard(depth)?;
        if self.eat("+") {
            return self.unary(depth + 1);
        }
        if self.eat("-") {
            let a = self.unary(depth + 1)?;
            return self.node(Expr::Negate(Box::new(a)));
        }
        let a = self.atom(depth)?;
        if self.eat("^") {
            let n = self.integer()?;
            self.node(Expr::Power(Box::new(a), n))
        } else {
            Ok(a)
        }
    }
    fn integer(&mut self) -> Result<usize, SeriesError> {
        self.space();
        let start = self.at;
        while self
            .text
            .as_bytes()
            .get(self.at)
            .is_some_and(u8::is_ascii_digit)
        {
            self.at += 1;
        }
        let n = self.text[start..self.at].parse::<usize>().map_err(|_| {
            SeriesError(format!(
                "at byte {start}: expected bounded nonnegative integer"
            ))
        })?;
        if n > self.limits.exponent {
            return self.error("exponent/truncation degree limit exceeded");
        }
        Ok(n)
    }
    fn atom(&mut self, depth: usize) -> Result<Expr, SeriesError> {
        if self.eat("(") {
            let a = self.sum(depth + 1)?;
            self.require(")")?;
            return Ok(a);
        }
        self.space();
        let start = self.at;
        if self.text[self.at..].starts_with("coeffs:") {
            self.at += 7;
            let exponential = self.eat("exponential:");
            if !exponential {
                self.eat("ordinary:");
            }
            self.require("[")?;
            let mut p = Vec::new();
            if !self.eat("]") {
                loop {
                    self.space();
                    let start = self.at;
                    while self.at < self.text.len()
                        && !matches!(self.text.as_bytes()[self.at], b',' | b']')
                    {
                        self.at += 1;
                    }
                    let token = self.text[start..self.at].trim();
                    let raw = if token.starts_with('"') || token.ends_with('"') {
                        token
                            .strip_prefix('"')
                            .and_then(|s| s.strip_suffix('"'))
                            .ok_or_else(|| {
                                SeriesError(format!(
                                    "at byte {start}: expected paired rational quotes"
                                ))
                            })?
                    } else {
                        token
                    };
                    let mut q = parse_rational_coeff(raw)
                        .map_err(|e| SeriesError(format!("at byte {start}: {e}")))?;
                    check_bits(std::slice::from_ref(&q), self.limits)?;
                    if exponential {
                        q /= Q::from_integer(factorial_bigint(p.len()));
                    }
                    p.push(q);
                    if p.len() > self.limits.exponent.saturating_add(1) {
                        return self.error("coefficient list length limit exceeded");
                    }
                    if self.eat("]") {
                        break;
                    }
                    self.require(",")?;
                }
            }
            return self.node(Expr::Polynomial(p));
        }
        while self
            .text
            .as_bytes()
            .get(self.at)
            .is_some_and(u8::is_ascii_digit)
        {
            self.at += 1;
        }
        if self.at > start {
            let q = parse_rational_coeff(&self.text[start..self.at])
                .map_err(|e| SeriesError(format!("at byte {start}: {e}")))?;
            check_bits(std::slice::from_ref(&q), self.limits)?;
            return self.node(Expr::Polynomial(vec![q]));
        }
        while self
            .text
            .as_bytes()
            .get(self.at)
            .is_some_and(|c| c.is_ascii_alphabetic() || *c == b'_')
        {
            self.at += 1;
        }
        let name = &self.text[start..self.at];
        if name == "z" || name == "x" {
            return self.node(Expr::Variable);
        }
        if !matches!(
            name,
            "exp" | "log" | "inv" | "compose" | "trunc" | "derivative"
        ) {
            return self.error("expected number, z, x, coeffs:[...], or a supported function");
        }
        self.require("(")?;
        let a = Box::new(self.sum(depth + 1)?);
        let expr = match name {
            "compose" => {
                self.require(",")?;
                let b = Box::new(self.sum(depth + 1)?);
                Expr::Compose(a, b)
            }
            "trunc" => {
                self.require(",")?;
                Expr::Truncate(a, self.integer()?)
            }
            "exp" => Expr::Exp(a),
            "log" => Expr::Log(a),
            "inv" => Expr::Inverse(a),
            _ => Expr::Derivative(a),
        };
        self.require(")")?;
        self.node(expr)
    }
}
