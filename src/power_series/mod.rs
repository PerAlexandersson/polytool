//! Exact rational series modulo a specified power of the variable.
//!
//! [`TruncatedSeries`](crate::power_series::TruncatedSeries) stores ordinary coefficients, including known zeros.
//! Precision is a coefficient count: `order=4` means modulo `z^4`.
//! Binary operations use the smaller precision and never invent missing terms.
//! [`SeriesExpression`](crate::power_series::SeriesExpression) parses bounded expressions and evaluates only the
//! requested prefix (with extra terms when differentiation requires them).
//!
//! ```
//! use polytool::power_series::SeriesExpression;
//! let series = SeriesExpression::parse("exp(z)-1-z")?.evaluate(6)?;
//! assert_eq!(series.order(), 7);
//! assert_eq!(series.coefficients()[2].to_string(), "1/2");
//! # Ok::<(), polytool::power_series::SeriesError>(())
//! ```

mod expression;
#[cfg(test)]
mod tests;
pub(crate) use expression::Expr;
pub use expression::{SeriesExpression, SeriesExpressionLimits};

use num_rational::BigRational as Q;
use num_traits::{One, Zero};

/// Invalid series operation, expression, or resource bound.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SeriesError(pub String);
impl std::fmt::Display for SeriesError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for SeriesError {}

/// Ordinary rational coefficients modulo `z^order`, in ascending degree order.
/// Empty coefficients mean no known terms. Trailing known zeros are preserved.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TruncatedSeries {
    coefficients: Vec<Q>,
}
impl TruncatedSeries {
    /// Store exactly these known ordinary coefficients, including constant term.
    pub fn new(coefficients: Vec<Q>) -> Self {
        Self { coefficients }
    }
    /// Known coefficient count; the series is specified modulo `z^order`.
    pub fn order(&self) -> usize {
        self.coefficients.len()
    }
    /// Known ordinary coefficients, ascending from degree zero.
    pub fn coefficients(&self) -> &[Q] {
        &self.coefficients
    }
    /// Reduce precision to at most `order`; never extend it.
    pub fn truncate(&self, order: usize) -> Self {
        Self::new(self.coefficients.iter().take(order).cloned().collect())
    }
    /// Add at the smaller input precision.
    pub fn add(&self, other: &Self) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .zip(&other.coefficients)
                .map(|(a, b)| a + b)
                .collect(),
        )
    }
    /// Subtract at the smaller input precision.
    pub fn subtract(&self, other: &Self) -> Self {
        self.add(&other.scale(&-Q::one()))
    }
    /// Multiply every known coefficient by an exact rational scalar.
    pub fn scale(&self, scalar: &Q) -> Self {
        let mut coefficients = crate::vec_poly::scale(&self.coefficients, scalar);
        coefficients.resize(self.order(), Q::zero());
        coefficients.truncate(self.order());
        Self::new(coefficients)
    }
    /// Multiply modulo the smaller input order, discarding higher terms early.
    pub fn multiply(&self, other: &Self) -> Self {
        let order = self.order().min(other.order());
        let mut out = vec![Q::zero(); order];
        for (i, a) in self
            .coefficients
            .iter()
            .take(order)
            .enumerate()
            .filter(|(_, a)| !a.is_zero())
        {
            for (j, b) in other
                .coefficients
                .iter()
                .take(order - i)
                .enumerate()
                .filter(|(_, b)| !b.is_zero())
            {
                out[i + j] += a * b;
            }
        }
        Self::new(out)
    }
    fn constant(value: Q, order: usize) -> Self {
        let mut out = vec![Q::zero(); order];
        if let Some(c) = out.first_mut() {
            *c = value;
        }
        Self::new(out)
    }
    /// Nonnegative integer power, at the original precision (`f^0=1`).
    pub fn power(&self, mut exponent: usize) -> Self {
        let mut result = Self::constant(Q::one(), self.order());
        let mut base = self.clone();
        while exponent != 0 {
            if exponent % 2 == 1 {
                result = result.multiply(&base);
            }
            exponent /= 2;
            if exponent != 0 {
                base = base.multiply(&base);
            }
        }
        result
    }
    /// Formal derivative; loses one known coefficient of precision.
    pub fn derivative(&self) -> Self {
        Self::new(
            self.coefficients
                .iter()
                .enumerate()
                .skip(1)
                .map(|(i, c)| c * Q::from_integer(i.into()))
                .collect(),
        )
    }
    /// Multiplicative inverse, requiring a known nonzero constant coefficient.
    pub fn inverse(&self) -> Result<Self, SeriesError> {
        let c = self
            .coefficients
            .first()
            .filter(|c| !c.is_zero())
            .ok_or_else(|| SeriesError("inverse requires a known nonzero constant term".into()))?;
        let mut out = vec![Q::zero(); self.order()];
        out[0] = c.recip();
        for n in 1..self.order() {
            out[n] = -(1..=n)
                .map(|i| &self.coefficients[i] * &out[n - i])
                .sum::<Q>()
                / c;
        }
        Ok(Self::new(out))
    }
    /// Exact formal exponential, requiring constant term zero (`exp(0)=1`).
    /// A nonzero rational constant generally has a nonrational exponential.
    pub fn exponential(&self) -> Result<Self, SeriesError> {
        if self.coefficients.first() != Some(&Q::zero()) {
            return Err(SeriesError(
                "exp requires a known zero constant term".into(),
            ));
        }
        let mut out = vec![Q::zero(); self.order()];
        out[0] = Q::one();
        for n in 1..self.order() {
            out[n] = (1..=n)
                .map(|i| Q::from_integer(i.into()) * &self.coefficients[i] * &out[n - i])
                .sum::<Q>()
                / Q::from_integer(n.into());
        }
        Ok(Self::new(out))
    }
    /// Formal logarithm, requiring constant term one (`log(1)=0`).
    pub fn logarithm(&self) -> Result<Self, SeriesError> {
        if self.coefficients.first() != Some(&Q::one()) {
            return Err(SeriesError(
                "log requires a known constant term of one".into(),
            ));
        }
        let derivative = self.derivative().multiply(&self.inverse()?);
        let mut out = vec![Q::zero(); self.order()];
        for (i, c) in derivative.coefficients.iter().enumerate() {
            out[i + 1] = c / Q::from_integer((i + 1).into());
        }
        Ok(Self::new(out))
    }
    /// Compose modulo the smaller input order; inner constant must be zero.
    pub fn compose(&self, inner: &Self) -> Result<Self, SeriesError> {
        if inner.coefficients.first() != Some(&Q::zero()) {
            return Err(SeriesError(
                "composition requires a known zero inner constant term".into(),
            ));
        }
        let order = self.order().min(inner.order());
        let mut out = Self::constant(Q::zero(), order);
        for c in self.coefficients.iter().take(order).rev() {
            out = out.multiply(inner);
            if let Some(first) = out.coefficients.first_mut() {
                *first += c;
            }
        }
        Ok(out)
    }
}

impl std::ops::Add for &TruncatedSeries {
    type Output = TruncatedSeries;
    fn add(self, other: Self) -> Self::Output {
        TruncatedSeries::add(self, other)
    }
}
impl std::ops::Sub for &TruncatedSeries {
    type Output = TruncatedSeries;
    fn sub(self, other: Self) -> Self::Output {
        self.subtract(other)
    }
}
impl std::ops::Mul for &TruncatedSeries {
    type Output = TruncatedSeries;
    fn mul(self, other: Self) -> Self::Output {
        self.multiply(other)
    }
}
impl std::ops::Mul<&Q> for &TruncatedSeries {
    type Output = TruncatedSeries;
    fn mul(self, scalar: &Q) -> Self::Output {
        self.scale(scalar)
    }
}
impl std::ops::Neg for &TruncatedSeries {
    type Output = TruncatedSeries;
    fn neg(self) -> Self::Output {
        self.scale(&-Q::one())
    }
}
