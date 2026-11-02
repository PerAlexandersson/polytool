//! Bivariate polynomials in the index `n` and a variable `y`.

use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Zero};

/// A bivariate polynomial `sum c_(i,j) n^i y^j` over the rationals:
/// `rows()[i][j]` is the coefficient of `n^i y^j`.
///
/// The representation is kept trimmed (no trailing zeros in a row, no
/// trailing empty rows; zero is the empty vector), so equality is
/// structural.  The variable `y` is the space variable `x`, the Euler
/// operator `theta = x D` or a coefficient index `k`, depending on context.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct BivariatePolynomial(Vec<Vec<Rational>>);

impl BivariatePolynomial {
    /// The polynomial with `rows[i][j]` the coefficient of `n^i y^j`
    /// (ragged rows allowed), trimmed.
    pub fn new(mut rows: Vec<Vec<Rational>>) -> Self {
        for row in &mut rows {
            poly::strip_trailing_zeros(row);
        }
        while rows.last().is_some_and(Vec::is_empty) {
            rows.pop();
        }
        Self(rows)
    }

    /// The trimmed rows: `rows()[i][j]` is the coefficient of `n^i y^j`.
    pub fn rows(&self) -> &[Vec<Rational>] {
        &self.0
    }

    /// The zero polynomial.
    pub fn zero() -> Self {
        Self(Vec::new())
    }

    /// The constant `c`.
    pub fn constant(c: Rational) -> Self {
        Self::new(vec![vec![c]])
    }

    /// The constant one.
    pub fn one() -> Self {
        Self::constant(Rational::one())
    }

    /// `c n^index_power y^variable_power`.
    pub fn monomial(index_power: usize, variable_power: usize, c: Rational) -> Self {
        let mut rows = vec![Vec::new(); index_power + 1];
        let mut row = vec![Rational::zero(); variable_power + 1];
        row[variable_power] = c;
        rows[index_power] = row;
        Self::new(rows)
    }

    /// A polynomial in `y` only (ascending coefficients).
    pub fn from_variable_polynomial(p: &[Rational]) -> Self {
        Self::new(vec![p.to_vec()])
    }

    /// A polynomial in `n` only (ascending coefficients).
    pub fn from_index_polynomial(p: &[Rational]) -> Self {
        Self::new(p.iter().map(|c| vec![c.clone()]).collect())
    }

    /// Whether this is the zero polynomial.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The constant value, when the polynomial is a constant (zero included).
    pub fn as_constant(&self) -> Option<Rational> {
        match self.0.as_slice() {
            [] => Some(Rational::zero()),
            [row] if row.len() == 1 => Some(row[0].clone()),
            _ => None,
        }
    }

    /// Whether the polynomial does not involve `n`.
    pub fn is_index_free(&self) -> bool {
        self.0.len() <= 1
    }

    /// Whether the polynomial does not involve `y`.
    pub fn is_variable_free(&self) -> bool {
        self.0.iter().all(|row| row.len() <= 1)
    }

    /// The sum `self + other`.
    pub fn add(&self, other: &Self) -> Self {
        let len = self.0.len().max(other.0.len());
        Self::new(
            (0..len)
                .map(|i| {
                    let a = self.0.get(i).map_or(&[][..], Vec::as_slice);
                    let b = other.0.get(i).map_or(&[][..], Vec::as_slice);
                    poly::add(a, b)
                })
                .collect(),
        )
    }

    /// The difference `self - other`.
    pub fn subtract(&self, other: &Self) -> Self {
        self.add(&other.scale(&-Rational::one()))
    }

    /// The scalar multiple `c self`.
    pub fn scale(&self, c: &Rational) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        Self::new(
            self.0
                .iter()
                .map(|row| row.iter().map(|x| x * c).collect())
                .collect(),
        )
    }

    /// The product `self other`.
    pub fn multiply(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut rows = vec![Vec::new(); self.0.len() + other.0.len() - 1];
        for (i, left) in self.0.iter().enumerate() {
            if left.is_empty() {
                continue;
            }
            for (j, right) in other.0.iter().enumerate() {
                if right.is_empty() {
                    continue;
                }
                rows[i + j] = poly::add(&rows[i + j], &poly::mul(left, right));
            }
        }
        Self::new(rows)
    }

    /// The power `self^exponent` (one for exponent zero).
    pub fn power(&self, exponent: usize) -> Self {
        (0..exponent).fold(Self::one(), |acc, _| acc.multiply(self))
    }

    /// The translate `p(n + index_shift, y + variable_shift)`.
    pub fn shift(&self, index_shift: &Rational, variable_shift: &Rational) -> Self {
        let n_linear = [index_shift.clone(), Rational::one()];
        let mut n_power = vec![Rational::one()];
        let mut result = Self::zero();
        for (i, row) in self.0.iter().enumerate() {
            if i > 0 {
                n_power = poly::mul(&n_power, &n_linear);
            }
            if row.is_empty() {
                continue;
            }
            let shifted = if variable_shift.is_zero() {
                row.clone()
            } else {
                poly::compose(row, &[variable_shift.clone(), Rational::one()])
            };
            let term = Self::new(n_power.iter().map(|c| poly::scale(&shifted, c)).collect());
            result = result.add(&term);
        }
        result
    }

    /// The index translate `p(n + index_shift, y)`.
    pub fn shift_index(&self, index_shift: i64) -> Self {
        self.shift(
            &Rational::from_integer(index_shift.into()),
            &Rational::zero(),
        )
    }

    /// The partial derivative `d/dy`.
    pub fn variable_derivative(&self) -> Self {
        Self::new(self.0.iter().map(|row| poly::derivative(row)).collect())
    }

    /// The polynomial in `y` at `n = value` (canonical: zero is `[0]`).
    pub fn evaluate_index(&self, value: &Rational) -> Vec<Rational> {
        let mut result = vec![Rational::zero()];
        for row in self.0.iter().rev() {
            result = poly::scale(&result, value);
            result = poly::add(&result, row);
        }
        poly::trim(&result)
    }

    /// The polynomial in `n` at `y = value` (canonical: zero is `[0]`).
    pub fn evaluate_variable(&self, value: &Rational) -> Vec<Rational> {
        poly::trim(
            &self
                .0
                .iter()
                .map(|row| poly::evaluate(row, value))
                .collect::<Vec<_>>(),
        )
    }

    /// The coefficient of `y^j` as a polynomial in `n` (canonical).
    pub fn variable_coefficient(&self, j: usize) -> Vec<Rational> {
        poly::trim(
            &self
                .0
                .iter()
                .map(|row| row.get(j).cloned().unwrap_or_else(Rational::zero))
                .collect::<Vec<_>>(),
        )
    }

    /// The degree in `y`, `None` for zero.
    pub fn variable_degree(&self) -> Option<usize> {
        self.0.iter().filter_map(|row| poly::degree(row)).max()
    }

    /// The degree in `n`, `None` for zero.
    pub fn index_degree(&self) -> Option<usize> {
        self.0.len().checked_sub(1)
    }

    /// The smallest power of `y` with a nonzero coefficient, `None` for zero.
    pub fn variable_valuation(&self) -> Option<usize> {
        self.0
            .iter()
            .filter_map(|row| row.iter().position(|c| !c.is_zero()))
            .min()
    }

    /// Divide by `y^e`, dropping the coefficients of `y^j`, `j < e` (the
    /// caller guarantees divisibility).
    pub fn divide_by_variable_power(&self, e: usize) -> Self {
        Self::new(
            self.0
                .iter()
                .map(|row| row.get(e..).map(<[Rational]>::to_vec).unwrap_or_default())
                .collect(),
        )
    }

    /// Multiply by `y^e`.
    pub fn multiply_by_variable_power(&self, e: usize) -> Self {
        Self::new(
            self.0
                .iter()
                .map(|row| {
                    if row.is_empty() {
                        Vec::new()
                    } else {
                        poly::shift_by(row, e)
                    }
                })
                .collect(),
        )
    }

    /// The lexicographically leading coefficient: highest power of `n`,
    /// then of `y`; `None` for zero.
    pub fn leading(&self) -> Option<&Rational> {
        self.0.last().and_then(|row| row.last())
    }

    /// `(lead, p / lead)`; the zero polynomial is returned unchanged with
    /// lead one.
    pub fn monic(&self) -> (Rational, Self) {
        match self.leading() {
            Some(lead) => {
                let lead = lead.clone();
                (lead.clone(), self.scale(&lead.recip()))
            }
            None => (Rational::one(), Self::zero()),
        }
    }
}
