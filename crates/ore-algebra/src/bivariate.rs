//! Bivariate polynomials in the index `n` and a variable `y`.

use crate::univariate as poly;
use crate::Rational;
use num_bigint::BigInt;
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
        let mut sum = self.clone();
        sum.add_assign(other);
        sum
    }

    /// `self += other`, in place (keeps the trimmed representation).
    pub(crate) fn add_assign(&mut self, other: &Self) {
        self.combine_assign(other, poly::add_rational_into);
    }

    /// `self -= other`, in place (keeps the trimmed representation).
    pub(crate) fn subtract_assign(&mut self, other: &Self) {
        self.combine_assign(other, poly::subtract_rational_from);
    }

    fn combine_assign(&mut self, other: &Self, combine: fn(&mut Rational, &Rational)) {
        if self.0.len() < other.0.len() {
            self.0.resize(other.0.len(), Vec::new());
        }
        for (row, addend) in self.0.iter_mut().zip(&other.0) {
            if addend.is_empty() {
                continue;
            }
            if row.len() < addend.len() {
                row.resize(addend.len(), Rational::zero());
            }
            for (x, y) in row.iter_mut().zip(addend) {
                combine(x, y);
            }
            poly::strip_trailing_zeros(row);
        }
        while self.0.last().is_some_and(Vec::is_empty) {
            self.0.pop();
        }
    }

    /// The difference `self - other`.
    pub fn subtract(&self, other: &Self) -> Self {
        let mut difference = self.clone();
        difference.subtract_assign(other);
        difference
    }

    /// The scalar multiple `c self`.
    pub fn scale(&self, c: &Rational) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        if c.is_one() {
            return self.clone();
        }
        Self::new(
            self.0
                .iter()
                .map(|row| row.iter().map(|x| poly::rational_product(x, c)).collect())
                .collect(),
        )
    }

    /// The product `self other`.
    pub fn multiply(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        if let Some(c) = other.as_constant() {
            return self.scale(&c);
        }
        if let Some(c) = self.as_constant() {
            return other.scale(&c);
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
                let target: &mut Vec<Rational> = &mut rows[i + j];
                let width = left.len() + right.len() - 1;
                if target.len() < width {
                    target.resize(width, Rational::zero());
                }
                for (a, x) in left.iter().enumerate() {
                    if x.is_zero() {
                        continue;
                    }
                    for (b, y) in right.iter().enumerate() {
                        if !y.is_zero() {
                            poly::add_product_into(&mut target[a + b], x, y);
                        }
                    }
                }
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
        if index_shift.is_zero() && variable_shift.is_zero() {
            return self.clone();
        }
        // p(n + h, y + v) = sum_k n^k sum_(i >= k) binom(i, k) h^(i - k) p_i(y + v).
        let len = self.0.len();
        let mut powers = Vec::with_capacity(len);
        if len > 0 {
            powers.push(Rational::one());
        }
        for i in 1..len {
            powers.push(poly::rational_product(&powers[i - 1], index_shift));
        }
        let mut binomials: Vec<BigInt> = Vec::with_capacity(len);
        let mut rows: Vec<Vec<Rational>> = vec![Vec::new(); len];
        for (i, row) in self.0.iter().enumerate() {
            // binomials[k] = binom(i, k).
            binomials.push(BigInt::one());
            for k in (1..i).rev() {
                let previous = binomials[k - 1].clone();
                binomials[k] += previous;
            }
            if row.is_empty() {
                continue;
            }
            let shifted = if variable_shift.is_zero() {
                row.clone()
            } else {
                poly::translate(row, variable_shift)
            };
            for (k, binomial) in binomials.iter().enumerate() {
                let power = &powers[i - k];
                if power.is_zero() {
                    continue;
                }
                let factor =
                    poly::rational_product(&Rational::from_integer(binomial.clone()), power);
                let target = &mut rows[k];
                if target.len() < shifted.len() {
                    target.resize(shifted.len(), Rational::zero());
                }
                for (value, c) in target.iter_mut().zip(&shifted) {
                    if !c.is_zero() {
                        poly::add_product_into(value, c, &factor);
                    }
                }
            }
        }
        Self::new(rows)
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
        let width = self.0.iter().map(Vec::len).max().unwrap_or(0);
        let mut result = vec![Rational::zero(); width.max(1)];
        for row in self.0.iter().rev() {
            for x in result.iter_mut() {
                if !x.is_zero() {
                    *x = poly::rational_product(x, value);
                }
            }
            for (x, c) in result.iter_mut().zip(row) {
                poly::add_rational_into(x, c);
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn q(numerator: i64, denominator: i64) -> Rational {
        Rational::new(numerator.into(), denominator.into())
    }

    fn sample(seed: u64, height: usize, width: usize) -> BivariatePolynomial {
        let mut state = seed;
        let rows = (0..height)
            .map(|_| {
                (0..width)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        let value = ((state >> 33) % 9) as i64 - 4;
                        // Some entries non-integral, some zero.
                        if state % 5 == 0 {
                            q(value, 3)
                        } else {
                            q(value, 1)
                        }
                    })
                    .collect()
            })
            .collect();
        BivariatePolynomial::new(rows)
    }

    /// The translate by the expanded binomial formula equals the sum of the
    /// products `(n + h)^i p_i(y + v)`.
    #[test]
    fn shift_agrees_with_the_naive_expansion() {
        for seed in 0..12 {
            let p = sample(seed, 5, 4);
            for (h, v) in [
                (q(0, 1), q(0, 1)),
                (q(1, 1), q(0, 1)),
                (q(-3, 1), q(0, 1)),
                (q(2, 7), q(0, 1)),
                (q(0, 1), q(1, 1)),
                (q(-1, 2), q(5, 3)),
            ] {
                let linear = BivariatePolynomial::new(vec![vec![h.clone()], vec![q(1, 1)]]);
                let mut expected = BivariatePolynomial::zero();
                for (i, row) in p.rows().iter().enumerate() {
                    let shifted = poly::compose(row, &[v.clone(), q(1, 1)]);
                    let term = linear
                        .power(i)
                        .multiply(&BivariatePolynomial::from_variable_polynomial(&shifted));
                    expected = expected.add(&term);
                }
                assert_eq!(p.shift(&h, &v), expected, "seed {seed}, h {h}, v {v}");
            }
        }
    }

    /// Products, sums and evaluations against coefficientwise references.
    #[test]
    fn arithmetic_agrees_with_coefficientwise_references() {
        for seed in 0..12 {
            let a = sample(seed, 4, 3);
            let b = sample(seed + 100, 3, 4);
            let mut product = vec![vec![q(0, 1); 6]; 6];
            for (i, left) in a.rows().iter().enumerate() {
                for (j, right) in b.rows().iter().enumerate() {
                    for (k, x) in left.iter().enumerate() {
                        for (l, y) in right.iter().enumerate() {
                            product[i + j][k + l] += x * y;
                        }
                    }
                }
            }
            assert_eq!(a.multiply(&b), BivariatePolynomial::new(product));
            let mut sum = vec![vec![q(0, 1); 4]; 4];
            for (i, row) in a.rows().iter().enumerate() {
                for (k, x) in row.iter().enumerate() {
                    sum[i][k] += x;
                }
            }
            for (i, row) in b.rows().iter().enumerate() {
                for (k, x) in row.iter().enumerate() {
                    sum[i][k] -= x;
                }
            }
            assert_eq!(a.subtract(&b), BivariatePolynomial::new(sum));
            assert!(a.subtract(&a).is_zero() && a.subtract(&a).rows().is_empty());
            for value in [q(0, 1), q(3, 1), q(-2, 5)] {
                let mut expected = vec![q(0, 1)];
                for row in a.rows().iter().rev() {
                    expected = poly::add(&poly::scale(&expected, &value), row);
                }
                assert_eq!(a.evaluate_index(&value), poly::trim(&expected));
            }
        }
    }
}
