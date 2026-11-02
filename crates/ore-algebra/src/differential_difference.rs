//! Normal-ordered operators in the Ore algebra of the index shift, the
//! derivative and the shift of the variable.
//!
//! The variable shift `E` acts on polynomials by `(E f)(x) = f(x + 1)`, and
//! the forward difference is `Delta = E - 1`.  Together with the index shift
//! `S`, the derivative `D = d/dx` and multiplication by polynomials in `n`
//! and `x`, it generates the Ore algebra with the rules
//!
//! ```text
//! S n = (n + 1) S,   D x = x D + 1,   E x = (x + 1) E,
//! ```
//!
//! and `S`, `D`, `E` pairwise commuting (`E n = n E`, `D n = n D`,
//! `S x = x S`).  Since `E` is invertible on polynomials, integer powers
//! `E^e`, `e < 0` included, are allowed (`E^(-1) x = (x - 1) E^(-1)`).
//! Every element has a unique normal-ordered form
//! `sum c_(i,d,e)(n, x) D^d E^e S^i` ([`DifferentialDifferenceOperator`]).
//!
//! # Representation
//!
//! This is a separate type, not a generalisation of [`OreOperator`]: the
//! existing modules (pseudo-division, row maps, row annihilators, image
//! annihilators) rely on `OreOperator`'s two-index keys and on the absence
//! of `E` (the reversal `x -> 1/x` and the coefficient multipliers do not
//! commute with `E` in any useful way), so a third key in that type would
//! change every caller.  The embedding
//! [`DifferentialDifferenceOperator::from_ore_operator`] and the partial
//! inverse [`DifferentialDifferenceOperator::to_ore_operator`] (defined
//! exactly when no `E` occurs) connect the two, and the composition of
//! embedded operators is the composition in `OreOperator`.
//!
//! # Composition
//!
//! Moving `D^a E^e S^i` past a coefficient `c(n, x)`:
//! `S^i c = c(n + i, x) S^i`, `E^e c = c(n, x + e) E^e` and
//! `D^a c = sum_t binom(a, t) (d/dx)^t c D^(a - t)`, so
//!
//! ```text
//! (c D^a E^e S^i) o (c' D^b E^f S^j)
//!     = sum_(t = 0..a) binom(a, t) c (d/dx)^t c'(n + i, x + e) D^(a - t + b) E^(e + f) S^(i + j).
//! ```
//!
//! # Action on rows
//!
//! `(L F)_n(x) = sum c_(i,d,e)(n, x) F_(n+i)^((d))(x + e)`; the three
//! operators commute, so the order of `D^d`, `E^e` and `S^i` in a term does
//! not matter.  The action is a left module action: `(A o B) F = A (B F)`.

use crate::bivariate::BivariatePolynomial;
use crate::operator::{binomial, rational, OreOperator};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// The key `(shift i, derivative order d, variable shift e)` of a term
/// `c D^d E^e S^i`.
pub type DifferentialDifferenceKey = (usize, usize, i64);

/// A normal-ordered operator `sum c_(i,d,e)(n, x) D^d E^e S^i` with
/// polynomial coefficients (module documentation), keyed by
/// `(shift i, derivative order d, variable shift e)`.
///
/// The field is the term map itself.  Its invariant, kept by every method,
/// is that no zero coefficient is stored; code that writes the map directly
/// must keep it (use [`DifferentialDifferenceOperator::add_term`]).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DifferentialDifferenceOperator(
    pub BTreeMap<DifferentialDifferenceKey, BivariatePolynomial>,
);

/// `p(x + shift)` for an ascending coefficient vector.
pub(crate) fn translate_polynomial(p: &[Rational], shift: &Rational) -> Vec<Rational> {
    if shift.is_zero() {
        return poly::trim(p);
    }
    poly::compose(p, &[shift.clone(), Rational::one()])
}

impl DifferentialDifferenceOperator {
    /// Add `coefficient D^d E^e S^i` for `key = (i, d, e)`.
    pub fn add_term(&mut self, key: DifferentialDifferenceKey, coefficient: BivariatePolynomial) {
        if coefficient.is_zero() {
            return;
        }
        let entry = self.0.entry(key).or_default();
        *entry = entry.add(&coefficient);
        if entry.is_zero() {
            self.0.remove(&key);
        }
    }

    /// The single term `coefficient D^derivative E^variable_shift S^shift`.
    pub fn single(
        shift: usize,
        derivative: usize,
        variable_shift: i64,
        coefficient: BivariatePolynomial,
    ) -> Self {
        let mut op = Self::default();
        op.add_term((shift, derivative, variable_shift), coefficient);
        op
    }

    /// The identity operator.
    pub fn identity() -> Self {
        Self::single(0, 0, 0, BivariatePolynomial::one())
    }

    /// The index shift `S^shift`.
    pub fn index_shift(shift: usize) -> Self {
        Self::single(shift, 0, 0, BivariatePolynomial::one())
    }

    /// The derivative `D^order`.
    pub fn derivative(order: usize) -> Self {
        Self::single(0, order, 0, BivariatePolynomial::one())
    }

    /// The variable shift `E^power`, `(E^power f)(x) = f(x + power)`.
    pub fn variable_shift(power: i64) -> Self {
        Self::single(0, 0, power, BivariatePolynomial::one())
    }

    /// Multiplication by `coefficient(n, x)`.
    pub fn multiplication(coefficient: BivariatePolynomial) -> Self {
        Self::single(0, 0, 0, coefficient)
    }

    /// The forward difference `Delta^order = (E - 1)^order`.
    pub fn forward_difference(order: usize) -> Self {
        let mut out = Self::default();
        for t in 0..=order {
            let sign = if (order - t) % 2 == 0 {
                Rational::one()
            } else {
                -Rational::one()
            };
            out.add_term(
                (0, 0, t as i64),
                BivariatePolynomial::constant(binomial(order, t) * sign),
            );
        }
        out
    }

    /// The backward difference `nabla^order = (1 - E^(-1))^order`,
    /// `(nabla f)(x) = f(x) - f(x - 1)`.
    pub fn backward_difference(order: usize) -> Self {
        let mut out = Self::default();
        for t in 0..=order {
            let sign = if t % 2 == 0 {
                Rational::one()
            } else {
                -Rational::one()
            };
            out.add_term(
                (0, 0, -(t as i64)),
                BivariatePolynomial::constant(binomial(order, t) * sign),
            );
        }
        out
    }

    /// The embedding of an operator without `E`.
    pub fn from_ore_operator(op: &OreOperator) -> Self {
        let mut out = Self::default();
        for ((i, d), c) in &op.0 {
            out.add_term((*i, *d, 0), c.clone());
        }
        out
    }

    /// The operator as an [`OreOperator`] when no power of `E` other than
    /// `E^0` occurs; `None` otherwise.
    pub fn to_ore_operator(&self) -> Option<OreOperator> {
        let mut out = OreOperator::default();
        for ((i, d, e), c) in &self.0 {
            if *e != 0 {
                return None;
            }
            out.add_term((*i, *d), c.clone());
        }
        Some(out)
    }

    /// Whether this is the zero operator.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The number of stored terms.
    pub fn term_count(&self) -> usize {
        self.0.len()
    }

    /// Whether no power of `E` other than `E^0` occurs.
    pub fn is_variable_shift_free(&self) -> bool {
        self.0.keys().all(|(_, _, e)| *e == 0)
    }

    /// Whether no derivative occurs.
    pub fn is_derivative_free(&self) -> bool {
        self.0.keys().all(|(_, d, _)| *d == 0)
    }

    /// The sum `self + other`.
    pub fn add(&self, other: &Self) -> Self {
        let mut sum = self.clone();
        for (key, c) in &other.0 {
            sum.add_term(*key, c.clone());
        }
        sum
    }

    /// The scalar multiple `c self`.
    pub fn scale(&self, c: &Rational) -> Self {
        let mut out = Self::default();
        for (key, coefficient) in &self.0 {
            out.add_term(*key, coefficient.scale(c));
        }
        out
    }

    /// The difference `self - other`.
    pub fn subtract(&self, other: &Self) -> Self {
        self.add(&other.scale(&-Rational::one()))
    }

    /// `factor o self` for a coefficient `factor(n, x)` (multiplication on
    /// the left, which is coefficientwise).
    pub fn left_multiply(&self, factor: &BivariatePolynomial) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, factor.multiply(c));
        }
        out
    }

    /// The largest shift (zero for the zero operator).
    pub fn max_shift(&self) -> usize {
        self.0.keys().map(|(i, _, _)| *i).max().unwrap_or(0)
    }

    /// The largest derivative order (zero for the zero operator).
    pub fn max_derivative(&self) -> usize {
        self.0.keys().map(|(_, d, _)| *d).max().unwrap_or(0)
    }

    /// The smallest power of `E` (zero for the zero operator).
    pub fn min_variable_shift(&self) -> i64 {
        self.0.keys().map(|(_, _, e)| *e).min().unwrap_or(0)
    }

    /// The largest power of `E` (zero for the zero operator).
    pub fn max_variable_shift(&self) -> i64 {
        self.0.keys().map(|(_, _, e)| *e).max().unwrap_or(0)
    }

    /// The largest degree in `x` of a coefficient (zero for the zero
    /// operator).
    pub fn max_variable_degree(&self) -> usize {
        self.0
            .values()
            .filter_map(BivariatePolynomial::variable_degree)
            .max()
            .unwrap_or(0)
    }

    /// The terms at shift `shift`, as a shift-free operator.
    pub fn column(&self, shift: usize) -> Self {
        let mut out = Self::default();
        for ((i, d, e), c) in &self.0 {
            if *i == shift {
                out.add_term((0, *d, *e), c.clone());
            }
        }
        out
    }

    /// The terms with `E^variable_shift`, as an operator without `E`.
    pub fn variable_shift_component(&self, variable_shift: i64) -> OreOperator {
        let mut out = OreOperator::default();
        for ((i, d, e), c) in &self.0 {
            if *e == variable_shift {
                out.add_term((*i, *d), c.clone());
            }
        }
        out
    }

    /// The composition `self o other` (module documentation).
    pub fn compose(&self, other: &Self) -> Self {
        let mut out = Self::default();
        for ((i, a, e), c) in &self.0 {
            let dn = rational(*i);
            let de = Rational::from_integer((*e).into());
            for ((j, b, f), coefficient) in &other.0 {
                let mut derivative = coefficient.shift(&dn, &de);
                for t in 0..=*a {
                    if derivative.is_zero() {
                        break;
                    }
                    let term = c.multiply(&derivative).scale(&binomial(*a, t));
                    out.add_term((i + j, a - t + b, e + f), term);
                    derivative = derivative.variable_derivative();
                }
            }
        }
        out
    }

    /// `self^exponent` (the identity for exponent zero).
    pub fn power(&self, exponent: usize) -> Self {
        (0..exponent).fold(Self::identity(), |acc, _| acc.compose(self))
    }

    /// `S^shift o self`.
    pub fn shifted(&self, shift: usize) -> Self {
        let mut out = Self::default();
        let dn = rational(shift);
        for ((i, d, e), c) in &self.0 {
            out.add_term((i + shift, *d, *e), c.shift(&dn, &Rational::zero()));
        }
        out
    }

    /// `E^power o self`: every coefficient is translated, `x -> x + power`.
    pub fn variable_shifted(&self, power: i64) -> Self {
        let mut out = Self::default();
        let dx = Rational::from_integer(power.into());
        for ((i, d, e), c) in &self.0 {
            out.add_term((*i, *d, e + power), c.shift(&Rational::zero(), &dx));
        }
        out
    }

    /// The conjugate `E^power o self o E^(-power)`: `x -> x + power` in every
    /// coefficient.
    pub fn conjugate_by_variable_shift(&self, power: i64) -> Self {
        let mut out = Self::default();
        let dx = Rational::from_integer(power.into());
        for (key, c) in &self.0 {
            out.add_term(*key, c.shift(&Rational::zero(), &dx));
        }
        out
    }

    /// The operator with `n` replaced by `n + index_shift` in every
    /// coefficient (the relation read at base `n + index_shift`).
    pub fn shift_index(&self, index_shift: i64) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.shift_index(index_shift));
        }
        out
    }

    /// The leading coefficient: the largest key, then the lexicographically
    /// leading coefficient of its polynomial.
    pub fn leading(&self) -> Option<Rational> {
        self.0
            .values()
            .next_back()
            .and_then(|c| c.leading().cloned())
    }

    /// One term `c(n, x) F^((d))(x + e)` at the index value `value`.
    fn term_value(
        c: &BivariatePolynomial,
        d: usize,
        e: i64,
        row: &[Rational],
        value: &Rational,
    ) -> Vec<Rational> {
        let mut derivative = poly::trim(row);
        for _ in 0..d {
            derivative = poly::derivative(&derivative);
        }
        let shifted = translate_polynomial(&derivative, &Rational::from_integer(e.into()));
        poly::mul(&c.evaluate_index(value), &shifted)
    }

    /// `(self F)_n` when `rows[j]` is `F_(origin + j)` and the rows cover
    /// every referenced index; `None` otherwise.  Canonical (zero is `[0]`).
    pub fn apply(&self, rows: &[Vec<Rational>], origin: usize, n: usize) -> Option<Vec<Rational>> {
        let value = rational(n);
        let mut sum = vec![Rational::zero()];
        for ((i, d, e), c) in &self.0 {
            let index = n.checked_add(*i)?.checked_sub(origin)?;
            let row = rows.get(index)?;
            sum = poly::add(&sum, &Self::term_value(c, *d, *e, row, &value));
        }
        Some(poly::trim(&sum))
    }

    /// Apply a shift-free operator to one row at index `n` (shifts are
    /// ignored).  Canonical (zero is `[0]`).
    pub fn apply_row(&self, row: &[Rational], n: usize) -> Vec<Rational> {
        let value = rational(n);
        let mut sum = vec![Rational::zero()];
        for ((_, d, e), c) in &self.0 {
            sum = poly::add(&sum, &Self::term_value(c, *d, *e, row, &value));
        }
        poly::trim(&sum)
    }
}
