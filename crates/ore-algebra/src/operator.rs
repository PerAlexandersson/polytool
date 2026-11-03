//! Normal-ordered operators in the Ore algebra of shifts and derivatives.

use crate::bivariate::BivariatePolynomial;
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Zero};
use std::collections::btree_map::Entry;
use std::collections::BTreeMap;

/// `binom(a, t)` as an exact rational.
pub(crate) fn binomial(a: usize, t: usize) -> Rational {
    if t > a {
        return Rational::zero();
    }
    let mut value = num_bigint::BigInt::one();
    for i in 0..t {
        value *= num_bigint::BigInt::from(a - i);
        value /= num_bigint::BigInt::from(i + 1);
    }
    Rational::from_integer(value)
}

pub(crate) fn rational(value: usize) -> Rational {
    Rational::from_integer(value.into())
}

/// A normal-ordered operator `sum c_(i,d)(n, x) D^d S^i` with polynomial
/// coefficients, keyed by `(shift i, derivative order d)`.
///
/// The operator acts on a sequence of polynomials `F_n(x)` by
/// `(L F)_n = sum c_(i,d)(n, x) D^d F_(n+i)`; composition follows the Ore
/// rules `S n = (n + 1) S` and `D x = x D + 1` (and `S D = D S`,
/// `S x = x S`, `D n = n D`).
///
/// The field is the term map itself.  Its invariant, kept by every method,
/// is that no zero coefficient is stored; code that writes the map directly
/// must keep it (use [`OreOperator::add_term`]).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OreOperator(pub BTreeMap<(usize, usize), BivariatePolynomial>);

impl OreOperator {
    /// Add `coefficient D^d S^i` for `key = (i, d)`.
    pub fn add_term(&mut self, key: (usize, usize), coefficient: BivariatePolynomial) {
        if coefficient.is_zero() {
            return;
        }
        match self.0.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient);
            }
            Entry::Occupied(mut entry) => {
                entry.get_mut().add_assign(&coefficient);
                if entry.get().is_zero() {
                    entry.remove();
                }
            }
        }
    }

    /// Add `coefficient D^d S^i` for `key = (i, d)`, borrowing it.
    pub(crate) fn add_term_ref(&mut self, key: (usize, usize), coefficient: &BivariatePolynomial) {
        if coefficient.is_zero() {
            return;
        }
        match self.0.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient.clone());
            }
            Entry::Occupied(mut entry) => {
                entry.get_mut().add_assign(coefficient);
                if entry.get().is_zero() {
                    entry.remove();
                }
            }
        }
    }

    /// Subtract `coefficient D^d S^i` for `key = (i, d)`.
    pub(crate) fn subtract_term(&mut self, key: (usize, usize), coefficient: &BivariatePolynomial) {
        if coefficient.is_zero() {
            return;
        }
        match self.0.entry(key) {
            Entry::Vacant(entry) => {
                entry.insert(coefficient.scale(&-Rational::one()));
            }
            Entry::Occupied(mut entry) => {
                entry.get_mut().subtract_assign(coefficient);
                if entry.get().is_zero() {
                    entry.remove();
                }
            }
        }
    }

    /// `self += other`, in place.
    pub(crate) fn add_assign(&mut self, other: &Self) {
        for (key, c) in &other.0 {
            self.add_term_ref(*key, c);
        }
    }

    /// `self -= other`, in place.
    pub(crate) fn subtract_assign(&mut self, other: &Self) {
        for (key, c) in &other.0 {
            self.subtract_term(*key, c);
        }
    }

    /// The single term `coefficient D^derivative S^shift`.
    pub fn single(shift: usize, derivative: usize, coefficient: BivariatePolynomial) -> Self {
        let mut op = Self::default();
        op.add_term((shift, derivative), coefficient);
        op
    }

    /// The identity operator.
    pub fn identity() -> Self {
        Self::single(0, 0, BivariatePolynomial::one())
    }

    /// Whether this is the zero operator.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The number of stored terms.
    pub fn term_count(&self) -> usize {
        self.0.len()
    }

    /// The sum `self + other`.
    pub fn add(&self, other: &Self) -> Self {
        let mut sum = self.clone();
        sum.add_assign(other);
        sum
    }

    /// The scalar multiple `c self`.
    pub fn scale(&self, c: &Rational) -> Self {
        if c.is_one() {
            return self.clone();
        }
        let mut out = Self::default();
        for (key, coefficient) in &self.0 {
            out.add_term(*key, coefficient.scale(c));
        }
        out
    }

    /// The difference `self - other`.
    pub fn subtract(&self, other: &Self) -> Self {
        let mut difference = self.clone();
        difference.subtract_assign(other);
        difference
    }

    /// `factor o self` for a coefficient `factor(n, x)` (multiplication on
    /// the left, which is coefficientwise).
    pub fn left_multiply(&self, factor: &BivariatePolynomial) -> Self {
        if factor.as_constant().is_some_and(|c| c.is_one()) {
            return self.clone();
        }
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, factor.multiply(c));
        }
        out
    }

    /// The largest shift (zero for the zero operator).
    pub fn max_shift(&self) -> usize {
        self.0.keys().map(|(i, _)| *i).max().unwrap_or(0)
    }

    /// The largest derivative order (zero for the zero operator).
    pub fn max_derivative(&self) -> usize {
        self.0.keys().map(|(_, d)| *d).max().unwrap_or(0)
    }

    /// The terms at shift `shift`, as a shift-free operator.
    pub fn column(&self, shift: usize) -> Self {
        let mut out = Self::default();
        for ((i, d), c) in &self.0 {
            if *i == shift {
                out.add_term((0, *d), c.clone());
            }
        }
        out
    }

    /// The composition `self o other`.
    pub fn compose(&self, other: &Self) -> Self {
        let mut out = Self::default();
        // For the current shift i of `self`: every coefficient e of `other`
        // translated to e(n + i, x), with its x-derivatives (computed once
        // per shift and extended on demand).
        let mut cached_shift = None;
        let mut translated: Vec<Vec<BivariatePolynomial>> = Vec::new();
        for ((i, a), c) in &self.0 {
            if cached_shift != Some(*i) {
                let dn = rational(*i);
                translated = other
                    .0
                    .values()
                    .map(|e| vec![e.shift(&dn, &Rational::zero())])
                    .collect();
                cached_shift = Some(*i);
            }
            for (((j, b), _), derivatives) in other.0.iter().zip(translated.iter_mut()) {
                for t in 0..=*a {
                    if derivatives.len() <= t {
                        let next = derivatives.last().expect("nonempty").variable_derivative();
                        derivatives.push(next);
                    }
                    let derivative = &derivatives[t];
                    if derivative.is_zero() {
                        break;
                    }
                    let coefficient = c.multiply(derivative).scale(&binomial(*a, t));
                    out.add_term((i + j, a - t + b), coefficient);
                }
            }
        }
        out
    }

    /// `S^shift o self`.
    pub fn shifted(&self, shift: usize) -> Self {
        let mut out = Self::default();
        let dn = rational(shift);
        for ((i, d), c) in &self.0 {
            out.add_term((i + shift, *d), c.shift(&dn, &Rational::zero()));
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

    /// The smallest power of `x` dividing every coefficient.
    pub fn variable_valuation(&self) -> usize {
        self.0
            .values()
            .filter_map(BivariatePolynomial::variable_valuation)
            .min()
            .unwrap_or(0)
    }

    /// Divide every coefficient by `x^e` (the caller guarantees
    /// divisibility).
    pub fn divide_by_variable_power(&self, e: usize) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.divide_by_variable_power(e));
        }
        out
    }

    /// Multiply every coefficient by `x^e` (on the left).
    pub fn multiply_by_variable_power(&self, e: usize) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.multiply_by_variable_power(e));
        }
        out
    }

    /// The leading coefficient: the largest `(shift, derivative)` key, then
    /// the lexicographically leading coefficient of its polynomial.
    pub fn leading(&self) -> Option<Rational> {
        self.0
            .values()
            .next_back()
            .and_then(|c| c.leading().cloned())
    }

    /// `(lead, e, N)` with `self = lead x^e N`, `N` having no common power
    /// of `x` and leading coefficient one.
    pub fn normalize(&self) -> (Rational, usize, Self) {
        let e = self.variable_valuation();
        let stripped = self.divide_by_variable_power(e);
        let lead = stripped.leading().unwrap_or_else(Rational::one);
        let normalized = stripped.scale(&lead.recip());
        (lead, e, normalized)
    }

    /// `(self F)_n` when `rows[j]` is `F_(origin + j)` and the rows cover
    /// every referenced index; `None` otherwise.  Canonical (zero is `[0]`).
    pub fn apply(&self, rows: &[Vec<Rational>], origin: usize, n: usize) -> Option<Vec<Rational>> {
        let value = rational(n);
        let mut sum = vec![Rational::zero()];
        for ((i, d), c) in &self.0 {
            let index = n.checked_add(*i)?.checked_sub(origin)?;
            let row = rows.get(index)?;
            let mut derivative = row.clone();
            for _ in 0..*d {
                derivative = poly::derivative(&derivative);
            }
            sum = poly::add(&sum, &poly::mul(&c.evaluate_index(&value), &derivative));
        }
        Some(poly::trim(&sum))
    }

    /// Apply a shift-free operator to one row at index `n` (shifts are
    /// ignored).  Canonical (zero is `[0]`).
    pub fn apply_row(&self, row: &[Rational], n: usize) -> Vec<Rational> {
        let value = rational(n);
        let mut sum = vec![Rational::zero()];
        for ((_, d), c) in &self.0 {
            let mut derivative = row.to_vec();
            for _ in 0..*d {
                derivative = poly::derivative(&derivative);
            }
            sum = poly::add(&sum, &poly::mul(&c.evaluate_index(&value), &derivative));
        }
        poly::trim(&sum)
    }
}
