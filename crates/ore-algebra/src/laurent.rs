//! Localization at `x`: operators with coefficients in `Q[n, x, x^(-1)]`.
//!
//! The Ore algebra of [`OreOperator`]s extends to Laurent coefficients,
//! because the derivative of a Laurent polynomial is one
//! (`D x^(-1) = x^(-1) D - x^(-2)`) and the index shift does not touch `x`.
//! (The variable shift `E` does not extend: `E x^(-1) = (x + 1)^(-1) E`
//! leaves the Laurent polynomials, so [`crate::DifferentialDifferenceOperator`]s
//! keep polynomial coefficients.)  Every element has a unique normal-ordered
//! form `sum c_(i,d)(n, x) D^d S^i` with Laurent coefficients
//! ([`LaurentOreOperator`]); it acts on sequences of Laurent polynomials,
//! and multiplication by a nonzero Laurent polynomial is injective there.
//!
//! # Clearing
//!
//! `x^N L` is an [`OreOperator`] once `N` is at least minus the smallest
//! power of `x` in a coefficient ([`LaurentOreOperator::clear_denominators`],
//! [`LaurentOreOperator::clear_with_power`]).  Since `x^N` is injective,
//! `(L F)_n = 0` exactly when `(x^N L F)_n = 0`, so a relation with Laurent
//! coefficients and its cleared form annihilate the same rows, and the
//! cleared form is a polynomial certificate.
//!
//! # Native conjugations
//!
//! For a map `C` of sequences that is invertible on Laurent polynomials,
//! `C^(-1) L C` is computed on the generators (it is an algebra
//! automorphism), and `(L C F)_n = (C (C^(-1) L C) F)_n` exactly:
//!
//! - the **reversal** `(C F)_n = x^(d_n) F_n(1/x)`, `d_n = d_0 + h n` with
//!   an integer `h` of either sign: `c(n, x) -> c(n, 1/x)`,
//!   `D -> rho(n) = d_n x - x^2 D` (from
//!   `D (x^d K(1/x)) = x^d (d x^(-1) K - x^(-2) K')(1/x)`) and
//!   `S -> x^(-h) S` (from `x^(d_(n+1)) = x^(d_n) x^h`); `C` is an
//!   involution, so push-forward and pullback agree;
//! - the **power gauge** `(C F)_n = x^(w_n) F_n`, `w_n = w_0 + h n` with an
//!   integer `h` (dividing rows by powers of `x` is `w_0 < 0`):
//!   `c -> c`, `D -> D + w_n x^(-1)` and `S -> x^h S`; the push-forward is
//!   the pullback through `x^(-w_n)`;
//! - the **dilation** `(C F)_n(x) = F_n(a x)`, `a != 0`: `c(n, x) -> c(n, x/a)`,
//!   `D -> a D`, `S -> S`.
//!
//! The reversal conjugate of `M = sum m_(i,d) D^d S^i` is
//! `sum m_(i,d)(n, 1/x) rho(n)^d x^(-h i) S^i
//!  = sum x^(-h i) m_(i,d)(n, 1/x) rho(n + i)^d S^i`
//! (as `rho(n) x^(-h) = x^(-h) rho(n + 1)`), so `x^N` times it, with
//! `N = b + h R` for `b >= deg_x m` and `R >=` the largest shift, is the
//! operator that [`crate::pull_back`] computes for a
//! [`crate::RowMapStep::Reversal`] (`h >= 0`); here `h < 0` is allowed too.

use crate::bivariate::BivariatePolynomial;
use crate::operator::{binomial, rational, OreOperator};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, ToPrimitive, Zero};
use std::collections::BTreeMap;
use std::fmt;

/// A Laurent polynomial `x^valuation p(n, x)` in `x` with polynomial
/// coefficients in `n`, kept with `p` not divisible by `x` (and valuation
/// zero for the zero polynomial), so equality is structural.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct LaurentBivariatePolynomial {
    valuation: i64,
    polynomial: BivariatePolynomial,
}

/// Why a Laurent conjugation cannot be formed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LaurentError {
    /// A degree or power `d_0 + h n` is not of degree at most one in `n`
    /// with an integer slope `h` (and, for a power gauge, an integer
    /// `d_0`).
    InvalidExponent,
    /// A dilation factor is zero.
    ZeroDilation,
}

impl fmt::Display for LaurentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::InvalidExponent => {
                "an exponent must be d_0 + h n with an integer h (and an integer d_0 for a gauge)"
            }
            Self::ZeroDilation => "a dilation factor must be nonzero",
        })
    }
}

impl std::error::Error for LaurentError {}

impl LaurentBivariatePolynomial {
    /// `x^valuation p`, normalised.
    pub fn new(valuation: i64, polynomial: BivariatePolynomial) -> Self {
        match polynomial.variable_valuation() {
            None => Self::default(),
            Some(e) => Self {
                valuation: valuation + e as i64,
                polynomial: polynomial.divide_by_variable_power(e),
            },
        }
    }

    /// The polynomial `p` itself.
    pub fn from_polynomial(polynomial: &BivariatePolynomial) -> Self {
        Self::new(0, polynomial.clone())
    }

    /// The monomial `c n^index_power x^variable_power`.
    pub fn monomial(index_power: usize, variable_power: i64, c: Rational) -> Self {
        Self::new(
            variable_power,
            BivariatePolynomial::monomial(index_power, 0, c),
        )
    }

    /// The smallest power of `x` (zero for the zero polynomial).
    pub fn valuation(&self) -> i64 {
        self.valuation
    }

    /// `p` with `self = x^valuation p`, not divisible by `x`.
    pub fn stripped_polynomial(&self) -> &BivariatePolynomial {
        &self.polynomial
    }

    /// Whether this is zero.
    pub fn is_zero(&self) -> bool {
        self.polynomial.is_zero()
    }

    /// The polynomial `x^power self` when it has no negative powers of `x`.
    pub fn to_polynomial_times_power(&self, power: usize) -> Option<BivariatePolynomial> {
        if self.is_zero() {
            return Some(BivariatePolynomial::zero());
        }
        let total = self.valuation + power as i64;
        let total = usize::try_from(total).ok()?;
        Some(self.polynomial.multiply_by_variable_power(total))
    }

    /// The sum.
    pub fn add(&self, other: &Self) -> Self {
        if self.is_zero() {
            return other.clone();
        }
        if other.is_zero() {
            return self.clone();
        }
        let low = self.valuation.min(other.valuation);
        let a = self
            .polynomial
            .multiply_by_variable_power((self.valuation - low) as usize);
        let b = other
            .polynomial
            .multiply_by_variable_power((other.valuation - low) as usize);
        Self::new(low, a.add(&b))
    }

    /// The scalar multiple.
    pub fn scale(&self, c: &Rational) -> Self {
        Self::new(self.valuation, self.polynomial.scale(c))
    }

    /// The product.
    pub fn multiply(&self, other: &Self) -> Self {
        Self::new(
            self.valuation + other.valuation,
            self.polynomial.multiply(&other.polynomial),
        )
    }

    /// `x^power self`.
    pub fn multiply_by_variable_power(&self, power: i64) -> Self {
        Self::new(self.valuation + power, self.polynomial.clone())
    }

    /// `p(n + index_shift, x)`.
    pub fn shift_index(&self, index_shift: &Rational) -> Self {
        Self::new(
            self.valuation,
            self.polynomial.shift(index_shift, &Rational::zero()),
        )
    }

    /// `d/dx`.
    pub fn variable_derivative(&self) -> Self {
        // d/dx (x^v p) = x^(v-1) (v p + x p').
        let v = Rational::from_integer(self.valuation.into());
        let inner = self.polynomial.scale(&v).add(
            &self
                .polynomial
                .variable_derivative()
                .multiply_by_variable_power(1),
        );
        Self::new(self.valuation - 1, inner)
    }

    /// `self(n, 1/x)`.
    pub fn reciprocal_substitution(&self) -> Self {
        let Some(m) = self.polynomial.variable_degree() else {
            return Self::default();
        };
        let rows = self
            .polynomial
            .rows()
            .iter()
            .map(|row| {
                let mut reversed = vec![Rational::zero(); m + 1];
                for (j, c) in row.iter().enumerate() {
                    reversed[m - j] = c.clone();
                }
                reversed
            })
            .collect();
        Self::new(-self.valuation - m as i64, BivariatePolynomial::new(rows))
    }

    /// `self(n, factor x)`, `factor != 0`.
    pub fn dilate(&self, factor: &Rational) -> Self {
        let rows = self
            .polynomial
            .rows()
            .iter()
            .map(|row| {
                let mut power = Rational::one();
                row.iter()
                    .map(|c| {
                        let value = c * &power;
                        power *= factor;
                        value
                    })
                    .collect()
            })
            .collect();
        let mut outer = Rational::one();
        let base = if self.valuation >= 0 {
            factor.clone()
        } else {
            factor.recip()
        };
        for _ in 0..self.valuation.unsigned_abs() {
            outer *= &base;
        }
        Self::new(self.valuation, BivariatePolynomial::new(rows).scale(&outer))
    }

    /// The Laurent polynomial in `x` at `n = value`: `(valuation, ascending
    /// coefficients)`, canonical (zero is `(0, [0])`).
    pub fn evaluate_index(&self, value: &Rational) -> (i64, Vec<Rational>) {
        normalize_laurent_row(self.valuation, self.polynomial.evaluate_index(value))
    }
}

/// `(v, p)` with `p(0) != 0` (or `(0, [0])`).
fn normalize_laurent_row(valuation: i64, row: Vec<Rational>) -> (i64, Vec<Rational>) {
    let row = poly::trim(&row);
    match row.iter().position(|c| !c.is_zero()) {
        None => (0, vec![Rational::zero()]),
        Some(e) => (valuation + e as i64, row[e..].to_vec()),
    }
}

/// The sum of two Laurent rows.
fn add_laurent_rows(a: &(i64, Vec<Rational>), b: &(i64, Vec<Rational>)) -> (i64, Vec<Rational>) {
    if poly::is_zero(&a.1) {
        return b.clone();
    }
    if poly::is_zero(&b.1) {
        return a.clone();
    }
    let low = a.0.min(b.0);
    let left = poly::shift_by(&a.1, (a.0 - low) as usize);
    let right = poly::shift_by(&b.1, (b.0 - low) as usize);
    normalize_laurent_row(low, poly::add(&left, &right))
}

/// `(d_0, h)` for an exponent `d_0 + h n` of degree at most one with an
/// integer `h`; `d_0` must be an integer too when `integer_constant`.
fn affine_exponent(exponent: &[Rational], integer_constant: bool) -> Option<(Rational, i64)> {
    let exponent = poly::trim(exponent);
    if exponent.len() > 2 {
        return None;
    }
    let d0 = exponent[0].clone();
    let h = exponent.get(1).cloned().unwrap_or_else(Rational::zero);
    if !h.is_integer() || (integer_constant && !d0.is_integer()) {
        return None;
    }
    Some((d0, h.to_integer().to_i64()?))
}

/// A normal-ordered operator `sum c_(i,d)(n, x) D^d S^i` with Laurent
/// coefficients (module documentation), keyed by `(shift i, derivative
/// order d)`.  No zero coefficient is stored.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LaurentOreOperator(pub BTreeMap<(usize, usize), LaurentBivariatePolynomial>);

impl LaurentOreOperator {
    /// Add `coefficient D^d S^i` for `key = (i, d)`.
    pub fn add_term(&mut self, key: (usize, usize), coefficient: LaurentBivariatePolynomial) {
        if coefficient.is_zero() {
            return;
        }
        let entry = self.0.entry(key).or_default();
        *entry = entry.add(&coefficient);
        if entry.is_zero() {
            self.0.remove(&key);
        }
    }

    /// The single term `coefficient D^derivative S^shift`.
    pub fn single(
        shift: usize,
        derivative: usize,
        coefficient: LaurentBivariatePolynomial,
    ) -> Self {
        let mut op = Self::default();
        op.add_term((shift, derivative), coefficient);
        op
    }

    /// The identity operator.
    pub fn identity() -> Self {
        Self::single(
            0,
            0,
            LaurentBivariatePolynomial::from_polynomial(&BivariatePolynomial::one()),
        )
    }

    /// The embedding of a polynomial operator.
    pub fn from_ore_operator(op: &OreOperator) -> Self {
        let mut out = Self::default();
        for (key, c) in &op.0 {
            out.add_term(*key, LaurentBivariatePolynomial::from_polynomial(c));
        }
        out
    }

    /// Whether this is the zero operator.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The sum.
    pub fn add(&self, other: &Self) -> Self {
        let mut sum = self.clone();
        for (key, c) in &other.0 {
            sum.add_term(*key, c.clone());
        }
        sum
    }

    /// The scalar multiple.
    pub fn scale(&self, c: &Rational) -> Self {
        let mut out = Self::default();
        for (key, coefficient) in &self.0 {
            out.add_term(*key, coefficient.scale(c));
        }
        out
    }

    /// The difference.
    pub fn subtract(&self, other: &Self) -> Self {
        self.add(&other.scale(&-Rational::one()))
    }

    /// `factor o self` (coefficientwise).
    pub fn left_multiply(&self, factor: &LaurentBivariatePolynomial) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, factor.multiply(c));
        }
        out
    }

    /// `x^power o self`.
    pub fn multiply_by_variable_power(&self, power: i64) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.multiply_by_variable_power(power));
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

    /// The smallest power of `x` in a coefficient (zero for the zero
    /// operator).
    pub fn variable_valuation(&self) -> i64 {
        self.0
            .values()
            .map(LaurentBivariatePolynomial::valuation)
            .min()
            .unwrap_or(0)
    }

    /// The composition `self o other` (the rules of [`OreOperator`]).
    pub fn compose(&self, other: &Self) -> Self {
        let mut out = Self::default();
        for ((i, a), c) in &self.0 {
            let dn = rational(*i);
            for ((j, b), e) in &other.0 {
                let mut derivative = e.shift_index(&dn);
                for t in 0..=*a {
                    if derivative.is_zero() {
                        break;
                    }
                    let term = c.multiply(&derivative).scale(&binomial(*a, t));
                    out.add_term((i + j, a - t + b), term);
                    derivative = derivative.variable_derivative();
                }
            }
        }
        out
    }

    /// The operator as an [`OreOperator`] when no negative power of `x`
    /// occurs.
    pub fn to_ore_operator(&self) -> Option<OreOperator> {
        self.clear_with_power(0)
    }

    /// `x^power self` as an [`OreOperator`], or `None` when a negative power
    /// of `x` remains.
    pub fn clear_with_power(&self, power: usize) -> Option<OreOperator> {
        let mut out = OreOperator::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.to_polynomial_times_power(power)?);
        }
        Some(out)
    }

    /// `(N, x^N self)` with the smallest `N >= 0` making every coefficient a
    /// polynomial (module documentation).
    pub fn clear_denominators(&self) -> (usize, OreOperator) {
        let power = (-self.variable_valuation()).max(0) as usize;
        let cleared = self
            .clear_with_power(power)
            .expect("the power clears every coefficient");
        (power, cleared)
    }

    /// The conjugate `C^(-1) self C` by the reversal
    /// `(C F)_n = x^(d_n) F_n(1/x)`, `degree = [d_0, h]` with an integer `h`
    /// of either sign (module documentation).
    ///
    /// # Errors
    ///
    /// [`LaurentError::InvalidExponent`] for a degree that is not `d_0 + h n`
    /// with an integer `h`.
    pub fn pull_back_through_reversal(&self, degree: &[Rational]) -> Result<Self, LaurentError> {
        let (_, h) = affine_exponent(degree, false).ok_or(LaurentError::InvalidExponent)?;
        let degree_coefficient = LaurentBivariatePolynomial::from_polynomial(
            &BivariatePolynomial::from_index_polynomial(degree),
        );
        // rho(n) = d_n x - x^2 D.
        let mut rho = Self::single(
            0,
            1,
            LaurentBivariatePolynomial::monomial(0, 2, -Rational::one()),
        );
        rho.add_term((0, 0), degree_coefficient.multiply_by_variable_power(1));
        let mut powers = vec![Self::identity()];
        let mut out = Self::default();
        for ((i, d), c) in &self.0 {
            while powers.len() <= *d {
                let next = powers.last().expect("nonempty").compose(&rho);
                powers.push(next);
            }
            let shift = Self::single(
                *i,
                0,
                LaurentBivariatePolynomial::monomial(0, -h * *i as i64, Rational::one()),
            );
            let term = powers[*d]
                .compose(&shift)
                .left_multiply(&c.reciprocal_substitution());
            out = out.add(&term);
        }
        Ok(out)
    }

    /// The conjugate `C^(-1) self C` by the power gauge
    /// `(C F)_n = x^(w_n) F_n`, `power = [w_0, h]` integers of either sign
    /// (module documentation).
    ///
    /// # Errors
    ///
    /// [`LaurentError::InvalidExponent`] when `power` is not `w_0 + h n` with
    /// integers `w_0`, `h`.
    pub fn pull_back_through_power_gauge(&self, power: &[Rational]) -> Result<Self, LaurentError> {
        let (_, h) = affine_exponent(power, true).ok_or(LaurentError::InvalidExponent)?;
        let weight = LaurentBivariatePolynomial::from_polynomial(
            &BivariatePolynomial::from_index_polynomial(power),
        );
        // D + w_n x^(-1).
        let mut derivative = Self::single(
            0,
            1,
            LaurentBivariatePolynomial::from_polynomial(&BivariatePolynomial::one()),
        );
        derivative.add_term((0, 0), weight.multiply_by_variable_power(-1));
        let mut powers = vec![Self::identity()];
        let mut out = Self::default();
        for ((i, d), c) in &self.0 {
            while powers.len() <= *d {
                let next = powers.last().expect("nonempty").compose(&derivative);
                powers.push(next);
            }
            let shift = Self::single(
                *i,
                0,
                LaurentBivariatePolynomial::monomial(0, h * *i as i64, Rational::one()),
            );
            out = out.add(&powers[*d].compose(&shift).left_multiply(c));
        }
        Ok(out)
    }

    /// The conjugate `C^(-1) self C` by the dilation
    /// `(C F)_n(x) = F_n(factor x)`: `x -> x/factor` in every coefficient and
    /// `D^d -> factor^d D^d`.
    ///
    /// # Errors
    ///
    /// [`LaurentError::ZeroDilation`] for a zero factor.
    pub fn pull_back_through_dilation(&self, factor: &Rational) -> Result<Self, LaurentError> {
        if factor.is_zero() {
            return Err(LaurentError::ZeroDilation);
        }
        let inverse = factor.recip();
        let mut out = Self::default();
        for ((i, d), c) in &self.0 {
            let mut power = Rational::one();
            for _ in 0..*d {
                power *= factor;
            }
            out.add_term((*i, *d), c.dilate(&inverse).scale(&power));
        }
        Ok(out)
    }

    /// `(self F)_n` as a Laurent polynomial `(valuation, ascending
    /// coefficients)` when `rows[j]` is the Laurent polynomial `F_(origin+j)`
    /// given the same way and the rows cover every referenced index; `None`
    /// otherwise.  Canonical (zero is `(0, [0])`).
    pub fn apply(
        &self,
        rows: &[(i64, Vec<Rational>)],
        origin: usize,
        n: usize,
    ) -> Option<(i64, Vec<Rational>)> {
        let value = rational(n);
        let mut sum = (0, vec![Rational::zero()]);
        for ((i, d), c) in &self.0 {
            let index = n.checked_add(*i)?.checked_sub(origin)?;
            let (valuation, row) = rows.get(index)?;
            let mut derivative = LaurentBivariatePolynomial::new(
                *valuation,
                BivariatePolynomial::from_variable_polynomial(row),
            );
            for _ in 0..*d {
                derivative = derivative.variable_derivative();
            }
            let product = c.multiply(&derivative).evaluate_index(&value);
            sum = add_laurent_rows(&sum, &product);
        }
        Some(sum)
    }
}
