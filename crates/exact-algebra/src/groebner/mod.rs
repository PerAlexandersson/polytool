//! Sparse commutative Gröbner bases over exact coefficient fields.
//!
//! [`SparsePolynomial`] borrows a consumer's term map, avoiding conversions in
//! normal-form loops. [`Polynomial`] supplies an owned representation and
//! [`RationalField`] supplies arbitrary-precision rational coefficients.
//! Variables are zero-indexed; exponent vectors have fixed length and terms
//! are sparse. [`reduced_groebner_basis`] returns descending leading terms.
//! Computation is unbounded and intended for small systems. Matrix reduction
//! uses shared symbolic reducer rows, followed by scalar verification/fallback;
//! it is not a full F4 implementation.

use num_rational::BigRational;
use num_traits::{One, Zero};
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::ops::{Add, Mul, Sub};

mod buchberger;
mod division;
mod monomial_order;

pub use buchberger::*;
pub use division::*;
pub use monomial_order::*;

/// Arithmetic provider for a coefficient ring, independent of consumer traits.
pub trait CoefficientRing: Clone + Debug + Eq {
    /// Exact coefficient type. Bounded types may panic on overflow.
    type Element: Clone
        + Debug
        + Eq
        + Add<Output = Self::Element>
        + Sub<Output = Self::Element>
        + Mul<Output = Self::Element>;
    /// Additive identity.
    fn zero() -> Self::Element;
    /// Multiplicative identity.
    fn one() -> Self::Element;
    /// Whether a coefficient is zero.
    fn is_zero(value: &Self::Element) -> bool;
}

/// Exact division by nonzero coefficients, satisfying the field laws.
///
/// Implementors must not use integer truncation or approximate division.
pub trait ExactField: CoefficientRing {
    /// Divide two coefficients exactly.
    ///
    /// # Panics
    ///
    /// May panic for zero denominator or overflow of bounded coefficient types.
    fn divide(numerator: Self::Element, denominator: Self::Element) -> Self::Element;
}

/// Coefficient type associated with a sparse polynomial representation.
pub type Coefficient<P> = <<P as SparsePolynomial>::Field as CoefficientRing>::Element;

/// Representation contract for sparse polynomials without arithmetic adapters.
///
/// Maps contain only nonzero coefficients and fixed-length exponent vectors;
/// map iteration is ascending lexicographic order. Implementations of `Add` and
/// `Sub` must preserve those invariants and reject mismatched variable counts.
pub trait SparsePolynomial: Clone + Debug + Eq + Add<Output = Self> + Sub<Output = Self> {
    /// Coefficient arithmetic provider; field algorithms additionally require [`ExactField`].
    type Field: CoefficientRing;
    /// Number of variables, including for zero polynomials.
    fn num_vars(&self) -> usize;
    /// Borrow the nonzero terms, indexed by exponent vectors.
    fn terms(&self) -> &BTreeMap<Vec<u32>, Coefficient<Self>>;
    /// Construct zero in a specified ambient ring.
    fn zero(num_vars: usize) -> Self;
    /// Construct from sparse terms, removing zeros.
    ///
    /// # Panics
    ///
    /// Panics if an exponent vector has length different from `num_vars`.
    fn from_terms(num_vars: usize, terms: BTreeMap<Vec<u32>, Coefficient<Self>>) -> Self;
    /// Construct a single term.
    ///
    /// # Panics
    ///
    /// Panics if the exponent vector has length different from `num_vars`.
    fn monomial(num_vars: usize, exponents: Vec<u32>, coefficient: Coefficient<Self>) -> Self {
        Self::from_terms(num_vars, BTreeMap::from([(exponents, coefficient)]))
    }
    /// Whether there are no nonzero terms.
    fn is_zero(&self) -> bool {
        self.terms().is_empty()
    }
    /// Multiply every coefficient by a scalar.
    fn scale(&self, coefficient: &Coefficient<Self>) -> Self {
        Self::from_terms(
            self.num_vars(),
            self.terms()
                .iter()
                .map(|(m, c)| (m.clone(), c.clone() * coefficient.clone()))
                .collect(),
        )
    }
}

/// Arbitrary-precision rational coefficient arithmetic.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RationalField;
impl CoefficientRing for RationalField {
    type Element = BigRational;
    fn zero() -> BigRational {
        BigRational::zero()
    }
    fn one() -> BigRational {
        BigRational::one()
    }
    fn is_zero(value: &BigRational) -> bool {
        value.is_zero()
    }
}
impl ExactField for RationalField {
    fn divide(numerator: BigRational, denominator: BigRational) -> BigRational {
        numerator / denominator
    }
}

/// An owned sparse polynomial over a coefficient provider.
///
/// Use the [`SparsePolynomial`] constructors and accessors. For example:
/// ```
/// use exact_algebra::groebner::{Polynomial, RationalField, SparsePolynomial,
///     MonomialOrder, reduced_groebner_basis};
/// let x = Polynomial::<RationalField>::monomial(1, vec![1], num_rational::BigRational::from_integer(1.into()));
/// assert_eq!(reduced_groebner_basis(&[x.clone(), x.clone()], MonomialOrder::Lex), vec![x]);
/// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial<F: CoefficientRing> {
    num_vars: usize,
    terms: BTreeMap<Vec<u32>, F::Element>,
}
impl<F: CoefficientRing> SparsePolynomial for Polynomial<F> {
    type Field = F;
    fn num_vars(&self) -> usize {
        self.num_vars
    }
    fn terms(&self) -> &BTreeMap<Vec<u32>, F::Element> {
        &self.terms
    }
    fn zero(num_vars: usize) -> Self {
        Self {
            num_vars,
            terms: BTreeMap::new(),
        }
    }
    fn from_terms(num_vars: usize, mut terms: BTreeMap<Vec<u32>, F::Element>) -> Self {
        assert!(
            terms.keys().all(|e| e.len() == num_vars),
            "all exponent vectors must have length num_vars"
        );
        terms.retain(|_, c| !F::is_zero(c));
        Self { num_vars, terms }
    }
}
impl<F: CoefficientRing> Add for Polynomial<F> {
    type Output = Self;
    fn add(mut self, rhs: Self) -> Self {
        assert_eq!(
            self.num_vars, rhs.num_vars,
            "polynomials have different rings"
        );
        for (m, c) in rhs.terms {
            let coefficient = self.terms.get(&m).cloned().unwrap_or_else(F::zero) + c;
            if F::is_zero(&coefficient) {
                self.terms.remove(&m);
            } else {
                self.terms.insert(m, coefficient);
            }
        }
        self
    }
}
impl<F: CoefficientRing> Sub for Polynomial<F> {
    type Output = Self;
    fn sub(mut self, rhs: Self) -> Self {
        assert_eq!(
            self.num_vars, rhs.num_vars,
            "polynomials have different rings"
        );
        for (m, c) in rhs.terms {
            let coefficient = self.terms.get(&m).cloned().unwrap_or_else(F::zero) - c;
            if F::is_zero(&coefficient) {
                self.terms.remove(&m);
            } else {
                self.terms.insert(m, coefficient);
            }
        }
        self
    }
}

#[cfg(test)]
mod tests;
