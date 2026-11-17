//! Monomial orders, leading terms and exponent-vector divisibility.

use std::cmp::Ordering;
use std::fmt;
use std::str::FromStr;

use super::{Coefficient, SparsePolynomial};

/// A monomial order, with the convention `x_1 > x_2 > ... > x_n`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MonomialOrder {
    /// Lexicographic order with `x_1 > x_2 > ... > x_n`.
    Lex,
    /// Total degree, then lexicographic order.
    GrLex,
    /// Total degree, then graded reverse lexicographic order.
    GrevLex,
}

/// The largest term of a polynomial with respect to a [`MonomialOrder`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LeadingTerm<C> {
    /// Exponent vector of the leading monomial.
    pub exponents: Vec<u32>,
    /// Its (nonzero) coefficient.
    pub coefficient: C,
}

/// Error from parsing an unknown monomial-order name (`MonomialOrder::from_str`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParseMonomialOrderError {
    input: String,
}

impl ParseMonomialOrderError {
    /// The rejected input string.
    pub fn input(&self) -> &str {
        &self.input
    }
}

impl fmt::Display for ParseMonomialOrderError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "unknown monomial order `{}`", self.input)
    }
}

impl std::error::Error for ParseMonomialOrderError {}

impl MonomialOrder {
    /// All supported orders, in the order `Lex`, `GrLex`, `GrevLex`.
    pub const STANDARD_ORDERS: [MonomialOrder; 3] = [
        MonomialOrder::Lex,
        MonomialOrder::GrLex,
        MonomialOrder::GrevLex,
    ];

    /// Compare two monomials; `Greater` means `a` is the larger monomial.
    ///
    /// # Panics
    ///
    /// Panics if the exponent vectors have different lengths, or (for `GrLex` and `GrevLex`) if a
    /// total degree overflows `u32`.
    #[inline]
    pub fn compare(&self, a: &[u32], b: &[u32]) -> Ordering {
        assert_eq!(a.len(), b.len(), "monomials have different lengths");
        match self {
            MonomialOrder::Lex => compare_lex(a, b),
            MonomialOrder::GrLex => compare_total_degree(a, b).then_with(|| compare_lex(a, b)),
            MonomialOrder::GrevLex => {
                compare_total_degree(a, b).then_with(|| compare_grevlex_same_degree(a, b))
            }
        }
    }
}

impl fmt::Display for MonomialOrder {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            MonomialOrder::Lex => "lex",
            MonomialOrder::GrLex => "grlex",
            MonomialOrder::GrevLex => "grevlex",
        };
        f.write_str(name)
    }
}

impl FromStr for MonomialOrder {
    type Err = ParseMonomialOrderError;

    fn from_str(input: &str) -> Result<Self, Self::Err> {
        match input.trim().to_ascii_lowercase().as_str() {
            "lex" | "lexicographic" => Ok(MonomialOrder::Lex),
            "grlex" | "graded-lex" | "graded_lex" | "gradedlex" => Ok(MonomialOrder::GrLex),
            "grevlex" | "graded-revlex" | "graded_revlex" | "gradedrevlex"
            | "graded-reverse-lex" | "graded_reverse_lex" | "gradedreverselex" => {
                Ok(MonomialOrder::GrevLex)
            }
            _ => Err(ParseMonomialOrderError {
                input: input.to_string(),
            }),
        }
    }
}

/// Leading term of `polynomial` under `order`, or `None` if it is zero.
///
/// # Panics
///
/// Panics if a graded comparison overflows the total-degree `u32` range.
pub fn leading_term<P: SparsePolynomial>(
    polynomial: &P,
    order: MonomialOrder,
) -> Option<LeadingTerm<Coefficient<P>>> {
    let term = match order {
        MonomialOrder::Lex => polynomial.terms().iter().next_back(),
        MonomialOrder::GrLex | MonomialOrder::GrevLex => polynomial
            .terms()
            .iter()
            .max_by(|(a, _), (b, _)| order.compare(a, b)),
    };
    term.map(|(exponents, coefficient)| LeadingTerm {
        exponents: exponents.clone(),
        coefficient: coefficient.clone(),
    })
}

/// Whether `x^divisor` divides `x^dividend` (componentwise `<=`).
///
/// # Panics
///
/// Panics if the exponent vectors have different lengths.
#[inline]
pub fn monomial_divides(divisor: &[u32], dividend: &[u32]) -> bool {
    assert_eq!(
        divisor.len(),
        dividend.len(),
        "monomials have different lengths"
    );
    divisor.iter().zip(dividend.iter()).all(|(&a, &b)| a <= b)
}

/// Exponent vector of `x^dividend / x^divisor`, or `None` if it does not divide.
///
/// # Panics
///
/// Panics if the exponent vectors have different lengths.
#[inline]
pub fn monomial_quotient(dividend: &[u32], divisor: &[u32]) -> Option<Vec<u32>> {
    if !monomial_divides(divisor, dividend) {
        return None;
    }
    Some(
        dividend
            .iter()
            .zip(divisor.iter())
            .map(|(&a, &b)| a - b)
            .collect(),
    )
}

#[inline]
fn compare_total_degree(a: &[u32], b: &[u32]) -> Ordering {
    checked_total_degree(a).cmp(&checked_total_degree(b))
}

#[inline]
fn compare_lex(a: &[u32], b: &[u32]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()) {
        match x.cmp(&y) {
            Ordering::Equal => continue,
            ordering => return ordering,
        }
    }
    Ordering::Equal
}

#[inline]
fn compare_grevlex_same_degree(a: &[u32], b: &[u32]) -> Ordering {
    for (&x, &y) in a.iter().zip(b.iter()).rev() {
        match x.cmp(&y) {
            Ordering::Equal => continue,
            Ordering::Less => return Ordering::Greater,
            Ordering::Greater => return Ordering::Less,
        }
    }
    Ordering::Equal
}

#[inline]
fn checked_total_degree(exponents: &[u32]) -> u32 {
    exponents
        .iter()
        .try_fold(0u32, |sum, &e| sum.checked_add(e))
        .expect("total degree overflow")
}
