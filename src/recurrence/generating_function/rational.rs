//! Construction, exact cancellation and Riordan recognition.

use super::super::{parse_rational_coeff, Recurrence, RecurrenceError};
use crate::{vec_poly, BigRational as Q, Polynomial};
use num_traits::{One, Zero};
use ore_algebra::BivariatePolynomial as B;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Stable schema of [`RationalGeneratingFunction`].
pub const RATIONAL_GENERATING_FUNCTION_SCHEMA: &str = "polytool.generating-function.v1";

/// Reduced exact `N(x,z)/D(x,z)` with `D(x,0) = 1`.
///
/// Coefficient arrays are ragged and ascending: `[i][j]` is the coefficient
/// of `z^i x^j`. Zero is `[]`; coefficients are canonical rational strings.
/// Public fields support serialization; recognition validates them before use.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RationalGeneratingFunction {
    /// Always [`RATIONAL_GENERATING_FUNCTION_SCHEMA`].
    pub schema: String,
    /// `z^0` represents this index, rather than the absolute index zero.
    pub first_index: usize,
    /// The reduced numerator, indexed `[z_power][x_power]`.
    pub numerator: Vec<Vec<String>>,
    /// The reduced denominator, indexed `[z_power][x_power]`.
    pub denominator: Vec<Vec<String>>,
}

/// A reduced univariate rational function regular at zero.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RationalGeneratingSeries {
    /// Numerator, ascending in `z`; zero is `[]`.
    pub numerator: Vec<String>,
    /// Denominator, ascending in `z`, with constant coefficient one.
    pub denominator: Vec<String>,
}

impl RationalGeneratingSeries {
    /// Cancel the gcd and normalize the denominator at zero.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] for a denominator vanishing at zero.
    pub fn new(numerator: &[Q], denominator: &[Q]) -> Result<Self, RecurrenceError> {
        let n = Polynomial::new(numerator.to_vec());
        let d = Polynomial::new(denominator.to_vec());
        if d.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "zero series denominator".into(),
            ));
        }
        let gcd = n.gcd(&d);
        let n = n.exact_div(&gcd);
        let d = d.exact_div(&gcd);
        let constant = d
            .coeffs()
            .first()
            .filter(|c| !c.is_zero())
            .ok_or_else(|| RecurrenceError::InvalidInput("series has a pole at zero".into()))?;
        Ok(Self {
            numerator: encode(&vec_poly::scale(n.coeffs(), &constant.recip())),
            denominator: encode(&vec_poly::scale(d.coeffs(), &constant.recip())),
        })
    }

    /// Parse the ascending exact numerator and denominator.
    ///
    /// # Errors
    ///
    /// Invalid rational strings, zero denominator or a pole at zero.
    pub fn coefficients(&self) -> Result<(Vec<Q>, Vec<Q>), RecurrenceError> {
        let n = decode(&self.numerator)?;
        let d = decode(&self.denominator)?;
        Self::new(&n, &d)?;
        Ok((n, d))
    }
}

/// The identity `F = d/(1 - x h)`, with `h = z d K`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RiordanDecomposition {
    /// The first-column generating function `N/D_0`.
    pub d: RationalGeneratingSeries,
    /// The column multiplier `L/D_0`, with zero constant coefficient.
    pub h: RationalGeneratingSeries,
    /// The left-production factor `L/(z N)`.
    pub k: RationalGeneratingSeries,
}

fn encode(p: &[Q]) -> Vec<String> {
    Polynomial::new(p.to_vec())
        .coeffs()
        .iter()
        .map(ToString::to_string)
        .collect()
}

fn decode(p: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    p.iter()
        .map(|s| parse_rational_coeff(s).map_err(Into::into))
        .collect()
}

pub(in crate::recurrence) fn encode_bivariate(p: &B) -> Vec<Vec<String>> {
    p.rows().iter().map(|r| encode(r)).collect()
}

pub(in crate::recurrence) fn decode_bivariate(p: &[Vec<String>]) -> Result<B, RecurrenceError> {
    Ok(B::new(
        p.iter().map(|r| decode(r)).collect::<Result<_, _>>()?,
    ))
}

// Primitive polynomial PRS over Q[z][x]. Gauss' lemma reduces a bivariate
// gcd to content gcds in Q[z] and pseudo-remainders in x. Each outer step
// strictly decreases the x degree; each inner step cancels its leading term.
fn content(p: &B) -> Polynomial<Q> {
    (0..=p.variable_degree().unwrap_or(0)).fold(Polynomial::zero(), |g, j| {
        g.gcd(&Polynomial::new(p.variable_coefficient(j)))
    })
}

fn divide_content(p: &B, c: &Polynomial<Q>) -> B {
    let mut out = B::zero();
    for j in 0..=p.variable_degree().unwrap_or(0) {
        let q = Polynomial::new(p.variable_coefficient(j)).exact_div(c);
        out = out.add(&B::from_index_polynomial(q.coeffs()).multiply_by_variable_power(j));
    }
    out
}

fn primitive(p: &B) -> B {
    if p.is_zero() {
        p.clone()
    } else {
        divide_content(p, &content(p)).monic().1
    }
}

// Exact coprimality certificate over Q[x][z]. A specialization preserving
// both z-degrees preserves the z-degree of every common divisor. Thus a
// coprime specialized pair excludes every divisor of positive z-degree.
// The gcd of all coefficient rows separately excludes divisors in Q[x].
// Degree preservation is essential: at x=0, 1+xz would become a unit.
fn certified_coprime(a: &B, b: &B) -> bool {
    let coefficient_content = a
        .rows()
        .iter()
        .chain(b.rows())
        .fold(Polynomial::zero(), |g, row| {
            g.gcd(&Polynomial::new(row.clone()))
        });
    if coefficient_content.degree() != Some(0) {
        return false;
    }
    [0, 1, -1, 2].into_iter().any(|point| {
        let point = Q::from_integer(point.into());
        let specialized_a = Polynomial::new(a.evaluate_variable(&point));
        let specialized_b = Polynomial::new(b.evaluate_variable(&point));
        specialized_a.degree() == a.index_degree()
            && specialized_b.degree() == b.index_degree()
            && specialized_a.gcd(&specialized_b).degree() == Some(0)
    })
}

fn gcd_bivariate(a: &B, b: &B) -> B {
    if a.is_zero() {
        return b.monic().1;
    }
    if b.is_zero() {
        return a.monic().1;
    }
    if certified_coprime(a, b) {
        return B::one();
    }
    let c = content(a).gcd(&content(b));
    let mut a = primitive(a);
    let mut b = primitive(b);
    while !b.is_zero() {
        let degree = b.variable_degree().unwrap();
        let lead = B::from_index_polynomial(&b.variable_coefficient(degree));
        while !a.is_zero() && a.variable_degree().unwrap() >= degree {
            let shift = a.variable_degree().unwrap() - degree;
            let coefficient = B::from_index_polynomial(&a.variable_coefficient(degree + shift));
            a = a
                .multiply(&lead)
                .subtract(&b.multiply(&coefficient).multiply_by_variable_power(shift));
            a = primitive(&a);
        }
        (a, b) = (b, a);
    }
    a.multiply(&B::from_index_polynomial(c.coeffs())).monic().1
}

fn exact_divide(a: &B, b: &B) -> Result<B, RecurrenceError> {
    let mut remainder = a.clone();
    let mut result = B::zero();
    let degree = b
        .variable_degree()
        .ok_or_else(|| RecurrenceError::InvalidInput("zero divisor".into()))?;
    let lead = Polynomial::new(b.variable_coefficient(degree));
    while !remainder.is_zero() && remainder.variable_degree().unwrap() >= degree {
        let shift = remainder.variable_degree().unwrap() - degree;
        let (coefficient, rest) =
            Polynomial::new(remainder.variable_coefficient(degree + shift)).div_rem(&lead);
        if !rest.is_zero() {
            break;
        }
        let term = B::from_index_polynomial(coefficient.coeffs()).multiply_by_variable_power(shift);
        remainder = remainder.subtract(&term.multiply(b));
        result = result.add(&term);
    }
    if !remainder.is_zero() {
        return Err(RecurrenceError::VerificationFailed(
            "bivariate gcd did not divide exactly".into(),
        ));
    }
    Ok(result)
}

fn reduced_pair(n: &B, d: &B) -> Result<(B, B), RecurrenceError> {
    if d.is_zero() {
        return Err(RecurrenceError::InvalidInput(
            "zero generating-function denominator".into(),
        ));
    }
    let gcd = gcd_bivariate(n, d);
    let n = exact_divide(n, &gcd)?;
    let d = exact_divide(d, &gcd)?;
    let constant = B::from_variable_polynomial(&d.evaluate_index(&Q::zero()))
        .as_constant()
        .filter(|c| !c.is_zero())
        .ok_or_else(|| {
            RecurrenceError::Unsupported(
                "reduced generating-function denominator at z=0 must be a nonzero scalar".into(),
            )
        })?;
    Ok((n.scale(&constant.recip()), d.scale(&constant.recip())))
}

/// Build the exact rational generating function from a homogeneous rule
/// free of the index and of derivatives, and **every** prescribed seed row.
///
/// The rule starts after the last seed; `z^0` is the row at `first_index`.
/// With `D = f(x) - sum a_j(x) z^j`, the numerator is the part of
/// `D sum_(m < seeds.len()) P_(first_index+m)(x) z^m` below `z^seeds.len()`.
/// Its degree can exceed the recurrence order. The bivariate gcd is cancelled
/// exactly, then the scalar is normalized so `D(x,0)=1`.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] for index dependence, derivatives,
/// alternating signs, inhomogeneity, or a reduced denominator whose value at
/// `z=0` is not a scalar (the requested polynomial normalization is impossible).
/// [`RecurrenceError::InvalidInput`] for zero left factor, offset zero, too
/// few seeds or overflowing index range.
pub fn rational_generating_function(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<RationalGeneratingFunction, RecurrenceError> {
    let unsupported = || {
        RecurrenceError::Unsupported(
            "rational generating functions require a homogeneous index-free rule \
             without derivatives or alternating signs"
                .into(),
        )
    };
    if rule.inhomogeneous.is_some() {
        return Err(unsupported());
    }
    first_index
        .checked_add(seeds.len())
        .ok_or_else(|| RecurrenceError::InvalidInput("index overflow".into()))?;
    let coefficient = |p: &super::super::BivarPoly| -> Result<Vec<Q>, RecurrenceError> {
        let p = B::new(p.coeffs.clone());
        if !p.is_index_free() {
            return Err(unsupported());
        }
        Ok(p.evaluate_index(&Q::zero()))
    };
    let lhs = rule
        .denominator
        .as_ref()
        .map(coefficient)
        .transpose()?
        .unwrap_or_else(|| vec![Q::one()]);
    if vec_poly::is_zero(&lhs) {
        return Err(RecurrenceError::InvalidInput("zero left factor".into()));
    }
    let mut d = B::from_variable_polynomial(&lhs);
    let mut order = 0;
    for term in &rule.terms {
        if term.deriv_order != 0 || term.is_alternating() {
            return Err(unsupported());
        }
        if term.offset == 0 {
            return Err(RecurrenceError::InvalidInput(
                "recurrence offset must be positive".into(),
            ));
        }
        let c = coefficient(&term.coeff)?;
        if !vec_poly::is_zero(&c) {
            order = order.max(term.offset);
            let mut rows = vec![Vec::new(); term.offset];
            rows.push(c);
            d = d.subtract(&B::new(rows));
        }
    }
    if seeds.len() < order {
        return Err(RecurrenceError::InvalidInput(
            "fewer seeds than the recurrence order".into(),
        ));
    }
    let product = d.multiply(&B::new(seeds.to_vec()));
    let n = B::new(product.rows().iter().take(seeds.len()).cloned().collect());
    let (n, d) = reduced_pair(&n, &d)?;
    Ok(RationalGeneratingFunction {
        schema: RATIONAL_GENERATING_FUNCTION_SCHEMA.into(),
        first_index,
        numerator: encode_bivariate(&n),
        denominator: encode_bivariate(&d),
    })
}

/// Recognize `F=d/(1-xh)` with `h=z d K`, returning reduced exact series.
///
/// Returns `None` for malformed data, a numerator involving `x`, a denominator
/// of degree greater than one in `x`, zero numerator, or a pole at zero in
/// any returned series. The input is reduced again before recognition.
pub fn riordan_decomposition(gf: &RationalGeneratingFunction) -> Option<RiordanDecomposition> {
    if gf.schema != RATIONAL_GENERATING_FUNCTION_SCHEMA {
        return None;
    }
    let (n, d) = reduced_pair(
        &decode_bivariate(&gf.numerator).ok()?,
        &decode_bivariate(&gf.denominator).ok()?,
    )
    .ok()?;
    if n.is_zero() || !n.is_variable_free() || d.variable_degree()? > 1 {
        return None;
    }
    let n = n.variable_coefficient(0);
    let d0 = d.variable_coefficient(0);
    let l = vec_poly::neg(&d.variable_coefficient(1));
    let zn = vec_poly::shift_by(&n, 1);
    Some(RiordanDecomposition {
        d: RationalGeneratingSeries::new(&n, &d0).ok()?,
        h: RationalGeneratingSeries::new(&l, &d0).ok()?,
        k: RationalGeneratingSeries::new(&l, &zn).ok()?,
    })
}

fn polynomial_text(p: &B) -> String {
    let mut terms = Vec::new();
    for (i, row) in p.rows().iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let mut term = format!("({c})");
            if j > 0 {
                term.push_str(&format!("*x^{j}"));
            }
            if i > 0 {
                term.push_str(&format!("*z^{i}"));
            }
            terms.push(term);
        }
    }
    if terms.is_empty() {
        "0".into()
    } else {
        terms.join(" + ")
    }
}

impl fmt::Display for RationalGeneratingFunction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match (
            decode_bivariate(&self.numerator),
            decode_bivariate(&self.denominator),
        ) {
            (Ok(n), Ok(d)) => write!(f, "({})/({})", polynomial_text(&n), polynomial_text(&d)),
            _ => f.write_str("invalid rational generating function"),
        }
    }
}

impl fmt::Display for RationalGeneratingSeries {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.coefficients() {
            Ok((n, d)) => write!(
                f,
                "({})/({})",
                polynomial_text(&B::from_index_polynomial(&n)),
                polynomial_text(&B::from_index_polynomial(&d))
            ),
            Err(_) => f.write_str("invalid rational series"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn coprimality_specialization_checks_degrees_and_coefficient_content() {
        let q = |n: i64| Q::from_integer(n.into());
        let common = B::new(vec![vec![q(1)], vec![q(0), q(1)]]); // 1+xz
        assert!(!certified_coprime(&common, &common.multiply(&common)));
        let x = B::from_variable_polynomial(&[q(0), q(1)]);
        assert!(!certified_coprime(&x, &x.multiply(&common)));
        let a = B::new(vec![vec![q(1)], vec![q(0), q(0), q(1)]]);
        let b = B::new(vec![vec![q(1)], vec![q(0), q(0), q(2)]]);
        assert!(certified_coprime(&a, &b));
        assert_eq!(gcd_bivariate(&a, &b), B::one());
        assert_eq!(
            gcd_bivariate(&common, &common.multiply(&common)),
            common.monic().1
        );
    }
}
