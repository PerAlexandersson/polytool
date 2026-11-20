//! Exact homogeneous Möbius substitution in O(n²) Gaussian operations.

use super::{trim_gaussian, StabilityError};
use crate::GaussianRational as G;

/// A nonsingular map `M(w) = (a*w + b)/(c*w + d)`.
///
/// The fields are private so every instance has nonzero determinant.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MobiusTransform {
    a: G,
    b: G,
    c: G,
    d: G,
}

/// Result of the homogeneous substitution `(c*w+d)^n p(M(w))`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TransformedPolynomial {
    /// Ascending coefficients, trimmed; empty for the zero polynomial.
    pub coefficients: Vec<G>,
    /// Degree before substitution, or `None` for the zero polynomial.
    pub original_degree: Option<usize>,
    /// Multiplicity of the original root at `M(∞) = a/c`.
    ///
    /// This is the degree drop, counted on the region boundary. It is zero
    /// when `c = 0`, including for constant and zero polynomials.
    pub boundary_at_infinity: usize,
}

impl MobiusTransform {
    /// Construct a map, rejecting `a*d - b*c = 0` exactly.
    pub fn new(a: G, b: G, c: G, d: G) -> Result<Self, StabilityError> {
        if a.multiply(&d).subtract(&b.multiply(&c)).is_zero() {
            Err(StabilityError::SingularMobiusTransform)
        } else {
            Ok(Self { a, b, c, d })
        }
    }

    /// Coefficients `[a, b, c, d]` in the defining fractional linear map.
    pub fn coefficients(&self) -> [&G; 4] {
        [&self.a, &self.b, &self.c, &self.d]
    }

    /// The identity map.
    pub fn identity() -> Self {
        Self {
            a: G::one(),
            b: G::zero(),
            c: G::zero(),
            d: G::one(),
        }
    }

    /// The inverse map, with matrix `[d, -b; -c, a]`.
    pub fn inverse(&self) -> Self {
        Self {
            a: self.d.clone(),
            b: G::zero().subtract(&self.b),
            c: G::zero().subtract(&self.c),
            d: self.a.clone(),
        }
    }

    /// Compose maps: the result sends `w` to `self(other(w))`.
    pub fn compose(&self, other: &Self) -> Self {
        Self {
            a: self.a.multiply(&other.a).add(&self.b.multiply(&other.c)),
            b: self.a.multiply(&other.b).add(&self.b.multiply(&other.d)),
            c: self.c.multiply(&other.a).add(&self.d.multiply(&other.c)),
            d: self.c.multiply(&other.b).add(&self.d.multiply(&other.d)),
        }
    }

    /// Evaluate at a finite point; `None` means the image is infinity.
    pub fn evaluate(&self, point: &G) -> Option<G> {
        self.a
            .multiply(point)
            .add(&self.b)
            .checked_divide(&self.c.multiply(point).add(&self.d))
    }

    /// Pull a polynomial back to the upper half-plane, tracking degree loss.
    ///
    /// Homogeneous Horner evaluation costs O(n²) field operations and O(n)
    /// working storage; no division by `c*w+d` occurs.
    pub fn transform_polynomial(&self, coefficients: &[G]) -> TransformedPolynomial {
        self.transform_polynomial_with_checkpoint(coefficients, &|| false)
            .expect("unlimited substitution")
    }

    pub(super) fn transform_polynomial_with_checkpoint(
        &self,
        coefficients: &[G],
        expired: &impl Fn() -> bool,
    ) -> Result<TransformedPolynomial, StabilityError> {
        super::checkpoint(expired)?;
        let mut p = coefficients.to_vec();
        trim_gaussian(&mut p);
        let Some(n) = p.len().checked_sub(1) else {
            return Ok(TransformedPolynomial {
                coefficients: p,
                original_degree: None,
                boundary_at_infinity: 0,
            });
        };
        let mut q = vec![p[n].clone()];
        let mut denominator_power = vec![G::one()];
        for coefficient in p[..n].iter().rev() {
            super::checkpoint(expired)?;
            q = multiply_linear(&q, &self.a, &self.b);
            denominator_power = multiply_linear(&denominator_power, &self.c, &self.d);
            q.resize(denominator_power.len(), G::zero());
            for (entry, power) in q.iter_mut().zip(&denominator_power) {
                *entry = entry.add(&coefficient.multiply(power));
            }
        }
        trim_gaussian(&mut q);
        // Nonsingular changes of homogeneous coordinates preserve nonzeroness.
        let degree = q.len() - 1;
        Ok(TransformedPolynomial {
            coefficients: q,
            original_degree: Some(n),
            boundary_at_infinity: n - degree,
        })
    }
}

fn multiply_linear(p: &[G], a: &G, b: &G) -> Vec<G> {
    let mut q = vec![G::zero(); p.len() + 1];
    for (k, coefficient) in p.iter().enumerate() {
        q[k] = q[k].add(&coefficient.multiply(b));
        q[k + 1] = q[k + 1].add(&coefficient.multiply(a));
    }
    q
}
