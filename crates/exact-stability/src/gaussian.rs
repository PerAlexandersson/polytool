//! Gaussian-rational coefficients, with exact elementary arithmetic.

use exact_algebra::BigRational as Q;
use num_bigint::BigInt;
use num_traits::{One, Zero};

/// An exact number `real + imaginary * i` in ℚ(i).
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct GaussianRational {
    /// Real part.
    pub real: Q,
    /// Imaginary part.
    pub imaginary: Q,
}

impl GaussianRational {
    /// Construct from exact real and imaginary parts.
    pub fn new(real: Q, imaginary: Q) -> Self {
        Self { real, imaginary }
    }

    /// The additive identity.
    pub fn zero() -> Self {
        Self::default()
    }

    /// The multiplicative identity.
    pub fn one() -> Self {
        Self::from(Q::one())
    }

    /// The imaginary unit.
    pub fn imaginary_unit() -> Self {
        Self::new(Q::zero(), Q::one())
    }

    /// Whether both parts are zero.
    pub fn is_zero(&self) -> bool {
        self.real.is_zero() && self.imaginary.is_zero()
    }

    /// Complex conjugation.
    pub fn conjugate(&self) -> Self {
        Self::new(self.real.clone(), -&self.imaginary)
    }

    /// Squared complex modulus.
    pub fn norm_squared(&self) -> Q {
        &self.real * &self.real + &self.imaginary * &self.imaginary
    }

    /// Exact sum.
    pub fn add(&self, other: &Self) -> Self {
        Self::new(&self.real + &other.real, &self.imaginary + &other.imaginary)
    }

    /// Exact difference.
    pub fn subtract(&self, other: &Self) -> Self {
        Self::new(&self.real - &other.real, &self.imaginary - &other.imaginary)
    }

    /// Exact product.
    pub fn multiply(&self, other: &Self) -> Self {
        Self::new(
            &self.real * &other.real - &self.imaginary * &other.imaginary,
            &self.real * &other.imaginary + &self.imaginary * &other.real,
        )
    }

    /// Scale by a real rational number.
    pub fn scale(&self, value: &Q) -> Self {
        Self::new(&self.real * value, &self.imaginary * value)
    }

    /// Exact quotient, or `None` when the divisor is zero.
    pub fn checked_divide(&self, other: &Self) -> Option<Self> {
        let norm = other.norm_squared();
        if norm.is_zero() {
            None
        } else {
            Some(self.multiply(&other.conjugate()).scale(&norm.recip()))
        }
    }
}

impl From<Q> for GaussianRational {
    fn from(real: Q) -> Self {
        Self::new(real, Q::zero())
    }
}

impl From<BigInt> for GaussianRational {
    fn from(real: BigInt) -> Self {
        Self::from(Q::from_integer(real))
    }
}

impl From<i64> for GaussianRational {
    fn from(real: i64) -> Self {
        Self::from(BigInt::from(real))
    }
}
