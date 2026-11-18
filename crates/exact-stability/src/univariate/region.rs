//! Circular regions as exact images of the open upper half-plane.

use super::{mobius::MobiusTransform, StabilityConvention, StabilityError};
use crate::{BigRational as Q, GaussianRational as G};
use num_traits::Signed;

/// An open region of the complex plane; its boundary is counted separately.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum StabilityRegion {
    /// `Im(z) > 0`; the default convention is avoidance.
    Upper,
    /// `Im(z) < 0`; the default convention is avoidance.
    Lower,
    /// `Re(z) < 0`; the default convention is all roots inside (Hurwitz).
    Left,
    /// `Re(z) > 0`; the default convention is all roots inside.
    Right,
    /// `|z| < 1`; the default convention is all roots inside (Schur).
    Disk,
    /// `|z| > 1`; the default convention is all roots inside.
    Exterior,
    /// `|z-center| < radius`, with a strictly positive rational radius.
    CenteredDisk {
        /// Gaussian-rational center.
        center: G,
        /// Strictly positive radius, checked before computation.
        radius: Q,
    },
    /// Image of the open upper half-plane under a nonsingular map.
    /// The default convention is avoidance.
    Mobius(Box<MobiusTransform>),
}

impl From<MobiusTransform> for StabilityRegion {
    fn from(map: MobiusTransform) -> Self {
        Self::Mobius(Box::new(map))
    }
}

impl StabilityRegion {
    /// The classical convention associated with this region.
    pub fn default_convention(&self) -> StabilityConvention {
        match self {
            Self::Upper | Self::Lower | Self::Mobius(_) => StabilityConvention::Avoid,
            _ => StabilityConvention::Inside,
        }
    }

    /// A map from the upper half-plane onto this open region.
    ///
    /// A nonpositive centered-disk radius returns a typed error.
    pub fn mobius_transform(&self) -> Result<MobiusTransform, StabilityError> {
        let zero = G::zero();
        let one = G::one();
        let i = G::imaginary_unit();
        let minus_i = zero.subtract(&i);
        match self {
            Self::Upper => Ok(MobiusTransform::identity()),
            Self::Lower => MobiusTransform::new(G::from(-1), zero.clone(), zero, one),
            Self::Left => MobiusTransform::new(i, zero.clone(), zero, one),
            Self::Right => MobiusTransform::new(minus_i, zero.clone(), zero, one),
            Self::Disk => MobiusTransform::new(one.clone(), minus_i, one, i),
            Self::Exterior => MobiusTransform::new(one.clone(), i, one, minus_i),
            Self::CenteredDisk { center, radius } => {
                if !radius.is_positive() {
                    return Err(StabilityError::NonpositiveRadius);
                }
                let r = G::from(radius.clone());
                MobiusTransform::new(
                    center.add(&r),
                    i.multiply(&center.subtract(&r)),
                    G::one(),
                    i,
                )
            }
            Self::Mobius(map) => Ok(map.as_ref().clone()),
        }
    }

    /// Exact membership sign: positive inside, zero on the boundary, negative outside.
    ///
    /// Finite poles of the inverse map are boundary points, including `M(∞)`.
    pub fn membership_sign(&self, point: &G) -> Result<i8, StabilityError> {
        let inverse = self.mobius_transform()?.inverse();
        let Some(preimage) = inverse.evaluate(point) else {
            return Ok(0);
        };
        Ok(if preimage.imaginary.is_positive() {
            1
        } else if preimage.imaginary.is_negative() {
            -1
        } else {
            0
        })
    }
}
