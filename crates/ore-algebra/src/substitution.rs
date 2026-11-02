//! Substitutions of the variable as first-class operations: affine maps,
//! reversals, Möbius maps, the `f`-polynomial map and power gauges, each
//! with its exact row map, inverse, pullback and forward rule.
//!
//! A [`Substitution`] `C` sends a sequence `F_n(x)` to `G_n = C_n(F_n)`:
//!
//! - `Affine`: `G_n(x) = F_n(a x + b)`, `a != 0`;
//! - `Reversal`: `G_n(x) = x^(d_n) F_n(1/x)`;
//! - `Mobius`: `G_n(x) = (s (c x + e))^(d_n) F_n((a x + b)/(c x + e))`,
//!   `a e - b c != 0`, `s != 0` (the scale `s` makes the family closed under
//!   inverses);
//! - `FPolynomial`: `G_n(x) = (1 + x)^(d_n) F_n(x/(1 + x))`, the Möbius map
//!   `(a, b, c, e) = (1, 0, 1, 1)`;
//! - `PowerGauge`: `G_n = x^(w_n) F_n` (with `w_n < 0` this divides the
//!   rows by a power of `x`).
//!
//! Exponents are `d_n = d_0 + h n` and `w_n = w_0 + h n` with integers
//! `d_0`, `h` (of either sign), given as ascending coefficient vectors
//! `[d_0, h]`.  The row map [`Substitution::apply_row`] is exact and returns
//! `None` where the image is not a polynomial (a row of degree above `d_n`,
//! a negative `d_n`, or a power of `x` that does not divide the row); the
//! operator identities below hold on Laurent polynomials and rational
//! functions regardless.
//!
//! # Pullbacks
//!
//! For an operator `M` acting on `G`, [`pull_back_substitution`] returns an
//! [`OreOperator`] `P` acting on `F` with `P = q o C^(-1) M C` for a nonzero
//! left factor `q` (a power of a linear form in `x`, times a nonzero
//! constant), so
//!
//! ```text
//! (M G)_n = 0   exactly when   (P F)_n = 0,
//! ```
//!
//! at every `n` (no exceptional indices: every factor is injective on
//! rational functions).  The conjugations are the Laurent ones of
//! [`crate::LaurentOreOperator`]:
//!
//! - an affine map substitutes `x -> (y - b)/a` and multiplies `D^d` by
//!   `a^d` (the [`crate::RowMapStep::Affine`] pullback);
//! - a reversal is [`crate::LaurentOreOperator::pull_back_through_reversal`],
//!   cleared with `x^N`, `N = B + max(h, 0) R` (`B` the bound on the degree
//!   in `x`, `R` on the shift), which is [`crate::pull_back`]'s
//!   [`crate::RowMapStep::Reversal`] operator for `h >= 0`;
//! - a power gauge is [`crate::LaurentOreOperator::pull_back_through_power_gauge`],
//!   cleared with `x^N`, `N = J + max(-h, 0) R` (`J` the bound on the
//!   derivative order), as `(D + w x^(-1))^d` has powers down to `x^(-d)`;
//! - the scale `s^(d_n) = s^(d_0) (s^h)^n` multiplies the term of shift `i`
//!   by `s^(h i)`;
//! - a Möbius map with `c != 0` is the composition (applied first to last)
//!   of `F(y) -> F(a/c + kappa y)`, `kappa = (b c - a e)/c`, the reversal
//!   at `d_n` and `K(u) -> K(c x + e)` (because
//!   `a/c + kappa/(c x + e) = (a x + b)/(c x + e)`), pulled back last to
//!   first with the reversal cleared before the second affine map (an
//!   affine substitution does not keep Laurent polynomials Laurent); with
//!   `c = 0` it is the constant scale `(s e)^(d_n)` and the affine map
//!   `x -> (a x + b)/e`.
//!
//! For fixed bounds ([`crate::PullbackBounds`]) the pullback is linear in
//! `M`, so it can be used inside an ansatz.
//!
//! # Forward rules
//!
//! The inverse of a substitution is a substitution
//! ([`Substitution::inverse`]): `x -> (x - b)/a` for an affine map; the
//! reversal is an involution; the power gauge `w_n -> -w_n`; and a Möbius
//! map `(a, b, c, e, s)` has the inverse `(e, -b, -c, a, 1/(s (a e - b c)))`
//! (from `c x + e = (a e - b c)/(a - c y)` at `x = (e y - b)/(a - c y)`).
//! Since a pullback is a conjugation up to a left factor,
//! `pull_back(L, C^(-1)) = q C L C^(-1)` annihilates `C F` when `L`
//! annihilates `F`: that is the forward rule [`push_forward_substitution`].
//! [`crate::DifferenceRowMapStep::Substitution`] puts substitutions into the
//! certified searches of [`crate::find_difference_image_annihilator`], whose
//! certificates replay the identity `w(n) P = H o L` on their own.

use crate::laurent::{LaurentError, LaurentOreOperator};
use crate::operator::OreOperator;
use crate::row_map::{pull_back, PullbackBounds, RowMapStep};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, ToPrimitive, Zero};
use std::fmt;

/// A substitution of the variable (module documentation).  Exponents are
/// ascending `[d_0, h]` with integer entries.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Substitution {
    /// `G_n(x) = F_n(slope x + intercept)`.
    Affine {
        /// `a != 0`.
        slope: Rational,
        /// `b`.
        intercept: Rational,
    },
    /// `G_n(x) = x^(d_n) F_n(1/x)`.
    Reversal {
        /// `d_n = d_0 + h n`.
        degree: Vec<Rational>,
    },
    /// `G_n(x) = (scale (c x + e))^(d_n) F_n((a x + b)/(c x + e))`.
    Mobius(Box<MobiusMap>),
    /// `G_n(x) = (1 + x)^(d_n) F_n(x/(1 + x))`.
    FPolynomial {
        /// `d_n = d_0 + h n`.
        degree: Vec<Rational>,
    },
    /// `G_n = x^(w_n) F_n`.
    PowerGauge {
        /// `w_n = w_0 + h n`.
        power: Vec<Rational>,
    },
}

/// The parameters of a weighted Möbius map
/// `G_n(x) = (scale (c x + e))^(d_n) F_n((a x + b)/(c x + e))`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MobiusMap {
    /// `a`.
    pub a: Rational,
    /// `b`.
    pub b: Rational,
    /// `c`.
    pub c: Rational,
    /// `e`.
    pub e: Rational,
    /// `s != 0`.
    pub scale: Rational,
    /// `d_n = d_0 + h n`.
    pub degree: Vec<Rational>,
}

/// Why a [`Substitution`] cannot be used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SubstitutionError {
    /// An affine slope or a Möbius scale is zero.
    ZeroParameter,
    /// A Möbius map is singular (`a e - b c = 0`).
    SingularMobius,
    /// An exponent is not `d_0 + h n` with integers `d_0`, `h`.
    InvalidExponent,
}

impl fmt::Display for SubstitutionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::ZeroParameter => "an affine slope or a Möbius scale must be nonzero",
            Self::SingularMobius => "a Möbius map needs a e - b c != 0",
            Self::InvalidExponent => "an exponent must be d_0 + h n with integers d_0 and h",
        })
    }
}

impl std::error::Error for SubstitutionError {}

impl From<LaurentError> for SubstitutionError {
    fn from(_: LaurentError) -> Self {
        Self::InvalidExponent
    }
}

/// `(d_0, h)` of an exponent with integer entries of degree at most one.
fn integer_exponent(exponent: &[Rational]) -> Option<(i64, i64)> {
    let exponent = poly::trim(exponent);
    if exponent.len() > 2 || !exponent.iter().all(Rational::is_integer) {
        return None;
    }
    let d0 = exponent[0].to_integer().to_i64()?;
    let h = exponent
        .get(1)
        .map_or(Some(0), |h| h.to_integer().to_i64())?;
    Some((d0, h))
}

/// The value of an exponent at `n`.
fn exponent_at(exponent: &[Rational], n: usize) -> Option<i64> {
    let (d0, h) = integer_exponent(exponent)?;
    d0.checked_add(h.checked_mul(i64::try_from(n).ok()?)?)
}

fn power(base: &Rational, exponent: i64) -> Rational {
    let factor = if exponent >= 0 {
        base.clone()
    } else {
        base.recip()
    };
    (0..exponent.unsigned_abs()).fold(Rational::one(), |acc, _| acc * &factor)
}

/// The Möbius parameters `(a, b, c, e, s, degree)`.
type MobiusParameters = (
    Rational,
    Rational,
    Rational,
    Rational,
    Rational,
    Vec<Rational>,
);

impl Substitution {
    /// `f(x) -> f(x + amount)`.
    pub fn translation(amount: Rational) -> Self {
        Self::Affine {
            slope: Rational::one(),
            intercept: amount,
        }
    }

    /// Check the parameters (module documentation).
    pub fn validate(&self) -> Result<(), SubstitutionError> {
        match self {
            Self::Affine { slope, .. } => {
                if slope.is_zero() {
                    return Err(SubstitutionError::ZeroParameter);
                }
            }
            Self::PowerGauge { power: exponent }
            | Self::Reversal { degree: exponent }
            | Self::FPolynomial { degree: exponent } => {
                integer_exponent(exponent).ok_or(SubstitutionError::InvalidExponent)?;
            }
            Self::Mobius(map) => {
                let MobiusMap {
                    a,
                    b,
                    c,
                    e,
                    scale,
                    degree,
                } = map.as_ref();
                if scale.is_zero() {
                    return Err(SubstitutionError::ZeroParameter);
                }
                if (a * e - b * c).is_zero() {
                    return Err(SubstitutionError::SingularMobius);
                }
                integer_exponent(degree).ok_or(SubstitutionError::InvalidExponent)?;
            }
        }
        Ok(())
    }

    /// The Möbius parameters of a reversal, an `f`-polynomial map or a
    /// Möbius map.
    fn mobius_parameters(&self) -> Option<MobiusParameters> {
        let (zero, one) = (Rational::zero(), Rational::one());
        Some(match self {
            Self::Reversal { degree } => (
                zero.clone(),
                one.clone(),
                one.clone(),
                zero,
                one,
                degree.clone(),
            ),
            Self::FPolynomial { degree } => (
                one.clone(),
                zero.clone(),
                one.clone(),
                one.clone(),
                one,
                degree.clone(),
            ),
            Self::Mobius(map) => (
                map.a.clone(),
                map.b.clone(),
                map.c.clone(),
                map.e.clone(),
                map.scale.clone(),
                map.degree.clone(),
            ),
            Self::Affine { .. } | Self::PowerGauge { .. } => return None,
        })
    }

    /// The inverse substitution (module documentation).
    pub fn inverse(&self) -> Self {
        match self {
            Self::Affine { slope, intercept } => Self::Affine {
                slope: slope.recip(),
                intercept: -(intercept / slope),
            },
            Self::Reversal { .. } => self.clone(),
            Self::PowerGauge { power } => Self::PowerGauge {
                power: power.iter().map(|c| -c).collect(),
            },
            Self::FPolynomial { .. } | Self::Mobius { .. } => {
                let (a, b, c, e, s, degree) = self.mobius_parameters().expect("Möbius type");
                let determinant = &a * &e - &b * &c;
                Self::Mobius(Box::new(MobiusMap {
                    a: e,
                    b: -b,
                    c: -c,
                    e: a,
                    scale: (s * determinant).recip(),
                    degree,
                }))
            }
        }
    }

    /// The image `G_n` of one row `F_n` (exact; `None` where it is not a
    /// polynomial or the parameters are invalid).  Canonical (zero is
    /// `[0]`).
    pub fn apply_row(&self, row: &[Rational], n: usize) -> Option<Vec<Rational>> {
        self.validate().ok()?;
        let row = poly::trim(row);
        Some(match self {
            Self::Affine { slope, intercept } => {
                poly::compose(&row, &[intercept.clone(), slope.clone()])
            }
            Self::PowerGauge { power } => {
                let w = exponent_at(power, n)?;
                if poly::is_zero(&row) {
                    return Some(row);
                }
                if w >= 0 {
                    poly::shift_by(&row, usize::try_from(w).ok()?)
                } else {
                    let k = usize::try_from(-w).ok()?;
                    if row.len() <= k || !row[..k].iter().all(Zero::is_zero) {
                        return None;
                    }
                    row[k..].to_vec()
                }
            }
            _ => {
                let (a, b, c, e, s, degree) = self.mobius_parameters()?;
                let d = usize::try_from(exponent_at(&degree, n)?).ok()?;
                if poly::degree(&row).is_some_and(|k| k > d) {
                    return None;
                }
                // Horner in the homogeneous form: with m = deg F,
                // H_m = f_m, H_k = H_(k+1) (a x + b) + f_k (c x + e)^(m - k),
                // and G = s^d (c x + e)^(d - m) H_0.
                let numerator = [b, a];
                let denominator = [e, c];
                let m = poly::degree(&row).unwrap_or(0);
                let mut denominator_power = vec![Rational::one()];
                let mut out = vec![row[m].clone()];
                for k in (0..m).rev() {
                    denominator_power = poly::mul(&denominator_power, &denominator);
                    out = poly::add(
                        &poly::mul(&out, &numerator),
                        &poly::scale(&denominator_power, &row[k]),
                    );
                }
                for _ in m..d {
                    out = poly::mul(&out, &denominator);
                }
                poly::scale(&out, &power(&s, d as i64))
            }
        })
    }
}

/// Multiply the term of shift `i` by `ratio^i` (the pullback through the
/// gauge `ratio^n`).
fn constant_gauge(op: &OreOperator, ratio: &Rational) -> OreOperator {
    let mut out = OreOperator::default();
    for ((i, d), c) in &op.0 {
        out.add_term((*i, *d), c.scale(&power(ratio, *i as i64)));
    }
    out
}

fn affine(op: &OreOperator, slope: Rational, intercept: Rational) -> OreOperator {
    let bounds = PullbackBounds::of(op);
    pull_back(op, &[RowMapStep::Affine { slope, intercept }], bounds)
        .expect("a valid affine step")
        .0
}

/// The reversal pullback cleared with `x^N`, `N = B + max(h, 0) R`.
fn reversal(
    op: &OreOperator,
    degree: &[Rational],
    bounds: PullbackBounds,
) -> Result<OreOperator, SubstitutionError> {
    let (_, h) = integer_exponent(degree).ok_or(SubstitutionError::InvalidExponent)?;
    let clearing = bounds.variable_degree + h.max(0) as usize * bounds.top_shift;
    let conjugate = LaurentOreOperator::from_ore_operator(op).pull_back_through_reversal(degree)?;
    Ok(conjugate
        .clear_with_power(clearing)
        .expect("the bounds clear every coefficient"))
}

/// The pullback of `op` through `substitution` with bounds `bounds` on `op`
/// (module documentation): `(M G)_n = 0` exactly when `(P F)_n = 0`.
///
/// # Errors
///
/// An invalid substitution ([`Substitution::validate`]).
///
/// # Panics
///
/// Panics when `bounds` is below the shape of `op`.
pub fn pull_back_substitution(
    op: &OreOperator,
    substitution: &Substitution,
    bounds: PullbackBounds,
) -> Result<OreOperator, SubstitutionError> {
    substitution.validate()?;
    let actual = PullbackBounds::of(op);
    assert!(
        op.is_zero()
            || (actual.top_shift <= bounds.top_shift
                && actual.derivative_order <= bounds.derivative_order
                && actual.variable_degree <= bounds.variable_degree),
        "pullback bounds below the operator's shape"
    );
    Ok(match substitution {
        Substitution::Affine { slope, intercept } => affine(op, slope.clone(), intercept.clone()),
        Substitution::PowerGauge { power } => {
            let (_, h) = integer_exponent(power).ok_or(SubstitutionError::InvalidExponent)?;
            let clearing = bounds.derivative_order + (-h).max(0) as usize * bounds.top_shift;
            LaurentOreOperator::from_ore_operator(op)
                .pull_back_through_power_gauge(power)?
                .clear_with_power(clearing)
                .expect("the bounds clear every coefficient")
        }
        Substitution::Reversal { degree } => reversal(op, degree, bounds)?,
        Substitution::FPolynomial { .. } | Substitution::Mobius { .. } => {
            let (a, b, c, e, s, degree) = substitution.mobius_parameters().expect("Möbius type");
            let (_, h) = integer_exponent(&degree).ok_or(SubstitutionError::InvalidExponent)?;
            if c.is_zero() {
                // (s e)^(d_n) F((a x + b)/e).
                let scaled = constant_gauge(op, &power(&(&s * &e), h));
                affine(&scaled, &a / &e, &b / &e)
            } else {
                let scaled = constant_gauge(op, &power(&s, h));
                let outer = affine(&scaled, c.clone(), e.clone());
                let reversed = reversal(&outer, &degree, bounds)?;
                let kappa = (&b * &c - &a * &e) / &c;
                affine(&reversed, kappa, &a / &c)
            }
        }
    })
}

/// The forward rule: an annihilator of `G = C F` from an annihilator `op` of
/// `F`, the pullback of `op` through the inverse substitution with the
/// exact bounds of `op` (module documentation).
///
/// # Errors
///
/// An invalid substitution.
pub fn push_forward_substitution(
    op: &OreOperator,
    substitution: &Substitution,
) -> Result<OreOperator, SubstitutionError> {
    substitution.validate()?;
    pull_back_substitution(op, &substitution.inverse(), PullbackBounds::of(op))
}

/// Upper bounds on the pullback of an operator within `bounds` (used to
/// chain steps).
pub(crate) fn bounds_after_substitution(
    bounds: PullbackBounds,
    substitution: &Substitution,
) -> PullbackBounds {
    let slope = |exponent: &[Rational]| {
        integer_exponent(exponent).map_or(0, |(_, h)| h.unsigned_abs() as usize)
    };
    let variable_degree = match substitution {
        Substitution::Affine { .. } => bounds.variable_degree,
        Substitution::PowerGauge { power } => {
            bounds.variable_degree + bounds.derivative_order + slope(power) * bounds.top_shift
        }
        Substitution::Mobius(map) if map.c.is_zero() => bounds.variable_degree,
        Substitution::Reversal { degree } | Substitution::FPolynomial { degree } => {
            bounds.variable_degree + 2 * bounds.derivative_order + slope(degree) * bounds.top_shift
        }
        Substitution::Mobius(map) => {
            bounds.variable_degree
                + 2 * bounds.derivative_order
                + slope(&map.degree) * bounds.top_shift
        }
    };
    PullbackBounds {
        variable_degree,
        ..bounds
    }
}
