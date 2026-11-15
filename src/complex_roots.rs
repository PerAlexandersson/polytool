//! Exact Rouche certificates for roots in Gaussian-rational discs.
//!
//! [`certify_root_disc`] uses an exact Taylor expansion at a centre in
//! `Q(i)`.  It proves that the linear Taylor term dominates all remaining
//! terms on the boundary of the disc, and therefore that the disc contains
//! exactly one root.  The optional [`propose_roots_f64`] helper is only a
//! numerical source of candidate centres; no floating-point value enters a
//! certificate.

use crate::{BigRational, Complex64};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

#[derive(Clone, Debug, PartialEq, Eq)]
struct GaussianRational {
    re: Q,
    im: Q,
}

impl GaussianRational {
    fn new(re: Q, im: Q) -> Self {
        Self { re, im }
    }

    fn zero() -> Self {
        Self::new(Q::zero(), Q::zero())
    }

    fn one() -> Self {
        Self::new(Q::one(), Q::zero())
    }

    fn add(&self, other: &Self) -> Self {
        Self::new(&self.re + &other.re, &self.im + &other.im)
    }

    fn scale(&self, factor: &Q) -> Self {
        Self::new(&self.re * factor, &self.im * factor)
    }

    fn mul(&self, other: &Self) -> Self {
        Self::new(
            &self.re * &other.re - &self.im * &other.im,
            &self.re * &other.im + &self.im * &other.re,
        )
    }

    fn norm_squared(&self) -> Q {
        &self.re * &self.re + &self.im * &self.im
    }
}

/// A Taylor coefficient in `Q(i)`, stored as exact rational strings.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct GaussianRationalCoefficient {
    /// Real part of the coefficient.
    pub real: String,
    /// Imaginary part of the coefficient.
    pub imaginary: String,
}

impl GaussianRationalCoefficient {
    fn from_exact(value: &GaussianRational) -> Self {
        Self {
            real: value.re.to_string(),
            imaginary: value.im.to_string(),
        }
    }
}

/// Exact, replayable proof that a disc contains exactly one root.
///
/// The input polynomial has real rational coefficients in ascending degree
/// order.  The centre and all bounds are rational; the proposer in
/// [`propose_roots_f64`] is deliberately not part of this proof.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct RootDiscCertificate {
    /// The Gaussian-rational centre's real part.
    pub center_re: String,
    /// The Gaussian-rational centre's imaginary part.
    pub center_im: String,
    /// Positive rational disc radius.
    pub radius: String,
    /// Exact Taylor coefficients at the centre, in ascending powers of the
    /// local variable `w = z - center`.
    pub taylor_coefficients: Vec<GaussianRationalCoefficient>,
    /// A rational upper bound for the modulus of the constant Taylor term.
    pub constant_modulus_upper_bound: String,
    /// A rational lower bound for the modulus of the linear Taylor term.
    pub linear_modulus_lower_bound: String,
    /// A rational upper bound for the sum of all Taylor terms of degree at
    /// least two on the disc boundary.
    pub remainder_modulus_upper_bound: String,
    /// The strict left-hand side `|c_1| r` certified by the inequality.
    pub linear_boundary_lower_bound: String,
    /// The strict right-hand side bound `|c_0| + sum_(k>=2)|c_k|r^k`.
    pub non_linear_boundary_upper_bound: String,
    /// An exact upper bound for `max{|z| : |z-center| <= radius}`.
    pub maximal_modulus_upper_bound: String,
}

impl RootDiscCertificate {
    /// Recompute the Taylor data and Rouche inequality from the polynomial.
    ///
    /// This accepts the same ascending coefficient convention as
    /// [`certify_root_disc`].  A certificate verifies only when every stored
    /// exact field, including its safe square-root bounds, is reproduced.
    pub fn verify(&self, coeffs: &[Q]) -> bool {
        let Some(center_re) = parse_q(&self.center_re) else {
            return false;
        };
        let Some(center_im) = parse_q(&self.center_im) else {
            return false;
        };
        let Some(radius) = parse_q(&self.radius) else {
            return false;
        };
        let Some(expected) = certify_root_disc(coeffs, &center_re, &center_im, &radius) else {
            return false;
        };
        expected == *self
    }
}

/// Certify exactly one root in the open Gaussian-rational disc.
///
/// The strict condition is
/// `|c_1| r > |c_0| + sum_(k>=2) |c_k| r^k`, where `c_k` are the exact
/// Taylor coefficients at the supplied centre.  Moduli on the right are
/// bounded above and the linear modulus is bounded below by adaptive,
/// relative-precision rational square-root bounds.  Thus `Some` proves by
/// Rouche's theorem that the open disc contains exactly one root, counted
/// with multiplicity.  `None` means only that this particular test did not
/// certify the disc.
pub fn certify_root_disc(
    coeffs: &[Q],
    center_re: &Q,
    center_im: &Q,
    radius: &Q,
) -> Option<RootDiscCertificate> {
    if radius <= &Q::zero() || coeffs.is_empty() || coeffs.iter().all(Zero::is_zero) {
        return None;
    }
    let center = GaussianRational::new(center_re.clone(), center_im.clone());
    let coefficients = taylor_coefficients(coeffs, &center);
    if coefficients.len() < 2 || coefficients[1].norm_squared().is_zero() {
        return None;
    }

    let c0_upper = sqrt_upper(&coefficients[0].norm_squared());
    let c1_lower = sqrt_lower(&coefficients[1].norm_squared());
    if c1_lower.is_zero() {
        return None;
    }
    let mut remainder = Q::zero();
    let mut radius_power = radius.clone() * radius;
    for coefficient in coefficients.iter().skip(2) {
        remainder += sqrt_upper(&coefficient.norm_squared()) * &radius_power;
        radius_power *= radius;
    }
    let right = &c0_upper + &remainder;
    let left = &c1_lower * radius;
    if left <= right {
        return None;
    }

    let maximal_modulus_upper_bound = sqrt_upper(&center.norm_squared()) + radius;
    Some(RootDiscCertificate {
        center_re: center_re.to_string(),
        center_im: center_im.to_string(),
        radius: radius.to_string(),
        taylor_coefficients: coefficients
            .iter()
            .map(GaussianRationalCoefficient::from_exact)
            .collect(),
        constant_modulus_upper_bound: c0_upper.to_string(),
        linear_modulus_lower_bound: c1_lower.to_string(),
        remainder_modulus_upper_bound: remainder.to_string(),
        linear_boundary_lower_bound: left.to_string(),
        non_linear_boundary_upper_bound: right.to_string(),
        maximal_modulus_upper_bound: maximal_modulus_upper_bound.to_string(),
    })
}

fn taylor_coefficients(coeffs: &[Q], center: &GaussianRational) -> Vec<GaussianRational> {
    let degree = coeffs.iter().rposition(|c| !c.is_zero()).unwrap_or(0);
    let mut powers = vec![GaussianRational::one(); degree + 1];
    for power in 1..=degree {
        powers[power] = powers[power - 1].mul(center);
    }
    let mut result = Vec::with_capacity(degree + 1);
    for k in 0..=degree {
        let mut coefficient = GaussianRational::zero();
        for j in k..=degree {
            let binomial = binomial(j, k);
            coefficient = coefficient.add(&powers[j - k].scale(&(&coeffs[j] * binomial)));
        }
        result.push(coefficient);
    }
    result
}

fn binomial(n: usize, k: usize) -> Q {
    let mut value = BigInt::one();
    for j in 1..=k.min(n - k) {
        value = value * BigInt::from(n + 1 - j) / BigInt::from(j);
    }
    Q::from_integer(value)
}

fn parse_q(value: &str) -> Option<Q> {
    value.parse().ok()
}

fn sqrt_lower(value: &Q) -> Q {
    sqrt_scaled(value, false)
}

fn sqrt_upper(value: &Q) -> Q {
    sqrt_scaled(value, true)
}

fn sqrt_scaled(value: &Q, upper: bool) -> Q {
    debug_assert!(!value.is_negative());
    if value.is_zero() {
        return Q::zero();
    }
    let numerator_bits = value.numer().bits() as i64;
    let denominator_bits = value.denom().bits() as i64;
    // The precision is relative to the rational's bit scale, so a value such
    // as 10^-150 still receives hundreds of effective fractional bits.  A
    // fixed absolute precision would incorrectly turn such values into zero.
    let scale_bits = (64 + (denominator_bits - numerator_bits) / 2).max(0) as usize;
    let scaled_numerator = value.numer() << (2 * scale_bits);
    let floor = (&scaled_numerator / value.denom()).sqrt();
    let rounded = if upper { floor + BigInt::one() } else { floor };
    Q::new(rounded, BigInt::one() << scale_bits)
}

/// Refine a numerical complex root of a real rational polynomial.
///
/// The input coefficients are in ascending degree order.  Newton steps are
/// performed with dyadic rationals, rounding every fixed-point operation to
/// the current common `2^-precision` scale.  The working precision starts at
/// 64 bits and doubles until `bits`; the returned pair is therefore suitable
/// for an exact certificate, while the initial `f64` values are only seeds.
/// `None` means that the derivative vanished, Newton escaped, or the
/// requested convergence criterion was not reached.
pub fn refine_root(coeffs: &[Q], start_re: f64, start_im: f64, bits: u32) -> Option<(Q, Q)> {
    if coeffs.len() < 2
        || coeffs.iter().all(Zero::is_zero)
        || !start_re.is_finite()
        || !start_im.is_finite()
        || bits == 0
    {
        return None;
    }
    let target_bits = bits.max(64);
    let mut point = FixedGaussian::new(
        round_to_dyadic(&f64_to_rational(start_re)?, 64),
        round_to_dyadic(&f64_to_rational(start_im)?, 64),
    );
    let mut precision = 64u32;
    loop {
        point = point.rounded(precision);
        let rounded_coefficients: Vec<Q> = coeffs
            .iter()
            .map(|coefficient| round_to_dyadic(coefficient, precision))
            .collect();
        let derivative: Vec<Q> = rounded_coefficients
            .iter()
            .enumerate()
            .skip(1)
            .map(|(degree, coefficient)| {
                round_to_dyadic(
                    &(coefficient * Q::from_integer(BigInt::from(degree))),
                    precision,
                )
            })
            .collect();
        let mut converged = false;
        let mut previous_step: Option<Q> = None;
        for _ in 0..160 {
            let value = evaluate_fixed(&rounded_coefficients, &point, precision);
            let slope = evaluate_fixed(&derivative, &point, precision);
            let step = value.div(&slope, precision)?;
            let step_size = step.norm_squared();
            let point_size = point.norm_squared().max(Q::one());
            let threshold_exponent = precision.saturating_sub(8) as usize;
            let threshold_squared =
                Q::new(BigInt::one(), BigInt::one() << (2 * threshold_exponent));
            point = point.sub(&step, precision);
            if !point.is_bounded() {
                return None;
            }
            if step_size <= threshold_squared * point_size {
                converged = true;
                break;
            }
            if previous_step
                .as_ref()
                .is_some_and(|previous| step_size > previous * Q::from_integer(256.into()))
            {
                return None;
            }
            previous_step = Some(step_size);
        }
        if !converged {
            return None;
        }
        if precision >= target_bits {
            return Some((point.re, point.im));
        }
        precision = (precision.saturating_mul(2)).min(target_bits);
    }
}

#[derive(Clone, Debug)]
struct FixedGaussian {
    re: Q,
    im: Q,
}

impl FixedGaussian {
    fn new(re: Q, im: Q) -> Self {
        Self { re, im }
    }

    fn rounded(&self, bits: u32) -> Self {
        Self::new(
            round_to_dyadic(&self.re, bits),
            round_to_dyadic(&self.im, bits),
        )
    }

    fn add(&self, other: &Self, bits: u32) -> Self {
        Self::new(
            round_to_dyadic(&(&self.re + &other.re), bits),
            round_to_dyadic(&(&self.im + &other.im), bits),
        )
    }

    fn sub(&self, other: &Self, bits: u32) -> Self {
        Self::new(
            round_to_dyadic(&(&self.re - &other.re), bits),
            round_to_dyadic(&(&self.im - &other.im), bits),
        )
    }

    fn mul(&self, other: &Self, bits: u32) -> Self {
        let real = round_to_dyadic(&(&self.re * &other.re - &self.im * &other.im), bits);
        let imaginary = round_to_dyadic(&(&self.re * &other.im + &self.im * &other.re), bits);
        Self::new(real, imaginary)
    }

    fn div(&self, other: &Self, bits: u32) -> Option<Self> {
        let denominator = round_to_dyadic(&(&other.re * &other.re + &other.im * &other.im), bits);
        if denominator.is_zero() {
            return None;
        }
        Some(Self::new(
            round_to_dyadic(
                &((&self.re * &other.re + &self.im * &other.im) / &denominator),
                bits,
            ),
            round_to_dyadic(
                &((&self.im * &other.re - &self.re * &other.im) / &denominator),
                bits,
            ),
        ))
    }

    fn norm_squared(&self) -> Q {
        &self.re * &self.re + &self.im * &self.im
    }

    fn is_bounded(&self) -> bool {
        let limit = Q::from_integer(BigInt::one() << 1024usize);
        self.norm_squared() < limit
    }
}

fn evaluate_fixed(coeffs: &[Q], point: &FixedGaussian, bits: u32) -> FixedGaussian {
    let mut value = FixedGaussian::new(Q::zero(), Q::zero());
    for coefficient in coeffs.iter().rev() {
        value = value.mul(point, bits).add(
            &FixedGaussian::new(round_to_dyadic(coefficient, bits), Q::zero()),
            bits,
        );
    }
    value
}

fn round_to_dyadic(value: &Q, bits: u32) -> Q {
    if value.is_zero() {
        return Q::zero();
    }
    let scale = BigInt::one() << bits as usize;
    let numerator = value.numer() * &scale;
    let denominator = value.denom();
    let (quotient, remainder) = numerator.div_rem(denominator);
    let twice_remainder = remainder.abs() * 2;
    let rounded = if twice_remainder >= *denominator {
        quotient
            + if numerator.is_negative() {
                -BigInt::one()
            } else {
                BigInt::one()
            }
    } else {
        quotient
    };
    Q::new(rounded, scale)
}

pub(crate) fn f64_to_rational(value: f64) -> Option<Q> {
    if !value.is_finite() {
        return None;
    }
    if value == 0.0 {
        return Some(Q::zero());
    }
    let bits = value.to_bits();
    let negative = bits >> 63 != 0;
    let exponent = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1u64 << 52) - 1);
    let (mantissa, power) = if exponent == 0 {
        (fraction, -1074)
    } else {
        ((1u64 << 52) | fraction, exponent - 1023 - 52)
    };
    let integer = BigInt::from(mantissa);
    let value = if power >= 0 {
        Q::from_integer(integer << power as usize)
    } else {
        Q::new(integer, BigInt::one() << (-power) as usize)
    };
    Some(if negative { -value } else { value })
}

/// Propose roots using a floating-point Durand--Kerner iteration.
///
/// The returned centres are numerical guesses only.  They carry no
/// certification and must be passed through [`certify_root_disc`] before any
/// mathematical conclusion is drawn.  Coefficients are still accepted in
/// ascending exact rational order so that only the proposer performs the
/// lossy conversion.
pub fn propose_roots_f64(coeffs: &[Q]) -> Vec<Complex64> {
    let degree = coeffs.iter().rposition(|c| !c.is_zero()).unwrap_or(0);
    if degree == 0 {
        return Vec::new();
    }
    let values: Vec<f64> = coeffs
        .iter()
        .take(degree + 1)
        .map(|c| {
            c.to_f64().unwrap_or(if c.is_negative() {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            })
        })
        .collect();
    if values.iter().any(|c| !c.is_finite()) {
        return Vec::new();
    }
    let leading = values[degree];
    if leading == 0.0 {
        return Vec::new();
    }
    if degree == 1 {
        return vec![Complex64::new(-values[0] / leading, 0.0)];
    }
    // Balance the variable before forming the monic polynomial.  Bulk
    // symbols can have a tiny leading coefficient while all their roots are
    // close to one; the unscaled Cauchy radius would then overflow a
    // Durand--Kerner initialisation even though the roots themselves are
    // modest.  This changes only the numerical proposal variable w=z/scale.
    let root_scale = if values[0] != 0.0 {
        (values[0].abs() / leading.abs())
            .powf(1.0 / degree as f64)
            .max(1e-12)
    } else {
        1.0
    };
    let mut scale_power = 1.0;
    let scaled: Vec<f64> = values
        .iter()
        .map(|coefficient| {
            let result = coefficient * scale_power;
            scale_power *= root_scale;
            result
        })
        .collect();
    let scaled_leading = scaled[degree];
    if !scaled_leading.is_finite() || scaled_leading == 0.0 {
        return Vec::new();
    }
    let monic: Vec<f64> = scaled.iter().map(|c| c / scaled_leading).collect();
    let radius = 1.0 + monic[..degree].iter().map(|c| c.abs()).fold(0.0, f64::max);
    if !radius.is_finite() {
        return Vec::new();
    }
    let mut roots: Vec<Complex64> = (0..degree)
        .map(|k| {
            Complex64::from_polar(
                radius,
                0.37 + 2.0 * std::f64::consts::PI * k as f64 / degree as f64,
            )
        })
        .collect();
    for _ in 0..2000 {
        let old = roots.clone();
        let mut maximum_delta: f64 = 0.0;
        for i in 0..degree {
            let mut denominator = Complex64::ONE;
            for j in 0..degree {
                if i != j {
                    denominator = denominator * (old[i] - old[j]);
                }
            }
            if denominator.abs() < 1e-30 {
                denominator = denominator + Complex64::new(1e-14, 1e-14);
            }
            let delta = evaluate_complex(&monic, old[i]) / denominator;
            roots[i] = old[i] - delta;
            maximum_delta = maximum_delta.max(delta.abs());
        }
        if maximum_delta <= 1e-13 {
            break;
        }
    }
    if roots
        .iter()
        .any(|root| !root.re.is_finite() || !root.im.is_finite())
    {
        return newton_fallback(&monic, root_scale, radius, degree);
    }
    roots.into_iter().map(|root| root * root_scale).collect()
}

fn newton_fallback(monic: &[f64], root_scale: f64, radius: f64, degree: usize) -> Vec<Complex64> {
    let derivative: Vec<f64> = monic
        .iter()
        .enumerate()
        .skip(1)
        .map(|(index, coefficient)| index as f64 * coefficient)
        .collect();
    let count = (degree * 2).max(24);
    let mut seeds = (0..count)
        .map(|index| {
            Complex64::from_polar(
                radius * 0.8,
                0.19 + 2.0 * std::f64::consts::PI * index as f64 / count as f64,
            )
        })
        .collect::<Vec<_>>();
    for real_index in -20..=20 {
        for imaginary_index in -12..=12 {
            seeds.push(Complex64::new(
                real_index as f64 * 0.125,
                imaginary_index as f64 * 0.125,
            ));
        }
    }
    // A repeated or nearly repeated real cluster often defeats a circular
    // Durand--Kerner start.  These local perturbations are still merely
    // numerical seeds; they help Newton reach a nearby nonreal pair.
    for real_index in -8..=8 {
        for imaginary_index in -8..=8 {
            seeds.push(Complex64::new(
                -1.0 + real_index as f64 * 0.025,
                imaginary_index as f64 * 0.025,
            ));
        }
    }
    // The fallback seeds live in the balanced variable w=z/root_scale.  A
    // root near -1 in the original variable is therefore near
    // -1/root_scale, not necessarily near -1 in w.  This matters for the
    // repeated bulk factor in the degree-40 fibre fixture.
    for real_center in [-1.1, -1.05, -1.0, -0.95, -0.9] {
        for real_index in -8..=8 {
            for imaginary_index in -4..=4 {
                seeds.push(Complex64::new(
                    (real_center + real_index as f64 * 0.025) / root_scale,
                    imaginary_index as f64 * 0.025 / root_scale,
                ));
            }
        }
    }
    let mut roots = Vec::new();
    for mut root in seeds {
        let mut converged = false;
        for _ in 0..160 {
            let value = evaluate_complex(monic, root);
            let slope = evaluate_complex(&derivative, root);
            if !value.re.is_finite() || !value.im.is_finite() || slope.abs() < 1e-14 {
                break;
            }
            let correction = value / slope;
            root = root - correction;
            if correction.abs() < 1e-12 {
                converged = true;
                break;
            }
        }
        let residual = evaluate_complex(monic, root).abs();
        if root.re.is_finite()
            && root.im.is_finite()
            // The residual is evaluated in a high-degree monic scaling and
            // can lose several digits near a bulk cluster.  The subsequent
            // exact Rouche test is the acceptance gate, so this is only a
            // permissive filter for Newton's numerical proposals.
            && (residual < 1e-3 || converged)
            && !roots
                .iter()
                .any(|old: &Complex64| (*old - root).abs() < 1e-5)
        {
            roots.push(root * root_scale);
            if roots.len() == degree {
                break;
            }
        }
    }
    roots
}

fn evaluate_complex(coeffs: &[f64], point: Complex64) -> Complex64 {
    let mut value = Complex64::ZERO;
    for coefficient in coeffs.iter().rev() {
        value = value * point + Complex64::new(*coefficient, 0.0);
    }
    value
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(n.into())
    }

    fn frac(n: i64, d: i64) -> Q {
        Q::new(n.into(), d.into())
    }

    #[test]
    fn certifies_x_squared_plus_one_at_i() {
        let certificate = certify_root_disc(&[q(1), q(0), q(1)], &q(0), &q(1), &frac(1, 2))
            .expect("Rouche disc around i");
        assert!(certificate.verify(&[q(1), q(0), q(1)]));
        assert!(
            certificate
                .maximal_modulus_upper_bound
                .parse::<Q>()
                .unwrap()
                > Q::new(3.into(), 2.into())
        );
    }

    #[test]
    fn rejects_a_disc_that_is_too_large() {
        assert!(certify_root_disc(&[q(1), q(0), q(1)], &q(0), &q(0), &q(1)).is_none());
    }

    #[test]
    fn certifies_a_tightly_clustered_nonreal_pair() {
        let coefficients = vec![Q::new(101.into(), 100.into()), q(-2), q(1)];
        let certificate = certify_root_disc(&coefficients, &q(1), &frac(1, 10), &frac(1, 50))
            .expect("disc around one member of a close conjugate pair");
        assert!(certificate.verify(&coefficients));
    }

    #[test]
    fn scales_to_tiny_coefficients() {
        let scale = Q::new(BigInt::one(), BigInt::from(10).pow(150));
        let certificate = certify_root_disc(
            &[scale.clone(), q(0), scale.clone()],
            &q(0),
            &q(1),
            &frac(1, 2),
        )
        .expect("relative square-root bounds must retain tiny values");
        assert!(certificate.verify(&[scale.clone(), q(0), scale]));
    }

    #[test]
    fn proposes_a_cluster_without_claiming_certification() {
        let roots = propose_roots_f64(&[
            Q::from_integer(1.into()),
            Q::from_integer((-6).into()),
            Q::from_integer(11.into()),
            Q::from_integer((-6).into()),
            Q::from_integer(1.into()),
        ]);
        assert_eq!(roots.len(), 4);
        assert!(roots
            .iter()
            .all(|root| root.re.is_finite() && root.im.is_finite()));
    }

    #[test]
    fn refines_a_complex_root_in_fixed_point() {
        let root = refine_root(&[q(1), q(0), q(1)], 0.07, 0.93, 128)
            .expect("Newton should refine the seed near i");
        let tolerance = Q::new(BigInt::one(), BigInt::one() << 90usize);
        assert!(root.0.abs() < tolerance);
        assert!((root.1 - q(1)).abs() < tolerance);
    }
}
