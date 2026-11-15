//! Exact endpoint and fibre screens for real Hessenberg limit spectra.
//!
//! For `b(z)=z^-1+b_0+sum b_j z^j`, set `F(z)=z b(z)` and
//! `N(z)=z F'(z)-F(z)`.  The endpoint screen is Proposition 1.1 of the
//! limiting-spectrum notes: a real limit spectrum forces a critical point on
//! each side of zero, within the exact bound
//! `|R|^(s+1)|b_s| <= 1`.  The fibre screen is Corollary 1.2: between zero and
//! the nearest critical point, `F(z)-b(l)z` must have no root of smaller
//! modulus.  Candidate roots are numerical only; every accepted fibre root is
//! enclosed by an exact [`crate::complex_roots::RootDiscCertificate`].

use super::*;
use crate::bkw::Complex64;
use crate::complex_roots::{
    certify_root_disc, f64_to_rational, propose_roots_f64, refine_root, RootDiscCertificate,
};
use crate::sturm::SturmChain;
use num_bigint::BigInt;
use num_traits::{One, Signed, ToPrimitive, Zero};
use serde::{Deserialize, Serialize};

/// Provenance shared by endpoint and fibre screen certificates.
pub const HESSENBERG_SCREEN_PROVENANCE: &str =
    "endpoint theorem for real Hessenberg limit spectra (unpublished manuscript); Schmidt-Spitzer 1960";

/// The side of zero on which a screen is applied.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScreenSide {
    /// The interval to the right of zero.
    Right,
    /// The interval to the left of zero.
    Left,
}

/// Exact status of one endpoint or fibre test.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ScreenVerdict {
    /// The test supplies a necessary-condition violation.
    Refuted,
    /// The tested necessary condition was not violated; this is not a proof.
    Passes,
    /// No critical point on this side was available for the fibre test.
    NoCriticalPoint,
    /// The exact data did not settle this test.
    Inconclusive,
}

/// Replayable endpoint-bound result for one side.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EndpointScreen {
    /// Right or left side of zero.
    pub side: ScreenSide,
    /// Endpoint theorem verdict.
    pub verdict: ScreenVerdict,
    /// Rational `u` with `u^d |b_s| >= 1` used by the test.
    pub bound_u: Option<String>,
    /// An exact isolating interval for the nearest critical point, when found.
    pub critical_interval: Option<[String; 2]>,
}

/// Replayable fibre result for one side.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FiberScreen {
    /// Right or left side of zero.
    pub side: ScreenSide,
    /// Fibre theorem verdict.
    pub verdict: ScreenVerdict,
    /// Rational `l` between zero and the nearest critical point.
    pub l: Option<String>,
    /// Exact value `b(l)` when `l` was selected.
    pub b_at_l: Option<String>,
    /// Certified root-disc witness, if the fibre test refuted the condition.
    pub root_disc: Option<RootDiscCertificate>,
}

/// Exact witness used to refute a real limiting spectrum.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HessenbergScreenWitness {
    /// Side on which the necessary condition fails.
    pub side: ScreenSide,
    /// The endpoint or fibre test producing the witness.
    pub test: String,
    /// Endpoint `u`, present for an endpoint witness.
    pub bound_u: Option<String>,
    /// Fibre endpoint `l`, present for a fibre witness.
    pub l: Option<String>,
    /// Exact `b(l)` for a fibre witness.
    pub b_at_l: Option<String>,
    /// Exact root-disc certificate for a fibre witness.
    pub root_disc: Option<RootDiscCertificate>,
}

/// Replayable exact endpoint/fibre screen for a Hessenberg symbol.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HessenbergLimitSetScreenReport {
    /// Versioned report schema.
    pub schema: String,
    /// The exact theorem provenance, with no numerical attribution.
    pub provenance: String,
    /// Constant coefficient `b_0` of the symbol.
    pub b0: String,
    /// Exact tails `b_1,...,b_s` in ascending symbol order.
    pub tails: Vec<String>,
    /// Exact ascending coefficients of `F(z)=z b(z)`.
    pub f_coefficients: Vec<String>,
    /// Exact ascending coefficients of `N(z)=zF'(z)-F(z)`.
    pub critical_polynomial: Vec<String>,
    /// Optional caller-supplied Gaussian-rational centre seeds, retained so
    /// that a report produced with better numerical approximations replays.
    #[serde(default)]
    pub centers: Vec<[String; 2]>,
    /// The two side-by-side endpoint and fibre results.
    pub sides: Vec<HessenbergScreenSideReport>,
    /// The first selected exact refutation witness, if any.
    pub witness: Option<HessenbergScreenWitness>,
}

/// Endpoint and fibre results grouped by side.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct HessenbergScreenSideReport {
    /// Side of zero.
    pub side: ScreenSide,
    /// Endpoint-bound result.
    pub endpoint: EndpointScreen,
    /// Least-root fibre result.
    pub fiber: FiberScreen,
}

impl HessenbergLimitSetScreenReport {
    /// Whether this screen contains an exact refutation.
    pub fn is_refuted(&self) -> bool {
        self.witness.is_some()
    }

    /// A short stable name for the refuting theorem screen.
    pub fn refutation_name(&self) -> Option<&'static str> {
        self.witness.as_ref().map(|w| match w.test.as_str() {
            "endpoint" => "Hessenberg endpoint screen (Proposition 1.1)",
            "fiber" => "Hessenberg fiber screen (Corollary 1.2)",
            _ => "Hessenberg limit-spectrum screen",
        })
    }

    /// Recompute both exact screens and compare the complete report.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let centers: Vec<(String, String)> = self
            .centers
            .iter()
            .map(|center| (center[0].clone(), center[1].clone()))
            .collect();
        let expected =
            certify_hessenberg_limit_screens_with_centers(&self.b0, &self.tails, &centers)?
                .ok_or_else(|| {
                    RecurrenceError::VerificationFailed("screen has no Hessenberg tails".into())
                })?;
        if expected != *self {
            return Err(RecurrenceError::VerificationFailed(
                "Hessenberg endpoint/fiber screen does not replay".into(),
            ));
        }
        Ok(())
    }
}

/// Certify the endpoint and fibre screens from a Hessenberg symbol.
///
/// `b0` is `b_0`; `tails[j]` is `b_(j+1)`.  A passing test is only the absence
/// of this necessary-condition refutation and never proves a real limit set.
/// The returned report is self-contained and can be replayed with
/// [`HessenbergLimitSetScreenReport::verify`].
pub fn certify_hessenberg_limit_screens(
    b0: &str,
    tails: &[String],
) -> Result<Option<HessenbergLimitSetScreenReport>, RecurrenceError> {
    certify_hessenberg_limit_screens_with_centers(b0, tails, &[])
}

/// Certify the endpoint and fibre screens, adding exact caller-supplied
/// Gaussian-rational centre seeds to the numerical proposer.
///
/// The centres are `(real, imaginary)` rational or decimal strings.  They are
/// only seeds unless they themselves pass the exact root-disc certificate;
/// accepted reports retain them and replay them exactly.
pub fn certify_hessenberg_limit_screens_with_centers(
    b0: &str,
    tails: &[String],
    centers: &[(String, String)],
) -> Result<Option<HessenbergLimitSetScreenReport>, RecurrenceError> {
    let b0 = parse_screen_rational(b0)?;
    let tails = tails
        .iter()
        .map(|tail| parse_screen_rational(tail))
        .collect::<Result<Vec<_>, _>>()?;
    let centers = centers
        .iter()
        .map(|(real, imaginary)| {
            Ok((
                parse_screen_rational(real)?,
                parse_screen_rational(imaginary)?,
            ))
        })
        .collect::<Result<Vec<_>, RecurrenceError>>()?;
    certify_symbol(&b0, &tails, &centers)
}

fn parse_screen_rational(value: &str) -> Result<Q, RecurrenceError> {
    if let Ok(value) = parse_rational_coeff(value) {
        return Ok(value);
    }
    let trimmed = value.trim();
    let (mantissa, exponent) = trimmed
        .split_once(['e', 'E'])
        .map_or((trimmed, 0i32), |(mantissa, exponent)| {
            (mantissa, exponent.parse::<i32>().unwrap_or(i32::MIN))
        });
    if exponent == i32::MIN {
        return Err(parse_rational_coeff(value).unwrap_err().into());
    }
    let negative = mantissa.starts_with('-');
    let digits = mantissa.trim_start_matches(['+', '-']);
    let Some((whole, fraction)) = digits.split_once('.') else {
        return Err(parse_rational_coeff(value).unwrap_err().into());
    };
    if whole.is_empty() && fraction.is_empty()
        || !whole.chars().all(|c| c.is_ascii_digit())
        || !fraction.chars().all(|c| c.is_ascii_digit())
    {
        return Err(parse_rational_coeff(value).unwrap_err().into());
    }
    let mut numerator = BigInt::parse_bytes(format!("{}{}", whole, fraction).as_bytes(), 10)
        .unwrap_or_else(BigInt::zero);
    if negative {
        numerator = -numerator;
    }
    let decimal_places = fraction.len() as i32 - exponent;
    if decimal_places >= 0 {
        Ok(Q::new(
            numerator,
            BigInt::from(10).pow(decimal_places as u32),
        ))
    } else {
        Ok(Q::from_integer(
            numerator * BigInt::from(10).pow((-decimal_places) as u32),
        ))
    }
}

pub(super) fn certify_from_denominator(
    denominator: &[P],
) -> Result<Option<HessenbergLimitSetScreenReport>, RecurrenceError> {
    if denominator.len() < 3
        || denominator[1].degree() != Some(1)
        || denominator[1].coeff(1) != -Q::one()
        || denominator
            .iter()
            .skip(2)
            .any(|polynomial| polynomial.degree().is_some_and(|degree| degree > 0))
    {
        return Ok(None);
    }
    let tails: Vec<Q> = denominator
        .iter()
        .skip(2)
        .map(|polynomial| polynomial.coeff(0))
        .collect();
    let last_nonzero = tails.iter().rposition(|tail| !tail.is_zero());
    let Some(last_nonzero) = last_nonzero else {
        return Ok(None);
    };
    let tails = tails[..=last_nonzero].to_vec();
    certify_symbol(&denominator[1].coeff(0), &tails, &[])
}

fn certify_symbol(
    b0: &Q,
    tails: &[Q],
    centers: &[(Q, Q)],
) -> Result<Option<HessenbergLimitSetScreenReport>, RecurrenceError> {
    if tails.is_empty() || tails.iter().all(Zero::is_zero) {
        return Ok(None);
    }
    let mut f = vec![Q::one(), b0.clone()];
    f.extend(tails.iter().cloned());
    while f.last().is_some_and(Zero::is_zero) {
        f.pop();
    }
    let degree = f.len().saturating_sub(1);
    if degree < 2 {
        return Ok(None);
    }
    let mut critical = vec![Q::zero(); f.len()];
    critical[0] = -Q::one();
    for (index, coefficient) in tails.iter().enumerate() {
        critical[index + 2] = coefficient * Q::from_integer(BigInt::from(index + 1));
    }
    let integer_critical = vec_poly::clear_denominators(&critical);
    let sturm = SturmChain::from_bigint_coeffs(&integer_critical);
    let mut sides = Vec::new();
    for side in [ScreenSide::Right, ScreenSide::Left] {
        let bound_u = threshold_bound(tails.last().expect("degree >= 2"), degree);
        let (a, b) = match side {
            ScreenSide::Right => (Q::zero(), bound_u.clone()),
            ScreenSide::Left => (-bound_u.clone(), Q::zero()),
        };
        let endpoint_verdict = if closed_root_count(&critical, &sturm, &a, &b) == 0 {
            ScreenVerdict::Refuted
        } else {
            ScreenVerdict::Passes
        };
        let interval = nearest_interval(&critical, &sturm, &bound_u, side);
        let endpoint = EndpointScreen {
            side,
            verdict: endpoint_verdict,
            bound_u: Some(bound_u.to_string()),
            critical_interval: interval
                .as_ref()
                .map(|(lo, hi)| [lo.to_string(), hi.to_string()]),
        };
        let fiber = fiber_screen(&f, &critical, &sturm, side, interval.as_ref(), centers)?;
        sides.push(HessenbergScreenSideReport {
            side,
            endpoint,
            fiber,
        });
    }

    let witness = sides.iter().find_map(|side| {
        // Prefer the fibre witness when the endpoint bound passes: this is
        // the informative degree-40 case.  Otherwise retain the cheaper
        // endpoint witness.
        if side.endpoint.verdict == ScreenVerdict::Refuted
            && side.fiber.verdict != ScreenVerdict::Refuted
        {
            return Some(HessenbergScreenWitness {
                side: side.side,
                test: "endpoint".into(),
                bound_u: side.endpoint.bound_u.clone(),
                l: None,
                b_at_l: None,
                root_disc: None,
            });
        }
        if side.fiber.verdict == ScreenVerdict::Refuted {
            return Some(HessenbergScreenWitness {
                side: side.side,
                test: "fiber".into(),
                bound_u: None,
                l: side.fiber.l.clone(),
                b_at_l: side.fiber.b_at_l.clone(),
                root_disc: side.fiber.root_disc.clone(),
            });
        }
        if side.endpoint.verdict == ScreenVerdict::Refuted {
            Some(HessenbergScreenWitness {
                side: side.side,
                test: "endpoint".into(),
                bound_u: side.endpoint.bound_u.clone(),
                l: None,
                b_at_l: None,
                root_disc: None,
            })
        } else {
            None
        }
    });
    Ok(Some(HessenbergLimitSetScreenReport {
        schema: "polytool.hessenberg-limit-screen.v1".into(),
        provenance: HESSENBERG_SCREEN_PROVENANCE.into(),
        b0: b0.to_string(),
        tails: tails.iter().map(ToString::to_string).collect(),
        f_coefficients: f.iter().map(ToString::to_string).collect(),
        critical_polynomial: critical.iter().map(ToString::to_string).collect(),
        centers: centers
            .iter()
            .map(|(real, imaginary)| [real.to_string(), imaginary.to_string()])
            .collect(),
        sides,
        witness,
    }))
}

fn threshold_bound(leading_tail: &Q, degree: usize) -> Q {
    let absolute = leading_tail.abs();
    let mut low = Q::zero();
    let mut high = Q::one();
    while high.pow(degree as i32) * &absolute < Q::one() {
        high *= Q::from_integer(2.into());
    }
    for _ in 0..48 {
        let middle = (&low + &high) / Q::from_integer(2.into());
        if middle.pow(degree as i32) * &absolute < Q::one() {
            low = middle;
        } else {
            high = middle;
        }
    }
    high
}

fn evaluate(coefficients: &[Q], point: &Q) -> Q {
    coefficients
        .iter()
        .rev()
        .fold(Q::zero(), |value, coefficient| value * point + coefficient)
}

/// Distinct roots in the closed interval between `a` and `b`, in either order.
///
/// `SturmChain::count_roots_in` saturates to zero on a reversed interval, which
/// would make the right-side fibre guard `[0, l]` vacuous.
fn closed_root_count(coefficients: &[Q], sturm: &SturmChain, a: &Q, b: &Q) -> usize {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    sturm.count_roots_in(a, b)
        + usize::from(evaluate(coefficients, a).is_zero())
        + usize::from(evaluate(coefficients, b).is_zero())
}

fn nearest_interval(
    critical: &[Q],
    sturm: &SturmChain,
    bound: &Q,
    side: ScreenSide,
) -> Option<(Q, Q)> {
    let zero = Q::zero();
    let (left, right) = match side {
        ScreenSide::Right => (zero.clone(), bound.clone()),
        ScreenSide::Left => (-bound.clone(), zero.clone()),
    };
    if closed_root_count(critical, sturm, &left, &right) == 0 {
        return None;
    }
    // Locate only the first root from zero.  Unlike full isolation this does
    // not bisect the other 39 roots of the degree-40 critical polynomial.
    let mut low = Q::zero();
    let mut high = bound.clone();
    for _ in 0..12 {
        let middle = (&low + &high) / Q::from_integer(2.into());
        let (a, b) = match side {
            ScreenSide::Right => (Q::zero(), middle.clone()),
            ScreenSide::Left => (-middle.clone(), Q::zero()),
        };
        if closed_root_count(critical, sturm, &a, &b) > 0 {
            high = middle;
        } else {
            low = middle;
        }
    }
    Some(match side {
        ScreenSide::Right => (low, high),
        ScreenSide::Left => (-high, -low),
    })
}

fn choose_l(interval: &(Q, Q), side: ScreenSide) -> Q {
    let (a, b) = interval;
    match side {
        ScreenSide::Right => {
            if a.is_positive() {
                if a == b {
                    a / Q::from_integer(2.into())
                } else {
                    a.clone()
                }
            } else {
                b / Q::from_integer(2.into())
            }
        }
        ScreenSide::Left => {
            if b.is_negative() {
                if a == b {
                    b / Q::from_integer(2.into())
                } else {
                    b.clone()
                }
            } else {
                a / Q::from_integer(2.into())
            }
        }
    }
}

fn fiber_screen(
    f: &[Q],
    critical: &[Q],
    sturm: &SturmChain,
    side: ScreenSide,
    interval: Option<&(Q, Q)>,
    supplied_centers: &[(Q, Q)],
) -> Result<FiberScreen, RecurrenceError> {
    let Some(interval) = interval else {
        return Ok(FiberScreen {
            side,
            verdict: ScreenVerdict::NoCriticalPoint,
            l: None,
            b_at_l: None,
            root_disc: None,
        });
    };
    let l = choose_l(interval, side);
    let zero = Q::zero();
    if l.is_zero() || closed_root_count(critical, sturm, &l, &zero) != 0 {
        return Ok(FiberScreen {
            side,
            verdict: ScreenVerdict::Inconclusive,
            l: Some(l.to_string()),
            b_at_l: None,
            root_disc: None,
        });
    }
    let b_at_l = evaluate(f, &l) / &l;
    let mut polynomial = f.to_vec();
    if polynomial.len() > 1 {
        polynomial[1] -= &b_at_l;
    }
    let l_abs = l.abs();
    let l_float = l_abs.to_f64().unwrap_or(f64::INFINITY);
    let mut candidates = propose_roots_f64(&polynomial);
    // In tails-only input the small factor of F contains the useful initial
    // locations, while the fibre polynomial can have a poorly conditioned
    // cluster.  Keep both proposal sets; neither is a proof.
    candidates.extend(propose_roots_f64(f));
    // The bulk examples retain a conjugate pair just inside the critical
    // radius.  Supplement a failed global proposer with a small deterministic
    // local Newton stencil; any resulting point still has to pass the exact
    // disc certificate below.
    let sign = match side {
        ScreenSide::Right => 1.0,
        ScreenSide::Left => -1.0,
    };
    for real_fraction in [0.75, 0.85, 0.9, 0.95] {
        for imaginary_fraction in [-0.08, -0.04, -0.02, 0.02, 0.04, 0.08] {
            candidates.push(Complex64::new(
                sign * l_float * real_fraction,
                l_float * imaginary_fraction,
            ));
        }
    }
    // The fibre obstruction is just inside |l|, so prioritize the largest
    // admissible proposals.  A small-modulus proposal can be a perfectly
    // valid root but is irrelevant to this necessary-condition screen.
    candidates.sort_by(|left, right| right.abs().total_cmp(&left.abs()));
    let precision_levels = [256u32, 512, 1024, 2048, 4096, 8192];

    // A client may have a substantially better approximation than an f64
    // proposer.  Try these exact supplied centres first, but still apply the
    // same exact disc and modulus checks.
    for (center_re, center_im) in supplied_centers {
        if let Some(certificate) =
            try_fiber_center(&polynomial, center_re, center_im, &l_abs, l_float)
        {
            return Ok(FiberScreen {
                side,
                verdict: ScreenVerdict::Refuted,
                l: Some(l.to_string()),
                b_at_l: Some(b_at_l.to_string()),
                root_disc: Some(certificate),
            });
        }
    }

    // Refine the most promising proposals at increasing fixed-point
    // precision.  The radius list below is deliberately dyadic; as the
    // centre's residual shrinks, a smaller radius can make the c0-dominated
    // Rouche inequality strict while retaining the exact |z| < |l| check.
    for candidate in candidates
        .into_iter()
        .filter(|candidate| {
            candidate.re.is_finite()
                && candidate.im.is_finite()
                && candidate.abs() < l_float * (1.0 - 1e-10)
        })
        .take(24)
    {
        let candidate = refine_seed_f64(&polynomial, candidate);
        for bits in precision_levels {
            let Some((center_re, center_im)) =
                refine_root(&polynomial, candidate.re, candidate.im, bits)
            else {
                break;
            };
            if let Some(certificate) =
                try_fiber_center(&polynomial, &center_re, &center_im, &l_abs, l_float)
            {
                return Ok(FiberScreen {
                    side,
                    verdict: ScreenVerdict::Refuted,
                    l: Some(l.to_string()),
                    b_at_l: Some(b_at_l.to_string()),
                    root_disc: Some(certificate),
                });
            }
        }
    }
    Ok(FiberScreen {
        side,
        verdict: ScreenVerdict::Passes,
        l: Some(l.to_string()),
        b_at_l: Some(b_at_l.to_string()),
        root_disc: None,
    })
}

fn try_fiber_center(
    polynomial: &[Q],
    center_re: &Q,
    center_im: &Q,
    l_abs: &Q,
    l_float: f64,
) -> Option<RootDiscCertificate> {
    let center_norm_squared = center_re * center_re + center_im * center_im;
    if center_norm_squared >= l_abs * l_abs {
        return None;
    }
    let center_abs = center_norm_squared.to_f64()?.sqrt();
    let gap = l_float - center_abs;
    if !(gap.is_finite() && gap > 0.0) {
        return None;
    }
    let gap = f64_to_rational(gap / 4.0)?;
    for exponent in [
        8u32, 12, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128, 160, 192, 256,
    ] {
        let dyadic_radius = Q::new(BigInt::one(), BigInt::one() << exponent as usize);
        let radius = if gap < dyadic_radius {
            gap.clone()
        } else {
            dyadic_radius
        };
        let Some(certificate) = certify_root_disc(polynomial, center_re, center_im, &radius) else {
            continue;
        };
        if certificate
            .maximal_modulus_upper_bound
            .parse::<Q>()
            .map(|upper| upper < *l_abs)
            .unwrap_or(false)
        {
            return Some(certificate);
        }
    }
    None
}

fn refine_seed_f64(coefficients: &[Q], mut root: Complex64) -> Complex64 {
    let scale = coefficients
        .iter()
        .map(|coefficient| coefficient.to_f64().unwrap_or(0.0).abs())
        .fold(0.0, f64::max);
    if !(scale.is_finite() && scale > 0.0) {
        return root;
    }
    let values: Vec<f64> = coefficients
        .iter()
        .map(|coefficient| coefficient.to_f64().unwrap_or(0.0) / scale)
        .collect();
    let derivative: Vec<f64> = values
        .iter()
        .enumerate()
        .skip(1)
        .map(|(degree, coefficient)| degree as f64 * coefficient)
        .collect();
    for _ in 0..120 {
        let value = evaluate_seed_f64(&values, root);
        let slope = evaluate_seed_f64(&derivative, root);
        if !value.re.is_finite() || !value.im.is_finite() || slope.abs() < 1e-30 {
            break;
        }
        let correction = value / slope;
        root = root - correction;
        if correction.abs() < 1e-14 {
            break;
        }
    }
    root
}

fn evaluate_seed_f64(coefficients: &[f64], point: Complex64) -> Complex64 {
    coefficients
        .iter()
        .rev()
        .fold(Complex64::ZERO, |value, coefficient| {
            value * point + Complex64::new(*coefficient, 0.0)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> String {
        n.to_string()
    }

    #[test]
    fn endpoint_screen_refutes_the_hand_checked_symbol() {
        let report =
            certify_hessenberg_limit_screens("1", &[q(5) + "/9", q(1) + "/3", q(2) + "/27"])
                .unwrap()
                .unwrap();
        assert!(report.is_refuted());
        assert_eq!(report.witness.as_ref().unwrap().side, ScreenSide::Left);
        report.verify().unwrap();
    }

    #[test]
    fn closed_root_count_accepts_either_order() {
        let half = Q::new(1.into(), 2.into());
        let critical = vec![-half, Q::one()];
        let sturm = SturmChain::from_bigint_coeffs(&vec_poly::clear_denominators(&critical));
        let (zero, one) = (Q::zero(), Q::one());
        assert_eq!(closed_root_count(&critical, &sturm, &zero, &one), 1);
        assert_eq!(closed_root_count(&critical, &sturm, &one, &zero), 1);
    }

    #[test]
    fn tridiagonal_boundary_is_not_refuted() {
        let report = certify_hessenberg_limit_screens("0", &["1".into()])
            .unwrap()
            .unwrap();
        assert!(!report.is_refuted());
        assert!(report.sides.iter().all(|side| {
            side.endpoint.verdict == ScreenVerdict::Passes
                && side.fiber.verdict == ScreenVerdict::Passes
        }));
    }
}
