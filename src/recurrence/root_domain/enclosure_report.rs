//! Exact interval interface and fixed-cap diagonal multiplier tests.
//!
//! [`root_enclosures`] answers a [`RootEnclosureRequest`]: for each supplied
//! polynomial (ascending rational strings) it returns rational isolating
//! enclosures with multiplicities, the degree, the valuation and the exact
//! real-rootedness verdict, optionally refined to a maximal width and tested
//! against an interval (constant, index-dependent or a ray).  Enclosures come
//! from [`crate::real_rootedness::isolate_real_roots_rational_coeffs`]; their midpoints
//! are representatives only, never certified bounds or roots.
//!
//! [`jensen_fixed_cap_check`] tests a diagonal map `x^k -> lambda_k x^k` at a
//! fixed degree cap `N` through the Jensen polynomial
//! `J_N(x) = sum_(k <= N) C(N, k) lambda_k x^k = T[(1 + x)^N]`.  This is
//! evidence at the stated cap only, kept separate from a proof for every
//! degree (a multiplier sequence needs every cap).

use super::*;
use crate::real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval};
use crate::vec_poly;
use num_traits::{Signed, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Schema tag of [`RootEnclosureReport`].
pub const ROOT_ENCLOSURES_SCHEMA: &str = "polytool.root-enclosures.v1";

/// Evidence label of [`JensenFixedCapReport`]: a statement for degrees up to the cap only.
pub const FIXED_CAP_EVIDENCE: &str = "fixed_cap";

/// Largest number of polynomials in a [`RootEnclosureRequest`].
pub const ROOT_ENCLOSURES_MAX_POLYNOMIALS: usize = 256;

/// Largest degree accepted by [`root_enclosures`] and [`jensen_fixed_cap_check`].
pub const ROOT_ENCLOSURES_MAX_DEGREE: usize = 1024;

/// Request for [`root_enclosures`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEnclosureRequest {
    /// Polynomials as ascending rational strings (`"3/2"`).
    pub polynomials: Vec<Vec<String>>,
    /// Refine every enclosure to at most this positive rational width.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_width: Option<String>,
    /// Interval to test; index-dependent endpoints use `first_index + position`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interval: Option<RootInterval>,
    /// Index of the first polynomial (default 0).
    #[serde(default)]
    pub first_index: usize,
}

/// Exact root data of one polynomial.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEnclosureItem {
    /// Index `first_index + position`.
    pub index: usize,
    /// Degree; omitted for the zero polynomial.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    /// Multiplicity of the root 0.
    pub valuation: usize,
    /// Every root real (true for zero and nonzero constants).
    pub real_rooted: bool,
    /// Real roots counted with multiplicity.
    pub real_root_count: usize,
    /// Distinct real roots, increasing, with exact enclosures.
    pub roots: Vec<RootEnclosureJson>,
    /// The tested interval at this index.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interval: Option<String>,
    /// Every real root lies in the interval (decided exactly).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub real_roots_in_interval: Option<bool>,
    /// Real-rooted with every root in the interval (decided exactly).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub real_rooted_in_interval: Option<bool>,
}

/// Versioned report of [`root_enclosures`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEnclosureReport {
    /// `polytool.root-enclosures.v1`.
    pub schema: String,
    /// The request, retained for replay.
    pub request: RootEnclosureRequest,
    /// One item per polynomial.
    pub items: Vec<RootEnclosureItem>,
}

/// Exact rational root enclosures with multiplicities, and exact interval
/// tests, for a batch of polynomials.
///
/// Each verdict is exact for the supplied polynomial; it says nothing about
/// any other member of a family.
pub fn root_enclosures(request: &RootEnclosureRequest) -> Result<RootEnclosureReport, String> {
    if request.polynomials.len() > ROOT_ENCLOSURES_MAX_POLYNOMIALS {
        return Err(format!(
            "at most {ROOT_ENCLOSURES_MAX_POLYNOMIALS} polynomials are accepted"
        ));
    }
    let max_width = request
        .max_width
        .as_deref()
        .map(parse_rational_coeff)
        .transpose()?;
    if max_width.as_ref().is_some_and(|w| !w.is_positive()) {
        return Err("max_width must be positive".into());
    }
    let mut items = Vec::with_capacity(request.polynomials.len());
    for (position, text) in request.polynomials.iter().enumerate() {
        let index = request
            .first_index
            .checked_add(position)
            .ok_or("index overflow")?;
        let polynomial = text
            .iter()
            .map(|c| parse_rational_coeff(c))
            .collect::<Result<Vec<_>, _>>()?;
        if vec_poly::degree(&polynomial).is_some_and(|d| d > ROOT_ENCLOSURES_MAX_DEGREE) {
            return Err(format!(
                "polynomial {index} exceeds degree {ROOT_ENCLOSURES_MAX_DEGREE}"
            ));
        }
        let mut isolation = isolate_real_roots_rational_coeffs(&polynomial);
        if let Some(width) = &max_width {
            isolation.refine(width);
        }
        let interval = request.interval.as_ref().map(|i| i.at(index)).transpose()?;
        items.push(RootEnclosureItem {
            index,
            degree: isolation.degree,
            valuation: isolation.valuation,
            real_rooted: isolation.is_real_rooted(),
            real_root_count: isolation.real_root_count,
            roots: RootEnclosureJson::from_isolation(&isolation),
            interval: interval.as_ref().map(ToString::to_string),
            real_roots_in_interval: interval.as_ref().map(|i| isolation.all_real_roots_in(i)),
            real_rooted_in_interval: interval.as_ref().map(|i| isolation.real_rooted_in(i)),
        });
    }
    Ok(RootEnclosureReport {
        schema: ROOT_ENCLOSURES_SCHEMA.to_string(),
        request: request.clone(),
        items,
    })
}

/// Verdict of [`jensen_fixed_cap_check`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum JensenCapStatus {
    /// `J_N` has degree `N` and only real roots of one sign: the diagonal map
    /// sends every real-rooted polynomial of degree at most `N` to a
    /// real-rooted polynomial or zero (Schur--Malo composition theorem).
    PreservesUpToCap,
    /// `J_N = T[(1 + x)^N]` is not real-rooted: a witness that the map does not
    /// preserve real-rootedness in degree `N`.
    FailsAtCap,
    /// `J_N` is real-rooted but has roots of both signs or degree below `N`;
    /// no conclusion is drawn.
    Inconclusive,
}

/// Fixed-cap Jensen-polynomial test of a diagonal map `x^k -> lambda_k x^k`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct JensenFixedCapReport {
    /// Always `fixed_cap`: a statement for degrees up to `cap` only.
    pub evidence: String,
    /// The degree cap `N`.
    pub cap: usize,
    /// `lambda_0, ..., lambda_N` as rational strings.
    pub multipliers: Vec<String>,
    /// `J_N(x) = sum_k C(N, k) lambda_k x^k`, ascending rational strings.
    pub jensen_polynomial: Vec<String>,
    /// Exact root enclosures of `J_N`.
    pub roots: Vec<RootEnclosureJson>,
    /// Whether `J_N` is real-rooted.
    pub real_rooted: bool,
    /// The verdict.
    pub status: JensenCapStatus,
    /// What the verdict proves, and what it does not.
    pub statement: String,
}

/// Test the diagonal map `x^k -> lambda_k x^k` on polynomials of degree at
/// most `cap` with the Jensen polynomial `J_N = T[(1 + x)^N]`, `N = cap`.
///
/// `PreservesUpToCap` is a proof for degrees at most `cap` (Schur--Malo:
/// if `sum C(N,k) a_k x^k` is real-rooted and `sum C(N,k) b_k x^k` has degree
/// `N` and only real roots of one sign, then `sum C(N,k) a_k b_k x^k` is
/// real-rooted; lower degrees follow by Hurwitz's theorem).  It is not a proof
/// for larger degrees.  `FailsAtCap` exhibits `(1 + x)^N` as a real-rooted
/// polynomial with non-real-rooted image.
///
/// ```
/// use num_bigint::BigInt;
/// use num_rational::Ratio;
/// use polytool::recurrence::{jensen_fixed_cap_check, JensenCapStatus};
///
/// // P = (x Q)': lambda_k = k + 1 preserves real-rootedness at every cap ...
/// let up: Vec<_> = (1..=6).map(|k| Ratio::from_integer(BigInt::from(k))).collect();
/// assert_eq!(jensen_fixed_cap_check(&up, 5).unwrap().status, JensenCapStatus::PreservesUpToCap);
/// // ... but its inverse lambda_k = 1/(k + 1) fails at cap 2: J_2 = 1 + x + x^2/3.
/// let down: Vec<_> = (1..=3).map(|k| Ratio::new(BigInt::from(1), BigInt::from(k))).collect();
/// assert_eq!(jensen_fixed_cap_check(&down, 2).unwrap().status, JensenCapStatus::FailsAtCap);
/// ```
pub fn jensen_fixed_cap_check(
    multipliers: &[BigRational],
    cap: usize,
) -> Result<JensenFixedCapReport, String> {
    if cap > ROOT_ENCLOSURES_MAX_DEGREE {
        return Err(format!("cap exceeds {ROOT_ENCLOSURES_MAX_DEGREE}"));
    }
    if multipliers.len() <= cap {
        return Err(format!(
            "need lambda_0..lambda_{cap} ({} multipliers given)",
            multipliers.len()
        ));
    }
    let lambdas = &multipliers[..=cap];
    let jensen: Vec<Q> = lambdas
        .iter()
        .enumerate()
        .map(|(k, lambda)| Q::from_integer(crate::numeric::binomial_bigint(cap, k)) * lambda)
        .collect();
    let isolation = isolate_real_roots_rational_coeffs(&jensen);
    let one_signed = isolation.all_real_roots_in(&RealInterval::at_most(Q::zero()))
        || isolation.all_real_roots_in(&RealInterval::at_least(Q::zero()));
    let status = if !isolation.is_real_rooted() {
        JensenCapStatus::FailsAtCap
    } else if isolation.degree == Some(cap) && one_signed {
        JensenCapStatus::PreservesUpToCap
    } else {
        JensenCapStatus::Inconclusive
    };
    let statement = match status {
        JensenCapStatus::PreservesUpToCap => format!(
            "every real-rooted polynomial of degree at most {cap} has a real-rooted (or zero) \
             image; this is not a proof for degree {} or more",
            cap + 1
        ),
        JensenCapStatus::FailsAtCap => format!(
            "(1 + x)^{cap} is real-rooted and its image J_{cap} is not, so the map does not \
             preserve real-rootedness in degree {cap}"
        ),
        JensenCapStatus::Inconclusive => format!(
            "J_{cap} is real-rooted but has degree below {cap} or roots of both signs; no \
             conclusion is drawn"
        ),
    };
    Ok(JensenFixedCapReport {
        evidence: FIXED_CAP_EVIDENCE.to_string(),
        cap,
        multipliers: lambdas.iter().map(format_rational_coeff).collect(),
        jensen_polynomial: vec_poly::trim(&jensen)
            .iter()
            .map(format_rational_coeff)
            .collect(),
        roots: RootEnclosureJson::from_isolation(&isolation),
        real_rooted: isolation.is_real_rooted(),
        status,
        statement,
    })
}
