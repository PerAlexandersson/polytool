//! Optional root-domain information attached to transformed recurrence
//! results.
//!
//! Algebraic discovery ([`transform_recurrence`], transport, independent
//! search) never needs an interval or a real-rootedness proof.  When a
//! request carries a [`RootDomainRequest`] (field `root_domain`), the report
//! gains a [`RootDomainReport`] that records, separately,
//!
//! - **implications** ([`RootImplication`]): theorem adapters for the map,
//!   each with its direction, hypotheses (with their status), degree bounds,
//!   coefficient-sign hypotheses, endpoint conventions, zero cases, index
//!   domain and, when it fails, the converse;
//! - **claims** ([`RootDomainClaim`]) about the source rows `P_n` or the
//!   transformed rows, each labelled by its [`RootEvidence`]: a finite exact
//!   observation, a supplied assumption, a proposed interval, a completed
//!   recurrence induction, or a certified family invariant;
//! - a **recurrence theorem match** ([`RecurrenceRootTheorem`]) of the
//!   transported rule against Rolle's theorem for derivatives of products,
//!   with every interval attempt;
//! - per-row **observations** ([`RowRootObservation`]): exact rational root
//!   enclosures with multiplicities from
//!   [`crate::real_rootedness::isolate_real_roots_rational_coeffs`].
//!
//! # Theorem adapters
//!
//! An invertible map need not preserve real-rootedness in both directions,
//! so every implication states its direction.
//!
//! - `monomial_factor` (`multiply_power`, `divide_power`): `Q_n = x^r P_n`
//!   is real-rooted iff `P_n` is; only the root `0` changes multiplicity.
//! - `real_mobius_equivalence` (`reversal`, `affine`, `mobius`,
//!   `f_polynomial`): for `Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))`
//!   with real parameters, `a e - b c != 0` and `deg P_n <= d_n`, the roots of
//!   `Q_n` are `psi(rho) = (e rho - b)/(a - c rho)` for the roots `rho` of
//!   `P_n` (a root `rho = a/c` lowers `deg Q_n` instead), together with the
//!   pole `-e/c` of multiplicity `d_n - deg P_n`.  A real Möbius map is a
//!   bijection of the extended real line that maps non-real points to
//!   non-real points, so `P_n` is real-rooted iff `Q_n` is.  Intervals are
//!   mapped exactly (rays included) when the pole is not interior.
//! - `nonnegative_coefficient_window`: a nonzero `P_n` with nonnegative
//!   coefficients has every real root in `(-inf, 0]`, so every real root of
//!   `Q_n` lies in the image window.  For the `f`-polynomial
//!   `Q = (1 + x)^d P(x/(1 + x)) = sum_j p_j x^j (1 + x)^(d - j)` this is
//!   `[-1, 0]` (also visible from `(-1)^d Q(-1 - x) = sum_j p_j (1 + x)^j x^(d - j)`):
//!   `-1` records missing top degree and `0` the valuation.  The window is a
//!   *proposed interval*: it does not prove that `Q_n` splits (for
//!   `P = 1 + z + z^2` the window holds and `Q` has no real roots).
//! - `gamma_quarter_bound` (`gamma`): for `P = (1 + x)^d G(x/(1 + x)^2)` with
//!   `deg G <= d/2`, `P` is real-rooted iff every root of `G` is real and at
//!   most `1/4`; each root `y` contributes `x - y (1 + x)^2`, whose
//!   discriminant is `1 - 4 y`.  Real-rootedness of `G` alone is
//!   insufficient: `G = 1 - z` at `d = 2` gives `1 + x + x^2`.
//! - `hermite_kakeya_obreschkoff_symmetric_pair` (`symmetric_decomposition`):
//!   if `a_n` and `x b_n` are real-rooted and interlace, then
//!   `P_n = a_n + x b_n` is real-rooted.  Real-rooted components alone do not
//!   suffice (`a = (1 + x)^3`, `b = x` at `d = 3` give `1 + 3x + 4x^2 + x^3`),
//!   and the converse fails (`P = (1 + 2x)^2` at `d = 2` has `a = 1 + x + x^2`).
//! - `rolle_derivative_of_product` (recurrence theorem): if every root of `g`
//!   and of `Q_(m-1)` is real and lies in an interval `I`, then
//!   `s(m) D^j (g Q_(m-1))` is zero or real-rooted with every root in `I`.
//!   Matching the transported rule `Q_m = s(m) D^j(g Q_(m-1))` and checking
//!   `g` and the seeds exactly closes an induction over the whole index
//!   domain.  The canonical Eulerian rows A008292 under the `f`-polynomial map
//!   give `Q_(n+1) = (x(1 + x) Q_n)'`, `Q_1 = 1`, closed on `[-1, 0]`; the
//!   Möbius equivalence then makes every Eulerian polynomial real-rooted.
//!   The factor may also be `g(x) (c(m) + l(m) x)` with an index-dependent
//!   linear part, and candidate intervals may have endpoints affine in `n`;
//!   the step `I_(m-1) ⊆ I_m` and the root `-c(m)/l(m)` in `I_m` are then
//!   checked symbolically (exact integer-range sign checks from the first
//!   target), and attempts that cannot be certified are reported as failed.
//!   In gamma coordinates a completed induction for `G_n` is transferred to
//!   `P_n` only when every interval lies in `(-inf, 1/4]`.
//!
//! Index-dependent intervals map exactly through affine maps and monomial
//! factors; through a Möbius map with a pole (whose image endpoints are
//! rational in `n`) their hull over the index domain is mapped instead, and
//! the claim says so.  Fixed-cap diagonal multiplier tests (Jensen
//! polynomials) are in [`jensen_fixed_cap_check`]; their evidence is kept
//! separate from a proof for every degree.  Diagonal maps and derivative
//! antecedents (implicit transports) have their own adapters in
//! [`analyze_implicit_root_domain`].
//!
//! # Evidence
//!
//! Exact root isolation of a fixed polynomial proves statements about that
//! polynomial only.  A claim about every index needs a completed recurrence
//! induction or a transfer through a theorem whose conditions are certified
//! for every index; finitely many exact checks are never promoted to a proof.
//! Supplied assumptions stay visible as such (with exact checks on the rows),
//! and everything is replayed by [`verify_recurrence_transform`].

use super::*;
use crate::integer_range_sign::{
    decide_index_polynomial_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::real_rootedness::{
    interlaces_either_order_bigint_coeffs, isolate_real_roots_rational_coeffs, RealInterval,
    RealRootIsolation,
};
use crate::vec_poly;
use crate::Polynomial;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

pub(super) type Q = BigRational;

/// Largest `max_candidate_intervals`, `max_observed_rows`, number of
/// candidate intervals and number of assumptions accepted in a request.
pub const ROOT_DOMAIN_MAX_ITEMS: usize = 256;

/// Largest `max_observed_degree` accepted in a request.
pub const ROOT_DOMAIN_MAX_DEGREE: usize = 1024;

pub(super) const ZERO_CASES: &str =
    "the zero polynomial satisfies every root-domain property vacuously \
                          (P_n = 0 iff the transformed row is zero); nonzero constants are \
                          real-rooted with no roots; roots are counted with multiplicity";

// ---------------------------------------------------------------------------
// Contract types
// ---------------------------------------------------------------------------

/// How a root-domain claim was obtained.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum RootEvidence {
    /// Exact root isolation of finitely many listed rows; nothing about other indices.
    FiniteExactObservation,
    /// Supplied in the request (or derived from a supplied assumption); checked on rows only.
    SuppliedAssumption,
    /// A candidate window (for example from nonnegative coefficients), not proved.
    ProposedInterval,
    /// A closed induction over the transported rule: every index of the domain.
    CompletedRecurrenceInduction,
    /// Transferred to every index by a theorem whose conditions are certified.
    CertifiedFamilyInvariant,
}

/// Which rows a claim is about.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum RootSubject {
    /// The source rows `P_n`.
    Source,
    /// A component of the transformed rows (`Q_n`, or `a_n`, `b_n`).
    Transformed,
}

/// A root-domain property of a polynomial `F`.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum RootProperty {
    /// Every root of `F` is real.
    RealRooted,
    /// Every root of `F` is real and lies in the interval.
    RealRootedIn,
    /// Every real root of `F` lies in the interval (non-real roots allowed).
    RealRootsIn,
}

/// Direction of a root implication between source and transformed rows.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RootImplicationDirection {
    /// The source property implies the target property, not conversely.
    SourceToTarget,
    /// The target property implies the source property, not conversely.
    TargetToSource,
    /// The two properties are equivalent.
    Equivalent,
}

/// Status of one hypothesis of an implication.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum HypothesisStatus {
    /// Proved for every index of the domain (from the source rule).
    CertifiedAllIndices,
    /// Verified exactly on every supplied (or observed) row only.
    CheckedSuppliedRows,
    /// Violated on at least one row (see the detail).
    FailedOnRows,
    /// A fixed polynomial checked exactly: a proof for that polynomial.
    ExactFixedPolynomial,
    /// Supplied in the request and not proved.
    Assumed,
    /// Not checked.
    NotChecked,
    /// Shown false exactly (for example a fixed polynomial with a non-real
    /// root); the implication does not apply.
    Violated,
}

/// One endpoint `value + per_index * n` of a [`RootInterval`].
///
/// # Exact quadratic endpoints
///
/// An index-independent endpoint may also carry its exact value
/// `a + b sqrt(d)` ([`RootEndpoint::exact`], an irrational
/// [`ExactQuadraticEndpoint`]; omitted for a rational endpoint, whose exact
/// value is `value` itself).  The rational `value` is then an **outer bound**
/// of the exact one, strictly on the outer side: below it for a lower
/// endpoint, above it for an upper endpoint.  `closed` refers to the exact
/// endpoint.  Since the rational bound lies strictly outside, the rational
/// interval with the same closedness contains the exact interval, so a
/// consumer that reads only `value` (and `per_index`, `closed`) still works
/// with a correct, slightly larger interval.  The invariant is checked by
/// [`RootInterval::check_exact_endpoints`].  Exact forms are produced by
/// the library (the basis-map lifts, whose target intervals have endpoints
/// in a real quadratic field) and are not accepted in request intervals,
/// whose checks decide the rational interval only.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEndpoint {
    /// Rational constant term (integer or `num/den` string); with an exact
    /// form, a rational outer bound of it.
    pub value: String,
    /// Rational coefficient of the original index `n`; omitted means `0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub per_index: Option<String>,
    /// Whether the endpoint belongs to the interval (default `true`); with an
    /// exact form, whether the exact endpoint does.
    #[serde(default = "default_true")]
    pub closed: bool,
    /// The exact irrational value `a + b sqrt(d)` of an index-independent
    /// endpoint, when known; omitted for a rational endpoint (and in records
    /// written before the field existed).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exact: Option<ExactQuadraticEndpoint>,
}

fn default_true() -> bool {
    true
}

impl RootEndpoint {
    /// The rational endpoint `value` (index-independent, no exact form).
    pub fn rational(value: &Q, closed: bool) -> Self {
        Self {
            value: format_rational_coeff(value),
            per_index: None,
            closed,
            exact: None,
        }
    }

    /// The endpoint as plain ASCII text: the exact value when it is
    /// recorded (`-3-2*sqrt(2)`, see [`ExactQuadraticEndpoint::to_plain_text`]),
    /// else `value`, followed by ` + (per_index) n` for an index-dependent
    /// endpoint.
    pub fn to_exact_text(&self) -> String {
        let constant = match &self.exact {
            Some(exact) => exact.to_plain_text(),
            None => self.value.clone(),
        };
        match self.per_index.as_deref() {
            None | Some("0") => constant,
            Some(p) => format!("{constant} + ({p}) n"),
        }
    }
}

/// The exact value `a + b sqrt(d)` of an irrational endpoint of a
/// [`RootInterval`] (see [`RootEndpoint`], "Exact quadratic endpoints").
///
/// The three parts are rational strings (integer or `num/den`), `d > 0` is
/// not the square of a rational and `b != 0`.  The JSON form is
/// `{"rational_part": "-3", "sqrt_coefficient": "-2", "radicand": "2"}`
/// for `-3 - 2 sqrt(2)`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ExactQuadraticEndpoint {
    /// The rational part `a`.
    pub rational_part: String,
    /// The coefficient `b` of `sqrt(d)`.
    pub sqrt_coefficient: String,
    /// The radicand `d`.
    pub radicand: String,
}

impl ExactQuadraticEndpoint {
    /// The exact form of `value`, or `None` when `value` is rational.
    pub fn from_quadratic_number(value: &QuadraticNumber) -> Option<Self> {
        (!value.is_rational()).then(|| Self {
            rational_part: format_rational_coeff(value.rational_part()),
            sqrt_coefficient: format_rational_coeff(value.sqrt_coefficient()),
            radicand: format_rational_coeff(value.radicand()),
        })
    }

    /// The exact number.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] when a part is not a rational, the
    /// radicand is negative, or the number is rational (a rational endpoint
    /// carries no exact form).
    pub fn to_quadratic_number(&self) -> Result<QuadraticNumber, RecurrenceError> {
        let number = QuadraticNumber::new(
            parse_rational_coeff(&self.rational_part)?,
            parse_rational_coeff(&self.sqrt_coefficient)?,
            parse_rational_coeff(&self.radicand)?,
        )?;
        if number.is_rational() {
            return Err(RecurrenceError::InvalidInput(format!(
                "the exact endpoint {} is rational (record it as the value)",
                number
            )));
        }
        Ok(number)
    }

    /// Plain ASCII text without spaces: `-3-2*sqrt(2)`, `1+sqrt(2)`,
    /// `-sqrt(2)`, `1/2+3/2*sqrt(5)`.  The parts are printed as recorded.
    pub fn to_plain_text(&self) -> String {
        let a = self.rational_part.trim();
        let b = self.sqrt_coefficient.trim();
        let root = format!("sqrt({})", self.radicand.trim());
        let (negative, magnitude) = match b.strip_prefix('-') {
            Some(rest) => (true, rest),
            None => (false, b.strip_prefix('+').unwrap_or(b)),
        };
        let term = if magnitude == "1" {
            root
        } else {
            format!("{magnitude}*{root}")
        };
        let sign = if negative { "-" } else { "+" };
        if a == "0" {
            if negative {
                format!("-{term}")
            } else {
                term
            }
        } else {
            format!("{a}{sign}{term}")
        }
    }
}

/// The interval as plain ASCII text with exact endpoints where they are
/// recorded (`[-3-2*sqrt(2), -3+2*sqrt(2)]`) and the rational values
/// elsewhere; the [`fmt::Display`] of [`RootInterval`] prints the rational
/// values only (the canonical key of the interval).
pub fn format_root_interval_exact(interval: &RootInterval) -> String {
    let mut text = String::new();
    match &interval.lower {
        None => text.push_str("(-inf"),
        Some(e) => {
            text.push(if e.closed { '[' } else { '(' });
            text.push_str(&e.to_exact_text());
        }
    }
    text.push_str(", ");
    match &interval.upper {
        None => text.push_str("+inf)"),
        Some(e) => {
            text.push_str(&e.to_exact_text());
            text.push(if e.closed { ']' } else { ')' });
        }
    }
    text
}

/// An interval of the real line with rational, index-dependent or infinite
/// endpoints; an omitted endpoint is infinite (a ray, or the whole line).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootInterval {
    /// Lower endpoint; omitted means `-infinity`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower: Option<RootEndpoint>,
    /// Upper endpoint; omitted means `+infinity`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub upper: Option<RootEndpoint>,
}

impl RootInterval {
    /// The index-independent interval with the endpoints of `interval`.
    pub fn from_real(interval: &RealInterval) -> Self {
        let end = |value: &Option<Q>, closed: bool| {
            value.as_ref().map(|v| RootEndpoint::rational(v, closed))
        };
        Self {
            lower: end(&interval.lower, interval.lower_closed),
            upper: end(&interval.upper, interval.upper_closed),
        }
    }

    /// Whether no endpoint depends on the index.
    pub fn is_index_independent(&self) -> bool {
        [&self.lower, &self.upper].into_iter().all(|end| {
            end.as_ref()
                .and_then(|e| e.per_index.as_ref())
                .is_none_or(|p| parse_rational_coeff(p).is_ok_and(|q| q.is_zero()))
        })
    }

    /// The concrete interval at the original index `index`.
    pub fn at(&self, index: usize) -> Result<RealInterval, RecurrenceError> {
        let n = Q::from_integer(BigInt::from(index));
        let value = |end: &RootEndpoint| -> Result<Q, RecurrenceError> {
            let constant = parse_rational_coeff(&end.value)?;
            let slope = match &end.per_index {
                Some(p) => parse_rational_coeff(p)?,
                None => Q::zero(),
            };
            Ok(constant + slope * &n)
        };
        let lower = self.lower.as_ref().map(value).transpose()?;
        let upper = self.upper.as_ref().map(value).transpose()?;
        Ok(RealInterval {
            lower_closed: lower.is_some() && self.lower.as_ref().is_some_and(|e| e.closed),
            lower,
            upper_closed: upper.is_some() && self.upper.as_ref().is_some_and(|e| e.closed),
            upper,
        })
    }

    fn validate(&self) -> Result<(), RecurrenceError> {
        self.at(0)?;
        self.check_request_interval()
    }

    /// Requests (candidate intervals, assumptions, enclosure intervals) take
    /// rational endpoints only: their checks decide the rational interval,
    /// so an exact form there would be claimed without proof.
    pub(super) fn check_request_interval(&self) -> Result<(), RecurrenceError> {
        if self.has_exact_endpoints() {
            return Err(RecurrenceError::InvalidInput(
                "request intervals take rational endpoints only (no exact form)".into(),
            ));
        }
        Ok(())
    }

    /// Check the exact forms of the endpoints (see [`RootEndpoint`], "Exact
    /// quadratic endpoints"): each is a valid irrational quadratic number on
    /// an index-independent endpoint, and the rational `value` lies strictly
    /// below it at the lower end and strictly above it at the upper end.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] naming the violated condition.
    pub fn check_exact_endpoints(&self) -> Result<(), RecurrenceError> {
        for (end, is_lower) in [(&self.lower, true), (&self.upper, false)] {
            let Some(end) = end else { continue };
            let Some(exact) = &end.exact else { continue };
            let side = if is_lower { "lower" } else { "upper" };
            let slope = match &end.per_index {
                Some(p) => parse_rational_coeff(p)?,
                None => Q::zero(),
            };
            if !slope.is_zero() {
                return Err(RecurrenceError::InvalidInput(format!(
                    "the {side} endpoint has an exact form and depends on the index"
                )));
            }
            let number = exact.to_quadratic_number()?;
            let bound = parse_rational_coeff(&end.value)?;
            let outer = match number.cmp_rational(&bound) {
                std::cmp::Ordering::Greater => is_lower,
                std::cmp::Ordering::Less => !is_lower,
                std::cmp::Ordering::Equal => false,
            };
            if !outer {
                return Err(RecurrenceError::InvalidInput(format!(
                    "the rational value {} of the {side} endpoint is not a strict outer bound \
                     of its exact value {number}",
                    end.value
                )));
            }
        }
        Ok(())
    }

    /// Whether some endpoint records an exact form.
    pub fn has_exact_endpoints(&self) -> bool {
        [&self.lower, &self.upper]
            .into_iter()
            .any(|end| end.as_ref().is_some_and(|e| e.exact.is_some()))
    }

    /// The interval without exact forms (the rational outer interval every
    /// consumer of `value` uses).
    pub fn without_exact_endpoints(&self) -> Self {
        let strip =
            |end: &Option<RootEndpoint>| end.clone().map(|e| RootEndpoint { exact: None, ..e });
        Self {
            lower: strip(&self.lower),
            upper: strip(&self.upper),
        }
    }

    /// Whether a recorded interval `self` agrees with its recomputation
    /// `recomputed`, allowing a record written before exact endpoints were
    /// recorded: equal, or equal after dropping each exact form that
    /// `self` lacks, where such a legacy endpoint may also be closed when
    /// the recomputed exact endpoint is open (the rational outer bound lies
    /// strictly outside, so the closed legacy endpoint is a correct, weaker
    /// statement).  A recorded exact form must equal the recomputed one.
    pub fn agrees_with_recomputed(&self, recomputed: &Self) -> bool {
        let end_agrees =
            |recorded: &Option<RootEndpoint>, fresh: &Option<RootEndpoint>| match (recorded, fresh)
            {
                (None, None) => true,
                (Some(r), Some(f)) => {
                    r.value == f.value
                        && r.per_index == f.per_index
                        && match (&r.exact, &f.exact) {
                            (None, Some(_)) => r.closed || !f.closed,
                            (recorded_exact, fresh_exact) => {
                                recorded_exact == fresh_exact && r.closed == f.closed
                            }
                        }
                }
                _ => false,
            };
        end_agrees(&self.lower, &recomputed.lower) && end_agrees(&self.upper, &recomputed.upper)
    }
}

impl fmt::Display for RootInterval {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = |end: &RootEndpoint| match end.per_index.as_deref() {
            None | Some("0") => end.value.clone(),
            Some(p) => format!("{} + ({p}) n", end.value),
        };
        match &self.lower {
            None => write!(f, "(-inf")?,
            Some(e) => write!(f, "{}{}", if e.closed { '[' } else { '(' }, text(e))?,
        }
        match &self.upper {
            None => write!(f, ", +inf)"),
            Some(e) => write!(f, ", {}{}", text(e), if e.closed { ']' } else { ')' }),
        }
    }
}

/// A supplied root-domain assumption about the source or transformed rows.
///
/// It is checked exactly on the observed rows and transferred through the
/// map's implication when the direction allows; it is never proved.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootAssumption {
    /// `source` (`P_n`) or `transformed`.
    pub subject: RootSubject,
    /// Component of a transformed state (0 for scalar maps and for the source).
    #[serde(default)]
    pub component: usize,
    /// `real_rooted`, `real_rooted_in` or `real_roots_in`.
    pub property: RootProperty,
    /// The interval; required for `real_rooted_in` and `real_roots_in` only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interval: Option<RootInterval>,
    /// First original index the assumption covers; omitted means the first row.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub first_index: Option<usize>,
    /// Last original index it covers; omitted means unbounded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_index: Option<usize>,
    /// Free-text provenance of the assumption.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

/// Opt-in root-domain analysis of a transformation request
/// (`RecurrenceTransformRequest::root_domain`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct RootDomainRequest {
    /// Candidate intervals for the transformed rows, tried (in order, before
    /// the default window and the root hull) in the recurrence theorem match.
    pub candidate_intervals: Vec<RootInterval>,
    /// Supplied assumptions, reported as such.
    pub assumptions: Vec<RootAssumption>,
    /// Propose the nonnegative-coefficient window when every source row has
    /// nonnegative coefficients (default `true`).
    pub default_window: bool,
    /// Largest number of interval attempts in the recurrence theorem match.
    pub max_candidate_intervals: usize,
    /// Largest number of rows whose roots are isolated (from the first row).
    pub max_observed_rows: usize,
    /// Rows of larger degree (source or transformed) are not isolated.
    pub max_observed_degree: usize,
}

impl Default for RootDomainRequest {
    fn default() -> Self {
        Self {
            candidate_intervals: Vec::new(),
            assumptions: Vec::new(),
            default_window: true,
            max_candidate_intervals: 8,
            max_observed_rows: 32,
            max_observed_degree: 128,
        }
    }
}

impl RootDomainRequest {
    /// Reject malformed intervals and limits above the library ceilings
    /// [`ROOT_DOMAIN_MAX_ITEMS`] and [`ROOT_DOMAIN_MAX_DEGREE`].
    pub fn validate(&self) -> Result<(), RecurrenceError> {
        for (name, value) in [
            ("max_candidate_intervals", self.max_candidate_intervals),
            ("max_observed_rows", self.max_observed_rows),
            ("candidate_intervals", self.candidate_intervals.len()),
            ("assumptions", self.assumptions.len()),
        ] {
            if value > ROOT_DOMAIN_MAX_ITEMS {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "root_domain.{name} = {value} exceeds {ROOT_DOMAIN_MAX_ITEMS}"
                )));
            }
        }
        if self.max_observed_degree > ROOT_DOMAIN_MAX_DEGREE {
            return Err(RecurrenceError::LimitExceeded(format!(
                "root_domain.max_observed_degree exceeds {ROOT_DOMAIN_MAX_DEGREE}"
            )));
        }
        for interval in &self.candidate_intervals {
            interval.validate()?;
        }
        for (k, assumption) in self.assumptions.iter().enumerate() {
            match (&assumption.property, &assumption.interval) {
                (RootProperty::RealRooted, Some(_)) => {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "assumption {k}: real_rooted takes no interval (use real_rooted_in)"
                    )))
                }
                (RootProperty::RealRootedIn | RootProperty::RealRootsIn, None) => {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "assumption {k}: the property needs an interval"
                    )))
                }
                (_, Some(interval)) => interval.validate()?,
                _ => {}
            }
            if assumption.subject == RootSubject::Source && assumption.component != 0 {
                return Err(RecurrenceError::InvalidInput(format!(
                    "assumption {k}: the source has one component"
                )));
            }
            if let (Some(first), Some(last)) = (assumption.first_index, assumption.last_index) {
                if last < first {
                    return Err(RecurrenceError::InvalidInput(format!(
                        "assumption {k}: last_index precedes first_index"
                    )));
                }
            }
        }
        Ok(())
    }
}

/// One hypothesis of an implication and its status.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootHypothesis {
    /// The hypothesis.
    pub condition: String,
    /// Its status.
    pub status: HypothesisStatus,
    /// Why this status was assigned.
    pub detail: String,
}

/// A theorem adapter: a root implication between source and transformed rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootImplication {
    /// Stable theorem identifier (for example `real_mobius_equivalence`).
    pub theorem: String,
    /// Direction of the implication.
    pub direction: RootImplicationDirection,
    /// The property of the source rows `P_n`.
    pub source_property: String,
    /// The property of the transformed rows.
    pub target_property: String,
    /// The precise statement, including how roots and intervals correspond.
    pub statement: String,
    /// Hypotheses with their status for this request.
    pub hypotheses: Vec<RootHypothesis>,
    /// Degree bounds used by the map.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree_bounds: Option<String>,
    /// Coefficient-sign hypotheses.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coefficient_signs: Option<String>,
    /// Endpoint conventions (missing degree, valuation, poles).
    pub endpoint_conventions: Vec<String>,
    /// Treatment of zero rows and constants.
    pub zero_cases: String,
    /// Index domain of the map.
    pub index_domain: String,
    /// Why the converse fails, for one-directional implications.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub converse: Option<String>,
}

/// A root-domain claim with its evidence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootDomainClaim {
    /// `source` or `transformed`.
    pub subject: RootSubject,
    /// Component of a transformed state (0 otherwise).
    pub component: usize,
    /// The property.
    pub property: RootProperty,
    /// The interval of `real_rooted_in` / `real_roots_in`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interval: Option<RootInterval>,
    /// First original index covered.
    pub first_index: usize,
    /// Last original index covered; omitted means unbounded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_index: Option<usize>,
    /// How the claim was obtained.
    pub evidence: RootEvidence,
    /// Unproved hypotheses the claim depends on; empty means none.
    pub conditional_on: Vec<String>,
    /// How the claim follows.
    pub derivation: String,
}

/// One candidate interval tried in the recurrence theorem match.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IntervalAttempt {
    /// The candidate interval.
    pub interval: RootInterval,
    /// `supplied`, `default_window` or `root_hull`.
    pub origin: String,
    /// Whether the induction closed on this interval.
    pub success: bool,
    /// The exact checks performed, or why the attempt failed.
    pub detail: String,
}

/// Result of matching the transported rule against a recurrence theorem.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceRootTheorem {
    /// `rolle_derivative_of_product`.
    pub theorem: String,
    /// The theorem statement.
    pub statement: String,
    /// Whether the rule has the form `Q_m = s(m) D^j (g Q_(m-1))`.
    pub matched: bool,
    /// Why the rule does not match.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    /// The matched form, in text.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub form: Option<String>,
    /// The factor `g(x)`, ascending rational strings.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub factor: Option<Vec<String>>,
    /// The derivative order `j`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derivative_order: Option<usize>,
    /// The scalar `s(m)` in the target index, ascending rational strings in `m`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub index_scalar: Option<Vec<String>>,
    /// The index-dependent linear part of the factor, when the rule is
    /// `Q_m = D^j (g(x) (c(m) + l(m) x) Q_(m-1))` (then `index_scalar` is `1`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub index_linear_factor: Option<IndexLinearFactorJson>,
    /// Exact degree formula of the generated rows, when `s(m)` provably never
    /// vanishes and `deg g >= j`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree_formula: Option<String>,
    /// Interval attempts, in order (bounded by `max_candidate_intervals`).
    pub attempts: Vec<IntervalAttempt>,
}

/// The index-dependent linear part `c(m) + l(m) x` of a Rolle factor.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct IndexLinearFactorJson {
    /// `c(m)`, ascending rational strings in the target index `m`.
    pub constant: Vec<String>,
    /// `l(m)`, ascending rational strings in the target index `m`.
    pub slope: Vec<String>,
}

/// An exact enclosure of one distinct real root.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEnclosureJson {
    /// Lower endpoint (rational string).
    pub lower: String,
    /// Upper endpoint; equal to `lower` iff the root is this rational, otherwise
    /// the root lies in the open interval and neither endpoint is a root.
    pub upper: String,
    /// Multiplicity of the root.
    pub multiplicity: usize,
}

impl RootEnclosureJson {
    /// Encode every distinct root of an isolation.
    pub fn from_isolation(isolation: &RealRootIsolation) -> Vec<Self> {
        isolation
            .roots
            .iter()
            .map(|root| Self {
                lower: format_rational_coeff(&root.lower),
                upper: format_rational_coeff(&root.upper),
                multiplicity: root.multiplicity,
            })
            .collect()
    }
}

/// Exact root isolation of one observed row.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RowRootObservation {
    /// `source` or `transformed`.
    pub subject: RootSubject,
    /// Component of a transformed state (0 otherwise).
    pub component: usize,
    /// Original index of the row.
    pub index: usize,
    /// Degree; omitted for the zero polynomial.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree: Option<usize>,
    /// Multiplicity of the root 0.
    pub valuation: usize,
    /// Whether every root is real (true for zero and constants).
    pub real_rooted: bool,
    /// Number of real roots with multiplicity.
    pub real_root_count: usize,
    /// Distinct real roots in increasing order.
    pub roots: Vec<RootEnclosureJson>,
}

/// Exact check of a supplied assumption on the observed rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct AssumptionCheck {
    /// Position of the assumption in the request.
    pub assumption: usize,
    /// `consistent_with_rows`, `refuted_on_rows` or `not_checked`.
    pub status: String,
    /// Observed indices at which the assumption was checked.
    pub checked_indices: Vec<usize>,
    /// First index at which it fails, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub refuted_at: Option<usize>,
}

/// Root-domain section of a transformation report.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootDomainReport {
    /// Conventions for properties, intervals, zero rows and evidence.
    pub conventions: Vec<String>,
    /// Theorem adapters for the map, with hypothesis statuses.
    pub implications: Vec<RootImplication>,
    /// Claims, each with its evidence.
    pub claims: Vec<RootDomainClaim>,
    /// The recurrence theorem match of the transported scalar rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recurrence_theorem: Option<RecurrenceRootTheorem>,
    /// Exact checks of the supplied assumptions.
    pub assumption_checks: Vec<AssumptionCheck>,
    /// Exact root isolation of the observed rows.
    pub observations: Vec<RowRootObservation>,
    /// Limits reached and transfers not performed.
    pub notes: Vec<String>,
}

// ---------------------------------------------------------------------------
// Recurrence theorem: Rolle's theorem for derivatives of products
// ---------------------------------------------------------------------------

/// Trimmed rows of a bivariate polynomial (`[n power][x power]`, zero is empty).
pub(super) fn bivar_rows(p: &BivarPoly) -> Vec<Vec<Q>> {
    let mut rows: Vec<Vec<Q>> = p
        .coeffs
        .iter()
        .map(|row| {
            let mut row = row.clone();
            vec_poly::strip_trailing_zeros(&mut row);
            row
        })
        .collect();
    while rows.last().is_some_and(Vec::is_empty) {
        rows.pop();
    }
    rows
}

pub(super) fn binomial(n: usize, k: usize) -> Q {
    Q::from_integer(crate::numeric::binomial_bigint(n, k))
}

/// Whether the polynomial `s(m)` provably never vanishes for `m >= first_target`:
/// one-signed coefficients, and a nonzero constant term or `first_target >= 1`.
pub(super) fn scalar_nonvanishing(scalar: &[Q], first_target: usize) -> bool {
    let Some(constant) = scalar.first() else {
        return false;
    };
    let one_signed =
        scalar.iter().all(|c| !c.is_negative()) || scalar.iter().all(|c| !c.is_positive());
    one_signed && (!constant.is_zero() || first_target >= 1)
}

/// The sign of a polynomial `p(m)` (ascending in `m`) at every integer
/// `m >= start`, decided exactly by [`crate::decide_index_polynomial_sign`];
/// `None` when `p` vanishes or takes both signs at integers `m >= start`, or
/// when the checker's budget is exceeded (which proves nothing).
pub(crate) fn index_polynomial_sign_from(p: &[Q], start: usize) -> Option<std::cmp::Ordering> {
    let range = IntegerRange::starting_at(start);
    match decide_index_polynomial_sign(p, &range, IndexSignCondition::Positive) {
        IndexSignDecision::Holds { .. } => Some(std::cmp::Ordering::Greater),
        // Negative at the first index: the only way p can be negative throughout.
        IndexSignDecision::FailsAt {
            index,
            value: Some(value),
        } if index == range.first && value.is_negative() => {
            decide_index_polynomial_sign(p, &range, IndexSignCondition::Negative)
                .holds()
                .then_some(std::cmp::Ordering::Less)
        }
        _ => None,
    }
}

/// Whether `p(m) >= 0` (or `> 0` when `strict`) at every integer `m >= start`,
/// decided exactly by [`crate::decide_index_polynomial_sign`] (the zero
/// polynomial is nonnegative; `false` also when the checker's budget is
/// exceeded).
pub(super) fn index_polynomial_nonnegative_from(p: &[Q], start: usize, strict: bool) -> bool {
    let condition = if strict {
        IndexSignCondition::Positive
    } else {
        IndexSignCondition::Nonnegative
    };
    decide_index_polynomial_sign(p, &IntegerRange::starting_at(start), condition).holds()
}

/// Whether `f(m, x)` is a nonzero polynomial in `x` at every integer
/// `m >= start`: `f(m, .)` vanishes exactly at the common integer roots of
/// its `x`-coefficients, that is at the integer roots of their gcd, and the
/// gcd is decided nonzero on the domain by
/// [`crate::decide_index_polynomial_sign`] (`false` also over the checker's
/// budget).
pub(crate) fn bivariate_leader_nonvanishing(f: &BivarPoly, start: usize) -> bool {
    let rows = bivar_rows(f);
    let width = rows.iter().map(Vec::len).max().unwrap_or(0);
    let mut gcd = Polynomial::<Q>::zero();
    for j in 0..width {
        let column: Vec<Q> = rows
            .iter()
            .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
            .collect();
        gcd = gcd.gcd(&Polynomial::new(column));
    }
    !gcd.is_zero()
        && decide_index_polynomial_sign(
            gcd.coeffs(),
            &IntegerRange::starting_at(start),
            IndexSignCondition::Nonzero,
        )
        .holds()
}

pub(super) fn poly_text(p: &[Q], variable: &str) -> String {
    let mut parts = Vec::new();
    for (k, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let coefficient = format_rational_coeff(c);
        parts.push(match k {
            0 => coefficient,
            1 => format!("({coefficient}) {variable}"),
            _ => format!("({coefficient}) {variable}^{k}"),
        });
    }
    if parts.is_empty() {
        "0".into()
    } else {
        parts.join(" + ")
    }
}

pub(super) fn encode(row: &[Q]) -> Vec<String> {
    vec_poly::trim(row)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

// ---------------------------------------------------------------------------
// Report assembly
// ---------------------------------------------------------------------------

/// The transported scalar rule handed to the analysis.
pub(crate) struct TransportedInput<'a> {
    pub(crate) rule: &'a TransportedRecurrence,
    /// First index of the source seeds.
    pub(crate) first_index: usize,
    /// First target index of the rule (`first_index + seeds`).
    pub(crate) valid_from: usize,
}

/// Everything the analysis reads from a transformation report.
pub(crate) struct RootDomainInputs<'a> {
    pub(crate) request: &'a RootDomainRequest,
    pub(crate) transformation: &'a RecurrenceTransformation,
    pub(crate) origin: usize,
    pub(crate) rows: &'a [Vec<Q>],
    pub(crate) states: &'a [Vec<Vec<Q>>],
    pub(crate) conditions: &'a [TransformConditionStatus],
    pub(crate) transported: Option<TransportedInput<'a>>,
}

#[cfg(test)]
mod tests;

mod geometry;
#[cfg(test)]
use geometry::mobius_interval_image;
pub use geometry::root_implication;
pub(super) use geometry::{condition_hypotheses, uncertified};
use geometry::{conventions, Geometry};
mod checking;
pub(crate) use checking::analyze_root_domain;
pub(super) use checking::{check_assumption, observation_claims, observe_rows};
