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
//!   [`crate::real_rootedness::isolate_real_roots_exact`].
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
//!
//! Fixed-cap diagonal multiplier tests (Jensen polynomials) are in
//! [`jensen_fixed_cap_check`]; their evidence is kept separate from a proof
//! for every degree.
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
use crate::real_rootedness::{
    interlaces_either_order_bigint_coeffs, isolate_real_roots_exact, RealInterval,
    RealRootIsolation,
};
use crate::vec_poly;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

type Q = BigRational;

/// Largest `max_candidate_intervals`, `max_observed_rows`, number of
/// candidate intervals and number of assumptions accepted in a request.
pub const ROOT_DOMAIN_MAX_ITEMS: usize = 256;

/// Largest `max_observed_degree` accepted in a request.
pub const ROOT_DOMAIN_MAX_DEGREE: usize = 1024;

const ZERO_CASES: &str = "the zero polynomial satisfies every root-domain property vacuously \
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
}

/// One endpoint `value + per_index * n` of a [`RootInterval`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RootEndpoint {
    /// Rational constant term (integer or `num/den` string).
    pub value: String,
    /// Rational coefficient of the original index `n`; omitted means `0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub per_index: Option<String>,
    /// Whether the endpoint belongs to the interval (default `true`).
    #[serde(default = "default_true")]
    pub closed: bool,
}

fn default_true() -> bool {
    true
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
            value.as_ref().map(|v| RootEndpoint {
                value: format_rational_coeff(v),
                per_index: None,
                closed,
            })
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
    pub fn at(&self, index: usize) -> Result<RealInterval, String> {
        let n = Q::from_integer(BigInt::from(index));
        let value = |end: &RootEndpoint| -> Result<Q, String> {
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

    fn validate(&self) -> Result<(), String> {
        self.at(0).map(|_| ())
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
    pub fn validate(&self) -> Result<(), String> {
        for (name, value) in [
            ("max_candidate_intervals", self.max_candidate_intervals),
            ("max_observed_rows", self.max_observed_rows),
            ("candidate_intervals", self.candidate_intervals.len()),
            ("assumptions", self.assumptions.len()),
        ] {
            if value > ROOT_DOMAIN_MAX_ITEMS {
                return Err(format!(
                    "root_domain.{name} = {value} exceeds {ROOT_DOMAIN_MAX_ITEMS}"
                ));
            }
        }
        if self.max_observed_degree > ROOT_DOMAIN_MAX_DEGREE {
            return Err(format!(
                "root_domain.max_observed_degree exceeds {ROOT_DOMAIN_MAX_DEGREE}"
            ));
        }
        for interval in &self.candidate_intervals {
            interval.validate()?;
        }
        for (k, assumption) in self.assumptions.iter().enumerate() {
            match (&assumption.property, &assumption.interval) {
                (RootProperty::RealRooted, Some(_)) => {
                    return Err(format!(
                        "assumption {k}: real_rooted takes no interval (use real_rooted_in)"
                    ))
                }
                (RootProperty::RealRootedIn | RootProperty::RealRootsIn, None) => {
                    return Err(format!("assumption {k}: the property needs an interval"))
                }
                (_, Some(interval)) => interval.validate()?,
                _ => {}
            }
            if assumption.subject == RootSubject::Source && assumption.component != 0 {
                return Err(format!("assumption {k}: the source has one component"));
            }
            if let (Some(first), Some(last)) = (assumption.first_index, assumption.last_index) {
                if last < first {
                    return Err(format!("assumption {k}: last_index precedes first_index"));
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
    /// Exact degree formula of the generated rows, when `s(m)` provably never
    /// vanishes and `deg g >= j`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub degree_formula: Option<String>,
    /// Interval attempts, in order (bounded by `max_candidate_intervals`).
    pub attempts: Vec<IntervalAttempt>,
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

fn conventions() -> Vec<String> {
    [
        "rows are ascending coefficient vectors; P_n are the source rows and the transformed \
         rows are named by the transformation's components; n is the original index",
        "real_rooted: every root is real; real_rooted_in I: every root is real and lies in I; \
         real_roots_in I: every real root lies in I (non-real roots allowed)",
        ZERO_CASES,
        "intervals have rational, infinite (omitted) or index-dependent endpoints \
         value + per_index * n, each closed or open",
        "evidence: finite_exact_observation (exact on the listed rows only), \
         supplied_assumption (from the request, checked on rows only), proposed_interval \
         (candidate window, not proved), completed_recurrence_induction (closed induction over \
         the transported rule for every index of the domain), certified_family_invariant \
         (transferred to every index by a map theorem whose conditions are certified)",
        "exact root isolation proves statements about the isolated polynomials only, never \
         about the family",
    ]
    .into_iter()
    .map(String::from)
    .collect()
}

// ---------------------------------------------------------------------------
// Map geometry: exact interval images
// ---------------------------------------------------------------------------

/// The empty interval `(0, 0)`.
fn empty_interval() -> RealInterval {
    RealInterval {
        lower: Some(Q::zero()),
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: false,
    }
}

/// Image of `interval` under `y -> (p y + q)/(r y + s)` (`p s - q r != 0`),
/// or `None` when the pole lies in the interior (the image is not an interval).
/// The pole itself maps to infinity and is dropped.
pub(crate) fn mobius_interval_image(
    interval: &RealInterval,
    p: &Q,
    q: &Q,
    r: &Q,
    s: &Q,
) -> Option<RealInterval> {
    if interval.is_empty() {
        return Some(interval.clone());
    }
    if !r.is_zero() {
        let pole = -(s / r);
        let is_lower = interval.lower.as_ref() == Some(&pole);
        let is_upper = interval.upper.as_ref() == Some(&pole);
        if is_lower && is_upper {
            return Some(empty_interval());
        }
        if !is_lower && !is_upper && interval.contains(&pole) {
            return None;
        }
    }
    let image = |end: &Option<Q>, closed: bool| -> (Option<Q>, bool) {
        match end {
            None if r.is_zero() => (None, false),
            None => (Some(p / r), false),
            Some(y) => {
                let denominator = r * y + s;
                if denominator.is_zero() {
                    (None, false)
                } else {
                    (Some((p * y + q) / denominator), closed)
                }
            }
        }
    };
    let low = image(&interval.lower, interval.lower_closed);
    let high = image(&interval.upper, interval.upper_closed);
    let increasing = (p * s - q * r).is_positive();
    let (lower, upper) = if increasing { (low, high) } else { (high, low) };
    Some(RealInterval {
        lower_closed: lower.0.is_some() && lower.1,
        lower: lower.0,
        upper_closed: upper.0.is_some() && upper.1,
        upper: upper.0,
    })
}

// One short-lived value per analysis; boxing the Möbius parameters buys nothing.
#[allow(clippy::large_enum_variant)]
enum Geometry {
    Monomial {
        power: usize,
        multiply: bool,
    },
    Mobius {
        a: Q,
        b: Q,
        c: Q,
        e: Q,
        theorem: &'static str,
    },
    Gamma,
    Symmetric,
    Unsupported,
}

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

impl Geometry {
    fn of(transformation: &RecurrenceTransformation) -> Result<Self, String> {
        let parse = |text: &str| parse_rational_coeff(text);
        Ok(match transformation {
            RecurrenceTransformation::MultiplyPower { power } => Self::Monomial {
                power: *power,
                multiply: true,
            },
            RecurrenceTransformation::DividePower { power } => Self::Monomial {
                power: *power,
                multiply: false,
            },
            RecurrenceTransformation::Reversal { .. } => Self::Mobius {
                a: q(0),
                b: q(1),
                c: q(1),
                e: q(0),
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::FPolynomial { .. } => Self::Mobius {
                a: q(1),
                b: q(0),
                c: q(1),
                e: q(1),
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::Affine { a, b } => Self::Mobius {
                a: parse(a)?,
                b: parse(b)?,
                c: q(0),
                e: q(1),
                theorem: "real_affine_equivalence",
            },
            RecurrenceTransformation::Mobius { a, b, c, e, .. } => Self::Mobius {
                a: parse(a)?,
                b: parse(b)?,
                c: parse(c)?,
                e: parse(e)?,
                theorem: "real_mobius_equivalence",
            },
            RecurrenceTransformation::Gamma { .. } => Self::Gamma,
            RecurrenceTransformation::SymmetricDecomposition { .. } => Self::Symmetric,
            #[allow(unreachable_patterns)]
            _ => Self::Unsupported,
        })
    }

    /// Whether real-rootedness is equivalent and intervals map exactly.
    fn is_scalar_equivalence(&self) -> bool {
        matches!(self, Self::Monomial { .. } | Self::Mobius { .. })
    }

    /// Interval containing the roots of the transformed row when the roots of
    /// the source row lie in `interval` (scalar equivalences only).
    fn forward(&self, interval: &RealInterval) -> Option<RealInterval> {
        match self {
            Self::Monomial { power, multiply } => Some(if *multiply && *power > 0 {
                interval.hull(&RealInterval::point(Q::zero()))
            } else {
                interval.clone()
            }),
            Self::Mobius { a, b, c, e, .. } => {
                let image = mobius_interval_image(interval, e, &-b.clone(), &-c.clone(), a)?;
                Some(if c.is_zero() {
                    image
                } else {
                    // The pole -e/c of phi is a root of multiplicity d_n - deg P_n.
                    image.hull(&RealInterval::point(-(e / c)))
                })
            }
            _ => None,
        }
    }

    /// Interval containing the roots of the source row when the roots of the
    /// transformed row lie in `interval`; `full_degree` means `deg Q_n = d_n`.
    fn backward(&self, interval: &RealInterval, full_degree: bool) -> Option<RealInterval> {
        match self {
            Self::Monomial { power, multiply } => Some(if !*multiply && *power > 0 {
                interval.hull(&RealInterval::point(Q::zero()))
            } else {
                interval.clone()
            }),
            Self::Mobius { a, b, c, e, .. } => {
                let image = mobius_interval_image(interval, a, b, c, e)?;
                Some(if c.is_zero() || full_degree {
                    image
                } else {
                    // deg Q_n < d_n gives roots of P_n at a/c.
                    image.hull(&RealInterval::point(a / c))
                })
            }
            _ => None,
        }
    }

    fn theorem(&self) -> &'static str {
        match self {
            Self::Monomial { .. } => "monomial_factor",
            Self::Mobius { theorem, .. } => theorem,
            Self::Gamma => "gamma_quarter_bound",
            Self::Symmetric => "hermite_kakeya_obreschkoff_symmetric_pair",
            Self::Unsupported => "none",
        }
    }
}

// ---------------------------------------------------------------------------
// Implications
// ---------------------------------------------------------------------------

fn condition_hypotheses(conditions: &[TransformConditionStatus]) -> Vec<RootHypothesis> {
    conditions
        .iter()
        .map(|c| RootHypothesis {
            condition: c.condition.clone(),
            status: if c.status == CERTIFIED_ALL_INDICES {
                HypothesisStatus::CertifiedAllIndices
            } else {
                HypothesisStatus::CheckedSuppliedRows
            },
            detail: c.detail.clone(),
        })
        .collect()
}

fn uncertified(conditions: &[TransformConditionStatus]) -> Vec<String> {
    conditions
        .iter()
        .filter(|c| c.status != CERTIFIED_ALL_INDICES)
        .map(|c| format!("{} (checked on the supplied rows only)", c.condition))
        .collect()
}

/// The theorem adapter of a transformation, with every hypothesis
/// [`HypothesisStatus::NotChecked`].  Reports fill in the statuses.
///
/// ```
/// use polytool::recurrence::*;
///
/// let f = RecurrenceTransformation::FPolynomial {
///     schedule: DegreeSchedule::new(1, 0, 1),
/// };
/// let implication = root_implication(&f).unwrap();
/// assert_eq!(implication.theorem, "real_mobius_equivalence");
/// assert_eq!(implication.direction, RootImplicationDirection::Equivalent);
/// ```
pub fn root_implication(
    transformation: &RecurrenceTransformation,
) -> Result<RootImplication, String> {
    transformation.validate().map_err(|e| e.to_string())?;
    let geometry = Geometry::of(transformation)?;
    let record = transformation.record();
    let hypotheses = record
        .image_conditions
        .iter()
        .map(|condition| RootHypothesis {
            condition: condition.clone(),
            status: HypothesisStatus::NotChecked,
            detail: "not checked".into(),
        })
        .collect();
    let mut implication = RootImplication {
        theorem: geometry.theorem().to_string(),
        direction: RootImplicationDirection::Equivalent,
        source_property: "P_n real-rooted".into(),
        target_property: "Q_n real-rooted".into(),
        statement: String::new(),
        hypotheses,
        degree_bounds: record.degree_schedule.clone(),
        coefficient_signs: None,
        endpoint_conventions: Vec::new(),
        zero_cases: ZERO_CASES.to_string(),
        index_domain: record.index_domain.clone(),
        converse: None,
    };
    match &geometry {
        Geometry::Monomial { power, multiply } => {
            implication.statement = format!(
                "{}: the roots of Q_n are those of P_n with the multiplicity of 0 {} by {power}; \
                 real-rootedness and every interval containing 0 are preserved both ways",
                record.row_map,
                if *multiply { "raised" } else { "lowered" }
            );
            implication.endpoint_conventions = vec![format!(
                "the root 0 of {} has multiplicity at least {power}",
                if *multiply { "Q_n" } else { "P_n" }
            )];
        }
        Geometry::Mobius { a, b, c, e, .. } => {
            let pole = (!c.is_zero()).then(|| format_rational_coeff(&-(e / c)));
            let infinity_image = (!c.is_zero()).then(|| format_rational_coeff(&(a / c)));
            let (a, b, c, e) = (
                format_rational_coeff(a),
                format_rational_coeff(b),
                format_rational_coeff(c),
                format_rational_coeff(e),
            );
            implication.statement = format!(
                "{}: a real Möbius map with a e - b c != 0 is a bijection of the extended real \
                 line mapping non-real points to non-real points; the roots of Q_n are \
                 psi(rho) = (({e}) rho - ({b}))/(({a}) - ({c}) rho) for the roots rho of P_n, so \
                 P_n is real-rooted iff Q_n is, and intervals map by psi",
                record.row_map
            );
            if let (Some(pole), Some(infinity_image)) = (pole, infinity_image) {
                implication.endpoint_conventions = vec![
                    format!(
                        "the pole x = {pole} of the substitution is a root of Q_n of \
                         multiplicity d_n - deg P_n (missing top degree of P_n)"
                    ),
                    format!(
                        "a root rho = {infinity_image} of P_n lowers deg Q_n below d_n; \
                         conversely deg Q_n < d_n gives roots of P_n at {infinity_image}"
                    ),
                ];
            }
            if transformation.kind_name() == "f_polynomial" {
                implication
                    .endpoint_conventions
                    .push("the root 0 of Q_n has the multiplicity of the root 0 of P_n".into());
            }
        }
        Geometry::Gamma => {
            implication.target_property =
                "Q_n (gamma coordinates) real-rooted with every root <= 1/4".into();
            implication.statement = "P_n = (1 + x)^(d_n) Q_n(x/(1 + x)^2) with \
                                     deg Q_n <= floor(d_n/2); each root y of Q_n contributes the \
                                     factor x - y (1 + x)^2 of discriminant 1 - 4y, so P_n is \
                                     real-rooted iff every root of Q_n is real and at most 1/4"
                .into();
            implication.endpoint_conventions = vec![
                "a root 1/4 of Q_n gives the double root x = 1 of P_n".into(),
                "a root 0 of Q_n gives the root x = 0 of P_n and lowers deg P_n".into(),
                "(1 + x)^(d_n - 2 deg Q_n) contributes the root -1".into(),
            ];
            implication.converse = Some(
                "real-rootedness of Q_n alone is insufficient: Q = 1 - z at d = 2 gives \
                 P = 1 + x + x^2"
                    .into(),
            );
        }
        Geometry::Symmetric => {
            implication.direction = RootImplicationDirection::TargetToSource;
            implication.source_property = "P_n = a_n + x b_n real-rooted".into();
            implication.target_property =
                "a_n and x b_n real-rooted and interlacing (in some order)".into();
            implication.statement = "Hermite-Kakeya-Obreschkoff: if a_n and x b_n are \
                                     real-rooted and interlace, every real combination, in \
                                     particular a_n + x b_n, is real-rooted (with nonnegative \
                                     coefficients this is the Brändén-Solus condition that b_n \
                                     interlaces a_n)"
                .into();
            implication.hypotheses.push(RootHypothesis {
                condition: "a_n and x b_n are real-rooted and interlace".into(),
                status: HypothesisStatus::NotChecked,
                detail: "not checked".into(),
            });
            implication.converse = Some(
                "fails: P = (1 + 2x)^2 at d = 2 is real-rooted but a = 1 + x + x^2 is not; and \
                 real-rooted components alone do not suffice: a = (1 + x)^3, b = x at d = 3 give \
                 P = 1 + 3x + 4x^2 + x^3, which is not real-rooted"
                    .into(),
            );
        }
        Geometry::Unsupported => {
            implication.statement = "no root implication is recorded for this map".into();
            implication.direction = RootImplicationDirection::SourceToTarget;
            implication.target_property = "none".into();
        }
    }
    Ok(implication)
}

// ---------------------------------------------------------------------------
// Recurrence theorem: Rolle's theorem for derivatives of products
// ---------------------------------------------------------------------------

/// `Q_m = s(m) D^j (g Q_(m-1))`.
struct RolleForm {
    factor: Vec<Q>,
    order: usize,
    scalar: Vec<Q>,
}

const ROLLE_STATEMENT: &str =
    "if every root of g and of Q_(m-1) is real and lies in an interval I, then \
     s(m) D^j (g Q_(m-1)) is zero or real-rooted with every root in I (Rolle's theorem; the \
     scalar s(m) does not move roots); with the seeds checked exactly this closes an induction";

/// Trimmed rows of a bivariate polynomial (`[n power][x power]`, zero is empty).
fn bivar_rows(p: &BivarPoly) -> Vec<Vec<Q>> {
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

fn binomial(n: usize, k: usize) -> Q {
    Q::from_integer(crate::numeric::binomial_bigint(n, k))
}

fn match_rolle_form(rule: &Recurrence) -> Result<RolleForm, String> {
    let rule = canonical_recurrence(rule);
    if rule.denominator.is_some() {
        return Err("the left factor f(n,x) is not a nonzero constant".into());
    }
    if rule
        .inhomogeneous
        .as_ref()
        .is_some_and(|g| !bivar_rows(g).is_empty())
    {
        return Err("the rule has an additive term".into());
    }
    if rule.terms.is_empty() {
        return Err("the rule has no terms".into());
    }
    if rule.terms.iter().any(|t| t.offset != 1) {
        return Err("only first-order rules (every offset 1) are matched".into());
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err("(-1)^n factors are not matched".into());
    }
    let order = rule.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0);
    let coefficient = |k: usize| -> Vec<Vec<Q>> {
        rule.terms
            .iter()
            .find(|t| t.deriv_order == k)
            .map(|t| bivar_rows(&t.coeff))
            .unwrap_or_default()
    };
    let top = coefficient(order);
    let factor = top
        .iter()
        .find(|row| !row.is_empty())
        .cloned()
        .ok_or("the top coefficient vanishes")?;
    let lead = factor.last().expect("nonzero factor").clone();
    let mut scalar = Vec::with_capacity(top.len());
    for row in &top {
        let s = row.last().map_or_else(Q::zero, |c| c / &lead);
        if vec_poly::trim(row) != vec_poly::trim(&vec_poly::scale(&factor, &s)) {
            return Err("the top coefficient is not of the form s(n) g(x)".into());
        }
        scalar.push(s);
    }
    vec_poly::strip_trailing_zeros(&mut scalar);
    let mut derivative = factor.clone();
    for k in (0..order).rev() {
        derivative = vec_poly::derivative(&derivative);
        let base = vec_poly::scale(&derivative, &binomial(order, k));
        let mut expected: Vec<Vec<Q>> = scalar
            .iter()
            .map(|s| {
                let mut row = vec_poly::scale(&base, s);
                vec_poly::strip_trailing_zeros(&mut row);
                row
            })
            .collect();
        while expected.last().is_some_and(Vec::is_empty) {
            expected.pop();
        }
        if coefficient(k) != expected {
            return Err(format!(
                "the coefficient of D^{k} is not s(n) C({order},{k}) g^({})(x)",
                order - k
            ));
        }
    }
    Ok(RolleForm {
        factor,
        order,
        scalar,
    })
}

/// Whether the polynomial `s(m)` provably never vanishes for `m >= first_target`:
/// one-signed coefficients, and a nonzero constant term or `first_target >= 1`.
fn scalar_nonvanishing(scalar: &[Q], first_target: usize) -> bool {
    let Some(constant) = scalar.first() else {
        return false;
    };
    let one_signed =
        scalar.iter().all(|c| !c.is_negative()) || scalar.iter().all(|c| !c.is_positive());
    one_signed && (!constant.is_zero() || first_target >= 1)
}

fn poly_text(p: &[Q], variable: &str) -> String {
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

fn encode(row: &[Q]) -> Vec<String> {
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

type IsolationKey = (RootSubject, usize, usize);

fn property_holds(
    isolation: &RealRootIsolation,
    property: RootProperty,
    interval: Option<&RealInterval>,
) -> bool {
    match (property, interval) {
        (RootProperty::RealRooted, _) => isolation.is_real_rooted(),
        (RootProperty::RealRootedIn, Some(i)) => isolation.real_rooted_in(i),
        (RootProperty::RealRootsIn, Some(i)) => isolation.all_real_roots_in(i),
        _ => false,
    }
}

fn subject_name(subject: RootSubject, component: usize, names: &[String]) -> String {
    match subject {
        RootSubject::Source => "P".into(),
        RootSubject::Transformed => names.get(component).cloned().unwrap_or_else(|| "Q".into()),
    }
}

fn property_text(property: RootProperty, interval: Option<&RootInterval>) -> String {
    match (property, interval) {
        (RootProperty::RealRooted, _) => "real-rooted".into(),
        (RootProperty::RealRootedIn, Some(i)) => format!("real-rooted in {i}"),
        (RootProperty::RealRootsIn, Some(i)) => format!("every real root in {i}"),
        _ => "malformed".into(),
    }
}

/// Transfer `property` on `interval` through a scalar equivalence.
fn transfer(
    geometry: &Geometry,
    to_target: bool,
    property: RootProperty,
    interval: Option<&RootInterval>,
    full_degree: bool,
) -> Option<(RootProperty, Option<RootInterval>)> {
    if !geometry.is_scalar_equivalence() {
        return None;
    }
    let image = interval
        .filter(|i| i.is_index_independent())
        .and_then(|i| i.at(0).ok())
        .and_then(|i| {
            if to_target {
                geometry.forward(&i)
            } else {
                geometry.backward(&i, full_degree)
            }
        });
    match (property, image) {
        (RootProperty::RealRooted, _) => Some((RootProperty::RealRooted, None)),
        (p, Some(image)) => Some((p, Some(RootInterval::from_real(&image)))),
        (RootProperty::RealRootedIn, None) => Some((RootProperty::RealRooted, None)),
        (RootProperty::RealRootsIn, None) => None,
    }
}

pub(crate) fn analyze_root_domain(
    inputs: &RootDomainInputs<'_>,
) -> Result<RootDomainReport, String> {
    let request = inputs.request;
    request.validate()?;
    let transformation = inputs.transformation;
    let geometry = Geometry::of(transformation)?;
    let record = transformation.record();
    let names: Vec<String> = record.components.iter().map(|c| c.name.clone()).collect();
    let components = transformation.component_count();
    let last_index = transformation.schedule().and_then(|s| s.last_index);
    let domain_first = transformation
        .schedule()
        .map_or(inputs.origin, |s| s.first_index.max(inputs.origin));
    let mut notes = Vec::new();
    let mut claims = Vec::new();

    // Exact isolation of the observed rows.
    let mut isolations: BTreeMap<IsolationKey, RealRootIsolation> = BTreeMap::new();
    let mut observations = Vec::new();
    let mut skipped = Vec::new();
    for (i, row) in inputs
        .rows
        .iter()
        .enumerate()
        .take(request.max_observed_rows)
    {
        let index = inputs.origin + i;
        let state = &inputs.states[i];
        let too_large = std::iter::once(row)
            .chain(state.iter())
            .any(|p| vec_poly::degree(p).is_some_and(|d| d > request.max_observed_degree));
        if too_large {
            skipped.push(index);
            continue;
        }
        let rows = std::iter::once((RootSubject::Source, 0, row)).chain(
            state
                .iter()
                .enumerate()
                .map(|(c, p)| (RootSubject::Transformed, c, p)),
        );
        for (subject, component, polynomial) in rows {
            let isolation = isolate_real_roots_exact(polynomial);
            observations.push(RowRootObservation {
                subject,
                component,
                index,
                degree: isolation.degree,
                valuation: isolation.valuation,
                real_rooted: isolation.is_real_rooted(),
                real_root_count: isolation.real_root_count,
                roots: RootEnclosureJson::from_isolation(&isolation),
            });
            isolations.insert((subject, component, index), isolation);
        }
    }
    if inputs.rows.len() > request.max_observed_rows {
        notes.push(format!(
            "roots isolated on the first {} of {} rows (max_observed_rows)",
            request.max_observed_rows,
            inputs.rows.len()
        ));
    }
    if !skipped.is_empty() {
        notes.push(format!(
            "rows {skipped:?} exceed max_observed_degree = {} and were not isolated",
            request.max_observed_degree
        ));
    }

    // Finite exact observation claims, one per subject and component.
    let subjects = std::iter::once((RootSubject::Source, 0))
        .chain((0..components).map(|c| (RootSubject::Transformed, c)));
    for (subject, component) in subjects {
        let observed: Vec<(&usize, &RealRootIsolation)> = isolations
            .iter()
            .filter(|((s, c, _), _)| *s == subject && *c == component)
            .map(|((_, _, index), isolation)| (index, isolation))
            .collect();
        let (Some(first), Some(last)) = (observed.first(), observed.last()) else {
            continue;
        };
        let failing: Vec<usize> = observed
            .iter()
            .filter(|(_, iso)| !iso.is_real_rooted())
            .map(|(index, _)| **index)
            .collect();
        let hull = observed
            .iter()
            .filter_map(|(_, iso)| iso.enclosure_hull())
            .reduce(|a, b| a.hull(&b));
        let name = subject_name(subject, component, &names);
        let (property, interval) = match (&hull, failing.is_empty()) {
            (Some(h), true) => (RootProperty::RealRootedIn, Some(RootInterval::from_real(h))),
            (None, true) => (RootProperty::RealRooted, None),
            (Some(h), false) => (RootProperty::RealRootsIn, Some(RootInterval::from_real(h))),
            (None, false) => {
                notes.push(format!(
                    "{name}_n has no real roots and is not real-rooted at the observed indices \
                     {failing:?}"
                ));
                continue;
            }
        };
        let mut derivation = format!(
            "exact root isolation of {name}_n at the {} observed indices; the interval is the \
             closed hull of the rational isolating enclosures",
            observed.len()
        );
        if !failing.is_empty() {
            derivation.push_str(&format!("; not real-rooted at indices {failing:?}"));
        }
        claims.push(RootDomainClaim {
            subject,
            component,
            property,
            interval,
            first_index: *first.0,
            last_index: Some(*last.0),
            evidence: RootEvidence::FiniteExactObservation,
            conditional_on: Vec::new(),
            derivation,
        });
    }

    // The map's implication with hypothesis statuses.
    let mut implication = root_implication(transformation)?;
    let base_hypotheses = condition_hypotheses(inputs.conditions);
    let extra = implication
        .hypotheses
        .split_off(record.image_conditions.len());
    implication.hypotheses = base_hypotheses;
    for mut hypothesis in extra {
        if matches!(geometry, Geometry::Symmetric) {
            let mut checked = 0;
            let mut failed = Vec::new();
            let observed = inputs
                .states
                .iter()
                .enumerate()
                .take(request.max_observed_rows)
                .filter(|(_, state)| {
                    state.iter().all(|p| {
                        vec_poly::degree(p).is_none_or(|d| d <= request.max_observed_degree)
                    })
                });
            for (i, state) in observed {
                let a = integer_row(&state[0]);
                let xb = integer_row(&vec_poly::shift(&state[1]));
                checked += 1;
                if !interlaces_either_order_bigint_coeffs(&a, &xb) {
                    failed.push(inputs.origin + i);
                }
            }
            (hypothesis.status, hypothesis.detail) = if failed.is_empty() {
                (
                    HypothesisStatus::CheckedSuppliedRows,
                    format!("exact on the {checked} observed rows; finite evidence only"),
                )
            } else {
                (
                    HypothesisStatus::FailedOnRows,
                    format!("fails at indices {failed:?}"),
                )
            };
        }
        implication.hypotheses.push(hypothesis);
    }
    let conditional = uncertified(inputs.conditions);
    let all_certified = conditional.is_empty();
    let mut implications = vec![implication];

    // The nonnegative-coefficient window.
    let mut candidates: Vec<(RealInterval, String)> = Vec::new();
    for interval in &request.candidate_intervals {
        match interval.is_index_independent() {
            true => candidates.push((interval.at(0)?, "supplied".into())),
            false => notes.push(format!(
                "candidate interval {interval} is index-dependent and is not used in the \
                 induction"
            )),
        }
    }
    if request.default_window {
        if let Some(window) = geometry.forward(&RealInterval::at_most(Q::zero())) {
            let negative: Vec<usize> = inputs
                .rows
                .iter()
                .enumerate()
                .filter(|(_, row)| row.iter().any(Signed::is_negative))
                .map(|(i, _)| inputs.origin + i)
                .collect();
            let window_json = RootInterval::from_real(&window);
            let (status, detail) = if negative.is_empty() {
                (
                    HypothesisStatus::CheckedSuppliedRows,
                    format!(
                        "exact on all {} supplied rows; finite evidence only",
                        inputs.rows.len()
                    ),
                )
            } else {
                (
                    HypothesisStatus::FailedOnRows,
                    format!("negative coefficients at indices {negative:?}"),
                )
            };
            let mut hypotheses = condition_hypotheses(inputs.conditions);
            hypotheses.push(RootHypothesis {
                condition: "P_n has nonnegative coefficients".into(),
                status,
                detail,
            });
            let f_polynomial = matches!(
                &geometry,
                Geometry::Mobius { a, b, c, e, .. }
                    if a.is_one() && b.is_zero() && c.is_one() && e.is_one()
            );
            implications.push(RootImplication {
                theorem: "nonnegative_coefficient_window".into(),
                direction: RootImplicationDirection::SourceToTarget,
                source_property: "P_n nonzero with nonnegative coefficients (every real root \
                                  in (-inf, 0])"
                    .into(),
                target_property: format!("every real root of Q_n in {window_json}"),
                statement: if f_polynomial {
                    "Q = (1 + x)^d P(x/(1 + x)) = sum_j p_j x^j (1 + x)^(d - j) and \
                     (-1)^d Q(-1 - x) = sum_j p_j (1 + x)^j x^(d - j) have nonnegative \
                     coefficients and deg Q <= d, so every real root of Q lies in [-1, 0]; \
                     P is real-rooted iff Q is.  The window does not prove that Q splits."
                        .into()
                } else {
                    format!(
                        "the real roots of P_n lie in (-inf, 0] and the map sends them into \
                         {window_json}; the window does not prove real-rootedness"
                    )
                },
                hypotheses,
                degree_bounds: record.degree_schedule.clone(),
                coefficient_signs: Some("p_(n,j) >= 0 for every j".into()),
                endpoint_conventions: implications[0].endpoint_conventions.clone(),
                zero_cases: ZERO_CASES.to_string(),
                index_domain: record.index_domain.clone(),
                converse: Some(
                    "fails: P = 1 + z + z^2 has nonnegative coefficients, the window holds \
                     vacuously for Q = 1 + 3x + 3x^2, and neither is real-rooted"
                        .into(),
                ),
            });
            if negative.is_empty() && transformation.component_count() == 1 {
                let mut conditional_on = vec![
                    "P_n has nonnegative coefficients for every index (checked on the supplied \
                     rows only)"
                        .to_string(),
                ];
                conditional_on.extend(conditional.iter().cloned());
                claims.push(RootDomainClaim {
                    subject: RootSubject::Transformed,
                    component: 0,
                    property: RootProperty::RealRootsIn,
                    interval: Some(window_json.clone()),
                    first_index: domain_first,
                    last_index,
                    evidence: RootEvidence::ProposedInterval,
                    conditional_on,
                    derivation: "nonnegative_coefficient_window; a candidate window for \
                                 proofs, not a proof of real-rootedness"
                        .into(),
                });
                candidates.push((window, "default_window".into()));
            }
        }
    }

    // Supplied assumptions.
    let mut assumption_checks = Vec::new();
    for (k, assumption) in request.assumptions.iter().enumerate() {
        if assumption.subject == RootSubject::Transformed && assumption.component >= components {
            return Err(format!(
                "assumption {k}: component {} does not exist",
                assumption.component
            ));
        }
        let first = assumption.first_index.unwrap_or(domain_first);
        let mut checked = Vec::new();
        let mut refuted_at = None;
        for ((subject, component, index), isolation) in &isolations {
            if *subject != assumption.subject
                || *component != assumption.component
                || *index < first
                || assumption.last_index.is_some_and(|last| *index > last)
            {
                continue;
            }
            let interval = assumption
                .interval
                .as_ref()
                .map(|i| i.at(*index))
                .transpose()?;
            checked.push(*index);
            if refuted_at.is_none()
                && !property_holds(isolation, assumption.property, interval.as_ref())
            {
                refuted_at = Some(*index);
            }
        }
        let status = match (checked.is_empty(), refuted_at) {
            (_, Some(_)) => "refuted_on_rows",
            (true, None) => "not_checked",
            (false, None) => "consistent_with_rows",
        };
        assumption_checks.push(AssumptionCheck {
            assumption: k,
            status: status.into(),
            checked_indices: checked,
            refuted_at,
        });
        let name = subject_name(assumption.subject, assumption.component, &names);
        let label = format!(
            "assumption {k}: {name}_n {} for n >= {first}{}",
            property_text(assumption.property, assumption.interval.as_ref()),
            assumption
                .last_index
                .map_or(String::new(), |last| format!(" and n <= {last}"))
        );
        claims.push(RootDomainClaim {
            subject: assumption.subject,
            component: assumption.component,
            property: assumption.property,
            interval: assumption.interval.clone(),
            first_index: first,
            last_index: assumption.last_index,
            evidence: RootEvidence::SuppliedAssumption,
            conditional_on: vec![label.clone()],
            derivation: format!(
                "supplied in the request{}",
                assumption
                    .note
                    .as_ref()
                    .map_or(String::new(), |note| format!(" ({note})"))
            ),
        });
        if refuted_at.is_some() {
            notes.push(format!(
                "{label} is refuted on the rows; nothing is transferred"
            ));
            continue;
        }
        let to_target = assumption.subject == RootSubject::Source;
        if let Some((property, interval)) = transfer(
            &geometry,
            to_target,
            assumption.property,
            assumption.interval.as_ref(),
            false,
        ) {
            let mut conditional_on = vec![label];
            conditional_on.extend(conditional.iter().cloned());
            claims.push(RootDomainClaim {
                subject: if to_target {
                    RootSubject::Transformed
                } else {
                    RootSubject::Source
                },
                component: 0,
                property,
                interval,
                first_index: first,
                last_index: assumption.last_index,
                evidence: RootEvidence::SuppliedAssumption,
                conditional_on,
                derivation: format!("transferred from assumption {k} by {}", geometry.theorem()),
            });
        } else if matches!(geometry, Geometry::Gamma) {
            let quarter = RealInterval::at_most(Q::new(BigInt::one(), BigInt::from(4)));
            let transferred = match (to_target, assumption.property) {
                (true, RootProperty::RealRooted | RootProperty::RealRootedIn) => Some((
                    RootSubject::Transformed,
                    RootProperty::RealRootedIn,
                    Some(RootInterval::from_real(&quarter)),
                )),
                (false, RootProperty::RealRootedIn)
                    if assumption.interval.as_ref().is_some_and(|i| {
                        i.is_index_independent() && i.at(0).is_ok_and(|i| i.is_subset_of(&quarter))
                    }) =>
                {
                    Some((RootSubject::Source, RootProperty::RealRooted, None))
                }
                _ => None,
            };
            if let Some((subject, property, interval)) = transferred {
                let mut conditional_on = vec![label];
                conditional_on.extend(conditional.iter().cloned());
                claims.push(RootDomainClaim {
                    subject,
                    component: 0,
                    property,
                    interval,
                    first_index: first,
                    last_index: assumption.last_index,
                    evidence: RootEvidence::SuppliedAssumption,
                    conditional_on,
                    derivation: format!("transferred from assumption {k} by gamma_quarter_bound"),
                });
            }
        }
    }

    // The recurrence theorem match of the transported rule.
    let recurrence_theorem = inputs.transported.as_ref().map(|transported| {
        rolle_induction(
            inputs,
            transported,
            &geometry,
            &candidates,
            all_certified,
            &conditional,
            &mut claims,
            &mut notes,
        )
    });

    Ok(RootDomainReport {
        conventions: conventions(),
        implications,
        claims,
        recurrence_theorem,
        assumption_checks,
        observations,
        notes,
    })
}

fn integer_row(row: &[Q]) -> Vec<BigInt> {
    let scale = row.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    row.iter()
        .map(|c| (c * Q::from_integer(scale.clone())).to_integer())
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn rolle_induction(
    inputs: &RootDomainInputs<'_>,
    transported: &TransportedInput<'_>,
    geometry: &Geometry,
    supplied: &[(RealInterval, String)],
    all_certified: bool,
    conditional: &[String],
    claims: &mut Vec<RootDomainClaim>,
    notes: &mut Vec<String>,
) -> RecurrenceRootTheorem {
    let mut report = RecurrenceRootTheorem {
        theorem: "rolle_derivative_of_product".into(),
        statement: ROLLE_STATEMENT.into(),
        matched: false,
        reason: None,
        form: None,
        factor: None,
        derivative_order: None,
        index_scalar: None,
        degree_formula: None,
        attempts: Vec::new(),
    };
    let rule = match transported.rule {
        TransportedRecurrence::Scalar(rule) => rule,
        _ => {
            report.reason = Some("only scalar transported rules are matched".into());
            return report;
        }
    };
    let form = match match_rolle_form(rule) {
        Ok(form) => form,
        Err(reason) => {
            report.reason = Some(reason);
            return report;
        }
    };
    let names: Vec<String> = inputs
        .transformation
        .record()
        .components
        .into_iter()
        .map(|c| c.name)
        .collect();
    let name = names.first().cloned().unwrap_or_else(|| "Q".into());
    let first = transported.first_index;
    let valid_from = transported.valid_from;
    let last_index = inputs.transformation.schedule().and_then(|s| s.last_index);
    report.matched = true;
    report.form = Some(format!(
        "{name}_m = s(m) D^{} (g {name}_(m-1)) for m >= {valid_from}, s(m) = {}, g = {}",
        form.order,
        poly_text(&form.scalar, "m"),
        poly_text(&form.factor, "x")
    ));
    report.factor = Some(encode(&form.factor));
    report.derivative_order = Some(form.order);
    report.index_scalar = Some(encode(&form.scalar));

    let seeds: Vec<(usize, &Vec<Q>)> = (first..valid_from)
        .filter_map(|k| {
            k.checked_sub(inputs.origin)
                .and_then(|i| inputs.states.get(i))
                .map(|state| (k, &state[0]))
        })
        .collect();
    if seeds.len() != valid_from - first {
        report.reason = Some("the transformed seeds are not available".into());
        report.matched = false;
        return report;
    }

    // Exact degree tracking.
    let factor_degree = form.factor.len() - 1;
    let (last_seed_index, last_seed) = *seeds.last().expect("at least one seed");
    let mut full_degree = false;
    if factor_degree >= form.order
        && scalar_nonvanishing(&form.scalar, valid_from)
        && !vec_poly::is_zero(last_seed)
    {
        let start = vec_poly::degree(last_seed).expect("nonzero seed");
        let slope = factor_degree - form.order;
        report.degree_formula = Some(format!(
            "deg {name}_n = {start} + {slope} (n - {last_seed_index}) for n >= {last_seed_index}"
        ));
        full_degree = inputs.transformation.schedule().is_some_and(|schedule| {
            i64::try_from(slope).is_ok_and(|s| s == schedule.step)
                && seeds.iter().all(|(k, seed)| {
                    schedule
                        .degree_at(*k)
                        .is_ok_and(|d| vec_poly::degree(seed) == Some(d))
                })
        });
    }

    // Candidate intervals: supplied, default window, then the root hull.
    let factor_isolation = isolate_real_roots_exact(&form.factor);
    let seed_isolations: Vec<(usize, RealRootIsolation)> = seeds
        .iter()
        .map(|(k, seed)| (*k, isolate_real_roots_exact(seed)))
        .collect();
    let mut candidates: Vec<(RealInterval, String)> = supplied.to_vec();
    let hull = std::iter::once(&factor_isolation)
        .chain(seed_isolations.iter().map(|(_, iso)| iso))
        .filter_map(RealRootIsolation::enclosure_hull)
        .reduce(|a, b| a.hull(&b));
    if let Some(hull) = hull {
        candidates.push((hull, "root_hull".into()));
    }
    let limit = inputs.request.max_candidate_intervals;
    if candidates.len() > limit {
        notes.push(format!(
            "{} candidate intervals were not tried (max_candidate_intervals = {limit})",
            candidates.len() - limit
        ));
        candidates.truncate(limit);
    }
    let mut success = None;
    for (interval, origin) in candidates {
        let detail = if !factor_isolation.real_rooted_in(&interval) {
            format!(
                "g = {} is not real-rooted in {interval}",
                poly_text(&form.factor, "x")
            )
        } else if let Some((k, _)) = seed_isolations
            .iter()
            .find(|(_, iso)| !iso.real_rooted_in(&interval))
        {
            format!("the seed {name}_{k} is not real-rooted in {interval}")
        } else {
            String::new()
        };
        let ok = detail.is_empty();
        report.attempts.push(IntervalAttempt {
            interval: RootInterval::from_real(&interval),
            origin,
            success: ok,
            detail: if ok {
                format!(
                    "exact: g and the seeds {name}_{first}..{name}_{} are real-rooted in \
                     {interval}",
                    valid_from - 1
                )
            } else {
                detail
            },
        });
        if ok {
            success = Some(interval);
            break;
        }
    }
    let Some(interval) = success else {
        notes.push("no candidate interval closes the Rolle induction".into());
        return report;
    };
    let interval_json = RootInterval::from_real(&interval);
    claims.push(RootDomainClaim {
        subject: RootSubject::Transformed,
        component: 0,
        property: RootProperty::RealRootedIn,
        interval: Some(interval_json.clone()),
        first_index: first,
        last_index,
        evidence: RootEvidence::CompletedRecurrenceInduction,
        conditional_on: conditional.to_vec(),
        derivation: format!(
            "rolle_derivative_of_product on the transported rule (an exact consequence of the \
             defining source rule) with g and the seeds checked exactly on {interval}; the rows \
             it generates are the transformed rows"
        ),
    });
    if !all_certified {
        notes.push(
            "the induction is not transferred to P_n: some map conditions are only checked on \
             the supplied rows"
                .into(),
        );
        return report;
    }
    match transfer(
        geometry,
        false,
        RootProperty::RealRootedIn,
        Some(&interval_json),
        full_degree,
    ) {
        Some((property, source_interval)) => claims.push(RootDomainClaim {
            subject: RootSubject::Source,
            component: 0,
            property,
            interval: source_interval,
            first_index: first,
            last_index,
            evidence: RootEvidence::CertifiedFamilyInvariant,
            conditional_on: Vec::new(),
            derivation: format!(
                "the completed induction for {name}_n and {} (target to source; every map \
                 condition certified for all indices){}",
                geometry.theorem(),
                if full_degree {
                    "; deg Q_n = d_n by the degree formula, so no root of P_n is added at the \
                     image of infinity"
                } else {
                    ""
                }
            ),
        }),
        None => notes.push(format!(
            "the induction is not transferred to P_n: {} gives no scalar equivalence",
            geometry.theorem()
        )),
    }
    report
}

#[cfg(test)]
#[path = "root_domain_tests.rs"]
mod tests;
