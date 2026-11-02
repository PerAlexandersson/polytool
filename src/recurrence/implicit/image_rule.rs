//! Exact image recurrences: a runnable rule for `Q_n = lambda(n) T_n(P_(n+s))`
//! from a runnable rule of `P`, with seeds and a replayable certificate
//! ([`image_recurrence_of_operator`]).
//!
//! # Input
//!
//! The parent is a scalar rule `f(m) P_m = sum c(m, x) D^d P_(m-r)` with
//! seeds ([`RecurrenceJson`]): homogeneous, without `(-1)^n` factors, with a
//! leader `f` free of `x`, decided nonzero at every generated target.  The
//! image is described by an [`ImageOperatorSpecification`]: an index shift
//! `s`, a row map `T_n` ([`ImageRowMap`]: a Weyl operator
//! `(1/delta(n)) sum c_d(n, x) D^d`, an affine map, a reversal, a Möbius map
//! or a coefficient multiplier `x^k -> u(k) x^k`) and a gauge `lambda(n)`
//! ([`ImageGauge`]), for every `n >= n0` (the image's first index,
//! `n0 + s >= p0`).
//!
//! # Method
//!
//! The rows `P'_n = P_(n+s)` satisfy the shifted parent relation `L'`
//! (`n -> n + s`).  The map is a list of steps of the `ore-algebra` crate
//! (`RowMapStep`), the gauge (with `1/delta`) last.  The crate first pushes
//! `L'` forward through the steps (coordinate maps, gauges and multipliers,
//! same order); when that is not runnable (always for a Weyl operator) it
//! searches a common left multiple `M` of an increasing ansatz shape whose
//! pullback lies in the left ideal of `L'`.  When auxiliary annihilators
//! are enabled ([`ImageRecurrenceOptions::row_annihilators`]), two kinds are
//! fitted to the parent rows first and kept only with a proof (step 2
//! below): a runnable rule `B` of lower order than `L'` (for example the
//! first-order contiguity relation, with a derivative, of a family given by
//! a three-term recurrence), which then replaces `L'`, and a shift-free row
//! annihilator `A`, which enlarges the ideal (every rule found without it is
//! found with it).  Some images need them: the binomial rows `(1 + x)^n`
//! need `(1 + x) D - n` for their image under `1/k!`, and with the
//! contiguity relation of the Narayana rows `theta (theta + 1) N_n` has a
//! first-order rule.  A push-forward is certified without `A` when possible.
//! The result is either a runnable `M` with a certificate or the typed
//! error [`RecurrenceError::BudgetExhausted`] (nothing is guessed).
//!
//! # Soundness
//!
//! Write `M = sum m_(i,d)(n, x) D^d S^i` (base `n`, top shift `R`, top
//! column `m_R(n, x) S^R` without derivatives) and `Q_n` for the image rows.
//!
//! 1. **Identity.**  The certificate records `w(n, x) E = H o L + Z o A`
//!    with `(E, b) = pull_back(M, steps)` (the crate's
//!    `ImageAnnihilatorCertificate`, recomputed exactly), where `L` is `B`
//!    when a reduced parent is recorded and `L'` otherwise.  The pullback
//!    contract gives `(M Q)_n = 0` whenever `(E P')_n = 0` and `b(n) != 0`
//!    (for the gauge, `b(n) = prod_(t=1..R) beta(n + t)` and `lambda`,
//!    `delta` are decided finite and nonzero at every `n >= n0`; for a
//!    multiplier `x^k -> u(k) x^k` with `u(k + 1)/u(k) = A(k)/B(k)`, `A` and
//!    `B` are decided to have no integer root `>= 0`, and `A` none in
//!    `[-S, -1]`).  The identity gives `(E P')_n = 0` at every base
//!    `n >= max(n0, b_L)` where `w(n, .)` is not the zero polynomial
//!    (`b_L` the first base of `L`: that of `L'`, or `n0` for `B`; and `A`
//!    holds from `n0` on).
//! 2. **Auxiliary annihilators.**  `B` (relative to `L'`) and `A` (relative
//!    to `L`) are proved by invariance identities
//!    `W (S^j o C) = Y o L_0 + Z_C o C` (`j` the order of `L_0` minus that
//!    of `C`; recomputed) and exact checks of `(C P')_m = 0` at the bases
//!    `m = n0, ..., max(n0, b_(L_0)) + j - 1` and at `n + j` for every base
//!    `n` where `W(n, .)` vanishes; by induction `C` holds at every base
//!    `m >= n0` (the crate's `RowAnnihilatorCertificate`).
//! 3. **Exceptional bases.**  The bases `n >= n0` not covered by 1 (below
//!    `b_L`, integer roots of the content of `w` in `n` or of `b`) are
//!    finitely many; `(M Q)_n = 0` is checked there exactly on the rows
//!    `Q_n` computed from the parent rule.  So `M` holds at every base
//!    `n >= n0`.
//! 4. **Runnable normalisation.**  The leader `m_R(n, x)` vanishes
//!    identically only at the integer roots of its content in `n`; the rule
//!    keeps every row up to the last such target as a seed (at least `R`
//!    seeds), and the leader is decided nonzero at every later target by
//!    the sign checker.  The seeds are computed exactly as `lambda(n)
//!    T_n(P'_n)`.  By induction on the target `t`, the generated row and
//!    `Q_t` satisfy `m_R(t - R, x) Y_t = (the same right-hand side)`, and
//!    `Q[x]` has no zero divisors, so every generated row equals `Q_t`.
//!    (For a reversal, `Q_t` is a Laurent polynomial when
//!    `deg P'_t > d_t`; the identity holds in Laurent polynomials, so a row
//!    that the rule generates is `Q_t`, and the generator stops with an
//!    error where `Q_t` is not a polynomial.)
//!
//! A replay of the rule against freshly computed `Q_n` on further rows is
//! evidence only.  [`ImageRecurrenceCertificate::verify`] recomputes every
//! step above from the recorded parent, specification, rule and identity.

use super::{biv_from_json, biv_json, BigRational as Q, ImplicitRelation, ImplicitRelationJson};
use crate::integer_range_sign::{
    decide_index_polynomial_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::recurrence::implicit::algebra::{
    content_and_primitive, roots_in_domain, Biv, BivarConversion, DOp,
};
use crate::recurrence::root_domain::bivariate_leader_nonvanishing;
use crate::recurrence::{
    parse_rational_coeff, BivarPoly, BivarPolyJson, DegreeSchedule, OperatorTermJson, Recurrence,
    RecurrenceError, RecurrenceJson, RecurrenceTerm,
};
use crate::vec_poly;
use num_traits::{One, Zero};
use ore_algebra::{
    find_annihilator_system, find_image_annihilator, AnnihilatorSystemBudget,
    AnnihilatorSystemCertificate, AnnihilatorSystemEntry, ImageAnnihilatorBudget,
    ImageAnnihilatorCertificate, ImageAnnihilatorMethod, ImageAnnihilatorSearch, PullbackBounds,
    RowMapStep,
};
use serde::{Deserialize, Serialize};

/// Schema of [`ImageRecurrence`].
pub const IMAGE_RECURRENCE_SCHEMA: &str = "polytool.image_recurrence.v1";

/// Largest index shift `|s|`.
pub const MAX_IMAGE_RECURRENCE_INDEX_SHIFT: usize = 64;

/// Largest number of ansatz unknowns ([`ImageRecurrenceOptions::max_unknowns`]).
pub const MAX_IMAGE_RECURRENCE_UNKNOWNS: usize = 2000;

/// Largest number of ansatz shapes ([`ImageRecurrenceOptions::max_shapes`]).
pub const MAX_IMAGE_RECURRENCE_SHAPES: usize = 500;

/// Parent rows used to fit auxiliary annihilators.
const AUXILIARY_FIT_ROWS: usize = 32;

/// Largest number of rows, after the first index, that the construction
/// and its verification compute (seeds, exceptional bases, replay).
pub const MAX_IMAGE_RECURRENCE_ROWS: usize = 512;

/// The gauge `lambda(n)` of an image `Q_n = lambda(n) T_n(P_(n+s))`.
///
/// Polynomials in `n` are ascending exact rational strings.  Every form must
/// be finite and nonzero at every `n >= n0`; this is decided exactly.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ImageGauge {
    /// `lambda(n) = 1`.
    #[default]
    None,
    /// `lambda(n0) = initial_value` and
    /// `lambda(n)/lambda(n-1) = numerator(n)/denominator(n)` for `n > n0`
    /// (both nonzero at every integer `n > n0`).
    Ratio {
        /// `a(n)`, ascending coefficients.
        numerator: Vec<String>,
        /// `b(n)`, ascending coefficients.
        denominator: Vec<String>,
        /// `lambda(n0)`, nonzero.
        initial_value: String,
    },
    /// `lambda(n) = numerator(n)/denominator(n)` (both nonzero at every
    /// integer `n >= n0`).
    Closed {
        /// Ascending coefficients.
        numerator: Vec<String>,
        /// Ascending coefficients.
        denominator: Vec<String>,
    },
}

fn one_polynomial() -> Vec<String> {
    vec!["1".into()]
}

/// The row map `T_n` of an image (module documentation).  Index-dependent
/// quantities are functions of the image index `n`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum ImageRowMap {
    /// `T_n = 1` (only the index shift and the gauge act).
    Identity,
    /// `T_n = (1/denominator(n)) sum_d c_d(n, x) D^d`.
    WeylOperator {
        /// The terms `c_d(n, x) D^d`.
        terms: Vec<OperatorTermJson>,
        /// `delta(n)`, ascending; nonzero at every `n >= n0`.  Default `1`.
        #[serde(default = "one_polynomial")]
        denominator: Vec<String>,
    },
    /// `P -> P(a x + b)`, `a != 0`.
    Affine {
        /// `a`.
        slope: String,
        /// `b`.
        intercept: String,
    },
    /// `P -> x^(d_n) P(1/x)`; the schedule is unbounded, with a nonnegative
    /// step, and starts at or before `n0`.
    Reversal {
        /// `d_n`.
        schedule: DegreeSchedule,
    },
    /// `P -> (c x + e)^(d_n) P((a x + b)/(c x + e))` with `c != 0` and
    /// `a e - b c != 0`; the schedule is as for [`ImageRowMap::Reversal`].
    Mobius {
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `c`.
        c: String,
        /// `e`.
        e: String,
        /// `d_n`.
        schedule: DegreeSchedule,
    },
    /// `x^k -> u(k) x^k` with `u(0) = 1` and
    /// `u(k + 1)/u(k) = numerator(k)/denominator(k)`; both polynomials must
    /// have no integer root `k >= 0` (so `u` is finite and nonzero), and the
    /// numerator none in `[-S, -1]` for the Euler bound `S` of the
    /// certificate.  For example `1/k!` is `numerator = [1]`,
    /// `denominator = [1, 1]`, and `1/(alpha)_k` is `denominator = [alpha, 1]`.
    CoefficientMultiplier {
        /// `A(k)`, ascending.
        numerator: Vec<String>,
        /// `B(k)`, ascending.
        denominator: Vec<String>,
    },
}

/// The image `Q_n = lambda(n) T_n(P_(n+s))`, `n >= n0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageOperatorSpecification {
    /// `s`, `|s| <=` [`MAX_IMAGE_RECURRENCE_INDEX_SHIFT`].  Default 0.
    #[serde(default)]
    pub index_shift: i64,
    /// `lambda`.  Default none.
    #[serde(default)]
    pub gauge: ImageGauge,
    /// `T_n`.
    pub map: ImageRowMap,
}

impl ImageOperatorSpecification {
    /// The map `map` without a shift or a gauge.
    pub fn new(map: ImageRowMap) -> Self {
        Self {
            index_shift: 0,
            gauge: ImageGauge::None,
            map,
        }
    }
}

/// Options of [`image_recurrence_of_operator`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ImageRecurrenceOptions {
    /// Largest order above the parent's order.  Default 2.
    pub max_order_increase: usize,
    /// Largest derivative order of the ansatz.  Default 4.
    pub max_derivative_order: usize,
    /// Largest coefficient degree in `n` and in `x`.  Default 6.
    pub max_coefficient_degree: usize,
    /// Largest number of unknowns of one ansatz, at most
    /// [`MAX_IMAGE_RECURRENCE_UNKNOWNS`].  Default 400.
    pub max_unknowns: usize,
    /// Largest number of ansatz shapes solved per search, at most
    /// [`MAX_IMAGE_RECURRENCE_SHAPES`].  Default 60.
    pub max_shapes: usize,
    /// Fit and prove a shift-free annihilator of the parent rows when the
    /// recurrence alone gives nothing.  Default true.
    pub row_annihilators: bool,
    /// Rows of the rule replayed against the image after the seeds
    /// (evidence only).  Default 8.
    pub replay_rows: usize,
}

impl Default for ImageRecurrenceOptions {
    fn default() -> Self {
        Self {
            max_order_increase: 2,
            max_derivative_order: 4,
            max_coefficient_degree: 6,
            max_unknowns: 400,
            max_shapes: 60,
            row_annihilators: true,
            replay_rows: 8,
        }
    }
}

/// How the rule was found.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum ImageRecurrenceMethod {
    /// The parent rule pushed forward through the map (same order).
    PushForward,
    /// A common left multiple from an ansatz of this shape.
    LeftMultiple {
        /// The order `R`.
        order: usize,
        /// The derivative order below the top shift.
        derivative_order: usize,
        /// The coefficient degree in `n`.
        index_degree: usize,
        /// The coefficient degree in `x`.
        variable_degree: usize,
    },
}

/// The bounds that fix the pullback of the rule (the crate's
/// `PullbackBounds`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageRecurrenceBounds {
    /// At least the order of the rule.
    pub top_shift: usize,
    /// At least its derivative order.
    pub derivative_order: usize,
    /// At least the degree in `x` of its coefficients.
    pub variable_degree: usize,
}

/// One member `C` of an annihilator system with its invariance identity
/// `W (S^j o C) = Y o L' + sum_l Z_l o C_l` (module documentation, step 2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageAnnihilatorSystemMember {
    /// `C` (base form): shift-free (a row annihilator) or a runnable rule of
    /// lower order than `L'`.
    pub annihilator: ImplicitRelationJson,
    /// `W(n, x)`.
    pub multiplier: BivarPolyJson,
    /// `Y`.
    pub recurrence_quotient: ImplicitRelationJson,
    /// `Z_l`, one per member, each of shift below the step.
    pub member_quotients: Vec<ImplicitRelationJson>,
    /// The bases `m` (image indexing) where `(C P')_m = 0` was checked
    /// exactly.
    pub checked_bases: Vec<usize>,
}

/// Additional annihilators of the shifted parent rows `P'`, jointly
/// invariant under `L'` (module documentation, step 2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageAnnihilatorSystem {
    /// The step `j`.
    pub step: usize,
    /// The members.
    pub members: Vec<ImageAnnihilatorSystemMember>,
}

/// The certificate of an [`ImageRecurrence`] (module documentation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageRecurrenceCertificate {
    /// The parent rule with its first index `p0` and seeds.
    pub parent: RecurrenceJson,
    /// The image's first index `n0`.
    pub first_index: usize,
    /// The image.
    pub specification: ImageOperatorSpecification,
    /// The runnable image rule with its seeds (first index `n0`).
    pub rule: RecurrenceJson,
    /// The pullback bounds.
    pub bounds: ImageRecurrenceBounds,
    /// `w(n, x)` of `w E = H o L + Z o A` (`L` is `B` when a reduced parent
    /// is recorded, else `L'`).
    pub multiplier: BivarPolyJson,
    /// `H` (terms with shifts).
    pub quotient: ImplicitRelationJson,
    /// Additional annihilators of `P'` with their joint invariance, when
    /// used: a member `B` of positive order replaces `L'` in the identity,
    /// a shift-free member is the row annihilator `A`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub annihilator_system: Option<ImageAnnihilatorSystem>,
    /// `Z` (empty without `A`).
    pub row_quotient: ImplicitRelationJson,
    /// The bases `n >= n0` not covered by the identity, where `(M Q)_n = 0`
    /// was checked exactly.
    pub exceptional_bases: Vec<usize>,
    /// Rows replayed after the seeds (evidence).
    pub replay_rows: usize,
}

/// A runnable image rule with its certificate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImageRecurrence {
    /// [`IMAGE_RECURRENCE_SCHEMA`].
    pub schema: String,
    /// The runnable rule of `Q` with its seeds (first index `n0`).
    pub rule: RecurrenceJson,
    /// The order `R` of the rule.
    pub order: usize,
    /// How the rule was found.
    pub method: ImageRecurrenceMethod,
    /// The certificate.
    pub certificate: ImageRecurrenceCertificate,
}

impl ImageRecurrence {
    /// Verify the certificate and that it records this rule.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != IMAGE_RECURRENCE_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported image recurrence schema `{}`",
                self.schema
            )));
        }
        if self.rule != self.certificate.rule {
            return Err(RecurrenceError::VerificationFailed(
                "the rule differs from the certified rule".into(),
            ));
        }
        let order = self.certificate.verify()?;
        if order != self.order {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the recorded order {} differs from the rule's order {order}",
                self.order
            )));
        }
        Ok(())
    }

    /// The first `count` rows of the rule (seeds included), exact.
    pub fn rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let (rule, first_index, seeds) = self.rule.to_recurrence_parts()?;
        rule.generate_rows_rational(&seeds, first_index, count.max(seeds.len()))
            .map(|rows| rows.into_iter().take(count).collect())
            .map_err(|error| {
                RecurrenceError::VerificationFailed(format!(
                    "the image rule cannot be run: {error}"
                ))
            })
    }
}

// ---------------------------------------------------------------------------
// Parent
// ---------------------------------------------------------------------------

/// A validated parent: its relation `L` (base form), first index, seeds and
/// the first base from which `L` holds.
struct Parent {
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
    relation: DOp,
    first_base: usize,
}

impl Parent {
    fn new(json: &RecurrenceJson) -> Result<Self, RecurrenceError> {
        let (rule, first_index, rows) = json.to_recurrence_parts()?;
        if rule.terms.is_empty() {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has no terms".into(),
            ));
        }
        if rule.terms.iter().any(|t| t.offset == 0) {
            return Err(RecurrenceError::InvalidInput(
                "parent recurrence offsets must be positive".into(),
            ));
        }
        if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has a (-1)^n factor".into(),
            ));
        }
        if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
            return Err(RecurrenceError::Unsupported(
                "the parent recurrence has an inhomogeneous term".into(),
            ));
        }
        let leader = rule.denominator.clone().unwrap_or(BivarPoly {
            coeffs: vec![vec![Q::one()]],
        });
        if leader.is_zero() {
            return Err(RecurrenceError::InvalidInput(
                "the parent recurrence denominator is zero".into(),
            ));
        }
        if !Biv::from_bivar(&leader).is_y_free() {
            return Err(RecurrenceError::Unsupported(
                "the parent leader f(n, x) must be free of x".into(),
            ));
        }
        let order = rule.max_offset();
        if rows.len() < order {
            return Err(RecurrenceError::InvalidInput(format!(
                "the parent recurrence has order {order} but only {} initial rows",
                rows.len()
            )));
        }
        let count = rule.generation_initial_count(first_index, rows.len());
        if count > MAX_IMAGE_RECURRENCE_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the parent needs {count} seeds, more than {MAX_IMAGE_RECURRENCE_ROWS}"
            )));
        }
        let rows: Vec<Vec<Q>> = rows.iter().map(|row| vec_poly::trim(row)).collect();
        let seeds = rows[..count].to_vec();
        let generated = rule
            .generate_rows_rational(&seeds, first_index, rows.len())
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!(
                    "the parent recurrence cannot be replayed: {error}"
                ))
            })?;
        for (offset, (supplied, generated)) in rows.iter().zip(&generated).enumerate() {
            if supplied != &vec_poly::trim(generated) {
                return Err(RecurrenceError::InvalidInput(format!(
                    "supplied parent row {} disagrees with its recurrence",
                    first_index + offset
                )));
            }
        }
        if !bivariate_leader_nonvanishing(&leader, first_index + count) {
            return Err(RecurrenceError::Unsupported(format!(
                "the parent leader is not proved nonzero at every target index >= {}",
                first_index + count
            )));
        }
        let (relation, _) = ImplicitRelation::from_recurrence(&rule)?;
        Ok(Self {
            first_base: first_index + count - order,
            rule,
            first_index,
            seeds,
            relation: relation.op(),
        })
    }

    /// Rows `P_(start), ..., P_(start + count - 1)`, `start >= p0`.
    fn rows(&self, start: usize, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let offset = start - self.first_index;
        let total = (offset + count).max(self.seeds.len());
        let rows = self
            .rule
            .generate_rows_rational(&self.seeds, self.first_index, total)
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!("the parent cannot be replayed: {error}"))
            })?;
        Ok(rows[offset..offset + count]
            .iter()
            .map(|row| vec_poly::trim(row))
            .collect())
    }
}

// ---------------------------------------------------------------------------
// Specification
// ---------------------------------------------------------------------------

fn parse_value(text: &str, what: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(text).map_err(|error| RecurrenceError::from(error).with_context(what))
}

fn parse_polynomial(values: &[String], what: &str) -> Result<Vec<Q>, RecurrenceError> {
    if values.len() > 64 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {what} has more than 64 coefficients"
        )));
    }
    let parsed = values
        .iter()
        .map(|value| parse_value(value, what))
        .collect::<Result<Vec<_>, _>>()?;
    let trimmed = vec_poly::trim(&parsed);
    if vec_poly::is_zero(&trimmed) {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {what} is the zero polynomial"
        )));
    }
    Ok(trimmed)
}

/// `p(n) != 0` for every integer `n >= from`, decided by the sign checker.
fn require_nonzero_from(p: &[Q], from: usize, what: &str) -> Result<(), RecurrenceError> {
    match decide_index_polynomial_sign(
        p,
        &IntegerRange::starting_at(from),
        IndexSignCondition::Nonzero,
    ) {
        IndexSignDecision::Holds { .. } => Ok(()),
        IndexSignDecision::FailsAt { index, .. } => Err(RecurrenceError::InvalidInput(format!(
            "the {what} vanishes at n = {index}; it must be nonzero for every n >= {from}"
        ))),
        IndexSignDecision::Undecided { reason } => Err(RecurrenceError::BudgetExhausted(format!(
            "the {what} could not be decided nonzero for every n >= {from}: {reason}"
        ))),
    }
}

/// `p(k) != 0` for every integer `k >= from` (`from` may be negative),
/// decided by exact integer roots.
fn require_no_integer_roots_from(p: &[Q], from: i64, what: &str) -> Result<(), RecurrenceError> {
    let shifted = vec_poly::compose(p, &[Q::from_integer(from.into()), Q::one()]);
    if let Some(root) = roots_in_domain(&shifted, 0, None).first() {
        return Err(RecurrenceError::Unsupported(format!(
            "the {what} vanishes at k = {}; it must be nonzero for every integer k >= {from}",
            *root as i64 + from
        )));
    }
    Ok(())
}

/// The specification as crate steps with its gauge values.
struct PreparedImage {
    shift: i64,
    /// The map steps (without the gauge).
    map_steps: Vec<RowMapStep>,
    /// `mu(n)/mu(n-1) = alpha(n)/beta(n)`, `mu = lambda / delta`.
    alpha: Vec<Q>,
    beta: Vec<Q>,
    /// `mu(n0)`.
    initial_gauge: Q,
    /// The coefficient multiplier `(A, B)`, if any.
    multiplier: Option<(Vec<Q>, Vec<Q>)>,
}

impl PreparedImage {
    fn new(
        specification: &ImageOperatorSpecification,
        parent: &Parent,
        n0: usize,
    ) -> Result<Self, RecurrenceError> {
        let shift = specification.index_shift;
        if shift.unsigned_abs() > MAX_IMAGE_RECURRENCE_INDEX_SHIFT as u64 {
            return Err(RecurrenceError::LimitExceeded(format!(
                "|index_shift| exceeds {MAX_IMAGE_RECURRENCE_INDEX_SHIFT}"
            )));
        }
        let start = i64::try_from(n0)
            .ok()
            .and_then(|n0| n0.checked_add(shift))
            .and_then(|start| usize::try_from(start).ok());
        if start.is_none_or(|start| start < parent.first_index) {
            return Err(RecurrenceError::InvalidInput(format!(
                "the image uses parent row n0 + s = {n0} + {shift}, before the parent's first \
                 index {}",
                parent.first_index
            )));
        }
        let (map_steps, delta, multiplier) = map_steps(&specification.map, n0)?;
        require_nonzero_from(&delta, n0, "Weyl-operator denominator")?;
        let at = |p: &[Q], n: usize| vec_poly::evaluate(p, &Q::from_integer(n.into()));
        // lambda in ratio form.
        let (a, b, lambda0) = match &specification.gauge {
            ImageGauge::None => (vec![Q::one()], vec![Q::one()], Q::one()),
            ImageGauge::Ratio {
                numerator,
                denominator,
                initial_value,
            } => {
                let a = parse_polynomial(numerator, "gauge numerator")?;
                let b = parse_polynomial(denominator, "gauge denominator")?;
                let initial = parse_value(initial_value, "gauge initial value")?;
                if initial.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "the gauge initial value lambda(n0) must be nonzero".into(),
                    ));
                }
                require_nonzero_from(&a, n0 + 1, "gauge ratio numerator")?;
                require_nonzero_from(&b, n0 + 1, "gauge ratio denominator")?;
                (a, b, initial)
            }
            ImageGauge::Closed {
                numerator,
                denominator,
            } => {
                let p = parse_polynomial(numerator, "gauge numerator")?;
                let q = parse_polynomial(denominator, "gauge denominator")?;
                require_nonzero_from(&p, n0, "gauge numerator")?;
                require_nonzero_from(&q, n0, "gauge denominator")?;
                let previous = |f: &[Q]| vec_poly::compose(f, &[-Q::one(), Q::one()]);
                (
                    vec_poly::mul(&p, &previous(&q)),
                    vec_poly::mul(&q, &previous(&p)),
                    at(&p, n0) / at(&q, n0),
                )
            }
        };
        // mu = lambda / delta: alpha = a(n) delta(n - 1), beta = b(n) delta(n).
        let previous_delta = vec_poly::compose(&delta, &[-Q::one(), Q::one()]);
        Ok(Self {
            shift,
            map_steps,
            alpha: vec_poly::trim(&vec_poly::mul(&a, &previous_delta)),
            beta: vec_poly::trim(&vec_poly::mul(&b, &delta)),
            initial_gauge: lambda0 / at(&delta, n0),
            multiplier,
        })
    }

    /// Every step, the gauge last (omitted when `mu` is constant).
    fn steps(&self) -> Vec<RowMapStep> {
        let mut steps = self.map_steps.clone();
        if self.alpha != self.beta {
            steps.push(RowMapStep::Gauge {
                numerator: self.alpha.clone(),
                denominator: self.beta.clone(),
            });
        }
        steps
    }

    /// The rows `Q_n0, ..., Q_(n0 + count - 1)`.
    fn image_rows(
        &self,
        parent: &Parent,
        n0: usize,
        count: usize,
    ) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let start = (n0 as i64 + self.shift) as usize;
        let parent_rows = parent.rows(start, count)?;
        let mut mu = self.initial_gauge.clone();
        let mut out = Vec::with_capacity(count);
        for (offset, row) in parent_rows.iter().enumerate() {
            let n = n0 + offset;
            if offset > 0 {
                let index = Q::from_integer(n.into());
                mu = mu * vec_poly::evaluate(&self.alpha, &index)
                    / vec_poly::evaluate(&self.beta, &index);
            }
            let mut current = row.clone();
            for step in &self.map_steps {
                current = step.apply_row(&current, n).ok_or_else(|| {
                    RecurrenceError::Unsupported(format!(
                        "the row map is undefined at n = {n} (for a reversal: the parent row \
                         has degree above d_n)"
                    ))
                })?;
            }
            out.push(vec_poly::trim(&vec_poly::scale(&current, &mu)));
        }
        Ok(out)
    }

    /// The shifted parent rows `P'_n0, ...`.
    fn parent_rows(
        &self,
        parent: &Parent,
        n0: usize,
        count: usize,
    ) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        parent.rows((n0 as i64 + self.shift) as usize, count)
    }
}

/// The map steps, the Weyl denominator `delta` and the multiplier `(A, B)`.
#[allow(clippy::type_complexity)]
fn map_steps(
    map: &ImageRowMap,
    n0: usize,
) -> Result<(Vec<RowMapStep>, Vec<Q>, Option<(Vec<Q>, Vec<Q>)>), RecurrenceError> {
    let reversal_degree = |schedule: &DegreeSchedule| -> Result<Vec<Q>, RecurrenceError> {
        schedule.validate()?;
        if schedule.last_index.is_some() || schedule.step < 0 || schedule.first_index > n0 {
            return Err(RecurrenceError::InvalidInput(
                "the reversal schedule must be unbounded, with a nonnegative step, and start at \
                 or before n0"
                    .into(),
            ));
        }
        let step = Q::from_integer(schedule.step.into());
        let first = Q::from_integer(schedule.first_index.into());
        let degree = Q::from_integer(schedule.first_degree.into());
        Ok(vec_poly::trim(&[degree - &step * first, step]))
    };
    let one = vec![Q::one()];
    Ok(match map {
        ImageRowMap::Identity => (Vec::new(), one, None),
        ImageRowMap::WeylOperator { terms, denominator } => {
            if terms.is_empty() || terms.len() > 64 {
                return Err(RecurrenceError::InvalidInput(
                    "the Weyl operator needs between 1 and 64 terms".into(),
                ));
            }
            let mut op = DOp::default();
            for term in terms {
                if term.deriv_order > 16 {
                    return Err(RecurrenceError::LimitExceeded(
                        "a Weyl-operator derivative order exceeds 16".into(),
                    ));
                }
                op.add_term((0, term.deriv_order), biv_from_json(&term.coeff)?);
            }
            if op.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "the Weyl operator is zero".into(),
                ));
            }
            let delta = parse_polynomial(denominator, "Weyl-operator denominator")?;
            (vec![RowMapStep::WeylOperator(op)], delta, None)
        }
        ImageRowMap::Affine { slope, intercept } => {
            let slope = parse_value(slope, "affine slope")?;
            if slope.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "the affine slope must be nonzero".into(),
                ));
            }
            let intercept = parse_value(intercept, "affine intercept")?;
            (vec![RowMapStep::Affine { slope, intercept }], one, None)
        }
        ImageRowMap::Reversal { schedule } => (
            vec![RowMapStep::Reversal {
                degree: reversal_degree(schedule)?,
            }],
            one,
            None,
        ),
        ImageRowMap::Mobius {
            a,
            b,
            c,
            e,
            schedule,
        } => {
            let a = parse_value(a, "Mobius a")?;
            let b = parse_value(b, "Mobius b")?;
            let c = parse_value(c, "Mobius c")?;
            let e = parse_value(e, "Mobius e")?;
            if c.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "a Mobius map needs c != 0 (use the affine map otherwise)".into(),
                ));
            }
            let determinant = &a * &e - &b * &c;
            if determinant.is_zero() {
                return Err(RecurrenceError::InvalidInput(
                    "a Mobius map needs a e - b c != 0".into(),
                ));
            }
            // (c x + e)^d P((a x + b)/(c x + e)) = A2(R(A1(P))) with
            // A1: F(y) -> F(a/c + kappa y), kappa = (b c - a e)/c,
            // R: G(u) -> u^d G(1/u), A2: K(u) -> K(c x + e).
            let kappa = -determinant / &c;
            (
                vec![
                    RowMapStep::Affine {
                        slope: kappa,
                        intercept: a / &c,
                    },
                    RowMapStep::Reversal {
                        degree: reversal_degree(schedule)?,
                    },
                    RowMapStep::Affine {
                        slope: c,
                        intercept: e,
                    },
                ],
                one,
                None,
            )
        }
        ImageRowMap::CoefficientMultiplier {
            numerator,
            denominator,
        } => {
            let a = parse_polynomial(numerator, "multiplier numerator")?;
            let b = parse_polynomial(denominator, "multiplier denominator")?;
            require_no_integer_roots_from(&a, 0, "multiplier numerator A(k)")?;
            require_no_integer_roots_from(&b, 0, "multiplier denominator B(k)")?;
            (
                vec![RowMapStep::CoefficientMultiplier {
                    numerator: a.clone(),
                    denominator: b.clone(),
                }],
                one,
                Some((a, b)),
            )
        }
    })
}

// ---------------------------------------------------------------------------
// JSON of operators
// ---------------------------------------------------------------------------

fn operator_json(op: &DOp) -> ImplicitRelationJson {
    ImplicitRelation::from_op(op).to_json()
}

fn operator_from_json(json: &ImplicitRelationJson) -> Result<DOp, RecurrenceError> {
    Ok(ImplicitRelation::from_json(json)?.op())
}

/// The runnable rule of `M` (target form) and its order.
fn successor_rule(m: &DOp) -> Result<(Recurrence, usize), RecurrenceError> {
    ImplicitRelation::from_op(m)
        .to_successor_rule()
        .ok_or_else(|| RecurrenceError::VerificationFailed("the image rule is not runnable".into()))
}

// ---------------------------------------------------------------------------
// The construction and its checks
// ---------------------------------------------------------------------------

/// Bases `n >= from` where `p(n, .)` is the zero polynomial.
fn zero_bases(p: &Biv, from: usize) -> Vec<usize> {
    let (content, _) = content_and_primitive(p);
    roots_in_domain(&content, from, None)
}

/// The bases `m` where each member must be checked (module documentation,
/// step 2), image indexing; `first_base` is the first base of `L'`.
fn system_checks(
    certificate: &AnnihilatorSystemCertificate,
    n0: usize,
    first_base: usize,
) -> Vec<Vec<usize>> {
    let j = certificate.step;
    let start = n0.max(first_base);
    certificate
        .entries
        .iter()
        .map(|entry| {
            let mut bases: Vec<usize> = (n0..start + j).collect();
            for base in zero_bases(&entry.multiplier, n0) {
                bases.push(base + j);
            }
            bases.sort_unstable();
            bases.dedup();
            bases
        })
        .collect()
}

/// Check `(C P')_m = 0` for every member at its bases, exactly.
fn system_holds(
    certificate: &AnnihilatorSystemCertificate,
    checks: &[Vec<usize>],
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
) -> Result<bool, RecurrenceError> {
    let needed = certificate
        .members
        .iter()
        .zip(checks)
        .filter_map(|(c, bases)| bases.last().map(|m| m + c.max_shift() + 1 - n0))
        .max()
        .unwrap_or(0);
    require_row_budget(needed)?;
    let rows = prepared.parent_rows(parent, n0, needed)?;
    Ok(certificate.members.iter().zip(checks).all(|(c, bases)| {
        bases.iter().all(|m| {
            c.apply(&rows, n0, *m)
                .is_some_and(|value| vec_poly::is_zero(&value))
        })
    }))
}

fn system_json(
    certificate: &AnnihilatorSystemCertificate,
    checks: Vec<Vec<usize>>,
) -> ImageAnnihilatorSystem {
    ImageAnnihilatorSystem {
        step: certificate.step,
        members: certificate
            .members
            .iter()
            .zip(&certificate.entries)
            .zip(checks)
            .map(|((c, entry), bases)| ImageAnnihilatorSystemMember {
                annihilator: operator_json(c),
                multiplier: biv_json(&entry.multiplier),
                recurrence_quotient: operator_json(&entry.recurrence_quotient),
                member_quotients: entry.member_quotients.iter().map(operator_json).collect(),
                checked_bases: bases,
            })
            .collect(),
    }
}

/// Replay a recorded annihilator system relative to `L'`.
fn replay_system(
    json: &ImageAnnihilatorSystem,
    recurrence: &DOp,
    first_base: usize,
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
) -> Result<AnnihilatorSystemCertificate, RecurrenceError> {
    let members = json
        .members
        .iter()
        .map(|m| operator_from_json(&m.annihilator))
        .collect::<Result<Vec<_>, _>>()?;
    let entries = json
        .members
        .iter()
        .map(|m| {
            Ok(AnnihilatorSystemEntry {
                multiplier: biv_from_json(&m.multiplier)?,
                recurrence_quotient: operator_from_json(&m.recurrence_quotient)?,
                member_quotients: m
                    .member_quotients
                    .iter()
                    .map(operator_from_json)
                    .collect::<Result<Vec<_>, RecurrenceError>>()?,
            })
        })
        .collect::<Result<Vec<_>, RecurrenceError>>()?;
    let certificate = AnnihilatorSystemCertificate {
        recurrence: recurrence.clone(),
        members,
        step: json.step,
        entries,
    };
    certificate.verify().map_err(|error| {
        RecurrenceError::VerificationFailed(format!("annihilator system: {error}"))
    })?;
    if certificate
        .members
        .iter()
        .filter(|c| c.max_shift() > 0)
        .count()
        > 1
        || certificate
            .members
            .iter()
            .filter(|c| c.max_shift() == 0)
            .count()
            > 1
    {
        return Err(RecurrenceError::VerificationFailed(
            "the annihilator system has at most one member of each kind".into(),
        ));
    }
    let checks = system_checks(&certificate, n0, first_base);
    let recorded: Vec<Vec<usize>> = json
        .members
        .iter()
        .map(|m| m.checked_bases.clone())
        .collect();
    if checks != recorded {
        return Err(RecurrenceError::VerificationFailed(format!(
            "annihilator system: recorded checks {recorded:?}, recomputed {checks:?}"
        )));
    }
    if !system_holds(&certificate, &checks, prepared, parent, n0)? {
        return Err(RecurrenceError::VerificationFailed(
            "annihilator system: a base row check fails".into(),
        ));
    }
    Ok(certificate)
}

/// The relation in the identity and its first base, and the row
/// annihilator, from an annihilator system (module documentation).
fn working_relation(
    system: Option<&AnnihilatorSystemCertificate>,
    relation: &DOp,
    first_base: usize,
    n0: usize,
) -> (DOp, usize, Option<DOp>) {
    match system {
        None => (relation.clone(), first_base, None),
        Some(system) => {
            let row = system.row_annihilator().cloned();
            match system.lower_order_member() {
                Some(b) => (b.clone(), n0, row),
                None => (relation.clone(), first_base, row),
            }
        }
    }
}

/// The seed count of a runnable `M` with first index `n0`: at least `R`,
/// and past every target where the leader vanishes.
fn seed_count(m: &DOp, n0: usize) -> usize {
    let order = m.max_shift();
    let leader = m.0.get(&(order, 0)).cloned().unwrap_or_default();
    zero_bases(&leader, n0)
        .into_iter()
        .map(|base| base + order + 1 - n0)
        .max()
        .unwrap_or(0)
        .max(order)
}

/// The exceptional bases of the identity (module documentation, step 3).
fn exceptional_bases(
    multiplier: &Biv,
    base_factor: &[Q],
    n0: usize,
    first_base: usize,
) -> Vec<usize> {
    let mut bases: Vec<usize> = (n0..first_base.max(n0)).collect();
    bases.extend(zero_bases(multiplier, n0));
    bases.extend(roots_in_domain(base_factor, n0, None));
    bases.sort_unstable();
    bases.dedup();
    bases
}

fn require_row_budget(count: usize) -> Result<(), RecurrenceError> {
    if count > MAX_IMAGE_RECURRENCE_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the image certificate needs {count} rows, more than {MAX_IMAGE_RECURRENCE_ROWS}"
        )));
    }
    Ok(())
}

/// The multiplier side condition on `[-S, -1]` (module documentation).
fn check_multiplier_bound(
    prepared: &PreparedImage,
    bounds: PullbackBounds,
) -> Result<(), RecurrenceError> {
    if let Some((a, _)) = &prepared.multiplier {
        let s = bounds.multiplier_euler_bound() as i64;
        require_no_integer_roots_from(a, -s, "multiplier numerator A(k)")?;
    }
    Ok(())
}

/// The shifted parent relation `L'` and its first base `b_L'` (saturated
/// at zero).
fn shifted_relation(parent: &Parent, shift: i64) -> (DOp, usize) {
    let first_base = (parent.first_base as i64 - shift).max(0) as usize;
    (parent.relation.shift_index(shift), first_base)
}

/// The row checks of a candidate at its exceptional bases: returns
/// `(image rows, exceptional bases)` or the failure.
fn check_rows(
    certificate: &ImageAnnihilatorCertificate,
    base_factor: &[Q],
    prepared: &PreparedImage,
    parent: &Parent,
    n0: usize,
    first_base: usize,
    replay_rows: usize,
) -> Result<(Vec<Vec<Q>>, Vec<usize>), RecurrenceError> {
    let m = &certificate.annihilator;
    let order = m.max_shift();
    let seeds = seed_count(m, n0);
    let exceptional = exceptional_bases(&certificate.multiplier, base_factor, n0, first_base);
    let needed = (seeds + replay_rows).max(exceptional.last().map_or(0, |b| b + order + 1 - n0));
    require_row_budget(needed)?;
    let images = prepared.image_rows(parent, n0, needed)?;
    for base in &exceptional {
        let value = m.apply(&images, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at base {base}"))
        })?;
        if !vec_poly::is_zero(&value) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the image rule fails at the exceptional base {base}"
            )));
        }
    }
    Ok((images, exceptional))
}

/// The runnable rule of `Q_n = lambda(n) T_n(P_(n+s))`, `n >= first_index`,
/// with seeds and a certificate (module documentation).
///
/// Returns [`RecurrenceError::BudgetExhausted`] when no runnable rule is
/// found within the budget, [`RecurrenceError::Unsupported`] for a parent or
/// map outside the method (a parent leader involving `x`, a side condition
/// that fails) and [`RecurrenceError::InvalidInput`] for malformed input.
pub fn image_recurrence_of_operator(
    parent: &RecurrenceJson,
    first_index: usize,
    specification: &ImageOperatorSpecification,
    options: &ImageRecurrenceOptions,
) -> Result<ImageRecurrence, RecurrenceError> {
    if options.max_unknowns > MAX_IMAGE_RECURRENCE_UNKNOWNS
        || options.max_shapes > MAX_IMAGE_RECURRENCE_SHAPES
        || options.replay_rows > MAX_IMAGE_RECURRENCE_ROWS
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "max_unknowns <= {MAX_IMAGE_RECURRENCE_UNKNOWNS}, max_shapes <= \
             {MAX_IMAGE_RECURRENCE_SHAPES} and replay_rows <= {MAX_IMAGE_RECURRENCE_ROWS} are \
             required"
        )));
    }
    let n0 = first_index;
    let validated = Parent::new(parent)?;
    let prepared = PreparedImage::new(specification, &validated, n0)?;
    let (relation, first_base) = shifted_relation(&validated, prepared.shift);
    let steps = prepared.steps();
    let budget = ImageAnnihilatorBudget {
        max_order_increase: options.max_order_increase,
        max_derivative_order: options.max_derivative_order,
        max_index_degree: options.max_coefficient_degree,
        max_variable_degree: options.max_coefficient_degree,
        max_unknowns: options.max_unknowns,
        max_shapes: options.max_shapes,
    };
    let crate_error =
        |error: ore_algebra::ImageAnnihilatorError| RecurrenceError::Unsupported(error.to_string());
    // Additional annihilators proved from the rows: a runnable rule of lower
    // order replaces L' (it holds from n0 on and determines the same rows),
    // and a row annihilator enlarges the ideal (every rule found without it
    // is found with it).
    let mut system = None;
    if options.row_annihilators {
        let fit_rows = prepared.parent_rows(&validated, n0, AUXILIARY_FIT_ROWS)?;
        if let Some(candidate) = find_annihilator_system(
            &relation,
            &fit_rows,
            n0,
            &AnnihilatorSystemBudget::default(),
        ) {
            let checks = system_checks(&candidate, n0, first_base);
            if system_holds(&candidate, &checks, &prepared, &validated, n0)? {
                system = Some((candidate, checks));
            }
        }
    }
    let (working, working_base, row) =
        working_relation(system.as_ref().map(|(c, _)| c), &relation, first_base, n0);
    let search =
        find_image_annihilator(&working, row.as_ref(), &steps, &budget).map_err(crate_error)?;
    let (certificate, method) = match search {
        ImageAnnihilatorSearch::Found {
            certificate,
            method,
        } => (*certificate, method),
        ImageAnnihilatorSearch::NotFoundWithinBudget {
            shapes_tried,
            largest_shape,
        } => {
            return Err(RecurrenceError::BudgetExhausted(format!(
                "no runnable image rule within the budget ({shapes_tried} ansatz shapes solved, \
                 the largest {largest_shape:?}{})",
                if system.is_some() {
                    ", also with an annihilator system of the parent rows"
                } else {
                    ""
                }
            )))
        }
    };
    // The system is recorded when the identity uses one of its members (the
    // row annihilator, or the lower-order member as the relation).
    if certificate.row_annihilator.is_none()
        && system
            .as_ref()
            .is_some_and(|(c, _)| c.lower_order_member().is_none())
    {
        system = None;
    }
    check_multiplier_bound(&prepared, certificate.bounds)?;
    let base_factor = certificate
        .verify()
        .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
    let (images, exceptional) = check_rows(
        &certificate,
        &base_factor,
        &prepared,
        &validated,
        n0,
        working_base,
        options.replay_rows,
    )?;
    let m = &certificate.annihilator;
    let (rule, order) = successor_rule(m)?;
    let seeds = images[..seed_count(m, n0)].to_vec();
    let rule = RecurrenceJson::from_recurrence_rational(&rule, n0, &seeds, None);
    let image = ImageRecurrence {
        schema: IMAGE_RECURRENCE_SCHEMA.into(),
        rule: rule.clone(),
        order,
        method: match method {
            ImageAnnihilatorMethod::PushForward => ImageRecurrenceMethod::PushForward,
            ImageAnnihilatorMethod::LeftMultiple(shape) => ImageRecurrenceMethod::LeftMultiple {
                order: shape.order,
                derivative_order: shape.derivative_order,
                index_degree: shape.index_degree,
                variable_degree: shape.variable_degree,
            },
        },
        certificate: ImageRecurrenceCertificate {
            parent: parent.clone(),
            first_index: n0,
            specification: specification.clone(),
            rule,
            bounds: ImageRecurrenceBounds {
                top_shift: certificate.bounds.top_shift,
                derivative_order: certificate.bounds.derivative_order,
                variable_degree: certificate.bounds.variable_degree,
            },
            multiplier: biv_json(&certificate.multiplier),
            quotient: operator_json(&certificate.quotient),
            annihilator_system: system.map(|(c, checks)| system_json(&c, checks)),
            row_quotient: operator_json(&certificate.row_quotient),
            exceptional_bases: exceptional,
            replay_rows: options.replay_rows,
        },
    };
    image.verify()?;
    Ok(image)
}

impl ImageRecurrenceCertificate {
    /// Recompute everything (module documentation): the parent and its
    /// rows, the steps and side conditions, the identity
    /// `w E = H o L' + Z o A` and the invariance of `A`, the row checks at
    /// the exceptional bases and of `A`, the seeds as `lambda(n) T_n(P'_n)`,
    /// the leader's nonvanishing after the seeds, and the replay.  Returns
    /// the order of the rule.
    pub fn verify(&self) -> Result<usize, RecurrenceError> {
        let n0 = self.first_index;
        if self.replay_rows > MAX_IMAGE_RECURRENCE_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "replay_rows exceeds {MAX_IMAGE_RECURRENCE_ROWS}"
            )));
        }
        let parent = Parent::new(&self.parent)?;
        let prepared = PreparedImage::new(&self.specification, &parent, n0)?;
        let (relation, first_base) = shifted_relation(&parent, prepared.shift);
        let (rule, rule_first, seeds) = self.rule.to_recurrence_parts()?;
        if rule_first != n0 {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the rule starts at {rule_first}, not at n0 = {n0}"
            )));
        }
        let (m, order) = ImplicitRelation::from_recurrence(&rule)?;
        let m = m.op();
        let failed = |what: &str, error: &dyn std::fmt::Display| {
            RecurrenceError::VerificationFailed(format!("{what}: {error}"))
        };
        let system = self
            .annihilator_system
            .as_ref()
            .map(|json| replay_system(json, &relation, first_base, &prepared, &parent, n0))
            .transpose()?;
        let (working, working_base, row) =
            working_relation(system.as_ref(), &relation, first_base, n0);
        let bounds = PullbackBounds {
            top_shift: self.bounds.top_shift,
            derivative_order: self.bounds.derivative_order,
            variable_degree: self.bounds.variable_degree,
        };
        check_multiplier_bound(&prepared, bounds)?;
        let certificate = ImageAnnihilatorCertificate {
            annihilator: m.clone(),
            parent_annihilator: working,
            steps: prepared.steps(),
            bounds,
            row_annihilator: row,
            multiplier: biv_from_json(&self.multiplier)?,
            quotient: operator_from_json(&self.quotient)?,
            row_quotient: operator_from_json(&self.row_quotient)?,
        };
        let base_factor = certificate
            .verify()
            .map_err(|error| failed("image identity", &error))?;
        let (images, exceptional) = check_rows(
            &certificate,
            &base_factor,
            &prepared,
            &parent,
            n0,
            working_base,
            self.replay_rows,
        )?;
        if exceptional != self.exceptional_bases {
            return Err(RecurrenceError::VerificationFailed(format!(
                "recorded exceptional bases {:?}, recomputed {exceptional:?}",
                self.exceptional_bases
            )));
        }
        let count = seed_count(&m, n0);
        let seeds: Vec<Vec<Q>> = seeds.iter().map(|row| vec_poly::trim(row)).collect();
        if seeds.as_slice() != &images[..count] {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the rule's {} seeds are not the {count} image rows lambda(n) T_n(P'_n)",
                seeds.len()
            )));
        }
        let leader = rule.denominator.clone().unwrap_or(BivarPoly {
            coeffs: vec![vec![Q::one()]],
        });
        if !bivariate_leader_nonvanishing(&leader, n0 + count) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the leader is not proved nonzero at every target >= {}",
                n0 + count
            )));
        }
        let total = count + self.replay_rows;
        let generated = rule
            .generate_rows_rational(&seeds, n0, total)
            .map_err(|error| failed("replay", &error))?;
        for (offset, (generated, image)) in generated.iter().zip(&images).enumerate() {
            if &vec_poly::trim(generated) != image {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "the rule's row {} differs from lambda(n) T_n(P'_n)",
                    n0 + offset
                )));
            }
        }
        Ok(order)
    }
}
