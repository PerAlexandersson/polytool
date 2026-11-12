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
    content_and_primitive, q_usize, roots_in_domain, Biv, BivarConversion, DOp,
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

mod preparation;
use preparation::*;
mod replay;
pub(in crate::recurrence) use replay::runnable_record_from_annihilator;
use replay::*;
