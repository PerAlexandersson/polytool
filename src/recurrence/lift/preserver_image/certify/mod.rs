//! The candidate, the certificates and the certification of one
//! preserver-image lift `Q_n = lambda(n) T_n(P_(n+s))` (module
//! documentation of `preserver_image`).

use super::interlacing::{interlacing_claim, PreserverInterlacingClaim};
use super::interval::{find_bounded_root_region_report, wants_interval_parent};
use super::operator::{dividend, OperatorForm, PreserverOperator};
use super::summary::{preserver_image_summary, PreserverImageSummaryParts};
use super::transfer::operator_transfer;
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::recurrence::implicit::algebra::{roots_in_domain, Biv, DOp};
use crate::recurrence::implicit::{biv_from_json, biv_json};
use crate::recurrence::lift::hadamard_factor::annihilator::vanishing_bases;
use crate::recurrence::lift::hadamard_factor::certify::{
    identity_holds, system_from_json, system_json, verify_invariance, verify_system, ProductParent,
    RowAnnihilation, RowAnnihilationKind, RowAnnihilatorInvariance,
};
use crate::recurrence::lift::hadamard_factor::prepare::{
    prepare_lift_rule, shift_free_from_json, shift_free_json, PreparedLiftRule,
};
use crate::recurrence::lift::hadamard_factor::search::{certify_candidate, remove_remainder};
use crate::recurrence::lift::parent_proof::{prove_parent_with_policy, ParentRecurrenceRecord};
use crate::recurrence::lift::section_certify::{
    n_coefficients, n_polynomial, pseudo_divide_scalar, x_free_leader, RatioGauge,
};
use crate::recurrence::{
    check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem,
    check_three_term_recurrence_theorem, decode_poly, encode_poly, image_recurrence_of_operator,
    BivarPolyJson, HadamardAnnihilatorSystem, HadamardFactor, HadamardFactorCertificate,
    HadamardFactorIdentityStatus, HadamardFactorLiftOptions, ImageGauge,
    ImageOperatorSpecification, ImageRecurrence, ImageRecurrenceOptions, ImageRowMap,
    ImplicitRelation, ImplicitRelationJson, LiftSummary, OperatorTermJson,
    RealRootednessTheoremOptions, RealRootednessTheoremReport, Recurrence, RecurrenceError,
    RecurrenceJson, VeroneseParentAttestation, VeroneseParentProof, VeroneseParentProofRoute,
    VeroneseParentProofStatus, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// The pairs `(A, Z)` of the terms `Z o A` of an annihilation.
type AnnihilatorTerms = Vec<(DOp, DOp)>;

/// JSON schema tag of a serialized [`PreserverImageCertification`].
pub const PRESERVER_IMAGE_CERTIFICATION_SCHEMA: &str = "polytool.recurrence-preserver-image.v1";

/// Largest `|s|` of the index shift of a candidate.
pub const MAX_PRESERVER_IMAGE_INDEX_SHIFT: usize = 8;

/// Largest number of rows compared exactly (also the bound on the rows that
/// the exceptional bases of a certificate may require).
pub const MAX_PRESERVER_IMAGE_CHECKED_ROWS: usize = 64;

/// Largest number of operator terms in a pseudo-division.
pub const MAX_PRESERVER_IMAGE_OPERATOR_TERMS: usize = 20_000;

/// Largest number of terms, derivative order, degree in `x` and coefficient
/// bit length of a division remainder handed to the row-annihilator search
/// (whose exact right gcds of operators grow quickly with these); a larger
/// remainder leaves the identity supported by rows.  A deterministic budget.
const MAX_ANNIHILATOR_REMAINDER: (usize, usize, usize, u64) = (32, 4, 8, 16);

/// Whether a remainder is within [`MAX_ANNIHILATOR_REMAINDER`].
fn remainder_is_small(remainder: &DOp) -> bool {
    let (max_terms, max_order, max_degree, max_bits) = MAX_ANNIHILATOR_REMAINDER;
    remainder.0.len() <= max_terms
        && remainder.max_derivative() <= max_order
        && remainder.0.values().all(|c| {
            c.variable_degree().unwrap_or(0) <= max_degree
                && c.rows()
                    .iter()
                    .flatten()
                    .all(|v| v.numer().bits() <= max_bits && v.denom().bits() <= max_bits)
        })
}

/// Largest distance `n0 + s - p0` between the first used parent row and the
/// parent's first index.
const MAX_PARENT_OFFSET: usize = 256;

/// A candidate `Q_n = lambda(n) T_n(P_(n+s))` for every `n >= n0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverImageCandidate {
    /// A short name of the parent (for reports).
    #[serde(default)]
    pub parent_name: String,
    /// The parent rule with its first index `p0` and seeds
    /// (`polytool.recurrence.v1`): homogeneous, without `(-1)^n` factors,
    /// with a leader `g(n)` free of `x` and of derivatives.
    pub parent: RecurrenceJson,
    /// The index shift `s`, `|s| <=` [`MAX_PRESERVER_IMAGE_INDEX_SHIFT`] and
    /// `n0 + s >= p0`.  Default 0.
    #[serde(default)]
    pub index_shift: i64,
    /// The normalisation `lambda`.  Default none (`lambda = 1`).
    #[serde(default = "default_gauge")]
    pub gauge: VeroneseSectionGauge,
    /// The operator `T_n`.
    pub operator: PreserverOperator,
}

fn default_gauge() -> VeroneseSectionGauge {
    VeroneseSectionGauge::None
}

impl PreserverImageCandidate {
    /// The candidate as the image `Q_n = lambda(n) T_n(P_(n+s))` of
    /// [`crate::recurrence::image_recurrence_of_operator`]: the same index
    /// shift and gauge, and the operator's row map
    /// (`ImageRowMap::try_from`).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::Unsupported`] for a multiplier that depends on `n`,
    /// and the validation errors of the operator.
    pub fn image_operator_specification(
        &self,
    ) -> Result<ImageOperatorSpecification, RecurrenceError> {
        let gauge = match &self.gauge {
            VeroneseSectionGauge::None => ImageGauge::None,
            VeroneseSectionGauge::Ratio {
                numerator,
                denominator,
                initial_value,
            } => ImageGauge::Ratio {
                numerator: numerator.clone(),
                denominator: denominator.clone(),
                initial_value: initial_value.clone(),
            },
            VeroneseSectionGauge::Closed {
                numerator,
                denominator,
            } => ImageGauge::Closed {
                numerator: numerator.clone(),
                denominator: denominator.clone(),
            },
        };
        Ok(ImageOperatorSpecification {
            index_shift: self.index_shift,
            gauge,
            map: ImageRowMap::try_from(&self.operator)?,
        })
    }

    /// A runnable rule, with seeds and a certificate, for the candidate's
    /// image `lambda(n) T_n(P_(n+s))`, `n >= first_index`, computed from the
    /// parent rule ([`crate::recurrence::image_recurrence_of_operator`]),
    /// not fitted.
    ///
    /// # Errors
    ///
    /// Those of [`PreserverImageCandidate::image_operator_specification`]
    /// and of `image_recurrence_of_operator`
    /// ([`RecurrenceError::BudgetExhausted`] when no runnable rule is found
    /// within the budget).
    pub fn image_recurrence(
        &self,
        first_index: usize,
        options: &ImageRecurrenceOptions,
    ) -> Result<ImageRecurrence, RecurrenceError> {
        image_recurrence_of_operator(
            &self.parent,
            first_index,
            &self.image_operator_specification()?,
            options,
        )
    }
}

/// Options of [`certify_preserver_image_lift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PreserverImageOptions {
    /// Rows `Q_(n0), ...` compared exactly with the candidate (at least the
    /// identifying rows; more when the exceptional bases need them).  At most
    /// [`MAX_PRESERVER_IMAGE_CHECKED_ROWS`].  Default 16.
    pub checked_rows: usize,
    /// Operator-term budget, at most [`MAX_PRESERVER_IMAGE_OPERATOR_TERMS`].
    pub max_operator_terms: usize,
    /// Caller-supplied parent attestations (the attestation route of the
    /// parent provers; a lift is then conditional, never proved).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain parent route.  Default true.
    pub try_root_domain: bool,
    /// Permit caller and bundled-catalog attestations. Default true; structural discovery disables both.
    pub allow_attestations: bool,
    /// When set, a certified identity also gets the computed runnable rule
    /// of the image ([`PreserverImageCertification::image_rule`]) with these
    /// options.  Default `None`.
    pub image_rule: Option<ImageRecurrenceOptions>,
}

impl Default for PreserverImageOptions {
    fn default() -> Self {
        Self {
            checked_rows: 16,
            max_operator_terms: MAX_PRESERVER_IMAGE_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
            allow_attestations: true,
            image_rule: None,
        }
    }
}

/// One quotient `H_nu` of a [`PreserverDivisionCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverDivisionQuotient {
    /// The shift `nu` of `H_nu o S^nu o L'`.
    pub shift: usize,
    /// `H_nu = sum_d h_d(n, x) D^d`.
    pub differentiations: Vec<OperatorTermJson>,
}

/// How a nonzero division remainder `r` is removed by row annihilators of
/// the shifted parent `P'` (the `hadamard_factor` module documentation,
/// "Row annihilators"): `p(x) r = Z o A` for one annihilator `A`, or
/// `mu_2(n, x) r = sum_l Z_l o A^(l)` for an invariant system (after the
/// index content of the multiplier is folded into `w`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PreserverRowAnnihilation {
    /// One annihilator: that of a product parent (no invariance fields,
    /// recomputed from the parent by the replay), or an invariant
    /// annihilator of a first-order parent with
    /// `nu(n, x) A(n + 1) o T(n) = B(n) o A(n)`.
    Single {
        /// `p(x)`, ascending exact strings.
        row_multiplier: Vec<String>,
        /// `A(n)` (shift-free).
        annihilator: Vec<OperatorTermJson>,
        /// `Z(n)` (shift-free).
        quotient: Vec<OperatorTermJson>,
        /// `nu(n, x)` of the invariance identity (absent for a product
        /// parent).
        #[serde(default, skip_serializing_if = "Option::is_none")]
        step_multiplier: Option<BivarPolyJson>,
        /// `B(n)` of the invariance identity.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        step_quotient: Vec<OperatorTermJson>,
        /// The first base `s` of the annihilator induction.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        start_base: Option<usize>,
        /// The bases `m >= s` where `nu(m, .)` vanishes.
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        restart_bases: Vec<usize>,
    },
    /// An invariant row-annihilator system with its quotients and its
    /// multiplier `mu_2(n, x)`.
    System {
        /// The system.
        system: HadamardAnnihilatorSystem,
    },
}

/// The exact identity behind a certified lift whose operator is a Weyl
/// operator or a coordinate map (module documentation, "The certificate"):
/// `w(n) E = sum_nu H_nu S^nu L'` with remainder zero, where `E` is the
/// input relation gauged by `mu = lambda / delta` and composed with (or
/// conjugated through) the operator, and `L'` the parent relation shifted
/// to `P'_n = P_(n+s)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverDivisionCertificate {
    /// The input: rule `M`, first index `n0` and the `K` identifying seeds.
    pub input: RecurrenceJson,
    /// The parent: rule, first index `p0` and identifying seeds.
    pub parent: RecurrenceJson,
    /// The index shift `s`.
    pub index_shift: i64,
    /// The normalisation `lambda`.
    pub gauge: VeroneseSectionGauge,
    /// The operator `T_n`.
    pub operator: PreserverOperator,
    /// `L'`: the relation of `P'_n = P_(n+s)` at base index `n`; its leader
    /// `g(n + s) S^k` is free of `x`.
    pub shifted_parent_relation: ImplicitRelationJson,
    /// `k`.
    pub parent_leader_shift: usize,
    /// `E`, the dividend.
    pub dividend_relation: ImplicitRelationJson,
    /// `w(n)`, ascending exact strings in `n`.
    pub multiplier: Vec<String>,
    /// The nonzero quotients `H_nu` (multiplied by the row multiplier of
    /// the annihilation, when there is one).
    pub quotients: Vec<PreserverDivisionQuotient>,
    /// The removal of a nonzero division remainder by row annihilators of
    /// `P'`; the identity is then
    /// `p(x) mu_2(n, x) w(n) E = sum_nu H_nu S^nu L' + sum Z A`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_annihilation: Option<PreserverRowAnnihilation>,
    /// First base at which `M` holds for the input rows (`n0 + K - R`).
    pub input_valid_from_base: usize,
    /// First base at which `L'` holds for the shifted parent rows.
    pub parent_valid_from_base: usize,
    /// Bases from `input_valid_from_base` on where the identity gives
    /// nothing (before `parent_valid_from_base`, integer roots of `w`); `M`
    /// is evaluated on the candidate rows there and vanishes.
    pub exceptional_bases: Vec<usize>,
    /// The identifying input indices `n0, ..., n0 + K - 1`.
    pub identifying_indices: Vec<usize>,
    /// What the certificate proves.
    pub statement: String,
}

/// A certified lift with a hypergeometric multiplier: the
/// [`HadamardFactorCertificate`] of the Hadamard conjugation against the
/// shifted parent `P'_n = P_(n+s)`, with the records it is replayed from.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverHadamardCertificate {
    /// The input: rule `M`, first index `n0` and the identifying seeds.
    pub input: RecurrenceJson,
    /// The parent: rule, first index `p0` and identifying seeds.
    pub parent: RecurrenceJson,
    /// The index shift `s`.
    pub index_shift: i64,
    /// The factor `u`.
    pub factor: HadamardFactor,
    /// The certificate of the Hadamard conjugation; its parent is the
    /// shifted parent `P'` starting at `n0`.
    pub hadamard: HadamardFactorCertificate,
}

/// The certificate of a certified identity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "method", content = "certificate", rename_all = "snake_case")]
pub enum PreserverImageCertificate {
    /// Pseudo-division of the composed or conjugated input relation.
    OperatorDivision(Box<PreserverDivisionCertificate>),
    /// The Hadamard conjugation of a hypergeometric multiplier.
    HadamardConjugation(Box<PreserverHadamardCertificate>),
}

/// Status of the identity `Q_n = lambda(n) T_n(P_(n+s))`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum PreserverImageIdentityStatus {
    /// Proved for every `n >= n0`.
    Certified {
        /// The replayable certificate.
        certificate: Box<PreserverImageCertificate>,
    },
    /// The compared rows agree, but the identity is not proved: finite
    /// evidence only.
    SupportedByRows {
        /// Number of rows compared, from `n0`.
        rows: usize,
        /// Why the certificate does not apply.
        reason: String,
    },
    /// The identity fails exactly.
    Refuted {
        /// The absolute index of the first failing row.
        index: usize,
        /// What fails.
        reason: String,
    },
}

/// `polytool.recurrence-preserver-image.v1`: the certification of one
/// candidate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverImageCertification {
    /// Always [`PRESERVER_IMAGE_CERTIFICATION_SCHEMA`].
    pub schema: String,
    /// The candidate as supplied.
    pub candidate: PreserverImageCandidate,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number `K` of identifying input rows.
    pub identifying_rows: usize,
    /// Number of rows `Q_(n0), ...` compared exactly with the candidate.
    pub rows_compared: usize,
    /// The identity status, with the certificate when certified.
    pub identity: PreserverImageIdentityStatus,
    /// The parent statement "every parent row from `p0` on is zero or
    /// real-rooted with coefficients of one sign" and its proof attempts.
    pub parent_proof: VeroneseParentProof,
    /// When that statement is not proved and the operator's theorem needs
    /// only real-rootedness: the three-term theorem without a root region
    /// (real-rootedness on the whole line).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub whole_line_parent_proof: Option<RealRootednessTheoremReport>,
    /// A theorem report proving every parent row real-rooted with every root
    /// in a bounded interval (the Ma--Wang interval check, or the three-term
    /// theorem with a root region), computed when the operator's theorem
    /// gives a region from an interval or the parent statement above is not
    /// proved, and present when proved; replayed by `verify()`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interval_parent_proof: Option<RealRootednessTheoremReport>,
    /// The parent theorem report (replayed by `verify()`) whose conclusion
    /// carries the parent's interlacing: the whole-line three-term report,
    /// or the three-term / Liu--Wang check behind a proved parent statement.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_interlacing_report: Option<RealRootednessTheoremReport>,
    /// Whether consecutive input rows interlace (the `preserver_image`
    /// interlacing rule): claimed only for an operator independent of `n`
    /// from the parent's proved interlacing, never from finitely many rows.
    pub interlacing: PreserverInterlacingClaim,
    /// The row route of the Lean library RealRooted
    /// (<https://github.com/PerAlexandersson/RealRooted>) matching the
    /// operator's theorem, when one exists.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lean_route: Option<String>,
    /// RealRooted (or Mathlib) declarations behind the operator's theorem.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lean_declarations: Vec<String>,
    /// The uniform summary; its verdict is the verdict of the lift.
    pub lift_summary: LiftSummary,
    /// What the verdict means for the input rows.
    pub statement: String,
    /// Further remarks (division remainder, limits).
    pub notes: Vec<String>,
    /// The runnable rule of the image `lambda(n) T_n(P_(n+s))` computed from
    /// the parent rule ([`PreserverImageCandidate::image_recurrence`]), with
    /// its own certificate; present when requested
    /// ([`PreserverImageOptions::image_rule`]), the identity is certified
    /// and a rule was found.  Since the identity is certified, the rule
    /// generates the input rows `Q_n`; `verify()` replays it and compares
    /// its first rows with the input's.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub image_rule: Option<ImageRecurrence>,
}

// ---------------------------------------------------------------------------
// Shared preparation
// ---------------------------------------------------------------------------

/// `n0 + s`, which must be at least `p0`.
fn parent_start(n0: usize, shift: i64, p0: usize) -> Result<usize, RecurrenceError> {
    if shift.unsigned_abs() > MAX_PRESERVER_IMAGE_INDEX_SHIFT as u64 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "|index_shift| = {} exceeds {MAX_PRESERVER_IMAGE_INDEX_SHIFT}",
            shift.unsigned_abs()
        )));
    }
    let start = i64::try_from(n0)
        .ok()
        .and_then(|n0| n0.checked_add(shift))
        .and_then(|start| usize::try_from(start).ok())
        .filter(|start| *start >= p0)
        .ok_or_else(|| {
            RecurrenceError::InvalidInput(format!(
                "the candidate uses parent row n0 + s = {n0} + {shift}, before the parent's \
                 first index {p0}"
            ))
        })?;
    if start - p0 > MAX_PARENT_OFFSET {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the first used parent row {start} is more than {MAX_PARENT_OFFSET} rows after the \
             parent's first index {p0}"
        )));
    }
    Ok(start)
}

/// The shifted parent `P'_n = P_(n+s)`, `n >= n0`, as a prepared rule with
/// first index `n0`: the parent rule with `n -> n + s`, and parent rows as
/// seeds (enough of them that the shifted leader is nonzero at every later
/// target).
pub(super) fn shifted_parent(
    parent: &PreparedLiftRule,
    n0: usize,
    shift: i64,
) -> Result<PreparedLiftRule, RecurrenceError> {
    let start = parent_start(n0, shift, parent.first_index)?;
    let order = parent.leader_shift;
    let parent_end = parent.first_index + parent.seeds.len();
    let needed = parent_end.saturating_sub(start).max(order) + order + 1;
    let offset = start - parent.first_index;
    let rows = parent.rows(offset + needed)?;
    let rule = parent.rule.shift_index(shift);
    prepare_lift_rule(&rule, &rows[offset..], n0, "shifted parent")
}

/// The canonical record of a prepared rule (rule, first index, identifying
/// seeds).
fn record(rule: &PreparedLiftRule) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(&rule.rule, rule.first_index, &rule.seeds, None)
}

fn prepare_record(json: &RecurrenceJson, name: &str) -> Result<PreparedLiftRule, RecurrenceError> {
    let (rule, first_index, rows) = json.to_recurrence_parts()?;
    prepare_lift_rule(&rule, &rows, first_index, name)
}

/// The ratio form `(alpha, beta)` of `mu = lambda / delta`:
/// `mu(n)/mu(n-1) = a(n) delta(n-1) / (b(n) delta(n))`.
fn effective_ratio(gauge: &RatioGauge, delta: &[Q]) -> (Vec<Q>, Vec<Q>) {
    let previous = vec_poly::compose(delta, &[-Q::one(), Q::one()]);
    (
        vec_poly::trim(&vec_poly::mul(&gauge.numerator, &previous)),
        vec_poly::trim(&vec_poly::mul(&gauge.denominator, delta)),
    )
}

/// `delta(n) != 0` for every `n >= n0`, decided exactly.
fn denominator_decision(form: &OperatorForm, n0: usize) -> Option<IndexSignDecision> {
    form.index_denominator().map(|delta| {
        decide_index_rational_function_sign(
            delta,
            &[Q::one()],
            &IntegerRange::starting_at(n0),
            IndexSignCondition::Nonzero,
        )
    })
}

/// The candidate rows `H_n = lambda(n) T_n(P'_n)`, `n = n0, ..., n0 +
/// count - 1`; `None` where the operator is undefined.
fn candidate_rows(
    shifted: &PreparedLiftRule,
    form: &OperatorForm,
    gauge: &RatioGauge,
    n0: usize,
    count: usize,
) -> Result<Vec<Option<Vec<Q>>>, RecurrenceError> {
    let parent_rows = shifted.rows(count)?;
    let lambdas = gauge.values(n0, count);
    Ok(parent_rows
        .iter()
        .zip(&lambdas)
        .enumerate()
        .map(|(offset, (row, lambda))| {
            form.apply_row(row, n0 + offset)
                .map(|image| vec_poly::trim(&vec_poly::scale(&image, lambda)))
        })
        .collect())
}

fn exceptional_bases(multiplier: &Biv, input_from: usize, parent_from: usize) -> Vec<usize> {
    let mut bases: Vec<usize> = (input_from..parent_from).collect();
    bases.extend(roots_in_domain(
        &n_coefficients(multiplier),
        input_from.max(parent_from),
        None,
    ));
    bases.sort_unstable();
    bases.dedup();
    bases
}

/// Rows `n0, ...` needed to evaluate the input relation at every listed base.
fn rows_needed(bases: &[usize], leader_shift: usize, n0: usize, identifying: usize) -> usize {
    bases
        .iter()
        .map(|base| base + leader_shift + 1 - n0)
        .max()
        .unwrap_or(0)
        .max(identifying)
}

/// The identifying rows agree and `M` vanishes on the candidate rows at
/// every exceptional base (exact).
fn check_rows(
    input: &PreparedLiftRule,
    rows: &[Option<Vec<Q>>],
    exceptional: &[usize],
) -> Result<(), RecurrenceError> {
    let failed = |what: String| Err(RecurrenceError::VerificationFailed(what));
    let n0 = input.first_index;
    let mut defined = Vec::with_capacity(rows.len());
    for (offset, row) in rows.iter().enumerate() {
        match row {
            Some(row) => defined.push(row.clone()),
            None => return failed(format!("the candidate row {} is undefined", n0 + offset)),
        }
    }
    for (offset, seed) in input.seeds.iter().enumerate() {
        if defined.get(offset) != Some(&vec_poly::trim(seed)) {
            return failed(format!(
                "the identifying row {} differs from the candidate row",
                n0 + offset
            ));
        }
    }
    for base in exceptional {
        if *base < n0 {
            return failed(format!(
                "the exceptional base {base} precedes the first index"
            ));
        }
        let value = input.relation.apply(&defined, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at exceptional base {base}"))
        })?;
        if !value.iter().all(Zero::is_zero) {
            return failed(format!(
                "the input relation does not vanish on the candidate rows at base {base}"
            ));
        }
    }
    Ok(())
}

fn check_schedule_domain(operator: &PreserverOperator, n0: usize) -> Result<(), RecurrenceError> {
    if let Some(schedule) = operator.schedule() {
        if schedule.first_index > n0 {
            return Err(RecurrenceError::InvalidInput(format!(
                "the operator schedule starts at {}, after the input's first index {n0}",
                schedule.first_index
            )));
        }
    }
    Ok(())
}

/// Certify one candidate preserver-image lift `Q_n = lambda(n) T_n(P_(n+s))`
/// of the input rows (`input_rows[0]` is `Q_(first_index)`; at least the
/// identifying rows) and run the parent provers.
///
/// The identity is certified by the argument of the `preserver_image`
/// module documentation (pseudo-division of the gauged input relation
/// composed with, or conjugated through, `T`, by the shifted parent
/// relation; the Hadamard conjugation for hypergeometric multipliers),
/// the side conditions of the operator's theorem are decided exactly for
/// every `n >= n0`, and the parent statement is settled by the shared parent
/// provers (with the three-term theorem on the whole line as a fallback for
/// operators that need only real-rootedness).  The verdict is that of the
/// [`PreserverImageCertification::lift_summary`]; a verdict other than
/// `proved` says nothing against real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError`] for malformed input (rules that do not reproduce their
/// rows, a leader not proved nonzero, an invalid operator, schedule or gauge,
/// a shift before the parent's first index) and over the limits.
///
/// ```
/// use polytool::recurrence::*;
///
/// let q = |v: i64| BigRational::from_integer(v.into());
/// let bivariate = |rows: Vec<Vec<i64>>| {
///     BivarPoly::from_coefficients_normalized(
///         rows.into_iter().map(|row| row.into_iter().map(q).collect()).collect(),
///     )
/// };
/// // Q_n = A_n(2x + 1) for the Eulerian polynomials A_n, n >= 1:
/// // Q_n = (n - 2x + 2n x) Q_(n-1) - (x + 2x^2) Q'_(n-1), Q_1 = 1.
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, bivariate(vec![vec![0, -2], vec![1, 2]])),
///         RecurrenceTerm::new(1, 1, bivariate(vec![vec![0, -1, -2]])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let rows = rule.generate_rows_rational(&[vec![q(1)]], 1, 12).expect("rows");
/// let candidate = PreserverImageCandidate {
///     parent_name: "eulerian".into(),
///     parent: GammaLiftParent::eulerian().recurrence,
///     index_shift: 0,
///     gauge: VeroneseSectionGauge::None,
///     operator: PreserverOperator::Affine { slope: "2".into(), intercept: "1".into() },
/// };
/// let certification = certify_preserver_image_lift(
///     &rule,
///     &rows,
///     1,
///     &candidate,
///     &PreserverImageOptions::default(),
/// )?;
/// certification.verify()?;
/// assert!(certification.is_proved());
/// assert_eq!(
///     certification.lift_summary.root_transfer.input_root_region.to_string(),
///     "(-inf, -1/2]"
/// );
/// # Ok::<(), RecurrenceError>(())
/// ```
pub fn certify_preserver_image_lift(
    input: &Recurrence,
    input_rows: &[Vec<BigRational>],
    first_index: usize,
    candidate: &PreserverImageCandidate,
    options: &PreserverImageOptions,
) -> Result<PreserverImageCertification, RecurrenceError> {
    if options.checked_rows > MAX_PRESERVER_IMAGE_CHECKED_ROWS
        || options.max_operator_terms > MAX_PRESERVER_IMAGE_OPERATOR_TERMS
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "checked_rows <= {MAX_PRESERVER_IMAGE_CHECKED_ROWS} and max_operator_terms <= \
             {MAX_PRESERVER_IMAGE_OPERATOR_TERMS} are required"
        )));
    }
    let n0 = first_index;
    let input = prepare_lift_rule(input, input_rows, n0, "input")?;
    let parent = prepare_record(&candidate.parent, "parent")?;
    let shifted = shifted_parent(&parent, n0, candidate.index_shift)?;
    check_schedule_domain(&candidate.operator, n0)?;
    let form = candidate.operator.form()?;
    let gauge = candidate.gauge.ratio_form(n0)?;
    let outcome = certify_identity(&input, &parent, &shifted, candidate, &form, &gauge, options)?;
    let transfer = operator_transfer(&candidate.operator, n0);
    let needed_from = parent_start(n0, candidate.index_shift, parent.first_index)?;
    let (parent_proof, whole_line_parent_proof) =
        prove_parents(&parent, needed_from, options, transfer.from_real_line)?;
    let interval_parent_proof = if wants_interval_parent(
        &candidate.operator,
        parent_proof.status == VeroneseParentProofStatus::Proved,
    ) {
        find_bounded_root_region_report(&parent.rule, parent.first_index, &parent.seeds)?
    } else {
        None
    };
    let parent_interlacing_report = parent_interlacing_report(
        &parent,
        &parent_proof,
        whole_line_parent_proof.as_ref(),
        interval_parent_proof.as_ref(),
    )?;
    let interlacing = interlacing_claim(
        candidate,
        n0,
        &outcome.identity,
        parent_interlacing_report.as_ref(),
    );
    let lift_summary = preserver_image_summary(PreserverImageSummaryParts {
        candidate,
        input_first_index: n0,
        identity: &outcome.identity,
        parent_proof: &parent_proof,
        whole_line_parent_proof: whole_line_parent_proof.as_ref(),
        interval_parent_proof: interval_parent_proof.as_ref(),
        interlacing: &interlacing,
    });
    let mut notes = outcome.notes;
    let image_rule =
        certified_image_rule(candidate, &input, &outcome.identity, options, &mut notes)?;
    let s = candidate.index_shift;
    let lift = format!(
        "lambda(n) T_n(P_(n{s:+})) with T = {}",
        candidate.operator.text()
    );
    let statement = if lift_summary.is_proved() {
        format!(
            "proved: Q_n = {lift} for every n >= {n0}, the side conditions hold for every n, \
             and {}; so every Q_n, n >= {n0}, is zero or real-rooted with roots in {}",
            lift_summary.parent_proof.statement(),
            lift_summary.root_transfer.input_root_region
        )
    } else {
        match &outcome.identity {
            PreserverImageIdentityStatus::Certified { .. } => format!(
                "Q_n = {lift} for every n >= {n0} is certified; the lift is {:?} ({})",
                lift_summary.verdict,
                lift_summary.parent_proof.statement()
            ),
            PreserverImageIdentityStatus::SupportedByRows { rows, reason } => format!(
                "Q_n = {lift} holds on the {rows} compared rows only (not proved: {reason})"
            ),
            PreserverImageIdentityStatus::Refuted { index, reason } => {
                format!("refuted at n = {index}: {reason}")
            }
        }
    };
    Ok(PreserverImageCertification {
        schema: PRESERVER_IMAGE_CERTIFICATION_SCHEMA.into(),
        candidate: candidate.clone(),
        input_first_index: n0,
        identifying_rows: input.seeds.len(),
        rows_compared: outcome.rows_compared,
        identity: outcome.identity,
        parent_proof,
        whole_line_parent_proof,
        interval_parent_proof,
        parent_interlacing_report,
        interlacing,
        lean_route: transfer.lean_route,
        lean_declarations: transfer.lean_declarations,
        lift_summary,
        statement,
        notes,
        image_rule,
    })
}

/// Whether the first `count` rows of the image rule are those of `input`.
fn image_rows_match(
    image: &ImageRecurrence,
    input: &PreparedLiftRule,
    count: usize,
) -> Result<bool, RecurrenceError> {
    let generated: Vec<Vec<Q>> = image
        .rows(count)?
        .iter()
        .map(|row| vec_poly::trim(row))
        .collect();
    Ok(generated == input.rows(count)?)
}

/// The computed image rule of a certified identity (when requested), with
/// its rows checked against the input's; a rule that cannot be computed is
/// a note.
fn certified_image_rule(
    candidate: &PreserverImageCandidate,
    input: &PreparedLiftRule,
    identity: &PreserverImageIdentityStatus,
    options: &PreserverImageOptions,
    notes: &mut Vec<String>,
) -> Result<Option<ImageRecurrence>, RecurrenceError> {
    let Some(rule_options) = &options.image_rule else {
        return Ok(None);
    };
    if !matches!(identity, PreserverImageIdentityStatus::Certified { .. }) {
        return Ok(None);
    }
    match candidate.image_recurrence(input.first_index, rule_options) {
        Ok(image) => {
            let count = input.seeds.len() + image.order + 1;
            if !image_rows_match(&image, input, count)? {
                return Err(RecurrenceError::VerificationFailed(
                    "the computed image rule disagrees with the certified input rows".into(),
                ));
            }
            Ok(Some(image))
        }
        Err(error) => {
            notes.push(format!("no image rule was computed: {error}"));
            Ok(None)
        }
    }
}

mod replay;
use replay::*;
mod identity;
use identity::*;
