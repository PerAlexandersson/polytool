//! The candidate, the certificates and the certification of one
//! preserver-image lift `Q_n = lambda(n) T_n(P_(n+s))` (module
//! documentation of `preserver_image`).

use super::interlacing::{interlacing_claim, PreserverInterlacingClaim};
use super::interval::{interval_parent_report, wants_interval_parent};
use super::operator::{dividend, gauged_input, OperatorForm, PreserverOperator};
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
    prepare_rule, shift_free_from_json, shift_free_json, PreparedRule,
};
use crate::recurrence::lift::hadamard_factor::search::{certify_candidate, remove_remainder};
use crate::recurrence::lift::parent_proof::{prove_parent, ParentRecord};
use crate::recurrence::lift::section_certify::{
    n_coefficients, n_polynomial, pseudo_divide_row, x_free_leader, RatioGauge,
};
use crate::recurrence::lift::section_operators::OperatorMatrix;
use crate::recurrence::{
    check_liu_wang_derivative_lag_theorem, check_liu_wang_theorem,
    check_three_term_recurrence_theorem, decode_poly, encode_poly, BivarPolyJson,
    HadamardAnnihilatorSystem, HadamardFactor, HadamardFactorCertificate,
    HadamardFactorIdentityStatus, HadamardFactorLiftOptions, ImplicitRelation,
    ImplicitRelationJson, LiftSummary, OperatorTermJson, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, Recurrence, RecurrenceError, RecurrenceJson,
    VeroneseParentAttestation, VeroneseParentProof, VeroneseParentProofRoute,
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
            c.y_degree().unwrap_or(0) <= max_degree
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
}

impl Default for PreserverImageOptions {
    fn default() -> Self {
        Self {
            checked_rows: 16,
            max_operator_terms: MAX_PRESERVER_IMAGE_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
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
    parent: &PreparedRule,
    n0: usize,
    shift: i64,
) -> Result<PreparedRule, RecurrenceError> {
    let start = parent_start(n0, shift, parent.first_index)?;
    let order = parent.leader_shift;
    let parent_end = parent.first_index + parent.seeds.len();
    let needed = parent_end.saturating_sub(start).max(order) + order + 1;
    let offset = start - parent.first_index;
    let rows = parent.rows(offset + needed)?;
    let rule = parent.rule.shift_index(shift);
    prepare_rule(&rule, &rows[offset..], n0, "shifted parent")
}

/// The canonical record of a prepared rule (rule, first index, identifying
/// seeds).
fn record(rule: &PreparedRule) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(&rule.rule, rule.first_index, &rule.seeds, None)
}

fn prepare_record(json: &RecurrenceJson, name: &str) -> Result<PreparedRule, RecurrenceError> {
    let (rule, first_index, rows) = json.to_recurrence_parts()?;
    prepare_rule(&rule, &rows, first_index, name)
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
    shifted: &PreparedRule,
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
    input: &PreparedRule,
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

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

/// Replay an annihilation against the shifted parent from base `from`:
/// returns the row multiplier `p(x) mu_2(n, x)`, the pairs `(A, Z)` and the
/// bases where `mu_2` vanishes (exceptional).
fn replay_annihilation(
    annihilation: &PreserverRowAnnihilation,
    shifted: &PreparedRule,
    from: usize,
) -> Result<(Biv, AnnihilatorTerms, Vec<usize>), RecurrenceError> {
    let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    match annihilation {
        PreserverRowAnnihilation::Single {
            row_multiplier,
            annihilator,
            quotient,
            step_multiplier,
            step_quotient,
            start_base,
            restart_bases,
        } => {
            let p = Biv::from_y(&decode_poly(row_multiplier)?);
            if p.is_zero() {
                return failed("the row multiplier is zero");
            }
            let a = shift_free_from_json(annihilator)?;
            match (step_multiplier, start_base) {
                (Some(nu), Some(start)) => {
                    let invariance = RowAnnihilatorInvariance {
                        step_multiplier: biv_from_json(nu)?,
                        step_quotient: shift_free_from_json(step_quotient)?,
                        start_base: *start,
                        restart_bases: restart_bases.clone(),
                    };
                    verify_invariance(&a, &invariance, shifted, from)?;
                }
                (None, None) if step_quotient.is_empty() && restart_bases.is_empty() => {
                    let Some(product) = ProductParent::recognise(shifted) else {
                        return failed("the shifted parent is not a product parent");
                    };
                    let expected = product.annihilator();
                    if expected != a || !product.annihilator_is_invariant(&a) {
                        return failed("the row annihilator does not match the parent");
                    }
                }
                _ => return failed("incomplete row-annihilator invariance fields"),
            }
            Ok((p, vec![(a, shift_free_from_json(quotient)?)], Vec::new()))
        }
        PreserverRowAnnihilation::System { system } => {
            let (decoded, terms) = system_from_json(system)?;
            verify_system(
                &decoded,
                system.start_base,
                &system.restart_bases,
                shifted,
                from,
            )?;
            let second = biv_from_json(&system.row_multiplier)?;
            if second.is_zero() {
                return failed("the row-annihilator system multiplier is zero");
            }
            let bases = vanishing_bases(&second, from);
            Ok((second, terms, bases))
        }
    }
}

/// The JSON form of an annihilation found by the Hadamard machinery, with
/// the full row multiplier and the pairs `(A, Z)`.
fn annihilation_record(
    annihilation: RowAnnihilation,
) -> Result<(PreserverRowAnnihilation, Biv, AnnihilatorTerms), RecurrenceError> {
    let p = Biv::from_y(&annihilation.row_multiplier);
    match annihilation.kind {
        RowAnnihilationKind::Single {
            annihilator,
            quotient,
            invariance,
        } => {
            let record = PreserverRowAnnihilation::Single {
                row_multiplier: encode_poly(&annihilation.row_multiplier),
                annihilator: shift_free_json(&annihilator),
                quotient: shift_free_json(&quotient),
                step_multiplier: invariance.as_ref().map(|i| biv_json(&i.step_multiplier)),
                step_quotient: invariance
                    .as_ref()
                    .map_or_else(Vec::new, |i| shift_free_json(&i.step_quotient)),
                start_base: invariance.as_ref().map(|i| i.start_base),
                restart_bases: invariance.map_or_else(Vec::new, |i| i.restart_bases),
            };
            Ok((record, p, vec![(annihilator, quotient)]))
        }
        RowAnnihilationKind::System {
            system,
            row_multiplier,
            quotients,
            start_base,
            restart_bases,
        } => {
            let mut json = system_json(&system, &quotients, start_base, restart_bases);
            json.row_multiplier = biv_json(&row_multiplier);
            let (_, terms) = system_from_json(&json)?;
            Ok((
                PreserverRowAnnihilation::System { system: json },
                p.mul(&row_multiplier),
                terms,
            ))
        }
    }
}

impl PreserverDivisionCertificate {
    /// Replay the certificate from its recorded input, parent, shift, gauge
    /// and operator: regenerate the input and parent rows from the recorded
    /// rules and seeds (and re-validate both rules), recompute `L'` and `E`,
    /// check `w E = sum_nu H_nu S^nu L'` exactly, check that every integer
    /// root of `w` from the valid bases on is listed as exceptional, and
    /// re-check the identifying rows and the input relation at every
    /// exceptional base on freshly computed candidate rows.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let input = prepare_record(&self.input, "input")?;
        let n0 = input.first_index;
        if self.identifying_indices != (n0..n0 + input.seeds.len()).collect::<Vec<_>>() {
            return failed("the identifying indices do not match the input");
        }
        let parent = prepare_record(&self.parent, "parent")?;
        let shifted = shifted_parent(&parent, n0, self.index_shift)?;
        if ImplicitRelation::from_op(&shifted.relation).to_json() != self.shifted_parent_relation {
            return failed("the shifted parent relation does not match its recomputation");
        }
        let (k, leader) = x_free_leader(&shifted.relation)?;
        if k != self.parent_leader_shift {
            return failed("the parent leader shift does not match");
        }
        check_schedule_domain(&self.operator, n0)?;
        let form = self.operator.form()?;
        if matches!(form, OperatorForm::Hadamard { .. }) {
            return failed("a Hadamard multiplier needs the Hadamard conjugation certificate");
        }
        if let Some(decision) = denominator_decision(&form, n0) {
            if !decision.holds() {
                return failed("the operator denominator is not proved nonzero on the domain");
            }
        }
        let gauge = self.gauge.ratio_form(n0)?;
        let delta = form.index_denominator().unwrap_or(&[]).to_vec();
        let delta = if delta.is_empty() {
            vec![Q::one()]
        } else {
            delta
        };
        let (alpha, beta) = effective_ratio(&gauge, &delta);
        let gauged = gauged_input(&input.relation, input.leader_shift, &alpha, &beta);
        let Some(e) = dividend(&gauged, &form) else {
            return failed("the operator has no dividend");
        };
        if ImplicitRelation::from_op(&e).to_json() != self.dividend_relation {
            return failed("the dividend does not match its recomputation");
        }
        let multiplier = n_polynomial(&decode_poly(&self.multiplier)?);
        if multiplier.is_zero() {
            return failed("the multiplier is zero");
        }
        let mut quotients: BTreeMap<usize, DOp> = BTreeMap::new();
        for entry in &self.quotients {
            let h = shift_free_from_json(&entry.differentiations)?;
            let slot = quotients.entry(entry.shift).or_default();
            *slot = slot.add(&h);
        }
        let _ = leader;
        let input_from = input.valid_from_base();
        let parent_from = shifted.valid_from_base().max(n0);
        if input_from != self.input_valid_from_base || parent_from != self.parent_valid_from_base {
            return failed("the recorded valid bases do not match");
        }
        let from = input_from.max(parent_from);
        let (row_multiplier, annihilators, mut expected) = match &self.row_annihilation {
            None => (Biv::one(), Vec::new(), Vec::new()),
            Some(annihilation) => replay_annihilation(annihilation, &shifted, from)?,
        };
        if !identity_holds(
            &e,
            &shifted.relation,
            &row_multiplier.mul(&multiplier),
            &quotients,
            &annihilators,
        ) {
            return failed("the identity w E = sum H_nu S^nu L' (+ Z A) fails");
        }
        expected.extend(exceptional_bases(&multiplier, input_from, parent_from));
        if expected
            .iter()
            .any(|base| !self.exceptional_bases.contains(base))
        {
            return failed("an exceptional base is missing");
        }
        let needed = rows_needed(
            &self.exceptional_bases,
            input.leader_shift,
            n0,
            input.seeds.len(),
        );
        if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the certificate needs {needed} candidate rows, more than \
                 {MAX_PRESERVER_IMAGE_CHECKED_ROWS}"
            )));
        }
        let rows = candidate_rows(&shifted, &form, &gauge, n0, needed)?;
        check_rows(&input, &rows, &self.exceptional_bases)
    }
}

impl PreserverHadamardCertificate {
    /// Replay: [`HadamardFactorCertificate::verify`], the agreement of its
    /// factor, gauge, input relation and parent with the records here (the
    /// parent being the shifted parent `P'` recomputed from the parent
    /// record), and the identifying rows and exceptional bases re-checked
    /// on candidate rows computed from the regenerated parent rows.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        self.hadamard.verify()?;
        let input = prepare_record(&self.input, "input")?;
        let n0 = input.first_index;
        let parent = prepare_record(&self.parent, "parent")?;
        let shifted = shifted_parent(&parent, n0, self.index_shift)?;
        if self.hadamard.factor != self.factor {
            return failed("the Hadamard certificate has another factor");
        }
        if self.hadamard.parent != record(&shifted) {
            return failed("the Hadamard certificate's parent is not the shifted parent");
        }
        if self.hadamard.input_relation != ImplicitRelation::from_op(&input.relation).to_json()
            || self.hadamard.input_leader_shift != input.leader_shift
        {
            return failed("the Hadamard certificate's input relation differs from the input");
        }
        if self.hadamard.identifying_indices != (n0..n0 + input.seeds.len()).collect::<Vec<_>>() {
            return failed("the identifying indices do not match the input");
        }
        if self.hadamard.input_valid_from_base != input.valid_from_base()
            || self.hadamard.parent_valid_from_base != shifted.valid_from_base().max(n0)
        {
            return failed("the recorded valid bases do not match");
        }
        let gauge = self.hadamard.gauge.ratio_form(n0)?;
        let form = OperatorForm::Hadamard {
            factor: self.factor.clone(),
        };
        let needed = rows_needed(
            &self.hadamard.exceptional_bases,
            input.leader_shift,
            n0,
            input.seeds.len(),
        );
        if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the certificate needs {needed} candidate rows"
            )));
        }
        let rows = candidate_rows(&shifted, &form, &gauge, n0, needed)?;
        check_rows(&input, &rows, &self.hadamard.exceptional_bases)
    }
}

impl PreserverImageCertificate {
    /// Replay the certificate ([`PreserverDivisionCertificate::verify`] or
    /// [`PreserverHadamardCertificate::verify`]).
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        match self {
            Self::OperatorDivision(certificate) => certificate.verify(),
            Self::HadamardConjugation(certificate) => certificate.verify(),
        }
    }

    /// The recorded input, parent, shift, gauge and operator.
    fn parts(
        &self,
    ) -> (
        &RecurrenceJson,
        &RecurrenceJson,
        i64,
        &VeroneseSectionGauge,
        PreserverOperator,
    ) {
        match self {
            Self::OperatorDivision(c) => (
                &c.input,
                &c.parent,
                c.index_shift,
                &c.gauge,
                c.operator.clone(),
            ),
            Self::HadamardConjugation(c) => (
                &c.input,
                &c.parent,
                c.index_shift,
                &c.hadamard.gauge,
                PreserverOperator::HadamardMultiplier {
                    factor: c.factor.clone(),
                },
            ),
        }
    }
}

impl PreserverImageCertification {
    /// Whether every `Q_n`, `n >= n0`, is proved zero or real-rooted in the
    /// summary's input root region.
    pub fn is_proved(&self) -> bool {
        self.lift_summary.is_proved()
    }

    /// The [`LiftSummary`] computed from the record.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        preserver_image_summary(PreserverImageSummaryParts {
            candidate: &self.candidate,
            input_first_index: self.input_first_index,
            identity: &self.identity,
            parent_proof: &self.parent_proof,
            whole_line_parent_proof: self.whole_line_parent_proof.as_ref(),
            interval_parent_proof: self.interval_parent_proof.as_ref(),
        })
    }

    /// Replay the certificate (when certified), check that it belongs to the
    /// candidate, replay the whole-line parent theorem report (when present),
    /// and check that the summary equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if self.schema != PRESERVER_IMAGE_CERTIFICATION_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        if let PreserverImageIdentityStatus::Certified { certificate } = &self.identity {
            certificate.verify()?;
            let (input, parent, shift, gauge, operator) = certificate.parts();
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if *parent != candidate_parent
                || shift != self.candidate.index_shift
                || *gauge != self.candidate.gauge
                || operator != self.candidate.operator
                || input.first_index != self.input_first_index
            {
                return failed("the certificate does not belong to the candidate");
            }
        }
        if let Some(report) = &self.whole_line_parent_proof {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if report.recurrence != candidate_parent || report.requested_root_region.is_some() {
                return failed("the whole-line parent report does not belong to the parent");
            }
        }
        if let Some(report) = &self.interval_parent_proof {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            let bounded = report
                .conclusion
                .as_ref()
                .and_then(|c| c.root_region.as_ref())
                .is_some_and(|r| r.lower.is_some() && r.upper.is_some());
            if report.recurrence != candidate_parent || !report.is_proved() || !bounded {
                return failed("the interval parent report does not belong to the parent");
            }
        }
        if let Some(report) = &self.parent_interlacing_report {
            report.verify()?;
            let candidate_parent = record(&prepare_record(&self.candidate.parent, "parent")?);
            if report.recurrence != candidate_parent || !report.is_proved() {
                return failed("the parent interlacing report does not belong to the parent");
            }
        }
        let claim = interlacing_claim(
            &self.candidate,
            self.input_first_index,
            &self.identity,
            self.parent_interlacing_report.as_ref(),
        );
        if claim != self.interlacing {
            return failed("the interlacing claim does not match the one recomputed");
        }
        let summary = self.recompute_lift_summary();
        if summary != self.lift_summary || !summary.is_consistent() {
            return failed("the lift summary does not match the one recomputed from the record");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// The outcome of the identity part of a certification.
struct IdentityOutcome {
    identity: PreserverImageIdentityStatus,
    rows_compared: usize,
    notes: Vec<String>,
}

fn supported(rows: usize, reason: String) -> PreserverImageIdentityStatus {
    PreserverImageIdentityStatus::SupportedByRows { rows, reason }
}

/// The first input row that differs from the candidate (an undefined
/// candidate row counts as a difference).
fn first_mismatch(input: &[Vec<Q>], candidate: &[Option<Vec<Q>>]) -> Option<(usize, bool)> {
    input
        .iter()
        .zip(candidate)
        .enumerate()
        .find_map(|(offset, (a, b))| match b {
            None => Some((offset, true)),
            Some(b) if vec_poly::trim(a) != *b => Some((offset, false)),
            Some(_) => None,
        })
}

fn refuted(n0: usize, offset: usize, undefined: bool) -> PreserverImageIdentityStatus {
    let index = n0 + offset;
    PreserverImageIdentityStatus::Refuted {
        index,
        reason: if undefined {
            format!("the operator is undefined on the parent row used for n = {index}")
        } else {
            format!("Q_{index} differs from the candidate lambda(n) T_n(P_(n+s))")
        },
    }
}

#[allow(clippy::too_many_arguments)]
fn certify_identity(
    input: &PreparedRule,
    parent: &PreparedRule,
    shifted: &PreparedRule,
    candidate: &PreserverImageCandidate,
    form: &OperatorForm,
    gauge: &RatioGauge,
    options: &PreserverImageOptions,
) -> Result<IdentityOutcome, RecurrenceError> {
    let n0 = input.first_index;
    let mut notes = Vec::new();
    let mut count = options.checked_rows.max(input.seeds.len());
    let compare =
        |count: usize| -> Result<Option<PreserverImageIdentityStatus>, RecurrenceError> {
            let inputs = input.rows(count)?;
            let rows = candidate_rows(shifted, form, gauge, n0, count)?;
            Ok(first_mismatch(&inputs, &rows)
                .map(|(offset, undefined)| refuted(n0, offset, undefined)))
        };
    if let Some(status) = compare(count)? {
        return Ok(IdentityOutcome {
            identity: status,
            rows_compared: count,
            notes,
        });
    }
    if let Some(decision) = denominator_decision(form, n0) {
        if !decision.holds() {
            return Ok(IdentityOutcome {
                identity: supported(
                    count,
                    format!(
                        "the operator coefficients are not proved finite at every n >= {n0}: {}",
                        decision.describe("delta(n)")
                    ),
                ),
                rows_compared: count,
                notes,
            });
        }
    }

    // Hypergeometric multipliers: the Hadamard conjugation.
    if let OperatorForm::Hadamard { factor } = form {
        let hadamard_options = HadamardFactorLiftOptions {
            checked_rows: options
                .checked_rows
                .clamp(1, MAX_PRESERVER_IMAGE_CHECKED_ROWS),
            max_operator_terms: options
                .max_operator_terms
                .min(crate::recurrence::MAX_HADAMARD_OPERATOR_TERMS),
            ..HadamardFactorLiftOptions::default()
        };
        let (status, rows) =
            certify_candidate(input, shifted, factor, &candidate.gauge, &hadamard_options)?;
        count = count.max(rows);
        let identity = match status {
            HadamardFactorIdentityStatus::Certified { certificate } => {
                let wrapped = PreserverHadamardCertificate {
                    input: record(input),
                    parent: record(parent),
                    index_shift: candidate.index_shift,
                    factor: factor.clone(),
                    hadamard: *certificate,
                };
                wrapped.verify()?;
                PreserverImageIdentityStatus::Certified {
                    certificate: Box::new(PreserverImageCertificate::HadamardConjugation(
                        Box::new(wrapped),
                    )),
                }
            }
            HadamardFactorIdentityStatus::SupportedByRows { rows, reason } => supported(
                rows,
                format!("the Hadamard conjugation does not apply: {reason}"),
            ),
            HadamardFactorIdentityStatus::Refuted { index } => {
                PreserverImageIdentityStatus::Refuted {
                    index,
                    reason: format!("Q_{index} differs from the candidate Hadamard image"),
                }
            }
        };
        return Ok(IdentityOutcome {
            identity,
            rows_compared: count,
            notes,
        });
    }

    // Weyl operators and coordinate maps: pseudo-division of the dividend.
    let (k, leader) = match x_free_leader(&shifted.relation) {
        Ok(found) => found,
        Err(error) => {
            return Ok(IdentityOutcome {
                identity: supported(count, format!("the parent relation is not usable: {error}")),
                rows_compared: count,
                notes,
            })
        }
    };
    let delta = form
        .index_denominator()
        .map(<[Q]>::to_vec)
        .unwrap_or_else(|| vec![Q::one()]);
    let (alpha, beta) = effective_ratio(gauge, &delta);
    let gauged = gauged_input(&input.relation, input.leader_shift, &alpha, &beta);
    let e = dividend(&gauged, form).ok_or_else(|| {
        RecurrenceError::VerificationFailed("the operator form has no dividend".into())
    })?;
    let system = OperatorMatrix {
        entries: vec![vec![shifted.relation.clone()]],
    };
    let division = match pseudo_divide_row(
        std::slice::from_ref(&e),
        &system,
        k,
        &leader,
        options.max_operator_terms,
    ) {
        Ok(division) => division,
        Err(RecurrenceError::BudgetExhausted(reason) | RecurrenceError::LimitExceeded(reason)) => {
            return Ok(IdentityOutcome {
                identity: supported(count, format!("the division ran out of budget: {reason}")),
                rows_compared: count,
                notes,
            })
        }
        Err(error) => return Err(error),
    };
    let input_from = input.valid_from_base();
    let parent_from = shifted.valid_from_base().max(n0);
    let from = input_from.max(parent_from);
    let mut multiplier = division.multiplier.clone();
    let mut row_multiplier = Biv::one();
    let mut annihilators: Vec<(DOp, DOp)> = Vec::new();
    let mut annihilation_json = None;
    let mut extra_exceptional = Vec::new();
    if !division.remainder_is_zero() {
        let hadamard_options = HadamardFactorLiftOptions {
            max_operator_terms: options
                .max_operator_terms
                .min(crate::recurrence::MAX_HADAMARD_OPERATOR_TERMS),
            ..HadamardFactorLiftOptions::default()
        };
        if !remainder_is_small(&division.remainder[0]) {
            return Ok(IdentityOutcome {
                identity: supported(
                    count,
                    "the division leaves a nonzero remainder too large for the row-annihilator \
                     search"
                        .into(),
                ),
                rows_compared: count,
                notes,
            });
        }
        let annihilation =
            match remove_remainder(shifted, &division.remainder[0], from, &hadamard_options)? {
                Ok(annihilation) => annihilation,
                Err(reason) => {
                    notes.push(format!(
                        "the pseudo-division of the dividend by the shifted parent relation \
                         leaves a nonzero remainder of shift < {k}"
                    ));
                    return Ok(IdentityOutcome {
                        identity: supported(count, reason),
                        rows_compared: count,
                        notes,
                    });
                }
            };
        multiplier = multiplier.mul(&n_polynomial(&annihilation.index_multiplier));
        let (record, full, pairs) = annihilation_record(annihilation)?;
        if let PreserverRowAnnihilation::System { system } = &record {
            extra_exceptional = vanishing_bases(&biv_from_json(&system.row_multiplier)?, from);
        }
        row_multiplier = full;
        annihilators = pairs;
        annihilation_json = Some(record);
        notes.push(
            "the division remainder is removed by row annihilators of the shifted parent".into(),
        );
    }
    let quotients: BTreeMap<usize, DOp> = division
        .quotients
        .iter()
        .map(|(nu, row)| {
            let mut h = DOp::default();
            for (key, c) in &row[0].0 {
                h.add_term(*key, row_multiplier.mul(c));
            }
            (*nu, h)
        })
        .filter(|(_, h)| !h.is_zero())
        .collect();
    if !identity_holds(
        &e,
        &shifted.relation,
        &row_multiplier.mul(&multiplier),
        &quotients,
        &annihilators,
    ) {
        return Err(RecurrenceError::VerificationFailed(
            "the pseudo-division identity does not hold".into(),
        ));
    }
    let mut exceptional = exceptional_bases(&multiplier, input_from, parent_from);
    exceptional.extend(extra_exceptional);
    exceptional.sort_unstable();
    exceptional.dedup();
    let needed = rows_needed(&exceptional, input.leader_shift, n0, input.seeds.len());
    if needed > MAX_PRESERVER_IMAGE_CHECKED_ROWS {
        return Ok(IdentityOutcome {
            identity: supported(
                count,
                format!("the exceptional bases need {needed} compared rows"),
            ),
            rows_compared: count,
            notes,
        });
    }
    if needed > count {
        count = needed;
        if let Some(status) = compare(count)? {
            return Ok(IdentityOutcome {
                identity: status,
                rows_compared: count,
                notes,
            });
        }
    }
    let rows = candidate_rows(shifted, form, gauge, n0, needed)?;
    if let Err(RecurrenceError::VerificationFailed(reason)) = check_rows(input, &rows, &exceptional)
    {
        return Ok(IdentityOutcome {
            identity: supported(count, reason),
            rows_compared: count,
            notes,
        });
    }
    let identifying = input.seeds.len();
    let s = candidate.index_shift;
    let certificate = PreserverDivisionCertificate {
        input: record(input),
        parent: record(parent),
        index_shift: s,
        gauge: candidate.gauge.clone(),
        operator: candidate.operator.clone(),
        shifted_parent_relation: ImplicitRelation::from_op(&shifted.relation).to_json(),
        parent_leader_shift: k,
        dividend_relation: ImplicitRelation::from_op(&e).to_json(),
        multiplier: encode_poly(&n_coefficients(&multiplier)),
        quotients: quotients
            .iter()
            .map(|(nu, h)| PreserverDivisionQuotient {
                shift: *nu,
                differentiations: shift_free_json(h),
            })
            .collect(),
        row_annihilation: annihilation_json,
        input_valid_from_base: input_from,
        parent_valid_from_base: parent_from,
        exceptional_bases: exceptional,
        identifying_indices: (n0..n0 + identifying).collect(),
        statement: format!(
            "w(n) E = sum_nu H_nu S^nu L' (with any remainder removed by row annihilators of \
             the parent), where L' is the parent rule at \
             P'_n = P_(n{s:+}) (valid from base {parent_from}) and E the input rule gauged by \
             lambda(n)/delta(n) and composed with (or conjugated through) T = {}; so the \
             candidate rows satisfy the input rule at every base n >= {input_from} (the \
             exceptional bases checked exactly), and agreement on the {identifying} \
             identifying rows gives Q_n = lambda(n) T_n(P_(n{s:+})) for every n >= {n0}",
            candidate.operator.text()
        ),
    };
    certificate.verify()?;
    Ok(IdentityOutcome {
        identity: PreserverImageIdentityStatus::Certified {
            certificate: Box::new(PreserverImageCertificate::OperatorDivision(Box::new(
                certificate,
            ))),
        },
        rows_compared: count,
        notes,
    })
}

/// The parent provers: `prove_parent` (rows zero or real-rooted with
/// coefficients of one sign) and, when that is not proved and the
/// operator's theorem needs only real-rootedness, the three-term theorem
/// without a root region.
fn prove_parents(
    parent: &PreparedRule,
    needed_from: usize,
    options: &PreserverImageOptions,
    whole_line_suffices: bool,
) -> Result<(VeroneseParentProof, Option<RealRootednessTheoremReport>), RecurrenceError> {
    let proof = prove_parent(
        &ParentRecord {
            rule: &parent.rule,
            seeds: &parent.seeds,
            first_index: parent.first_index,
        },
        needed_from,
        &options.attestations,
        options.try_root_domain,
    )?;
    if proof.status == VeroneseParentProofStatus::Proved || !whole_line_suffices {
        return Ok((proof, None));
    }
    let report = check_three_term_recurrence_theorem(
        &parent.rule,
        parent.first_index,
        &parent.seeds,
        &RealRootednessTheoremOptions {
            root_region: None,
            ..RealRootednessTheoremOptions::default()
        },
    )?;
    Ok((proof, report.is_proved().then_some(report)))
}

/// The parent theorem report carrying the parent's interlacing: the
/// whole-line report, else (for a parent proved by the parent provers) the
/// first of the three-term and Liu--Wang checks (roots in `(-inf, 0)`) and
/// the derivative-lag check that proves the parent with an interlacing
/// conclusion, else the interval report.
fn parent_interlacing_report(
    parent: &PreparedRule,
    proof: &VeroneseParentProof,
    whole_line: Option<&RealRootednessTheoremReport>,
    interval: Option<&RealRootednessTheoremReport>,
) -> Result<Option<RealRootednessTheoremReport>, RecurrenceError> {
    if let Some(report) = whole_line {
        return Ok(Some(report.clone()));
    }
    let interval = interval
        .filter(|report| {
            report.is_proved()
                && report
                    .conclusion
                    .as_ref()
                    .is_some_and(|c| c.interlacing.is_some())
        })
        .cloned();
    if proof.status != VeroneseParentProofStatus::Proved {
        return Ok(interval);
    }
    let negative = RealRootednessTheoremOptions {
        root_region: Some(crate::real_rootedness::RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let (rule, p0, seeds) = (&parent.rule, parent.first_index, &parent.seeds);
    if proof.route == Some(VeroneseParentProofRoute::Attestation) {
        return Ok(interval);
    }
    let reports = [
        check_three_term_recurrence_theorem(rule, p0, seeds, &negative)?,
        check_liu_wang_theorem(rule, p0, seeds, &negative)?,
        check_liu_wang_derivative_lag_theorem(
            rule,
            p0,
            seeds,
            &RealRootednessTheoremOptions::default(),
        )?,
    ];
    Ok(reports
        .into_iter()
        .find(|report| {
            report.is_proved()
                && report
                    .conclusion
                    .as_ref()
                    .is_some_and(|c| c.interlacing.is_some())
        })
        .or(interval))
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
    let input = prepare_rule(input, input_rows, n0, "input")?;
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
        interval_parent_report(&parent.rule, parent.first_index, &parent.seeds)?
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
    });
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
        notes: outcome.notes,
    })
}
