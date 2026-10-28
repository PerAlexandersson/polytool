//! [`certify_veronese_section_lift`]: certify one candidate Veronese-section
//! lift `Q_n = lambda(n) V_(r, j)(P_(n+s))` of a recurrence-defined sequence
//! with index-dependent and differential terms, prove its parent, and report
//! the result as `polytool.recurrence-veronese-certificate.v1`.
//!
//! The kernel is in `section_operators` (the transport `phi`),
//! `section_certify` (gauge, pseudo-division, the forward certificate) and
//! `parent_proof` (the parent provers).

use super::super::implicit::algebra::{roots_in_domain, Biv, DOp};
use super::parent_proof::{prove_parent, ParentRecord};
use super::section_certify::{
    gauge_parent_relation, n_coefficients, pseudo_divide_row, residue_row, x_free_leader,
};
use super::section_operators::{section_row, transport_to_sections};
use super::*;
use crate::vec_poly;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// JSON schema tag of a serialized [`VeroneseSectionCertification`].
pub const VERONESE_SECTION_CERTIFICATE_SCHEMA: &str = "polytool.recurrence-veronese-certificate.v1";

/// Largest modulus `r` accepted by [`certify_veronese_section_lift`].
pub const MAX_VERONESE_CERTIFY_MODULUS: usize = 4;

/// Largest `|s|` of the index shift of a candidate.
pub const MAX_VERONESE_CERTIFY_INDEX_SHIFT: usize = 4;

/// Largest derivative order of the input and the parent rules.
pub const MAX_VERONESE_CERTIFY_DERIVATIVE_ORDER: usize = 6;

/// Largest number of operator terms in a transported system, a remainder or
/// the quotients of a pseudo-division.
pub const MAX_VERONESE_OPERATOR_TERMS: usize = 10_000;

/// Largest number of rows compared exactly (also the bound on the rows that
/// the exceptional bases of a certificate may require).
pub const MAX_VERONESE_CHECKED_ROWS: usize = 64;

/// Largest distance `n0 + s - p0` between the first used parent row and the
/// parent's first index (the parent rows in between are generated).
pub const MAX_VERONESE_PARENT_OFFSET: usize = 1024;

/// Value of [`VeroneseSectionCertification::lean_replay`] in this version.
pub const LEAN_REPLAY_NOT_ATTEMPTED: &str = "not_attempted";

fn default_gauge() -> VeroneseSectionGauge {
    VeroneseSectionGauge::None
}

/// A candidate lift `Q_n = lambda(n) V_(r, j)(P_(n+s))` for every `n >= n0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseSectionCandidate {
    /// The modulus `r`, `2 <= r <=` [`MAX_VERONESE_CERTIFY_MODULUS`].
    pub modulus: usize,
    /// The residue `j`, `0 <= j < r`.
    pub residue: usize,
    /// The parent `P_n(x)`: a `polytool.recurrence.v1` record (rule, first
    /// index `p0` and seeds).  Its rule must be homogeneous, without `(-1)^n`
    /// factors, with a leader (denominator) free of `x`.
    pub parent: RecurrenceJson,
    /// The index shift `s`, `|s| <=` [`MAX_VERONESE_CERTIFY_INDEX_SHIFT`] and
    /// `n0 + s >= p0`.  Default 0.
    #[serde(default)]
    pub index_shift: i64,
    /// The normalisation `lambda`.  Default none (`lambda = 1`).
    #[serde(default = "default_gauge")]
    pub gauge: VeroneseSectionGauge,
}

/// Options of [`certify_veronese_section_lift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VeroneseSectionCertifyOptions {
    /// Rows `Q_(n0), ...` compared exactly with the candidate (at least the
    /// identifying rows; more when the exceptional bases need them).  At most
    /// [`MAX_VERONESE_CHECKED_ROWS`].  Default 16.
    pub checked_rows: usize,
    /// Operator-term budget, at most [`MAX_VERONESE_OPERATOR_TERMS`].
    /// Default [`MAX_VERONESE_OPERATOR_TERMS`].
    pub max_operator_terms: usize,
    /// Caller-supplied parent attestations (route 4 of the parent provers).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain parent route (`transform_recurrence`).  Default true.
    pub try_root_domain: bool,
}

impl Default for VeroneseSectionCertifyOptions {
    fn default() -> Self {
        Self {
            checked_rows: 16,
            max_operator_terms: MAX_VERONESE_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
        }
    }
}

/// How an identity was certified.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseIdentityCertificateKind {
    /// Pseudo-division of `e_j^T M` by the transported parent system with
    /// remainder zero ([`VeroneseForwardDivisionCertificate`]).
    ForwardDivision,
}

/// Status of the identity `Q_n = lambda(n) V_(r, j)(P_(n+s))`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum VeroneseSectionIdentityStatus {
    /// Proved for every `n >= n0`.
    Certified {
        /// The certificate kind.
        kind: VeroneseIdentityCertificateKind,
        /// The replayable certificate.
        certificate: Box<VeroneseForwardDivisionCertificate>,
    },
    /// The first `rows` rows agree exactly, but the division left a nonzero
    /// remainder: finite evidence only, never a proof.
    SupportedByRows {
        /// Number of rows compared, from `n0`.
        rows: usize,
    },
    /// `Q_index != lambda(index) V_(r, j)(P_(index+s))` exactly.
    Refuted {
        /// The absolute index of the first differing row.
        index: usize,
    },
}

/// The overall verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseSectionVerdict {
    /// The identity is certified and the parent is proved (periodic product,
    /// paper theorem or root domain): every `Q_n`, `n >= n0`, is zero or
    /// real-rooted with nonpositive roots.
    Proved,
    /// The identity is certified and the parent statement is attested by
    /// the caller: the conclusion holds if the attestation does.
    ConditionalOnAttestation,
    /// The identity is certified; the parent statement is open.
    IdentityCertifiedParentOpen,
    /// The identity is certified, but a parent row is nonzero and not
    /// real-rooted with one-signed coefficients: this lift gives no proof.
    IdentityCertifiedParentRefuted,
    /// The identity is supported by the compared rows only.
    IdentitySupportedByRows,
    /// The identity is false.
    IdentityRefuted,
}

impl VeroneseSectionVerdict {
    /// Whether the verdict proves that every `Q_n` is zero or real-rooted.
    pub fn is_proved(self) -> bool {
        self == Self::Proved
    }
}

/// `polytool.recurrence-veronese-certificate.v1`: the certification of one
/// candidate lift.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseSectionCertification {
    /// Always [`VERONESE_SECTION_CERTIFICATE_SCHEMA`].
    pub schema: String,
    /// The candidate as supplied.
    pub candidate: VeroneseSectionCandidate,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number `K` of identifying input rows (`Recurrence::generation_initial_count`).
    pub identifying_rows: usize,
    /// Number of rows `Q_(n0), ...` compared exactly with the candidate.
    pub rows_compared: usize,
    /// The identity status, with the certificate when certified.
    pub identity: VeroneseSectionIdentityStatus,
    /// The parent statement and its proof attempts.
    pub parent_proof: VeroneseParentProof,
    /// The overall verdict.
    pub verdict: VeroneseSectionVerdict,
    /// What the verdict means for the input rows.
    pub statement: String,
    /// Always [`LEAN_REPLAY_NOT_ATTEMPTED`] in this version.
    pub lean_replay: String,
    /// Further remarks (division remainder, limits).
    pub notes: Vec<String>,
    /// The uniform [`LiftSummary`] of this certification
    /// ([`VeroneseSectionCertification::recompute_lift_summary`]); always
    /// present in reports of this version, optional in JSON for records
    /// written before the field existed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lift_summary: Option<LiftSummary>,
}

impl VeroneseSectionCertification {
    /// Replay the certificate (when certified), check that the verdict is
    /// consistent with the identity and parent statuses, and that a recorded
    /// `lift_summary` equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != VERONESE_SECTION_CERTIFICATE_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        let certified = match &self.identity {
            VeroneseSectionIdentityStatus::Certified { certificate, .. } => {
                certificate.verify()?;
                if certificate.modulus != self.candidate.modulus
                    || certificate.residue != self.candidate.residue
                    || certificate.index_shift != self.candidate.index_shift
                    || certificate.gauge != self.candidate.gauge
                {
                    return Err(RecurrenceError::VerificationFailed(
                        "the certificate does not belong to the candidate".into(),
                    ));
                }
                true
            }
            _ => false,
        };
        let expected = verdict(&self.identity, &self.parent_proof);
        if expected != self.verdict || (self.verdict.is_proved() && !certified) {
            return Err(RecurrenceError::VerificationFailed(
                "the verdict is inconsistent with the identity and parent statuses".into(),
            ));
        }
        if let Some(summary) = &self.lift_summary {
            if summary != &self.recompute_lift_summary() || !summary.is_consistent() {
                return Err(RecurrenceError::VerificationFailed(
                    "the lift summary does not match the one recomputed from the record".into(),
                ));
            }
        }
        Ok(())
    }
}

fn verdict(
    identity: &VeroneseSectionIdentityStatus,
    parent: &VeroneseParentProof,
) -> VeroneseSectionVerdict {
    match identity {
        VeroneseSectionIdentityStatus::Refuted { .. } => VeroneseSectionVerdict::IdentityRefuted,
        VeroneseSectionIdentityStatus::SupportedByRows { .. } => {
            VeroneseSectionVerdict::IdentitySupportedByRows
        }
        VeroneseSectionIdentityStatus::Certified { .. } => match parent.status {
            VeroneseParentProofStatus::Proved
                if parent
                    .route
                    .is_some_and(|r| r != VeroneseParentProofRoute::Attestation) =>
            {
                VeroneseSectionVerdict::Proved
            }
            VeroneseParentProofStatus::ConditionalOnAttestation => {
                VeroneseSectionVerdict::ConditionalOnAttestation
            }
            VeroneseParentProofStatus::Refuted => {
                VeroneseSectionVerdict::IdentityCertifiedParentRefuted
            }
            _ => VeroneseSectionVerdict::IdentityCertifiedParentOpen,
        },
    }
}

/// A validated rule with its identifying seeds and relation.
struct PreparedRule {
    rule: Recurrence,
    first_index: usize,
    /// The `K` identifying seeds.
    seeds: Vec<Vec<Q>>,
    /// The rule at base index `n`, leader shift `leader_shift`.
    relation: DOp,
    leader_shift: usize,
}

impl PreparedRule {
    fn rows(&self, count: usize, name: &str) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        self.rule
            .generate_rows_rational(&self.seeds, self.first_index, count)
            .map(|rows| rows.iter().map(|row| vec_poly::trim(row)).collect())
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!(
                    "the {name} recurrence cannot be replayed: {error}"
                ))
            })
    }
}

/// Targets `m >= from` where the leader `f(m, .)` is the zero polynomial.
fn leader_vanishing_targets(denominator: Option<&BivarPoly>, from: usize) -> Vec<usize> {
    let Some(denominator) = denominator else {
        return Vec::new();
    };
    let leader = Biv::from_bivar(denominator);
    let Some(degree) = leader.y_degree() else {
        return vec![from];
    };
    let columns: Vec<Vec<Q>> = (0..=degree).map(|j| leader.y_coefficient(j)).collect();
    let Some(first) = columns.iter().find(|c| !c.is_empty()) else {
        return vec![from];
    };
    roots_in_domain(first, from, None)
        .into_iter()
        .filter(|m| {
            let value = Q::from_integer((*m).into());
            columns
                .iter()
                .all(|c| vec_poly::evaluate(c, &value).is_zero())
        })
        .collect()
}

fn prepare_rule(
    rule: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    name: &str,
) -> Result<PreparedRule, RecurrenceError> {
    if rule.terms.is_empty() {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has no terms"
        )));
    }
    if let Some(position) = rule.terms.iter().position(|t| t.offset == 0) {
        return Err(RecurrenceError::InvalidInput(format!(
            "{name} recurrence term {position} has offset 0; offsets must be positive"
        )));
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has a (-1)^n factor"
        )));
    }
    if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
        return Err(RecurrenceError::Unsupported(format!(
            "the {name} recurrence has an inhomogeneous term"
        )));
    }
    if let Some(term) = rule
        .terms
        .iter()
        .find(|t| t.deriv_order > MAX_VERONESE_CERTIFY_DERIVATIVE_ORDER)
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {name} recurrence has derivative order {} > {MAX_VERONESE_CERTIFY_DERIVATIVE_ORDER}",
            term.deriv_order
        )));
    }
    if rule.denominator.as_ref().is_some_and(BivarPoly::is_zero) {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {name} recurrence denominator is zero"
        )));
    }
    let order = rule.max_offset();
    if rows.len() < order {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {name} recurrence has order {order} but only {} initial rows were supplied",
            rows.len()
        )));
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    if count > MAX_VERONESE_CHECKED_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {name} recurrence needs {count} identifying rows, more than \
             {MAX_VERONESE_CHECKED_ROWS}"
        )));
    }
    let normalized: Vec<Vec<Q>> = rows.iter().map(|row| vec_poly::trim(row)).collect();
    let seeds = normalized[..count].to_vec();
    let generated = rule
        .generate_rows_rational(&seeds, first_index, normalized.len())
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the {name} recurrence cannot be replayed: {error}"
            ))
        })?;
    for (offset, (supplied, generated)) in normalized.iter().zip(&generated).enumerate() {
        if supplied != &vec_poly::trim(generated) {
            return Err(RecurrenceError::InvalidInput(format!(
                "supplied {name} row {} disagrees with the row generated by its recurrence",
                first_index + offset
            )));
        }
    }
    if let Some(target) =
        leader_vanishing_targets(rule.denominator.as_ref(), first_index + count).first()
    {
        return Err(RecurrenceError::InvalidInput(format!(
            "the {name} leader vanishes identically at index {target}, so the recurrence does \
             not define that row; supply the rows through index {target}"
        )));
    }
    let (relation, leader_shift) = ImplicitRelation::from_recurrence(rule)?;
    Ok(PreparedRule {
        rule: rule.clone(),
        first_index,
        seeds,
        relation: relation.op(),
        leader_shift,
    })
}

fn validate(
    candidate: &VeroneseSectionCandidate,
    options: &VeroneseSectionCertifyOptions,
) -> Result<(), RecurrenceError> {
    if !(2..=MAX_VERONESE_CERTIFY_MODULUS).contains(&candidate.modulus) {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the modulus must be between 2 and {MAX_VERONESE_CERTIFY_MODULUS}, got {}",
            candidate.modulus
        )));
    }
    if candidate.residue >= candidate.modulus {
        return Err(RecurrenceError::InvalidInput(format!(
            "the residue {} is not below the modulus {}",
            candidate.residue, candidate.modulus
        )));
    }
    if candidate.index_shift.unsigned_abs() > MAX_VERONESE_CERTIFY_INDEX_SHIFT as u64 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "|index_shift| must be at most {MAX_VERONESE_CERTIFY_INDEX_SHIFT}, got {}",
            candidate.index_shift
        )));
    }
    if options.checked_rows == 0 || options.checked_rows > MAX_VERONESE_CHECKED_ROWS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "checked_rows must be between 1 and {MAX_VERONESE_CHECKED_ROWS}, got {}",
            options.checked_rows
        )));
    }
    if options.max_operator_terms > MAX_VERONESE_OPERATOR_TERMS {
        return Err(RecurrenceError::LimitExceeded(format!(
            "max_operator_terms must be at most {MAX_VERONESE_OPERATOR_TERMS}, got {}",
            options.max_operator_terms
        )));
    }
    Ok(())
}

/// The input rows and the candidate's section rows `lambda(n) V_j(P_(n+s))`
/// for `n0 <= n < n0 + count`.
struct ComparedRows {
    input: Vec<Vec<Q>>,
    candidate: Vec<Vec<Q>>,
}

impl ComparedRows {
    fn first_mismatch(&self) -> Option<usize> {
        self.input
            .iter()
            .zip(&self.candidate)
            .position(|(a, b)| a != b)
    }
}

/// Certify the candidate Veronese-section lift `Q_n = lambda(n) V_(r, j)(P_(n+s))`
/// of the sequence `Q` defined by `input` and `input_rows` (`input_rows[0]`
/// is `Q_(first_index)`, ascending coefficients in `y`), and try to prove
/// the parent statement.
///
/// The input may have index-dependent coefficients, derivatives on earlier
/// rows and a denominator `f(m, y)`; its leading row carries no derivative
/// (a [`Recurrence`] always has this form).  The parent must have a leader
/// free of `x`.
///
/// # The identity
///
/// Write `n0 = first_index`, `K` for the identifying rows
/// ([`Recurrence::generation_initial_count`]), `R` for the input order, `M`
/// for the input rule as a relation at base index `n`
/// ([`ImplicitRelation::from_recurrence`]: `Q` satisfies `M` at every base
/// `n >= n0 + K - R`), `p0`, `Kp` and `k` for the same data of the parent,
/// and `P~_n = lambda(n) P_(n+s)`, `F~` for its section vector.  The gauged
/// parent relation `L~` (leader `g~(n) S^k`, free of `x`) holds at every
/// base `n >= b = max(n0, p0 + Kp - k - s)`.
///
/// 1. The first `checked_rows` rows (at least `K`) are compared exactly; a
///    difference is [`VeroneseSectionIdentityStatus::Refuted`] with its index.
/// 2. `e_j^T M` is pseudo-divided by `phi(L~)`:
///    `w(n) e_j^T M = sum_nu H_nu o S^nu o phi(L~) + R`.  A nonzero remainder
///    gives [`VeroneseSectionIdentityStatus::SupportedByRows`] (evidence
///    only).
/// 3. With `R = 0`: applying the identity to `F~` gives
///    `w(n) (M V_j(P~))_n = sum_nu H_nu (phi(L~) F~)_(n+nu) = 0` for every
///    `n >= b`, because `phi(L~) F~` are the sections of `L~ P~` and
///    sections are unique.  `w` is a nonzero polynomial in `n` only, so
///    `(M V_j(P~))_n = 0` for every base `n >= max(b, n0 + K - R)` except the
///    finitely many integer roots of `w`, found exactly; at those roots and
///    at the bases `n0 + K - R <= n < b`, `M V_j(P~)` is evaluated exactly on
///    the rows (all of which are compared in step 1).  Hence `V_j(P~)`
///    satisfies `M` at every base `n >= n0 + K - R`, like `Q`.
/// 4. `S_n = Q_n - V_j(P~_n)` then satisfies `f(m, y) S_m = (terms in
///    S_(m-1), ..., S_(m-R))` for every target `m >= n0 + K`, where the leader
///    `f(m, .)` is not the zero polynomial (its common integer roots with
///    `m >= n0 + K` are excluded exactly when the input is validated), and
///    `S` vanishes at the `K` identifying indices.  By induction `S = 0`:
///    `Q_n = lambda(n) V_j(P_(n+s))` for every `n >= n0`
///    ([`VeroneseSectionIdentityStatus::Certified`]).
///
/// The gauge `lambda` is checked to be finite and nonzero at every
/// `n >= n0`, and the parent leader `g(m)` to be nonzero at every target
/// `m >= p0 + Kp` (otherwise the parent is not defined by its record).
///
/// # The parent and the verdict
///
/// The parent provers (`parent_proof`: periodic product, paper theorem,
/// root domain, attestation) decide whether every `P_m`, `m >= n0 + s`, is
/// zero or real-rooted with coefficients of one sign.  Then every section is
/// zero or real-rooted with nonpositive roots (the transfer theorem of
/// [`find_veronese_section_lifts`]), and so is `Q_n`, a nonzero multiple of
/// one.  The verdict is [`VeroneseSectionVerdict::Proved`] only for a
/// certified identity and a parent proved by one of the first three routes;
/// an attestation gives [`VeroneseSectionVerdict::ConditionalOnAttestation`].
/// Lean replay is not attempted.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for malformed or inconsistent input
/// (rows disagreeing with their rule, a vanishing gauge or leader,
/// `n0 + s < p0`); [`RecurrenceError::Unsupported`] for a rule outside the
/// class (alternating or inhomogeneous terms, a parent leader depending on
/// `x`); [`RecurrenceError::LimitExceeded`] for options or candidates beyond
/// the documented bounds (including [`MAX_VERONESE_PARENT_OFFSET`]);
/// [`RecurrenceError::BudgetExhausted`] when the operator-term or checked-row
/// budget runs out;
/// [`RecurrenceError::VerificationFailed`] if an internal self-check fails.
///
/// # Example
///
/// ```
/// use polytool::recurrence::{
///     certify_veronese_section_lift, BigRational, BivarPoly, Recurrence, RecurrenceJson,
///     RecurrenceTerm, VeroneseSectionCandidate, VeroneseSectionCertifyOptions,
///     VeroneseSectionGauge, VeroneseSectionVerdict,
/// };
///
/// let q = |n: i64| BigRational::from_integer(n.into());
/// // Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n, Q_0 = Q_1 = 1: the even sections of (1 + x)^n.
/// let input = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(2)]] }),
///         RecurrenceTerm::new(2, 0, BivarPoly { coeffs: vec![vec![q(-1), q(1)]] }),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// // The parent P_n = (1 + x) P_(n-1), P_0 = 1.
/// let parent_rule = Recurrence {
///     terms: vec![RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(1), q(1)]] })],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let candidate = VeroneseSectionCandidate {
///     modulus: 2,
///     residue: 0,
///     parent: RecurrenceJson::from_recurrence_rational(&parent_rule, 0, &[vec![q(1)]], None),
///     index_shift: 0,
///     gauge: VeroneseSectionGauge::None,
/// };
/// let certification = certify_veronese_section_lift(
///     &input,
///     &[vec![q(1)], vec![q(1)]],
///     0,
///     &candidate,
///     &VeroneseSectionCertifyOptions::default(),
/// )?;
/// assert_eq!(certification.verdict, VeroneseSectionVerdict::Proved);
/// certification.verify()?;
/// # Ok::<(), polytool::recurrence::RecurrenceError>(())
/// ```
pub fn certify_veronese_section_lift(
    input: &Recurrence,
    input_rows: &[Vec<BigRational>],
    first_index: usize,
    candidate: &VeroneseSectionCandidate,
    options: &VeroneseSectionCertifyOptions,
) -> Result<VeroneseSectionCertification, RecurrenceError> {
    validate(candidate, options)?;
    let (modulus, residue, shift) = (candidate.modulus, candidate.residue, candidate.index_shift);
    let input = prepare_rule(input, input_rows, first_index, "input")?;
    let n0 = first_index;
    let identifying = input.seeds.len();

    let (parent_rule, parent_first, parent_rows) = candidate.parent.to_recurrence_parts()?;
    let parent = prepare_rule(&parent_rule, &parent_rows, parent_first, "parent")?;
    let (parent_k, _) = x_free_leader(&parent.relation)?;
    let parent_start = i64::try_from(n0)
        .ok()
        .and_then(|n0| n0.checked_add(shift))
        .filter(|start| *start >= 0 && usize::try_from(*start).is_ok_and(|s| s >= parent_first))
        .ok_or_else(|| {
            RecurrenceError::InvalidInput(format!(
                "the candidate uses parent row n0 + s = {n0} + {shift}, before the parent's \
                 first index {parent_first}"
            ))
        })? as usize;
    if parent_start - parent_first > MAX_VERONESE_PARENT_OFFSET {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the first used parent row {parent_start} is more than \
             {MAX_VERONESE_PARENT_OFFSET} rows after the parent's first index {parent_first}"
        )));
    }
    let gauge = candidate.gauge.ratio_form(n0)?;

    let compare = |count: usize| -> Result<ComparedRows, RecurrenceError> {
        let input_rows = input.rows(count, "input")?;
        let parent_rows = parent.rows(parent_start + count - parent_first, "parent")?;
        let lambdas = gauge.values(n0, count);
        let candidate_rows = (0..count)
            .map(|offset| {
                let section = section_row(
                    &parent_rows[parent_start - parent_first + offset],
                    modulus,
                    residue,
                );
                vec_poly::trim(&vec_poly::scale(&section, &lambdas[offset]))
            })
            .collect();
        Ok(ComparedRows {
            input: input_rows,
            candidate: candidate_rows,
        })
    };

    let mut notes = Vec::new();
    let mut rows_compared = options.checked_rows.max(identifying);
    let mut rows = compare(rows_compared)?;
    if let Some(offset) = rows.first_mismatch() {
        return finish(
            candidate,
            options,
            &parent,
            (n0, identifying, rows_compared),
            VeroneseSectionIdentityStatus::Refuted { index: n0 + offset },
            notes,
        );
    }

    // Step 2: pseudo-division of e_j^T M by phi(L~).
    let gauged = gauge_parent_relation(&parent.relation, shift, &gauge)?;
    let (k, leader) = x_free_leader(&gauged)?;
    let system = transport_to_sections(&gauged, modulus, options.max_operator_terms)?;
    let row = residue_row(&input.relation, modulus, residue);
    let division = pseudo_divide_row(&row, &system, k, &leader, options.max_operator_terms)?;
    if !super::section_certify::division_identity_holds(&row, &system, &division) {
        return Err(RecurrenceError::VerificationFailed(
            "the pseudo-division identity does not hold".into(),
        ));
    }
    if !division.remainder_is_zero() {
        notes.push(format!(
            "the pseudo-division left a nonzero remainder ({} terms of shift < {k}); the \
             identity is supported by the compared rows only",
            division
                .remainder
                .iter()
                .map(|op| op.0.len())
                .sum::<usize>()
        ));
        return finish(
            candidate,
            options,
            &parent,
            (n0, identifying, rows_compared),
            VeroneseSectionIdentityStatus::SupportedByRows {
                rows: rows_compared,
            },
            notes,
        );
    }

    // Step 3: the exceptional bases.
    let input_from = n0 + identifying - input.leader_shift;
    let parent_from = {
        let base = i64::try_from(parent.first_index + parent.seeds.len() - parent_k)
            .map_err(|_| RecurrenceError::LimitExceeded("index overflow".into()))?
            - shift;
        usize::try_from(base.max(0)).unwrap_or(0).max(n0)
    };
    let multiplier = n_coefficients(&division.multiplier);
    let mut exceptional: Vec<usize> = (input_from..parent_from).collect();
    exceptional.extend(roots_in_domain(
        &multiplier,
        input_from.max(parent_from),
        None,
    ));
    exceptional.sort_unstable();
    exceptional.dedup();
    if let Some(last) = exceptional.last() {
        let needed = last + input.leader_shift + 1 - n0;
        if needed > MAX_VERONESE_CHECKED_ROWS {
            return Err(RecurrenceError::BudgetExhausted(format!(
                "the exceptional base {last} needs {needed} compared rows, more than \
                 {MAX_VERONESE_CHECKED_ROWS}"
            )));
        }
        if needed > rows_compared {
            rows_compared = needed;
            rows = compare(rows_compared)?;
            if let Some(offset) = rows.first_mismatch() {
                return finish(
                    candidate,
                    options,
                    &parent,
                    (n0, identifying, rows_compared),
                    VeroneseSectionIdentityStatus::Refuted { index: n0 + offset },
                    notes,
                );
            }
        }
    }
    for base in &exceptional {
        let value = input
            .relation
            .apply(&rows.candidate, n0, *base)
            .ok_or_else(|| {
                RecurrenceError::VerificationFailed(format!(
                    "rows missing at exceptional base {base}"
                ))
            })?;
        if !value.iter().all(Zero::is_zero) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the input relation does not vanish on the candidate rows at base {base}"
            )));
        }
    }

    let certificate = VeroneseForwardDivisionCertificate {
        modulus,
        residue,
        input_relation: ImplicitRelation::from_op(&input.relation).to_json(),
        input_leader_shift: input.leader_shift,
        parent_relation: ImplicitRelation::from_op(&parent.relation).to_json(),
        index_shift: shift,
        gauge: candidate.gauge.clone(),
        gauged_parent_relation: ImplicitRelation::from_op(&gauged).to_json(),
        parent_leader_shift: k,
        multiplier: encode_poly(&multiplier),
        quotients: VeroneseForwardDivisionCertificate::quotient_entries(&division),
        input_valid_from_base: input_from,
        parent_valid_from_base: parent_from,
        exceptional_bases: exceptional,
        identifying_indices: (n0..n0 + identifying).collect(),
        statement: format!(
            "w(n) e_{residue}^T M = sum_nu H_nu S^nu phi(L~) with remainder zero, so \
             w(n) (M V_({modulus},{residue})(P~))_n = 0 for every base n >= {parent_from}; with \
             the exceptional bases checked exactly, V_({modulus},{residue})(P~) satisfies the \
             input rule at every base n >= {input_from}, and agreement on the {identifying} \
             identifying rows gives Q_n = lambda(n) V_({modulus},{residue})(P_(n{shift:+})) \
             for every n >= {n0}"
        ),
    };
    certificate.verify()?;
    finish(
        candidate,
        options,
        &parent,
        (n0, identifying, rows_compared),
        VeroneseSectionIdentityStatus::Certified {
            kind: VeroneseIdentityCertificateKind::ForwardDivision,
            certificate: Box::new(certificate),
        },
        notes,
    )
}

fn finish(
    candidate: &VeroneseSectionCandidate,
    options: &VeroneseSectionCertifyOptions,
    parent: &PreparedRule,
    (n0, identifying, rows_compared): (usize, usize, usize),
    identity: VeroneseSectionIdentityStatus,
    notes: Vec<String>,
) -> Result<VeroneseSectionCertification, RecurrenceError> {
    let needed_from = (n0 as i64 + candidate.index_shift) as usize;
    let parent_proof = prove_parent(
        &ParentRecord {
            rule: &parent.rule,
            seeds: &parent.seeds,
            first_index: parent.first_index,
        },
        needed_from,
        &options.attestations,
        options.try_root_domain,
    )?;
    let verdict = verdict(&identity, &parent_proof);
    let (r, j, s) = (candidate.modulus, candidate.residue, candidate.index_shift);
    let section = format!("lambda(n) V_({r},{j})(P_(n{s:+}))");
    let statement = match verdict {
        VeroneseSectionVerdict::Proved => format!(
            "proved: Q_n = {section} for every n >= {n0} and {}; so every Q_n, n >= {n0}, is \
             zero or real-rooted with nonpositive roots",
            parent_proof.statement
        ),
        VeroneseSectionVerdict::ConditionalOnAttestation => format!(
            "conditional: Q_n = {section} for every n >= {n0} is certified; if the attestation \
             {:?} holds, every Q_n is zero or real-rooted with nonpositive roots",
            parent_proof.conditional_on
        ),
        VeroneseSectionVerdict::IdentityCertifiedParentOpen => format!(
            "Q_n = {section} for every n >= {n0} is certified; the parent statement is open"
        ),
        VeroneseSectionVerdict::IdentityCertifiedParentRefuted => format!(
            "Q_n = {section} for every n >= {n0} is certified, but {}",
            parent_proof.statement
        ),
        VeroneseSectionVerdict::IdentitySupportedByRows => {
            format!("Q_n = {section} holds on the {rows_compared} compared rows only (not proved)")
        }
        VeroneseSectionVerdict::IdentityRefuted => match &identity {
            VeroneseSectionIdentityStatus::Refuted { index } => {
                format!("refuted: Q_{index} differs from the candidate {section} at n = {index}")
            }
            _ => unreachable!("the verdict follows the identity status"),
        },
    };
    let mut certification = VeroneseSectionCertification {
        schema: VERONESE_SECTION_CERTIFICATE_SCHEMA.into(),
        candidate: candidate.clone(),
        input_first_index: n0,
        identifying_rows: identifying,
        rows_compared,
        identity,
        parent_proof,
        verdict,
        statement,
        lean_replay: LEAN_REPLAY_NOT_ATTEMPTED.into(),
        notes,
        lift_summary: None,
    };
    certification.lift_summary = Some(certification.recompute_lift_summary());
    Ok(certification)
}
