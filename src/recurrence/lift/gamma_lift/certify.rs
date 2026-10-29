//! The candidate, the certificate and the certification of one gamma lift
//! `Q_n = lambda(n) Gamma_(n+s)(c z)` (module documentation of `gamma_lift`).

use crate::recurrence::implicit::algebra::{roots_in_domain, Biv, DOp};
use crate::recurrence::lift::hadamard_factor::prepare::{
    left_multiply, prepare_rule, shift_free_from_json, shift_free_json, PreparedRule,
};
use crate::recurrence::lift::lift_summary::{gamma_lift_summary, GammaLiftSummaryParts};
use crate::recurrence::lift::parent_proof::{prove_parent, ParentRecord};
use crate::recurrence::lift::section_certify::{
    gauge_parent_relation, n_coefficients, n_polynomial, pseudo_divide_row, x_free_leader,
    RatioGauge,
};
use crate::recurrence::lift::section_operators::OperatorMatrix;
use crate::recurrence::root_domain::bivariate_leader_nonvanishing;
use crate::recurrence::{
    canonical_recurrence, decode_poly, encode_poly, format_rational_coeff, parse_rational_coeff,
    transport_gamma, BivarPoly, DegreeSchedule, ImplicitRelation, ImplicitRelationJson,
    LiftSummary, OperatorTermJson, Recurrence, RecurrenceError, RecurrenceJson,
    RecurrenceTransformation, TransportedRecurrence, VeroneseParentAttestation,
    VeroneseParentProof, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// JSON schema tag of a serialized [`GammaLiftCertification`].
pub const GAMMA_LIFT_CERTIFICATION_SCHEMA: &str = "polytool.recurrence-gamma-lift.v1";

/// Largest `|s|` of the index shift of a candidate.
pub const MAX_GAMMA_LIFT_INDEX_SHIFT: usize = 4;

/// Largest number of rows compared exactly (also the bound on the rows that
/// the exceptional bases of a certificate may require).
pub const MAX_GAMMA_LIFT_CHECKED_ROWS: usize = 64;

/// Largest number of operator terms in a pseudo-division.
pub const MAX_GAMMA_LIFT_OPERATOR_TERMS: usize = 10_000;

/// Largest distance `n0 + s - p0` between the first used parent row and the
/// parent's first index.
const MAX_GAMMA_LIFT_PARENT_OFFSET: usize = 1024;

/// A parent family `P_m`, symmetric at index `d_m` for every `m >= p0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftParent {
    /// A short name (for reports).
    pub name: String,
    /// The parent rule with its first index `p0` and seeds
    /// (`polytool.recurrence.v1`): homogeneous, without `(-1)^n` factors.
    pub recurrence: RecurrenceJson,
    /// The symmetry indices `d_m` (`x^(d_m) P_m(1/x) = P_m`); unbounded, with
    /// `first_index <= p0`.
    pub symmetry_schedule: DegreeSchedule,
}

/// A candidate lift `Q_n = lambda(n) Gamma_(n+s)(c z)` for every `n >= n0`,
/// with `Gamma_m` the gamma polynomial of `P_m` at index `d_m`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftCandidate {
    /// The parent.
    pub parent: GammaLiftParent,
    /// The index shift `s`, `|s| <=` [`MAX_GAMMA_LIFT_INDEX_SHIFT`] and
    /// `n0 + s >= p0`.  Default 0.
    #[serde(default)]
    pub index_shift: i64,
    /// The scale `c != 0` (exact rational string).
    pub scale: String,
    /// The normalisation `lambda`.  Default none (`lambda = 1`).
    #[serde(default = "default_gauge")]
    pub gauge: VeroneseSectionGauge,
}

fn default_gauge() -> VeroneseSectionGauge {
    VeroneseSectionGauge::None
}

/// Options of [`certify_gamma_lift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GammaLiftOptions {
    /// Rows `Q_(n0), ...` compared exactly with the candidate (at least the
    /// identifying rows; more when the exceptional bases need them).  At most
    /// [`MAX_GAMMA_LIFT_CHECKED_ROWS`].  Default 16.
    pub checked_rows: usize,
    /// Operator-term budget, at most [`MAX_GAMMA_LIFT_OPERATOR_TERMS`].
    pub max_operator_terms: usize,
    /// Caller-supplied parent attestations (the attestation route of the
    /// parent provers).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain parent route.  Default true.
    pub try_root_domain: bool,
}

impl Default for GammaLiftOptions {
    fn default() -> Self {
        Self {
            checked_rows: 16,
            max_operator_terms: MAX_GAMMA_LIFT_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
        }
    }
}

/// Whether every parent row `P_m`, `m >= p0`, is symmetric at `d_m`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum GammaParentSymmetryStatus {
    /// Proved for every `m >= p0` by the reversal induction (module
    /// documentation, step 1).
    Certified {
        /// The argument with its index ranges.
        detail: String,
    },
    /// The compared rows are symmetric, but the induction does not apply.
    SupportedByRows {
        /// Number of parent rows compared, from `p0`.
        rows: usize,
        /// Why the induction does not apply.
        reason: String,
    },
    /// A parent row is not symmetric at its index.
    Refuted {
        /// The absolute index of the first asymmetric row.
        index: usize,
    },
}

impl GammaParentSymmetryStatus {
    /// Whether the symmetry is proved for every parent row.
    pub fn is_certified(&self) -> bool {
        matches!(self, Self::Certified { .. })
    }
}

/// One quotient `H_nu` of a [`GammaLiftCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaDivisionQuotient {
    /// The shift `nu` of `H_nu o S^nu o L`.
    pub shift: usize,
    /// `H_nu = sum_d h_d(n, z) D^d`.
    pub differentiations: Vec<OperatorTermJson>,
}

/// The exact identity behind a certified gamma lift (module documentation):
/// the parent symmetry induction, the gamma relation `T`, the gauged and
/// scaled relation `L` of the candidate rows, and
/// `w(n) M = sum_nu H_nu S^nu L` with remainder zero.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftCertificate {
    /// The input: rule `M`, first index `n0` and the `K` identifying seeds.
    pub input: RecurrenceJson,
    /// The parent: rule, first index `p0` and the `K_p` identifying seeds.
    pub parent: RecurrenceJson,
    /// The symmetry indices `d_m` of the parent rows.
    pub symmetry_schedule: DegreeSchedule,
    /// The relation `T` of the gamma rows `G_m` at base index (from
    /// [`transport_gamma`]).
    pub gamma_relation: ImplicitRelationJson,
    /// The scale `c`.
    pub scale: String,
    /// The index shift `s`.
    pub index_shift: i64,
    /// The normalisation `lambda`.
    pub gauge: VeroneseSectionGauge,
    /// `L`: the relation of `H_n = lambda(n) G_(n+s)(c z)`, leader `g(n) S^k`.
    pub candidate_relation: ImplicitRelationJson,
    /// `k`.
    pub candidate_leader_shift: usize,
    /// `w(n)`, ascending exact strings in `n`.
    pub multiplier: Vec<String>,
    /// The nonzero quotients `H_nu`.
    pub quotients: Vec<GammaDivisionQuotient>,
    /// First base at which `M` holds for the input rows (`n0 + K - R`).
    pub input_valid_from_base: usize,
    /// First base at which `L` holds for the candidate rows.
    pub candidate_valid_from_base: usize,
    /// Bases from `input_valid_from_base` on where the identity gives nothing
    /// (before `candidate_valid_from_base`, integer roots of `w`); `M` is
    /// evaluated on the candidate rows there and vanishes.
    pub exceptional_bases: Vec<usize>,
    /// The identifying input indices `n0, ..., n0 + K - 1`.
    pub identifying_indices: Vec<usize>,
    /// What the certificate proves.
    pub statement: String,
}

/// Status of the identity `Q_n = lambda(n) Gamma_(n+s)(c z)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum GammaLiftIdentityStatus {
    /// Proved for every `n >= n0`.
    Certified {
        /// The replayable certificate.
        certificate: Box<GammaLiftCertificate>,
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

/// `polytool.recurrence-gamma-lift.v1`: the certification of one candidate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GammaLiftCertification {
    /// Always [`GAMMA_LIFT_CERTIFICATION_SCHEMA`].
    pub schema: String,
    /// The candidate as supplied.
    pub candidate: GammaLiftCandidate,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number `K` of identifying input rows.
    pub identifying_rows: usize,
    /// Number of rows `Q_(n0), ...` compared exactly with the candidate.
    pub rows_compared: usize,
    /// The symmetry of the parent rows.
    pub parent_symmetry: GammaParentSymmetryStatus,
    /// The transported rule of the gamma rows `G_m` (variable `t` is `z`),
    /// when the gamma transport succeeded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gamma_rule: Option<String>,
    /// The identity status, with the certificate when certified.
    pub identity: GammaLiftIdentityStatus,
    /// The parent statement (every parent row used is zero or real-rooted
    /// with coefficients of one sign) and its proof attempts.
    pub parent_proof: VeroneseParentProof,
    /// The uniform summary; its verdict is the verdict of the lift.
    pub lift_summary: LiftSummary,
    /// What the verdict means for the input rows.
    pub statement: String,
    /// Further remarks (division remainder, limits).
    pub notes: Vec<String>,
}

impl GammaLiftCertification {
    /// Whether every `Q_n`, `n >= n0`, is proved zero or real-rooted in the
    /// summary's input root region.
    pub fn is_proved(&self) -> bool {
        self.lift_summary.is_proved()
    }

    /// The [`LiftSummary`] computed from the record.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        gamma_lift_summary(GammaLiftSummaryParts {
            candidate: &self.candidate,
            input_first_index: self.input_first_index,
            parent_symmetry: &self.parent_symmetry,
            identity: &self.identity,
            parent_proof: &self.parent_proof,
        })
    }

    /// Replay the certificate (when certified), check that it belongs to the
    /// candidate, and that the summary equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if self.schema != GAMMA_LIFT_CERTIFICATION_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        if let GammaLiftIdentityStatus::Certified { certificate } = &self.identity {
            certificate.verify()?;
            let candidate = &self.candidate;
            if certificate.parent != canonical_parent_record(&candidate.parent.recurrence)?
                || certificate.symmetry_schedule != candidate.parent.symmetry_schedule
                || certificate.index_shift != candidate.index_shift
                || certificate.gauge != candidate.gauge
                || parse_scale(&certificate.scale)? != parse_scale(&candidate.scale)?
                || certificate.identifying_indices.first() != Some(&self.input_first_index)
            {
                return failed("the certificate does not belong to the candidate");
            }
            if !self.parent_symmetry.is_certified() {
                return failed("a certified identity needs certified parent symmetry");
            }
        }
        let summary = self.recompute_lift_summary();
        if summary != self.lift_summary || !summary.is_consistent() {
            return failed("the lift summary does not match the one recomputed from the record");
        }
        Ok(())
    }
}

fn parse_scale(text: &str) -> Result<Q, RecurrenceError> {
    let scale = parse_rational_coeff(text)
        .map_err(|error| RecurrenceError::from(error).with_context("scale"))?;
    if scale.is_zero() {
        return Err(RecurrenceError::InvalidInput(
            "the scale c must be nonzero".into(),
        ));
    }
    Ok(scale)
}

/// The parent record reduced to its identifying seeds.
fn canonical_parent_record(record: &RecurrenceJson) -> Result<RecurrenceJson, RecurrenceError> {
    let (rule, first_index, rows) = record.to_recurrence_parts()?;
    let parent = prepare_rule(&rule, &rows, first_index, "parent")?;
    Ok(RecurrenceJson::from_recurrence_rational(
        &parent.rule,
        parent.first_index,
        &parent.seeds,
        None,
    ))
}

// ---------------------------------------------------------------------------
// The algebra shared by certification and replay
// ---------------------------------------------------------------------------

/// The parent symmetry induction (module documentation, step 1) from the
/// rule and the seeds alone; `Err` gives the reason it does not apply.
fn symmetry_induction(parent: &PreparedRule, schedule: &DegreeSchedule) -> Result<String, String> {
    let p0 = parent.first_index;
    let valid_from = p0 + parent.seeds.len();
    let reversal = RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    };
    let reversed = match reversal.transport(&parent.rule) {
        Ok(TransportedRecurrence::Scalar(rule)) => rule,
        Ok(TransportedRecurrence::Vector(_)) => {
            return Err("the reversal transport is not scalar".into())
        }
        Err(error) => return Err(format!("the reversal transport failed: {error}")),
    };
    if !crate::recurrence::transform::rules_proportional(&reversed, &parent.rule) {
        return Err(
            "the reversal-transported rule does not define the same row map as the parent rule"
                .into(),
        );
    }
    let one = BivarPoly {
        coeffs: vec![vec![Q::one()]],
    };
    for (what, rule) in [("parent", &parent.rule), ("reversed", &reversed)] {
        let leader = canonical_recurrence(rule)
            .denominator
            .unwrap_or(one.clone());
        if !bivariate_leader_nonvanishing(&leader, valid_from) {
            return Err(format!(
                "the {what} leader is not proved nonzero at every target >= {valid_from}"
            ));
        }
    }
    for (offset, seed) in parent.seeds.iter().enumerate() {
        let index = p0 + offset;
        if !row_is_symmetric(seed, schedule, index) {
            return Err(format!("the seed P_{index} is not symmetric at d_{index}"));
        }
    }
    Ok(format!(
        "the reversal transport x^(d_m) P_m(1/x) of the parent rule defines the same row map, \
         both left factors are nonzero polynomials at every target m >= {valid_from}, and the \
         seeds P_{p0}..P_{} are symmetric; so the reversed rows coincide with the rows by \
         induction and every P_m, m >= {p0}, is symmetric at d_m",
        valid_from - 1
    ))
}

fn row_is_symmetric(row: &[Q], schedule: &DegreeSchedule, index: usize) -> bool {
    let row = vec_poly::trim(row);
    if vec_poly::is_zero(&row) {
        return true;
    }
    match schedule.degree_at(index) {
        Ok(d) => vec_poly::reverse(&row, d).is_some_and(|reversed| reversed == row),
        Err(_) => false,
    }
}

/// The relation `T` of the gamma rows: the transported scalar rule, else a
/// part of the transported identity with a nonzero leader.
fn gamma_rule(
    parent: &Recurrence,
    schedule: &DegreeSchedule,
) -> Result<Recurrence, RecurrenceError> {
    let transport = transport_gamma(parent, schedule)?;
    if let Some(scalar) = transport.scalar {
        return Ok(scalar);
    }
    for part in [&transport.even, &transport.odd] {
        if !part.leader_vanishes {
            return Ok(part.relation.clone());
        }
    }
    Err(RecurrenceError::Unsupported(
        "no part of the gamma-transported identity has a nonzero leader".into(),
    ))
}

/// `T`, then `L` (the scaled, shifted and gauged relation of the candidate
/// rows) with its leader shift.
struct CandidateRelation {
    gamma: DOp,
    gamma_rule: Recurrence,
    candidate: DOp,
    leader_shift: usize,
}

fn candidate_relation(
    parent: &Recurrence,
    schedule: &DegreeSchedule,
    scale: &Q,
    index_shift: i64,
    gauge: &RatioGauge,
) -> Result<CandidateRelation, RecurrenceError> {
    let gamma_rule = gamma_rule(parent, schedule)?;
    let (gamma, _) = ImplicitRelation::from_recurrence(&gamma_rule)?;
    let affine = RecurrenceTransformation::Affine {
        a: format_rational_coeff(scale),
        b: "0".into(),
    };
    let scaled = match affine.transport(&gamma_rule)? {
        TransportedRecurrence::Scalar(rule) => rule,
        TransportedRecurrence::Vector(_) => {
            return Err(RecurrenceError::Unsupported(
                "the scaled gamma rule is not scalar".into(),
            ))
        }
    };
    let (scaled, _) = ImplicitRelation::from_recurrence(&scaled)?;
    let candidate = gauge_parent_relation(&scaled.op(), index_shift, gauge)?;
    let (leader_shift, _) = x_free_leader(&candidate)?;
    Ok(CandidateRelation {
        gamma: gamma.op(),
        gamma_rule,
        candidate,
        leader_shift,
    })
}

/// `w M - sum_nu H_nu o S^nu o L` vanishes.
fn division_identity_holds(
    input: &DOp,
    candidate: &DOp,
    multiplier: &Biv,
    quotients: &BTreeMap<usize, DOp>,
) -> bool {
    let mut difference = left_multiply(multiplier, input);
    for (nu, h) in quotients {
        difference = difference.sub(&h.compose(&candidate.shifted(*nu)));
    }
    difference.is_zero()
}

/// The candidate rows `H_n = lambda(n) G_(n+s)(c z)`, `n = n0, ..., n0 + count - 1`;
/// `Err(m)` when the parent row `P_m` is not symmetric at `d_m`.
struct CandidateRows<'a> {
    parent: &'a PreparedRule,
    schedule: &'a DegreeSchedule,
    scale: &'a Q,
    gauge: &'a RatioGauge,
    first_index: usize,
    parent_start: usize,
}

enum RowsError {
    Asymmetric(usize),
    Failed(RecurrenceError),
}

impl CandidateRows<'_> {
    fn rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RowsError> {
        let offset = self.parent_start - self.parent.first_index;
        let parent_rows = self
            .parent
            .rows(offset + count)
            .map_err(RowsError::Failed)?;
        let gamma = RecurrenceTransformation::Gamma {
            schedule: self.schedule.clone(),
        };
        let lambdas = self.gauge.values(self.first_index, count);
        (0..count)
            .map(|i| {
                let m = self.parent_start + i;
                let row = &parent_rows[offset + i];
                if !row_is_symmetric(row, self.schedule, m) {
                    return Err(RowsError::Asymmetric(m));
                }
                let mut g = gamma
                    .map_row(m, row)
                    .map_err(|error| RowsError::Failed(error.into()))?
                    .swap_remove(0);
                let mut power = lambdas[i].clone();
                for coefficient in &mut g {
                    *coefficient *= &power;
                    power *= self.scale;
                }
                Ok(vec_poly::trim(&g))
            })
            .collect()
    }
}

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

impl GammaLiftCertificate {
    /// Replay the certificate from its recorded input, parent, schedule,
    /// scale, shift and gauge: the parent symmetry induction, the gamma
    /// relation, the candidate relation, the division identity, the
    /// exceptional bases (every integer root of `w` from the valid bases on
    /// is listed) with `M` evaluated on the candidate rows there, and the
    /// identifying rows.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let (input_rule, n0, input_rows) = self.input.to_recurrence_parts()?;
        let input = prepare_rule(&input_rule, &input_rows, n0, "input")?;
        let (parent_rule, p0, parent_rows) = self.parent.to_recurrence_parts()?;
        let parent = prepare_rule(&parent_rule, &parent_rows, p0, "parent")?;
        validate_schedule(&self.symmetry_schedule, p0)?;
        if let Err(reason) = symmetry_induction(&parent, &self.symmetry_schedule) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the parent symmetry induction fails: {reason}"
            )));
        }
        let scale = parse_scale(&self.scale)?;
        let gauge = self.gauge.ratio_form(n0)?;
        let parent_start = parent_start(n0, self.index_shift, p0)?;
        let relation = candidate_relation(
            &parent.rule,
            &self.symmetry_schedule,
            &scale,
            self.index_shift,
            &gauge,
        )?;
        if ImplicitRelation::from_op(&relation.gamma).to_json() != self.gamma_relation {
            return failed("the gamma relation does not match its recomputation");
        }
        if ImplicitRelation::from_op(&relation.candidate).to_json() != self.candidate_relation
            || relation.leader_shift != self.candidate_leader_shift
        {
            return failed("the candidate relation does not match its recomputation");
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
        if !division_identity_holds(
            &input.relation,
            &relation.candidate,
            &multiplier,
            &quotients,
        ) {
            return failed("the identity w M = sum H_nu S^nu L fails");
        }
        let input_from = input.valid_from_base();
        let candidate_from =
            candidate_valid_from(n0, &parent, relation.leader_shift, self.index_shift);
        if input_from != self.input_valid_from_base
            || candidate_from != self.candidate_valid_from_base
        {
            return failed("the recorded valid bases do not match");
        }
        let expected = exceptional_bases(&multiplier, input_from, candidate_from);
        if expected
            .iter()
            .any(|base| !self.exceptional_bases.contains(base))
        {
            return failed("an exceptional base is missing");
        }
        if self.identifying_indices != (n0..n0 + input.seeds.len()).collect::<Vec<_>>() {
            return failed("the identifying indices do not match the input");
        }
        if self.exceptional_bases.iter().any(|base| *base < n0) {
            return failed("an exceptional base precedes the input's first index");
        }
        let needed = self
            .exceptional_bases
            .iter()
            .map(|base| base + input.leader_shift + 1 - n0)
            .max()
            .unwrap_or(0)
            .max(input.seeds.len());
        if needed > MAX_GAMMA_LIFT_CHECKED_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "the certificate needs {needed} candidate rows, more than \
                 {MAX_GAMMA_LIFT_CHECKED_ROWS}"
            )));
        }
        let rows = CandidateRows {
            parent: &parent,
            schedule: &self.symmetry_schedule,
            scale: &scale,
            gauge: &gauge,
            first_index: n0,
            parent_start,
        }
        .rows(needed)
        .map_err(|error| match error {
            RowsError::Asymmetric(m) => RecurrenceError::VerificationFailed(format!(
                "the parent row P_{m} is not symmetric"
            )),
            RowsError::Failed(error) => error,
        })?;
        for (offset, seed) in input.seeds.iter().enumerate() {
            if vec_poly::trim(seed) != rows[offset] {
                return failed("an identifying row differs from the candidate row");
            }
        }
        for base in &self.exceptional_bases {
            let value = input.relation.apply(&rows, n0, *base).ok_or_else(|| {
                RecurrenceError::VerificationFailed(format!(
                    "rows missing at exceptional base {base}"
                ))
            })?;
            if !value.iter().all(Zero::is_zero) {
                return failed("the input relation does not vanish at an exceptional base");
            }
        }
        Ok(())
    }
}

fn validate_schedule(
    schedule: &DegreeSchedule,
    parent_first: usize,
) -> Result<(), RecurrenceError> {
    schedule.validate()?;
    if schedule.last_index.is_some() || schedule.first_index > parent_first {
        return Err(RecurrenceError::InvalidInput(format!(
            "the symmetry schedule must be unbounded and start at or before the parent's first \
             index {parent_first}"
        )));
    }
    Ok(())
}

/// `n0 + s`, which must be at least `p0`.
fn parent_start(n0: usize, shift: i64, p0: usize) -> Result<usize, RecurrenceError> {
    if shift.unsigned_abs() > MAX_GAMMA_LIFT_INDEX_SHIFT as u64 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "|index_shift| = {} exceeds {MAX_GAMMA_LIFT_INDEX_SHIFT}",
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
    if start - p0 > MAX_GAMMA_LIFT_PARENT_OFFSET {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the first used parent row {start} is more than {MAX_GAMMA_LIFT_PARENT_OFFSET} rows \
             after the parent's first index {p0}"
        )));
    }
    Ok(start)
}

/// `max(n0, p0 + K_p - k - s)`: from there on `L` holds for the candidate
/// rows (module documentation, step 3).
fn candidate_valid_from(n0: usize, parent: &PreparedRule, k: usize, shift: i64) -> usize {
    let base = (parent.first_index + parent.seeds.len()) as i64 - k as i64 - shift;
    usize::try_from(base.max(0)).unwrap_or(0).max(n0)
}

fn exceptional_bases(multiplier: &Biv, input_from: usize, candidate_from: usize) -> Vec<usize> {
    let mut bases: Vec<usize> = (input_from..candidate_from).collect();
    bases.extend(roots_in_domain(
        &n_coefficients(multiplier),
        input_from.max(candidate_from),
        None,
    ));
    bases.sort_unstable();
    bases.dedup();
    bases
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// Certify one candidate gamma lift `Q_n = lambda(n) Gamma_(n+s)(c z)` of the
/// input rows (`input_rows[0]` is `Q_(first_index)`; at least the
/// identifying rows) and run the parent provers.
///
/// The identity is certified by the argument of the `gamma_lift` module
/// documentation (parent symmetry induction, gamma transport, pseudo-division
/// of the input rule by the candidate relation); the parent statement (every
/// used parent row zero or real-rooted with coefficients of one sign) by the
/// shared parent provers.  The verdict is that of the
/// [`GammaLiftCertification::lift_summary`]: `proved` exactly when both hold
/// and `lambda` is finite and nonzero on the domain; then every `Q_n`,
/// `n >= n0`, is zero or real-rooted with roots in `(-inf, 0]` (`c > 0`) or
/// `[0, inf)` (`c < 0`).  A verdict other than `proved` says nothing against
/// real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError`] for malformed input (rules that do not reproduce their
/// rows, a leader not proved nonzero, a zero scale, an invalid schedule or
/// gauge, a shift before the parent's first index) and over the limits.
///
/// ```
/// use polytool::recurrence::*;
/// use polytool::sequences::literature::peak_polynomials_recurrence;
///
/// // The peak polynomials W_n = 2^(n-1) Gamma(A_n)(z/4), n >= 1.
/// let (rule, seeds, first_index) = peak_polynomials_recurrence();
/// let candidate = GammaLiftCandidate {
///     parent: GammaLiftParent::eulerian(),
///     index_shift: 0,
///     scale: "1/4".into(),
///     gauge: VeroneseSectionGauge::Ratio {
///         numerator: vec!["2".into()],
///         denominator: vec!["1".into()],
///         initial_value: "1".into(),
///     },
/// };
/// let certification =
///     certify_gamma_lift(&rule, &seeds, first_index, &candidate, &GammaLiftOptions::default())?;
/// assert!(certification.is_proved());
/// certification.verify()?;
/// # Ok::<(), RecurrenceError>(())
/// ```
pub fn certify_gamma_lift(
    input: &Recurrence,
    input_rows: &[Vec<BigRational>],
    first_index: usize,
    candidate: &GammaLiftCandidate,
    options: &GammaLiftOptions,
) -> Result<GammaLiftCertification, RecurrenceError> {
    if options.checked_rows > MAX_GAMMA_LIFT_CHECKED_ROWS
        || options.max_operator_terms > MAX_GAMMA_LIFT_OPERATOR_TERMS
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "checked_rows <= {MAX_GAMMA_LIFT_CHECKED_ROWS} and max_operator_terms <= \
             {MAX_GAMMA_LIFT_OPERATOR_TERMS} are required"
        )));
    }
    let n0 = first_index;
    let input = prepare_rule(input, input_rows, n0, "input")?;
    let identifying = input.seeds.len();
    let (parent_rule, p0, parent_rows) = candidate.parent.recurrence.to_recurrence_parts()?;
    let parent = prepare_rule(&parent_rule, &parent_rows, p0, "parent")?;
    let schedule = &candidate.parent.symmetry_schedule;
    validate_schedule(schedule, p0)?;
    let scale = parse_scale(&candidate.scale)?;
    let gauge = candidate.gauge.ratio_form(n0)?;
    let shift = candidate.index_shift;
    let start = parent_start(n0, shift, p0)?;
    let candidate_rows = CandidateRows {
        parent: &parent,
        schedule,
        scale: &scale,
        gauge: &gauge,
        first_index: n0,
        parent_start: start,
    };

    let mut notes = Vec::new();
    let mut rows_compared = options.checked_rows.max(identifying);
    let mut parent_symmetry = match symmetry_induction(&parent, schedule) {
        Ok(detail) => GammaParentSymmetryStatus::Certified { detail },
        Err(reason) => GammaParentSymmetryStatus::SupportedByRows {
            rows: rows_compared,
            reason,
        },
    };
    let gamma_rule = gamma_rule(&parent.rule, schedule)
        .ok()
        .map(|rule| rule.to_string());
    let context = FinishContext {
        candidate,
        options,
        parent: &parent,
        n0,
        identifying,
    };

    // Step 1: compare rows exactly.
    let compare = |count: usize| -> Result<Compared, RecurrenceError> {
        let inputs = input.rows(count)?;
        match candidate_rows.rows(count) {
            Ok(rows) => Ok(match first_mismatch(&inputs, &rows) {
                Some(offset) => Compared::Mismatch(n0 + offset),
                None => Compared::Rows(rows),
            }),
            Err(RowsError::Asymmetric(m)) => Ok(Compared::Asymmetric(m)),
            Err(RowsError::Failed(error)) => Err(error),
        }
    };
    let refuted = |outcome: Compared, symmetry: &mut GammaParentSymmetryStatus| match outcome {
        Compared::Mismatch(index) => GammaLiftIdentityStatus::Refuted {
            index,
            reason: format!("Q_{index} differs from the candidate lambda(n) Gamma_(n+s)(c z)"),
        },
        Compared::Asymmetric(m) => {
            *symmetry = GammaParentSymmetryStatus::Refuted { index: m };
            GammaLiftIdentityStatus::Refuted {
                index: (m as i64 - shift) as usize,
                reason: format!(
                    "the parent row P_{m} is not symmetric at d_{m}, so it has no gamma \
                     polynomial"
                ),
            }
        }
        Compared::Rows(_) => unreachable!("rows are not a refutation"),
    };
    let mut rows = match compare(rows_compared)? {
        Compared::Rows(rows) => rows,
        outcome => {
            let identity = refuted(outcome, &mut parent_symmetry);
            return finish(
                &context,
                parent_symmetry,
                gamma_rule,
                identity,
                rows_compared,
                notes,
            );
        }
    };
    if !parent_symmetry.is_certified() {
        let reason = match &parent_symmetry {
            GammaParentSymmetryStatus::SupportedByRows { reason, .. } => reason.clone(),
            _ => String::new(),
        };
        notes.push(format!(
            "the parent symmetry is not proved for every row: {reason}"
        ));
        let identity = GammaLiftIdentityStatus::SupportedByRows {
            rows: rows_compared,
            reason: "the parent symmetry induction does not apply".into(),
        };
        return finish(
            &context,
            parent_symmetry,
            gamma_rule,
            identity,
            rows_compared,
            notes,
        );
    }

    // Steps 2-4: the candidate relation and the division.
    let relation = match candidate_relation(&parent.rule, schedule, &scale, shift, &gauge) {
        Ok(relation) => relation,
        Err(error) => {
            let identity = GammaLiftIdentityStatus::SupportedByRows {
                rows: rows_compared,
                reason: format!("no candidate relation: {error}"),
            };
            return finish(
                &context,
                parent_symmetry,
                gamma_rule,
                identity,
                rows_compared,
                notes,
            );
        }
    };
    let k = relation.leader_shift;
    let (_, leader) = x_free_leader(&relation.candidate)?;
    let system = OperatorMatrix {
        entries: vec![vec![relation.candidate.clone()]],
    };
    let division = pseudo_divide_row(
        std::slice::from_ref(&input.relation),
        &system,
        k,
        &leader,
        options.max_operator_terms,
    )?;
    if !division.remainder_is_zero() {
        let identity = GammaLiftIdentityStatus::SupportedByRows {
            rows: rows_compared,
            reason: format!(
                "the pseudo-division of the input rule by the candidate relation left a nonzero \
                 remainder of shift < {k}"
            ),
        };
        return finish(
            &context,
            parent_symmetry,
            gamma_rule,
            identity,
            rows_compared,
            notes,
        );
    }
    let quotients: BTreeMap<usize, DOp> = division
        .quotients
        .iter()
        .map(|(nu, row)| (*nu, row[0].clone()))
        .filter(|(_, h)| !h.is_zero())
        .collect();
    if !division_identity_holds(
        &input.relation,
        &relation.candidate,
        &division.multiplier,
        &quotients,
    ) {
        return Err(RecurrenceError::VerificationFailed(
            "the pseudo-division identity does not hold".into(),
        ));
    }

    let input_from = input.valid_from_base();
    let candidate_from = candidate_valid_from(n0, &parent, k, shift);
    let exceptional = exceptional_bases(&division.multiplier, input_from, candidate_from);
    if let Some(last) = exceptional.last() {
        let needed = last + input.leader_shift + 1 - n0;
        if needed > MAX_GAMMA_LIFT_CHECKED_ROWS {
            return Err(RecurrenceError::BudgetExhausted(format!(
                "the exceptional base {last} needs {needed} compared rows, more than \
                 {MAX_GAMMA_LIFT_CHECKED_ROWS}"
            )));
        }
        if needed > rows_compared {
            rows_compared = needed;
            rows = match compare(rows_compared)? {
                Compared::Rows(rows) => rows,
                outcome => {
                    let identity = refuted(outcome, &mut parent_symmetry);
                    return finish(
                        &context,
                        parent_symmetry,
                        gamma_rule,
                        identity,
                        rows_compared,
                        notes,
                    );
                }
            };
        }
    }
    for base in &exceptional {
        let value = input.relation.apply(&rows, n0, *base).ok_or_else(|| {
            RecurrenceError::VerificationFailed(format!("rows missing at exceptional base {base}"))
        })?;
        if !value.iter().all(Zero::is_zero) {
            let identity = GammaLiftIdentityStatus::SupportedByRows {
                rows: rows_compared,
                reason: format!(
                    "the input rule does not vanish on the candidate rows at the exceptional \
                     base {base}"
                ),
            };
            return finish(
                &context,
                parent_symmetry,
                gamma_rule,
                identity,
                rows_compared,
                notes,
            );
        }
    }

    let multiplier = n_coefficients(&division.multiplier);
    let certificate = GammaLiftCertificate {
        input: RecurrenceJson::from_recurrence_rational(&input.rule, n0, &input.seeds, None),
        parent: RecurrenceJson::from_recurrence_rational(&parent.rule, p0, &parent.seeds, None),
        symmetry_schedule: schedule.clone(),
        gamma_relation: ImplicitRelation::from_op(&relation.gamma).to_json(),
        scale: format_rational_coeff(&scale),
        index_shift: shift,
        gauge: candidate.gauge.clone(),
        candidate_relation: ImplicitRelation::from_op(&relation.candidate).to_json(),
        candidate_leader_shift: k,
        multiplier: encode_poly(&multiplier),
        quotients: quotients
            .iter()
            .map(|(nu, h)| GammaDivisionQuotient {
                shift: *nu,
                differentiations: shift_free_json(h),
            })
            .collect(),
        input_valid_from_base: input_from,
        candidate_valid_from_base: candidate_from,
        exceptional_bases: exceptional,
        identifying_indices: (n0..n0 + identifying).collect(),
        statement: format!(
            "every parent row P_m, m >= {p0}, is symmetric at d_m; the gamma rows satisfy the \
             transported rule {}; w(n) M = sum_nu H_nu S^nu L with remainder zero for the \
             relation L of H_n = lambda(n) Gamma_(n{shift:+})(c z), c = {}, valid from base \
             {candidate_from}; with the exceptional bases checked exactly H satisfies the input \
             rule at every base n >= {input_from}, and agreement on the {identifying} \
             identifying rows gives Q_n = H_n for every n >= {n0}",
            relation.gamma_rule,
            format_rational_coeff(&scale)
        ),
    };
    certificate.verify()?;
    let identity = GammaLiftIdentityStatus::Certified {
        certificate: Box::new(certificate),
    };
    finish(
        &context,
        parent_symmetry,
        gamma_rule,
        identity,
        rows_compared,
        notes,
    )
}

/// The outcome of an exact comparison of input and candidate rows.
enum Compared {
    Rows(Vec<Vec<Q>>),
    /// The first differing input index.
    Mismatch(usize),
    /// The first parent index whose row is not symmetric.
    Asymmetric(usize),
}

fn first_mismatch(input: &[Vec<Q>], candidate: &[Vec<Q>]) -> Option<usize> {
    input
        .iter()
        .zip(candidate)
        .position(|(a, b)| vec_poly::trim(a) != vec_poly::trim(b))
}

struct FinishContext<'a> {
    candidate: &'a GammaLiftCandidate,
    options: &'a GammaLiftOptions,
    parent: &'a PreparedRule,
    n0: usize,
    identifying: usize,
}

fn finish(
    context: &FinishContext<'_>,
    parent_symmetry: GammaParentSymmetryStatus,
    gamma_rule: Option<String>,
    identity: GammaLiftIdentityStatus,
    rows_compared: usize,
    notes: Vec<String>,
) -> Result<GammaLiftCertification, RecurrenceError> {
    let candidate = context.candidate;
    let n0 = context.n0;
    let needed_from = (n0 as i64 + candidate.index_shift) as usize;
    let parent_proof = prove_parent(
        &ParentRecord {
            rule: &context.parent.rule,
            seeds: &context.parent.seeds,
            first_index: context.parent.first_index,
        },
        needed_from,
        &context.options.attestations,
        context.options.try_root_domain,
    )?;
    let lift_summary = gamma_lift_summary(GammaLiftSummaryParts {
        candidate,
        input_first_index: n0,
        parent_symmetry: &parent_symmetry,
        identity: &identity,
        parent_proof: &parent_proof,
    });
    let s = candidate.index_shift;
    let lift = format!("lambda(n) Gamma_(n{s:+})(c z), c = {}", candidate.scale);
    let statement = if lift_summary.is_proved() {
        format!(
            "proved: Q_n = {lift} for every n >= {n0} and {}; so every Q_n, n >= {n0}, is zero \
             or real-rooted with roots in {}",
            parent_proof.statement, lift_summary.root_transfer.input_root_region
        )
    } else {
        match &identity {
            GammaLiftIdentityStatus::Certified { .. } => format!(
                "Q_n = {lift} for every n >= {n0} is certified; the lift is not proved ({})",
                parent_proof.statement
            ),
            GammaLiftIdentityStatus::SupportedByRows { rows, reason } => format!(
                "Q_n = {lift} holds on the {rows} compared rows only (not proved: {reason})"
            ),
            GammaLiftIdentityStatus::Refuted { index, reason } => {
                format!("refuted at n = {index}: {reason}")
            }
        }
    };
    Ok(GammaLiftCertification {
        schema: GAMMA_LIFT_CERTIFICATION_SCHEMA.into(),
        candidate: candidate.clone(),
        input_first_index: n0,
        identifying_rows: context.identifying,
        rows_compared,
        parent_symmetry,
        gamma_rule,
        identity,
        parent_proof,
        lift_summary,
        statement,
        notes,
    })
}
