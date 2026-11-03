//! Certification of one basis-map lift `Q_n = T_P(R_n)`: the division
//! certificate of the identity, the root location of the coefficient
//! family and the [`LiftSummary`] (module documentation of `basis_map`,
//! "The lift certificate").

use super::super::hadamard_factor::prepare::{
    prepare_lift_rule, shift_free_from_json, shift_free_json, PreparedLiftRule,
};
use super::super::lift_summary::{fixed_preserver_interlacing, proved_parent_interlacing};
use super::super::parent_proof::{prove_parent_real_rootedness, ParentRecurrenceRecord};
use super::super::section_certify::{
    division_identity_holds, n_coefficients, pseudo_divide_row, x_free_leader, ForwardDivision,
};
use super::super::section_operators::OperatorMatrix;
use super::map::{apply_with_rows, BasisMap};
use super::quadratic::AlgebraicInterval;
use super::theorems::{
    check_basis_map_preserves_pf, BasisMapConclusion, BasisMapTheoremOptions, BasisMapTheoremReport,
};
use super::transport::transported_system;
use crate::real_rootedness::{isolate_real_roots_rational_coeffs, RealInterval};
use crate::recurrence::implicit::algebra::{roots_in_domain, Biv, DOp};
use crate::recurrence::{
    format_rational_coeff, parse_rational_coeff, LiftHypothesis, LiftHypothesisStatus,
    LiftIdentityStatus, LiftInterlacing, LiftInterlacingClaim, LiftOperation, LiftParent,
    LiftParentProof, LiftRootTransfer, LiftSummary, OperatorTermJson, Recurrence, RecurrenceError,
    RecurrenceJson, RecurrenceTransformation, RootInterval, TransportedRecurrence,
    VeroneseParentAttestation, VeroneseParentProof, VeroneseParentProofStatus,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// JSON schema tag of a serialized [`BasisMapLiftCertification`].
pub const BASIS_MAP_LIFT_SCHEMA: &str = "polytool.recurrence-basis-map-lift.v1";

/// Largest number of rows compared exactly (also the bound on the rows the
/// exceptional bases may require).
pub const MAX_BASIS_MAP_LIFT_CHECKED_ROWS: usize = 64;

/// Largest number of operator terms in a transported system or a division.
pub const MAX_BASIS_MAP_OPERATOR_TERMS: usize = 20_000;

/// Bits of the rational inner approximation of an irrational source end.
const INNER_BITS: u32 = 32;

/// Options of [`certify_basis_map_lift`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BasisMapLiftOptions {
    /// Rows `Q_(n0), ...` compared exactly with the candidate (at least the
    /// identifying rows), at most [`MAX_BASIS_MAP_LIFT_CHECKED_ROWS`].
    /// Default 16.
    pub checked_rows: usize,
    /// Operator-term budget, at most [`MAX_BASIS_MAP_OPERATOR_TERMS`].
    pub max_operator_terms: usize,
    /// Caller-supplied attestations for the coefficient family (the
    /// attestation route of the parent provers; a lift is then conditional,
    /// never proved).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain route of the parent provers.  Default true.
    pub try_root_domain: bool,
    /// Options of the map theorem check.
    pub theorem: BasisMapTheoremOptions,
}

impl Default for BasisMapLiftOptions {
    fn default() -> Self {
        Self {
            checked_rows: 16,
            max_operator_terms: MAX_BASIS_MAP_OPERATOR_TERMS,
            attestations: Vec::new(),
            try_root_domain: true,
            theorem: BasisMapTheoremOptions::default(),
        }
    }
}

/// Which rule is divided by which in a [`BasisMapDivisionCertificate`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BasisMapDivisionDirection {
    /// The input row `e_0^T M` is divided by the transported system `Phi`:
    /// the candidate rows `H_n = T(R_n)` satisfy the input rule `M`.
    InputByTransported,
    /// (Maps of order one.)  The transported rule `M' = Phi` is divided by
    /// the input rule `M`: the input rows satisfy `M'`.
    TransportedByInput,
}

/// One quotient `H_nu` (a row of shift-free operators, one per state
/// component) of a [`BasisMapDivisionCertificate`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapDivisionQuotient {
    /// The shift `nu` of `H_nu o S^nu o divisor`.
    pub shift: usize,
    /// `H_nu`: for each column, the terms `h_d(n, t) D^d`.
    pub row: Vec<Vec<OperatorTermJson>>,
}

/// The exact identity behind a certified basis-map lift `Q_n = T_P(R_n)`
/// for every `n >= n0` (module documentation of `basis_map`, "The lift
/// certificate"): `w(n) rho = sum_nu H_nu S^nu Delta` with remainder zero,
/// where for [`BasisMapDivisionDirection::InputByTransported`] `rho =
/// e_0^T M` and `Delta = Phi`, and for
/// [`BasisMapDivisionDirection::TransportedByInput`] `rho = Phi` and
/// `Delta = M`; plus the exceptional bases where the identity does not
/// apply, which are checked exactly on rows.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapDivisionCertificate {
    /// The input: rule `M`, first index `n0` and the identifying seeds.
    pub input: RecurrenceJson,
    /// The basis map.
    pub basis: BasisMap,
    /// The coefficient family `R`: rule, first index `n0` and identifying
    /// seeds.
    pub coefficient_family: RecurrenceJson,
    /// The direction of the division.
    pub direction: BasisMapDivisionDirection,
    /// `w(n)`, ascending exact coefficients in `n`.
    pub multiplier: Vec<String>,
    /// The quotients `H_nu`.
    pub quotients: Vec<BasisMapDivisionQuotient>,
    /// The bases `n` (absolute) where the divided rule is evaluated
    /// exactly on rows instead (below the first base of the divisor, or a
    /// root of `w`).
    pub exceptional_bases: Vec<usize>,
}

/// Status of the identity `Q_n = T_P(R_n)`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum BasisMapIdentityStatus {
    /// Proved for every `n >= n0`.
    Certified {
        /// The certificate.
        certificate: Box<BasisMapDivisionCertificate>,
    },
    /// The compared rows agree, but the identity is not proved.
    SupportedByRows {
        /// Rows compared from `n0`.
        rows: usize,
        /// Why the identity is not proved.
        reason: String,
    },
    /// A compared row differs: this candidate is wrong (nothing is said
    /// about the input).
    Refuted {
        /// The absolute index of the first differing row.
        index: usize,
    },
    /// The candidate leaves the supported class.
    Unsupported {
        /// Why.
        reason: String,
    },
}

/// The root location of the coefficient family: every `R_n` is zero or
/// real-rooted with every root in `proved_interval`, a rational interval
/// inside the map's source interval.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapCoefficientRootProof {
    /// The source interval of the map conclusion used.
    pub source_interval: AlgebraicInterval,
    /// The rational interval `[L', U']` the proof aims at (inside the
    /// source).
    pub target_interval: RootInterval,
    /// The parent provers on `R_n(v + U')` (every root `<= U'`).
    pub upper_bound: VeroneseParentProof,
    /// The parent provers on `R_n(L' - v)` (every root `>= L'`), for a
    /// finite lower end.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lower_bound: Option<VeroneseParentProof>,
    /// The parent statement and how it was settled.
    pub parent_proof: LiftParentProof,
}

/// A candidate `Q_n = T_P(R_n)`: the map and the coefficient family.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapLiftCandidate {
    /// The basis map.
    pub basis: BasisMap,
    /// The coefficient family `R` (rule in `u`, first index `n0`, seeds).
    pub coefficient_family: RecurrenceJson,
}

/// The certification of one candidate (schema
/// `polytool.recurrence-basis-map-lift.v1`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapLiftCertification {
    /// [`BASIS_MAP_LIFT_SCHEMA`].
    pub schema: String,
    /// The input: rule, first index `n0` and identifying seeds.
    pub input: RecurrenceJson,
    /// The candidate.
    pub candidate: BasisMapLiftCandidate,
    /// The identity status.
    pub identity: BasisMapIdentityStatus,
    /// The map theorem.
    pub map_theorem: BasisMapTheoremReport,
    /// The conclusion of the map theorem used (source and target).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conclusion: Option<BasisMapConclusion>,
    /// The root location of the coefficient family, when a conclusion was
    /// chosen.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coefficient_roots: Option<BasisMapCoefficientRootProof>,
    /// The uniform summary.
    pub lift_summary: LiftSummary,
}

impl BasisMapLiftCertification {
    /// Whether the summary verdict is proved.
    pub fn is_proved(&self) -> bool {
        self.lift_summary.is_proved()
    }

    /// Replay the certification: the schema and the summary's consistency,
    /// the division certificate ([`BasisMapDivisionCertificate::verify`])
    /// with the recorded input and candidate, the map theorem
    /// ([`BasisMapTheoremReport::verify`]), and a full recomputation with
    /// [`certify_basis_map_lift`] and default options from the recorded
    /// input (rows regenerated from its seeds) and candidate, whose identity
    /// status and summary must equal the recorded ones.  The summary's input
    /// root region may be a record without exact endpoints
    /// ([`RootInterval::agrees_with_recomputed`]); a recorded exact endpoint
    /// must be valid ([`RootInterval::check_exact_endpoints`]) and equal the
    /// recomputed one.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] for any mismatch, and the
    /// errors of the replays.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let fail = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if self.schema != BASIS_MAP_LIFT_SCHEMA {
            return fail("unsupported basis-map lift schema");
        }
        if !self.lift_summary.is_consistent() {
            return fail("the lift summary verdict is not that of its statuses");
        }
        if let BasisMapIdentityStatus::Certified { certificate } = &self.identity {
            if certificate.input != self.input
                || certificate.basis != self.candidate.basis
                || certificate.coefficient_family != self.candidate.coefficient_family
            {
                return fail("the certificate is about another input or candidate");
            }
            certificate.verify()?;
        }
        self.map_theorem.verify()?;
        let (rule, n0, seeds) = self.input.to_recurrence_parts()?;
        let rows = rule
            .generate_rows_rational(&seeds, n0, seeds.len() + 2)
            .map_err(|error| {
                RecurrenceError::VerificationFailed(format!(
                    "the input cannot be replayed: {error}"
                ))
            })?;
        let recomputed = certify_basis_map_lift(
            &rule,
            &rows,
            n0,
            &self.candidate,
            &BasisMapLiftOptions::default(),
        )?;
        // A record written before exact endpoints were recorded carries the
        // rational outer region only: compare it in that form (a recorded
        // exact form must equal the recomputed one).
        let recorded_region = &self.lift_summary.root_transfer.input_root_region;
        recorded_region
            .check_exact_endpoints()
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        let mut recomputed_summary = recomputed.lift_summary;
        let fresh_region = &mut recomputed_summary.root_transfer.input_root_region;
        if recorded_region.agrees_with_recomputed(fresh_region) {
            *fresh_region = recorded_region.clone();
        }
        if recomputed.identity != self.identity || recomputed_summary != self.lift_summary {
            return fail("the recomputed certification differs from the recorded one");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// The identity
// ---------------------------------------------------------------------------

/// The state row with `op` in column 0.
fn first_component_row(op: &DOp, size: usize) -> Vec<DOp> {
    let mut row = vec![DOp::default(); size];
    row[0] = op.clone();
    row
}

/// Exceptional bases: `from <= n < divisor_from`, and the integer roots of
/// `w` from `max(from, divisor_from)` on.
fn exceptional_bases(multiplier: &Biv, from: usize, divisor_from: usize) -> Vec<usize> {
    let mut bases: Vec<usize> = (from..divisor_from).collect();
    let w = n_coefficients(multiplier);
    if vec_poly::degree(&w).unwrap_or(0) > 0 {
        bases.extend(roots_in_domain(&w, from.max(divisor_from), None));
    }
    bases.sort_unstable();
    bases.dedup();
    bases
}

/// Candidate rows `H_n = T(R_n)` for `n0 <= n < n0 + count`.
fn candidate_rows(
    map: &BasisMap,
    coefficient: &PreparedLiftRule,
    count: usize,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    let r_rows = coefficient.rows(count)?;
    let degree = r_rows.iter().map(|r| r.len()).max().unwrap_or(1);
    let basis = map.basis_rows(degree)?;
    Ok(r_rows.iter().map(|r| apply_with_rows(r, &basis)).collect())
}

/// The internal outcome of the identity step.
struct IdentityOutcome {
    status: BasisMapIdentityStatus,
}

/// Certify `Q = T(R)` (module documentation).
fn certify_identity(
    input: &PreparedLiftRule,
    map: &BasisMap,
    coefficient: &PreparedLiftRule,
    options: &BasisMapLiftOptions,
) -> Result<IdentityOutcome, RecurrenceError> {
    let n0 = input.first_index;
    let checked = options
        .checked_rows
        .max(input.seeds.len())
        .max(coefficient.seeds.len());
    let q_rows = input.rows(checked)?;
    let h_rows = candidate_rows(map, coefficient, checked)?;
    if let Some(offset) = q_rows.iter().zip(&h_rows).position(|(q, h)| q != h) {
        return Ok(IdentityOutcome {
            status: BasisMapIdentityStatus::Refuted { index: n0 + offset },
        });
    }
    let supported = |reason: String| IdentityOutcome {
        status: BasisMapIdentityStatus::SupportedByRows {
            rows: checked,
            reason,
        },
    };
    let system = match transported_system(&coefficient.relation, map, options.max_operator_terms) {
        Ok(system) => system,
        Err(RecurrenceError::Unsupported(reason)) => {
            return Ok(IdentityOutcome {
                status: BasisMapIdentityStatus::Unsupported { reason },
            })
        }
        Err(RecurrenceError::LimitExceeded(reason)) => return Ok(supported(reason)),
        Err(error) => return Err(error),
    };
    let (k, leader) = match x_free_leader(&coefficient.relation) {
        Ok(found) => found,
        Err(RecurrenceError::Unsupported(reason)) => {
            return Ok(IdentityOutcome {
                status: BasisMapIdentityStatus::Unsupported { reason },
            })
        }
        Err(error) => return Err(error),
    };
    let mut reasons = Vec::new();
    // Direction A.
    let rho = first_component_row(&input.relation, system.size());
    match pseudo_divide_row(&rho, &system, k, &leader, options.max_operator_terms) {
        Ok(division) if division.remainder_is_zero() => {
            let certificate = build_certificate(
                input,
                map,
                coefficient,
                BasisMapDivisionDirection::InputByTransported,
                &division,
                input.valid_from_base(),
                coefficient.valid_from_base(),
            );
            match check_certificate(&certificate, input, map, coefficient, &system) {
                Ok(()) => {
                    return Ok(IdentityOutcome {
                        status: BasisMapIdentityStatus::Certified {
                            certificate: Box::new(certificate),
                        },
                    })
                }
                Err(RecurrenceError::VerificationFailed(reason)) => reasons.push(reason),
                Err(error) => return Err(error),
            }
        }
        Ok(_) => reasons.push(
            "the input rule is not a left multiple of the transported system (nonzero \
             remainder)"
                .to_string(),
        ),
        Err(RecurrenceError::LimitExceeded(reason))
        | Err(RecurrenceError::VerificationFailed(reason)) => reasons.push(reason),
        Err(error) => return Err(error),
    }
    // Direction B (maps of order one, input leader free of t).
    if system.size() == 1 {
        if let Ok((input_order, input_leader)) = x_free_leader(&input.relation) {
            let divisor = OperatorMatrix {
                entries: vec![vec![input.relation.clone()]],
            };
            let rho = vec![system.entries[0][0].clone()];
            match pseudo_divide_row(
                &rho,
                &divisor,
                input_order,
                &input_leader,
                options.max_operator_terms,
            ) {
                Ok(division) if division.remainder_is_zero() => {
                    let certificate = build_certificate(
                        input,
                        map,
                        coefficient,
                        BasisMapDivisionDirection::TransportedByInput,
                        &division,
                        coefficient.valid_from_base(),
                        input.valid_from_base(),
                    );
                    match check_certificate(&certificate, input, map, coefficient, &system) {
                        Ok(()) => {
                            return Ok(IdentityOutcome {
                                status: BasisMapIdentityStatus::Certified {
                                    certificate: Box::new(certificate),
                                },
                            })
                        }
                        Err(RecurrenceError::VerificationFailed(reason)) => reasons.push(reason),
                        Err(error) => return Err(error),
                    }
                }
                Ok(_) => reasons.push(
                    "the transported rule is not a left multiple of the input rule either"
                        .to_string(),
                ),
                Err(RecurrenceError::LimitExceeded(reason))
                | Err(RecurrenceError::VerificationFailed(reason)) => reasons.push(reason),
                Err(error) => return Err(error),
            }
        }
    }
    Ok(supported(reasons.join("; ")))
}

fn build_certificate(
    input: &PreparedLiftRule,
    map: &BasisMap,
    coefficient: &PreparedLiftRule,
    direction: BasisMapDivisionDirection,
    division: &ForwardDivision,
    from: usize,
    divisor_from: usize,
) -> BasisMapDivisionCertificate {
    BasisMapDivisionCertificate {
        input: RecurrenceJson::from_recurrence_rational(
            &input.rule,
            input.first_index,
            &input.seeds,
            None,
        ),
        basis: map.clone(),
        coefficient_family: RecurrenceJson::from_recurrence_rational(
            &coefficient.rule,
            coefficient.first_index,
            &coefficient.seeds,
            None,
        ),
        direction,
        multiplier: n_coefficients(&division.multiplier)
            .iter()
            .map(format_rational_coeff)
            .collect(),
        quotients: division
            .quotients
            .iter()
            .map(|(shift, row)| BasisMapDivisionQuotient {
                shift: *shift,
                row: row.iter().map(shift_free_json).collect(),
            })
            .collect(),
        exceptional_bases: exceptional_bases(&division.multiplier, from, divisor_from),
    }
}

/// The division identity, the exceptional bases and the identifying rows
/// of `certificate`, checked against the prepared rules.
fn check_certificate(
    certificate: &BasisMapDivisionCertificate,
    input: &PreparedLiftRule,
    map: &BasisMap,
    coefficient: &PreparedLiftRule,
    system: &OperatorMatrix,
) -> Result<(), RecurrenceError> {
    let fail = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    let multiplier: Vec<Q> = certificate
        .multiplier
        .iter()
        .map(|c| parse_rational_coeff(c))
        .collect::<Result<_, _>>()?;
    let multiplier = Biv::new(multiplier.into_iter().map(|c| vec![c]).collect());
    let mut quotients = BTreeMap::new();
    for quotient in &certificate.quotients {
        let row = quotient
            .row
            .iter()
            .map(|terms| shift_free_from_json(terms))
            .collect::<Result<Vec<_>, _>>()?;
        quotients.insert(quotient.shift, row);
    }
    let (rho, divisor, from, divisor_from) = match certificate.direction {
        BasisMapDivisionDirection::InputByTransported => (
            first_component_row(&input.relation, system.size()),
            system.clone(),
            input.valid_from_base(),
            coefficient.valid_from_base(),
        ),
        BasisMapDivisionDirection::TransportedByInput => {
            if system.size() != 1 {
                return fail("the reverse division needs a map of order one");
            }
            (
                vec![system.entries[0][0].clone()],
                OperatorMatrix {
                    entries: vec![vec![input.relation.clone()]],
                },
                coefficient.valid_from_base(),
                input.valid_from_base(),
            )
        }
    };
    if quotients
        .values()
        .any(|row: &Vec<DOp>| row.len() != rho.len())
    {
        return fail("a quotient row has the wrong length");
    }
    let division = ForwardDivision {
        multiplier: multiplier.clone(),
        quotients,
        remainder: vec![DOp::default(); rho.len()],
    };
    if !division_identity_holds(&rho, &divisor, &division) {
        return fail("the division identity w rho = sum H S^nu Delta does not hold");
    }
    if exceptional_bases(&multiplier, from, divisor_from) != certificate.exceptional_bases {
        return fail("the recorded exceptional bases differ from the recomputed ones");
    }
    // Rows: the divided rule must hold at the exceptional bases.
    let top = certificate
        .exceptional_bases
        .iter()
        .max()
        .copied()
        .unwrap_or(0);
    let order = rho.iter().map(DOp::max_shift).max().unwrap_or(0);
    let n0 = input.first_index;
    let count = (top + order + 1)
        .saturating_sub(n0)
        .max(input.seeds.len().max(coefficient.seeds.len()));
    if count > MAX_BASIS_MAP_LIFT_CHECKED_ROWS {
        return fail("the exceptional bases need more rows than the checked-row bound");
    }
    let q_rows = input.rows(count)?;
    let h_rows = candidate_rows(map, coefficient, count)?;
    let (rows, relation) = match certificate.direction {
        BasisMapDivisionDirection::InputByTransported => (&h_rows, &input.relation),
        BasisMapDivisionDirection::TransportedByInput => (&q_rows, &system.entries[0][0]),
    };
    for &n in &certificate.exceptional_bases {
        let value = relation
            .apply(rows, n0, n)
            .ok_or_else(|| RecurrenceError::VerificationFailed("missing row".into()))?;
        if !vec_poly::is_zero(&value) {
            return fail(&format!(
                "the divided rule fails at the exceptional base {n}"
            ));
        }
    }
    // Identifying rows.
    let identifying = match certificate.direction {
        BasisMapDivisionDirection::InputByTransported => input.seeds.len(),
        BasisMapDivisionDirection::TransportedByInput => coefficient.seeds.len(),
    };
    if q_rows[..identifying] != h_rows[..identifying] {
        return fail("the identifying rows differ");
    }
    Ok(())
}

impl BasisMapDivisionCertificate {
    /// Replay the certificate from the recorded rules and seeds: re-validate
    /// the input and the coefficient family (rows regenerated, leaders
    /// nonvanishing), recompute the transported system, check the division
    /// identity exactly, recompute the exceptional bases and evaluate the
    /// divided rule there on freshly computed rows, and compare the
    /// identifying rows.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] for any mismatch, and the
    /// errors of the replays.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let input = prepare_recorded(&self.input, "input")?;
        let coefficient = prepare_recorded(&self.coefficient_family, "coefficient family")?;
        if coefficient.first_index != input.first_index {
            return Err(RecurrenceError::VerificationFailed(
                "the coefficient family and the input start at different indices".into(),
            ));
        }
        let system = transported_system(
            &coefficient.relation,
            &self.basis,
            MAX_BASIS_MAP_OPERATOR_TERMS,
        )?;
        x_free_leader(&coefficient.relation)?;
        check_certificate(self, &input, &self.basis, &coefficient, &system)
    }
}

/// Prepare a recorded rule from its seeds (rows regenerated).
fn prepare_recorded(
    record: &RecurrenceJson,
    name: &str,
) -> Result<PreparedLiftRule, RecurrenceError> {
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    let count = seeds.len() + 2;
    let rows = rule
        .generate_rows_rational(&seeds, first_index, count)
        .map_err(|error| {
            RecurrenceError::VerificationFailed(format!("the {name} cannot be replayed: {error}"))
        })?;
    prepare_lift_rule(&rule, &rows, first_index, name)
}

// ---------------------------------------------------------------------------
// Root location of the coefficient family
// ---------------------------------------------------------------------------

/// Transport `rule` and `seeds` through `R(a v + b)`.
fn affine_family(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    a: &Q,
    b: &Q,
) -> Result<(Recurrence, Vec<Vec<Q>>), RecurrenceError> {
    let map = RecurrenceTransformation::Affine {
        a: format_rational_coeff(a),
        b: format_rational_coeff(b),
    };
    let TransportedRecurrence::Scalar(transported) = map.transport(rule)? else {
        return Err(RecurrenceError::Unsupported(
            "an affine transport did not give a scalar rule".into(),
        ));
    };
    let mapped = seeds
        .iter()
        .enumerate()
        .map(|(offset, row)| {
            map.map_row(first_index + offset, row)
                .map(|mut components| components.remove(0))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((transported, mapped))
}

fn prove_shifted(
    coefficient: &PreparedLiftRule,
    a: &Q,
    b: &Q,
    options: &BasisMapLiftOptions,
) -> Result<VeroneseParentProof, RecurrenceError> {
    let (rule, seeds) = affine_family(
        &coefficient.rule,
        coefficient.first_index,
        &coefficient.seeds,
        a,
        b,
    )?;
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: coefficient.first_index,
    };
    prove_parent_real_rootedness(
        &record,
        coefficient.first_index,
        &options.attestations,
        options.try_root_domain,
    )
}

/// Prove that every `R_n` is zero or real-rooted in a rational interval
/// inside `source` (module documentation, "Root location").
fn prove_coefficient_roots(
    coefficient: &PreparedLiftRule,
    r_rows: &[Vec<Q>],
    source: &AlgebraicInterval,
    options: &BasisMapLiftOptions,
) -> Result<BasisMapCoefficientRootProof, RecurrenceError> {
    let inner: RealInterval = source.rational_inner(INNER_BITS);
    let one = Q::from_integer(1.into());
    // The hull of the roots of the checked rows (rational enclosure ends).
    let mut smallest: Option<Q> = None;
    let mut largest: Option<Q> = None;
    for row in r_rows.iter().filter(|row| !vec_poly::is_zero(row)) {
        let isolation = isolate_real_roots_rational_coeffs(row);
        if let (Some(first), Some(last)) = (isolation.roots.first(), isolation.roots.last()) {
            if smallest.as_ref().is_none_or(|s| first.lower < *s) {
                smallest = Some(first.lower.clone());
            }
            if largest.as_ref().is_none_or(|l| last.upper > *l) {
                largest = Some(last.upper.clone());
            }
        }
    }
    // Upper end: 0 when the source reaches 0 (or beyond) and the checked
    // rows have no positive root, else a simple rational between the
    // largest root and the rational inner end U'.
    let rows_nonpositive = largest.as_ref().is_none_or(|l| !l.is_positive())
        || r_rows.iter().all(|row| {
            vec_poly::is_zero(row)
                || AlgebraicInterval::nonpositive_ray().contains_all_roots_of(row)
        });
    let upper_candidates: Vec<Q> = match &inner.upper {
        Some(u) if u.is_negative() || !rows_nonpositive => {
            let mut c = Vec::new();
            if let Some(l) = largest.as_ref().filter(|l| *l <= u) {
                c.push(simplest_rational_between(l, u));
            }
            c.push(u.clone());
            c
        }
        _ => vec![Q::zero()],
    };
    let lower_candidates: Vec<Q> = match &inner.lower {
        Some(l) => {
            let mut c = Vec::new();
            if let Some(s) = smallest.as_ref().filter(|s| *s >= l) {
                c.push(simplest_rational_between(l, s));
            }
            c.push(l.clone());
            c
        }
        None => Vec::new(),
    };
    let (upper, upper_bound) = first_proved(&upper_candidates, |b| {
        prove_shifted(coefficient, &one, b, options)
    })?;
    let (lower, lower_bound) = if lower_candidates.is_empty() {
        (None, None)
    } else {
        let (l, proof) = first_proved(&lower_candidates, |b| {
            prove_shifted(coefficient, &-one.clone(), b, options)
        })?;
        (Some(l), Some(proof))
    };
    let target = RealInterval {
        lower_closed: lower.is_some(),
        lower,
        upper: Some(upper),
        upper_closed: true,
    };
    let statement = format!(
        "every coefficient row R_n, n >= {}, is zero or real-rooted with every root in {}",
        coefficient.first_index, target
    );
    let proofs: Vec<&VeroneseParentProof> = std::iter::once(&upper_bound)
        .chain(lower_bound.as_ref())
        .collect();
    let parent_proof = if let Some(failed) = proofs
        .iter()
        .find(|p| p.status == VeroneseParentProofStatus::Refuted)
    {
        LiftParentProof::Failed {
            index: failed.failing_row_index,
            statement: format!(
                "{statement} fails: a shifted coefficient row is not real-rooted with one-signed \
                 coefficients ({})",
                failed.statement
            ),
        }
    } else if proofs
        .iter()
        .all(|p| p.status == VeroneseParentProofStatus::Proved)
    {
        let routes: Vec<String> = proofs
            .iter()
            .map(|p| {
                p.route
                    .and_then(|r| {
                        serde_json::to_value(r)
                            .ok()
                            .and_then(|v| v.as_str().map(String::from))
                    })
                    .unwrap_or_else(|| "unknown".into())
            })
            .collect();
        LiftParentProof::Proved {
            route: format!("basis_map_coefficient_roots({})", routes.join(", ")),
            statement,
        }
    } else if proofs.iter().all(|p| {
        matches!(
            p.status,
            VeroneseParentProofStatus::Proved | VeroneseParentProofStatus::ConditionalOnAttestation
        )
    }) {
        LiftParentProof::Attested {
            labels: proofs
                .iter()
                .flat_map(|p| p.conditional_on.iter().cloned())
                .collect(),
            statement,
        }
    } else {
        LiftParentProof::Open { statement }
    };
    Ok(BasisMapCoefficientRootProof {
        source_interval: source.clone(),
        target_interval: RootInterval::from_real(&target),
        upper_bound,
        lower_bound,
        parent_proof,
    })
}

/// Try the candidate bounds in order; the first proved one wins, else the
/// last attempt is returned.
fn first_proved(
    candidates: &[Q],
    mut prove: impl FnMut(&Q) -> Result<VeroneseParentProof, RecurrenceError>,
) -> Result<(Q, VeroneseParentProof), RecurrenceError> {
    let mut last = None;
    for (i, candidate) in candidates.iter().enumerate() {
        if i > 0 && candidates[..i].contains(candidate) {
            continue;
        }
        let proof = prove(candidate)?;
        let proved = proof.status == VeroneseParentProofStatus::Proved;
        last = Some((candidate.clone(), proof));
        if proved {
            break;
        }
    }
    last.ok_or_else(|| RecurrenceError::InvalidInput("no candidate bound".into()))
}

/// The rational with the smallest denominator, then the smallest absolute
/// numerator, in `[a, b]` (`a <= b`; continued fractions).
fn simplest_rational_between(a: &Q, b: &Q) -> Q {
    if a > b {
        return simplest_rational_between(b, a);
    }
    if !a.is_positive() && !b.is_negative() {
        return Q::zero();
    }
    if b.is_negative() {
        return -simplest_rational_between(&-b.clone(), &-a.clone());
    }
    let floor = a.floor();
    if floor == *a {
        return a.clone();
    }
    let ceil = &floor + Q::from_integer(1.into());
    if ceil <= *b {
        return ceil;
    }
    // a, b in (floor, floor + 1): invert the fractional parts.
    let inner = simplest_rational_between(&(b - &floor).recip(), &(a - &floor).recip());
    floor + inner.recip()
}

/// The conclusions in the order of preference: PF statements, then
/// bounded targets, then the rest.
fn ordered_conclusions(report: &BasisMapTheoremReport) -> Vec<BasisMapConclusion> {
    let mut conclusions = report.conclusions.clone();
    conclusions.sort_by_key(|c| {
        (
            !c.is_pf_conclusion(),
            c.output_interval.lower.is_none() || c.output_interval.upper.is_none(),
        )
    });
    conclusions
}

// ---------------------------------------------------------------------------
// Certification
// ---------------------------------------------------------------------------

/// Certify one candidate `Q_n = T_P(R_n)` for every `n >= n0`, where `Q`
/// is the input (`recurrence` with `rows`, `rows[0] = Q_(n0)`), `P` the
/// basis of `candidate.basis` and `R` the coefficient family (module
/// documentation of `basis_map`, "The lift certificate" and "Root
/// location").
///
/// The verdict of the returned [`LiftSummary`] is `proved` exactly when the
/// identity is certified, the map theorem gives a conclusion whose source
/// contains the proved root interval of every `R_n`, and that root
/// interval is proved by the parent provers; then every `Q_n` is zero or
/// real-rooted with all roots in the conclusion's target.
///
/// # Errors
///
/// [`RecurrenceError`] for an input or coefficient rule that does not
/// reproduce its rows or whose leader is not proved nonvanishing, for
/// families starting at different indices, and for options beyond their
/// bounds.
pub fn certify_basis_map_lift(
    recurrence: &Recurrence,
    rows: &[Vec<Q>],
    n0: usize,
    candidate: &BasisMapLiftCandidate,
    options: &BasisMapLiftOptions,
) -> Result<BasisMapLiftCertification, RecurrenceError> {
    if !(1..=MAX_BASIS_MAP_LIFT_CHECKED_ROWS).contains(&options.checked_rows)
        || options.max_operator_terms > MAX_BASIS_MAP_OPERATOR_TERMS
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "checked_rows must be between 1 and {MAX_BASIS_MAP_LIFT_CHECKED_ROWS} and \
             max_operator_terms at most {MAX_BASIS_MAP_OPERATOR_TERMS}"
        )));
    }
    let input = prepare_lift_rule(recurrence, rows, n0, "input")?;
    let (r_rule, r_first, r_seeds) = candidate.coefficient_family.to_recurrence_parts()?;
    if r_first != n0 {
        return Err(RecurrenceError::InvalidInput(format!(
            "the coefficient family starts at {r_first}, the input at {n0}"
        )));
    }
    let r_rows = r_rule
        .generate_rows_rational(&r_seeds, n0, r_seeds.len() + 2)
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the coefficient family cannot be replayed: {error}"
            ))
        })?;
    let coefficient = prepare_lift_rule(&r_rule, &r_rows, n0, "coefficient family")?;
    certify_prepared(&input, &coefficient, candidate, options)
}

pub(super) fn certify_prepared(
    input: &PreparedLiftRule,
    coefficient: &PreparedLiftRule,
    candidate: &BasisMapLiftCandidate,
    options: &BasisMapLiftOptions,
) -> Result<BasisMapLiftCertification, RecurrenceError> {
    let map = &candidate.basis;
    let n0 = input.first_index;
    let identity = certify_identity(input, map, coefficient, options)?.status;
    let map_theorem = check_basis_map_preserves_pf(map, &options.theorem)?;
    // Choose the conclusion whose source contains the checked rows' roots.
    let checked = options.checked_rows.max(coefficient.seeds.len());
    let r_rows = coefficient.rows(checked)?;
    let mut chosen = None;
    let mut root_failure = None;
    for conclusion in ordered_conclusions(&map_theorem) {
        let failing = r_rows.iter().position(|row| {
            !vec_poly::is_zero(row) && !conclusion.input_interval.contains_all_roots_of(row)
        });
        match failing {
            None => {
                chosen = Some(conclusion);
                break;
            }
            Some(offset) => {
                if root_failure.is_none() {
                    root_failure = Some((n0 + offset, conclusion.input_interval.clone()));
                }
            }
        }
    }
    let coefficient_roots = match &chosen {
        Some(conclusion) => Some(prove_coefficient_roots(
            coefficient,
            &r_rows,
            &conclusion.input_interval,
            options,
        )?),
        None => None,
    };
    let parent_proof = match (&coefficient_roots, &root_failure) {
        (Some(proof), _) => proof.parent_proof.clone(),
        (None, Some((index, source))) => LiftParentProof::Failed {
            index: Some(*index),
            statement: format!(
                "the coefficient row R_{index} has a root outside the source interval {source} \
                 of the map theorem (checked exactly); this route gives no proof"
            ),
        },
        (None, None) => LiftParentProof::Open {
            statement: "the map theorem gives no conclusion for this basis".into(),
        },
    };
    let summary = summary(
        input,
        coefficient,
        map,
        &identity,
        &map_theorem,
        chosen.as_ref(),
        coefficient_roots.as_ref(),
        parent_proof,
    );
    Ok(BasisMapLiftCertification {
        schema: BASIS_MAP_LIFT_SCHEMA.into(),
        input: RecurrenceJson::from_recurrence_rational(&input.rule, n0, &input.seeds, None),
        candidate: candidate.clone(),
        identity,
        map_theorem,
        conclusion: chosen,
        coefficient_roots,
        lift_summary: summary,
    })
}

#[allow(clippy::too_many_arguments)]
fn summary(
    input: &PreparedLiftRule,
    coefficient: &PreparedLiftRule,
    map: &BasisMap,
    identity: &BasisMapIdentityStatus,
    report: &BasisMapTheoremReport,
    conclusion: Option<&BasisMapConclusion>,
    coefficient_roots: Option<&BasisMapCoefficientRootProof>,
    parent_proof: LiftParentProof,
) -> LiftSummary {
    let parent_region = coefficient_roots.map(|p| p.target_interval.clone());
    let operation = LiftOperation::BasisMap {
        basis: map.name().to_string(),
        basis_family: map.recurrence_json(),
        source_interval: conclusion
            .map(|c| c.input_interval.to_string())
            .unwrap_or_default(),
        target_interval: conclusion
            .map(|c| c.output_interval.to_string())
            .unwrap_or_default(),
    };
    let parent = LiftParent {
        recurrence: RecurrenceJson::from_recurrence_rational(
            &coefficient.rule,
            coefficient.first_index,
            &coefficient.seeds,
            None,
        ),
        relation: format!(
            "Q_n = T_P(R_n), T_P(u^k) = P_({} + k), P the basis family {}",
            map.first_index(),
            map.name()
        ),
        gauge: None,
        index_shift: 0,
    };
    let identity = match identity {
        BasisMapIdentityStatus::Certified { certificate } => LiftIdentityStatus::Certified {
            method: match certificate.direction {
                BasisMapDivisionDirection::InputByTransported => {
                    "basis_map_division_input_by_transported".into()
                }
                BasisMapDivisionDirection::TransportedByInput => {
                    "basis_map_division_transported_by_input".into()
                }
            },
        },
        BasisMapIdentityStatus::SupportedByRows { rows, .. } => {
            LiftIdentityStatus::SupportedByRows { rows: *rows }
        }
        BasisMapIdentityStatus::Refuted { index } => LiftIdentityStatus::Refuted {
            index: Some(*index),
            reason: "a candidate row T_P(R_n) differs from the input row".into(),
        },
        BasisMapIdentityStatus::Unsupported { reason } => LiftIdentityStatus::Unsupported {
            reason: reason.clone(),
        },
    };
    let hypotheses: Vec<LiftHypothesis> = match conclusion {
        Some(c) => report
            .hypotheses
            .iter()
            .filter(|h| c.hypotheses.contains(&h.name))
            .cloned()
            .collect(),
        None => Vec::new(),
    };
    let root_transfer = LiftRootTransfer {
        theorem: report.theorem.clone(),
        parent_hypothesis: "every coefficient row R_n is zero or real-rooted with every root in \
                            the source interval of the map theorem"
            .into(),
        parent_root_region: match &parent_proof {
            LiftParentProof::Proved { .. } | LiftParentProof::Attested { .. } => parent_region,
            _ => None,
        },
        input_root_region: conclusion
            .map(|c| c.output_interval.to_root_interval())
            .unwrap_or_else(|| AlgebraicInterval::whole_line().to_root_interval()),
        input_conclusion: match conclusion {
            Some(c) => format!(
                "every input row Q_n, n >= {}, is zero or real-rooted with every root in {}{}",
                input.first_index,
                c.output_interval,
                if c.nonnegative_coefficients {
                    " (nonnegative coefficients when the coefficient row has a positive leading \
                     coefficient)"
                } else {
                    ""
                }
            ),
            None => "no conclusion".into(),
        },
    };
    let summary = LiftSummary::new(
        operation,
        input.first_index,
        parent,
        identity,
        root_transfer,
        hypotheses,
        parent_proof,
    );
    let interlacing = basis_map_interlacing(&summary, map, report, coefficient_roots);
    summary.with_input_interlacing(interlacing)
}

/// The interlacing rule of the basis-map lifts: the fixed-preserver
/// (Hermite--Kakeya--Obreschkoff) rule, which applies because `T_P` does
/// not depend on `n`, when the map theorem proves that `T_P` maps every
/// real-rooted polynomial (roots anywhere on the line) to a real-rooted
/// polynomial or zero.  The coefficient rows' interlacing comes from the
/// parent provers on the translates `R_n(v + U')`, which have the zeros of
/// `R_n` moved by `-U'` (same kind and orientation).  A map proved only on a
/// smaller source interval gets no claim: a combination
/// `alpha R_(n-1) + beta R_n` may have zeros outside it.
fn basis_map_interlacing(
    summary: &LiftSummary,
    map: &BasisMap,
    report: &BasisMapTheoremReport,
    coefficient_roots: Option<&BasisMapCoefficientRootProof>,
) -> LiftInterlacingClaim {
    let whole_line = report.is_proved()
        && report.conclusions.iter().any(|c| {
            c.input_interval == AlgebraicInterval::whole_line()
                && c.hypotheses.iter().all(|name| {
                    report
                        .hypotheses
                        .iter()
                        .any(|h| &h.name == name && h.status == LiftHypothesisStatus::Verified)
                })
        });
    if !whole_line {
        return LiftInterlacingClaim::not_proved(
            "the map theorem does not prove that T_P preserves real-rootedness on the whole \
             line",
        );
    }
    let Some(roots) = coefficient_roots else {
        return LiftInterlacingClaim::not_proved("the coefficient rows are not located");
    };
    let translated = proved_parent_interlacing(&roots.upper_bound).map(|shifted| {
        let mut interlacing: LiftInterlacing = shifted.clone();
        interlacing.justification.premises.push(format!(
            "the rows R_n(v + U'), U' = {}, have the zeros of R_n moved by -U', so R_(m-1), R_m \
             interlace with the same kind and orientation",
            roots
                .target_interval
                .upper
                .as_ref()
                .map_or_else(|| "?".to_string(), |end| end.value.clone())
        ));
        interlacing
    });
    fixed_preserver_interlacing(
        summary,
        translated.as_ref().map_err(Clone::clone),
        format!(
            "T_P (T_P(u^k) = P_({} + k), P = {}) does not depend on n and maps every real-rooted \
             polynomial to a real-rooted polynomial or zero (map theorem: {})",
            map.first_index(),
            map.name(),
            report.theorem
        ),
    )
}
