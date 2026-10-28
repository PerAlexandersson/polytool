//! Detection of row factorizations `Q_n = F^(e(n)) S_n`, the cofactor rule,
//! the report types (`polytool.recurrence-row-factor-lifts.v1`) and
//! [`find_row_factor_lifts`].

use super::super::hadamard_factor::prepare::{prepare_rule, PreparedRule};
use super::super::lift_summary::{nonpositive_ray, parent_proof_summary, real_line};
use super::super::veronese_algebra::classify_row;
use super::certify::{divide_operator, substitute, RowFactorCertificate};
use crate::factorization::irreducible_factorization_rational_coeffs_with_budget;
use crate::recurrence::lift::parent_proof::{prove_parent, ParentRecord};
use crate::recurrence::{
    encode_poly, ImplicitRelation, LiftHypothesis, LiftHypothesisStatus, LiftIdentityStatus,
    LiftOperation, LiftParent, LiftRootTransfer, LiftSummary, Recurrence, RecurrenceError,
    RecurrenceJson, VeroneseParentAttestation, VeroneseParentProof,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::Zero;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// JSON schema tag of a serialized [`RowFactorLiftReport`].
pub const ROW_FACTOR_LIFTS_SCHEMA: &str = "polytool.recurrence-row-factor-lifts.v1";

/// Largest number of rows factored or compared.
pub const MAX_ROW_FACTOR_ROWS: usize = 64;

/// Bounds of [`find_row_factor_lifts`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RowFactorLiftOptions {
    /// Rows factored exactly to detect common factors, between 2 and
    /// [`MAX_ROW_FACTOR_ROWS`].  Default 8.
    pub detection_rows: usize,
    /// Rows `Q_(n0), ...` checked for divisibility by `F^(e(n))` (at least the
    /// identifying rows), at most [`MAX_ROW_FACTOR_ROWS`].  Default 16.
    pub checked_rows: usize,
    /// Recombination budget of each row factorization.  Default 4096.
    pub max_recombination_subsets: usize,
    /// Caller-supplied attestations for the cofactor (the attestation route
    /// of the parent provers; a lift is then conditional, never proved).
    pub attestations: Vec<VeroneseParentAttestation>,
    /// Try the root-domain parent route.  Default true.
    pub try_root_domain: bool,
}

impl Default for RowFactorLiftOptions {
    fn default() -> Self {
        Self {
            detection_rows: 8,
            checked_rows: 16,
            max_recombination_subsets: 4096,
            attestations: Vec::new(),
            try_root_domain: true,
        }
    }
}

/// Status of the identity `Q_n = F^(e(n)) S_n` for every `n >= n0`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
pub enum RowFactorIdentityStatus {
    /// Proved for every `n >= n0` (`row_factor` module documentation).
    Certified {
        /// The replayable certificate.
        certificate: Box<RowFactorCertificate>,
    },
    /// The compared rows agree, but the certificate is incomplete: finite
    /// evidence only.
    SupportedByRows {
        /// Number of rows compared from `n0`.
        rows: usize,
        /// Why.
        reason: String,
    },
}

/// One row-factor lift `Q_n = F(x)^(e(n)) S_n`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowFactorLift {
    /// `F`, ascending exact coefficients (primitive, positive leading
    /// coefficient).
    pub factor: Vec<String>,
    /// The slope `alpha >= 0` of `e(n) = alpha n + beta`.
    pub exponent_slope: usize,
    /// The offset `beta`.
    pub exponent_offset: i64,
    /// Whether `F` is real-rooted (exact).
    pub factor_real_rooted: bool,
    /// Whether `F` is real-rooted with every root in `(-inf, 0]`.
    pub factor_roots_nonpositive: bool,
    /// The cofactor `S` (`polytool.recurrence.v1`, first index `n0`).
    pub cofactor: RecurrenceJson,
    /// Number of rows compared exactly.
    pub rows_compared: usize,
    /// The identity status.
    pub identity: RowFactorIdentityStatus,
    /// The cofactor statement (every `S_n` zero or real-rooted with
    /// coefficients of one sign) and the provers tried.
    pub parent_proof: VeroneseParentProof,
    /// The uniform [`LiftSummary`] (always equal to
    /// [`RowFactorLift::recompute_lift_summary`]).
    pub lift_summary: LiftSummary,
}

/// A detected factor schedule that gave no lift.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowFactorRejection {
    /// `F`, ascending exact coefficients.
    pub factor: Vec<String>,
    /// The slope of `e(n)`.
    pub exponent_slope: usize,
    /// The offset of `e(n)`.
    pub exponent_offset: i64,
    /// The first absolute index `n` with `F^(e(n))` not dividing `Q_n`, when
    /// that is the reason.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub failing_index: Option<usize>,
    /// Why.
    pub detail: String,
}

/// `polytool.recurrence-row-factor-lifts.v1`: the result of
/// [`find_row_factor_lifts`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowFactorLiftReport {
    /// Always [`ROW_FACTOR_LIFTS_SCHEMA`].
    pub schema: String,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Number of rows factored.
    pub rows_factored: usize,
    /// The lifts found.
    pub lifts: Vec<RowFactorLift>,
    /// Detected schedules that gave no lift (near misses).
    pub rejected: Vec<RowFactorRejection>,
    /// Further remarks (budgets, zero rows).
    pub notes: Vec<String>,
}

impl RowFactorLift {
    /// The [`LiftSummary`] of this lift, computed from the record.
    pub fn recompute_lift_summary(&self) -> LiftSummary {
        row_factor_summary(RowFactorParts {
            factor: &self.factor,
            slope: self.exponent_slope,
            offset: self.exponent_offset,
            real_rooted: self.factor_real_rooted,
            nonpositive: self.factor_roots_nonpositive,
            cofactor: &self.cofactor,
            identity: &self.identity,
            parent_proof: &self.parent_proof,
        })
    }

    /// Replay the certificate (when certified), check that it belongs to
    /// this lift, and that the recorded summary equals the recomputed one.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if let RowFactorIdentityStatus::Certified { certificate } = &self.identity {
            certificate.verify()?;
            if certificate.factor != self.factor
                || certificate.exponent_slope != self.exponent_slope
                || certificate.exponent_offset != self.exponent_offset
                || certificate.cofactor != self.cofactor
            {
                return Err(RecurrenceError::VerificationFailed(
                    "the certificate does not belong to the lift".into(),
                ));
            }
        }
        let f = crate::recurrence::decode_poly(&self.factor)?;
        let shape = classify_row(&f);
        if shape.is_real_rooted != self.factor_real_rooted
            || shape.has_nonpositive_real_roots() != self.factor_roots_nonpositive
            || self.lift_summary != self.recompute_lift_summary()
            || !self.lift_summary.is_consistent()
        {
            return Err(RecurrenceError::VerificationFailed(
                "the factor checks or the lift summary do not match the record".into(),
            ));
        }
        Ok(())
    }
}

impl RowFactorLiftReport {
    /// Check the schema and replay every lift ([`RowFactorLift::verify`]).
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != ROW_FACTOR_LIFTS_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        self.lifts.iter().try_for_each(RowFactorLift::verify)
    }

    /// The lifts whose summary is proved.
    pub fn proved_lifts(&self) -> impl Iterator<Item = &RowFactorLift> {
        self.lifts
            .iter()
            .filter(|lift| lift.lift_summary.is_proved())
    }
}

fn factor_text(factor: &[String]) -> String {
    let mut parts = Vec::new();
    for (k, c) in factor.iter().enumerate() {
        if c == "0" {
            continue;
        }
        parts.push(match k {
            0 => c.clone(),
            1 => format!("{c} x"),
            _ => format!("{c} x^{k}"),
        });
    }
    parts.join(" + ")
}

fn validate_options(options: &RowFactorLiftOptions) -> Result<(), RecurrenceError> {
    if !(2..=MAX_ROW_FACTOR_ROWS).contains(&options.detection_rows) {
        return Err(RecurrenceError::LimitExceeded(format!(
            "detection_rows must be between 2 and {MAX_ROW_FACTOR_ROWS}"
        )));
    }
    if !(1..=MAX_ROW_FACTOR_ROWS).contains(&options.checked_rows) {
        return Err(RecurrenceError::LimitExceeded(format!(
            "checked_rows must be between 1 and {MAX_ROW_FACTOR_ROWS}"
        )));
    }
    Ok(())
}

/// The multiplicity of `f` in `row` (nonzero), by repeated exact division.
fn multiplicity(row: &[Q], f: &[Q]) -> usize {
    let divisor = crate::polynomial::Polynomial::new(f.to_vec());
    let mut current = crate::polynomial::Polynomial::new(row.to_vec());
    let mut count = 0;
    loop {
        let (quotient, remainder) = current.div_rem(&divisor);
        if !remainder.is_zero() || quotient.is_zero() {
            return count;
        }
        current = quotient;
        count += 1;
    }
}

/// `e(n) = e(n0) + slope (n - n0)` with the largest slope `>= 0` such that
/// `e(n) <= m(n)` on every row; `None` when `e` vanishes identically.
fn fit_schedule(multiplicities: &[usize], n0: usize) -> Option<(usize, i64)> {
    let first = *multiplicities.first()?;
    let slope = multiplicities
        .iter()
        .enumerate()
        .skip(1)
        .map(|(offset, m)| m.saturating_sub(first) / offset)
        .min()
        .unwrap_or(0);
    let fits = multiplicities
        .iter()
        .enumerate()
        .all(|(offset, m)| first + slope * offset <= *m);
    if !fits || (first == 0 && slope == 0) {
        return None;
    }
    let offset = first as i64 - (slope * n0) as i64;
    Some((slope, offset))
}

/// Find row-factor lifts `Q_n = F(x)^(e(n)) S_n` of the sequence defined by
/// `recurrence` and `rows` (`rows[0]` is `Q_(first_index)`): common
/// irreducible factors of the first rows, an affine exponent schedule, and
/// the cofactor rule obtained by substitution.
///
/// The detection, the certificate and its soundness argument and the root
/// transfer are documented in the `row_factor` module (see the
/// [`crate::recurrence`] module list).  A missing lift says nothing against
/// real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError::InvalidInput`] for rows disagreeing with the rule;
/// [`RecurrenceError::Unsupported`] for a rule outside the class;
/// [`RecurrenceError::LimitExceeded`] for options beyond their bounds;
/// [`RecurrenceError::VerificationFailed`] if an internal self-check fails.
/// A factorization over its recombination budget is reported in the notes.
///
/// # Example
///
/// ```
/// use polytool::recurrence::{
///     find_row_factor_lifts, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RowFactorLiftOptions,
/// };
///
/// let q = |n: i64| BigRational::from_integer(n.into());
/// // Q_n = (1 + x)^2 Q_(n-1), Q_0 = 1 + x: Q_n = (1 + x)^(2 n + 1).
/// let rule = Recurrence {
///     terms: vec![RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(1), q(2), q(1)]] })],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let report = find_row_factor_lifts(&rule, &[vec![q(1), q(1)]], 0, &RowFactorLiftOptions::default())?;
/// let lift = report.proved_lifts().next().expect("a proved lift");
/// assert_eq!((lift.exponent_slope, lift.exponent_offset), (2, 1));
/// report.verify()?;
/// # Ok::<(), polytool::recurrence::RecurrenceError>(())
/// ```
pub fn find_row_factor_lifts(
    recurrence: &Recurrence,
    rows: &[Vec<BigRational>],
    first_index: usize,
    options: &RowFactorLiftOptions,
) -> Result<RowFactorLiftReport, RecurrenceError> {
    validate_options(options)?;
    let input = prepare_rule(recurrence, rows, first_index, "input")?;
    let n0 = first_index;
    let detection_rows = input.rows(options.detection_rows)?;
    let mut report = RowFactorLiftReport {
        schema: ROW_FACTOR_LIFTS_SCHEMA.into(),
        input_first_index: n0,
        rows_factored: detection_rows.len(),
        lifts: Vec::new(),
        rejected: Vec::new(),
        notes: Vec::new(),
    };
    if detection_rows
        .iter()
        .any(|row| row.iter().all(Zero::is_zero))
    {
        report
            .notes
            .push("a detection row is zero; no common factor is detected".into());
        return Ok(report);
    }
    // Irreducible factors of the detection rows, with their multiplicities
    // in every detection row.
    let mut irreducibles = std::collections::BTreeSet::new();
    for (offset, row) in detection_rows.iter().enumerate() {
        match irreducible_factorization_rational_coeffs_with_budget(
            row,
            options.max_recombination_subsets,
        ) {
            Ok(factorization) => {
                irreducibles.extend(factorization.factors.into_iter().map(|(f, _)| f));
            }
            Err(error) => report.notes.push(format!(
                "the factorization of row {} failed: {error}",
                n0 + offset
            )),
        }
    }
    let mut schedules: BTreeMap<(usize, i64), Vec<Q>> = BTreeMap::new();
    for coefficients in &irreducibles {
        let g: Vec<Q> = coefficients.iter().cloned().map(Q::from_integer).collect();
        let multiplicities: Vec<usize> = detection_rows
            .iter()
            .map(|row| multiplicity(row, &g))
            .collect();
        if let Some(schedule) = fit_schedule(&multiplicities, n0) {
            let entry = schedules
                .entry(schedule)
                .or_insert_with(|| vec![Q::from_integer(1.into())]);
            *entry = vec_poly::mul(entry, &g);
        }
    }
    for ((slope, offset), f) in schedules {
        match certify_schedule(&input, &f, slope, offset, options)? {
            Ok(lift) => report.lifts.push(*lift),
            Err(rejection) => report.rejected.push(rejection),
        }
    }
    Ok(report)
}

fn certify_schedule(
    input: &PreparedRule,
    f: &[Q],
    slope: usize,
    offset: i64,
    options: &RowFactorLiftOptions,
) -> Result<Result<Box<RowFactorLift>, RowFactorRejection>, RecurrenceError> {
    let n0 = input.first_index;
    let factor = encode_poly(f);
    let reject = |failing_index: Option<usize>, detail: String| {
        Ok(Err(RowFactorRejection {
            factor: factor.clone(),
            exponent_slope: slope,
            exponent_offset: offset,
            failing_index,
            detail,
        }))
    };
    let exponent = |n: usize| (slope * n) as i64 + offset;
    let r = input.leader_shift;
    let j = input.relation.max_derivative();
    let substituted = substitute(&input.relation, f, slope, offset);
    let h = slope * r + j;
    let Some(cofactor_relation) = divide_operator(&substituted, f, h) else {
        return reject(
            None,
            format!(
                "the substituted relation is not divisible by F^{h}, so the cofactor rule would \
                 have a leader depending on x"
            ),
        );
    };
    let leader = cofactor_relation.column(r);
    if leader.0.len() != 1 || !leader.0.get(&(0, 0)).is_some_and(|c| c.is_y_free()) {
        return reject(None, "the cofactor leader depends on x or D".into());
    }
    let Some((cofactor_rule, _)) =
        ImplicitRelation::from_op(&cofactor_relation).to_successor_rule()
    else {
        return reject(None, "the cofactor relation is not a successor rule".into());
    };
    // Compare F^(e(n)) S_n with Q_n on the checked rows.
    let count = options.checked_rows.max(input.seeds.len()).max(r);
    let input_rows = input.rows(count)?;
    let mut cofactor_rows = Vec::with_capacity(count);
    for (k, row) in input_rows.iter().enumerate() {
        let power = vec_poly::pow(f, exponent(n0 + k) as usize);
        let (quotient, remainder) = crate::polynomial::Polynomial::new(row.clone())
            .div_rem(&crate::polynomial::Polynomial::new(power));
        if !remainder.is_zero() {
            return reject(
                Some(n0 + k),
                format!("F^(e(n)) does not divide Q_n at n = {}", n0 + k),
            );
        }
        cofactor_rows.push(vec_poly::trim(quotient.coeffs()));
    }
    let count_s = cofactor_rule.generation_initial_count(n0, count);
    let cofactor = match prepare_rule(
        &cofactor_rule,
        &cofactor_rows[..count_s.max(r)],
        n0,
        "cofactor",
    ) {
        Ok(cofactor) => cofactor,
        Err(error) => return reject(None, format!("the cofactor rule is not usable: {error}")),
    };
    let generated = cofactor.rows(count)?;
    if let Some(offset_index) = generated
        .iter()
        .zip(&cofactor_rows)
        .position(|(a, b)| a != b)
    {
        return reject(
            Some(n0 + offset_index),
            format!(
                "the cofactor rule disagrees with Q_n / F^(e(n)) at n = {}",
                n0 + offset_index
            ),
        );
    }
    // Exceptional bases: M holds for Q from input_from, M'_S for S from
    // cofactor_from; in between M is evaluated on F^e S (= Q on these rows).
    let input_from = input.valid_from_base();
    let cofactor_from = cofactor.valid_from_base();
    let exceptional: Vec<usize> = (input_from..cofactor_from).collect();
    let mut rows_compared = count;
    let identity = if let Some(last) = exceptional.last() {
        let needed = last + r + 1 - n0;
        if needed > count {
            rows_compared = needed;
        }
        if needed > MAX_ROW_FACTOR_ROWS {
            RowFactorIdentityStatus::SupportedByRows {
                rows: count,
                reason: format!("the exceptional base {last} needs {needed} rows"),
            }
        } else {
            let candidate: Vec<Vec<Q>> = cofactor
                .rows(needed.max(count))?
                .iter()
                .enumerate()
                .map(|(k, s)| {
                    vec_poly::trim(&vec_poly::mul(
                        &vec_poly::pow(f, exponent(n0 + k) as usize),
                        s,
                    ))
                })
                .collect();
            for base in &exceptional {
                let value = input.relation.apply(&candidate, n0, *base).ok_or_else(|| {
                    RecurrenceError::VerificationFailed(format!("rows missing at base {base}"))
                })?;
                if !value.iter().all(Zero::is_zero) {
                    return Err(RecurrenceError::VerificationFailed(format!(
                        "the input relation does not vanish on the candidate rows at base {base}"
                    )));
                }
            }
            certified(
                input,
                &cofactor,
                &factor,
                (slope, offset),
                &cofactor_relation,
                exceptional,
            )?
        }
    } else {
        certified(
            input,
            &cofactor,
            &factor,
            (slope, offset),
            &cofactor_relation,
            exceptional,
        )?
    };
    let parent_proof = prove_parent(
        &ParentRecord {
            rule: &cofactor.rule,
            seeds: &cofactor.seeds,
            first_index: n0,
        },
        n0,
        &options.attestations,
        options.try_root_domain,
    )?;
    let shape = classify_row(f);
    let cofactor_json =
        RecurrenceJson::from_recurrence_rational(&cofactor.rule, n0, &cofactor.seeds, None);
    let lift_summary = row_factor_summary(RowFactorParts {
        factor: &factor,
        slope,
        offset,
        real_rooted: shape.is_real_rooted,
        nonpositive: shape.has_nonpositive_real_roots(),
        cofactor: &cofactor_json,
        identity: &identity,
        parent_proof: &parent_proof,
    });
    Ok(Ok(Box::new(RowFactorLift {
        factor,
        exponent_slope: slope,
        exponent_offset: offset,
        factor_real_rooted: shape.is_real_rooted,
        factor_roots_nonpositive: shape.has_nonpositive_real_roots(),
        cofactor: cofactor_json,
        rows_compared,
        identity,
        parent_proof,
        lift_summary,
    })))
}

fn certified(
    input: &PreparedRule,
    cofactor: &PreparedRule,
    factor: &[String],
    (slope, offset): (usize, i64),
    cofactor_relation: &crate::recurrence::implicit::algebra::DOp,
    exceptional: Vec<usize>,
) -> Result<RowFactorIdentityStatus, RecurrenceError> {
    let n0 = input.first_index;
    let r = input.leader_shift;
    let j = input.relation.max_derivative();
    let identifying = input.seeds.len();
    let certificate = RowFactorCertificate {
        factor: factor.to_vec(),
        exponent_slope: slope,
        exponent_offset: offset,
        input_relation: ImplicitRelation::from_op(&input.relation).to_json(),
        input_leader_shift: r,
        derivative_order: j,
        removed_factor_power: slope * r + j,
        cofactor_relation: ImplicitRelation::from_op(cofactor_relation).to_json(),
        cofactor: RecurrenceJson::from_recurrence_rational(
            &cofactor.rule,
            n0,
            &cofactor.seeds,
            None,
        ),
        input_valid_from_base: input.valid_from_base(),
        cofactor_valid_from_base: cofactor.valid_from_base(),
        exceptional_bases: exceptional,
        identifying_indices: (n0..n0 + identifying).collect(),
        statement: format!(
            "M_S = F^h M'_S with an x-free leader, so the cofactor rows S_n are polynomials; \
             F^(e(n)) S_n satisfies the input rule at every base n >= {} and agrees with Q_n on \
             the {identifying} identifying rows, so Q_n = F^(e(n)) S_n for every n >= {n0}",
            input.valid_from_base()
        ),
    };
    certificate.verify()?;
    Ok(RowFactorIdentityStatus::Certified {
        certificate: Box::new(certificate),
    })
}

/// The parts of a [`RowFactorLift`] its summary depends on.
struct RowFactorParts<'a> {
    factor: &'a [String],
    slope: usize,
    offset: i64,
    real_rooted: bool,
    nonpositive: bool,
    cofactor: &'a RecurrenceJson,
    identity: &'a RowFactorIdentityStatus,
    parent_proof: &'a VeroneseParentProof,
}

fn row_factor_summary(parts: RowFactorParts<'_>) -> LiftSummary {
    let n0 = parts.cofactor.first_index;
    let f = factor_text(parts.factor);
    let (alpha, beta) = (parts.slope, parts.offset);
    let identity = match parts.identity {
        RowFactorIdentityStatus::Certified { .. } => LiftIdentityStatus::Certified {
            method: "substitution_division".into(),
        },
        RowFactorIdentityStatus::SupportedByRows { rows, .. } => {
            LiftIdentityStatus::SupportedByRows { rows: *rows }
        }
    };
    let hypothesis = LiftHypothesis {
        name: "factor_real_rooted".into(),
        statement: format!("F = {f} is real-rooted"),
        status: if parts.real_rooted {
            LiftHypothesisStatus::Verified
        } else {
            LiftHypothesisStatus::Failed { witness: None }
        },
        evidence: if parts.real_rooted {
            "exact real-root count (is_real_rooted_bigint_coeffs)".into()
        } else {
            "F has a non-real root; every Q_n with e(n) >= 1 has one too".into()
        },
        sign_decision: None,
    };
    let (input_region, conclusion) = if parts.nonpositive {
        (
            nonpositive_ray(),
            format!("every Q_n, n >= {n0}, is zero or real-rooted with roots in (-inf, 0]"),
        )
    } else {
        (
            real_line(),
            format!("every Q_n, n >= {n0}, is zero or real-rooted (roots of F and of S_n)"),
        )
    };
    LiftSummary::new(
        LiftOperation::RowFactor {
            factor: parts.factor.to_vec(),
            exponent_slope: alpha,
            exponent_offset: beta,
        },
        n0,
        LiftParent {
            recurrence: parts.cofactor.clone(),
            relation: format!("Q_n = ({f})^({alpha} n {beta:+}) S_n for every n >= {n0}"),
            gauge: None,
            index_shift: 0,
        },
        identity,
        LiftRootTransfer {
            theorem: "a product of polynomials is real-rooted iff every factor is; its roots \
                          are the union of theirs"
                .into(),
            parent_hypothesis: "every cofactor row S_n is zero or real-rooted with all roots \
                                    in (-inf, 0] (coefficients of one sign)"
                .into(),
            parent_root_region: Some(nonpositive_ray()),
            input_root_region: input_region,
            input_conclusion: conclusion,
        },
        vec![hypothesis],
        parent_proof_summary(parts.parent_proof),
    )
}
