//! Provers for the parent statement of a Veronese-section lift: every parent
//! row `P_n` (from a given index on) is zero or real-rooted with all
//! coefficients of one sign, equivalently zero or with every root real and
//! `<= 0`.  By the transfer theorem (see [`find_veronese_section_lifts`])
//! every section of such a row is zero or real-rooted with nonpositive roots.
//!
//! [`prove_parent`] tries, in order:
//!
//! 1. **periodic product**: the rule is `P_(n+k) = c(x) P_n` with `c` free of
//!    `n` (a nonzero constant denominator is allowed), and `c` and the `k`
//!    seeds are nonzero, real-rooted and of one sign; then
//!    `P_(n0 + k m + i) = c^m P_(n0 + i)` (shared with the v1 search);
//! 2. **paper theorem**: [`match_recurrence_theorems`] with boundary status
//!    `canonical` or `canonical_after_relabel` and, for the predecessor
//!    theorem, an all-index certificate scope; both theorems conclude that
//!    every row has nonnegative coefficients and is real-rooted;
//! 3. **root domain**: [`transform_recurrence`] with a root-domain request
//!    over the maps identity (`affine 1, 0`), `f_polynomial`, `reversal` (the
//!    last two with the degree schedule `d_n` fitted to the first rows, which
//!    the transform certifies or reports as a condition) and the inverse
//!    `f`-map `mobius (1, 0, -1, 1)`.  Only a claim about the source rows
//!    with property `real_rooted_in` an index-independent interval whose
//!    upper end is `<= 0`, evidence `completed_recurrence_induction` or
//!    `certified_family_invariant`, no `conditional_on` entries, no last
//!    index and a first index `<=` the parent's first index is accepted;
//! 4. **attestation**: a caller-supplied [`VeroneseParentAttestation`] whose
//!    record equals the parent exactly (same first index, the same relation
//!    term by term, the same identifying seeds); the verdict is then
//!    conditional on the attestation, never proved.
//!
//! A few parent rows are also checked exactly (nonzero rows of one sign and
//! real-rooted).  A failing row refutes the parent statement; no route can
//! then succeed and an attestation is not applied.

use super::veronese_algebra::classify_row;
use super::*;
use crate::renewal::{match_recurrence_theorems, CertificateScope, TheoremBoundaryStatus};
use crate::vec_poly;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Number of parent rows checked exactly for being zero or real-rooted with
/// coefficients of one sign (evidence only; a failure refutes the parent).
const PARENT_ROWS_CHECKED: usize = 10;

/// Rows generated for the root-domain transform of a parent.
const ROOT_DOMAIN_ROWS: usize = 8;

/// A caller-supplied statement that every parent row `P_n`, `n >= first_index`,
/// is zero or real-rooted with coefficients of one sign (for example a
/// published theorem).  It is applied only to a parent whose record equals
/// `recurrence` exactly; the verdict is then conditional on it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentAttestation {
    /// Short label, repeated in the verdict.
    pub label: String,
    /// The attested parent: rule, first index and seeds.
    pub recurrence: RecurrenceJson,
    /// First parent index covered by the statement; it must be at most the
    /// first parent index the lift uses.
    pub first_index: usize,
    /// The attested statement, in words.
    pub statement: String,
    /// Where the statement comes from (reference, theorem, URL).
    pub source: String,
}

/// Outcome of the parent provers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseParentProofStatus {
    /// Proved by a periodic product, a paper theorem or the root-domain route.
    Proved,
    /// Not proved; a matching caller-supplied attestation asserts it.
    ConditionalOnAttestation,
    /// No route succeeded.
    NotProved,
    /// A checked parent row is not zero and not real-rooted with one-signed
    /// coefficients: the parent statement is false.
    Refuted,
}

/// The route that settled the parent statement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum VeroneseParentProofRoute {
    /// `P_(n+k) = c(x) P_n` with `c` and the seeds real-rooted of one sign.
    PeriodicProduct,
    /// A matched paper theorem (`match_recurrence_theorems`).
    PaperTheorem,
    /// A root-domain claim of `transform_recurrence`.
    RootDomain,
    /// A caller-supplied attestation.
    Attestation,
}

/// One tried route.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentProofAttempt {
    /// The route.
    pub route: VeroneseParentProofRoute,
    /// Theorem identifier, transformation kind or attestation label.
    pub variant: String,
    /// Whether it settled the parent statement.
    pub success: bool,
    /// What was checked, or why it did not apply.
    pub detail: String,
}

/// The parent statement and how it was settled.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VeroneseParentProof {
    /// The outcome.
    pub status: VeroneseParentProofStatus,
    /// The successful route, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub route: Option<VeroneseParentProofRoute>,
    /// The parent statement as proved, attested or refuted, or why it is open.
    pub statement: String,
    /// Attestation labels the status depends on (empty unless conditional).
    pub conditional_on: Vec<String>,
    /// Every route tried, in order.
    pub attempts: Vec<VeroneseParentProofAttempt>,
    /// Absolute index of the first checked parent row that is nonzero and not
    /// real-rooted with coefficients of one sign, if any.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub failing_row_index: Option<usize>,
    /// Number of parent rows checked exactly (evidence only).
    pub rows_checked: usize,
}

/// A parent rule with its identifying seeds and first index.
pub(super) struct ParentRecord<'a> {
    pub(super) rule: &'a Recurrence,
    /// The `K` identifying seeds (`generation_initial_count`).
    pub(super) seeds: &'a [Vec<Q>],
    pub(super) first_index: usize,
}

/// Route 1: `P_(n+k) = c(x) P_n` with `c` and the `k` seeds nonzero,
/// real-rooted and of one sign.  Returns the reason when it does not apply.
pub(super) fn periodic_product_parent(rule: &Recurrence, seeds: &[Vec<Q>]) -> Result<(), String> {
    let (relation, k) = ImplicitRelation::from_recurrence(rule).map_err(|e| e.to_string())?;
    let op = relation.op();
    let leader =
        op.0.get(&(k, 0))
            .and_then(|c| c.as_constant())
            .filter(|c| !c.is_zero())
            .ok_or("the leader is not a nonzero constant")?;
    let mut factor = None;
    for ((shift, order), coefficient) in &op.0 {
        if *shift == k {
            continue;
        }
        if *shift != 0 || *order != 0 || !coefficient.is_n_free() {
            return Err("the rule is not P_(n+k) = c(x) P_n with c free of n".into());
        }
        factor = coefficient.rows().first().cloned();
    }
    let c: Vec<Q> = factor
        .unwrap_or_default()
        .iter()
        .map(|value| -value / &leader)
        .collect();
    if !classify_row(&vec_poly::trim(&c)).has_nonpositive_real_roots() {
        return Err("the factor c(x) is not nonzero, real-rooted and of one sign".into());
    }
    if seeds.len() != k {
        return Err(format!(
            "the rule needs {k} seeds, {} identify it",
            seeds.len()
        ));
    }
    if let Some(position) = seeds
        .iter()
        .position(|seed| !classify_row(seed).has_nonpositive_real_roots())
    {
        return Err(format!(
            "seed {position} is not nonzero, real-rooted and of one sign"
        ));
    }
    Ok(())
}

fn theorem_route(record: &ParentRecord<'_>, rows: &[Vec<Q>]) -> Vec<VeroneseParentProofAttempt> {
    let matches = match_recurrence_theorems(record.rule, record.first_index, rows);
    let mut attempts = Vec::new();
    for theorem_match in matches.theorem_matches {
        let boundary_ok = matches!(
            theorem_match.boundary.status,
            TheoremBoundaryStatus::Canonical | TheoremBoundaryStatus::CanonicalAfterRelabel
        );
        let scope_ok = theorem_match
            .variable_predecessor
            .as_ref()
            .is_none_or(|certificate| certificate.scope == CertificateScope::AllIndices);
        let success = boundary_ok && scope_ok;
        let detail = if success {
            format!(
                "{} ({}): {}; the identifying seeds equal the canonical rows ({:?}), so every \
                 parent row is a canonical row",
                theorem_match.family,
                theorem_match.proof_status,
                theorem_match.conclusion,
                theorem_match.boundary.status
            )
        } else if !boundary_ok {
            format!(
                "the parent rows are not the canonical boundary rows ({:?})",
                theorem_match.boundary.status
            )
        } else {
            "the theorem certificate covers finitely many indices only".into()
        };
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: theorem_match.theorem.clone(),
            success,
            detail,
        });
    }
    if attempts.is_empty() {
        let reasons: Vec<String> = matches
            .rejected_theorems
            .iter()
            .map(|r| format!("{}: {}", r.theorem, r.reason))
            .collect();
        attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::PaperTheorem,
            variant: "match_recurrence_theorems".into(),
            success: false,
            detail: format!("no theorem matched ({})", reasons.join("; ")),
        });
    }
    attempts
}

/// `d_n = d_(n0) + step (n - n0)` fitted to the degrees of the first rows,
/// when every observed row is nonzero and the degrees are affine.
fn fitted_schedule(first_index: usize, rows: &[Vec<Q>]) -> Option<DegreeSchedule> {
    let degrees: Vec<i64> = rows
        .iter()
        .map(|row| vec_poly::degree(row).and_then(|d| i64::try_from(d).ok()))
        .collect::<Option<_>>()?;
    let first = *degrees.first()?;
    let step = degrees.get(1).map_or(0, |second| second - first);
    let affine = degrees
        .iter()
        .enumerate()
        .all(|(i, d)| i64::try_from(i).is_ok_and(|i| *d == first + step * i));
    affine.then(|| DegreeSchedule::new(first_index, first as usize, step))
}

fn accepted_source_claim(claim: &RootDomainClaim, first_index: usize) -> bool {
    let evidence_ok = matches!(
        claim.evidence,
        RootEvidence::CompletedRecurrenceInduction | RootEvidence::CertifiedFamilyInvariant
    );
    let upper_ok = claim.interval.as_ref().is_some_and(|interval| {
        interval.is_index_independent()
            && interval.upper.as_ref().is_some_and(|upper| {
                parse_rational_coeff(&upper.value).is_ok_and(|value| value <= Q::zero())
            })
    });
    claim.subject == RootSubject::Source
        && claim.component == 0
        && claim.property == RootProperty::RealRootedIn
        && upper_ok
        && evidence_ok
        && claim.conditional_on.is_empty()
        && claim.last_index.is_none()
        && claim.first_index <= first_index
}

fn root_domain_route(
    record: &ParentRecord<'_>,
    rows: &[Vec<Q>],
) -> Vec<VeroneseParentProofAttempt> {
    let source = RecurrenceJson::from_recurrence_rational(
        record.rule,
        record.first_index,
        record.seeds,
        None,
    );
    let schedule = fitted_schedule(record.first_index, rows);
    let mut menu = vec![(
        "identity (affine 1, 0)".to_string(),
        Some(RecurrenceTransformation::Affine {
            a: "1".into(),
            b: "0".into(),
        }),
    )];
    for (name, make) in [
        (
            "f_polynomial",
            (|s: DegreeSchedule| RecurrenceTransformation::FPolynomial { schedule: s })
                as fn(DegreeSchedule) -> RecurrenceTransformation,
        ),
        ("reversal", |s| RecurrenceTransformation::Reversal {
            schedule: s,
        }),
        ("mobius (1, 0, -1, 1)", |s| {
            RecurrenceTransformation::Mobius {
                a: "1".into(),
                b: "0".into(),
                c: "-1".into(),
                e: "1".into(),
                schedule: s,
            }
        }),
    ] {
        menu.push((name.to_string(), schedule.clone().map(make)));
    }
    let mut attempts = Vec::new();
    for (name, transformation) in menu {
        let Some(transformation) = transformation else {
            attempts.push(VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: name,
                success: false,
                detail: "no affine degree schedule fits the first parent rows".into(),
            });
            continue;
        };
        let mut request = RecurrenceTransformRequest::from_source(
            transformation,
            source.clone(),
            ROOT_DOMAIN_ROWS,
        );
        request.root_domain = Some(RootDomainRequest::default());
        let attempt = match transform_recurrence(&request) {
            Err(error) => VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: name,
                success: false,
                detail: format!("the transform failed: {error}"),
            },
            Ok(report) => {
                let claim = report.root_domain.as_ref().and_then(|root| {
                    root.claims
                        .iter()
                        .find(|claim| accepted_source_claim(claim, record.first_index))
                });
                match claim {
                    Some(claim) => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: name,
                        success: true,
                        detail: format!(
                            "source claim: every P_n, n >= {}, is real-rooted in {} ({:?}): {}",
                            claim.first_index,
                            claim
                                .interval
                                .as_ref()
                                .map_or_else(String::new, ToString::to_string),
                            claim.evidence,
                            claim.derivation
                        ),
                    },
                    None => VeroneseParentProofAttempt {
                        route: VeroneseParentProofRoute::RootDomain,
                        variant: name,
                        success: false,
                        detail: "no unconditional all-index source claim with roots in \
                                 (-inf, 0]"
                            .into(),
                    },
                }
            }
        };
        let success = attempt.success;
        attempts.push(attempt);
        if success {
            break;
        }
    }
    attempts
}

/// Whether two records define the same parent: same first index, the same
/// relation term by term, the same number of identifying seeds and equal
/// seeds (the supplied rows of `other` must agree with its rule).
fn same_record(record: &ParentRecord<'_>, other: &RecurrenceJson) -> Result<bool, RecurrenceError> {
    let (rule, first_index, rows) = other.to_recurrence_parts()?;
    if first_index != record.first_index || rows.len() < rule.max_offset() {
        return Ok(false);
    }
    let (left, left_shift) = ImplicitRelation::from_recurrence(record.rule)?;
    let (right, right_shift) = ImplicitRelation::from_recurrence(&rule)?;
    if left_shift != right_shift || left.op() != right.op() {
        return Ok(false);
    }
    let count = rule.generation_initial_count(first_index, rows.len());
    if count != record.seeds.len() {
        return Ok(false);
    }
    let generated = match rule.generate_rows_rational(&rows[..count], first_index, rows.len()) {
        Ok(generated) => generated,
        Err(_) => return Ok(false),
    };
    let equal = |a: &[Q], b: &[Q]| vec_poly::trim(a) == vec_poly::trim(b);
    Ok(rows.iter().zip(&generated).all(|(a, b)| equal(a, b))
        && rows[..count]
            .iter()
            .zip(record.seeds)
            .all(|(a, b)| equal(a, b)))
}

/// First checked row that is nonzero and not real-rooted of one sign.
fn first_failing_row(rows: &[Vec<Q>], first_index: usize) -> Option<usize> {
    rows.iter()
        .position(|row| {
            let shape = classify_row(row);
            !shape.is_zero && !(shape.has_one_sign && shape.is_real_rooted)
        })
        .map(|offset| first_index + offset)
}

/// Run the parent provers (module documentation).  `needed_from` is the
/// first parent index the lift uses (at least `first_index`);
/// `try_root_domain` enables route 3.
pub(super) fn prove_parent(
    record: &ParentRecord<'_>,
    needed_from: usize,
    attestations: &[VeroneseParentAttestation],
    try_root_domain: bool,
) -> Result<VeroneseParentProof, RecurrenceError> {
    let rows = record
        .rule
        .generate_rows_rational(
            record.seeds,
            record.first_index,
            PARENT_ROWS_CHECKED.max(record.seeds.len()),
        )
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!(
                "the parent recurrence cannot be replayed: {error}"
            ))
        })?;
    let mut proof = VeroneseParentProof {
        status: VeroneseParentProofStatus::NotProved,
        route: None,
        statement: String::new(),
        conditional_on: Vec::new(),
        attempts: Vec::new(),
        failing_row_index: first_failing_row(&rows, record.first_index),
        rows_checked: rows.len(),
    };
    let conclusion = format!(
        "every parent row P_n, n >= {}, is zero or real-rooted with coefficients of one sign",
        record.first_index
    );
    if let Some(index) = proof.failing_row_index {
        proof.status = VeroneseParentProofStatus::Refuted;
        proof.statement = format!(
            "parent row P_{index} is nonzero and not real-rooted with coefficients of one sign; \
             this parent gives no proof"
        );
        return Ok(proof);
    }

    // Route 1.
    let periodic = periodic_product_parent(record.rule, record.seeds);
    proof.attempts.push(VeroneseParentProofAttempt {
        route: VeroneseParentProofRoute::PeriodicProduct,
        variant: "P_(n+k) = c(x) P_n".into(),
        success: periodic.is_ok(),
        detail: match &periodic {
            Ok(()) => "c and the seeds are nonzero, real-rooted and of one sign, so \
                       P_(n0 + k m + i) = c^m P_(n0 + i) is too"
                .into(),
            Err(reason) => reason.clone(),
        },
    });
    let mut proved = periodic
        .is_ok()
        .then_some(VeroneseParentProofRoute::PeriodicProduct);

    // Route 2.
    if proved.is_none() {
        let attempts = theorem_route(record, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::PaperTheorem);
        }
        proof.attempts.extend(attempts);
    }

    // Route 3.
    if proved.is_none() && try_root_domain {
        let attempts = root_domain_route(record, &rows);
        if attempts.iter().any(|a| a.success) {
            proved = Some(VeroneseParentProofRoute::RootDomain);
        }
        proof.attempts.extend(attempts);
    }

    if let Some(route) = proved {
        proof.status = VeroneseParentProofStatus::Proved;
        proof.route = Some(route);
        proof.statement = conclusion;
        return Ok(proof);
    }

    // Route 4.
    for attestation in attestations {
        let matches = same_record(record, &attestation.recurrence)?;
        let covers = attestation.first_index <= needed_from;
        let success = matches && covers;
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: if success {
                format!("{} (source: {})", attestation.statement, attestation.source)
            } else if !matches {
                "the attested record differs from the parent".into()
            } else {
                format!(
                    "the attestation starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                )
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!("{conclusion} (attested, not proved here)");
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    proof.statement = format!(
        "open: no route proved that {conclusion}; the first {} parent rows pass the exact \
         check (evidence only)",
        proof.rows_checked
    );
    Ok(proof)
}
