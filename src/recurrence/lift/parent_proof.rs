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
//!    index and a first index `<=` the parent's first index is accepted.
//!    When no map gives such a claim, the menu is repeated with the parent's
//!    first index moved `1, 2, 3` rows later (the moved seeds are generated
//!    parent rows); a claim from the moved first index on suffices, because
//!    the skipped rows are among the rows checked exactly below (each zero or
//!    real-rooted with coefficients of one sign, which is the parent
//!    statement for those rows);
//! 4. **three-term recurrence**: [`check_three_term_recurrence_theorem`] with
//!    the root region `(-inf, 0)` (orthogonal-polynomial recurrences, also
//!    after the coordinate change `y = 1/(x + v_0)`, as for the Narayana
//!    polynomials), proved only when every hypothesis is verified;
//! 5. **Liu--Wang**: [`check_liu_wang_theorem`] (first-order rules
//!    `P_n = a_n(x) P_(n-1) + b_n(x) P'_(n-1)` with every root kept in
//!    `(-inf, 0)`, as for the type-B Eulerian polynomials), proved only when
//!    every hypothesis is verified;
//! 6. **Liu--Wang derivative lag**: [`check_liu_wang_derivative_lag_theorem`]
//!    (second-order rules `P_n = A_n(x) P_(n-1) + V_n(x) P'_(n-1) +
//!    W_n(x) P_(n-2)`, the RealRooted routes `halfGrowth` and `sign`, as for
//!    the Delannoy and deco rows; also after dividing every row by a power
//!    of `x`, as for the derangement excedance rows), proved only when every
//!    hypothesis is verified;
//! 7. **Ma--Wang interval**: [`check_ma_wang_interval_theorem`] (first-order
//!    rules `P_n = (alpha_n x + beta_n) P_(n-1) + kappa_n (x - L)(x - U)
//!    P'_(n-1)` with every root kept in `[L, U]` by the exact two-endpoint
//!    criterion, as for the rows `k! S(n, k)` of the ordered set partitions),
//!    proved when every hypothesis is verified and the certified interval
//!    lies in `(-inf, 0]` (`U <= 0`);
//! 8. **attestation**: a caller-supplied [`VeroneseParentAttestation`] whose
//!    record equals the parent exactly (same first index, the same relation
//!    term by term, the same identifying seeds); the verdict is then
//!    conditional on the attestation, never proved;
//! 9. **catalog flag**: a bundled OEIS catalog entry flagged
//!    [`OeisRealRootednessStatus::LeanProved`] (every row zero or
//!    real-rooted by a Lean proof made outside polytool) whose rule, first
//!    index and initial rows equal the parent exactly, matched as in route 8
//!    with the attestation label `oeis_catalog:<A-number>:lean_proved`.  The
//!    flag gives real-rootedness only, so the route also needs
//!    coefficients of one sign for every row, which polytool proves here
//!    (see "One-signed rows" below); the verdict is conditional on the flag,
//!    never proved.
//!
//! Routes 4 to 7 run whether or not the root-domain route is enabled.
//! Routes 4 and 5 conclude that every parent row from the first index on is
//! real-rooted with every root in `(-inf, 0)`; route 6 that every parent row
//! is zero (only seed rows) or real-rooted with every root in `(-inf, 0]`
//! (`(-inf, 0)` without the power of `x`), so it runs without a requested
//! root region; route 7 that every parent row is real-rooted with every root
//! in `[L, U]`, `U <= 0`.  All four move their induction start up to three rows past
//! the parent seeds when the theorem fails at the seeds (the default
//! [`RealRootednessTheoremOptions::maximum_start_delay`]; the skipped rows
//! are checked exactly by the theorem), and the attempt names the delay
//! used.
//!
//! A few parent rows are also checked exactly (nonzero rows of one sign and
//! real-rooted).  A failing row refutes the parent statement; no route can
//! then succeed and an attestation is not applied.
//!
//! # One-signed rows (route 9)
//!
//! Write the rule as `f(n) P_n = sum_terms c(n,x) D^d P_(n-r) + g(n,x)` for
//! `n >= N1` (the first computed index, after the `K` seeds).  Every row has
//! nonnegative coefficients (after one global sign change of a homogeneous
//! rule) when:
//!
//! - every seed has nonnegative coefficients;
//! - `f` does not depend on `x` and `f(n) > 0` for every integer `n >= N1`,
//!   no term carries the factor `(-1)^n`, and every coefficient of `g` (a
//!   polynomial in `n`) is `>= 0` for `n >= N1` (exact sign checker);
//! - for every offset `r` and shift `e`, the combined weight
//!   `W_(r,e)(n,s) = sum over the terms at offset r and the monomials c_j(n) x^j
//!   with j - d = e of c_j(n) s (s-1) ... (s-d+1)`, which multiplies the
//!   coefficient of `x^s` in `P_(n-r)` in the coefficient of `x^(s+e)` in
//!   `f(n) P_n`, is `>= 0` for every integer `n >= N1` and every
//!   `0 <= s <= deg P_(n-r)`.
//!
//! Then `f(n) P_n` is a sum of nonnegative multiples of nonnegative
//! coefficients, and induction on `n` gives the claim.  The weight condition
//! is shown by one of two sufficient tests: the polynomial
//! `W(N1 + v, s)` or `W(N1 + v, D(N1 + v - r) - u)` has nonnegative
//! coefficients in the basis `v(v-1)...(v-i+1) w(w-1)...(w-j+1)` (`w = s` or
//! `u`), which is nonnegative at all integers `v, w >= 0`.  Here
//! `D(m) = alpha + beta m` is an upper bound for `deg P_m`, `m >= ` the first
//! index: `beta >= max(0, (deg_x c - d) / r)` over the terms and `alpha`
//! covers the seeds and `g`, so that `deg P_n <= max(deg_x c - d + D(n - r))
//! <= D(n)` by induction.  Every test is exact; a failure only means the
//! route does not apply.

use super::veronese_algebra::classify_row;
use super::*;
use crate::integer_range_sign::{decide_index_polynomial_sign, IndexSignCondition, IntegerRange};
use crate::oeis::{entries_with_real_rootedness_flag, OeisRealRootednessStatus};
use crate::real_rootedness::RealInterval;
use crate::renewal::{match_recurrence_theorems, CertificateScope, TheoremBoundaryStatus};
use crate::vec_poly;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Number of parent rows checked exactly for being zero or real-rooted with
/// coefficients of one sign (evidence only; a failure refutes the parent).
const PARENT_ROWS_CHECKED: usize = 10;

/// Rows generated for the root-domain transform of a parent.
const ROOT_DOMAIN_ROWS: usize = 8;

/// Largest number of rows by which the root-domain route moves the parent's
/// first index (smaller than [`PARENT_ROWS_CHECKED`], so that every skipped
/// row is among the rows checked exactly).
const MAXIMUM_ROOT_DOMAIN_START_DELAY: usize = 3;

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
    /// Proved by a periodic product, a paper theorem, the root-domain route,
    /// the three-term recurrence theorem, the Liu--Wang theorem, the
    /// Liu--Wang derivative-lag theorem or the Ma--Wang interval theorem.
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
    /// [`check_three_term_recurrence_theorem`] with roots in `(-inf, 0)`.
    ThreeTermRecurrence,
    /// [`check_liu_wang_theorem`].
    LiuWang,
    /// [`check_liu_wang_derivative_lag_theorem`].
    LiuWangDerivativeLag,
    /// [`check_ma_wang_interval_theorem`] with an interval `[L, U]`,
    /// `U <= 0`.
    MaWangInterval,
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
    let mut attempts = Vec::new();
    let seed_count = record.seeds.len();
    for delay in 0..=MAXIMUM_ROOT_DOMAIN_START_DELAY {
        // The moved record needs its seeds among the exactly checked rows.
        let Some(seeds) = rows.get(delay..delay + seed_count) else {
            break;
        };
        let moved = ParentRecord {
            rule: record.rule,
            seeds,
            first_index: record.first_index + delay,
        };
        let (moved_attempts, success) = root_domain_menu(&moved, &rows[delay..], delay);
        attempts.extend(moved_attempts);
        if success {
            break;
        }
    }
    attempts
}

/// The root-domain menu for `record` (moved `delay` rows past the parent's
/// own seeds; `rows` start at the record's first index).
fn root_domain_menu(
    record: &ParentRecord<'_>,
    rows: &[Vec<Q>],
    delay: usize,
) -> (Vec<VeroneseParentProofAttempt>, bool) {
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
    let variant = |name: String| {
        if delay == 0 {
            name
        } else {
            format!("{name}, induction start moved {delay} row(s)")
        }
    };
    let moved_note = if delay == 0 {
        String::new()
    } else {
        format!(
            " (source moved {delay} row(s) past the parent seeds; the earlier parent rows are \
             among the rows checked exactly to be zero or real-rooted with coefficients of one \
             sign)"
        )
    };
    let mut attempts = Vec::new();
    for (name, transformation) in menu {
        let Some(transformation) = transformation else {
            attempts.push(VeroneseParentProofAttempt {
                route: VeroneseParentProofRoute::RootDomain,
                variant: variant(name),
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
                variant: variant(name),
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
                        variant: variant(name),
                        success: true,
                        detail: format!(
                            "source claim: every P_n, n >= {}, is real-rooted in {} ({:?}): \
                             {}{moved_note}",
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
                        variant: variant(name),
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
            return (attempts, true);
        }
    }
    (attempts, false)
}

/// Routes 4 to 7: a real-rootedness theorem with every root in
/// `(-inf, 0)` (routes 4 and 5) or, for route 6, in `(-inf, 0]` (every
/// conclusion of that theorem, so it runs without a requested region); route
/// 7 (Ma--Wang) runs without a requested region and succeeds when its
/// interval `[L, U]` has `U <= 0`.
fn theorem_attempt(
    record: &ParentRecord<'_>,
    route: VeroneseParentProofRoute,
    variant: &str,
) -> Result<VeroneseParentProofAttempt, RecurrenceError> {
    let options = RealRootednessTheoremOptions {
        root_region: Some(RealInterval {
            lower: None,
            lower_closed: false,
            upper: Some(Q::zero()),
            upper_closed: false,
        }),
        ..RealRootednessTheoremOptions::default()
    };
    let report = match route {
        VeroneseParentProofRoute::ThreeTermRecurrence => check_three_term_recurrence_theorem(
            record.rule,
            record.first_index,
            record.seeds,
            &options,
        )?,
        VeroneseParentProofRoute::LiuWangDerivativeLag => check_liu_wang_derivative_lag_theorem(
            record.rule,
            record.first_index,
            record.seeds,
            &RealRootednessTheoremOptions::default(),
        )?,
        VeroneseParentProofRoute::MaWangInterval => check_ma_wang_interval_theorem(
            record.rule,
            record.first_index,
            record.seeds,
            &RealRootednessTheoremOptions::default(),
        )?,
        _ => check_liu_wang_theorem(record.rule, record.first_index, record.seeds, &options)?,
    };
    // The Ma--Wang interval settles the parent statement only when it lies
    // in (-inf, 0].
    let outside_ray = route == VeroneseParentProofRoute::MaWangInterval
        && report.is_proved()
        && !report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref())
            .and_then(|region| region.at(0).ok())
            .and_then(|region| region.upper)
            .is_some_and(|upper| !upper.is_positive());
    let detail = match (&report.status, &report.conclusion) {
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if outside_ray => format!(
            "{}: {}; the interval is not contained in (-inf, 0], so the rows are real-rooted \
             but this route does not give coefficients of one sign",
            report.normal_form, conclusion.statement
        ),
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) if report.start_delay > 0 => {
            format!(
                "induction start moved {} row(s) past the seeds: {}: {}",
                report.start_delay, report.normal_form, conclusion.statement
            )
        }
        (RealRootednessTheoremStatus::Proved, Some(conclusion)) => {
            format!("{}: {}", report.normal_form, conclusion.statement)
        }
        (RealRootednessTheoremStatus::NotApplicable, _) => {
            format!("not applicable: {}", report.normal_form)
        }
        (status, _) => format!(
            "{status:?}: {}; open hypotheses: {}",
            report.normal_form,
            report.open_hypotheses().join("; ")
        ),
    };
    Ok(VeroneseParentProofAttempt {
        route,
        variant: variant.into(),
        success: report.is_proved() && !outside_ray,
        detail,
    })
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

/// A bivariate polynomial `sum a[i][j] v^i w^j` (dense, exact).
type Bivariate = Vec<Vec<Q>>;

fn bivariate_add(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let rows = left.len().max(right.len());
    (0..rows)
        .map(|i| {
            let a = left.get(i).map_or(&[][..], Vec::as_slice);
            let b = right.get(i).map_or(&[][..], Vec::as_slice);
            (0..a.len().max(b.len()))
                .map(|j| {
                    a.get(j).cloned().unwrap_or_default() + b.get(j).cloned().unwrap_or_default()
                })
                .collect()
        })
        .collect()
}

fn bivariate_mul(left: &Bivariate, right: &Bivariate) -> Bivariate {
    let mut out: Bivariate = Vec::new();
    for (i, row) in left.iter().enumerate() {
        for (j, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (k, other) in right.iter().enumerate() {
                for (l, b) in other.iter().enumerate() {
                    if b.is_zero() {
                        continue;
                    }
                    if out.len() <= i + k {
                        out.resize(i + k + 1, Vec::new());
                    }
                    if out[i + k].len() <= j + l {
                        out[i + k].resize(j + l + 1, Q::zero());
                    }
                    out[i + k][j + l] += a * b;
                }
            }
        }
    }
    out
}

/// `p(a + b v + c w)` for a univariate `p`, as a bivariate polynomial in `(v, w)`.
fn bivariate_affine_powers(a: &Q, b: &Q, c: &Q, degree: usize) -> Vec<Bivariate> {
    let affine: Bivariate = vec![vec![a.clone(), c.clone()], vec![b.clone()]];
    let mut powers = vec![vec![vec![Q::one()]]];
    for _ in 0..degree {
        let next = bivariate_mul(powers.last().expect("nonempty"), &affine);
        powers.push(next);
    }
    powers
}

/// Stirling numbers of the second kind `S(m, i)`, `m, i <= max`.
fn stirling_second_kind_table(max: usize) -> Vec<Vec<Q>> {
    let mut table = vec![vec![Q::zero(); max + 1]; max + 1];
    table[0][0] = Q::one();
    for m in 1..=max {
        for i in 1..=m {
            table[m][i] = &table[m - 1][i - 1] + Q::from_integer(i.into()) * &table[m - 1][i];
        }
    }
    table
}

/// Whether every coefficient of `p` in the basis `v^(falling i) w^(falling j)`
/// is `>= 0`; then `p(v, w) >= 0` at all integers `v, w >= 0`
/// (`x^m = sum_i S(m, i) x^(falling i)`).
fn has_nonnegative_falling_factorial_coefficients(p: &Bivariate) -> bool {
    let max = p
        .len()
        .max(p.iter().map(Vec::len).max().unwrap_or(0))
        .max(1);
    let stirling = stirling_second_kind_table(max);
    let mut falling = vec![vec![Q::zero(); max]; max];
    for (m, row) in p.iter().enumerate() {
        for (l, a) in row.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for i in 0..=m {
                for j in 0..=l {
                    let weight = &stirling[m][i] * &stirling[l][j];
                    if !weight.is_zero() {
                        falling[i][j] += a * weight;
                    }
                }
            }
        }
    }
    falling.iter().flatten().all(|value| !value.is_negative())
}

/// `x (x - 1) ... (x - d + 1)` in ascending powers.
fn falling_factorial_polynomial(d: usize) -> Vec<Q> {
    (0..d).fold(vec![Q::one()], |acc, k| {
        vec_poly::mul(&acc, &[-Q::from_integer(k.into()), Q::one()])
    })
}

/// Degree in `x` of `c(n, x)` (`None` for the zero polynomial).
fn x_degree(poly: &BivarPoly) -> Option<usize> {
    poly.coeffs
        .iter()
        .filter_map(|row| row.iter().rposition(|value| !value.is_zero()))
        .max()
}

/// The coefficient of `x^j` of `c(n, x)` as a polynomial in `n`.
fn x_coefficient(poly: &BivarPoly, j: usize) -> Vec<Q> {
    poly.coeffs
        .iter()
        .map(|row| row.get(j).cloned().unwrap_or_default())
        .collect()
}

fn integer_q(value: i64) -> Q {
    Q::from_integer(value.into())
}

/// Route 9: an exact proof that every row from `first_index` on has
/// coefficients of one sign (module documentation, "One-signed rows").
/// Returns what was shown, or why the sufficient test does not apply.
fn one_signed_rows_for_all_indices(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<String, String> {
    let start = first_index + seeds.len();
    let sign = seeds
        .iter()
        .flat_map(|seed| seed.iter())
        .find(|value| !value.is_zero())
        .map_or(Q::one(), |value| {
            if value.is_negative() {
                -Q::one()
            } else {
                Q::one()
            }
        });
    let inhomogeneous = rule
        .inhomogeneous
        .as_ref()
        .filter(|g| g.coeffs.iter().flatten().any(|value| !value.is_zero()));
    if sign.is_negative() && inhomogeneous.is_some() {
        return Err("the seeds are nonpositive but the rule is inhomogeneous".into());
    }
    if let Some(position) = seeds
        .iter()
        .position(|seed| seed.iter().any(|value| (value * &sign).is_negative()))
    {
        return Err(format!("seed {position} has coefficients of both signs"));
    }
    let positive_from_start = |polynomial: &[Q], condition: IndexSignCondition| {
        decide_index_polynomial_sign(
            polynomial,
            &IntegerRange::starting_at(BigInt::from(start)),
            condition,
        )
        .holds()
    };
    if let Some(leader) = &rule.denominator {
        if x_degree(leader).is_some_and(|degree| degree > 0) {
            return Err("the leader depends on x".into());
        }
        if !positive_from_start(&x_coefficient(leader, 0), IndexSignCondition::Positive) {
            return Err(format!("the leader is not positive for every n >= {start}"));
        }
    }
    if rule
        .terms
        .iter()
        .any(|term| term.sign == RecurrenceSign::AlternatingN)
    {
        return Err("a term carries the factor (-1)^n".into());
    }
    if let Some(g) = inhomogeneous {
        for j in 0..=x_degree(g).unwrap_or(0) {
            if !positive_from_start(&x_coefficient(g, j), IndexSignCondition::Nonnegative) {
                return Err(format!(
                    "the coefficient of x^{j} of the inhomogeneous term is not >= 0 for every \
                     n >= {start}"
                ));
            }
        }
    }

    // Degree bound D(m) = alpha + beta m for every row m >= first_index.
    let mut beta = 0i64;
    for term in &rule.terms {
        let Some(degree) = x_degree(&term.coeff) else {
            continue;
        };
        let excess = i64::try_from(degree).unwrap_or(i64::MAX)
            - i64::try_from(term.deriv_order).unwrap_or(i64::MAX);
        let offset = i64::try_from(term.offset).unwrap_or(i64::MAX).max(1);
        if excess > 0 {
            beta = beta.max((excess + offset - 1) / offset);
        }
    }
    let as_i64 = |value: usize| i64::try_from(value).unwrap_or(i64::MAX);
    let mut alpha: Option<i64> = None;
    for (position, seed) in seeds.iter().enumerate() {
        if let Some(degree) = vec_poly::degree(seed) {
            let needed = as_i64(degree) - beta * as_i64(first_index + position);
            alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
        }
    }
    if let Some(degree) = inhomogeneous.and_then(x_degree) {
        let needed = as_i64(degree) - beta * as_i64(start);
        alpha = Some(alpha.map_or(needed, |alpha| alpha.max(needed)));
    }
    let Some(alpha) = alpha else {
        return Ok("every seed is zero and the rule is homogeneous, so every row is zero".into());
    };

    // Combined weights W_(r,e)(n, s).
    let mut weights: std::collections::BTreeMap<(usize, i64), Bivariate> =
        std::collections::BTreeMap::new();
    for term in &rule.terms {
        let falling = falling_factorial_polynomial(term.deriv_order);
        for j in 0..=x_degree(&term.coeff).unwrap_or(0) {
            let c = x_coefficient(&term.coeff, j);
            if c.iter().all(Zero::is_zero) {
                continue;
            }
            let shift = as_i64(j) - as_i64(term.deriv_order);
            let product: Bivariate = c
                .iter()
                .map(|cn| falling.iter().map(|value| cn * value).collect())
                .collect();
            let entry = weights.entry((term.offset, shift)).or_default();
            *entry = bivariate_add(entry, &product);
        }
    }
    let start_q = integer_q(as_i64(start));
    for ((offset, shift), weight) in &weights {
        if weight.iter().flatten().all(Zero::is_zero) {
            continue;
        }
        let n_degree = weight.len().saturating_sub(1);
        let s_degree = weight
            .iter()
            .map(Vec::len)
            .max()
            .unwrap_or(1)
            .saturating_sub(1);
        let n_powers = bivariate_affine_powers(&start_q, &Q::one(), &Q::zero(), n_degree);
        let substitute = |s_powers: &[Bivariate]| {
            let mut total: Bivariate = Vec::new();
            for (a, row) in weight.iter().enumerate() {
                for (b, value) in row.iter().enumerate() {
                    if value.is_zero() {
                        continue;
                    }
                    let monomial = bivariate_mul(&n_powers[a], &s_powers[b]);
                    let scaled: Bivariate = monomial
                        .iter()
                        .map(|row| row.iter().map(|x| x * value).collect())
                        .collect();
                    total = bivariate_add(&total, &scaled);
                }
            }
            total
        };
        // (a) s >= 0 free.
        let free = bivariate_affine_powers(&Q::zero(), &Q::zero(), &Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&free)) {
            continue;
        }
        // (b) s = D(n - r) - u with n = start + v.
        let base = alpha + beta * (as_i64(start) - as_i64(*offset));
        let bounded =
            bivariate_affine_powers(&integer_q(base), &integer_q(beta), &-Q::one(), s_degree);
        if has_nonnegative_falling_factorial_coefficients(&substitute(&bounded)) {
            continue;
        }
        return Err(format!(
            "the weight of the coefficient of x^s of P_(n-{offset}) in the coefficient of \
             x^(s{shift:+}) is not shown to be >= 0 for n >= {start}"
        ));
    }
    Ok(format!(
        "every row has coefficients of one sign: the seeds do, the leader is positive and every \
         combined weight is >= 0 for n >= {start}, with deg P_m <= {alpha} + {beta} m"
    ))
}

/// Attestations from bundled catalog entries flagged `lean_proved` whose
/// first index and leading initial rows equal the parent's (a cheap filter;
/// [`same_record`] decides the exact match).
fn catalog_flag_attestations(record: &ParentRecord<'_>) -> Vec<VeroneseParentAttestation> {
    let wanted: Vec<Vec<String>> = record
        .seeds
        .iter()
        .map(|seed| {
            vec_poly::trim(seed)
                .iter()
                .map(format_rational_coeff)
                .collect()
        })
        .collect();
    let trimmed = |row: &[&str]| -> Vec<String> {
        let end = row
            .iter()
            .rposition(|value| *value != "0")
            .map_or(0, |i| i + 1);
        row[..end]
            .iter()
            .map(|value| (*value).to_string())
            .collect()
    };
    let zero_row = vec!["0".to_string()];
    let mut attestations = Vec::new();
    for (entry, flag) in entries_with_real_rootedness_flag() {
        if flag.status != OeisRealRootednessStatus::LeanProved
            || entry.recurrence_first_index != record.first_index
            || entry.initial_rows.len() < wanted.len()
        {
            continue;
        }
        let seeds_agree = entry.initial_rows.iter().zip(&wanted).all(|(row, seed)| {
            let row = trimmed(row);
            row == *seed || (row.is_empty() && *seed == zero_row)
        });
        if !seeds_agree {
            continue;
        }
        let Ok((rule, initial_rows)) = entry.recurrence_parts() else {
            continue;
        };
        let as_of = flag
            .as_of
            .map_or_else(String::new, |date| format!(", as of {date}"));
        attestations.push(VeroneseParentAttestation {
            label: format!("oeis_catalog:{}:lean_proved", entry.id),
            recurrence: RecurrenceJson::from_recurrence_rational(
                &rule,
                entry.recurrence_first_index,
                &initial_rows,
                None,
            ),
            first_index: entry.recurrence_first_index,
            statement: format!(
                "every row of the bundled OEIS catalog entry {} is zero or real-rooted (flag \
                 lean_proved{as_of}: a machine-checked Lean proof made outside polytool)",
                entry.id
            ),
            source: format!(
                "polytool OEIS catalog, real_rootedness flag of {}",
                entry.id
            ),
        });
    }
    attestations
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

    // Routes 4 to 7.
    if proved.is_none() {
        for (route, variant) in [
            (
                VeroneseParentProofRoute::ThreeTermRecurrence,
                "three_term_recurrence",
            ),
            (VeroneseParentProofRoute::LiuWang, "liu_wang"),
            (
                VeroneseParentProofRoute::LiuWangDerivativeLag,
                "liu_wang_derivative_lag",
            ),
            (VeroneseParentProofRoute::MaWangInterval, "ma_wang_interval"),
        ] {
            let attempt = theorem_attempt(record, route, variant)?;
            let success = attempt.success;
            proof.attempts.push(attempt);
            if success {
                proved = Some(route);
                break;
            }
        }
    }

    if let Some(route) = proved {
        proof.status = VeroneseParentProofStatus::Proved;
        proof.route = Some(route);
        proof.statement = conclusion;
        return Ok(proof);
    }

    // Route 8.
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
    // Route 9.
    let candidates = catalog_flag_attestations(record);
    let mut differing = Vec::new();
    for attestation in candidates {
        if !same_record(record, &attestation.recurrence)? {
            differing.push(attestation.label);
            continue;
        }
        let covers = attestation.first_index <= needed_from;
        let one_signed = covers.then(|| {
            one_signed_rows_for_all_indices(record.rule, record.first_index, record.seeds)
        });
        let success = matches!(one_signed, Some(Ok(_)));
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: attestation.label.clone(),
            success,
            detail: match &one_signed {
                Some(Ok(signs)) => format!("{}; {signs}", attestation.statement),
                Some(Err(reason)) => format!(
                    "the parent equals the flagged entry, but the flag gives real-rootedness \
                     only and one-signed coefficients are not shown: {reason}"
                ),
                None => format!(
                    "the flagged entry starts at index {}, after the first used parent index \
                     {needed_from}",
                    attestation.first_index
                ),
            },
        });
        if success {
            proof.status = VeroneseParentProofStatus::ConditionalOnAttestation;
            proof.route = Some(VeroneseParentProofRoute::Attestation);
            proof.statement = format!(
                "{conclusion} (real-rootedness from the catalog flag, not proved here; \
                 one-signed coefficients proved here)"
            );
            proof.conditional_on = vec![attestation.label.clone()];
            return Ok(proof);
        }
    }
    if !differing.is_empty() {
        proof.attempts.push(VeroneseParentProofAttempt {
            route: VeroneseParentProofRoute::Attestation,
            variant: "oeis_catalog".into(),
            success: false,
            detail: format!(
                "no lean_proved catalog entry equals the parent; {} with the same first index \
                 and seeds differ in the rule or in later initial rows ({})",
                differing.len(),
                differing.join(", ")
            ),
        });
    }
    proof.statement = format!(
        "open: no route proved that {conclusion}; the first {} parent rows pass the exact \
         check (evidence only)",
        proof.rows_checked
    );
    Ok(proof)
}

#[cfg(test)]
mod catalog_flag_tests {
    use super::*;

    fn q(value: i64) -> Q {
        Q::from_integer(value.into())
    }

    /// The catalog rule, identifying seeds and first index of a bundled entry.
    fn catalog_parent(id: &str) -> (Recurrence, Vec<Vec<Q>>, usize) {
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
        (rule, rows[..count].to_vec(), entry.recurrence_first_index)
    }

    fn prove(rule: &Recurrence, seeds: &[Vec<Q>], first_index: usize) -> VeroneseParentProof {
        let record = ParentRecord {
            rule,
            seeds,
            first_index,
        };
        prove_parent(&record, first_index, &[], false).unwrap()
    }

    fn catalog_attempts(proof: &VeroneseParentProof) -> Vec<&VeroneseParentProofAttempt> {
        proof
            .attempts
            .iter()
            .filter(|attempt| attempt.variant.starts_with("oeis_catalog:"))
            .collect()
    }

    /// The Lucas polynomials A034807 (`P_n = P_(n-1) + x P_(n-2)`, `P_0 = 2`,
    /// `P_1 = 1`) are flagged `lean_proved`; no polytool route proves them,
    /// so the verdict is conditional on the flag, with one-signed rows proved.
    #[test]
    fn a_flagged_catalog_parent_is_conditional_on_the_flag() {
        let (rule, seeds, first_index) = catalog_parent("A034807");
        assert_eq!(seeds, vec![vec![q(2)], vec![q(1)]]);
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        assert_eq!(proof.route, Some(VeroneseParentProofRoute::Attestation));
        assert_eq!(proof.conditional_on, ["oeis_catalog:A034807:lean_proved"]);
        assert!(proof.statement.contains("not proved here"));
        let attempt = proof.attempts.last().unwrap();
        assert!(attempt.success);
        assert!(attempt.detail.contains("lean_proved, as of 2026-10-03"));
        assert!(attempt.detail.contains("coefficients of one sign"));
        // A caller attestation of the same record is applied first.
        let attestation = VeroneseParentAttestation {
            label: "caller".into(),
            recurrence: RecurrenceJson::from_recurrence_rational(&rule, first_index, &seeds, None),
            first_index,
            statement: "attested".into(),
            source: "test".into(),
        };
        let record = ParentRecord {
            rule: &rule,
            seeds: &seeds,
            first_index,
        };
        let proof = prove_parent(&record, first_index, &[attestation], false).unwrap();
        assert_eq!(proof.conditional_on, ["caller"]);
    }

    /// A157153 is not flagged (its Lean proof is open), so no catalog
    /// attestation applies.
    #[test]
    fn an_unflagged_catalog_parent_gets_no_catalog_attestation() {
        assert!(crate::oeis::by_id("A157153")
            .unwrap()
            .real_rootedness()
            .is_none());
        let (rule, seeds, first_index) = catalog_parent("A157153");
        let proof = prove(&rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation,
            "{proof:#?}"
        );
        // Flagged entries with the same seeds differ in the rule.
        assert!(catalog_attempts(&proof).is_empty(), "{proof:#?}");
        assert!(proof
            .attempts
            .iter()
            .filter(|attempt| attempt.route == VeroneseParentProofRoute::Attestation)
            .all(|attempt| !attempt.success && !attempt.detail.contains("A157153")));
    }

    /// Matching is exact: one changed rule coefficient or seed and the flag no
    /// longer applies.
    #[test]
    fn a_parent_differing_in_one_coefficient_or_seed_does_not_match() {
        let (rule, seeds, first_index) = catalog_parent("A034807");
        let mut changed_rule = rule.clone();
        changed_rule.terms[1].coeff.coeffs[0][1] = q(2);
        let proof = prove(&changed_rule, &seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        let attempts = catalog_attempts(&proof);
        assert!(attempts.is_empty(), "{proof:#?}");
        let summary = proof.attempts.last().unwrap();
        assert_eq!(summary.variant, "oeis_catalog");
        assert!(!summary.success);
        assert!(summary.detail.contains("oeis_catalog:A034807:lean_proved"));

        let mut changed_seeds = seeds.clone();
        changed_seeds[1] = vec![q(2)];
        let proof = prove(&rule, &changed_seeds, first_index);
        assert_ne!(
            proof.status,
            VeroneseParentProofStatus::ConditionalOnAttestation
        );
        assert!(catalog_attempts(&proof).is_empty());

        // The same rule and seeds at another first index do not match either.
        let proof = prove(&rule, &seeds, first_index + 1);
        assert!(catalog_attempts(&proof).is_empty());
    }

    /// The flag gives real-rootedness only: a flagged parent whose rule
    /// cancels (A003506, `P_n = (2 + 2x) P_(n-1) - (1 + x)^2 P_(n-2)`) is not
    /// shown one-signed, so the catalog route fails.
    #[test]
    fn a_flagged_parent_without_a_one_sign_proof_stays_open() {
        let (rule, seeds, first_index) = catalog_parent("A003506");
        assert!(crate::oeis::by_id("A003506")
            .unwrap()
            .real_rootedness()
            .is_some());
        let proof = prove(&rule, &seeds, first_index);
        assert_eq!(proof.status, VeroneseParentProofStatus::NotProved);
        let attempts = catalog_attempts(&proof);
        assert_eq!(attempts.len(), 1);
        assert!(attempts[0]
            .detail
            .contains("one-signed coefficients are not shown: the weight"));
    }

    fn term(offset: usize, derivative: usize, coefficients: Vec<Vec<Q>>) -> RecurrenceTerm {
        RecurrenceTerm::new(
            offset,
            derivative,
            BivarPoly::from_coefficients_normalized(coefficients),
        )
    }

    #[test]
    fn one_signed_rows_need_nonnegative_combined_weights() {
        // Eulerian: P_n = (1 - x + n x) P_(n-1) + (x - x^2) P'_(n-1), P_1 = 1.
        // The weight n - 1 - s of x^s into x^(s+1) needs the degree bound
        // deg P_m <= m - 1.
        let eulerian = Recurrence {
            terms: vec![
                term(1, 0, vec![vec![q(1), q(-1)], vec![q(0), q(1)]]),
                term(1, 1, vec![vec![q(0), q(1), q(-1)]]),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        let shown = one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1)]]).unwrap();
        assert!(shown.ends_with("deg P_m <= -1 + 1 m"), "{shown}");
        // Nonpositive seeds of a homogeneous rule: the rows are nonpositive.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(-1)]]).is_ok());
        assert!(
            one_signed_rows_for_all_indices(&eulerian, 1, &[vec![q(1), q(-1)]])
                .unwrap_err()
                .contains("both signs")
        );
        // A coefficient x - 2 x^2 makes the weight n - 2 - 2 s, negative at
        // s = n - 2 (the top coefficient of P_(n-1)).
        let mut too_negative = eulerian.clone();
        too_negative.terms[1] = term(1, 1, vec![vec![q(0), q(1), q(-2)]]);
        assert!(
            one_signed_rows_for_all_indices(&too_negative, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("not shown to be >= 0")
        );
        // A leader 3 - n vanishes at n = 3.
        let mut leader = eulerian.clone();
        leader.denominator = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(3)],
            vec![q(-1)],
        ]));
        assert!(one_signed_rows_for_all_indices(&leader, 1, &[vec![q(1)]])
            .unwrap_err()
            .contains("leader is not positive"));
        // (-1)^n factors and negative inhomogeneous parts are refused.
        let mut alternating = eulerian.clone();
        alternating.terms[0].sign = RecurrenceSign::AlternatingN;
        assert!(one_signed_rows_for_all_indices(&alternating, 1, &[vec![q(1)]]).is_err());
        let mut inhomogeneous = eulerian.clone();
        inhomogeneous.inhomogeneous = Some(BivarPoly::from_coefficients_normalized(vec![
            vec![q(5)],
            vec![q(-1)],
        ]));
        assert!(
            one_signed_rows_for_all_indices(&inhomogeneous, 1, &[vec![q(1)]])
                .unwrap_err()
                .contains("inhomogeneous")
        );
        // Zero seeds of a homogeneous rule give zero rows.
        assert!(one_signed_rows_for_all_indices(&eulerian, 1, &[vec![]]).is_ok());
    }
}
