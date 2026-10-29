//! Conditional theorem matches for found or supplied recurrences.
//!
//! [`match_recurrence_theorems`] checks a [`Recurrence`] against the
//! interlacing theorems of this module family and returns explicit
//! *conditional* matches: the hypotheses are verified exactly, the conclusion
//! and interlacing window are stated, and supplied rows are compared with the
//! canonical boundary rows.  A match never upgrades finite fit evidence: it
//! does not show that the recurrence governs the intended sequence for all
//! `n`, which stays listed in `remaining_obligations`.
//!
//! Theorems (paper proofs in
//! <https://github.com/PerAlexandersson/RealRooted/issues/2037>, since
//! formalized in the RealRooted Lean library:
//! `RealRooted/Interlacing/VariableDelay.lean` and
//! `RealRooted/CombinatorialExamples/WeightedRenewal.lean`; the matches
//! themselves are not replayed in Lean, so `proof_status` stays
//! `paper_proof_not_lean_checked`):
//!
//! - `weighted_renewal_window`: `P_n = sum_j r_j t^j P_(n-1-ell*j)`, see
//!   [`super::recognize_weighted_renewal`];
//! - `fixed_delay_predecessor_window`: `P_n = P_(n-1) + w(n) t P_(n-L)` with a
//!   polynomial weight, see [`super::recognize_fixed_delay_predecessor`].
//!
//! [`RecurrenceTheoremReport`] (schema `polytool.recurrence-theorems.v1`)
//! wraps the matches with the recurrence record and its finite evidence; it
//! is what the CLI (`recurrence --theorem-matches`, `recurrence-theorems`) and
//! the MCP server emit.

use super::predecessor::{
    recognize_fixed_delay_predecessor, VariablePredecessorCertificate,
    VariablePredecessorRecurrence, WeightCertificate, PAPER_PROOF_STATUS,
    VARIABLE_PREDECESSOR_REFERENCE,
};
use super::{recognize_weighted_renewal, RenewalError, WeightedRenewal};
use crate::recurrence::{
    format_rational_coeff, AdaptiveSearchResult, BigRational, Recurrence, RecurrenceJson,
};
use num_traits::Zero;
use serde::{Deserialize, Serialize};

/// JSON schema tag of [`RecurrenceTheoremReport`].
pub const RECURRENCE_THEOREM_SCHEMA: &str = "polytool.recurrence-theorems.v1";

/// Theorem identifier of the stationary weighted renewal window.
pub const WEIGHTED_RENEWAL_THEOREM: &str = "weighted_renewal_window";

/// Theorem identifier of the fixed-delay variable-predecessor window.
pub const FIXED_DELAY_PREDECESSOR_THEOREM: &str = "fixed_delay_predecessor_window";

/// Status string of every theorem match.
pub const CONDITIONAL_THEOREM_MATCH: &str = "conditional_theorem_match";

/// Boundary rows are compared only for supplied labels up to this index.
const MAX_BOUNDARY_INDEX: usize = 1 << 12;

/// How the supplied rows relate to the canonical boundary rows.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TheoremBoundaryStatus {
    /// Every supplied row equals the canonical row with the same index.
    Canonical,
    /// Every supplied row equals the canonical row after relabeling the
    /// first supplied row as index 0 (supplied index = canonical index +
    /// `index_shift`).  For an index-dependent weight the relabeled family
    /// has weight `w(n + index_shift)`, which inherits nonnegativity.
    CanonicalAfterRelabel,
    /// Some supplied row differs from the canonical rows.
    Mismatch,
    /// Fewer rows than the recurrence lag were supplied.
    InsufficientRows,
}

/// Finite comparison of supplied rows with the canonical boundary contract.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TheoremBoundaryCheck {
    /// The canonical boundary/initial-data contract of the theorem.
    pub contract: String,
    /// Result of the comparison.
    pub status: TheoremBoundaryStatus,
    /// True for `canonical` and `canonical_after_relabel`.
    pub rows_match_canonical: bool,
    /// Index of the first supplied row (as labeled by the record).
    pub supplied_first_index: usize,
    /// Number of supplied rows compared.
    pub supplied_rows: usize,
    /// Rows needed to determine the recurrence (its maximal lag).
    pub required_rows: usize,
    /// For `canonical_after_relabel`: supplied index minus canonical index.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub index_shift: Option<usize>,
    /// For `mismatch`: supplied index of the first differing row (relabeled
    /// comparison when `supplied_first_index > 0`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub first_mismatch_index: Option<usize>,
    /// Always `finite_supplied_rows`: only the supplied rows were compared.
    pub scope: String,
}

/// One hypothesis of a matched theorem; every listed hypothesis holds.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TheoremHypothesis {
    /// Short identifier.
    pub name: String,
    /// The hypothesis as instantiated for this recurrence.
    pub statement: String,
    /// How it was verified: `exact_structural`, `exact_rational`,
    /// `exact_root_count` or `exact_symbolic_certificate`.
    pub check: String,
}

/// The interlacing window of a matched theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TheoremWindow {
    /// The window inequality with its parameters, e.g. `2*s <= 2*1+1`.
    pub inequality: String,
    /// Largest index step `s` with `P_n << P_(n+s)`; `None` means every step
    /// (constant family).
    pub max_step: Option<usize>,
    /// Plain-text window statement.
    pub statement: String,
}

/// Parameters of a weighted renewal match.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct WeightedRenewalParameters {
    /// Delay `ell >= 1`.
    pub ell: usize,
    /// Degree `d` of the seed `R(z)`.
    pub degree: usize,
    /// Seed coefficients `r_0 = 1, r_1, ..., r_d` as exact strings.
    pub seed: Vec<String>,
}

/// An explicit conditional theorem match.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecurrenceTheoremMatch {
    /// `weighted_renewal_window` or `fixed_delay_predecessor_window`.
    pub theorem: String,
    /// Always `conditional_theorem_match`.
    pub status: String,
    /// The canonical family the theorem is about, with its boundary.
    pub family: String,
    /// Exactly checked hypotheses.
    pub hypotheses: Vec<TheoremHypothesis>,
    /// Conclusion for the canonical family.
    pub conclusion: String,
    /// Interlacing window.
    pub window: TheoremWindow,
    /// Finite comparison of the supplied rows with the canonical rows.
    pub boundary: TheoremBoundaryCheck,
    /// What must still be shown before the conclusion applies to the
    /// intended sequence.  Never empty.
    pub remaining_obligations: Vec<String>,
    /// Always `paper_proof_not_lean_checked`.
    pub proof_status: String,
    /// Source of the paper proof.
    pub reference: String,
    /// Present for `weighted_renewal_window`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub weighted_renewal: Option<WeightedRenewalParameters>,
    /// Present for `fixed_delay_predecessor_window`: the replayable
    /// hypotheses certificate.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub variable_predecessor: Option<VariablePredecessorCertificate>,
}

/// A theorem that did not match, with the reason.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecurrenceTheoremRejection {
    /// Theorem identifier.
    pub theorem: String,
    /// Why the recognizer or a hypothesis check failed.
    pub reason: String,
}

/// Matches and rejections for one recurrence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecurrenceTheoremMatches {
    /// Conditional matches (possibly several).
    pub theorem_matches: Vec<RecurrenceTheoremMatch>,
    /// Theorems that did not apply.
    pub rejected_theorems: Vec<RecurrenceTheoremRejection>,
}

/// Where the recurrence came from and what finite evidence supports it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecurrenceTheoremEvidence {
    /// `adaptive_search_fit` or `supplied_recurrence`.
    pub source: String,
    /// Always `finite_evidence_only`.
    pub status: String,
    /// Index label of the first compared row.
    pub first_index: usize,
    /// Number of rows compared with the canonical boundary rows.
    pub rows_compared: usize,
    /// Rows used to fit the recurrence (search only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fit_polynomials: Option<usize>,
    /// Held-out rows verified after fitting (search only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub verification_polynomials: Option<usize>,
    /// Plain-text caveat.
    pub note: String,
}

/// `polytool.recurrence-theorems.v1`: a recurrence record, its finite
/// evidence, and conditional theorem matches, kept in separate fields.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RecurrenceTheoremReport {
    /// Always `polytool.recurrence-theorems.v1`.
    pub schema: String,
    /// The recurrence as text.
    pub recurrence_text: String,
    /// The analysed recurrence (`polytool.recurrence.v1`), reusable for row
    /// generation.
    pub recurrence_json: RecurrenceJson,
    /// Finite fit or supply evidence; never a proof of the recurrence.
    pub finite_evidence: RecurrenceTheoremEvidence,
    /// Conditional theorem matches.
    pub theorem_matches: Vec<RecurrenceTheoremMatch>,
    /// Theorems that did not apply, with reasons.
    pub rejected_theorems: Vec<RecurrenceTheoremRejection>,
}

fn strip(row: &[BigRational]) -> &[BigRational] {
    let mut end = row.len();
    while end > 0 && row[end - 1].is_zero() {
        end -= 1;
    }
    &row[..end]
}

/// First `k` with `rows[k] != expected[k]` (trailing zeros ignored).
fn first_mismatch(expected: &[Vec<BigRational>], rows: &[Vec<BigRational>]) -> Option<usize> {
    rows.iter()
        .zip(expected)
        .position(|(actual, expected)| strip(actual) != strip(expected))
}

/// Compare supplied rows with canonical rows at their labels and, if that
/// fails and `first_index > 0`, after relabeling the first row as index 0.
fn boundary_check(
    contract: String,
    required_rows: usize,
    first_index: usize,
    rows: &[Vec<BigRational>],
    canonical: impl Fn(usize) -> Option<Vec<Vec<BigRational>>>,
    relabeled: impl Fn(usize) -> Option<Vec<Vec<BigRational>>>,
) -> TheoremBoundaryCheck {
    let mut check = TheoremBoundaryCheck {
        contract,
        status: TheoremBoundaryStatus::InsufficientRows,
        rows_match_canonical: false,
        supplied_first_index: first_index,
        supplied_rows: rows.len(),
        required_rows,
        index_shift: None,
        first_mismatch_index: None,
        scope: "finite_supplied_rows".into(),
    };
    if rows.len() < required_rows {
        return check;
    }
    let exact = (first_index <= MAX_BOUNDARY_INDEX)
        .then(|| canonical(first_index + rows.len()))
        .flatten()
        .map(|expected| first_mismatch(&expected[first_index..], rows));
    if exact == Some(None) {
        check.status = TheoremBoundaryStatus::Canonical;
        check.rows_match_canonical = true;
        return check;
    }
    let mismatch = if first_index > 0 {
        match relabeled(rows.len()).map(|expected| first_mismatch(&expected, rows)) {
            Some(None) => {
                check.status = TheoremBoundaryStatus::CanonicalAfterRelabel;
                check.rows_match_canonical = true;
                check.index_shift = Some(first_index);
                return check;
            }
            Some(Some(k)) => Some(k),
            None => exact.flatten(),
        }
    } else {
        exact.flatten()
    };
    check.status = TheoremBoundaryStatus::Mismatch;
    check.first_mismatch_index = mismatch.map(|k| first_index + k);
    check
}

fn obligations(boundary: &TheoremBoundaryCheck, relabel_note: &str) -> Vec<String> {
    let mut out = vec![
        "Show that the intended sequence satisfies this recurrence for every n; agreement \
         with finitely many rows (a fit) is evidence, not proof."
            .to_string(),
    ];
    match boundary.status {
        TheoremBoundaryStatus::Canonical => {}
        TheoremBoundaryStatus::CanonicalAfterRelabel => out.push(format!(
            "The conclusion is about the relabeled sequence: supplied index = canonical \
             index + {}.{relabel_note}",
            boundary.index_shift.unwrap_or(0)
        )),
        TheoremBoundaryStatus::Mismatch => out.push(
            "The supplied rows differ from the canonical boundary rows, so the theorem \
             applies to the canonical family, not to these rows."
                .into(),
        ),
        TheoremBoundaryStatus::InsufficientRows => out.push(format!(
            "Supply at least {} consecutive initial rows to compare with the canonical \
             boundary rows.",
            boundary.required_rows
        )),
    }
    out
}

fn monomial(c: &BigRational, j: usize, row: &str) -> String {
    let coefficient = if c == &BigRational::from_integer(1.into()) {
        String::new()
    } else {
        format!("{} ", format_rational_coeff(c))
    };
    match j {
        0 => format!("{coefficient}{row}"),
        1 => format!("{coefficient}t {row}"),
        _ => format!("{coefficient}t^{j} {row}"),
    }
}

fn renewal_match(
    family: &WeightedRenewal,
    first_index: usize,
    rows: &[Vec<BigRational>],
) -> RecurrenceTheoremMatch {
    let ell = family.ell();
    let seed = family.seed();
    let d = seed.len() - 1;
    let rhs = seed
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(j, c)| monomial(c, j, &format!("P_(n-{})", 1 + ell * j)))
        .collect::<Vec<_>>()
        .join(" + ");
    let seed_text = seed
        .iter()
        .enumerate()
        .filter(|(_, c)| !c.is_zero())
        .map(|(j, c)| match j {
            0 => format_rational_coeff(c),
            1 => format!("{} z", format_rational_coeff(c)),
            _ => format!("{} z^{j}", format_rational_coeff(c)),
        })
        .collect::<Vec<_>>()
        .join(" + ");
    let contract = "P_0 = 1 and P_m = 0 for m < 0; the recurrence holds for every n >= 1".into();
    let required = 1 + ell * d;
    let boundary = boundary_check(
        contract,
        required,
        first_index,
        rows,
        |count| Some(family.rows(count)),
        |count| Some(family.rows(count)),
    );
    let (conclusion, window) = match family.max_interlacing_step() {
        None => (
            "P_n = 1 for every n >= 0 (constant family); any two rows interlace trivially."
                .to_string(),
            TheoremWindow {
                inequality: "R = 1: no restriction".into(),
                max_step: None,
                statement: "P_n << P_(n+s) for every s >= 0".into(),
            },
        ),
        Some(step) => (
            format!(
                "For all n, s >= 0 with {d}*s <= {d}*{ell}+1: P_n << P_(n+s) (oriented weak \
                 interlacing; shared and repeated zeros allowed). Every P_n has nonnegative \
                 coefficients and constant term 1 and is real-rooted."
            ),
            TheoremWindow {
                inequality: format!("{d}*s <= {d}*{ell}+1 (d*s <= d*ell+1)"),
                max_step: Some(step),
                statement: format!("P_n << P_(n+s) for every n >= 0 and 0 <= s <= {step}"),
            },
        ),
    };
    let remaining_obligations = obligations(
        &boundary,
        " The recurrence is index-independent, so relabeling preserves it.",
    );
    RecurrenceTheoremMatch {
        theorem: WEIGHTED_RENEWAL_THEOREM.into(),
        status: CONDITIONAL_THEOREM_MATCH.into(),
        family: format!(
            "P_n = {rhs}; P_0 = 1, P_m = 0 for m < 0; generating function 1/(1 - x R(t x^{ell}))"
        ),
        hypotheses: vec![
            TheoremHypothesis {
                name: "renewal_structure".into(),
                statement: format!(
                    "homogeneous, index-independent, no derivatives or sign factors; every \
                     t^j term has lag 1 + {ell}*j with ell = {ell} >= 1"
                ),
                check: "exact_structural".into(),
            },
            TheoremHypothesis {
                name: "seed_nonnegative".into(),
                statement: format!("R(z) = {seed_text} has R(0) = 1 and nonnegative coefficients"),
                check: "exact_rational".into(),
            },
            TheoremHypothesis {
                name: "seed_real_rooted".into(),
                statement: "R(z) has only real zeros".into(),
                check: "exact_root_count".into(),
            },
        ],
        conclusion,
        window,
        boundary,
        remaining_obligations,
        proof_status: PAPER_PROOF_STATUS.into(),
        reference: VARIABLE_PREDECESSOR_REFERENCE.into(),
        weighted_renewal: Some(WeightedRenewalParameters {
            ell,
            degree: d,
            seed: seed.iter().map(format_rational_coeff).collect(),
        }),
        variable_predecessor: None,
    }
}

fn predecessor_match(
    family: &VariablePredecessorRecurrence,
    certificate: VariablePredecessorCertificate,
    first_index: usize,
    rows: &[Vec<BigRational>],
) -> RecurrenceTheoremMatch {
    let delay = family.interlacing_window().unwrap_or(0);
    let contract = format!(
        "Q_0 = 1 and Q_n = Q_(n-1) + w(n) t Q_(max(0, n-{delay})) for n >= 1: for n < {delay} \
         the predecessor is Q_0, not a zero row"
    );
    let boundary = boundary_check(
        contract,
        delay,
        first_index,
        rows,
        |count| family.rows(count).ok(),
        |count| family.with_index_shift(first_index)?.rows(count).ok(),
    );
    let remaining_obligations = obligations(
        &boundary,
        " The relabeled family has weight w(n + shift), nonnegative for n >= 1.",
    );
    let weight_statement = match &certificate.weight {
        WeightCertificate::Polynomial {
            shift,
            shifted_coefficients,
            ..
        } => format!(
            "w(n) >= 0 for every integer n >= 1: w({shift}+m) has nonnegative coefficients \
             [{}] in m{}",
            shifted_coefficients.join(", "),
            if *shift > 1 {
                format!(
                    ", and w(1), ..., w({}) >= 0 were evaluated exactly",
                    shift - 1
                )
            } else {
                String::new()
            }
        ),
        _ => "w(n) >= 0 for every integer n >= 1 (see variable_predecessor.weight)".into(),
    };
    RecurrenceTheoremMatch {
        theorem: FIXED_DELAY_PREDECESSOR_THEOREM.into(),
        status: CONDITIONAL_THEOREM_MATCH.into(),
        family: family.describe(),
        hypotheses: vec![
            TheoremHypothesis {
                name: "predecessor".into(),
                statement: format!(
                    "p(n) = max(0, n-{delay}) with {delay} >= 1 is nondecreasing and \
                     0 <= p(n) < n for n >= 1"
                ),
                check: "exact_structural".into(),
            },
            TheoremHypothesis {
                name: "weight_nonnegative".into(),
                statement: weight_statement,
                check: "exact_symbolic_certificate".into(),
            },
        ],
        conclusion: certificate.conclusion.clone(),
        window: TheoremWindow {
            inequality: format!("B - A <= {delay}"),
            max_step: Some(delay),
            statement: format!("Q_A << Q_B for 0 <= A <= B with B - A <= {delay}"),
        },
        boundary,
        remaining_obligations,
        proof_status: PAPER_PROOF_STATUS.into(),
        reference: VARIABLE_PREDECESSOR_REFERENCE.into(),
        weighted_renewal: None,
        variable_predecessor: Some(certificate),
    }
}

/// Check `recurrence` against every supported theorem.
///
/// `rows` are supplied rows labeled `first_index, first_index+1, ...`; they
/// are compared with the canonical boundary rows (finite check only).  The
/// hypotheses are verified exactly.  Every match is conditional: see
/// [`RecurrenceTheoremMatch::remaining_obligations`].
pub fn match_recurrence_theorems(
    recurrence: &Recurrence,
    first_index: usize,
    rows: &[Vec<BigRational>],
) -> RecurrenceTheoremMatches {
    let mut out = RecurrenceTheoremMatches {
        theorem_matches: Vec::new(),
        rejected_theorems: Vec::new(),
    };
    let reject = |theorem: &str, reason: String| RecurrenceTheoremRejection {
        theorem: theorem.into(),
        reason,
    };
    match recognize_weighted_renewal(recurrence) {
        Ok(family) => out
            .theorem_matches
            .push(renewal_match(&family, first_index, rows)),
        Err(reason) => out
            .rejected_theorems
            .push(reject(WEIGHTED_RENEWAL_THEOREM, reason.to_string())),
    }
    match recognize_fixed_delay_predecessor(recurrence)
        .and_then(|family| family.certify().map(|cert| (family, cert)))
    {
        Ok((family, certificate)) => {
            out.theorem_matches
                .push(predecessor_match(&family, certificate, first_index, rows))
        }
        Err(reason) => out
            .rejected_theorems
            .push(reject(FIXED_DELAY_PREDECESSOR_THEOREM, reason.to_string())),
    }
    out
}

impl RecurrenceTheoremMatch {
    /// One-line summary, e.g. for CLI diagnostics.
    pub fn summary(&self) -> String {
        let boundary = match self.boundary.status {
            TheoremBoundaryStatus::Canonical => "canonical".to_string(),
            TheoremBoundaryStatus::CanonicalAfterRelabel => format!(
                "canonical after relabel (shift {})",
                self.boundary.index_shift.unwrap_or(0)
            ),
            TheoremBoundaryStatus::Mismatch => "mismatch".to_string(),
            TheoremBoundaryStatus::InsufficientRows => "insufficient rows".to_string(),
        };
        format!(
            "conditional theorem match {}: window {} ({}); supplied boundary rows: {}; \
             {} (all-n validity of the recurrence is not established)",
            self.theorem,
            self.window.inequality,
            self.window.statement,
            boundary,
            self.proof_status
        )
    }
}

impl RecurrenceTheoremReport {
    /// Report for a caller-supplied recurrence record; its initial rows are
    /// compared with the canonical boundary rows.
    pub fn for_supplied(record: &RecurrenceJson) -> Result<Self, RenewalError> {
        let (recurrence, first_index, rows) = record.to_recurrence_parts()?;
        let matches = match_recurrence_theorems(&recurrence, first_index, &rows);
        Ok(Self {
            schema: RECURRENCE_THEOREM_SCHEMA.into(),
            recurrence_text: recurrence.to_string(),
            recurrence_json: record.clone(),
            finite_evidence: RecurrenceTheoremEvidence {
                source: "supplied_recurrence".into(),
                status: "finite_evidence_only".into(),
                first_index,
                rows_compared: rows.len(),
                fit_polynomials: None,
                verification_polynomials: None,
                note: "Supplied by the caller; only the record's initial rows were compared. \
                       That the recurrence governs the intended family is the caller's \
                       obligation."
                    .into(),
            },
            theorem_matches: matches.theorem_matches,
            rejected_theorems: matches.rejected_theorems,
        })
    }

    /// Report for an adaptive-search result.  `searched_rows` are the rows
    /// after `skip_prefix`, labeled from `recurrence_json.first_index`.
    pub fn for_search(
        result: &AdaptiveSearchResult,
        recurrence_json: RecurrenceJson,
        searched_rows: &[Vec<BigRational>],
    ) -> Self {
        let first_index = recurrence_json.first_index;
        let matches = match_recurrence_theorems(&result.recurrence, first_index, searched_rows);
        Self {
            schema: RECURRENCE_THEOREM_SCHEMA.into(),
            recurrence_text: result.recurrence.to_string(),
            recurrence_json,
            finite_evidence: RecurrenceTheoremEvidence {
                source: "adaptive_search_fit".into(),
                status: "finite_evidence_only".into(),
                first_index,
                rows_compared: searched_rows.len(),
                fit_polynomials: Some(result.fit_polynomials),
                verification_polynomials: Some(result.verification_polynomials),
                note: "The recurrence was fitted to finitely many rows and checked on held-out \
                       rows. This does not prove it for all n."
                    .into(),
            },
            theorem_matches: matches.theorem_matches,
            rejected_theorems: matches.rejected_theorems,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::{BivarPoly, RecurrenceTerm};

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(n.into())
    }
    fn renewal(ell: usize, seed: &[i64]) -> Recurrence {
        Recurrence {
            denominator: None,
            inhomogeneous: None,
            terms: seed
                .iter()
                .enumerate()
                .map(|(j, c)| {
                    let mut row = vec![q(0); j + 1];
                    row[j] = q(*c);
                    RecurrenceTerm::new(1 + ell * j, 0, BivarPoly { coeffs: vec![row] })
                })
                .collect(),
        }
    }

    #[test]
    fn renewal_match_reports_window_and_boundary() {
        let rec = renewal(1, &[1, 3, 2]);
        let rows = WeightedRenewal::new(1, vec![q(1), q(3), q(2)])
            .unwrap()
            .rows(8);
        let matches = match_recurrence_theorems(&rec, 0, &rows);
        assert_eq!(matches.theorem_matches.len(), 1);
        let m = &matches.theorem_matches[0];
        assert_eq!(m.theorem, WEIGHTED_RENEWAL_THEOREM);
        assert_eq!(m.status, CONDITIONAL_THEOREM_MATCH);
        assert_eq!(m.window.max_step, Some(1));
        assert_eq!(m.window.inequality, "2*s <= 2*1+1 (d*s <= d*ell+1)");
        assert_eq!(m.boundary.status, TheoremBoundaryStatus::Canonical);
        assert_eq!(m.boundary.required_rows, 3);
        assert!(m
            .family
            .starts_with("P_n = P_(n-1) + 3 t P_(n-2) + 2 t^2 P_(n-3)"));
        assert_eq!(m.remaining_obligations.len(), 1);
        assert_eq!(m.proof_status, PAPER_PROOF_STATUS);
        assert_eq!(
            matches.rejected_theorems[0].theorem,
            FIXED_DELAY_PREDECESSOR_THEOREM
        );

        // Search output labels the first row 1: relabeled boundary match.
        let m = &match_recurrence_theorems(&rec, 1, &rows).theorem_matches[0];
        assert_eq!(
            m.boundary.status,
            TheoremBoundaryStatus::CanonicalAfterRelabel
        );
        assert_eq!(m.boundary.index_shift, Some(1));
        assert_eq!(m.remaining_obligations.len(), 2);

        let mut wrong = rows.clone();
        wrong[4][1] += q(1);
        let m = &match_recurrence_theorems(&rec, 0, &wrong).theorem_matches[0];
        assert_eq!(m.boundary.status, TheoremBoundaryStatus::Mismatch);
        assert_eq!(m.boundary.first_mismatch_index, Some(4));
        assert!(!m.boundary.rows_match_canonical);
        let m = &match_recurrence_theorems(&rec, 0, &rows[..2]).theorem_matches[0];
        assert_eq!(m.boundary.status, TheoremBoundaryStatus::InsufficientRows);
    }

    #[test]
    fn non_renewal_recurrences_are_rejected_with_reasons() {
        let fibonacci_like = renewal(1, &[1, 1, 1]); // R = 1 + z + z^2 is not real-rooted
        let matches = match_recurrence_theorems(&fibonacci_like, 0, &[]);
        assert!(matches.theorem_matches.is_empty());
        assert_eq!(matches.rejected_theorems.len(), 2);
        assert!(matches.rejected_theorems[0].reason.contains("real-rooted"));
    }

    #[test]
    fn constant_weight_matches_both_theorems_with_their_own_boundaries() {
        // P_n = P_(n-1) + 2 t P_(n-2): renewal (ell = 1, d = 1) and fixed delay L = 2.
        let rec = renewal(1, &[1, 2]);
        let renewal_rows = WeightedRenewal::new(1, vec![q(1), q(2)]).unwrap().rows(6);
        let matches = match_recurrence_theorems(&rec, 0, &renewal_rows);
        assert_eq!(matches.theorem_matches.len(), 2);
        assert_eq!(
            matches.theorem_matches[0].boundary.status,
            TheoremBoundaryStatus::Canonical
        );
        assert_eq!(matches.theorem_matches[0].window.max_step, Some(2));
        let predecessor = &matches.theorem_matches[1];
        assert_eq!(predecessor.theorem, FIXED_DELAY_PREDECESSOR_THEOREM);
        assert_eq!(predecessor.window.max_step, Some(2));
        assert_eq!(predecessor.boundary.status, TheoremBoundaryStatus::Mismatch);
        assert_eq!(predecessor.boundary.first_mismatch_index, Some(1));
    }

    #[test]
    fn index_dependent_weight_relabels_with_shifted_weight() {
        // P_n = P_(n-1) + n t P_(n-2), rows Q_0.. supplied with first label 1.
        let rec = Recurrence {
            denominator: None,
            inhomogeneous: None,
            terms: vec![
                RecurrenceTerm::new(
                    1,
                    0,
                    BivarPoly {
                        coeffs: vec![vec![q(1)]],
                    },
                ),
                RecurrenceTerm::new(
                    2,
                    0,
                    BivarPoly {
                        coeffs: vec![vec![q(0), q(0)], vec![q(0), q(1)]],
                    },
                ),
            ],
        };
        let family = recognize_fixed_delay_predecessor(&rec).unwrap();
        let rows = family.rows(7).unwrap();
        let matches = match_recurrence_theorems(&rec, 0, &rows);
        assert_eq!(matches.theorem_matches.len(), 1);
        let m = &matches.theorem_matches[0];
        assert_eq!(m.boundary.status, TheoremBoundaryStatus::Canonical);
        assert!(m.variable_predecessor.is_some());
        assert!(matches.rejected_theorems[0]
            .reason
            .contains("index-independent"));
        // Exact labels: supplied rows labeled 3.. are Q_3, Q_4, ...
        let m = &match_recurrence_theorems(&rec, 3, &rows[3..]).theorem_matches[0];
        assert_eq!(m.boundary.status, TheoremBoundaryStatus::Canonical);
        // Relabel with w(n+1) = n + 1.
        let shifted = family.with_index_shift(1).unwrap().rows(6).unwrap();
        let m = &match_recurrence_theorems(&rec, 1, &shifted).theorem_matches[0];
        assert_eq!(
            m.boundary.status,
            TheoremBoundaryStatus::CanonicalAfterRelabel
        );
        // A weight negative at some n >= 1 is a rejection, not a match.
        let mut negative = rec;
        negative.terms[1].coeff.coeffs[0][1] = q(-3);
        let matches = match_recurrence_theorems(&negative, 0, &rows);
        assert!(matches.theorem_matches.is_empty());
        assert!(matches.rejected_theorems[1].reason.contains("w(1) = -2"));
    }

    #[test]
    fn supplied_record_report_round_trips() {
        let rec = renewal(2, &[1, 1]);
        let rows = WeightedRenewal::new(2, vec![q(1), q(1)]).unwrap().rows(3);
        let record = RecurrenceJson::from_recurrence_rational(&rec, 0, &rows, None);
        let report = RecurrenceTheoremReport::for_supplied(&record).unwrap();
        assert_eq!(report.schema, RECURRENCE_THEOREM_SCHEMA);
        assert_eq!(report.finite_evidence.status, "finite_evidence_only");
        assert_eq!(report.theorem_matches.len(), 2);
        assert_eq!(report.theorem_matches[0].window.max_step, Some(3));
        let encoded = serde_json::to_string(&report).unwrap();
        let decoded: RecurrenceTheoremReport = serde_json::from_str(&encoded).unwrap();
        assert_eq!(decoded, report);
        assert!(report.theorem_matches[0]
            .summary()
            .contains("not established"));
    }
}
