//! Implicit differential relations of polynomial sequences, generating
//! descriptions with reconstruction data, uniqueness arguments and
//! reconstruction certificates.
//!
//! # Relations
//!
//! An [`ImplicitRelation`] is a homogeneous linear relation
//!
//! ```text
//! sum_terms c(n,x) D^d P_(n+i) = 0        (base index n)
//! ```
//!
//! with derivatives allowed on every row, the latest included. A relation
//! whose terms all have shift zero is a row differential equation; one whose
//! leading row (largest shift) carries no derivative is a successor rule in
//! the sense of [`Recurrence`]. The same type, restricted to shift zero,
//! represents a differential operator `T = sum t_d(n,x) D^d` acting on one row
//! (for example the map of a reconstruction certificate).
//!
//! A [`CoefficientCondition`] is a linear condition
//! `sum_j c_j(n) [x^(k_j)] P_(n+i_j) = 0` on finitely many coefficients, valid
//! for base indices `n >= first_base`. Exceptional-index conditions produced
//! by a clearing multiplier and coefficient boundary conditions (such as
//! `P_(n+1)(0) = P_n(0)`) share this form but are recorded separately.
//!
//! # Generating descriptions
//!
//! A [`GeneratingDescription`] stores seeds, relations with their index
//! domains and evidence, exceptional-index conditions and boundary
//! conditions. [`assess_generating_description`] keeps two questions apart:
//!
//! - annihilation: does the sequence satisfy every component? Each component
//!   carries its evidence (`transported_identity`, `certified_by_reconstruction`,
//!   `checked_supplied_rows`, `finite_fit`, ...);
//! - uniqueness: do the components determine the sequence from the seeds? This
//!   is a property of the constraint system alone ([`UniquenessArgument`]).
//!
//! The uniqueness argument picks a relation of positive leader shift `R` and
//! compares coefficients of the difference of two solutions. With the leader
//! operator `sum_d c_(R,d)(n,x) D^d`, every monomial `n^u x^j D^d` moves the
//! degree by `j - d`. In the degree grading (largest `j - d = v`), the top
//! coefficient of a nonzero difference `y` of degree `k` is multiplied by
//! `psi(n,k) = sum_(j-d=v) c_(u,j) n^u k(k-1)...(k-d+1)`; if `psi(n,k) != 0`
//! for all `k >= 0` and every base index used, the leader is injective. In the
//! valuation grading (smallest `j - d`) the lowest coefficient is multiplied by
//! the analogous `phi(n,k)`, and every exceptional index `k` with
//! `phi(n,k) = 0` must be pinned by a condition whose terms on the new row
//! involve only `[x^k]` with a nonvanishing coefficient and lower indices. The
//! exceptional set must be certified uniformly in `n`: an index-free
//! polynomial (exact roots), or a one-signed shifted polynomial (no zeros).
//! Anything else stays uncertified; prefixes are never tested.
//!
//! # Reconstruction certificates
//!
//! A [`ReconstructionCertificate`] records an operator identity
//! `L o T = sum_nu H_nu o S^nu o M`, where `L` is a relation of the source
//! sequence `P`, `T` the map with `P_n = T Q_n`, `M` a relation of `Q`, and
//! each `H_nu = sum_d h_(nu,d)(n,x) D^d` records its differentiations.
//! [`ReconstructionCertificate::verify_identity`] checks the identity in the
//! Ore algebra; [`ReconstructionCertificate::replay_on_rows`] first derives the
//! differentiated identities `D^d (M Q)_(n+nu)` and then forms the linear
//! combination, comparing it with `L (T Q)` at every base index. Source
//! uniqueness, initial agreement and exceptional indices complete an
//! equivalence proof; see `implicit_transform_recurrence`.
//!
//! Related modules: `diagonal` (diagonal transport and derivative
//! antecedents), `transform` (the coordinate-change registry).

use super::implicit_algebra::{
    analyze_zeros, biv_text, falling_factorial, poly_text, roots_in_domain, Biv, DOp,
};
use super::*;
use crate::vec_poly;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

type Q = BigRational;

/// Evidence: the defining rule of the source sequence.
pub const DEFINING_EVIDENCE: &str = "defining";

/// Evidence: a fitted relation proved for every index by a reconstruction
/// certificate, source uniqueness and initial agreement.
pub const CERTIFIED_BY_RECONSTRUCTION_EVIDENCE: &str = "certified_by_reconstruction";

/// Evidence: an exact consequence of a supplied source relation that is
/// itself verified only on the supplied rows.
pub const TRANSPORTED_SUPPLIED_EVIDENCE: &str = "transported_supplied_relation";

/// Description status: certified annihilation and certified uniqueness.
pub const EQUIVALENT_DESCRIPTION: &str = "equivalent_description";

/// Description status: certified annihilation without a uniqueness certificate.
pub const CERTIFIED_CONSEQUENCE: &str = "certified_consequence";

/// Description status: some component is supported by finite evidence only.
pub const FINITE_EVIDENCE: &str = "finite_evidence";

/// Description status: the seeds violate a component.
pub const INCONSISTENT_SEEDS: &str = "inconsistent_seeds";

/// Relation role: the leading row carries no derivative.
pub const SUCCESSOR_ROLE: &str = "successor";

/// Relation role: the leading row carries derivatives (an implicit rule).
pub const IMPLICIT_SUCCESSOR_ROLE: &str = "implicit_successor";

/// Relation role: every term has shift zero (a row differential equation).
pub const ROW_EQUATION_ROLE: &str = "row_equation";

/// Whether an evidence label proves annihilation for every index of its domain.
pub fn evidence_is_certified(evidence: &str) -> bool {
    [
        DEFINING_EVIDENCE,
        TRANSPORTED_IDENTITY_EVIDENCE,
        CERTIFIED_BY_RECONSTRUCTION_EVIDENCE,
    ]
    .contains(&evidence)
}

// ---------------------------------------------------------------------------
// JSON helpers
// ---------------------------------------------------------------------------

pub(super) fn biv_json(p: &Biv) -> BivarPolyJson {
    BivarPolyJson {
        coeffs: p
            .to_bivar()
            .normalized_coefficients()
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect(),
    }
}

pub(super) fn biv_from_json(p: &BivarPolyJson) -> Result<Biv, String> {
    p.coeffs
        .iter()
        .map(|row| row.iter().map(|c| parse_rational_coeff(c)).collect())
        .collect::<Result<Vec<Vec<Q>>, String>>()
        .map(Biv::new)
}

pub(super) fn encode_poly(p: &[Q]) -> Vec<String> {
    vec_poly::trim(p)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

pub(super) fn decode_poly(p: &[String]) -> Result<Vec<Q>, String> {
    p.iter()
        .map(|c| parse_rational_coeff(c))
        .collect::<Result<Vec<_>, _>>()
        .map(|r| vec_poly::trim(&r))
}

// ---------------------------------------------------------------------------
// Relations
// ---------------------------------------------------------------------------

/// One term `coeff(n,x) D^deriv_order P_(n+shift)` of an [`ImplicitRelation`].
#[derive(Debug, Clone)]
pub struct ImplicitTerm {
    /// Forward shift `i` of the row `P_(n+i)`.
    pub shift: usize,
    /// Derivative order `d` in `x`.
    pub deriv_order: usize,
    /// Coefficient `c(n,x)`.
    pub coeff: BivarPoly,
}

/// A homogeneous relation `sum c(n,x) D^d P_(n+i) = 0` at base index `n`.
///
/// Terms with equal shift and derivative order may repeat; they are summed.
#[derive(Debug, Clone, Default)]
pub struct ImplicitRelation {
    /// The terms.
    pub terms: Vec<ImplicitTerm>,
}

/// JSON form of one [`ImplicitTerm`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitTermJson {
    /// Forward shift `i` of the row `P_(n+i)`.
    pub shift: usize,
    /// Derivative order `d`.
    pub deriv_order: usize,
    /// Coefficient `c(n,x)`: `coeffs[a][b]` is the coefficient of `n^a x^b`.
    pub coeff: BivarPolyJson,
}

/// JSON form of an [`ImplicitRelation`]: `sum c(n,x) D^d P_(n+shift) = 0`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitRelationJson {
    /// Terms, in canonical order (shift, then derivative order) on output.
    pub terms: Vec<ImplicitTermJson>,
}

impl ImplicitRelation {
    pub(super) fn op(&self) -> DOp {
        let mut op = DOp::default();
        for term in &self.terms {
            op.add_term((term.shift, term.deriv_order), Biv::from_bivar(&term.coeff));
        }
        op
    }

    pub(super) fn from_op(op: &DOp) -> Self {
        Self {
            terms: op
                .0
                .iter()
                .map(|((shift, deriv_order), c)| ImplicitTerm {
                    shift: *shift,
                    deriv_order: *deriv_order,
                    coeff: c.to_bivar(),
                })
                .collect(),
        }
    }

    /// Convert a scalar rule `f(m,x) P_m = sum c(m,x) D^d P_(m-r)` (target index
    /// `m`) to the forward form with base `n = m - R`, `R` the largest offset.
    /// Returns the relation and `R`.
    ///
    /// Rules with `(-1)^n` factors or an additive term are rejected: implicit
    /// relations are homogeneous and sign free.
    pub fn from_recurrence(rule: &Recurrence) -> Result<(Self, usize), String> {
        if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
            return Err("implicit relations carry no (-1)^n factors".into());
        }
        if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
            return Err("implicit relations are homogeneous; the rule has an additive term".into());
        }
        let r = rule.max_offset();
        let shift = i64::try_from(r).map_err(|_| "offset overflow")?;
        let mut op = DOp::default();
        let leader = rule
            .denominator
            .as_ref()
            .map(Biv::from_bivar)
            .unwrap_or_else(Biv::one);
        op.add_term((r, 0), leader.shift_n(shift));
        for term in &rule.terms {
            op.add_term(
                (r - term.offset, term.deriv_order),
                Biv::from_bivar(&term.coeff)
                    .shift_n(shift)
                    .scale(&-Q::one()),
            );
        }
        Ok((Self::from_op(&op), r))
    }

    /// Encode as JSON (canonical order, terms merged).
    pub fn to_json(&self) -> ImplicitRelationJson {
        ImplicitRelationJson {
            terms: self
                .op()
                .0
                .iter()
                .map(|((shift, deriv_order), c)| ImplicitTermJson {
                    shift: *shift,
                    deriv_order: *deriv_order,
                    coeff: biv_json(c),
                })
                .collect(),
        }
    }

    /// Decode from JSON.
    pub fn from_json(json: &ImplicitRelationJson) -> Result<Self, String> {
        let mut op = DOp::default();
        for term in &json.terms {
            op.add_term((term.shift, term.deriv_order), biv_from_json(&term.coeff)?);
        }
        Ok(Self::from_op(&op))
    }

    /// Largest shift (0 for a row equation or the zero relation).
    pub fn max_shift(&self) -> usize {
        self.op().max_shift()
    }

    /// Largest derivative order.
    pub fn max_derivative(&self) -> usize {
        self.op().max_derivative()
    }

    /// Whether every coefficient vanishes.
    pub fn is_zero(&self) -> bool {
        self.op().is_zero()
    }

    /// `successor`, `implicit_successor` or `row_equation`.
    pub fn role(&self) -> &'static str {
        let op = self.op();
        let r = op.max_shift();
        if r == 0 {
            ROW_EQUATION_ROLE
        } else if op.column(r).max_derivative() == 0 {
            SUCCESSOR_ROLE
        } else {
            IMPLICIT_SUCCESSOR_ROLE
        }
    }

    /// The canonical form: common powers of `x` removed and the leading
    /// coefficient (largest shift and derivative, then highest `n`, `x`
    /// powers) scaled to one. Two relations are the same up to a nonzero
    /// factor `c x^e` exactly when their canonical forms agree.
    pub fn canonical(&self) -> Self {
        Self::from_op(&self.op().normalize().2)
    }

    /// Whether two relations agree up to a factor `c x^e`, `c != 0`.
    pub fn same_relation(&self, other: &Self) -> bool {
        self.op().normalize().2 == other.op().normalize().2
    }

    /// `(L P)_base` when `rows[j]` is `P_(origin + j)` and every referenced row
    /// is present.
    pub fn residual(&self, rows: &[Vec<Q>], origin: usize, base: usize) -> Option<Vec<Q>> {
        self.op().apply(rows, origin, base)
    }

    /// Whether the relation holds exactly at `base`; `None` if a row is missing.
    pub fn holds_at(&self, rows: &[Vec<Q>], origin: usize, base: usize) -> Option<bool> {
        self.residual(rows, origin, base)
            .map(|r| vec_poly::is_zero(&r))
    }

    /// Apply a shift-free relation (an operator on one row) to the row at index `n`.
    pub fn apply_to_row(&self, row: &[Q], n: usize) -> Vec<Q> {
        self.op().apply_row(row, n)
    }

    /// Human-readable form `sum (c) D^d Q[n+i] = 0` with `symbol` as the row name.
    pub fn text(&self, symbol: &str) -> String {
        op_text(&self.op(), symbol)
    }
}

/// `""` for 1, `"-"` for -1, `"m "` for a monomial `m`, `"(p) "` otherwise.
fn coefficient_prefix(text: &str, monomial: bool) -> String {
    match text {
        "1" => String::new(),
        "-1" => "-".into(),
        _ if monomial => format!("{text} "),
        _ => format!("({text}) "),
    }
}

/// Join signed parts with ` + ` / ` - `.
fn join_signed(parts: &[String]) -> String {
    let mut text = parts[0].clone();
    for part in &parts[1..] {
        match part.strip_prefix('-') {
            Some(rest) => {
                text.push_str(" - ");
                text.push_str(rest);
            }
            None => {
                text.push_str(" + ");
                text.push_str(part);
            }
        }
    }
    text
}

pub(super) fn op_text(op: &DOp, symbol: &str) -> String {
    if op.is_zero() {
        return "0 = 0".into();
    }
    let parts =
        op.0.iter()
            .rev()
            .map(|((i, d), c)| {
                let row = if *i == 0 {
                    format!("{symbol}[n]")
                } else {
                    format!("{symbol}[n+{i}]")
                };
                let derivative = match d {
                    0 => String::new(),
                    1 => "D ".into(),
                    _ => format!("D^{d} "),
                };
                let monomial = c.rows().iter().flatten().filter(|v| !v.is_zero()).count() == 1;
                let factor = coefficient_prefix(&biv_text(c, "x"), monomial);
                format!("{factor}{derivative}{row}")
            })
            .collect::<Vec<_>>();
    format!("{} = 0", join_signed(&parts))
}

// ---------------------------------------------------------------------------
// Coefficient conditions
// ---------------------------------------------------------------------------

/// One term `coeff(n) [x^coefficient_index] P_(n+shift)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConditionTerm {
    /// Forward shift of the row.
    pub shift: usize,
    /// The coefficient index `k`.
    pub coefficient_index: usize,
    /// Polynomial in `n`, ascending.
    pub coeff: Vec<Q>,
}

/// A linear condition `sum_j c_j(n) [x^(k_j)] P_(n+i_j) = 0` at base index `n`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct CoefficientCondition {
    /// The terms; equal `(shift, coefficient_index)` pairs are summed.
    pub terms: Vec<ConditionTerm>,
}

/// JSON form of one [`ConditionTerm`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ConditionTermJson {
    /// Forward shift of the row.
    pub shift: usize,
    /// The coefficient index `k` (power of `x`).
    pub coefficient_index: usize,
    /// Coefficient polynomial in `n`, ascending rational strings.
    pub coeff: Vec<String>,
}

/// JSON form of a [`CoefficientCondition`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct CoefficientConditionJson {
    /// Terms in canonical order (shift, then coefficient index) on output.
    pub terms: Vec<ConditionTermJson>,
}

impl CoefficientCondition {
    /// Merge equal `(shift, coefficient_index)` terms and drop zero ones.
    pub fn canonical(&self) -> Self {
        let mut merged: BTreeMap<(usize, usize), Vec<Q>> = BTreeMap::new();
        for term in &self.terms {
            let entry = merged
                .entry((term.shift, term.coefficient_index))
                .or_insert_with(|| vec![Q::zero()]);
            *entry = vec_poly::add(entry, &term.coeff);
        }
        Self {
            terms: merged
                .into_iter()
                .filter(|(_, c)| !vec_poly::is_zero(c))
                .map(|((shift, coefficient_index), coeff)| ConditionTerm {
                    shift,
                    coefficient_index,
                    coeff: vec_poly::trim(&coeff),
                })
                .collect(),
        }
    }

    /// Whether the canonical condition has no terms (it holds vacuously).
    pub fn is_vacuous(&self) -> bool {
        self.canonical().terms.is_empty()
    }

    /// Largest shift.
    pub fn max_shift(&self) -> usize {
        self.terms.iter().map(|t| t.shift).max().unwrap_or(0)
    }

    /// The value of the left side at `base`; `None` if a row is missing.
    pub fn residual(&self, rows: &[Vec<Q>], origin: usize, base: usize) -> Option<Q> {
        let n = Q::from_integer(BigInt::from(base));
        let mut sum = Q::zero();
        for term in &self.terms {
            let row = rows.get(base.checked_add(term.shift)?.checked_sub(origin)?)?;
            let coefficient = row
                .get(term.coefficient_index)
                .cloned()
                .unwrap_or_else(Q::zero);
            sum += vec_poly::evaluate(&term.coeff, &n) * coefficient;
        }
        Some(sum)
    }

    /// Whether the condition holds exactly at `base`; `None` if a row is missing.
    pub fn holds_at(&self, rows: &[Vec<Q>], origin: usize, base: usize) -> Option<bool> {
        self.residual(rows, origin, base).map(|r| r.is_zero())
    }

    /// Encode as JSON (canonical form).
    pub fn to_json(&self) -> CoefficientConditionJson {
        CoefficientConditionJson {
            terms: self
                .canonical()
                .terms
                .iter()
                .map(|t| ConditionTermJson {
                    shift: t.shift,
                    coefficient_index: t.coefficient_index,
                    coeff: encode_poly(&t.coeff),
                })
                .collect(),
        }
    }

    /// Decode from JSON.
    pub fn from_json(json: &CoefficientConditionJson) -> Result<Self, String> {
        Ok(Self {
            terms: json
                .terms
                .iter()
                .map(|t| {
                    Ok(ConditionTerm {
                        shift: t.shift,
                        coefficient_index: t.coefficient_index,
                        coeff: decode_poly(&t.coeff)?,
                    })
                })
                .collect::<Result<Vec<_>, String>>()?,
        }
        .canonical())
    }

    /// Human-readable form `sum (c(n)) [x^k] Q[n+i] = 0`.
    pub fn text(&self, symbol: &str) -> String {
        let canonical = self.canonical();
        if canonical.terms.is_empty() {
            return "0 = 0".into();
        }
        let parts = canonical
            .terms
            .iter()
            .rev()
            .map(|t| {
                let row = if t.shift == 0 {
                    format!("{symbol}[n]")
                } else {
                    format!("{symbol}[n+{}]", t.shift)
                };
                let monomial = t.coeff.iter().filter(|v| !v.is_zero()).count() == 1;
                format!(
                    "{}[x^{}] {row}",
                    coefficient_prefix(&poly_text(&t.coeff, "n"), monomial),
                    t.coefficient_index
                )
            })
            .collect::<Vec<_>>();
        format!("{} = 0", join_signed(&parts))
    }
}

// ---------------------------------------------------------------------------
// Generating descriptions
// ---------------------------------------------------------------------------

/// A relation of a description with its domain and evidence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DescribedRelation {
    /// The relation.
    pub relation: ImplicitRelationJson,
    /// It holds for base indices `n >= first_base`.
    pub first_base: usize,
    /// `successor`, `implicit_successor` or `row_equation`.
    pub role: String,
    /// Evidence label (see [`evidence_is_certified`]).
    pub evidence: String,
    /// Human-readable form.
    pub text: String,
    /// Provenance or further explanation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

impl DescribedRelation {
    /// Describe a relation; role and text are derived from it.
    pub fn new(
        relation: &ImplicitRelation,
        first_base: usize,
        evidence: &str,
        symbol: &str,
    ) -> Self {
        Self {
            relation: relation.to_json(),
            first_base,
            role: relation.role().to_string(),
            evidence: evidence.to_string(),
            text: relation.text(symbol),
            note: None,
        }
    }

    /// Return this record with a note.
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.note = Some(note.into());
        self
    }
}

/// A coefficient condition of a description with its domain and evidence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DescribedCondition {
    /// The condition.
    pub condition: CoefficientConditionJson,
    /// It holds for base indices `n >= first_base`.
    pub first_base: usize,
    /// Evidence label.
    pub evidence: String,
    /// The exceptional coefficient index it pins, for exceptional-index conditions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub coefficient_index: Option<usize>,
    /// Human-readable form.
    pub text: String,
    /// Provenance or further explanation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

impl DescribedCondition {
    /// Describe a condition; the text is derived from it.
    pub fn new(
        condition: &CoefficientCondition,
        first_base: usize,
        evidence: &str,
        coefficient_index: Option<usize>,
        symbol: &str,
    ) -> Self {
        Self {
            condition: condition.to_json(),
            first_base,
            evidence: evidence.to_string(),
            coefficient_index,
            text: condition.text(symbol),
            note: None,
        }
    }

    /// Return this record with a note.
    pub fn with_note(mut self, note: impl Into<String>) -> Self {
        self.note = Some(note.into());
        self
    }
}

/// Seeds, relations and coefficient conditions describing a sequence
/// `Q_first_index, Q_(first_index+1), ...`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct GeneratingDescription {
    /// Index of `seeds[0]`.
    pub first_index: usize,
    /// Seed rows, ascending rational strings.
    pub seeds: Vec<Vec<String>>,
    /// Relations (successor rules, implicit rules, row equations).
    pub relations: Vec<DescribedRelation>,
    /// Conditions at exceptional coefficient indices of a clearing multiplier.
    pub exceptional_index_conditions: Vec<DescribedCondition>,
    /// Coefficient boundary conditions such as `Q_(n+1)(0) = Q_n(0)`.
    pub boundary_conditions: Vec<DescribedCondition>,
}

/// Why the components of a description determine the sequence, or why that
/// is not certified. Independent of the evidence for annihilation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct UniquenessArgument {
    /// `certified` or `not_certified`.
    pub status: String,
    /// Index into `relations` of the determining relation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub relation: Option<usize>,
    /// Leader shift `R` of the determining relation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub leader_shift: Option<usize>,
    /// `degree` or `valuation`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub grading: Option<String>,
    /// The indicial polynomial `phi(n,k)`: `coeffs[a][b]` multiplies `n^a k^b`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub indicial_polynomial: Option<BivarPolyJson>,
    /// Certified exceptional coefficient indices of the grading.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub exceptional_indices: Vec<usize>,
    /// The certificate of the exceptional set.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exceptional_certificate: Option<String>,
    /// Conditions pinning the exceptional indices, as `list[i]` references.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub covering_conditions: Vec<String>,
    /// The argument, when certified.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub statement: Option<String>,
    /// Every failed attempt, with its reason.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reasons: Vec<String>,
}

/// Result of [`assess_generating_description`].
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct DescriptionAssessment {
    /// `equivalent_description`, `certified_consequence`, `finite_evidence`
    /// or `inconsistent_seeds`.
    pub status: String,
    /// Whether every relation and condition instance inside the seeds holds.
    pub seeds_consistent: bool,
    /// Components whose evidence is not certified, as `list[i]: evidence`.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub uncertified_components: Vec<String>,
    /// The uniqueness argument.
    pub uniqueness: UniquenessArgument,
    /// Explanations (seed violations, missing certificates).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub reasons: Vec<String>,
}

struct ParsedCondition {
    label: String,
    condition: CoefficientCondition,
    first_base: usize,
}

/// Assess a generating description: seed consistency, annihilation evidence
/// and a uniqueness argument, kept separately.
///
/// The status is `equivalent_description` only when every component carries
/// certified evidence and uniqueness is certified; removing a condition that
/// the argument needs (for example a lost constant term) removes the claim.
pub fn assess_generating_description(
    description: &GeneratingDescription,
) -> Result<DescriptionAssessment, String> {
    let seeds = description
        .seeds
        .iter()
        .map(|row| decode_poly(row))
        .collect::<Result<Vec<_>, _>>()?;
    if seeds.is_empty() {
        return Err("a generating description needs at least one seed row".into());
    }
    let first = description.first_index;
    let relations = description
        .relations
        .iter()
        .map(|r| Ok((ImplicitRelation::from_json(&r.relation)?, r.first_base)))
        .collect::<Result<Vec<_>, String>>()?;
    let mut conditions = Vec::new();
    for (list, records) in [
        (
            "exceptional_index_conditions",
            &description.exceptional_index_conditions,
        ),
        ("boundary_conditions", &description.boundary_conditions),
    ] {
        for (i, record) in records.iter().enumerate() {
            conditions.push(ParsedCondition {
                label: format!("{list}[{i}]"),
                condition: CoefficientCondition::from_json(&record.condition)?,
                first_base: record.first_base,
            });
        }
    }

    let mut reasons = Vec::new();
    let seed_end = first + seeds.len();
    for (i, (relation, first_base)) in relations.iter().enumerate() {
        let r = relation.max_shift();
        for base in (*first_base).max(first)..seed_end.saturating_sub(r) {
            if relation.holds_at(&seeds, first, base) == Some(false) {
                reasons.push(format!(
                    "relations[{i}] fails on the seeds at base n = {base}"
                ));
            }
        }
    }
    for c in &conditions {
        let r = c.condition.max_shift();
        for base in c.first_base.max(first)..seed_end.saturating_sub(r) {
            if c.condition.holds_at(&seeds, first, base) == Some(false) {
                reasons.push(format!("{} fails on the seeds at base n = {base}", c.label));
            }
        }
    }
    let seeds_consistent = reasons.is_empty();

    let mut uncertified = Vec::new();
    for (i, r) in description.relations.iter().enumerate() {
        if !evidence_is_certified(&r.evidence) {
            uncertified.push(format!("relations[{i}]: {}", r.evidence));
        }
    }
    for (list, records) in [
        (
            "exceptional_index_conditions",
            &description.exceptional_index_conditions,
        ),
        ("boundary_conditions", &description.boundary_conditions),
    ] {
        for (i, record) in records.iter().enumerate() {
            if !evidence_is_certified(&record.evidence) {
                uncertified.push(format!("{list}[{i}]: {}", record.evidence));
            }
        }
    }

    let uniqueness = uniqueness_argument(first, seeds.len(), &relations, &conditions);
    let unique = uniqueness.status == "certified";
    let status = if !seeds_consistent {
        INCONSISTENT_SEEDS
    } else if !uncertified.is_empty() {
        FINITE_EVIDENCE
    } else if unique {
        EQUIVALENT_DESCRIPTION
    } else {
        CERTIFIED_CONSEQUENCE
    };
    if !uncertified.is_empty() {
        reasons.push(
            "some components are supported by finite evidence only; no equivalence is claimed"
                .into(),
        );
    }
    if !unique {
        reasons.push("uniqueness is not certified; the components are consequences only".into());
    }
    Ok(DescriptionAssessment {
        status: status.to_string(),
        seeds_consistent,
        uncertified_components: uncertified,
        uniqueness,
        reasons,
    })
}

/// The indicial polynomial of a shift-free leader in a grading:
/// `(v, sum_(j-d=v) c_(u,j) n^u k^(d falling))`.
fn indicial(leader: &DOp, top: bool) -> Option<(i64, Biv)> {
    let mut shifts = Vec::new();
    for ((_, d), c) in &leader.0 {
        for row in c.rows() {
            for (j, value) in row.iter().enumerate() {
                if !value.is_zero() {
                    shifts.push(j as i64 - *d as i64);
                }
            }
        }
    }
    let v = if top {
        shifts.into_iter().max()?
    } else {
        shifts.into_iter().min()?
    };
    let mut phi = Biv::zero();
    for ((_, d), c) in &leader.0 {
        let falling = falling_factorial(*d);
        for (u, row) in c.rows().iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if value.is_zero() || j as i64 - *d as i64 != v {
                    continue;
                }
                let mut rows = vec![Vec::new(); u + 1];
                rows[u] = vec_poly::scale(&falling, value);
                phi = phi.add(&Biv::new(rows));
            }
        }
    }
    Some((v, phi))
}

fn uniqueness_argument(
    first: usize,
    seed_count: usize,
    relations: &[(ImplicitRelation, usize)],
    conditions: &[ParsedCondition],
) -> UniquenessArgument {
    let seed_end = first + seed_count;
    let mut reasons = Vec::new();
    for (index, (relation, first_base)) in relations.iter().enumerate() {
        let op = relation.op();
        let r = op.max_shift();
        if r == 0 {
            reasons.push(format!(
                "relations[{index}] is a row equation and does not determine later rows"
            ));
            continue;
        }
        if *first_base < first {
            reasons.push(format!(
                "relations[{index}] starts at base {first_base}, before the first seed {first}"
            ));
            continue;
        }
        if first_base + r > seed_end {
            reasons.push(format!(
                "relations[{index}] determines rows from index {} on, but the seeds end at {}",
                first_base + r,
                seed_end - 1
            ));
            continue;
        }
        let start = (seed_end - r).max(*first_base);
        let leader = op.column(r);
        for top in [true, false] {
            let grading = if top { "degree" } else { "valuation" };
            let Some((_, phi)) = indicial(&leader, top) else {
                continue;
            };
            let analysis = analyze_zeros(&phi, start, None);
            if !analysis.singular_indices.is_empty() {
                reasons.push(format!(
                    "relations[{index}]: the {grading} part of the leader vanishes identically \
                     at base n = {:?}, so row n + {r} is not determined (singular leader)",
                    analysis.singular_indices
                ));
                continue;
            }
            if !analysis.certified {
                reasons.push(format!(
                    "relations[{index}]: {grading} grading: {}",
                    analysis.certificate
                ));
                continue;
            }
            let mut argument = UniquenessArgument {
                status: "certified".into(),
                relation: Some(index),
                leader_shift: Some(r),
                grading: Some(grading.into()),
                indicial_polynomial: Some(biv_json(&phi)),
                exceptional_indices: analysis.uniform_indices.clone(),
                exceptional_certificate: Some(analysis.certificate.clone()),
                covering_conditions: Vec::new(),
                statement: None,
                reasons: Vec::new(),
            };
            if top {
                if !analysis.uniform_indices.is_empty() {
                    reasons.push(format!(
                        "relations[{index}]: degree grading has exceptional indices {:?}",
                        analysis.uniform_indices
                    ));
                    continue;
                }
                argument.statement = Some(format!(
                    "relations[{index}] determines row n + {r} for every base n >= {start}: the \
                     top coefficient of a nonzero difference of degree k is multiplied by \
                     psi(n,k) = {}, which has no zero with n >= {start}, k >= 0 ({}). With \
                     the seeds at indices {first}..{} every row is determined.",
                    biv_text(&phi, "k"),
                    analysis.certificate,
                    seed_end - 1
                ));
                return argument;
            }
            let mut missing = Vec::new();
            for &k in &analysis.uniform_indices {
                match conditions
                    .iter()
                    .find(|c| condition_covers(&c.condition, c.first_base, k, r, start, first))
                {
                    Some(c) => argument
                        .covering_conditions
                        .push(format!("{} pins [x^{k}] of row n + {r}", c.label)),
                    None => missing.push(k),
                }
            }
            if !missing.is_empty() {
                reasons.push(format!(
                    "relations[{index}]: valuation grading has exceptional indices {:?}; no \
                     condition pins [x^k] of the new row for k in {missing:?}",
                    analysis.uniform_indices
                ));
                continue;
            }
            argument.statement = Some(format!(
                "relations[{index}] determines row n + {r} for every base n >= {start}: \
                 comparing coefficients of a difference in ascending order, the lowest one \
                 at index k is multiplied by phi(n,k) = {}, nonzero except at the certified \
                 exceptional indices {:?} ({}), which the listed conditions pin. With the \
                 seeds at indices {first}..{} every row is determined.",
                biv_text(&phi, "k"),
                analysis.uniform_indices,
                analysis.certificate,
                seed_end - 1
            ));
            return argument;
        }
    }
    if relations.is_empty() {
        reasons.push("the description has no relations".into());
    }
    UniquenessArgument {
        status: "not_certified".into(),
        relation: None,
        leader_shift: None,
        grading: None,
        indicial_polynomial: None,
        exceptional_indices: Vec::new(),
        exceptional_certificate: None,
        covering_conditions: Vec::new(),
        statement: None,
        reasons,
    }
}

/// Whether `condition`, aligned so that its largest shift reaches the new row
/// `n + r` (`n >= start`), pins `[x^k]` of that row.
fn condition_covers(
    condition: &CoefficientCondition,
    first_base: usize,
    k: usize,
    r: usize,
    start: usize,
    first: usize,
) -> bool {
    let canonical = condition.canonical();
    let rc = canonical.max_shift();
    if canonical.terms.is_empty() || start + r < rc {
        return false;
    }
    let aligned = start + r - rc;
    if aligned < first_base || aligned < first {
        return false;
    }
    let top = canonical
        .terms
        .iter()
        .filter(|t| t.shift == rc)
        .collect::<Vec<_>>();
    if top.iter().any(|t| t.coefficient_index > k) {
        return false;
    }
    top.iter()
        .find(|t| t.coefficient_index == k)
        .is_some_and(|t| {
            !vec_poly::is_zero(&t.coeff) && roots_in_domain(&t.coeff, aligned, None).is_empty()
        })
}

// ---------------------------------------------------------------------------
// Reconstruction certificates
// ---------------------------------------------------------------------------

/// One term `coeff(n,x) D^deriv_order` of a shift-free operator.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct OperatorTermJson {
    /// Derivative order `d`.
    pub deriv_order: usize,
    /// Coefficient `c(n,x)`: `coeffs[a][b]` multiplies `n^a x^b`.
    pub coeff: BivarPolyJson,
}

/// The multiplier `H_nu = sum_d h_d(n,x) D^d` of `S^nu M`: the identity
/// `(M Q)_(n+nu) = 0` is differentiated `d` times and multiplied by `h_d`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ReconstructionMultiplier {
    /// The shift `nu`.
    pub shift: usize,
    /// The recorded differentiations `h_d(n,x) D^d`.
    pub differentiations: Vec<OperatorTermJson>,
}

/// An exact operator identity `L o T = sum_nu H_nu o S^nu o M`.
///
/// `L` is a relation of the source rows `P`, `T` the shift-free map with
/// `P_n = T Q_n`, and `M` a relation of `Q`. If `Q` satisfies `M`, then
/// `T Q` satisfies `L`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ReconstructionCertificate {
    /// `L`, in source coordinates.
    pub source_relation: ImplicitRelationJson,
    /// `T`, with `P_n = T Q_n`.
    pub map: Vec<OperatorTermJson>,
    /// `M`, in the new coordinates.
    pub target_relation: ImplicitRelationJson,
    /// The multipliers `H_nu`.
    pub multipliers: Vec<ReconstructionMultiplier>,
    /// What the identity proves.
    pub statement: String,
}

/// `(L, T, M, [(nu, H_nu)])` of a certificate.
type CertificateParts = (DOp, DOp, DOp, Vec<(usize, DOp)>);

fn shift_free_from_json(terms: &[OperatorTermJson]) -> Result<DOp, String> {
    let mut op = DOp::default();
    for term in terms {
        op.add_term((0, term.deriv_order), biv_from_json(&term.coeff)?);
    }
    Ok(op)
}

fn shift_free_to_json(op: &DOp) -> Vec<OperatorTermJson> {
    op.0.iter()
        .map(|((_, d), c)| OperatorTermJson {
            deriv_order: *d,
            coeff: biv_json(c),
        })
        .collect()
}

impl ReconstructionCertificate {
    fn parts(&self) -> Result<CertificateParts, String> {
        let l = ImplicitRelation::from_json(&self.source_relation)?.op();
        let t = shift_free_from_json(&self.map)?;
        let m = ImplicitRelation::from_json(&self.target_relation)?.op();
        let h = self
            .multipliers
            .iter()
            .map(|mu| Ok((mu.shift, shift_free_from_json(&mu.differentiations)?)))
            .collect::<Result<Vec<_>, String>>()?;
        Ok((l, t, m, h))
    }

    /// Check `L o T = sum_nu H_nu o S^nu o M` exactly in the Ore algebra
    /// (`S n = (n + 1) S`, `D x = x D + 1`).
    pub fn verify_identity(&self) -> Result<(), String> {
        let (l, t, m, h) = self.parts()?;
        let left = l.compose(&t);
        let right = h.iter().fold(DOp::default(), |acc, (nu, op)| {
            acc.add(&op.compose(&m.shifted(*nu)))
        });
        if left == right {
            Ok(())
        } else {
            Err("reconstruction identity L o T = sum H_nu S^nu M fails".into())
        }
    }

    /// Replay on rows `rows[j] = Q_(origin + j)`: derive the identities
    /// `(M Q)_(n+nu)`, differentiate them as recorded, form the linear
    /// combination and compare with `L (T Q)` at every base index `n` whose
    /// rows are available. The comparison holds for any rows when the
    /// identity is valid. Returns the number of base indices replayed.
    pub fn replay_on_rows(&self, rows: &[Vec<Q>], origin: usize) -> Result<usize, String> {
        let (l, t, m, h) = self.parts()?;
        let mapped = rows
            .iter()
            .enumerate()
            .map(|(j, row)| t.apply_row(row, origin + j))
            .collect::<Vec<_>>();
        let mut checked = 0;
        for base in origin..origin + rows.len() {
            let Some(left) = l.apply(&mapped, origin, base) else {
                continue;
            };
            let mut right = vec![Q::zero()];
            let mut available = true;
            for (nu, multiplier) in &h {
                let Some(identity) = m.apply(rows, origin, base + nu) else {
                    available = false;
                    break;
                };
                let n = Q::from_integer(BigInt::from(base));
                for ((_, d), coefficient) in &multiplier.0 {
                    let mut derived = identity.clone();
                    for _ in 0..*d {
                        derived = vec_poly::derivative(&derived);
                    }
                    right =
                        vec_poly::add(&right, &vec_poly::mul(&coefficient.eval_n(&n), &derived));
                }
            }
            if !available {
                continue;
            }
            if vec_poly::trim(&left) != vec_poly::trim(&right) {
                return Err(format!("replay fails at base n = {base}"));
            }
            checked += 1;
        }
        Ok(checked)
    }
}

/// Find `H_nu` with `L o T = sum_nu H_nu o S^nu o M`.
///
/// Two exact constructions are tried: `M` proportional to `L o T` up to a
/// factor `c x^e` (then `H_0 = c x^e`), and right division by `M` when its
/// leading row carries a nonzero constant coefficient and no derivative
/// (each step removes the largest shift; the remainder must vanish).
/// `None` when neither applies; no search over multiplier ansatzes is made.
pub fn find_reconstruction_certificate(
    source: &ImplicitRelation,
    map: &ImplicitRelation,
    target: &ImplicitRelation,
) -> Option<ReconstructionCertificate> {
    let l = source.op();
    let t = map.op();
    if t.max_shift() != 0 {
        return None;
    }
    let m = target.op();
    if m.is_zero() {
        return None;
    }
    let product = l.compose(&t);
    let (lead_lt, e_lt, normal_lt) = product.normalize();
    let (lead_m, e_m, normal_m) = m.normalize();
    let mut multipliers: Vec<(usize, DOp)> = Vec::new();
    let construction;
    if normal_lt == normal_m && e_lt >= e_m {
        let factor = Biv::monomial(0, e_lt - e_m, lead_lt / lead_m);
        multipliers.push((0, DOp::single(0, 0, factor)));
        construction = "M is proportional to L o T";
    } else {
        let rm = m.max_shift();
        let leader = m.column(rm);
        let kappa = match leader.0.iter().collect::<Vec<_>>().as_slice() {
            [((0, 0), c)] => c.as_constant().filter(|c| !c.is_zero())?,
            _ => return None,
        };
        let mut remainder = product;
        let mut by_shift: BTreeMap<usize, DOp> = BTreeMap::new();
        while !remainder.is_zero() && remainder.max_shift() >= rm {
            let r = remainder.max_shift();
            let nu = r - rm;
            let h = remainder.column(r).scale(&kappa.recip());
            remainder = remainder.sub(&h.compose(&m.shifted(nu)));
            let entry = by_shift.entry(nu).or_default();
            *entry = entry.add(&h);
        }
        if !remainder.is_zero() {
            return None;
        }
        multipliers = by_shift.into_iter().filter(|(_, h)| !h.is_zero()).collect();
        construction = "right division by the monic leader of M";
    }
    let certificate = ReconstructionCertificate {
        source_relation: source.to_json(),
        map: shift_free_to_json(&t),
        target_relation: target.to_json(),
        multipliers: multipliers
            .iter()
            .map(|(nu, h)| ReconstructionMultiplier {
                shift: *nu,
                differentiations: shift_free_to_json(h),
            })
            .collect(),
        statement: format!(
            "L o T = sum_nu H_nu o S^nu o M ({construction}): every sequence Q satisfying M at \
             the bases n + nu yields rows T Q satisfying L at base n"
        ),
    };
    certificate.verify_identity().ok()?;
    Some(certificate)
}
