//! Certified rowwise products of two rule-defined polynomial families.
//!
//! Hadamard closure uses invariant coefficient relations of the two rules,
//! `p(n,k) a(n,k) + q(n,k) a(n,k-1) = 0`. Their product gives the coefficient
//! relation of the image (with weight `1`, `k!`, or `1/k!`). An invariant
//! annihilator system proves that relation for a proposed image rule. A
//! constant coefficient and a degree envelope make the relation unique.
//! Ordinary products use the tensor product of shift-only companion systems.
//! The image rule is then compared with the requested rule by certified Ore
//! recurrence equality. A fit supplies candidates only; it never proves one.
//!
//! The bounded recognizer can return `Unsupported`; it makes no claim about
//! families outside these exactly decidable classes. In particular a zero of
//! a coefficient-relation leader cannot be ignored. No interlacing is claimed.

use super::hadamard_factor::prepare::{prepare_lift_rule, PreparedLiftRule};
use crate::recurrence::implicit::algebra::{Biv, DOp};
use crate::recurrence::*;
use crate::{vec_poly, BigRational};
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};

mod coefficient;
mod ordinary;
#[cfg(test)]
mod tests;

type Q = BigRational;

/// Schema of [`PolynomialFamilyProductCertificate`].
pub const POLYNOMIAL_FAMILY_PRODUCT_SCHEMA: &str = "polytool.polynomial-family-product.v1";

/// The rowwise operation on the two parent families.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PolynomialFamilyProductOperation {
    /// `[x^k] Q_n = [x^k] A_n [x^k] B_n` (Maló).
    Hadamard,
    /// `[x^k] Q_n = k! [x^k] A_n [x^k] B_n` (Schur).
    FactorialHadamard,
    /// `[x^k] Q_n = [x^k] A_n [x^k] B_n / k!` (Maló and division by factorials).
    InverseFactorialHadamard,
    /// `Q_n = A_n B_n` (multiplication of split polynomials).
    Ordinary,
}

/// Two defining records and their index alignment with the image.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyProductCandidate {
    /// Definition of `A`.
    pub left: RecurrenceJson,
    /// Definition of `B`.
    pub right: RecurrenceJson,
    /// The image uses `A_(n + left_index_shift)`.
    pub left_index_shift: i64,
    /// The image uses `B_(n + right_index_shift)`.
    pub right_index_shift: i64,
    /// Operation applied at each index.
    pub operation: PolynomialFamilyProductOperation,
}

/// Bounded construction options. These bounds affect applicability, never soundness.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyProductOptions {
    /// Rows used to find a proposed image rule and coefficient relations (8..=40).
    pub search_rows: usize,
    /// Candidate budget for the image recurrence search (1..=2000).
    pub search_candidates: usize,
}

impl Default for PolynomialFamilyProductOptions {
    fn default() -> Self {
        Self {
            search_rows: 18,
            search_candidates: 600,
        }
    }
}

/// Replayable proof of a parent's root hypotheses.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyProductParentProof {
    /// All-index real-rootedness with coefficients of one sign, when proved.
    pub nonpositive: VeroneseParentProof,
    /// A whole-line theorem, used when the nonpositive-root statement is unavailable.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub whole_line: Option<RealRootednessTheoremReport>,
    /// A fixed positive constant coefficient, proved by induction on the rule;
    /// combined with `nonpositive`, this proves PF, including its sign condition.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub positive_constant_coefficient: Option<String>,
}

impl PolynomialFamilyProductParentProof {
    fn is_pf(&self) -> bool {
        self.nonpositive.status == VeroneseParentProofStatus::Proved
            && self.positive_constant_coefficient.is_some()
    }

    fn is_real_rooted(&self) -> bool {
        self.nonpositive.status == VeroneseParentProofStatus::Proved
            || self
                .whole_line
                .as_ref()
                .is_some_and(RealRootednessTheoremReport::is_proved)
    }
}

/// The complete all-index identity and published root transfer for two families.
/// The statement permits zero image rows; no output nonvanishing is inferred
/// from a finite sample.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PolynomialFamilyProductCertificate {
    /// [`POLYNOMIAL_FAMILY_PRODUCT_SCHEMA`].
    pub schema: String,
    /// Both source definitions, operation and index shifts.
    pub candidate: PolynomialFamilyProductCandidate,
    /// Definition of the requested target family.
    pub target: RecurrenceJson,
    /// Certified definition of the forward image.
    pub image: RecurrenceJson,
    /// Construction and replay bounds.
    pub options: PolynomialFamilyProductOptions,
    /// Proof of the left parent, at the aligned indices actually used.
    pub left_proof: PolynomialFamilyProductParentProof,
    /// Proof of the right parent, at the aligned indices actually used.
    pub right_proof: PolynomialFamilyProductParentProof,
    /// Annihilator closure, with the identifying and exceptional row checks.
    pub closure: PolynomialFamilyProductClosure,
    /// Certified all-index equality of the forward image and the target rule.
    pub equality: RecurrenceEqualityCertificate,
    /// Published theorem and its exact hypotheses.
    pub theorem: String,
    /// `pf` or `zero_or_real_rooted`; no interlacing or strictness claim.
    pub conclusion: String,
    /// Source-inspected Lean declarations; a correspondence, not a Lean replay.
    pub lean_declarations: Vec<String>,
    /// `not_attempted`: no formal proof was generated or checked.
    pub lean_replay: String,
}

/// The algebraic proof of the forward image identity.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum PolynomialFamilyProductClosure {
    /// Hadamard closure of certified coefficient relations of the two parents.
    CoefficientRelations {
        /// Left coefficient relation and its invariant annihilator proof.
        left: Box<coefficient::PolynomialFamilyCoefficientRelation>,
        /// Right coefficient relation and its invariant annihilator proof.
        right: Box<coefficient::PolynomialFamilyCoefficientRelation>,
        /// Product relation, certified against the image rule.
        image: Box<coefficient::PolynomialFamilyCoefficientRelation>,
    },
    /// Tensor product of shift-only companion systems and component elimination.
    Ordinary {
        /// Exact identity and finite exceptional-base checks.
        certificate: Box<ordinary::PolynomialFamilyOrdinaryClosure>,
    },
    /// A rule-defined parent proved equal to binomial rows, followed by an
    /// independently certified Hadamard-factor conjugation. This also handles
    /// a non-hypergeometric other parent such as the Eulerian family.
    BinomialFactor {
        /// Whether the binomial parent is the right input.
        right_is_binomial: bool,
        /// All-index equality with the binomial rule.
        binomial_equality: Box<RecurrenceEqualityCertificate>,
        /// Conjugation certificate for the image and the other parent.
        certificate: Box<HadamardFactorCertificate>,
    },
    /// A certified coefficient relation identifies one rule-defined source
    /// with a hypergeometric multiplier; conjugation proves the product.
    CoefficientMultiplier {
        /// Whether the represented source is the right input.
        right_is_multiplier: bool,
        /// All-index relation, support, nonsingular leader and anchor of that source.
        source: Box<coefficient::PolynomialFamilyCoefficientRelation>,
        /// The rational identity between its coefficient ratio and the
        /// multiplier's exact binomial ratio, checked symbolically on replay.
        ratio_identity: String,
        /// Exact conjugation certificate, including the normalization gauge.
        certificate: Box<HadamardFactorCertificate>,
    },
}

pub use coefficient::PolynomialFamilyCoefficientRelation;
pub use ordinary::PolynomialFamilyOrdinaryClosure;

fn unsupported(reason: impl Into<String>) -> RecurrenceError {
    RecurrenceError::Unsupported(reason.into())
}

fn prepare_record(
    record: &RecurrenceJson,
    name: &str,
) -> Result<PreparedLiftRule, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    prepare_lift_rule(&rule, &seeds, first, name)
}

fn record(parent: &PreparedLiftRule) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(&parent.rule, parent.first_index, &parent.seeds, None)
}

/// Align by taking rows of the *original* definition, then substituting n+s
/// in the rule. An exceptional supplied prefix remains part of the source;
/// the requested tail must start late enough to admit a runnable definition.
fn aligned(
    record: &RecurrenceJson,
    shift: i64,
    first: usize,
) -> Result<PreparedLiftRule, RecurrenceError> {
    if shift.unsigned_abs() > 64 || first > 256 {
        return Err(unsupported("index alignment exceeds the product bounds"));
    }
    let (rule, source_first, seeds) = record.to_recurrence_parts()?;
    let index = i64::try_from(first)
        .ok()
        .and_then(|n| n.checked_add(shift))
        .and_then(|n| usize::try_from(n).ok())
        .filter(|n| *n >= source_first)
        .ok_or_else(|| unsupported("a product uses a row before its source definition"))?;
    let skip = index - source_first;
    if skip > 256 {
        return Err(unsupported("source alignment needs too many rows"));
    }
    let count = seeds.len() + rule.max_offset() + 2;
    let rows = rule
        .generate_rows_rational(&seeds, source_first, skip + count)
        .map_err(|error| unsupported(format!("source rows cannot be generated: {error}")))?;
    prepare_lift_rule(
        &rule.shift_index(shift),
        &rows[skip..],
        first,
        "aligned parent",
    )
}

fn constant_coefficient(parent: &PreparedLiftRule) -> Option<Q> {
    let constant = parent
        .seeds
        .first()?
        .first()
        .cloned()
        .unwrap_or_else(Q::zero);
    if parent
        .seeds
        .iter()
        .any(|row| row.first().cloned().unwrap_or_else(Q::zero) != constant)
    {
        return None;
    }
    let mut sum = Biv::zero();
    for term in &parent.rule.terms {
        let coefficient = Biv::new(term.coeff.coeffs.clone());
        let column = coefficient.variable_coefficient(0);
        if term.deriv_order > 0 {
            if column.iter().any(|c| !c.is_zero()) {
                return None;
            }
        } else {
            sum = sum.add(&Biv::from_index_polynomial(&column));
        }
    }
    let leader = parent
        .rule
        .denominator
        .as_ref()
        .map_or_else(Biv::one, |p| Biv::new(p.coeffs.clone()));
    (sum == leader).then_some(constant)
}

fn parent_proof(
    parent: &PreparedLiftRule,
) -> Result<PolynomialFamilyProductParentProof, RecurrenceError> {
    let nonpositive = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &parent.rule,
            seeds: &parent.seeds,
            first_index: parent.first_index,
        },
        parent.first_index,
        &[],
        true,
    )?;
    let whole_line = if nonpositive.status != VeroneseParentProofStatus::Proved {
        let proof = check_three_term_recurrence_theorem(
            &parent.rule,
            parent.first_index,
            &parent.seeds,
            &RealRootednessTheoremOptions {
                root_region: None,
                ..Default::default()
            },
        )?;
        proof.is_proved().then_some(proof)
    } else {
        None
    };
    let positive_constant_coefficient = constant_coefficient(parent)
        .filter(Signed::is_positive)
        .map(|value| format_rational_coeff(&value));
    Ok(PolynomialFamilyProductParentProof {
        nonpositive,
        whole_line,
        positive_constant_coefficient,
    })
}

fn transfer(
    operation: PolynomialFamilyProductOperation,
    left: &PolynomialFamilyProductParentProof,
    right: &PolynomialFamilyProductParentProof,
) -> Result<(String, String, Vec<String>), RecurrenceError> {
    use PolynomialFamilyProductOperation::*;
    let (theorem, conclusion, declarations): (&str, &str, &[&str]) = match operation {
        Ordinary if left.is_real_rooted() && right.is_real_rooted() => (
            "The product of two split real polynomials is split (Polynomial.Splits.mul); both parents are proved real-rooted at every used index.",
            "zero_or_real_rooted", &["Polynomial.Splits.mul", "RealRooted.productSequence_splits"]),
        Hadamard if left.is_pf() && right.is_pf() => (
            "Malo; Garloff-Wagner (1996), Theorem 4(a): the coefficientwise product of two PF polynomials is PF. Each parent has proved nonpositive real roots and a proved positive constant coefficient, hence nonnegative coefficients.",
            "pf", &["RealRooted.IsPFPolynomial.hadamardProduct", "RealRooted.Tactic.hadamardProduct_sequence_pf"]),
        FactorialHadamard if (left.is_real_rooted() && right.is_pf()) || (right.is_real_rooted() && left.is_pf()) => (
            "Schur (1914): sum_k k! a_k b_k x^k is zero or split when one input splits and the other is PF. Commutativity permits either input to supply the PF hypothesis.",
            if left.is_pf() && right.is_pf() { "pf" } else { "zero_or_real_rooted" },
            &["RealRooted.factorialHadamardProduct_eq_zero_or_splits_of_isPFPolynomial", "RealRooted.IsPFPolynomial.factorialHadamardProduct"]),
        InverseFactorialHadamard if left.is_pf() && right.is_pf() => (
            "Malo's PF Hadamard closure followed by PF-preserving division of coefficient k by k! (Polya-Schur). Both parents are proved PF; no arbitrary real-rooted input is admitted.",
            "pf", &["RealRooted.IsPFPolynomial.hadamardProduct", "RealRooted.IsPFPolynomial.divFactorial"]),
        _ => return Err(unsupported("the published product theorem's all-index parent hypotheses are not proved")),
    };
    Ok((
        theorem.into(),
        conclusion.into(),
        declarations.iter().map(|s| (*s).into()).collect(),
    ))
}

fn product_row(left: &[Q], right: &[Q], operation: PolynomialFamilyProductOperation) -> Vec<Q> {
    if operation == PolynomialFamilyProductOperation::Ordinary {
        return vec_poly::trim(&vec_poly::mul(left, right));
    }
    let mut factorial = Q::one();
    vec_poly::trim(
        &left
            .iter()
            .zip(right)
            .enumerate()
            .map(|(k, (a, b))| {
                if k > 0 {
                    factorial *= Q::from_integer(k.into());
                }
                match operation {
                    PolynomialFamilyProductOperation::FactorialHadamard => a * b * &factorial,
                    PolynomialFamilyProductOperation::InverseFactorialHadamard => {
                        a * b / &factorial
                    }
                    _ => a * b,
                }
            })
            .collect::<Vec<_>>(),
    )
}

fn image_rows(
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
    operation: PolynomialFamilyProductOperation,
    count: usize,
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    Ok(left
        .rows(count)?
        .iter()
        .zip(right.rows(count)?)
        .map(|(a, b)| product_row(a, &b, operation))
        .collect())
}

fn fit_image(
    left: &PreparedLiftRule,
    right: &PreparedLiftRule,
    operation: PolynomialFamilyProductOperation,
    options: &PolynomialFamilyProductOptions,
) -> Result<PreparedLiftRule, RecurrenceError> {
    let rows = image_rows(left, right, operation, options.search_rows)?;
    let search = AdaptiveSearchOptions {
        max_rec_len: 4,
        max_var_deg: 4,
        max_idx_deg: 3,
        max_diff_deg: 2,
        try_denominator: true,
        max_denom_var_deg: 0,
        max_denom_idx_deg: 3,
        ..Default::default()
    };
    let fallback_shapes = [
        (1, 2, 2, 1, 1),
        (2, 3, 3, 2, 2),
        (2, 3, 3, 1, 2),
        (3, 3, 3, 1, 2),
    ];
    let reserved = fallback_shapes
        .len()
        .min(options.search_candidates.saturating_sub(1));
    if let AdaptiveSearchOutcome::Found(found) = find_recurrence_adaptive_rational_with_budget(
        &rows,
        &search,
        AdaptiveSearchBudget::limited(options.search_candidates - reserved),
    ) {
        let rule = found.recurrence.shift_index(1 - left.first_index as i64);
        if let Ok(prepared) = prepare_lift_rule(&rule, &rows, left.first_index, "product image") {
            return Ok(prepared);
        }
    }
    // Reserve a few larger differential shapes so a long list of small
    // impossible shapes cannot consume the entire deterministic budget.
    // Every fit still needs the independent all-index closure proof.
    for &(rec_len, var_deg, idx_deg, diff_deg, denom_idx_deg) in
        fallback_shapes.iter().take(reserved)
    {
        let shape = RecurrenceOptions {
            rec_len,
            var_deg,
            idx_deg,
            diff_deg,
            denom_idx_deg,
            ..Default::default()
        };
        if let Some(rule) = find_polynomial_recurrence_rational(&rows, &shape) {
            let rule = rule.shift_index(1 - left.first_index as i64);
            if let Ok(prepared) = prepare_lift_rule(&rule, &rows, left.first_index, "product image")
            {
                return Ok(prepared);
            }
        }
    }
    Err(unsupported(
        "no product image rule found within the search budget",
    ))
}

fn validate_options(options: &PolynomialFamilyProductOptions) -> Result<(), RecurrenceError> {
    if !(8..=40).contains(&options.search_rows) || !(1..=2000).contains(&options.search_candidates)
    {
        return Err(RecurrenceError::LimitExceeded(
            "invalid product search bounds".into(),
        ));
    }
    Ok(())
}

/// Certify a requested target as a rowwise product of two proved families.
/// Parent records are never accepted on attestation. Index shifts are explicit;
/// all identifying rows, exceptional indices and theorem hypotheses are exact.
///
/// # Errors
/// Returns [`RecurrenceError::Unsupported`] if a hypothesis, annihilator closure,
/// or recurrence equality cannot be proved within the documented class/bounds;
/// malformed definitions and invalid bounds return their typed errors.
pub fn certify_polynomial_family_product(
    target: &RecurrenceJson,
    candidate: &PolynomialFamilyProductCandidate,
    options: &PolynomialFamilyProductOptions,
) -> Result<PolynomialFamilyProductCertificate, RecurrenceError> {
    validate_options(options)?;
    let target_prepared = prepare_record(target, "target")?;
    let first = target_prepared.first_index;
    let left = aligned(&candidate.left, candidate.left_index_shift, first)?;
    let right = aligned(&candidate.right, candidate.right_index_shift, first)?;
    let left_proof = parent_proof(&left)?;
    let right_proof = parent_proof(&right)?;
    let (theorem, conclusion, lean_declarations) =
        transfer(candidate.operation, &left_proof, &right_proof)?;
    let count = target_prepared.seeds.len().max(options.search_rows);
    if target_prepared.rows(count)? != image_rows(&left, &right, candidate.operation, count)? {
        return Err(unsupported(
            "the target differs from the proposed product on identifying/search rows",
        ));
    }
    // Prefer the supplied target rule as the image candidate: it too must pass
    // the closure proof. This avoids making certification depend on a fit.
    let closure = if candidate.operation == PolynomialFamilyProductOperation::Ordinary {
        PolynomialFamilyProductClosure::Ordinary {
            certificate: Box::new(ordinary::certify(&target_prepared, &left, &right)?),
        }
    } else {
        match coefficient::binomial_closure(&target_prepared, &left, &right, candidate.operation) {
            Ok(closure) => closure,
            Err(first_error) => coefficient::certify(
                &target_prepared,
                &left,
                &right,
                candidate.operation,
                options,
            )
            .map_err(|second_error| {
                unsupported(format!("{second_error}; binomial closure: {first_error}"))
            })?,
        }
    };
    let image = record(&target_prepared);
    let equality = certify_recurrence_equality(target, &image)?;
    if !equality.is_certified() {
        return Err(unsupported(
            "the target and product image rules were not certified equal",
        ));
    }
    Ok(PolynomialFamilyProductCertificate {
        schema: POLYNOMIAL_FAMILY_PRODUCT_SCHEMA.into(),
        candidate: candidate.clone(),
        target: target.clone(),
        image,
        options: options.clone(),
        left_proof,
        right_proof,
        closure,
        equality,
        theorem,
        conclusion,
        lean_declarations,
        lean_replay: "not_attempted".into(),
    })
}

/// Construct a forward product family from two defining records. The returned
/// recurrence is certified by annihilator closure; fitted rows alone cannot
/// produce a certificate. `first_index` selects the first common image row.
///
/// # Errors
/// Errors of [`certify_polynomial_family_product`], and `Unsupported` when no
/// candidate image rule is found within the bounded search.
pub fn apply_polynomial_family_product(
    candidate: &PolynomialFamilyProductCandidate,
    first_index: usize,
    options: &PolynomialFamilyProductOptions,
) -> Result<PolynomialFamilyProductCertificate, RecurrenceError> {
    validate_options(options)?;
    let left = aligned(&candidate.left, candidate.left_index_shift, first_index)?;
    let right = aligned(&candidate.right, candidate.right_index_shift, first_index)?;
    // Fail before fitting if there is no published root transfer.
    transfer(
        candidate.operation,
        &parent_proof(&left)?,
        &parent_proof(&right)?,
    )?;
    let image = fit_image(&left, &right, candidate.operation, options)?;
    certify_polynomial_family_product(&record(&image), candidate, options)
}

impl PolynomialFamilyProductCertificate {
    /// Replay the complete derivation from the two parent definitions, target
    /// rule and recorded bounds, rejecting altered algebra, shifts, hypotheses,
    /// conclusions or identifying rows.
    ///
    /// # Errors
    /// [`RecurrenceError::VerificationFailed`] for any changed certificate, or
    /// errors when a recorded definition no longer admits the claimed proof.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != POLYNOMIAL_FAMILY_PRODUCT_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(
                "unknown product certificate schema".into(),
            ));
        }
        let recomputed =
            certify_polynomial_family_product(&self.target, &self.candidate, &self.options)?;
        if &recomputed != self {
            return Err(RecurrenceError::VerificationFailed(
                "the product certificate differs from exact replay".into(),
            ));
        }
        Ok(())
    }
}
