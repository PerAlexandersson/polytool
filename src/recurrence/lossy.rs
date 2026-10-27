//! Lossy coordinate maps: plain differentiation `Q_n = D^m P_n`, its lost
//! data, the transport of relations and the description of the source
//! through the image.
//!
//! # The map and its lost data
//!
//! The registry map [`RecurrenceTransformation::Derivative`] sends `P_n` to
//! `Q_n = D^m P_n`. It forgets the coefficients `[x^0..x^(m-1)] P_n`, the
//! *lost data* `L_n = P_n mod x^m` ([`RecurrenceTransformation::lost_data`]).
//! With them the map is inverted exactly,
//! `P_n = I^m Q_n + L_n` with `[x^(k+m)] I^m Q_n = k!/(k+m)! [x^k] Q_n`
//! ([`RecurrenceTransformation::inverse_row_with_lost_data`]); without them
//! [`RecurrenceTransformation::inverse_row`] rejects the map.
//!
//! # Transport
//!
//! Let `L = sum c_(i,d)(n,x) D^d S^i` be a relation of `P` and
//! `r = max(0, max_terms (deg_x c - d))`. By the Leibniz rule
//! `D^(m+r) c(x) D^d = sum_t binom(m+r,t) c^(t)(x) D^(m+r-t+d)`, and every
//! nonzero term has `t <= deg_x c <= r + d`, so its derivative order is at
//! least `m`: the operator `D^(m+r) o L` is right divisible by `D^m`. With
//! `M = D^(m+r) L D^(-m)` (normal order, `D^m` removed from the right of each
//! term) the identity `M o D^m = D^(m+r) o L` holds in the Ore algebra, so
//! `M Q_n = D^(m+r) (L P)_n = 0` at every base where `L` holds. The identity
//! is recorded as a [`ReconstructionCertificate`] (source relation `M`, map
//! `D^m` with `Q = D^m P`, target relation `L`, one multiplier `D^(m+r)`) and
//! replayed by [`ReconstructionCertificate::verify_identity`]. The transported
//! relation is a consequence of `L`; it never determines the lost data.
//!
//! # Pullback and reconstruction
//!
//! A relation `M` of `Q` pulls back to the relation `M o D^m` of `P`, and a
//! coefficient condition of `Q` to one of `P` by
//! `[x^k] Q_n = (k+m)!/k! [x^(k+m)] P_n`. The leader of `M o D^m` kills
//! `x^0..x^(m-1)`, so a uniqueness argument for `P` needs conditions that pin
//! the lost coefficients. [`lossy_reconstruction`] assembles the pulled-back
//! description of `P` with the supplied conditions of `P` (typically
//! `[x^k] (L P)_n = 0` for a certified relation `L`,
//! [`relation_coefficient_condition`]) and again without every condition
//! that refers to a lost coefficient, and assesses both with
//! [`assess_generating_description`]: a description of `Q` plus the lost data
//! can be equivalent, without the lost data it is a consequence only.
//!
//! Related modules: `transform` (the registry), `diagonal` (implicit
//! transport of the bijective implicit maps), `implicit` (relations,
//! descriptions, certificates), `representation_search` (lossy steps).

use super::implicit::{biv_json, encode_poly};
use super::implicit_algebra::{factorial, Biv, DOp};
use super::*;
use crate::vec_poly;
use num_traits::Zero;
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// `n!/(n-d)!` as an exact rational: the factor of `[x^n] P` in
/// `[x^(n-d)] D^d P`.
fn falling(n: usize, d: usize) -> Q {
    if d > n {
        return Q::zero();
    }
    factorial(n) / factorial(n - d)
}

/// `[x^0]` or `[x^0..x^(m-1)]`: the coefficients lost by `D^m`.
pub(super) fn lost_range(order: usize) -> String {
    if order <= 1 {
        "[x^0]".into()
    } else {
        format!("[x^0..x^{}]", order - 1)
    }
}

pub(super) fn derivative_order_check(order: usize) -> Result<(), TransformError> {
    if order == 0 || order > MAX_IMPLICIT_MAP_ORDER {
        return Err(TransformError::InvalidParameter(format!(
            "derivative order must lie in 1..={MAX_IMPLICIT_MAP_ORDER}"
        )));
    }
    Ok(())
}

/// `D^order` of an ascending coefficient vector (zero as `[0]`).
pub(super) fn derivative_row(order: usize, row: &[Q]) -> Vec<Q> {
    let mut out = vec_poly::trim(row);
    for _ in 0..order {
        out = vec_poly::derivative(&out);
    }
    vec_poly::trim(&out)
}

pub(super) fn lossy_inverse_error(index: usize, order: usize) -> TransformError {
    TransformError::ImageCondition {
        index,
        reason: format!(
            "the derivative map is lossy: reconstruction needs the lost coefficients {} P_n \
             (inverse_row_with_lost_data)",
            lost_range(order)
        ),
    }
}

/// Record parts (row map, inverse, conditions, transport) of `derivative`.
pub(super) fn lossy_map_record(order: usize) -> (String, String, Vec<String>, String) {
    (
        format!(
            "Q_n = D^{order} P_n (lossy: the coefficients {} P_n are lost)",
            lost_range(order)
        ),
        format!(
            "P_n = I^{order} Q_n + L_n, [x^(k+{order})] I^{order} Q_n = k!/(k+{order})! [x^k] Q_n, \
             where L_n = P_n mod x^{order} is the recorded lost data"
        ),
        vec![format!(
            "inverse: the lost coefficients {} P_n must be supplied as boundary data",
            lost_range(order)
        )],
        "implicit: M = D^(m+r) L D^(-m) with the identity M o D^m = D^(m+r) o L (a \
         consequence); a description of Q plus conditions pinning the lost coefficients \
         describes P (lossy_reconstruction)"
            .to_string(),
    )
}

impl RecurrenceTransformation {
    /// Whether the map forgets data (`derivative`).
    pub fn is_lossy(&self) -> bool {
        matches!(self, Self::Derivative { .. })
    }

    /// The data a lossy map forgets on the row `P_index`: for `derivative`
    /// of order `m`, the coefficients `[x^0..x^(m-1)] P_index` (ascending,
    /// trimmed, zero as `[0]`). `None` for the bijective maps.
    ///
    /// ```
    /// use polytool::recurrence::RecurrenceTransformation;
    /// use polytool::recurrence::BigRational as Q;
    ///
    /// let d = RecurrenceTransformation::Derivative { order: 1 };
    /// let row: Vec<Q> = [1, 3, 1].iter().map(|&c| Q::from_integer(c.into())).collect();
    /// let image = d.map_row(2, &row).unwrap();
    /// let lost = d.lost_data(2, &row).unwrap().unwrap();
    /// assert_eq!(lost, vec![Q::from_integer(1.into())]);
    /// assert!(d.inverse_row(2, &image).is_err());
    /// assert_eq!(d.inverse_row_with_lost_data(2, &image, Some(&lost)).unwrap(), row);
    /// ```
    pub fn lost_data(&self, index: usize, row: &[Q]) -> Result<Option<Vec<Q>>, TransformError> {
        let _ = index;
        self.validate()?;
        match self {
            Self::Derivative { order } => {
                let row = vec_poly::trim(row);
                let low = row.iter().take(*order).cloned().collect::<Vec<_>>();
                Ok(Some(vec_poly::trim(&low)))
            }
            _ => Ok(None),
        }
    }

    /// Reconstruct `P_index` from its image and, for a lossy map, the lost
    /// data returned by [`RecurrenceTransformation::lost_data`]. For the
    /// bijective maps `lost` must be `None` and this is
    /// [`RecurrenceTransformation::inverse_row`].
    pub fn inverse_row_with_lost_data(
        &self,
        index: usize,
        components: &[Vec<Q>],
        lost: Option<&[Q]>,
    ) -> Result<Vec<Q>, TransformError> {
        self.validate()?;
        match (self, lost) {
            (Self::Derivative { order }, Some(lost)) => {
                if components.len() != 1 {
                    return Err(TransformError::ComponentCount {
                        expected: 1,
                        found: components.len(),
                    });
                }
                let lost = vec_poly::trim(lost);
                if vec_poly::degree(&lost).is_some_and(|d| d >= *order) {
                    return Err(TransformError::ImageCondition {
                        index,
                        reason: format!("the lost data must have degree below {order}"),
                    });
                }
                let image = vec_poly::trim(&components[0]);
                let mut row = vec![Q::zero(); image.len() + order];
                for (k, c) in image.iter().enumerate() {
                    row[k + order] = c / falling(k + order, *order);
                }
                for (k, c) in lost.iter().enumerate() {
                    row[k] = c.clone();
                }
                Ok(vec_poly::trim(&row))
            }
            (Self::Derivative { order }, None) => Err(lossy_inverse_error(index, *order)),
            (_, None) => self.inverse_row(index, components),
            (_, Some(_)) => Err(TransformError::InvalidParameter(format!(
                "{} is not lossy; no lost data is expected",
                self.kind_name()
            ))),
        }
    }
}

/// `D^order` as a shift-free operator.
fn derivative_op(order: usize) -> DOp {
    DOp::single(0, order, Biv::one())
}

/// Transport a relation of `P` through `Q = D^order P` (see the module notes).
pub(super) fn transport_derivative(
    relation: &ImplicitRelation,
    order: usize,
) -> Result<ImplicitTransport, TransformError> {
    derivative_order_check(order)?;
    let l = relation.op();
    if l.is_zero() {
        return Err(TransformError::TransportUnsupported(
            "the zero relation".into(),
        ));
    }
    let r =
        l.0.iter()
            .map(|((_, d), c)| c.y_degree().unwrap_or(0) as i64 - *d as i64)
            .max()
            .unwrap_or(0)
            .max(0) as usize;
    if order + r > 4 * MAX_IMPLICIT_MAP_ORDER {
        return Err(TransformError::Limit(
            "derivative transport order too large".into(),
        ));
    }
    let multiplier = derivative_op(order + r);
    let product = multiplier.compose(&l);
    let mut quotient = DOp::default();
    for ((i, d), c) in &product.0 {
        let reduced = d.checked_sub(order).ok_or_else(|| {
            TransformError::Limit("internal: D^(m+r) L is not right divisible by D^m".into())
        })?;
        quotient.add_term((*i, reduced), c.clone());
    }
    let (_, _, normal) = quotient.normalize();
    let target = ImplicitRelation::from_op(&normal);
    let unit = biv_json(&Biv::one());
    let certificate = ReconstructionCertificate {
        source_relation: ImplicitRelation::from_op(&quotient).to_json(),
        map: vec![OperatorTermJson {
            deriv_order: order,
            coeff: unit.clone(),
        }],
        target_relation: relation.to_json(),
        multipliers: vec![ReconstructionMultiplier {
            shift: 0,
            differentiations: vec![OperatorTermJson {
                deriv_order: order + r,
                coeff: unit,
            }],
        }],
        statement: format!(
            "M o D^{order} = D^{} o L with M = D^{} L D^-{order} (lossy transport): every P \
             satisfying L at base n gives Q = D^{order} P satisfying M at base n; the \
             transported relation is c x^e M for its canonical form",
            order + r,
            order + r
        ),
    };
    certificate
        .verify_identity()
        .map_err(|e| TransformError::Limit(format!("internal: {e}")))?;
    Ok(ImplicitTransport {
        relation: target,
        exceptional_conditions: Vec::new(),
        clearing: None,
        reconstruction: Some(certificate),
    })
}

/// The relation `M o D^order` of `P` for a relation `M` of `Q = D^order P`.
pub fn pull_back_relation_through_derivative(
    relation: &ImplicitRelation,
    order: usize,
) -> ImplicitRelation {
    ImplicitRelation::from_op(&relation.op().compose(&derivative_op(order)))
}

/// The condition of `P` equivalent to a condition of `Q = D^order P`, by
/// `[x^k] Q_n = (k+order)!/k! [x^(k+order)] P_n`.
pub fn pull_back_condition_through_derivative(
    condition: &CoefficientCondition,
    order: usize,
) -> CoefficientCondition {
    CoefficientCondition {
        terms: condition
            .terms
            .iter()
            .map(|t| ConditionTerm {
                shift: t.shift,
                coefficient_index: t.coefficient_index + order,
                coeff: vec_poly::scale(&t.coeff, &falling(t.coefficient_index + order, order)),
            })
            .collect(),
    }
    .canonical()
}

/// The condition of `Q = D^order P` equivalent to a condition of `P`, or
/// `None` when the condition refers to a lost coefficient `[x^k] P`,
/// `k < order`.
pub fn push_condition_through_derivative(
    condition: &CoefficientCondition,
    order: usize,
) -> Option<CoefficientCondition> {
    let canonical = condition.canonical();
    if canonical.terms.iter().any(|t| t.coefficient_index < order) {
        return None;
    }
    Some(
        CoefficientCondition {
            terms: canonical
                .terms
                .iter()
                .map(|t| ConditionTerm {
                    shift: t.shift,
                    coefficient_index: t.coefficient_index - order,
                    coeff: vec_poly::scale(&t.coeff, &falling(t.coefficient_index, order).recip()),
                })
                .collect(),
        }
        .canonical(),
    )
}

/// Whether a condition refers to a coefficient `[x^k]` with `k < order`.
pub fn condition_refers_below(condition: &CoefficientCondition, order: usize) -> bool {
    condition
        .canonical()
        .terms
        .iter()
        .any(|t| t.coefficient_index < order)
}

/// The coefficient condition `[x^k] (L P)_n = 0` of a relation `L`:
/// `sum c_(u,j) n^u (k-j+d)!/(k-j)! [x^(k-j+d)] P_(n+i) = 0`, summed over the
/// terms `c(n,x) D^d P_(n+i)` of `L` and the powers `x^j` with `j <= k`. It
/// holds wherever `L` holds.
pub fn relation_coefficient_condition(
    relation: &ImplicitRelation,
    k: usize,
) -> CoefficientCondition {
    let mut terms = Vec::new();
    for ((i, d), c) in &relation.op().0 {
        let width = c.y_degree().map_or(0, |deg| deg + 1);
        for j in 0..width.min(k + 1) {
            let coefficient = c.y_coefficient(j);
            if vec_poly::is_zero(&coefficient) {
                continue;
            }
            let index = k - j + d;
            terms.push(ConditionTerm {
                shift: *i,
                coefficient_index: index,
                coeff: vec_poly::scale(&coefficient, &falling(index, *d)),
            });
        }
    }
    CoefficientCondition { terms }.canonical()
}

/// A description of the source `P` assembled from a description of
/// `Q = D^m P`, with and without the conditions that pin the lost data.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct LossyReconstruction {
    /// The derivative order `m`.
    pub order: usize,
    /// Description of `P`: the relations of `Q` pulled back (`M o D^m`), the
    /// conditions of `Q` pulled back (`[x^k] Q = (k+m)!/k! [x^(k+m)] P`) and
    /// the supplied conditions of `P` (boundary conditions).
    pub with_lost_data: GeneratingDescription,
    /// Assessment of `with_lost_data`.
    pub with_lost_data_assessment: DescriptionAssessment,
    /// `with_lost_data` without every condition that refers to a lost
    /// coefficient `[x^k] P`, `k < m`.
    pub without_lost_data: GeneratingDescription,
    /// Assessment of `without_lost_data`.
    pub without_lost_data_assessment: DescriptionAssessment,
    /// The conditions of `with_lost_data` that refer to lost coefficients, as
    /// `list[i]` references.
    pub lost_data_conditions: Vec<String>,
    /// What the two assessments show.
    pub statement: String,
}

/// Assemble and assess the description of `P` through `Q = D^order P`.
///
/// `q_description` describes `Q` (its relations and conditions; its seeds
/// are not used). `source_rows[j]` is `P_(first_index + j)` and supplies the
/// seeds of `P`; `source_conditions` are conditions of `P` with their
/// evidence, typically `[x^k] (L P)_n = 0` for a certified relation `L`.
pub fn lossy_reconstruction(
    order: usize,
    q_description: &GeneratingDescription,
    first_index: usize,
    source_rows: &[Vec<Q>],
    source_conditions: &[DescribedCondition],
) -> Result<LossyReconstruction, String> {
    derivative_order_check(order).map_err(|e| e.to_string())?;
    let symbol = "P";
    let mut relations = Vec::new();
    let mut seed_count = q_description.seeds.len().max(1);
    for record in &q_description.relations {
        let relation = ImplicitRelation::from_json(&record.relation)?;
        let pulled = pull_back_relation_through_derivative(&relation, order);
        seed_count = seed_count.max(
            (record.first_base + pulled.max_shift())
                .saturating_sub(first_index)
                .max(1),
        );
        relations.push(
            DescribedRelation::new(&pulled, record.first_base, &record.evidence, symbol).with_note(
                format!("pulled back from {} through D^{order}", record.text),
            ),
        );
    }
    let pull = |records: &[DescribedCondition]| -> Result<Vec<DescribedCondition>, String> {
        records
            .iter()
            .map(|record| {
                let condition = CoefficientCondition::from_json(&record.condition)?;
                let pulled = pull_back_condition_through_derivative(&condition, order);
                Ok(DescribedCondition::new(
                    &pulled,
                    record.first_base,
                    &record.evidence,
                    record.coefficient_index.map(|k| k + order),
                    symbol,
                )
                .with_note(format!(
                    "pulled back from {} through D^{order}",
                    record.text
                )))
            })
            .collect()
    };
    let exceptional = pull(&q_description.exceptional_index_conditions)?;
    let mut boundary = pull(&q_description.boundary_conditions)?;
    for record in source_conditions {
        if !exceptional
            .iter()
            .chain(&boundary)
            .any(|c| c.condition == record.condition)
        {
            boundary.push(record.clone());
        }
    }
    if source_rows.len() < seed_count {
        return Err(format!(
            "lossy reconstruction needs {seed_count} source rows for the seeds"
        ));
    }
    let with_lost_data = GeneratingDescription {
        first_index,
        seeds: source_rows[..seed_count]
            .iter()
            .map(|r| encode_poly(r))
            .collect(),
        relations,
        exceptional_index_conditions: exceptional,
        boundary_conditions: boundary,
    };
    let lost = |records: &[DescribedCondition], list: &str| -> Result<Vec<String>, String> {
        let mut labels = Vec::new();
        for (i, record) in records.iter().enumerate() {
            if condition_refers_below(&CoefficientCondition::from_json(&record.condition)?, order) {
                labels.push(format!("{list}[{i}]"));
            }
        }
        Ok(labels)
    };
    let mut lost_data_conditions = lost(
        &with_lost_data.exceptional_index_conditions,
        "exceptional_index_conditions",
    )?;
    lost_data_conditions.extend(lost(
        &with_lost_data.boundary_conditions,
        "boundary_conditions",
    )?);
    let keep = |records: &[DescribedCondition]| -> Result<Vec<DescribedCondition>, String> {
        let mut kept = Vec::new();
        for record in records {
            if !condition_refers_below(&CoefficientCondition::from_json(&record.condition)?, order)
            {
                kept.push(record.clone());
            }
        }
        Ok(kept)
    };
    let without_lost_data = GeneratingDescription {
        exceptional_index_conditions: keep(&with_lost_data.exceptional_index_conditions)?,
        boundary_conditions: keep(&with_lost_data.boundary_conditions)?,
        ..with_lost_data.clone()
    };
    let with_assessment = assess_generating_description(&with_lost_data)?;
    let without_assessment = assess_generating_description(&without_lost_data)?;
    let statement = format!(
        "Q_n = D^{order} P_n forgets {} P_n. The relations and conditions of Q pulled \
         back to P, with the {} condition(s) that refer to lost coefficients, give a description \
         of P assessed as {}; without them it is assessed as {}.",
        lost_range(order),
        lost_data_conditions.len(),
        with_assessment.status,
        without_assessment.status
    );
    Ok(LossyReconstruction {
        order,
        with_lost_data,
        with_lost_data_assessment: with_assessment,
        without_lost_data,
        without_lost_data_assessment: without_assessment,
        lost_data_conditions,
        statement,
    })
}

/// The conditions `[x^k] (L P)_n = 0` (`k` in `indices`) of a described
/// relation `L` of `P`, with its domain and evidence; vacuous conditions are
/// skipped.
pub(super) fn coefficient_conditions_of(
    relation: &DescribedRelation,
    indices: &[usize],
    symbol: &str,
) -> Result<Vec<DescribedCondition>, String> {
    let parsed = ImplicitRelation::from_json(&relation.relation)?;
    let mut out = Vec::new();
    for &k in indices {
        let condition = relation_coefficient_condition(&parsed, k);
        if condition.is_vacuous() {
            continue;
        }
        out.push(
            DescribedCondition::new(
                &condition,
                relation.first_base,
                &relation.evidence,
                Some(k),
                symbol,
            )
            .with_note(format!("[x^{k}] of the relation {}", relation.text)),
        );
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(n: i64) -> Q {
        Q::from_integer(n.into())
    }

    #[test]
    fn derivative_transport_of_pascal_rows() {
        // P_(n+1) = (1 + x) P_n, Q = D P: Q'_(n+1) - (1 + x) Q'_n - 2 Q_n = 0.
        let mut op = DOp::default();
        op.add_term((1, 0), Biv::one());
        op.add_term((0, 0), Biv::from_y(&[q(-1), q(-1)]));
        let l = ImplicitRelation::from_op(&op);
        let transport = transport_derivative(&l, 1).unwrap();
        assert_eq!(
            transport.relation.text("Q"),
            "D Q[n+1] + (-1 - x) D Q[n] - 2 Q[n] = 0"
        );
        let certificate = transport.reconstruction.unwrap();
        certificate.verify_identity().unwrap();
        // Rows Q_n = n (1 + x)^(n-1).
        let rows = (0..6)
            .map(|n| {
                let p = (0..n).fold(vec![q(1)], |acc, _| vec_poly::mul(&acc, &[q(1), q(1)]));
                derivative_row(1, &p)
            })
            .collect::<Vec<_>>();
        for base in 0..5 {
            assert_eq!(transport.relation.holds_at(&rows, 0, base), Some(true));
        }
        // [x^0] (L P)_n = [x^0] P_(n+1) - [x^0] P_n.
        let c = relation_coefficient_condition(&l, 0);
        assert_eq!(c.text("P"), "[x^0] P[n+1] - [x^0] P[n] = 0");
        assert!(condition_refers_below(&c, 1));
        assert!(push_condition_through_derivative(&c, 1).is_none());
        let pulled = pull_back_condition_through_derivative(
            &CoefficientCondition {
                terms: vec![ConditionTerm {
                    shift: 1,
                    coefficient_index: 1,
                    coeff: vec![q(1)],
                }],
            },
            2,
        );
        // [x^1] Q = 3!/1! [x^3] P.
        assert_eq!(pulled.terms[0].coefficient_index, 3);
        assert_eq!(pulled.terms[0].coeff, vec![q(6)]);
        assert_eq!(
            push_condition_through_derivative(&pulled, 2).unwrap(),
            CoefficientCondition {
                terms: vec![ConditionTerm {
                    shift: 1,
                    coefficient_index: 1,
                    coeff: vec![q(1)],
                }],
            }
        );
    }
}
