//! Even-derivative maps `T_c f = sum_k c_k x^k f^(2k)`.
//!
//! Ratios and finite supports propose a basis, never a preserver by themselves.
//! The hypergeometric row certificate identifies its recurrence for every index;
//! ordinary proof discovery (without catalog lifts) must prove every row real-rooted.
//! Finite supports use the fixed signed symbol and Hermite--Poulain instead.
//! Conversely, if every finite-support test polynomial is real-rooted, its
//! reversal scaled by `N` has coefficients `(-1)^k c_k (N)_(2k)/N^(2k)`.
//! These polynomials converge coefficientwise to `C(-w^2)`. Real-rootedness
//! is closed in this fixed finite degree, and the limit has constant term
//! one. Thus the signed-symbol criterion is necessary as well as sufficient.

use super::map::BasisMap;
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IndexSignDecision, IntegerRange,
};
use crate::recurrence::*;
use crate::{vec_poly, BigRational};
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::{Mutex, OnceLock};

type Q = BigRational;
const RECOGNITION_ROWS: usize = 18;
const RECOVERY_FINGERPRINT_ROWS: usize = 6;
const MAX_PARENT_SHIFT: usize = 4;

/// Exact description of `c`, with `c_0 = 1`.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum EvenDerivativeCoefficients {
    /// `c_(k+1) / c_k = numerator(k) / denominator(k)`, ascending coefficients.
    /// This route requires a strictly positive ratio on every integer `k >= 0`;
    /// terminating sequences can be supplied with `Finite`.
    RationalRatio {
        /// Numerator, ascending in `k`.
        numerator: Vec<String>,
        /// Denominator, ascending in `k`.
        denominator: Vec<String>,
    },
    /// A finite list, followed by zeros (at most eight entries).
    Finite {
        /// `c_0, ..., c_M`, exact rational strings, with `c_0 = 1`.
        coefficients: Vec<String>,
    },
}

/// Replayable all-degree hypotheses of the even-derivative theorem.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EvenDerivativeCertificate {
    /// The coefficient sequence.
    pub coefficients: EvenDerivativeCoefficients,
    /// All-integer positivity decision for a rational ratio.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ratio_sign: Option<IndexSignDecision>,
    /// Whether all coefficients are nonnegative (exact).
    pub nonnegative_coefficients: bool,
    /// Hypergeometric coefficient identity for all basis rows.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub triangle: Option<Box<HypergeometricRecordCertificate>>,
    /// Discovery on the canonical basis recurrence, with catalog lifts disabled.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_proof: Option<Box<RecurrenceProofDiscoveryReport>>,
    /// Even and odd reversed sections, when the direct row proof does not close.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub row_sections: Vec<EvenDerivativeSectionCertificate>,
    /// `C(-w^2)`, ascending, for finite support.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signed_symbol: Option<Vec<String>>,
    /// Every canonical basis row is PF, proved for all indices. False means
    /// not proved; a counterexample belongs to the enclosing theorem report.
    pub all_rows_pf: bool,
    /// Equality of the supplied basis recurrence with the canonical one.
    pub basis_identity: Box<RecurrenceEqualityCertificate>,
}

/// A proof for `J_m(x)=sum_k c_k (2m+parity)_(2k) x^k`.
/// Reversal and multiplication by a power of `x` recover `H_(2m+parity)`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EvenDerivativeSectionCertificate {
    /// `0` for even degrees, `1` for odd degrees.
    pub parity: usize,
    /// Exact coefficient identity for every section row.
    pub triangle: HypergeometricRecordCertificate,
    /// All-index real-rootedness proof, with catalog lifts disabled.
    pub row_proof: RecurrenceProofDiscoveryReport,
}

impl EvenDerivativeCertificate {
    /// Replay all coefficient, row and basis-identity hypotheses by exact
    /// recomputation from the recorded basis and coefficient description.
    /// A valid replay can also record failed hypotheses; inspect `all_rows_pf`
    /// and `nonnegative_coefficients` before using the preserver theorem.
    ///
    /// # Errors
    /// Malformed records, exhausted proof budgets, or a changed certificate.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let (rule, first, seeds) = self.basis_identity.alternative.to_recurrence_parts()?;
        let map = BasisMap::new("even derivative replay", rule, first, seeds)?;
        if certificate(&map, &self.coefficients)? != *self {
            return Err(RecurrenceError::VerificationFailed(
                "even-derivative hypotheses differ from their replay".into(),
            ));
        }
        Ok(())
    }
}

fn parse(values: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    values
        .iter()
        .map(|s| parse_rational_coeff(s).map_err(RecurrenceError::from))
        .collect()
}
fn text(values: &[Q]) -> Vec<String> {
    values.iter().map(format_rational_coeff).collect()
}
fn q(value: usize) -> Q {
    Q::from_integer(value.into())
}

impl EvenDerivativeCoefficients {
    fn values(&self, count: usize) -> Result<Vec<Q>, RecurrenceError> {
        match self {
            Self::Finite { coefficients } => {
                let mut values = parse(coefficients)?;
                if values.is_empty() || values.len() > 8 || values[0] != Q::one() {
                    return Err(RecurrenceError::InvalidInput(
                        "finite c must have c_0=1 and at most eight entries".into(),
                    ));
                }
                values.resize(count, Q::zero());
                values.truncate(count);
                Ok(values)
            }
            Self::RationalRatio {
                numerator,
                denominator,
            } => {
                let (a, b) = (parse(numerator)?, parse(denominator)?);
                let mut values = vec![Q::one()];
                for k in 0..count.saturating_sub(1) {
                    let divisor = vec_poly::evaluate(&b, &q(k));
                    if divisor.is_zero() {
                        return Err(RecurrenceError::InvalidInput(
                            "coefficient ratio has a pole".into(),
                        ));
                    }
                    values.push(values[k].clone() * vec_poly::evaluate(&a, &q(k)) / divisor);
                }
                values.truncate(count);
                Ok(values)
            }
        }
    }

    /// Generate `H_N = sum_k c_k (N)_(2k) x^(N-k)`, `0 <= N < count`.
    /// Coefficients are ascending and use exact rational arithmetic.
    ///
    /// # Errors
    /// A malformed coefficient, a pole, or more than 512 rows.
    pub fn basis_rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        if count > super::MAX_BASIS_MAP_ROWS {
            return Err(RecurrenceError::LimitExceeded(
                "too many even-derivative rows".into(),
            ));
        }
        let c = self.values(count / 2 + 1)?;
        Ok((0..count)
            .map(|n| {
                let mut row = vec![Q::zero(); n + 1];
                let mut falling = Q::one();
                for k in 0..=n / 2 {
                    if k > 0 {
                        falling *= q(n - 2 * k + 2) * q(n - 2 * k + 1);
                    }
                    row[n - k] = &c[k] * &falling;
                }
                row
            })
            .collect())
    }
}

fn polynomial_expression(values: &[String], variable: &str) -> Result<String, RecurrenceError> {
    let factors = crate::factorization::irreducible_factorization_rational_coeffs_with_budget(
        &parse(values)?,
        256,
    )
    .map_err(|e| RecurrenceError::Unsupported(e.to_string()))?;
    let mut expression = format!("({})", format_rational_coeff(&factors.constant));
    for (factor, multiplicity) in factors.factors {
        if factor.len() != 2 {
            return Err(RecurrenceError::Unsupported(
                "the triangle API requires rational linear factors in the ratio".into(),
            ));
        }
        expression.push_str(&format!(
            "*({}*({variable})+({}))^{multiplicity}",
            factor[1], factor[0]
        ));
    }
    Ok(expression)
}
fn triangle(
    coefficients: &EvenDerivativeCoefficients,
) -> Result<HypergeometricTriangle, RecurrenceError> {
    let EvenDerivativeCoefficients::RationalRatio {
        numerator,
        denominator,
    } = coefficients
    else {
        return Err(RecurrenceError::Unsupported(
            "finite coefficients have no hypergeometric triangle".into(),
        ));
    };
    Ok(HypergeometricTriangle::single(
        HypergeometricTermDefinition::new(
            1,
            &format!(
                "({})/(({})*(2*k-n+1)*(2*k-n+2))",
                polynomial_expression(denominator, "n-k-1")?,
                polynomial_expression(numerator, "n-k-1")?
            ),
            &format!(
                "({})*(n+1)*(2*k-n)/({})",
                polynomial_expression(numerator, "n-k")?,
                polynomial_expression(denominator, "n-k")?
            ),
            "n/2",
            "n",
            "1",
        ),
    ))
}

// For finite support M, (S-x)^(2M+1) annihilates H_N: after factoring
// x^N, every coefficient is a polynomial in N of degree at most 2M.
fn finite_map(coefficients: &EvenDerivativeCoefficients) -> Result<BasisMap, RecurrenceError> {
    let EvenDerivativeCoefficients::Finite { coefficients: c } = coefficients else {
        unreachable!()
    };
    if c.is_empty() || c.len() > 8 || parse(c)?[0] != Q::one() {
        return Err(RecurrenceError::InvalidInput(
            "finite c must have c_0=1 and at most eight entries".into(),
        ));
    }
    let order = 2 * (c.len() - 1) + 1;
    let terms = (1..=order)
        .map(|j| {
            let mut p = vec![Q::zero(); j + 1];
            p[j] = Q::from_integer(crate::numeric::binomial_bigint(order, j));
            if j % 2 == 0 {
                p[j] = -p[j].clone();
            }
            RecurrenceTerm::new(j, 0, BivarPoly::from_coefficients_normalized(vec![p]))
        })
        .collect();
    BasisMap::new(
        "even derivative (finite support)",
        Recurrence {
            terms,
            denominator: None,
            inhomogeneous: None,
        },
        0,
        coefficients.basis_rows(order)?,
    )
}

fn canonical_map(coefficients: &EvenDerivativeCoefficients) -> Result<BasisMap, RecurrenceError> {
    static CACHE: OnceLock<Mutex<BTreeMap<EvenDerivativeCoefficients, BasisMap>>> = OnceLock::new();
    let cache = CACHE.get_or_init(Mutex::default);
    if let Some(map) = cache.lock().expect("cache").get(coefficients).cloned() {
        return Ok(map);
    }
    let map = match coefficients {
        EvenDerivativeCoefficients::Finite { .. } => finite_map(coefficients)?,
        EvenDerivativeCoefficients::RationalRatio {
            numerator,
            denominator,
        } => {
            let sign = decide_index_rational_function_sign(
                &parse(numerator)?,
                &parse(denominator)?,
                &IntegerRange::starting_at(0),
                IndexSignCondition::Positive,
            );
            if !matches!(sign, IndexSignDecision::Holds { .. }) {
                return Err(RecurrenceError::Unsupported("ratio not proved strictly positive for every k>=0; use finite support for a terminating sequence".into()));
            }
            let triangle = triangle(coefficients)?;
            let numerator = vec_poly::trim(&parse(numerator)?);
            let denominator = vec_poly::trim(&parse(denominator)?);
            if numerator.len() == 1 && denominator.len() == 2 && denominator[0] == denominator[1] {
                // The Hermite row recurrence is a cheap proposal. Certify it
                // against the triangle before using it for any proof.
                let scale = &numerator[0] / &denominator[0];
                let two_scale = q(2) * scale;
                let rule = Recurrence {
                    terms: vec![
                        RecurrenceTerm::new(
                            1,
                            0,
                            BivarPoly::from_coefficients_normalized(vec![vec![
                                Q::zero(),
                                Q::one(),
                            ]]),
                        ),
                        RecurrenceTerm::new(
                            2,
                            0,
                            BivarPoly::from_coefficients_normalized(vec![
                                vec![Q::zero(), -&two_scale],
                                vec![Q::zero(), two_scale],
                            ]),
                        ),
                    ],
                    denominator: None,
                    inhomogeneous: None,
                };
                let map = BasisMap::new(
                    "even derivative (rational ratio)",
                    rule,
                    0,
                    coefficients.basis_rows(2)?,
                )?;
                certify_hypergeometric_row_record(&triangle, &map.tail(1)?.recurrence_json())?;
                map
            } else {
                let derived = derive_hypergeometric_row_recurrence(
                    &triangle,
                    &HypergeometricRowRecurrenceBudget {
                        max_order: 2,
                        max_variable_power: 2,
                        max_euler_power: 2,
                        max_index_degree: 3,
                        max_unknowns: 120,
                        max_checked_row: 64,
                    },
                )?;
                let (rule, first, seeds) = derived.record.to_recurrence_parts()?;
                debug_assert_eq!(first, 1);
                let mut full_seeds = vec![vec![Q::one()]];
                full_seeds.extend(seeds);
                BasisMap::new("even derivative (rational ratio)", rule, 0, full_seeds)?
            }
        }
    };
    cache
        .lock()
        .expect("cache")
        .insert(coefficients.clone(), map.clone());
    Ok(map)
}

impl BasisMap {
    /// Construct the basis `u^N -> sum_k c_k (N)_(2k) x^(N-k)` from exact
    /// coefficients. Construction certifies the row recurrence; preservation
    /// still requires [`check_basis_map_preserves_pf`].
    ///
    /// # Errors
    /// Invalid coefficients, a ratio not proved positive, or a row recurrence
    /// outside the bounded hypergeometric search.
    pub fn even_derivative(
        coefficients: &EvenDerivativeCoefficients,
    ) -> Result<Self, RecurrenceError> {
        canonical_map(coefficients)
    }
}

/// Fit a rational ratio from coefficient values (proposal only).
/// Numerator/denominator degrees are at most two, and the denominator is
/// normalized to constant term one. At least two extra equations are required.
fn fit_ratio(c: &[Q]) -> Option<EvenDerivativeCoefficients> {
    if c.iter().any(|x| !x.is_positive()) {
        return None;
    }
    for total in 1..=4 {
        for a_degree in 0..=2 {
            if total < a_degree {
                continue;
            }
            let b_degree = total - a_degree;
            if b_degree > 2 || c.len() < total + 4 {
                continue;
            }
            let mut matrix = Vec::new();
            let mut rhs = Vec::new();
            for k in 0..c.len() - 1 {
                let x = q(k);
                let r = &c[k + 1] / &c[k];
                let mut row = (0..=a_degree).map(|i| x.pow(i as i32)).collect::<Vec<_>>();
                row.extend((1..=b_degree).map(|i| -&r * x.pow(i as i32)));
                matrix.push(row);
                rhs.push(r);
            }
            let Some(solution) = crate::linalg::solve_linear_system(&matrix, &rhs) else {
                continue;
            };
            let a = vec_poly::trim(&solution[..=a_degree]);
            let mut b = vec![Q::one()];
            b.extend_from_slice(&solution[a_degree + 1..]);
            let candidate = EvenDerivativeCoefficients::RationalRatio {
                numerator: text(&a),
                denominator: text(&vec_poly::trim(&b)),
            };
            if candidate.values(c.len()).ok().as_ref() == Some(&c.to_vec()) {
                return Some(candidate);
            }
        }
    }
    None
}

/// Recover coefficients from one image/parent pair. Row scalars are not
/// silently discarded: an all-index identity must use the actual rows.
fn recover_row(image: &[Q], parent: &[Q]) -> Option<Vec<Q>> {
    let n = parent.len().checked_sub(1)?;
    if image.len() != parent.len() || image[n] != parent[n] || parent[n].is_zero() {
        return None;
    }
    let mut c = vec![Q::one()];
    for k in 1..=n / 2 {
        let mut value = image[n - k].clone();
        for (j, cj) in c.iter().enumerate() {
            let m = n - k + j;
            if m >= 2 * j {
                value -= cj
                    * &parent[m]
                    * Q::from_integer(crate::numeric::falling_factorial_bigint(&m.into(), 2 * j));
            }
        }
        c.push(
            value
                / (&parent[n]
                    * Q::from_integer(crate::numeric::falling_factorial_bigint(&n.into(), 2 * k))),
        );
    }
    // Check the whole row, including coefficients below the recovery band.
    let mut actual = vec![Q::zero(); n + 1];
    for (m, pm) in parent.iter().enumerate() {
        for k in 0..=m / 2 {
            actual[m - k] += pm
                * &c[k]
                * Q::from_integer(crate::numeric::falling_factorial_bigint(&m.into(), 2 * k));
        }
    }
    (actual == image).then_some(c)
}

fn recover_values(image: &[Vec<Q>], parent: &[Vec<Q>]) -> Option<Vec<Q>> {
    if image.len() != parent.len() {
        return None;
    }
    let mut coefficients = vec![Q::one()];
    for (q, p) in image.iter().zip(parent) {
        let c = recover_row(q, p)?;
        if c.iter().zip(&coefficients).any(|(a, b)| a != b) {
            return None;
        }
        if c.len() > coefficients.len() {
            coefficients = c;
        }
    }
    Some(coefficients)
}

fn recover_rows(image: &[Vec<Q>], parent: &[Vec<Q>]) -> Option<EvenDerivativeCoefficients> {
    let coefficients = recover_values(image, parent)?;
    if coefficients.len() < 4 || coefficients[1..].iter().all(Zero::is_zero) {
        return None;
    }
    if let Some(last) = coefficients.iter().rposition(|c| !c.is_zero()) {
        if last + 1 < coefficients.len() {
            return Some(EvenDerivativeCoefficients::Finite {
                coefficients: text(&coefficients[..=last]),
            });
        }
    }
    fit_ratio(&coefficients)
}

/// Recognize the row shape, then certify its all-index identity and hypotheses.
pub(super) fn coefficient_candidate(
    map: &BasisMap,
) -> Result<Option<EvenDerivativeCoefficients>, RecurrenceError> {
    let prefix = map.basis_rows(3)?;
    if prefix[0] != vec![Q::one()] || prefix[1] != vec![Q::zero(), Q::one()] {
        return Ok(None);
    }
    let rows = map.basis_rows(RECOGNITION_ROWS)?;
    if rows
        .iter()
        .enumerate()
        .any(|(n, r)| r.len() != n + 1 || !r[n].is_one())
    {
        return Ok(None);
    }
    let monomials = (0..RECOGNITION_ROWS)
        .map(|n| {
            let mut p = vec![Q::zero(); n + 1];
            p[n] = Q::one();
            p
        })
        .collect::<Vec<_>>();
    Ok(recover_rows(&rows, &monomials))
}

pub(super) fn recognize(
    map: &BasisMap,
) -> Result<Option<EvenDerivativeCertificate>, RecurrenceError> {
    let Some(coefficients) = coefficient_candidate(map)? else {
        return Ok(None);
    };
    certificate(map, &coefficients).map(Some)
}

fn prove_record(
    record: &RecurrenceJson,
) -> Result<RecurrenceProofDiscoveryReport, RecurrenceError> {
    discover_recurrence_proofs_with_stages(
        record,
        &RealRootednessTheoremOptions::default(),
        &AlternativeRepresentationOptions {
            max_candidates: 0,
            max_alternatives: 0,
            ..AlternativeRepresentationOptions::default()
        },
        &CatalogLiftOptions::disabled(),
    )
}

fn section_proofs(
    coefficients: &EvenDerivativeCoefficients,
) -> Result<Vec<EvenDerivativeSectionCertificate>, RecurrenceError> {
    let EvenDerivativeCoefficients::RationalRatio {
        numerator,
        denominator,
    } = coefficients
    else {
        return Ok(Vec::new());
    };
    let mut sections = Vec::new();
    for parity in 0..=1 {
        let triangle = HypergeometricTriangle::single(HypergeometricTermDefinition::new(
            0,
            &format!(
                "({})*(2*n+{parity}-2*k)*(2*n+{parity}-2*k-1)/({})",
                polynomial_expression(numerator, "k")?,
                polynomial_expression(denominator, "k")?
            ),
            &format!(
                "(2*n+{parity}+1)*(2*n+{parity}+2)/((2*n+{parity}-2*k+1)*(2*n+{parity}-2*k+2))"
            ),
            "0",
            "n",
            "1",
        ));
        let derived = derive_hypergeometric_row_recurrence(
            &triangle,
            &HypergeometricRowRecurrenceBudget {
                max_order: 2,
                max_variable_power: 2,
                max_euler_power: 2,
                max_index_degree: 3,
                max_unknowns: 120,
                max_checked_row: 64,
            },
        )?;
        let identity = certify_hypergeometric_row_record(&triangle, &derived.record)?;
        let row_proof = prove_record(&derived.record)?;
        sections.push(EvenDerivativeSectionCertificate {
            parity,
            triangle: identity,
            row_proof,
        });
    }
    Ok(sections)
}

pub(super) fn certificate(
    map: &BasisMap,
    coefficients: &EvenDerivativeCoefficients,
) -> Result<EvenDerivativeCertificate, RecurrenceError> {
    let canonical = canonical_map(coefficients)?;
    let shifted = canonical.rule().shift_index(-(map.first_index() as i64));
    let record = RecurrenceJson::from_recurrence_rational(
        &shifted,
        map.first_index(),
        canonical.seeds(),
        None,
    );
    let basis_identity = Box::new(certify_recurrence_equality(
        &record,
        &map.recurrence_json(),
    )?);
    let (ratio_sign, nonnegative_coefficients, tri, row_proof, signed_symbol, all_rows_pf) =
        match coefficients {
            EvenDerivativeCoefficients::Finite { coefficients: c } => {
                let c = parse(c)?;
                let nonnegative = c.iter().all(|x| !x.is_negative());
                let mut symbol = vec![Q::zero(); 2 * c.len() - 1];
                for (k, ck) in c.iter().enumerate() {
                    symbol[2 * k] = if k % 2 == 0 { ck.clone() } else { -ck };
                }
                let pf = nonnegative
                    && super::AlgebraicInterval::whole_line().contains_all_roots_of(&symbol);
                (None, nonnegative, None, None, Some(text(&symbol)), pf)
            }
            EvenDerivativeCoefficients::RationalRatio {
                numerator,
                denominator,
            } => {
                let sign = decide_index_rational_function_sign(
                    &parse(numerator)?,
                    &parse(denominator)?,
                    &IntegerRange::starting_at(0),
                    IndexSignCondition::Positive,
                );
                let positive = matches!(sign, IndexSignDecision::Holds { .. });
                let tri = certify_hypergeometric_row_record(
                    &triangle(coefficients)?,
                    &canonical.tail(1)?.recurrence_json(),
                )?;
                let proof = prove_record(&canonical.recurrence_json())?;
                let pf =
                    positive && proof.status == "proved" && proof.real_rooted_from_index == Some(0);
                (
                    Some(sign),
                    positive,
                    Some(Box::new(tri)),
                    Some(Box::new(proof)),
                    None,
                    pf,
                )
            }
        };
    let row_sections = if !all_rows_pf && tri.is_some() {
        match section_proofs(coefficients) {
            Ok(sections) => sections,
            Err(RecurrenceError::BudgetExhausted(_)) | Err(RecurrenceError::Unsupported(_)) => {
                Vec::new()
            }
            Err(error) => return Err(error),
        }
    } else {
        Vec::new()
    };
    let all_rows_pf = all_rows_pf
        || (nonnegative_coefficients
            && row_sections.len() == 2
            && row_sections.iter().all(|s| {
                s.row_proof.status == "proved" && s.row_proof.real_rooted_from_index == Some(0)
            }));
    Ok(EvenDerivativeCertificate {
        coefficients: coefficients.clone(),
        ratio_sign,
        nonnegative_coefficients,
        triangle: tri,
        row_proof,
        row_sections,
        signed_symbol,
        all_rows_pf,
        basis_identity,
    })
}

/// An exponential Riordan transport. For the Touchard coefficient basis,
/// R(z,u)=exp(u*(exp(z)-1)); normal ordering gives
/// H(z,x)=exp(x*(f+a*f^2)), f=exp(z)-1. With b=2a-1 and
/// K=S-x-2xD, K H=b*x*f*H. Since f'=1+f and f+a*f^2=D H/H,
/// [a*K^2+(b*x+a)*K-b^2*x^2*D-a*b*x] H=0.
/// The parent and basis identities are certified before using this formula.
pub(super) fn exponential_riordan_transport(
    coefficient: &crate::recurrence::lift::hadamard_factor::prepare::PreparedLiftRule,
    map: &BasisMap,
) -> Result<Option<super::transport::Transport>, RecurrenceError> {
    use crate::recurrence::implicit::algebra::{Biv, DOp};
    use crate::recurrence::lift::section_operators::OperatorMatrix;
    let prefix = map.basis_rows(3)?;
    if prefix[0] != vec![Q::one()] || prefix[1] != vec![Q::zero(), Q::one()] {
        return Ok(None);
    }
    let rows = map.basis_rows(RECOGNITION_ROWS)?;
    let a = rows[2].get(1).cloned().unwrap_or_else(Q::zero) / q(2);
    if !a.is_positive() {
        return Ok(None);
    }
    let coefficients = EvenDerivativeCoefficients::RationalRatio {
        numerator: text(std::slice::from_ref(&a)),
        denominator: vec!["1".into(), "1".into()],
    };
    if coefficients.basis_rows(RECOGNITION_ROWS)? != rows {
        return Ok(None);
    }
    let canonical = canonical_map(&coefficients)?;
    let record = RecurrenceJson::from_recurrence_rational(
        &canonical.rule().shift_index(-(map.first_index() as i64)),
        map.first_index(),
        canonical.seeds(),
        None,
    );
    if !certify_recurrence_equality(&record, &map.recurrence_json())?.is_certified() {
        return Ok(None);
    }
    let touchard = Recurrence {
        terms: vec![
            RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![Q::zero(), Q::one()]]),
            ),
            RecurrenceTerm::new(
                1,
                1,
                BivarPoly::from_coefficients_normalized(vec![vec![Q::zero(), Q::one()]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let parent = RecurrenceJson::from_recurrence_rational(
        &touchard,
        coefficient.first_index,
        &[vec![Q::one()]],
        None,
    );
    let actual = RecurrenceJson::from_recurrence_rational(
        &coefficient.rule,
        coefficient.first_index,
        &coefficient.seeds,
        None,
    );
    // A cheap fingerprint is followed by the all-index equality engine.
    if coefficient.rows(6)?
        != touchard
            .generate_rows_rational(&[vec![Q::one()]], coefficient.first_index, 6)
            .map_err(|e| RecurrenceError::InvalidInput(e.to_string()))?
    {
        return Ok(None);
    }
    if !certify_recurrence_equality(&parent, &actual)?.is_certified() {
        return Ok(None);
    }
    let x = Biv::new(vec![vec![Q::zero(), Q::one()]]);
    let mut k = DOp::single(1, 0, Biv::one());
    k.add_term((0, 0), x.scale(&(-Q::one())));
    k.add_term((0, 1), x.scale(&(-q(2))));
    let b = q(2) * &a - Q::one();
    let affine = Biv::new(vec![vec![a.clone(), b.clone()]]);
    let relation = k
        .compose(&k)
        .scale(&a)
        .add(&DOp::single(0, 0, affine).compose(&k))
        .subtract(&DOp::single(0, 1, x.multiply(&x).scale(&(&b * &b))))
        .subtract(&DOp::single(0, 0, x.scale(&(&a * &b))));
    Ok(Some(super::transport::Transport {
        system: OperatorMatrix {
            entries: vec![vec![relation]],
        },
        order: 2,
        leader: Biv::constant(a),
        valid_from_base: coefficient.first_index,
        identifying_rows: 2,
        kind: super::transport::TransportKind::Riordan,
    }))
}

/// Recovered candidates from catalog parents, with bounded row and ratio
/// searches. Neither a finite match nor the catalog's flags is a proof.
pub(super) fn recovered_candidates(
    rows: &[Vec<Q>],
    first: usize,
    limit: usize,
) -> Vec<BasisMapLiftCandidate> {
    struct Parent {
        first: usize,
        rows: Vec<Vec<Q>>,
        seeds: Vec<Vec<Q>>,
        rule: Recurrence,
    }
    static PARENTS: OnceLock<Vec<Parent>> = OnceLock::new();
    let parents = PARENTS.get_or_init(|| {
        crate::oeis::catalog()
            .iter()
            .filter_map(|entry| {
                let (rule, seeds, first) = entry.recurrence_export_parts().ok()?;
                let first = usize::try_from(first).ok()?;
                let rows = rule
                    .generate_rows_rational(
                        &seeds,
                        first,
                        RECOVERY_FINGERPRINT_ROWS + MAX_PARENT_SHIFT,
                    )
                    .ok()?;
                Some(Parent {
                    first,
                    rows: rows.iter().map(|r| vec_poly::trim(r)).collect(),
                    seeds,
                    rule,
                })
            })
            .collect()
    });
    let mut candidates = Vec::new();
    let fingerprint_count = rows.len().min(RECOVERY_FINGERPRINT_ROWS);
    for parent in parents {
        for shift in 0..=MAX_PARENT_SHIFT {
            let Some(prefix) = parent.rows.get(shift..shift + fingerprint_count) else {
                continue;
            };
            // A full identity must pass these exact prefix equations. Generate
            // the larger rows only for survivors, keeping the candidate order.
            if rows[..fingerprint_count]
                .iter()
                .zip(prefix)
                .any(|(a, b)| a.len() != b.len() || a.last() != b.last())
                || recover_values(&rows[..fingerprint_count], prefix).is_none()
            {
                continue;
            }
            let Ok(full) =
                parent
                    .rule
                    .generate_rows_rational(&parent.seeds, parent.first, shift + rows.len())
            else {
                continue;
            };
            let full: Vec<_> = full.iter().map(|row| vec_poly::trim(row)).collect();
            let Some(window) = full.get(shift..shift + rows.len()) else {
                continue;
            };
            let Some(coefficients) = recover_rows(rows, window) else {
                continue;
            };
            if coefficients.values(8).is_err() {
                continue;
            }
            let Ok(map) = canonical_map(&coefficients) else {
                continue;
            };
            let rule = parent
                .rule
                .shift_index((parent.first + shift) as i64 - first as i64);
            let Ok(prepared) = crate::recurrence::lift::hadamard_factor::prepare::prepare_lift_rule(
                &rule,
                window,
                first,
                "even-derivative coefficient family",
            ) else {
                continue;
            };
            let candidate = BasisMapLiftCandidate {
                basis: map,
                coefficient_family: RecurrenceJson::from_recurrence_rational(
                    &prepared.rule,
                    first,
                    &prepared.seeds,
                    None,
                ),
            };
            if !candidates.contains(&candidate) {
                candidates.push(candidate);
            }
            if candidates.len() >= limit {
                return candidates;
            }
        }
    }
    candidates
}
