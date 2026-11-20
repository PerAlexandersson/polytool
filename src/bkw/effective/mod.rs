//! Effective dominant-pair certificates for rational generating functions.
//!
//! Unlike the limit-set scout, this module checks signs over entire parameter
//! intervals, a phase lift, a degree budget and every row below an exact cutoff.
//! [`certify_effective_bkw`] accepts ascending rational-string coefficients;
//! [`EffectiveBkwReport::verify`] replays all arithmetic from the embedded request.
//! No floating-point computations or caller attestations enter a proof.

use crate::recurrence::{parse_rational_coeff, RationalGeneratingFunction, RecurrenceError};
use crate::{vec_poly, BigRational as Q, Polynomial};
use num_traits::{One, Signed, Zero};
use ore_algebra::BivariatePolynomial as B;
use serde::{Deserialize, Serialize};

mod algebra;
mod compact;
mod cubic;
mod signs;
#[cfg(test)]
mod tests;
use algebra::*;
pub use compact::uniform_tail_cutoff;
use signs::*;
pub use signs::{EffectiveBkwBox, EffectiveBkwCellSample, EffectiveBkwSignCertificate};
type P = Polynomial<Q>;

/// Version of the effective-BKW request and replay protocol.
pub const EFFECTIVE_BKW_SCHEMA: &str = "polytool.effective-bkw.v1";

/// A real endpoint, represented exactly, including either infinity.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum EffectiveBkwEndpoint {
    /// Negative infinity.
    NegativeInfinity,
    /// Positive infinity.
    PositiveInfinity,
    /// A rational endpoint.
    Rational {
        /// Canonical rational string.
        value: String,
    },
    /// The unique root of a polynomial in a closed rational enclosure.
    Algebraic {
        /// Ascending rational coefficients of the defining polynomial.
        polynomial: Vec<String>,
        /// Rational isolating bounds; equal bounds denote an exact root.
        enclosure: [String; 2],
    },
}

/// An open real parameter interval; algebraic endpoints are supported.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwInterval {
    /// Left endpoint, excluded.
    pub lower: EffectiveBkwEndpoint,
    /// Right endpoint, excluded.
    pub upper: EffectiveBkwEndpoint,
}

/// The count `slope*n + offsets[n % offsets.len()]`, for the local GF index n.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwZeroCount {
    /// Linear growth per local row index.
    pub slope: usize,
    /// Periodic offsets. A proof checks the complete period symbolically.
    pub offsets: Vec<i64>,
}

/// Uniform compact-interval tail constants (all exact rational strings).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwUniformConstants {
    /// Upper bound on the remainder-root modulus ratio, in `(0,1)`.
    pub modulus_ratio: String,
    /// Positive lower bound for `|N/H|` at the reciprocal pair root.
    pub amplitude_lower_bound: String,
    /// Bounds for `|N_j|/|lambda|^j`, including j=0.
    pub numerator_weights: Vec<String>,
    /// Optional bounds `|H_j|/|lambda|^j <= gamma_j`, j=1,...,d-2.
    /// For a supplied higher-order factorization these must satisfy
    /// `sum gamma_j/q^j < 1`, a uniform Cauchy root bound.
    #[serde(default)]
    pub cofactor_weights: Option<Vec<String>>,
}

/// Supplied rational-polynomial pair factorization for higher-order inputs.
/// The identity is checked over `Q[t,z]`; this is data, not an attestation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwPairFactorization {
    /// Ascending t coefficients of Re(lambda).
    pub real_part: Vec<String>,
    /// Ascending t coefficients of |lambda|^2.
    pub norm_squared: Vec<String>,
    /// Ascending z coefficients of H(t,z), each ascending in t, with H_0=1.
    pub cofactor: Vec<Vec<String>>,
}

/// Exact input to the effective-BKW checker.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwRequest {
    /// Seed-dependent rational GF; z^n denotes local index n.
    pub generating_function: RationalGeneratingFunction,
    /// Requested pair interval; omitted selects a discriminant component.
    #[serde(default)]
    pub interval: Option<EffectiveBkwInterval>,
    /// Requested total real-zero count; omitted uses the certified degree budget.
    #[serde(default)]
    pub target_zero_count: Option<EffectiveBkwZeroCount>,
    /// Compact template constants; omitted quartic bounds use an exact grid
    /// search. Required with a supplied higher-order factorization.
    #[serde(default)]
    pub uniform_constants: Option<EffectiveBkwUniformConstants>,
    /// Optional exact pair factorization; required beyond automatic order four.
    #[serde(default)]
    pub pair_factorization: Option<EffectiveBkwPairFactorization>,
}

/// A finite exact row check below the cutoff.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwFiniteRow {
    /// Absolute row index (GF first_index plus local index).
    pub index: usize,
    /// Ascending rational coefficients.
    pub coefficients: Vec<String>,
    /// Ordinary polynomial degree.
    pub degree: usize,
    /// Number of distinct real roots, equal to the degree for a simple row.
    pub distinct_real_roots: usize,
}

/// A rational expression in the branch coordinate and parameter t.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwBranchRational {
    /// Numerator coefficients `[branch_power][t_power]`, ascending.
    pub numerator: Vec<Vec<String>>,
    /// Nonzero denominator in the same coefficient convention.
    pub denominator: Vec<Vec<String>>,
}

/// Machine-readable spectral factorization feeding a future pair emitter.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwExpansionData {
    /// `Re(lambda)` in the algebraic branch coordinate.
    pub pair_real_part: EffectiveBkwBranchRational,
    /// `|lambda|^2` in the algebraic branch coordinate.
    pub pair_norm_squared: EffectiveBkwBranchRational,
    /// Ascending GF cofactor H in `D=(1-2*a*z+b*z^2)H`.
    pub cofactor: Vec<EffectiveBkwBranchRational>,
    /// Single real remainder eigenvalue; absent for the other suppliers.
    pub remainder_root: Option<EffectiveBkwBranchRational>,
    /// Its amplitude C in `rest_n=C*r^n`.
    pub remainder_amplitude: Option<EffectiveBkwBranchRational>,
}

/// Limits of the chosen cosine phase `n*arg(lambda)+arg(A)`, in units of pi.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwPhaseLimits {
    /// Left eigenvalue argument limit, 1 for the supported orientation.
    pub left_eigenvalue_argument: String,
    /// Right eigenvalue argument limit, 0 for the supported orientation.
    pub right_eigenvalue_argument: String,
    /// Limit of the continuous amplitude argument lift at the left edge.
    pub left_amplitude_argument: String,
    /// Limit of the same lift at the right edge.
    pub right_amplitude_argument: String,
}

/// Exact degree and leading-coefficient evidence for exterior signs.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwDegreeCertificate {
    /// All-index upper bound, with the same local indexing as the zero count.
    pub upper_bound: EffectiveBkwZeroCount,
    /// Reduced GF for `[t^n]P_n`, ascending rational coefficients in z.
    pub leading_numerator: Vec<String>,
    /// Denominator of that GF, normalized to constant coefficient one.
    pub leading_denominator: Vec<String>,
    /// The top coefficient is positive exactly in these residue classes
    /// modulo `upper_bound.offsets.len()`; elsewhere it is zero.
    pub positive_top_residues: Vec<usize>,
}

/// Exact coordinate transport applied before the displayed certificate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EffectiveBkwCoordinateTransport {
    /// The chart family is the input family.
    Identity,
    /// The chart family is `Q_n(u)=(-1)^n P_n(-u)`, with local GF index n.
    /// This reflects a positive-infinity edge to negative infinity and
    /// preserves degrees, simplicity and real-rootedness in both directions.
    ParameterReflectionWithAlternatingRows,
}

/// Certified location of all roots of the input family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum EffectiveBkwRootRegion {
    /// Every root is real, with no additional certified location.
    WholeRealLine,
    /// Every root lies in `(-infinity,0)`.
    NegativeRay,
    /// Every root lies in `(0,infinity)`.
    PositiveRay,
}

/// Exact expansion, domination, phase, exterior and finite-prefix evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwCertificate {
    /// Exact transformation from the request to the chart family.
    pub coordinate_transport: EffectiveBkwCoordinateTransport,
    /// Certified region in the input variable, after coordinate transport.
    pub root_region: EffectiveBkwRootRegion,
    /// Checked open pair interval in the chart coordinate.
    pub interval: EffectiveBkwInterval,
    /// Ascending characteristic coefficients in lambda, each ascending in t.
    pub characteristic: Vec<Vec<String>>,
    /// Algebraic branch equation, outer variable the real root (cubic) or
    /// pair product (quartic); a constant branch is used for order two.
    pub branch_equation: Vec<Vec<String>>,
    /// How the continuous branch is selected on the whole open interval.
    pub branch_selection: String,
    /// Exact expansion formulas in the branch coordinate.
    pub expansion: Vec<String>,
    /// Exact pair factors, independent of the displayed formulas.
    pub expansion_data: EffectiveBkwExpansionData,
    /// Full-interval resultant/sign certificates, with off-branch boxes.
    pub signs: Vec<EffectiveBkwSignCertificate>,
    /// Domination supplier, independent of phase/degree closure.
    pub domination: String,
    /// Certified local cutoff for strict domination and the phase count.
    pub cutoff: usize,
    /// Exact tail bound at cutoff; zero for the no-remainder supplier.
    pub tail_at_cutoff: String,
    /// Certified endpoint and argument-lift data, in multiples of pi.
    pub phase: Vec<String>,
    /// Exact phase-limit constants for later emission.
    pub phase_limits: EffectiveBkwPhaseLimits,
    /// Exterior signs and the all-index degree argument.
    pub exterior: Vec<String>,
    /// Structured degree and leading-sign evidence for exterior closure.
    pub degree_certificate: EffectiveBkwDegreeCertificate,
    /// Total distinct real-root count, closing the degree budget.
    pub zero_count: EffectiveBkwZeroCount,
    /// All local rows 0 <= n < cutoff, checked exactly.
    pub finite_rows: Vec<EffectiveBkwFiniteRow>,
    /// Lean API declarations this evidence is intended to feed; no Lean replay
    /// is claimed by this arithmetic certificate.
    pub lean_declarations: Vec<String>,
}

/// A successful all-index proof or a precise failed sufficient hypothesis.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwReport {
    /// Always [`EFFECTIVE_BKW_SCHEMA`].
    pub schema: String,
    /// Complete replay input.
    pub request: EffectiveBkwRequest,
    /// `proved` or `not_proved`; the latter is not a counterexample verdict.
    pub status: String,
    /// The proved conclusion or the precise failed/unsupported obligation.
    pub reason: String,
    /// Present exactly when every proof obligation passed.
    pub certificate: Option<EffectiveBkwCertificate>,
}

impl EffectiveBkwReport {
    /// Replay the exact computation, rejecting any modified input or witness.
    ///
    /// # Errors
    /// Returns [`RecurrenceError::VerificationFailed`] on any mismatch.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let expected = certify_effective_bkw(&self.request)?;
        if self != &expected {
            return Err(RecurrenceError::VerificationFailed(
                "effective-BKW replay mismatch".into(),
            ));
        }
        Ok(())
    }
}

/// Construct the seed-dependent GF and try the automatic pair suppliers.
///
/// The quartic supplier searches a bounded grid of simple rationals for
/// modulus, amplitude and numerator bounds, checking each over the interval.
/// Use [`certify_effective_bkw`] to supply other constants or a factorization.
/// Extra prescribed seeds are retained; an improper numerator is unsupported
/// rather than silently shifting the claimed starting index.
///
/// # Errors
/// Returns [`RecurrenceError`] when the recurrence cannot be converted to
/// an index-free rational generating function or is malformed.
pub fn certify_effective_bkw_recurrence(
    record: &crate::recurrence::RecurrenceJson,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let gf = crate::recurrence::rational_generating_function(&rule, first, &seeds)?;
    certify_effective_bkw(&EffectiveBkwRequest {
        generating_function: gf,
        interval: None,
        target_zero_count: None,
        uniform_constants: None,
        pair_factorization: None,
    })
}

fn failure(message: impl Into<String>) -> RecurrenceError {
    RecurrenceError::Unsupported(message.into())
}

// Keep chart selection separate from diagnostic text and the public error wire.
#[derive(Debug)]
enum ChartError {
    NoNegativeCubicCollision,
    Recurrence(RecurrenceError),
}
impl From<RecurrenceError> for ChartError {
    fn from(error: RecurrenceError) -> Self {
        Self::Recurrence(error)
    }
}
impl From<ChartError> for RecurrenceError {
    fn from(error: ChartError) -> Self {
        match error {
            ChartError::NoNegativeCubicCollision => failure("no negative cubic collision endpoint"),
            ChartError::Recurrence(error) => error,
        }
    }
}

fn parse_uniform_constants(
    constants: &EffectiveBkwUniformConstants,
    numerator_count: usize,
) -> Result<(Q, Q, Vec<Q>), RecurrenceError> {
    let ratio = rational(&constants.modulus_ratio)?;
    let amplitude = rational(&constants.amplitude_lower_bound)?;
    let weights = constants
        .numerator_weights
        .iter()
        .map(|s| rational(s))
        .collect::<Result<Vec<_>, _>>()?;
    if ratio <= Q::zero()
        || ratio >= Q::one()
        || amplitude <= Q::zero()
        || weights.len() != numerator_count
        || weights.iter().any(|w| w.is_negative())
    {
        return Err(RecurrenceError::InvalidInput(
            "invalid compact constants or numerator weight count".into(),
        ));
    }
    Ok((ratio, amplitude, weights))
}
fn q(n: i64) -> Q {
    Q::from_integer(n.into())
}
fn encode(p: &P) -> Vec<String> {
    p.coeffs().iter().map(ToString::to_string).collect()
}
fn decode(p: &[String]) -> Result<P, RecurrenceError> {
    Ok(P::new(
        p.iter()
            .map(|s| parse_rational_coeff(s).map_err(Into::into))
            .collect::<Result<_, RecurrenceError>>()?,
    ))
}
fn encode_b(p: &B) -> Vec<Vec<String>> {
    p.rows()
        .iter()
        .map(|r| encode(&P::new(r.clone())))
        .collect()
}
fn decode_b(p: &[Vec<String>]) -> Result<B, RecurrenceError> {
    Ok(B::new(
        p.iter()
            .map(|r| decode(r).map(|p| p.coeffs().to_vec()))
            .collect::<Result<_, _>>()?,
    ))
}

/// Check an effective-BKW proof, using exact arithmetic throughout.
///
/// Unsupported geometry and failed sufficient inequalities produce a
/// replayable `not_proved` report, never an all-index assertion. Malformed
/// rational strings or invalid input shapes return [`RecurrenceError`].
pub fn certify_effective_bkw(
    request: &EffectiveBkwRequest,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    let mut report = EffectiveBkwReport {
        schema: EFFECTIVE_BKW_SCHEMA.into(),
        request: request.clone(),
        status: "not_proved".into(),
        reason: String::new(),
        certificate: None,
    };
    match certify(request) {
        Ok(certificate) => {
            report.status = "proved".into();
            report.reason = "every row is nonzero and has only simple real zeros; exact phase count, degree budget and finite prefix checked".into();
            report.certificate = Some(certificate);
        }
        Err(RecurrenceError::Unsupported(reason)) => report.reason = reason,
        Err(e) => return Err(e),
    }
    Ok(report)
}

fn certify(request: &EffectiveBkwRequest) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    if request
        .interval
        .as_ref()
        .is_some_and(|i| matches!(i.upper, EffectiveBkwEndpoint::PositiveInfinity))
    {
        return reflected_certificate(request);
    }
    match certify_in_chart(request) {
        Err(ChartError::NoNegativeCubicCollision) if request.interval.is_none() => {
            reflected_certificate(request)
        }
        result => result.map_err(Into::into),
    }
}

fn reflected_certificate(
    request: &EffectiveBkwRequest,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    if request.pair_factorization.is_some() {
        return Err(failure(
            "reflection of supplied higher-order factors is not implemented",
        ));
    }
    let reflected = reflect_request(request)?;
    let mut certificate = certify_in_chart(&reflected).map_err(RecurrenceError::from)?;
    certificate.coordinate_transport =
        EffectiveBkwCoordinateTransport::ParameterReflectionWithAlternatingRows;
    certificate.root_region = match certificate.root_region {
        EffectiveBkwRootRegion::NegativeRay => EffectiveBkwRootRegion::PositiveRay,
        EffectiveBkwRootRegion::PositiveRay => EffectiveBkwRootRegion::NegativeRay,
        EffectiveBkwRootRegion::WholeRealLine => EffectiveBkwRootRegion::WholeRealLine,
    };
    certificate.exterior.push("exact transport Q_n(u)=(-1)^n P_n(-u), obtained by N(t,z),D(t,z) -> N(-u,-z),D(-u,-z); nonzero scalar and real reflection preserve every root and multiplicity".into());
    Ok(certificate)
}

fn reflect_request(request: &EffectiveBkwRequest) -> Result<EffectiveBkwRequest, RecurrenceError> {
    let mut reflected = request.clone();
    let transform = |rows: &[Vec<String>]| -> Result<Vec<Vec<String>>, RecurrenceError> {
        rows.iter()
            .enumerate()
            .map(|(j, row)| {
                row.iter()
                    .enumerate()
                    .map(|(k, s)| {
                        rational(s).map(|v| {
                            if (j + k) % 2 == 0 {
                                v.to_string()
                            } else {
                                (-v).to_string()
                            }
                        })
                    })
                    .collect()
            })
            .collect()
    };
    reflected.generating_function.numerator = transform(&request.generating_function.numerator)?;
    reflected.generating_function.denominator =
        transform(&request.generating_function.denominator)?;
    let endpoint = |e: &EffectiveBkwEndpoint| -> Result<EffectiveBkwEndpoint, RecurrenceError> {
        Ok(match e {
            EffectiveBkwEndpoint::NegativeInfinity => EffectiveBkwEndpoint::PositiveInfinity,
            EffectiveBkwEndpoint::PositiveInfinity => EffectiveBkwEndpoint::NegativeInfinity,
            EffectiveBkwEndpoint::Rational { value } => EffectiveBkwEndpoint::Rational {
                value: (-rational(value)?).to_string(),
            },
            EffectiveBkwEndpoint::Algebraic {
                polynomial,
                enclosure,
            } => EffectiveBkwEndpoint::Algebraic {
                polynomial: polynomial
                    .iter()
                    .enumerate()
                    .map(|(j, s)| {
                        rational(s).map(|v| {
                            if j % 2 == 0 {
                                v.to_string()
                            } else {
                                (-v).to_string()
                            }
                        })
                    })
                    .collect::<Result<_, _>>()?,
                enclosure: [
                    (-rational(&enclosure[1])?).to_string(),
                    (-rational(&enclosure[0])?).to_string(),
                ],
            },
        })
    };
    reflected.interval = request
        .interval
        .as_ref()
        .map(|i| {
            Ok::<_, RecurrenceError>(EffectiveBkwInterval {
                lower: endpoint(&i.upper)?,
                upper: endpoint(&i.lower)?,
            })
        })
        .transpose()?;
    Ok(reflected)
}

fn certify_in_chart(request: &EffectiveBkwRequest) -> Result<EffectiveBkwCertificate, ChartError> {
    let gf = &request.generating_function;
    if gf.schema != crate::recurrence::RATIONAL_GENERATING_FUNCTION_SCHEMA {
        return Err(
            RecurrenceError::InvalidInput("invalid generating-function schema".into()).into(),
        );
    }
    let numerator = decode_b(&gf.numerator)?;
    let denominator = decode_b(&gf.denominator)?;
    if let Some(constants) = &request.uniform_constants {
        parse_uniform_constants(constants, numerator.rows().len())?;
    }
    let d: Vec<P> = denominator
        .rows()
        .iter()
        .map(|r| P::new(r.clone()))
        .collect();
    if d.first() != Some(&P::one()) {
        return Err(RecurrenceError::InvalidInput("effective BKW requires D(t,0)=1".into()).into());
    }
    let order = d.len() - 1;
    if !(2..=8).contains(&order) || (order > 4 && request.pair_factorization.is_none()) {
        return Err(failure("automatic pair factorization supports orders two, three and four; orders up to eight accept a supplied rational-polynomial pair factorization").into());
    }
    if numerator.is_zero() || numerator.index_degree().unwrap_or(0) >= order {
        return Err(failure(
            "requires a nonzero proper numerator (polynomial parts are not yet supported)",
        )
        .into());
    }
    if d.iter().any(|p| p.degree().unwrap_or(0) > 8) {
        return Err(failure("parameter degree exceeds the phase-A elimination guard (8)").into());
    }
    let degree = degree_budget(&numerator, &denominator)?;
    if request
        .target_zero_count
        .as_ref()
        .is_some_and(|c| c != &degree.count)
    {
        return Err(failure(
            "requested zero count does not match the exact all-index degree budget",
        )
        .into());
    }
    let mut certificate = if order == 3 && request.pair_factorization.is_none() {
        cubic::certify(&numerator, &d, request, &degree)?
    } else {
        compact::certify(&numerator, &d, request, &degree)?
    };
    let mut rows: Vec<P> = Vec::new();
    let mut negative_prefix = certificate.root_region == EffectiveBkwRootRegion::NegativeRay;
    for n in 0..certificate.cutoff {
        let mut row = P::new(numerator.rows().get(n).cloned().unwrap_or_default());
        for j in 1..=order.min(n) {
            row = row - d[j].clone() * rows[n - j].clone();
        }
        let integer_row = vec_poly::clear_denominators(row.coeffs());
        let real_roots = crate::root_count::count_real_roots_prs_bigint_coeffs(&integer_row);
        let degree = row
            .degree()
            .ok_or_else(|| failure(format!("finite row {n} is zero")))?;
        if real_roots != degree {
            return Err(failure(format!("finite row {n} is not simple real-rooted")).into());
        }
        if negative_prefix {
            let constant = row.coeff(0);
            negative_prefix = if constant.is_positive() {
                vec_poly::has_nonnegative_coefficients(row.coeffs())
            } else if constant.is_negative() {
                vec_poly::has_nonnegative_coefficients(&vec_poly::neg(row.coeffs()))
            } else {
                false
            };
        }
        certificate.finite_rows.push(EffectiveBkwFiniteRow {
            index: gf
                .first_index
                .checked_add(n)
                .ok_or_else(|| failure("absolute index overflow"))?,
            coefficients: encode(&row),
            degree,
            distinct_real_roots: real_roots,
        });
        rows.push(row);
    }
    if !negative_prefix {
        certificate.root_region = EffectiveBkwRootRegion::WholeRealLine;
    }
    Ok(certificate)
}
