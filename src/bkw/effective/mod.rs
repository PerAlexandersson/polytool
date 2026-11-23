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
mod degree;
mod edges;
mod higher;
mod pointwise;
mod quadratic;
mod quartic_ray;
mod replay;
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

/// The count `(slope*n + offsets[n % offsets.len()])/divisor`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwZeroCount {
    /// Numerator of the linear growth per local row index, divided by `divisor`.
    pub slope: usize,
    /// Periodic offsets. A proof checks the complete period symbolically.
    pub offsets: Vec<i64>,
    /// Positive denominator of a quasi-linear count; omitted means one.
    #[serde(default = "one_usize", skip_serializing_if = "is_one")]
    pub divisor: usize,
}

fn one_usize() -> usize {
    1
}
fn is_one(value: &usize) -> bool {
    *value == 1
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
    /// Requested real-zero count law from the cutoff; omitted uses the
    /// certified degree budget, with earlier rows checked individually.
    #[serde(default)]
    pub target_zero_count: Option<EffectiveBkwZeroCount>,
    /// Compact template constants; omitted bounds use exact rational grids.
    #[serde(default)]
    pub uniform_constants: Option<EffectiveBkwUniformConstants>,
    /// Optional exact pair factorization, bypassing automatic discovery.
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
    /// Left eigenvalue argument limit in units of pi, possibly fractional.
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
    /// Degree law in local GF indices; a Newton law starts at `newton.start`,
    /// with every earlier degree in `newton.prefix_degrees`.
    pub upper_bound: EffectiveBkwZeroCount,
    /// Reduced leading-coefficient GF for Newton laws; legacy diagonal
    /// templates store the GF for `[t^n]P_n`. Ascending coefficients in z.
    pub leading_numerator: Vec<String>,
    /// Denominator of that GF, normalized to constant coefficient one.
    pub leading_denominator: Vec<String>,
    /// Positive leading-coefficient classes from the Newton law's start.
    /// For legacy diagonal templates, these are the positive diagonal top
    /// classes and the remaining classes have zero diagonal coefficient.
    pub positive_top_residues: Vec<usize>,
    /// Newton weight and exact coefficient layers for a quasi-linear law.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub newton: Option<EffectiveBkwNewtonDegree>,
}

/// Weighted generating-function proof of the exact eventual degree law.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwNewtonDegree {
    /// Substitution `t=u^(-weight_denominator), z=z*u^weight_numerator`.
    pub weight_numerator: usize,
    /// Denominator of the Newton slope, relatively prime to the numerator.
    pub weight_denominator: usize,
    /// Multiply the substituted numerator by this power of u.
    pub numerator_shift: usize,
    /// Local start of the exact degree law; earlier rows are checked finitely.
    pub start: usize,
    /// Coefficient of each u power as a reduced rational function in z.
    pub layers: Vec<crate::recurrence::RationalGeneratingSeries>,
    /// First eventual nonzero layer for each index residue. Earlier layers
    /// can have finite support, covered by `prefix_degrees`.
    pub first_layers: Vec<usize>,
    /// Local degrees before `start`, absent for zero rows.
    pub prefix_degrees: Vec<Option<usize>>,
    /// Positive support seeds modulo the smallest Newton-face step, grouped
    /// by degree-law residue. Adding that step preserves strict positivity.
    pub reachability_seeds: Vec<Vec<usize>>,
    /// Exact sufficient positivity argument used for the coefficient layers.
    pub positivity_supplier: String,
    /// Optional all-index signs of polynomial coefficient formulas for a
    /// binomial Newton face, in the residue index k.
    pub leading_sign_certificates: Vec<crate::integer_range_sign::IndexSignCertificate>,
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
    /// pair product (higher orders); rational pair charts use a constant branch.
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
    /// Strict domination holds for every local n at least this value.
    /// Omitted in older certificates, whose joint cutoff supplies the bound.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub domination_cutoff: Option<usize>,
    /// Certified local cutoff for strict domination and the phase count.
    pub cutoff: usize,
    /// Exact tail bound at cutoff; zero for the no-remainder supplier.
    pub tail_at_cutoff: String,
    /// Automatically found compact/Cauchy constants for higher-order pairs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub automatic_uniform_constants: Option<EffectiveBkwUniformConstants>,
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
    /// Original complete recurrence when a prefix was removed or its validation failed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_recurrence: Option<crate::recurrence::RecurrenceJson>,
    /// Exactly checked prefix omitted from the nonzero simple-root conclusion.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub skipped_prefix: Vec<EffectiveBkwPrefixRow>,
}

/// A removed zero or simple real-rooted row, bound to the original recurrence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwPrefixRow {
    /// Absolute index in the original record.
    pub index: usize,
    /// Ascending exact rational coefficients; the zero row is empty.
    pub coefficients: Vec<String>,
    /// Degree, absent exactly for the zero polynomial.
    pub degree: Option<usize>,
    /// Distinct real roots, equal to degree for every nonzero prefix row.
    pub distinct_real_roots: usize,
}

impl EffectiveBkwReport {
    /// Check stored factorization, sign resultants, samples, boxes and finite rows.
    /// Phase endpoints, domination and Newton closure additionally regenerate
    /// the current supplier and compare the complete report. Legacy proofs
    /// without `domination_cutoff` are rejected and must be regenerated.
    ///
    /// # Errors
    /// Returns [`RecurrenceError::VerificationFailed`] on any mismatch.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        replay::witnesses(self)?;
        // Phase endpoints, domination inequalities and Newton closure still
        // use deterministic supplier regeneration, always with current guards.
        let expected = if let Some(source) = &self.source_recurrence {
            certify_effective_bkw_recurrence(source)?
        } else {
            certify_effective_bkw(&self.request)?
        };
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
/// A zero prefix or an improper numerator is removed exactly when possible;
/// the report binds the shifted GF and checked prefix to the original record.
///
/// # Errors
/// Returns [`RecurrenceError`] when the recurrence cannot be converted to
/// an index-free rational generating function or is malformed.
pub fn certify_effective_bkw_recurrence(
    record: &crate::recurrence::RecurrenceJson,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    certify_recurrence_with_mode(record, false)
}

pub(crate) fn certify_recurrence_with_mode(
    record: &crate::recurrence::RecurrenceJson,
    legacy: bool,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let gf = crate::recurrence::rational_generating_function(&rule, first, &seeds)?;
    match certify_recurrence_from_gf(record, gf.clone(), seeds.len(), legacy) {
        Err(RecurrenceError::Unsupported(reason)) => Ok(EffectiveBkwReport {
            schema: EFFECTIVE_BKW_SCHEMA.into(),
            request: EffectiveBkwRequest {
                generating_function: gf,
                interval: None,
                target_zero_count: None,
                uniform_constants: None,
                pair_factorization: None,
            },
            status: "not_proved".into(),
            reason,
            certificate: None,
            source_recurrence: Some(record.clone()),
            skipped_prefix: Vec::new(),
        }),
        result => result,
    }
}

fn certify_recurrence_from_gf(
    record: &crate::recurrence::RecurrenceJson,
    mut gf: RationalGeneratingFunction,
    seed_count: usize,
    legacy: bool,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    let mut numerator = decode_b(&gf.numerator)?;
    let denominator = decode_b(&gf.denominator)?;
    let order = denominator.index_degree().unwrap_or(0);
    let mut skipped = Vec::new();
    // The complete seed-dependent GF is shifted, never rebuilt from a
    // truncated seed list. This retains all prescribed exceptions to the rule.
    while skipped.len() <= seed_count
        && (numerator
            .rows()
            .first()
            .is_none_or(|r| vec_poly::is_zero(r))
            || numerator.index_degree().is_some_and(|d| d >= order))
    {
        if denominator.rows().first() != Some(&vec![Q::one()]) {
            return Err(failure("prefix shifting requires D(t,0)=1"));
        }
        if numerator.is_zero() || order < 2 {
            break;
        }
        let row = P::new(numerator.rows().first().cloned().unwrap_or_default());
        let degree = row.degree();
        let roots = crate::root_count::count_real_roots_prs_bigint_coeffs(
            &vec_poly::clear_denominators(row.coeffs()),
        );
        if degree.is_some_and(|d| d != roots) {
            return Err(failure(
                "removed finite prefix contains a row that is not simple real-rooted",
            ));
        }
        skipped.push(EffectiveBkwPrefixRow {
            index: gf.first_index,
            coefficients: encode(&row),
            degree,
            distinct_real_roots: roots,
        });
        numerator = numerator.subtract(&denominator.multiply(&t(&row)));
        if numerator
            .rows()
            .first()
            .is_some_and(|r| !vec_poly::is_zero(r))
        {
            return Err(failure(
                "prefix subtraction did not cancel its constant row",
            ));
        }
        numerator = B::new(numerator.rows().iter().skip(1).cloned().collect());
        gf.first_index = gf
            .first_index
            .checked_add(1)
            .ok_or_else(|| failure("shift index overflow"))?;
    }
    gf.numerator = encode_b(&numerator);
    let mut report = certify_with_mode(
        &EffectiveBkwRequest {
            generating_function: gf,
            interval: None,
            target_zero_count: None,
            uniform_constants: None,
            pair_factorization: None,
        },
        legacy,
    )?;
    if !skipped.is_empty() {
        report.source_recurrence = Some(record.clone());
        report.skipped_prefix = skipped;
        if report.certificate.is_some() {
            report.reason = format!("every row from index {} is nonzero and has only simple real zeros; removed prefix checked exactly and bound to the original recurrence", report.request.generating_function.first_index);
        }
    }
    Ok(report)
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
    certify_with_mode(request, false)
}

fn certify_with_mode(
    request: &EffectiveBkwRequest,
    legacy: bool,
) -> Result<EffectiveBkwReport, RecurrenceError> {
    let mut report = EffectiveBkwReport {
        schema: EFFECTIVE_BKW_SCHEMA.into(),
        request: request.clone(),
        status: "not_proved".into(),
        reason: String::new(),
        certificate: None,
        source_recurrence: None,
        skipped_prefix: Vec::new(),
    };
    match certify(request, legacy) {
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

fn certify(
    request: &EffectiveBkwRequest,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    if request
        .interval
        .as_ref()
        .is_some_and(|i| matches!(i.upper, EffectiveBkwEndpoint::PositiveInfinity))
    {
        return reflected_certificate(request, legacy);
    }
    match certify_in_chart(request, legacy) {
        Err(ChartError::NoNegativeCubicCollision) if request.interval.is_none() => {
            reflected_certificate(request, legacy)
        }
        result => result.map_err(Into::into),
    }
}

fn reflected_certificate(
    request: &EffectiveBkwRequest,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    if request.pair_factorization.is_some() {
        return Err(failure(
            "reflection of supplied higher-order factors is not implemented",
        ));
    }
    let reflected = reflect_request(request)?;
    let mut certificate = certify_in_chart(&reflected, legacy).map_err(RecurrenceError::from)?;
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

fn certify_in_chart(
    request: &EffectiveBkwRequest,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, ChartError> {
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
    if order == 2 && request.uniform_constants.is_some() && !legacy {
        return Err(RecurrenceError::InvalidInput(
            "uniform_constants are unsupported for order-two no-remainder proofs".into(),
        )
        .into());
    }
    if !(2..=8).contains(&order) {
        return Err(failure("effective pair suppliers support orders two through eight").into());
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
    if !legacy {
        // An exact early rejection, independent of the asymptotic degree or
        // phase suppliers. A later zero row is outside this route's nonzero
        // simple-root conclusion, even when the initial numerator is proper.
        let first = P::new(numerator.rows().first().cloned().unwrap_or_default());
        let next =
            P::new(numerator.rows().get(1).cloned().unwrap_or_default()) - d[1].clone() * first;
        if next.is_zero() {
            if let Some(k) = degree::zero_edge_obstruction(&numerator, &denominator) {
                return Err(failure(format!("local row 1 is zero; moreover the zero-end Newton face has {k} equal limiting modes: at n divisible by {k}, remainder/pair amplitude tends ({k}-2)/2>1, excluding global pair domination near t=0" )).into());
            }
            return Err(failure("local row 1 is zero; the nonzero simple-root conclusion requires a different tail or a zero/multiplicity-aware supplier").into());
        }
    }
    // Two finite collisions with theta limits pi,0 supply n+O(1) roots.
    // A sublinear Newton cap cannot close that count, regardless of the
    // numerator's fixed intercept. Reject before any expensive elimination.
    if (order >= if legacy { 4 } else { 5 } || request.pair_factorization.is_some())
        && d.iter()
            .enumerate()
            .skip(1)
            .all(|(j, p)| p.degree().is_none_or(|degree| degree < j))
    {
        if !legacy
            && request
                .interval
                .as_ref()
                .is_none_or(|i| matches!(i.lower, EffectiveBkwEndpoint::NegativeInfinity))
        {
            if let Some(k) = degree::equimodular_obstruction(&numerator, &denominator) {
                return Err(failure(format!("negative-infinity Newton face has {k} equal limiting modes and constant numerator: at every n divisible by {k}, remainder/pair amplitude tends ({k}-2)/2>1; no global domination cutoff exists for this pair scheme")).into());
            }
        }
        return Err(failure(
            if legacy {"compact pair phase has slope one but the Newton degree cap is sublinear"} else {"higher-order fractional ray requires a multi-mode remainder and Newton edge supplier; the compact phase has slope one but the Newton degree cap is sublinear"},
        )
        .into());
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
    let mut certificate =
        if order == 2 && degree.data.newton.is_some() && request.pair_factorization.is_none() {
            match compact::certify(&numerator, &d, request, &degree, legacy) {
                Ok(certificate) => certificate,
                Err(RecurrenceError::Unsupported(_)) => {
                    quadratic::certify(&numerator, &d, request, &degree, legacy)?
                }
                Err(error) => return Err(error.into()),
            }
        } else if order == 3 && request.pair_factorization.is_none() {
            cubic::certify(&numerator, &d, request, &degree, legacy)?
        } else {
            compact::certify(&numerator, &d, request, &degree, legacy)?
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
