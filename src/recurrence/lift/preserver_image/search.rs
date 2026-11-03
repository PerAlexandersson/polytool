//! Detection of preserver-image lifts ([`find_preserver_image_lifts`]): for
//! every library parent, index shift and operator shape, fit the operator
//! parameters exactly from a few rows and certify every match.

use super::certify::{
    certify_preserver_image_lift, PreserverImageCandidate, PreserverImageCertification,
    PreserverImageOptions, MAX_PRESERVER_IMAGE_CHECKED_ROWS,
};
use super::library::{PreserverParent, PreserverParentLibrary};
use super::operator::{IndexRationalFunction, PreserverOperator};
use crate::linalg::solve_linear_system_with_rank;
use crate::recurrence::lift::hadamard_factor::prepare::prepare_lift_rule;
use crate::recurrence::{
    format_rational_coeff, DegreeSchedule, HadamardFactor, HadamardFactorTerm, Recurrence,
    RecurrenceError, VeroneseSectionGauge,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::sync::Arc;

type Q = BigRational;

/// JSON schema tag of a serialized [`PreserverImageReport`].
pub const PRESERVER_IMAGE_LIFTS_SCHEMA: &str = "polytool.recurrence-preserver-image-lifts.v1";

/// One operator shape of the detection menu.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PreserverOperatorShape {
    /// `f(n, theta)`, fitted coefficient by coefficient.
    EulerPolynomial,
    /// `phi(n, D)` of degree at least one, fitted coefficient by coefficient.
    DerivativePolynomial,
    /// `D^j o x^m`.
    DerivativeOfPowerMultiple,
    /// `x^m D^j` with `m >= 1`.
    PowerTimesDerivative,
    /// `P(a x + b)`.
    Affine,
    /// `x^(d_n) P(1/x)`.
    Reversal,
    /// The `f`-polynomial `(1 + x)^(d_n) P(x/(1 + x))`.
    FPolynomial,
    /// A hypergeometric multiplier (`1/k!`, `1/(1/2)_k`, `1/(n - k)!`,
    /// `binom(n + a, k + b)` with `0 <= a <= 2`, `0 <= b <= 1`).
    HadamardMultiplier,
    /// The Euler step `kappa (theta + a)(theta + b) + x (u + v theta)`,
    /// fitted coefficient by coefficient (with a raising term).
    EulerStep,
}

impl PreserverOperatorShape {
    /// Every shape, in the order tried.
    pub fn all() -> Vec<Self> {
        vec![
            Self::EulerPolynomial,
            Self::DerivativePolynomial,
            Self::DerivativeOfPowerMultiple,
            Self::PowerTimesDerivative,
            Self::Affine,
            Self::Reversal,
            Self::FPolynomial,
            Self::HadamardMultiplier,
            Self::EulerStep,
        ]
    }
}

/// Bounds of [`find_preserver_image_lifts`].
#[derive(Debug, Clone)]
pub struct PreserverImageSearchOptions {
    /// The parents; `None` builds [`PreserverParentLibrary::named_families`]
    /// for this call (share a built library across calls to cache the
    /// parent proofs).
    pub library: Option<Arc<PreserverParentLibrary>>,
    /// The shapes tried, in order.  Default [`PreserverOperatorShape::all`].
    pub shapes: Vec<PreserverOperatorShape>,
    /// Index shifts `s` with `|s| <=` this bound are tried (`0, 1, -1, 2,
    /// ...`).  Default 2.
    pub max_index_shift: usize,
    /// Largest degree of `f(n, theta)` and `phi(n, D)`, of the power `m`
    /// and of the derivative order `j`.  Default 2.
    pub max_operator_degree: usize,
    /// Largest total degree of the rational functions of `n` fitted to the
    /// operator coefficients and the gauge.  Default 2.
    pub max_coefficient_degree: usize,
    /// Input rows used for detection, between 4 and
    /// [`MAX_PRESERVER_IMAGE_CHECKED_ROWS`].  Default 10.
    pub detection_rows: usize,
    /// Largest number of candidates certified; reaching it sets
    /// [`PreserverImageReport::exhausted_budget`].  Default 64.
    pub max_candidates: usize,
    /// Stop at the first proved lift.  Default true.
    pub stop_at_first_proof: bool,
    /// The options of each certification.
    pub certification: PreserverImageOptions,
}

impl Default for PreserverImageSearchOptions {
    fn default() -> Self {
        Self {
            library: None,
            shapes: PreserverOperatorShape::all(),
            max_index_shift: 2,
            max_operator_degree: 2,
            max_coefficient_degree: 2,
            detection_rows: 10,
            max_candidates: 64,
            stop_at_first_proof: true,
            certification: PreserverImageOptions::default(),
        }
    }
}

/// `polytool.recurrence-preserver-image-lifts.v1`: the result of
/// [`find_preserver_image_lifts`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PreserverImageReport {
    /// Always [`PRESERVER_IMAGE_LIFTS_SCHEMA`].
    pub schema: String,
    /// The input's first index `n0`.
    pub input_first_index: usize,
    /// Input rows used for detection.
    pub rows_examined: usize,
    /// Library parents tried.
    pub parents_tried: usize,
    /// Candidates certified.
    pub candidates_tried: usize,
    /// Every certified candidate (proved or not), in the order found.
    pub lifts: Vec<PreserverImageCertification>,
    /// Whether the candidate budget was reached.
    pub exhausted_budget: bool,
    /// Further remarks.
    pub notes: Vec<String>,
}

impl PreserverImageReport {
    /// The proved lifts.
    pub fn proved_lifts(&self) -> impl Iterator<Item = &PreserverImageCertification> {
        self.lifts.iter().filter(|lift| lift.is_proved())
    }

    /// Replay every certification.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != PRESERVER_IMAGE_LIFTS_SCHEMA {
            return Err(RecurrenceError::InvalidInput(format!(
                "unsupported schema `{}`",
                self.schema
            )));
        }
        for lift in &self.lifts {
            lift.verify()?;
            if lift.input_first_index != self.input_first_index {
                return Err(RecurrenceError::VerificationFailed(
                    "a lift has another first index".into(),
                ));
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Exact fitting
// ---------------------------------------------------------------------------

fn q_usize(value: usize) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// A rational function `p(n)/q(n)` of total degree at most `max_total`
/// through every sample `(n, v)` (with at least one sample to spare), the
/// smallest degrees first; `q` monic.
fn fit_rational(samples: &[(usize, Q)], max_total: usize) -> Option<IndexRationalFunction> {
    let first = &samples.first()?.1;
    if samples.iter().all(|(_, v)| v == first) {
        return Some(IndexRationalFunction::constant(first));
    }
    for total in 1..=max_total {
        for den in 0..=total {
            let num = total - den;
            let unknowns = num + 1 + den;
            if samples.len() < unknowns + 1 {
                continue;
            }
            // p(n) - v (q_0 + ... + q_(den-1) n^(den-1)) = v n^den.
            let mut matrix = Vec::with_capacity(samples.len());
            let mut rhs = Vec::with_capacity(samples.len());
            for (n, v) in samples {
                let x = q_usize(*n);
                let mut row = Vec::with_capacity(unknowns);
                let mut power = Q::one();
                for _ in 0..=num {
                    row.push(power.clone());
                    power *= &x;
                }
                let mut power = Q::one();
                for _ in 0..den {
                    row.push(-(v * &power));
                    power *= &x;
                }
                rhs.push(v * power);
                matrix.push(row);
            }
            let Some(solution) = solve_linear_system_with_rank(&matrix, &rhs) else {
                continue;
            };
            if solution.nullity != 0 {
                continue;
            }
            let numerator = solution.particular[..=num].to_vec();
            let mut denominator = solution.particular[num + 1..].to_vec();
            denominator.push(Q::one());
            if samples
                .iter()
                .any(|(n, _)| vec_poly::evaluate(&denominator, &q_usize(*n)).is_zero())
            {
                continue;
            }
            return Some(IndexRationalFunction::from_polynomials(
                &numerator,
                &denominator,
            ));
        }
    }
    None
}

/// The gauge `lambda` through the samples `(n, lambda(n))` (consecutive
/// from `n0`, nonzero): `1`, a closed rational function, or a rational
/// ratio `lambda(n)/lambda(n-1)`.
fn fit_gauge(values: &[(usize, Q)], n0: usize, max_total: usize) -> Option<VeroneseSectionGauge> {
    if values.first().map(|(n, _)| *n) != Some(n0) || values.iter().any(|(_, v)| v.is_zero()) {
        return None;
    }
    if values.iter().all(|(_, v)| v.is_one()) {
        return Some(VeroneseSectionGauge::None);
    }
    let valid = |gauge: VeroneseSectionGauge| gauge.ratio_form(n0).is_ok().then_some(gauge);
    if let Some(closed) = fit_rational(values, max_total) {
        if let Some(gauge) = valid(VeroneseSectionGauge::Closed {
            numerator: closed.numerator,
            denominator: closed.denominator,
        }) {
            return Some(gauge);
        }
    }
    let ratios: Vec<(usize, Q)> = values
        .windows(2)
        .filter(|w| w[1].0 == w[0].0 + 1)
        .map(|w| (w[1].0, &w[1].1 / &w[0].1))
        .collect();
    if ratios.len() + 1 != values.len() {
        return None;
    }
    let ratio = fit_rational(&ratios, max_total)?;
    valid(VeroneseSectionGauge::Ratio {
        numerator: ratio.numerator,
        denominator: ratio.denominator,
        initial_value: format_rational_coeff(&values[0].1),
    })
}

/// `lambda` with `target = lambda image` exactly, or `None`.
fn scalar_ratio(target: &[Q], image: &[Q]) -> Option<Q> {
    let target = vec_poly::trim(target);
    let image = vec_poly::trim(image);
    let degree = vec_poly::degree(&image)?;
    if vec_poly::degree(&target) != Some(degree) {
        return None;
    }
    let lambda = &target[degree] / &image[degree];
    (vec_poly::scale(&image, &lambda) == target).then_some(lambda)
}

/// The gauge for images `I_n` of the parent rows with `Q_n = lambda(n) I_n`
/// on every detection row.
fn gauge_for_images(
    inputs: &[Vec<Q>],
    images: &[Option<Vec<Q>>],
    n0: usize,
    max_total: usize,
) -> Option<VeroneseSectionGauge> {
    let mut values = Vec::with_capacity(inputs.len());
    for (offset, (target, image)) in inputs.iter().zip(images).enumerate() {
        values.push((n0 + offset, scalar_ratio(target, image.as_ref()?)?));
    }
    fit_gauge(&values, n0, max_total)
}

fn theta_image(row: &[Q], power: usize) -> Vec<Q> {
    row.iter()
        .enumerate()
        .map(|(k, c)| {
            let mut weight = Q::one();
            for _ in 0..power {
                weight *= q_usize(k);
            }
            c * weight
        })
        .collect()
}

/// The basis `P, theta P, theta^2 P, x P, x theta P` of the Euler step.
fn euler_step_image(row: &[Q], index: usize) -> Vec<Q> {
    match index {
        0..=2 => theta_image(row, index),
        _ => {
            let mut shifted = vec![Q::zero()];
            shifted.extend(theta_image(row, index - 3));
            shifted
        }
    }
}

fn derivative_image(row: &[Q], order: usize) -> Vec<Q> {
    let mut out = vec_poly::trim(row);
    for _ in 0..order {
        out = vec_poly::derivative(&out);
    }
    out
}

/// Fits of a linear shape `Q_n = sum_j y_j(n) B_j(P'_n)`: the coefficients
/// as rational functions (no gauge), or normalised by the first nonzero
/// coefficient with a fitted gauge.
fn fit_linear_shape(
    inputs: &[Vec<Q>],
    parents: &[Vec<Q>],
    n0: usize,
    degree: usize,
    basis: impl Fn(&[Q], usize) -> Vec<Q>,
    max_total: usize,
) -> Vec<(Vec<IndexRationalFunction>, VeroneseSectionGauge)> {
    let mut samples: Vec<(usize, Vec<Q>)> = Vec::new();
    for (offset, (target, parent)) in inputs.iter().zip(parents).enumerate() {
        let images: Vec<Vec<Q>> = (0..=degree).map(|j| basis(parent, j)).collect();
        let width = images
            .iter()
            .map(Vec::len)
            .chain(std::iter::once(target.len()))
            .max()
            .unwrap_or(0);
        let matrix: Vec<Vec<Q>> = (0..width)
            .map(|k| {
                images
                    .iter()
                    .map(|image| vec_poly::coeff(image, k))
                    .collect()
            })
            .collect();
        let rhs: Vec<Q> = (0..width).map(|k| vec_poly::coeff(target, k)).collect();
        let Some(solution) = solve_linear_system_with_rank(&matrix, &rhs) else {
            return Vec::new();
        };
        if solution.nullity == 0 {
            samples.push((n0 + offset, solution.particular));
        }
    }
    if samples.len() < 4 {
        return Vec::new();
    }
    let mut used = degree;
    while used > 0 && samples.iter().all(|(_, y)| y[used].is_zero()) {
        used -= 1;
    }
    let mut out = Vec::new();
    // Coefficients as rational functions of n.
    let direct: Option<Vec<IndexRationalFunction>> = (0..=used)
        .map(|j| {
            let values: Vec<(usize, Q)> = samples.iter().map(|(n, y)| (*n, y[j].clone())).collect();
            fit_rational(&values, max_total)
        })
        .collect();
    if let Some(coefficients) = direct {
        out.push((coefficients, VeroneseSectionGauge::None));
    }
    // Normalised by a pivot coefficient, with a gauge.
    if let Some(pivot) = (0..=used).find(|j| samples.iter().all(|(_, y)| !y[*j].is_zero())) {
        let lambdas: Vec<(usize, Q)> = samples
            .iter()
            .map(|(n, y)| (*n, y[pivot].clone()))
            .collect();
        let gauge = fit_gauge(&lambdas, n0, max_total);
        let coefficients: Option<Vec<IndexRationalFunction>> = (0..=used)
            .map(|j| {
                let values: Vec<(usize, Q)> = samples
                    .iter()
                    .map(|(n, y)| (*n, &y[j] / &y[pivot]))
                    .collect();
                fit_rational(&values, max_total)
            })
            .collect();
        if let (Some(gauge), Some(coefficients)) = (gauge, coefficients) {
            if gauge != VeroneseSectionGauge::None {
                out.push((coefficients, gauge));
            }
        }
    }
    out
}

/// Exact square root of a nonnegative rational, when it is one.
fn rational_sqrt(value: &Q) -> Option<Q> {
    if value.is_negative() {
        return None;
    }
    let root = |v: &BigInt| {
        let r = v.sqrt();
        (&r * &r == *v).then_some(r)
    };
    Some(Q::new(root(value.numer())?, root(value.denom())?))
}

/// `(e, s1, V)` of a row of degree `e >= 2`: the root sum `s1` and the
/// spread `V = e sum r^2 - s1^2`.
fn root_moments(row: &[Q]) -> Option<(usize, Q, Q)> {
    let row = vec_poly::trim(row);
    let e = vec_poly::degree(&row)?;
    if e < 2 {
        return None;
    }
    let lead = &row[e];
    let s1 = -(&row[e - 1] / lead);
    let e2 = &row[e - 2] / lead;
    let power_sum = &s1 * &s1 - q_usize(2) * e2;
    let spread = q_usize(e) * power_sum - &s1 * &s1;
    Some((e, s1, spread))
}

/// Candidate affine maps `(a, b)` with `Q_n = lambda P'_n(a x + b)`, from the
/// root moments of the detection rows.
fn affine_candidates(inputs: &[Vec<Q>], parents: &[Vec<Q>]) -> Vec<(Q, Q)> {
    let mut candidates: Option<Vec<(Q, Q)>> = None;
    for (target, parent) in inputs.iter().zip(parents) {
        let (Some((e, s1q, vq)), Some((ep, s1p, vp))) =
            (root_moments(target), root_moments(parent))
        else {
            continue;
        };
        if e != ep || vq.is_zero() || vp.is_zero() {
            continue;
        }
        // V(Q) = V(P)/a^2 and s1(Q) = (s1(P) - e b)/a.
        let Some(a) = rational_sqrt(&(&vp / &vq)) else {
            return Vec::new();
        };
        let row_candidates: Vec<(Q, Q)> = [a.clone(), -a]
            .into_iter()
            .map(|a| {
                let b = (&s1p - &a * &s1q) / q_usize(e);
                (a, b)
            })
            .collect();
        candidates = Some(match candidates {
            None => row_candidates,
            Some(previous) => previous
                .into_iter()
                .filter(|c| row_candidates.contains(c))
                .collect(),
        });
    }
    candidates
        .unwrap_or_default()
        .into_iter()
        .filter(|(a, b)| !(a.is_one() && b.is_zero()))
        .collect()
}

/// An affine schedule `d_n` through the values (consecutive from `n0`).
fn affine_schedule(values: &[i64], n0: usize) -> Option<DegreeSchedule> {
    let first = *values.first()?;
    let step = values.get(1).map_or(0, |second| second - first);
    if first < 0 || step < 0 {
        return None;
    }
    let consistent = values
        .iter()
        .enumerate()
        .all(|(i, v)| *v == first + step * i as i64);
    consistent.then(|| DegreeSchedule::new(n0, first as usize, step))
}

fn hadamard_menu() -> Vec<HadamardFactor> {
    let mut terms = vec![
        HadamardFactorTerm::InverseFactorial,
        HadamardFactorTerm::InverseRisingFactorial {
            parameter: "1/2".into(),
        },
        HadamardFactorTerm::InverseComplementFactorial,
    ];
    for a in 0..=2 {
        for b in 0..=1 {
            terms.push(HadamardFactorTerm::Binomial {
                index_offset: a,
                coefficient_offset: b,
            });
        }
    }
    terms.into_iter().map(HadamardFactor::single).collect()
}

/// The operator candidates (operator, gauge) of one shape for one parent
/// and shift, fitted on the detection rows.
fn shape_candidates(
    shape: PreserverOperatorShape,
    inputs: &[Vec<Q>],
    parents: &[Vec<Q>],
    n0: usize,
    options: &PreserverImageSearchOptions,
) -> Vec<(PreserverOperator, VeroneseSectionGauge)> {
    let degree = options.max_operator_degree;
    let total = options.max_coefficient_degree;
    let with_images =
        |operator: PreserverOperator| -> Option<(PreserverOperator, VeroneseSectionGauge)> {
            let form = operator.form().ok()?;
            let images: Vec<Option<Vec<Q>>> = parents
                .iter()
                .enumerate()
                .map(|(offset, row)| form.apply_row(row, n0 + offset))
                .collect();
            let gauge = gauge_for_images(inputs, &images, n0, total)?;
            Some((operator, gauge))
        };
    match shape {
        PreserverOperatorShape::EulerPolynomial => {
            fit_linear_shape(inputs, parents, n0, degree, theta_image, total)
                .into_iter()
                .map(|(coefficients, gauge)| {
                    (PreserverOperator::EulerPolynomial { coefficients }, gauge)
                })
                .collect()
        }
        PreserverOperatorShape::DerivativePolynomial => {
            fit_linear_shape(inputs, parents, n0, degree, derivative_image, total)
                .into_iter()
                .filter(|(coefficients, _)| coefficients.len() >= 2)
                .map(|(coefficients, gauge)| {
                    (
                        PreserverOperator::DerivativePolynomial { coefficients },
                        gauge,
                    )
                })
                .collect()
        }
        PreserverOperatorShape::DerivativeOfPowerMultiple => {
            let mut out = Vec::new();
            for derivative_order in 1..=degree {
                for power in 0..=degree {
                    out.extend(with_images(PreserverOperator::DerivativeOfPowerMultiple {
                        power,
                        derivative_order,
                    }));
                }
            }
            out
        }
        PreserverOperatorShape::PowerTimesDerivative => {
            let mut out = Vec::new();
            for power in 1..=degree {
                for derivative_order in 1..=degree {
                    out.extend(with_images(PreserverOperator::PowerTimesDerivative {
                        power,
                        derivative_order,
                    }));
                }
            }
            out
        }
        PreserverOperatorShape::Affine => affine_candidates(inputs, parents)
            .into_iter()
            .filter_map(|(a, b)| {
                with_images(PreserverOperator::Affine {
                    slope: format_rational_coeff(&a),
                    intercept: format_rational_coeff(&b),
                })
            })
            .collect(),
        PreserverOperatorShape::Reversal => {
            let mut degrees = Vec::with_capacity(inputs.len());
            for (target, parent) in inputs.iter().zip(parents) {
                let (Some(dq), Some(valuation)) = (
                    vec_poly::degree(target),
                    parent.iter().position(|c| !c.is_zero()),
                ) else {
                    return Vec::new();
                };
                degrees.push((dq + valuation) as i64);
            }
            affine_schedule(&degrees, n0)
                .and_then(|schedule| with_images(PreserverOperator::Reversal { schedule }))
                .into_iter()
                .collect()
        }
        PreserverOperatorShape::FPolynomial => {
            let Some(degrees) = parents
                .iter()
                .map(|row| vec_poly::degree(row).map(|d| d as i64))
                .collect::<Option<Vec<i64>>>()
            else {
                return Vec::new();
            };
            (0..=1)
                .filter_map(|extra| {
                    let shifted: Vec<i64> = degrees.iter().map(|d| d + extra).collect();
                    affine_schedule(&shifted, n0)
                })
                .filter_map(|schedule| with_images(PreserverOperator::f_polynomial(schedule)))
                .collect()
        }
        PreserverOperatorShape::EulerStep => {
            fit_linear_shape(inputs, parents, n0, 4, euler_step_image, total)
                .into_iter()
                .filter(|(coefficients, _)| coefficients.len() >= 4)
                .map(|(mut coefficients, gauge)| {
                    coefficients.resize(5, IndexRationalFunction::constant(&Q::zero()));
                    let raising = coefficients.split_off(3);
                    (
                        PreserverOperator::EulerStep {
                            diagonal: coefficients,
                            raising,
                        },
                        gauge,
                    )
                })
                .collect()
        }
        PreserverOperatorShape::HadamardMultiplier => hadamard_menu()
            .into_iter()
            .filter(|factor| factor.validate(n0).is_ok())
            .filter_map(|factor| with_images(PreserverOperator::HadamardMultiplier { factor }))
            .collect(),
    }
}

/// Index shifts in the order tried: `0, 1, -1, 2, -2, ...`.
fn shifts(bound: usize) -> Vec<i64> {
    let mut out = vec![0];
    for s in 1..=bound as i64 {
        out.push(s);
        out.push(-s);
    }
    out
}

/// Find preserver-image lifts of the sequence defined by `recurrence` and
/// `rows` (`rows[0]` is `Q_(first_index)`, ascending coefficients):
/// `Q_n = lambda(n) T_n(P_(n+s))` for a library parent `P` and an operator
/// `T_n` of the menu (module documentation of `preserver_image`).
///
/// For each parent (library order), index shift (`0, 1, -1, ...`) and
/// shape ([`PreserverOperatorShape`] order) the parameters are fitted
/// exactly on the detection rows: the polynomial shapes by solving
/// `Q_n = sum_j y_j(n) B_j(P_(n+s))` row by row and fitting each `y_j` (or
/// `y_j / y_pivot` with a gauge) as a rational function of `n`; the affine
/// map from the root moments; the schedules from the degrees; the gauge
/// from exact row ratios.  Every match is certified with
/// [`certify_preserver_image_lift`].  Bounded and deterministic; budgets
/// are reported, never a panic.  A candidate whose certification fails
/// (invalid parameters, a budget, a failed internal self-check) is
/// rejected with a note in [`PreserverImageReport::notes`] and the search
/// goes on; it never aborts the search.  A missing lift says nothing
/// against real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError`] for an input rule that does not reproduce its rows,
/// is not homogeneous or has a leader not proved nonvanishing, and for
/// options beyond their bounds.
pub fn find_preserver_image_lifts(
    recurrence: &Recurrence,
    rows: &[Vec<BigRational>],
    first_index: usize,
    options: &PreserverImageSearchOptions,
) -> Result<PreserverImageReport, RecurrenceError> {
    if !(4..=MAX_PRESERVER_IMAGE_CHECKED_ROWS).contains(&options.detection_rows)
        || options.max_operator_degree > super::operator::MAX_PRESERVER_OPERATOR_DEGREE
        || options.max_coefficient_degree > 4
    {
        return Err(RecurrenceError::LimitExceeded(format!(
            "detection_rows must be between 4 and {MAX_PRESERVER_IMAGE_CHECKED_ROWS}, \
             max_operator_degree at most {} and max_coefficient_degree at most 4",
            super::operator::MAX_PRESERVER_OPERATOR_DEGREE
        )));
    }
    let n0 = first_index;
    let input = prepare_lift_rule(recurrence, rows, n0, "input")?;
    let detection = options.detection_rows;
    let inputs = input.rows(detection)?;
    let library = match &options.library {
        Some(library) => library.clone(),
        None => Arc::new(PreserverParentLibrary::named_families()?),
    };
    let mut report = PreserverImageReport {
        schema: PRESERVER_IMAGE_LIFTS_SCHEMA.into(),
        input_first_index: n0,
        rows_examined: detection,
        parents_tried: 0,
        candidates_tried: 0,
        lifts: Vec::new(),
        exhausted_budget: false,
        notes: Vec::new(),
    };
    for parent in &library.parents {
        report.parents_tried += 1;
        for shift in shifts(options.max_index_shift) {
            let Some(parent_rows) = shifted_rows(parent, n0, shift, detection) else {
                continue;
            };
            for shape in &options.shapes {
                for (operator, gauge) in
                    shape_candidates(*shape, &inputs, &parent_rows, n0, options)
                {
                    if report.candidates_tried >= options.max_candidates {
                        report.exhausted_budget = true;
                        report.notes.push(format!(
                            "the candidate budget of {} was reached",
                            options.max_candidates
                        ));
                        return Ok(report);
                    }
                    report.candidates_tried += 1;
                    let candidate = PreserverImageCandidate {
                        parent_name: parent.name.clone(),
                        parent: parent.recurrence.clone(),
                        index_shift: shift,
                        gauge,
                        operator,
                    };
                    let certification = match certify_preserver_image_lift(
                        recurrence,
                        rows,
                        n0,
                        &candidate,
                        &options.certification,
                    ) {
                        Ok(certification) => certification,
                        // Every failure of one candidate (invalid parameters,
                        // a budget, or a self-check of its certificate that
                        // fails) rejects that candidate only: it is recorded
                        // as a note and proves nothing.
                        Err(error) => {
                            report.notes.push(format!(
                                "{} (shift {shift}, {}): rejected: {error}",
                                parent.name,
                                candidate.operator.text()
                            ));
                            continue;
                        }
                    };
                    let proved = certification.is_proved();
                    report.lifts.push(certification);
                    if proved && options.stop_at_first_proof {
                        return Ok(report);
                    }
                }
            }
        }
    }
    Ok(report)
}

/// The parent rows `P_(n+s)`, `n = n0, ..., n0 + count - 1`, from the
/// library's precomputed rows.
fn shifted_rows(
    parent: &PreserverParent,
    n0: usize,
    shift: i64,
    count: usize,
) -> Option<Vec<Vec<Q>>> {
    let start = usize::try_from(n0 as i64 + shift).ok()?;
    let offset = start.checked_sub(parent.first_index())?;
    let rows = parent.rows();
    (offset + count <= rows.len()).then(|| rows[offset..offset + count].to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;
    use crate::recurrence::PreserverImageIdentityStatus;

    #[test]
    fn rational_functions_are_fitted_exactly() {
        // (2n + 1)/(n + 3).
        let samples: Vec<(usize, Q)> = (1..8)
            .map(|n| (n, Q::new((2 * n + 1).into(), (n + 3).into())))
            .collect();
        let fitted = fit_rational(&samples, 2).unwrap();
        assert_eq!(
            fitted.parse().unwrap(),
            (vec![q(1), q(2)], vec![q(3), q(1)])
        );
        // n! needs a ratio gauge.
        let mut factorial = Q::one();
        let values: Vec<(usize, Q)> = (1..8)
            .map(|n| {
                factorial *= q(n as i64);
                (n, factorial.clone())
            })
            .collect();
        assert!(matches!(
            fit_gauge(&values, 1, 2),
            Some(VeroneseSectionGauge::Ratio { .. })
        ));
    }

    #[test]
    fn affine_maps_are_read_from_root_moments() {
        // P = (x + 1)(x + 2)(x + 4); Q = P(2x + 1).
        let p = vec_poly::mul(&vec_poly::mul(&[q(1), q(1)], &[q(2), q(1)]), &[q(4), q(1)]);
        let target = vec_poly::scale(&vec_poly::compose(&p, &[q(1), q(2)]), &q(3));
        let candidates = affine_candidates(&[target], &[p]);
        assert!(candidates.contains(&(q(2), q(1))));
    }

    #[test]
    fn the_euler_step_is_fitted_from_independent_images() {
        // T = (theta + 1)(theta + 2) + x (1 + theta) on the Eulerian rows
        // A_1, ..., A_10 (whose five basis images are independent from A_4
        // on; the binomial rows are not, since x P = (theta P + x theta P)/n).
        let parents: Vec<Vec<Q>> = crate::sequences::eulerian_polynomials_bigint(10)[1..]
            .iter()
            .map(|row| {
                vec_poly::trim(&row.iter().cloned().map(Q::from_integer).collect::<Vec<_>>())
            })
            .collect();
        let inputs: Vec<Vec<Q>> = parents
            .iter()
            .map(|p| {
                let diagonal = [q(2), q(3), q(1)]
                    .iter()
                    .enumerate()
                    .fold(vec![q(0)], |acc, (j, c)| {
                        vec_poly::add(&acc, &vec_poly::scale(&euler_step_image(p, j), c))
                    });
                let raising = vec_poly::add(&euler_step_image(p, 3), &euler_step_image(p, 4));
                vec_poly::trim(&vec_poly::add(&diagonal, &raising))
            })
            .collect();
        let options = PreserverImageSearchOptions::default();
        let candidates = shape_candidates(
            PreserverOperatorShape::EulerStep,
            &inputs,
            &parents,
            1,
            &options,
        );
        let constant = |v: i64| IndexRationalFunction::constant(&q(v));
        let expected = PreserverOperator::EulerStep {
            diagonal: vec![constant(2), constant(3), constant(1)],
            raising: vec![constant(1), constant(1)],
        };
        assert!(
            candidates
                .iter()
                .any(|(operator, gauge)| *operator == expected
                    && *gauge == VeroneseSectionGauge::None),
            "{candidates:?}"
        );
        // The operator reproduces the rows.
        let form = expected.form().unwrap();
        for (offset, (input, parent)) in inputs.iter().zip(&parents).enumerate() {
            assert_eq!(
                form.apply_row(parent, 1 + offset)
                    .map(|r| vec_poly::trim(&r)),
                Some(input.clone())
            );
        }
    }

    /// Whether the nonzero `row` is real-rooted with every root in `region`
    /// (exact: Sturm counts on the integer row, signs of the coefficients
    /// of the row shifted to each endpoint, the endpoint itself excluded
    /// for an open end).
    fn row_in_region(row: &[Q], region: &crate::real_rootedness::RealInterval) -> bool {
        let one_sign = |p: &[Q]| {
            let p = vec_poly::trim(p);
            !(p.iter().any(Signed::is_positive) && p.iter().any(Signed::is_negative))
        };
        let common = row.iter().fold(BigInt::one(), |acc, c| {
            num_integer::Integer::lcm(&acc, c.denom())
        });
        let integers: Vec<BigInt> = row
            .iter()
            .map(|c| (c * Q::from_integer(common.clone())).to_integer())
            .collect();
        if !crate::is_real_rooted_bigint_coeffs(&integers) {
            return false;
        }
        let lower_ok = region.lower.as_ref().is_none_or(|l| {
            one_sign(&vec_poly::compose(row, &[l.clone(), -Q::one()]))
                && (region.lower_closed || !vec_poly::evaluate(row, l).is_zero())
        });
        let upper_ok = region.upper.as_ref().is_none_or(|u| {
            one_sign(&vec_poly::compose(row, &[u.clone(), Q::one()]))
                && (region.upper_closed || !vec_poly::evaluate(row, u).is_zero())
        });
        lower_ok && upper_ok
    }

    /// The search on a catalog entry with a library of the given catalog
    /// parents completes, rejects no candidate by a failed self-check, and
    /// every proved lift replays and has its first 40 rows real-rooted in
    /// its claimed region.
    fn assert_search_completes(id: &str, parents: &[&str]) -> PreserverImageReport {
        let library = PreserverParentLibrary::from_catalog_entries(
            parents,
            &super::super::PreserverParentLibraryOptions {
                include_oeis_catalog: true,
                max_catalog_order: 3,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            library.parents.len(),
            parents.len(),
            "{:?}",
            library.skipped
        );
        let entry = crate::oeis::by_id(id).unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let n0 = entry.recurrence_first_index;
        let options = PreserverImageSearchOptions {
            library: Some(Arc::new(library)),
            stop_at_first_proof: false,
            ..Default::default()
        };
        let report =
            find_preserver_image_lifts(&rule, &rows, n0, &options).expect("the search completes");
        report.verify().unwrap();
        assert!(
            !report
                .notes
                .iter()
                .any(|note| note.contains("rejected: the")),
            "{:#?}",
            report.notes
        );
        let count = rule.generation_initial_count(n0, rows.len());
        let generated = rule.generate_rows_rational(&rows[..count], n0, 40).unwrap();
        for lift in report.proved_lifts() {
            let region = &lift.lift_summary.root_transfer.input_root_region;
            for (offset, row) in generated.iter().enumerate() {
                let row = vec_poly::trim(row);
                if vec_poly::is_zero(&row) {
                    continue;
                }
                assert!(
                    row_in_region(&row, &region.at(n0 + offset).unwrap()),
                    "{id}: row {} is not real-rooted in {region}",
                    n0 + offset
                );
            }
        }
        report
    }

    /// Regression: the Euler-step candidates from A011973 and A052553 used
    /// to abort the search on the Lucas triangle A061896 with "the
    /// pseudo-division identity does not hold" (the quotients missed the
    /// index content of the row-annihilator multiplier).
    #[test]
    fn the_search_on_a061896_completes() {
        let report = assert_search_completes("A061896", &["A011973", "A052553"]);
        // Both Euler-step identities are now certified (their side
        // conditions fail, so nothing is proved).
        for parent in ["A011973", "A052553"] {
            assert!(
                report
                    .lifts
                    .iter()
                    .any(|lift| lift.candidate.parent_name == parent
                        && matches!(lift.candidate.operator, PreserverOperator::EulerStep { .. })
                        && matches!(
                            lift.identity,
                            PreserverImageIdentityStatus::Certified { .. }
                        )),
                "{parent}"
            );
        }
    }

    /// Regression: the same for A199478 and the Euler-step candidate from
    /// A128099 at shift 0.
    #[test]
    fn the_search_on_a199478_completes() {
        let report = assert_search_completes("A199478", &["A128099"]);
        assert!(report
            .lifts
            .iter()
            .any(|lift| lift.candidate.index_shift == 0
                && matches!(lift.candidate.operator, PreserverOperator::EulerStep { .. })
                && matches!(
                    lift.identity,
                    PreserverImageIdentityStatus::Certified { .. }
                )));
    }
}
