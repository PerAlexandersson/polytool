//! The operator menu of the preserver-image lifts ([`PreserverOperator`]),
//! its internal normal forms, its action on rows and the conjugation of an
//! input relation through it (module documentation of `preserver_image`).

use crate::recurrence::implicit::algebra::{from_theta, Biv, DOp, ThetaOp};
use crate::recurrence::lift::hadamard_factor::prepare::{shift_free_from_json, shift_free_json};
use crate::recurrence::lift::section_certify::n_polynomial;
use crate::recurrence::{
    format_rational_coeff, parse_rational_coeff, DegreeSchedule, HadamardFactor,
    HadamardFactorTerm, ImageRowMap, OperatorTermJson, RecurrenceError,
};
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};
use ore_algebra::{pull_back, PullbackBounds, RowMapStep};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Largest degree of an operator polynomial `f(n, theta)` or `phi(n, D)`, of
/// a power `x^m` and of a derivative order `D^j` in [`PreserverOperator`].
pub const MAX_PRESERVER_OPERATOR_DEGREE: usize = 6;

/// A rational function `numerator(n) / denominator(n)` of the index `n`
/// (ascending exact rational strings; the denominator defaults to `1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IndexRationalFunction {
    /// Ascending coefficients of the numerator in `n`.
    pub numerator: Vec<String>,
    /// Ascending coefficients of the denominator in `n` (nonzero polynomial).
    #[serde(default = "one_polynomial")]
    pub denominator: Vec<String>,
}

fn one_polynomial() -> Vec<String> {
    vec!["1".into()]
}

impl IndexRationalFunction {
    /// The constant `value`.
    pub fn constant(value: &BigRational) -> Self {
        Self {
            numerator: vec![format_rational_coeff(value)],
            denominator: one_polynomial(),
        }
    }

    /// `numerator(n) / denominator(n)` from exact coefficients.
    pub fn from_polynomials(numerator: &[BigRational], denominator: &[BigRational]) -> Self {
        let text = |p: &[Q]| {
            vec_poly::trim(p)
                .iter()
                .map(format_rational_coeff)
                .collect()
        };
        Self {
            numerator: text(numerator),
            denominator: text(denominator),
        }
    }

    /// The parsed `(numerator, denominator)`, trimmed; the denominator must
    /// be a nonzero polynomial.
    pub fn parse(&self) -> Result<(Vec<BigRational>, Vec<BigRational>), RecurrenceError> {
        let parse = |values: &[String], what: &str| -> Result<Vec<Q>, RecurrenceError> {
            let parsed = values
                .iter()
                .map(|value| {
                    parse_rational_coeff(value)
                        .map_err(|error| RecurrenceError::from(error).with_context(what))
                })
                .collect::<Result<Vec<_>, _>>()?;
            if parsed.len() > MAX_PRESERVER_OPERATOR_DEGREE * 4 + 1 {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "the {what} has more than {} coefficients",
                    MAX_PRESERVER_OPERATOR_DEGREE * 4 + 1
                )));
            }
            Ok(vec_poly::trim(&parsed))
        };
        let numerator = parse(&self.numerator, "coefficient numerator")?;
        let denominator = parse(&self.denominator, "coefficient denominator")?;
        if vec_poly::is_zero(&denominator) {
            return Err(RecurrenceError::InvalidInput(
                "an operator coefficient has the zero denominator".into(),
            ));
        }
        Ok((numerator, denominator))
    }

    /// Plain-text form, for example `(2n+1)/(n+1)`.
    pub fn text(&self) -> String {
        let polynomial = |values: &[String]| -> String {
            let mut terms = Vec::new();
            for (power, value) in values.iter().enumerate() {
                if value == "0" {
                    continue;
                }
                terms.push(match power {
                    0 => value.clone(),
                    1 if value == "1" => "n".into(),
                    1 => format!("{value} n"),
                    _ if value == "1" => format!("n^{power}"),
                    _ => format!("{value} n^{power}"),
                });
            }
            if terms.is_empty() {
                "0".into()
            } else {
                terms.join(" + ")
            }
        };
        let numerator = polynomial(&self.numerator);
        if self.denominator.len() <= 1 && self.denominator.first().is_none_or(|d| d == "1") {
            numerator
        } else {
            format!("({numerator})/({})", polynomial(&self.denominator))
        }
    }
}

/// A linear operator `T_n` acting on the parent rows (module documentation
/// of `preserver_image`: the menu, its transfer theorems and side
/// conditions).  Every index-dependent quantity is a function of the
/// *input* index `n`: the candidate is `Q_n = lambda(n) T_n(P_(n+s))`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "shape", rename_all = "snake_case", deny_unknown_fields)]
#[non_exhaustive]
pub enum PreserverOperator {
    /// The diagonal operator `x^k -> f(n, k) x^k`, that is `f(n, theta)` with
    /// the Euler operator `theta = x D` and
    /// `f(n, theta) = sum_j coefficients[j](n) theta^j` (a multiplier
    /// sequence `k -> f(n, k)` for each `n` when `f(n, .)` is real-rooted with
    /// every root `<= 0`).
    EulerPolynomial {
        /// `c_j(n)`, `j = 0, 1, ...`.
        coefficients: Vec<IndexRationalFunction>,
    },
    /// The differential operator `phi(n, D) = sum_j coefficients[j](n) D^j`
    /// (a Laguerre--Pólya operator for each `n` when `phi(n, .)` is
    /// real-rooted), for example `P + c(n) P'`.
    DerivativePolynomial {
        /// `c_j(n)`, `j = 0, 1, ...`.
        coefficients: Vec<IndexRationalFunction>,
    },
    /// `P -> D^j (x^m P)`.
    DerivativeOfPowerMultiple {
        /// `m`.
        power: usize,
        /// `j`.
        derivative_order: usize,
    },
    /// `P -> x^m D^j P`.
    PowerTimesDerivative {
        /// `m`.
        power: usize,
        /// `j`.
        derivative_order: usize,
    },
    /// The affine substitution `P -> P(a x + b)`, `a != 0`.
    Affine {
        /// `a`, an exact rational string.
        slope: String,
        /// `b`, an exact rational string.
        intercept: String,
    },
    /// The reversal `P -> x^(d_n) P(1/x)`; the schedule is indexed by the
    /// input index `n`, unbounded and starts at or before `n0`.
    Reversal {
        /// `d_n`.
        schedule: DegreeSchedule,
    },
    /// The Möbius map `P -> (c x + e)^(d_n) P((a x + b)/(c x + e))` with
    /// `c != 0` and `a e - b c != 0` (for example the `f`-polynomial,
    /// `(a, b, c, e) = (1, 0, 1, 1)`); the schedule is as for
    /// [`PreserverOperator::Reversal`].
    Mobius {
        /// `a`.
        a: String,
        /// `b`.
        b: String,
        /// `c`.
        c: String,
        /// `e`.
        e: String,
        /// `d_n`.
        schedule: DegreeSchedule,
    },
    /// The diagonal operator `x^k -> u(n, k) x^k` for a hypergeometric
    /// factor of the Hadamard-factor menu (`1/k!`, `1/(alpha)_k`,
    /// `binom(n + a, k + b)`, ...); certified by the Hadamard conjugation.
    HadamardMultiplier {
        /// The factor `u`.
        factor: HadamardFactor,
    },
    /// The Euler step `kappa (theta + a)(theta + b) + x (u + v theta)` with
    /// index-dependent parameters, given by its expanded coefficients:
    /// `diagonal = [c_0, c_1, c_2]` is the polynomial
    /// `c_0 + c_1 t + c_2 t^2 = kappa (t + a)(t + b)` in `t = theta`, and
    /// `raising = [u, v]`.  It maps polynomials with nonnegative coefficients
    /// and only real roots to such polynomials under the side conditions of
    /// RealRooted `IsPFPolynomial.generalStep_of_natDegree_le` (`kappa, a,
    /// b > 0`, `u + v k > 0`, a negative comparison defect).
    EulerStep {
        /// `[kappa a b, kappa (a + b), kappa]`.
        diagonal: Vec<IndexRationalFunction>,
        /// `[u, v]`.
        raising: Vec<IndexRationalFunction>,
    },
    /// A general element `(1/denominator(n)) sum_d c_d(n, x) D^d` of the
    /// Weyl algebra.  No transfer theorem applies: the identity can be
    /// certified, the lift stays open.
    WeylOperator {
        /// The terms `c_d(n, x) D^d`.
        terms: Vec<OperatorTermJson>,
        /// The common denominator in `n` (ascending; default `1`).
        #[serde(default = "one_polynomial")]
        denominator: Vec<String>,
    },
}

/// One step `F -> G` of a coordinate map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum CoordinateStep {
    /// `G(y) = F(a y + b)`.
    Affine { a: Q, b: Q },
    /// `G_n(y) = y^(d_n) F_n(1/y)` with `d_n = degree(n)` (affine in `n`)
    /// and increment `step >= 0`.
    Reversal { degree: Vec<Q>, step: i64 },
}

/// The internal normal form of an operator.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum OperatorForm {
    /// `T_n = numerator / denominator(n)`, `numerator` a shift-free
    /// operator with coefficients in `(n, x)`.
    Weyl { numerator: DOp, denominator: Vec<Q> },
    /// `T_n = step_last o ... o step_first` (applied first to last).
    Coordinate { steps: Vec<CoordinateStep> },
    /// `x^k -> u(n, k) x^k`.
    Hadamard { factor: HadamardFactor },
}

fn parse_value(text: &str, what: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(text).map_err(|error| RecurrenceError::from(error).with_context(what))
}

/// The coefficients with a common denominator: `(numerators, delta)` with
/// `c_j = numerators[j] / delta` and `delta = prod_j q_j`.
pub(super) fn common_denominator(
    coefficients: &[IndexRationalFunction],
) -> Result<(Vec<Vec<Q>>, Vec<Q>), RecurrenceError> {
    if coefficients.is_empty() || coefficients.len() > MAX_PRESERVER_OPERATOR_DEGREE + 1 {
        return Err(RecurrenceError::LimitExceeded(format!(
            "an operator polynomial needs between 1 and {} coefficients",
            MAX_PRESERVER_OPERATOR_DEGREE + 1
        )));
    }
    let parsed = coefficients
        .iter()
        .map(IndexRationalFunction::parse)
        .collect::<Result<Vec<_>, _>>()?;
    if parsed.iter().all(|(p, _)| vec_poly::is_zero(p)) {
        return Err(RecurrenceError::InvalidInput(
            "the operator polynomial is zero".into(),
        ));
    }
    let delta = parsed
        .iter()
        .fold(vec![Q::one()], |acc, (_, q)| vec_poly::mul(&acc, q));
    let numerators = parsed
        .iter()
        .enumerate()
        .map(|(j, (p, _))| {
            parsed
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != j)
                .fold(p.clone(), |acc, (_, (_, q))| vec_poly::mul(&acc, q))
        })
        .collect();
    Ok((numerators, delta))
}

fn check_degree(value: usize, what: &str) -> Result<(), RecurrenceError> {
    if value > MAX_PRESERVER_OPERATOR_DEGREE {
        return Err(RecurrenceError::LimitExceeded(format!(
            "the {what} {value} exceeds {MAX_PRESERVER_OPERATOR_DEGREE}"
        )));
    }
    Ok(())
}

fn validate_schedule(schedule: &DegreeSchedule) -> Result<Vec<Q>, RecurrenceError> {
    schedule.validate()?;
    if schedule.last_index.is_some() || schedule.step < 0 {
        return Err(RecurrenceError::InvalidInput(
            "the operator schedule must be unbounded with a nonnegative step".into(),
        ));
    }
    let step = Q::from_integer(schedule.step.into());
    let first = Q::from_integer(schedule.first_index.into());
    let degree = Q::from_integer(schedule.first_degree.into());
    Ok(vec_poly::trim(&[degree - &step * first, step]))
}

impl PreserverOperator {
    /// `f(n, theta)` from exact constant coefficients.
    pub fn euler_polynomial(coefficients: &[BigRational]) -> Self {
        Self::EulerPolynomial {
            coefficients: coefficients
                .iter()
                .map(IndexRationalFunction::constant)
                .collect(),
        }
    }

    /// `phi(n, D)` from exact constant coefficients.
    pub fn derivative_polynomial(coefficients: &[BigRational]) -> Self {
        Self::DerivativePolynomial {
            coefficients: coefficients
                .iter()
                .map(IndexRationalFunction::constant)
                .collect(),
        }
    }

    /// The `f`-polynomial map `P -> (1 + x)^(d_n) P(x/(1 + x))`.
    pub fn f_polynomial(schedule: DegreeSchedule) -> Self {
        Self::Mobius {
            a: "1".into(),
            b: "0".into(),
            c: "1".into(),
            e: "1".into(),
            schedule,
        }
    }

    /// Plain-text form, for example `phi(n, D) = 1 + (1/2) D`.
    pub fn text(&self) -> String {
        let polynomial = |coefficients: &[IndexRationalFunction], variable: &str| {
            let terms: Vec<String> = coefficients
                .iter()
                .enumerate()
                .filter(|(_, c)| c.numerator.iter().any(|v| v != "0"))
                .map(|(j, c)| match j {
                    0 => c.text(),
                    1 => format!("({}) {variable}", c.text()),
                    _ => format!("({}) {variable}^{j}", c.text()),
                })
                .collect();
            terms.join(" + ")
        };
        match self {
            Self::EulerPolynomial { coefficients } => {
                format!("f(n, theta) = {}", polynomial(coefficients, "theta"))
            }
            Self::DerivativePolynomial { coefficients } => {
                format!("phi(n, D) = {}", polynomial(coefficients, "D"))
            }
            Self::DerivativeOfPowerMultiple {
                power,
                derivative_order,
            } => format!("D^{derivative_order} o x^{power}"),
            Self::PowerTimesDerivative {
                power,
                derivative_order,
            } => format!("x^{power} D^{derivative_order}"),
            Self::Affine { slope, intercept } => format!("P(({slope}) x + ({intercept}))"),
            Self::Reversal { schedule } => format!(
                "x^(d_n) P(1/x), d_n = {} + {} (n - {})",
                schedule.first_degree, schedule.step, schedule.first_index
            ),
            Self::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => format!(
                "(({c}) x + ({e}))^(d_n) P((({a}) x + ({b}))/(({c}) x + ({e}))), d_n = {} + {} \
                 (n - {})",
                schedule.first_degree, schedule.step, schedule.first_index
            ),
            Self::HadamardMultiplier { factor } => format!("x^k -> {} x^k", factor.text()),
            Self::EulerStep { diagonal, raising } => format!(
                "kappa (theta + a)(theta + b) + x (u + v theta) with kappa (t + a)(t + b) = {} \
                 and u + v theta = {}",
                polynomial(diagonal, "t"),
                polynomial(raising, "theta")
            ),
            Self::WeylOperator { terms, .. } => {
                format!("a Weyl-algebra operator with {} terms", terms.len())
            }
        }
    }

    /// The internal normal form, after validating the parameters.
    pub(super) fn form(&self) -> Result<OperatorForm, RecurrenceError> {
        match self {
            Self::EulerPolynomial { coefficients } => {
                let (numerators, delta) = common_denominator(coefficients)?;
                // f'(n, theta) = sum_j numerators[j](n) theta^j as one Biv in (n, theta).
                let mut biv = Biv::zero();
                for (j, p) in numerators.iter().enumerate() {
                    biv = biv.add(&n_polynomial(p).multiply(&Biv::monomial(0, j, Q::one())));
                }
                let mut theta = ThetaOp::new();
                theta.insert((0, 0), biv);
                Ok(OperatorForm::Weyl {
                    numerator: from_theta(&theta),
                    denominator: delta,
                })
            }
            Self::DerivativePolynomial { coefficients } => {
                let (numerators, delta) = common_denominator(coefficients)?;
                let mut op = DOp::default();
                for (j, p) in numerators.iter().enumerate() {
                    op.add_term((0, j), n_polynomial(p));
                }
                Ok(OperatorForm::Weyl {
                    numerator: op,
                    denominator: delta,
                })
            }
            Self::DerivativeOfPowerMultiple {
                power,
                derivative_order,
            } => {
                check_degree(*power, "power")?;
                check_degree(*derivative_order, "derivative order")?;
                let d = DOp::single(0, *derivative_order, Biv::one());
                let x = DOp::single(0, 0, Biv::monomial(0, *power, Q::one()));
                Ok(OperatorForm::Weyl {
                    numerator: d.compose(&x),
                    denominator: vec![Q::one()],
                })
            }
            Self::PowerTimesDerivative {
                power,
                derivative_order,
            } => {
                check_degree(*power, "power")?;
                check_degree(*derivative_order, "derivative order")?;
                Ok(OperatorForm::Weyl {
                    numerator: DOp::single(
                        0,
                        *derivative_order,
                        Biv::monomial(0, *power, Q::one()),
                    ),
                    denominator: vec![Q::one()],
                })
            }
            Self::Affine { slope, intercept } => {
                let a = parse_value(slope, "affine slope")?;
                let b = parse_value(intercept, "affine intercept")?;
                if a.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "the affine slope a must be nonzero".into(),
                    ));
                }
                Ok(OperatorForm::Coordinate {
                    steps: vec![CoordinateStep::Affine { a, b }],
                })
            }
            Self::Reversal { schedule } => Ok(OperatorForm::Coordinate {
                steps: vec![CoordinateStep::Reversal {
                    degree: validate_schedule(schedule)?,
                    step: schedule.step,
                }],
            }),
            Self::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => {
                let a = parse_value(a, "Mobius a")?;
                let b = parse_value(b, "Mobius b")?;
                let c = parse_value(c, "Mobius c")?;
                let e = parse_value(e, "Mobius e")?;
                if c.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "a Mobius map needs c != 0 (use the affine operator otherwise)".into(),
                    ));
                }
                let determinant = &a * &e - &b * &c;
                if determinant.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "a Mobius map needs a e - b c != 0".into(),
                    ));
                }
                // (c x + e)^d P((a x + b)/(c x + e)) = A2(R(A1(P))) with
                // A1: F(y) -> F(a/c + kappa y), kappa = (b c - a e)/c,
                // R: G(u) -> u^d G(1/u), A2: K(u) -> K(c x + e).
                let kappa = -determinant / &c;
                Ok(OperatorForm::Coordinate {
                    steps: vec![
                        CoordinateStep::Affine {
                            a: kappa,
                            b: a / &c,
                        },
                        CoordinateStep::Reversal {
                            degree: validate_schedule(schedule)?,
                            step: schedule.step,
                        },
                        CoordinateStep::Affine { a: c, b: e },
                    ],
                })
            }
            Self::HadamardMultiplier { factor } => {
                if factor.terms.is_empty() || factor.terms.len() > 2 {
                    return Err(RecurrenceError::InvalidInput(
                        "a Hadamard multiplier has one or two terms".into(),
                    ));
                }
                Ok(OperatorForm::Hadamard {
                    factor: factor.clone(),
                })
            }
            Self::EulerStep { diagonal, raising } => {
                if diagonal.len() != 3 || raising.len() != 2 {
                    return Err(RecurrenceError::InvalidInput(
                        "an Euler step has three diagonal and two raising coefficients".into(),
                    ));
                }
                let all: Vec<IndexRationalFunction> =
                    diagonal.iter().chain(raising).cloned().collect();
                let (numerators, delta) = common_denominator(&all)?;
                let polynomial_in_theta = |coefficients: &[Vec<Q>]| {
                    let mut biv = Biv::zero();
                    for (j, p) in coefficients.iter().enumerate() {
                        biv = biv.add(&n_polynomial(p).multiply(&Biv::monomial(0, j, Q::one())));
                    }
                    biv
                };
                let mut theta = ThetaOp::new();
                let diagonal_part = polynomial_in_theta(&numerators[..3]);
                if !diagonal_part.is_zero() {
                    theta.insert((0, 0), diagonal_part);
                }
                let raising_part = polynomial_in_theta(&numerators[3..]);
                if !raising_part.is_zero() {
                    theta.insert((0, 1), raising_part);
                }
                Ok(OperatorForm::Weyl {
                    numerator: from_theta(&theta),
                    denominator: delta,
                })
            }
            Self::WeylOperator { terms, denominator } => {
                if terms.is_empty() {
                    return Err(RecurrenceError::InvalidInput(
                        "the Weyl operator has no terms".into(),
                    ));
                }
                if terms
                    .iter()
                    .any(|t| t.deriv_order > MAX_PRESERVER_OPERATOR_DEGREE)
                {
                    return Err(RecurrenceError::LimitExceeded(format!(
                        "a Weyl-operator derivative order exceeds {MAX_PRESERVER_OPERATOR_DEGREE}"
                    )));
                }
                let numerator = shift_free_from_json(terms)?;
                if numerator.is_zero() {
                    return Err(RecurrenceError::InvalidInput(
                        "the Weyl operator is zero".into(),
                    ));
                }
                let (_, denominator) = IndexRationalFunction {
                    numerator: vec!["1".into()],
                    denominator: denominator.clone(),
                }
                .parse()?;
                Ok(OperatorForm::Weyl {
                    numerator,
                    denominator,
                })
            }
        }
    }

    /// The schedule of a reversal or Möbius operator.
    pub(super) fn schedule(&self) -> Option<&DegreeSchedule> {
        match self {
            Self::Reversal { schedule } | Self::Mobius { schedule, .. } => Some(schedule),
            _ => None,
        }
    }
}

/// The ratio `u(k + 1)/u(k) = A(k)/B(k)` of a factor `u(k)` that does not
/// depend on `n`, with `u(0) = 1`: `k!`, `1/k!`, `(alpha)_k`, `1/(alpha)_k`
/// and their products.  `None` for a factor that depends on `n`.
fn index_free_factor_ratio(factor: &HadamardFactor) -> Option<(Vec<Q>, Vec<Q>)> {
    let mut numerator = vec![Q::one()];
    let mut denominator = vec![Q::one()];
    for term in &factor.terms {
        let (a, b) = match term {
            HadamardFactorTerm::Factorial => (vec![Q::one(), Q::one()], vec![Q::one()]),
            HadamardFactorTerm::InverseFactorial => (vec![Q::one()], vec![Q::one(), Q::one()]),
            HadamardFactorTerm::RisingFactorial { parameter } => (
                vec![parse_rational_coeff(parameter).ok()?, Q::one()],
                vec![Q::one()],
            ),
            HadamardFactorTerm::InverseRisingFactorial { parameter } => (
                vec![Q::one()],
                vec![parse_rational_coeff(parameter).ok()?, Q::one()],
            ),
            HadamardFactorTerm::Binomial { .. }
            | HadamardFactorTerm::ComplementFactorial
            | HadamardFactorTerm::InverseComplementFactorial => return None,
        };
        numerator = vec_poly::mul(&numerator, &a);
        denominator = vec_poly::mul(&denominator, &b);
    }
    Some((numerator, denominator))
}

fn coefficient_strings(p: &[Q]) -> Vec<String> {
    vec_poly::trim(p)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

/// The row map of [`crate::recurrence::image_recurrence_of_operator`] with
/// the same action on the parent rows: every Weyl form (`f(n, theta)`,
/// `phi(n, D)`, the Rolle operators, the Euler step, a general Weyl
/// operator) as [`ImageRowMap::WeylOperator`] with the same numerator and
/// denominator, the coordinate maps as themselves, and a hypergeometric
/// multiplier that does not depend on `n` (`k!`, `1/k!`, `(alpha)_k`,
/// `1/(alpha)_k` and their products) as
/// [`ImageRowMap::CoefficientMultiplier`].
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] for a multiplier that depends on `n`
/// (`binom(n + a, k + b)`, `(n - k)!`, `1/(n - k)!`), and the validation
/// errors of the operator's parameters.
impl TryFrom<&PreserverOperator> for ImageRowMap {
    type Error = RecurrenceError;

    fn try_from(operator: &PreserverOperator) -> Result<Self, Self::Error> {
        Ok(match operator {
            PreserverOperator::Affine { slope, intercept } => {
                operator.form()?;
                ImageRowMap::Affine {
                    slope: slope.clone(),
                    intercept: intercept.clone(),
                }
            }
            PreserverOperator::Reversal { schedule } => {
                operator.form()?;
                ImageRowMap::Reversal {
                    schedule: schedule.clone(),
                }
            }
            PreserverOperator::Mobius {
                a,
                b,
                c,
                e,
                schedule,
            } => {
                operator.form()?;
                ImageRowMap::Mobius {
                    a: a.clone(),
                    b: b.clone(),
                    c: c.clone(),
                    e: e.clone(),
                    schedule: schedule.clone(),
                }
            }
            PreserverOperator::HadamardMultiplier { factor } => {
                operator.form()?;
                let (numerator, denominator) =
                    index_free_factor_ratio(factor).ok_or_else(|| {
                        RecurrenceError::Unsupported(format!(
                            "the multiplier {} depends on n; image rules support multipliers \
                             u(k) only",
                            factor.text()
                        ))
                    })?;
                ImageRowMap::CoefficientMultiplier {
                    numerator: coefficient_strings(&numerator),
                    denominator: coefficient_strings(&denominator),
                }
            }
            _ => match operator.form()? {
                OperatorForm::Weyl {
                    numerator,
                    denominator,
                } => ImageRowMap::WeylOperator {
                    terms: shift_free_json(&numerator),
                    denominator: coefficient_strings(&denominator),
                },
                _ => {
                    return Err(RecurrenceError::Unsupported(
                        "the operator has no image row map".into(),
                    ))
                }
            },
        })
    }
}

impl OperatorForm {
    /// The polynomials in `n` that must not vanish at any `n >= n0` for the
    /// operator to be defined there: the denominator of a Weyl form.
    pub(super) fn index_denominator(&self) -> Option<&[Q]> {
        match self {
            Self::Weyl { denominator, .. } => Some(denominator),
            _ => None,
        }
    }

    /// `T_n(row)`; `None` where it is undefined (a vanishing denominator, a
    /// reversal of a row of too high degree, an undefined factor value).
    pub(super) fn apply_row(&self, row: &[Q], n: usize) -> Option<Vec<Q>> {
        let index = Q::from_integer(n.into());
        match self {
            Self::Weyl {
                numerator,
                denominator,
            } => {
                let delta = vec_poly::evaluate(denominator, &index);
                if delta.is_zero() {
                    return None;
                }
                let image = numerator.apply_row(row, n);
                Some(vec_poly::trim(&vec_poly::scale(&image, &delta.recip())))
            }
            Self::Coordinate { steps } => {
                let mut current = vec_poly::trim(row);
                for step in steps {
                    current = match step {
                        CoordinateStep::Affine { a, b } => {
                            vec_poly::trim(&vec_poly::compose(&current, &[b.clone(), a.clone()]))
                        }
                        CoordinateStep::Reversal { degree, .. } => {
                            let d = vec_poly::evaluate(degree, &index);
                            if !d.is_integer() || d < Q::zero() {
                                return None;
                            }
                            let d = usize::try_from(d.to_integer()).ok()?;
                            vec_poly::trim(&vec_poly::reverse(&current, d)?)
                        }
                    };
                }
                Some(current)
            }
            Self::Hadamard { factor } => {
                let mut out = Vec::with_capacity(row.len());
                for (k, p) in row.iter().enumerate() {
                    if p.is_zero() {
                        out.push(Q::zero());
                    } else {
                        out.push(factor.value(n, k)? * p);
                    }
                }
                Some(vec_poly::trim(&out))
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Conjugation of the input relation
// ---------------------------------------------------------------------------

impl OperatorForm {
    /// The map `P' -> T'(P')` as steps of the `ore-algebra` crate (applied
    /// first to last): the Weyl numerator `T'`, or the coordinate steps.
    /// `None` for a Hadamard form (certified by the Hadamard conjugation).
    pub(super) fn row_map_steps(&self) -> Option<Vec<RowMapStep>> {
        match self {
            Self::Weyl { numerator, .. } => Some(vec![RowMapStep::WeylOperator(numerator.clone())]),
            Self::Coordinate { steps } => Some(
                steps
                    .iter()
                    .map(|step| match step {
                        CoordinateStep::Affine { a, b } => RowMapStep::Affine {
                            slope: a.clone(),
                            intercept: b.clone(),
                        },
                        CoordinateStep::Reversal { degree, .. } => RowMapStep::Reversal {
                            degree: degree.clone(),
                        },
                    })
                    .collect(),
            ),
            Self::Hadamard { .. } => None,
        }
    }
}

/// The dividend `E` of the certificate (module documentation of
/// `preserver_image`, steps 1 and 2): the pullback
/// (`ore_algebra::pull_back`) of the input relation `M` (leader shift `R`)
/// through `T'` followed by the gauge `mu(n)/mu(n-1) = alpha(n)/beta(n)`,
/// with the bounds of `M` and top shift `R`, divided by the largest power
/// of `x` dividing every coefficient.  So `E = M_mu o T'` for a Weyl form
/// and the conjugation of `M_mu` through the coordinate steps (last step
/// first, with `y^N`, `N = b + h R`, for a reversal) otherwise; the gauge
/// factor of the pullback is `prod_(t=1..R) beta(n + t)`.  Dividing by
/// `x^v` keeps `(E P')_n = 0` equivalent to the undivided identity, since
/// `Q[x]` (and the Laurent polynomials) have no zero divisors.  `None` for
/// a Hadamard form or an invalid step.
pub(super) fn dividend(
    input: &DOp,
    leader_shift: usize,
    alpha: &[Q],
    beta: &[Q],
    form: &OperatorForm,
) -> Option<DOp> {
    let mut steps = form.row_map_steps()?;
    steps.push(RowMapStep::Gauge {
        numerator: alpha.to_vec(),
        denominator: beta.to_vec(),
    });
    let bounds = PullbackBounds {
        top_shift: leader_shift.max(input.max_shift()),
        ..PullbackBounds::of(input)
    };
    let (e, _) = pull_back(input, &steps, bounds).ok()?;
    Some(e.divide_by_variable_power(e.variable_valuation()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    fn rows(seed: u64) -> Vec<Vec<Q>> {
        let mut state = seed;
        (0..9)
            .map(|n| {
                (0..n + 2)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6364136223846793005)
                            .wrapping_add(1442695040888963407);
                        q(((state >> 33) % 7) as i64 - 3)
                    })
                    .collect()
            })
            .collect()
    }

    /// The conjugation identity on arbitrary rows: `(M H)_n` vanishes
    /// exactly when `(E F)_n` does, and is its image (up to `mu`, `beta`,
    /// `y^N` and the coordinate map).
    #[test]
    fn conjugation_matches_direct_application() {
        // M = (n + 2) S^2 - (1 + x) D S + x^2 (n + 1).
        let mut m = DOp::default();
        m.add_term((2, 0), Biv::new(vec![vec![q(2)], vec![q(1)]]));
        m.add_term((1, 1), Biv::from_variable_polynomial(&[q(-1), q(-1)]));
        m.add_term(
            (0, 0),
            Biv::new(vec![vec![q(0), q(0), q(1)], vec![q(0), q(0), q(1)]]),
        );
        let n0 = 1;
        let parent = rows(11);
        let operators = [
            PreserverOperator::Affine {
                slope: "2".into(),
                intercept: "1".into(),
            },
            PreserverOperator::Reversal {
                schedule: DegreeSchedule::new(0, 3, 1),
            },
            PreserverOperator::f_polynomial(DegreeSchedule::new(0, 2, 1)),
            PreserverOperator::Mobius {
                a: "2".into(),
                b: "-1".into(),
                c: "3".into(),
                e: "1/2".into(),
                schedule: DegreeSchedule::new(1, 3, 2),
            },
            PreserverOperator::EulerPolynomial {
                coefficients: vec![
                    IndexRationalFunction::from_polynomials(&[q(1)], &[q(1), q(1)]),
                    IndexRationalFunction::constant(&q(2)),
                    IndexRationalFunction::constant(&q(1)),
                ],
            },
            PreserverOperator::DerivativePolynomial {
                coefficients: vec![
                    IndexRationalFunction::constant(&q(1)),
                    IndexRationalFunction::from_polynomials(&[q(1), q(1)], &[q(2)]),
                ],
            },
            PreserverOperator::DerivativeOfPowerMultiple {
                power: 2,
                derivative_order: 1,
            },
        ];
        let alpha = [q(1), q(2)];
        let beta = [q(3)];
        let product_beta = |n: usize| {
            (1..=2).fold(Q::one(), |acc, t| {
                acc * vec_poly::evaluate(&beta, &q((n + t) as i64))
            })
        };
        for operator in operators {
            let form = operator.form().unwrap();
            let e = dividend(&m, 2, &alpha, &beta, &form).unwrap();
            // T' (the Weyl numerator, or the coordinate map) on the rows.
            let apply = |row: &[Q], n: usize| match &form {
                OperatorForm::Weyl { numerator, .. } => numerator.apply_row(row, n),
                _ => form.apply_row(row, n).unwrap(),
            };
            // mu(n0) = 1, mu(n)/mu(n-1) = alpha(n)/beta(n).
            let mut mus = vec![Q::one()];
            for offset in 1..parent.len() {
                let value = q((n0 + offset) as i64);
                let next = mus[offset - 1].clone() * vec_poly::evaluate(&alpha, &value)
                    / vec_poly::evaluate(&beta, &value);
                mus.push(next);
            }
            let images: Vec<Vec<Q>> = parent
                .iter()
                .enumerate()
                .map(|(offset, row)| vec_poly::scale(&apply(row, n0 + offset), &mus[offset]))
                .collect();
            for base in n0..n0 + parent.len() - 2 {
                let direct = m.apply(&images, n0, base).unwrap();
                let conjugated = e.apply(&parent, n0, base).unwrap();
                assert!(!vec_poly::is_zero(&direct), "the test rows are generic");
                assert!(
                    !vec_poly::is_zero(&conjugated),
                    "{operator:?} at base {base}"
                );
                let lhs = vec_poly::scale(&direct, &product_beta(base));
                match &form {
                    // (M H)_n prod beta(n + t) = mu(n) (E F)_n exactly.
                    OperatorForm::Weyl { .. } => assert_eq!(
                        lhs,
                        vec_poly::scale(&conjugated, &mus[base - n0]),
                        "{operator:?} at base {base}"
                    ),
                    // ... = mu(n) (E F)_n (a x + b) for an affine map.
                    OperatorForm::Coordinate { steps } if steps.len() == 1 => {
                        if let CoordinateStep::Affine { a, b } = &steps[0] {
                            let image = vec_poly::compose(&conjugated, &[b.clone(), a.clone()]);
                            assert_eq!(lhs, vec_poly::scale(&image, &mus[base - n0]));
                        }
                    }
                    _ => {}
                }
            }
        }
    }
}
