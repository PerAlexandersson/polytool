//! The transfer theorems of the operator menu, their side conditions
//! (decided exactly for every index with [`crate::integer_range_sign`]) and
//! the root regions they give (module documentation of `preserver_image`,
//! "The preserver menu").

use super::operator::{IndexRationalFunction, PreserverOperator};
use crate::integer_range_sign::{
    decide_index_rational_function_sign, IndexSignCondition, IntegerRange,
};
use crate::recurrence::lift::hadamard_factor::factor::POLYA_SCHUR_THEOREM;
use crate::recurrence::lift::lift_summary::{nonpositive_ray, real_line};
use crate::recurrence::lift::veronese_algebra::classify_row;
use crate::recurrence::{
    format_rational_coeff, parse_rational_coeff, LiftHypothesis, LiftHypothesisStatus,
    RootEndpoint, RootInterval,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

type Q = BigRational;

/// What the menu theorem of one operator gives.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct OperatorTransfer {
    /// The theorem with its reference.
    pub(super) theorem: String,
    /// The side conditions of the theorem for every `n >= n0`.
    pub(super) hypotheses: Vec<LiftHypothesis>,
    /// The input region when every used parent row is zero or real-rooted
    /// with every root in `(-inf, 0]`; `None` when no theorem applies.
    pub(super) from_nonpositive: Option<RootInterval>,
    /// Whether the theorem also gives real-rootedness of the input rows
    /// from parent rows that are real-rooted anywhere on the line.
    pub(super) from_real_line: bool,
    /// The row route of the Lean library RealRooted, when one matches.
    pub(super) lean_route: Option<String>,
    /// RealRooted (or Mathlib) declarations behind the theorem.
    pub(super) lean_declarations: Vec<String>,
}

fn strings(values: &[&str]) -> Vec<String> {
    values.iter().map(|value| value.to_string()).collect()
}

fn verified(name: &str, statement: &str, evidence: String) -> LiftHypothesis {
    LiftHypothesis {
        name: name.into(),
        statement: statement.into(),
        status: LiftHypothesisStatus::Verified,
        evidence,
        sign_decision: None,
    }
}

fn unresolved(name: &str, statement: &str, evidence: String) -> LiftHypothesis {
    LiftHypothesis {
        name: name.into(),
        statement: statement.into(),
        status: LiftHypothesisStatus::Unresolved,
        evidence,
        sign_decision: None,
    }
}

fn failed(name: &str, statement: &str, n0: usize, evidence: String) -> LiftHypothesis {
    LiftHypothesis {
        name: name.into(),
        statement: statement.into(),
        status: LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(n0)),
        },
        evidence,
        sign_decision: None,
    }
}

/// `numerator(n) / denominator(n)` `condition` for every `n >= n0`, decided
/// exactly.
fn sign_hypothesis(
    name: &str,
    statement: &str,
    numerator: &[Q],
    denominator: &[Q],
    n0: usize,
    condition: IndexSignCondition,
    function: &str,
) -> LiftHypothesis {
    LiftHypothesis::from_sign_decision(
        name,
        statement,
        decide_index_rational_function_sign(
            numerator,
            denominator,
            &IntegerRange::starting_at(n0),
            condition,
        ),
        function,
    )
}

/// Parsed coefficients `c_j = p_j / q_j`, `j = 0..=degree` (trailing zero
/// numerators dropped).
type Coefficients = Vec<(Vec<Q>, Vec<Q>)>;

fn parse_coefficients(coefficients: &[IndexRationalFunction]) -> Option<Coefficients> {
    let mut parsed: Coefficients = coefficients
        .iter()
        .map(IndexRationalFunction::parse)
        .collect::<Result<_, _>>()
        .ok()?;
    while parsed
        .last()
        .is_some_and(|(p, _)| vec_poly::is_zero(p) && parsed.len() > 1)
    {
        parsed.pop();
    }
    Some(parsed)
}

/// The constant polynomial `sum_j c_j t^j` when every coefficient is free
/// of `n`.
fn constant_polynomial(coefficients: &Coefficients) -> Option<Vec<Q>> {
    coefficients
        .iter()
        .map(|(p, q)| {
            let constant = |v: &[Q]| vec_poly::degree(v).is_none_or(|d| d == 0);
            (constant(p) && constant(q)).then(|| vec_poly::coeff(p, 0) / vec_poly::coeff(q, 0))
        })
        .collect()
}

/// `c_i c_j >= 0` for every `n >= n0` (the coefficients of a quadratic or
/// linear polynomial have one sign), one hypothesis per pair.
fn one_sign_hypotheses(
    prefix: &str,
    variable: &str,
    coefficients: &Coefficients,
    n0: usize,
) -> Vec<LiftHypothesis> {
    let mut out = Vec::new();
    for i in 0..coefficients.len() {
        for j in i + 1..coefficients.len() {
            let (pi, qi) = &coefficients[i];
            let (pj, qj) = &coefficients[j];
            if vec_poly::is_zero(pi) || vec_poly::is_zero(pj) {
                continue;
            }
            out.push(sign_hypothesis(
                &format!("{prefix}_coefficients_{i}{j}_same_sign"),
                &format!(
                    "the coefficients of {variable}^{i} and {variable}^{j} have the same sign \
                     (or one vanishes) for every n >= n0"
                ),
                &vec_poly::mul(pi, pj),
                &vec_poly::mul(qi, qj),
                n0,
                IndexSignCondition::Nonnegative,
                &format!("c_{i}(n) c_{j}(n)"),
            ));
        }
    }
    out
}

/// `c_1^2 - 4 c_0 c_2 >= 0` for every `n >= n0`.
fn discriminant_hypothesis(
    prefix: &str,
    variable: &str,
    coefficients: &Coefficients,
    n0: usize,
) -> LiftHypothesis {
    let (p0, q0) = &coefficients[0];
    let (p1, q1) = &coefficients[1];
    let (p2, q2) = &coefficients[2];
    // (p1^2 q0 q2 - 4 p0 p2 q1^2) / (q1^2 q0 q2).
    let q1_squared = vec_poly::mul(q1, q1);
    let numerator = vec_poly::sub(
        &vec_poly::mul(&vec_poly::mul(p1, p1), &vec_poly::mul(q0, q2)),
        &vec_poly::scale(
            &vec_poly::mul(&vec_poly::mul(p0, p2), &q1_squared),
            &Q::from_integer(4.into()),
        ),
    );
    let denominator = vec_poly::mul(&q1_squared, &vec_poly::mul(q0, q2));
    sign_hypothesis(
        &format!("{prefix}_discriminant_nonnegative"),
        &format!(
            "the quadratic polynomial in {variable} has a nonnegative discriminant (real roots) \
             for every n >= n0"
        ),
        &vec_poly::trim(&numerator),
        &vec_poly::trim(&denominator),
        n0,
        IndexSignCondition::Nonnegative,
        "c_1(n)^2 - 4 c_0(n) c_2(n)",
    )
}

/// The side conditions of `f(n, theta)` (`require_nonpositive`: real roots,
/// all `<= 0`) or `phi(n, D)` (real roots); returns the hypotheses and
/// whether the coefficients are proved to have one sign for every `n`.
fn polynomial_hypotheses(
    prefix: &str,
    variable: &str,
    coefficients: &[IndexRationalFunction],
    n0: usize,
    require_nonpositive: bool,
) -> (Vec<LiftHypothesis>, bool) {
    let name = if require_nonpositive {
        format!("{prefix}_roots_nonpositive")
    } else {
        format!("{prefix}_real_rooted")
    };
    let statement = if require_nonpositive {
        format!(
            "the polynomial in {variable} is real-rooted with every root <= 0 for every n >= n0"
        )
    } else {
        format!("the polynomial in {variable} is real-rooted for every n >= n0")
    };
    let Some(coefficients) = parse_coefficients(coefficients) else {
        return (
            vec![unresolved(
                &name,
                &statement,
                "the coefficients do not parse".into(),
            )],
            false,
        );
    };
    if let Some(constant) = constant_polynomial(&coefficients) {
        let shape = classify_row(&constant);
        let text = constant
            .iter()
            .map(format_rational_coeff)
            .collect::<Vec<_>>()
            .join(", ");
        let real_rooted = shape.is_zero || shape.is_real_rooted || constant.len() <= 1;
        let holds = real_rooted && (!require_nonpositive || shape.has_one_sign);
        let evidence = format!(
            "the polynomial with ascending coefficients [{text}] is fixed; exact check: \
             real-rooted = {real_rooted}, coefficients of one sign = {}",
            shape.has_one_sign
        );
        let hypothesis = if holds {
            verified(&name, &statement, evidence)
        } else {
            failed(&name, &statement, n0, evidence)
        };
        return (vec![hypothesis], shape.has_one_sign);
    }
    let degree = coefficients.len() - 1;
    match degree {
        0 => (
            vec![verified(
                &name,
                &statement,
                format!("the polynomial has degree 0 in {variable} for every n"),
            )],
            true,
        ),
        1 | 2 => {
            let mut hypotheses = Vec::new();
            if degree == 2 {
                hypotheses.push(discriminant_hypothesis(prefix, variable, &coefficients, n0));
            } else {
                hypotheses.push(verified(
                    &name,
                    &statement,
                    format!("the polynomial has degree <= 1 in {variable}, so it is real-rooted"),
                ));
            }
            let signs = one_sign_hypotheses(prefix, variable, &coefficients, n0);
            let one_sign = signs
                .iter()
                .all(|h| h.status == LiftHypothesisStatus::Verified);
            if require_nonpositive {
                // A real-rooted polynomial has every root <= 0 iff its
                // coefficients have one sign.
                hypotheses.extend(signs);
            }
            (hypotheses, one_sign)
        }
        _ => (
            vec![unresolved(
                &name,
                &statement,
                format!(
                    "the polynomial has degree {degree} in {variable} with index-dependent \
                     coefficients; only fixed polynomials and degrees <= 2 are decided"
                ),
            )],
            false,
        ),
    }
}

/// `operator_coefficients_defined`: the common denominator `delta(n)` of
/// the operator coefficients is nonzero for every `n >= n0`.
pub(super) fn denominator_hypothesis(delta: &[Q], n0: usize) -> LiftHypothesis {
    sign_hypothesis(
        "operator_coefficients_defined",
        "the operator coefficients are finite at every n >= n0 (their denominators do not \
         vanish)",
        delta,
        &[Q::one()],
        n0,
        IndexSignCondition::Nonzero,
        "delta(n)",
    )
}

fn endpoint(value: &Q) -> RootEndpoint {
    RootEndpoint {
        value: format_rational_coeff(value),
        per_index: None,
        closed: true,
    }
}

fn parse(text: &str) -> Option<Q> {
    parse_rational_coeff(text).ok()
}

/// The image of `(-inf, 0]` under the inverse of the Möbius map
/// `x -> (a x + b)/(c x + e)`, together with the point `-e/c` (the roots
/// of the factor `(c x + e)^(d - deg P)`), as one closed interval when the
/// pole `a/c` of the inverse lies outside `(-inf, 0]`.
fn mobius_region(a: &Q, b: &Q, c: &Q, e: &Q) -> RootInterval {
    if c.is_zero() || a.is_zero() {
        return real_line();
    }
    let pole = a / c;
    if !pole.is_positive() {
        return real_line();
    }
    let at_infinity = -(e / c);
    let at_zero = -(b / a);
    let (low, high) = if at_infinity <= at_zero {
        (at_infinity, at_zero)
    } else {
        (at_zero, at_infinity)
    };
    RootInterval {
        lower: Some(endpoint(&low)),
        upper: Some(endpoint(&high)),
    }
}

const LAGUERRE_MULTIPLIER: &str = "Laguerre's theorem on multiplier sequences: if f is a real \
     polynomial whose roots are all real and <= 0, then (f(k))_(k >= 0) is a multiplier \
     sequence, so x^k -> f(k) x^k maps every real-rooted polynomial to a real-rooted polynomial \
     or zero (Laguerre; see T. Craven and G. Csordas, Multiplier sequences for fields, Illinois \
     J. Math. 21 (1977); Polya--Schur, Theorem 40 of P. Branden, arXiv:1410.6601: \
     sum_k f(k) x^k / k! = f(theta) e^x = C e^x x^m prod (1 + x/a_j), a_j > 0); f(k) has one \
     sign for k >= 0, so coefficients of one sign stay of one sign; applied to f(n, .) for \
     each n";

const LAGUERRE_DERIVATIVE: &str = "Laguerre's theorem on differential operators: if phi is a \
     real polynomial with only real roots, then phi(D) maps every real-rooted polynomial to a \
     real-rooted polynomial or zero (write phi(D) = C D^j prod (1 + c_i D); D preserves \
     real-rootedness by Rolle's theorem and p + c p' is real-rooted for every real c by Rolle's \
     theorem applied to e^(x/c) p; Polya--Schur 1914); when the coefficients of phi have one \
     sign, coefficients of one sign stay of one sign; applied to phi(n, .) for each n";

const ROLLE: &str = "Rolle's theorem: x^m p is real-rooted when p is, and every derivative of \
     a real-rooted polynomial is real-rooted or zero, with its roots in the convex hull of the \
     roots; coefficients of one sign stay of one sign";

const EULER_STEP: &str = "The Euler step T = kappa (theta + a)(theta + b) + x (u + v theta), \
     theta = x D, with kappa, a, b > 0, u + v k > 0 for every k <= d and a negative comparison \
     defect kappa (a b - (a + b + 1)^2/4) + (u - (a + b + 1) v/2) s - v^2 s^2/(4 kappa) < 0 for \
     every s <= 0, maps every polynomial of degree <= d with nonnegative coefficients and only \
     real roots to such a polynomial (RealRooted IsPFPolynomial.generalStep_of_natDegree_le, \
     the degree-local strict-to-weak closure of P. Branden, Linear transformations preserving \
     the Polya frequency property, Trans. Amer. Math. Soc. 358 (2006), Corollary 3.5); here \
     u > 0 and v >= 0 give u + v k > 0 for every degree; by linearity a row with \
     nonpositive coefficients is handled through its negative; applied to the parameters at \
     each n";

const SUBSTITUTION: &str = "a real affine substitution x -> a x + b (a != 0), the reversal \
     x^d p(1/x) (d >= deg p) and a real Mobius map (c x + e)^d p((a x + b)/(c x + e)) \
     (a e - b c != 0, d >= deg p) map real roots to real roots: a root r of p gives the root \
     (r - b)/a, 1/r (and 0) or (e r - b)/(a - c r) (and -e/c); for an input row that is a \
     polynomial the degree condition holds automatically";

/// A rational function of `n` as `(numerator, denominator)`, ascending.
type IndexFraction = (Vec<Q>, Vec<Q>);

fn fraction_mul(a: &IndexFraction, b: &IndexFraction) -> IndexFraction {
    (vec_poly::mul(&a.0, &b.0), vec_poly::mul(&a.1, &b.1))
}

fn fraction_add(a: &IndexFraction, b: &IndexFraction) -> IndexFraction {
    (
        vec_poly::add(&vec_poly::mul(&a.0, &b.1), &vec_poly::mul(&b.0, &a.1)),
        vec_poly::mul(&a.1, &b.1),
    )
}

fn fraction_scale(a: &IndexFraction, factor: &Q) -> IndexFraction {
    (vec_poly::scale(&a.0, factor), a.1.clone())
}

fn fraction_divide(a: &IndexFraction, b: &IndexFraction) -> IndexFraction {
    (vec_poly::mul(&a.0, &b.1), vec_poly::mul(&a.1, &b.0))
}

/// The value of `fraction` at `n`, when defined.
fn fraction_at(fraction: &IndexFraction, n: &BigInt) -> Option<Q> {
    let n = Q::from_integer(n.clone());
    let denominator = vec_poly::evaluate(&fraction.1, &n);
    (!denominator.is_zero()).then(|| vec_poly::evaluate(&fraction.0, &n) / denominator)
}

/// The side conditions of the Euler step (the operator's documentation and
/// [`EULER_STEP`]), each for every `n >= n0`.
fn euler_step_hypotheses(
    diagonal: &[IndexRationalFunction],
    raising: &[IndexRationalFunction],
    n0: usize,
) -> Vec<LiftHypothesis> {
    let parsed: Option<Vec<IndexFraction>> = diagonal
        .iter()
        .chain(raising)
        .map(|c| c.parse().ok())
        .collect();
    let Some(parsed) = parsed.filter(|p| p.len() == 5) else {
        return vec![unresolved(
            "euler_step_parameters",
            "the Euler step has three diagonal and two raising coefficients",
            "the coefficients do not parse".into(),
        )];
    };
    let (c0, c1, c2, u, v) = (&parsed[0], &parsed[1], &parsed[2], &parsed[3], &parsed[4]);
    let decide = |name: &str, statement: &str, f: &IndexFraction, condition, text: &str| {
        sign_hypothesis(
            name,
            statement,
            &vec_poly::trim(&f.0),
            &vec_poly::trim(&f.1),
            n0,
            condition,
            text,
        )
    };
    let four = Q::from_integer(BigInt::from(4));
    let discriminant = fraction_add(
        &fraction_mul(c1, c1),
        &fraction_scale(&fraction_mul(c2, c0), &-&four),
    );
    let mut hypotheses = vec![
        decide(
            "euler_step_kappa_positive",
            "kappa(n) = c_2(n) > 0 for every n >= n0",
            c2,
            IndexSignCondition::Positive,
            "c_2(n)",
        ),
        decide(
            "euler_step_diagonal_sum_positive",
            "kappa (a + b) = c_1(n) > 0 for every n >= n0 (with kappa > 0: a + b > 0)",
            c1,
            IndexSignCondition::Positive,
            "c_1(n)",
        ),
        decide(
            "euler_step_diagonal_product_positive",
            "kappa a b = c_0(n) > 0 for every n >= n0 (with kappa > 0: a b > 0)",
            c0,
            IndexSignCondition::Positive,
            "c_0(n)",
        ),
        decide(
            "euler_step_diagonal_real",
            "c_1(n)^2 - 4 c_2(n) c_0(n) >= 0 for every n >= n0 (a and b real, hence both > 0)",
            &discriminant,
            IndexSignCondition::Nonnegative,
            "c_1(n)^2 - 4 c_2(n) c_0(n)",
        ),
        decide(
            "euler_step_raising_constant_positive",
            "u(n) > 0 for every n >= n0",
            u,
            IndexSignCondition::Positive,
            "u(n)",
        ),
        decide(
            "euler_step_raising_slope_nonnegative",
            "v(n) >= 0 for every n >= n0 (with u > 0: u + v k > 0 for every k >= 0)",
            v,
            IndexSignCondition::Nonnegative,
            "v(n)",
        ),
    ];
    // B = u - (a + b + 1) v / 2 with a + b + 1 = (c_1 + c_2)/c_2, and
    // H = (a + b + 1)^2/4 - a b.
    let one_plus_sum = fraction_divide(&fraction_add(c1, c2), c2);
    let half = Q::from_integer(BigInt::from(2)).recip();
    let b = fraction_add(u, &fraction_scale(&fraction_mul(&one_plus_sum, v), &-&half));
    let h = fraction_add(
        &fraction_scale(&fraction_mul(&one_plus_sum, &one_plus_sum), &four.recip()),
        &fraction_scale(&fraction_divide(c0, c2), &-Q::one()),
    );
    let margin = fraction_add(
        &fraction_mul(&fraction_mul(v, v), &h),
        &fraction_scale(&fraction_mul(&b, &b), &-Q::one()),
    );
    let name = "euler_step_comparison_defect_negative";
    let statement = "the comparison defect kappa (a b - (a + b + 1)^2/4) + B s - v^2 s^2/(4 kappa),                      B = u - (a + b + 1) v/2, is < 0 for every s <= 0 and every n >= n0: its                      value at 0 is < 0, so this holds when B >= 0, and otherwise exactly when                      B^2 < v^2 ((a + b + 1)^2/4 - a b)";
    let range = IntegerRange::starting_at(n0);
    let linear = decide_index_rational_function_sign(
        &vec_poly::trim(&b.0),
        &vec_poly::trim(&b.1),
        &range,
        IndexSignCondition::Nonnegative,
    );
    let quadratic = decide_index_rational_function_sign(
        &vec_poly::trim(&margin.0),
        &vec_poly::trim(&margin.1),
        &range,
        IndexSignCondition::Positive,
    );
    hypotheses.push(if linear.holds() {
        LiftHypothesis::from_sign_decision(name, statement, linear, "B(n)")
    } else if quadratic.holds() {
        LiftHypothesis::from_sign_decision(
            name,
            statement,
            quadratic,
            "v(n)^2 ((a + b + 1)^2/4 - a b) - B(n)^2",
        )
    } else {
        // Failing at an index where both tests fail (B < 0 and the margin
        // <= 0 there) is exact; otherwise undecided.
        let both_fail = [linear.failing_index(), quadratic.failing_index()]
            .into_iter()
            .flatten()
            .find(|index| {
                fraction_at(&b, index).is_some_and(|value| value.is_negative())
                    && fraction_at(&margin, index).is_none_or(|value| !value.is_positive())
            });
        match both_fail {
            Some(index) => LiftHypothesis {
                name: name.into(),
                statement: statement.into(),
                status: LiftHypothesisStatus::Failed {
                    witness: Some(index.clone()),
                },
                evidence: format!(
                    "at n = {index}: B < 0 and B^2 >= v^2 ((a + b + 1)^2/4 - a b) (exact), so \
                     the defect is >= 0 at some s <= 0"
                ),
                sign_decision: None,
            },
            None => unresolved(
                name,
                statement,
                "neither B >= 0 nor B^2 < v^2 ((a + b + 1)^2/4 - a b) is proved for every n, and \
                 no index fails both"
                    .into(),
            ),
        }
    });
    hypotheses
}

/// The transfer theorem, side conditions, regions and Lean names of
/// `operator` on the index domain `n >= n0`.
pub(super) fn operator_transfer(operator: &PreserverOperator, n0: usize) -> OperatorTransfer {
    match operator {
        PreserverOperator::EulerPolynomial { coefficients } => {
            let (hypotheses, _) =
                polynomial_hypotheses("euler_polynomial", "theta", coefficients, n0, true);
            OperatorTransfer {
                theorem: LAGUERRE_MULTIPLIER.into(),
                hypotheses,
                from_nonpositive: Some(nonpositive_ray()),
                from_real_line: true,
                lean_route: None,
                lean_declarations: strings(&[
                    "RealRooted.isMultiplierSequence_eval_of_roots_nonpos",
                    "RealRooted.isPFMultiplierSequence_eval_of_roots_nonpos",
                    "RealRooted.isMultiplierSequence_iff_preserves",
                ]),
            }
        }
        PreserverOperator::DerivativePolynomial { coefficients } => {
            let (hypotheses, one_sign) =
                polynomial_hypotheses("derivative_polynomial", "D", coefficients, n0, false);
            OperatorTransfer {
                theorem: LAGUERRE_DERIVATIVE.into(),
                hypotheses,
                from_nonpositive: Some(if one_sign {
                    nonpositive_ray()
                } else {
                    real_line()
                }),
                from_real_line: true,
                lean_route: None,
                lean_declarations: strings(&[
                    "RealRooted.splits_add_C_mul_derivative_all",
                    "RealRooted.splits_C_mul_add_derivative",
                    "Polynomial.Splits.derivative",
                ]),
            }
        }
        PreserverOperator::DerivativeOfPowerMultiple {
            power,
            derivative_order,
        }
        | PreserverOperator::PowerTimesDerivative {
            power,
            derivative_order,
        } => OperatorTransfer {
            theorem: ROLLE.into(),
            hypotheses: Vec::new(),
            from_nonpositive: Some(nonpositive_ray()),
            from_real_line: true,
            lean_route: (*power == 1
                && *derivative_order == 1
                && matches!(
                    operator,
                    PreserverOperator::DerivativeOfPowerMultiple { .. }
                ))
            .then(|| "derivativeXMul".into()),
            lean_declarations: strings(&[
                "RealRooted.DegreeDropReversal.splits_X_pow_mul_iff",
                "Polynomial.Splits.iterate_derivative",
            ]),
        },
        PreserverOperator::Affine { slope, intercept } => {
            let region = match (parse(slope), parse(intercept)) {
                (Some(a), Some(b)) if !a.is_zero() => {
                    let end = endpoint(&-(&b / &a));
                    if a.is_positive() {
                        RootInterval {
                            lower: None,
                            upper: Some(end),
                        }
                    } else {
                        RootInterval {
                            lower: Some(end),
                            upper: None,
                        }
                    }
                }
                _ => real_line(),
            };
            OperatorTransfer {
                theorem: SUBSTITUTION.into(),
                hypotheses: Vec::new(),
                from_nonpositive: Some(region),
                from_real_line: true,
                lean_route: None,
                lean_declarations: strings(&["Polynomial.Splits.comp_of_degree_le_one"]),
            }
        }
        PreserverOperator::Reversal { .. } => OperatorTransfer {
            theorem: SUBSTITUTION.into(),
            hypotheses: Vec::new(),
            from_nonpositive: Some(nonpositive_ray()),
            from_real_line: true,
            lean_route: Some("reversed".into()),
            lean_declarations: strings(&[
                "RealRooted.X_pow_mul_reverse_eq_zero_or_splits",
                "RealRooted.reflect_eq_zero_or_splits",
            ]),
        },
        PreserverOperator::Mobius { a, b, c, e, .. } => {
            let parsed = (parse(a), parse(b), parse(c), parse(e));
            let region = match &parsed {
                (Some(a), Some(b), Some(c), Some(e)) => mobius_region(a, b, c, e),
                _ => real_line(),
            };
            let f_polynomial = matches!(
                &parsed,
                (Some(a), Some(b), Some(c), Some(e))
                    if a.is_one() && b.is_zero() && c.is_one() && e.is_one()
            );
            OperatorTransfer {
                theorem: SUBSTITUTION.into(),
                hypotheses: Vec::new(),
                from_nonpositive: Some(region),
                from_real_line: true,
                lean_route: f_polynomial.then(|| "fTransform".into()),
                lean_declarations: if f_polynomial {
                    strings(&["RealRooted.splits_fPolynomial"])
                } else {
                    Vec::new()
                },
            }
        }
        PreserverOperator::HadamardMultiplier { factor } => {
            let name = "factor_multiplier_sequence";
            let statement = "k -> u(n, k) is a nonnegative multiplier sequence for every n >= n0";
            let (hypothesis, justified) = match factor.multiplier_justification() {
                Some(justification) => (verified(name, statement, justification), true),
                None => (
                    unresolved(
                        name,
                        statement,
                        format!(
                            "no multiplier-sequence theorem is recorded for {}",
                            factor.text()
                        ),
                    ),
                    false,
                ),
            };
            let validity = match factor.validate(n0) {
                Ok(()) => verified(
                    "factor_defined",
                    "the factor u(n, k) is defined on the index domain",
                    format!("{} is defined for every n >= {n0}", factor.text()),
                ),
                Err(error) => failed(
                    "factor_defined",
                    "the factor u(n, k) is defined on the index domain",
                    n0,
                    error.to_string(),
                ),
            };
            let declarations = factor
                .terms
                .iter()
                .filter_map(|term| match term {
                    crate::recurrence::HadamardFactorTerm::InverseFactorial => {
                        Some("RealRooted.isMultiplierSequence_inv_factorial")
                    }
                    crate::recurrence::HadamardFactorTerm::InverseRisingFactorial { .. } => {
                        Some("RealRooted.isPFMultiplierSequence_inv_ascPochhammer")
                    }
                    _ => None,
                })
                .map(String::from)
                .collect();
            OperatorTransfer {
                theorem: POLYA_SCHUR_THEOREM.into(),
                hypotheses: vec![validity, hypothesis],
                from_nonpositive: justified.then(nonpositive_ray),
                from_real_line: justified,
                lean_route: None,
                lean_declarations: declarations,
            }
        }
        PreserverOperator::EulerStep { diagonal, raising } => OperatorTransfer {
            theorem: EULER_STEP.into(),
            hypotheses: euler_step_hypotheses(diagonal, raising, n0),
            from_nonpositive: Some(nonpositive_ray()),
            from_real_line: false,
            lean_route: None,
            lean_declarations: strings(&[
                "RealRooted.IsPFPolynomial.generalStep_of_natDegree_le",
                "RealRooted.EulerBidiagonal.generalStep",
                "RealRooted.EulerBidiagonal.comparisonDefect",
            ]),
        },
        PreserverOperator::WeylOperator { .. } => OperatorTransfer {
            theorem: "none: a general Weyl-algebra operator need not preserve real-rootedness"
                .into(),
            hypotheses: vec![unresolved(
                "operator_preserves_real_rootedness",
                "the operator maps every parent row to a real-rooted polynomial or zero",
                "no transfer theorem is recorded for a general Weyl-algebra operator; finite \
                 row checks are not a proof"
                    .into(),
            )],
            from_nonpositive: None,
            from_real_line: false,
            lean_route: None,
            lean_declarations: Vec::new(),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    fn status(hypotheses: &[LiftHypothesis]) -> Vec<LiftHypothesisStatus> {
        hypotheses.iter().map(|h| h.status.clone()).collect()
    }

    #[test]
    fn fixed_polynomials_are_decided_exactly() {
        // theta (theta + 1): roots 0, -1.
        let good = PreserverOperator::euler_polynomial(&[q(0), q(1), q(1)]);
        let transfer = operator_transfer(&good, 0);
        assert_eq!(
            status(&transfer.hypotheses),
            vec![LiftHypothesisStatus::Verified]
        );
        // theta - 2: a positive root.
        let bad = PreserverOperator::euler_polynomial(&[q(-2), q(1)]);
        assert!(matches!(
            operator_transfer(&bad, 0).hypotheses[0].status,
            LiftHypothesisStatus::Failed { .. }
        ));
        // 1 + D^2 is not real-rooted.
        let bad = PreserverOperator::derivative_polynomial(&[q(1), q(0), q(1)]);
        assert!(matches!(
            operator_transfer(&bad, 0).hypotheses[0].status,
            LiftHypothesisStatus::Failed { .. }
        ));
        // 1 - 3 D + D^2 is real-rooted with mixed signs: the region is the line.
        let mixed = PreserverOperator::derivative_polynomial(&[q(1), q(-3), q(1)]);
        let transfer = operator_transfer(&mixed, 0);
        assert_eq!(
            transfer.hypotheses[0].status,
            LiftHypothesisStatus::Verified
        );
        assert_eq!(transfer.from_nonpositive, Some(real_line()));
    }

    #[test]
    fn index_dependent_polynomials_use_the_sign_checker() {
        let c = |num: &[i64], den: &[i64]| {
            IndexRationalFunction::from_polynomials(
                &num.iter().map(|v| q(*v)).collect::<Vec<_>>(),
                &den.iter().map(|v| q(*v)).collect::<Vec<_>>(),
            )
        };
        // theta + (n - 3): the root 3 - n is positive for n < 3.
        let op = PreserverOperator::EulerPolynomial {
            coefficients: vec![c(&[-3, 1], &[1]), c(&[1], &[1])],
        };
        let transfer = operator_transfer(&op, 0);
        let failing = transfer
            .hypotheses
            .iter()
            .find(|h| matches!(h.status, LiftHypothesisStatus::Failed { .. }))
            .expect("a failing hypothesis");
        assert_eq!(
            failing.status,
            LiftHypothesisStatus::Failed {
                witness: Some(BigInt::from(0))
            }
        );
        // From n0 = 3 on it holds.
        assert!(operator_transfer(&op, 3)
            .hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified));
        // 1 + (n/(n - 2)) D: real-rooted (linear) everywhere; the denominator
        // vanishes at n = 2.
        assert!(matches!(
            denominator_hypothesis(&[q(-2), q(1)], 0).status,
            LiftHypothesisStatus::Failed { .. }
        ));
        // theta^2 + n theta + n: discriminant n^2 - 4n < 0 at n = 1, 2, 3.
        let op = PreserverOperator::EulerPolynomial {
            coefficients: vec![c(&[0, 1], &[1]), c(&[0, 1], &[1]), c(&[1], &[1])],
        };
        let hypotheses = operator_transfer(&op, 1).hypotheses;
        assert_eq!(
            hypotheses[0].status,
            LiftHypothesisStatus::Failed {
                witness: Some(BigInt::from(1))
            }
        );
        assert!(operator_transfer(&op, 4)
            .hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified));
    }

    #[test]
    fn mobius_regions() {
        // The f-polynomial maps (-inf, 0] to [-1, 0].
        let region = mobius_region(&q(1), &q(0), &q(1), &q(1));
        assert_eq!(region.to_string(), "[-1, 0]");
        // A pole of the inverse inside (-inf, 0]: the whole line.
        assert_eq!(mobius_region(&q(-1), &q(0), &q(1), &q(1)), real_line());
    }
}
