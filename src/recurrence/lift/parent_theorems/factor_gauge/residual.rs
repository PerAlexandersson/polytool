//! Factor gauges, bounded residuals and palindromic hypotheses.

use super::*;

// ---------------------------------------------------------------------------
// The gauged rule and the scalar recurrence
// ---------------------------------------------------------------------------

/// The rule as `f(n) P_n = sum_j kappa_j(n) q(x)^j P_(n-j)`.
pub(super) struct GaugedRule {
    /// Monic, nonzero.
    pub(super) q: Vec<Q>,
    /// Nonzero, trimmed.
    pub(super) f: Vec<Q>,
    /// `kappa[j - 1]`, ascending in `n` (empty or `[0]` for a missing lag).
    pub(super) kappa: Vec<Vec<Q>>,
}

impl GaugedRule {
    pub(super) fn order(&self) -> usize {
        self.kappa.len()
    }

    /// Whether `f(n) pi(n) - sum_j kappa_j(n) pi(n - j)` is the zero
    /// polynomial in `n`.
    pub(super) fn annihilates(&self, pi: &[Q]) -> bool {
        let mut total = vec_poly::mul(&self.f, pi);
        for (lag, kappa) in (1i64..).zip(&self.kappa) {
            let shifted = vec_poly::compose(pi, &[Q::from_integer(BigInt::from(-lag)), Q::one()]);
            total = vec_poly::sub(&total, &vec_poly::mul(kappa, &shifted));
        }
        vec_poly::is_zero(&total)
    }

    /// The rows `S_(start), ..., S_(start + count - 1)` from the window
    /// `window = (S_start, ..., S_(start + L - 1))`; `None` where `f(n) = 0`.
    pub(super) fn residual_rows(
        &self,
        window: &[Vec<Q>],
        start: usize,
        count: usize,
    ) -> Option<Vec<Vec<Q>>> {
        let mut rows = window.to_vec();
        while rows.len() < count {
            let n = index_value(start + rows.len());
            let f = vec_poly::evaluate(&self.f, &n);
            if f.is_zero() {
                return None;
            }
            let mut sum = vec![Q::zero()];
            for (lag, kappa) in (1usize..).zip(&self.kappa) {
                let scalar = vec_poly::evaluate(kappa, &n);
                if !scalar.is_zero() {
                    sum = vec_poly::add(&sum, &vec_poly::scale(&rows[rows.len() - lag], &scalar));
                }
            }
            rows.push(vec_poly::trim(&vec_poly::scale(&sum, &f.recip())));
        }
        Some(rows)
    }

    /// The largest degree in `n` of `f` and the `kappa_j`.
    fn index_degree(&self) -> usize {
        std::iter::once(&self.f)
            .chain(&self.kappa)
            .filter_map(|p| vec_poly::degree(p))
            .max()
            .unwrap_or(0)
    }
}

pub(super) fn index_value(n: usize) -> Q {
    Q::from_integer(BigInt::from(n))
}

pub(super) fn negate_bivariate(b: &Bivariate) -> Bivariate {
    b.iter().map(|row| vec_poly::neg(row)).collect()
}

pub(super) fn parse_polynomial(coefficients: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    coefficients
        .iter()
        .map(|c| {
            parse_rational_coeff(c).map_err(|error| {
                RecurrenceError::VerificationFailed(format!("factor gauge: {error}"))
            })
        })
        .collect()
}

pub(super) fn format_polynomial(p: &[Q]) -> Vec<String> {
    vec_poly::trim(p)
        .iter()
        .map(format_rational_coeff)
        .collect()
}

/// The monic `q` with `q^power = f` for a monic `f`, when it exists:
/// `x^deg f f(1/x)` has constant term `1`, so its `power`-th root as a power
/// series is determined coefficient by coefficient; the candidate is
/// checked by raising it to the power.
pub(super) fn monic_polynomial_root(f: &[Q], power: usize) -> Option<Vec<Q>> {
    let f = vec_poly::trim(f);
    let degree = vec_poly::degree(&f)?;
    if power == 0 || degree % power != 0 || !f[degree].is_one() {
        return None;
    }
    let reversed: Vec<Q> = f.iter().rev().cloned().collect();
    let root_degree = degree / power;
    let exponent = index_value(power);
    let mut root = vec![Q::one()];
    for (k, target) in reversed.iter().enumerate().take(root_degree + 1).skip(1) {
        let mut product = vec![Q::one()];
        for _ in 0..power {
            product = vec_poly::mul(&product, &root);
            product.truncate(k + 1);
        }
        let current = vec_poly::coeff(&product, k);
        root.push((target.clone() - current) / &exponent);
    }
    let candidate: Vec<Q> = root.into_iter().rev().collect();
    (vec_poly::trim(&vec_poly::pow(&candidate, power)) == f).then_some(candidate)
}

/// Match the rule (module documentation), or say why it does not match.
pub(super) fn match_gauged_rule(shape: &RecurrenceShape) -> Result<GaugedRule, String> {
    if shape.terms.is_empty() {
        return Err("the rule has no nonzero right-hand term".into());
    }
    if shape.terms.keys().any(|&(_, derivative)| derivative != 0) {
        return Err("the rule has a derivative term".into());
    }
    let order = shape.terms.keys().map(|&(lag, _)| lag).max().unwrap_or(0);
    if order < 2 {
        return Err("a first-order rule is the product route's case (order 1)".into());
    }
    let mut separated = BTreeMap::new();
    for (&(lag, _), coefficients) in &shape.terms {
        let (gamma, w) = separate(coefficients).ok_or_else(|| {
            format!("the multiplier of P_(n-{lag}) is not a product gamma(n) w(x)")
        })?;
        separated.insert(lag, (gamma, w));
    }
    let (&smallest, (_, w)) = separated.iter().next().expect("a nonzero term");
    let q = monic_polynomial_root(&monic(w), smallest).ok_or_else(|| {
        format!(
            "the multiplier of P_(n-{smallest}) is not a constant times a {smallest}-th power of a \
             polynomial in x"
        )
    })?;
    let mut kappa = vec![Vec::new(); order];
    for (lag, (gamma, w)) in separated {
        if monic(&w) != vec_poly::trim(&vec_poly::pow(&q, lag)) {
            return Err(format!(
                "the multipliers are not kappa_j(n) q(x)^j for one q (q = {} from lag \
                 {smallest}, lag {lag} differs)",
                format_index_polynomial(&q, "x")
            ));
        }
        let leading = vec_poly::trim(&w).last().cloned().expect("nonzero");
        kappa[lag - 1] = vec_poly::scale(&gamma, &leading);
    }
    Ok(GaugedRule {
        q,
        f: vec_poly::trim(&shape.denominator),
        kappa,
    })
}

/// The shape of the space spanned by `window`, or `None` when it is not
/// one of the decided shapes (the zero polynomial is of every shape).
pub(super) fn residual_shape(window: &[Vec<Q>]) -> Option<FactorGaugeResidualShape> {
    let degree = window.iter().filter_map(|t| vec_poly::degree(t)).max()?;
    Some(match degree {
        0 => FactorGaugeResidualShape::Constant,
        1 => FactorGaugeResidualShape::Linear,
        2 => FactorGaugeResidualShape::Quadratic,
        3 => FactorGaugeResidualShape::Cubic,
        4 if window.iter().all(|t| {
            let c = |k| vec_poly::coeff(t, k);
            c(0) == c(4) && c(1) == c(3)
        }) =>
        {
            FactorGaugeResidualShape::PalindromicQuartic
        }
        _ => return None,
    })
}

/// The coefficient polynomials `pi_k(n)`, `k = 0, ..., degree`, of the
/// residual rows (module documentation), or `None`.
pub(super) fn polynomial_coordinates(
    rule: &GaugedRule,
    window: &[Vec<Q>],
    start: usize,
    degree: usize,
) -> Option<Vec<Vec<Q>>> {
    let order = rule.order();
    let bound = (order - 1 + 2 * rule.index_degree() + 2).min(MAX_COORDINATE_DEGREE);
    let count = bound + order + 2;
    let rows = rule.residual_rows(window, start, count)?;
    let points: Vec<Q> = (0..count).map(|t| index_value(start + t)).collect();
    let mut coefficients = Vec::with_capacity(degree + 1);
    for k in 0..=degree {
        let values: Vec<Q> = rows.iter().map(|row| vec_poly::coeff(row, k)).collect();
        let pi = (0..=bound).find_map(|d| {
            let candidate = Polynomial::lagrange_interpolation(&points[..=d], &values[..=d])
                .coeffs()
                .to_vec();
            points
                .iter()
                .zip(&values)
                .all(|(point, value)| &vec_poly::evaluate(&candidate, point) == value)
                .then_some(candidate)
        })?;
        if !rule.annihilates(&pi) {
            return None;
        }
        coefficients.push(vec_poly::trim(&pi));
    }
    Some(coefficients)
}

// ---------------------------------------------------------------------------
// Root sides
// ---------------------------------------------------------------------------

/// Where the roots of a real-rooted polynomial lie, read from its
/// coefficient signs: all coefficients nonzero of one sign (every root
/// negative), the nonzero ones of one sign (every root `<= 0`), or neither.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub(super) enum RootSide {
    Negative,
    NonPositive,
    Other,
}

pub(super) fn fixed_root_side(p: &[Q]) -> RootSide {
    let p = vec_poly::trim(p);
    let positive = p.iter().all(|c| c.is_positive());
    let negative = p.iter().all(|c| c.is_negative());
    if positive || negative {
        return RootSide::Negative;
    }
    let nonzero: Vec<&Q> = p.iter().filter(|c| !c.is_zero()).collect();
    if nonzero.iter().all(|c| c.is_positive()) || nonzero.iter().all(|c| c.is_negative()) {
        RootSide::NonPositive
    } else {
        RootSide::Other
    }
}

/// The root side of `S_n` for every `n >= from` from its coefficient
/// polynomials, with the (verified) sign hypotheses that give it.
pub(super) fn residual_root_side(
    coefficients: &[Vec<Q>],
    from: usize,
) -> (RootSide, Vec<LiftHypothesis>) {
    let attempts = [
        (IndexSignCondition::Positive, "positive", RootSide::Negative),
        (IndexSignCondition::Negative, "negative", RootSide::Negative),
        (
            IndexSignCondition::Nonnegative,
            "nonnegative",
            RootSide::NonPositive,
        ),
        (
            IndexSignCondition::Nonpositive,
            "nonpositive",
            RootSide::NonPositive,
        ),
    ];
    for (condition, word, side) in attempts {
        let mut hypotheses = Vec::new();
        for (k, pi) in coefficients.iter().enumerate() {
            let function = IndexFunction::polynomial(pi.clone());
            let hypothesis = sign_hypothesis(
                &format!("residual_coefficient_{k}_{word}"),
                &format!("[x^{k}] S_n is {word} for every n >= {from}"),
                &format!("[x^{k}] S_n"),
                &function,
                from,
                condition,
            );
            if hypothesis.status != LiftHypothesisStatus::Verified {
                break;
            }
            hypotheses.push(hypothesis);
        }
        if hypotheses.len() == coefficients.len() {
            return (side, hypotheses);
        }
    }
    (RootSide::Other, Vec::new())
}

// ---------------------------------------------------------------------------
// The report
// ---------------------------------------------------------------------------

/// Which palindromic case proved the quartic (for the Lean declarations).
#[derive(Clone, Copy, PartialEq, Eq)]
pub(super) enum PalindromicCase {
    BelowMinusTwo,
    AboveTwo,
    Straddle,
}

/// The real-rootedness hypotheses of `S_n` for the shape (module
/// documentation), and the palindromic case used.
pub(super) fn residual_hypotheses(
    shape: FactorGaugeResidualShape,
    coefficients: Option<&[Vec<Q>]>,
    from: usize,
) -> (Vec<LiftHypothesis>, Option<PalindromicCase>) {
    let p = |pi: &[Q], qi: &[Q]| vec_poly::mul(pi, qi);
    let s = |c: i64, pi: &[Q]| vec_poly::scale(pi, &Q::from_integer(BigInt::from(c)));
    let sum = |parts: &[Vec<Q>]| {
        parts
            .iter()
            .fold(vec![Q::zero()], |total, part| vec_poly::add(&total, part))
    };
    let degree_at_most_one = || {
        vec![exact_hypothesis(
            "residual_degree_at_most_one",
            "every element of the residual space has degree at most 1, so every nonzero S_n is \
             real-rooted",
            Ok("exact: every residual window row has degree at most 1".into()),
        )]
    };
    let c = match (shape, coefficients) {
        (FactorGaugeResidualShape::Constant | FactorGaugeResidualShape::Linear, _) => {
            return (degree_at_most_one(), None);
        }
        (_, None) => {
            return (
                vec![unresolved_hypothesis(
                    "residual_real_rooted",
                    "S_n is real-rooted for every n >= N",
                    "the coefficients of S_n are not polynomials in n".into(),
                )],
                None,
            );
        }
        (_, Some(c)) => c,
    };
    match shape {
        FactorGaugeResidualShape::Constant | FactorGaugeResidualShape::Linear => {
            (degree_at_most_one(), None)
        }
        FactorGaugeResidualShape::Quadratic => {
            // B^2 - 4 A C with S_n = A x^2 + B x + C.
            let discriminant = sum(&[p(&c[1], &c[1]), s(-4, &p(&c[2], &c[0]))]);
            (
                vec![sign_hypothesis(
                    "residual_discriminant_nonnegative",
                    &format!(
                        "the discriminant B(n)^2 - 4 A(n) C(n) of S_n = A(n) x^2 + B(n) x + C(n) \
                         is >= 0 for every n >= {from} (for A(n) = 0 it is B(n)^2 and S_n has \
                         degree at most 1)"
                    ),
                    "B(n)^2 - 4 A(n) C(n)",
                    &IndexFunction::polynomial(discriminant),
                    from,
                    IndexSignCondition::Nonnegative,
                )],
                None,
            )
        }
        FactorGaugeResidualShape::Cubic => {
            let (a, b, cc, d) = (&c[3], &c[2], &c[1], &c[0]);
            let discriminant = sum(&[
                s(18, &p(&p(a, b), &p(cc, d))),
                s(-4, &p(&p(&p(b, b), b), d)),
                p(&p(b, b), &p(cc, cc)),
                s(-4, &p(&p(a, cc), &p(cc, cc))),
                s(-27, &p(&p(a, a), &p(d, d))),
            ]);
            (
                vec![sign_hypothesis(
                    "residual_cubic_discriminant_nonnegative",
                    &format!(
                        "the formal discriminant 18abcd - 4b^3 d + b^2 c^2 - 4ac^3 - 27a^2 d^2 of \
                         S_n = a x^3 + b x^2 + c x + d is >= 0 for every n >= {from} (for a = 0 \
                         it is b^2 (c^2 - 4bd))"
                    ),
                    "the cubic discriminant",
                    &IndexFunction::polynomial(discriminant),
                    from,
                    IndexSignCondition::Nonnegative,
                )],
                None,
            )
        }
        FactorGaugeResidualShape::PalindromicQuartic => {
            palindromic_hypotheses(&c[0], &c[1], &c[2], from)
        }
    }
}

fn palindromic_hypotheses(
    a: &[Q],
    b: &[Q],
    c: &[Q],
    from: usize,
) -> (Vec<LiftHypothesis>, Option<PalindromicCase>) {
    let outer = |condition, word: &str| {
        sign_hypothesis(
            "residual_outer_coefficient_sign",
            &format!(
                "a(n) = [x^0] S_n = [x^4] S_n is {word} for every n >= {from} (S_n = a x^4 + b \
                 x^3 + c x^2 + b x + a)"
            ),
            "a(n)",
            &IndexFunction::polynomial(a.to_vec()),
            from,
            condition,
        )
    };
    let positive = outer(IndexSignCondition::Positive, "positive");
    let (sign, outer_hypothesis) = if positive.status == LiftHypothesisStatus::Verified {
        (1, positive)
    } else {
        let negative = outer(IndexSignCondition::Negative, "negative");
        if negative.status == LiftHypothesisStatus::Verified {
            (-1, negative)
        } else {
            return (vec![positive], None);
        }
    };
    let scale = |p: &[Q], factor: i64| vec_poly::scale(p, &Q::from_integer(BigInt::from(factor)));
    let (a, b, c) = (scale(a, sign), scale(b, sign), scale(c, sign));
    let combination = |terms: &[(i64, &Vec<Q>)]| {
        terms.iter().fold(vec![Q::zero()], |total, (factor, p)| {
            vec_poly::add(&total, &scale(p, *factor))
        })
    };
    // b^2 - 4 a (c - 2a), h(2) = 2a + 2b + c, h(-2) = 2a - 2b + c.
    let discriminant = combination(&[
        (1, &vec_poly::mul(&b, &b)),
        (-4, &vec_poly::mul(&a, &c)),
        (8, &vec_poly::mul(&a, &a)),
    ]);
    let at_two = combination(&[(2, &a), (2, &b), (1, &c)]);
    let at_minus_two = combination(&[(2, &a), (-2, &b), (1, &c)]);
    let normalized = if sign == 1 {
        ""
    } else {
        " (after S_n -> -S_n)"
    };
    let vertex_below = combination(&[(1, &b), (-4, &a)]);
    let vertex_above = combination(&[(-1, &b), (-4, &a)]);
    let cases: [(PalindromicCase, Vec<PalindromicCondition<'_>>); 3] = [
        (
            PalindromicCase::BelowMinusTwo,
            vec![
                (
                    "residual_u_vertex_at_most_minus_two",
                    "b - 4a >= 0 (the vertex of h is at most -2)",
                    &vertex_below,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_value_at_minus_two_nonnegative",
                    "h(-2) = 2a - 2b + c >= 0",
                    &at_minus_two,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_discriminant_nonnegative",
                    "b^2 - 4a(c - 2a) >= 0 (h has real roots)",
                    &discriminant,
                    IndexSignCondition::Nonnegative,
                ),
            ],
        ),
        (
            PalindromicCase::AboveTwo,
            vec![
                (
                    "residual_u_vertex_at_least_two",
                    "-b - 4a >= 0 (the vertex of h is at least 2)",
                    &vertex_above,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_value_at_two_nonnegative",
                    "h(2) = 2a + 2b + c >= 0",
                    &at_two,
                    IndexSignCondition::Nonnegative,
                ),
                (
                    "residual_u_discriminant_nonnegative",
                    "b^2 - 4a(c - 2a) >= 0 (h has real roots)",
                    &discriminant,
                    IndexSignCondition::Nonnegative,
                ),
            ],
        ),
        (
            PalindromicCase::Straddle,
            vec![
                (
                    "residual_u_value_at_two_nonpositive",
                    "h(2) = 2a + 2b + c <= 0",
                    &at_two,
                    IndexSignCondition::Nonpositive,
                ),
                (
                    "residual_u_value_at_minus_two_nonpositive",
                    "h(-2) = 2a - 2b + c <= 0",
                    &at_minus_two,
                    IndexSignCondition::Nonpositive,
                ),
            ],
        ),
    ];
    let mut evaluated = Vec::new();
    for (case, conditions) in &cases {
        let hypotheses: Vec<LiftHypothesis> = conditions
            .iter()
            .map(|(name, statement, polynomial, condition)| {
                sign_hypothesis(
                    name,
                    &format!(
                        "{statement} for every n >= {from}, with h(u) = a u^2 + b u + (c - 2a), \
                         u = x + 1/x{normalized}"
                    ),
                    name,
                    &IndexFunction::polynomial((*polynomial).clone()),
                    from,
                    *condition,
                )
            })
            .collect();
        if hypotheses
            .iter()
            .all(|h| h.status == LiftHypothesisStatus::Verified)
        {
            let mut all = vec![outer_hypothesis];
            all.extend(hypotheses);
            return (all, Some(*case));
        }
        evaluated.push(hypotheses);
    }
    // No case holds: report the case on the side of the vertex of h at the
    // first index (b >= 0: both roots of h at most -2).
    let b_first = IndexFunction::polynomial(b.clone())
        .at(from)
        .unwrap_or_else(Q::zero);
    let reported = usize::from(b_first.is_negative());
    let mut all = vec![outer_hypothesis];
    all.extend(evaluated.swap_remove(reported));
    (all, None)
}

/// One sign condition of a palindromic case: name, statement, polynomial
/// in `n` and condition.
type PalindromicCondition<'a> = (&'static str, &'static str, &'a Vec<Q>, IndexSignCondition);
