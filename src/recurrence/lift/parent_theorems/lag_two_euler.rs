//! Positive lag-two Euler recurrences, with the simultaneous adjacent and
//! auxiliary interlacing induction of RealRooted's `LagTwoRecurrence`.
//!
//! Put `m = n - first_index` and `theta = x D`. The matched operators are
//! `(a + b theta) + x(a + b(m-1) - b theta)` at lag one and
//! `c x (theta+1)(m-1-theta)` at lag two. We convert `x^d D^d` to the
//! falling factorial in `theta` and compare every coefficient in `Q[n]`.
//! The scalar functions are divided by the rule's original denominator;
//! its nonvanishing and the signs `a > 0`, `b-a >= 0`, `c >= 0` are decided
//! on the entire integer rule range, never sampled.
//!
//! The initial rows are a positive constant and a positive linear row.
//! Additional seeds (including a moved start) must themselves be positive
//! Euler steps with nonnegative lag correction. Their scalars are solved
//! exactly, independently of the rule before its start. This supplies the
//! hypotheses of the same Lean recurrence theorem at every seed index.

use super::*;
use crate::recurrence::IndexRationalFunction;

const REFERENCE: &str = "Positive lag-two Euler recurrence: a_n > 0, b_n >= a_n, c_n >= 0, \
    a positive constant and positive linear seed give simple negative roots, adjacent and \
    auxiliary strict interlacing (RealRooted/EulerianTriangles/LagTwoRecurrence.lean, \
    lagTwoEuler_recurrence_certificate)";

/// One additional seed as an exact step of the lag-two Euler recurrence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LagTwoEulerSeedStep {
    /// Absolute index of the output row (relative degree `index - first_index`).
    pub index: usize,
    /// Positive successor scalar `a`, as an exact rational string.
    pub a: String,
    /// Successor scalar `b >= a`, as an exact rational string.
    pub b: String,
    /// Lag scalar `c >= 0`, as an exact rational string.
    pub c: String,
}

/// Replayable operator identity and exact initial rows for
/// [`check_lag_two_euler_theorem`]. Infinite-range sign certificates are
/// kept in the enclosing report's hypotheses.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct LagTwoEulerCertificate {
    /// Normalized rational function `a(n)` on the rule range.
    pub a: IndexRationalFunction,
    /// Normalized rational function `b(n)` on the rule range.
    pub b: IndexRationalFunction,
    /// Normalized rational function `c(n)` on the rule range.
    pub c: IndexRationalFunction,
    /// Positive constant initial row `w`.
    pub first_row: String,
    /// Positive constant coefficient `u` of the linear initial row.
    pub linear_constant: String,
    /// Positive coefficient `v` of `x` in the linear initial row.
    pub linear_slope: String,
    /// Every seed after the first two, in increasing absolute index order.
    pub seed_steps: Vec<LagTwoEulerSeedStep>,
}

impl LagTwoEulerCertificate {
    /// Verify the normalized operator identity in `Q(n)[x, D]`, using the
    /// degree origin `first_index`. Does not assert the coefficient signs.
    ///
    /// # Errors
    ///
    /// Invalid rational functions or an operator mismatch.
    pub fn verify_operator_identity(
        &self,
        rule: &Recurrence,
        first_index: usize,
    ) -> Result<(), RecurrenceError> {
        let (rule, _) = rule_with_x_free_left_hand_factor(rule)?;
        let expected =
            recognize(&rule, first_index).map_err(RecurrenceError::VerificationFailed)?;
        for (recorded, expected) in [&self.a, &self.b, &self.c].into_iter().zip(expected) {
            // The representation is shared with preservers, whose parser
            // bounds the degree for their search. This theorem places no
            // such bound on its scalar functions of the index.
            let parse = |coefficients: &[String]| {
                coefficients
                    .iter()
                    .map(|c| parse_rational_coeff(c).map_err(RecurrenceError::from))
                    .collect::<Result<Vec<_>, _>>()
            };
            let numerator = parse(&recorded.numerator)?;
            let denominator = parse(&recorded.denominator)?;
            if vec_poly::is_zero(&denominator) {
                return Err(RecurrenceError::VerificationFailed(
                    "lag-two Euler coefficient has a zero denominator".into(),
                ));
            }
            if vec_poly::trim(&vec_poly::mul(&numerator, &expected.denominator))
                != vec_poly::trim(&vec_poly::mul(&expected.numerator, &denominator))
            {
                return Err(RecurrenceError::VerificationFailed(
                    "lag-two Euler operator coefficients differ".into(),
                ));
            }
        }
        Ok(())
    }

    /// Replay the positive constant/linear seeds and every recorded later
    /// seed step, including all three scalar inequalities.
    ///
    /// # Errors
    ///
    /// Invalid scalars, wrong indices, missing steps or unequal rows.
    pub fn verify_seed_steps(
        &self,
        first_index: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        let fail = || RecurrenceError::VerificationFailed("lag-two Euler seed mismatch".into());
        let w = parse_rational_coeff(&self.first_row)?;
        let u = parse_rational_coeff(&self.linear_constant)?;
        let v = parse_rational_coeff(&self.linear_slope)?;
        if seeds.len() < 2
            || self.seed_steps.len() != seeds.len() - 2
            || !w.is_positive()
            || !u.is_positive()
            || !v.is_positive()
            || vec_poly::trim(&seeds[0]) != [w]
            || vec_poly::trim(&seeds[1]) != [u, v]
        {
            return Err(fail());
        }
        for (m, step) in (2..seeds.len()).zip(&self.seed_steps) {
            let a = parse_rational_coeff(&step.a)?;
            let b = parse_rational_coeff(&step.b)?;
            let c = parse_rational_coeff(&step.c)?;
            if step.index != first_index + m
                || !a.is_positive()
                || b < a
                || c.is_negative()
                || vec_poly::trim(&seeds[m])
                    != seed_successor(&seeds[m - 1], &seeds[m - 2], m, &a, &b, &c)
            {
                return Err(fail());
            }
        }
        Ok(())
    }
}

/// Check RealRooted's positive lag-two Euler recurrence theorem.
///
/// With `m = n - first_index` and `theta = x D`, the homogeneous rule must
/// be `P_n = [(a+b theta)+x(a+b(m-1)-b theta)]P_(n-1) +
/// c x(theta+1)(m-1-theta)P_(n-2)`. For every integer on the rule range,
/// the scalar denominator must be nonzero, `a > 0`, `b >= a`, `c >= 0`.
/// Two initial rows must be `w` and `u+v x` with `w,u,v > 0`. Extra seeds
/// must be exact steps with these inequalities (their scalars may differ
/// from the rule before its start). The usual maximum start delay applies.
///
/// Every row has its full relative degree of simple negative roots;
/// consecutive rows strictly interlace, as do the auxiliary rows
/// `(theta+1)(m-theta)P_(n-1)` and `P_n`. A requested region must contain
/// `(-inf, 0)`. The negative lag correction is excluded.
///
/// # Errors
///
/// Invalid recurrence records, seed normalization or moved row generation.
pub fn check_lag_two_euler_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::LagTwoEuler,
        rule,
        first_index,
        seeds,
        options,
    )
}

pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::LagTwoEuler,
    check_at_supplied_seeds: lag_two_euler_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: Some("lag-two Euler"),
    discovery: DiscoveryPolicy::WhenApplicable,
    in_residue_class_menu: true,
    in_index_section_menu: true,
    parent: ParentPolicy {
        stage: ParentStage::WhileUnproved,
        region: ParentRegion::OpenNonpositiveRay,
        record_only_if_applicable: true,
    },
};

/// Theta coefficients keyed by `(lag, x power, theta power)`, polynomial in n.
type ThetaOperator = BTreeMap<(usize, usize, usize), Vec<Q>>;

fn put(operator: &mut ThetaOperator, key: (usize, usize, usize), p: Vec<Q>) {
    let combined = vec_poly::trim(&vec_poly::add(operator.entry(key).or_default(), &p));
    if vec_poly::is_zero(&combined) {
        operator.remove(&key);
    } else {
        operator.insert(key, combined);
    }
}

fn recognize(rule: &Recurrence, first: usize) -> Result<[IndexFunction; 3], String> {
    let shape = recurrence_shape(rule)?;
    let form = "the operators are not the lag-two Euler theta form";
    let mut actual = ThetaOperator::new();
    // x^0 D^0 = 1, x D = theta, x^2 D^2 = theta^2 - theta.
    // Higher derivatives cannot occur in an operator of this form.
    for (&(lag, derivative), coefficient) in &shape.terms {
        if !(1..=2).contains(&lag) || derivative > if lag == 1 { 1 } else { 2 } {
            return Err(form.into());
        }
        for power in 0..x_width(coefficient) {
            let p = x_column(coefficient, power);
            if p.is_empty() {
                continue;
            }
            if power < derivative {
                return Err(form.into());
            }
            if derivative == 2 {
                put(&mut actual, (lag, power - 2, 1), vec_poly::neg(&p));
            }
            put(&mut actual, (lag, power - derivative, derivative), p);
        }
    }
    let get = |key| actual.get(&key).cloned().unwrap_or_default();
    let a = get((1, 0, 0));
    let b = get((1, 0, 1));
    let c = vec_poly::neg(&get((2, 1, 2)));
    let degree = vec![-Q::from_integer(BigInt::from(first) + 1), Q::one()];
    let mut expected = ThetaOperator::new();
    put(&mut expected, (1, 0, 0), a.clone());
    put(&mut expected, (1, 0, 1), b.clone());
    put(
        &mut expected,
        (1, 1, 0),
        vec_poly::add(&a, &vec_poly::mul(&b, &degree)),
    );
    put(&mut expected, (1, 1, 1), vec_poly::neg(&b));
    put(&mut expected, (2, 1, 0), vec_poly::mul(&c, &degree));
    put(
        &mut expected,
        (2, 1, 1),
        vec_poly::mul(&c, &vec_poly::sub(&degree, &[Q::one()])),
    );
    put(&mut expected, (2, 1, 2), vec_poly::neg(&c));
    if actual != expected {
        return Err(form.into());
    }
    Ok([a, b, c].map(|p| IndexFunction::new(p, shape.denominator.clone())))
}

fn seed_successor(previous: &[Q], delayed: &[Q], m: usize, a: &Q, b: &Q, c: &Q) -> Vec<Q> {
    let mut result = vec![Q::zero(); previous.len().max(delayed.len()) + 1];
    let degree = Q::from_integer(BigInt::from(m) - 1);
    for (index, coefficient) in previous.iter().enumerate() {
        let k = Q::from_integer(BigInt::from(index));
        result[index] += (a + b * &k) * coefficient;
        result[index + 1] += (a + b * (&degree - k)) * coefficient;
    }
    for (k, coefficient) in delayed.iter().enumerate() {
        let index = Q::from_integer(BigInt::from(k));
        result[k + 1] += c * (&index + Q::one()) * (&degree - index) * coefficient;
    }
    vec_poly::trim(&result)
}

/// Solve `target = s left + c right`, `s,c >= 0`, in dimension at most two.
fn nonnegative_pair(target: &[Q], left: &[Q], right: &[Q]) -> Option<(Q, Q)> {
    let width = target.len().max(left.len()).max(right.len());
    let at = |p: &[Q], i| p.get(i).cloned().unwrap_or_else(Q::zero);
    let accepts = |s: Q, c: Q| {
        (!s.is_negative()
            && !c.is_negative()
            && vec_poly::trim(target)
                == vec_poly::trim(&vec_poly::add(
                    &vec_poly::scale(left, &s),
                    &vec_poly::scale(right, &c),
                )))
        .then_some((s, c))
    };
    let pivot = (0..width).find(|&i| !at(left, i).is_zero() || !at(right, i).is_zero());
    let Some(i) = pivot else {
        return accepts(Q::zero(), Q::zero());
    };
    for j in i + 1..width {
        let determinant = at(left, i) * at(right, j) - at(left, j) * at(right, i);
        if !determinant.is_zero() {
            let s = (at(target, i) * at(right, j) - at(target, j) * at(right, i)) / &determinant;
            let c = (at(left, i) * at(target, j) - at(left, j) * at(target, i)) / determinant;
            return accepts(s, c);
        }
    }
    if !at(left, i).is_zero() {
        if let Some(pair) = accepts(at(target, i) / at(left, i), Q::zero()) {
            return Some(pair);
        }
    }
    if !at(right, i).is_zero() {
        return accepts(Q::zero(), at(target, i) / at(right, i));
    }
    None
}

fn seed_step(seed: &[Q], previous: &[Q], delayed: &[Q], m: usize) -> Option<(Q, Q, Q)> {
    let previous_constant = previous.first().filter(|c| c.is_positive())?;
    let a = seed.first()?.clone() / previous_constant;
    if !a.is_positive() {
        return None;
    }
    let base = seed_successor(previous, delayed, m, &a, &a, &Q::zero());
    let left = seed_successor(previous, delayed, m, &Q::zero(), &Q::one(), &Q::zero());
    let right = seed_successor(previous, delayed, m, &Q::zero(), &Q::zero(), &Q::one());
    let (surplus, c) = nonnegative_pair(&vec_poly::sub(seed, &base), &left, &right)?;
    Some((a.clone(), a + surplus, c))
}

pub(super) fn verify_certificate_in_report(
    certificate: &LagTwoEulerCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    let (rule, first, seeds) = report.recurrence.to_recurrence_parts()?;
    certificate.verify_operator_identity(&rule, first)?;
    if report.is_proved() {
        let seeds = rule
            .generate_rows_rational(&seeds, first, seeds.len() + report.start_delay)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        certificate.verify_seed_steps(first, &seeds)?;
    }
    Ok(())
}

fn lag_two_euler_report(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::LagTwoEuler,
        REFERENCE,
        rule,
        first,
        seeds,
        options,
    );
    let [a, b, c] = match recognize(rule, first) {
        Ok(parameters) => parameters,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let start = first + seeds.len();
    let shape = recurrence_shape(rule).expect("recognition checked the shape");
    report.normal_form = format!("m = n - {first}; P_n = [(a+b theta)+x(a+b(m-1)-b theta)] P_(n-1) + c x(theta+1)(m-1-theta) P_(n-2), n >= {start}; a(n) = {}, b(n) = {}, c(n) = {}", a.text(), b.text(), c.text());
    report.lean_route = Some("lagTwoEuler".into());
    report.lean_declarations = [
        "RealRooted.lagTwoEuler_recurrence_certificate",
        "RealRooted.lagTwoEuler_recurrence_simpleNegRooted",
    ]
    .map(str::to_string)
    .to_vec();
    let mut hypotheses: Vec<_> = denominator_hypothesis(&shape.denominator, start)
        .into_iter()
        .collect();
    hypotheses.push(sign_hypothesis(
        "successor_constant_positive",
        &format!("a(n) > 0 for every n >= {start}"),
        "a(n)",
        &a,
        start,
        IndexSignCondition::Positive,
    ));
    hypotheses.push(sign_hypothesis(
        "successor_theta_dominates_constant",
        &format!("b(n) - a(n) >= 0 for every n >= {start}"),
        "b(n) - a(n)",
        &b.sub(&a),
        start,
        IndexSignCondition::Nonnegative,
    ));
    hypotheses.push(sign_hypothesis(
        "lag_coefficient_nonnegative",
        &format!("c(n) >= 0 for every n >= {start}"),
        "c(n)",
        &c,
        start,
        IndexSignCondition::Nonnegative,
    ));
    let initial = seeds.first().map(|p| vec_poly::trim(p));
    let linear = seeds.get(1).map(|p| vec_poly::trim(p));
    let w = initial
        .as_ref()
        .filter(|p| p.len() == 1 && p[0].is_positive())
        .map(|p| p[0].clone());
    let uv = linear
        .as_ref()
        .filter(|p| p.len() == 2 && p.iter().all(Signed::is_positive))
        .map(|p| (p[0].clone(), p[1].clone()));
    hypotheses.push(exact_hypothesis(
        "positive_constant_seed",
        "the first row is w > 0",
        w.as_ref()
            .map(|w| format!("w = {}", format_rational_coeff(w)))
            .ok_or((
                Some(first),
                "the first row is not a positive constant".into(),
            )),
    ));
    hypotheses.push(exact_hypothesis(
        "positive_linear_seed",
        "the second row is u + v x with u,v > 0",
        uv.as_ref()
            .map(|(u, v)| {
                format!(
                    "u = {}, v = {}",
                    format_rational_coeff(u),
                    format_rational_coeff(v)
                )
            })
            .ok_or((
                Some(first + 1),
                "the second row is not positive linear (or is missing)".into(),
            )),
    ));
    let mut steps = Vec::new();
    let mut failed = None;
    for m in 2..seeds.len() {
        match seed_step(&seeds[m], &seeds[m - 1], &seeds[m - 2], m) {
            Some((a, b, c)) => steps.push(LagTwoEulerSeedStep {
                index: first + m,
                a: format_rational_coeff(&a),
                b: format_rational_coeff(&b),
                c: format_rational_coeff(&c),
            }),
            None => {
                failed = Some(first + m);
                break;
            }
        }
    }
    if seeds.len() > 2 {
        hypotheses.push(exact_hypothesis(
            "seed_euler_steps",
            "every later seed is an exact Euler step with a > 0, b >= a, c >= 0",
            match failed {
                None => Ok(format!(
                    "{} exact seed steps, scalars recorded in the certificate",
                    steps.len()
                )),
                Some(index) => Err((
                    Some(index),
                    format!("P_{index} has no Euler seed step with a > 0, b >= a, c >= 0"),
                )),
            },
        ));
    }
    let negative_axis = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: false,
    };
    if let Some(region) = &options.root_region {
        let contains = region.lower.is_none()
            && region
                .upper
                .as_ref()
                .is_none_or(|upper| !upper.is_negative());
        hypotheses.push(if contains {
            exact_hypothesis(
                "requested_root_region",
                "the requested region contains (-inf, 0)",
                Ok("exact".into()),
            )
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                "the requested region contains (-inf, 0)",
                "the theorem certifies (-inf, 0) only".into(),
            )
        });
    }
    report.hypotheses = hypotheses;
    if let (Some(w), Some((u, v))) = (w, uv) {
        let encode = |f: &IndexFunction| {
            IndexRationalFunction::from_polynomials(&f.numerator, &f.denominator)
        };
        report.certificate = Some(RouteCertificate::LagTwoEuler(LagTwoEulerCertificate {
            a: encode(&a),
            b: encode(&b),
            c: encode(&c),
            first_row: format_rational_coeff(&w),
            linear_constant: format_rational_coeff(&u),
            linear_slope: format_rational_coeff(&v),
            seed_steps: steps,
        }));
    }
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: first,
        zero_residue_classes: None,
        interlacing_from: Some(first + 1),
        interlacing: Some(RowInterlacing::Strict),
        simple_roots: true,
        root_region: Some(open_root_interval(&negative_axis)),
        statement: format!("every P_n, n >= {first}, has n - {first} simple negative roots and positive coefficients; P_(n-1) strictly interlaces P_n for n > {first}, and (theta+1)(n-{first}-theta) P_(n-1) StrictInterl P_n (auxiliary interlacing)"),
    }))
}
