//! Fixed-delay recurrences `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)`,
//! `L >= 2`: the RealRooted row route `variableDelay`.
//!
//! # The theorem
//!
//! RealRooted (`RealRooted/Interlacing/VariableDelay.lean`), after Fisk,
//! *Polynomials, roots, and interlacing*, Example 16.10, with the cone rules
//! of Wagner (1992), Section 3: let `Q_0 = 1` and, for every `m >= 0`,
//!
//! ```text
//! Q_(m+1) = a(m+1) Q_m + b(m+1) x Q_(max(m+1-L, 0))
//! ```
//!
//! with `a > 0` and `b >= 0` (`L >= 1`).  Then every `Q_m` has nonnegative
//! coefficients and constant term `1`
//! (`RealRooted.hasNonnegCoeffs_of_scaled_fixedDelay_recurrence`), is
//! real-rooted, and `Q_A << Q_B` whenever `A <= B <= A + L`
//! (`RealRooted.strictInterl_of_scaled_fixedDelay_recurrence`; for `a = 1`,
//! `RealRooted.strictInterl_of_fixedDelay_recurrence`).  The relation
//! `f << g` (`RealRooted.StrictInterl`) means: both are nonzero and
//! real-rooted, `deg g` is `deg f` or `deg f + 1`, and the roots alternate
//! weakly with the largest root belonging to `g`, which is
//! [`RowInterlacing::Weak`].  Constant term `1` and nonnegative coefficients
//! put every root in `(-inf, 0)`.
//!
//! # Matching a polytool family
//!
//! The family is `P_n`, `n >= n_0` (the first index), with `K` seed rows; the
//! rule applies from `N = n_0 + K` on.  Put `Q_m = P_(n_0 + m) / P_(n_0)`.
//! The scalars of Lean's recurrence are free functions of the index, so the
//! hypotheses split into
//!
//! - the rule range `n >= N`: `a(n) = A(n) / f(n) > 0` and
//!   `b(n) = B(n) / f(n) >= 0`, decided exactly for every integer `n >= N` by
//!   [`crate::decide_index_rational_function_sign`] (with `f(n) != 0`), and
//!   `K >= L`, so that the delayed index `n - L >= n_0` is never truncated;
//! - the seed range: `P_(n_0)` is a nonzero constant, and every later seed
//!   row is a step of the window, `P_(n_0+m) = a_m P_(n_0+m-1) +
//!   b_m x P_(n_0 + max(m-L, 0))` with rational `a_m > 0`, `b_m >= 0`, solved
//!   exactly and recorded in [`VariableDelayCertificate`].
//!
//! `L = 1` is excluded: the rule is then the product
//! `P_n = (a(n) + b(n) x) P_(n-1)`, the product route's case.
//!
//! Related: the paper-theorem route `fixed_delay_predecessor_window` of
//! [`crate::renewal::match_recurrence_theorems`] (`a = 1`, a polynomial
//! weight, canonical seed rows only).

use super::*;

const REFERENCE: &str = "Fixed-delay window (Fisk, Polynomials, roots, and interlacing, \
     Example 16.10; Wagner 1992, Section 3): Q_0 = 1 and Q_(m+1) = a(m+1) Q_m + b(m+1) x \
     Q_(max(m+1-L,0)) with a > 0, b >= 0 give nonnegative coefficients and Q_A << Q_B for \
     A <= B <= A + L; RealRooted route variableDelay \
     (RealRooted/Interlacing/VariableDelay.lean)";

/// The Lean declarations of the scaled form (any `a > 0`).
const SCALED_DECLARATIONS: [&str; 2] = [
    "RealRooted.strictInterl_of_scaled_fixedDelay_recurrence",
    "RealRooted.hasNonnegCoeffs_of_scaled_fixedDelay_recurrence",
];

/// The Lean declarations of the unscaled form (`a = 1`).
const UNSCALED_DECLARATIONS: [&str; 2] = [
    "RealRooted.strictInterl_of_fixedDelay_recurrence",
    "RealRooted.hasNonnegCoeffs_of_scaled_fixedDelay_recurrence",
];

/// One seed row written as a step of the fixed-delay window:
/// `P_index = previous_scalar P_(index-1) + delayed_scalar x P_delayed_index`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VariableDelaySeedStep {
    /// The absolute index of the seed row.
    pub index: usize,
    /// `a_m > 0`, an exact rational as a string.
    pub previous_scalar: String,
    /// `b_m >= 0`, an exact rational as a string.
    pub delayed_scalar: String,
    /// The absolute index `n_0 + max(m - L, 0)` of the delayed row.
    pub delayed_index: usize,
}

/// The exact data of a [`check_variable_delay_theorem`] report: the delay
/// and the seed rows as steps of the window (module documentation of the
/// `variable_delay` submodule, summarized at
/// [`check_variable_delay_theorem`]).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct VariableDelayCertificate {
    /// The delay `L >= 2` of the `x`-term.
    pub delay: usize,
    /// The constant first row `P_(n_0)`, an exact rational as a string.
    pub first_row: String,
    /// Whether `a(n) = 1` on the rule range and every seed step has
    /// `a_m = 1` (Lean's unscaled form).
    pub unscaled: bool,
    /// The seed rows `P_(n_0+1), ..., P_(N-1)` as steps of the window (empty
    /// for a single seed row); a row that is not a step ends the list.
    pub seed_steps: Vec<VariableDelaySeedStep>,
}

impl VariableDelayCertificate {
    /// Replay the recorded seed steps on `seeds` (`seeds[i]` is the row
    /// `first_index + i`): the first row is the recorded nonzero constant,
    /// every step has `a_m > 0`, `b_m >= 0`, the recorded delayed index and
    /// holds exactly, and every seed row after the first has a step.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] on any mismatch.
    pub fn verify_seed_steps(
        &self,
        first_index: usize,
        seeds: &[Vec<Q>],
    ) -> Result<(), RecurrenceError> {
        let fail = |what: String| Err(RecurrenceError::VerificationFailed(what));
        let first = vec_poly::trim(seeds.first().map_or(&[][..], Vec::as_slice));
        let constant = parse_rational_coeff(&self.first_row)?;
        if constant.is_zero() || first != [constant] {
            return fail("the first row is not the recorded nonzero constant".into());
        }
        if self.delay < 2 || self.seed_steps.len() + 1 != seeds.len() {
            return fail("the certificate does not cover every seed row".into());
        }
        for (m, step) in (1usize..).zip(&self.seed_steps) {
            let delayed = m.saturating_sub(self.delay);
            if step.index != first_index + m || step.delayed_index != first_index + delayed {
                return fail(format!("seed step {m} has the wrong indices"));
            }
            let a = parse_rational_coeff(&step.previous_scalar)?;
            let b = parse_rational_coeff(&step.delayed_scalar)?;
            if !a.is_positive() || b.is_negative() {
                return fail(format!("seed step {m} has a scalar of the wrong sign"));
            }
            let right = vec_poly::add(
                &vec_poly::scale(&seeds[m - 1], &a),
                &vec_poly::scale(&times_x(&seeds[delayed]), &b),
            );
            if vec_poly::trim(&right) != vec_poly::trim(&seeds[m]) {
                return fail(format!("seed step {m} does not hold"));
            }
        }
        Ok(())
    }
}

/// Check the fixed-delay window theorem (RealRooted route `variableDelay`)
/// for the family given by `rule`, `first_index` and the seed rows `seeds`
/// (the rule applies from `N = first_index + seeds.len()` on).
///
/// The rule must be `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)` with `L >= 2`,
/// `f`, `A` and `B` free of `x` (two terms, no derivative, no `(-1)^n`
/// factor, no inhomogeneous term).  The hypotheses, each decided exactly:
///
/// - `f(n) != 0` for `n >= N` (sign checker);
/// - `a(n) = A(n) / f(n) > 0` and `b(n) = B(n) / f(n) >= 0` for every
///   `n >= N` (sign checker);
/// - `seeds.len() >= L`;
/// - `P_(first_index)` is a nonzero constant;
/// - every later seed row is `a_m P_(n-1) + b_m x P_(first_index + max(m-L,
///   0))` with rational `a_m > 0`, `b_m >= 0` (solved exactly, recorded in
///   [`RealRootednessTheoremReport::variable_delay`]);
/// - a requested root region must contain `(-inf, 0)`.
///
/// When all hold, every `P_n`, `n >= first_index`, is nonzero and
/// real-rooted with every root in `(-inf, 0)` and coefficients of the sign
/// of `P_(first_index)`, and `P_A << P_B` (weakly, in the sense of
/// [`RowInterlacing::Weak`]) whenever `first_index <= A <= B <= A + L`; in
/// particular `P_(n-1)` interlaces `P_n` for every `n > first_index`.  The
/// induction start is moved past the seeds as for the other theorems
/// ([`RealRootednessTheoremOptions::maximum_start_delay`]); the moved seed
/// rows are then window steps too.
///
/// # Errors
///
/// Errors of the record itself (seed normalization, row generation for a
/// moved start).
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_variable_delay_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // P_n = (n + 2) P_(n-1) + (n + 1)(n + 2)/2 x P_(n-3), P_0 = 1, P_1 = 3, P_2 = 12
/// // (OEIS A201013).
/// let half = BigRational::new(BigInt::from(1), BigInt::from(2));
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![
///             vec![q(2)],
///             vec![q(1)],
///         ])),
///         RecurrenceTerm::new(3, 0, BivarPoly::from_coefficients_normalized(vec![
///             vec![q(0), q(1)],
///             vec![q(0), q(3) * &half],
///             vec![q(0), half.clone()],
///         ])),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = vec![vec![q(1)], vec![q(3)], vec![q(12)]];
/// let report =
///     check_variable_delay_theorem(&rule, 0, &seeds, &RealRootednessTheoremOptions::default())
///         .unwrap();
/// assert!(report.is_proved());
/// report.verify().unwrap();
/// ```
pub fn check_variable_delay_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    check_real_rootedness_theorem(
        RealRootednessTheorem::VariableDelay,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// `x p`.
fn times_x(p: &[Q]) -> Vec<Q> {
    let mut shifted = vec![Q::zero()];
    shifted.extend_from_slice(p);
    shifted
}

/// `seed = a previous + b x delayed` with `a > 0`, `b >= 0`, when it holds
/// (the constant term of `previous` is nonzero).
fn window_step(seed: &[Q], previous: &[Q], delayed: &[Q]) -> Option<(Q, Q)> {
    let previous_constant = previous.first().filter(|c| !c.is_zero())?;
    let delayed_constant = delayed.first().filter(|c| !c.is_zero())?;
    let a = seed.first().cloned().unwrap_or_else(Q::zero) / previous_constant;
    let remainder = vec_poly::sub(seed, &vec_poly::scale(previous, &a));
    let b = remainder.get(1).cloned().unwrap_or_else(Q::zero) / delayed_constant;
    let expected = vec_poly::scale(&times_x(delayed), &b);
    (vec_poly::trim(&remainder) == vec_poly::trim(&expected) && a.is_positive() && !b.is_negative())
        .then_some((a, b))
}

/// The registry entry of the fixed-delay theorem ([`super::registry::ROUTES`]).
pub(super) const ROUTE: TheoremRoute = TheoremRoute {
    theorem: RealRootednessTheorem::VariableDelay,
    check_at_supplied_seeds: variable_delay_report,
    start_policy: StartPolicy::MovesStart,
    is_preferred: any_proof,
    certificate_name: Some("fixed-delay"),
    discovery: DiscoveryPolicy::WhenApplicable {
        tolerated_when_missing: false,
    },
    in_residue_class_menu: true,
    in_index_section_menu: true,
};

/// Replay the certificate of a recomputed fixed-delay report: the seed steps
/// of a proved report on the checked (possibly moved) seed rows.
pub(super) fn verify_certificate_in_report(
    certificate: &VariableDelayCertificate,
    report: &RealRootednessTheoremReport,
) -> Result<(), RecurrenceError> {
    if !report.is_proved() {
        return Ok(());
    }
    let (rule, first_index, seeds) = report.recurrence.to_recurrence_parts()?;
    let seeds = if report.start_delay == 0 {
        seeds
    } else {
        let generated = rule
            .generate_rows_rational(&seeds, first_index, seeds.len() + report.start_delay)
            .map_err(|error| RecurrenceError::VerificationFailed(error.to_string()))?;
        normalized_seed_rows(&rule, first_index, &generated)?
    };
    certificate.verify_seed_steps(first_index, &seeds)
}

fn variable_delay_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::VariableDelay,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let form = "the rule is not f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L) with L >= 2";
    let keys: Vec<(usize, usize)> = shape.terms.keys().copied().collect();
    let delay = match keys.as_slice() {
        [(1, 0), (delay, 0)] if *delay >= 2 => *delay,
        _ => return Ok(report.not_applicable(form)),
    };
    let previous = &shape.terms[&(1, 0)];
    let delayed = &shape.terms[&(delay, 0)];
    if x_degree(previous) != Some(0)
        || x_degree(delayed) != Some(1)
        || !x_column(delayed, 0).is_empty()
    {
        return Ok(report.not_applicable(form));
    }
    if seeds.len() < delay {
        return Ok(report.not_applicable(format!(
            "the rule has delay {delay} but only {} seed row(s)",
            seeds.len()
        )));
    }
    let f = &shape.denominator;
    let over_f = |p: Vec<Q>| IndexFunction::new(p, f.clone());
    let a = over_f(x_column(previous, 0));
    let b = over_f(x_column(delayed, 1));
    let n_first = first_index + seeds.len();
    report.normal_form = format!(
        "P_n = a(n) P_(n-1) + b(n) x P_(n-{delay}) for n >= {n_first}, with a(n) = {} and \
         b(n) = {}",
        a.text(),
        b.text()
    );
    report.lean_route = Some("variableDelay".into());

    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    hypotheses.push(sign_hypothesis(
        "previous_coefficient_positive",
        &format!("a(n) > 0 for every n >= {n_first}"),
        "a(n)",
        &a,
        n_first,
        IndexSignCondition::Positive,
    ));
    hypotheses.push(sign_hypothesis(
        "delayed_coefficient_nonnegative",
        &format!("b(n) >= 0 for every n >= {n_first}"),
        "b(n)",
        &b,
        n_first,
        IndexSignCondition::Nonnegative,
    ));
    let first_row = vec_poly::trim(&seeds[0]);
    let constant = match first_row.as_slice() {
        [c] if !c.is_zero() => Some(c.clone()),
        _ => None,
    };
    hypotheses.push(exact_hypothesis(
        "first_row_constant",
        &format!("P_{first_index} is a nonzero constant"),
        match &constant {
            Some(c) => Ok(format!("P_{first_index} = {}", format_rational_coeff(c))),
            None => Err((
                Some(first_index),
                format!("P_{first_index} is zero or not constant (exact)"),
            )),
        },
    ));
    let mut seed_steps = Vec::new();
    let mut unit_seed_scalars = true;
    let mut failing_step = None;
    if constant.is_some() {
        for m in 1..seeds.len() {
            let delayed_offset = m.saturating_sub(delay);
            match window_step(&seeds[m], &seeds[m - 1], &seeds[delayed_offset]) {
                Some((a_m, b_m)) => {
                    unit_seed_scalars &= a_m.is_one();
                    seed_steps.push(VariableDelaySeedStep {
                        index: first_index + m,
                        previous_scalar: format_rational_coeff(&a_m),
                        delayed_scalar: format_rational_coeff(&b_m),
                        delayed_index: first_index + delayed_offset,
                    });
                }
                None => {
                    failing_step = Some(first_index + m);
                    break;
                }
            }
        }
    }
    if seeds.len() > 1 {
        let last = n_first - 1;
        let statement = format!(
            "every seed row P_m, {} <= m <= {last}, is a_m P_(m-1) + b_m x P_(max(m-{delay}, \
             {first_index})) with rational a_m > 0, b_m >= 0",
            first_index + 1
        );
        hypotheses.push(exact_hypothesis(
            "seed_window_steps",
            &statement,
            match (constant.is_some(), failing_step) {
                (true, None) => Ok(format!(
                    "solved exactly for {} seed row(s); the scalars are recorded in the \
                     certificate",
                    seed_steps.len()
                )),
                (true, Some(index)) => Err((
                    Some(index),
                    format!("P_{index} is not such a step of the earlier rows (exact)"),
                )),
                (false, _) => Err((
                    Some(first_index),
                    "the window needs a nonzero constant first row".into(),
                )),
            },
        ));
    }
    if let Some(region) = &options.root_region {
        let contains_negative_axis = region.lower.is_none()
            && region
                .upper
                .as_ref()
                .is_none_or(|upper| !upper.is_negative());
        let statement = format!(
            "the requested region {} contains (-inf, 0), the region this theorem certifies",
            open_root_interval(region)
        );
        hypotheses.push(if contains_negative_axis {
            exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
        } else {
            unresolved_hypothesis(
                "requested_root_region",
                &statement,
                "this theorem certifies the region (-inf, 0) only".into(),
            )
        });
    }
    report.hypotheses = hypotheses;
    let unscaled = a == IndexFunction::constant(Q::one()) && unit_seed_scalars;
    report.lean_declarations = if unscaled {
        UNSCALED_DECLARATIONS
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        SCALED_DECLARATIONS.iter().map(|s| s.to_string()).collect()
    };
    report.certificate = constant.map(|c| {
        RouteCertificate::VariableDelay(VariableDelayCertificate {
            delay,
            first_row: format_rational_coeff(&c),
            unscaled,
            seed_steps,
        })
    });
    let negative_axis = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: false,
    };
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        zero_residue_classes: None,
        interlacing_from: Some(first_index + 1),
        interlacing: Some(RowInterlacing::Weak),
        simple_roots: false,
        root_region: Some(open_root_interval(&negative_axis)),
        statement: format!(
            "every P_n, n >= {first_index}, is nonzero and real-rooted with every root in \
             (-inf, 0) (coefficients of the sign of P_{first_index}), and P_A << P_B whenever \
             {first_index} <= A <= B <= A + {delay}; in particular P_(n-1) interlaces P_n for \
             every n > {first_index} (fixed-delay window, delay {delay})"
        ),
    }))
}
