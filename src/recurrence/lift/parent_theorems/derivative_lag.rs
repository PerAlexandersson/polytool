//! Second-order derivative-lag recurrences
//! `P_n = A_n(x) P_(n-1) + V_n(x) P'_(n-1) + W_n(x) P_(n-2)`: strict
//! interlacing of consecutive rows and every root in `(-inf, 0)`, by the
//! Liu--Wang sign argument with the derivative and the lag row as
//! interlacers, optionally after dividing every row by a power of `x`.
//!
//! # The theorem
//!
//! For the rule `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1) +
//! c_2(n, x) P_(n-2)` (applied for `n >= N`) write `A_n = c_0 / f`,
//! `V_n = c_1 / f`, `W_n = c_2 / f`, and require `deg_x c_0 <= 1`,
//! `deg_x c_1 <= 2`, `deg_x c_2 <= 2`.  Let `d = deg P_(N-1)` and, for
//! `n >= N`,
//!
//! ```text
//! lambda_n = [x^1] A_n + (d + n - N) [x^2] V_n,      mu_n = [x^2] W_n.
//! ```
//!
//! Assume, each for every `n >= N`:
//!
//! 1. `f(n) != 0`;
//! 2. `P_(N-2)` and `P_(N-1)` are nonzero and real-rooted with every root
//!    `< 0`, `deg P_(N-1) = deg P_(N-2) + 1`, their leading coefficients have
//!    one sign, and `P_(N-2)` strictly interlaces `P_(N-1)` (simple roots, no
//!    common root; the exact Bézout test);
//! 3. degree law: `lambda_n > 0` when `mu = 0`, otherwise `lambda_n >= 0`,
//!    `mu_n >= 0` and `lambda_n + mu_n > 0`;
//! 4. `V_n(0) = 0` and `W_n(0) = 0`, and `A_n(0) > 0`;
//! 5. `V_n(x) <= 0` and `W_n(x) <= 0` for every `x < 0`;
//! 6. for every `x < 0`, `V_n(x) < 0` or `W_n(x) < 0` (decided as
//!    `V_n < 0` on `(-inf, 0)` for every `n`, or `W_n < 0` there for every
//!    `n`).
//!
//! Then for every `n >= N - 1`, `P_n` has `deg P_n = d + n - N + 1`, simple
//! real roots, every root `< 0`, and `P_(n-1)` strictly interlaces `P_n`.
//!
//! # Proof
//!
//! Induction on `n >= N`; replacing every row by `-P_n` (the rule is linear)
//! we may assume positive leading coefficients.  Suppose `P_(n-2)`,
//! `P_(n-1)` have positive leading coefficients, degrees `D - 1`, `D`
//! (`D >= 1`), every root `< 0`, and `P_(n-2)` strictly interlaces
//! `P_(n-1)`; let `r_1 > ... > r_D` be the (simple) roots of `P_(n-1)`.
//!
//! - *Degree.*  The coefficient of `x^(D+1)` on the right is
//!   `lambda_n lead P_(n-1) + mu_n lead P_(n-2)` (from `[x^1] A_n`,
//!   `[x^2] V_n` times `lead P'_(n-1) = D lead P_(n-1)`, and `[x^2] W_n`), and
//!   no higher power occurs; by 3 it is positive.  So `deg P_n = D + 1` with
//!   a positive leading coefficient.
//! - *Value at `0`.*  By 4, `P_n(0) = A_n(0) P_(n-1)(0) > 0`
//!   (`P_(n-1)(0) > 0` because its roots are negative and its leading
//!   coefficient positive).
//! - *Signs at the roots of `P_(n-1)`.*  `P'_(n-1)(r_i)` and `P_(n-2)(r_i)`
//!   both have the sign `(-1)^(i-1)` (simple roots; strict interlacing puts
//!   exactly `i - 1` roots of `P_(n-2)` right of `r_i`).  So
//!   `P_n(r_i) = V_n(r_i) P'_(n-1)(r_i) + W_n(r_i) P_(n-2)(r_i)` is a sum of
//!   two terms of sign `(-1)^i` or zero (by 5, as `r_i < 0`), not both zero
//!   (by 6): `P_n(r_i)` has the sign `(-1)^i`.
//! - *Counting.*  `P_n` has a root in each `(r_(i+1), r_i)`; one in
//!   `(r_1, 0)` (`P_n(r_1) < 0 < P_n(0)`); and one in `(-inf, r_D)` (the sign
//!   `(-1)^D` at `r_D` against `(-1)^(D+1)` at `-inf`).  These are
//!   `D + 1 = deg P_n` distinct roots, all `< 0`, strictly interlaced by the
//!   roots of `P_(n-1)` (none of which is a root of `P_n`).
//!
//! # Nonnegative coefficients
//!
//! The Lean statements (below) and Liu--Wang's Corollary 2.4 assume that
//! every row has nonnegative coefficients and (Lean) a positive value at
//! `0`.  For a real-rooted row with a positive leading coefficient, "every
//! root `< 0`" is equivalent to "nonnegative coefficients and a positive
//! value at `0`" (`prod (x - r_i)` with `r_i < 0` has positive coefficients;
//! conversely nonnegative coefficients exclude positive roots and
//! `P(0) > 0` excludes `0`).  The induction above proves this root location
//! together with real-rootedness instead of assuming the coefficients, and
//! hypothesis 4 is what Lean's `derivLag_eval_zero_pos` needs for
//! `P_n(0) > 0` (there `U(0) >= 0`, `V(0) >= 0`, `W(0) >= 0` with `V <= 0`,
//! `W <= 0` left of `0` force `V(0) = W(0) = 0`).
//!
//! # Gauge by a power of `x`
//!
//! When `x^e` (`e >= 1`) divides both `P_(N-2)` and `P_(N-1)` and `V_n(0) = 0`
//! (hypothesis 4), `x^e` divides every `P_n`, `n >= N - 2`: if it divides
//! `P_(n-1)` and `P_(n-2)`, then `x^(e-1)` divides `P'_(n-1)` and `x` divides
//! `V_n`, so `x^e` divides the right-hand side.  With `P_n = x^e Q_n`,
//! `P'_(n-1) = e x^(e-1) Q_(n-1) + x^e Q'_(n-1)`, so
//!
//! ```text
//! Q_n = (A_n + e V_n / x) Q_(n-1) + V_n Q'_(n-1) + W_n Q_(n-2),
//! ```
//!
//! and the theorem is applied to this rule with the seeds `Q_(N-2)`,
//! `Q_(N-1)` (`e` is the largest such power, so hypothesis 2 then asks for
//! `Q(0) != 0`).  Every `P_n = x^e Q_n`, `n >= N - 2`, is then real-rooted
//! with every root in `(-inf, 0]` (`0` with multiplicity exactly `e`), and
//! `P_(n-1)` interlaces `P_n` in the weak sense of Liu and Wang (the common
//! root `0` aside, the roots strictly alternate).  Example: the derangement
//! excedance rows (OEIS A271697) `P_n = (n-2) x P_(n-1) + x (1-x) P'_(n-1) +
//! (n-2) x P_(n-2)` are `x Q_n` from `P_3 = x` on, and the `Q_n` are the rows
//! of A046739.
//!
//! # Earlier rows
//!
//! Seed rows before `P_(N-2)` are checked exactly: each must be zero or
//! real-rooted with every root in the certified region (`(-inf, 0)`, or
//! `(-inf, 0]` with a gauge).  Zero rows are named in the conclusion, whose
//! [`RealRootednessTheoremConclusion::first_index`] is then the index after
//! the last zero row.
//!
//! # Lean correspondence
//!
//! The hypotheses mirror the row routes of the public Lean library
//! RealRooted (<https://github.com/PerAlexandersson/RealRooted>), recorded in
//! [`RealRootednessTheoremReport::lean_route`] and
//! [`RealRootednessTheoremReport::lean_declarations`]:
//!
//! - `V_n < 0` on `(-inf, 0)`: route `halfGrowth`,
//!   `RealRooted.derivLag_strictInterl_of_nonnegCoeffs` (the degree law
//!   `D0 + n`), whose induction `derivLag_strictInterl_and_noCommonRoot` is
//!   the argument above, and `derivLag_eval_zero_pos`;
//! - `V = 0` (three-term rules such as the Delannoy rows): route `sign`,
//!   `RealRooted.threeTerm_interlaces_of_nonnegCoeffs`, with the strict form
//!   `RealRooted.LiuWang.strictInterl_nonpos_lag_sequence`;
//! - `V_n <= 0` and `W_n < 0`: the named tactic `LiuWang`, step theorem
//!   `RealRooted.LiuWang.strictInterl_derivative_lag_of_eval_nonpos`.
//!
//! The Lean theorems take nonnegative coefficients, positive leading
//! coefficients, `P_n(0) > 0` and the degree law as hypotheses for every
//! row; here they are conclusions of the induction (see "Nonnegative
//! coefficients"), so a Lean replay must establish them separately.  The
//! gauge has no Lean route of its own (the Lean proofs divide by the power
//! of `x` by hand).  The half-growth degree law `D0 + (n + e)/2` of
//! `derivLag_strictInterl_of_nonnegCoeffs_half` (degree plateaus) is not
//! covered: hypothesis 3 asks for growth one at every step.
//!
//! # Deciding the hypotheses
//!
//! Hypotheses 1, 3 and `A_n(0) > 0` are sign conditions on rational
//! functions of `n` ([`crate::decide_index_rational_function_sign`]);
//! `V_n(0) = W_n(0) = 0` is checked on the coefficient polynomials; 2 and the
//! earlier rows are exact row computations.  A condition "`c(n, x) <= 0`
//! (or `< 0`) for every `x < 0`" is decided for a separable
//! `c = gamma(n) w(x)` with `w` without roots in `(-inf, 0)` (so of one sign
//! there) by the sign of `gamma / f` (exact; a failure is a disproof at that
//! `n`), and otherwise by the sufficient test after `x = -t`, `t > 0`: every
//! `(-1)^k [x^k] c / f <= 0`, and for `< 0` one of them `< 0` for every `n`
//! (a failure is not a disproof).

use super::*;
use crate::real_rootedness::isolate_real_roots_rational_coeffs;
use std::cmp::Ordering;

const REFERENCE: &str = "L. L. Liu and Y. Wang, A unified approach to polynomial sequences with \
     only real zeros, Adv. in Appl. Math. 38 (2007) 542-560, Theorem 2.1 with the interlacers \
     f' and the lag row (strict form, as in the derivative-lag theorems of the RealRooted \
     library, https://github.com/PerAlexandersson/RealRooted), with nonnegativity of the \
     coefficients replaced by an exact induction keeping every root in (-inf, 0)";

/// Lean route of the derivative-lag rows with `V_n < 0` (RealRooted).
const LEAN_ROUTE_HALF_GROWTH: &str = "halfGrowth";
/// Lean route of the three-term rows (`V = 0`).
const LEAN_ROUTE_SIGN: &str = "sign";
/// Lean named tactic of the derivative-lag rows with `W_n < 0`.
const LEAN_ROUTE_LIU_WANG: &str = "LiuWang";

/// Check the Liu--Wang derivative-lag theorem (module documentation of the
/// `derivative_lag` submodule, summarized here) for the family given by
/// `rule`, `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on; at least two seed rows).
///
/// The rule must be `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1) +
/// c_2(n, x) P_(n-2)` with `f` free of `x`, `deg_x c_0 <= 1`,
/// `deg_x c_1 <= 2`, `deg_x c_2 <= 2`, and a derivative or a lag term.  With
/// `A_n = c_0 / f`, `V_n = c_1 / f`, `W_n = c_2 / f`, the hypotheses, each
/// decided exactly for every `n >= N`, are:
///
/// - `f(n) != 0` (sign checker);
/// - `P_(N-2)`, `P_(N-1)` nonzero, real-rooted with every root `< 0`, of
///   degrees `d - 1`, `d` with leading coefficients of one sign, `P_(N-2)`
///   strictly interlacing `P_(N-1)` (exact); every earlier seed row zero or
///   real-rooted with every root `< 0` (exact);
/// - the degree law `deg P_n = deg P_(n-1) + 1` with leading coefficients of
///   one sign: `lambda_n = [x^1] A_n + (d + n - N) [x^2] V_n > 0` when
///   `deg_x c_2 <= 1`, otherwise `lambda_n >= 0`, `mu_n = [x^2] W_n >= 0`
///   and `lambda_n + mu_n > 0` (sign checker);
/// - `V_n(0) = W_n(0) = 0` (coefficient polynomials) and `A_n(0) > 0` (sign
///   checker);
/// - `V_n(x) <= 0` and `W_n(x) <= 0` for `x < 0`, and `V_n < 0` there for
///   every `n` or `W_n < 0` there for every `n` (separable form and sign
///   checker, or the sufficient monomial test);
/// - a requested root region must contain the certified region.
///
/// When all hold, every `P_n`, `n >= first_index`, is zero (only seed rows,
/// listed in the conclusion) or real-rooted with every root in `(-inf, 0)`
/// (so its coefficients have one sign), and for `n >= N - 1` the rows have
/// simple roots and `P_(n-1)` strictly interlaces `P_n`.
///
/// When `x^e` divides both `P_(N-2)` and `P_(N-1)` (for example
/// `P_n(0) = 0`), the check runs on `Q_n = P_n / x^e` (the rule becomes
/// `Q_n = (A_n + e V_n / x) Q_(n-1) + V_n Q'_(n-1) + W_n Q_(n-2)`, which
/// needs `V_n(0) = 0`), and concludes every root in `(-inf, 0]` with weak
/// (Liu--Wang) interlacing of consecutive rows.  The induction start is moved
/// as in [`check_liu_wang_theorem`] when the check fails at the supplied
/// seeds ([`RealRootednessTheoremOptions::maximum_start_delay`]).  The report
/// records the mirrored RealRooted route and declarations
/// ([`RealRootednessTheoremReport::lean_route`]).
///
/// # Errors
///
/// Errors of the seed record itself (as for [`check_liu_wang_theorem`]); a
/// rule that cannot generate the rows of a moved induction start ends the
/// retries instead of returning an error.
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_liu_wang_derivative_lag_theorem, BigRational, BivarPoly, Recurrence,
///     RecurrenceTerm, RealRootednessTheoremOptions, RowInterlacing,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // P_n = (1 + n x) P_(n-1) + x (1 - x) P'_(n-1) + x P_(n-2), P_0 = 1, P_1 = 1 + x
/// // (OEIS A144438).
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(
///             1,
///             0,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(0)], vec![q(0), q(1)]]),
///         ),
///         RecurrenceTerm::new(
///             1,
///             1,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1), q(-1)]]),
///         ),
///         RecurrenceTerm::new(
///             2,
///             0,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1)]]),
///         ),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seeds = [vec![q(1)], vec![q(1), q(1)]];
/// let report = check_liu_wang_derivative_lag_theorem(
///     &rule,
///     0,
///     &seeds,
///     &RealRootednessTheoremOptions::default(),
/// )
/// .unwrap();
/// assert!(report.is_proved());
/// assert_eq!(report.lean_route.as_deref(), Some("halfGrowth"));
/// let conclusion = report.conclusion.as_ref().unwrap();
/// assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
/// assert_eq!(conclusion.root_region.as_ref().unwrap().to_string(), "(-inf, 0)");
/// report.verify().unwrap();
/// ```
pub fn check_liu_wang_derivative_lag_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    super::check_with_first_proving_start(
        RealRootednessTheorem::LiuWangDerivativeLag,
        rule,
        first_index,
        seeds,
        options,
    )
}

/// The check with the induction at the given (normalized) seeds.
pub(super) fn derivative_lag_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::LiuWangDerivativeLag,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.len() < 2 {
        return Ok(report.not_applicable("the theorem needs two seed rows"));
    }
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    let allowed = [(1, 0), (1, 1), (2, 0)];
    if shape.terms.keys().any(|key| !allowed.contains(key))
        || !(shape.terms.contains_key(&(1, 1)) || shape.terms.contains_key(&(2, 0)))
    {
        return Ok(report.not_applicable(
            "the rule is not f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1) + c_2(n, x) \
             P_(n-2) with a derivative or a lag term",
        ));
    }
    let empty = Bivariate::new();
    let c0 = shape.terms.get(&(1, 0)).unwrap_or(&empty);
    let c1 = shape.terms.get(&(1, 1)).unwrap_or(&empty);
    let c2 = shape.terms.get(&(2, 0)).unwrap_or(&empty);
    if x_degree(c0).is_some_and(|d| d > 1)
        || x_degree(c1).is_some_and(|d| d > 2)
        || x_degree(c2).is_some_and(|d| d > 2)
    {
        return Ok(report.not_applicable(
            "deg_x c_0 > 1, deg_x c_1 > 2 or deg_x c_2 > 2: the degree may grow by more than one",
        ));
    }
    let f = &shape.denominator;
    let n_first = first_index + seeds.len();
    let previous_index = n_first - 2;
    let current_index = n_first - 1;
    let previous = vec_poly::trim(&seeds[seeds.len() - 2]);
    let current = vec_poly::trim(&seeds[seeds.len() - 1]);

    // The gauge x^e, valid when x divides c_1 (V_n(0) = 0).
    let derivative_vanishes_at_zero = x_column(c1, 0).is_empty();
    let gauge = if derivative_vanishes_at_zero {
        match (order_at_zero(&previous), order_at_zero(&current)) {
            (Some(a), Some(b)) => a.min(b),
            _ => 0,
        }
    } else {
        0
    };
    // c_0 of the rule for Q_n = P_n / x^gauge (c_1 and c_2 are unchanged).
    let gauged_current = if gauge == 0 {
        c0.clone()
    } else {
        add_bivariate(c0, &scale_bivariate(&divide_by_x(c1), gauge))
    };
    let previous_quotient = previous[gauge..].to_vec();
    let current_quotient = current[gauge..].to_vec();
    let row = |index: usize| {
        if gauge == 0 {
            format!("P_{index}")
        } else {
            format!("P_{index} / x^{gauge}")
        }
    };
    let over_f = |p: Vec<Q>| IndexFunction::new(p, f.clone());
    let multipliers = format!(
        "A_n(x) = {}, V_n(x) = {}, W_n(x) = {}",
        x_polynomial_text(&gauged_current, f),
        x_polynomial_text(c1, f),
        x_polynomial_text(c2, f)
    );
    report.normal_form = if gauge == 0 {
        format!(
            "P_n = A_n(x) P_(n-1) + V_n(x) P'_(n-1) + W_n(x) P_(n-2) for n >= {n_first}, with \
             {multipliers}"
        )
    } else {
        format!(
            "P_n = x^{gauge} Q_n for n >= {previous_index} (x^{gauge} divides P_{previous_index} \
             and P_{current_index}, and V_n(0) = 0), with Q_n = A_n(x) Q_(n-1) + V_n(x) \
             Q'_(n-1) + W_n(x) Q_(n-2) for n >= {n_first}, {multipliers} (A_n = a_n + \
             {gauge} V_n / x for the rule's a_n(x) = {})",
            x_polynomial_text(c0, f)
        )
    };

    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    // Hypothesis 2: the two seed rows of the induction.
    let previous_ok = has_real_roots_below(&previous_quotient, &Q::zero());
    let current_ok = has_real_roots_below(&current_quotient, &Q::zero());
    hypotheses.push(exact_hypothesis(
        "seed_roots_negative",
        &format!(
            "{} and {} are nonzero and real-rooted with every root < 0",
            row(previous_index),
            row(current_index)
        ),
        match (previous_ok, current_ok) {
            (true, true) => Ok("checked exactly".into()),
            (false, _) => Err((
                Some(previous_index),
                format!(
                    "{} is zero, not real-rooted, or has a root >= 0 (exact)",
                    row(previous_index)
                ),
            )),
            (true, false) => Err((
                Some(current_index),
                format!(
                    "{} is zero, not real-rooted, or has a root >= 0 (exact)",
                    row(current_index)
                ),
            )),
        },
    ));
    let previous_degree = vec_poly::degree(&previous_quotient);
    let current_degree = vec_poly::degree(&current_quotient);
    let degree_step = match (previous_degree, current_degree) {
        (Some(a), Some(b)) => {
            let same_sign = previous_quotient[a].is_positive() == current_quotient[b].is_positive();
            if b == a + 1 && same_sign {
                Ok(format!(
                    "degrees {a} and {b}, leading coefficients of one sign"
                ))
            } else {
                Err((
                    Some(current_index),
                    format!(
                        "degrees {a} and {b}{} (exact)",
                        if same_sign {
                            ""
                        } else {
                            ", leading coefficients of opposite signs"
                        }
                    ),
                ))
            }
        }
        _ => Err((
            Some(if previous_degree.is_none() {
                previous_index
            } else {
                current_index
            }),
            "a seed row of the induction is zero".into(),
        )),
    };
    hypotheses.push(exact_hypothesis(
        "seed_degree_step",
        &format!(
            "deg {} = deg {} + 1, with leading coefficients of one sign",
            row(current_index),
            row(previous_index)
        ),
        degree_step,
    ));
    let interlaces = crate::check_interlacing_bigint_coeffs(
        &integer_coefficients(&previous_quotient),
        &integer_coefficients(&current_quotient),
    ) == Some(true)
        && !vec_poly::is_zero(&previous_quotient);
    hypotheses.push(exact_hypothesis(
        "seed_strict_interlacing",
        &format!(
            "{} strictly interlaces {} (simple real roots, no common root)",
            row(previous_index),
            row(current_index)
        ),
        if interlaces {
            Ok("Bezout test, exact".into())
        } else {
            Err((
                Some(current_index),
                "the Bezout test fails or the degrees do not differ by one (exact)".into(),
            ))
        },
    ));
    let (earlier, zero_rows) =
        earlier_rows_hypothesis(&seeds[..seeds.len() - 2], first_index, gauge > 0);
    hypotheses.extend(earlier);
    // Hypothesis 3: the degree law.
    let d = current_degree.unwrap_or(0);
    let previous_row_degree = vec![
        Q::from_integer(BigInt::from(d as i64 - n_first as i64)),
        Q::one(),
    ];
    let lambda = over_f(vec_poly::add(
        &vec_poly::trim(&x_column(&gauged_current, 1)),
        &vec_poly::mul(&previous_row_degree, &vec_poly::trim(&x_column(c1, 2))),
    ));
    let lag_leading = x_column(c2, 2);
    let degree_statement = "(deg P_n = deg P_(n-1) + 1 with leading coefficients of one sign)";
    if lag_leading.is_empty() {
        hypotheses.push(sign_hypothesis(
            "degree_step_positive",
            &format!(
                "lambda_n = [x^1] A_n + (deg P_(n-1)) [x^2] V_n > 0 for every n >= {n_first} \
                 {degree_statement}"
            ),
            "lambda_n",
            &lambda,
            n_first,
            IndexSignCondition::Positive,
        ));
    } else {
        let mu = over_f(lag_leading);
        hypotheses.push(sign_hypothesis(
            "derivative_degree_step_nonnegative",
            &format!(
                "lambda_n = [x^1] A_n + (deg P_(n-1)) [x^2] V_n >= 0 for every n >= {n_first}"
            ),
            "lambda_n",
            &lambda,
            n_first,
            IndexSignCondition::Nonnegative,
        ));
        hypotheses.push(sign_hypothesis(
            "lag_degree_step_nonnegative",
            &format!("mu_n = [x^2] W_n >= 0 for every n >= {n_first}"),
            "mu_n",
            &mu,
            n_first,
            IndexSignCondition::Nonnegative,
        ));
        hypotheses.push(sign_hypothesis(
            "degree_step_positive",
            &format!("lambda_n + mu_n > 0 for every n >= {n_first} {degree_statement}"),
            "lambda_n + mu_n",
            &lambda.add(&mu),
            n_first,
            IndexSignCondition::Positive,
        ));
    }
    // Hypothesis 4.
    hypotheses.push(vanishing_at_zero_hypothesis(c1, c2, n_first));
    hypotheses.push(sign_hypothesis(
        "constant_term_positive",
        &format!("A_n(0) > 0 for every n >= {n_first}"),
        "A_n(0)",
        &over_f(x_column(&gauged_current, 0)),
        n_first,
        IndexSignCondition::Positive,
    ));
    // Hypotheses 5 and 6.
    hypotheses.push(negative_axis_sign_hypothesis(
        "derivative_coefficient_nonpositive",
        &format!("V_n(x) <= 0 for every x < 0 and every n >= {n_first}"),
        "V_n",
        c1,
        f,
        n_first,
        false,
    ));
    hypotheses.push(negative_axis_sign_hypothesis(
        "lag_coefficient_nonpositive",
        &format!("W_n(x) <= 0 for every x < 0 and every n >= {n_first}"),
        "W_n",
        c2,
        f,
        n_first,
        false,
    ));
    let (strict, lean_route) = strict_sign_hypothesis(c1, c2, f, n_first);
    hypotheses.push(strict);
    // The requested region.
    let certified = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: gauge > 0,
    };
    if let Some(region) = &options.root_region {
        hypotheses.push(requested_region_hypothesis(region, gauge, current_index));
    }
    report.hypotheses = hypotheses;
    report.lean_route = Some(lean_route.into());
    report.lean_declarations = lean_declarations(lean_route);

    let root_region = RootInterval::from_real(&certified);
    let conclusion_first = zero_rows.last().map_or(first_index, |last| last + 1);
    let zero_text = if zero_rows.is_empty() {
        String::new()
    } else {
        let list: Vec<String> = zero_rows.iter().map(|j| format!("P_{j}")).collect();
        format!(
            "{} {} zero; ",
            list.join(", "),
            if zero_rows.len() == 1 { "is" } else { "are" }
        )
    };
    let statement = if gauge == 0 {
        format!(
            "{zero_text}every {}P_n, n >= {first_index}, is nonzero and real-rooted with every \
             root in {root_region} (so its coefficients have one sign); for n >= \
             {current_index}, P_n has simple roots and P_(n-1) strictly interlaces P_n (no \
             common root), and deg P_n = deg P_(n-1) + 1 for n >= {n_first}",
            if zero_rows.is_empty() { "" } else { "other " }
        )
    } else {
        format!(
            "{zero_text}every {}P_n, n >= {first_index}, is nonzero and real-rooted with every \
             root in {root_region} (so its coefficients have one sign); P_n = x^{gauge} Q_n for \
             n >= {previous_index} with Q_n(0) != 0, and for n >= {current_index} the Q_n have \
             simple roots, Q_(n-1) strictly interlaces Q_n, and P_(n-1) interlaces P_n \
             (Liu--Wang, common root 0)",
            if zero_rows.is_empty() { "" } else { "other " }
        )
    };
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index: conclusion_first,
        interlacing_from: Some(current_index),
        interlacing: Some(if gauge == 0 {
            RowInterlacing::Strict
        } else {
            RowInterlacing::Weak
        }),
        simple_roots: gauge <= 1,
        root_region: Some(root_region),
        statement,
    }))
}

/// The multiplicity of the root `0` of a nonzero `p`.
fn order_at_zero(p: &[Q]) -> Option<usize> {
    p.iter().position(|c| !c.is_zero())
}

/// `c(n, x) / x` for a `c` divisible by `x`.
fn divide_by_x(c: &Bivariate) -> Bivariate {
    c.iter()
        .map(|row| row.iter().skip(1).cloned().collect())
        .collect()
}

fn scale_bivariate(c: &Bivariate, factor: usize) -> Bivariate {
    let factor = Q::from_integer(BigInt::from(factor));
    c.iter()
        .map(|row| row.iter().map(|value| value * &factor).collect())
        .collect()
}

/// `c(n, x) / f(n)` as `(..) + (..) x + (..) x^2`, `0` when zero.
fn x_polynomial_text(c: &Bivariate, f: &[Q]) -> String {
    let terms: Vec<String> = (0..x_width(c))
        .filter_map(|k| {
            let column = x_column(c, k);
            (!column.is_empty()).then(|| {
                let coefficient = IndexFunction::new(column, f.to_vec()).text();
                match k {
                    0 => format!("({coefficient})"),
                    1 => format!("({coefficient}) x"),
                    _ => format!("({coefficient}) x^{k}"),
                }
            })
        })
        .collect();
    if terms.is_empty() {
        "0".into()
    } else {
        terms.join(" + ")
    }
}

/// The seed rows before the induction, each zero or real-rooted with every
/// root `< 0` (`<= 0` when `closed_at_zero`), and the indices of the zero
/// rows.
fn earlier_rows_hypothesis(
    rows: &[Vec<Q>],
    first_index: usize,
    closed_at_zero: bool,
) -> (Option<LiftHypothesis>, Vec<usize>) {
    if rows.is_empty() {
        return (None, Vec::new());
    }
    let last = first_index + rows.len() - 1;
    let region = if closed_at_zero {
        "(-inf, 0]"
    } else {
        "(-inf, 0)"
    };
    let mut zero_rows = Vec::new();
    let mut failing = None;
    for (offset, row) in rows.iter().enumerate() {
        let row = vec_poly::trim(row);
        if vec_poly::is_zero(&row) {
            zero_rows.push(first_index + offset);
            continue;
        }
        let part = if closed_at_zero {
            row[order_at_zero(&row).unwrap_or(0)..].to_vec()
        } else {
            row
        };
        if failing.is_none() && !has_real_roots_below(&part, &Q::zero()) {
            failing = Some(first_index + offset);
        }
    }
    let outcome = match failing {
        None => Ok(format!(
            "checked exactly for P_{first_index}, ..., P_{last}{}",
            if zero_rows.is_empty() {
                String::new()
            } else {
                format!(" ({} zero)", zero_rows.len())
            }
        )),
        Some(index) => Err((
            Some(index),
            format!("P_{index} is not real-rooted with every root in {region} (exact)"),
        )),
    };
    (
        Some(exact_hypothesis(
            "earlier_seeds",
            &format!(
                "the seed rows P_{first_index}, ..., P_{last} are each zero or real-rooted with \
                 every root in {region}"
            ),
            outcome,
        )),
        zero_rows,
    )
}

/// `V_n(0) = 0` and `W_n(0) = 0` for `n >= first`: the `x^0` coefficient
/// polynomials vanish.
fn vanishing_at_zero_hypothesis(c1: &Bivariate, c2: &Bivariate, first: usize) -> LiftHypothesis {
    let name = "multipliers_vanish_at_zero";
    let statement = format!("V_n(0) = 0 and W_n(0) = 0 for every n >= {first}");
    let mut failures = Vec::new();
    let mut witness: Option<usize> = None;
    for (label, c) in [("c_1", c1), ("c_2", c2)] {
        let column = x_column(c, 0);
        if column.is_empty() {
            continue;
        }
        // A nonzero polynomial in n has at most deg roots.
        let at = (first..=first + column.len())
            .find(|&n| !vec_poly::evaluate(&column, &Q::from_integer(BigInt::from(n))).is_zero());
        witness = match (witness, at) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        };
        failures.push(format!(
            "[x^0] {label}(n, x) = {} is not the zero polynomial",
            format_index_polynomial(&column, "n")
        ));
    }
    exact_hypothesis(
        name,
        &statement,
        if failures.is_empty() {
            Ok("[x^0] c_1 and [x^0] c_2 are zero polynomials".into())
        } else {
            Err((witness, failures.join("; ")))
        },
    )
}

/// `c(n, x) / f(n) <= 0` (`< 0` when `strict`) for every `x < 0` and
/// `n >= first` (module documentation, "Deciding the hypotheses").
fn negative_axis_sign_hypothesis(
    name: &str,
    statement: &str,
    label: &str,
    c: &Bivariate,
    f: &[Q],
    first: usize,
    strict: bool,
) -> LiftHypothesis {
    let relation = if strict { "< 0" } else { "<= 0" };
    let Some(degree) = x_degree(c) else {
        return if strict {
            exact_hypothesis(
                name,
                statement,
                Err((Some(first), format!("{label} = 0 is not < 0"))),
            )
        } else {
            exact_hypothesis(name, statement, Ok(format!("{label} = 0")))
        };
    };
    if let Some((gamma, w)) = separate(c) {
        let isolation = isolate_real_roots_rational_coeffs(&w);
        let zero = Q::zero();
        let negative_root = isolation
            .roots
            .iter()
            .any(|root| isolation.compare_root(root, &zero) == Ordering::Less);
        if !negative_root {
            // w has one sign on (-inf, 0); it is the sign at -1.
            let w_positive = vec_poly::evaluate(&w, &-Q::one()).is_positive();
            let condition = match (w_positive, strict) {
                (true, false) => IndexSignCondition::Nonpositive,
                (true, true) => IndexSignCondition::Negative,
                (false, false) => IndexSignCondition::Nonnegative,
                (false, true) => IndexSignCondition::Positive,
            };
            let mut hypothesis = sign_hypothesis(
                name,
                statement,
                "gamma(n)/f(n)",
                &IndexFunction::new(gamma, f.to_vec()),
                first,
                condition,
            );
            hypothesis.evidence = format!(
                "{label}(x) = gamma(n) w(x) / f(n) with w(x) = {}, which is {} on (-inf, 0) (no \
                 root there, exact): {}",
                format_index_polynomial(&w, "x"),
                if w_positive { "positive" } else { "negative" },
                hypothesis.evidence
            );
            return hypothesis;
        }
    }
    // Sufficient test after x = -t: (-1)^k [x^k] c / f <= 0 for every k, and
    // for the strict form < 0 for some k.
    let mut evidence = Vec::new();
    let mut strict_term = None;
    for k in 0..=degree {
        let column = x_column(c, k);
        if column.is_empty() {
            continue;
        }
        let function = IndexFunction::new(column, f.to_vec());
        let (weak, strong) = if k % 2 == 0 {
            (
                IndexSignCondition::Nonpositive,
                IndexSignCondition::Negative,
            )
        } else {
            (
                IndexSignCondition::Nonnegative,
                IndexSignCondition::Positive,
            )
        };
        let coefficient_label = format!("[x^{k}] {label}");
        let decision = function.decide(first, weak);
        if !decision.holds() {
            return unresolved_hypothesis(
                name,
                statement,
                format!(
                    "{label} is not separable with a factor of one sign on (-inf, 0), and the \
                     sufficient monomial test (-1)^k [x^k] {label} <= 0 fails for k = {k} (not a \
                     disproof): {}",
                    decision.describe(&coefficient_label)
                ),
            );
        }
        evidence.push(decision.describe(&coefficient_label));
        if strict && strict_term.is_none() {
            let decision = function.decide(first, strong);
            if decision.holds() {
                strict_term = Some(decision.describe(&coefficient_label));
            }
        }
    }
    if strict {
        match strict_term {
            Some(term) => evidence.push(format!("strict: {term}")),
            None => {
                return unresolved_hypothesis(
                    name,
                    statement,
                    format!(
                        "the monomial test gives {label} <= 0 on (-inf, 0) but no term \
                         (-1)^k [x^k] {label} is < 0 for every n (not a disproof): {}",
                        evidence.join("; ")
                    ),
                )
            }
        }
    }
    exact_hypothesis(
        name,
        statement,
        Ok(format!(
            "monomial test after x = -t, every (-1)^k [x^k] {label} <= 0 ({label} {relation}): {}",
            evidence.join("; ")
        )),
    )
}

/// Hypothesis 6 and the Lean route it selects: `V_n < 0` on `(-inf, 0)` for
/// every `n` (`halfGrowth`), or `W_n < 0` there (`sign` when `V = 0`, the
/// `LiuWang` step otherwise).
fn strict_sign_hypothesis(
    c1: &Bivariate,
    c2: &Bivariate,
    f: &[Q],
    first: usize,
) -> (LiftHypothesis, &'static str) {
    let name = "derivative_or_lag_coefficient_negative";
    let statement = format!(
        "for every n >= {first} and x < 0, V_n(x) < 0 or W_n(x) < 0 (V_n < 0 on (-inf, 0) for \
         every n, or W_n < 0 there for every n)"
    );
    let derivative_present = x_degree(c1).is_some();
    let lag_present = x_degree(c2).is_some();
    let derivative = derivative_present
        .then(|| negative_axis_sign_hypothesis(name, &statement, "V_n", c1, f, first, true));
    if let Some(hypothesis) = &derivative {
        if hypothesis.status == LiftHypothesisStatus::Verified || !lag_present {
            return (hypothesis.clone(), LEAN_ROUTE_HALF_GROWTH);
        }
    }
    let lag = negative_axis_sign_hypothesis(name, &statement, "W_n", c2, f, first, true);
    match derivative {
        None => (lag, LEAN_ROUTE_SIGN),
        Some(_) if lag.status == LiftHypothesisStatus::Verified => (lag, LEAN_ROUTE_LIU_WANG),
        Some(derivative) => (
            unresolved_hypothesis(
                name,
                &statement,
                format!(
                    "neither holds uniformly (not a disproof): V_n: {}; W_n: {}",
                    derivative.evidence, lag.evidence
                ),
            ),
            LEAN_ROUTE_HALF_GROWTH,
        ),
    }
}

/// The RealRooted declarations mirrored by `route`.
fn lean_declarations(route: &str) -> Vec<String> {
    let names: &[&str] = match route {
        LEAN_ROUTE_SIGN => &[
            "RealRooted.threeTerm_interlaces_of_nonnegCoeffs",
            "RealRooted.LiuWang.strictInterl_nonpos_lag_sequence",
        ],
        LEAN_ROUTE_LIU_WANG => &["RealRooted.LiuWang.strictInterl_derivative_lag_of_eval_nonpos"],
        _ => &[
            "RealRooted.derivLag_strictInterl_of_nonnegCoeffs",
            "RealRooted.derivLag_strictInterl_and_noCommonRoot",
            "RealRooted.derivLag_eval_zero_pos",
        ],
    };
    names.iter().map(|name| (*name).to_string()).collect()
}

/// The requested region must contain `(-inf, 0)` (`(-inf, 0]` with a gauge).
fn requested_region_hypothesis(
    region: &RealInterval,
    gauge: usize,
    current_index: usize,
) -> LiftHypothesis {
    let certified = if gauge == 0 { "(-inf, 0)" } else { "(-inf, 0]" };
    let statement = format!(
        "the requested region {} contains {certified}, the region this theorem certifies",
        open_root_interval(region)
    );
    let upper_ok = region.upper.as_ref().is_none_or(|upper| {
        if gauge == 0 {
            !upper.is_negative()
        } else {
            upper.is_positive()
        }
    });
    if region.lower.is_none() && upper_ok {
        exact_hypothesis("requested_root_region", &statement, Ok("exact".into()))
    } else if gauge > 0 && region.lower.is_none() {
        exact_hypothesis(
            "requested_root_region",
            &statement,
            Err((
                Some(current_index),
                format!(
                    "x^{gauge} divides P_{current_index}, so 0 is a root outside the open region \
                     (exact)"
                ),
            )),
        )
    } else {
        unresolved_hypothesis(
            "requested_root_region",
            &statement,
            format!("this theorem certifies the region {certified} only"),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::parent_proof::{prove_parent, ParentRecord};
    use super::*;

    /// `prove_parent` settles the catalog derivative-lag parents with the new
    /// route (the earlier routes do not apply to a rule with a derivative and
    /// a lag term).
    #[test]
    fn prove_parent_uses_the_derivative_lag_route() {
        for id in ["A144438", "A046739", "A271697"] {
            let entry = crate::oeis::by_id(id).expect("catalog entry");
            let (rule, seeds) = entry.recurrence_parts().expect("catalog recurrence");
            let record = ParentRecord {
                rule: &rule,
                seeds: &seeds,
                first_index: entry.recurrence_first_index,
            };
            let proof = prove_parent(&record, entry.recurrence_first_index, &[], false)
                .expect("the provers run");
            assert_eq!(
                proof.route,
                Some(VeroneseParentProofRoute::LiuWangDerivativeLag),
                "{id}: {proof:#?}"
            );
            assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
        }
    }
}
