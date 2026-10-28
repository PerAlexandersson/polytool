//! First-order derivative recurrences `P_n = a_n(x) P_(n-1) + b_n(x) P'_(n-1)`
//! by the interlacing criterion of Liu and Wang, with an exact root-location
//! induction in place of their nonnegativity hypothesis.
//!
//! # The theorem
//!
//! L. L. Liu and Y. Wang, *A unified approach to polynomial sequences with
//! only real zeros*, Adv. in Appl. Math. 38 (2007) 542--560, Theorem 2.1: let
//! `F = a(x) f + b(x) g` with `deg F` equal to `deg f` or `deg f + 1`, `f` and
//! `g` real-rooted with `g ⪯ f`, and `F`, `g` with leading coefficients of
//! one sign.  If `b(r) <= 0` at every zero `r` of `f`, then `F` is
//! real-rooted and `f ⪯ F`.  (`g ⪯ f` means that the zeros alternate with the
//! largest zero belonging to `f`: `... <= s_2 <= r_2 <= s_1 <= r_1`.)  Their
//! Corollary 2.4 applies it with `g = f' ⪯ f` to sequences with nonnegative
//! coefficients; here nonnegativity is replaced by the following induction,
//! whose every step is an exact sign condition on rational functions of `n`.
//!
//! # The induction
//!
//! Let `P_(N-1)` be nonzero and real-rooted with every root `< 0`, of degree
//! `d`, and for `n >= N` write `a_n(x) = a_(n,0) + a_(n,1) x`,
//! `b_n(x) = b_(n,0) + b_(n,1) x + b_(n,2) x^2` (the rule
//! `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) D P_(n-1)` divided by `f(n)`) and
//! `d_(n-1) = d + n - N`.  Assume for every `n >= N`:
//!
//! 1. `lambda_n = a_(n,1) + d_(n-1) b_(n,2) > 0`; then `deg P_n = d_(n-1) + 1`
//!    and `lead P_n = lambda_n lead P_(n-1)` has the sign of
//!    `lead P'_(n-1)` (for `d_(n-1) >= 1`);
//! 2. `b_(n,0) = 0`, that is `b_n(0) = 0`;
//! 3. `b_n(x) <= 0` for every `x < 0`;
//! 4. `a_n(0) > 0`.
//!
//! If every root of `P_(n-1)` is `< 0`, Theorem 2.1 with `f = P_(n-1)`,
//! `g = P'_(n-1)` (Rolle gives `g ⪯ f`) and condition 3 at the zeros of `f`
//! makes `P_n` real-rooted with `P_(n-1) ⪯ P_n`; for a constant `P_(n-1)`,
//! `P_n = a_n P_(n-1)` is linear.  Interlacing leaves at most the largest
//! root `t_1` of `P_n` at or right of `0`, and by condition 2
//!
//! ```text
//! P_n(0) / lead P_n = (a_n(0) / lambda_n) * P_(n-1)(0) / lead P_(n-1) > 0,
//! ```
//!
//! since the last factor is `prod (-s_i) > 0`; with one root `t_1 >= 0` the
//! left side would be `<= 0`.  So every root of `P_n` is `< 0`.
//!
//! Conditions 1 and 4 are sign conditions on rational functions of `n`
//! ([`crate::decide_index_rational_function_sign`]); condition 2 is checked
//! on the coefficient polynomial; condition 3 is decided for a separable
//! `b_n(x) = gamma(n) w(x) / f(n)` (`w` without roots in `(-inf, 0)`, so of
//! one sign there, and the sign of `gamma / f` by the sign checker; a failure
//! is a disproof at that `n`) and otherwise by the sufficient monomial test
//! `(-1)^k b_(n,k) <= 0` for every `k` (a failure is not a disproof).

use super::*;
use crate::real_rootedness::isolate_real_roots_rational_coeffs;
use std::cmp::Ordering;

const REFERENCE: &str = "L. L. Liu and Y. Wang, A unified approach to polynomial sequences with \
     only real zeros, Adv. in Appl. Math. 38 (2007) 542-560, Theorem 2.1 with g = f' (as in \
     their Corollary 2.4), with nonnegativity of the coefficients replaced by an exact \
     induction keeping every root in (-inf, 0)";

/// Check the Liu--Wang theorem (module documentation of the `liu_wang`
/// submodule, summarized here) for the family given by `rule`,
/// `first_index` and the seed rows `seeds` (the rule applies from
/// `N = first_index + seeds.len()` on).
///
/// The rule must be `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1)` with
/// `f` free of `x`, `deg_x c_0 <= 1` and `deg_x c_1 <= 2` (either term may be
/// absent).  With `a_n = c_0 / f`, `b_n = c_1 / f`, `d = deg P_(N-1)`, the
/// hypotheses, each decided exactly, are:
///
/// - `f(n) != 0` for `n >= N` (sign checker);
/// - `P_(N-1)` is nonzero and real-rooted with every root `< 0`, and so is
///   every earlier seed row (exact);
/// - `lambda_n = [x^1] a_n + (d + n - N) [x^2] b_n > 0` for `n >= N` (sign
///   checker);
/// - `b_n(0) = 0` for `n >= N` (the polynomial `[x^0] c_1` is zero);
/// - `b_n(x) <= 0` for `x < 0` and `n >= N` (separable form and sign checker,
///   or the sufficient monomial test);
/// - `a_n(0) > 0` for `n >= N` (sign checker);
/// - a requested root region must contain `(-inf, 0)`.
///
/// When all hold, every `P_n`, `n >= first_index`, is real-rooted with every
/// root in `(-inf, 0)` (so its coefficients have one sign),
/// `deg P_n = d + n - N + 1` for `n >= N`, and `P_(n-1) ⪯ P_n` for `n >= N`.
///
/// # Errors
///
/// None at present; the `Result` is kept for symmetry with
/// [`check_three_term_recurrence_theorem`].
///
/// # Example
///
/// ```
/// use num_bigint::BigInt;
/// use polytool::recurrence::{
///     check_liu_wang_theorem, BigRational, BivarPoly, Recurrence, RecurrenceTerm,
///     RealRootednessTheoremOptions,
/// };
///
/// let q = |v: i64| BigRational::from_integer(BigInt::from(v));
/// // Type-B Eulerian polynomials: B_n = (1 + (2n - 1) x) B_(n-1) + 2 x (1 - x) B'_(n-1).
/// let rule = Recurrence {
///     terms: vec![
///         RecurrenceTerm::new(
///             1,
///             0,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(-1)], vec![q(0), q(2)]]),
///         ),
///         RecurrenceTerm::new(
///             1,
///             1,
///             BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(2), q(-2)]]),
///         ),
///     ],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let report =
///     check_liu_wang_theorem(&rule, 0, &[vec![q(1)]], &RealRootednessTheoremOptions::default())
///         .unwrap();
/// assert!(report.is_proved());
/// report.verify().unwrap();
/// ```
pub fn check_liu_wang_theorem(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<BigRational>],
    options: &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError> {
    let mut report = RealRootednessTheoremReport::new(
        RealRootednessTheorem::LiuWang,
        REFERENCE,
        rule,
        first_index,
        seeds,
        options,
    );
    if seeds.is_empty() {
        return Ok(report.not_applicable("the theorem needs a seed row"));
    }
    let shape = match recurrence_shape(rule) {
        Ok(shape) => shape,
        Err(reason) => return Ok(report.not_applicable(reason)),
    };
    if shape.terms.is_empty()
        || shape
            .terms
            .keys()
            .any(|key| *key != (1, 0) && *key != (1, 1))
    {
        return Ok(report
            .not_applicable("the rule is not f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1)"));
    }
    let empty = Bivariate::new();
    let c0 = shape.terms.get(&(1, 0)).unwrap_or(&empty);
    let c1 = shape.terms.get(&(1, 1)).unwrap_or(&empty);
    if x_degree(c0).is_some_and(|d| d > 1) || x_degree(c1).is_some_and(|d| d > 2) {
        return Ok(report.not_applicable(
            "deg_x c_0 > 1 or deg_x c_1 > 2: the degree may grow by more than one",
        ));
    }
    let n_first = first_index + seeds.len();
    let start = n_first - 1;
    let seed = &seeds[seeds.len() - 1];
    let d = vec_poly::degree(seed).unwrap_or(0);
    let f = &shape.denominator;
    let over_f = |p: Vec<Q>| IndexFunction::new(p, f.clone());
    report.normal_form = format!(
        "P_n = a_n(x) P_(n-1) + b_n(x) P'_(n-1) for n >= {n_first}, with a_n(x) = ({}) + ({}) x \
         and b_n(x) = ({}) + ({}) x + ({}) x^2",
        over_f(x_column(c0, 0)).text(),
        over_f(x_column(c0, 1)).text(),
        over_f(x_column(c1, 0)).text(),
        over_f(x_column(c1, 1)).text(),
        over_f(x_column(c1, 2)).text(),
    );

    let negative_axis = RealInterval {
        lower: None,
        lower_closed: false,
        upper: Some(Q::zero()),
        upper_closed: false,
    };
    let mut hypotheses: Vec<LiftHypothesis> = Vec::new();
    hypotheses.extend(denominator_hypothesis(f, n_first));
    hypotheses.push(exact_hypothesis(
        "seed_roots_negative",
        &format!("P_{start} is nonzero and real-rooted with every root < 0"),
        if has_real_roots_below(seed, &Q::zero()) {
            Ok(format!("P_{start} has degree {d}; checked exactly"))
        } else {
            Err((
                Some(start),
                format!("P_{start} is zero, not real-rooted, or has a root >= 0 (exact)"),
            ))
        },
    ));
    hypotheses.extend(earlier_seeds_hypothesis(
        &seeds[..seeds.len() - 1],
        first_index,
        Some(&negative_axis),
    ));
    // lambda_n = a_(n,1) + (d + n - N) b_(n,2), over f.
    let previous_degree = vec![
        Q::from_integer(BigInt::from(d as i64 - n_first as i64)),
        Q::one(),
    ];
    let lambda = over_f(vec_poly::add(
        &vec_poly::trim(&x_column(c0, 1)),
        &vec_poly::mul(&previous_degree, &vec_poly::trim(&x_column(c1, 2))),
    ));
    hypotheses.push(sign_hypothesis(
        "degree_step_positive",
        &format!(
            "lambda_n = [x^1] a_n + (deg P_(n-1)) [x^2] b_n > 0 for every n >= {n_first} \
             (deg P_n = deg P_(n-1) + 1 with leading coefficients of one sign)"
        ),
        "lambda_n",
        &lambda,
        n_first,
        IndexSignCondition::Positive,
    ));
    let b_at_zero = x_column(c1, 0);
    let statement = format!("b_n(0) = 0 for every n >= {n_first}");
    hypotheses.push(if b_at_zero.is_empty() {
        exact_hypothesis(
            "derivative_coefficient_vanishes_at_zero",
            &statement,
            Ok("[x^0] c_1(n, x) is the zero polynomial".into()),
        )
    } else {
        // A nonzero polynomial in n has at most deg roots.
        let witness = (n_first..=n_first + b_at_zero.len()).find(|&n| {
            !vec_poly::evaluate(&b_at_zero, &Q::from_integer(BigInt::from(n))).is_zero()
        });
        exact_hypothesis(
            "derivative_coefficient_vanishes_at_zero",
            &statement,
            Err((
                witness,
                format!(
                    "[x^0] c_1(n, x) = {} is not the zero polynomial",
                    format_index_polynomial(&b_at_zero, "n")
                ),
            )),
        )
    });
    hypotheses.push(derivative_coefficient_hypothesis(c1, f, n_first));
    hypotheses.push(sign_hypothesis(
        "constant_term_positive",
        &format!("a_n(0) > 0 for every n >= {n_first}"),
        "a_n(0)",
        &over_f(x_column(c0, 0)),
        n_first,
        IndexSignCondition::Positive,
    ));
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
    let root_region = open_root_interval(&negative_axis);
    Ok(report.finish(RealRootednessTheoremConclusion {
        first_index,
        interlacing_from: Some(n_first),
        interlacing: Some(RowInterlacing::Weak),
        simple_roots: false,
        statement: format!(
            "every P_n, n >= {first_index}, is nonzero and real-rooted with every root in \
             {root_region} (so its coefficients have one sign); deg P_n = deg P_(n-1) + 1 and \
             P_(n-1) interlaces P_n (Liu--Wang) for n >= {n_first}"
        ),
        root_region: Some(root_region),
    }))
}

/// `b_n(x) <= 0` for every `x < 0` and `n >= first`.
fn derivative_coefficient_hypothesis(c1: &Bivariate, f: &[Q], first: usize) -> LiftHypothesis {
    let name = "derivative_coefficient_nonpositive";
    let statement = format!("b_n(x) <= 0 for every x < 0 and every n >= {first}");
    if x_degree(c1).is_none() {
        return exact_hypothesis(name, &statement, Ok("b_n = 0".into()));
    }
    if let Some((gamma, w)) = separate(c1) {
        let isolation = isolate_real_roots_rational_coeffs(&w);
        let zero = Q::zero();
        let negative_root = isolation
            .roots
            .iter()
            .any(|root| isolation.compare_root(root, &zero) == Ordering::Less);
        if !negative_root {
            // w has one sign on (-inf, 0); it is the sign at -1.
            let w_sign = vec_poly::evaluate(&w, &-Q::one());
            let condition = if w_sign.is_positive() {
                IndexSignCondition::Nonpositive
            } else {
                IndexSignCondition::Nonnegative
            };
            let mut hypothesis = sign_hypothesis(
                name,
                &statement,
                "gamma(n)/f(n)",
                &IndexFunction::new(gamma, f.to_vec()),
                first,
                condition,
            );
            hypothesis.evidence = format!(
                "b_n(x) = gamma(n) w(x) / f(n) with w(x) = {}, which is {} on (-inf, 0) (no \
                 root there, exact): {}",
                format_index_polynomial(&w, "x"),
                if w_sign.is_positive() {
                    "positive"
                } else {
                    "negative"
                },
                hypothesis.evidence
            );
            return hypothesis;
        }
    }
    // Sufficient monomial test: x^k has the sign (-1)^k on (-inf, 0).
    let mut evidence = Vec::new();
    for k in 0..=x_degree(c1).unwrap_or(0) {
        let column = IndexFunction::new(x_column(c1, k), f.to_vec());
        let condition = if k % 2 == 0 {
            IndexSignCondition::Nonpositive
        } else {
            IndexSignCondition::Nonnegative
        };
        let decision = column.decide(first, condition);
        let label = format!("[x^{k}] b_n");
        if !decision.holds() {
            return unresolved_hypothesis(
                name,
                &statement,
                format!(
                    "b_n is not separable with a factor of one sign on (-inf, 0), and the \
                     sufficient monomial test (-1)^k [x^k] b_n <= 0 fails for k = {k} (not a \
                     disproof): {}",
                    decision.describe(&label)
                ),
            );
        }
        evidence.push(decision.describe(&label));
    }
    exact_hypothesis(
        name,
        &statement,
        Ok(format!(
            "monomial test, every (-1)^k [x^k] b_n <= 0: {}",
            evidence.join("; ")
        )),
    )
}
