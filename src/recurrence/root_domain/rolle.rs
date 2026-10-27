//! Rolle's theorem for derivatives of products as a recurrence induction,
//! with index-dependent factors and index-dependent intervals.
//!
//! A first-order rule `Q_m = D^j (G(m, x) Q_(m-1))` is matched when the
//! coefficient of `D^k Q_(m-1)` is `C(j,k) d^(j-k)/dx^(j-k) G(m, x)` for every
//! `k`.  The factor is either separable, `G = s(m) g(x)`, or
//! `G = g(x) (c(m) + l(m) x)` with an index-dependent linear part (found by
//! taking the `x`-content `g` of the rows of `G`).
//!
//! An interval candidate `I_n = [l_0 + l_1 n, u_0 + u_1 n]` (ends optional,
//! open or closed) closes the induction when
//!
//! 1. `I_(m-1) ⊆ I_m` for every `m` (`l_1 <= 0` and `u_1 >= 0`; automatic for
//!    index-independent intervals), checked symbolically;
//! 2. `g` is real-rooted in `I_(m_0)`, `m_0` the first target (hence in every
//!    later `I_m`), checked by exact isolation;
//! 3. for a linear part, `l(m)` has a certified sign for `m >= m_0` and the
//!    root `-c(m)/l(m)` lies in `I_m` for every `m >= m_0`: both are polynomial
//!    inequalities in `m`, certified when every Taylor coefficient of the
//!    polynomial at `m_0` has the right sign (a sufficient criterion; a
//!    failure only means "not certified");
//! 4. every seed `Q_k` is real-rooted in `I_k`, checked exactly.
//!
//! Then `G(m, x) Q_(m-1)` is zero or real-rooted with every root in `I_m`, and
//! by Rolle's theorem so is `D^j (G(m, x) Q_(m-1))`.  Steps that cannot be
//! checked symbolically are reported as failed attempts, never as proofs.

use super::*;
use crate::polynomial::Polynomial;
use crate::real_rootedness::{isolate_real_roots_exact, RealInterval, RealRootIsolation};
use crate::vec_poly;
use num_traits::{One, Signed, Zero};
use std::cmp::Ordering;

pub(super) const ROLLE_STATEMENT: &str =
    "if every root of g and of Q_(m-1) is real and lies in an interval I, then \
     s(m) D^j (g Q_(m-1)) is zero or real-rooted with every root in I (Rolle's theorem; the \
     scalar s(m) does not move roots); with the seeds checked exactly this closes an \
     induction; for an index-dependent factor g(m,x) and interval I_m the step needs \
     I_(m-1) ⊆ I_m and the roots of g(m,x) in I_m for every m";

/// The index-dependent linear part `c(m) + l(m) x` of a Rolle factor.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct LinearIndexFactor {
    /// `c(m)`, ascending in the target index `m`.
    pub(super) constant: Vec<Q>,
    /// `l(m)`, ascending in the target index `m`.
    pub(super) slope: Vec<Q>,
}

/// `Q_m = s(m) D^j (g Q_(m-1))`, or `Q_m = D^j (g (c(m) + l(m) x) Q_(m-1))`
/// when `linear` is present (then `scalar` is `[1]`).
#[derive(Debug, Clone)]
pub(super) struct RolleForm {
    pub(super) factor: Vec<Q>,
    pub(super) order: usize,
    pub(super) scalar: Vec<Q>,
    pub(super) linear: Option<LinearIndexFactor>,
}

fn trim_rows(mut rows: Vec<Vec<Q>>) -> Vec<Vec<Q>> {
    for row in &mut rows {
        vec_poly::strip_trailing_zeros(row);
    }
    while rows.last().is_some_and(Vec::is_empty) {
        rows.pop();
    }
    rows
}

pub(super) fn match_rolle_form(rule: &Recurrence) -> Result<RolleForm, String> {
    let rule = canonical_recurrence(rule);
    if rule.denominator.is_some() {
        return Err("the left factor f(n,x) is not a nonzero constant".into());
    }
    if rule
        .inhomogeneous
        .as_ref()
        .is_some_and(|g| !bivar_rows(g).is_empty())
    {
        return Err("the rule has an additive term".into());
    }
    if rule.terms.is_empty() {
        return Err("the rule has no terms".into());
    }
    if rule.terms.iter().any(|t| t.offset != 1) {
        return Err("only first-order rules (every offset 1) are matched".into());
    }
    if rule.terms.iter().any(RecurrenceTerm::is_alternating) {
        return Err("(-1)^n factors are not matched".into());
    }
    let order = rule.terms.iter().map(|t| t.deriv_order).max().unwrap_or(0);
    let coefficient = |k: usize| -> Vec<Vec<Q>> {
        rule.terms
            .iter()
            .find(|t| t.deriv_order == k)
            .map(|t| bivar_rows(&t.coeff))
            .unwrap_or_default()
    };
    let top = coefficient(order);
    if top.is_empty() {
        return Err("the top coefficient vanishes".into());
    }
    // Leibniz: the coefficient of D^k is C(j,k) d^(j-k)/dx^(j-k) G(m,x).
    let mut derivative = top.clone();
    for k in (0..order).rev() {
        derivative = trim_rows(derivative.iter().map(|r| vec_poly::derivative(r)).collect());
        let expected = trim_rows(
            derivative
                .iter()
                .map(|r| vec_poly::scale(r, &binomial(order, k)))
                .collect(),
        );
        if coefficient(k) != expected {
            return Err(format!(
                "the coefficient of D^{k} is not C({order},{k}) times the {}-th x-derivative of \
                 the top coefficient g(n,x)",
                order - k
            ));
        }
    }
    // Separable: G = s(m) g(x).
    let factor = top
        .iter()
        .find(|row| !row.is_empty())
        .cloned()
        .expect("nonzero top coefficient");
    let lead = factor.last().expect("nonzero factor").clone();
    let mut scalar = Vec::with_capacity(top.len());
    let separable = top.iter().all(|row| {
        let s = row.last().map_or_else(Q::zero, |c| c / &lead);
        let ok = vec_poly::trim(row) == vec_poly::trim(&vec_poly::scale(&factor, &s));
        scalar.push(s);
        ok
    });
    if separable {
        vec_poly::strip_trailing_zeros(&mut scalar);
        return Ok(RolleForm {
            factor,
            order,
            scalar,
            linear: None,
        });
    }
    // G = g(x) (c(m) + l(m) x) with g the x-content of the rows.
    let content = top
        .iter()
        .filter(|row| !row.is_empty())
        .fold(Polynomial::<Q>::zero(), |acc, row| {
            acc.gcd(&Polynomial::new(row.clone()))
        });
    let mut constant = Vec::with_capacity(top.len());
    let mut slope = Vec::with_capacity(top.len());
    for row in &top {
        let (quotient, remainder) = Polynomial::new(row.clone()).div_rem(&content);
        debug_assert!(remainder.is_zero());
        let quotient = vec_poly::trim(quotient.coeffs());
        if vec_poly::degree(&quotient).is_some_and(|d| d > 1) {
            return Err(
                "the top coefficient g(n,x) is neither s(n) g(x) nor g(x) times a \
                        factor linear in x with index-dependent coefficients"
                    .into(),
            );
        }
        constant.push(vec_poly::coeff(&quotient, 0));
        slope.push(vec_poly::coeff(&quotient, 1));
    }
    vec_poly::strip_trailing_zeros(&mut constant);
    vec_poly::strip_trailing_zeros(&mut slope);
    Ok(RolleForm {
        factor: vec_poly::trim(content.coeffs()),
        order,
        scalar: vec![Q::one()],
        linear: Some(LinearIndexFactor { constant, slope }),
    })
}

// ---------------------------------------------------------------------------
// Index-dependent intervals
// ---------------------------------------------------------------------------

/// `(value, per_index)` of an endpoint `value + per_index * n`.
fn endpoint_parts(end: &RootEndpoint) -> Result<(Q, Q), String> {
    let value = parse_rational_coeff(&end.value)?;
    let slope = match &end.per_index {
        Some(p) => parse_rational_coeff(p)?,
        None => Q::zero(),
    };
    Ok((value, slope))
}

/// The maximum (`max`) or minimum of `value + slope n` over
/// `first <= n (<= last)`; `None` when it is unbounded.
fn affine_extreme(value: &Q, slope: &Q, first: usize, last: Option<usize>, max: bool) -> Option<Q> {
    let at = |n: usize| value + slope * Q::from_integer(BigInt::from(n));
    let grows = if max {
        slope.is_positive()
    } else {
        slope.is_negative()
    };
    if grows {
        last.map(at)
    } else {
        Some(at(first))
    }
}

/// The smallest index-independent interval containing every `I_n` with
/// `first <= n (<= last)` (a superset of their union).
pub(super) fn domain_hull(
    interval: &RootInterval,
    first: usize,
    last: Option<usize>,
) -> Result<RealInterval, String> {
    let end = |end: &Option<RootEndpoint>, max: bool| -> Result<(Option<Q>, bool), String> {
        match end {
            None => Ok((None, false)),
            Some(e) => {
                let (value, slope) = endpoint_parts(e)?;
                let extreme = affine_extreme(&value, &slope, first, last, max);
                let closed = extreme.is_some() && e.closed;
                Ok((extreme, closed))
            }
        }
    };
    let (lower, lower_closed) = end(&interval.lower, false)?;
    let (upper, upper_closed) = end(&interval.upper, true)?;
    Ok(RealInterval {
        lower,
        lower_closed,
        upper,
        upper_closed,
    })
}

/// Whether `I_(n-1) ⊆ I_n` for every `n`: the lower end does not increase and
/// the upper end does not decrease.
pub(super) fn grows_with_index(interval: &RootInterval) -> Result<bool, String> {
    let lower = match &interval.lower {
        Some(e) => !endpoint_parts(e)?.1.is_positive(),
        None => true,
    };
    let upper = match &interval.upper {
        Some(e) => !endpoint_parts(e)?.1.is_negative(),
        None => true,
    };
    Ok(lower && upper)
}

/// Whether `point` lies in every `I_n` with `first <= n (<= last)`.
pub(super) fn contains_point_throughout(
    interval: &RootInterval,
    point: &Q,
    first: usize,
    last: Option<usize>,
) -> Result<bool, String> {
    let lower_ok = match &interval.lower {
        None => true,
        Some(e) => {
            let (value, slope) = endpoint_parts(e)?;
            match affine_extreme(&value, &slope, first, last, true) {
                None => false,
                Some(high) => high < *point || (high == *point && e.closed),
            }
        }
    };
    let upper_ok = match &interval.upper {
        None => true,
        Some(e) => {
            let (value, slope) = endpoint_parts(e)?;
            match affine_extreme(&value, &slope, first, last, false) {
                None => false,
                Some(low) => low > *point || (low == *point && e.closed),
            }
        }
    };
    Ok(lower_ok && upper_ok)
}

/// The exact image of `I_n` under `y -> alpha y + beta` (`alpha != 0`):
/// affine endpoints stay affine in `n`.
pub(super) fn affine_interval_image(
    interval: &RootInterval,
    alpha: &Q,
    beta: &Q,
) -> Result<RootInterval, String> {
    let map = |end: &Option<RootEndpoint>| -> Result<Option<RootEndpoint>, String> {
        end.as_ref()
            .map(|e| {
                let (value, slope) = endpoint_parts(e)?;
                let slope = alpha * slope;
                Ok(RootEndpoint {
                    value: format_rational_coeff(&(alpha * value + beta)),
                    per_index: (!slope.is_zero()).then(|| format_rational_coeff(&slope)),
                    closed: e.closed,
                })
            })
            .transpose()
    };
    let (lower, upper) = (map(&interval.lower)?, map(&interval.upper)?);
    Ok(if alpha.is_negative() {
        RootInterval {
            lower: upper,
            upper: lower,
        }
    } else {
        RootInterval { lower, upper }
    })
}

// ---------------------------------------------------------------------------
// The induction
// ---------------------------------------------------------------------------

/// Certify that the root `-c(m)/l(m)` of the linear part lies in `I_m` for
/// every `m >= start`.
fn certify_linear_root(
    linear: &LinearIndexFactor,
    interval: &RootInterval,
    start: usize,
) -> Result<String, String> {
    let (constant, slope) = match index_polynomial_sign_from(&linear.slope, start) {
        Some(Ordering::Greater) => (linear.constant.clone(), linear.slope.clone()),
        Some(Ordering::Less) => (
            vec_poly::neg(&linear.constant),
            vec_poly::neg(&linear.slope),
        ),
        _ => {
            return Err(format!(
                "the coefficient l(m) = {} of x in the index-dependent factor is not certified \
                 nonzero for m >= {start}",
                poly_text(&linear.slope, "m")
            ))
        }
    };
    let root = format!(
        "-({})/({})",
        poly_text(&linear.constant, "m"),
        poly_text(&linear.slope, "m")
    );
    if let Some(lower) = &interval.lower {
        let (value, per_index) = endpoint_parts(lower)?;
        // -c(m) - lower(m) l(m) >= 0 (> 0 for an open end), l(m) > 0.
        let excess = vec_poly::sub(
            &vec_poly::neg(&constant),
            &vec_poly::mul(&[value, per_index], &slope),
        );
        if !index_polynomial_nonnegative_from(&excess, start, !lower.closed) {
            return Err(format!(
                "the root {root} of the index-dependent factor is not certified to lie above \
                 the lower end of I_m for every m >= {start} (Taylor-shift criterion)"
            ));
        }
    }
    if let Some(upper) = &interval.upper {
        let (value, per_index) = endpoint_parts(upper)?;
        // upper(m) l(m) + c(m) >= 0 (> 0 for an open end).
        let excess = vec_poly::add(&vec_poly::mul(&[value, per_index], &slope), &constant);
        if !index_polynomial_nonnegative_from(&excess, start, !upper.closed) {
            return Err(format!(
                "the root {root} of the index-dependent factor is not certified to lie below \
                 the upper end of I_m for every m >= {start} (Taylor-shift criterion)"
            ));
        }
    }
    Ok(format!(
        "the root {root} of the index-dependent factor lies in I_m for every m >= {start} \
         (Taylor-shift certificates)"
    ))
}

/// The inputs of one Rolle induction.
pub(super) struct RolleSetup<'a> {
    pub(super) rule: &'a Recurrence,
    /// Row name (`Q`, `P`, ...).
    pub(super) name: String,
    /// First seed index.
    pub(super) first: usize,
    /// First target index.
    pub(super) valid_from: usize,
    /// Seeds `(index, row)` for `first..valid_from`.
    pub(super) seeds: Vec<(usize, Vec<Q>)>,
    /// Supplied and default candidates, tried before the root hull.
    pub(super) candidates: Vec<(RootInterval, String)>,
    /// Largest number of attempts.
    pub(super) limit: usize,
    /// Degree schedule of the rows, for the degree formula.
    pub(super) schedule: Option<&'a DegreeSchedule>,
}

/// The theorem report, the proven interval (if any) and whether the degree
/// formula gives `deg = d_n` on the schedule.
pub(super) struct RolleOutcome {
    pub(super) theorem: RecurrenceRootTheorem,
    pub(super) proven: Option<RootInterval>,
    pub(super) full_degree: bool,
}

fn check_attempt(
    form: &RolleForm,
    factor_isolation: &RealRootIsolation,
    seed_isolations: &[(usize, RealRootIsolation)],
    interval: &RootInterval,
    valid_from: usize,
    first: usize,
    name: &str,
) -> Result<String, String> {
    let dependent = !interval.is_index_independent();
    let mut steps = Vec::new();
    if dependent {
        if !grows_with_index(interval)? {
            return Err(format!(
                "{interval} does not grow with the index (I_(m-1) is not contained in I_m), so \
                 the Rolle step is not certified"
            ));
        }
        steps.push("I_(m-1) ⊆ I_m for every m (lower slope <= 0, upper slope >= 0)".to_string());
    }
    let start = interval.at(valid_from)?;
    if !factor_isolation.real_rooted_in(&start) {
        let place = if dependent {
            format!("I_{valid_from} = {start}")
        } else {
            interval.to_string()
        };
        return Err(format!(
            "g = {} is not real-rooted in {place}",
            poly_text(&form.factor, "x")
        ));
    }
    if dependent {
        steps.push(format!(
            "g is real-rooted in I_{valid_from} and hence in every later I_m"
        ));
    }
    if let Some(linear) = &form.linear {
        steps.push(certify_linear_root(linear, interval, valid_from)?);
    }
    for (k, isolation) in seed_isolations {
        let at = interval.at(*k)?;
        if !isolation.real_rooted_in(&at) {
            return Err(format!("the seed {name}_{k} is not real-rooted in {at}"));
        }
    }
    let mut detail = format!(
        "exact: g and the seeds {name}_{first}..{name}_{} are real-rooted in {interval}",
        valid_from - 1
    );
    for step in steps {
        detail.push_str("; ");
        detail.push_str(&step);
    }
    Ok(detail)
}

/// A theorem report for a rule that was not matched.
pub(super) fn unmatched_theorem(reason: &str) -> RecurrenceRootTheorem {
    RecurrenceRootTheorem {
        theorem: "rolle_derivative_of_product".into(),
        statement: ROLLE_STATEMENT.into(),
        matched: false,
        reason: Some(reason.into()),
        form: None,
        factor: None,
        derivative_order: None,
        index_scalar: None,
        index_linear_factor: None,
        degree_formula: None,
        attempts: Vec::new(),
    }
}

/// Match `setup.rule` against Rolle's theorem and try the candidate
/// intervals, then the root hull of the factor and the seeds.
pub(super) fn rolle_core(setup: &RolleSetup<'_>, notes: &mut Vec<String>) -> RolleOutcome {
    let mut report = unmatched_theorem("");
    report.reason = None;
    let failed = |report| RolleOutcome {
        theorem: report,
        proven: None,
        full_degree: false,
    };
    let form = match match_rolle_form(setup.rule) {
        Ok(form) => form,
        Err(reason) => {
            report.reason = Some(reason);
            return failed(report);
        }
    };
    let name = setup.name.as_str();
    let (first, valid_from) = (setup.first, setup.valid_from);
    report.matched = true;
    report.form = Some(match &form.linear {
        None => format!(
            "{name}_m = s(m) D^{} (g {name}_(m-1)) for m >= {valid_from}, s(m) = {}, g = {}",
            form.order,
            poly_text(&form.scalar, "m"),
            poly_text(&form.factor, "x")
        ),
        Some(linear) => format!(
            "{name}_m = D^{} (g(m,x) {name}_(m-1)) for m >= {valid_from}, g(m,x) = ({}) ({} + \
             ({}) x)",
            form.order,
            poly_text(&form.factor, "x"),
            poly_text(&linear.constant, "m"),
            poly_text(&linear.slope, "m")
        ),
    });
    report.factor = Some(encode(&form.factor));
    report.derivative_order = Some(form.order);
    report.index_scalar = Some(encode(&form.scalar));
    report.index_linear_factor = form.linear.as_ref().map(|l| IndexLinearFactorJson {
        constant: encode(&l.constant),
        slope: encode(&l.slope),
    });

    if setup.seeds.is_empty() || setup.seeds.len() != valid_from - first {
        report.reason = Some("the seeds are not available".into());
        report.matched = false;
        return failed(report);
    }

    // Exact degree tracking (separable factors only).
    let factor_degree = form.factor.len() - 1;
    let (last_seed_index, last_seed) = setup.seeds.last().expect("at least one seed");
    let mut full_degree = false;
    if form.linear.is_none()
        && factor_degree >= form.order
        && scalar_nonvanishing(&form.scalar, valid_from)
        && !vec_poly::is_zero(last_seed)
    {
        let start = vec_poly::degree(last_seed).expect("nonzero seed");
        let slope = factor_degree - form.order;
        report.degree_formula = Some(format!(
            "deg {name}_n = {start} + {slope} (n - {last_seed_index}) for n >= {last_seed_index}"
        ));
        full_degree = setup.schedule.is_some_and(|schedule| {
            i64::try_from(slope).is_ok_and(|s| s == schedule.step)
                && setup.seeds.iter().all(|(k, seed)| {
                    schedule
                        .degree_at(*k)
                        .is_ok_and(|d| vec_poly::degree(seed) == Some(d))
                })
        });
    }

    let factor_isolation = isolate_real_roots_exact(&form.factor);
    let seed_isolations: Vec<(usize, RealRootIsolation)> = setup
        .seeds
        .iter()
        .map(|(k, seed)| (*k, isolate_real_roots_exact(seed)))
        .collect();
    let mut candidates = setup.candidates.clone();
    let hull = std::iter::once(&factor_isolation)
        .chain(seed_isolations.iter().map(|(_, iso)| iso))
        .filter_map(RealRootIsolation::enclosure_hull)
        .reduce(|a, b| a.hull(&b));
    if let Some(hull) = hull {
        candidates.push((RootInterval::from_real(&hull), "root_hull".into()));
    }
    if candidates.len() > setup.limit {
        notes.push(format!(
            "{} candidate intervals were not tried (max_candidate_intervals = {})",
            candidates.len() - setup.limit,
            setup.limit
        ));
        candidates.truncate(setup.limit);
    }
    let mut proven = None;
    for (interval, origin) in candidates {
        let outcome = check_attempt(
            &form,
            &factor_isolation,
            &seed_isolations,
            &interval,
            valid_from,
            first,
            name,
        );
        let success = outcome.is_ok();
        report.attempts.push(IntervalAttempt {
            interval: interval.clone(),
            origin,
            success,
            detail: outcome.unwrap_or_else(|reason| reason),
        });
        if success {
            proven = Some(interval);
            break;
        }
    }
    if proven.is_none() {
        notes.push("no candidate interval closes the Rolle induction".into());
    }
    RolleOutcome {
        theorem: report,
        proven,
        full_degree,
    }
}

#[cfg(test)]
mod tests;
