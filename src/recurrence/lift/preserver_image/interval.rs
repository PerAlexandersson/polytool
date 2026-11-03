//! Interval regions of the preserver-image lifts: parents proved real-rooted
//! with every root in a bounded interval `[L, U]` (rational ends), and the
//! image of such an interval under the operators of the menu whose theorem
//! gives one.
//!
//! # Interval parents
//!
//! [`find_bounded_root_region_report`] tries, in order:
//!
//! 1. the Ma--Wang interval check
//!    ([`crate::recurrence::check_ma_wang_interval_theorem`]) without a
//!    requested region; its interval `[L, U]` is read off the rule;
//! 2. the three-term theorem
//!    ([`crate::recurrence::check_three_term_recurrence_theorem`]) with a
//!    requested open region `(l, u)`, `l < u` integers: `u` is the smallest
//!    integer in `[-64, 64]` above every root of the first rows and `l` the
//!    largest integer in that range below them (exact Descartes tests on
//!    the shifted rows).  This only chooses the request; the region is
//!    certified by the theorem's positive-ratio induction for every row.
//!
//! Each report is replayed by
//! [`crate::recurrence::RealRootednessTheoremReport::verify`] when a
//! certification is verified.
//!
//! # Interval images
//!
//! For a parent rooted in `[L, U]` (each end open or closed) the input rows
//! `Q_n = lambda(n) T(P_(n+s))` are rooted in:
//!
//! - `P(a x + b)`, `a != 0`: `[(L - b)/a, (U - b)/a]` (ends exchanged for
//!   `a < 0`), since the roots are `(r - b)/a`;
//! - `D^j o x^m`, `x^m D^j`: the convex hull of `[L, U]` and `0` (only
//!   `[L, U]` for `m = 0`), since `x^m P` has the roots of `P` and `0`, and
//!   every derivative of a real-rooted polynomial is zero or real-rooted
//!   with its roots between the smallest and the largest root (Rolle); the
//!   end `0` is closed;
//! - `x^(d_n) P(1/x)`: for `L > 0` the roots `1/r` lie in `[1/U, 1/L]` and the
//!   extra roots `0` (degree defect) give `[0, 1/L]`; for `U < 0` likewise
//!   `[1/U, 0]`; for `L = 0` the ray `[0, +inf)`, for `U = 0` the ray
//!   `(-inf, 0]`, otherwise the line;
//! - `(c x + e)^(d_n) P((a x + b)/(c x + e))`: a root `r` gives
//!   `phi(r) = (e r - b)/(a - c r)`, which is monotone on `[L, U]` when the
//!   pole `a/c` lies outside it, so the roots lie between `phi(L)` and
//!   `phi(U)`; the extra roots `-e/c` widen this to the hull with `-e/c`
//!   (closed).  With the pole in `[L, U]`, the line; for `c = 0` the affine
//!   map `x -> (a x + b)/e`;
//! - `f(n, theta)`, `phi(n, D)`, hypergeometric multipliers and the Euler
//!   step: only through
//!   the ray theorems, when `[L, U]` lies in `(-inf, 0]` (no interval theorem
//!   is used for them);
//! - a general Weyl operator: none.
//!
//! The substitution maps are injective on roots, and the side conditions of
//! the operators (none for these maps) are unchanged; so a proved interval
//! parent gives the image region for every input row.

use super::operator::PreserverOperator;
use crate::real_rootedness::RealInterval;
use crate::recurrence::lift::lift_summary::real_line;
use crate::recurrence::{
    check_ma_wang_interval_theorem, check_three_term_recurrence_theorem, format_rational_coeff,
    parse_rational_coeff, RealRootednessTheoremOptions, RealRootednessTheoremReport, Recurrence,
    RecurrenceError, RootEndpoint, RootInterval,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

type Q = BigRational;

/// Rows used to choose the requested region of the three-term theorem.
const CANDIDATE_ROWS: usize = 12;

/// Largest absolute value of an integer end tried for the three-term region.
const CANDIDATE_BOUND: i64 = 64;

/// The integer polynomial with the same roots.
fn integer_coefficients(p: &[Q]) -> Vec<BigInt> {
    let common = p
        .iter()
        .fold(BigInt::one(), |accumulator, c| accumulator.lcm(c.denom()));
    p.iter()
        .map(|c| c.numer() * (&common / c.denom()))
        .collect()
}

/// Whether every coefficient of `p(point + sign z)` is nonzero and of one
/// sign (for a real-rooted `p`: every root `< point` for `sign = 1`, every
/// root `> point` for `sign = -1`).
fn shifted_one_sign(p: &[Q], point: &Q, sign: i64) -> bool {
    let shifted = vec_poly::compose(p, &[point.clone(), Q::from_integer(BigInt::from(sign))]);
    shifted.iter().all(Signed::is_positive) || shifted.iter().all(Signed::is_negative)
}

/// The integer interval `(l, u)` chosen for the three-term request
/// (module documentation), or `None` when a row is not real-rooted or no
/// integer in range bounds the roots.
fn three_term_candidate(rows: &[Vec<Q>]) -> Option<(Q, Q)> {
    let rows: Vec<Vec<Q>> = rows
        .iter()
        .map(|row| vec_poly::trim(row))
        .filter(|row| vec_poly::degree(row).is_some_and(|d| d > 0))
        .collect();
    if rows.is_empty()
        || rows
            .iter()
            .any(|row| !crate::is_real_rooted_bigint_coeffs(&integer_coefficients(row)))
    {
        return None;
    }
    let q = |v: i64| Q::from_integer(BigInt::from(v));
    let upper = (-CANDIDATE_BOUND..=CANDIDATE_BOUND)
        .find(|&u| rows.iter().all(|row| shifted_one_sign(row, &q(u), 1)))?;
    let lower = (-CANDIDATE_BOUND..=CANDIDATE_BOUND)
        .rev()
        .find(|&l| rows.iter().all(|row| shifted_one_sign(row, &q(l), -1)))?;
    (lower < upper).then(|| (q(lower), q(upper)))
}

/// Whether a proved report certifies a bounded root region.
fn bounded_region(report: &RealRootednessTheoremReport) -> bool {
    report.is_proved()
        && report
            .conclusion
            .as_ref()
            .and_then(|c| c.root_region.as_ref())
            .is_some_and(|r| r.lower.is_some() && r.upper.is_some() && r.is_index_independent())
}

/// A replayable theorem report proving every row `P_n`, `n >= first_index`,
/// of the family `rule` with identifying `seeds` real-rooted with every root
/// in a bounded interval `[L, U]` with rational ends (module documentation,
/// "Interval parents"), or `None` when neither the Ma--Wang interval check
/// nor the three-term theorem with a chosen bounded region proves one.  The
/// report replays with [`RealRootednessTheoremReport::verify`].
///
/// # Errors
///
/// Errors of the theorem checks (malformed records).
pub fn find_bounded_root_region_report(
    rule: &Recurrence,
    first_index: usize,
    seeds: &[Vec<Q>],
) -> Result<Option<RealRootednessTheoremReport>, RecurrenceError> {
    let ma_wang = check_ma_wang_interval_theorem(
        rule,
        first_index,
        seeds,
        &RealRootednessTheoremOptions::default(),
    )?;
    if bounded_region(&ma_wang) {
        return Ok(Some(ma_wang));
    }
    let whole_line = check_three_term_recurrence_theorem(
        rule,
        first_index,
        seeds,
        &RealRootednessTheoremOptions::default(),
    )?;
    if !whole_line.is_proved() {
        return Ok(None);
    }
    let Ok(rows) = rule.generate_rows_rational(seeds, first_index, CANDIDATE_ROWS) else {
        return Ok(None);
    };
    let Some((lower, upper)) = three_term_candidate(&rows) else {
        return Ok(None);
    };
    let report = check_three_term_recurrence_theorem(
        rule,
        first_index,
        seeds,
        &RealRootednessTheoremOptions {
            root_region: Some(RealInterval {
                lower: Some(lower),
                lower_closed: false,
                upper: Some(upper),
                upper_closed: false,
            }),
            ..RealRootednessTheoremOptions::default()
        },
    )?;
    Ok(bounded_region(&report).then_some(report))
}

/// Whether a certification computes the interval parent report: for the
/// operators whose theorem gives a region from an interval (substitutions
/// and Rolle operators), and for the others (not a general Weyl operator)
/// when the parent statement with roots in `(-inf, 0]` is not proved.
pub(super) fn wants_interval_parent(operator: &PreserverOperator, parent_proved: bool) -> bool {
    match operator {
        PreserverOperator::Affine { .. }
        | PreserverOperator::Reversal { .. }
        | PreserverOperator::Mobius { .. }
        | PreserverOperator::DerivativeOfPowerMultiple { .. }
        | PreserverOperator::PowerTimesDerivative { .. } => true,
        PreserverOperator::EulerPolynomial { .. }
        | PreserverOperator::DerivativePolynomial { .. }
        | PreserverOperator::HadamardMultiplier { .. }
        | PreserverOperator::EulerStep { .. } => !parent_proved,
        PreserverOperator::WeylOperator { .. } => false,
    }
}

/// The route label of an interval parent report.
pub(super) fn interval_route_label(report: &RealRootednessTheoremReport) -> &'static str {
    match report.theorem {
        crate::recurrence::RealRootednessTheorem::MaWangInterval => "ma_wang_interval",
        _ => "three_term_interval",
    }
}

fn endpoint(value: &Q, closed: bool) -> RootEndpoint {
    RootEndpoint {
        value: format_rational_coeff(value),
        per_index: None,
        closed,
    }
}

/// The interval with the given ends (`None`: infinite).
fn interval(lower: Option<(Q, bool)>, upper: Option<(Q, bool)>) -> RootInterval {
    RootInterval {
        lower: lower.map(|(v, closed)| endpoint(&v, closed)),
        upper: upper.map(|(v, closed)| endpoint(&v, closed)),
    }
}

/// The smallest interval containing `[lower, upper]` and the point `p`
/// (`p` closed).
fn hull_with_point(lower: (Q, bool), upper: (Q, bool), p: &Q) -> RootInterval {
    let lower = if *p <= lower.0 {
        (p.clone(), true)
    } else {
        lower
    };
    let upper = if *p >= upper.0 {
        (p.clone(), true)
    } else {
        upper
    };
    interval(Some(lower), Some(upper))
}

fn parse(text: &str) -> Option<Q> {
    parse_rational_coeff(text).ok()
}

/// The roots of `P(slope x + intercept)` for `P` rooted in the bounded
/// interval `[lower, upper]` (each end with its closedness): a root `r` of
/// `P` gives the root `(r - intercept)/slope`.
fn affine_image(lower: &(Q, bool), upper: &(Q, bool), slope: &Q, intercept: &Q) -> RootInterval {
    let map = |end: &(Q, bool)| ((&end.0 - intercept) / slope, end.1);
    if slope.is_positive() {
        interval(Some(map(lower)), Some(map(upper)))
    } else {
        interval(Some(map(upper)), Some(map(lower)))
    }
}

/// The root region of the input rows when every used parent row is
/// real-rooted with every root in the bounded `parent` interval (module
/// documentation, "Interval images"); `None` when the operator's theorem
/// gives no region from an interval.  `from_nonpositive` is the region the
/// operator's ray theorem gives from `(-inf, 0]`.
pub(super) fn interval_image(
    operator: &PreserverOperator,
    parent: &RealInterval,
    from_nonpositive: Option<&RootInterval>,
) -> Option<RootInterval> {
    let lower = (parent.lower.clone()?, parent.lower_closed);
    let upper = (parent.upper.clone()?, parent.upper_closed);
    let zero = Q::zero();
    match operator {
        PreserverOperator::Affine { slope, intercept } => {
            let (a, b) = (parse(slope)?, parse(intercept)?);
            (!a.is_zero()).then(|| affine_image(&lower, &upper, &a, &b))
        }
        PreserverOperator::DerivativeOfPowerMultiple { power, .. }
        | PreserverOperator::PowerTimesDerivative { power, .. } => Some(if *power == 0 {
            interval(Some(lower), Some(upper))
        } else {
            hull_with_point(lower, upper, &zero)
        }),
        PreserverOperator::Reversal { .. } => Some(if lower.0.is_positive() {
            interval(Some((zero, true)), Some((lower.0.recip(), lower.1)))
        } else if upper.0.is_negative() {
            interval(Some((upper.0.recip(), upper.1)), Some((zero, true)))
        } else if lower.0.is_zero() {
            interval(Some((zero, true)), None)
        } else if upper.0.is_zero() {
            interval(None, Some((zero, true)))
        } else {
            real_line()
        }),
        PreserverOperator::Mobius { a, b, c, e, .. } => {
            let (a, b, c, e) = (parse(a)?, parse(b)?, parse(c)?, parse(e)?);
            if (&a * &e - &b * &c).is_zero() {
                return None;
            }
            if c.is_zero() {
                return Some(affine_image(&lower, &upper, &(&a / &e), &(&b / &e)));
            }
            let pole = &a / &c;
            if lower.0 <= pole && pole <= upper.0 {
                return Some(real_line());
            }
            let phi = |r: &Q| (&e * r - &b) / (&a - &c * r);
            let (at_lower, at_upper) = ((phi(&lower.0), lower.1), (phi(&upper.0), upper.1));
            let (low, high) = if at_lower.0 <= at_upper.0 {
                (at_lower, at_upper)
            } else {
                (at_upper, at_lower)
            };
            Some(hull_with_point(low, high, &-(&e / &c)))
        }
        PreserverOperator::EulerPolynomial { .. }
        | PreserverOperator::DerivativePolynomial { .. }
        | PreserverOperator::HadamardMultiplier { .. }
        | PreserverOperator::EulerStep { .. } => (!upper.0.is_positive())
            .then(|| from_nonpositive.cloned())
            .flatten(),
        PreserverOperator::WeylOperator { .. } => None,
    }
}

/// The concrete ends of an index-independent root interval.
fn concrete(region: &RootInterval) -> Option<RealInterval> {
    region
        .is_index_independent()
        .then(|| region.at(0).ok())
        .flatten()
}

/// Whether `inner` is contained in `outer` (both index-independent; `false`
/// otherwise).
pub(super) fn region_contains(outer: &RootInterval, inner: &RootInterval) -> bool {
    let (Some(outer), Some(inner)) = (concrete(outer), concrete(inner)) else {
        return false;
    };
    let lower_ok = match (&outer.lower, &inner.lower) {
        (None, _) => true,
        (Some(_), None) => false,
        (Some(o), Some(i)) => o < i || (o == i && (outer.lower_closed || !inner.lower_closed)),
    };
    let upper_ok = match (&outer.upper, &inner.upper) {
        (None, _) => true,
        (Some(_), None) => false,
        (Some(o), Some(i)) => o > i || (o == i && (outer.upper_closed || !inner.upper_closed)),
    };
    lower_ok && upper_ok
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;
    use crate::recurrence::DegreeSchedule;

    fn closed(lower: i64, upper: i64) -> RealInterval {
        RealInterval::closed(q(lower), q(upper))
    }

    fn schedule() -> DegreeSchedule {
        DegreeSchedule {
            first_index: 0,
            first_degree: 0,
            step: 1,
            last_index: None,
        }
    }

    #[test]
    fn substitution_images_of_an_interval() {
        let parent = closed(-1, 0);
        let image = |operator: PreserverOperator| {
            interval_image(&operator, &parent, None).map(|r| r.to_string())
        };
        // P(2x + 1): roots (r - 1)/2.
        assert_eq!(
            image(PreserverOperator::Affine {
                slope: "2".into(),
                intercept: "1".into()
            }),
            Some("[-1, -1/2]".into())
        );
        // P(-x): roots -r.
        assert_eq!(
            image(PreserverOperator::Affine {
                slope: "-1".into(),
                intercept: "0".into()
            }),
            Some("[0, 1]".into())
        );
        // The f-polynomial (x + 1)^d P(x/(x + 1)): r/(1 - r) and -1.
        assert_eq!(
            image(PreserverOperator::Mobius {
                a: "1".into(),
                b: "0".into(),
                c: "1".into(),
                e: "1".into(),
                schedule: schedule(),
            }),
            Some("[-1, 0]".into())
        );
        // Reversal of [-1, 0] (end 0 attained): (-inf, 0].
        assert_eq!(
            image(PreserverOperator::Reversal {
                schedule: schedule()
            }),
            Some("(-inf, 0]".into())
        );
        // Reversal of (1, 2]: [0, 1).
        let open_lower = RealInterval {
            lower: Some(q(1)),
            lower_closed: false,
            upper: Some(q(2)),
            upper_closed: true,
        };
        assert_eq!(
            interval_image(
                &PreserverOperator::Reversal {
                    schedule: schedule()
                },
                &open_lower,
                None
            )
            .map(|r| r.to_string()),
            Some("[0, 1)".into())
        );
        // Reversal of [-2, -1]: the roots 1/r lie in [-1, -1/2], with 0.
        assert_eq!(
            interval_image(
                &PreserverOperator::Reversal {
                    schedule: schedule()
                },
                &closed(-2, -1),
                None
            )
            .map(|r| r.to_string()),
            Some("[-1, 0]".into())
        );
        // Rolle with x^1 on [1, 2]: the hull with 0.
        assert_eq!(
            interval_image(
                &PreserverOperator::DerivativeOfPowerMultiple {
                    power: 1,
                    derivative_order: 1
                },
                &closed(1, 2),
                None
            )
            .map(|r| r.to_string()),
            Some("[0, 2]".into())
        );
        // A multiplier on an interval with a positive end: no interval theorem.
        assert_eq!(
            interval_image(
                &PreserverOperator::EulerPolynomial {
                    coefficients: Vec::new()
                },
                &closed(0, 1),
                Some(&real_line())
            ),
            None
        );
    }

    #[test]
    fn containment_respects_open_ends() {
        let region = |lower: Option<(i64, bool)>, upper: Option<(i64, bool)>| {
            interval(lower.map(|(v, c)| (q(v), c)), upper.map(|(v, c)| (q(v), c)))
        };
        let ray = region(None, Some((0, true)));
        assert!(region_contains(
            &ray,
            &region(Some((-1, true)), Some((0, true)))
        ));
        assert!(!region_contains(
            &region(None, Some((0, false))),
            &region(Some((-1, true)), Some((0, true)))
        ));
        assert!(region_contains(
            &region(None, Some((0, false))),
            &region(Some((-1, true)), Some((0, false)))
        ));
        assert!(!region_contains(
            &region(Some((-1, true)), Some((0, true))),
            &ray
        ));
    }

    #[test]
    fn three_term_candidates_bound_the_roots() {
        // Chebyshev U_2 = 4x^2 - 1 and U_3 = 8x^3 - 4x: roots in (-1, 1).
        let rows = vec![vec![q(-1), q(0), q(4)], vec![q(0), q(-4), q(0), q(8)]];
        assert_eq!(three_term_candidate(&rows), Some((q(-1), q(1))));
        // x^2 + 1 is not real-rooted.
        assert_eq!(three_term_candidate(&[vec![q(1), q(0), q(1)]]), None);
    }
}
