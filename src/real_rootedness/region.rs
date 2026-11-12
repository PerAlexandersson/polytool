//! Exact root-region predicates for nonzero rational polynomials.
//!
//! The shifted-coefficient tests use Descartes' rule after exact
//! real-rootedness testing; finite endpoints can be open or closed.

use super::{is_real_rooted_bigint_coeffs, RealInterval, Q};
use crate::vec_poly;
use num_traits::{One, Signed, Zero};

fn is_real_rooted_rational(p: &[Q]) -> bool {
    is_real_rooted_bigint_coeffs(&vec_poly::clear_denominators(p))
}

/// Whether the nonzero `p` is real-rooted with every root `< point`: then
/// `p(point + z)` has all its coefficients nonzero and of one sign, and
/// conversely (Descartes' rule of signs).
pub(crate) fn has_real_roots_below(p: &[Q], point: &Q) -> bool {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) || !is_real_rooted_rational(&p) {
        return false;
    }
    let shifted = vec_poly::compose(&p, &[point.clone(), Q::one()]);
    shifted.iter().all(Signed::is_positive) || shifted.iter().all(Signed::is_negative)
}

/// Whether the nonzero real-rooted `p` has every root in the open
/// `region`.
pub(crate) fn has_real_roots_in(p: &[Q], region: &RealInterval) -> bool {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) || !is_real_rooted_rational(&p) {
        return false;
    }
    let upper_ok = region
        .upper
        .as_ref()
        .is_none_or(|upper| has_real_roots_below(&p, upper));
    let lower_ok = region.lower.as_ref().is_none_or(|lower| {
        let reflected: Vec<Q> = p
            .iter()
            .enumerate()
            .map(|(k, c)| if k % 2 == 1 { -c } else { c.clone() })
            .collect();
        has_real_roots_below(&reflected, &-lower)
    });
    upper_ok && lower_ok
}

/// Whether the nonzero `p` is real-rooted with every root in `region`
/// (closed or open ends).
///
/// Coefficients are ascending and trailing zeros are allowed. The zero
/// polynomial returns `false`; a nonzero constant returns `true`.
/// Unlike [`super::RealRootIsolation::real_rooted_in`], this predicate
/// excludes the zero polynomial and does not isolate roots.
pub fn polynomial_has_all_real_roots_in(p: &[Q], region: &RealInterval) -> bool {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) || !is_real_rooted_rational(&p) {
        return false;
    }
    let one_signed =
        |q: &[Q]| q.iter().all(|c| !c.is_negative()) || q.iter().all(|c| !c.is_positive());
    let upper_ok = region.upper.as_ref().is_none_or(|upper| {
        if region.upper_closed {
            one_signed(&vec_poly::compose(&p, &[upper.clone(), Q::one()]))
        } else {
            has_real_roots_below(&p, upper)
        }
    });
    let lower_ok = region.lower.as_ref().is_none_or(|lower| {
        let reflected = vec_poly::compose(&p, &[Q::zero(), -Q::one()]);
        if region.lower_closed {
            one_signed(&vec_poly::compose(&p, &[lower.clone(), -Q::one()]))
        } else {
            has_real_roots_below(&reflected, &-lower.clone())
        }
    });
    upper_ok && lower_ok
}

/// Whether the nonzero `p` is real-rooted with every root in `region`,
/// roots at a closed finite end allowed (they are divided out, then the open
/// region is tested).
pub(crate) fn has_real_roots_in_closed_region(p: &[Q], region: &RealInterval) -> bool {
    let mut q = vec_poly::trim(p);
    for (end, closed) in [
        (&region.lower, region.lower_closed),
        (&region.upper, region.upper_closed),
    ] {
        if let (Some(value), true) = (end, closed) {
            while q.len() > 1 {
                match vec_poly::exact_div_by_t_minus(&q, value) {
                    Some(quotient) => q = vec_poly::trim(&quotient),
                    None => break,
                }
            }
        }
    }
    !vec_poly::is_zero(&q) && (q.len() <= 1 || has_real_roots_in(&q, region))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::real_rootedness::isolate_real_roots_rational_coeffs;

    #[test]
    fn shifted_coefficients_agree_with_exact_isolation_and_endpoint_deflation() {
        let endpoints = [
            None,
            Some(Q::from_integer((-1).into())),
            Some(Q::zero()),
            Some(Q::one()),
        ];
        for code in 0..81 {
            let mut digits = code;
            let p: Vec<Q> = (0..4)
                .map(|_| {
                    let coefficient = Q::from_integer((digits % 3 - 1).into());
                    digits /= 3;
                    coefficient
                })
                .collect();
            let isolation = isolate_real_roots_rational_coeffs(&p);
            for lower in &endpoints {
                for upper in &endpoints {
                    for lower_closed in [false, true] {
                        for upper_closed in [false, true] {
                            let region = RealInterval {
                                lower: lower.clone(),
                                lower_closed,
                                upper: upper.clone(),
                                upper_closed,
                            };
                            if region.is_empty() {
                                continue;
                            }
                            let expected =
                                !vec_poly::is_zero(&p) && isolation.real_rooted_in(&region);
                            assert_eq!(
                                polynomial_has_all_real_roots_in(&p, &region),
                                expected,
                                "p={p:?}, region={region}"
                            );
                            assert_eq!(
                                has_real_roots_in_closed_region(&p, &region),
                                expected,
                                "deflation: p={p:?}, region={region}"
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn containment_keeps_open_ends_and_literal_empty_bounds_distinct() {
        let zero = Q::zero();
        let point = RealInterval::point(zero.clone());
        let open = RealInterval {
            upper_closed: false,
            ..RealInterval::at_most(zero.clone())
        };
        assert!(!open.contains_interval_bounds(&point));
        assert!(RealInterval::at_most(zero.clone()).contains_interval_bounds(&point));
        let empty = RealInterval {
            lower_closed: false,
            ..point.clone()
        };
        let outside = RealInterval::at_least(Q::one());
        assert!(empty.is_subset_of(&outside));
        assert!(!outside.contains_interval_bounds(&empty));
        assert!(polynomial_has_all_real_roots_in(
            &[zero.clone(), Q::one()],
            &point
        ));
        assert!(!has_real_roots_in(&[zero, Q::one()], &point));
    }
}
