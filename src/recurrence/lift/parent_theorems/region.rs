//! Root-region helpers shared by the theorem checks of `parent_theorems`.
//!
//! - exact row tests: [`has_real_roots_below`] (every root `< point`),
//!   [`has_real_roots_in`] (every root in an open region),
//!   [`has_real_roots_in_region`] (closed or open ends, by Descartes' rule on
//!   the shifted row) and [`has_real_roots_in_closed_region`] (closed ends,
//!   by dividing out the roots at a closed end);
//! - regions of reports: [`open_root_interval`] (the requested region as
//!   recorded), [`proved_region`] (the index-independent region of a proved
//!   report) and [`requested_region_contains`] (whether a certified region
//!   satisfies an open request).
//!
//! Every test is exact; a zero row is in no region.

use super::{is_real_rooted_rational, RealRootednessTheoremReport, RootInterval, Q};
use crate::real_rootedness::RealInterval;
use crate::vec_poly;
use num_traits::{One, Signed, Zero};

/// Whether the nonzero `p` is real-rooted with every root `< point`: then
/// `p(point + z)` has all its coefficients nonzero and of one sign, and
/// conversely (Descartes' rule of signs).
pub(super) fn has_real_roots_below(p: &[Q], point: &Q) -> bool {
    let p = vec_poly::trim(p);
    if vec_poly::is_zero(&p) || !is_real_rooted_rational(&p) {
        return false;
    }
    let shifted = vec_poly::compose(&p, &[point.clone(), Q::one()]);
    shifted.iter().all(Signed::is_positive) || shifted.iter().all(Signed::is_negative)
}

/// Whether the nonzero real-rooted `p` has every root in the open
/// `region`.
pub(super) fn has_real_roots_in(p: &[Q], region: &RealInterval) -> bool {
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
pub(super) fn has_real_roots_in_region(p: &[Q], region: &RealInterval) -> bool {
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
pub(super) fn has_real_roots_in_closed_region(p: &[Q], region: &RealInterval) -> bool {
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

/// The interval with the same finite ends, both open.
pub(super) fn open_root_interval(interval: &RealInterval) -> RootInterval {
    RootInterval::from_real(&RealInterval {
        lower: interval.lower.clone(),
        lower_closed: false,
        upper: interval.upper.clone(),
        upper_closed: false,
    })
}

/// The index-independent root region of a proved report, if any.
pub(super) fn proved_region(report: &RealRootednessTheoremReport) -> Option<RealInterval> {
    let region = report.conclusion.as_ref()?.root_region.as_ref()?;
    if !region.is_index_independent() {
        return None;
    }
    region.at(0).ok()
}

/// Whether the requested (open) region contains `certified`.
pub(super) fn requested_region_contains(requested: &RealInterval, certified: &RealInterval) -> bool {
    let lower_ok = match (&requested.lower, &certified.lower) {
        (None, _) => true,
        (Some(_), None) => false,
        (Some(a), Some(l)) => a < l || (a == l && !certified.lower_closed),
    };
    let upper_ok = match (&requested.upper, &certified.upper) {
        (None, _) => true,
        (Some(_), None) => false,
        (Some(b), Some(u)) => b > u || (b == u && !certified.upper_closed),
    };
    lower_ok && upper_ok
}
