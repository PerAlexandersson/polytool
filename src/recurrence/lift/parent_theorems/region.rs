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

use super::{RealRootednessTheoremReport, RootInterval};
use crate::real_rootedness::RealInterval;
pub(super) use crate::real_rootedness::{
    has_real_roots_below, has_real_roots_in, has_real_roots_in_closed_region,
    polynomial_has_all_real_roots_in as has_real_roots_in_region,
};

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
pub(super) fn requested_region_contains(
    requested: &RealInterval,
    certified: &RealInterval,
) -> bool {
    let open = RealInterval {
        lower_closed: false,
        upper_closed: false,
        ..requested.clone()
    };
    open.contains_interval_bounds(certified)
}
