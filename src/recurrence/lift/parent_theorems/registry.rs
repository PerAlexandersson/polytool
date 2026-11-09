//! The registry of the theorem routes: [`ROUTES`] holds one [`TheoremRoute`]
//! per [`RealRootednessTheorem`], declared by the route's own module (its
//! `ROUTE` constant), in the order of [`RealRootednessTheorem::ALL`].
//!
//! The order is policy: it is the order of the theorem stages of the parent
//! cascade (`parent_proof`), of the attempts of [`discover_recurrence_proofs`]
//! and of the residue-class and index-section menus.  Every dispatch on a
//! theorem goes through [`route`].

use super::*;

/// A theorem check with the induction start at the supplied seeds (module
/// documentation of `parent_theorems`, "Moving the induction start").
pub(super) type TheoremCheck = fn(
    &Recurrence,
    usize,
    &[Vec<Q>],
    &RealRootednessTheoremOptions,
) -> Result<RealRootednessTheoremReport, RecurrenceError>;

/// Whether the check of a theorem may move its induction start (module
/// documentation of `parent_theorems`, "Moving the induction start").
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum StartPolicy {
    /// Unless the check at the supplied seeds gives a preferred proof or does
    /// not apply, retry with the start moved `1, ..., maximum_start_delay`
    /// rows later.
    MovesStart,
    /// The report at the supplied seeds is final: the check passes the
    /// start-delay bound to its sub-checks (the residue classes, the parent
    /// family of an index section), which move their own starts.
    OwnStarts,
}

/// One theorem route: its identity, how it is checked, and how its
/// induction start moves.
pub(in crate::recurrence) struct TheoremRoute {
    /// The theorem.
    pub(in crate::recurrence) theorem: RealRootednessTheorem,
    /// The check at the supplied seeds.
    pub(super) check_at_supplied_seeds: TheoremCheck,
    /// Whether the induction start may move.
    pub(super) start_policy: StartPolicy,
    /// Whether a proved report is the preferred proof; a proof that is not
    /// is kept only until a later start gives a preferred one.
    pub(super) is_preferred: fn(&RealRootednessTheoremReport) -> bool,
}

/// [`TheoremRoute::is_preferred`] of a route without weaker proofs.
pub(super) fn any_proof(_: &RealRootednessTheoremReport) -> bool {
    true
}

/// Every theorem route, in the order of [`RealRootednessTheorem::ALL`].
pub(in crate::recurrence) static ROUTES: [TheoremRoute; 12] = [
    product::ROUTE,
    three_term::ROUTE,
    liu_wang::ROUTE,
    derivative_lag::ROUTE,
    ma_wang::ROUTE,
    liu_wang_three_term::ROUTE,
    iterated_preserver::ROUTE,
    hypergeometric_row::ROUTE,
    variable_delay::ROUTE,
    residue_class::ROUTE,
    index_section::ROUTE,
    factor_gauge::ROUTE,
];

/// The route of `theorem`.
pub(in crate::recurrence) fn route(theorem: RealRootednessTheorem) -> &'static TheoremRoute {
    ROUTES
        .iter()
        .find(|route| route.theorem == theorem)
        .expect("every theorem has a route")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn theorems(routes: &[TheoremRoute]) -> Vec<RealRootednessTheorem> {
        routes.iter().map(|route| route.theorem).collect()
    }

    #[test]
    fn every_theorem_has_one_route_in_the_order_of_all() {
        assert_eq!(theorems(&ROUTES), RealRootednessTheorem::ALL);
        for theorem in RealRootednessTheorem::ALL {
            assert_eq!(route(theorem).theorem, theorem);
        }
    }

    #[test]
    fn the_discovery_lists_follow_the_route_order() {
        let listed: Vec<RealRootednessTheorem> = DISCOVERY_THEOREMS
            .iter()
            .chain(ROW_ROUTE_THEOREMS.iter())
            .copied()
            .collect();
        assert_eq!(listed, theorems(&ROUTES)[..listed.len()]);
    }
}
