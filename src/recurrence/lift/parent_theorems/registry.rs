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

/// When [`discover_recurrence_proofs`] records the attempt of a theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DiscoveryPolicy {
    /// Recorded for every record ([`DISCOVERY_THEOREMS`]).  With
    /// `first_version`, the theorem is one of the six of the first version
    /// of the discovery schema, which every saved report records.
    Always {
        /// One of the theorems of the first schema version.
        first_version: bool,
    },
    /// Recorded only when the check applies (its status is not
    /// `not_applicable`), so that the reports of other rules are unchanged.
    /// With `tolerated_when_missing`, a saved tail or alternative attempt
    /// list may omit it (it was saved before the theorem existed); without,
    /// it is one of [`ROW_ROUTE_THEOREMS`].
    WhenApplicable {
        /// A saved attempt list without the theorem still replays.
        tolerated_when_missing: bool,
    },
}

impl DiscoveryPolicy {
    /// Whether the attempt is recorded.
    pub(super) fn records(self, attempt: &RealRootednessTheoremReport) -> bool {
        match self {
            DiscoveryPolicy::Always { .. } => true,
            DiscoveryPolicy::WhenApplicable { .. } => {
                attempt.status != RealRootednessTheoremStatus::NotApplicable
            }
        }
    }

    /// Whether every saved report records the attempt.
    pub(super) fn in_first_version(self) -> bool {
        matches!(
            self,
            DiscoveryPolicy::Always {
                first_version: true
            }
        )
    }

    /// Whether a saved tail or alternative attempt list may omit the
    /// attempt although the current search records it.
    pub(super) fn tolerated_when_missing(self) -> bool {
        matches!(
            self,
            DiscoveryPolicy::WhenApplicable {
                tolerated_when_missing: true
            }
        )
    }
}

/// When the parent cascade (`prove_parent_real_rootedness`) runs a theorem
/// route; within a stage the routes run in the registry order and the first
/// proof ends the stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::recurrence) enum ParentStage {
    /// After the periodic product and before the paper-theorem and
    /// root-domain routes, while the parent is unproved.  A proof settles
    /// real-rootedness; a route of [`ParentStage::ReplacesEarly`] may still
    /// replace it.
    Early,
    /// After the root-domain route, while the parent is unproved or proved
    /// only by an [`ParentStage::Early`] route: a proof here replaces that
    /// proof (a theorem may conclude the interlacing of consecutive rows).
    ReplacesEarly,
    /// After the [`ParentStage::ReplacesEarly`] routes, only while the
    /// parent is unproved.
    WhileUnproved,
}

/// The root region a parent route requests from its theorem.  Whatever the
/// request, the route proves the parent statement only when the certified
/// region lies in `(-inf, 0]` (every row then has coefficients of one sign).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::recurrence) enum ParentRegion {
    /// The open ray `(-inf, 0)`.
    OpenNonpositiveRay,
    /// No region: the theorem certifies its own (for example a closed end
    /// at `0`, or an interval `[L, U]`).
    Unrestricted,
    /// The open ray `(-inf, 0)` and, when the theorem applies but is not
    /// proved with it, no region (sub-proofs may have roots at `0`).
    OpenNonpositiveRayThenUnrestricted,
}

/// How the parent cascade runs a theorem route.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::recurrence) struct ParentPolicy {
    /// When it runs.
    pub(in crate::recurrence) stage: ParentStage,
    /// The requested root region.
    pub(in crate::recurrence) region: ParentRegion,
    /// Whether the attempt is recorded only when the theorem applies (its
    /// status is not `not_applicable`).
    pub(in crate::recurrence) record_only_if_applicable: bool,
}

/// One theorem route: its identity, how it is checked, how its induction
/// start moves, when discovery records it, which menus offer it and how the
/// parent cascade runs it.
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
    /// The name of the route certificate ([`RouteCertificate`]) that a
    /// proved report carries, as named in replay errors (for example
    /// `"fixed-delay"`); `None` for a theorem without one.
    pub(super) certificate_name: Option<&'static str>,
    /// When [`discover_recurrence_proofs`] records the attempt.
    pub(super) discovery: DiscoveryPolicy,
    /// Whether [`check_residue_class_subsequence_theorem`] tries the theorem
    /// on a residue class.
    pub(super) in_residue_class_menu: bool,
    /// Whether [`check_index_section_lift_theorem`] tries the theorem on the
    /// parent family.
    pub(super) in_index_section_menu: bool,
    /// How the parent cascade runs the route.
    pub(in crate::recurrence) parent: ParentPolicy,
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

    fn filtered(predicate: impl Fn(&TheoremRoute) -> bool) -> Vec<RealRootednessTheorem> {
        ROUTES
            .iter()
            .filter(|route| predicate(route))
            .map(|route| route.theorem)
            .collect()
    }

    #[test]
    fn the_public_discovery_lists_are_the_discovery_policies() {
        assert_eq!(
            filtered(|route| matches!(route.discovery, DiscoveryPolicy::Always { .. })),
            DISCOVERY_THEOREMS
        );
        assert_eq!(
            filtered(|route| route.discovery
                == DiscoveryPolicy::WhenApplicable {
                    tolerated_when_missing: false
                }),
            ROW_ROUTE_THEOREMS
        );
        // The theorems of the first schema version: the first six.
        assert_eq!(
            filtered(|route| route.discovery.in_first_version()),
            DISCOVERY_THEOREMS[..6]
        );
        assert_eq!(
            filtered(|route| route.discovery.tolerated_when_missing()),
            [
                RealRootednessTheorem::IndexSectionLift,
                RealRootednessTheorem::FactorGauge
            ]
        );
    }

    #[test]
    fn the_menus_are_prefixes_of_the_route_order() {
        let all = theorems(&ROUTES);
        // A residue class is not split again; the parent of an index section
        // may itself be proved by residue classes.
        assert_eq!(filtered(|route| route.in_residue_class_menu), all[..9]);
        assert_eq!(filtered(|route| route.in_index_section_menu), all[..10]);
    }

    #[test]
    fn the_start_policies_and_certificates_are_those_of_the_routes() {
        assert_eq!(
            filtered(|route| route.start_policy == StartPolicy::OwnStarts),
            [
                RealRootednessTheorem::ResidueClassSubsequences,
                RealRootednessTheorem::IndexSectionLift
            ]
        );
        assert_eq!(
            filtered(|route| route.certificate_name.is_some()),
            [
                RealRootednessTheorem::IteratedPreserver,
                RealRootednessTheorem::VariableDelay,
                RealRootednessTheorem::ResidueClassSubsequences,
                RealRootednessTheorem::IndexSectionLift,
                RealRootednessTheorem::FactorGauge
            ]
        );
    }
}
