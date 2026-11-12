//! The registry of the theorem routes: [`ROUTES`] holds one [`TheoremRoute`]
//! per [`RealRootednessTheorem`]. The shared `theorem_routes!` table supplies
//! identity, names, summaries and exports; each module's `ROUTE` supplies
//! its implementation, base reference and policies.
//!
//! The order is policy: it is the order of the theorem stages of the parent
//! cascade (`parent_proof`), of the attempts of [`discover_recurrence_proofs`]
//! and of the residue-class and index-section menus.  Every dispatch on a
//! theorem goes through [`route`] or iterates [`ROUTES`].  Adding a route:
//! one module with its `ROUTE` policy, one table row, and one typed
//! certificate row when needed (see `doc/development.md`).

use super::*;

// One row enrolls a route in both public enums and their mappings, the
// ordered registry, its module declaration and its public exports.
macro_rules! theorem_routes {
    ($consumer:ident) => {
        $consumer! {
            /// Product induction with uniformly real-rooted multipliers.
            ProductRecurrence, product, "product_recurrence",
            "Product induction with uniformly real-rooted multipliers.", {
                check_product_recurrence_theorem,
                is_trivially_real_rooted,
            };
            /// [`check_three_term_recurrence_theorem`].
            ThreeTermRecurrence, three_term, "three_term_recurrence",
            "Three-term recurrence with strict interlacing and a certified root region.", {
                check_three_term_recurrence_theorem,
            };
            /// [`check_liu_wang_theorem`].
            LiuWang, liu_wang, "liu_wang",
            "Liu--Wang first-order derivative recurrence.", {
                check_liu_wang_theorem,
            };
            /// [`check_liu_wang_derivative_lag_theorem`].
            LiuWangDerivativeLag, derivative_lag, "liu_wang_derivative_lag",
            "Liu--Wang derivative-lag recurrence.", {
                check_liu_wang_derivative_lag_theorem,
            };
            /// [`check_ma_wang_interval_theorem`].
            MaWangInterval, ma_wang, "ma_wang_interval",
            "Ma--Wang derivative recurrence on a certified interval.", {
                check_ma_wang_interval_theorem,
            };
            /// [`check_liu_wang_three_term_theorem`].
            LiuWangThreeTerm, liu_wang_three_term, "liu_wang_three_term",
            "Liu--Wang three-term recurrence with polynomial coefficients.", {
                check_liu_wang_three_term_theorem,
            };
            /// [`check_iterated_preserver_theorem`].
            IteratedPreserver, iterated_preserver, "iterated_preserver",
            "Induction by iterated real-rootedness preservers.", {
                check_iterated_preserver_theorem,
                check_weyl_step_preserver,
                IteratedPreserverCertificate,
                IteratedPreserverCoordinates,
                IteratedPreserverFactor,
                WeylStepPreserverVerdict,
            };
            /// [`check_hypergeometric_row_theorem`].
            HypergeometricRow, hypergeometric_row, "hypergeometric_row",
            "Exact row-hypergeometric real-rootedness criterion.", {
                check_hypergeometric_row_theorem,
            };
            /// [`check_variable_delay_theorem`] (RealRooted route `variableDelay`).
            VariableDelay, variable_delay, "variable_delay",
            "Fixed-delay recurrence with a certified interlacing window.", {
                check_variable_delay_theorem,
                VariableDelayCertificate,
                VariableDelaySeedStep,
            };
            /// [`check_residue_class_subsequence_theorem`] (RealRooted route
            /// `subseq`).
            ResidueClassSubsequences, residue_class, "residue_class_subsequences",
            "Proof by independent residue-class subsequences.", {
                check_residue_class_subsequence_theorem,
                ResidueClassCertificate,
                ResidueClassProof,
            };
            /// [`check_index_section_lift_theorem`] (RealRooted route `pathPower`).
            IndexSectionLift, index_section, "index_section_lift",
            "Proof by an index section of a certified parent family.", {
                check_index_section_lift_theorem,
                IndexSectionCertificate,
                IndexSectionTransform,
                MAX_INDEX_SECTION_PERIOD,
            };
            /// [`check_factor_gauge_theorem`] (RealRooted routes `powerForm`,
            /// `closedForm`, `factorPower`).
            FactorGauge, factor_gauge, "factor_gauge",
            "Factor gauge with a bounded residual family.", {
                check_factor_gauge_theorem,
                FactorGaugeCertificate,
                FactorGaugeDominantRoot,
                FactorGaugeResidualShape,
            };
            /// [`check_lag_two_euler_theorem`] (RealRooted positive lag-two Euler recurrence).
            LagTwoEuler, lag_two_euler, "lag_two_euler",
            "Positive lag-two Euler recurrence with strict interlacing.", {
                check_lag_two_euler_theorem,
                LagTwoEulerCertificate,
                LagTwoEulerSeedStep,
            };
            /// [`check_binomial_perturbed_euler_orbit_theorem`] (RealRooted issue #2513).
            BinomialPerturbedEulerOrbit, binomial_perturbed_euler_orbit, "binomial_perturbed_euler_orbit",
            "Binomial perturbation of an Euler-operator orbit.", {
                check_binomial_perturbed_euler_orbit_theorem,
                BinomialPerturbedEulerOrbitCertificate,
            };
            /// [`check_rational_pf_riordan_theorem`].
            RationalPfRiordan, rational_pf_riordan, "rational_pf_riordan",
            "Exact rational PF Riordan certificate for all rows.", {
                check_rational_pf_riordan_theorem,
                check_rational_polya_frequency,
                PolyaFrequencyRootLayer,
                RationalPfRiordanCertificate,
                RationalPolyaFrequencyEvidence,
            };
            /// [`check_degree_law_reversal_theorem`], a shared closure of the theorem menu.
            DegreeLawReversal, degree_law_reversal, "degree_law_reversal",
            "Degree-law reversal followed by a theorem on the transformed family.", {
                check_degree_law_reversal_theorem,
                DegreeLawReversalCertificate,
                DegreeLawTransform,
                MAXIMUM_REVERSAL_DEGREE_SLOPE,
            };
        }
    };
}
pub(in crate::recurrence) use theorem_routes;

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
    /// The report at the supplied seeds is final: either it proves an
    /// all-index identity unaffected by extending the seed prefix, or it
    /// passes the start-delay bound to sub-checks that move their own starts.
    OwnStarts,
}

/// When [`discover_recurrence_proofs`] records the attempt of a theorem.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DiscoveryPolicy {
    /// Recorded for every record.
    Always,
    /// Recorded only when the check applies (its status is not
    /// `not_applicable`), so that the reports of other rules are unchanged.
    WhenApplicable,
}

impl DiscoveryPolicy {
    /// Whether the attempt is recorded.
    pub(super) fn records(self, attempt: &RealRootednessTheoremReport) -> bool {
        match self {
            DiscoveryPolicy::Always => true,
            DiscoveryPolicy::WhenApplicable => {
                attempt.status != RealRootednessTheoremStatus::NotApplicable
            }
        }
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

/// When discovery runs a theorem check.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DiscoveryPhase {
    /// On the supplied family, its tails and direct representations.
    Direct,
    /// After all direct routes, lifts and affine coordinates fail.
    AfterDirectLifts,
}

/// Membership in the existing discovery replay protocol.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum DiscoveryReplayPolicy {
    /// Required in exactly the current direct-attempt membership and order.
    Direct,
    /// Replayed separately when recorded as the last attempt. The message
    /// preserves the route's rejection after an earlier direct proof or lift.
    TrailingClosure { premature_error: &'static str },
}

impl DiscoveryReplayPolicy {
    /// The rejection message of a trailing closure, absent for direct routes.
    pub(super) fn premature_closure_error(self) -> Option<&'static str> {
        match self {
            Self::Direct => None,
            Self::TrailingClosure { premature_error } => Some(premature_error),
        }
    }
}

/// Replay membership of the parent theorem routes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::recurrence) enum ParentReplayPolicy {
    /// Replay the recorded route alone, even if another route now wins.
    RecordedRoute,
}

/// A late route's eligibility check against records already tried directly.
type AfterDirectEligibility =
    fn(&Recurrence, usize, &[Vec<Q>], &[&RecurrenceJson]) -> Result<bool, RecurrenceError>;

/// One theorem route: its identity, how it is checked, how its induction
/// start moves, when discovery records it, which menus offer it and how the
/// parent cascade runs it.
pub(in crate::recurrence) struct TheoremRoute {
    /// The identity generated by the shared route table.
    pub(in crate::recurrence) theorem: RealRootednessTheorem,
    /// The route-owned implementation and policies.
    pub(in crate::recurrence) policy: RoutePolicy,
    /// The declarations provided by the route's matched report, including
    /// branch-specific and supporting declarations of composed checks.
    pub(super) lean_declarations: fn(&RealRootednessTheoremReport) -> &[String],
}

/// The route implementation and policies, supplied by its owning module.
pub(in crate::recurrence) struct RoutePolicy {
    /// The base reference used to initialize a theorem report.
    pub(super) reference: &'static str,
    /// The discovery stage that runs this route.
    pub(super) phase: DiscoveryPhase,
    /// Optional eligibility check after direct lifts (including duplicate detection).
    pub(super) after_direct_eligible: Option<AfterDirectEligibility>,
    /// How a recorded parent proof is replayed.
    pub(in crate::recurrence) replay: ParentReplayPolicy,
    /// How discovery includes this route in replay membership.
    pub(super) discovery_replay: DiscoveryReplayPolicy,
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
    /// Whether degree-law reversal may use the theorem on its transformed family.
    pub(super) in_degree_law_reversal_menu: bool,
    /// How the parent cascade runs the route.
    pub(in crate::recurrence) parent: ParentPolicy,
}

/// [`RoutePolicy::is_preferred`] of a route without weaker proofs.
pub(super) fn any_proof(_: &RealRootednessTheoremReport) -> bool {
    true
}

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

    #[test]
    fn late_discovery_and_recorded_route_replay_membership_are_explicit() {
        assert_eq!(
            filtered(|route| route.policy.phase == DiscoveryPhase::AfterDirectLifts),
            [RealRootednessTheorem::DegreeLawReversal]
        );
        assert_eq!(
            filtered(|route| route.policy.after_direct_eligible.is_some()),
            [RealRootednessTheorem::DegreeLawReversal]
        );
        assert_eq!(
            filtered(|route| route
                .policy
                .discovery_replay
                .premature_closure_error()
                .is_some()),
            [RealRootednessTheorem::DegreeLawReversal]
        );
        for route in &ROUTES {
            assert_eq!(route.policy.replay, ParentReplayPolicy::RecordedRoute);
        }
    }

    #[test]
    fn parent_region_recording_stage_and_menu_policies_are_frozen() {
        use DiscoveryPolicy::{Always, WhenApplicable};
        use ParentRegion::{
            OpenNonpositiveRay as Ray, OpenNonpositiveRayThenUnrestricted as Fallback, Unrestricted,
        };
        use ParentStage::{Early, ReplacesEarly, WhileUnproved};
        use RealRootednessTheorem::*;
        // Independent expectations: never derive these from ALL or table slices.
        let expected = [
            (ProductRecurrence, Early, Ray, false, Always, true, true),
            (
                ThreeTermRecurrence,
                ReplacesEarly,
                Ray,
                false,
                Always,
                true,
                true,
            ),
            (LiuWang, ReplacesEarly, Ray, false, Always, true, true),
            (
                LiuWangDerivativeLag,
                ReplacesEarly,
                Unrestricted,
                false,
                Always,
                true,
                true,
            ),
            (
                MaWangInterval,
                ReplacesEarly,
                Unrestricted,
                false,
                Always,
                true,
                true,
            ),
            (
                LiuWangThreeTerm,
                ReplacesEarly,
                Unrestricted,
                false,
                Always,
                true,
                true,
            ),
            (
                IteratedPreserver,
                ReplacesEarly,
                Unrestricted,
                true,
                Always,
                true,
                true,
            ),
            (
                HypergeometricRow,
                ReplacesEarly,
                Unrestricted,
                true,
                Always,
                true,
                true,
            ),
            (
                VariableDelay,
                WhileUnproved,
                Ray,
                true,
                WhenApplicable,
                true,
                true,
            ),
            (
                ResidueClassSubsequences,
                WhileUnproved,
                Fallback,
                true,
                WhenApplicable,
                false,
                true,
            ),
            (
                IndexSectionLift,
                WhileUnproved,
                Fallback,
                true,
                WhenApplicable,
                false,
                false,
            ),
            (
                FactorGauge,
                WhileUnproved,
                Fallback,
                true,
                WhenApplicable,
                false,
                false,
            ),
            (
                LagTwoEuler,
                WhileUnproved,
                Ray,
                true,
                WhenApplicable,
                true,
                true,
            ),
            (
                BinomialPerturbedEulerOrbit,
                WhileUnproved,
                Unrestricted,
                true,
                WhenApplicable,
                true,
                true,
            ),
            (
                RationalPfRiordan,
                WhileUnproved,
                Unrestricted,
                true,
                WhenApplicable,
                true,
                true,
            ),
            (
                DegreeLawReversal,
                WhileUnproved,
                Unrestricted,
                true,
                WhenApplicable,
                false,
                false,
            ),
        ];
        let actual: Vec<_> = ROUTES
            .iter()
            .map(|r| {
                (
                    r.theorem,
                    r.policy.parent.stage,
                    r.policy.parent.region,
                    r.policy.parent.record_only_if_applicable,
                    r.policy.discovery,
                    r.policy.in_residue_class_menu,
                    r.policy.in_index_section_menu,
                )
            })
            .collect();
        assert_eq!(actual, expected);
    }

    #[test]
    fn recording_distinguishes_not_applicable_from_all_other_outcomes() {
        let entry = crate::oeis::by_id("A008292").unwrap();
        let (rule, rows) = entry.recurrence_parts().unwrap();
        let mut report = check_liu_wang_theorem(
            &rule,
            entry.recurrence_first_index,
            &rows[..1],
            &Default::default(),
        )
        .unwrap();
        for status in [
            RealRootednessTheoremStatus::Proved,
            RealRootednessTheoremStatus::HypothesisFailed,
            RealRootednessTheoremStatus::Unresolved,
            RealRootednessTheoremStatus::NotApplicable,
        ] {
            report.status = status;
            assert!(DiscoveryPolicy::Always.records(&report));
            assert_eq!(
                DiscoveryPolicy::WhenApplicable.records(&report),
                status != RealRootednessTheoremStatus::NotApplicable
            );
        }
    }

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
    fn the_discovery_policies_record_the_first_eight_always() {
        let all = theorems(&ROUTES);
        assert_eq!(
            filtered(|route| route.policy.discovery == DiscoveryPolicy::Always),
            all[..8]
        );
    }

    #[test]
    fn every_theorem_has_its_parent_route_and_the_stages_follow_the_order() {
        for route in &ROUTES {
            let parent_route = VeroneseParentProofRoute::of_theorem(route.theorem);
            assert_eq!(parent_route.theorem(), Some(route.theorem));
            assert_eq!(parent_route.as_str(), route.theorem.as_str());
        }
        let stages: Vec<ParentStage> = ROUTES
            .iter()
            .map(|route| route.policy.parent.stage)
            .collect();
        let rank = |stage: &ParentStage| match stage {
            ParentStage::Early => 0,
            ParentStage::ReplacesEarly => 1,
            ParentStage::WhileUnproved => 2,
        };
        assert!(stages
            .windows(2)
            .all(|pair| rank(&pair[0]) <= rank(&pair[1])));
        assert_eq!(
            filtered(|route| route.policy.parent.stage == ParentStage::Early),
            [RealRootednessTheorem::ProductRecurrence]
        );
    }

    #[test]
    fn the_menus_follow_the_route_order() {
        let all = theorems(&ROUTES);
        // A residue class is not split again; the parent of an index section
        // may itself be proved by residue classes.
        assert_eq!(
            filtered(|route| route.policy.in_residue_class_menu),
            [&all[..9], &all[12..15]].concat()
        );
        assert_eq!(
            filtered(|route| route.policy.in_index_section_menu),
            [&all[..10], &all[12..15]].concat()
        );
    }

    #[test]
    fn the_start_policies_and_certificates_are_those_of_the_routes() {
        assert_eq!(
            filtered(|route| route.policy.start_policy == StartPolicy::OwnStarts),
            [
                RealRootednessTheorem::ResidueClassSubsequences,
                RealRootednessTheorem::IndexSectionLift,
                RealRootednessTheorem::RationalPfRiordan,
                RealRootednessTheorem::DegreeLawReversal
            ]
        );
        assert_eq!(
            filtered(|route| route.policy.certificate_name.is_some()),
            [
                RealRootednessTheorem::IteratedPreserver,
                RealRootednessTheorem::VariableDelay,
                RealRootednessTheorem::ResidueClassSubsequences,
                RealRootednessTheorem::IndexSectionLift,
                RealRootednessTheorem::FactorGauge,
                RealRootednessTheorem::LagTwoEuler,
                RealRootednessTheorem::BinomialPerturbedEulerOrbit,
                RealRootednessTheorem::RationalPfRiordan,
                RealRootednessTheorem::DegreeLawReversal
            ]
        );
    }
}
