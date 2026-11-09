//! The registry of real-rootedness theorem routes: names, the frozen route
//! order, and the parent routes that run a theorem.

use polytool::recurrence::{RealRootednessTheorem, VeroneseParentProofRoute};

/// The serialized name of a value of a unit-variant enum.
fn serialized<T: serde::Serialize>(value: T) -> String {
    serde_json::to_value(value)
        .unwrap()
        .as_str()
        .unwrap()
        .to_string()
}

/// Every parent route (the enum has no `ALL`).
const PARENT_ROUTES: [VeroneseParentProofRoute; 16] = [
    VeroneseParentProofRoute::PeriodicProduct,
    VeroneseParentProofRoute::ProductRecurrence,
    VeroneseParentProofRoute::PaperTheorem,
    VeroneseParentProofRoute::RootDomain,
    VeroneseParentProofRoute::ThreeTermRecurrence,
    VeroneseParentProofRoute::LiuWang,
    VeroneseParentProofRoute::LiuWangDerivativeLag,
    VeroneseParentProofRoute::MaWangInterval,
    VeroneseParentProofRoute::LiuWangThreeTerm,
    VeroneseParentProofRoute::IteratedPreserver,
    VeroneseParentProofRoute::HypergeometricRow,
    VeroneseParentProofRoute::VariableDelay,
    VeroneseParentProofRoute::ResidueClassSubsequences,
    VeroneseParentProofRoute::IndexSectionLift,
    VeroneseParentProofRoute::FactorGauge,
    VeroneseParentProofRoute::Attestation,
];

/// Fails to compile when a parent route is added without extending
/// [`PARENT_ROUTES`].
#[allow(dead_code)]
fn parent_routes_are_listed(route: VeroneseParentProofRoute) {
    use VeroneseParentProofRoute as R;
    match route {
        R::PeriodicProduct
        | R::ProductRecurrence
        | R::PaperTheorem
        | R::RootDomain
        | R::ThreeTermRecurrence
        | R::LiuWang
        | R::LiuWangDerivativeLag
        | R::MaWangInterval
        | R::LiuWangThreeTerm
        | R::IteratedPreserver
        | R::HypergeometricRow
        | R::VariableDelay
        | R::ResidueClassSubsequences
        | R::IndexSectionLift
        | R::FactorGauge
        | R::Attestation => {}
    }
}

#[test]
fn theorem_names_are_their_serialized_names() {
    for theorem in RealRootednessTheorem::ALL {
        assert_eq!(theorem.as_str(), serialized(theorem));
        let parsed: RealRootednessTheorem =
            serde_json::from_value(serde_json::Value::from(theorem.as_str())).unwrap();
        assert_eq!(parsed, theorem);
    }
    for route in PARENT_ROUTES {
        assert_eq!(route.as_str(), serialized(route));
    }
}

/// The route order is policy: it is the order of the parent cascade, of the
/// discovery attempts and of the residue-class and index-section menus.
#[test]
fn the_theorem_order_is_frozen() {
    let names: Vec<&str> = RealRootednessTheorem::ALL
        .iter()
        .map(|theorem| theorem.as_str())
        .collect();
    assert_eq!(
        names,
        [
            "product_recurrence",
            "three_term_recurrence",
            "liu_wang",
            "liu_wang_derivative_lag",
            "ma_wang_interval",
            "liu_wang_three_term",
            "iterated_preserver",
            "hypergeometric_row",
            "variable_delay",
            "residue_class_subsequences",
            "index_section_lift",
            "factor_gauge",
        ]
    );
}

#[test]
fn every_theorem_has_exactly_one_parent_route_of_the_same_name() {
    for theorem in RealRootednessTheorem::ALL {
        let routes: Vec<VeroneseParentProofRoute> = PARENT_ROUTES
            .into_iter()
            .filter(|route| route.theorem() == Some(theorem))
            .collect();
        assert_eq!(routes.len(), 1, "{theorem:?}");
        assert_eq!(routes[0].as_str(), theorem.as_str());
    }
    let without_theorem: Vec<&str> = PARENT_ROUTES
        .into_iter()
        .filter(|route| route.theorem().is_none())
        .map(VeroneseParentProofRoute::as_str)
        .collect();
    assert_eq!(
        without_theorem,
        ["periodic_product", "paper_theorem", "root_domain", "attestation"]
    );
}
