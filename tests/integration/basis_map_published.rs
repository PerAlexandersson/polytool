//! Published basis maps with index-dependent rules (the Narayana, Eulerian
//! and derangement transformations): theorem recognition for every index,
//! the index-side transport through Riordan coefficient families, forward
//! and backward lifts to OEIS catalog targets, replay, tampering and
//! negative cases; and row-factor lifts with exceptional leading rows
//! (A008518).

use crate::common::q;
use polytool::oeis;
use polytool::recurrence::{
    apply_basis_map_to_family, certify_recurrence_equality, check_basis_map_preserves_pf,
    find_basis_map_lifts, find_row_factor_lifts, BasisMap, BasisMapIdentityStatus,
    BasisMapImageConstruction, BasisMapImageOptions, BasisMapLibrary, BasisMapLibraryOptions,
    BasisMapLiftCertification, BasisMapLiftSearchOptions, BasisMapShape, BasisMapTheoremOptions,
    BasisMapTheoremReport, BasisMapTheoremStatus, BigRational, BivarPoly, LiftVerdict, Recurrence,
    RecurrenceJson, RecurrenceTerm, RowFactorIdentityStatus, RowFactorLiftOptions,
};
use std::sync::Arc;

type Q = BigRational;

/// `(rule, initial rows, first row)` of a catalog entry at its displayed
/// indices.
fn catalog(id: &str) -> (Recurrence, Vec<Vec<Q>>, usize) {
    let entry = oeis::by_id(id).expect("catalog entry");
    let (rule, rows, first) = entry.recurrence_export_parts().expect("catalog recurrence");
    (
        rule,
        rows,
        usize::try_from(first).expect("nonnegative first row"),
    )
}

fn catalog_json(id: &str) -> RecurrenceJson {
    let (rule, rows, first) = catalog(id);
    RecurrenceJson::from_recurrence_rational(&rule, first, &rows, None)
}

/// `R_n = (u + c)^n`, `n >= 0`.
fn power_family(c: i64) -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q(c), q(1)]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

fn theorem(map: &BasisMap) -> BasisMapTheoremReport {
    check_basis_map_preserves_pf(map, &BasisMapTheoremOptions::default()).unwrap()
}

fn library(maps: Vec<BasisMap>) -> Arc<BasisMapLibrary> {
    Arc::new(
        BasisMapLibrary::build(&BasisMapLibraryOptions {
            maps,
            ..BasisMapLibraryOptions::default()
        })
        .unwrap(),
    )
}

/// The backward lifts of a catalog entry with the given maps.
fn backward(id: &str, maps: Vec<BasisMap>) -> Vec<BasisMapLiftCertification> {
    let (rule, rows, first) = catalog(id);
    find_basis_map_lifts(
        &rule,
        &rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(library(maps)),
            ..BasisMapLiftSearchOptions::default()
        },
    )
    .unwrap()
    .lifts
}

fn proved_backward(id: &str, map: BasisMap) -> BasisMapLiftCertification {
    let lifts = backward(id, vec![map]);
    let lift = lifts
        .into_iter()
        .find(BasisMapLiftCertification::is_proved)
        .unwrap_or_else(|| panic!("{id}: no proved lift"));
    lift.verify().unwrap();
    lift
}

#[test]
fn published_maps_are_recognised_for_every_index() {
    for (map, shape) in [
        (
            BasisMap::from_oeis("A001263").unwrap(),
            BasisMapShape::Narayana { parameter: 1 },
        ),
        (
            BasisMap::from_oeis("A008459").unwrap(),
            BasisMapShape::Narayana { parameter: 0 },
        ),
        (
            BasisMap::narayana_transformation(2),
            BasisMapShape::Narayana { parameter: 2 },
        ),
        (
            BasisMap::from_oeis("A123125").unwrap(),
            BasisMapShape::EulerianTransformation,
        ),
        (
            BasisMap::derangement_transformation(),
            BasisMapShape::DerangementTransformation,
        ),
    ] {
        let report = theorem(&map);
        assert_eq!(report.shape, shape, "{}", map.name());
        assert!(report.is_proved(), "{report:#?}");
        report.verify().unwrap();
        let pf = report.pf_conclusion().expect("a PF conclusion");
        let expected_source = match shape {
            BasisMapShape::Narayana { .. } => "(-inf, 0]",
            _ => "[-1, 0]",
        };
        assert_eq!(pf.input_interval.to_string(), expected_source);
        assert_eq!(pf.output_interval.to_string(), "(-inf, 0]");
        let lean = match shape {
            BasisMapShape::Narayana { .. } => "RealRooted.narayanaTransformPreservesPF",
            BasisMapShape::EulerianTransformation => {
                "RealRooted.EulerianTransform.transform_magicExpansion_isPF"
            }
            _ => "RealRooted.DerangementTransform.transform_magicExpansion_isPF",
        };
        assert!(report.lean_declarations.iter().any(|d| d == lean));
    }
}

#[test]
fn a_family_that_only_starts_like_a_published_one_gets_no_conclusion() {
    // P_n = (n + (n-1)(n-2)...(n-6)) t P_(n-1) + t (1 - t) P'_(n-1): the
    // Eulerian-transformation rows up to P_6, then different.
    let product = [720, -1764, 1624, -735, 175, -21, 1];
    let mut coeffs: Vec<Vec<Q>> = product.iter().map(|&c| vec![q(0), q(c)]).collect();
    coeffs[1][1] += q(1);
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(coeffs)),
            RecurrenceTerm::new(
                1,
                1,
                BivarPoly::from_coefficients_normalized(vec![vec![q(0), q(1), q(-1)]]),
            ),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let map = BasisMap::new("lookalike", rule, 0, vec![vec![q(1)]]).unwrap();
    let rows = map.basis_rows(8).unwrap();
    let eulerian = BasisMap::eulerian_transformation().basis_rows(8).unwrap();
    assert_eq!(rows[..7], eulerian[..7]);
    assert_ne!(rows[7], eulerian[7]);
    let report = theorem(&map);
    assert_eq!(report.shape, BasisMapShape::EulerianTransformation);
    assert_eq!(report.status, BasisMapTheoremStatus::HypothesisFailed);
    assert!(report.conclusions.is_empty());
}

#[test]
fn a046802_is_the_eulerian_transformation_of_the_binomial_rows() {
    // Forward: A°((1 + u)^n), certified by the index-side transport.
    let image = apply_basis_map_to_family(
        &BasisMap::eulerian_transformation(),
        &power_family(1),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert_eq!(image.construction, BasisMapImageConstruction::Transported);
    assert!(image.is_proved());
    image.certification.verify().unwrap();
    let equality = certify_recurrence_equality(&catalog_json("A046802"), &image.image).unwrap();
    assert!(equality.is_certified(), "{equality:#?}");
    // Backward: the expansion of A046802 in the basis 1, t A_k(t).
    let lift = proved_backward("A046802", BasisMap::eulerian_transformation());
    let (_, _, coefficient_seeds) = lift
        .candidate
        .coefficient_family
        .to_recurrence_parts()
        .unwrap();
    assert_eq!(coefficient_seeds[0], vec![q(1)]);
    assert_eq!(lift.lift_summary.verdict, LiftVerdict::Proved);
}

#[test]
fn narayana_transformations_reach_a162303_and_a257142() {
    // Forward: N_1 applied to the proved A026729 (Riordan: 1 / (1 - u z (1 + z))).
    let image = apply_basis_map_to_family(
        &BasisMap::narayana_transformation(1),
        &catalog_json("A026729"),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert!(image.is_proved());
    let equality = certify_recurrence_equality(&catalog_json("A162303"), &image.image).unwrap();
    assert!(equality.is_certified(), "{equality:#?}");
    // Backward lifts from the catalog rules.
    for id in ["A162303", "A257142"] {
        let lift = proved_backward(id, BasisMap::narayana_transformation(1));
        assert!(matches!(
            lift.identity,
            BasisMapIdentityStatus::Certified { .. }
        ));
    }
}

#[test]
fn non_riordan_coefficient_families_are_not_certified() {
    // A318397 = N_0(A008459), A095801 = N_1(A001263): the coefficient
    // families have index-dependent rules, so neither transport applies.
    for (id, map) in [
        ("A318397", BasisMap::narayana_transformation(0)),
        ("A095801", BasisMap::narayana_transformation(1)),
    ] {
        let lifts = backward(id, vec![map]);
        assert!(!lifts.is_empty(), "{id}");
        for lift in &lifts {
            assert!(!lift.is_proved(), "{id}");
            assert!(
                matches!(lift.identity, BasisMapIdentityStatus::Unsupported { .. }),
                "{id}: {:?}",
                lift.identity
            );
        }
    }
}

#[test]
fn inputs_outside_the_hypotheses_are_not_proved() {
    // A° needs roots in [-1, 0]: (u + 2)^n has the root -2.
    let image = apply_basis_map_to_family(
        &BasisMap::eulerian_transformation(),
        &power_family(2),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert!(matches!(
        image.certification.identity,
        BasisMapIdentityStatus::Certified { .. }
    ));
    assert!(!image.is_proved());
    // Narayana needs roots <= 0: (u - 1)^n has the root 1.
    let image = apply_basis_map_to_family(
        &BasisMap::narayana_transformation(0),
        &power_family(-1),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert!(!image.is_proved());
}

#[test]
fn the_derangement_transformation_is_a_forward_map() {
    // D((1 + u)^n) = A_n(t), the Eulerian polynomials of A173018.
    let image = apply_basis_map_to_family(
        &BasisMap::derangement_transformation(),
        &power_family(1),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert!(image.is_proved());
    image.certification.verify().unwrap();
    let equality = certify_recurrence_equality(&catalog_json("A173018"), &image.image).unwrap();
    assert!(equality.is_certified(), "{equality:#?}");
    // Not injective (d_1 = 0): never a library map for the backward lift.
    let mut library = BasisMapLibrary::default();
    assert!(!library
        .admit(
            BasisMap::derangement_transformation(),
            &BasisMapTheoremOptions::default()
        )
        .unwrap());
}

#[test]
fn the_catalog_library_ends_with_the_published_maps() {
    let library = BasisMapLibrary::oeis_catalog().unwrap();
    let names: Vec<&str> = library.entries().iter().map(|e| e.map.name()).collect();
    let tail = &names[names.len() - 4..];
    assert_eq!(
        tail,
        [
            "Narayana transformation N_(k,0)",
            "Narayana transformation N_(k,1)",
            "Narayana transformation N_(k,2)",
            "Eulerian transformation A°",
        ]
    );
    assert!(!names.contains(&"A001263"));
}

#[test]
fn tampered_published_certificates_are_rejected() {
    let lift = proved_backward("A046802", BasisMap::eulerian_transformation());
    let BasisMapIdentityStatus::Certified { certificate } = &lift.identity else {
        panic!("certified");
    };
    let mut tampered = (**certificate).clone();
    tampered.multiplier.push("1".into());
    assert!(tampered.verify().is_err());
    let mut tampered = (**certificate).clone();
    tampered.basis = BasisMap::narayana_transformation(1);
    assert!(tampered.verify().is_err());
    let mut tampered = (**certificate).clone();
    tampered.coefficient_family = power_family(2);
    assert!(tampered.verify().is_err());
    let mut report = lift.map_theorem.clone();
    report.shape = BasisMapShape::DerangementTransformation;
    assert!(report.verify().is_err());
    // The serialized lift replays.
    let json = serde_json::to_string(&lift).unwrap();
    let back: BasisMapLiftCertification = serde_json::from_str(&json).unwrap();
    back.verify().unwrap();
}

#[test]
fn a008518_is_a_row_factor_lift_with_an_exceptional_first_row() {
    let (rule, rows, first) = catalog("A008518");
    let report =
        find_row_factor_lifts(&rule, &rows, first, &RowFactorLiftOptions::default()).unwrap();
    report.verify().unwrap();
    let lift = report.proved_lifts().next().expect("a proved lift");
    assert_eq!(lift.factor, vec!["1".to_string(), "1".to_string()]);
    assert_eq!((lift.exponent_slope, lift.exponent_offset), (0, 1));
    let prefix = lift
        .exceptional_prefix
        .as_ref()
        .expect("an exceptional row");
    assert_eq!(prefix.rows, vec![vec!["1".to_string()]]);
    assert!(prefix.roots_nonpositive);
    assert_eq!(lift.cofactor.first_index, 1);
    assert_eq!(lift.lift_summary.input_first_index, 0);
    let RowFactorIdentityStatus::CertifiedByDivision { certificate } = &lift.identity else {
        panic!("certified by division: {:?}", lift.identity);
    };
    certificate.verify().unwrap();
    // Tampering.
    let mut tampered = (**certificate).clone();
    tampered.multiplier = vec!["2".into(), "1".into()];
    assert!(tampered.verify().is_err());
    let mut tampered = (**certificate).clone();
    tampered.exponent_offset = 2;
    assert!(tampered.verify().is_err());
    let mut tampered = lift.clone();
    tampered.exceptional_prefix.as_mut().unwrap().rows =
        vec![vec!["1".into(), "1".into(), "1".into()]];
    assert!(tampered.verify().is_err());
    // Without exceptional rows the input rule is rejected as before.
    let strict = RowFactorLiftOptions {
        max_exceptional_prefix_rows: 0,
        ..RowFactorLiftOptions::default()
    };
    assert!(find_row_factor_lifts(&rule, &rows, first, &strict).is_err());
}

#[test]
fn a_non_real_rooted_exceptional_row_fails_the_route() {
    // The A008518 rule with the first row replaced by 1 + x + x^2.
    let (rule, mut rows, first) = catalog("A008518");
    rows[0] = vec![q(1), q(1), q(1)];
    let report =
        find_row_factor_lifts(&rule, &rows, first, &RowFactorLiftOptions::default()).unwrap();
    report.verify().unwrap();
    assert_eq!(report.proved_lifts().count(), 0);
    let lift = report.lifts.first().expect("a lift record");
    assert!(!lift.exceptional_prefix.as_ref().unwrap().real_rooted);
    assert_eq!(lift.lift_summary.verdict, LiftVerdict::RouteFailed);
}

/// Every catalog entry: backward lifts with the published maps only, and
/// row-factor lifts that need an exceptional prefix or the division route.
/// Prints the proved ids.  `cargo test --release --test
/// basis_map_published -- --ignored --nocapture`.
#[test]
#[ignore = "catalog sweep, a few minutes in release; run with --ignored"]
fn catalog_sweep_of_the_published_maps() {
    let published = library(vec![
        BasisMap::narayana_transformation(0),
        BasisMap::narayana_transformation(1),
        BasisMap::narayana_transformation(2),
        BasisMap::eulerian_transformation(),
    ]);
    let mut basis = Vec::new();
    let mut row_factor = Vec::new();
    for entry in oeis::catalog() {
        let Ok((rule, rows, first)) = entry.recurrence_export_parts() else {
            continue;
        };
        let Ok(first) = usize::try_from(first) else {
            continue;
        };
        if let Ok(report) = find_basis_map_lifts(
            &rule,
            &rows,
            first,
            &BasisMapLiftSearchOptions {
                library: Some(published.clone()),
                ..BasisMapLiftSearchOptions::default()
            },
        ) {
            for lift in report.lifts.iter().filter(|l| l.is_proved()) {
                lift.verify().unwrap();
                basis.push(format!("{} <- {}", entry.id, lift.candidate.basis.name()));
            }
        }
        if let Ok(report) =
            find_row_factor_lifts(&rule, &rows, first, &RowFactorLiftOptions::default())
        {
            for lift in report.proved_lifts() {
                if lift.exceptional_prefix.is_some()
                    || matches!(
                        lift.identity,
                        RowFactorIdentityStatus::CertifiedByDivision { .. }
                    )
                {
                    lift.verify().unwrap();
                    row_factor.push(format!(
                        "{} (prefix rows {})",
                        entry.id,
                        lift.exceptional_prefix.as_ref().map_or(0, |p| p.rows.len())
                    ));
                }
            }
        }
    }
    eprintln!("published basis-map lifts: {basis:#?}");
    eprintln!("new row-factor lifts: {row_factor:#?}");
    assert!(basis.iter().any(|s| s.starts_with("A046802")));
    assert!(row_factor.iter().any(|s| s.starts_with("A008518")));
}
