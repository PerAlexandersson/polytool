//! Even-derivative PF preservers: all-index certificates and exact failures.

use num_traits::{One, Zero};
use polytool::recurrence::*;
use polytool::BigRational as Q;

fn ratio(numerator: &[&str], denominator: &[&str]) -> EvenDerivativeCoefficients {
    EvenDerivativeCoefficients::RationalRatio {
        numerator: numerator.iter().map(|s| (*s).into()).collect(),
        denominator: denominator.iter().map(|s| (*s).into()).collect(),
    }
}
fn finite(c: &[&str]) -> EvenDerivativeCoefficients {
    EvenDerivativeCoefficients::Finite {
        coefficients: c.iter().map(|s| (*s).into()).collect(),
    }
}
fn catalog(id: &str) -> (Recurrence, usize, Vec<Vec<Q>>) {
    let (rule, seeds, first) = polytool::oeis::by_id(id)
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    (rule, first.try_into().unwrap(), seeds)
}

#[test]
fn even_derivative_fixed_symbols_and_degree_five_failure() {
    for (c, proved) in [
        (finite(&["1", "2", "1"]), true),
        (finite(&["1", "1", "1"]), false),
    ] {
        let map = BasisMap::even_derivative(&c).unwrap();
        assert_eq!(map.basis_rows(12).unwrap(), c.basis_rows(12).unwrap());
        let report =
            check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default()).unwrap();
        assert!(matches!(report.shape, BasisMapShape::EvenDerivative { .. }));
        assert_eq!(report.is_proved(), proved, "{report:#?}");
        if !proved {
            assert_eq!(report.status, BasisMapTheoremStatus::NotPreserving);
            let witness = report.counterexample.as_ref().unwrap();
            assert_eq!(witness.input.len(), 6);
            assert_eq!(witness.image, ["0", "0", "0", "120", "20", "1"]);
        }
        report.verify().unwrap();
        let mut tampered = report.clone();
        if let BasisMapShape::EvenDerivative { certificate } = &mut tampered.shape {
            certificate.all_rows_pf = !certificate.all_rows_pf;
        }
        assert!(tampered.verify().is_err());
        let mut tampered = report;
        tampered.lean_declarations.push("forged".into());
        assert!(tampered.verify().is_err());
    }
}

#[test]
fn even_derivative_hypergeometric_rows_and_pf_proofs() {
    for c in [
        ratio(&["1"], &["1", "1"]),
        ratio(&["1"], &["1", "2", "1"]),
        ratio(&["1"], &["2", "3", "1"]),
    ] {
        let map = BasisMap::even_derivative(&c).unwrap();
        assert_eq!(map.basis_rows(18).unwrap(), c.basis_rows(18).unwrap());
        let report =
            check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default()).unwrap();
        eprintln!("ratio {c:?}: {:?}", report.status);
        assert!(report.is_proved(), "{report:#?}");
        report.verify().unwrap();
        if let BasisMapShape::EvenDerivative { certificate } = &report.shape {
            certificate.verify().unwrap();
        }
        let mut tampered = report;
        if let BasisMapShape::EvenDerivative { certificate } = &mut tampered.shape {
            certificate
                .triangle
                .as_mut()
                .unwrap()
                .record
                .initial_polynomials[0][0] = "2".into();
        }
        assert!(tampered.verify().is_err());
    }
}

#[test]
fn even_derivative_ratio_late_sign_failure_and_poles() {
    // Positive through k=50, negative thereafter, with no terminating zero.
    assert!(BasisMap::even_derivative(&ratio(&["101", "-2"], &["1", "1"])).is_err());
    assert!(BasisMap::even_derivative(&ratio(&["1"], &["50", "-1"])).is_err());
    assert!(BasisMap::even_derivative(&finite(&["2", "1"])).is_err());
}

#[test]
fn even_derivative_stirling_image_and_replay() {
    let map = BasisMap::even_derivative(&ratio(&["1"], &["1", "1"])).unwrap();
    let (parent, first, seeds) = catalog("A048993");
    let parent_rows = parent.generate_rows_rational(&seeds, first, 16).unwrap();
    let candidate = BasisMapLiftCandidate {
        basis: map,
        coefficient_family: RecurrenceJson::from_recurrence_rational(&parent, first, &seeds, None),
    };
    let (rule, n0, seeds) = catalog("A143395");
    let rows = rule.generate_rows_rational(&seeds, n0, 16).unwrap();
    assert_eq!(
        parent_rows
            .iter()
            .map(|r| candidate.basis.apply(r).unwrap())
            .collect::<Vec<_>>(),
        rows
    );
    let result = certify_basis_map_lift(
        &rule,
        &rows,
        n0,
        &candidate,
        &BasisMapLiftOptions::default(),
    )
    .unwrap();
    assert!(result.is_proved(), "{result:#?}");
    result.verify().unwrap();
    let mut unrelated = result.clone();
    let other_map = BasisMap::even_derivative(&ratio(&["2"], &["1", "1"])).unwrap();
    unrelated.map_theorem =
        check_basis_map_preserves_pf(&other_map, &BasisMapTheoremOptions::default()).unwrap();
    assert!(unrelated.verify().is_err());
    let mut bad = result;
    bad.candidate.coefficient_family.initial_polynomials[0][0] = "2".into();
    assert!(bad.verify().is_err());
}

#[test]
fn even_derivative_catalog_lift_discovers_a143395() {
    let (rule, first, seeds) = catalog("A143395");
    let record = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
    let result =
        discover_recurrence_proofs(&record, &RealRootednessTheoremOptions::default()).unwrap();
    assert_eq!(result.status, "proved", "{result:#?}");
    let stage = result.catalog_lifts.as_ref().unwrap();
    assert!(stage.is_proved());
    result.verify().unwrap();
}

#[test]
fn even_derivative_rows_keep_the_monomial_leading_coefficient() {
    let c = ratio(&["1"], &["1", "1"]);
    for (n, row) in c.basis_rows(20).unwrap().iter().enumerate() {
        assert!(row[n].is_one());
        if n > 0 {
            assert!(row[0].is_zero());
        }
    }
}

#[test]
fn even_derivative_non_riordan_acceptance_candidates_report_identity_blockers() {
    for (id, parent, c, reversed) in [
        ("A318397", "A104684", ratio(&["1"], &["1", "2", "1"]), false),
        ("A367177", "A104684", ratio(&["1"], &["1", "2", "1"]), true),
        ("A095801", "A126216", ratio(&["1"], &["2", "3", "1"]), false),
    ] {
        let map = BasisMap::even_derivative(&c).unwrap();
        let (parent, first, seeds) = catalog(parent);
        let candidate = BasisMapLiftCandidate {
            basis: map,
            coefficient_family: RecurrenceJson::from_recurrence_rational(
                &parent, first, &seeds, None,
            ),
        };
        let (mut rule, n0, seeds) = catalog(id);
        let mut rows = rule.generate_rows_rational(&seeds, n0, 16).unwrap();
        if reversed {
            let reversal = RecurrenceTransformation::Reversal {
                schedule: DegreeSchedule::new(n0, rows[0].len() - 1, 1),
            };
            rows = rows
                .iter()
                .enumerate()
                .map(|(i, row)| reversal.map_row(n0 + i, row).unwrap().remove(0))
                .collect();
            let TransportedRecurrence::Scalar(reversed_rule) = reversal.transport(&rule).unwrap()
            else {
                panic!("scalar reversal")
            };
            rule = reversed_rule;
        }
        let result = certify_basis_map_lift(
            &rule,
            &rows,
            n0,
            &candidate,
            &BasisMapLiftOptions::default(),
        )
        .unwrap();
        assert!(
            result.map_theorem.is_proved(),
            "{id}: {:?}",
            result.map_theorem.status
        );
        assert!(
            matches!(result.identity, BasisMapIdentityStatus::Unsupported { .. }),
            "{id}: {:?}",
            result.identity
        );
        assert!(!result.is_proved());
        eprintln!("acceptance {id}: {:?}", result.identity);
        result.verify().unwrap();
    }
}

#[test]
fn even_derivative_stirling_forward_transport_and_library_admission() {
    let map = BasisMap::even_derivative(&ratio(&["2"], &["1", "1"])).unwrap();
    let mut library = BasisMapLibrary::default();
    assert!(library
        .admit(map.clone(), &BasisMapTheoremOptions::default())
        .unwrap());
    let (rule, first, seeds) = catalog("A048993");
    let parent = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
    let image = apply_basis_map_to_family(&map, &parent, &BasisMapImageOptions::default()).unwrap();
    assert!(image.is_proved());
    image.certification.verify().unwrap();
}

#[test]
fn even_derivative_a_long_row_match_does_not_certify_the_basis() {
    let map = BasisMap::even_derivative(&ratio(&["1"], &["1", "1"])).unwrap();
    let mut rule = map.rule().clone();
    let mut perturbation = vec![Q::one()];
    for i in 1..=18 {
        perturbation =
            polytool::vec_poly::mul(&perturbation, &[-Q::from_integer(i.into()), Q::one()]);
    }
    rule.terms.push(RecurrenceTerm::new(
        1,
        0,
        BivarPoly::from_coefficients_normalized(
            perturbation.into_iter().map(|c| vec![c]).collect(),
        ),
    ));
    let impostor = BasisMap::new("long prefix", rule, 0, map.seeds().to_vec()).unwrap();
    assert_eq!(
        map.basis_rows(18).unwrap(),
        impostor.basis_rows(18).unwrap()
    );
    let report =
        check_basis_map_preserves_pf(&impostor, &BasisMapTheoremOptions::default()).unwrap();
    assert!(!report.is_proved());
    assert!(report.conclusions.is_empty());
}

#[test]
fn even_derivative_positive_ratio_does_not_replace_the_row_pf_hypothesis() {
    let map = BasisMap::even_derivative(&ratio(&["1"], &["1"])).unwrap();
    let report = check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default()).unwrap();
    assert_eq!(report.status, BasisMapTheoremStatus::NotPreserving);
    assert!(report.conclusions.is_empty());
    assert_eq!(report.counterexample.as_ref().unwrap().input.len(), 6);
    report.verify().unwrap();
}

#[test]
fn even_derivative_missing_row_proof_is_unresolved_without_a_witness() {
    let map = BasisMap::even_derivative(&ratio(&["1"], &["1"])).unwrap();
    let report = check_basis_map_preserves_pf(
        &map,
        &BasisMapTheoremOptions {
            witness_search_limit: 0,
        },
    )
    .unwrap();
    assert!(report.conclusions.is_empty());
    assert!(report.counterexample.is_none());
    assert_eq!(
        report
            .hypotheses
            .iter()
            .find(|h| h.name == "even_derivative_all_rows_pf")
            .unwrap()
            .status,
        LiftHypothesisStatus::Unresolved
    );
}

#[test]
fn even_derivative_recovery_respects_explicit_libraries_and_budgets() {
    use std::sync::Arc;
    let (rule, first, seeds) = catalog("A143395");
    let rows = rule.generate_rows_rational(&seeds, first, 16).unwrap();
    let empty = find_basis_map_lifts(
        &rule,
        &rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(Arc::new(BasisMapLibrary::default())),
            ..BasisMapLiftSearchOptions::default()
        },
    )
    .unwrap();
    assert!(empty.lifts.is_empty());
    assert_eq!(empty.maps_tried, 0);

    let library = Arc::new(
        BasisMapLibrary::build(&BasisMapLibraryOptions {
            maps: vec![BasisMap::even_derivative(&ratio(&["1"], &["1", "1"])).unwrap()],
            ..BasisMapLibraryOptions::default()
        })
        .unwrap(),
    );
    let exhausted = find_basis_map_lifts(
        &rule,
        &rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(library.clone()),
            max_maps: 0,
            ..BasisMapLiftSearchOptions::default()
        },
    )
    .unwrap();
    assert!(exhausted.exhausted_budget);
    assert!(exhausted.lifts.is_empty());
    let found = find_basis_map_lifts(
        &rule,
        &rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(library),
            ..BasisMapLiftSearchOptions::default()
        },
    )
    .unwrap();
    assert!(found.proved().is_some());
    found.verify().unwrap();
}

#[test]
fn even_derivative_scalar_transport_uses_recurrence_equality() {
    for (id, c) in [
        ("A008459", ratio(&["1"], &["1", "2", "1"])),
        ("A001263", ratio(&["1"], &["2", "3", "1"])),
    ] {
        let (rule, first, seeds) = catalog(id);
        let one = Q::one();
        let parent = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![one.clone(), one.clone()]]),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        let candidate = BasisMapLiftCandidate {
            basis: BasisMap::even_derivative(&c).unwrap(),
            coefficient_family: RecurrenceJson::from_recurrence_rational(
                &parent,
                first,
                &[vec![one]],
                None,
            ),
        };
        let rows = rule.generate_rows_rational(&seeds, first, 16).unwrap();
        let mut result = certify_basis_map_lift(
            &rule,
            &rows,
            first,
            &candidate,
            &BasisMapLiftOptions::default(),
        )
        .unwrap();
        assert!(result.is_proved());
        assert!(matches!(
            result.identity,
            BasisMapIdentityStatus::RecurrenceEquality { .. }
        ));
        result.verify().unwrap();
        if let BasisMapIdentityStatus::RecurrenceEquality { certificate } = &mut result.identity {
            certificate.alternative.initial_polynomials[0][0] = "2".into();
        }
        assert!(result.verify().is_err());
    }
}

#[test]
fn even_derivative_finite_symbol_disproves_without_a_small_witness() {
    let map = BasisMap::even_derivative(&finite(&["1", "1", "1"])).unwrap();
    let report = check_basis_map_preserves_pf(
        &map,
        &BasisMapTheoremOptions {
            witness_search_limit: 0,
        },
    )
    .unwrap();
    assert_eq!(report.status, BasisMapTheoremStatus::NotPreserving);
    assert!(report.counterexample.is_none());
    report.verify().unwrap();
}
