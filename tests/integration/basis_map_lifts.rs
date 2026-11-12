//! Basis maps `T_P(u^k) = P_k(t)`: the theorem checks on worked examples
//! (OEIS catalog families), exact consistency checks of their conclusions,
//! negative results with exact counterexamples, the forward move, the
//! backward lift with a near miss, and replay and tampering of the
//! certificates.

use num_bigint::BigInt;
use polytool::recurrence::{
    apply_basis_map_to_family, certify_basis_map_lift, check_basis_map_preserves_pf,
    find_basis_map_lifts, AlgebraicInterval, BasisMap, BasisMapIdentityStatus,
    BasisMapImageConstruction, BasisMapImageOptions, BasisMapLibrary, BasisMapLibraryOptions,
    BasisMapLiftCandidate, BasisMapLiftCertification, BasisMapLiftOptions, BasisMapLiftRejection,
    BasisMapLiftSearchOptions, BasisMapShape, BasisMapTheoremOptions, BasisMapTheoremReport,
    BasisMapTheoremStatus, BigRational, BivarPoly, LiftParentProof, LiftVerdict, QuadraticNumber,
    Recurrence, RecurrenceJson, RecurrenceTerm,
};
use polytool::vec_poly;
use std::sync::Arc;

type Q = BigRational;

use crate::common::q;

fn qq(numerator: i64, denominator: i64) -> Q {
    Q::new(BigInt::from(numerator), BigInt::from(denominator))
}

fn poly(coefficients: &[i64]) -> Vec<Q> {
    coefficients.iter().map(|&c| q(c)).collect()
}

/// `sum_(i, j) c n^i x^j` from `(i, j, c)` triples.
fn bivariate(monomials: &[(usize, usize, Q)]) -> BivarPoly {
    let rows = monomials.iter().map(|m| m.0).max().unwrap_or(0) + 1;
    let width = monomials.iter().map(|m| m.1).max().unwrap_or(0) + 1;
    let mut coeffs = vec![vec![q(0); width]; rows];
    for (i, j, c) in monomials {
        coeffs[*i][*j] += c.clone();
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

fn theorem(map: &BasisMap) -> BasisMapTheoremReport {
    check_basis_map_preserves_pf(map, &BasisMapTheoremOptions::default()).unwrap()
}

fn rational(value: Q) -> QuadraticNumber {
    QuadraticNumber::rational(value)
}

/// `R_n = (u + c)^n`, `n >= 0`.
fn power_family(c: Q) -> RecurrenceJson {
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            bivariate(&[(0, 0, c), (0, 1, q(1))]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    RecurrenceJson::from_recurrence_rational(&rule, 0, &[vec![q(1)]], None)
}

fn is_real_rooted(p: &[Q]) -> bool {
    AlgebraicInterval::whole_line().contains_all_roots_of(p)
}

#[test]
fn a035607_has_the_sharp_pf_radius_one_plus_sqrt_two() {
    // P_0 = 1, P_1 = 1 + 2t, P_(n+2) = (1 + t) P_(n+1) + t P_n.
    let map = BasisMap::from_oeis("A035607").unwrap();
    assert_eq!(
        map.basis_rows(3).unwrap(),
        vec![poly(&[1]), poly(&[1, 2]), poly(&[1, 4, 2])]
    );
    let report = theorem(&map);
    assert_eq!(report.status, BasisMapTheoremStatus::Proved);
    assert!(matches!(report.shape, BasisMapShape::AffineCurrent { .. }));
    let pf = report.pf_conclusion().unwrap();
    let radius = QuadraticNumber::new(q(1), q(1), q(2)).unwrap();
    assert_eq!(
        pf.input_interval,
        AlgebraicInterval::closed(radius.negated(), rational(q(0)))
    );
    assert!(pf.nonnegative_coefficients);
    // The target interval I = [-3 - 2 sqrt(2), -3 + 2 sqrt(2)] holds for E.
    let target = AlgebraicInterval::closed(
        QuadraticNumber::new(q(-3), q(-2), q(2)).unwrap(),
        QuadraticNumber::new(q(-3), q(2), q(2)).unwrap(),
    );
    assert!(report
        .conclusions
        .iter()
        .any(|c| c.output_interval == target));
    // The survey's failure on the full ray: T((u + 4)^3) is not real-rooted,
    // and -4 lies outside [-1 - sqrt(2), 0].
    let image = map.apply(&poly(&[64, 48, 12, 1])).unwrap();
    assert_eq!(image, poly(&[125, 150, 32, 2]));
    assert!(!is_real_rooted(&image));
    assert!(!pf.input_interval.contains(&rational(q(-4))));
    // Consistency (finite evidence): inputs with roots in the source.
    for r in [q(0), q(1), q(2), qq(12, 5)] {
        for n in 0..=7 {
            let input = vec_poly::pow(&[r.clone(), q(1)], n);
            let image = map.apply(&input).unwrap();
            assert!(
                pf.output_interval.contains_all_roots_of(&image),
                "r = {r}, n = {n}"
            );
            assert!(image.iter().all(|c| *c >= q(0)));
        }
    }
    report.verify().unwrap();
}

#[test]
fn the_survey_corollaries_hold() {
    // A099040 preserves [-2, 0]; A099097 sends [-2/3, 0] into [-4/9, 0].
    for (id, source, target) in [("A099040", q(-2), q(-2)), ("A099097", qq(-2, 3), qq(-4, 9))] {
        let report = theorem(&BasisMap::from_oeis(id).unwrap());
        let pf = report.pf_conclusion().unwrap();
        assert_eq!(
            pf.input_interval,
            AlgebraicInterval::closed(rational(source), rational(q(0))),
            "{id}"
        );
        assert_eq!(
            pf.output_interval,
            AlgebraicInterval::closed(rational(target), rational(q(0))),
            "{id}"
        );
    }
    // A056242: real polynomials supported in [0, 1] go to [-1, 0]; no PF
    // statement (the radius is zero).
    let report = theorem(&BasisMap::from_oeis("A056242").unwrap());
    assert!(report.pf_conclusion().is_none());
    let unit = AlgebraicInterval::closed(rational(q(0)), rational(q(1)));
    let minus_unit = AlgebraicInterval::closed(rational(q(-1)), rational(q(0)));
    assert!(report
        .conclusions
        .iter()
        .any(|c| c.input_interval == unit && c.output_interval == minus_unit));
    // T((u - 1/2)^3): a root inside [0, 1], image in [-1, 0].
    let map = BasisMap::from_oeis("A056242").unwrap();
    let image = map.apply(&[qq(-1, 8), qq(3, 4), qq(-3, 2), q(1)]).unwrap();
    assert!(minus_unit.contains_all_roots_of(&image));
}

#[test]
fn constant_current_maps_follow_the_scalar_criterion() {
    // A034839: (a, b, c, s) = (2, 1, -1, 1), PF with a negative lag constant.
    let report = theorem(&BasisMap::from_oeis("A034839").unwrap());
    assert!(matches!(
        report.shape,
        BasisMapShape::ConstantCurrent { .. }
    ));
    let pf = report.pf_conclusion().unwrap();
    assert_eq!(pf.input_interval, AlgebraicInterval::nonpositive_ray());
    // a = s = b = 1, c = -1/2: P_4 = t^2 + 2t - 1/4 is not PF, real-rootedness
    // with roots <= 1/2 remains.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivariate(&[(0, 0, q(1))])),
            RecurrenceTerm::new(2, 0, bivariate(&[(0, 0, qq(-1, 2)), (0, 1, q(1))])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let map = BasisMap::new("c = -1/2", rule, 0, vec![poly(&[1]), poly(&[1])]).unwrap();
    assert_eq!(map.basis_rows(5).unwrap()[4], vec![qq(-1, 4), q(2), q(1)]);
    let report = theorem(&map);
    assert_eq!(report.status, BasisMapTheoremStatus::Proved);
    assert!(report.pf_conclusion().is_none());
    let counterexample = report.counterexample.as_ref().unwrap();
    assert_eq!(counterexample.input.len(), 5);
    assert!(report
        .conclusions
        .iter()
        .any(|c| c.output_interval == AlgebraicInterval::at_most(rational(qq(1, 2)))));
}

#[test]
fn geometric_and_differential_maps() {
    // Binomial basis (1 + t)^k: roots <= 1 go to roots <= 0.
    let report = theorem(&BasisMap::from_oeis("A007318").unwrap());
    assert!(matches!(report.shape, BasisMapShape::Geometric { .. }));
    assert_eq!(
        report.pf_conclusion().unwrap().input_interval,
        AlgebraicInterval::at_most(rational(q(1)))
    );
    // Ordered Bell rows F_(n+1) = t F_n + t (t + 1) F_n': [-1, 0] -> [-1, 0].
    let report = theorem(&BasisMap::from_oeis("A131689").unwrap());
    assert!(matches!(
        report.shape,
        BasisMapShape::IntervalDifferential { .. }
    ));
    let unit = AlgebraicInterval::closed(rational(q(-1)), rational(q(0)));
    let conclusion = report.pf_conclusion().unwrap();
    assert_eq!(conclusion.input_interval, unit);
    assert_eq!(conclusion.output_interval, unit);
    assert!(conclusion.sharp);
    // Bell rows B_(n+1) = t B_n + t B_n': a ray map.
    let report = theorem(&BasisMap::from_oeis("A048993").unwrap());
    assert!(matches!(
        report.shape,
        BasisMapShape::RayDifferential { .. }
    ));
    assert!(report.pf_conclusion().is_some());
    // Hermite H_(n+1) = 2t H_n - H_n': every real-rooted input.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivariate(&[(0, 1, q(2))])),
            RecurrenceTerm::new(1, 1, bivariate(&[(0, 0, q(-1))])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let map = BasisMap::new("hermite", rule, 0, vec![poly(&[1])]).unwrap();
    let report = theorem(&map);
    assert!(matches!(
        report.shape,
        BasisMapShape::HermiteDifferential { .. }
    ));
    assert_eq!(
        report.conclusions[0].input_interval,
        AlgebraicInterval::whole_line()
    );
    // Rising factorials t (t + 1) ... (t + n - 1): an index-dependent scalar.
    let report = theorem(&BasisMap::from_oeis("A132393").unwrap());
    assert!(matches!(report.shape, BasisMapShape::ScalarIndex { .. }));
    assert_eq!(
        report.pf_conclusion().unwrap().input_interval,
        AlgebraicInterval::nonpositive_ray()
    );
}

#[test]
fn negative_results_carry_exact_counterexamples() {
    // Catalog indexing of A071919: rows 1, 1, then (1 + t)^(n-1); an
    // eventually geometric family with a prefix discrepancy.
    let map = BasisMap::from_oeis("A071919").unwrap();
    let report = theorem(&map);
    assert!(matches!(
        report.shape,
        BasisMapShape::EventuallyGeometric { .. }
    ));
    assert_eq!(report.status, BasisMapTheoremStatus::NotPreserving);
    let counterexample = report.counterexample.as_ref().unwrap();
    let image: Vec<Q> = counterexample
        .image
        .iter()
        .map(|c| polytool::recurrence::parse_rational_coeff(c).unwrap())
        .collect();
    assert!(!is_real_rooted(&image));
    // A Newton basis with decreasing nodes (A039683).
    let report = theorem(&BasisMap::from_oeis("A039683").unwrap());
    assert_eq!(report.status, BasisMapTheoremStatus::NotPreserving);
    assert!(report.counterexample.is_some());
}

fn a035607_library() -> Arc<BasisMapLibrary> {
    Arc::new(
        BasisMapLibrary::build(&BasisMapLibraryOptions {
            maps: vec![BasisMap::from_oeis("A035607").unwrap()],
            ..BasisMapLibraryOptions::default()
        })
        .unwrap(),
    )
}

#[test]
fn the_forward_move_certifies_images() {
    let map = BasisMap::from_oeis("A035607").unwrap();
    for c in [q(1), q(2)] {
        let image = apply_basis_map_to_family(
            &map,
            &power_family(c.clone()),
            &BasisMapImageOptions::default(),
        )
        .unwrap();
        assert_eq!(image.construction, BasisMapImageConstruction::Eliminated);
        assert!(image.is_proved(), "c = {c}");
        let (rule, n0, seeds) = image.image.to_recurrence_parts().unwrap();
        let rows = rule.generate_rows_rational(&seeds, n0, 8).unwrap();
        for (n, row) in rows.iter().enumerate() {
            let input = vec_poly::pow(&[c.clone(), q(1)], n);
            assert_eq!(vec_poly::trim(row), map.apply(&input).unwrap());
        }
        let BasisMapIdentityStatus::Certified { certificate } = &image.certification.identity
        else {
            panic!("not certified");
        };
        certificate.verify().unwrap();
        let json = serde_json::to_string(&image).unwrap();
        let back: polytool::recurrence::BasisMapImage = serde_json::from_str(&json).unwrap();
        assert_eq!(back, image);
    }
    // c = 4 lies outside [-1 - sqrt(2), 0] (as a root -4): the identity is
    // certified, but the route fails at the root check.
    let image =
        apply_basis_map_to_family(&map, &power_family(q(4)), &BasisMapImageOptions::default())
            .unwrap();
    assert_eq!(
        image.certification.lift_summary.verdict,
        LiftVerdict::RouteFailed
    );
    assert!(matches!(
        image.certification.lift_summary.parent_proof,
        LiftParentProof::Failed { .. }
    ));
    // A map of order one gives the transported rule directly.
    let binomial = BasisMap::from_oeis("A007318").unwrap();
    let image = apply_basis_map_to_family(
        &binomial,
        &power_family(q(1)),
        &BasisMapImageOptions::default(),
    )
    .unwrap();
    assert_eq!(image.construction, BasisMapImageConstruction::Transported);
    assert!(image.is_proved());
}

/// The image of `(u + c)^n` under A035607 as an input rule with rows.
fn image_input(c: i64) -> (Recurrence, Vec<Vec<Q>>, usize) {
    let map = BasisMap::from_oeis("A035607").unwrap();
    let image =
        apply_basis_map_to_family(&map, &power_family(q(c)), &BasisMapImageOptions::default())
            .unwrap();
    let (rule, n0, seeds) = image.image.to_recurrence_parts().unwrap();
    let rows = rule.generate_rows_rational(&seeds, n0, 16).unwrap();
    (rule, rows, n0)
}

#[test]
fn the_backward_lift_recovers_the_coefficient_family() {
    let (rule, rows, n0) = image_input(2);
    let options = BasisMapLiftSearchOptions {
        library: Some(a035607_library()),
        ..BasisMapLiftSearchOptions::default()
    };
    let report = find_basis_map_lifts(&rule, &rows, n0, &options).unwrap();
    let lift = report.proved().expect("a proved lift");
    assert_eq!(lift.candidate.basis.name(), "A035607");
    let (r_rule, r0, r_seeds) = lift
        .candidate
        .coefficient_family
        .to_recurrence_parts()
        .unwrap();
    let r_rows = r_rule.generate_rows_rational(&r_seeds, r0, 8).unwrap();
    for (n, row) in r_rows.iter().enumerate() {
        assert_eq!(vec_poly::trim(row), vec_poly::pow(&[q(2), q(1)], n));
    }
    assert!(lift.lift_summary.is_consistent());
    report.verify().unwrap();
    // A certification whose coefficient family was swapped does not replay.
    let mut tampered = lift.clone();
    tampered.candidate.coefficient_family = power_family(q(1));
    assert!(tampered.verify().is_err());
    let json = serde_json::to_string(&report).unwrap();
    let back: polytool::recurrence::BasisMapLiftReport = serde_json::from_str(&json).unwrap();
    assert_eq!(back, report);
}

#[test]
fn a_near_miss_is_rejected() {
    // T((u + 4)^n): the coefficient rows have the root -4, outside the source.
    let (rule, rows, n0) = image_input(4);
    let options = BasisMapLiftSearchOptions {
        library: Some(a035607_library()),
        ..BasisMapLiftSearchOptions::default()
    };
    let report = find_basis_map_lifts(&rule, &rows, n0, &options).unwrap();
    assert!(report.proved().is_none());
    assert!(matches!(
        report.rejections.as_slice(),
        [BasisMapLiftRejection::RootOutsideSource { index: 1, .. }]
    ));
    // Certifying the candidate directly gives a failed route, not a proof.
    let candidate = BasisMapLiftCandidate {
        basis: BasisMap::from_oeis("A035607").unwrap(),
        coefficient_family: power_family(q(4)),
    };
    let certification = certify_basis_map_lift(
        &rule,
        &rows,
        n0,
        &candidate,
        &BasisMapLiftOptions::default(),
    )
    .unwrap();
    assert!(matches!(
        certification.identity,
        BasisMapIdentityStatus::Certified { .. }
    ));
    assert_eq!(certification.lift_summary.verdict, LiftVerdict::RouteFailed);
    // A wrong coefficient family is refuted at a row.
    let wrong = BasisMapLiftCandidate {
        coefficient_family: power_family(q(3)),
        ..candidate
    };
    let certification =
        certify_basis_map_lift(&rule, &rows, n0, &wrong, &BasisMapLiftOptions::default()).unwrap();
    assert!(matches!(
        certification.identity,
        BasisMapIdentityStatus::Refuted { index: 1 }
    ));
}

fn certified(c: i64) -> BasisMapLiftCertification {
    let (rule, rows, n0) = image_input(c);
    let candidate = BasisMapLiftCandidate {
        basis: BasisMap::from_oeis("A035607").unwrap(),
        coefficient_family: power_family(q(c)),
    };
    certify_basis_map_lift(
        &rule,
        &rows,
        n0,
        &candidate,
        &BasisMapLiftOptions::default(),
    )
    .unwrap()
}

#[test]
fn tampered_certificates_are_rejected() {
    let certification = certified(1);
    assert!(certification.is_proved());
    let BasisMapIdentityStatus::Certified { certificate } = certification.identity else {
        panic!("not certified");
    };
    certificate.verify().unwrap();

    let mut tampered = certificate.clone();
    tampered.multiplier = vec!["2".into()];
    assert!(tampered.verify().is_err());

    let mut tampered = certificate.clone();
    tampered.exceptional_bases.push(40);
    assert!(tampered.verify().is_err());

    let mut tampered = certificate.clone();
    tampered.coefficient_family = power_family(q(2));
    assert!(tampered.verify().is_err());

    let mut tampered = certificate.clone();
    tampered.input.initial_polynomials[0] = vec!["2".into()];
    assert!(tampered.verify().is_err());

    let mut tampered = certificate.clone();
    if let Some(quotient) = tampered.quotients.first_mut() {
        quotient.shift += 1;
    }
    assert!(tampered.verify().is_err());

    // A theorem report with an altered conclusion does not replay.
    let mut report = theorem(&BasisMap::from_oeis("A035607").unwrap());
    report.conclusions[0].input_interval = AlgebraicInterval::nonpositive_ray();
    assert!(report.verify().is_err());
}

#[test]
fn intertwining_relations_and_expansions_hold() {
    for id in [
        "A035607", "A131689", "A039683", "A008288", "A034839", "A071919",
    ] {
        let map = BasisMap::from_oeis(id).unwrap();
        map.intertwining_relation()
            .check_on_monomials(&map, 6)
            .unwrap_or_else(|e| panic!("{id}: {e}"));
    }
    for id in ["A035607", "A131689", "A007318"] {
        let map = BasisMap::from_oeis(id).unwrap();
        let p = poly(&[3, -1, 0, 2, 5]);
        assert_eq!(map.expand(&map.apply(&p).unwrap()).unwrap(), p, "{id}");
    }
}
