//! Exact checks of the Euler-orbit reduction and its catalog coordinates.

use num_bigint::BigInt;
use num_rational::BigRational as Q;
use num_traits::One;
use polytool::recurrence::*;
use polytool::vec_poly;

mod common;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn bivar(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().copied().map(q).collect())
            .collect(),
    )
}

fn orbit_rule(s: i64) -> Recurrence {
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[-3, -6], &[2, 4]])),
            RecurrenceTerm::new(1, 1, bivar(&[&[0, -2 * s, -2 * s], &[0, s, s]])),
            RecurrenceTerm::new(2, 0, bivar(&[&[1, 4, 4], &[-1, -4, -4]])),
            RecurrenceTerm::new(
                2,
                1,
                bivar(&[&[0, s, 3 * s, 2 * s], &[0, -s, -3 * s, -2 * s]]),
            ),
        ],
        denominator: Some(bivar(&[&[-2], &[1]])),
        inhomogeneous: None,
    }
}

fn orbit_seeds(kappa: i64) -> Vec<Vec<Q>> {
    vec![vec![q(1)], vec![q(1), q(2)], vec![q(1), q(kappa), q(kappa)]]
}

const TARGETS: [(&str, i64, i64); 5] = [
    ("A141689", 1, 5),
    ("A141690", 1, 8),
    ("A141696", 1, 10),
    ("A141697", 1, 18),
    ("A147565", 2, 6),
];

#[test]
fn issue_2513_reduction_is_an_exact_polynomial_identity() {
    for (s, kappa) in [(1, 5), (1, 8), (1, 10), (1, 18), (2, 6), (3, 7), (4, 13)] {
        let rows = orbit_rule(s)
            .generate_rows_rational(&orbit_seeds(kappa), 0, 13)
            .unwrap();
        let alpha = Q::new(s.into(), 2.into());
        let r = Q::new((kappa - 4).into(), (2 * s).into());
        assert!(r >= Q::new(1.into(), 2.into()));
        let mut p = vec![q(1)];
        for (n, row) in rows.iter().enumerate() {
            let transformed = vec_poly::compose(
                row,
                &[Q::new((-1).into(), 2.into()), Q::new(1.into(), 2.into())],
            );
            let expected = vec_poly::add(
                &vec_poly::scale(&p, &r),
                &vec_poly::scale(&vec_poly::pow(&[q(0), q(1)], n), &(Q::one() - &r)),
            );
            assert_eq!(
                vec_poly::trim(&transformed),
                vec_poly::trim(&expected),
                "s={s}, kappa={kappa}, n={n}"
            );
            p = vec_poly::add(
                &vec_poly::shift(&p),
                &vec_poly::scale(
                    &vec_poly::mul(&[q(-1), q(0), q(1)], &vec_poly::derivative(&p)),
                    &alpha,
                ),
            );
        }
    }
}

#[test]
fn reduction_recurrence_identity_holds_symbolically_for_every_index() {
    use ore_algebra::{BivariatePolynomial as B, OreOperator as O};
    // In the x coordinate, at target m=n+2, the recurrence operator is
    // [n S-(n+1)x] [S-L]. This annihilates the orbit p for every index.
    // The coefficients are affine in s; s=0 and s=1 verify that entire
    // parameter identity, without an index specialization.
    for s in [0, 1] {
        let map = RecurrenceTransformation::Affine {
            a: "1/2".into(),
            b: "-1/2".into(),
        };
        let TransportedRecurrence::Scalar(rule) = map.transport(&orbit_rule(s)).unwrap() else {
            panic!("scalar affine transport")
        };
        let mut actual = O::single(
            2,
            0,
            B::new(rule.denominator.unwrap().coeffs).shift_index(2),
        );
        for term in rule.terms {
            actual = actual.subtract(&O::single(
                2 - term.offset,
                term.deriv_order,
                B::new(term.coeff.coeffs).shift_index(2),
            ));
        }
        let alpha = Q::new(s.into(), 2.into());
        let mut orbit = O::single(1, 0, B::one());
        orbit = orbit.subtract(&O::single(0, 0, B::from_variable_polynomial(&[q(0), q(1)])));
        orbit = orbit.subtract(&O::single(
            0,
            1,
            B::from_variable_polynomial(&[-alpha.clone(), q(0), alpha.clone()]),
        ));
        let factor = O::single(1, 0, B::new(vec![vec![q(0)], vec![q(1)]])).subtract(&O::single(
            0,
            0,
            B::new(vec![vec![q(0), q(1)], vec![q(0), q(1)]]),
        ));
        assert_eq!(actual, factor.compose(&orbit));
        // For z_m=x^m, z_m-L z_(m-1) = -alpha*(m-1)*(x^2-1)*x^(m-2).
        // The remaining coefficient of (x^2-1)*x^n is identically zero:
        // n*(-alpha*(n+1))-(n+1)*(-alpha*n)=0 in Q[n].
        let n = vec![q(0), q(1)];
        let next = vec![q(1), q(1)];
        assert!(vec_poly::is_zero(&vec_poly::sub(
            &vec_poly::mul(&n, &vec_poly::scale(&next, &-alpha.clone())),
            &vec_poly::mul(&next, &vec_poly::scale(&n, &-alpha)),
        )));
    }
}

#[test]
fn catalog_rows_map_to_the_stated_parameters_exactly() {
    for (id, s, kappa) in TARGETS {
        let entry = polytool::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let rows = rule
            .generate_rows_rational(&seeds, entry.recurrence_first_index, 11)
            .unwrap();
        let first_degree = vec_poly::degree(&seeds[0]).unwrap();
        let map = RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule {
                first_index: entry.recurrence_first_index,
                first_degree,
                step: 1,
                last_index: None,
            },
        };
        let expected = orbit_rule(s)
            .generate_rows_rational(&orbit_seeds(kappa), 0, first_degree + rows.len())
            .unwrap();
        for (i, row) in rows.iter().enumerate() {
            let actual = map.map_row(entry.recurrence_first_index + i, row).unwrap();
            assert_eq!(
                vec_poly::trim(&actual[0]),
                vec_poly::trim(&expected[first_degree + i]),
                "{id}, degree {}",
                first_degree + i
            );
        }
    }
}

fn check(rule: &Recurrence, seeds: &[Vec<Q>]) -> RealRootednessTheoremReport {
    check_binomial_perturbed_euler_orbit_theorem(
        rule,
        0,
        seeds,
        &RealRootednessTheoremOptions {
            maximum_start_delay: 0,
            ..Default::default()
        },
    )
    .unwrap()
}

#[test]
fn direct_rules_catalog_rules_and_cli_prove_and_replay() {
    for (id, s, kappa) in TARGETS {
        let direct = check(&orbit_rule(s), &orbit_seeds(kappa));
        assert!(direct.is_proved(), "{id}: {:?}", direct.open_hypotheses());
        direct.verify().unwrap();
        let conclusion = direct.conclusion.as_ref().unwrap();
        assert_eq!(
            conclusion
                .root_region
                .as_ref()
                .unwrap()
                .at(0)
                .unwrap()
                .lower,
            Some(q(-1))
        );
        assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
        assert!(conclusion.simple_roots);
        assert_eq!(
            direct.lean_declarations[0],
            "RealRooted.strictInterl_of_orbitPerturbation_recurrence"
        );
        assert!(direct.reference.contains("RealRooted/issues/2513"));

        let entry = polytool::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let first = entry.recurrence_first_index;
        let report =
            check_binomial_perturbed_euler_orbit_theorem(&rule, first, &seeds, &Default::default())
                .unwrap();
        assert!(
            report.is_proved(),
            "{id}: {:?}; {}",
            report.open_hypotheses(),
            report.normal_form
        );
        report.verify().unwrap();
        let certificate = report.binomial_perturbed_euler_orbit_certificate().unwrap();
        assert_eq!(certificate.s, s.to_string());
        assert_eq!(certificate.kappa, kappa.to_string());
        assert!(certificate.coordinate.is_some());
        certificate.verify(&rule, first, &seeds).unwrap();
        let conclusion = report.conclusion.as_ref().unwrap();
        let region = conclusion.root_region.as_ref().unwrap().at(0).unwrap();
        assert!(region.lower.is_none());
        assert_eq!(region.upper, Some(q(0)));
        assert!(!region.upper_closed);
        assert_eq!(conclusion.interlacing_from, Some(first + 1));

        let record = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None);
        let output = common::run_polytool(
            &["recurrence-prove", "--json"],
            &serde_json::to_string(&record).unwrap(),
        );
        assert!(output.status.success(), "{id}: {}", common::stderr(&output));
        let discovery: RecurrenceProofDiscoveryReport =
            serde_json::from_slice(&output.stdout).unwrap();
        assert_eq!(discovery.status, "proved");
        assert!(discovery.attempts.iter().any(|r| r.theorem
            == RealRootednessTheorem::BinomialPerturbedEulerOrbit
            && r.is_proved()));
        discovery.verify().unwrap();
        let replay = common::run_polytool(
            &["recurrence-prove", "--verify", "--json"],
            &common::stdout(&output),
        );
        assert!(replay.status.success(), "{id}: {}", common::stderr(&replay));
        let parent = prove_parent_real_rootedness(
            &ParentRecurrenceRecord {
                rule: &rule,
                seeds: &seeds,
                first_index: first,
            },
            first,
            &[],
            true,
        )
        .unwrap();
        assert_eq!(
            parent.route,
            Some(VeroneseParentProofRoute::BinomialPerturbedEulerOrbit)
        );
    }
}

#[test]
fn exact_parameter_failures_bad_seeds_and_perturbed_operators_are_rejected() {
    for (s, kappa, name) in [
        (1, 4, "binomial_parameter_bound"),
        (2, 5, "binomial_parameter_bound"),
        (0, 5, "euler_parameter_at_least_one"),
    ] {
        let report = check(&orbit_rule(s), &orbit_seeds(kappa));
        assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
        assert!(
            report
                .hypotheses
                .iter()
                .any(|h| h.name == name && matches!(h.status, LiftHypothesisStatus::Failed { .. })),
            "{:?}",
            report.open_hypotheses()
        );
        report.verify().unwrap();
    }
    let mut wrong = orbit_rule(1);
    wrong
        .terms
        .push(RecurrenceTerm::new(2, 1, bivar(&[&[0, 1]])));
    let report = check(&wrong, &orbit_seeds(5));
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    report.verify().unwrap();
    let mut late = orbit_rule(1);
    let invisible_on_prefix = (3..=15).fold(vec![q(1)], |p, n| vec_poly::mul(&p, &[q(-n), q(1)]));
    late.terms.push(RecurrenceTerm::new(
        2,
        1,
        BivarPoly::from_coefficients_normalized(
            invisible_on_prefix
                .into_iter()
                .map(|v| vec![q(0), v])
                .collect(),
        ),
    ));
    assert_eq!(
        late.generate_rows_rational(&orbit_seeds(5), 0, 16).unwrap(),
        orbit_rule(1)
            .generate_rows_rational(&orbit_seeds(5), 0, 16)
            .unwrap(),
    );
    assert_eq!(
        check(&late, &orbit_seeds(5)).status,
        RealRootednessTheoremStatus::NotApplicable
    );
    for seed_index in 0..3 {
        let mut wrong = orbit_seeds(5);
        wrong[seed_index][0] += q(1);
        let report = check(&orbit_rule(1), &wrong);
        assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
        report.verify().unwrap();
    }
    let report = check(&orbit_rule(1), &orbit_seeds(5)[..2]);
    assert!(!report.is_proved());
    report.verify().unwrap();
}

#[test]
fn certificate_tampering_cannot_forge_parameters_maps_rules_seeds_or_conclusions() {
    let entry = polytool::oeis::by_id("A141689").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let first = entry.recurrence_first_index;
    let report =
        check_binomial_perturbed_euler_orbit_theorem(&rule, first, &seeds, &Default::default())
            .unwrap();
    let certificate = report.binomial_perturbed_euler_orbit_certificate().unwrap();
    for field in [
        "s",
        "kappa",
        "degree_offset",
        "coordinate",
        "transformed_recurrence",
    ] {
        let mut value = serde_json::to_value(certificate).unwrap();
        match field {
            "s" | "kappa" => value[field] = serde_json::json!("100"),
            "degree_offset" => value[field] = serde_json::json!(0),
            "coordinate" => value[field] = serde_json::Value::Null,
            _ => value[field]["initial_polynomials"][0][0] = serde_json::json!("2"),
        }
        let forged: BinomialPerturbedEulerOrbitCertificate = serde_json::from_value(value).unwrap();
        assert!(forged.verify(&rule, first, &seeds).is_err(), "{field}");
        let mut forged_report = report.clone();
        forged_report.certificate = Some(RouteCertificate::BinomialPerturbedEulerOrbit(forged));
        assert!(forged_report.verify().is_err(), "{field}");
    }
    let mut forged = report.clone();
    forged.certificate = None;
    assert!(forged.verify().is_err());
    let mut operator = serde_json::to_value(certificate).unwrap();
    operator["transformed_recurrence"]["recurrence"]["terms"][0]["coeff"]["coeffs"][0][0] =
        serde_json::json!("314");
    let operator: BinomialPerturbedEulerOrbitCertificate =
        serde_json::from_value(operator).unwrap();
    assert!(operator.verify(&rule, first, &seeds).is_err());
    let mut forged = report;
    forged.conclusion.as_mut().unwrap().interlacing_from = Some(first);
    assert!(forged.verify().is_err());
}

#[test]
fn rational_parameters_common_factors_shifted_starts_and_region_requests() {
    let mut rule = orbit_rule(3);
    // s=3/2, kappa=11/2 is the exact lower boundary.
    for term in rule.terms.iter_mut().filter(|t| t.deriv_order == 1) {
        term.coeff = BivarPoly::from_coefficients_normalized(
            term.coeff
                .coeffs
                .iter()
                .map(|r| vec_poly::scale(r, &Q::new(1.into(), 2.into())))
                .collect(),
        );
    }
    let mut seeds = orbit_seeds(11);
    seeds[2][1] /= q(2);
    seeds[2][2] /= q(2);
    let report = check(&rule, &seeds);
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    report.verify().unwrap();
    let rows = rule.generate_rows_rational(&seeds, 0, 9).unwrap();
    let report =
        check_binomial_perturbed_euler_orbit_theorem(&rule, 7, &rows[7..], &Default::default())
            .unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    report.verify().unwrap();
    for lower in [-2, -1, 0] {
        let report = check_binomial_perturbed_euler_orbit_theorem(
            &rule,
            0,
            &seeds,
            &RealRootednessTheoremOptions {
                root_region: Some(polytool::real_rootedness::RealInterval {
                    lower: Some(q(lower)),
                    lower_closed: false,
                    upper: Some(q(0)),
                    upper_closed: false,
                }),
                maximum_start_delay: 0,
            },
        )
        .unwrap();
        assert_eq!(report.is_proved(), lower <= -1);
        report.verify().unwrap();
    }
    // A common factor preserves the rational operator, but its zero at 100
    // still invalidates the defining rule, even though all early rows agree.
    for coefficient in rule
        .terms
        .iter_mut()
        .map(|t| &mut t.coeff)
        .chain(rule.denominator.iter_mut())
    {
        let source = coefficient.coeffs.clone();
        let width = source.iter().map(Vec::len).max().unwrap();
        let mut product = vec![vec![q(0); width]; source.len() + 1];
        for (i, row) in source.iter().enumerate() {
            for (j, v) in row.iter().enumerate() {
                product[i][j] -= q(100) * v;
                product[i + 1][j] += v;
            }
        }
        *coefficient = BivarPoly::from_coefficients_normalized(product);
    }
    let report = check(&rule, &seeds);
    assert!(!report.is_proved());
    assert!(report
        .hypotheses
        .iter()
        .any(|h| h.name == "left_factor_nonzero"
            && matches!(h.status, LiftHypothesisStatus::Failed { .. })));
    report.verify().unwrap();
}

#[test]
fn finite_exact_root_sanity_through_degree_sixty_is_not_the_proof() {
    for (id, s, kappa) in TARGETS {
        let entry = polytool::oeis::by_id(id).unwrap();
        // Include the authoritative prefix rows of the three catalog rules
        // whose recurrence seeds start at degree two: 61 checks per family.
        let rows = entry.generate_rows(61).unwrap();
        let expected = orbit_rule(s)
            .generate_rows_rational(&orbit_seeds(kappa), 0, 61)
            .unwrap();
        let coordinate = RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(0, 0, 1),
        };
        for (m, integer) in rows.iter().enumerate() {
            let row: Vec<Q> = integer.iter().cloned().map(Q::from_integer).collect();
            assert!(integer.iter().all(|v| v > &BigInt::from(0)));
            assert!(
                polytool::is_real_rooted_bigint_coeffs(integer),
                "{id}, m={m}"
            );
            // Positive coefficients and nonzero constant exclude [0,inf).
            // The exact inverse Mobius map carries every negative root to
            // (-1,0), and the degree and endpoint checks exclude added roots.
            let mapped = coordinate.map_row(m, &row).unwrap().remove(0);
            assert_eq!(mapped, expected[m], "{id}, m={m}");
            assert_eq!(vec_poly::degree(&mapped), Some(m));
            assert_eq!(vec_poly::evaluate(&mapped, &q(0)), q(1));
            assert_eq!(
                vec_poly::evaluate(&mapped, &q(-1)),
                q(if m % 2 == 0 { 1 } else { -1 })
            );
        }
    }
}
