//! Tests of the coordinate-transformation registry and recurrence transport.

use super::*;

fn r(n: i64, d: i64) -> Q {
    Q::new(BigInt::from(n), BigInt::from(d))
}

fn row(coefficients: &[i64]) -> Vec<Q> {
    vec_poly::trim(&coefficients.iter().map(|&c| q(c)).collect::<Vec<_>>())
}

/// `rows[i][j]` is the coefficient of `n^i x^j`.
fn bivar(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&c| q(c)).collect())
            .collect(),
    )
}

fn term(offset: usize, derivative: usize, coefficient: BivarPoly) -> RecurrenceTerm {
    RecurrenceTerm::new(offset, derivative, coefficient)
}

/// Canonical A008292 with `E_1 = 1`: `E_m = (1 + (m-1)x) E_(m-1) + x(1-x) E'_(m-1)`.
fn eulerian_rule() -> Recurrence {
    Recurrence {
        terms: vec![
            term(1, 0, bivar(&[&[1, -1], &[0, 1]])),
            term(1, 1, bivar(&[&[0, 1, -1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    }
}

fn scalar(transported: TransportedRecurrence) -> Recurrence {
    match transported {
        TransportedRecurrence::Scalar(rule) => rule,
        TransportedRecurrence::Vector(_) => panic!("expected a scalar rule"),
    }
}

fn generate(rule: &Recurrence, origin: usize, seeds: &[Vec<Q>], count: usize) -> Vec<Vec<Q>> {
    rule.generate_rows_rational(seeds, origin, count)
        .unwrap()
        .into_iter()
        .map(|r| vec_poly::trim(&r))
        .collect()
}

/// Map rows, transport the rule and check the transported identity on every
/// target, together with exact reconstruction of every row.
fn assert_transport_holds(
    transformation: &RecurrenceTransformation,
    rule: &Recurrence,
    origin: usize,
    rows: &[Vec<Q>],
    first_target: usize,
) -> Recurrence {
    let mapped = rows
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let state = transformation.map_row(origin + i, p).unwrap();
            assert_eq!(transformation.inverse_row(origin + i, &state).unwrap(), *p);
            state[0].clone()
        })
        .collect::<Vec<_>>();
    let transported = scalar(transformation.transport(rule).unwrap());
    for target in first_target..origin + rows.len() {
        assert_eq!(
            transformed_rule_holds_at(&transported, &mapped, origin, target),
            Some(true),
            "{transformation:?} fails at {target}: {transported}"
        );
    }
    // A perturbed row breaks the identity: the check is not vacuous.
    let mut broken = mapped.clone();
    let last = broken.len() - 1;
    broken[last] = vec_poly::add(&broken[last], &[q(1)]);
    assert_eq!(
        transformed_rule_holds_at(&transported, &broken, origin, origin + last),
        Some(false)
    );
    transported
}

#[test]
fn eulerian_f_polynomial_transport_matches_the_issue_example() {
    let transformation = RecurrenceTransformation::FPolynomial {
        schedule: DegreeSchedule::new(1, 0, 1),
    };
    let transported = scalar(transformation.transport(&eulerian_rule()).unwrap());
    let expected = Recurrence {
        terms: vec![
            term(1, 0, bivar(&[&[1, 2]])),
            term(1, 1, bivar(&[&[0, 1, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    assert!(transported.denominator.is_none());
    assert!(
        same_canonical_recurrence(&transported, &expected),
        "{transported}"
    );
    assert_eq!(
        transported.to_string(),
        "P(n) = (1 + 2t) P(n-1) + (t + t^2) P'(n-1)"
    );

    let rows = generate(&eulerian_rule(), 1, &[row(&[1])], 8);
    assert_eq!(rows[2], row(&[1, 4, 1]));
    let mapped = assert_transport_holds(&transformation, &eulerian_rule(), 1, &rows, 2);
    let q_rows = generate(&mapped, 1, &[row(&[1])], 4);
    assert_eq!(
        q_rows,
        vec![
            row(&[1]),
            row(&[1, 2]),
            row(&[1, 6, 6]),
            row(&[1, 14, 36, 24])
        ]
    );
}

#[test]
fn oeis_request_retains_transport_and_independent_search() {
    let request = RecurrenceTransformRequest {
        transformation: RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(1, 0, 1),
        },
        source_recurrence: None,
        oeis_id: Some("A008292".into()),
        rows: None,
        row_origin: None,
        row_count: Some(9),
        search: Some(TransformSearchBounds::default()),
        limits: TransformLimits::default(),
    };
    let report = transform_recurrence(&request).unwrap();
    assert!(report.reconstruction_verified);
    assert_eq!(report.row_origin, 1);
    assert_eq!(report.transformed_states[1], vec![vec!["1", "2"]]);

    let transported = report.transported.as_ref().unwrap();
    assert_eq!(transported.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert!(transported.holds_on_transformed_rows);
    assert_eq!(transported.valid_from_index, 2);
    assert_eq!(transported.checked_targets, 8);
    let scalar_json = transported.scalar.as_ref().unwrap();
    assert_eq!(scalar_json.first_index, 1);
    assert_eq!(scalar_json.initial_polynomials, vec![vec!["1"]]);
    let (rule, first, seeds) = scalar_json.to_recurrence_parts().unwrap();
    assert_eq!(
        generate(&rule, first, &seeds, 3),
        vec![row(&[1]), row(&[1, 2]), row(&[1, 6, 6])]
    );

    let search = report.independent_search.as_ref().unwrap();
    assert_eq!(search.evidence, FINITE_FIT_EVIDENCE);
    assert_eq!(search.status, "found");
    assert_eq!(search.holds_on_transformed_rows, Some(true));
    assert!(search.scalar.as_ref().unwrap().search.is_some());
    assert_ne!(search.evidence, transported.evidence);

    assert_eq!(report.conditions.len(), 1);
    assert_eq!(report.conditions[0].status, CERTIFIED_ALL_INDICES);
}

#[test]
fn report_round_trips_through_json_and_replays() {
    let request = RecurrenceTransformRequest::from_source(
        RecurrenceTransformation::SymmetricDecomposition {
            schedule: DegreeSchedule::new(1, 0, 1),
        },
        RecurrenceJson::from_recurrence_rational(&eulerian_rule(), 1, &[row(&[1])], None),
        7,
    );
    let report = transform_recurrence(&request).unwrap();
    let encoded = serde_json::to_string(&report).unwrap();
    let decoded: RecurrenceTransformReport = serde_json::from_str(&encoded).unwrap();
    assert_eq!(decoded, report);
    verify_recurrence_transform(&decoded).unwrap();

    let vector = report
        .transported
        .as_ref()
        .unwrap()
        .vector
        .as_ref()
        .unwrap();
    let (system, first, state) = vector.to_recurrence().unwrap();
    assert_eq!(first, 1);
    let next = system.evaluate_next_rational(&state, first).unwrap();
    assert_eq!(next, vec![row(&[1, 1]), row(&[0])]);

    let mut tampered = decoded.clone();
    tampered.transformed_states[2][0][0] = "2".into();
    assert!(verify_recurrence_transform(&tampered).is_err());
    let mut wrong_origin = decoded;
    wrong_origin.request.row_origin = Some(1);
    wrong_origin.row_origin = 0;
    assert!(verify_recurrence_transform(&wrong_origin).is_err());
}

#[test]
fn map_and_inverse_are_exact_on_basis_and_boundary_rows() {
    let schedule = DegreeSchedule::new(0, 4, 0);
    let transformations = vec![
        RecurrenceTransformation::MultiplyPower { power: 2 },
        RecurrenceTransformation::Reversal {
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::Affine {
            a: "-3/2".into(),
            b: "5".into(),
        },
        RecurrenceTransformation::Mobius {
            a: "2".into(),
            b: "-1/3".into(),
            c: "3/2".into(),
            e: "5".into(),
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::Mobius {
            a: "2".into(),
            b: "1".into(),
            c: "0".into(),
            e: "3".into(),
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::FPolynomial {
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::SymmetricDecomposition {
            schedule: schedule.clone(),
        },
    ];
    let mut rows = vec![row(&[]), row(&[0, 0, 0, 0, 0, 0])];
    for k in 0..=4 {
        let mut basis = vec![0; k + 1];
        basis[k] = 1;
        rows.push(row(&basis));
    }
    // Missing top coefficient, endpoint factors x^2 and (3x - 4) = (c y - a) * 2.
    rows.push(vec![r(1, 7), q(-2), q(0), r(5, 3)]);
    rows.push(row(&[0, 0, 1, 1]));
    rows.push(row(&[-4, 3]));
    rows.push(vec_poly::mul(&row(&[-4, 3]), &row(&[0, 1, 1])));
    for transformation in &transformations {
        for p in &rows {
            let state = transformation.map_row(0, p).unwrap();
            assert_eq!(state.len(), transformation.component_count());
            assert_eq!(
                transformation.inverse_row(0, &state).unwrap(),
                vec_poly::trim(p),
                "{transformation:?} on {p:?}"
            );
        }
    }
    // Endpoint factor (3x - 4) of P lowers the degree of the Möbius image.
    let state = transformations[3].map_row(0, &row(&[-4, 3])).unwrap();
    assert_eq!(vec_poly::degree(&state[0]), Some(3));
    // Reversal of a polynomial with a missing top coefficient has an x factor.
    let reversed = transformations[1].map_row(0, &row(&[1, 1])).unwrap();
    assert_eq!(reversed[0], row(&[0, 0, 0, 1, 1]));
    // Divisible rows round-trip through removal of x^2.
    let divide = RecurrenceTransformation::DividePower { power: 2 };
    for p in [row(&[]), row(&[0, 0, 1]), row(&[0, 0, 0, 2, 5])] {
        let state = divide.map_row(3, &p).unwrap();
        assert_eq!(divide.inverse_row(3, &state).unwrap(), p);
    }
}

#[test]
fn gamma_coordinates_handle_endpoint_factors_and_rationals() {
    let gamma = RecurrenceTransformation::Gamma {
        schedule: DegreeSchedule::new(0, 5, 0),
    };
    // x (1 + x)^3 + 2 x^2 (1 + x) is symmetric at index 5 but has degree 4.
    let p = vec_poly::add(
        &vec_poly::shift_by(&vec_poly::pow(&row(&[1, 1]), 3), 1),
        &vec_poly::scale(&vec_poly::shift_by(&row(&[1, 1]), 2), &q(2)),
    );
    let state = gamma.map_row(0, &p).unwrap();
    assert_eq!(state[0], row(&[0, 1, 2]));
    assert_eq!(gamma.inverse_row(0, &state).unwrap(), p);
    let half = vec_poly::scale(&p, &r(1, 2));
    assert_eq!(
        gamma.map_row(0, &half).unwrap()[0],
        vec![q(0), r(1, 2), q(1)]
    );
    assert_eq!(gamma.map_row(0, &row(&[])).unwrap()[0], row(&[]));
    for basis in 0..=2 {
        let mut g = vec![0; basis + 1];
        g[basis] = 1;
        let p = gamma.inverse_row(0, &[row(&g)]).unwrap();
        assert_eq!(gamma.map_row(0, &p).unwrap()[0], row(&g));
    }
}

#[test]
fn transported_residuals_hold_for_scalar_maps() {
    // P_m = (1 + m x) P_(m-1) + (2 - x^2) P'_(m-1) - m x^2 P''_(m-2), degree grows by one.
    let rule = Recurrence {
        terms: vec![
            term(1, 0, bivar(&[&[1], &[0, 1]])),
            term(1, 1, bivar(&[&[2, 0, -1]])),
            term(2, 2, bivar(&[&[0], &[0, 0, -1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let rows = generate(&rule, 2, &[row(&[1, 3]), row(&[2, -1, 1])], 9);
    let schedule = DegreeSchedule::new(2, 2, 1);
    let maps = vec![
        RecurrenceTransformation::MultiplyPower { power: 3 },
        RecurrenceTransformation::Reversal {
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::Affine {
            a: "2/3".into(),
            b: "-1".into(),
        },
        RecurrenceTransformation::Mobius {
            a: "1".into(),
            b: "2".into(),
            c: "-1/2".into(),
            e: "3".into(),
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::Mobius {
            a: "1".into(),
            b: "2".into(),
            c: "0".into(),
            e: "3".into(),
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(2, 5, 2),
        },
    ];
    for transformation in &maps {
        assert_transport_holds(transformation, &rule, 2, &rows, 4);
    }
}

#[test]
fn transported_residuals_hold_with_denominators_signs_and_forcing() {
    // (m + 1) P_m = (-1)^m x P'_(m-1) + (m + x^2) P_(m-1) + (1 + m x), degree <= 3.
    let rule = Recurrence {
        terms: vec![
            term(1, 1, bivar(&[&[0, 1]])).with_sign(RecurrenceSign::AlternatingN),
            term(1, 0, bivar(&[&[0, 0, 1], &[1]])),
        ],
        denominator: Some(bivar(&[&[1], &[1]])),
        inhomogeneous: Some(bivar(&[&[1], &[0, 1]])),
    };
    let seeds = [row(&[1])];
    let mut rows = vec![seeds[0].clone()];
    for m in 1..6 {
        rows.push(rule.evaluate_next_rational(&rows, 0).unwrap());
        assert!(vec_poly::degree(&rows[m]).unwrap_or(0) <= 2 * m);
    }
    let rows = rows
        .into_iter()
        .map(|p| vec_poly::trim(&p))
        .collect::<Vec<_>>();
    let constant = DegreeSchedule::new(0, 10, 0);
    for transformation in [
        RecurrenceTransformation::Reversal {
            schedule: constant.clone(),
        },
        RecurrenceTransformation::Mobius {
            a: "3".into(),
            b: "1".into(),
            c: "2".into(),
            e: "-1".into(),
            schedule: constant.clone(),
        },
        RecurrenceTransformation::MultiplyPower { power: 1 },
        RecurrenceTransformation::Affine {
            a: "1".into(),
            b: "1".into(),
        },
    ] {
        let transported = assert_transport_holds(&transformation, &rule, 0, &rows, 1);
        assert!(transported.inhomogeneous.is_some());
    }
    let growing = RecurrenceTransformation::Mobius {
        a: "3".into(),
        b: "1".into(),
        c: "2".into(),
        e: "-1".into(),
        schedule: DegreeSchedule::new(0, 10, 2),
    };
    assert!(matches!(
        growing.transport(&rule),
        Err(TransformError::InhomogeneousNotTransportable)
    ));
}

#[test]
fn zero_rows_and_divisibility_transport() {
    // P_m = x P'_(m-1) + m P_(m-1) keeps x^2 | P_m.
    let rule = Recurrence {
        terms: vec![
            term(1, 1, bivar(&[&[0, 1]])),
            term(1, 0, bivar(&[&[0], &[1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let rows = generate(&rule, 1, &[row(&[0, 0, 1, 1])], 6);
    let divide = RecurrenceTransformation::DividePower { power: 2 };
    let transported = assert_transport_holds(&divide, &rule, 1, &rows, 2);
    assert!(transported.denominator.is_none());
    let zeros = generate(&rule, 1, &[row(&[])], 4);
    assert!(zeros.iter().all(|p| vec_poly::is_zero(p)));
    for transformation in [
        divide,
        RecurrenceTransformation::FPolynomial {
            schedule: DegreeSchedule::new(1, 0, 1),
        },
    ] {
        let transported = scalar(transformation.transport(&rule).unwrap());
        for (i, p) in zeros.iter().enumerate() {
            let state = transformation.map_row(1 + i, p).unwrap();
            assert!(vec_poly::is_zero(&state[0]));
            assert!(vec_poly::is_zero(
                &transformation.inverse_row(1 + i, &state).unwrap()
            ));
        }
        let mapped = zeros.clone();
        for target in 2..5 {
            assert_eq!(
                transformed_rule_holds_at(&transported, &mapped, 1, target),
                Some(true)
            );
        }
    }
}

#[test]
fn involutions_and_inverse_maps_transport_back() {
    // A non-symmetric rule: P_m = (1 + 2x) P_(m-1) + x(1 + x) P'_(m-1).
    let rule = Recurrence {
        terms: vec![
            term(1, 0, bivar(&[&[1, 2]])),
            term(1, 1, bivar(&[&[0, 1, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let schedule = DegreeSchedule::new(1, 0, 1);
    let reversal = RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    };
    let once = scalar(reversal.transport(&rule).unwrap());
    assert!(!same_canonical_recurrence(&once, &rule));
    let twice = scalar(reversal.transport(&once).unwrap());
    assert!(same_canonical_recurrence(&twice, &rule), "{twice}");
    // A different degree bound gives a different rule.
    let shifted = RecurrenceTransformation::Reversal {
        schedule: DegreeSchedule::new(1, 1, 1),
    };
    assert!(!same_canonical_recurrence(
        &scalar(shifted.transport(&rule).unwrap()),
        &once
    ));

    let f = RecurrenceTransformation::FPolynomial {
        schedule: schedule.clone(),
    };
    let inverse_f = RecurrenceTransformation::Mobius {
        a: "1".into(),
        b: "0".into(),
        c: "-1".into(),
        e: "1".into(),
        schedule,
    };
    let forward = scalar(f.transport(&eulerian_rule()).unwrap());
    let back = scalar(inverse_f.transport(&forward).unwrap());
    assert!(same_canonical_recurrence(&back, &eulerian_rule()), "{back}");
}

#[test]
fn symmetric_decomposition_transports_to_the_complete_vector_system() {
    let schedule = DegreeSchedule::new(1, 0, 1);
    let transformation = RecurrenceTransformation::SymmetricDecomposition {
        schedule: schedule.clone(),
    };
    let rule = Recurrence {
        terms: vec![
            term(1, 0, bivar(&[&[1, 2]])),
            term(1, 1, bivar(&[&[0, 1, 1]])),
            term(1, 2, bivar(&[&[0], &[0, 0, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    for rule in [rule, eulerian_rule()] {
        let rows = generate(&rule, 1, &[row(&[1])], 7);
        let wide = RecurrenceTransformation::SymmetricDecomposition {
            schedule: DegreeSchedule::new(1, 0, 2),
        };
        for transformation in [&transformation, &wide] {
            let states = rows
                .iter()
                .enumerate()
                .map(|(i, p)| transformation.map_row(1 + i, p).unwrap())
                .collect::<Vec<_>>();
            let TransportedRecurrence::Vector(system) = transformation.transport(&rule).unwrap()
            else {
                panic!("expected a vector system");
            };
            assert_eq!(system.dimension(), 2);
            for n in 1..rows.len() {
                assert!(system.holds_transition_rational(&states[n - 1], &states[n], n));
            }
        }
    }
}

#[test]
fn gamma_symmetry_is_certified_from_the_source_rule() {
    let source = RecurrenceJson::from_recurrence_rational(&eulerian_rule(), 1, &[row(&[1])], None);
    let request = RecurrenceTransformRequest::from_source(
        RecurrenceTransformation::Gamma {
            schedule: DegreeSchedule::new(1, 0, 1),
        },
        source,
        8,
    );
    let report = transform_recurrence(&request).unwrap();
    assert_eq!(report.conditions[0].status, CERTIFIED_ALL_INDICES);
    assert_eq!(report.transformed_states[3][0], vec!["1", "8"]);
    assert!(report.transported.is_none());
    assert!(report
        .transport_error
        .as_ref()
        .unwrap()
        .contains("gamma rows directly"));

    let rows_only = RecurrenceTransformRequest {
        source_recurrence: None,
        rows: Some(report.source_rows.clone()),
        row_origin: Some(1),
        row_count: None,
        search: Some(TransformSearchBounds::default()),
        ..request
    };
    let finite = transform_recurrence(&rows_only).unwrap();
    assert_eq!(finite.conditions[0].status, CHECKED_SUPPLIED_ROWS);
    assert!(finite.transported.is_none() && finite.transport_error.is_none());
    assert_eq!(finite.independent_search.as_ref().unwrap().status, "found");
}

#[test]
fn invalid_schedules_and_maps_are_rejected() {
    let invalid = [
        DegreeSchedule::new(3, 2, -1),
        DegreeSchedule::new(3, 2, 1).with_last_index(2),
        DegreeSchedule::new(0, 2, -1).with_last_index(3),
    ];
    for schedule in invalid {
        assert!(matches!(
            schedule.validate(),
            Err(TransformError::InvalidSchedule(_))
        ));
        assert!(RecurrenceTransformation::Reversal { schedule }
            .map_row(0, &row(&[1]))
            .is_err());
    }
    assert!(DegreeSchedule::new(0, 2, -1)
        .with_last_index(2)
        .validate()
        .is_ok());
    let singular = RecurrenceTransformation::Mobius {
        a: "1".into(),
        b: "2".into(),
        c: "2".into(),
        e: "4".into(),
        schedule: DegreeSchedule::new(0, 2, 0),
    };
    assert_eq!(singular.validate(), Err(TransformError::SingularMobius));
    assert!(matches!(
        singular.transport(&eulerian_rule()),
        Err(TransformError::SingularMobius)
    ));
    let flat = RecurrenceTransformation::Affine {
        a: "0".into(),
        b: "1".into(),
    };
    assert_eq!(flat.validate(), Err(TransformError::SingularMobius));
    let malformed = RecurrenceTransformation::Affine {
        a: "1/0".into(),
        b: "1".into(),
    };
    assert!(matches!(
        malformed.validate(),
        Err(TransformError::InvalidParameter(_))
    ));
}

#[test]
fn uncertified_preconditions_and_images_are_rejected() {
    let schedule = DegreeSchedule::new(1, 2, 1);
    let gamma = RecurrenceTransformation::Gamma {
        schedule: schedule.clone(),
    };
    assert_eq!(
        gamma.map_row(1, &row(&[1, 2, 3])),
        Err(TransformError::NotSymmetric {
            index: 1,
            symmetry_index: 2
        })
    );
    // Symmetric at the actual degree 2 but not at the declared index 3.
    assert!(matches!(
        gamma.map_row(2, &row(&[1, 2, 1])),
        Err(TransformError::NotSymmetric { .. })
    ));
    assert!(matches!(
        gamma.inverse_row(1, &[row(&[1, 1, 1])]),
        Err(TransformError::ImageCondition { .. })
    ));
    assert_eq!(
        RecurrenceTransformation::DividePower { power: 1 }.map_row(4, &row(&[1, 1])),
        Err(TransformError::NotDivisible { index: 4, power: 1 })
    );
    assert!(matches!(
        RecurrenceTransformation::MultiplyPower { power: 2 }.inverse_row(0, &[row(&[0, 1])]),
        Err(TransformError::ImageCondition { .. })
    ));
    let reversal = RecurrenceTransformation::Reversal {
        schedule: schedule.clone(),
    };
    assert_eq!(
        reversal.map_row(1, &row(&[1, 0, 0, 1])),
        Err(TransformError::DegreeExceedsBound {
            index: 1,
            degree: 3,
            bound: 2
        })
    );
    assert_eq!(
        reversal.map_row(0, &row(&[1])),
        Err(TransformError::IndexOutsideDomain { index: 0 })
    );
    assert!(matches!(
        reversal.inverse_row(1, &[row(&[1, 0, 0, 1])]),
        Err(TransformError::ImageCondition { .. })
    ));
    let pair = RecurrenceTransformation::SymmetricDecomposition { schedule };
    assert!(matches!(
        pair.inverse_row(1, &[row(&[1, 1]), row(&[1, 1])]),
        Err(TransformError::ImageCondition { .. })
    ));
    assert!(matches!(
        pair.inverse_row(1, &[row(&[1, 0, 1]), row(&[1])]),
        Err(TransformError::ImageCondition { .. })
    ));
    assert!(matches!(
        pair.inverse_row(1, &[row(&[1])]),
        Err(TransformError::ComponentCount { .. })
    ));
    let second_order = Recurrence {
        terms: vec![term(2, 0, bivar(&[&[1]]))],
        denominator: None,
        inhomogeneous: None,
    };
    assert!(matches!(
        pair.transport(&second_order),
        Err(TransformError::TransportUnsupported(_))
    ));

    // A request with a row violating the precondition is rejected.
    let request = RecurrenceTransformRequest {
        transformation: RecurrenceTransformation::DividePower { power: 1 },
        source_recurrence: None,
        oeis_id: None,
        rows: Some(vec![vec!["0".into(), "1".into()], vec!["1".into()]]),
        row_origin: None,
        row_count: None,
        search: None,
        limits: TransformLimits::default(),
    };
    assert!(transform_recurrence(&request)
        .unwrap_err()
        .contains("does not divide row 1"));
}

#[test]
fn transformation_json_is_tagged_and_strict() {
    let json = serde_json::json!({
        "kind": "mobius", "a": "1", "b": "0", "c": "1", "e": "1",
        "schedule": {"first_index": 1, "first_degree": 0, "step": 1}
    });
    let parsed: RecurrenceTransformation = serde_json::from_value(json).unwrap();
    assert!(parsed.validate().is_ok());
    assert_eq!(parsed.record().components.len(), 1);
    let unknown = serde_json::json!({"kind": "multiply_power", "power": 1, "extra": 2});
    assert!(serde_json::from_value::<RecurrenceTransformation>(unknown).is_err());
    let record = RecurrenceTransformation::SymmetricDecomposition {
        schedule: DegreeSchedule::new(1, 0, 1),
    }
    .record();
    assert_eq!(
        record.components[1].symmetry_index.as_deref(),
        Some("d_n - 1")
    );
    assert_eq!(record.index_domain, "n >= 1");
}
