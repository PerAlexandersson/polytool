//! Tests of componentwise maps and vector-system transport.

use super::*;

fn eulerian_rows(count: usize) -> (Recurrence, Vec<Vec<Q>>) {
    let (rule, seeds, _) = crate::oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    let rows = rule
        .generate_rows_rational(&seeds, 1, count)
        .unwrap()
        .into_iter()
        .map(|r| vec_poly::trim(&r))
        .collect();
    (rule, rows)
}

/// The symmetric decomposition of the Eulerian rows at `d_n = n` (index 1 on).
fn decomposition() -> RecurrenceTransformation {
    RecurrenceTransformation::SymmetricDecomposition {
        schedule: DegreeSchedule::new(1, 1, 1),
    }
}

fn componentwise_maps() -> Vec<RecurrenceTransformation> {
    let schedule = DegreeSchedule::new(1, 1, 1);
    vec![
        RecurrenceTransformation::Reversal {
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::FPolynomial {
            schedule: schedule.clone(),
        },
        RecurrenceTransformation::Mobius {
            a: "1".into(),
            b: "2".into(),
            c: "1".into(),
            e: "1".into(),
            schedule: DegreeSchedule::new(1, 2, 1),
        },
        RecurrenceTransformation::Affine {
            a: "2".into(),
            b: "1".into(),
        },
        RecurrenceTransformation::MultiplyPower { power: 1 },
    ]
}

#[test]
fn vector_systems_are_transported_through_componentwise_maps() {
    let (rule, rows) = eulerian_rows(12);
    let split = decomposition();
    let TransportedRecurrence::Vector(system) = split.transport(&rule).unwrap() else {
        panic!("vector system");
    };
    let offsets = split.component_degree_offsets();
    assert_eq!(offsets, vec![0, -1]);
    let states = rows
        .iter()
        .enumerate()
        .map(|(i, row)| {
            split
                .map_state(1 + i, std::slice::from_ref(row), &[0])
                .unwrap()
        })
        .collect::<Vec<_>>();
    for (i, pair) in states.windows(2).enumerate() {
        assert!(system.holds_transition_rational(&pair[0], &pair[1], 1 + i));
    }
    for map in componentwise_maps() {
        let mapped = states
            .iter()
            .enumerate()
            .map(|(i, state)| map.map_state(1 + i, state, &offsets).unwrap())
            .collect::<Vec<_>>();
        // The inverse is exact, componentwise.
        for (i, (state, image)) in states.iter().zip(&mapped).enumerate() {
            assert_eq!(&map.inverse_state(1 + i, image, &offsets).unwrap(), state);
        }
        let transported = transport_vector_componentwise(&system, &map, &offsets).unwrap();
        for (i, pair) in mapped.windows(2).enumerate() {
            assert!(
                transported.holds_transition_rational(&pair[0], &pair[1], 1 + i),
                "{} at n = {}",
                map.kind_name(),
                1 + i
            );
        }
        // The transported system regenerates the mapped states.
        let mut state = mapped[0].clone();
        for (i, expected) in mapped.iter().enumerate().skip(1) {
            state = transported.evaluate_next_rational(&state, i).unwrap();
            assert_eq!(&state, expected, "{}", map.kind_name());
        }
    }
}

#[test]
fn degree_offsets_shift_schedules_per_component() {
    let gamma = RecurrenceTransformation::Gamma {
        schedule: DegreeSchedule::new(1, 1, 1),
    };
    // The gamma vectors of a_n (index n) and b_n (index n - 1).
    let (_, rows) = eulerian_rows(6);
    let split = decomposition();
    let state = split
        .map_state(4, std::slice::from_ref(&rows[3]), &[0])
        .unwrap();
    let image = gamma.map_state(4, &state, &[0, -1]).unwrap();
    assert_eq!(image.len(), 2);
    assert_eq!(gamma.inverse_state(4, &image, &[0, -1]).unwrap(), state);
    // Without the offset b_n is not symmetric at index n.
    assert!(gamma.map_state(4, &state, &[0, 0]).is_err());
    assert_eq!(
        gamma.with_degree_offset(-1).unwrap(),
        RecurrenceTransformation::Gamma {
            schedule: DegreeSchedule::new(1, 0, 1)
        }
    );
    assert!(gamma.with_degree_offset(-2).is_err());
    assert_eq!(
        RecurrenceTransformation::MultiplyPower { power: 1 }
            .with_degree_offset(-5)
            .unwrap(),
        RecurrenceTransformation::MultiplyPower { power: 1 }
    );
}

#[test]
fn vector_maps_and_lossy_maps_do_not_act_componentwise() {
    let state = vec![vec![q(1), q(1)], vec![q(1)]];
    for map in [
        decomposition(),
        RecurrenceTransformation::Derivative { order: 1 },
    ] {
        assert!(matches!(
            map.map_state(1, &state, &[0, -1]),
            Err(TransformError::TransportUnsupported(_))
        ));
    }
    // Component count mismatch.
    assert!(RecurrenceTransformation::MultiplyPower { power: 1 }
        .map_state(1, &state, &[0])
        .is_err());
    // Non-Möbius maps have no vector transport.
    let (rule, _) = eulerian_rows(4);
    let TransportedRecurrence::Vector(system) = decomposition().transport(&rule).unwrap() else {
        panic!("vector system");
    };
    let diagonal = RecurrenceTransformation::Diagonal {
        weight: DiagonalWeight::Binomial { m: 1 },
    };
    assert!(matches!(
        transport_vector_componentwise(&system, &diagonal, &[0, -1]),
        Err(TransformError::TransportUnsupported(_))
    ));
    assert!(transport_vector_componentwise(&system, &diagonal, &[0]).is_err());
}

#[test]
fn forcing_terms_are_transported_only_with_index_free_factors() {
    let int = |c: i64| BivarPoly::from_coefficients_normalized(vec![vec![q(c)]]);
    // F_(n+1) = F_n + 1 (one component).
    let system = SharedDenominatorVectorRecurrence {
        denominator: int(1),
        matrix: vec![vec![WeylOperator::new(vec![int(1)])]],
        forcing: Some(vec![int(1)]),
    };
    let affine = RecurrenceTransformation::Affine {
        a: "1".into(),
        b: "1".into(),
    };
    assert_eq!(
        transport_vector_componentwise(&system, &affine, &[0]).unwrap(),
        system
    );
    let reversal = RecurrenceTransformation::Reversal {
        schedule: DegreeSchedule::new(0, 0, 1),
    };
    assert_eq!(
        transport_vector_componentwise(&system, &reversal, &[0]).unwrap_err(),
        TransformError::InhomogeneousNotTransportable
    );
}
