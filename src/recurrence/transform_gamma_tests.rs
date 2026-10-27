//! Tests of gamma-coordinate transport.

use super::*;

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

fn catalog(id: &str) -> (Recurrence, Vec<Vec<Q>>, usize) {
    let (rule, seeds, first) = crate::oeis::by_id(id)
        .unwrap()
        .recurrence_export_parts()
        .unwrap();
    (rule, seeds, usize::try_from(first).unwrap())
}

fn generate(rule: &Recurrence, seeds: &[Vec<Q>], first: usize, count: usize) -> Vec<Vec<Q>> {
    rule.generate_rows_rational(seeds, first, count)
        .unwrap()
        .into_iter()
        .map(|r| vec_poly::trim(&r))
        .collect()
}

/// Map `rows` to gamma coordinates, transport `rule`, and check that the
/// transported rule regenerates every gamma row from the gamma seeds.
fn assert_gamma_rule_reproduces(
    rule: &Recurrence,
    seeds: usize,
    first: usize,
    rows: &[Vec<Q>],
    schedule: &DegreeSchedule,
) -> Recurrence {
    let map = RecurrenceTransformation::Gamma {
        schedule: schedule.clone(),
    };
    let gamma_rows = rows
        .iter()
        .enumerate()
        .map(|(i, p)| map.map_row(first + i, p).unwrap().remove(0))
        .collect::<Vec<_>>();
    let TransportedRecurrence::Scalar(transported) = map.transport(rule).unwrap() else {
        panic!("expected a scalar gamma rule");
    };
    let regenerated = generate(&transported, &gamma_rows[..seeds], first, rows.len());
    assert_eq!(regenerated, gamma_rows, "{transported}");
    for target in first + seeds..first + rows.len() {
        assert_eq!(
            transformed_rule_holds_at(&transported, &gamma_rows, first, target),
            Some(true)
        );
    }
    transported
}

/// Acceptance: Eulerian A008292, Narayana A001263 and type B Eulerian A060187
/// transport to scalar gamma rules that regenerate the gamma rows exactly.
#[test]
fn gamma_transport_of_eulerian_narayana_and_type_b() {
    let cases = [
        (
            "A008292",
            "P(n) = (1 - 4t + 2nt) P(n-1) + (t - 4t^2) P'(n-1)",
        ),
        (
            "A001263",
            "(1 + n) P(n) = (-1 + 2n) P(n-1) + (2 - n - 8t + 4nt) P(n-2)",
        ),
        (
            "A060187",
            "P(n) = (1 - 8t + 4nt) P(n-1) + (2t - 8t^2) P'(n-1)",
        ),
    ];
    for (id, expected) in cases {
        let (rule, seeds, first) = catalog(id);
        // Every family has row index 1 of degree 0 and symmetry index n - 1.
        assert_eq!(first, 1, "{id}");
        let schedule = DegreeSchedule::new(1, 0, 1);
        let gamma = transport_gamma(&rule, &schedule).unwrap();
        assert!(
            gamma.proportional && gamma.odd.vanishes,
            "{id}: {}",
            gamma.odd.text
        );
        assert!(gamma.obstruction.is_none());
        let rows = generate(&rule, &seeds, first, 14);
        let transported = assert_gamma_rule_reproduces(&rule, seeds.len(), first, &rows, &schedule);
        assert_eq!(transported.to_string(), expected, "{id}");
    }
    // Gamma vectors of the Eulerian row 5 and the Narayana row 5.
    let map = RecurrenceTransformation::Gamma {
        schedule: DegreeSchedule::new(1, 0, 1),
    };
    let (rule, seeds, _) = catalog("A008292");
    let rows = generate(&rule, &seeds, 1, 5);
    assert_eq!(map.map_row(5, &rows[4]).unwrap()[0], row(&[1, 22, 16]));
    let (rule, seeds, _) = catalog("A001263");
    let rows = generate(&rule, &seeds, 1, 5);
    assert_eq!(map.map_row(5, &rows[4]).unwrap()[0], row(&[1, 6, 2]));
}

/// A second-order derivative two rows back exercises the chain-rule recursion
/// `T_(k+1,i) = (d - k - 2i) T_(k,i) + u T_(k,i)' + (1 - x) T_(k,i-1)`.
#[test]
fn gamma_transport_of_a_second_derivative_rule() {
    // E_m = (1 + (m-1)x) E_(m-1) + x(1-x) E'_(m-1) with E'_(m-1) expanded by
    // differentiating the rule one step back.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[1, -1], &[0, 1]])),
            RecurrenceTerm::new(2, 0, bivar(&[&[0, -2, 2], &[0, 1, -1]])),
            RecurrenceTerm::new(2, 1, bivar(&[&[0, 2, -6, 4], &[0, 0, 1, -1]])),
            RecurrenceTerm::new(2, 2, bivar(&[&[0, 0, 1, -2, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let (eulerian, seeds, _) = catalog("A008292");
    let rows = generate(&eulerian, &seeds, 1, 12);
    assert_eq!(generate(&rule, &rows[..2], 1, 12), rows);
    let schedule = DegreeSchedule::new(1, 0, 1);
    let map = RecurrenceTransformation::Gamma {
        schedule: schedule.clone(),
    };
    let gamma_rows = rows
        .iter()
        .enumerate()
        .map(|(i, p)| map.map_row(1 + i, p).unwrap().remove(0))
        .collect::<Vec<_>>();
    // The expansion is not reversal invariant: both parts are needed, and
    // transport refuses to pick one although these rows are symmetric.
    let gamma = transport_gamma(&rule, &schedule).unwrap();
    assert!(!gamma.proportional && gamma.obstruction.is_some());
    assert!(matches!(
        map.transport(&rule),
        Err(TransformError::TransportUnsupported(_))
    ));
    // Each part is still an exact consequence on symmetric rows.
    assert!(!gamma.even.leader_vanishes && gamma.odd.leader_vanishes);
    assert!(gamma
        .even
        .relation
        .terms
        .iter()
        .any(|t| t.offset == 2 && t.deriv_order == 2));
    assert_eq!(
        generate(&gamma.even.relation, &gamma_rows[..2], 1, 12),
        gamma_rows
    );
    for target in 3..13 {
        assert_eq!(
            transformed_rule_holds_at(&gamma.odd.relation, &gamma_rows, 1, target),
            Some(true),
            "{}",
            gamma.odd.text
        );
    }
    let mut broken = gamma_rows.clone();
    broken[9] = vec_poly::add(&broken[9], &[q(0), q(1)]);
    assert_eq!(
        transformed_rule_holds_at(&gamma.odd.relation, &broken, 1, 11),
        Some(false)
    );
}

/// Negative controls: a rule that is not reversal invariant has independent
/// even and odd parts, and non-symmetric rows are rejected by the map.
#[test]
fn gamma_transport_rejects_non_symmetric_rules_and_rows() {
    // P_m = (1 + 2x) P_(m-1): the odd part forces G_(m-1) = 0.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 2]]))],
        denominator: None,
        inhomogeneous: None,
    };
    let schedule = DegreeSchedule::new(0, 0, 1);
    let gamma = transport_gamma(&rule, &schedule).unwrap();
    assert!(!gamma.proportional);
    assert!(gamma.scalar.is_none());
    assert_eq!(gamma.even.text, "P(n) = 3/2 P(n-1)");
    assert_eq!(gamma.odd.text, "0 = -1/2 P(n-1)");
    let map = RecurrenceTransformation::Gamma {
        schedule: schedule.clone(),
    };
    match map.transport(&rule) {
        Err(TransformError::TransportUnsupported(message)) => {
            assert!(message.contains("independent"), "{message}")
        }
        other => panic!("expected an obstruction, got {other:?}"),
    }
    assert_eq!(
        map.map_row(1, &row(&[1, 2])),
        Err(TransformError::NotSymmetric {
            index: 1,
            symmetry_index: 1
        })
    );
    // A report on these rows is rejected at the first non-symmetric row.
    let source = RecurrenceJson::from_recurrence_rational(&rule, 0, &[row(&[1])], None);
    let request = RecurrenceTransformRequest::from_source(map, source, 4);
    assert!(transform_recurrence(&request)
        .unwrap_err()
        .contains("not symmetric"));

    // An additive term picks up the index-dependent factor (1 + x)^(-d_m).
    let inhomogeneous = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 1]]))],
        denominator: None,
        inhomogeneous: Some(bivar(&[&[1]])),
    };
    assert_eq!(
        transport_gamma(&inhomogeneous, &schedule).unwrap_err(),
        TransformError::InhomogeneousNotTransportable
    );
}

/// The transport report certifies symmetry from reversal invariance (also
/// with a nonvanishing left factor) and replays exactly.
#[test]
fn gamma_transport_report_certifies_and_replays() {
    for id in ["A008292", "A001263"] {
        let request = RecurrenceTransformRequest {
            transformation: RecurrenceTransformation::Gamma {
                schedule: DegreeSchedule::new(1, 0, 1),
            },
            source_recurrence: None,
            oeis_id: Some(id.into()),
            rows: None,
            row_origin: None,
            row_count: Some(10),
            search: None,
            limits: TransformLimits::default(),
            root_domain: None,
        };
        let report = transform_recurrence(&request).unwrap();
        assert_eq!(report.conditions[0].status, CERTIFIED_ALL_INDICES, "{id}");
        let transported = report.transported.as_ref().unwrap();
        assert!(transported.holds_on_transformed_rows);
        assert!(report.transport_error.is_none());
        verify_recurrence_transform(&report).unwrap();
    }
}

#[test]
fn symmetric_gamma_and_split_rows() {
    // 1 + 4x + x^2 = (1 + x)^2 + 2x.
    assert_eq!(symmetric_gamma(&row(&[1, 4, 1]), 2), Some(row(&[1, 2])));
    // Endpoint factors: x (1 + x) at index 3 is x (1 + x)^(3 - 2).
    assert_eq!(symmetric_gamma(&row(&[0, 1, 1]), 3), Some(row(&[0, 1])));
    assert_eq!(symmetric_gamma(&row(&[1, 2]), 1), None);
    // 1 + 2x = (3/2)(1 + x) + (1 - x)(-1/2) at index 1.
    let (even, odd) = split_row(&row(&[1, 2]), 1);
    assert_eq!(even, vec![Q::new(3.into(), 2.into())]);
    assert_eq!(odd, vec![Q::new((-1).into(), 2.into())]);
}
