//! Tests of the direct search for implicit relations and conditions.

use super::*;
use num_traits::One;
use serde_json::json;

fn request(value: serde_json::Value) -> ImplicitRelationSearchRequest {
    serde_json::from_value(value).expect("request")
}

fn factorial_map() -> serde_json::Value {
    json!({"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}})
}

#[test]
fn pascal_rows_divided_by_factorials_recover_the_implicit_rule_and_its_condition() {
    let report = search_implicit_relations(&request(json!({
        "oeis_id": "A007318",
        "transformation": factorial_map(),
        "supports": [{"role": "implicit_successor", "max_shift": 1, "derivative_order": 1}]
    })))
    .unwrap();
    assert_eq!(report.source.constant_index_order, Some(1));
    let result = &report.results[0];
    assert_eq!(result.status, "relations_found");
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    // H = U + m(A+J) + 1 = 2 residuals for the order-1 source.
    assert_eq!(result.annihilator_power, Some(2));
    assert_eq!(result.residual_count, 2);
    assert_eq!(result.kernel_dimension, 1);
    let found = &result.representatives[0];
    assert_eq!(found.relation.text, "D Q[n+1] - D Q[n] - Q[n] = 0");
    assert_eq!(found.relation.role, IMPLICIT_SUCCESSOR_ROLE);
    assert_eq!(found.relation.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    // The constant condition: a complete condition space (r (k+1)(U+1) = 2 bases).
    let search = &found.condition_searches[0];
    assert_eq!(search.coefficient_index, 0);
    assert_eq!(search.status, "condition_found");
    assert_eq!(search.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    assert_eq!(search.residual_count, 2);
    let condition = search.condition.as_ref().unwrap();
    assert_eq!(condition.text, "[x^0] Q[n+1] - [x^0] Q[n] = 0");
    assert_eq!(found.description.seeds, vec![vec!["1"]]);
    assert_eq!(found.assessment.status, EQUIVALENT_DESCRIPTION);
    assert_eq!(
        found.assessment.uniqueness.grading.as_deref(),
        Some("valuation")
    );
    verify_implicit_relation_search(&report).unwrap();

    // The same relation as the transported one.
    let rule = crate::oeis::by_id("A007318")
        .unwrap()
        .recurrence_export_parts()
        .unwrap()
        .0;
    let (source, _) = ImplicitRelation::from_recurrence(&rule).unwrap();
    let transported = RecurrenceTransformation::Diagonal {
        weight: DiagonalWeight::FactorialPower { power: 1 },
    }
    .transport_implicit(&source, 0)
    .unwrap();
    let direct = ImplicitRelation::from_json(&found.relation.relation).unwrap();
    assert!(direct.same_relation(&transported.relation));

    // Without the condition only a consequence.
    let mut weaker = found.description.clone();
    weaker.boundary_conditions.clear();
    assert_eq!(
        assess_generating_description(&weaker).unwrap().status,
        CERTIFIED_CONSEQUENCE
    );

    // Tampering with a row, a kernel certificate or a label is caught.
    let mut tampered = report.clone();
    tampered.rows[2][1] = "3".into();
    assert!(verify_implicit_relation_search(&tampered).is_err());
    let mut tampered = report.clone();
    tampered.results[0].certificate.as_mut().unwrap().basis[0][0] = "7".into();
    assert!(verify_implicit_relation_search(&tampered).is_err());
    let mut tampered = report;
    tampered.results[0].representatives[0].relation.evidence = FINITE_FIT_EVIDENCE.into();
    assert!(verify_implicit_relation_search(&tampered).is_err());
}

#[test]
fn row_equation_is_found_directly_and_is_not_a_successor_rule() {
    let report = search_implicit_relations(&request(json!({
        "oeis_id": "A007318",
        "transformation": factorial_map(),
        "supports": [{"role": "row_equation", "derivative_order": 2,
                      "index_degree": 1, "variable_degree": 1}]
    })))
    .unwrap();
    let result = &report.results[0];
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    // H = 1 + 1 (1 + 2) + 1 = 5.
    assert_eq!(result.residual_count, 5);
    assert_eq!(result.kernel_dimension, 1);
    let found = &result.representatives[0];
    assert_eq!(
        found.relation.text,
        "x D^2 Q[n] + (1 + x) D Q[n] - n Q[n] = 0"
    );
    assert_eq!(found.relation.role, ROW_EQUATION_ROLE);
    // A row equation alone determines nothing beyond its seeds.
    assert_eq!(found.assessment.status, CERTIFIED_CONSEQUENCE);
    assert!(found.condition_searches.is_empty());
    verify_implicit_relation_search(&report).unwrap();
}

#[test]
fn derivative_of_pascal_rows_has_a_complete_successor_space() {
    let report = search_implicit_relations(&request(json!({
        "oeis_id": "A007318",
        "transformation": {"kind": "derivative", "order": 1},
        "supports": [{"role": "successor", "max_shift": 1, "index_degree": 1,
                      "variable_degree": 1}]
    })))
    .unwrap();
    assert_eq!(report.symbol, "Q");
    let result = &report.results[0];
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    // h = max(u + s) + m + 1 = 3.
    assert_eq!(result.annihilator_power, Some(3));
    let found = &result.representatives[0];
    // Q_n = n (1 + x)^(n-1): n Q_(n+1) = (n + 1)(1 + x) Q_n.
    assert_eq!(found.relation.role, SUCCESSOR_ROLE);
    let relation = ImplicitRelation::from_json(&found.relation.relation).unwrap();
    let rows = (0..8)
        .map(|n: usize| {
            let p = (0..n).fold(vec![Q::one()], |acc, _| {
                vec_poly::mul(&acc, &[Q::one(), Q::one()])
            });
            vec_poly::derivative(&p)
        })
        .collect::<Vec<_>>();
    for base in 0..7 {
        assert_eq!(relation.holds_at(&rows, 0, base), Some(true));
    }
    verify_implicit_relation_search(&report).unwrap();
}

#[test]
fn eulerian_derivative_relation_is_certified_by_division() {
    let report = search_implicit_relations(&request(json!({
        "oeis_id": "A008292",
        "row_count": 10,
        "transformation": {"kind": "derivative", "order": 1},
        "supports": [{"role": "implicit_successor", "max_shift": 1, "derivative_order": 2,
                      "index_degree": 1, "variable_degree": 2}]
    })))
    .unwrap();
    assert_eq!(report.source.constant_index_order, None);
    let result = &report.results[0];
    assert_eq!(result.evidence, FINITE_FIT_EVIDENCE);
    assert!(result.completeness.is_none());
    // The transported relation M = D^2 L D^-1 is among the representatives,
    // proved by right division by the Eulerian rule.
    let rule = crate::oeis::by_id("A008292")
        .unwrap()
        .recurrence_export_parts()
        .unwrap()
        .0;
    let (source, _) = ImplicitRelation::from_recurrence(&rule).unwrap();
    let transported = RecurrenceTransformation::Derivative { order: 1 }
        .transport_implicit(&source, 1)
        .unwrap();
    assert_eq!(
        transported.relation.text("Q"),
        "D Q[n+1] + (-x + x^2) D^2 Q[n] + (-3 + 4 x - n x) D Q[n] + (2 - 2 n) Q[n] = 0"
    );
    transported
        .reconstruction
        .as_ref()
        .unwrap()
        .verify_identity()
        .unwrap();
    let certified = result
        .representatives
        .iter()
        .filter(|r| r.relation.evidence == CERTIFIED_BY_DIVISION_EVIDENCE)
        .collect::<Vec<_>>();
    assert!(!certified.is_empty());
    assert!(certified.iter().all(|r| r.division_certificate.is_some()));
    assert!(result.representatives.iter().any(|r| {
        ImplicitRelation::from_json(&r.relation.relation)
            .unwrap()
            .same_relation(&transported.relation)
    }));
    verify_implicit_relation_search(&report).unwrap();
}

#[test]
fn rows_only_search_gives_finite_fits_and_certified_absence() {
    // Q_n = sum_k binom(n,k) x^k / k! as plain rows.
    let rows = (0..8usize)
        .map(|n| {
            (0..=n)
                .map(|k| {
                    let value = crate::recurrence::implicit::algebra::binomial(n, k)
                        / crate::recurrence::implicit::algebra::factorial(k);
                    format_rational_coeff(&value)
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let report = search_implicit_relations(&request(json!({
        "rows": rows,
        "supports": [
            {"role": "implicit_successor", "max_shift": 1, "derivative_order": 1},
            {"role": "successor", "max_shift": 1}
        ]
    })))
    .unwrap();
    assert_eq!(report.source.kind, "rows_only");
    let found = &report.results[0];
    assert_eq!(found.evidence, FINITE_FIT_EVIDENCE);
    let rep = &found.representatives[0];
    assert_eq!(rep.relation.text, "D P[n+1] - D P[n] - P[n] = 0");
    assert_eq!(rep.relation.evidence, FINITE_FIT_EVIDENCE);
    assert_eq!(rep.condition_searches[0].evidence, FINITE_FIT_EVIDENCE);
    assert_eq!(rep.assessment.status, FINITE_EVIDENCE);
    // No constant-coefficient first-order successor rule: certified absence.
    let absent = &report.results[1];
    assert_eq!(absent.status, "absent_in_support");
    assert_eq!(absent.evidence, "rows_only_obstruction");
    assert!(absent.representatives.is_empty());
    verify_implicit_relation_search(&report).unwrap();
}

#[test]
fn invalid_supports_and_maps_are_rejected() {
    for support in [
        json!({"role": "row_equation", "max_shift": 1}),
        json!({"role": "successor", "max_shift": 0}),
        json!({"role": "sideways", "max_shift": 1}),
    ] {
        assert!(search_implicit_relations(&request(json!({
            "oeis_id": "A007318", "supports": [support]
        })))
        .is_err());
    }
    let error = search_implicit_relations(&request(json!({
        "oeis_id": "A007318",
        "transformation": {"kind": "symmetric_decomposition",
                           "schedule": {"first_index": 0, "first_degree": 0, "step": 1}}
    })))
    .unwrap_err();
    assert!(error.to_string().contains("scalar"));
}

// ---------------------------------------------------------------------------
// Completeness after Möbius-class maps and polynomial diagonal weights
// ---------------------------------------------------------------------------

/// Fibonacci polynomials `P_(n+2) = P_(n+1) + x P_n`, `P_0 = P_1 = 1`
/// (degree `floor(n/2)`), a constant-index monic source of order 2.
fn fibonacci_source() -> serde_json::Value {
    json!({
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"], ["1"]],
        "recurrence": {"terms": [
            {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1"]]}},
            {"offset": 2, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["0", "1"]]}}
        ]}
    })
}

fn schedule(first_degree: usize, step: i64) -> serde_json::Value {
    json!({"first_index": 0, "first_degree": first_degree, "step": step})
}

fn mobius_class_maps() -> Vec<serde_json::Value> {
    vec![
        json!({"kind": "reversal", "schedule": schedule(0, 1)}),
        json!({"kind": "f_polynomial", "schedule": schedule(0, 1)}),
        json!({"kind": "mobius", "a": "1", "b": "2", "c": "1", "e": "1",
               "schedule": schedule(1, 1)}),
        json!({"kind": "affine", "a": "2", "b": "1"}),
        json!({"kind": "multiply_power", "power": 2}),
    ]
}

/// The transported source rule of `map` in implicit form.
fn transported_source(map: &serde_json::Value) -> ImplicitRelation {
    let t: RecurrenceTransformation = serde_json::from_value(map.clone()).unwrap();
    let json: RecurrenceJson = serde_json::from_value(fibonacci_source()).unwrap();
    let (rule, _, _) = json.to_recurrence_parts().unwrap();
    let TransportedRecurrence::Scalar(transported) = t.transport(&rule).unwrap() else {
        panic!("scalar transport");
    };
    ImplicitRelation::from_recurrence(&transported).unwrap().0
}

/// Exact rows `T(P_n)` for `n < count`, generated independently of the search.
fn mapped_rows(map: &serde_json::Value, count: usize) -> Vec<Vec<Q>> {
    let t: RecurrenceTransformation = serde_json::from_value(map.clone()).unwrap();
    let json: RecurrenceJson = serde_json::from_value(fibonacci_source()).unwrap();
    let (rule, first, seeds) = json.to_recurrence_parts().unwrap();
    rule.generate_rows_rational(&seeds, first, count)
        .unwrap()
        .iter()
        .enumerate()
        .map(|(n, row)| t.map_row(n, row).unwrap().remove(0))
        .collect()
}

#[test]
fn mobius_class_maps_have_complete_spaces_with_the_identity_count() {
    for map in mobius_class_maps() {
        let report = search_implicit_relations(&request(json!({
            "source_recurrence": fibonacci_source(),
            "transformation": map,
            "supports": [{"role": "successor", "max_shift": 2, "variable_degree": 2}]
        })))
        .unwrap();
        assert_eq!(report.source.constant_index_order, Some(2), "{map}");
        assert!(
            report
                .source
                .completeness_basis
                .contains("transported-source"),
            "{map}"
        );
        let result = &report.results[0];
        assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE, "{map}");
        // h = max(u + s) + 1 = 1: k h = 2 residuals, as for the identity.
        assert_eq!(result.annihilator_power, Some(1), "{map}");
        assert_eq!(result.residual_count, 2, "{map}");
        assert!(result
            .completeness
            .as_ref()
            .unwrap()
            .contains("Möbius-class"));
        let expected = transported_source(&map);
        assert!(
            result.representatives.iter().any(|r| {
                ImplicitRelation::from_json(&r.relation.relation)
                    .unwrap()
                    .same_relation(&expected)
            }),
            "{map}: {:?}",
            result
                .representatives
                .iter()
                .map(|r| &r.relation.text)
                .collect::<Vec<_>>()
        );
        // Every certified relation holds far beyond the rows used.
        let rows = mapped_rows(&map, 24);
        for rep in &result.representatives {
            let relation = ImplicitRelation::from_json(&rep.relation.relation).unwrap();
            for base in 0..22 {
                assert_eq!(relation.holds_at(&rows, 0, base), Some(true), "{map}");
            }
        }
        verify_implicit_relation_search(&report).unwrap();
    }
}

#[test]
fn the_residual_count_is_needed_and_enforced() {
    let map = json!({"kind": "reversal", "schedule": schedule(0, 1)});
    let support = json!({"role": "successor", "max_shift": 2, "index_degree": 1,
                         "variable_degree": 2});
    let report = search_implicit_relations(&request(json!({
        "source_recurrence": fibonacci_source(),
        "transformation": map,
        "supports": [support]
    })))
    .unwrap();
    let result = &report.results[0];
    // h = max(u + s) + 1 = 2: k h = 4 residuals.
    assert_eq!(result.residual_count, 4);
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    // Fewer residuals admit spurious relations: the bound is not vacuous here.
    let rows = mapped_rows(&map, 12)
        .into_iter()
        .map(Poly::new)
        .collect::<Vec<_>>();
    let limits = RecurrenceSpaceLimits::default();
    let dimension = |count: usize| {
        let matrix = homogeneous_matrix(&rows, 0, 0, count, &result.monomials, &limits).unwrap();
        rational_kernel_with_certificate(&matrix, result.monomials.len(), kernel_limits(&limits))
            .unwrap()
            .basis
            .len()
    };
    assert_eq!(dimension(4), result.kernel_dimension);
    assert!(dimension(1) > dimension(4));
    // Too few rows for the theorem: never certified on fewer residuals.
    let short = search_implicit_relations(&request(json!({
        "source_recurrence": fibonacci_source(),
        "transformation": map,
        "row_count": 4,
        "supports": [support],
        "limits": {"max_rows": 5}
    })))
    .unwrap();
    assert_eq!(short.results[0].status, "incomplete_budget");
    assert_eq!(short.results[0].evidence, "none");
    assert!(short.results[0].representatives.is_empty());
    // A lowered residual count or annihilator power does not replay.
    let mut tampered = report.clone();
    tampered.results[0].residual_count = 3;
    assert!(verify_implicit_relation_search(&tampered).is_err());
    let mut tampered = report.clone();
    tampered.results[0].annihilator_power = Some(1);
    assert!(verify_implicit_relation_search(&tampered).is_err());
    let mut tampered = report;
    tampered.source.completeness_basis = "identity".into();
    assert!(verify_implicit_relation_search(&tampered).is_err());
}

#[test]
fn bounded_schedules_and_uncovered_weights_stay_finite_fits() {
    for map in [
        json!({"kind": "reversal", "schedule":
               {"first_index": 0, "first_degree": 0, "step": 1, "last_index": 40}}),
        json!({"kind": "diagonal", "weight": {"kind": "inverse_binomial", "m": 1}}),
        json!({"kind": "diagonal", "weight": {"kind": "product", "scale": "1",
               "factors": [{"coeffs": {"coeffs": [["1", "1"]]}, "exponent": -1}]}}),
        json!({"kind": "derivative_antecedent", "g": ["0", "1"]}),
    ] {
        let report = search_implicit_relations(&request(json!({
            "source_recurrence": fibonacci_source(),
            "transformation": map,
            "supports": [{"role": "successor", "max_shift": 2, "variable_degree": 2}]
        })))
        .unwrap();
        assert!(
            report
                .source
                .completeness_basis
                .starts_with("no finite theorem"),
            "{map}"
        );
        let result = &report.results[0];
        assert_ne!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE, "{map}");
        assert!(result.annihilator_power.is_none() && result.completeness.is_none());
        verify_implicit_relation_search(&report).unwrap();
    }
}

#[test]
fn transported_coefficients_with_a_pole_at_zero_are_detected() {
    let source = |a0: &[&str]| -> RecurrenceJson {
        serde_json::from_value(json!({
            "schema": "polytool.recurrence.v1", "first_index": 0,
            "initial_polynomials": [["1"]],
            "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                                      "coeff": {"coeffs": [a0]}}]}
        }))
        .unwrap()
    };
    let reversal = RecurrenceTransformation::Reversal {
        schedule: DegreeSchedule::new(0, 0, 1),
    };
    // P_(n+1) = (1 + x) P_n: a~ = (1 + 1/x) x = 1 + x.
    assert_eq!(
        transported_regular_at_zero(&reversal, &source(&["1", "1"])),
        Some(true)
    );
    // P_(n+1) = x^2 P_n under d_n = n: a~ = x^-2 x = 1/x, a pole at zero.
    assert_eq!(
        transported_regular_at_zero(&reversal, &source(&["0", "0", "1"])),
        Some(false)
    );
    let theorem = map_theorem(Some(&reversal), &source(&["0", "0", "1"])).unwrap();
    assert!(
        theorem
            == Theorem::Transported {
                regular_at_zero: false
            }
    );
    assert_eq!(condition_theorem(theorem, 3), None);
    assert_eq!(
        condition_theorem(
            Theorem::Transported {
                regular_at_zero: true
            },
            3
        ),
        Some((3, 0))
    );
}

/// `(1 + x)^n` with origin 0.
fn pascal_source() -> serde_json::Value {
    json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                                  "coeff": {"coeffs": [["1", "1"]]}}]}
    })
}

#[test]
fn polynomial_diagonal_weights_have_complete_spaces() {
    // binom(k+1, 1) = k + 1: Q_n = D(x P_n) = (1 + (n+1) x)(1 + x)^(n-1).
    let binomial = json!({"kind": "diagonal", "weight": {"kind": "binomial", "m": 1}});
    // The same weight written as a product, and an index-dependent weight
    // 2 (n + k + 1).
    let product = json!({"kind": "diagonal", "weight": {"kind": "product", "scale": "1",
        "factors": [{"coeffs": {"coeffs": [["1", "1"]]}, "exponent": 1}]}});
    let indexed = json!({"kind": "diagonal", "weight": {"kind": "product", "scale": "2",
        "factors": [{"coeffs": {"coeffs": [["1", "1"], ["1"]]}, "exponent": 1}]}});
    let support = json!({"role": "successor", "max_shift": 1, "index_degree": 1,
                         "variable_degree": 2});
    let mut texts = Vec::new();
    for (map, index) in [(&binomial, 0), (&product, 0), (&indexed, 1)] {
        let report = search_implicit_relations(&request(json!({
            "source_recurrence": pascal_source(),
            "transformation": map,
            "supports": [support]
        })))
        .unwrap();
        assert!(report.source.completeness_basis.contains("weight theorem"));
        assert!(report
            .source
            .completeness_basis
            .contains(&format!("index degree {index}")));
        let result = &report.results[0];
        assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE, "{map}");
        // h = max(u + s) + w + 1 with w = 1.
        assert_eq!(result.annihilator_power, Some(3), "{map}");
        assert_eq!(result.residual_count, 3);
        assert_eq!(result.kernel_dimension, 1, "{map}");
        let rep = &result.representatives[0];
        assert_eq!(rep.assessment.status, EQUIVALENT_DESCRIPTION, "{map}");
        let relation = ImplicitRelation::from_json(&rep.relation.relation).unwrap();
        let t: RecurrenceTransformation = serde_json::from_value(map.clone()).unwrap();
        let rows = (0..20usize)
            .map(|n| {
                let p = (0..n).fold(vec![Q::one()], |acc, _| {
                    vec_poly::mul(&acc, &[Q::one(), Q::one()])
                });
                t.map_row(n, &p).unwrap().remove(0)
            })
            .collect::<Vec<_>>();
        for base in 0..19 {
            assert_eq!(relation.holds_at(&rows, 0, base), Some(true), "{map}");
        }
        verify_implicit_relation_search(&report).unwrap();
        texts.push(rep.relation.text.clone());
    }
    assert_eq!(texts[0], texts[1]);
    // A factorial weight with a polynomial factor: h = U + m(A+J) + w + 1.
    let mixed = json!({"kind": "diagonal", "weight": {"kind": "product", "scale": "1",
        "factorial_exponent": -1,
        "factors": [{"coeffs": {"coeffs": [["1", "1"]]}, "exponent": 1}]}});
    let report = search_implicit_relations(&request(json!({
        "source_recurrence": pascal_source(),
        "transformation": mixed,
        "supports": [{"role": "implicit_successor", "max_shift": 1, "derivative_order": 1,
                      "index_degree": 1, "variable_degree": 1}]
    })))
    .unwrap();
    let result = &report.results[0];
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    assert_eq!(result.annihilator_power, Some(1 + 2 + 1 + 1));
    verify_implicit_relation_search(&report).unwrap();
    assert_eq!(
        condition_theorem(
            Theorem::Weight {
                m: 0,
                total: 1,
                index: 1
            },
            2
        ),
        Some((2, 1))
    );
}

#[test]
fn divide_power_has_a_complete_space() {
    // P_n = x (1 + x)^n, so Q_n = P_n / x = (1 + x)^n.
    let source = json!({
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["0", "1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                                  "coeff": {"coeffs": [["1", "1"]]}}]}
    });
    let report = search_implicit_relations(&request(json!({
        "source_recurrence": source,
        "transformation": {"kind": "divide_power", "power": 1},
        "supports": [{"role": "successor", "max_shift": 1, "index_degree": 1,
                      "variable_degree": 1}]
    })))
    .unwrap();
    let result = &report.results[0];
    assert_eq!(result.evidence, CERTIFIED_COMPLETE_SPACE_EVIDENCE);
    assert_eq!(result.residual_count, 2);
    assert_eq!(report.rows[3], vec!["1", "3", "3", "1"]);
    verify_implicit_relation_search(&report).unwrap();
}
