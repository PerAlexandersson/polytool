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
                    let value = super::super::implicit_algebra::binomial(n, k)
                        / super::super::implicit_algebra::factorial(k);
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
    assert!(error.contains("scalar"));
}
