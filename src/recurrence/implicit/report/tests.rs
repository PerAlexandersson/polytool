//! Tests of implicit relations, diagonal transport, derivative antecedents,
//! generating descriptions and reconstruction certificates.

use super::*;
use crate::recurrence::implicit::algebra::Biv;
use serde_json::json;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn request(value: serde_json::Value) -> ImplicitTransformRequest {
    serde_json::from_value(value).expect("request")
}

fn binomial_factorial_request() -> serde_json::Value {
    json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "transformation": {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
    })
}

/// `x Q'' + (1 + x) Q' - n Q = 0`.
fn row_equation_json() -> serde_json::Value {
    json!({"terms": [
        {"shift": 0, "deriv_order": 2, "coeff": {"coeffs": [["0", "1"]]}},
        {"shift": 0, "deriv_order": 1, "coeff": {"coeffs": [["1", "1"]]}},
        {"shift": 0, "deriv_order": 0, "coeff": {"coeffs": [["0"], ["-1"]]}}
    ]})
}

#[test]
fn binomial_rows_divided_by_factorials_give_an_implicit_rule_and_a_lost_constant() {
    let report = implicit_transform_recurrence(&request(binomial_factorial_request())).unwrap();
    // Q_n = sum_k binom(n,k) x^k / k!.
    assert_eq!(report.transformed_rows[3], vec!["1", "3", "3/2", "1/6"]);
    let defining = &report.transports[0];
    assert_eq!(defining.relation.text, "D Q[n+1] - D Q[n] - Q[n] = 0");
    assert_eq!(defining.relation.role, IMPLICIT_SUCCESSOR_ROLE);
    assert_eq!(defining.relation.evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert!(defining.holds_on_transformed_rows);
    let clearing = defining.clearing.as_ref().unwrap();
    assert_eq!(clearing.multiplier_text, "(theta)");
    assert_eq!(clearing.exceptional.status, CERTIFIED_ALL_INDICES);
    assert_eq!(clearing.exceptional.indices, vec![0]);
    // [x^0] Q_(n+1) = [x^0] Q_n.
    let constant = &defining.exceptional_index_conditions[0];
    assert_eq!(constant.coefficient_index, Some(0));
    assert_eq!(constant.text, "[x^0] Q[n+1] - [x^0] Q[n] = 0");
    assert_eq!(report.description.seeds, vec![vec!["1"]]);
    assert_eq!(report.assessment.status, EQUIVALENT_DESCRIPTION);
    assert_eq!(
        report.assessment.uniqueness.grading.as_deref(),
        Some("valuation")
    );
    assert_eq!(report.assessment.uniqueness.exceptional_indices, vec![0]);
    verify_recurrence_implicit(&report).unwrap();

    // Removing the constant condition removes the equivalence claim.
    let mut weaker = report.description.clone();
    weaker.exceptional_index_conditions.clear();
    let assessment = assess_generating_description(&weaker).unwrap();
    assert_eq!(assessment.status, CERTIFIED_CONSEQUENCE);
    assert_eq!(assessment.uniqueness.status, "not_certified");
    assert!(assessment
        .uniqueness
        .reasons
        .iter()
        .any(|r| r.contains("[0]")));

    // A supplied boundary condition restores uniqueness, but only as finite evidence.
    let mut supplied = weaker.clone();
    supplied
        .boundary_conditions
        .push(report.description.exceptional_index_conditions[0].clone());
    supplied.boundary_conditions[0].evidence = CHECKED_SUPPLIED_ROWS.into();
    let assessment = assess_generating_description(&supplied).unwrap();
    assert_eq!(assessment.uniqueness.status, "certified");
    assert_eq!(assessment.status, FINITE_EVIDENCE);

    let mut tampered = report.clone();
    tampered.transformed_rows[2][1] = "3".into();
    assert!(verify_recurrence_implicit(&tampered).is_err());
}

#[test]
fn row_differential_equation_is_stored_without_being_a_successor_rule() {
    let mut value = binomial_factorial_request();
    value["target_relations"] = json!([{"relation": row_equation_json(), "first_base": 0}]);
    // The source row equation (1 + x) P' - n P = 0 transports to the same relation.
    value["source_relations"] = json!([{"relation": {"terms": [
        {"shift": 0, "deriv_order": 1, "coeff": {"coeffs": [["1", "1"]]}},
        {"shift": 0, "deriv_order": 0, "coeff": {"coeffs": [["0"], ["-1"]]}}
    ]}, "first_base": 0}]);
    let report = implicit_transform_recurrence(&request(value)).unwrap();
    assert_eq!(report.description.relations.len(), 1);
    assert_eq!(report.additional_relations.len(), 2);
    for relation in &report.additional_relations {
        assert_eq!(relation.role, ROW_EQUATION_ROLE);
        assert_ne!(relation.role, SUCCESSOR_ROLE);
    }
    assert_eq!(
        report.additional_relations[0].evidence,
        TRANSPORTED_SUPPLIED_EVIDENCE
    );
    assert_eq!(
        report.additional_relations[1].evidence,
        CHECKED_SUPPLIED_ROWS
    );
    let transported =
        ImplicitRelation::from_json(&report.additional_relations[0].relation).unwrap();
    let supplied = ImplicitRelation::from_json(&report.additional_relations[1].relation).unwrap();
    assert!(transported.same_relation(&supplied));
    assert_eq!(
        supplied.text("Q"),
        "x D^2 Q[n] + (1 + x) D Q[n] - n Q[n] = 0"
    );
    assert_eq!(report.assessment.status, EQUIVALENT_DESCRIPTION);
    verify_recurrence_implicit(&report).unwrap();

    // A false supplied relation is rejected.
    let mut wrong = binomial_factorial_request();
    wrong["target_relations"] = json!([{"relation": {"terms": [
        {"shift": 0, "deriv_order": 1, "coeff": {"coeffs": [["1"]]}}
    ]}, "first_base": 0}]);
    let error = implicit_transform_recurrence(&request(wrong)).unwrap_err();
    assert!(error.contains("target_relations[0] fails"));
}

fn eulerian_rows(count: usize) -> Vec<Vec<String>> {
    crate::oeis::by_id("A008292")
        .unwrap()
        .generate_rows(count)
        .unwrap()
        .into_iter()
        .map(|row| row.iter().map(|c| c.to_string()).collect())
        .collect()
}

#[test]
fn weighted_eulerian_divided_by_k_plus_one_recovers_eulerian_with_replay() {
    for transformation in [
        json!({"kind": "diagonal", "weight": {"kind": "inverse_binomial", "m": 1}}),
        json!({"kind": "derivative_antecedent", "g": ["0", "1"]}),
    ] {
        let report = implicit_transform_recurrence(&request(json!({
            "oeis_id": "A065826",
            "row_count": 9,
            "search": {},
            "transformation": transformation
        })))
        .unwrap();
        assert_eq!(report.transformed_rows, eulerian_rows(9));
        assert!(report.reconstruction_verified);
        let defining = &report.transports[0];
        let certificate = defining.reconstruction.as_ref().expect("certificate");
        certificate.verify_identity().unwrap();
        let rows = report
            .transformed_rows
            .iter()
            .map(|r| decode_poly(r).unwrap())
            .collect::<Vec<_>>();
        assert!(
            certificate
                .replay_on_rows(&rows, report.row_origin)
                .unwrap()
                > 0
        );
        assert_eq!(report.assessment.status, EQUIVALENT_DESCRIPTION);
        assert_eq!(
            report.assessment.uniqueness.grading.as_deref(),
            Some("degree")
        );

        // The fitted Eulerian rule is proved by a reconstruction certificate.
        let search = report.independent_search.as_ref().unwrap();
        assert_eq!(search.status, "found");
        assert_eq!(
            search.evidence, CERTIFIED_BY_RECONSTRUCTION_EVIDENCE,
            "{search:?}"
        );
        let fitted = search.reconstruction.as_ref().unwrap();
        fitted.verify_identity().unwrap();
        assert!(fitted.replay_on_rows(&rows, report.row_origin).unwrap() > 0);
        assert!(fitted
            .multipliers
            .iter()
            .any(|m| m.differentiations.iter().any(|d| d.deriv_order > 0)));
        assert_eq!(
            search.assessment.as_ref().unwrap().status,
            EQUIVALENT_DESCRIPTION
        );
        verify_recurrence_implicit(&report).unwrap();

        // A tampered multiplier breaks the identity.
        let mut broken = fitted.clone();
        broken.multipliers[0].differentiations[0].coeff.coeffs[0][0] = "17".into();
        assert!(broken.verify_identity().is_err());
    }
}

#[test]
fn factorial_rescaling_fit_stays_a_finite_fit() {
    let mut value = binomial_factorial_request();
    value["search"] = json!({});
    let report = implicit_transform_recurrence(&request(value)).unwrap();
    let search = report.independent_search.unwrap();
    assert_eq!(search.status, "found");
    assert_eq!(search.evidence, FINITE_FIT_EVIDENCE);
    assert!(search
        .reason
        .unwrap()
        .contains("not a differential operator"));
}

#[test]
fn zero_weights_are_rejected() {
    for factor in [json!([["-2", "1"]]), json!([["-3"], ["1"]])] {
        let value = json!({
            "oeis_id": "A007318",
            "transformation": {"kind": "diagonal", "weight": {
                "kind": "product", "scale": "1",
                "factors": [{"coeffs": {"coeffs": factor}, "exponent": 1}]}}
        });
        let error = implicit_transform_recurrence(&request(value)).unwrap_err();
        assert!(error.contains("vanishes"), "{error}");
    }
    let weight = DiagonalWeight::Product {
        scale: "1".into(),
        factorial_exponent: 0,
        factors: vec![DiagonalWeightFactor {
            coeffs: BivarPolyJson {
                coeffs: vec![vec!["-2".into(), "1".into()]],
            },
            exponent: -1,
        }],
    };
    assert!(weight.map_row(0, &[q(1), q(1), q(1)]).is_err());
    assert!(weight.map_row(0, &[q(1), q(1)]).is_ok());
}

#[test]
fn exceptional_differential_denominators_are_certified_or_visible() {
    // Q = W_1 P (lambda = k + 1): the denominator theta + 1 of x-columns gives
    // the certified exceptional index 0.
    let report = implicit_transform_recurrence(&request(json!({
        "oeis_id": "A007318",
        "row_count": 6,
        "transformation": {"kind": "diagonal", "weight": {"kind": "binomial", "m": 1}}
    })))
    .unwrap();
    let clearing = report.transports[0].clearing.as_ref().unwrap();
    assert_eq!(clearing.exceptional.status, CERTIFIED_ALL_INDICES);
    assert_eq!(clearing.exceptional.indices, vec![0]);
    verify_recurrence_implicit(&report).unwrap();

    // lambda = 1/(n + 2 - k): an index-dependent factor stays uncertified and visible.
    let report = implicit_transform_recurrence(&request(json!({
        "oeis_id": "A007318",
        "row_count": 8,
        "transformation": {"kind": "diagonal", "weight": {
            "kind": "product", "scale": "1",
            "factors": [{"coeffs": {"coeffs": [["2", "-1"], ["1"]]}, "exponent": -1}]}}
    })))
    .unwrap();
    assert_eq!(report.map_conditions[0].status, CHECKED_SUPPLIED_ROWS);
    let clearing = report.transports[0].clearing.as_ref().unwrap();
    assert_eq!(clearing.exceptional.status, CHECKED_SUPPLIED_ROWS);
    assert!(!clearing.exceptional.observed_zeros.is_empty());
    assert_eq!(report.assessment.status, FINITE_EVIDENCE);
    assert!(report.transports[0]
        .relation
        .note
        .as_deref()
        .unwrap()
        .contains("not certified"));
    verify_recurrence_implicit(&report).unwrap();

    // One common multiplier for every column of an ansatz.
    let factorial = DiagonalWeight::FactorialPower { power: 1 };
    let h = factorial
        .clearing_multiplier(&[(0, 1), (1, 1), (0, 2)])
        .unwrap();
    assert_eq!(h.len(), 2); // theta (theta - 1)
}

/// `(1 + x)^n` described by `Q_(n+1) = (1 + x) Q_n` with the given seeds.
fn pascal_description(leader_root: Option<i64>, seeds: &[&[&str]]) -> GeneratingDescription {
    let factor = leader_root.map_or(vec![vec![q(1)]], |root| vec![vec![q(-root)], vec![q(1)]]);
    let leader = BivarPoly::from_coefficients_normalized(factor);
    let lower = Biv::from_bivar(&leader)
        .mul(&Biv::from_y(&[q(1), q(1)]))
        .scale(&q(-1))
        .to_bivar();
    let relation = ImplicitRelation {
        terms: vec![
            ImplicitTerm {
                shift: 1,
                deriv_order: 0,
                coeff: leader,
            },
            ImplicitTerm {
                shift: 0,
                deriv_order: 0,
                coeff: lower,
            },
        ],
    };
    GeneratingDescription {
        first_index: 0,
        seeds: seeds
            .iter()
            .map(|row| row.iter().map(|c| c.to_string()).collect())
            .collect(),
        relations: vec![DescribedRelation::new(
            &relation,
            0,
            TRANSPORTED_IDENTITY_EVIDENCE,
            "Q",
        )],
        exceptional_index_conditions: Vec::new(),
        boundary_conditions: Vec::new(),
    }
}

#[test]
fn invalid_seeds_and_singular_leaders() {
    let plain = pascal_description(None, &[&["1"]]);
    let assessment = assess_generating_description(&plain).unwrap();
    assert_eq!(assessment.status, EQUIVALENT_DESCRIPTION);
    assert_eq!(assessment.uniqueness.grading.as_deref(), Some("degree"));

    let invalid = pascal_description(None, &[&["1"], &["1", "2"]]);
    let assessment = assess_generating_description(&invalid).unwrap();
    assert_eq!(assessment.status, INCONSISTENT_SEEDS);
    assert!(!assessment.seeds_consistent);

    // (n - 2) Q_(n+1) = (n - 2)(1 + x) Q_n: row 3 is not determined.
    let singular = pascal_description(Some(2), &[&["1"]]);
    let assessment = assess_generating_description(&singular).unwrap();
    assert_eq!(assessment.uniqueness.status, "not_certified");
    assert!(assessment
        .uniqueness
        .reasons
        .iter()
        .any(|r| r.contains("singular leader") && r.contains("[2]")));
    // Seeds through index 3 cover the singular base.
    let covered = pascal_description(
        Some(2),
        &[&["1"], &["1", "1"], &["1", "2", "1"], &["1", "3", "3", "1"]],
    );
    let assessment = assess_generating_description(&covered).unwrap();
    assert_eq!(assessment.status, EQUIVALENT_DESCRIPTION);

    let empty = GeneratingDescription {
        seeds: Vec::new(),
        ..plain
    };
    assert!(assess_generating_description(&empty).is_err());

    // A source rule without enough seeds is rejected.
    let bad = json!({
        "source_recurrence": {"schema": "polytool.recurrence.v1", "first_index": 0,
            "initial_polynomials": [],
            "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                "coeff": {"coeffs": [["1", "1"]]}}]}},
        "transformation": {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
    });
    assert!(implicit_transform_recurrence(&request(bad)).is_err());
}

#[test]
fn registry_maps_and_the_coordinate_report() {
    let diagonal = RecurrenceTransformation::Diagonal {
        weight: DiagonalWeight::InverseBinomial { m: 2 },
    };
    let row = vec![q(1), q(6), q(12)];
    let image = diagonal.map_row(4, &row).unwrap();
    assert_eq!(image[0], vec![q(1), q(2), q(2)]);
    assert_eq!(diagonal.inverse_row(4, &image).unwrap(), row);
    let antecedent = RecurrenceTransformation::DerivativeAntecedent {
        g: vec!["1".into(), "0".into(), "2".into()],
    };
    let image = antecedent.map_row(0, &row).unwrap();
    assert_eq!(antecedent.inverse_row(0, &image).unwrap(), row);
    assert!(RecurrenceTransformation::DerivativeAntecedent {
        g: vec!["0".into()]
    }
    .validate()
    .is_err());
    assert_eq!(diagonal.record().kind, "diagonal");

    // The coordinate-change report maps rows and points to implicit transport.
    let report = transform_recurrence(&RecurrenceTransformRequest::from_source(
        diagonal,
        RecurrenceJson::from_recurrence_rational(
            &crate::oeis::by_id("A007318")
                .unwrap()
                .recurrence_export_parts()
                .unwrap()
                .0,
            0,
            &[vec![q(1)]],
            None,
        ),
        5,
    ))
    .unwrap();
    assert!(report.reconstruction_verified);
    assert!(report
        .transport_error
        .unwrap()
        .contains("transport_implicit"));
    assert_eq!(report.transformed_states[2][0], vec!["1", "2/3", "1/6"]);
}
