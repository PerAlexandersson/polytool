use super::*;

fn q(n: i64) -> Q {
    Q::from_integer(BigInt::from(n))
}

fn constant_condition() -> CoefficientCondition {
    CoefficientCondition {
        terms: vec![
            ConditionTerm {
                shift: 1,
                coefficient_index: 1,
                coeff: vec![q(1)],
            },
            ConditionTerm {
                shift: 0,
                coefficient_index: 0,
                coeff: vec![q(0), q(-1)],
            },
        ],
    }
}

#[test]
fn conditions_are_carried_through_implicit_maps() {
    // [x^1] P[n+1] - n [x^0] P[n] = 0.
    let condition = constant_condition();
    // Q = T P with lambda = binom(k+1, 1): [x^k] P = [x^k] Q / (k + 1).
    let binomial = RecurrenceTransformation::Diagonal {
        weight: DiagonalWeight::Binomial { m: 1 },
    };
    assert_eq!(
        transport_condition(&binomial, &condition)
            .unwrap()
            .text("Q"),
        "[x^1] Q[n+1] - 2n [x^0] Q[n] = 0"
    );
    // P = D(x Q): [x^e] P = (e + 1) [x^e] Q.
    let antecedent = RecurrenceTransformation::DerivativeAntecedent {
        g: vec!["0".into(), "1".into()],
    };
    assert_eq!(
        transport_condition(&antecedent, &condition)
            .unwrap()
            .text("Q"),
        "[x^1] Q[n+1] - 1/2n [x^0] Q[n] = 0"
    );
    // Q = D P forgets [x^0] P.
    let derivative = RecurrenceTransformation::Derivative { order: 1 };
    assert!(transport_condition(&derivative, &condition).is_none());
    let reversal = RecurrenceTransformation::Reversal {
        schedule: DegreeSchedule::new(0, 0, 1),
    };
    assert!(transport_condition(&reversal, &condition).is_none());
}

fn pascal_source(rule: &Recurrence) -> RecurrenceJson {
    RecurrenceJson::from_recurrence_rational(rule, 0, &[vec![q(1)]], None)
}

fn bivar(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&c| q(c)).collect())
            .collect(),
    )
}

fn gamma_step(rule: &Recurrence, source: usize) -> RepresentationStep {
    let map = RecurrenceTransformation::Gamma {
        schedule: DegreeSchedule::new(0, 0, 1),
    };
    let request = RecurrenceTransformRequest::from_source(map, pascal_source(rule), 8);
    RepresentationStep {
        source_description: Some(source),
        transform: Some(transform_recurrence(&request).unwrap()),
        implicit: None,
        implicit_relation: None,
        vector_transport: None,
        error: None,
    }
}

#[test]
fn gamma_systems_are_certified_only_with_a_certified_symmetry() {
    // (1 + x)^n: an invariant rule and a rank-two rule for the same rows.
    let invariant = Recurrence {
        terms: vec![RecurrenceTerm::new(1, 0, bivar(&[&[1, 1]]))],
        denominator: None,
        inhomogeneous: None,
    };
    let rank_two = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, bivar(&[&[1, 2], &[0, -1]])),
            RecurrenceTerm::new(1, 1, bivar(&[&[0, 1, 1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let certifying = gamma_step(&invariant, 0);
    let coupled = gamma_step(&rank_two, 1);
    let report = coupled.transform.as_ref().unwrap();
    assert_eq!(report.conditions[0].status, CHECKED_SUPPLIED_ROWS);
    assert_eq!(
        certifying.transform.as_ref().unwrap().conditions[0].status,
        CERTIFIED_ALL_INDICES
    );
    let transported = report.transported.as_ref().unwrap();
    let item = || RepresentationDescription {
        provenance: "steps[1].transform.transported".into(),
        ..gamma_system_item(
            report,
            transported,
            transported.system.as_ref().unwrap(),
            false,
        )
        .unwrap()
    };
    let mut alone = vec![item()];
    upgrade_gamma_systems(std::slice::from_ref(&coupled), &mut alone);
    assert_eq!(alone[0].label, FINITE_FIT_LABEL);
    assert_eq!(alone[0].evidence, CHECKED_SUPPLIED_ROWS);
    // Another step of the edge certifies the symmetry of the rows.
    let steps = vec![certifying, coupled.clone()];
    assert_eq!(symmetry_certifying_step(&steps), Some(0));
    let mut items = vec![item()];
    upgrade_gamma_systems(&steps, &mut items);
    assert_eq!(items[0].evidence, TRANSPORTED_IDENTITY_EVIDENCE);
    assert_eq!(items[0].label, EQUIVALENT_DESCRIPTION_LABEL);
    assert!(items[0].note.as_ref().unwrap().contains("steps[0]"));
    let description = items[0].description.as_ref().unwrap();
    assert!(description
        .relations
        .iter()
        .all(|r| r.evidence == TRANSPORTED_IDENTITY_EVIDENCE));
    assert_eq!(
        items[0].assessment.as_ref().unwrap().status,
        EQUIVALENT_DESCRIPTION
    );
    // A step without a certified source certifies nothing.
    let mut uncertified = steps.clone();
    uncertified[0].source_description = None;
    assert_eq!(symmetry_certifying_step(&uncertified), None);
}

#[test]
fn degree_bounds_compare_on_the_whole_domain() {
    let bound = DegreeSchedule::new(1, 1, 1);
    assert!(bound_at_most(&bound, &DegreeSchedule::new(1, 1, 1), 1));
    assert!(bound_at_most(&bound, &DegreeSchedule::new(0, 0, 2), 1));
    // Smaller at the start, or a smaller slope on an unbounded domain.
    assert!(!bound_at_most(&bound, &DegreeSchedule::new(1, 0, 1), 1));
    assert!(!bound_at_most(&bound, &DegreeSchedule::new(1, 5, 0), 1));
    // ... but fine on a bounded domain that ends in time.
    assert!(bound_at_most(
        &bound,
        &DegreeSchedule::new(1, 5, 0).with_last_index(5),
        1
    ));
    // A bound known only from a later index proves nothing earlier.
    assert!(!bound_at_most(
        &DegreeSchedule::new(2, 0, 1),
        &DegreeSchedule::new(0, 9, 1),
        1
    ));
}
