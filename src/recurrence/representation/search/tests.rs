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

#[test]
fn from_oeis_request_equals_the_json_request_and_finds_the_pascal_rule() {
    // `from_oeis` is documented as the request for a bundled family with the given
    // maps and otherwise default fields: it must equal the deserialized minimal JSON
    // request, and running it finds the factorial-divided Pascal rows.
    let factorial = RecurrenceTransformation::Diagonal {
        weight: DiagonalWeight::FactorialPower { power: 1 },
    };
    let built = RepresentationSearchRequest::from_oeis("A007318", vec![factorial]);
    let parsed: RepresentationSearchRequest = serde_json::from_value(serde_json::json!({
        "oeis_id": "A007318",
        "transformations": [
            {"kind": "diagonal", "weight": {"kind": "factorial_power", "power": 1}}
        ]
    }))
    .unwrap();
    assert_eq!(built, parsed);
    assert_eq!(built.oeis_id.as_deref(), Some("A007318"));
    assert!(built.source_recurrence.is_none() && built.rows.is_none());
    assert!(built.fits && built.certified_spaces);

    let mut limited = built.clone();
    limited.limits.max_depth = 1;
    let report = search_representations(&limited).unwrap();
    verify_representation_search(&report).unwrap();
    assert_eq!(report.source.kind, "oeis");
    // The factorial-divided rows of Pascal's triangle (A007318): Q_3 = 1 + 3x + 3x^2/2 + x^3/6.
    let node = report
        .representations
        .iter()
        .find(|node| node.chain == [0usize])
        .expect("one-step representation");
    assert_eq!(node.states[3][0], vec!["1", "3", "3/2", "1/6"]);
    assert!(node.reconstruction_verified);
}

fn decoded_states(node: &RepresentationNode) -> Vec<Vec<Vec<Q>>> {
    node.states
        .iter()
        .map(|state| {
            state
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|c| parse_rational_coeff(c).unwrap())
                        .collect()
                })
                .collect()
        })
        .collect()
}

#[test]
fn vector_description_replay_failures_are_verification_failures() {
    let schedule = serde_json::json!({"first_index": 1, "first_degree": 1, "step": 1});
    let request: RepresentationSearchRequest = serde_json::from_value(serde_json::json!({
        "oeis_id": "A008292",
        "row_count": 6,
        "transformations": [{"kind": "symmetric_decomposition", "schedule": schedule}],
        "limits": {"max_depth": 1}
    }))
    .unwrap();
    let report = search_representations(&request).unwrap();
    let node = report
        .representations
        .iter()
        .find(|n| n.chain == [0])
        .unwrap();
    let item = &node.descriptions[0];
    let states = decoded_states(node);
    let origin = node.row_origin;
    check_vector_description(item, &states, origin, "label").unwrap();

    // Each way the replay can disagree with the rows is a failed verification.
    let failed = |item: &RepresentationDescription, states: &[Vec<Vec<Q>>], needle: &str| {
        let error = check_vector_description(item, states, origin, "label").unwrap_err();
        assert!(
            matches!(&error, RecurrenceError::VerificationFailed(message)
                if message.contains(needle)),
            "{error:?}"
        );
    };
    let mut changed = item.clone();
    changed.vector.as_mut().unwrap().initial_state[0][0] = "7".into();
    failed(&changed, &states, "the initial state differs from the rows");
    let mut changed = item.clone();
    changed.vector.as_mut().unwrap().first_index = origin + states.len();
    failed(&changed, &states, "the system starts outside the rows");
    let mut shifted = states.clone();
    let last = shifted.len() - 1;
    shifted[last][0][0] += q(1);
    failed(item, &shifted, "the system fails at n");
    let mut changed = item.clone();
    changed.vector_assessment = None;
    failed(&changed, &states, "the vector assessment does not replay");
    let mut changed = item.clone();
    changed.evidence = FINITE_FIT_EVIDENCE.into();
    failed(&changed, &states, "the label does not match the assessment");
}
