//! Unit tests of the Veronese-section lifts through the crate-internal entry
//! points (the candidate `G` or the factorization supplied directly).

use super::veronese_algebra::*;
use super::veronese_search::*;
use super::*;
use crate::factorization::MonicBivariateFactor;

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn row(values: &[i64]) -> Vec<BigRational> {
    values.iter().map(|&value| q(value)).collect()
}

fn rows(values: &[&[i64]]) -> Vec<Vec<BigRational>> {
    values.iter().map(|values| row(values)).collect()
}

/// `Q_n = sum_k c_k(y) Q_(n-k)` from `(k, ascending coefficients of c_k)`.
fn recurrence(terms: &[(usize, &[i64])]) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, coefficients)| {
                RecurrenceTerm::new(
                    *offset,
                    0,
                    BivarPoly::from_coefficients_normalized(vec![row(coefficients)]),
                )
            })
            .collect(),
        denominator: None,
        inhomogeneous: None,
    }
}

fn factor(coefficients: &[&[i64]]) -> MonicBivariateFactor {
    MonicBivariateFactor {
        coefficients: rows(coefficients),
        multiplicity: 1,
    }
}

fn candidate(coefficients: &[&[i64]]) -> BivariatePolynomial {
    bivariate_from_coefficients(&rows(coefficients))
}

/// `Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n`.
fn worked_recurrence() -> Recurrence {
    recurrence(&[(1, &[2]), (2, &[-1, 1])])
}

fn prepared(recurrence: &Recurrence, seeds: &[&[i64]]) -> PreparedVeroneseInput {
    prepare_veronese_input(recurrence, &rows(seeds), 0).unwrap()
}

fn search(
    prepared: &PreparedVeroneseInput,
    modulus: usize,
    factors: &[MonicBivariateFactor],
    options: &VeroneseSectionLiftOptions,
) -> VeroneseSectionLiftReport {
    let mut accumulator = LiftAccumulator::default();
    search_modulus_with_factors(prepared, modulus, factors, options, &mut accumulator).unwrap();
    accumulator.into_report(prepared)
}

fn find_lift(evaluation: &CandidateEvaluation, residue: usize) -> Option<&VeroneseSectionLift> {
    evaluation.lifts.iter().find(|lift| lift.residue == residue)
}

fn binomial(n: u64, k: u64) -> i64 {
    if k > n {
        return 0;
    }
    (0..k).fold(1i64, |acc, i| acc * (n - i) as i64 / (i + 1) as i64)
}

/// `V_(r, j)((a + b x)^n)` for `n` in `0..count`.
fn product_sections(
    a: i64,
    b: i64,
    modulus: u64,
    residue: u64,
    count: u64,
) -> Vec<Vec<BigRational>> {
    (0..count)
        .map(|n| {
            let mut coefficients = Vec::new();
            let mut power = residue;
            while power <= n {
                let value = binomial(n, power) * a.pow((n - power) as u32) * b.pow(power as u32);
                coefficients.push(q(value));
                power += modulus;
            }
            if coefficients.is_empty() {
                coefficients.push(q(0));
            }
            coefficients
        })
        .collect()
}

fn assert_generates(
    lift: &VeroneseSectionLift,
    recurrence: &Recurrence,
    seeds: &[Vec<BigRational>],
) {
    let count = 10;
    let parent = lift
        .parent_recurrence
        .generate_rows_rational(&lift.parent_initial_rows, lift.parent_first_index, count)
        .unwrap();
    let input = recurrence.generate_rows_rational(seeds, 0, count).unwrap();
    for (parent_row, input_row) in parent.iter().zip(&input) {
        assert_eq!(
            coefficient_section(parent_row, lift.modulus, lift.residue),
            crate::Polynomial::new(input_row.clone())
        );
    }
}

#[test]
fn worked_example_even_sections_of_one_plus_x() {
    let input = prepared(&worked_recurrence(), &[&[1], &[1]]);
    assert_eq!(
        bivariate_to_coefficients(&input.characteristic),
        rows(&[&[1, -1], &[-2], &[1]])
    );
    let evaluation =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, -1], &[1]]), 12)
            .unwrap();
    let lift = find_lift(&evaluation, 0).expect("residue 0 lift");
    assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
    assert_eq!(lift.status, VeroneseSectionLiftStatus::ProvedProductParent);
    assert_eq!(
        lift.section_characteristic_polynomial,
        rows(&[&[1, -1], &[-2], &[1]])
    );
    assert_eq!(lift.cofactor, rows(&[&[1]]));
    assert_eq!(lift.parent_recurrence.to_string(), "P(n) = (1 + t) P(n-1)");
    assert_eq!(lift.parent_row_check.rows_checked, 12);
    // Residue 1 with the parent x (1 + x)^n is a lift too.
    let odd = find_lift(&evaluation, 1).expect("residue 1 lift");
    assert_eq!(odd.parent_initial_rows, rows(&[&[0, 1]]));
    assert_eq!(odd.status, VeroneseSectionLiftStatus::ProvedProductParent);
    assert_generates(lift, &worked_recurrence(), &rows(&[&[1], &[1]]));

    // The conjugate L - 1 + x is listed with a failing status.
    let conjugate =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, 1], &[1]]), 12)
            .unwrap();
    let lift = find_lift(&conjugate, 0).expect("conjugate lift");
    assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
    assert_eq!(
        lift.status,
        VeroneseSectionLiftStatus::ParentRowMixedSigns { index: 1 }
    );
    assert_eq!(lift.parent_row_check.first_mixed_sign_row_index, Some(1));
    assert_eq!(lift.parent_row_check.first_not_real_rooted_row_index, None);
}

#[test]
fn worked_example_search_orders_and_serializes_lifts() {
    let input = prepared(&worked_recurrence(), &[&[1], &[1]]);
    // chi(L, x^2) = (L - 1 + x)(L - 1 - x).
    let factors = [factor(&[&[-1, 1], &[1]]), factor(&[&[-1, -1], &[1]])];
    let report = search(&input, 2, &factors, &VeroneseSectionLiftOptions::default());
    let summary: Vec<_> = report
        .lifts
        .iter()
        .map(|lift| {
            (
                format_characteristic_polynomial(&lift.parent_characteristic_polynomial, "x"),
                lift.residue,
                lift.status.label(),
            )
        })
        .collect();
    assert_eq!(
        summary,
        [
            ("L - 1 - x".to_string(), 0, "proved_product_parent"),
            ("L - 1 - x".to_string(), 1, "proved_product_parent"),
            ("L - 1 + x".to_string(), 0, "parent_row_mixed_signs"),
            ("L - 1 + x".to_string(), 1, "parent_row_mixed_signs"),
        ]
    );
    assert!(report.rejected_candidates.is_empty());
    assert!(!report.exhausted_budget);

    let json = serde_json::to_value(&report).unwrap();
    assert_eq!(json["schema"], VERONESE_SECTION_LIFT_SCHEMA);
    assert_eq!(
        json["input_characteristic_polynomial_text"],
        "L^2 - 2L + 1 - y"
    );
    assert_eq!(json["lifts"][0]["status"], "proved_product_parent");
    assert_eq!(
        json["lifts"][0]["parent_characteristic_polynomial"][0][1],
        "-1"
    );
    assert_eq!(
        json["lifts"][0]["parent_recurrence"]["schema"],
        "polytool.recurrence.v1"
    );
    assert_eq!(json["lifts"][2]["failing_row_index"], 1);
    assert_eq!(
        json["lifts"][0]["failing_row_index"],
        serde_json::Value::Null
    );
}

#[test]
fn odd_section_of_one_plus_x() {
    let input = prepared(&worked_recurrence(), &[&[0], &[1]]);
    let evaluation =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, -1], &[1]]), 12)
            .unwrap();
    let lift = find_lift(&evaluation, 1).expect("odd lift");
    assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
    assert_eq!(lift.status, VeroneseSectionLiftStatus::ProvedProductParent);
    assert_eq!(
        evaluation.rejections,
        [VeroneseSectionLiftRejection {
            modulus: 2,
            residue: Some(0),
            parent_characteristic_polynomial: rows(&[&[-1, -1], &[1]]),
            reason: VeroneseSectionLiftRejectionReason::SeedsNotPolynomial,
        }]
    );
}

#[test]
fn modulus_three_sections() {
    // Q_(n+3) = 3 Q_(n+2) - 3 Q_(n+1) + (1 + y) Q_n.
    let rule = recurrence(&[(1, &[3]), (2, &[-3]), (3, &[1, 1])]);
    let parent = candidate(&[&[-1, -1], &[1]]);
    for residue in 0..3u64 {
        let seeds = product_sections(1, 1, 3, residue, 3);
        let input = prepare_veronese_input(&rule, &seeds, 0).unwrap();
        let evaluation = evaluate_veronese_section_lift_candidate(&input, 3, &parent, 12).unwrap();
        let lift = find_lift(&evaluation, residue as usize).expect("modulus-3 lift");
        assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
        assert_eq!(lift.status, VeroneseSectionLiftStatus::ProvedProductParent);
        assert_eq!(
            lift.section_characteristic_polynomial,
            rows(&[&[-1, -1], &[3], &[-3], &[1]])
        );
        assert_generates(lift, &rule, &seeds);
    }
    // The full search with the factorization of (L - 1)^3 - x^3.
    let seeds = product_sections(1, 1, 3, 0, 3);
    let input = prepare_veronese_input(&rule, &seeds, 0).unwrap();
    let factors = [
        factor(&[&[-1, -1], &[1]]),
        factor(&[&[1, -1, 1], &[-2, 1], &[1]]),
    ];
    let report = search(&input, 3, &factors, &VeroneseSectionLiftOptions::default());
    assert_eq!(report.lifts[0].residue, 0);
    assert_eq!(report.lifts[0].parent_initial_rows, rows(&[&[1]]));
    assert!(report.lifts[0].status.is_proved());
}

#[test]
fn sections_of_two_plus_three_x() {
    // (L - 2)^2 - 9y: Q_(n+2) = 4 Q_(n+1) + (9y - 4) Q_n.
    let rule = recurrence(&[(1, &[4]), (2, &[-4, 9])]);
    let parent = candidate(&[&[-2, -3], &[1]]);
    for residue in 0..2u64 {
        let seeds = product_sections(2, 3, 2, residue, 2);
        let input = prepare_veronese_input(&rule, &seeds, 0).unwrap();
        let evaluation = evaluate_veronese_section_lift_candidate(&input, 2, &parent, 12).unwrap();
        let lift = find_lift(&evaluation, residue as usize).expect("lift");
        assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
        assert_eq!(lift.status, VeroneseSectionLiftStatus::ProvedProductParent);
        assert_generates(lift, &rule, &seeds);
    }
}

#[test]
fn modulus_four_and_an_order_two_parent_at_modulus_two() {
    // (L - 1)^4 - y.
    let rule = recurrence(&[(1, &[4]), (2, &[-6]), (3, &[4]), (4, &[-1, 1])]);
    let seeds = product_sections(1, 1, 4, 0, 4);
    let input = prepare_veronese_input(&rule, &seeds, 0).unwrap();
    let evaluation =
        evaluate_veronese_section_lift_candidate(&input, 4, &candidate(&[&[-1, -1], &[1]]), 12)
            .unwrap();
    let lift = find_lift(&evaluation, 0).expect("modulus-4 lift");
    assert_eq!(lift.parent_initial_rows, rows(&[&[1]]));
    assert!(lift.status.is_proved());

    // At r = 2, chi(L, x^2) = ((L - 1)^2 - x)((L - 1)^2 + x): the parent
    // V_(2,0)((1 + t)^n)(x) of order two is found, with an open status.
    let factors = [
        factor(&[&[1, -1], &[-2], &[1]]),
        factor(&[&[1, 1], &[-2], &[1]]),
    ];
    let report = search(&input, 2, &factors, &VeroneseSectionLiftOptions::default());
    let open = &report.lifts[0];
    assert_eq!(
        open.parent_characteristic_polynomial,
        rows(&[&[1, -1], &[-2], &[1]])
    );
    assert_eq!(open.residue, 0);
    assert_eq!(open.parent_initial_rows, rows(&[&[1], &[1]]));
    assert_eq!(
        open.status,
        VeroneseSectionLiftStatus::ParentRealRootednessOpen
    );
    assert_eq!(open.cofactor, rows(&[&[1]]));
    assert_generates(open, &rule, &seeds);
}

#[test]
fn even_sections_of_fibonacci_polynomials() {
    // chi = L^4 - 2L^3 + L^2 - y.
    let rule = recurrence(&[(1, &[2]), (2, &[-1]), (4, &[0, 1])]);
    let input = prepared(&rule, &[&[1], &[1], &[1], &[1]]);
    let generated = rule
        .generate_rows_rational(&rows(&[&[1], &[1], &[1], &[1]]), 0, 6)
        .unwrap();
    assert_eq!(generated[4], row(&[1, 1]));
    assert_eq!(generated[5], row(&[1, 3]));
    let evaluation = evaluate_veronese_section_lift_candidate(
        &input,
        2,
        &candidate(&[&[0, -1], &[-1], &[1]]),
        12,
    )
    .unwrap();
    let lift = find_lift(&evaluation, 0).expect("Fibonacci lift");
    assert_eq!(lift.parent_initial_rows, rows(&[&[1], &[1]]));
    assert_eq!(
        lift.status,
        VeroneseSectionLiftStatus::ParentRealRootednessOpen
    );
    assert_eq!(lift.parent_row_check.rows_checked, 12);
    assert_eq!(lift.parent_row_check.first_not_real_rooted_row_index, None);
    assert_eq!(lift.parent_row_check.first_mixed_sign_row_index, None);
    assert_eq!(lift.parent_row_check.first_zero_row_index, None);
    assert_eq!(
        format_characteristic_polynomial(&lift.section_characteristic_polynomial, "y"),
        "L^4 - 2L^3 + L^2 - y"
    );
    assert_generates(lift, &rule, &rows(&[&[1], &[1], &[1], &[1]]));

    let factors = [
        factor(&[&[0, -1], &[-1], &[1]]),
        factor(&[&[0, 1], &[-1], &[1]]),
    ];
    let report = search(&input, 2, &factors, &VeroneseSectionLiftOptions::default());
    assert_eq!(
        report.lifts[0].parent_characteristic_polynomial,
        rows(&[&[0, -1], &[-1], &[1]])
    );
    assert_eq!(
        report.lifts[0].status.label(),
        "parent_real_rootedness_open"
    );
    assert!(report.lifts.last().unwrap().status.is_failing());
}

#[test]
fn non_minimal_input_recurrence() {
    // chi = (L - 2)((L - 1)^2 - y) = L^3 - 4L^2 + (5 - y) L + 2y - 2.
    let rule = recurrence(&[(1, &[4]), (2, &[-5, 1]), (3, &[2, -2])]);
    let input = prepared(&rule, &[&[1], &[1], &[1, 1]]);
    let factors = [
        factor(&[&[-2], &[1]]),
        factor(&[&[-1, -1], &[1]]),
        factor(&[&[-1, 1], &[1]]),
    ];
    let report = search(&input, 2, &factors, &VeroneseSectionLiftOptions::default());
    let best = &report.lifts[0];
    assert_eq!(
        best.parent_characteristic_polynomial,
        rows(&[&[-1, -1], &[1]])
    );
    assert_eq!(best.residue, 0);
    assert_eq!(best.parent_initial_rows, rows(&[&[1]]));
    assert!(best.status.is_proved());
    assert_eq!(best.cofactor, rows(&[&[-2], &[1]]));
    // L - 2 alone fails the norm test.
    assert!(report.rejected_candidates.iter().any(|rejection| {
        rejection.parent_characteristic_polynomial == rows(&[&[-2], &[1]])
            && rejection.reason
                == VeroneseSectionLiftRejectionReason::NormDoesNotDivideCharacteristicPolynomial
    }));
}

#[test]
fn near_misses() {
    // Changed seeds Q_0 = 1, Q_1 = 2: residue 0 needs the seed section 1/y.
    let input = prepared(&worked_recurrence(), &[&[1], &[2]]);
    let parent = candidate(&[&[-1, -1], &[1]]);
    let evaluation = evaluate_veronese_section_lift_candidate(&input, 2, &parent, 12).unwrap();
    assert!(find_lift(&evaluation, 0).is_none());
    assert!(evaluation
        .rejections
        .iter()
        .any(|rejection| rejection.residue == Some(0)
            && rejection.reason == VeroneseSectionLiftRejectionReason::SeedsNotPolynomial));
    // These rows are the odd sections of (1 + x)^(n+1): a different parent instance.
    let odd = find_lift(&evaluation, 1).expect("odd lift of the shifted parent");
    assert_eq!(odd.parent_initial_rows, rows(&[&[1, 1]]));

    // Changed coefficient: (L - 1)^2 + 1 - y; the norm of L - 1 - x does not divide.
    let changed = recurrence(&[(1, &[2]), (2, &[-2, 1])]);
    let input = prepared(&changed, &[&[1], &[1]]);
    let evaluation = evaluate_veronese_section_lift_candidate(&input, 2, &parent, 12).unwrap();
    assert!(evaluation.lifts.is_empty());
    assert_eq!(
        evaluation.rejections[0].reason,
        VeroneseSectionLiftRejectionReason::NormDoesNotDivideCharacteristicPolynomial
    );
    assert_eq!(evaluation.rejections[0].residue, None);
    // Its x^2-substitution (L - 1)^2 + 1 - x^2 is irreducible of L-degree 2 > d / r.
    let report = search(
        &input,
        2,
        &[factor(&[&[2, 0, -1], &[-2], &[1]])],
        &VeroneseSectionLiftOptions::default(),
    );
    assert!(report.lifts.is_empty() && report.rejected_candidates.is_empty());

    // An exhausted candidate budget.
    let input = prepared(&worked_recurrence(), &[&[1], &[1]]);
    let factors = [factor(&[&[-1, 1], &[1]]), factor(&[&[-1, -1], &[1]])];
    let options = VeroneseSectionLiftOptions {
        max_parent_candidates: 0,
        ..Default::default()
    };
    let report = search(&input, 2, &factors, &options);
    assert!(report.exhausted_budget && report.lifts.is_empty());
    let options = VeroneseSectionLiftOptions {
        max_parent_candidates: 1,
        ..Default::default()
    };
    let report = search(&input, 2, &factors, &options);
    assert!(report.exhausted_budget);
    assert_eq!(report.lifts.len(), 2);
}

#[test]
fn row_mismatch_and_wrong_factorizations() {
    // Seeds with a supplied third row that the recurrence reproduces, but a
    // candidate whose overdetermined rows disagree: chi = (L - 2)((L-1)^2 - y)
    // with rows that are not a section of (1 + x)^n.
    let rule = recurrence(&[(1, &[4]), (2, &[-5, 1]), (3, &[2, -2])]);
    let input = prepared(&rule, &[&[1], &[1], &[2]]);
    let evaluation =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, -1], &[1]]), 12)
            .unwrap();
    assert!(evaluation
        .rejections
        .iter()
        .any(|rejection| rejection.residue == Some(0)
            && rejection.reason == VeroneseSectionLiftRejectionReason::RowMismatch { index: 2 }));

    // Factors that do not multiply back to chi(L, x^r) are an internal error.
    let input = prepared(&worked_recurrence(), &[&[1], &[1]]);
    let mut accumulator = LiftAccumulator::default();
    let error = search_modulus_with_factors(
        &input,
        2,
        &[factor(&[&[-1, -1], &[1]])],
        &VeroneseSectionLiftOptions::default(),
        &mut accumulator,
    )
    .unwrap_err();
    assert!(matches!(error, RecurrenceError::VerificationFailed(_)));
}

#[test]
fn unsupported_and_invalid_inputs() {
    let seeds = rows(&[&[1], &[1]]);
    let options = VeroneseSectionLiftOptions::default();
    let unsupported = |rule: Recurrence| {
        matches!(
            find_veronese_section_lifts(&rule, &seeds, 0, &options),
            Err(RecurrenceError::Unsupported(_))
        )
    };
    let mut rule = worked_recurrence();
    rule.terms[0].coeff = BivarPoly::from_coefficients_normalized(vec![row(&[2]), row(&[1])]);
    assert!(unsupported(rule), "n-dependent coefficient");
    let mut rule = worked_recurrence();
    rule.terms[0].deriv_order = 1;
    assert!(unsupported(rule), "derivative term");
    let mut rule = worked_recurrence();
    rule.inhomogeneous = Some(BivarPoly::from_coefficients_normalized(vec![row(&[1])]));
    assert!(unsupported(rule), "inhomogeneous term");
    let mut rule = worked_recurrence();
    rule.terms[0].sign = RecurrenceSign::AlternatingN;
    assert!(unsupported(rule), "alternating term");
    let mut rule = worked_recurrence();
    rule.denominator = Some(BivarPoly::from_coefficients_normalized(vec![row(&[1, 1])]));
    assert!(unsupported(rule), "denominator depending on t");
    assert!(unsupported(Recurrence {
        terms: Vec::new(),
        denominator: None,
        inhomogeneous: None,
    }));

    let invalid = |rule: &Recurrence, rows: &[Vec<BigRational>], options| {
        matches!(
            find_veronese_section_lifts(rule, rows, 0, &options),
            Err(RecurrenceError::InvalidInput(_))
        )
    };
    let rule = worked_recurrence();
    assert!(
        invalid(&rule, &rows(&[&[1]]), options.clone()),
        "too few rows"
    );
    assert!(
        invalid(&rule, &rows(&[&[1], &[1], &[1, 2]]), options.clone()),
        "supplied row disagrees"
    );
    for max_modulus in [1, MAX_VERONESE_SECTION_MODULUS + 1] {
        let options = VeroneseSectionLiftOptions {
            max_modulus,
            ..Default::default()
        };
        assert!(invalid(&rule, &seeds, options), "max_modulus {max_modulus}");
    }

    // A constant denominator is divided out.
    let mut scaled = recurrence(&[(1, &[4]), (2, &[-2, 2])]);
    scaled.denominator = Some(BivarPoly::from_coefficients_normalized(vec![row(&[2])]));
    let input = prepared(&scaled, &[&[1], &[1]]);
    assert_eq!(
        bivariate_to_coefficients(&input.characteristic),
        rows(&[&[1, -1], &[-2], &[1]])
    );
}

#[test]
fn parent_recurrence_json_round_trip() {
    let input = prepared(
        &recurrence(&[(1, &[2]), (2, &[-1]), (4, &[0, 1])]),
        &[&[1], &[1], &[1], &[1]],
    );
    let evaluation = evaluate_veronese_section_lift_candidate(
        &input,
        2,
        &candidate(&[&[0, -1], &[-1], &[1]]),
        6,
    )
    .unwrap();
    let lift = find_lift(&evaluation, 0).unwrap();
    let encoded = serde_json::to_string(&lift.parent_recurrence_json()).unwrap();
    let decoded: RecurrenceJson = serde_json::from_str(&encoded).unwrap();
    let (rule, first_index, seeds) = decoded.to_recurrence_parts().unwrap();
    let replayed = rule.generate_rows_rational(&seeds, first_index, 8).unwrap();
    let direct = lift
        .parent_recurrence
        .generate_rows_rational(&lift.parent_initial_rows, lift.parent_first_index, 8)
        .unwrap();
    assert_eq!(replayed, direct);
    // Fibonacci polynomials F_n = F_(n-1) + x F_(n-2).
    assert_eq!(direct[4], row(&[1, 3, 1]));
}

#[test]
fn lift_summaries_of_proved_open_and_failing_parents() {
    let input = prepared(&worked_recurrence(), &[&[1], &[1]]);
    let evaluation =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, -1], &[1]]), 12)
            .unwrap();
    let proved = find_lift(&evaluation, 0).unwrap().lift_summary();
    assert_eq!(proved.verdict, LiftVerdict::Proved);
    assert!(proved.is_consistent());
    assert_eq!(
        proved.operation,
        LiftOperation::VeroneseSection {
            modulus: 2,
            residue: 0
        }
    );
    assert_eq!(
        proved.identity,
        LiftIdentityStatus::Certified {
            method: "identifying_rows".into()
        }
    );
    assert!(matches!(
        &proved.parent_proof,
        LiftParentProof::Proved { route, .. } if route == "periodic_product"
    ));
    assert!(proved.hypotheses.is_empty());
    assert_eq!(
        proved
            .root_transfer
            .parent_root_region
            .as_ref()
            .unwrap()
            .to_string(),
        "(-inf, 0]"
    );
    assert_eq!(
        proved.parent.recurrence,
        find_lift(&evaluation, 0).unwrap().parent_recurrence_json()
    );

    // The conjugate parent fails at row 1: refuted, with the row index.
    let conjugate =
        evaluate_veronese_section_lift_candidate(&input, 2, &candidate(&[&[-1, 1], &[1]]), 12)
            .unwrap();
    let failing = find_lift(&conjugate, 0).unwrap().lift_summary();
    assert_eq!(failing.verdict, LiftVerdict::RouteFailed);
    assert!(matches!(
        failing.parent_proof,
        LiftParentProof::Failed { index: Some(1), .. }
    ));

    // The Fibonacci parent is open.
    let rule = recurrence(&[(1, &[2]), (2, &[-1]), (4, &[0, 1])]);
    let fibonacci = prepared(&rule, &[&[1], &[1], &[1], &[1]]);
    let evaluation = evaluate_veronese_section_lift_candidate(
        &fibonacci,
        2,
        &candidate(&[&[0, -1], &[-1], &[1]]),
        12,
    )
    .unwrap();
    let open = find_lift(&evaluation, 0).unwrap();
    let summary = open.lift_summary();
    assert_eq!(summary.verdict, LiftVerdict::Open);
    assert!(matches!(summary.parent_proof, LiftParentProof::Open { .. }));

    // The JSON report carries the summary next to the existing fields.
    let json = serde_json::to_value(open).unwrap();
    assert_eq!(json["status"], "parent_real_rootedness_open");
    assert_eq!(json["lift_summary"]["verdict"], "open");
    assert_eq!(
        json["lift_summary"]["operation"]["kind"],
        "veronese_section"
    );
    assert_eq!(json["lift_summary"]["lean_replay"], "not_attempted");
    let back: LiftSummary = serde_json::from_value(json["lift_summary"].clone()).unwrap();
    assert_eq!(back, summary);
}
