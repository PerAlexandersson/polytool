//! Exact recognition, seed, sign, discovery and replay checks for the
//! positive lag-two Euler theorem; negative lag families stay excluded.

use num_bigint::BigInt;
use num_traits::{One, Zero};
use polytool::recurrence::*;
use polytool::vec_poly;
use serde_json::Value;

type Q = BigRational;

fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

fn polynomial(values: &[i64]) -> Vec<Q> {
    values.iter().copied().map(q).collect()
}

fn at_power(p: &[Q], power: usize) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        p.iter()
            .map(|value| {
                let mut row = vec![q(0); power + 1];
                row[power] = value.clone();
                row
            })
            .collect(),
    )
}

/// Expanded derivative form, with common scalar denominator f.
fn rule(a: &[i64], b: &[i64], c: &[i64], f: &[i64], first: usize) -> Recurrence {
    let a = polynomial(a);
    let b = polynomial(b);
    let c = polynomial(c);
    let degree = vec![-Q::from_integer(BigInt::from(first) + 1), q(1)];
    Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, at_power(&a, 0)),
            RecurrenceTerm::new(
                1,
                0,
                at_power(&vec_poly::add(&a, &vec_poly::mul(&b, &degree)), 1),
            ),
            RecurrenceTerm::new(1, 1, at_power(&b, 1)),
            RecurrenceTerm::new(1, 1, at_power(&vec_poly::neg(&b), 2)),
            RecurrenceTerm::new(2, 0, at_power(&vec_poly::mul(&c, &degree), 1)),
            RecurrenceTerm::new(
                2,
                1,
                at_power(&vec_poly::mul(&c, &vec_poly::sub(&degree, &[q(2)])), 2),
            ),
            RecurrenceTerm::new(2, 2, at_power(&vec_poly::neg(&c), 3)),
        ],
        denominator: Some(at_power(&polynomial(f), 0)),
        inhomogeneous: None,
    }
}

fn seeds() -> Vec<Vec<Q>> {
    vec![vec![q(1)], vec![q(1), q(1)]]
}

fn check(rule: &Recurrence, seeds: &[Vec<Q>]) -> RealRootednessTheoremReport {
    check_lag_two_euler_theorem(
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

fn assert_failed(report: &RealRootednessTheoremReport, name: &str, index: usize) {
    assert_eq!(report.status, RealRootednessTheoremStatus::HypothesisFailed);
    let hypothesis = report.hypotheses.iter().find(|h| h.name == name).unwrap();
    assert_eq!(
        hypothesis.status,
        LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(index))
        }
    );
    report.verify().unwrap();
}

#[test]
fn catalog_targets_have_simple_negative_roots_and_strict_interlacing() {
    for id in [
        "A157147", "A157148", "A157149", "A157150", "A157151", "A157637",
    ] {
        let entry = polytool::oeis::by_id(id).unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let report = check_lag_two_euler_theorem(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            &Default::default(),
        )
        .unwrap();
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        assert_eq!(report.theorem, RealRootednessTheorem::LagTwoEuler);
        assert_eq!(report.start_delay, 0);
        let conclusion = report.conclusion.as_ref().unwrap();
        assert!(conclusion.simple_roots);
        assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
        assert_eq!(conclusion.interlacing_from, Some(1));
        let region = conclusion.root_region.as_ref().unwrap().at(0).unwrap();
        assert!(region.lower.is_none());
        assert_eq!(region.upper, Some(q(0)));
        assert!(!region.upper_closed);
        assert!(report
            .lean_declarations
            .iter()
            .any(|name| name == "RealRooted.lagTwoEuler_recurrence_certificate"));
        report.verify().unwrap();
        let certificate = report.lag_two_euler_certificate().unwrap();
        certificate.verify_operator_identity(&rule, 0).unwrap();
        certificate.verify_seed_steps(0, &seeds).unwrap();
        let record = RecurrenceJson::from_recurrence_rational(&rule, 0, &seeds, None);
        let discovery = discover_recurrence_proofs(&record, &Default::default()).unwrap();
        assert_eq!(discovery.status, "proved");
        assert!(discovery
            .attempts
            .iter()
            .any(|r| r.theorem == RealRootednessTheorem::LagTwoEuler && r.is_proved()));
        discovery.verify().unwrap();
        let parent = prove_parent_real_rootedness(
            &ParentRecurrenceRecord {
                rule: &rule,
                seeds: &seeds,
                first_index: 0,
            },
            0,
            &[],
            true,
        )
        .unwrap();
        assert_eq!(parent.route, Some(VeroneseParentProofRoute::LagTwoEuler));
    }
}

/// Independently generate the negative-sign coefficient triangle, exactly.
fn negative_triangle(m: i64, max_n: usize) -> Vec<Vec<Q>> {
    let mut rows = vec![vec![q(1)], vec![q(1), q(1)]];
    for n in 2..=max_n {
        let mut row = vec![q(0); n + 1];
        for (k, value) in row.iter_mut().enumerate() {
            let get = |n: usize, k: usize| rows[n].get(k).cloned().unwrap_or_else(Q::zero);
            if k > 0 {
                *value += q(m * (n - k) as i64 + 1) * get(n - 1, k - 1);
                *value -= q(m * k as i64 * (n - k) as i64) * get(n - 2, k - 1);
            }
            *value += q(m * k as i64 + 1) * get(n - 1, k);
        }
        rows.push(row);
    }
    rows
}

#[test]
fn negative_lag_families_are_rejected_and_a157152_row_eight_is_not_real_rooted() {
    for m in 1..=5 {
        let rule = rule(&[1], &[m], &[-m], &[1], 0);
        let rows = negative_triangle(m, 8);
        assert_eq!(rule.generate_rows_rational(&seeds(), 0, 9).unwrap(), rows);
        let report = check(&rule, &seeds());
        assert_failed(&report, "lag_coefficient_nonnegative", 2);
        assert!(
            !check_lag_two_euler_theorem(&rule, 0, &seeds(), &Default::default())
                .unwrap()
                .is_proved()
        );
        if m == 1 {
            let integer: Vec<BigInt> = rows[8]
                .iter()
                .map(|v| {
                    assert!(v.denom().is_one());
                    v.to_integer()
                })
                .collect();
            assert_eq!(
                vec_poly::evaluate(&integer, &BigInt::from(-1)),
                BigInt::from(-96)
            );
            assert!(!polytool::is_real_rooted_bigint_coeffs(&integer));
        }
    }
}

#[test]
fn signs_are_decided_on_the_entire_integer_range() {
    assert_failed(
        &check(&rule(&[2], &[1], &[1], &[1], 0), &seeds()),
        "successor_theta_dominates_constant",
        2,
    );
    assert_failed(
        &check(&rule(&[1], &[101, -1], &[1], &[1], 0), &seeds()),
        "successor_theta_dominates_constant",
        101,
    );
    assert_failed(
        &check(&rule(&[1], &[1], &[100, -1], &[1], 0), &seeds()),
        "lag_coefficient_nonnegative",
        101,
    );
    assert_failed(
        &check(&rule(&[100, -1], &[100, -1], &[0], &[1], 0), &seeds()),
        "successor_constant_positive",
        100,
    );
    let singular = check(&rule(&[1], &[1], &[1], &[-100, 1], 0), &seeds());
    assert!(!singular.is_proved());
    assert!(singular
        .hypotheses
        .iter()
        .any(|h| h.name == "left_factor_nonzero"
            && matches!(h.status, LiftHypothesisStatus::Failed { .. })));
    singular.verify().unwrap();
}

#[test]
fn rational_coefficients_zero_lag_boundary_and_nonzero_index_origin() {
    for rule in [
        rule(&[1, 1], &[3, 2], &[0, 1], &[2, 1], 0),
        rule(&[-1], &[-1], &[0], &[-1], 0),
    ] {
        let seeds = vec![vec![q(2)], vec![q(3), q(5)]];
        let report = check(&rule, &seeds);
        assert!(report.is_proved());
        report.verify().unwrap();
        let rows = rule.generate_rows_rational(&seeds, 0, 7).unwrap();
        let extra = check(&rule, &rows);
        assert!(extra.is_proved(), "{:?}", extra.open_hypotheses());
        extra.verify().unwrap();
    }
    let shifted = rule(&[1], &[2], &[1], &[1], 7);
    let report = check_lag_two_euler_theorem(&shifted, 7, &seeds(), &Default::default()).unwrap();
    assert!(report.is_proved());
    assert_eq!(report.conclusion.as_ref().unwrap().first_index, 7);
    report.verify().unwrap();
}

#[test]
fn requested_regions_cannot_strengthen_the_certified_negative_ray() {
    use polytool::real_rootedness::RealInterval;
    let rule = rule(&[1], &[2], &[1], &[1], 0);
    for (upper, expected) in [(q(0), true), (q(1), true), (q(-1), false)] {
        let report = check_lag_two_euler_theorem(
            &rule,
            0,
            &seeds(),
            &RealRootednessTheoremOptions {
                root_region: Some(RealInterval {
                    lower: None,
                    lower_closed: false,
                    upper: Some(upper),
                    upper_closed: false,
                }),
                maximum_start_delay: 0,
            },
        )
        .unwrap();
        assert_eq!(report.is_proved(), expected);
        report.verify().unwrap();
    }
}

#[test]
fn index_function_replay_has_no_preserver_search_degree_bound() {
    let mut rule = rule(&[1], &[2], &[1], &[1], 0);
    // Keep the same sequence, multiplying both sides by 1+n^260.
    // The unreduced rational functions retain this nonvanishing factor.
    const DEGREE: usize = 260;
    for coefficient in rule
        .terms
        .iter_mut()
        .map(|t| &mut t.coeff)
        .chain(rule.denominator.iter_mut())
    {
        let original = &coefficient.coeffs;
        let width = original.iter().map(Vec::len).max().unwrap();
        let mut scaled = vec![vec![q(0); width]; original.len() + DEGREE];
        for (i, row) in original.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                scaled[i][j] += value;
                scaled[i + DEGREE][j] += value;
            }
        }
        *coefficient = BivarPoly::from_coefficients_normalized(scaled);
    }
    let report = check(&rule, &seeds());
    // The sign engine has its own degree budget. An unresolved report must
    // still replay its exact operator certificate, without a preserver cap.
    assert_eq!(report.status, RealRootednessTheoremStatus::Unresolved);
    let certificate = report.lag_two_euler_certificate().unwrap();
    assert_eq!(certificate.a.numerator.len(), DEGREE + 1);
    certificate.verify_operator_identity(&rule, 0).unwrap();
    report.verify().unwrap();
    let mut forged = certificate.clone();
    forged.a.denominator = vec!["0".into()];
    assert!(forged.verify_operator_identity(&rule, 0).is_err());
}

#[test]
fn seed_steps_use_the_lean_inequalities_and_start_delay_replays() {
    let rule = rule(&[1], &[3], &[-3, 1], &[1], 0);
    assert_failed(&check(&rule, &seeds()), "lag_coefficient_nonnegative", 2);
    // At n=2 the generated row is also a valid positive Euler step with
    // b=2,c=0. The rule's own c(n) is nonnegative from n=3 on.
    let report = check_lag_two_euler_theorem(&rule, 0, &seeds(), &Default::default()).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    assert_eq!(report.start_delay, 1);
    assert_eq!(
        report.lag_two_euler_certificate().unwrap().seed_steps.len(),
        1
    );
    report.verify().unwrap();
    let mut wrong = rule.generate_rows_rational(&seeds(), 0, 4).unwrap();
    wrong[3][1] = q(-100);
    assert_failed(&check(&rule, &wrong), "seed_euler_steps", 3);
    for bad in [
        vec![vec![q(-1)], vec![q(1), q(1)]],
        vec![vec![q(1)], vec![q(0), q(1)]],
        vec![vec![q(1)], vec![q(1), q(-1)]],
        vec![vec![q(1)], vec![q(1), q(1), q(1)]],
    ] {
        let report = check(&rule, &bad);
        assert!(!report.is_proved());
        report.verify().unwrap();
    }
}

#[test]
fn nonmatching_operators_and_inhomogeneous_rules_are_not_applicable() {
    let valid = rule(&[1], &[2], &[1], &[1], 0);
    for term in [
        RecurrenceTerm::new(2, 1, at_power(&[q(1)], 2)),
        RecurrenceTerm::new(2, 2, at_power(&[q(1)], 1)),
        RecurrenceTerm::new(3, 0, at_power(&[q(1)], 0)),
        RecurrenceTerm::new(1, 3, at_power(&[q(1)], 3)),
    ] {
        let mut wrong = valid.clone();
        wrong.terms.push(term);
        assert_eq!(
            check(&wrong, &seeds()).status,
            RealRootednessTheoremStatus::NotApplicable
        );
    }
    let mut wrong = valid;
    wrong.inhomogeneous = Some(at_power(&[q(1)], 0));
    assert_eq!(
        check(&wrong, &seeds()).status,
        RealRootednessTheoremStatus::NotApplicable
    );
}

#[test]
fn common_left_monomial_normalization_and_residue_menu_replay() {
    let mut normalized = rule(&[1], &[2], &[1], &[1], 0);
    for coefficient in normalized
        .terms
        .iter_mut()
        .map(|t| &mut t.coeff)
        .chain(normalized.denominator.iter_mut())
    {
        let shifted = coefficient
            .coeffs
            .iter()
            .map(|row| {
                let mut shifted = vec![q(0); 2];
                shifted.extend_from_slice(row);
                shifted
            })
            .collect();
        *coefficient = BivarPoly::from_coefficients_normalized(shifted);
    }
    let report = check(&normalized, &seeds());
    assert!(report.is_proved());
    assert!(report.left_hand_factor_normalization.is_some());
    report.verify().unwrap();

    let mut even = rule(&[1], &[2], &[1], &[1], 0);
    for term in &mut even.terms {
        term.offset *= 2;
        let coefficients = term
            .coeff
            .coeffs
            .iter()
            .enumerate()
            .map(|(i, row)| {
                let scale = Q::new(BigInt::one(), BigInt::from(2).pow(i as u32));
                vec_poly::scale(row, &scale)
            })
            .collect();
        term.coeff = BivarPoly::from_coefficients_normalized(coefficients);
    }
    let seeds = vec![vec![q(1)], vec![q(0)], vec![q(1), q(1)], vec![q(0)]];
    let report =
        check_residue_class_subsequence_theorem(&even, 0, &seeds, &Default::default()).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    let value = serde_json::to_value(&report).unwrap();
    assert!(value.to_string().contains("lag_two_euler"));
    report.verify().unwrap();
}

#[test]
fn cli_discovers_and_replays_the_catalog_rule() {
    use std::io::Write;
    use std::process::{Command, Stdio};
    let run = |arguments: &[&str], input: &[u8]| {
        let mut child = Command::new(env!("CARGO_BIN_EXE_polytool"))
            .args(arguments)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .unwrap();
        child.stdin.take().unwrap().write_all(input).unwrap();
        child.wait_with_output().unwrap()
    };
    let entry = polytool::oeis::by_id("A157637").unwrap();
    let (rule, seeds) = entry.recurrence_parts().unwrap();
    let record = RecurrenceJson::from_recurrence_rational(&rule, 0, &seeds, None);
    let result = run(
        &["recurrence-prove", "--json"],
        &serde_json::to_vec(&record).unwrap(),
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let report: RecurrenceProofDiscoveryReport = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(report.status, "proved");
    assert!(report
        .attempts
        .iter()
        .any(|a| a.theorem == RealRootednessTheorem::LagTwoEuler && a.is_proved()));
    let replay = run(&["recurrence-prove", "--verify"], &result.stdout);
    assert!(
        replay.status.success(),
        "{}",
        String::from_utf8_lossy(&replay.stderr)
    );
}

#[test]
fn certificate_tampering_is_rejected_independently_and_in_report_replay() {
    let rule = rule(&[1], &[2], &[1], &[1], 0);
    let rows = rule.generate_rows_rational(&seeds(), 0, 5).unwrap();
    let report = check(&rule, &rows);
    assert!(report.is_proved());
    report.verify().unwrap();
    let mut certificate = report.lag_two_euler_certificate().unwrap().clone();
    certificate.c.numerator[0] = "-1".into();
    assert!(certificate.verify_operator_identity(&rule, 0).is_err());
    let mut certificate = report.lag_two_euler_certificate().unwrap().clone();
    certificate.seed_steps[0].c = "-1".into();
    assert!(certificate.verify_seed_steps(0, &rows).is_err());
    let value = serde_json::to_value(&report).unwrap();
    for field in ["a", "b", "c"] {
        let mut forged = value.clone();
        forged["certificate"]["lag_two_euler"][field]["numerator"][0] = Value::from("999");
        let forged: RealRootednessTheoremReport = serde_json::from_value(forged).unwrap();
        assert!(forged.verify().is_err());
    }
    let mut forged = value.clone();
    forged["certificate"] = Value::Null;
    assert!(
        serde_json::from_value::<RealRootednessTheoremReport>(forged)
            .unwrap()
            .verify()
            .is_err()
    );
    let mut forged = value;
    forged["lean_declarations"] = serde_json::json!(["RealRooted.wrong"]);
    assert!(
        serde_json::from_value::<RealRootednessTheoremReport>(forged)
            .unwrap()
            .verify()
            .is_err()
    );
}
