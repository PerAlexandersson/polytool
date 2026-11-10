//! The RealRooted row routes `variableDelay` and `subseq`: catalog families
//! proved upstream by these routes, negative cases, tampered certificates
//! and the parent provers.

use num_bigint::BigInt;
use polytool::recurrence::*;

mod common;
use common::{q, rows};

fn b(rows: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(
        rows.iter()
            .map(|r| r.iter().map(|&n| q(n)).collect())
            .collect(),
    )
}

fn options() -> RealRootednessTheoremOptions {
    RealRootednessTheoremOptions::default()
}

/// The catalog family: rule, identifying seeds and first index.
fn catalog(id: &str) -> (Recurrence, Vec<Vec<BigRational>>, usize) {
    let entry = polytool::oeis::by_id(id).unwrap();
    let (rule, rows) = entry.recurrence_parts().unwrap();
    let count = rule.generation_initial_count(entry.recurrence_first_index, rows.len());
    (rule, rows[..count].to_vec(), entry.recurrence_first_index)
}

fn integer(p: &[BigRational]) -> Vec<BigInt> {
    let mut p = p.to_vec();
    while p.len() > 1 && p.last().is_some_and(|c| c == &q(0)) {
        p.pop();
    }
    let common = p.iter().fold(BigInt::from(1), |l, c| {
        num_integer::Integer::lcm(&l, c.denom())
    });
    p.iter()
        .map(|c| c.numer() * (&common / c.denom()))
        .collect()
}

fn hypothesis<'a>(report: &'a RealRootednessTheoremReport, name: &str) -> &'a LiftHypothesis {
    report
        .hypotheses
        .iter()
        .find(|h| h.name == name)
        .unwrap_or_else(|| panic!("no hypothesis {name}"))
}

fn tampered(report: &RealRootednessTheoremReport, edit: impl Fn(&mut serde_json::Value)) -> bool {
    let mut value = serde_json::to_value(report).unwrap();
    edit(&mut value);
    let report: RealRootednessTheoremReport = serde_json::from_value(value).unwrap();
    report.verify().is_err()
}

// ---------------------------------------------------------------------------
// variableDelay
// ---------------------------------------------------------------------------

#[test]
fn a201013_is_a_scaled_fixed_delay_window() {
    let (rule, seeds, first) = catalog("A201013");
    let report = check_variable_delay_theorem(&rule, first, &seeds, &options()).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    assert_eq!(report.lean_route.as_deref(), Some("variableDelay"));
    assert_eq!(
        report.lean_declarations[0],
        "RealRooted.strictInterl_of_scaled_fixedDelay_recurrence"
    );
    let certificate = report.variable_delay_certificate().unwrap();
    assert_eq!(certificate.delay, 3);
    assert!(!certificate.unscaled);
    assert_eq!(certificate.seed_steps.len(), seeds.len() - 1);
    let conclusion = report.conclusion.as_ref().unwrap();
    assert_eq!(conclusion.interlacing, Some(RowInterlacing::Weak));
    assert_eq!(conclusion.interlacing_from, Some(first + 1));
    assert_eq!(
        conclusion.root_region.as_ref().unwrap().to_string(),
        "(-inf, 0)"
    );
    report.verify().unwrap();

    // The conclusion on the first rows, exactly: real-rooted, nonnegative,
    // and the window P_A << P_B for B - A <= 3.
    let generated = rule.generate_rows_rational(&seeds, first, 14).unwrap();
    for (a, row) in generated.iter().enumerate() {
        let row = integer(row);
        assert!(polytool::is_real_rooted_bigint_coeffs(&row));
        assert!(row.iter().all(|c| c >= &BigInt::from(0)));
        for step in 1..=3 {
            if let Some(later) = generated.get(a + step) {
                assert_eq!(
                    polytool::check_weak_interlacing_bigint_coeffs(&row, &integer(later)),
                    Some(true),
                    "P_{} << P_{}",
                    first + a,
                    first + a + step
                );
            }
        }
    }
}

#[test]
fn unscaled_catalog_windows_name_the_unscaled_declaration() {
    // A102547: P_n = P_(n-1) + x P_(n-3) (also the paper-theorem route
    // fixed_delay_predecessor_window).
    let (rule, seeds, first) = catalog("A102547");
    let report = check_variable_delay_theorem(&rule, first, &seeds, &options()).unwrap();
    assert!(report.is_proved(), "{:?}", report.open_hypotheses());
    assert!(report.variable_delay_certificate().unwrap().unscaled);
    assert_eq!(
        report.lean_declarations[0],
        "RealRooted.strictInterl_of_fixedDelay_recurrence"
    );
    report.verify().unwrap();
}

#[test]
fn variable_delay_rejects_other_shapes_and_failing_hypotheses() {
    let delay_rule = |a: &[&[i64]], b_: &[&[i64]], lag: usize| Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(a)),
            RecurrenceTerm::new(lag, 0, b(b_)),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = rows(&[&[1], &[1], &[1]]);
    let check = |rule: &Recurrence, seeds: &[Vec<BigRational>]| {
        let report = check_variable_delay_theorem(rule, 0, seeds, &options()).unwrap();
        report.verify().unwrap();
        report
    };
    // L = 1 is the product route; a constant term at the delayed lag, and
    // an x-term at the previous lag, are other shapes.
    for rule in [
        delay_rule(&[&[1]], &[&[0, 1]], 1),
        delay_rule(&[&[1]], &[&[1, 1]], 2),
        delay_rule(&[&[1, 1]], &[&[0, 1]], 2),
    ] {
        let report = check(&rule, &seeds);
        assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
        assert!(report.variable_delay_certificate().is_none() && report.hypotheses.is_empty());
    }
    // P_n = P_(n-1) - x P_(n-2): the delayed weight is negative (P_2 = 1 - x
    // is real-rooted, but the window fails).
    let negative = check(&delay_rule(&[&[1]], &[&[0, -1]], 2), &seeds[..2]);
    assert_eq!(
        negative.status,
        RealRootednessTheoremStatus::HypothesisFailed
    );
    assert!(matches!(
        hypothesis(&negative, "delayed_coefficient_nonnegative").status,
        LiftHypothesisStatus::Failed { .. }
    ));
    // a(n) = n - 5 changes sign at n = 5; moving the start by up to three
    // rows does not help.
    let sign_change = check(&delay_rule(&[&[-5], &[1]], &[&[0, 1]], 2), &seeds[..2]);
    assert!(!sign_change.is_proved());
    // A nonconstant first row, and a seed that is not a window step.
    let rule = delay_rule(&[&[1]], &[&[0, 1]], 2);
    let first_row = check(&rule, &rows(&[&[1, 1], &[1, 1]]));
    assert!(matches!(
        hypothesis(&first_row, "first_row_constant").status,
        LiftHypothesisStatus::Failed { witness: Some(_) }
    ));
    let bad_seed = check(&rule, &rows(&[&[1], &[1, 0, 1]]));
    assert!(!bad_seed.is_proved());
    assert!(matches!(
        hypothesis(&bad_seed, "seed_window_steps").status,
        LiftHypothesisStatus::Failed { .. }
    ));
}

#[test]
fn tampered_variable_delay_certificates_are_rejected() {
    let (rule, seeds, first) = catalog("A201013");
    let report = check_variable_delay_theorem(&rule, first, &seeds, &options()).unwrap();
    report.verify().unwrap();
    assert!(tampered(&report, |v| {
        v["certificate"]["variable_delay"]["seed_steps"][1]["previous_scalar"] = "5".into();
    }));
    assert!(tampered(&report, |v| v["certificate"]["variable_delay"]
        ["delay"] = 2.into()));
    assert!(tampered(&report, |v| {
        v["certificate"] = serde_json::Value::Null;
    }));
    assert!(tampered(&report, |v| {
        v["conclusion"]["interlacing_from"] = 0.into();
    }));
    assert!(tampered(&report, |v| {
        v["lean_declarations"][0] = "RealRooted.strictInterl_of_fixedDelay_recurrence".into();
    }));
    // The certificate alone replays on the seeds, and fails on other seeds.
    let certificate = report.variable_delay_certificate().unwrap();
    certificate.verify_seed_steps(first, &seeds).unwrap();
    let mut other = seeds.clone();
    other[2] = vec![q(13)];
    assert!(certificate.verify_seed_steps(first, &other).is_err());
}

// ---------------------------------------------------------------------------
// subseq
// ---------------------------------------------------------------------------

/// Upstream `rr_row_splits` proofs whose recurrences have lags `2` and `4`.
const SUBSEQ_IDS: [&str; 8] = [
    "A099605", "A132918", "A137374", "A171608", "A211956", "A223168", "A286100", "A003982",
];

#[test]
fn catalog_families_are_proved_one_residue_class_at_a_time() {
    for id in SUBSEQ_IDS {
        let (rule, seeds, first) = catalog(id);
        let report =
            check_residue_class_subsequence_theorem(&rule, first, &seeds, &options()).unwrap();
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        assert_eq!(report.lean_route.as_deref(), Some("subseq"), "{id}");
        assert_eq!(report.lean_declarations[0], "RealRooted.exists_mul_add_lt");
        let certificate = report.residue_class_certificate().unwrap();
        assert_eq!(certificate.period, 2, "{id}");
        assert_eq!(certificate.classes.len(), 2);
        report.verify().unwrap();

        // Each class family generates exactly the rows of its class.
        let generated = rule.generate_rows_rational(&seeds, first, 24).unwrap();
        for class in &certificate.classes {
            let (class_rule, class_first, class_seeds) =
                class.recurrence.to_recurrence_parts().unwrap();
            assert_eq!(class_first, 0);
            let class_rows = class_rule
                .generate_rows_rational(&class_seeds, 0, 12)
                .unwrap();
            for (m, row) in class_rows.iter().enumerate() {
                let n = class.residue + 2 * m;
                let trim = |p: &[BigRational]| polytool::vec_poly::trim(p);
                assert_eq!(trim(row), trim(&generated[n]), "{id}: class row {m}");
            }
            if class.zero {
                assert!(class.attempts.is_empty() && class.selected.is_none());
            } else {
                assert!(class.attempts[class.selected.unwrap()].is_proved());
            }
        }
    }
}

#[test]
fn zero_classes_are_named_in_the_conclusion() {
    // A132918: P_n = x^2 P_(n-2), P_0 = 1, P_1 = 0 (the catalog rule with
    // first index 0); the odd rows vanish, so no index has every later row
    // nonzero.
    let (rule, seeds, first) = catalog("A132918");
    assert_eq!(first, 0);
    let report = check_residue_class_subsequence_theorem(&rule, first, &seeds, &options()).unwrap();
    assert!(report.is_proved());
    let certificate = report.residue_class_certificate().unwrap();
    assert_eq!(
        certificate.classes.iter().filter(|c| c.zero).count(),
        1,
        "one zero class"
    );
    let conclusion = report.conclusion.as_ref().unwrap();
    assert!(
        conclusion.statement.contains("are zero"),
        "{}",
        conclusion.statement
    );
    // Class j = 1 is n = 0 + 1 + 2 m: the residue 1 of n modulo 2.
    assert_eq!(
        conclusion.zero_residue_classes,
        Some(ZeroResidueClasses {
            period: 2,
            residues: vec![1]
        })
    );
    assert!(conclusion.has_zero_residue_classes());
    assert_eq!(conclusion.first_index, 0);
    assert!(conclusion.proves_row_nonzero_and_real_rooted(4));
    assert!(!conclusion.proves_row_nonzero_and_real_rooted(3));
    // Relabelled to start at index 1 (P_1 = 1, P_2 = 0), the residues are
    // those of the index itself: the even rows vanish.
    let shifted = RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None)
        .reindexed(1)
        .unwrap();
    let (shifted_rule, shifted_first, shifted_seeds) = shifted.to_recurrence_parts().unwrap();
    let shifted = check_residue_class_subsequence_theorem(
        &shifted_rule,
        shifted_first,
        &shifted_seeds,
        &options(),
    )
    .unwrap();
    let shifted_conclusion = shifted.conclusion.as_ref().unwrap();
    assert_eq!(
        shifted_conclusion
            .zero_residue_classes
            .as_ref()
            .unwrap()
            .residues,
        vec![0]
    );
    assert!(shifted_conclusion.proves_row_nonzero_and_real_rooted(3));
    assert!(!shifted_conclusion.proves_row_nonzero_and_real_rooted(4));
    assert!(!shifted_conclusion.proves_row_nonzero_and_real_rooted(0));
    assert!(report
        .lean_declarations
        .contains(&"RealRooted.eq_zero_of_product_rec".to_string()));
    assert!(conclusion.interlacing.is_none());
    report.verify().unwrap();

    // The zero classes cannot be dropped or moved.
    assert!(tampered(&report, |v| {
        v["conclusion"]
            .as_object_mut()
            .unwrap()
            .remove("zero_residue_classes");
    }));
    assert!(tampered(&report, |v| {
        v["conclusion"]["zero_residue_classes"]["residues"] = serde_json::json!([0]);
    }));
    assert!(tampered(&report, |v| {
        v["conclusion"]["zero_residue_classes"] = serde_json::Value::Null;
    }));
}

/// The catalog families of the row routes whose residue-class proof has a
/// zero class, with the residues of `n` modulo 2 of their zero rows.
const ZERO_CLASS_IDS: [(&str, usize); 3] = [("A003982", 1), ("A132918", 1), ("A286100", 1)];

#[test]
fn zero_classes_are_exactly_the_zero_rows() {
    for id in NEWLY_PROVED {
        let (rule, seeds, first) = catalog(id);
        let report =
            check_residue_class_subsequence_theorem(&rule, first, &seeds, &options()).unwrap();
        let expected = ZERO_CLASS_IDS
            .iter()
            .find(|(zero_id, _)| *zero_id == id)
            .map(|&(_, residue)| ZeroResidueClasses {
                period: 2,
                residues: vec![residue],
            });
        let Some(conclusion) = report.conclusion.as_ref() else {
            assert!(expected.is_none(), "{id}");
            continue;
        };
        assert_eq!(conclusion.zero_residue_classes, expected, "{id}");
        // Without zero classes the JSON has no such key (earlier reports
        // keep their exact serialization).
        let json = serde_json::to_string(&report).unwrap();
        assert_eq!(
            json.contains("zero_residue_classes"),
            expected.is_some(),
            "{id}"
        );
        let generated = rule.generate_rows_rational(&seeds, first, 24).unwrap();
        for (offset, row) in generated.iter().enumerate() {
            let n = first + offset;
            let zero = polytool::vec_poly::is_zero(&polytool::vec_poly::trim(row));
            if conclusion.proves_row_nonzero_and_real_rooted(n) {
                assert!(!zero, "{id}: P_{n}");
                assert!(
                    polytool::is_real_rooted_bigint_coeffs(&integer(row)),
                    "{id}"
                );
            }
            if n >= conclusion.first_index
                && conclusion
                    .zero_residue_classes
                    .as_ref()
                    .is_some_and(|classes| classes.contains(n))
            {
                assert!(zero, "{id}: P_{n} lies in a zero class");
            }
        }
    }
}

#[test]
fn only_zero_classes_are_not_applicable() {
    // P_n = x^2 P_(n-2) with P_0 = P_1 = 0: every row is zero.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(2, 0, b(&[&[0, 0, 1]]))],
        denominator: None,
        inhomogeneous: None,
    };
    let report =
        check_residue_class_subsequence_theorem(&rule, 0, &rows(&[&[0], &[0]]), &options())
            .unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(report.conclusion.is_none() && report.hypotheses.is_empty());
    report.verify().unwrap();
}

#[test]
fn residue_classes_do_not_apply_or_fail_honestly() {
    // gcd 1: not applicable.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, b(&[&[0, 1]])),
            RecurrenceTerm::new(2, 0, b(&[&[1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let report =
        check_residue_class_subsequence_theorem(&rule, 0, &rows(&[&[1], &[0, 1]]), &options())
            .unwrap();
    assert_eq!(report.status, RealRootednessTheoremStatus::NotApplicable);
    assert!(report.residue_class_certificate().is_none());

    // P_n = x P_(n-2) + P_(n-4): the classes Q_m = x Q_(m-1) + Q_(m-2) with
    // Q_0 = 1, Q_1 = x give Q_2 = x^2 + 1, which is not real-rooted.
    let rule = Recurrence {
        terms: vec![
            RecurrenceTerm::new(2, 0, b(&[&[0, 1]])),
            RecurrenceTerm::new(4, 0, b(&[&[1]])),
        ],
        denominator: None,
        inhomogeneous: None,
    };
    let seeds = rows(&[&[1], &[1], &[0, 1], &[0, 1]]);
    let report = check_residue_class_subsequence_theorem(&rule, 0, &seeds, &options()).unwrap();
    assert!(!report.is_proved());
    assert!(report.conclusion.is_none());
    report.verify().unwrap();
    let certificate = report.residue_class_certificate().unwrap();
    assert!(certificate.classes.iter().all(|c| c.selected.is_none()));
    assert!(report
        .hypotheses
        .iter()
        .all(|h| h.status != LiftHypothesisStatus::Verified));
}

#[test]
fn tampered_residue_class_certificates_are_rejected() {
    let (rule, seeds, first) = catalog("A211956");
    let report = check_residue_class_subsequence_theorem(&rule, first, &seeds, &options()).unwrap();
    report.verify().unwrap();
    assert!(tampered(&report, |v| {
        v["certificate"]["residue_class_subsequences"]["classes"][0]["recurrence"]
            ["initial_polynomials"][0][0] = "2".into();
    }));
    assert!(tampered(&report, |v| {
        v["certificate"]["residue_class_subsequences"]["classes"][1]["selected"] = 7.into();
    }));
    assert!(tampered(&report, |v| {
        v["certificate"]["residue_class_subsequences"]["classes"][1]["zero"] = true.into();
    }));
    assert!(tampered(&report, |v| {
        v["certificate"]["residue_class_subsequences"]["class_maximum_start_delay"] = 17.into();
    }));
    assert!(tampered(&report, |v| {
        v["certificate"]["residue_class_subsequences"]["period"] = 4.into();
    }));
    assert!(tampered(&report, |v| {
        v["certificate"] = serde_json::Value::Null;
    }));
    assert!(tampered(&report, |v| {
        v["theorem"] = "variable_delay".into();
    }));
}

// ---------------------------------------------------------------------------
// Discovery and parent provers
// ---------------------------------------------------------------------------

#[test]
fn discovery_records_the_row_routes_only_where_they_apply() {
    let record = |id: &str| {
        let (rule, seeds, first) = catalog(id);
        RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None)
    };
    let options = RealRootednessTheoremOptions {
        maximum_start_delay: 0,
        ..Default::default()
    };
    for (id, theorem) in [
        ("A201013", RealRootednessTheorem::VariableDelay),
        ("A211956", RealRootednessTheorem::ResidueClassSubsequences),
    ] {
        let report = discover_recurrence_proofs(&record(id), &options).unwrap();
        assert_eq!(report.status, "proved", "{id}");
        // The eight discovery theorems, then the applicable row route.
        assert_eq!(report.attempts.len(), 8 + 1, "{id}");
        let attempt = report.attempts.last().unwrap();
        assert_eq!(attempt.theorem, theorem);
        assert!(attempt.is_proved());
        report.verify().unwrap();
    }
    // A residue-class proof with a zero class: no index has every later row
    // nonzero, so `real_rooted_from_index` is omitted and the zero classes
    // are reported.
    for (id, residue) in ZERO_CLASS_IDS {
        let record = record(id);
        let report = discover_recurrence_proofs(&record, &options).unwrap();
        assert_eq!(report.status, "proved", "{id}");
        assert_eq!(report.real_rooted_from_index, None, "{id}");
        let rows = report
            .real_rooted_outside_zero_residue_classes
            .as_ref()
            .unwrap();
        assert_eq!(rows.from_index, record.first_index, "{id}");
        assert_eq!(rows.zero_residue_classes.residues, vec![residue], "{id}");
        assert!(
            report.statement.contains("infinitely many zero rows"),
            "{id}"
        );
        assert!(report.tail.is_none() && report.alternative_representations.is_none());
        report.verify().unwrap();
        // A summary that hides the zero classes is rejected, and so is an
        // attempt whose conclusion drops them.
        let mut hidden = report.clone();
        hidden.real_rooted_from_index = Some(record.first_index);
        hidden.real_rooted_outside_zero_residue_classes = None;
        assert!(hidden.verify().is_err(), "{id}");
        let mut dropped = report.clone();
        let attempt = dropped.attempts.iter_mut().find(|a| a.is_proved()).unwrap();
        attempt.conclusion.as_mut().unwrap().zero_residue_classes = None;
        assert!(dropped.verify().is_err(), "{id}");
    }
    // Without zero classes the summary is unchanged.
    let report = discover_recurrence_proofs(&record("A099605"), &options).unwrap();
    assert_eq!(report.real_rooted_from_index, Some(0));
    assert!(report.real_rooted_outside_zero_residue_classes.is_none());
    // Eulerian rows (first order, lag 1): neither fixed-delay nor residue
    // splitting applies, but the lag-two Euler form includes c(n) = 0.
    let report = discover_recurrence_proofs(&record("A008292"), &options).unwrap();
    assert_eq!(report.attempts.len(), 9);
    let attempt = report.attempts.last().unwrap();
    assert_eq!(attempt.theorem, RealRootednessTheorem::LagTwoEuler);
    assert_eq!(
        attempt.status,
        RealRootednessTheoremStatus::HypothesisFailed
    );
    assert_eq!(
        hypothesis(attempt, "positive_linear_seed").status,
        LiftHypothesisStatus::Failed {
            witness: Some(BigInt::from(2))
        }
    );
    report.verify().unwrap();
    let report = discover_recurrence_proofs(
        &record("A008292"),
        &RealRootednessTheoremOptions {
            maximum_start_delay: 1,
            ..options
        },
    )
    .unwrap();
    let attempt = report.attempts.last().unwrap();
    assert_eq!(attempt.theorem, RealRootednessTheorem::LagTwoEuler);
    assert!(attempt.is_proved());
    assert_eq!(attempt.start_delay, 1);
    assert_eq!(
        attempt.lag_two_euler_certificate().unwrap().c.numerator,
        ["0"]
    );
    report.verify().unwrap();
}

fn parent(id: &str) -> VeroneseParentProof {
    let (rule, seeds, first) = catalog(id);
    let record = ParentRecurrenceRecord {
        rule: &rule,
        seeds: &seeds,
        first_index: first,
    };
    let proof = prove_parent_real_rootedness(&record, first, &[], false).unwrap();
    check_recorded_parent_proof(&record, first, &[], false, &proof).unwrap();
    proof
}

#[test]
fn parent_provers_use_the_row_routes_for_open_parents() {
    let proof = parent("A201013");
    assert_eq!(proof.status, VeroneseParentProofStatus::Proved);
    assert_eq!(proof.route, Some(VeroneseParentProofRoute::VariableDelay));
    let interlacing = proof.interlacing.as_ref().unwrap();
    assert_eq!(interlacing.kind, RowInterlacing::Weak);

    for id in ["A211956", "A171608", "A132918", "A003982", "A286100"] {
        let proof = parent(id);
        assert_eq!(proof.status, VeroneseParentProofStatus::Proved, "{id}");
        assert_eq!(
            proof.route,
            Some(VeroneseParentProofRoute::ResidueClassSubsequences),
            "{id}"
        );
        assert!(proof.interlacing.is_none());
        // The parent statement allows zero rows, so route 8d accepts zero
        // residue classes; its detail names them.
        assert!(proof.statement.contains("is zero or real-rooted"), "{id}");
        let attempt = proof.attempts.last().unwrap();
        assert!(attempt.success);
        let zero_class = ZERO_CLASS_IDS.iter().any(|(zero_id, _)| *zero_id == id);
        assert_eq!(
            attempt.detail.contains("residue classes j = [1]"),
            zero_class,
            "{id}: {}",
            attempt.detail
        );
    }
    // A parent proved by an earlier route records no row-route attempt.
    let proof = parent("A008292");
    assert!(proof.attempts.iter().all(|a| !matches!(
        a.route,
        VeroneseParentProofRoute::VariableDelay
            | VeroneseParentProofRoute::ResidueClassSubsequences
    )));
}

/// Every catalog entry that `recurrence-prove` proved only after the row
/// routes were added (catalog sweep, 60 s per entry).
const NEWLY_PROVED: [&str; 21] = [
    "A003982", "A099605", "A102547", "A109446", "A109447", "A118394", "A118931", "A118933",
    "A132918", "A137374", "A171608", "A201013", "A211956", "A223168", "A286100", "A305540",
    "A317496", "A317497", "A318772", "A318773", "A321434",
];

#[test]
fn newly_proved_catalog_rows_pass_the_exact_row_checks() {
    for id in NEWLY_PROVED {
        let (rule, seeds, first) = catalog(id);
        let delay = check_variable_delay_theorem(&rule, first, &seeds, &options()).unwrap();
        let classes =
            check_residue_class_subsequence_theorem(&rule, first, &seeds, &options()).unwrap();
        let report = if delay.is_proved() { delay } else { classes };
        assert!(report.is_proved(), "{id}");
        let conclusion = report.conclusion.as_ref().unwrap();
        let region = conclusion.root_region.as_ref().map(|r| r.at(0).unwrap());
        let generated = rule.generate_rows_rational(&seeds, first, 20).unwrap();
        for (offset, row) in generated.iter().enumerate() {
            let row = integer(row);
            if row.iter().all(|c| c == &BigInt::from(0)) {
                continue;
            }
            assert!(
                polytool::is_real_rooted_bigint_coeffs(&row),
                "{id}: P_{}",
                first + offset
            );
            if region
                .as_ref()
                .is_some_and(|r| r.upper.as_ref().is_some_and(|u| u <= &q(0)) && r.lower.is_none())
            {
                let zero = BigInt::from(0);
                assert!(
                    row.iter().all(|c| c >= &zero) || row.iter().all(|c| c <= &zero),
                    "{id}: P_{} has coefficients of both signs",
                    first + offset
                );
            }
        }
    }
}
