//! Proof discovery beyond the supplied rule: the full theorem list, tails
//! after rows that are not real-rooted, certified alternative
//! representations, proportional rows, inhomogeneous rules, replay of
//! first-version reports, and the Veronese-section tails and chaining.

use polytool::recurrence::*;

fn q(n: i64) -> BigRational {
    BigRational::from_integer(n.into())
}

/// The bundled catalog entry as the record `polytool oeis info --json`
/// prints (`recurrence_data`).
fn catalog_record(id: &str) -> RecurrenceJson {
    let entry = polytool::oeis::by_id(id).expect("bundled entry");
    let (rule, rows, first) = entry.recurrence_export_parts().unwrap();
    RecurrenceJson::from_recurrence_rational(&rule, first as usize, &rows, None)
}

fn discover(record: &RecurrenceJson) -> RecurrenceProofDiscoveryReport {
    let report = discover_recurrence_proofs(record, &RealRootednessTheoremOptions::default())
        .expect("discovery runs");
    report.verify().expect("the report replays");
    let replayed: RecurrenceProofDiscoveryReport =
        serde_json::from_str(&serde_json::to_string(&report).unwrap()).unwrap();
    replayed.verify().expect("the serialized report replays");
    report
}

fn theorems(attempts: &[RealRootednessTheoremReport]) -> Vec<RealRootednessTheorem> {
    attempts
        .iter()
        .filter(|a| a.is_proved())
        .map(|a| a.theorem)
        .collect()
}

#[test]
fn every_parent_cascade_theorem_is_attempted_in_order() {
    let record = catalog_record("A001263");
    let report = discover(&record);
    let order: Vec<_> = report.attempts.iter().map(|a| a.theorem).collect();
    assert_eq!(order, DISCOVERY_THEOREMS.to_vec());
    assert_eq!(
        &DISCOVERY_THEOREMS[6..],
        &[
            RealRootednessTheorem::IteratedPreserver,
            RealRootednessTheorem::HypergeometricRow
        ]
    );
    assert_eq!(report.status, "proved");
    assert_eq!(report.real_rooted_from_index, Some(record.first_index));
    assert!(report.tail.is_none() && report.alternative_representations.is_none());
}

#[test]
fn iterated_preserver_and_hypergeometric_rows_prove_through_discovery() {
    // A124860, A128540, A135065: hypergeometric rows; A129184, A172106 to
    // A172108, A216118, A219570: iterated preserver after zero seed rows.
    for id in ["A124860", "A128540", "A135065"] {
        let report = discover(&catalog_record(id));
        assert!(
            theorems(&report.attempts).contains(&RealRootednessTheorem::HypergeometricRow),
            "{id}"
        );
    }
    for (id, from) in [
        ("A073424", 1),
        ("A129184", 2),
        ("A172106", 2),
        ("A172107", 3),
        ("A172108", 4),
        ("A185778", 6),
        ("A216118", 4),
        ("A219570", 1),
    ] {
        let report = discover(&catalog_record(id));
        assert_eq!(report.status, "proved", "{id}");
        assert_eq!(report.real_rooted_from_index, Some(from), "{id}");
        assert!(
            theorems(&report.attempts).contains(&RealRootednessTheorem::IteratedPreserver),
            "{id}"
        );
        let proved = report.attempts.iter().find(|a| a.is_proved()).unwrap();
        assert!(proved.conclusion.as_ref().unwrap().first_index >= from);
    }
}

#[test]
fn a_tail_after_rows_that_are_not_real_rooted_is_proved_from_its_index() {
    // A122431: the prefix rows 3 + x and 3 + 3x + x^2 (not real-rooted)
    // precede x^k (1 + x)^3.
    let record = catalog_record("A122431");
    let report = discover(&record);
    assert_eq!(report.status, "proved");
    assert!(report.attempts.iter().all(|a| !a.is_proved()));
    let tail = report.tail.as_ref().expect("a tail");
    assert_eq!(tail.first_index, 3);
    assert_eq!(report.real_rooted_from_index, Some(3));
    assert!(tail.reason.contains("P_2 is not real-rooted"));
    assert!(theorems(&tail.attempts).contains(&RealRootednessTheorem::ProductRecurrence));
    assert_eq!(tail.recurrence.first_index, 3);

    // A synthetic family: P_0 = 1 + x^2, P_1 = 1, P_n = (1 + x) P_(n-1).
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q(1), q(1)]]),
        )],
        denominator: None,
        inhomogeneous: None,
    };
    let synthetic = RecurrenceJson::from_recurrence_rational(
        &rule,
        0,
        &[vec![q(1), q(0), q(1)], vec![q(1)]],
        None,
    );
    let report = discover(&synthetic);
    assert_eq!(report.real_rooted_from_index, Some(1));
    let tail = report.tail.clone().unwrap();
    assert_eq!(tail.first_index, 1);
    assert!(tail.reason.starts_with("P_0 is not real-rooted"));

    // Tampering with the tail is rejected.
    let mut moved = report.clone();
    moved.tail.as_mut().unwrap().first_index = 2;
    assert!(moved.verify().is_err());
    let mut seeds = report.clone();
    seeds.tail.as_mut().unwrap().recurrence.initial_polynomials[0] = vec!["2".into()];
    assert!(seeds.verify().is_err());
    let mut summary = report;
    summary.real_rooted_from_index = Some(0);
    assert!(summary.verify().is_err());
}

#[test]
fn alternative_representations_are_certified_before_they_are_used() {
    for (id, theorem) in [
        ("A030524", RealRootednessTheorem::MaWangInterval),
        ("A049324", RealRootednessTheorem::MaWangInterval),
        ("A049325", RealRootednessTheorem::MaWangInterval),
        ("A049326", RealRootednessTheorem::MaWangInterval),
        ("A383348", RealRootednessTheorem::MaWangInterval),
        ("A321966", RealRootednessTheorem::LiuWang),
        ("A271704", RealRootednessTheorem::LiuWangThreeTerm),
        ("A271705", RealRootednessTheorem::LiuWangThreeTerm),
        ("A180013", RealRootednessTheorem::ProductRecurrence),
    ] {
        let record = catalog_record(id);
        let report = discover(&record);
        assert_eq!(report.status, "proved", "{id}");
        assert!(report.attempts.iter().all(|a| !a.is_proved()), "{id}");
        assert_eq!(
            report.real_rooted_from_index,
            Some(record.first_index),
            "{id}"
        );
        let search = report.alternative_representations.as_ref().unwrap();
        let proved = search.alternatives.last().unwrap();
        assert!(proved.proved && proved.equality.is_certified(), "{id}");
        assert!(theorems(&proved.attempts).contains(&theorem), "{id}");
        proved.equality.verify().unwrap();
        assert!(proved.equality.statement.contains("all-index identity"));

        // A changed alternative rule no longer replays.
        let mut tampered = report.clone();
        let alternatives = &mut tampered
            .alternative_representations
            .as_mut()
            .unwrap()
            .alternatives;
        let last = alternatives.last_mut().unwrap();
        last.equality.alternative.initial_polynomials[0] = vec!["7".into()];
        assert!(tampered.verify().is_err(), "{id}");
    }
}

#[test]
fn the_equality_certificate_rejects_a_different_family() {
    let record = catalog_record("A030524");
    // n P_n = (6 + 3n + 2x) P_(n-1) + x (x + 9) P'_(n-1), P_1 = 1.
    let n_poly = |rows: Vec<Vec<i64>>| {
        BivarPoly::from_coefficients_normalized(
            rows.into_iter()
                .map(|row| row.into_iter().map(q).collect())
                .collect(),
        )
    };
    let alternative = Recurrence {
        terms: vec![
            RecurrenceTerm::new(1, 0, n_poly(vec![vec![6, 2], vec![3]])),
            RecurrenceTerm::new(1, 1, n_poly(vec![vec![0, 9, 1]])),
        ],
        denominator: Some(n_poly(vec![vec![0], vec![1]])),
        inhomogeneous: None,
    };
    let good = RecurrenceJson::from_recurrence_rational(&alternative, 1, &[vec![q(1)]], None);
    let certificate = certify_recurrence_equality(&record, &good).unwrap();
    assert!(certificate.is_certified(), "{:?}", certificate.reasons);
    assert_eq!(certificate.schema, RECURRENCE_EQUALITY_SCHEMA);
    certificate.verify().unwrap();
    let mut forged = certificate.clone();
    forged.checked_bases.push(99);
    assert!(forged.verify().is_err());

    // Another seed: the rows differ at once.
    let bad = RecurrenceJson::from_recurrence_rational(&alternative, 1, &[vec![q(2)]], None);
    let refused = certify_recurrence_equality(&record, &bad).unwrap();
    assert!(!refused.is_certified());
    assert!(refused.reasons[0].contains("row 1"));
    // Another coefficient: the rows agree on the seeds, but no identity.
    let mut other = alternative.clone();
    other.terms[1] = RecurrenceTerm::new(1, 1, n_poly(vec![vec![0, 8, 1]]));
    let other = RecurrenceJson::from_recurrence_rational(&other, 1, &[vec![q(1)]], None);
    assert!(!certify_recurrence_equality(&record, &other)
        .unwrap()
        .is_certified());
}

#[test]
fn proportional_rows_close_the_liu_wang_three_term_check() {
    // Rows c_n g(x)^k P_(N-2): B_n = b_n g^2 >= 0 defeats the general sign
    // condition, the scalar sequence c_n settles the theorem.
    for id in ["A124860", "A127647", "A128540", "A127733", "A135065"] {
        let record = catalog_record(id);
        let (rule, first, seeds) = record.to_recurrence_parts().unwrap();
        let report = check_liu_wang_three_term_theorem(
            &rule,
            first,
            &seeds,
            &RealRootednessTheoremOptions::default(),
        )
        .unwrap();
        assert!(report.is_proved(), "{id}: {:?}", report.open_hypotheses());
        assert!(report.normal_form.contains("proportional rows"), "{id}");
        report.verify().unwrap();
    }
    // A127733: P_n = (12 - 4n) x P_(n-1) + (4n - 3) x^2 P_(n-2), a_3 = 0, rows
    // n^2 x^(n-1): the hypergeometric closed form.
    let (rule, first, seeds) = catalog_record("A127733").to_recurrence_parts().unwrap();
    let report = check_liu_wang_three_term_theorem(
        &rule,
        first,
        &seeds,
        &RealRootednessTheoremOptions::default(),
    )
    .unwrap();
    let scalar = report
        .hypotheses
        .iter()
        .find(|h| h.name == "scalar_sequence_nonzero")
        .unwrap();
    assert!(scalar.evidence.contains("hypergeometric closed form"));
}

#[test]
fn inhomogeneous_rules_are_reported_not_applicable() {
    let report = discover(&catalog_record("A099759"));
    assert_eq!(report.status, "not_applicable");
    assert!(report
        .not_applicable_reason
        .as_deref()
        .unwrap()
        .contains("inhomogeneous"));
    assert!(report
        .attempts
        .iter()
        .all(|a| a.status == RealRootednessTheoremStatus::NotApplicable));
    assert!(report.alternative_representations.is_none());
}

#[test]
fn first_version_reports_keep_replaying() {
    let report = discover(&catalog_record("A008292"));
    assert_eq!(report.status, "proved");
    // A report of the first schema version: six attempts, no new fields.
    let mut legacy = serde_json::to_value(&report).unwrap();
    let object = legacy.as_object_mut().unwrap();
    object["attempts"].as_array_mut().unwrap().truncate(6);
    object.remove("real_rooted_from_index");
    object.insert(
        "statement".into(),
        serde_json::json!(
            "Each proved attempt establishes its conclusion for the recurrence-defined \
             sequence. Applying it to an independently defined sequence requires an all-index \
             recurrence identity and the stated seeds; finite fitting or held-out checks do not \
             establish that identity. No Lean replay was attempted."
        ),
    );
    let legacy: RecurrenceProofDiscoveryReport = serde_json::from_value(legacy).unwrap();
    legacy.verify().unwrap();
    // Dropping a first-version attempt is rejected.
    let mut dropped = legacy.clone();
    dropped.attempts.remove(2);
    assert!(dropped.verify().is_err());
    // A report missing only a newer theorem in the middle still replays.
    let mut without_preserver = report.clone();
    without_preserver.attempts.remove(6);
    without_preserver.verify().unwrap();
    // Reordering is rejected.
    let mut reordered = report.clone();
    reordered.attempts.swap(6, 7);
    assert!(reordered.verify().is_err());
    // An extra (repeated) attempt is rejected.
    let mut extra = report.clone();
    let first = extra.attempts[0].clone();
    extra.attempts.push(first);
    assert!(extra.verify().is_err());
    // An attempt flipped to proved is rejected.
    let mut flipped = report.clone();
    let open = flipped
        .attempts
        .iter_mut()
        .find(|a| a.status != RealRootednessTheoremStatus::Proved)
        .unwrap();
    open.status = RealRootednessTheoremStatus::Proved;
    assert!(flipped.verify().is_err());
    // A not_proved status with a proved attempt is rejected.
    let mut status = legacy;
    status.status = "not_proved".into();
    assert!(status.verify().is_err());
}

#[test]
fn veronese_lifts_use_tails_and_the_parent_provers() {
    let options = VeroneseSectionLiftOptions::default();
    // A171608, A286100: an open parent proved by the parent provers.
    for id in ["A171608", "A286100"] {
        let record = catalog_record(id);
        let (rule, first, rows) = record.to_recurrence_parts().unwrap();
        let report = find_veronese_section_lifts(&rule, &rows, first, &options).unwrap();
        let best = &report.lifts[0];
        assert!(
            matches!(
                best.status,
                VeroneseSectionLiftStatus::ProvedParentRoute { .. }
            ),
            "{id}: {:?}",
            best.status
        );
        assert!(best.status.is_proved());
        assert_eq!(best.status.label(), "proved_parent_route");
        assert!(best.lift_summary().is_proved());
    }
    // Rows that do not follow the recurrence: the search runs on the tail.
    for (id, from) in [
        ("A073424", 1),
        ("A122431", 3),
        ("A129184", 2),
        ("A185778", 6),
    ] {
        let record = catalog_record(id);
        let (rule, first, rows) = record.to_recurrence_parts().unwrap();
        let report = find_veronese_section_lifts(&rule, &rows, first, &options).unwrap();
        assert_eq!(report.tail_first_index, Some(from), "{id}");
    }
}
