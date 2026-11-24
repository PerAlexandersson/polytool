//! Tests for the enclosing module.

use super::*;
use crate::vec_poly::trim;

#[test]
fn catalog_is_sorted_unique_and_decodable() {
    assert_eq!(catalog().len(), 930);
    assert_eq!(
        catalog()
            .iter()
            .filter(|entry| entry.status == OeisSequenceStatus::Verified)
            .count(),
        179
    );
    assert_eq!(
        catalog()
            .iter()
            .filter(|entry| entry.status == OeisSequenceStatus::Validated)
            .count(),
        751
    );
    assert_eq!(
        catalog()
            .iter()
            .filter(|entry| entry.status == OeisSequenceStatus::Experimental)
            .count(),
        0
    );
    for pair in catalog().windows(2) {
        assert!(pair[0].id < pair[1].id);
    }
    for entry in catalog() {
        entry.recurrence_parts().expect(entry.id);
        if entry.status == OeisSequenceStatus::Verified {
            assert!(
                entry.verification_rows > 0,
                "{} lacks holdout rows",
                entry.id
            );
        }
    }
}

fn catalog_data_text() -> String {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("data/oeis-catalog.jsonl");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
}

#[test]
fn catalog_data_file_is_canonical_and_matches_the_static_catalog() {
    let text = catalog_data_text();
    let records = OeisCatalogRecord::parse_jsonl(&text).unwrap();
    for (line, record) in text.lines().zip(&records) {
        assert_eq!(
            line,
            record.to_jsonl_line(),
            "{} is not canonical",
            record.id
        );
    }
    assert_eq!(text.lines().count(), records.len());
    assert!(text.ends_with('\n'));
    assert_eq!(records.len(), catalog().len());
    for (record, entry) in records.iter().zip(catalog()) {
        assert_eq!(record.id, entry.id);
        assert_eq!(record.provenance.status, entry.status);
        assert_eq!(record.provenance.source_rows, entry.source_rows);
        assert_eq!(record.row_lengths, entry.row_lengths);
    }
}

#[test]
fn migrated_catalog_provenance_tiers_are_recorded() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    let count = |source| {
        records
            .iter()
            .filter(|record| record.provenance.source == source)
            .count()
    };
    assert_eq!(count(OeisProvenanceSource::BenchmarkFixture), 73);
    assert_eq!(count(OeisProvenanceSource::LeanDefinition), 630);
    assert_eq!(count(OeisProvenanceSource::CuratedQueue), 71);
    assert_eq!(count(OeisProvenanceSource::SuppliedRecurrence), 102);
    for record in &records {
        let provenance = &record.provenance;
        match provenance.source {
            OeisProvenanceSource::LeanDefinition | OeisProvenanceSource::SuppliedRecurrence => {
                assert!(provenance.validation_row.is_some(), "{}", record.id)
            }
            OeisProvenanceSource::BenchmarkFixture | OeisProvenanceSource::CuratedQueue => {
                assert!(provenance.validation_row.is_none(), "{}", record.id)
            }
            OeisProvenanceSource::PolytoolRecurrenceSearch => {}
        }
    }
}

#[test]
fn catalog_records_reject_inconsistent_provenance() {
    let text = catalog_data_text();
    let mut record: OeisCatalogRecord = serde_json::from_str(text.lines().next().unwrap()).unwrap();
    record.validate().unwrap();
    let mut verified_without_holdout = record.clone();
    verified_without_holdout.provenance.status = OeisSequenceStatus::Verified;
    assert_eq!(
        verified_without_holdout.validate().unwrap_err(),
        OeisCatalogError::VerifiedWithoutHeldOutRows {
            id: record.id.clone()
        }
    );
    let mut bad_id = record.clone();
    bad_id.id = "a001263".to_string();
    let error = bad_id.validate().unwrap_err();
    assert_eq!(
        error,
        OeisCatalogError::InvalidId {
            id: "a001263".to_string()
        }
    );
    assert_eq!(error.to_string(), "`a001263` is not an upper-case A-number");
    record.initial_rows[0][0] = "1/0".to_string();
    assert!(matches!(
        record.validate().unwrap_err(),
        OeisCatalogError::InvalidRowCoefficient { .. }
    ));
    let error = OeisCatalogRecord::parse_jsonl("{\"id\":\"A000001\"}").unwrap_err();
    assert!(matches!(
        &error,
        OeisCatalogError::InvalidJsonLine { line: 1, .. }
    ));
    assert!(error.to_string().starts_with("line 1: "));
    let records = OeisCatalogRecord::parse_jsonl(&text).unwrap();
    let swapped = vec![records[1].clone(), records[0].clone()];
    let error = OeisCatalogRecord::validate_all(&swapped).unwrap_err();
    assert!(error
        .to_string()
        .starts_with("records are not strictly sorted by A-number: "));
}

#[test]
fn explicit_table_row_lengths_round_trip_and_reject_invalid_shapes() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    let record = records
        .iter()
        .find(|record| record.id == "A266178")
        .unwrap();
    assert_eq!(
        record.row_lengths,
        Some(OeisCatalogRowLengths {
            first: 1,
            step: 2,
            every: 1
        })
    );
    assert_eq!(
        serde_json::from_str::<OeisCatalogRecord>(&record.to_jsonl_line()).unwrap(),
        *record
    );
    for (layout, lengths, width) in [
        (
            OeisLayout::RegularTriangle,
            OeisCatalogRowLengths {
                first: 1,
                step: 2,
                every: 1,
            },
            1,
        ),
        (
            OeisLayout::Table,
            OeisCatalogRowLengths {
                first: 0,
                step: 2,
                every: 1,
            },
            1,
        ),
        (
            OeisLayout::Table,
            OeisCatalogRowLengths {
                first: 1,
                step: 2,
                every: 0,
            },
            1,
        ),
        (
            OeisLayout::Table,
            OeisCatalogRowLengths {
                first: 1,
                step: 0,
                every: 1,
            },
            1,
        ),
        (
            OeisLayout::Table,
            OeisCatalogRowLengths {
                first: 1,
                step: usize::MAX,
                every: 1,
            },
            1,
        ),
        (
            OeisLayout::Table,
            OeisCatalogRowLengths {
                first: 1,
                step: 2,
                every: 1,
            },
            2,
        ),
    ] {
        let mut invalid = record.clone();
        invalid.layout = layout;
        invalid.row_lengths = Some(lengths);
        invalid.recurrence_first_width = width;
        assert_eq!(
            invalid.validate(),
            Err(OeisCatalogError::InvalidRowLengths {
                id: record.id.clone()
            })
        );
    }
    let mut unknown = serde_json::to_value(record).unwrap();
    unknown["row_lengths"]["phase"] = serde_json::json!(0);
    assert!(serde_json::from_value::<OeisCatalogRecord>(unknown).is_err());
    let lengths = OeisCatalogRowLengths {
        first: 1,
        step: usize::MAX,
        every: 1,
    };
    assert_eq!(lengths.checked_width(1), None);
    let lengths = OeisCatalogRowLengths {
        first: 1,
        step: usize::MAX / 2 + 1,
        every: 1,
    };
    assert_eq!(lengths.checked_width(2), None);
    // Absence remains absent in canonical serialization for legacy records.
    assert!(!records[0].to_jsonl_line().contains("row_lengths"));
}

#[test]
fn padded_cellular_automaton_rows_keep_every_oeis_coordinate() {
    for entry in [A266178(), A266298()] {
        let rows = entry.generate_rows(12).unwrap();
        for (n, coefficients) in rows.iter().enumerate() {
            assert_eq!(coefficients.len(), 2 * n + 1);
            assert_eq!(coefficients[0], BigInt::from(1));
            if n > 0 {
                let second = if entry.id == "A266178" { n % 2 } else { 1 };
                assert_eq!(coefficients[1], BigInt::from(second));
                assert!(coefficients[2..].iter().all(Zero::is_zero));
            }
        }
        assert_eq!(entry.generate_rows_from(7, 3).unwrap(), rows[7..10]);
        let bfile = flatten_rows_to_bfile(
            &rows,
            &BFileOptions {
                max_terms: Some(24),
                ..BFileOptions::default()
            },
        )
        .unwrap();
        assert_eq!(bfile.rows_included, 4);
        assert_eq!(bfile.terms_included, 16);
        assert_eq!(bfile.text.lines().last(), Some("15 0"));
        assert_eq!(
            entry.real_rootedness().unwrap().status,
            OeisRealRootednessStatus::LeanProved
        );
    }
}

#[test]
fn explicit_lengths_use_catalog_positions_and_pad_zero_rows() {
    let entry = OeisSequenceDefinition {
        first_row: 10,
        prefix_rows: &[&["7"]],
        recurrence_first_index: 7,
        recurrence_first_width: 1,
        initial_rows: &[&["0"]],
        row_lengths: Some(OeisCatalogRowLengths {
            first: 1,
            step: 2,
            every: 2,
        }),
        bfile_prefix_verified: false,
        ..*A266298()
    };
    let rows = entry.generate_rows(6).unwrap();
    assert_eq!(rows[0], row(&[7]));
    assert_eq!(
        rows.iter().map(Vec::len).collect::<Vec<_>>(),
        [1, 1, 3, 3, 5, 5]
    );
    assert!(rows[1..].iter().flatten().all(Zero::is_zero));
    assert_eq!(entry.generate_rows_from(14, 2).unwrap(), rows[4..6]);
}

#[test]
fn explicit_lengths_never_discard_nonzero_coefficients_and_errors_are_fused() {
    let entry = OeisSequenceDefinition {
        row_lengths: Some(OeisCatalogRowLengths {
            first: 1,
            step: 0,
            every: 1,
        }),
        ..*A266298()
    };
    let mut rows = entry.generate_rows_iter().unwrap();
    assert_eq!(rows.next().unwrap().unwrap(), row(&[1]));
    assert!(matches!(
        rows.next().unwrap(),
        Err(OeisSequenceError::InvalidDefinition { .. })
    ));
    assert!(rows.next().is_none());
    let invalid = OeisSequenceDefinition {
        row_lengths: Some(OeisCatalogRowLengths {
            first: 1,
            step: 2,
            every: 0,
        }),
        ..entry
    };
    assert!(matches!(
        invalid.generate_rows(1),
        Err(OeisSequenceError::InvalidDefinition { .. })
    ));
}

#[test]
fn real_rootedness_flags_parse_validate_and_round_trip() {
    let text = catalog_data_text();
    let line = text
        .lines()
        .find(|line| line.starts_with("{\"id\":\"A001263\""))
        .unwrap();
    let record: OeisCatalogRecord = serde_json::from_str(line).unwrap();
    let flag = record.real_rootedness.clone().unwrap();
    assert_eq!(flag.status, OeisRealRootednessStatus::LeanProved);
    assert_eq!(flag.status.as_str(), "lean_proved");
    assert_eq!(record.to_jsonl_line(), line);
    assert!(line
        .ends_with(",\"real_rootedness\":{\"status\":\"lean_proved\",\"as_of\":\"2026-10-03\"}}"));

    // Absent, without a date, and the second status all round-trip.
    let mut unknown = record.clone();
    unknown.real_rootedness = None;
    let unknown_line = unknown.to_jsonl_line();
    assert!(!unknown_line.contains("real_rootedness"));
    assert_eq!(
        serde_json::from_str::<OeisCatalogRecord>(&unknown_line).unwrap(),
        unknown
    );
    let mut refuted = record.clone();
    refuted.real_rootedness = Some(OeisCatalogRealRootedness {
        status: OeisRealRootednessStatus::KnownNotRealRooted,
        as_of: None,
        audit: None,
    });
    let refuted_line = refuted.to_jsonl_line();
    assert!(refuted_line.ends_with("\"real_rootedness\":{\"status\":\"known_not_real_rooted\"}}"));
    assert_eq!(
        serde_json::from_str::<OeisCatalogRecord>(&refuted_line).unwrap(),
        refuted
    );
    refuted.validate().unwrap();

    // Malformed dates, statuses and fields are rejected.
    for date in ["2026-10-3", "2026/10/03", "20261003", "2026-13-01"] {
        let mut bad = record.clone();
        bad.real_rootedness.as_mut().unwrap().as_of = Some(date.to_string());
        assert_eq!(
            bad.validate().unwrap_err(),
            OeisCatalogError::InvalidRealRootednessDate {
                id: record.id.clone()
            },
            "{date}"
        );
    }
    let wrong_status = line.replace("\"lean_proved\"", "\"proved\"");
    assert!(OeisCatalogRecord::parse_jsonl(&wrong_status).is_err());
    let extra_field = line.replace("\"as_of\"", "\"by\":\"x\",\"as_of\"");
    assert!(OeisCatalogRecord::parse_jsonl(&extra_field).is_err());
}

#[test]
fn supplied_recurrence_provenance_is_checked() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    let record = records
        .iter()
        .find(|record| record.provenance.source == OeisProvenanceSource::SuppliedRecurrence)
        .unwrap()
        .clone();
    record.validate().unwrap();
    assert_eq!(record.provenance.status, OeisSequenceStatus::Validated);
    assert_eq!(
        OeisProvenanceSource::SuppliedRecurrence.as_str(),
        "supplied_recurrence"
    );
    let mut without_terms = record.clone();
    without_terms.provenance.checked_terms = None;
    assert_eq!(
        without_terms.validate().unwrap_err(),
        OeisCatalogError::InvalidSuppliedProvenance {
            id: record.id.clone()
        }
    );
    let mut with_search = record.clone();
    with_search.provenance.fit_rows = Some(3);
    assert!(with_search.validate().is_err());
}

#[test]
fn bundled_real_rootedness_flags_are_exposed() {
    let flag = A001263().real_rootedness().unwrap();
    assert_eq!(flag.status, OeisRealRootednessStatus::LeanProved);
    assert_eq!(flag.as_of, Some("2026-10-03"));
    // Row 2 of A122431 is 3 + 3x + x^2.
    assert_eq!(
        A122431().real_rootedness().unwrap().status,
        OeisRealRootednessStatus::KnownNotRealRooted
    );
    assert_eq!(
        A122431().generate_rows(3).unwrap()[2],
        [3, 3, 1].map(BigInt::from)
    );
    // The Lean proof of A157153 is open, so it carries no flag.
    assert_eq!(A157153().real_rootedness(), None);
    let flagged = entries_with_real_rootedness_flag().collect::<Vec<_>>();
    let count = |status| {
        flagged
            .iter()
            .filter(|(_, flag)| flag.status == status)
            .count()
    };
    assert_eq!(count(OeisRealRootednessStatus::LeanProved), 814);
    assert_eq!(count(OeisRealRootednessStatus::KnownNotRealRooted), 1);
    assert_eq!(count(OeisRealRootednessStatus::TheoremProved), 5);
    for (entry, flag) in &flagged {
        assert_eq!(entry.real_rootedness(), Some(*flag));
    }
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    for (record, entry) in records.iter().zip(catalog()) {
        assert_eq!(
            record
                .real_rootedness
                .as_ref()
                .map(|flag| (flag.status, flag.as_of.as_deref())),
            entry
                .real_rootedness()
                .map(|flag| (flag.status, flag.as_of)),
            "{}",
            record.id
        );
    }
}

fn row(coefficients: &[i64]) -> Vec<BigInt> {
    coefficients.iter().copied().map(BigInt::from).collect()
}

#[test]
fn structural_theorem_attestations_round_trip_with_public_sources() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    for (id, source) in [
        ("A035469", "https://arxiv.org/abs/0803.1659v1"),
        ("A049029", "https://arxiv.org/abs/0803.1659v1"),
        ("A089503", "https://arxiv.org/abs/0803.1659v1"),
        ("A143494", "poissonPolynomial_isPFPolynomial"),
        ("A143495", "poissonPolynomial_isPFPolynomial"),
    ] {
        let record = records.iter().find(|r| r.id == id).unwrap();
        let flag = record.real_rootedness.as_ref().unwrap();
        assert_eq!(flag.status, OeisRealRootednessStatus::TheoremProved, "{id}");
        assert_eq!(flag.status.as_str(), "theorem_proved");
        assert_eq!(flag.as_of.as_deref(), Some("2026-11-24"));
        let audit = flag.audit.as_deref().unwrap();
        assert!(audit.contains(source), "{id}");
        assert!(!audit.contains("/workspace/"));
        assert!(!audit.contains("/cargo-target/"));
        assert_eq!(
            by_id(id).unwrap().real_rootedness().unwrap().audit,
            Some(audit)
        );
        assert_eq!(
            OeisCatalogRecord::parse_jsonl(&record.to_jsonl_line()).unwrap(),
            vec![record.clone()]
        );
    }
}

#[test]
fn structural_rook_attestations_match_catalog_rows_through_degree_40() {
    // Test-local board evaluator: Polytool remains standalone, so its tests
    // cannot depend on the workspace's foundational combinatorics crate.
    // In increasing column order, every previously occupied row is present
    // in the new column. Placing its rook leaves height-(k-1) choices.
    fn rook_numbers(columns: &[usize]) -> Vec<BigInt> {
        let mut numbers = vec![BigInt::from(1)];
        for &height in columns {
            let mut next = numbers.clone();
            next.push(BigInt::from(0));
            for (k, count) in numbers.iter().enumerate() {
                next[k + 1] += count * BigInt::from(height - k);
            }
            numbers = next;
        }
        numbers
    }

    for (id, step) in [("A035469", 4), ("A049029", 5), ("A089503", 0)] {
        let entry = by_id(id).unwrap();
        assert_eq!(entry.first_row, 1);
        let rows = entry.generate_rows(41).unwrap();
        for (m, actual) in rows.iter().enumerate() {
            let n = m + 1;
            let columns: Vec<_> = if step == 0 {
                vec![2 * n; n - 1]
            } else {
                (1..=n - 1).map(|j| step * j).collect()
            };
            let mut expected = rook_numbers(&columns);
            if step == 0 {
                // Independently check the rectangle formula as well.
                let mut falling = BigInt::from(1);
                for (k, coefficient) in expected.iter().enumerate() {
                    assert_eq!(
                        *coefficient,
                        crate::numeric::binomial_bigint(n - 1, k) * &falling
                    );
                    falling *= BigInt::from(2 * n - k);
                }
            } else {
                expected.reverse();
            }
            assert_eq!(*actual, expected, "{id}, displayed n={n}");
        }
    }
}

#[test]
fn structural_poisson_attestations_match_catalog_rows_through_degree_40() {
    for (id, offset) in [("A143494", 2), ("A143495", 3)] {
        let entry = by_id(id).unwrap();
        assert_eq!(entry.first_row, offset as i64);
        let rows = entry.generate_rows(41).unwrap();
        let mut expected = vec![BigInt::from(1)];
        for (m, actual) in rows.iter().enumerate() {
            assert_eq!(*actual, expected, "{id}, displayed n={}", m + offset);
            // Poisson step (t+a)Q+tQ': [t^k] becomes
            // (a+k)[t^k]Q + [t^(k-1)]Q, starting from the empty product 1.
            let mut next = vec![BigInt::from(0); expected.len() + 1];
            for (k, coefficient) in expected.iter().enumerate() {
                next[k] += coefficient * BigInt::from(offset + k);
                next[k + 1] += coefficient;
            }
            expected = next;
        }
    }
}

#[test]
fn evidence_row_properties_follow_the_documented_definitions() {
    use OeisEvidenceProperty::*;
    let holds = |coefficients: &[i64]| {
        OeisEvidenceProperty::ALL.map(|p| row_has_evidence_property(&row(coefficients), p))
    };
    // (1 + t)^3: real-rooted, so ultra-log-concave, log-concave and unimodal.
    assert_eq!(holds(&[1, 3, 3, 1]), [true; 4]);
    // The zero row and constants satisfy everything.
    assert_eq!(holds(&[]), [true; 4]);
    assert_eq!(holds(&[0, 0]), [true; 4]);
    assert_eq!(holds(&[5]), [true; 4]);
    // 3 + 3t + t^2: not real-rooted (discriminant 9 - 12 < 0) but log-concave.
    assert_eq!(holds(&[3, 3, 1]), [false, false, true, true]);
    // Mixed signs: t^2 - 1 is real-rooted; the coefficient properties fail.
    assert_eq!(holds(&[-1, 0, 1]), [true, false, false, false]);
    // An internal zero: 1 + t^3 has a non-real root and fails every property
    // although a_1^2 >= a_0 a_2 holds pointwise.
    assert!(crate::is_log_concave_bigint_coeffs(&row(&[1, 0, 0, 1])));
    assert_eq!(holds(&[1, 0, 0, 1]), [false; 4]);
    // Leading and trailing zeros are allowed: t^2 (1 + t)^2 and (1 + t) t^0.
    assert_eq!(holds(&[0, 0, 1, 2, 1, 0]), [true; 4]);
    // 1 + 3t + 3t^2 is not ultra-log-concave at its degree 2.
    assert_eq!(holds(&[1, 3, 3]), [false, false, true, true]);
    assert!(!row_has_evidence_property(&row(&[1]), NoProperty));
}

#[test]
fn evidence_accumulator_reports_strongest_and_first_failures() {
    let mut accumulator = OeisEvidenceAccumulator::new(7);
    assert!(accumulator.finish("2026-11-02").is_none());
    accumulator.push_row(&row(&[1]));
    accumulator.push_row(&row(&[1, 1]));
    let evidence = accumulator.finish("2026-11-02").unwrap();
    assert_eq!(evidence.strongest, OeisEvidenceProperty::RealRooted);
    assert_eq!(evidence.rows_checked, 2);
    assert_eq!(evidence.first_failure, OeisEvidenceFirstFailures::default());
    // Row 9 is 3 + 3t + t^2: real-rooted and ultra-log-concave fail there,
    // log-concavity and unimodality hold.
    accumulator.push_row(&row(&[3, 3, 1]));
    let evidence = accumulator.finish("2026-11-02").unwrap();
    assert_eq!(evidence.strongest, OeisEvidenceProperty::LogConcave);
    assert_eq!(evidence.rows_checked, 3);
    assert_eq!(evidence.first_failure.real_rooted, Some(9));
    assert_eq!(evidence.first_failure.ultra_log_concave, Some(9));
    assert_eq!(evidence.first_failure.log_concave, None);
    // Row 10 has mixed signs, so log-concavity and unimodality fail there
    // and no property holds on every row.
    accumulator.push_row(&row(&[-1, 0, 1]));
    let evidence = accumulator.finish("2026-11-02").unwrap();
    assert_eq!(evidence.strongest, OeisEvidenceProperty::NoProperty);
    assert_eq!(evidence.first_failure.real_rooted, Some(9));
    assert_eq!(evidence.first_failure.log_concave, Some(10));
    assert_eq!(evidence.first_failure.unimodal, Some(10));
    // Real-rooted rows with mixed signs keep `real_rooted` as strongest and
    // report no weaker failures.
    let mut signed = OeisEvidenceAccumulator::new(0);
    signed.push_row(&row(&[-1, 0, 1]));
    let evidence = signed.finish("2026-11-02").unwrap();
    assert_eq!(evidence.strongest, OeisEvidenceProperty::RealRooted);
    assert_eq!(evidence.first_failure, OeisEvidenceFirstFailures::default());
}

#[test]
fn evidence_records_parse_validate_and_round_trip() {
    let text = catalog_data_text();
    let line = text
        .lines()
        .find(|line| line.starts_with("{\"id\":\"A122431\""))
        .unwrap();
    let record: OeisCatalogRecord = serde_json::from_str(line).unwrap();
    assert!(record.evidence.is_none());
    assert_eq!(record.to_jsonl_line(), line);

    // Computed evidence is written after `real_rootedness`, at the end of the
    // line, and parses back to the same record.
    let computed = A122431()
        .compute_evidence(12, "2026-11-02")
        .unwrap()
        .unwrap();
    let mut with_computed = record.clone();
    with_computed.evidence = Some(computed.clone());
    with_computed.validate().unwrap();
    let computed_line = with_computed.to_jsonl_line();
    assert!(computed_line.starts_with(line.strip_suffix('}').unwrap()));
    assert!(computed_line.ends_with(&format!(
        ",\"evidence\":{}}}",
        serde_json::to_string(&computed).unwrap()
    )));
    assert!(
        computed_line.find("\"real_rootedness\"").unwrap()
            < computed_line.find("\"evidence\"").unwrap()
    );
    assert_eq!(
        serde_json::from_str::<OeisCatalogRecord>(&computed_line).unwrap(),
        with_computed
    );
    let without = record.clone();
    assert!(!without.to_jsonl_line().contains("evidence"));

    // A synthetic record covers every field and key order.
    let synthetic = OeisCatalogEvidence {
        strongest: OeisEvidenceProperty::NoProperty,
        rows_checked: 12,
        first_failure: OeisEvidenceFirstFailures {
            real_rooted: Some(3),
            ultra_log_concave: Some(4),
            log_concave: Some(4),
            unimodal: Some(5),
        },
        as_of: "2026-11-02".to_string(),
    };
    assert_eq!(
        serde_json::to_string(&synthetic).unwrap(),
        "{\"strongest\":\"none\",\"rows_checked\":12,\"first_failure\":\
         {\"real_rooted\":3,\"ultra_log_concave\":4,\"log_concave\":4,\"unimodal\":5},\
         \"as_of\":\"2026-11-02\"}"
    );
    let mut with = record.clone();
    with.evidence = Some(synthetic.clone());
    let with_line = with.to_jsonl_line();
    assert_eq!(
        serde_json::from_str::<OeisCatalogRecord>(&with_line).unwrap(),
        with
    );
    // Rows 3..=14 are covered when the first row is 3.
    with.first_row = 3;
    with.validate().unwrap();

    // Inconsistent or malformed records are rejected.
    let inconsistent = |edit: &dyn Fn(&mut OeisCatalogEvidence)| {
        let mut bad = with.clone();
        edit(bad.evidence.as_mut().unwrap());
        bad.validate().unwrap_err()
    };
    let error = OeisCatalogError::InconsistentEvidence {
        id: record.id.clone(),
    };
    assert_eq!(inconsistent(&|e| e.rows_checked = 0), error);
    assert_eq!(inconsistent(&|e| e.first_failure.unimodal = None), error);
    assert_eq!(
        inconsistent(&|e| e.first_failure.real_rooted = Some(15)),
        error
    );
    assert_eq!(
        inconsistent(&|e| e.first_failure.real_rooted = Some(2)),
        error
    );
    assert_eq!(
        inconsistent(&|e| {
            e.strongest = OeisEvidenceProperty::LogConcave;
        }),
        error
    );
    assert_eq!(
        inconsistent(&|e| e.as_of = "2026-1-02".to_string()),
        OeisCatalogError::InvalidEvidenceDate {
            id: record.id.clone()
        }
    );
    let unknown_property = computed_line.replace("\"log_concave\"", "\"strong\"");
    assert!(OeisCatalogRecord::parse_jsonl(&unknown_property).is_err());
    let extra_field = computed_line.replace("\"rows_checked\"", "\"by\":\"x\",\"rows_checked\"");
    assert!(OeisCatalogRecord::parse_jsonl(&extra_field).is_err());
}

#[test]
fn computed_evidence_of_known_entries_and_flag_consistency() {
    // No bundled entry carries a record unless the data file says so.
    for entry in catalog() {
        if let Some(evidence) = entry.evidence() {
            assert!(evidence.rows_checked > 0, "{}", entry.id);
        }
    }
    let a008292 = A008292()
        .compute_evidence(30, "2026-11-02")
        .unwrap()
        .unwrap();
    assert_eq!(a008292.strongest, OeisEvidenceProperty::RealRooted);
    assert_eq!(a008292.rows_checked, 30);
    assert_eq!(a008292.first_failure, OeisEvidenceFirstFailures::default());
    // Row 2 of A122431 is 3 + 3t + t^2: it fails real-rootedness and
    // ultra-log-concavity, but every checked row is log-concave.
    let a122431 = A122431()
        .compute_evidence(30, "2026-11-02")
        .unwrap()
        .unwrap();
    assert_eq!(a122431.strongest, OeisEvidenceProperty::LogConcave);
    assert_eq!(a122431.first_failure.real_rooted, Some(2));
    assert_eq!(a122431.first_failure.ultra_log_concave, Some(2));
    assert_eq!(a122431.first_failure.log_concave, None);
    assert_eq!(A122431().first_row, 0);
    assert!(A122431()
        .compute_evidence(0, "2026-11-02")
        .unwrap()
        .is_none());

    // A sample of Lean-proved entries is real-rooted on the checked rows.
    let sample = entries_with_real_rootedness_flag()
        .filter(|(_, flag)| flag.status == OeisRealRootednessStatus::LeanProved)
        .step_by(101);
    for (entry, _) in sample {
        let evidence = entry.compute_evidence(20, "2026-11-02").unwrap().unwrap();
        assert_eq!(
            evidence.strongest,
            OeisEvidenceProperty::RealRooted,
            "{}",
            entry.id
        );
    }
}

#[test]
fn direct_id_function_and_lookup_agree() {
    assert_eq!(A008292(), by_id("a008292").unwrap());
}

#[test]
fn a008292_begins_with_eulerian_triangle_rows() {
    let rows = A008292().generate_rows(4).unwrap();
    assert_eq!(
        rows,
        vec![
            vec![1.into()],
            vec![1.into(), 1.into()],
            vec![1.into(), 4.into(), 1.into()],
            vec![1.into(), 11.into(), 11.into(), 1.into()],
        ]
    );
}

#[test]
fn a046741_degree_eight_matches_its_canonical_recurrence() {
    let rows = A046741().generate_rows(9).unwrap();
    assert_eq!(
        rows[8],
        [1, 22, 191, 838, 1982, 2504, 1577, 424, 34]
            .into_iter()
            .map(BigInt::from)
            .collect::<Vec<_>>()
    );
}

#[test]
fn a166345_export_preserves_delayed_rows_and_displayed_indices() {
    let sequence = A166345();
    let (recurrence, initial_rows, first_index) = sequence.recurrence_export_parts().unwrap();
    assert_eq!(first_index, 1);
    assert_eq!(initial_rows.len(), 3);
    assert_eq!(
        recurrence.to_string(),
        "P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)"
    );

    let exported = recurrence
        .generate_rows_rational(&initial_rows, first_index as usize, 8)
        .unwrap();
    let expected = sequence
        .generate_rows(8)
        .unwrap()
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(BigRational::from_integer)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(exported, expected);
}

#[test]
fn oeis_export_includes_prefix_rows_before_recurrence_initial_rows() {
    let sequence = A166344();
    let (recurrence, initial_rows, first_index) = sequence.recurrence_export_parts().unwrap();
    assert_eq!(first_index, 1);
    assert_eq!(initial_rows.len(), 3);

    let exported = recurrence
        .generate_rows_rational(&initial_rows, first_index as usize, 8)
        .unwrap();
    let expected = sequence
        .generate_rows(8)
        .unwrap()
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(BigRational::from_integer)
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(exported, expected);
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "full test tier: cargo test -p polytool --features full-tests"
)]
fn every_oeis_export_reproduces_displayed_rows_past_its_initial_prefix() {
    for sequence in catalog() {
        let (recurrence, initial_rows, first_index) = sequence
            .recurrence_export_parts()
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        let row_count = initial_rows.len() + 2;
        let exported = recurrence
            .generate_rows_rational(&initial_rows, first_index as usize, row_count)
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        let expected = sequence
            .generate_rows(row_count)
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id))
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(BigRational::from_integer)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        assert_eq!(
            exported.iter().map(|row| trim(row)).collect::<Vec<_>>(),
            expected.iter().map(|row| trim(row)).collect::<Vec<_>>(),
            "{}",
            sequence.id
        );
    }
}

#[test]
fn bfile_offsets_account_for_complete_rows() {
    let sequence = A008292();
    assert!(sequence.bfile_prefix_verified);
    assert_eq!(sequence.flattened_offset_for_row(1).unwrap(), 1);
    assert_eq!(sequence.flattened_offset_for_row(4).unwrap(), 7);
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "full test tier: cargo test -p polytool --features full-tests"
)]
fn every_sparse_definition_reproduces_its_fixture_rows() {
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("fixtures/recurrence-benchmarks/rows");
    for sequence in catalog() {
        if sequence.fixture_slug.is_empty() {
            continue;
        }
        let path = fixtures.join(format!("{}.txt", sequence.fixture_slug));
        let expected = std::fs::read_to_string(&path)
            .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()))
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split(',')
                    .map(|value| value.trim().parse::<BigInt>().unwrap())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let generated = sequence
            .generate_rows(sequence.fixture_row_offset + expected.len())
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        let actual = generated[sequence.fixture_row_offset..]
            .iter()
            .map(|row| trim(row))
            .collect::<Vec<_>>();
        let expected = expected.iter().map(|row| trim(row)).collect::<Vec<_>>();
        assert_eq!(actual, expected, "{}", sequence.id);
    }
}

/// 64-bit FNV-1a; a stable, dependency-free fingerprint for regression
/// fixtures (not a cryptographic hash).
fn fnv1a64(text: &str) -> u64 {
    text.bytes().fold(0xcbf2_9ce4_8422_2325, |hash, byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

const FINGERPRINT_ROWS: usize = 20;

fn catalog_fingerprint(entry: &OeisSequenceDefinition) -> String {
    // Preserve the pre-extension fingerprints for absent optional row lengths.
    let definition = fnv1a64(&format!("{entry:?}").replace("row_lengths: None, ", ""));
    let rows = match entry.generate_rows(FINGERPRINT_ROWS) {
        Ok(rows) => rows
            .iter()
            .map(|row| {
                row.iter()
                    .map(ToString::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .collect::<Vec<_>>()
            .join(";"),
        Err(error) => format!("error: {error}"),
    };
    format!("{}\t{definition:016x}\t{:016x}", entry.id, fnv1a64(&rows))
}

/// Every catalog entry must keep its exact definition (all fields, as
/// printed by `Debug`) and its first rows.  The fixture was captured from
/// the catalog before its data moved to `data/oeis-catalog.jsonl`; entries
/// added later are not listed and are therefore not constrained.  Set
/// `POLYTOOL_BLESS_OEIS_FINGERPRINTS=1` to rewrite the fixture after an
/// intentional change to an existing entry.
#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "full test tier: cargo test -p polytool --features full-tests"
)]
fn catalog_entries_match_recorded_fingerprints() {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/oeis-catalog-fingerprints.tsv");
    if std::env::var_os("POLYTOOL_BLESS_OEIS_FINGERPRINTS").is_some() {
        let mut text = String::new();
        for entry in catalog() {
            text.push_str(&catalog_fingerprint(entry));
            text.push('\n');
        }
        std::fs::write(&path, text).unwrap();
        return;
    }
    let expected = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("failed to read {}: {error}", path.display()));
    let mut checked = 0;
    for line in expected.lines().filter(|line| !line.trim().is_empty()) {
        let id = line.split('\t').next().unwrap();
        let entry = by_id(id).unwrap_or_else(|| panic!("{id} disappeared from the catalog"));
        assert_eq!(catalog_fingerprint(entry), line, "{id}");
        checked += 1;
    }
    assert!(checked >= 817);
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "full test tier: cargo test -p polytool --features full-tests"
)]
fn every_imported_lean_definition_reproduces_its_validation_row() {
    for (id, row_index, expected) in OEIS_IMPORTED_VALIDATION_ROWS {
        let generated = by_id(id)
            .unwrap()
            .generate_rows(row_index + 1)
            .unwrap_or_else(|error| panic!("{id}: {error}"));
        let actual = &generated[*row_index];
        let expected = expected
            .iter()
            .map(|value| value.parse::<BigInt>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(actual, &expected, "{id} row {row_index}");
    }
}

#[test]
fn stretch_additions_preserve_source_rows_and_supplied_provenance() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    let cases: [(&str, usize, usize, &[&[i64]]); 9] = [
        (
            "A056940",
            1000,
            44,
            &[
                &[1],
                &[1, 1],
                &[1, 5, 1],
                &[1, 15, 15, 1],
                &[1, 35, 105, 35, 1],
            ],
        ),
        (
            "A140982",
            48,
            9,
            &[
                &[1],
                &[3, 1],
                &[6, 16, 1],
                &[10, 100, 50, 1],
                &[15, 400, 750, 120, 1],
            ],
        ),
        (
            "A142470",
            1000,
            44,
            &[
                &[1],
                &[1, 1],
                &[1, 8, 1],
                &[1, 30, 30, 1],
                &[1, 80, 300, 80, 1],
            ],
        ),
        (
            "A155497",
            1000,
            44,
            &[
                &[1],
                &[1, 1],
                &[1, 18, 1],
                &[1, 90, 90, 1],
                &[1, 280, 1400, 280, 1],
            ],
        ),
        (
            "A155516",
            1000,
            44,
            &[
                &[1],
                &[1, 1],
                &[1, 20, 1],
                &[1, 105, 105, 1],
                &[1, 336, 1764, 336, 1],
            ],
        ),
        (
            "A155537",
            1000,
            44,
            &[
                &[3],
                &[5, 5],
                &[9, 27, 9],
                &[17, 102, 102, 17],
                &[33, 330, 660, 330, 33],
            ],
        ),
        (
            "A202750",
            231,
            21,
            &[
                &[1],
                &[1, 1],
                &[1, 16, 1],
                &[1, 81, 81, 1],
                &[1, 256, 1296, 256, 1],
            ],
        ),
        (
            "A278880",
            1000,
            44,
            &[
                &[1],
                &[1, 1],
                &[1, 5, 1],
                &[1, 14, 14, 1],
                &[1, 30, 81, 30, 1],
            ],
        ),
        (
            "A370232",
            55,
            10,
            &[
                &[1],
                &[1, 1],
                &[1, 9, 1],
                &[1, 36, 25, 1],
                &[1, 100, 225, 49, 1],
            ],
        ),
    ];
    for (id, checked_terms, source_rows, expected) in cases {
        let entry = by_id(id).unwrap();
        let record = records.iter().find(|record| record.id == id).unwrap();
        assert_eq!(
            record.provenance.source,
            OeisProvenanceSource::SuppliedRecurrence,
            "{id}"
        );
        assert_eq!(
            record.provenance.status,
            OeisSequenceStatus::Validated,
            "{id}"
        );
        assert_eq!(record.provenance.checked_terms, Some(checked_terms), "{id}");
        assert_eq!(record.provenance.source_rows, source_rows, "{id}");
        assert_eq!(record.provenance.verification_rows, 0, "{id}");
        assert_eq!(
            record.provenance.terms_source.as_deref(),
            Some("bfile"),
            "{id}"
        );
        let flag = entry.real_rootedness().unwrap();
        assert_eq!(flag.status, OeisRealRootednessStatus::LeanProved, "{id}");
        assert_eq!(flag.as_of, Some("2026-11-22"), "{id}");
        let expected: Vec<Vec<BigInt>> = expected
            .iter()
            .map(|row| row.iter().copied().map(BigInt::from).collect())
            .collect();
        assert_eq!(
            entry.generate_rows(expected.len()).unwrap(),
            expected,
            "{id}"
        );
    }
    for id in [
        "A039811", "A039812", "A048743", "A055314", "A090657", "A090683", "A111636", "A129062",
        "A130191", "A141618", "A141686", "A155467", "A155493", "A175124", "A245111", "A320280",
        "A326600", "A383149", "A383206", "A384968",
    ] {
        assert!(by_id(id).is_none(), "{id} must remain outside the catalog");
    }
}

#[test]
#[cfg_attr(
    not(feature = "full-tests"),
    ignore = "full test tier: catalog provenance replay"
)]
fn recorded_term_digests_and_validation_rows_recompute_from_catalog_generators() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    let mut digests = 0;
    let mut validations = 0;
    for record in records {
        let entry = by_id(&record.id).unwrap();
        let count = record.provenance.source_rows.saturating_add(1).max(
            record
                .provenance
                .validation_row
                .as_ref()
                .map_or(0, |row| row.index + 1),
        );
        if record.provenance.checked_terms.is_none() && record.provenance.validation_row.is_none() {
            continue;
        }
        let rows = entry
            .generate_rows(count)
            .unwrap_or_else(|error| panic!("{}: {error}", record.id));
        if let Some(validation) = &record.provenance.validation_row {
            let expected: Vec<BigInt> = validation
                .terms
                .iter()
                .map(|c| c.parse().unwrap())
                .collect();
            assert_eq!(
                rows[validation.index], expected,
                "{} validation row",
                record.id
            );
            validations += 1;
        }
        if let Some(checked) = record.provenance.checked_terms {
            let terms: Vec<_> = rows
                .iter()
                .flatten()
                .take(checked)
                .map(ToString::to_string)
                .collect();
            assert_eq!(terms.len(), checked, "{} term count", record.id);
            assert_eq!(
                crate::digest::sha256_hex(terms.join(",").as_bytes()),
                record.provenance.rows_sha256,
                "{} rows_sha256",
                record.id
            );
            digests += 1;
        }
    }
    assert!(digests >= 11 && validations >= 11);
}

#[test]
fn eleven_additions_record_external_build_and_standard_axiom_audit() {
    let records = OeisCatalogRecord::parse_jsonl(&catalog_data_text()).unwrap();
    for id in [
        "A056940", "A140982", "A142470", "A155497", "A155516", "A155537", "A202750", "A278880",
        "A370232", "A266178", "A266298",
    ] {
        let flag = records
            .iter()
            .find(|r| r.id == id)
            .unwrap()
            .real_rootedness
            .as_ref()
            .unwrap();
        assert_eq!(flag.status, OeisRealRootednessStatus::LeanProved);
        let audit = flag.audit.as_deref().unwrap();
        assert!(
            audit.contains("external formalization corpus, audited 2026-11-23, standard axioms")
        );
        assert!(audit.contains("downstream build and axiom audit passed"));
        assert_eq!(
            by_id(id).unwrap().real_rootedness().unwrap().audit,
            Some(audit)
        );
    }
}
