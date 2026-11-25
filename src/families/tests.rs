//! Reproducibility and malformed-sidecar regressions.
use super::*;
fn packing() -> FamilyRecord {
    catalog()
        .unwrap()
        .into_iter()
        .find(|r| r.id == "PT-F0001")
        .unwrap()
}
#[test]
fn bundled_records_validate_and_preserve_canonical_attestations() {
    let records = catalog().unwrap();
    assert_eq!(records.len(), oeis::catalog().len() + 2);
    validate_catalog(&records).unwrap();
    for e in oeis::catalog() {
        let r = records.iter().find(|r| r.id == e.id).unwrap();
        if let Some(flag) = e.real_rootedness() {
            assert_eq!(
                serde_json::to_value(r.real_rootedness.status).unwrap(),
                flag.status.as_str()
            );
            assert_eq!(r.real_rootedness.as_of.as_deref(), flag.as_of);
            assert_eq!(r.real_rootedness.statement, flag.audit.unwrap_or_default());
        }
    }
}
#[test]
fn packing_generators_reproduce_independent_published_rows() {
    let r = packing();
    assert_eq!(r.generate_rows(9).unwrap(), r.initial_rows);
    assert_eq!(r.generate_rows(5).unwrap()[4], vec!["1", "53", "3"]);
    assert_eq!(r.real_rootedness.status, FamilyProofStatus::Conjectured);
    let other = catalog()
        .unwrap()
        .into_iter()
        .find(|r| r.id == "PT-F0002")
        .unwrap();
    assert_eq!(other.generate_rows(9).unwrap(), other.initial_rows);
    assert_eq!(other.generate_rows(5).unwrap()[4], vec!["1", "59", "3"]);
    assert!(r.generate_rows(513).is_err());
}
#[test]
fn malformed_schema_references_dates_and_rows_are_rejected() {
    let records = catalog().unwrap();
    let original = records.iter().position(|r| r.id == "PT-F0001").unwrap();
    let reject = |edit: fn(&mut FamilyRecord)| {
        let mut changed = records.clone();
        edit(&mut changed[original]);
        assert!(validate_catalog(&changed).is_err());
    };
    reject(|r| r.related[0].id = "PT-F9999".into());
    reject(|r| r.initial_rows[4][1] = "54".into());
    reject(|r| r.oeis.checked = Some("2026-02-30".into()));
    reject(|r| r.oeis.checked = Some("éé-12-05".into()));
    reject(|r| r.real_rootedness.theorem = Some("missing".into()));
    reject(|r| r.submission.state = "published".into());
    reject(|r| r.submission.state = "submitted".into());
    reject(|r| r.schema_version = 2);
    reject(|r| {
        r.generator = FamilyGenerator::Oeis {
            id: "A000000".into(),
        }
    });
    let mut json = serde_json::to_value(packing()).unwrap();
    json["typo"] = true.into();
    assert!(parse_family(&json.to_string()).is_err());
    json.as_object_mut().unwrap().remove("typo");
    json["real_rootedness"]["status"] = "proved".into();
    assert!(parse_family(&json.to_string()).is_err());
}
#[test]
fn public_lean_links_have_optional_staleness_checks() {
    let records = catalog().unwrap();
    assert!(validate_lean_declarations(&records, &[]).is_err());
    let mut declarations = vec![];
    for r in &records {
        if let Some(l) = &r.real_rootedness.lean {
            declarations.push(FamilyLeanDeclaration {
                family_id: r.id.clone(),
                declaration: l.declaration.clone(),
                property: "real_rootedness".into(),
                repository: l.repository.clone(),
                commit: l.commit.clone(),
            });
        }
    }
    validate_lean_declarations(&records, &declarations).unwrap();
    declarations[0].commit = "0".repeat(40);
    assert!(validate_lean_declarations(&records, &declarations).is_err());
}
#[test]
fn exports_are_exact_prefixes_with_bounded_data_lines_and_bfiles() {
    let record = packing();
    let rows = record.generate_rows(40).unwrap();
    let export = export_oeis(&record, 40).unwrap();
    let lines: Vec<_> = export
        .blocks
        .lines()
        .filter(|l| ["%S ", "%T ", "%U "].iter().any(|tag| l.starts_with(tag)))
        .collect();
    assert!(!lines.is_empty() && lines.len() <= 3);
    assert!(lines.iter().all(|l| l.len() <= 80));
    let data = lines
        .iter()
        .map(|l| l.splitn(3, ' ').nth(2).unwrap())
        .collect::<String>();
    let terms: Vec<_> = data.split(',').collect();
    assert!(rows
        .iter()
        .flatten()
        .map(String::as_str)
        .eq(export.bfile.lines().map(|l| l.split_once(' ').unwrap().1)));
    assert!(rows
        .iter()
        .flatten()
        .take(terms.len())
        .map(String::as_str)
        .eq(terms));
    assert!(export.proposal_toml.contains("[oeis_search]"));
    assert!(export.blocks.contains("conjectured"));
    assert!(export_oeis(&record, 0).is_err());
    let bounded = export_oeis(&record, 210).unwrap();
    assert!(bounded.bfile_terms <= 10_000 && bounded.bfile.len() <= 1024 * 1024);
    assert!(bounded.bfile_truncated);
}

#[test]
fn recurrence_generator_reuses_canonical_evaluator_and_rejects_fractions() {
    let mut record = packing();
    let source = oeis::by_id("A008292").unwrap();
    let (rule, seeds) = source.recurrence_parts().unwrap();
    let mut recurrence = crate::recurrence::RecurrenceJson::from_recurrence_rational(
        &rule,
        source.recurrence_first_index,
        &seeds,
        None,
    );
    record.generator = FamilyGenerator::Recurrence {
        record: Box::new(recurrence.clone()),
    };
    let expected: Vec<Vec<String>> = source
        .generate_rows(6)
        .unwrap()
        .into_iter()
        .map(|row| row.iter().map(ToString::to_string).collect())
        .collect();
    assert_eq!(record.generate_rows(6).unwrap(), expected);
    recurrence.initial_polynomials[0][0] = "1/2".into();
    record.generator = FamilyGenerator::Recurrence {
        record: Box::new(recurrence),
    };
    assert!(record.generate_rows(1).is_err());
    record = packing();
    record.first_row = i64::MAX;
    record.generator = FamilyGenerator::Oeis {
        id: "A008292".into(),
    };
    assert!(export_oeis(&record, 2).is_err());
}
