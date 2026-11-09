//! Tests of the public building blocks shared by the lifts:
//! `recurrence::prepare_lift_rule`, `recurrence::prove_parent_real_rootedness`,
//! `recurrence::find_bounded_root_region_report` and
//! `PreserverParentLibrary::add_parent_record`.

use polytool::recurrence::{
    find_bounded_root_region_report, prepare_lift_rule, prove_parent_real_rootedness, BivarPoly,
    ParentRecurrenceRecord, PreserverParentLibrary, PreserverParentLibraryOptions,
    PreserverParentRegion, Recurrence, RecurrenceError, RecurrenceJson, RecurrenceTerm,
    VeroneseParentProofStatus,
};

mod common;
use common::rows;

/// `coeffs[i][j]` = coefficient of `n^i x^j`.
fn bivar(coeffs: &[&[i64]]) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(rows(coeffs))
}

/// `f(n) P_n = sum c_k(n, x) D^d P_(n-k)` from `(k, d, c_k)`.
fn rule(terms: &[(usize, usize, &[&[i64]])], denominator: Option<&[&[i64]]>) -> Recurrence {
    Recurrence {
        terms: terms
            .iter()
            .map(|(offset, order, coeffs)| RecurrenceTerm::new(*offset, *order, bivar(coeffs)))
            .collect(),
        denominator: denominator.map(bivar),
        inhomogeneous: None,
    }
}

/// Eulerian: `A_n = (1 + (n - 1) x) A_(n-1) + x (1 - x) A'_(n-1)`.
fn eulerian() -> Recurrence {
    rule(
        &[(1, 0, &[&[1, -1], &[0, 1]]), (1, 1, &[&[0, 1, -1]])],
        None,
    )
}

/// Chebyshev `U`: `U_n = 2 x U_(n-1) - U_(n-2)`, `U_0 = 1`, `U_1 = 2x`.
fn chebyshev_second_kind() -> Recurrence {
    rule(&[(1, 0, &[&[0, 2]]), (2, 0, &[&[-1]])], None)
}

#[test]
fn prepared_eulerian_rule_has_one_identifying_seed_and_replays() -> Result<(), RecurrenceError> {
    let supplied = rows(&[&[1], &[1, 1], &[1, 4, 1], &[1, 11, 11, 1]]);
    let prepared = prepare_lift_rule(&eulerian(), &supplied, 1, "parent")?;
    assert_eq!(prepared.first_index, 1);
    assert_eq!(prepared.leader_shift, 1);
    assert_eq!(prepared.seeds, rows(&[&[1]]));
    assert_eq!(prepared.rows(4)?, supplied);
    Ok(())
}

#[test]
fn prepare_rejects_rows_the_rule_does_not_generate() {
    let wrong = rows(&[&[1], &[1, 1], &[1, 5, 1]]);
    let error = prepare_lift_rule(&eulerian(), &wrong, 1, "parent").unwrap_err();
    assert!(matches!(error, RecurrenceError::InvalidInput(_)), "{error}");
    assert!(error.to_string().contains("parent row 3"), "{error}");
}

#[test]
fn eulerian_parent_is_proved_for_every_index() -> Result<(), RecurrenceError> {
    let prepared = prepare_lift_rule(&eulerian(), &rows(&[&[1]]), 1, "parent")?;
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &prepared.rule,
            seeds: &prepared.seeds,
            first_index: prepared.first_index,
        },
        prepared.first_index,
        &[],
        false,
    )?;
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Proved,
        "{proof:#?}"
    );
    assert!(proof.route.is_some());
    assert!(proof.conditional_on.is_empty());
    Ok(())
}

#[test]
fn a_parent_with_a_positive_root_is_refuted() -> Result<(), RecurrenceError> {
    // P_n = (x - 1) P_(n-1), P_0 = 1: P_1 = x - 1 has the positive root 1.
    let shifted = rule(&[(1, 0, &[&[-1, 1]])], None);
    let proof = prove_parent_real_rootedness(
        &ParentRecurrenceRecord {
            rule: &shifted,
            seeds: &rows(&[&[1]]),
            first_index: 0,
        },
        0,
        &[],
        false,
    )?;
    assert_eq!(
        proof.status,
        VeroneseParentProofStatus::Refuted,
        "{proof:#?}"
    );
    assert_eq!(proof.failing_row_index, Some(1));
    Ok(())
}

#[test]
fn chebyshev_rows_have_a_bounded_replayable_root_region() -> Result<(), RecurrenceError> {
    let report =
        find_bounded_root_region_report(&chebyshev_second_kind(), 0, &rows(&[&[1], &[0, 2]]))?
            .expect("a bounded region");
    assert!(report.is_proved());
    let region = report
        .conclusion
        .as_ref()
        .and_then(|conclusion| conclusion.root_region.as_ref())
        .expect("a root region");
    assert!(region.lower.is_some() && region.upper.is_some());
    report.verify()?;
    Ok(())
}

#[test]
fn eulerian_rows_have_no_bounded_root_region() -> Result<(), RecurrenceError> {
    // The roots of the Eulerian polynomials tend to -infinity.
    assert!(find_bounded_root_region_report(&eulerian(), 1, &rows(&[&[1]]))?.is_none());
    Ok(())
}

#[test]
fn add_parent_record_proves_keeps_and_deduplicates() -> Result<(), RecurrenceError> {
    let options = PreserverParentLibraryOptions {
        include_named_families: false,
        ..PreserverParentLibraryOptions::default()
    };
    let record = RecurrenceJson::from_recurrence_rational(&eulerian(), 1, &rows(&[&[1]]), None);
    let mut library = PreserverParentLibrary::default();
    library.add_parent_record("eulerian".into(), &record, &options)?;
    assert_eq!(library.parents.len(), 1);
    assert_eq!(library.parents[0].name, "eulerian");
    assert_eq!(
        library.parents[0].region,
        PreserverParentRegion::Nonpositive
    );
    library.add_parent_record("again".into(), &record, &options)?;
    assert_eq!(library.parents.len(), 1);
    assert_eq!(library.skipped.len(), 1);
    assert_eq!(library.skipped[0].0, "again");
    Ok(())
}
