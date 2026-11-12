//! Exact replay boundaries, singularities and primitive integral presentation.

use polytool::recurrence::{BigRational, RecurrenceJson};
use serde_json::json;

fn record(leader: &[&[&str]], coefficient: &[&[&str]]) -> RecurrenceJson {
    serde_json::from_value(json!({
        "schema":"polytool.recurrence.v1", "first_index":0, "initial_polynomials":[["1"]],
        "recurrence": {"denominator":{"coeffs":leader}, "terms":[{
            "offset":1,"deriv_order":0,"sign":"none","coeff":{"coeffs":coefficient}}]}
    }))
    .unwrap()
}

fn rows(values: &[i64]) -> Vec<Vec<BigRational>> {
    values
        .iter()
        .map(|&n| vec![BigRational::from_integer(n.into())])
        .collect()
}

#[test]
fn replay_minimizes_the_prefix_but_keeps_failures_and_missing_lags() {
    let record = record(&[&["1"]], &[&["2"]]);
    let ready = record
        .clone()
        .with_checked_presentation(&rows(&[1, 2, 4, 8]))
        .unwrap();
    assert_eq!(ready.initial_polynomials.len(), 1);
    assert_eq!(ready.presentation.unwrap().earliest_valid_index, 1);
    let failure = record
        .with_checked_presentation(&rows(&[7, 3, 6, 12]))
        .unwrap();
    assert_eq!(failure.initial_polynomials.len(), 2);
    assert_eq!(failure.presentation.unwrap().start_index, Some(2));
}

#[test]
fn all_integer_zeros_include_future_repeated_large_and_zero_roots() {
    // n(n-2)^2(n-100000000000000000000), with a harmless t factor.
    let record = record(
        &[
            &["0", "0"],
            &["-400000000000000000000", "-400000000000000000000"],
            &["400000000000000000004", "400000000000000000004"],
            &["-100000000000000000004", "-100000000000000000004"],
            &["1", "1"],
        ],
        &[&["1"]],
    );
    let rule = record.to_recurrence_parts().unwrap().0;
    assert_eq!(
        rule.left_hand_nonnegative_integer_zeros()
            .unwrap()
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>(),
        ["0", "2", "100000000000000000000"]
    );
    let checked = record.with_checked_presentation(&rows(&[1, 1, 1])).unwrap();
    let report = checked.presentation.as_ref().unwrap();
    assert_eq!(report.status, "insufficient_initial_rows");
    assert_eq!(report.start_index, None);
    assert!(checked
        .to_latex()
        .unwrap()
        .contains("need rows through index 100000000000000000000"));
    assert!(!checked.to_latex().unwrap().contains("\\geq"));
}

#[test]
fn t_dependent_leader_requires_common_zeros_and_zero_leader_is_explicit() {
    let nonzero = record(&[&["0", "1"], &["1"]], &[&["1"]]);
    assert!(nonzero
        .to_recurrence_parts()
        .unwrap()
        .0
        .left_hand_nonnegative_integer_zeros()
        .unwrap()
        .is_empty());
    let zero = record(&[&["0"]], &[&["0"]])
        .with_checked_presentation(&rows(&[1, 2]))
        .unwrap();
    assert_eq!(zero.presentation.unwrap().status, "zero_left_hand_factor");
}

#[test]
fn zeros_require_initial_rows_through_the_last_zero() {
    let checked = record(&[&["-1"], &["1"]], &[&["-1"], &["1"]])
        .with_checked_presentation(&rows(&[1, 5, 5, 5]))
        .unwrap();
    assert_eq!(checked.initial_polynomials.len(), 2);
    let report = checked.presentation.unwrap();
    assert_eq!(report.left_hand_zeros, ["1"]);
    assert_eq!(report.start_index, Some(2));
}

#[test]
fn normalization_keeps_fitted_scaling_and_normalizes_both_sides() {
    let record = record(
        &[&["1/2"], &["-1/2"]],
        &[&["0", "-3/2"], &["-3/2", "-1/2"], &["-1/2"]],
    );
    let decorated = record
        .clone()
        .with_checked_presentation(&rows(&[0, 1]))
        .unwrap();
    assert_eq!(record.recurrence, decorated.recurrence);
    let latex = decorated.to_latex().unwrap();
    assert!(
        latex.contains("(n - 1) P_{n} = (n + 3)(n + t) P_{n-1}"),
        "{latex}"
    );
    let normalized = decorated.normalized_recurrence.unwrap();
    assert_eq!(
        normalized.denominator.unwrap().coeffs,
        vec![vec!["-1"], vec!["1"]]
    );
}

#[test]
fn checked_indices_follow_relabeling_and_reject_overflow() {
    let record = record(&[&["-1"], &["1"]], &[&["-1"], &["1"]]);
    let shifted = record
        .reindexed(5)
        .unwrap()
        .with_checked_presentation(&rows(&[1, 2, 2]))
        .unwrap();
    assert_eq!(
        shifted.presentation.as_ref().unwrap().left_hand_zeros,
        ["6"]
    );
    assert_eq!(shifted.presentation.unwrap().start_index, Some(7));
    assert!(record
        .to_recurrence_parts()
        .unwrap()
        .0
        .presentation_for_rows(&rows(&[1]), usize::MAX)
        .is_err());
}

#[test]
fn scalar_left_hand_factors_are_printed_without_parentheses() {
    let source = record(&[&["2"]], &[&["1"]]);
    let rule = source.to_recurrence_parts().unwrap().0;
    assert_eq!(rule.to_string(), "2 P(n) = P(n-1)");
    assert!(rule.to_latex().starts_with("2\\, P_{n}"));
    assert!(source.to_latex().unwrap().contains("$2 P_{n} = P_{n-1}"));
    let padded = record(&[&["2", "0"], &["0"]], &[&["1"]]);
    assert_eq!(
        padded.to_recurrence_parts().unwrap().0.to_string(),
        "2 P(n) = P(n-1)"
    );
    for (leader, expected) in [
        (vec![vec!["0"], vec!["1"]], "(n) P(n) = P(n-1)"),
        (vec![vec!["0", "1"]], "(t) P(n) = P(n-1)"),
    ] {
        let slices: Vec<_> = leader.iter().map(Vec::as_slice).collect();
        let source = record(&slices, &[&["1"]]);
        assert_eq!(
            source.to_recurrence_parts().unwrap().0.to_string(),
            expected
        );
    }
}
