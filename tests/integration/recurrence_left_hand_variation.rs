//! Left-hand side variation of the adaptive recurrence search
//! (`AdaptiveSearchOptions::try_laurent`): normalization charts of the
//! left-hand factor, `n`-dependent factors without `try_denominator`, and the
//! left-hand monomial `x^e` (`RecurrenceOptions::leader_x_power`), with the
//! historical results left unchanged.

use num_bigint::BigInt;
use polytool::recurrence::{
    find_polynomial_recurrence_rational, find_recurrence_adaptive_rational_with_budget,
    AdaptiveSearchBudget, AdaptiveSearchOptions, AdaptiveSearchOutcome, AdaptiveSearchResult,
    BigRational, BivarPoly, Recurrence, RecurrenceEvaluationError, RecurrenceJson,
    RecurrenceOptions, RecurrenceTerm,
};
use polytool::sequences;

type Q = BigRational;

use crate::common::q;

fn rational_rows(rows: &[Vec<BigInt>]) -> Vec<Vec<Q>> {
    rows.iter()
        .map(|row| row.iter().cloned().map(Q::from_integer).collect())
        .collect()
}

fn search(rows: &[Vec<Q>], options: &AdaptiveSearchOptions) -> Option<AdaptiveSearchResult> {
    match find_recurrence_adaptive_rational_with_budget(
        rows,
        options,
        AdaptiveSearchBudget::unbounded(),
    ) {
        AdaptiveSearchOutcome::Found(result) => Some(result),
        AdaptiveSearchOutcome::NoRecurrence(_) => None,
        AdaptiveSearchOutcome::BudgetExhausted(_) => unreachable!("unbounded"),
    }
}

/// A small search space: order at most 2, degrees at most 2 in `t` and 1
/// in `n`, no derivatives, and `try_denominator` left off.
fn small_space() -> AdaptiveSearchOptions {
    AdaptiveSearchOptions {
        max_rec_len: 2,
        max_var_deg: 2,
        max_idx_deg: 1,
        max_diff_deg: 0,
        ..AdaptiveSearchOptions::default()
    }
    .with_laurent(true)
}

/// The shifted Legendre polynomials from `n = 1`: their rule
/// `n P(n) = (2n - 1)(2t - 1) P(n-1) - (n - 1) P(n-2)` has the left-hand
/// factor `n`, whose constant term (in the 1-based search index) is zero.
/// The historical normalization cannot reach it; the `n` chart does, without
/// `try_denominator`.
#[test]
fn an_index_dependent_left_hand_factor_needs_a_chart() {
    let rows = rational_rows(&sequences::shifted_legendre_polynomials_bigint(14)[1..]);
    let options = small_space();
    assert!(!options.try_denominator);
    assert!(search(&rows, &options.clone().with_laurent(false)).is_none());
    let found = search(&rows, &options).expect("the n chart finds the rule");
    assert_eq!(
        found.recurrence.to_string(),
        "(n) P(n) = (1 - 2n - 2t + 4nt) P(n-1) + (1 - n) P(n-2)"
    );
    let denominator = found.recurrence.denominator.as_ref().unwrap();
    assert!(denominator.coeffs[0].iter().all(|c| *c == q(0)));
    // The rule regenerates the rows (exact division by n).
    let generated = found
        .recurrence
        .generate_rows_rational(&rows[..2], 1, rows.len())
        .unwrap();
    let trim = |row: &Vec<Q>| {
        let mut row = row.clone();
        while row.last() == Some(&q(0)) {
            row.pop();
        }
        row
    };
    assert_eq!(
        generated.iter().map(trim).collect::<Vec<_>>(),
        rows.iter().map(trim).collect::<Vec<_>>()
    );
}

/// The Narayana polynomials have the left-hand factor `n + 1`; it is found
/// without `try_denominator`, with the same result and candidate count as
/// the historical search.
#[test]
fn an_index_dependent_left_hand_factor_without_try_denominator() {
    let rows = rational_rows(&sequences::narayana_polynomials_bigint(14)[1..]);
    let options = small_space();
    let found = search(&rows, &options).expect("a rule");
    assert_eq!(
        found.recurrence.to_string(),
        "(1 + n) P(n) = (-1 + 2n - t + 2nt) P(n-1) + (2 - n - 4t + 2nt + 2t^2 - nt^2) \
         P(n-2)"
    );
    let historical = search(&rows, &options.with_laurent(false)).unwrap();
    assert_eq!(
        historical.recurrence.to_string(),
        found.recurrence.to_string()
    );
    assert_eq!(historical.candidates_tried, found.candidates_tried);
}

/// Results the historical search finds are unchanged, with the same number
/// of candidates tried.
#[test]
fn historical_results_are_unchanged() {
    let families = [
        rational_rows(&sequences::eulerian_polynomials_bigint(12)[1..]),
        rational_rows(&sequences::chebyshev_polynomials_t_bigint(12)),
        rational_rows(&sequences::hermite_polynomials_bigint(12)),
        rational_rows(&sequences::type_b_eulerian_polynomials_bigint(12)),
        rational_rows(&sequences::literature::delannoy_polynomials_bigint(12)),
    ];
    for rows in families {
        let options = AdaptiveSearchOptions::default().with_laurent(true);
        let with = search(&rows, &options).expect("a rule");
        let without = search(&rows, &options.with_laurent(false)).expect("a rule");
        assert_eq!(with.recurrence.to_string(), without.recurrence.to_string());
        assert_eq!(with.candidates_tried, without.candidates_tried);
        assert_eq!(with.opts.leader_x_power, 0);
    }
}

/// The fixed-shape search honours `leader_x_power`: the Eulerian rows
/// satisfy `x P(n) = (x + (n-1) x^2) P(n-1) + (x^2 - x^3) P'(n-1)`, and the
/// rule stores the left-hand factor `x` as its denominator.
#[test]
fn a_left_hand_monomial_in_the_fixed_shape_search() {
    let rows = rational_rows(&sequences::eulerian_polynomials_bigint(14)[1..]);
    let options = RecurrenceOptions {
        rec_len: 1,
        var_deg: 3,
        idx_deg: 1,
        diff_deg: 1,
        leader_x_power: 1,
        ..RecurrenceOptions::default()
    };
    let rule = find_polynomial_recurrence_rational(&rows, &options).expect("a rule");
    let denominator = rule.denominator.as_ref().expect("the factor x");
    assert_eq!(denominator.coeffs, vec![vec![q(0), q(1)]]);
    let generated = rule
        .generate_rows_rational(&rows[..1], 1, rows.len())
        .unwrap();
    assert_eq!(generated.len(), rows.len());
    for (a, b) in generated.iter().zip(&rows) {
        let mut a = a.clone();
        while a.last() == Some(&q(0)) {
            a.pop();
        }
        assert_eq!(&a, b);
    }
    // The JSON record of the rule carries the factor x in its denominator.
    let json = RecurrenceJson::from_recurrence_rational(&rule, 1, &rows[..1], None);
    let (back, _, _) = json.to_recurrence_parts().unwrap();
    assert_eq!(back.denominator, rule.denominator);
}

/// Row generation with a left-hand monomial divides exactly; a row that is
/// not divisible is a typed error, never truncated.
#[test]
fn a_non_divisible_row_is_a_typed_error() {
    // x P(n) = P(n-1) from P(0) = 1: 1/x is not a polynomial.
    let rule = Recurrence {
        terms: vec![RecurrenceTerm::new(
            1,
            0,
            BivarPoly::from_coefficients_normalized(vec![vec![q(1)]]),
        )],
        denominator: Some(BivarPoly::from_coefficients_normalized(vec![vec![
            q(0),
            q(1),
        ]])),
        inhomogeneous: None,
    };
    let error = rule
        .generate_rows_rational(&[vec![q(1)]], 0, 3)
        .unwrap_err();
    assert_eq!(error, RecurrenceEvaluationError::NonPolynomialQuotient);
}

/// Old serialized search options (without the new fields) still load, with
/// the library defaults (no left-hand variation, `x^e` up to `e = 2`).
#[test]
fn old_search_options_json_loads_with_defaults() {
    let mut value = serde_json::to_value(AdaptiveSearchOptions::default()).unwrap();
    let object = value.as_object_mut().unwrap();
    object.remove("try_laurent");
    object.remove("max_laurent_x_power");
    let options: AdaptiveSearchOptions = serde_json::from_value(value).unwrap();
    assert!(!options.try_laurent);
    assert_eq!(options.max_laurent_x_power, 2);
}
