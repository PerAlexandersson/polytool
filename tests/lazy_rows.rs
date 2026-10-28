//! Lazy recurrence row iterators: the shared kernel, the OEIS adapter and the
//! named polynomial families.
//!
//! The batch functions delegate to the iterators, so agreement between them
//! proves nothing by itself.  These tests instead check laziness through the
//! evaluation counters and through recurrences that fail at a known later
//! step, compare against a naive full-history reference loop written here, and
//! pin a few independently known values (OEIS A-numbers cited inline).

use num_bigint::BigInt;
use polytool::oeis::{
    by_id, catalog, OeisLayout, OeisSequenceDefinition, OeisSequenceError, SparseCoefficient,
    SparseRecurrence, SparseRecurrenceTerm,
};
use polytool::recurrence::{
    BigRational, BivarPoly, Recurrence, RecurrenceEvaluationError, RecurrenceSign, RecurrenceTerm,
};
use polytool::sequences::{self, literature};

fn q(value: i64) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn row(values: &[i64]) -> Vec<BigRational> {
    values.iter().map(|&value| q(value)).collect()
}

fn ints(values: &[i64]) -> Vec<BigInt> {
    values.iter().map(|&value| BigInt::from(value)).collect()
}

/// `sum coeff * n^i * t^j` from `(i, j, coeff)` triples.
fn bivar(monomials: &[(usize, usize, i64)]) -> BivarPoly {
    let n_size = monomials.iter().map(|m| m.0 + 1).max().unwrap_or(1);
    let t_size = monomials.iter().map(|m| m.1 + 1).max().unwrap_or(1);
    let mut coeffs = vec![vec![q(0); t_size]; n_size];
    for &(i, j, c) in monomials {
        coeffs[i][j] = q(c);
    }
    BivarPoly::from_coefficients_normalized(coeffs)
}

fn recurrence(
    terms: Vec<RecurrenceTerm>,
    denominator: Option<BivarPoly>,
    inhomogeneous: Option<BivarPoly>,
) -> Recurrence {
    Recurrence {
        terms,
        denominator,
        inhomogeneous,
    }
}

/// `P_n = n P_{n-1}`, so `P_n = n!` when seeded with `P_0 = 1`.
fn factorial_recurrence() -> Recurrence {
    recurrence(
        vec![RecurrenceTerm::new(1, 0, bivar(&[(1, 0, 1)]))],
        None,
        None,
    )
}

/// `P_n = t P_{n-1} + P_{n-2}`: the Fibonacci polynomials (OEIS A011973).
fn fibonacci_polynomial_recurrence() -> Recurrence {
    recurrence(
        vec![
            RecurrenceTerm::new(1, 0, bivar(&[(0, 1, 1)])),
            RecurrenceTerm::new(2, 0, bivar(&[(0, 0, 1)])),
        ],
        None,
        None,
    )
}

/// `(n - 5) P_n = P_{n-1}`: the step at `n = 5` has a zero denominator.
fn zero_denominator_at_five() -> Recurrence {
    recurrence(
        vec![RecurrenceTerm::new(1, 0, bivar(&[(0, 0, 1)]))],
        Some(bivar(&[(1, 0, 1), (0, 0, -5)])),
        None,
    )
}

/// The pre-iterator batch algorithm: full history, one evaluation per row.
fn naive_rows(
    recurrence: &Recurrence,
    seeds: &[Vec<BigRational>],
    first_index: usize,
    total: usize,
) -> Result<Vec<Vec<BigRational>>, RecurrenceEvaluationError> {
    if total <= seeds.len() {
        return Ok(seeds[..total].to_vec());
    }
    let required = recurrence.max_offset();
    if seeds.len() < required {
        return Err(RecurrenceEvaluationError::InsufficientInitialRows {
            required,
            provided: seeds.len(),
        });
    }
    let mut rows = seeds.to_vec();
    while rows.len() < total {
        let next = recurrence.evaluate_next_rational(&rows, first_index)?;
        rows.push(next);
    }
    Ok(rows)
}

// ---------------------------------------------------------------------------
// Shared kernel
// ---------------------------------------------------------------------------

#[test]
fn kernel_prefix_pause_and_continue_without_recomputation() {
    // Fibonacci polynomials F_1 = 1, F_2 = t, F_3 = t^2 + 1, F_4 = t^3 + 2t,
    // F_5 = t^4 + 3t^2 + 1, F_6 = t^5 + 4t^3 + 3t (A011973).
    let seeds = vec![row(&[1]), row(&[0, 1])];
    let mut rows = fibonacci_polynomial_recurrence().generate_rows_rational_iter(&seeds, 1);
    assert_eq!(rows.first_index(), 1);
    assert_eq!(rows.steps_evaluated(), 0, "construction evaluates nothing");

    let prefix: Vec<_> = rows.by_ref().take(4).collect::<Result<_, _>>().unwrap();
    assert_eq!(
        prefix,
        vec![row(&[1]), row(&[0, 1]), row(&[1, 0, 1]), row(&[0, 2, 0, 1])]
    );
    assert_eq!(rows.rows_yielded(), 4);
    assert_eq!(rows.steps_evaluated(), 2, "only the two computed rows");
    assert_eq!(rows.next_index(), Some(5));

    // Pause, then continue: no row repeated, no step recomputed.
    assert_eq!(rows.next().unwrap().unwrap(), row(&[1, 0, 3, 0, 1]));
    assert_eq!(rows.next().unwrap().unwrap(), row(&[0, 3, 0, 4, 0, 1]));
    assert_eq!(rows.steps_evaluated(), 4);
    assert_eq!(rows.next_index(), Some(7));
}

#[test]
fn kernel_uses_absolute_indices_for_nonzero_origins() {
    // P_n = n P_{n-1} seeded with P_3 = 6 gives n!, so the index n inside the
    // coefficients is the absolute one.
    let rows: Vec<_> = factorial_recurrence()
        .generate_rows_rational_iter(&[row(&[6])], 3)
        .take(4)
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(rows, vec![row(&[6]), row(&[24]), row(&[120]), row(&[720])]);
}

#[test]
fn kernel_is_lazy_and_an_error_at_a_later_step_keeps_the_earlier_prefix() {
    // P_0 = 1, then (n - 5) P_n = P_{n-1}: rows 0..=4 are fine, row 5 fails.
    let mut rows = zero_denominator_at_five().generate_rows_rational_iter(&[row(&[1])], 0);
    // Stopping after three rows evaluated exactly two steps, not the failing one.
    let early: Vec<_> = rows.by_ref().take(3).collect::<Result<_, _>>().unwrap();
    assert_eq!(early[0], row(&[1]));
    assert_eq!(
        early[1],
        vec![BigRational::new(BigInt::from(-1), BigInt::from(4))]
    );
    assert_eq!(rows.steps_evaluated(), 2);
    assert_eq!(rows.next_index(), Some(3));

    // The rest of the valid prefix is still available, then the error, once.
    assert!(rows.next().unwrap().is_ok()); // n = 3
    assert!(rows.next().unwrap().is_ok()); // n = 4
    assert_eq!(rows.steps_evaluated(), 4);
    assert_eq!(
        rows.next(),
        Some(Err(RecurrenceEvaluationError::ZeroDenominator))
    );
    assert_eq!(rows.next(), None);
    assert_eq!(rows.next(), None, "fused after the terminal error");
    assert_eq!(rows.retained_row_count(), 0);
    assert_eq!(rows.size_hint(), (0, Some(0)));

    // The batch form reports the same error; a shorter request does not hit it.
    let recurrence = zero_denominator_at_five();
    assert_eq!(
        recurrence.generate_rows_rational(&[row(&[1])], 0, 6),
        Err(RecurrenceEvaluationError::ZeroDenominator)
    );
    assert_eq!(
        recurrence
            .generate_rows_rational(&[row(&[1])], 0, 5)
            .unwrap()
            .len(),
        5
    );
}

#[test]
fn kernel_retains_only_a_bounded_history_window() {
    let recurrence = fibonacci_polynomial_recurrence();
    let mut rows = recurrence.generate_rows_rational_iter(&[row(&[1]), row(&[0, 1])], 1);
    assert_eq!(rows.history_limit(), 2);
    for _ in 0..200 {
        rows.next().unwrap().unwrap();
        assert!(rows.retained_row_count() <= 2);
    }
    assert_eq!(rows.retained_row_count(), 2);
    assert_eq!(rows.steps_evaluated(), 198);

    // Surplus seeds are held only until yielded: five seeds, order one.
    let seeds: Vec<_> = (1..=5).map(|k| row(&[k])).collect();
    let mut rows = factorial_recurrence().generate_rows_rational_iter(&seeds, 1);
    assert_eq!(rows.history_limit(), 1);
    let mut retained = Vec::new();
    for _ in 0..7 {
        rows.next().unwrap().unwrap();
        retained.push(rows.retained_row_count());
    }
    assert_eq!(retained, vec![5, 4, 3, 2, 1, 1, 1]);

    // An order-zero recurrence-free window keeps nothing.
    let constant = recurrence_with_inhomogeneous_term();
    let mut rows = constant.generate_rows_rational_iter(&[], 0);
    assert_eq!(rows.history_limit(), 0);
    for expected in 0..5 {
        assert_eq!(rows.next().unwrap().unwrap(), row(&[expected]));
        assert_eq!(rows.retained_row_count(), 0);
    }
}

/// No terms and `P_n = n`: an order-zero (history-free) recurrence.
fn recurrence_with_inhomogeneous_term() -> Recurrence {
    recurrence(Vec::new(), None, Some(bivar(&[(1, 0, 1)])))
}

#[test]
fn kernel_handles_empty_requests_short_prefixes_and_missing_seeds() {
    let recurrence = fibonacci_polynomial_recurrence();
    let seeds = vec![row(&[1]), row(&[0, 1])];

    // Zero rows and prefixes of the seeds evaluate nothing and are returned as is.
    assert_eq!(recurrence.generate_rows_rational(&seeds, 1, 0), Ok(vec![]));
    assert_eq!(
        recurrence.generate_rows_rational(&seeds, 1, 1),
        Ok(vec![row(&[1])])
    );
    assert_eq!(
        recurrence.generate_rows_rational(&seeds, 1, 2),
        Ok(seeds.clone())
    );
    let mut rows = recurrence.generate_rows_rational_iter(&seeds, 1);
    assert_eq!(rows.by_ref().take(0).count(), 0);
    assert_eq!(rows.rows_yielded(), 0);

    // A seed prefix is allowed even when the seeds alone are insufficient ...
    assert_eq!(
        recurrence.generate_rows_rational(&seeds[..1], 1, 1),
        Ok(vec![row(&[1])])
    );
    // ... but the first computed row is not: the seed is yielded, then one error.
    let mut rows = recurrence.generate_rows_rational_iter(&seeds[..1], 1);
    assert_eq!(rows.next(), Some(Ok(row(&[1]))));
    let missing = RecurrenceEvaluationError::InsufficientInitialRows {
        required: 2,
        provided: 1,
    };
    assert_eq!(rows.next(), Some(Err(missing.clone())));
    assert_eq!(rows.next(), None);
    assert_eq!(
        recurrence.generate_rows_rational(&seeds[..1], 1, 2),
        Err(missing)
    );
}

#[test]
fn kernel_reports_index_overflow_and_zero_offset_terms_like_the_full_history_loop() {
    // Index overflow at the first computed row.
    let seeds = [row(&[1])];
    let one_term = factorial_recurrence();
    let mut rows = one_term.generate_rows_rational_iter(&seeds, usize::MAX);
    assert_eq!(rows.next_index(), Some(usize::MAX));
    assert_eq!(rows.next(), Some(Ok(row(&[1]))));
    assert_eq!(rows.next_index(), None);
    assert_eq!(
        rows.next(),
        Some(Err(RecurrenceEvaluationError::IndexOverflow {
            first_index: usize::MAX,
            row_offset: 1,
        }))
    );
    assert_eq!(rows.next(), None);

    // A term with offset 0 refers to the row being computed.
    let self_referential = recurrence(
        vec![RecurrenceTerm::new(0, 0, bivar(&[(0, 0, 1)]))],
        None,
        None,
    );
    let error = self_referential.generate_rows_rational(&seeds, 4, 3);
    assert_eq!(error, naive_rows(&self_referential, &seeds, 4, 3));
    assert!(matches!(
        error,
        Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
            n: 5,
            offset: 0,
            first_index: 4
        })
    ));

    // A denominator that does not divide exactly.
    let non_polynomial = recurrence(
        vec![RecurrenceTerm::new(1, 0, bivar(&[(0, 0, 1)]))],
        Some(bivar(&[(0, 1, 1)])),
        None,
    );
    assert_eq!(
        non_polynomial.generate_rows_rational(&seeds, 0, 3),
        Err(RecurrenceEvaluationError::NonPolynomialQuotient)
    );
    assert_eq!(
        non_polynomial.generate_rows_rational(&seeds, 0, 3),
        naive_rows(&non_polynomial, &seeds, 0, 3)
    );
}

#[test]
fn kernel_matches_the_naive_full_history_loop_including_derivative_and_sign_terms() {
    // A008292-type Eulerian step with derivative, plus an alternating term and
    // an inhomogeneous term, so every branch of the evaluator is exercised.
    let mixed = recurrence(
        vec![
            RecurrenceTerm::new(1, 0, bivar(&[(0, 0, 1), (1, 1, 1), (0, 1, -1)])),
            RecurrenceTerm::new(1, 1, bivar(&[(0, 1, 1), (0, 2, -1)])),
            RecurrenceTerm::new(2, 0, bivar(&[(0, 1, 1)])).with_sign(RecurrenceSign::AlternatingN),
        ],
        None,
        Some(bivar(&[(1, 0, 1)])),
    );
    let seeds = vec![row(&[1]), row(&[1, 1])];
    for first_index in [0, 1, 7] {
        for total in [0, 1, 2, 3, 12] {
            assert_eq!(
                mixed.generate_rows_rational(&seeds, first_index, total),
                naive_rows(&mixed, &seeds, first_index, total),
                "first_index {first_index}, total {total}"
            );
        }
    }
}

#[test]
fn kernel_clone_continues_independently() {
    let mut rows = factorial_recurrence().generate_rows_rational_iter(&[row(&[1])], 0);
    rows.by_ref().take(3).for_each(drop);
    let mut copy = rows.clone();
    assert_eq!(rows.next().unwrap().unwrap(), row(&[6]));
    assert_eq!(copy.next().unwrap().unwrap(), row(&[6]));
    assert_eq!(copy.next().unwrap().unwrap(), row(&[24]));
    assert_eq!(rows.rows_yielded(), 4);
    assert_eq!(copy.rows_yielded(), 5);
}

// ---------------------------------------------------------------------------
// OEIS adapter
// ---------------------------------------------------------------------------

#[test]
fn oeis_rows_prefix_pause_and_continue_with_independent_values() {
    // A008292, the Eulerian triangle: 1; 1, 1; 1, 4, 1; 1, 11, 11, 1;
    // 1, 26, 66, 26, 1; 1, 57, 302, 302, 57, 1.
    let sequence = by_id("A008292").unwrap();
    let mut rows = sequence.generate_rows_iter().unwrap();
    assert_eq!(rows.next_displayed_row(), Some(1));
    let first: Vec<_> = rows.by_ref().take(3).collect::<Result<_, _>>().unwrap();
    assert_eq!(first, vec![ints(&[1]), ints(&[1, 1]), ints(&[1, 4, 1])]);
    assert_eq!(rows.next_displayed_row(), Some(4));
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[1, 11, 11, 1]));
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[1, 26, 66, 26, 1]));
    assert_eq!(
        rows.next().unwrap().unwrap(),
        ints(&[1, 57, 302, 302, 57, 1])
    );

    // The batch APIs return the same rows.
    assert_eq!(
        sequence.generate_rows(6).unwrap()[5],
        ints(&[1, 57, 302, 302, 57, 1])
    );
}

#[test]
fn oeis_rows_from_a_displayed_row_use_displayed_indices() {
    let sequence = by_id("A008292").unwrap();
    assert_eq!(sequence.first_row, 1);
    let mut rows = sequence.generate_rows_from_iter(4).unwrap();
    assert_eq!(rows.next_displayed_row(), Some(4));
    assert_eq!(rows.recurrence_steps_evaluated(), 0, "skipping is deferred");
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[1, 11, 11, 1]));
    assert_eq!(rows.next_displayed_row(), Some(5));
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[1, 26, 66, 26, 1]));
    assert_eq!(
        sequence.generate_rows_from(5, 2).unwrap(),
        vec![ints(&[1, 26, 66, 26, 1]), ints(&[1, 57, 302, 302, 57, 1])]
    );
    assert_eq!(
        sequence.generate_rows_from_iter(0).unwrap_err(),
        OeisSequenceError::RowBeforeStart {
            requested: 0,
            first: 1
        }
    );
    assert_eq!(
        sequence.generate_rows_from(0, 3).unwrap_err(),
        OeisSequenceError::RowBeforeStart {
            requested: 0,
            first: 1
        }
    );
    // An empty request is empty and a row that already exists starts anywhere.
    assert_eq!(
        sequence.generate_rows(0).unwrap(),
        Vec::<Vec<BigInt>>::new()
    );
    assert_eq!(sequence.generate_rows_from(7, 0).unwrap().len(), 0);
}

#[test]
fn oeis_rows_stop_generating_when_consumption_stops() {
    for sequence in catalog() {
        let prefix = sequence.prefix_rows.len();
        let seeds = sequence.initial_rows.len();
        let (recurrence, _) = sequence.recurrence_parts().unwrap();

        // Rows from the authoritative prefix and the initial rows need no step.
        let mut rows = sequence.generate_rows_iter().unwrap();
        rows.by_ref()
            .take(prefix + seeds)
            .for_each(|row| drop(row.unwrap()));
        assert_eq!(rows.recurrence_steps_evaluated(), 0, "{}", sequence.id);

        // Then exactly one step per further row, with bounded retention.
        for step in 1..=4 {
            rows.next()
                .unwrap()
                .unwrap_or_else(|e| panic!("{}: {e}", sequence.id));
            assert_eq!(rows.recurrence_steps_evaluated(), step, "{}", sequence.id);
            assert!(
                rows.retained_row_count() <= recurrence.max_offset().max(seeds),
                "{}",
                sequence.id
            );
        }
    }
}

#[test]
fn oeis_zero_and_short_requests_do_not_touch_the_recurrence() {
    // A definition whose recurrence cannot even be decoded.  Anything inside
    // the authoritative prefix must still work, and only the first row beyond
    // it fails, once, then the iterator ends.
    static BAD_COEFFICIENT: [SparseCoefficient; 1] = [SparseCoefficient::new(0, 0, "not-a-number")];
    static BAD_TERMS: [SparseRecurrenceTerm; 1] = [SparseRecurrenceTerm {
        offset: 1,
        derivative_order: 0,
        alternating_sign: false,
        coefficient: &BAD_COEFFICIENT,
    }];
    let corrupt = OeisSequenceDefinition {
        prefix_rows: &[&["1"], &["2", "3"]],
        initial_rows: &[&["1"]],
        recurrence_first_index: 0,
        recurrence: SparseRecurrence {
            terms: &BAD_TERMS,
            denominator: None,
            inhomogeneous: None,
        },
        ..*by_id("A008292").unwrap()
    };
    assert_eq!(corrupt.generate_rows(0).unwrap().len(), 0);
    assert_eq!(corrupt.generate_rows(1).unwrap(), vec![ints(&[1])]);
    assert_eq!(corrupt.generate_rows(2).unwrap().len(), 2);
    assert!(matches!(
        corrupt.generate_rows(3),
        Err(OeisSequenceError::InvalidDefinition { .. })
    ));

    let mut rows = corrupt.generate_rows_iter().unwrap();
    assert_eq!(rows.by_ref().take(2).count(), 2);
    assert!(matches!(
        rows.next(),
        Some(Err(OeisSequenceError::InvalidDefinition { .. }))
    ));
    assert!(rows.next().is_none());
    assert!(rows.next().is_none());
}

fn custom_definition(
    terms: &'static [SparseRecurrenceTerm],
    denominator: Option<&'static [SparseCoefficient]>,
) -> OeisSequenceDefinition {
    OeisSequenceDefinition {
        layout: OeisLayout::RegularTriangle,
        first_row: 10,
        bfile_prefix_verified: true,
        prefix_rows: &[&["7", "8"]],
        recurrence_first_index: 0,
        recurrence_first_width: 3,
        initial_rows: &[&["1"]],
        recurrence: SparseRecurrence {
            terms,
            denominator,
            inhomogeneous: None,
        },
        ..*by_id("A008292").unwrap()
    }
}

#[test]
fn oeis_rows_pad_regular_triangle_rows_but_not_the_prefix() {
    // P_n = P_{n-1}: the recurrence rows are all `1`, padded to the structural
    // width 3, 4, 5, ... of the triangle; the prefix row keeps its own length.
    static ONE: [SparseCoefficient; 1] = [SparseCoefficient::new(0, 0, "1")];
    static TERMS: [SparseRecurrenceTerm; 1] = [SparseRecurrenceTerm {
        offset: 1,
        derivative_order: 0,
        alternating_sign: false,
        coefficient: &ONE,
    }];
    let sequence = custom_definition(&TERMS, None);
    let rows: Vec<_> = sequence
        .generate_rows_iter()
        .unwrap()
        .take(4)
        .collect::<Result<_, _>>()
        .unwrap();
    assert_eq!(
        rows,
        vec![
            ints(&[7, 8]),
            ints(&[1, 0, 0]),
            ints(&[1, 0, 0, 0]),
            ints(&[1, 0, 0, 0, 0])
        ]
    );
    assert_eq!(rows, sequence.generate_rows(4).unwrap());
    // Without the verified b-file mapping there is no padding.
    let unpadded = OeisSequenceDefinition {
        bfile_prefix_verified: false,
        ..sequence
    };
    assert_eq!(unpadded.generate_rows(3).unwrap()[2], ints(&[1]));
}

#[test]
fn oeis_rows_report_nonintegral_coefficients_with_displayed_rows() {
    // P_n = P_{n-1} / 2 from P_0 = 1: row P_1 = 1/2 is displayed row 12
    // (first_row 10, one prefix row, so P_0 is displayed row 11).
    static HALF: [SparseCoefficient; 1] = [SparseCoefficient::new(0, 0, "1/2")];
    static TERMS: [SparseRecurrenceTerm; 1] = [SparseRecurrenceTerm {
        offset: 1,
        derivative_order: 0,
        alternating_sign: false,
        coefficient: &HALF,
    }];
    let sequence = custom_definition(&TERMS, None);
    let id = sequence.id;
    let mut rows = sequence.generate_rows_iter().unwrap();
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[7, 8]));
    assert_eq!(rows.next().unwrap().unwrap(), ints(&[1, 0, 0]));
    let error = OeisSequenceError::NonIntegralCoefficient { id, row: 12 };
    assert_eq!(rows.next(), Some(Err(error.clone())));
    assert_eq!(rows.next(), None, "ends after the typed error");
    assert_eq!(rows.next(), None);
    assert_eq!(rows.retained_row_count(), 0);
    assert_eq!(sequence.generate_rows(3), Err(error.clone()));
    assert_eq!(sequence.generate_rows(2).unwrap().len(), 2);
    // Skipped rows are still validated, as before.
    assert_eq!(sequence.generate_rows_from(13, 0), Err(error.clone()));
    assert_eq!(
        sequence.generate_rows_from_iter(13).unwrap().next(),
        Some(Err(error))
    );
}

#[test]
fn oeis_later_step_error_does_not_prevent_consuming_the_earlier_prefix() {
    // (n - 3) P_n = (n - 3) P_{n-1}: P_n = 1 until n = 3, where the
    // denominator vanishes.
    static COEFFICIENT: [SparseCoefficient; 2] = [
        SparseCoefficient::new(1, 0, "1"),
        SparseCoefficient::new(0, 0, "-3"),
    ];
    static TERMS: [SparseRecurrenceTerm; 1] = [SparseRecurrenceTerm {
        offset: 1,
        derivative_order: 0,
        alternating_sign: false,
        coefficient: &COEFFICIENT,
    }];
    let sequence = OeisSequenceDefinition {
        bfile_prefix_verified: false,
        ..custom_definition(&TERMS, Some(&COEFFICIENT))
    };
    let mut rows = sequence.generate_rows_iter().unwrap();
    let valid: Vec<_> = rows.by_ref().take(4).collect::<Result<_, _>>().unwrap();
    // prefix, P_0, P_1, P_2
    assert_eq!(
        valid,
        vec![ints(&[7, 8]), ints(&[1]), ints(&[1]), ints(&[1])]
    );
    assert_eq!(rows.recurrence_steps_evaluated(), 2);
    assert_eq!(
        rows.next(),
        Some(Err(OeisSequenceError::Evaluation(
            RecurrenceEvaluationError::ZeroDenominator
        )))
    );
    assert_eq!(rows.next(), None);
    assert_eq!(sequence.generate_rows(4).unwrap(), valid);
    assert_eq!(
        sequence.generate_rows(5),
        Err(OeisSequenceError::Evaluation(
            RecurrenceEvaluationError::ZeroDenominator
        ))
    );
}

#[test]
fn every_catalog_entry_matches_a_naive_full_history_evaluation() {
    // Independent of the iterator: replay each exported recurrence with the
    // old whole-history loop, convert, and compare with the lazy rows past the
    // authoritative prefix and the initial rows.
    for sequence in catalog() {
        let (recurrence, seeds, first_index) = sequence
            .recurrence_export_parts()
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        let count = seeds.len() + 6;
        let naive = naive_rows(&recurrence, &seeds, first_index as usize, count)
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        let lazy: Vec<_> = sequence
            .generate_rows_iter()
            .unwrap()
            .take(count)
            .collect::<Result<_, _>>()
            .unwrap_or_else(|error| panic!("{}: {error}", sequence.id));
        assert_eq!(lazy.len(), count, "{}", sequence.id);
        for (index, (naive_row, lazy_row)) in naive.iter().zip(&lazy).enumerate() {
            let mut naive_integers: Vec<BigInt> =
                naive_row.iter().map(|c| c.to_integer()).collect();
            let mut lazy_integers = lazy_row.clone();
            for trimmed in [&mut naive_integers, &mut lazy_integers] {
                while trimmed.last().is_some_and(|c| *c == BigInt::from(0)) {
                    trimmed.pop();
                }
            }
            assert_eq!(naive_integers, lazy_integers, "{} row {index}", sequence.id);
        }
    }
}

#[test]
fn every_catalog_triangle_row_has_its_structural_width() {
    let mut padded_somewhere = false;
    for sequence in catalog() {
        if !(sequence.bfile_prefix_verified && sequence.layout == OeisLayout::RegularTriangle) {
            continue;
        }
        let prefix = sequence.prefix_rows.len();
        let rows = sequence.generate_rows(prefix + 12).unwrap();
        for (index, row) in rows.iter().enumerate().skip(prefix) {
            assert_eq!(
                row.len(),
                sequence.recurrence_first_width + index - prefix,
                "{} row {index}",
                sequence.id
            );
            padded_somewhere |= row.last() == Some(&BigInt::from(0));
        }
    }
    assert!(padded_somewhere, "no catalog row exercises zero padding");
}

#[test]
fn bfile_offsets_do_not_retain_the_skipped_rows() {
    // A008292: complete rows 1..=3 have 1 + 2 + 3 coefficients.
    let sequence = by_id("A008292").unwrap();
    assert_eq!(sequence.flattened_offset_for_row(1).unwrap(), 1);
    assert_eq!(sequence.flattened_offset_for_row(4).unwrap(), 7);
    assert_eq!(sequence.flattened_offset_for_row(7).unwrap(), 22);
}

// ---------------------------------------------------------------------------
// Named families
// ---------------------------------------------------------------------------

#[test]
fn named_family_iterators_match_their_independent_integer_generators() {
    let n = 40;
    type Rows = Vec<Vec<BigInt>>;
    let cases: Vec<(&str, Rows, Rows)> = vec![
        (
            "eulerian",
            sequences::eulerian_polynomials_iter().take(n + 1).collect(),
            sequences::eulerian_polynomials_bigint(n),
        ),
        (
            "type B Eulerian",
            sequences::type_b_eulerian_polynomials_iter()
                .take(n + 1)
                .collect(),
            sequences::type_b_eulerian_polynomials_bigint(n),
        ),
        (
            "Chebyshev T",
            sequences::chebyshev_polynomials_t_iter()
                .take(n + 1)
                .collect(),
            sequences::chebyshev_polynomials_t_bigint(n),
        ),
        (
            "Chebyshev U",
            sequences::chebyshev_polynomials_u_iter()
                .take(n + 1)
                .collect(),
            sequences::chebyshev_polynomials_u_bigint(n),
        ),
        (
            "Hermite",
            sequences::hermite_polynomials_iter().take(n + 1).collect(),
            sequences::hermite_polynomials_bigint(n),
        ),
        (
            "Delannoy",
            literature::delannoy_polynomials_iter()
                .take(n + 1)
                .collect(),
            literature::delannoy_polynomials_bigint(n),
        ),
        (
            "Delannoy square",
            literature::delannoy_square_polynomials_iter()
                .take(n + 1)
                .collect(),
            literature::delannoy_square_polynomials_bigint(n),
        ),
    ];
    for (name, lazy, batch) in cases {
        assert_eq!(lazy, batch, "{name}");
    }
}

#[test]
fn named_family_iterators_have_independently_known_values() {
    // A008292 row 4 and A060187 row 3 (type B Eulerian), A053120 (T_4 =
    // 8t^4 - 8t^2 + 1), A053117 (U_3 = 8t^3 - 4t), A099174 (He_4 = t^4 - 6t^2 + 3),
    // A008288 (Delannoy, D_3 = 1 + 5x + 5x^2 + x^3).
    assert_eq!(
        sequences::eulerian_polynomials_iter().nth(4).unwrap(),
        ints(&[1, 11, 11, 1])
    );
    assert_eq!(
        sequences::type_b_eulerian_polynomials_iter()
            .nth(3)
            .unwrap(),
        ints(&[1, 23, 23, 1])
    );
    assert_eq!(
        sequences::chebyshev_polynomials_t_iter().nth(4).unwrap(),
        ints(&[1, 0, -8, 0, 8])
    );
    assert_eq!(
        sequences::chebyshev_polynomials_u_iter().nth(3).unwrap(),
        ints(&[0, -4, 0, 8])
    );
    assert_eq!(
        sequences::hermite_polynomials_iter().nth(4).unwrap(),
        ints(&[3, 0, -6, 0, 1])
    );
    assert_eq!(
        literature::delannoy_polynomials_iter().nth(3).unwrap(),
        ints(&[1, 5, 5, 1])
    );
    // The first computed Delannoy-square row, expanded by hand from the
    // four-term recurrence (Mao--Wang):
    // G_4 = (x+3) G_3 + (x-1) G_2 + (x-1) G_1 + x G_0 = 29 + 88x + 61x^2 + 14x^3 + x^4.
    assert_eq!(
        literature::delannoy_square_polynomials_iter()
            .nth(4)
            .unwrap(),
        ints(&[29, 88, 61, 14, 1])
    );
}

#[test]
fn named_family_iterators_continue_and_keep_a_bounded_window() {
    let mut rows = sequences::chebyshev_polynomials_t_iter();
    let first: Vec<_> = rows.by_ref().take(5).collect();
    assert_eq!(rows.rows_yielded(), 5);
    assert_eq!(
        rows.next().unwrap(),
        sequences::chebyshev_polynomials_t_bigint(5)[5]
    );
    assert_eq!(first, sequences::chebyshev_polynomials_t_bigint(4));
    for _ in 0..100 {
        rows.next().unwrap();
        assert!(rows.retained_row_count() <= 2);
    }
    let mut eulerian = sequences::eulerian_polynomials_iter();
    eulerian.by_ref().take(30).for_each(drop);
    assert_eq!(eulerian.retained_row_count(), 1);
}
