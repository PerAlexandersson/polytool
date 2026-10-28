//! Lazy, stateful row generation for a scalar recurrence.
//!
//! [`RecurrenceRowIterator`] is the pull-based counterpart of
//! [`Recurrence::generate_rows_rational`]: it produces one row `P_n` per call
//! to [`Iterator::next`], evaluating the recurrence step only when that row is
//! requested.  The batch function collects from this iterator, so both agree
//! row for row.

use super::*;
use std::collections::VecDeque;
use std::iter::FusedIterator;

/// A lazy iterator over the rows `P_{first_index}, P_{first_index + 1}, ...`
/// of a [`Recurrence`], with exact rational coefficients.
///
/// # Contract
///
/// * **Rows yielded.** The supplied seed rows (`initial_rows`) are yielded
///   first, in order, as the rows `P_{first_index}, ...`; every later row is
///   computed with [`Recurrence::evaluate_next_rational`].  Collecting the first
///   `k` items therefore gives exactly the `k` rows of
///   [`Recurrence::generate_rows_rational`].  The stream is infinite unless a
///   step fails.
/// * **Index.** The first item has absolute recurrence index `first_index`;
///   the index `n` substituted into the recurrence coefficients for a computed
///   row is its own absolute index (`first_index + seed_count + steps`).
/// * **Laziness.** Construction does no evaluation.  Each `next` evaluates at
///   most one recurrence step, and none while seeds remain.  Dropping the
///   iterator, or no longer calling `next`, stops all further evaluation.
/// * **Continuation.** The iterator keeps its position: `by_ref().take(k)`
///   followed by further calls continues with row `k` without recomputing
///   earlier rows.
/// * **Bounded history.** Only the last `max_offset()` rows (the
///   [`history_limit`](Self::history_limit)) are retained for evaluation, in
///   addition to the seeds that have not been yielded yet.  Polynomial sizes
///   themselves may grow.
/// * **Errors.** Construction cannot fail.  Every failure is a *step-time*
///   [`RecurrenceEvaluationError`], produced when the offending row is
///   requested (after any seeds): too few seeds
///   ([`InsufficientInitialRows`](RecurrenceEvaluationError::InsufficientInitialRows),
///   reported after the seeds were yielded), a vanishing or non-dividing
///   denominator, or index overflow.  The error is yielded exactly once; after
///   it the iterator is exhausted and keeps returning `None`
///   ([`FusedIterator`]).  Rows before the failing step remain consumable.
///
/// ```
/// use polytool::recurrence::{BigRational, BivarPoly, Recurrence, RecurrenceTerm};
///
/// // P_n = 2 P_{n-1}, with P_0 = 1.
/// let two = BigRational::from_integer(2.into());
/// let recurrence = Recurrence {
///     terms: vec![RecurrenceTerm::new(1, 0, BivarPoly::from_coefficients_normalized(vec![vec![two]]))],
///     denominator: None,
///     inhomogeneous: None,
/// };
/// let seed = vec![vec![BigRational::from_integer(1.into())]];
/// let mut rows = recurrence.generate_rows_rational_iter(&seed, 0);
/// let first_three: Vec<_> = rows.by_ref().take(3).collect::<Result<_, _>>().unwrap();
/// assert_eq!(first_three[2], vec![BigRational::from_integer(4.into())]);
/// // Continue with the same iterator: this is P_3 = 8, no earlier row is recomputed.
/// assert_eq!(rows.next().unwrap().unwrap(), vec![BigRational::from_integer(8.into())]);
/// assert_eq!(rows.steps_evaluated(), 3);
/// ```
#[derive(Debug, Clone)]
pub struct RecurrenceRowIterator {
    recurrence: Recurrence,
    first_index: usize,
    history_limit: usize,
    seed_count: usize,
    pending_seeds: VecDeque<Vec<BigRational>>,
    history: VecDeque<Vec<BigRational>>,
    rows_yielded: usize,
    steps_evaluated: usize,
    finished: bool,
}

impl RecurrenceRowIterator {
    /// Create an iterator over the rows of `recurrence` starting at absolute
    /// index `first_index`, whose first rows are `initial_rows`.
    ///
    /// Infallible and lazy; see the [type-level contract](Self).
    pub fn new(
        recurrence: Recurrence,
        initial_rows: Vec<Vec<BigRational>>,
        first_index: usize,
    ) -> Self {
        let history_limit = recurrence.max_offset();
        Self {
            recurrence,
            first_index,
            history_limit,
            seed_count: initial_rows.len(),
            pending_seeds: initial_rows.into(),
            history: VecDeque::new(),
            rows_yielded: 0,
            steps_evaluated: 0,
            finished: false,
        }
    }

    /// Absolute recurrence index of the first row (the index of the first seed).
    pub fn first_index(&self) -> usize {
        self.first_index
    }

    /// Absolute index of the row the next call to `next` produces, or `None`
    /// if it does not fit in `usize`.
    pub fn next_index(&self) -> Option<usize> {
        self.first_index.checked_add(self.rows_yielded)
    }

    /// Number of rows yielded so far (seeds included, the terminal error not).
    pub fn rows_yielded(&self) -> usize {
        self.rows_yielded
    }

    /// Number of recurrence steps evaluated so far, that is, the number of
    /// yielded rows that were computed rather than supplied as seeds.  Failed
    /// steps are not counted.
    pub fn steps_evaluated(&self) -> usize {
        self.steps_evaluated
    }

    /// The history bound: the largest recurrence offset, which is the number of
    /// previously yielded rows kept for evaluation.
    pub fn history_limit(&self) -> usize {
        self.history_limit
    }

    /// Number of rows currently retained: the history window (at most
    /// [`history_limit`](Self::history_limit)) plus seeds not yet yielded.
    /// Zero after a terminal error.
    pub fn retained_row_count(&self) -> usize {
        self.history.len() + self.pending_seeds.len()
    }

    fn remember(&mut self, row: &[BigRational]) {
        if self.history_limit == 0 {
            return;
        }
        if self.history.len() == self.history_limit {
            self.history.pop_front();
        }
        self.history.push_back(row.to_vec());
    }

    fn evaluate_step(&mut self) -> Result<Vec<BigRational>, RecurrenceEvaluationError> {
        if self.seed_count < self.history_limit {
            return Err(RecurrenceEvaluationError::InsufficientInitialRows {
                required: self.history_limit,
                provided: self.seed_count,
            });
        }
        let next_index = self
            .next_index()
            .ok_or(RecurrenceEvaluationError::IndexOverflow {
                first_index: self.first_index,
                row_offset: self.rows_yielded,
            })?;
        // `history` holds the rows `P_{next_index - len}, ..., P_{next_index - 1}`.
        let window_first_index = next_index - self.history.len();
        let window = self.history.make_contiguous();
        self.recurrence
            .evaluate_next_rational(window, window_first_index)
            .map_err(|error| match error {
                // Report indices relative to the whole stream, as the batch
                // evaluation over all rows does.
                RecurrenceEvaluationError::OffsetBeforeInitialRow { n, offset, .. } => {
                    RecurrenceEvaluationError::OffsetBeforeInitialRow {
                        n,
                        offset,
                        first_index: self.first_index,
                    }
                }
                other => other,
            })
    }
}

impl Iterator for RecurrenceRowIterator {
    type Item = Result<Vec<BigRational>, RecurrenceEvaluationError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.finished {
            return None;
        }
        if let Some(seed) = self.pending_seeds.pop_front() {
            self.remember(&seed);
            self.rows_yielded += 1;
            return Some(Ok(seed));
        }
        match self.evaluate_step() {
            Ok(row) => {
                self.remember(&row);
                self.rows_yielded += 1;
                self.steps_evaluated += 1;
                Some(Ok(row))
            }
            Err(error) => {
                self.finished = true;
                self.history.clear();
                Some(Err(error))
            }
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        if self.finished {
            (0, Some(0))
        } else {
            (self.pending_seeds.len(), None)
        }
    }
}

impl FusedIterator for RecurrenceRowIterator {}

impl Recurrence {
    /// Lazily iterate over the rows `P_{first_index}, P_{first_index + 1}, ...`.
    ///
    /// `initial_rows[0]` is `P(first_index)`; the seeds are yielded first and
    /// the recurrence is evaluated only for later rows.  This clones the
    /// recurrence and the seeds; see [`RecurrenceRowIterator`] for the laziness,
    /// bounded-history and error contract.  [`Self::generate_rows_rational`]
    /// collects a finite prefix of this iterator.
    pub fn generate_rows_rational_iter(
        &self,
        initial_rows: &[Vec<BigRational>],
        first_index: usize,
    ) -> RecurrenceRowIterator {
        RecurrenceRowIterator::new(self.clone(), initial_rows.to_vec(), first_index)
    }
}
