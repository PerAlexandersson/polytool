//! Descent and peak polynomials for permutations avoiding an increasing pattern.
//!
//! The generators use the Robinson--Schensted correspondence rather than
//! enumerating permutations.  Avoiding `12...m` is equivalent to the first row
//! of the RSK shape having length at most `m - 1`.  The descent set of a
//! permutation is the descent set of its recording tableau, so adjoining the
//! next largest tableau entry gives a finite-state recurrence.

use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

#[derive(Clone, Copy)]
enum TrackedStatistic {
    Descents,
    Peaks,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct TableauState {
    shape: Vec<usize>,
    maximum_row: usize,
    last_descent: bool,
}

/// Descent polynomials of permutations avoiding the increasing pattern
/// `12...m`, for sizes `0..=max_n`.
///
/// The coefficient of `t^k` in row `n` is the number of permutations in
/// `Av_n(12...m)` with `k` descents.  Coefficients are returned in ascending
/// degree and use arbitrary-precision unsigned integers.  The empty
/// permutation contributes row `[1]`.  For `m = 1`, all positive-size rows are
/// `[0]`.
///
/// This is an exact RSK/tableau recurrence.  A state consists of a partition
/// whose first row has length at most `m - 1`, together with the row containing
/// the current maximum.  On
/// adjoining the next maximum in row `s`, a descent is created precisely when
/// `s` is below the previous maximum's row.  If `p(n)` is the number of
/// partitions of `n` with largest part at most `m - 1`, the state count is
/// `O(n p(n))`; the polynomial updates add another factor linear in `n`.
///
/// # Panics
///
/// Panics if `m = 0`, since the empty word is not a permutation pattern.
pub fn increasing_pattern_avoiding_descent_polynomials(
    max_n: usize,
    m: usize,
) -> Vec<Vec<BigUint>> {
    increasing_pattern_avoiding_polynomials(max_n, m, TrackedStatistic::Descents)
}

/// Ordinary peak polynomials of permutations avoiding the increasing pattern
/// `12...m`, for sizes `0..=max_n`.
///
/// The coefficient of `t^k` in row `n` is the number of permutations in
/// `Av_n(12...m)` with `k` peaks, where a peak is a position
/// `i in {2, ..., n-1}` with `pi(i-1) < pi(i) > pi(i+1)`.  Coefficients are in
/// ascending degree and use arbitrary precision.  Boundary behavior and
/// complexity are as in [`increasing_pattern_avoiding_descent_polynomials`].
/// The state additionally records whether the preceding tableau position was
/// a descent; a new peak occurs exactly on a `0 -> 1` transition after the
/// first position.
///
/// # Panics
///
/// Panics if `m = 0`.
pub fn increasing_pattern_avoiding_peak_polynomials(max_n: usize, m: usize) -> Vec<Vec<BigUint>> {
    increasing_pattern_avoiding_polynomials(max_n, m, TrackedStatistic::Peaks)
}

fn increasing_pattern_avoiding_polynomials(
    max_n: usize,
    m: usize,
    statistic: TrackedStatistic,
) -> Vec<Vec<BigUint>> {
    assert!(m > 0, "the increasing pattern length must be positive");

    let mut rows = Vec::with_capacity(max_n + 1);
    rows.push(vec![BigUint::one()]);
    if max_n == 0 {
        return rows;
    }
    if m == 1 {
        rows.resize(max_n + 1, vec![BigUint::zero()]);
        return rows;
    }

    let initial = TableauState {
        shape: vec![1],
        maximum_row: 0,
        last_descent: false,
    };
    let mut states = BTreeMap::from([(initial, vec![BigUint::one()])]);
    rows.push(aggregate_rsk_pairs(&states));

    for size in 1..max_n {
        let mut next = BTreeMap::<TableauState, Vec<BigUint>>::new();
        for (state, polynomial) in &states {
            for new_row in 0..=state.shape.len() {
                if new_row < state.shape.len() && state.shape[new_row] == m - 1 {
                    continue;
                }
                if new_row < state.shape.len()
                    && new_row > 0
                    && state.shape[new_row] == state.shape[new_row - 1]
                {
                    continue;
                }
                let creates_descent = new_row > state.maximum_row;
                let shift = match statistic {
                    TrackedStatistic::Descents => usize::from(creates_descent),
                    TrackedStatistic::Peaks => {
                        usize::from(size >= 2 && !state.last_descent && creates_descent)
                    }
                };
                let mut next_shape = state.shape.clone();
                if new_row == next_shape.len() {
                    next_shape.push(1);
                } else {
                    next_shape[new_row] += 1;
                }
                let next_state = TableauState {
                    shape: next_shape,
                    maximum_row: new_row,
                    last_descent: creates_descent,
                };
                add_shifted(next.entry(next_state).or_default(), polynomial, shift);
            }
        }
        states = next;
        rows.push(aggregate_rsk_pairs(&states));
    }
    rows
}

fn add_shifted(target: &mut Vec<BigUint>, source: &[BigUint], shift: usize) {
    target.resize(target.len().max(source.len() + shift), BigUint::zero());
    for (degree, coefficient) in source.iter().enumerate() {
        target[degree + shift] += coefficient;
    }
}

fn add_assign(target: &mut Vec<BigUint>, source: &[BigUint]) {
    target.resize(target.len().max(source.len()), BigUint::zero());
    for (degree, coefficient) in source.iter().enumerate() {
        target[degree] += coefficient;
    }
}

fn aggregate_rsk_pairs(states: &BTreeMap<TableauState, Vec<BigUint>>) -> Vec<BigUint> {
    let mut recording_by_shape = BTreeMap::<&[usize], Vec<BigUint>>::new();
    for (state, polynomial) in states {
        add_assign(
            recording_by_shape.entry(&state.shape).or_default(),
            polynomial,
        );
    }

    let mut result = Vec::new();
    for recording_polynomial in recording_by_shape.values() {
        let insertion_tableaux: BigUint = recording_polynomial.iter().sum();
        let weighted = recording_polynomial
            .iter()
            .map(|coefficient| coefficient * &insertion_tableaux)
            .collect::<Vec<_>>();
        add_assign(&mut result, &weighted);
    }
    if result.is_empty() {
        vec![BigUint::zero()]
    } else {
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::permutation::avoiding_permutations;
    use crate::statistics::{compute, Stat};

    fn brute_polynomial(n: u8, pattern: &[u8], statistic: Stat) -> Vec<BigUint> {
        let mut polynomial = vec![BigUint::zero(); usize::from(n) + 1];
        for permutation in avoiding_permutations(n, &[pattern.to_vec()]) {
            polynomial[compute(&permutation, statistic)] += 1u8;
        }
        while polynomial.len() > 1 && polynomial.last().is_some_and(BigUint::is_zero) {
            polynomial.pop();
        }
        polynomial
    }

    #[test]
    fn increasing_patterns_of_lengths_four_and_five_match_direct_avoidance() {
        for pattern in [&[1, 2, 3, 4][..], &[1, 2, 3, 4, 5][..]] {
            let descents = increasing_pattern_avoiding_descent_polynomials(8, pattern.len());
            let peaks = increasing_pattern_avoiding_peak_polynomials(8, pattern.len());
            for n in 0..=8 {
                assert_eq!(descents[n], brute_polynomial(n as u8, pattern, Stat::Des));
                assert_eq!(peaks[n], brute_polynomial(n as u8, pattern, Stat::Peak));
            }
        }
    }

    #[test]
    fn total_counts_for_1234_avoidance_are_known() {
        let rows = increasing_pattern_avoiding_descent_polynomials(10, 4);
        let totals = rows
            .iter()
            .map(|row| row.iter().sum::<BigUint>())
            .collect::<Vec<_>>();
        let expected = [1u64, 1, 2, 6, 23, 103, 513, 2761, 15767, 94359, 586590];
        assert_eq!(
            totals,
            expected.into_iter().map(BigUint::from).collect::<Vec<_>>()
        );
    }

    #[test]
    fn increasing_pattern_length_two_leaves_only_decreasing_permutation() {
        let descents = increasing_pattern_avoiding_descent_polynomials(6, 2);
        let peaks = increasing_pattern_avoiding_peak_polynomials(6, 2);
        for n in 1..=6 {
            let mut expected_descents = vec![BigUint::zero(); n];
            expected_descents[n - 1] = BigUint::one();
            assert_eq!(descents[n], expected_descents);
            assert_eq!(peaks[n], vec![BigUint::one()]);
        }
    }

    #[test]
    fn increasing_pattern_length_one_has_no_nonempty_avoiders() {
        assert_eq!(
            increasing_pattern_avoiding_descent_polynomials(3, 1),
            vec![
                vec![BigUint::one()],
                vec![BigUint::zero()],
                vec![BigUint::zero()],
                vec![BigUint::zero()]
            ]
        );
        assert_eq!(
            increasing_pattern_avoiding_peak_polynomials(3, 1),
            vec![
                vec![BigUint::one()],
                vec![BigUint::zero()],
                vec![BigUint::zero()],
                vec![BigUint::zero()]
            ]
        );
    }

    #[test]
    #[should_panic(expected = "pattern length must be positive")]
    fn increasing_pattern_length_zero_is_rejected() {
        let _ = increasing_pattern_avoiding_descent_polynomials(3, 0);
    }
}
