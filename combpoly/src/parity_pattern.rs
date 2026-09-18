//! Exact counting for parity-oriented permutations avoiding `321` or `123`.
//!
//! The implementation is count-only: it memoizes boundary states and never
//! constructs an output permutation.  See [`count_parity_oriented_321_avoiding`]
//! for the state invariant and recurrence.

use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::HashMap;

/// The parity relation required between every position `i` and value `pi(i)`.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ParityOrientation {
    /// `pi(i)` and `i` have the same parity.
    Preserving,
    /// `pi(i)` and `i` have opposite parity.
    Reversing,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct FillState {
    parity_word: u64,
    len: u8,
    legal_low: u8,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct SubsequenceState {
    parity_word: u64,
    len: u8,
    next_parity: u8,
    length_parity: u8,
}

/// One aggregated way to delete an alternating subsequence.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct SubsequenceResult {
    selected_len: u8,
    skipped_word: u64,
}

#[derive(Default)]
struct Counter {
    fill_cache: HashMap<FillState, BigUint>,
    subsequence_cache: HashMap<SubsequenceState, Vec<(SubsequenceResult, u64)>>,
}

impl Counter {
    fn fill(&mut self, parity_word: u64, len: u8, legal_low: u8) -> BigUint {
        if len == 0 {
            return BigUint::one();
        }
        debug_assert!(legal_low < len);
        debug_assert!(parity_word < (1u64 << len));

        // If no value is legal to the right of the maximum, every remaining
        // value is forced into increasing order.  The same is true for one
        // remaining value.  The target positions begin odd and alternate.
        if legal_low == 0 || len == 1 {
            return if parity_word == alternating_word(len) {
                BigUint::one()
            } else {
                BigUint::zero()
            };
        }

        let state = FillState {
            parity_word,
            len,
            legal_low,
        };
        if let Some(value) = self.fill_cache.get(&state) {
            return value.clone();
        }

        let top_parity = ((parity_word >> (len - 1)) & 1) as u8;
        let top_bit = (top_parity as u64) << (len - 1);
        let required_suffix_parity = (len ^ top_parity) & 1;
        let mut total = BigUint::zero();

        // An empty suffix is possible exactly when its required length is even.
        if required_suffix_parity == 0 {
            total += self.fill(parity_word ^ top_bit, len - 1, legal_low.min(len - 2));
        }

        // Otherwise choose the least value in the increasing suffix, followed
        // by an alternating-parity subsequence of the still-legal interval.
        for first in 0..legal_low {
            if ((parity_word >> first) & 1) as u8 != (top_parity ^ 1) {
                continue;
            }

            let middle_len = legal_low - first - 1;
            let middle_word = low_bits(parity_word >> (first + 1), middle_len);
            let choices = self.alternating_subsequences(
                middle_word,
                middle_len,
                top_parity,
                required_suffix_parity ^ 1,
            );

            for (choice, multiplicity) in choices {
                let low_word = low_bits(parity_word, first);
                let high_word = (parity_word ^ top_bit) >> legal_low;
                let retained_middle = middle_len - choice.selected_len;
                let next_word = low_word
                    | (choice.skipped_word << first)
                    | (high_word << (first + retained_middle));
                let next_len = len - 2 - choice.selected_len;
                // If the suffix used every value above `first`, then all
                // values in the remaining prefix are legal; the next maximum
                // itself is not counted among its legal lower values.
                let next_legal = if next_len == 0 {
                    0
                } else {
                    first.min(next_len - 1)
                };
                let contribution = self.fill(next_word, next_len, next_legal);
                total += contribution * multiplicity;
            }
        }

        self.fill_cache.insert(state, total.clone());
        total
    }

    /// Return `(result, multiplicity)` pairs for parity-alternating
    /// subsequences with the requested first parity and length parity.
    fn alternating_subsequences(
        &mut self,
        parity_word: u64,
        len: u8,
        next_parity: u8,
        length_parity: u8,
    ) -> Vec<(SubsequenceResult, u64)> {
        if len == 0 {
            return if length_parity == 0 {
                vec![(
                    SubsequenceResult {
                        selected_len: 0,
                        skipped_word: 0,
                    },
                    1,
                )]
            } else {
                Vec::new()
            };
        }

        let state = SubsequenceState {
            parity_word,
            len,
            next_parity,
            length_parity,
        };
        if let Some(results) = self.subsequence_cache.get(&state) {
            return results.clone();
        }

        let low_parity = (parity_word & 1) as u8;
        let tail_word = parity_word >> 1;
        let mut aggregated: HashMap<SubsequenceResult, u64> = HashMap::new();

        // Skip the least value; it remains as the least bit of the skipped word.
        for (mut result, multiplicity) in
            self.alternating_subsequences(tail_word, len - 1, next_parity, length_parity)
        {
            result.skipped_word = (result.skipped_word << 1) | low_parity as u64;
            *aggregated.entry(result).or_default() += multiplicity;
        }

        // Select it when it is the required next parity.
        if low_parity == next_parity {
            for (mut result, multiplicity) in self.alternating_subsequences(
                tail_word,
                len - 1,
                next_parity ^ 1,
                length_parity ^ 1,
            ) {
                result.selected_len += 1;
                *aggregated.entry(result).or_default() += multiplicity;
            }
        }

        let results: Vec<_> = aggregated.into_iter().collect();
        // Caching very long result lists dominates memory and gives little
        // reuse.  The bound mirrors the practical split in the original
        // bitmask computation while keeping all arithmetic exact.
        if len <= 35 {
            self.subsequence_cache.insert(state, results.clone());
        }
        results
    }
}

/// Count `321`-avoiding permutations with a fixed parity orientation.
///
/// Permutations and positions are one-indexed.  A preserving permutation has
/// `pi(i) congruent to i (mod 2)`; a reversing permutation has the opposite
/// parity at every position.  The empty permutation has count one.
///
/// # Dynamic-programming invariant
///
/// A state `(w, r, m)` fills `r` consecutive target positions, beginning with
/// an odd position, using `r` remaining values.  The bits of `w`, from least to
/// most significant, are the parities of those values in increasing order.
/// The `m` smallest values are precisely the values still allowed to occur to
/// the right of the largest value.  In a `321`-avoider, the entries to the
/// right of the largest value must be increasing.  Hence a transition chooses
/// an increasing suffix from those `m` values whose parity word alternates as
/// the target positions do, removes that suffix and the maximum, and recurses
/// on the remaining prefix.  If the chosen suffix has least value of rank
/// `j + 1`, then exactly the `j` smaller values remain legal in the next state.
///
/// The method never constructs a permutation.  With `n <= 63`, there are at
/// most `O(n^2 2^n)` fill states and `O(n 2^n)` alternating-subsequence states;
/// the deliberately loose worst-case time bound is `O(n^2 3^n)` and memory is
/// `O(n^2 2^n)`.  Memoization and parity aggregation make the intended range
/// substantially smaller in practice.
///
/// # Panics
///
/// Panics when `n > 63`, the limit of the compact parity-word representation.
pub fn count_parity_oriented_321_avoiding(n: u8, orientation: ParityOrientation) -> BigUint {
    assert!(n <= 63, "parity-word DP supports n <= 63");
    let mut counter = Counter::default();
    counter.count_321(n, orientation)
}

impl Counter {
    fn count_321(&mut self, n: u8, orientation: ParityOrientation) -> BigUint {
        if n == 0 {
            return BigUint::one();
        }
        let mut parity_word = alternating_word(n);
        if orientation == ParityOrientation::Reversing {
            parity_word ^= low_mask(n);
        }
        self.fill(parity_word, n, n - 1)
    }
}

/// Count parity-preserving `123`-avoiding permutations (OEIS A355089).
///
/// Complementing every value sends `123`-avoidance to `321`-avoidance.  It
/// preserves value parity for odd `n` and reverses it for even `n`, which gives
/// the two orientations used here.
pub fn count_parity_preserving_123_avoiding(n: u8) -> BigUint {
    assert!(n <= 63, "parity-word DP supports n <= 63");
    let orientation = if n % 2 == 1 {
        ParityOrientation::Preserving
    } else {
        ParityOrientation::Reversing
    };
    count_parity_oriented_321_avoiding(n, orientation)
}

/// Return A355089 from index `0` through `max_n`, inclusive.
///
/// One memo table is shared across indices.  The returned entries are exact
/// [`BigUint`] values.
pub fn parity_preserving_123_avoiding_counts(max_n: u8) -> Vec<BigUint> {
    assert!(max_n <= 63, "parity-word DP supports n <= 63");
    let mut counter = Counter::default();
    (0..=max_n)
        .map(|n| {
            let orientation = if n % 2 == 1 {
                ParityOrientation::Preserving
            } else {
                ParityOrientation::Reversing
            };
            counter.count_321(n, orientation)
        })
        .collect()
}

fn alternating_word(len: u8) -> u64 {
    let mut word = 0u64;
    for bit in (0..len).step_by(2) {
        word |= 1u64 << bit;
    }
    word
}

fn low_mask(len: u8) -> u64 {
    if len == 0 {
        0
    } else {
        (1u64 << len) - 1
    }
}

fn low_bits(word: u64, len: u8) -> u64 {
    word & low_mask(len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::permutation::{all_permutations, contains_pattern};

    fn brute_a355089(n: u8) -> usize {
        all_permutations(n)
            .into_iter()
            .filter(|permutation| {
                permutation
                    .iter()
                    .enumerate()
                    .all(|(i, &value)| value as usize % 2 == (i + 1) % 2)
                    && !contains_pattern(permutation, &[1, 2, 3])
            })
            .count()
    }

    #[test]
    fn a355089_matches_independent_brute_force_through_nine() {
        for n in 0..=9 {
            assert_eq!(
                count_parity_preserving_123_avoiding(n),
                BigUint::from(brute_a355089(n)),
                "n={n}"
            );
        }
    }

    #[test]
    fn a355089_matches_all_published_terms() {
        let published = [
            1u64, 1, 1, 1, 3, 3, 10, 11, 37, 44, 146, 185, 603, 808, 2576, 3635, 11294, 16736,
            50545, 78466, 230012, 373203, 1061236, 1795611, 4953447, 8721086, 23350320, 42691298,
            111013825, 210379132, 531720722,
        ];
        let actual = parity_preserving_123_avoiding_counts(30);
        assert_eq!(
            actual,
            published.into_iter().map(BigUint::from).collect::<Vec<_>>()
        );
    }

    #[test]
    fn preserving_321_matches_a354208_prefix() {
        let published = [
            1u64, 1, 1, 1, 2, 3, 6, 11, 22, 44, 89, 185, 382, 808, 1702, 3635, 7779, 16736, 36229,
            78466, 171238,
        ];
        for (n, expected) in published.into_iter().enumerate() {
            assert_eq!(
                count_parity_oriented_321_avoiding(n as u8, ParityOrientation::Preserving),
                BigUint::from(expected),
                "n={n}"
            );
        }
    }
}
