//! Exact counting for parity-oriented permutations avoiding `321` or `123`.
//!
//! The implementation memoizes boundary states and never constructs an output
//! permutation. It supports exact counts and selected additive-statistic
//! distributions. See [`count_parity_oriented_321_avoiding`] for the state
//! invariant and recurrence.

use num_bigint::BigUint;
use num_traits::{One, Zero};
use std::collections::HashMap;
use std::marker::PhantomData;

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

trait DpWeight {
    type Value: Clone;

    fn zero() -> Self::Value;
    fn one() -> Self::Value;
    fn forced_increasing(len: u8) -> Self::Value;
    fn transition_shift(len: u8, suffix_len: u8) -> usize;
    fn add_transition(
        total: &mut Self::Value,
        child: &Self::Value,
        multiplicity: u64,
        shift: usize,
    );
}

struct CountWeight;
struct DescentWeight;
struct MajorIndexWeight;
struct LeftToRightMaximumWeight;

impl DpWeight for CountWeight {
    type Value = BigUint;

    fn zero() -> Self::Value {
        BigUint::zero()
    }

    fn one() -> Self::Value {
        BigUint::one()
    }

    fn forced_increasing(_len: u8) -> Self::Value {
        BigUint::one()
    }

    fn transition_shift(_len: u8, _suffix_len: u8) -> usize {
        0
    }

    fn add_transition(
        total: &mut Self::Value,
        child: &Self::Value,
        multiplicity: u64,
        shift: usize,
    ) {
        debug_assert_eq!(shift, 0);
        *total += child * multiplicity;
    }
}

fn monomial(degree: usize) -> Vec<BigUint> {
    let mut polynomial = vec![BigUint::zero(); degree + 1];
    polynomial[degree] = BigUint::one();
    polynomial
}

fn add_shifted_multiple(
    total: &mut Vec<BigUint>,
    child: &[BigUint],
    multiplicity: u64,
    shift: usize,
) {
    if child.is_empty() {
        return;
    }
    total.resize(total.len().max(child.len() + shift), BigUint::zero());
    for (degree, coefficient) in child.iter().enumerate() {
        total[degree + shift] += coefficient * multiplicity;
    }
}

macro_rules! impl_polynomial_weight {
    ($weight:ty, $forced:expr, $shift:expr) => {
        impl DpWeight for $weight {
            type Value = Vec<BigUint>;

            fn zero() -> Self::Value {
                Vec::new()
            }

            fn one() -> Self::Value {
                monomial(0)
            }

            fn forced_increasing(len: u8) -> Self::Value {
                monomial($forced(len))
            }

            fn transition_shift(len: u8, suffix_len: u8) -> usize {
                $shift(len, suffix_len)
            }

            fn add_transition(
                total: &mut Self::Value,
                child: &Self::Value,
                multiplicity: u64,
                shift: usize,
            ) {
                add_shifted_multiple(total, child, multiplicity, shift);
            }
        }
    };
}

impl_polynomial_weight!(DescentWeight, |_len: u8| 0, |_len: u8, suffix_len: u8| {
    usize::from(suffix_len > 0)
});
impl_polynomial_weight!(MajorIndexWeight, |_len: u8| 0, |len: u8, suffix_len: u8| {
    if suffix_len == 0 {
        0
    } else {
        usize::from(len - suffix_len)
    }
});
impl_polynomial_weight!(
    LeftToRightMaximumWeight,
    |len: u8| usize::from(len),
    |_len: u8, _suffix_len: u8| 1
);

struct Counter<W: DpWeight> {
    fill_cache: HashMap<FillState, W::Value>,
    subsequence_cache: HashMap<SubsequenceState, Vec<(SubsequenceResult, u64)>>,
    weight: PhantomData<W>,
}

impl<W: DpWeight> Default for Counter<W> {
    fn default() -> Self {
        Self {
            fill_cache: HashMap::new(),
            subsequence_cache: HashMap::new(),
            weight: PhantomData,
        }
    }
}

impl<W: DpWeight> Counter<W> {
    fn fill(&mut self, parity_word: u64, len: u8, legal_low: u8) -> W::Value {
        if len == 0 {
            return W::one();
        }
        debug_assert!(legal_low < len);
        debug_assert!(parity_word < (1u64 << len));

        // If no value is legal to the right of the maximum, every remaining
        // value is forced into increasing order.  The same is true for one
        // remaining value.  The target positions begin odd and alternate.
        if legal_low == 0 || len == 1 {
            return if parity_word == alternating_word(len) {
                W::forced_increasing(len)
            } else {
                W::zero()
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
        let mut total = W::zero();

        // An empty suffix is possible exactly when its required length is even.
        if required_suffix_parity == 0 {
            let child = self.fill(parity_word ^ top_bit, len - 1, legal_low.min(len - 2));
            W::add_transition(&mut total, &child, 1, W::transition_shift(len, 0));
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
                let suffix_len = 1 + choice.selected_len;
                W::add_transition(
                    &mut total,
                    &contribution,
                    multiplicity,
                    W::transition_shift(len, suffix_len),
                );
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
    let mut counter = Counter::<CountWeight>::default();
    counter.weighted_321(n, orientation)
}

impl<W: DpWeight> Counter<W> {
    fn weighted_321(&mut self, n: u8, orientation: ParityOrientation) -> W::Value {
        if n == 0 {
            return W::one();
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
    count_parity_oriented_321_avoiding(n, a355089_orientation(n))
}

/// Return A355089 from index `0` through `max_n`, inclusive.
///
/// One memo table is shared across indices.  The returned entries are exact
/// [`BigUint`] values.
pub fn parity_preserving_123_avoiding_counts(max_n: u8) -> Vec<BigUint> {
    assert!(max_n <= 63, "parity-word DP supports n <= 63");
    let mut counter = Counter::<CountWeight>::default();
    (0..=max_n)
        .map(|n| counter.weighted_321(n, a355089_orientation(n)))
        .collect()
}

/// A statistic whose distribution is computed on parity-preserving
/// `123`-avoiding permutations.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Parity123Statistic {
    /// Number of descents.
    Descents,
    /// Number of ascents.
    Ascents,
    /// Sum of descent positions.
    MajorIndex,
    /// Sum of ascent positions.
    AscentMajorIndex,
    /// Number of left-to-right minima.
    LeftToRightMinima,
}

/// Return exact statistic polynomials through `max_n` for A355089.
///
/// Row `n` is in ascending degree order: its coefficient of `t^k` counts
/// parity-preserving `123`-avoiders of length `n` with statistic value `k`.
/// Rows retain their natural full degree (`n-1` for descents or ascents,
/// `binomial(n,2)` for either position sum, and `n` for left-to-right minima),
/// including structural zero coefficients.
///
/// The same maximum-and-gap boundary states as the count-only DP suffice.
/// A nonempty increasing gap contributes one descent at position
/// `len-suffix_len`; each peeled maximum contributes one left-to-right maximum.
/// Value complement turns descents into ascents, major index into complementary
/// major index, and left-to-right maxima into left-to-right minima.
/// At even lengths the left-to-right-minimum rows are palindromic: reverse
/// complement preserves the class and exchanges `k` left-to-right minima with
/// `n-k` of them. See the project documentation for the disjoint-record proof.
///
/// # Panics
///
/// Panics when `max_n > 63`.
pub fn parity_preserving_123_statistic_polynomials(
    max_n: u8,
    statistic: Parity123Statistic,
) -> Vec<Vec<BigUint>> {
    assert!(max_n <= 63, "parity-word DP supports n <= 63");
    match statistic {
        Parity123Statistic::Descents => {
            let mut counter = Counter::<DescentWeight>::default();
            (0..=max_n)
                .map(|n| {
                    let polynomial = counter.weighted_321(n, a355089_orientation(n));
                    reverse_polynomial(polynomial, n.saturating_sub(1) as usize)
                })
                .collect()
        }
        Parity123Statistic::Ascents => {
            let mut counter = Counter::<DescentWeight>::default();
            (0..=max_n)
                .map(|n| {
                    let mut polynomial = counter.weighted_321(n, a355089_orientation(n));
                    polynomial.resize(usize::from(n.max(1)), BigUint::zero());
                    polynomial
                })
                .collect()
        }
        Parity123Statistic::MajorIndex => {
            let mut counter = Counter::<MajorIndexWeight>::default();
            (0..=max_n)
                .map(|n| {
                    let polynomial = counter.weighted_321(n, a355089_orientation(n));
                    let degree = usize::from(n) * usize::from(n.saturating_sub(1)) / 2;
                    reverse_polynomial(polynomial, degree)
                })
                .collect()
        }
        Parity123Statistic::AscentMajorIndex => {
            let mut counter = Counter::<MajorIndexWeight>::default();
            (0..=max_n)
                .map(|n| {
                    let mut polynomial = counter.weighted_321(n, a355089_orientation(n));
                    let degree = usize::from(n) * usize::from(n.saturating_sub(1)) / 2;
                    polynomial.resize(degree + 1, BigUint::zero());
                    polynomial
                })
                .collect()
        }
        Parity123Statistic::LeftToRightMinima => {
            let mut counter = Counter::<LeftToRightMaximumWeight>::default();
            (0..=max_n)
                .map(|n| {
                    let mut polynomial = counter.weighted_321(n, a355089_orientation(n));
                    polynomial.resize(usize::from(n) + 1, BigUint::zero());
                    polynomial
                })
                .collect()
        }
    }
}

fn a355089_orientation(n: u8) -> ParityOrientation {
    if n % 2 == 1 {
        ParityOrientation::Preserving
    } else {
        ParityOrientation::Reversing
    }
}

fn reverse_polynomial(polynomial: Vec<BigUint>, total_degree: usize) -> Vec<BigUint> {
    let mut reversed = vec![BigUint::zero(); total_degree + 1];
    for (degree, coefficient) in polynomial.into_iter().enumerate() {
        debug_assert!(degree <= total_degree);
        reversed[total_degree - degree] = coefficient;
    }
    reversed
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
    use crate::statistics::{compute, Stat};

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

    fn brute_statistic_polynomial(n: u8, statistic: Stat, degree: usize) -> Vec<BigUint> {
        let mut polynomial = vec![BigUint::zero(); degree + 1];
        for permutation in all_permutations(n).into_iter().filter(|permutation| {
            permutation
                .iter()
                .enumerate()
                .all(|(i, &value)| value as usize % 2 == (i + 1) % 2)
                && !contains_pattern(permutation, &[1, 2, 3])
        }) {
            polynomial[compute(&permutation, statistic)] += 1u8;
        }
        polynomial
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

    #[test]
    fn statistic_refinements_match_independent_brute_force() {
        for (statistic, brute_statistic) in [
            (Parity123Statistic::Descents, Stat::Des),
            (Parity123Statistic::Ascents, Stat::Asc),
            (Parity123Statistic::MajorIndex, Stat::Maj),
            (Parity123Statistic::LeftToRightMinima, Stat::Lrmin),
        ] {
            let actual = parity_preserving_123_statistic_polynomials(8, statistic);
            for n in 0..=8u8 {
                let degree = match statistic {
                    Parity123Statistic::Descents | Parity123Statistic::Ascents => {
                        usize::from(n.saturating_sub(1))
                    }
                    Parity123Statistic::MajorIndex | Parity123Statistic::AscentMajorIndex => {
                        usize::from(n) * usize::from(n.saturating_sub(1)) / 2
                    }
                    Parity123Statistic::LeftToRightMinima => usize::from(n),
                };
                assert_eq!(
                    actual[n as usize],
                    brute_statistic_polynomial(n, brute_statistic, degree),
                    "statistic={statistic:?}, n={n}"
                );
            }
        }

        let actual =
            parity_preserving_123_statistic_polynomials(8, Parity123Statistic::AscentMajorIndex);
        for n in 0..=8u8 {
            let degree = usize::from(n) * usize::from(n.saturating_sub(1)) / 2;
            let mut expected = vec![BigUint::zero(); degree + 1];
            for permutation in all_permutations(n).into_iter().filter(|permutation| {
                permutation
                    .iter()
                    .enumerate()
                    .all(|(i, &value)| value as usize % 2 == (i + 1) % 2)
                    && !contains_pattern(permutation, &[1, 2, 3])
            }) {
                let ascent_major = permutation
                    .windows(2)
                    .enumerate()
                    .filter(|(_, adjacent)| adjacent[0] < adjacent[1])
                    .map(|(i, _)| i + 1)
                    .sum::<usize>();
                expected[ascent_major] += 1u8;
            }
            assert_eq!(actual[n as usize], expected, "ascent major, n={n}");
        }
    }

    #[test]
    fn statistic_refinements_specialize_to_a355089() {
        let counts = parity_preserving_123_avoiding_counts(20);
        for statistic in [
            Parity123Statistic::Descents,
            Parity123Statistic::Ascents,
            Parity123Statistic::MajorIndex,
            Parity123Statistic::AscentMajorIndex,
            Parity123Statistic::LeftToRightMinima,
        ] {
            let rows = parity_preserving_123_statistic_polynomials(20, statistic);
            for (n, row) in rows.iter().enumerate() {
                assert_eq!(
                    row.iter().sum::<BigUint>(),
                    counts[n],
                    "{statistic:?}, n={n}"
                );
            }
        }
    }

    #[test]
    fn even_left_to_right_minimum_rows_are_palindromic() {
        let rows =
            parity_preserving_123_statistic_polynomials(20, Parity123Statistic::LeftToRightMinima);
        for n in (0..=20usize).step_by(2) {
            assert!(
                rows[n].iter().eq(rows[n].iter().rev()),
                "left-to-right minima, n={n}"
            );
        }
    }
}
