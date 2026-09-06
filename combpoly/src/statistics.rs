//! Permutation and word statistics (des, inv, maj, exc, etc.).
//!
//! Includes both scalar statistics ([`Stat`] / [`compute`]) and
//! set-valued statistics ([`SetStat`] / [`compute_set`]).

use std::collections::BTreeSet;
use std::fmt;

use clap::ValueEnum;

/// A combinatorial statistic on words or permutations.
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Stat {
    /// Number of descents
    Des,
    /// Number of ascents
    Asc,
    /// Number of adjacent ties
    Tie,
    /// Number of excedances (w_i > i)
    Exc,
    /// Number of strict peaks: w_{i-1} < w_i > w_{i+1}
    Peak,
    /// Number of weak-left peaks: w_{i-1} <= w_i > w_{i+1}
    WeakPeak,
    /// Number of inversions
    Inv,
    /// Major index (sum of descent positions)
    Maj,
    /// Comajor index (sum of n-i at descent positions)
    Comaj,
    /// Number of fixed points (w_i = i)
    Fix,
    /// Number of cycles (permutations only)
    Cyc,
    /// Number of valleys
    Valley,
    /// Number of left-to-right minima
    Lrmin,
    /// Number of left-to-right maxima
    Lrmax,
    /// Number of right-to-left minima
    Rlmin,
    /// Number of right-to-left maxima
    Rlmax,
    /// Charge (Lascoux-Schutzenberger, permutations only)
    Charge,
    /// Cocharge = C(n,2) - charge (permutations only)
    Cocharge,
    /// Number of coinversions
    Coinv,
    /// Length of longest increasing subsequence
    Lis,
    /// Length of longest decreasing subsequence
    Lds,
    /// Number of long swaps: |{i ∈ [n-1] : σ⁻¹(i) < σ⁻¹(i+1) - 1}|
    Swaps,
}

impl fmt::Display for Stat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Stat::Des => write!(f, "des"),
            Stat::Asc => write!(f, "asc"),
            Stat::Tie => write!(f, "tie"),
            Stat::Exc => write!(f, "exc"),
            Stat::Peak => write!(f, "peak"),
            Stat::WeakPeak => write!(f, "weak_peak"),
            Stat::Inv => write!(f, "inv"),
            Stat::Maj => write!(f, "maj"),
            Stat::Comaj => write!(f, "comaj"),
            Stat::Fix => write!(f, "fix"),
            Stat::Cyc => write!(f, "cyc"),
            Stat::Valley => write!(f, "valley"),
            Stat::Lrmin => write!(f, "lrmin"),
            Stat::Lrmax => write!(f, "lrmax"),
            Stat::Rlmin => write!(f, "rlmin"),
            Stat::Rlmax => write!(f, "rlmax"),
            Stat::Charge => write!(f, "charge"),
            Stat::Cocharge => write!(f, "cocharge"),
            Stat::Coinv => write!(f, "coinv"),
            Stat::Lis => write!(f, "lis"),
            Stat::Lds => write!(f, "lds"),
            Stat::Swaps => write!(f, "swaps"),
        }
    }
}

/// Evaluate the given statistic on a word or permutation.
pub fn compute(w: &[u8], stat: Stat) -> usize {
    match stat {
        Stat::Des => descents(w),
        Stat::Asc => ascents(w),
        Stat::Tie => ties(w),
        Stat::Exc => excedances(w),
        Stat::Peak => peaks(w),
        Stat::WeakPeak => weak_left_peaks(w),
        Stat::Inv => inversions(w),
        Stat::Maj => major_index(w),
        Stat::Comaj => comajor_index(w),
        Stat::Fix => fixed_points(w),
        Stat::Cyc => cycles(w),
        Stat::Valley => valleys(w),
        Stat::Lrmin => left_to_right_minima(w),
        Stat::Lrmax => left_to_right_maxima(w),
        Stat::Rlmin => right_to_left_minima(w),
        Stat::Rlmax => right_to_left_maxima(w),
        Stat::Charge => charge(w),
        Stat::Cocharge => cocharge(w),
        Stat::Coinv => coinversions(w),
        Stat::Lis => longest_increasing_subseq(w),
        Stat::Lds => longest_decreasing_subseq(w),
        Stat::Swaps => long_swaps(w),
    }
}

// ---------------------------------------------------------------------------
// Set-valued statistics
// ---------------------------------------------------------------------------

/// A set-valued combinatorial statistic on words or permutations.
///
/// Each variant returns a `BTreeSet<usize>` of 1-indexed positions (or values).
#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum SetStat {
    /// Descent set {i : w_i > w_{i+1}}, 1-indexed positions
    DesSet,
    /// Ascent set {i : w_i < w_{i+1}}, 1-indexed positions
    AscSet,
    /// Tie set {i : w_i = w_{i+1}}, 1-indexed positions
    TieSet,
    /// Descent bottom set {w_{i+1} : i in Des(w)}, values at bottom of descents
    DesBottomSet,
    /// Descent top set {w_i : i in Des(w)}, values at top of descents
    DesTopSet,
    /// Excedance set {i : w_i > i}, 1-indexed positions
    ExcSet,
    /// Fixed-point set {i : w_i = i}, 1-indexed positions
    FixSet,
    /// Strict peak set {i in {2,..,n-1} : w_{i-1} < w_i > w_{i+1}}
    PeakSet,
    /// Weak-left peak set {i in {2,..,n-1} : w_{i-1} <= w_i > w_{i+1}}
    WeakPeakSet,
    /// Valley set {i in {2,..,n-1} : w_{i-1} > w_i < w_{i+1}}, 1-indexed positions
    ValleySet,
    /// Left-to-right minima positions, 1-indexed
    LrminSet,
    /// Left-to-right maxima positions, 1-indexed
    LrmaxSet,
    /// Right-to-left minima positions, 1-indexed
    RlminSet,
    /// Right-to-left maxima positions, 1-indexed
    RlmaxSet,
    /// Long swaps set {i ∈ [n-1] : σ⁻¹(i) < σ⁻¹(i+1) - 1}, 1-indexed values
    SwapsSet,
}

impl fmt::Display for SetStat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            SetStat::DesSet => write!(f, "des_set"),
            SetStat::AscSet => write!(f, "asc_set"),
            SetStat::TieSet => write!(f, "tie_set"),
            SetStat::DesBottomSet => write!(f, "des_bottom_set"),
            SetStat::DesTopSet => write!(f, "des_top_set"),
            SetStat::ExcSet => write!(f, "exc_set"),
            SetStat::FixSet => write!(f, "fix_set"),
            SetStat::PeakSet => write!(f, "peak_set"),
            SetStat::WeakPeakSet => write!(f, "weak_peak_set"),
            SetStat::ValleySet => write!(f, "valley_set"),
            SetStat::LrminSet => write!(f, "lrmin_set"),
            SetStat::LrmaxSet => write!(f, "lrmax_set"),
            SetStat::RlminSet => write!(f, "rlmin_set"),
            SetStat::RlmaxSet => write!(f, "rlmax_set"),
            SetStat::SwapsSet => write!(f, "swaps_set"),
        }
    }
}

/// Evaluate the given set-valued statistic on a word or permutation.
///
/// Returns a `BTreeSet<usize>` of 1-indexed positions (or values, for
/// `DesBottomSet` and `DesTopSet`).
pub fn compute_set(w: &[u8], stat: SetStat) -> BTreeSet<usize> {
    match stat {
        SetStat::DesSet => descent_set(w),
        SetStat::AscSet => ascent_set(w),
        SetStat::TieSet => tie_set(w),
        SetStat::DesBottomSet => descent_bottom_set(w),
        SetStat::DesTopSet => descent_top_set(w),
        SetStat::ExcSet => excedance_set(w),
        SetStat::FixSet => fixed_point_set(w),
        SetStat::PeakSet => peak_set(w),
        SetStat::WeakPeakSet => weak_left_peak_set(w),
        SetStat::ValleySet => valley_set(w),
        SetStat::LrminSet => lrmin_set(w),
        SetStat::LrmaxSet => lrmax_set(w),
        SetStat::RlminSet => rlmin_set(w),
        SetStat::RlmaxSet => rlmax_set(w),
        SetStat::SwapsSet => long_swaps_set(w),
    }
}

fn descent_set(w: &[u8]) -> BTreeSet<usize> {
    (1..w.len()).filter(|&i| w[i - 1] > w[i]).collect()
}

fn ascent_set(w: &[u8]) -> BTreeSet<usize> {
    (1..w.len()).filter(|&i| w[i - 1] < w[i]).collect()
}

fn tie_set(w: &[u8]) -> BTreeSet<usize> {
    (1..w.len()).filter(|&i| w[i - 1] == w[i]).collect()
}

fn descent_bottom_set(w: &[u8]) -> BTreeSet<usize> {
    (1..w.len())
        .filter(|&i| w[i - 1] > w[i])
        .map(|i| w[i] as usize)
        .collect()
}

fn descent_top_set(w: &[u8]) -> BTreeSet<usize> {
    (1..w.len())
        .filter(|&i| w[i - 1] > w[i])
        .map(|i| w[i - 1] as usize)
        .collect()
}

fn excedance_set(w: &[u8]) -> BTreeSet<usize> {
    (0..w.len())
        .filter(|&i| w[i] as usize > i + 1)
        .map(|i| i + 1)
        .collect()
}

fn fixed_point_set(w: &[u8]) -> BTreeSet<usize> {
    (0..w.len())
        .filter(|&i| w[i] as usize == i + 1)
        .map(|i| i + 1)
        .collect()
}

fn peak_set(w: &[u8]) -> BTreeSet<usize> {
    if w.len() < 3 {
        return BTreeSet::new();
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .map(|i| i + 1) // 1-indexed
        .collect()
}

fn weak_left_peak_set(w: &[u8]) -> BTreeSet<usize> {
    if w.len() < 3 {
        return BTreeSet::new();
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] <= w[i] && w[i] > w[i + 1])
        .map(|i| i + 1)
        .collect()
}

fn valley_set(w: &[u8]) -> BTreeSet<usize> {
    if w.len() < 3 {
        return BTreeSet::new();
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] > w[i] && w[i] < w[i + 1])
        .map(|i| i + 1) // 1-indexed
        .collect()
}

fn lrmin_set(w: &[u8]) -> BTreeSet<usize> {
    let mut result = BTreeSet::new();
    let mut min_so_far = u8::MAX;
    for (i, &v) in w.iter().enumerate() {
        if v < min_so_far {
            result.insert(i + 1);
            min_so_far = v;
        }
    }
    result
}

fn lrmax_set(w: &[u8]) -> BTreeSet<usize> {
    let mut result = BTreeSet::new();
    let mut max_so_far = 0u8;
    for (i, &v) in w.iter().enumerate() {
        if v > max_so_far {
            result.insert(i + 1);
            max_so_far = v;
        }
    }
    result
}

fn rlmin_set(w: &[u8]) -> BTreeSet<usize> {
    let mut result = BTreeSet::new();
    let mut min_so_far = u8::MAX;
    for (i, &v) in w.iter().enumerate().rev() {
        if v < min_so_far {
            result.insert(i + 1);
            min_so_far = v;
        }
    }
    result
}

fn rlmax_set(w: &[u8]) -> BTreeSet<usize> {
    let mut result = BTreeSet::new();
    let mut max_so_far = 0u8;
    for (i, &v) in w.iter().enumerate().rev() {
        if v > max_so_far {
            result.insert(i + 1);
            max_so_far = v;
        }
    }
    result
}

// --- Long swaps ---

/// Compute the long swaps statistic.
///
/// `swaps(σ) = |{i ∈ [n-1] : σ⁻¹(i) < σ⁻¹(i+1) - 1}|`
///
/// Counts values i such that i appears to the left of i+1 in σ, but not adjacent.
fn long_swaps(w: &[u8]) -> usize {
    let n = w.len();
    if n <= 1 {
        return 0;
    }
    let mut inv = vec![0usize; n + 1];
    for (pos, &val) in w.iter().enumerate() {
        inv[val as usize] = pos;
    }
    let mut count = 0;
    for i in 1..n {
        if inv[i] < inv[i + 1] && inv[i + 1] - inv[i] > 1 {
            count += 1;
        }
    }
    count
}

/// Return the set of values i forming long swaps.
///
/// `{i ∈ [n-1] : σ⁻¹(i) < σ⁻¹(i+1) - 1}`
fn long_swaps_set(w: &[u8]) -> BTreeSet<usize> {
    let n = w.len();
    if n <= 1 {
        return BTreeSet::new();
    }
    let mut inv = vec![0usize; n + 1];
    for (pos, &val) in w.iter().enumerate() {
        inv[val as usize] = pos;
    }
    (1..n)
        .filter(|&i| inv[i] < inv[i + 1] && inv[i + 1] - inv[i] > 1)
        .collect()
}

/// Descent set as a bitmask: bit j is set iff position j+1 is a descent (1-indexed).
///
/// That is, bit j (0-indexed) is set iff `w[j] > w[j+1]`,
/// corresponding to a descent at 1-indexed position j+1.
///
/// This is efficient for grouping permutations by descent set.
///
/// # Panics
///
/// Panics if `w` has a descent beyond the 64 positions representable by a
/// `u64`. Use [`checked_descent_set_bitmask`] when this is not known in
/// advance.
pub fn descent_set_bitmask(w: &[u8]) -> u64 {
    checked_descent_set_bitmask(w)
        .expect("descent at position 65 or later cannot be represented by a u64 bitmask")
}

/// Descent set as a `u64` bitmask, returning `None` if it does not fit.
///
/// A word longer than 65 entries can still be represented when all of its
/// descents occur among the first 64 positions.
pub fn checked_descent_set_bitmask(w: &[u8]) -> Option<u64> {
    let mut s = 0u64;
    for i in 0..w.len().saturating_sub(1) {
        if w[i] > w[i + 1] {
            if i >= u64::BITS as usize {
                return None;
            }
            s |= 1 << i;
        }
    }
    Some(s)
}

// --- Basic statistics ---

fn descents(w: &[u8]) -> usize {
    (1..w.len()).filter(|&i| w[i - 1] > w[i]).count()
}

fn ascents(w: &[u8]) -> usize {
    (1..w.len()).filter(|&i| w[i - 1] < w[i]).count()
}

fn ties(w: &[u8]) -> usize {
    (1..w.len()).filter(|&i| w[i - 1] == w[i]).count()
}

fn excedances(w: &[u8]) -> usize {
    (0..w.len()).filter(|&i| w[i] as usize > i + 1).count()
}

/// Number of peaks in a word or permutation.
///
/// A peak is an index `i` with `1 < i < n` such that
/// `w[i - 1] < w[i] > w[i + 1]` (using 0-based indexing internally).
pub fn peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] < w[i] && w[i] > w[i + 1])
        .count()
}

/// Number of weak-left peaks in a word.
///
/// A weak-left peak is an index `i` with `1 < i < n` such that
/// `w[i - 1] <= w[i] > w[i + 1]`.  Unlike a strict peak, it is determined by
/// the descent set.
pub fn weak_left_peaks(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] <= w[i] && w[i] > w[i + 1])
        .count()
}

fn valleys(w: &[u8]) -> usize {
    if w.len() < 3 {
        return 0;
    }
    (1..w.len() - 1)
        .filter(|&i| w[i - 1] > w[i] && w[i] < w[i + 1])
        .count()
}

fn inversions(w: &[u8]) -> usize {
    let mut count = 0;
    for i in 0..w.len() {
        for j in i + 1..w.len() {
            if w[i] > w[j] {
                count += 1;
            }
        }
    }
    count
}

fn coinversions(w: &[u8]) -> usize {
    let n = w.len();
    n * (n - 1) / 2 - inversions(w)
}

fn major_index(w: &[u8]) -> usize {
    (0..w.len().saturating_sub(1))
        .filter(|&i| w[i] > w[i + 1])
        .map(|i| i + 1)
        .sum()
}

fn comajor_index(w: &[u8]) -> usize {
    let n = w.len();
    (0..n.saturating_sub(1))
        .filter(|&i| w[i] > w[i + 1])
        .map(|i| n - 1 - i)
        .sum()
}

fn fixed_points(w: &[u8]) -> usize {
    (0..w.len()).filter(|&i| w[i] as usize == i + 1).count()
}

// --- Cycle statistics (permutations only) ---

fn cycles(w: &[u8]) -> usize {
    let n = w.len();
    let mut visited = vec![false; n];
    let mut count = 0;
    for i in 0..n {
        if !visited[i] {
            count += 1;
            let mut j = i;
            loop {
                if j >= n || visited[j] {
                    break;
                }
                visited[j] = true;
                j = (w[j] as usize).wrapping_sub(1);
            }
        }
    }
    count
}

// --- Left-to-right / right-to-left extrema ---

fn left_to_right_minima(w: &[u8]) -> usize {
    let mut count = 0;
    let mut min_so_far = u8::MAX;
    for &v in w {
        if v < min_so_far {
            count += 1;
            min_so_far = v;
        }
    }
    count
}

fn left_to_right_maxima(w: &[u8]) -> usize {
    let mut count = 0;
    let mut max_so_far = 0u8;
    for &v in w {
        if v > max_so_far {
            count += 1;
            max_so_far = v;
        }
    }
    count
}

fn right_to_left_minima(w: &[u8]) -> usize {
    let mut count = 0;
    let mut min_so_far = u8::MAX;
    for &v in w.iter().rev() {
        if v < min_so_far {
            count += 1;
            min_so_far = v;
        }
    }
    count
}

fn right_to_left_maxima(w: &[u8]) -> usize {
    let mut count = 0;
    let mut max_so_far = 0u8;
    for &v in w.iter().rev() {
        if v > max_so_far {
            count += 1;
            max_so_far = v;
        }
    }
    count
}

// --- Charge/cocharge (Lascoux-Schutzenberger, permutations only) ---
// charge(w) = maj(reverse(inverse(w)))

fn charge(w: &[u8]) -> usize {
    let n = w.len();
    if n == 0 {
        return 0;
    }
    let mut inv_perm = vec![0u8; n];
    for (i, &w_i) in w.iter().enumerate() {
        let v = w_i as usize;
        if v == 0 || v > n {
            return 0;
        }
        inv_perm[v - 1] = (i + 1) as u8;
    }
    inv_perm.reverse();
    major_index(&inv_perm)
}

fn cocharge(w: &[u8]) -> usize {
    let n = w.len();
    n * (n - 1) / 2 - charge(w)
}

// --- Subsequence statistics ---

fn longest_increasing_subseq(w: &[u8]) -> usize {
    if w.is_empty() {
        return 0;
    }
    let n = w.len();
    let mut dp = vec![1usize; n];
    for i in 1..n {
        for j in 0..i {
            if w[j] < w[i] {
                dp[i] = dp[i].max(dp[j] + 1);
            }
        }
    }
    *dp.iter().max().unwrap()
}

fn longest_decreasing_subseq(w: &[u8]) -> usize {
    if w.is_empty() {
        return 0;
    }
    let n = w.len();
    let mut dp = vec![1usize; n];
    for i in 1..n {
        for j in 0..i {
            if w[j] > w[i] {
                dp[i] = dp[i].max(dp[j] + 1);
            }
        }
    }
    *dp.iter().max().unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_descents() {
        assert_eq!(compute(&[1, 2, 3], Stat::Des), 0);
        assert_eq!(compute(&[3, 2, 1], Stat::Des), 2);
        assert_eq!(compute(&[3, 1, 4, 2], Stat::Des), 2);
    }

    #[test]
    fn test_ascents() {
        assert_eq!(compute(&[1, 2, 3], Stat::Asc), 2);
        assert_eq!(compute(&[3, 2, 1], Stat::Asc), 0);
    }

    #[test]
    fn test_ties() {
        assert_eq!(compute(&[1, 1, 2, 2, 1], Stat::Tie), 2);
        assert_eq!(compute(&[1, 2, 1], Stat::Tie), 0);
        assert_eq!(compute_set(&[1, 1, 2, 2, 1], SetStat::TieSet), set(&[1, 3]));
    }

    #[test]
    fn test_inversions_and_coinversions() {
        assert_eq!(compute(&[1, 2, 3], Stat::Inv), 0);
        assert_eq!(compute(&[3, 2, 1], Stat::Inv), 3);
        assert_eq!(compute(&[2, 1, 3], Stat::Inv), 1);
        // coinv = C(n,2) - inv
        assert_eq!(compute(&[1, 2, 3], Stat::Coinv), 3);
        assert_eq!(compute(&[3, 2, 1], Stat::Coinv), 0);
    }

    #[test]
    fn test_major_index() {
        assert_eq!(compute(&[1, 2, 3], Stat::Maj), 0);
        assert_eq!(compute(&[3, 2, 1], Stat::Maj), 3); // des at 1-indexed pos 1,2 -> 1+2
                                                       // [2,3,1]: descent at 1-indexed position 2 (3>1), so maj = 2
        assert_eq!(compute(&[2, 3, 1], Stat::Maj), 2);
    }

    #[test]
    fn test_comajor_index() {
        // comaj(321) = (3-1-0) + (3-1-1) = 2+1 = 3
        assert_eq!(compute(&[3, 2, 1], Stat::Comaj), 3);
        assert_eq!(compute(&[1, 2, 3], Stat::Comaj), 0);
    }

    #[test]
    fn test_excedances() {
        assert_eq!(compute(&[1, 2, 3], Stat::Exc), 0);
        assert_eq!(compute(&[2, 3, 1], Stat::Exc), 2);
    }

    #[test]
    fn test_peaks_and_valleys() {
        assert_eq!(compute(&[1, 3, 2], Stat::Peak), 1);
        assert_eq!(compute(&[1, 2, 3], Stat::Peak), 0);
        assert_eq!(compute(&[3, 1, 2], Stat::Valley), 1);
        assert_eq!(compute(&[1, 2, 3], Stat::Valley), 0);

        // Ties distinguish strict peaks from weak-left peaks.
        assert_eq!(compute(&[1, 2, 2, 1], Stat::Peak), 0);
        assert_eq!(compute(&[1, 2, 2, 1], Stat::WeakPeak), 1);
        assert_eq!(compute(&[1, 1, 2, 1], Stat::Peak), 1);
        assert_eq!(compute(&[1, 1, 2, 1], Stat::WeakPeak), 1);
    }

    #[test]
    fn test_fixed_points() {
        assert_eq!(compute(&[1, 2, 3], Stat::Fix), 3);
        assert_eq!(compute(&[2, 1, 3], Stat::Fix), 1);
        assert_eq!(compute(&[2, 3, 1], Stat::Fix), 0);
    }

    #[test]
    fn test_cycles() {
        assert_eq!(compute(&[1, 2, 3], Stat::Cyc), 3); // 3 fixed points
        assert_eq!(compute(&[2, 3, 1], Stat::Cyc), 1); // one 3-cycle
        assert_eq!(compute(&[2, 1, 3], Stat::Cyc), 2); // one 2-cycle + one fix
    }

    #[test]
    fn test_lr_rl_extrema() {
        assert_eq!(compute(&[3, 2, 1], Stat::Lrmin), 3);
        assert_eq!(compute(&[1, 2, 3], Stat::Lrmin), 1);
        assert_eq!(compute(&[3, 1, 2], Stat::Lrmin), 2);

        assert_eq!(compute(&[1, 2, 3], Stat::Lrmax), 3);
        assert_eq!(compute(&[3, 2, 1], Stat::Lrmax), 1);

        assert_eq!(compute(&[3, 2, 1], Stat::Rlmin), 1);
        assert_eq!(compute(&[1, 2, 3], Stat::Rlmin), 3);

        assert_eq!(compute(&[1, 2, 3], Stat::Rlmax), 1);
        assert_eq!(compute(&[3, 2, 1], Stat::Rlmax), 3);
    }

    #[test]
    fn test_charge_cocharge() {
        // charge(123) = maj(reverse(inverse(123))) = maj(321) = 3
        assert_eq!(compute(&[1, 2, 3], Stat::Charge), 3);
        // charge(321) = maj(reverse(inverse(321))) = maj(123) = 0
        assert_eq!(compute(&[3, 2, 1], Stat::Charge), 0);
        // cocharge = C(n,2) - charge
        assert_eq!(compute(&[1, 2, 3], Stat::Cocharge), 0);
        assert_eq!(compute(&[3, 2, 1], Stat::Cocharge), 3);
        // charge(3,5,1,4,2) = 6
        assert_eq!(compute(&[3, 5, 1, 4, 2], Stat::Charge), 6);
    }

    #[test]
    fn test_lis_lds() {
        assert_eq!(compute(&[1, 2, 3], Stat::Lis), 3);
        assert_eq!(compute(&[3, 2, 1], Stat::Lis), 1);
        assert_eq!(compute(&[2, 1, 4, 3], Stat::Lis), 2);

        assert_eq!(compute(&[1, 2, 3], Stat::Lds), 1);
        assert_eq!(compute(&[3, 2, 1], Stat::Lds), 3);
        assert_eq!(compute(&[2, 1, 4, 3], Stat::Lds), 2);
    }

    // --- Set-valued statistic tests ---

    fn set(v: &[usize]) -> BTreeSet<usize> {
        v.iter().copied().collect()
    }

    #[test]
    fn test_descent_set() {
        // [3,1,4,2]: descents at positions 1 (3>1) and 3 (4>2)
        assert_eq!(compute_set(&[3, 1, 4, 2], SetStat::DesSet), set(&[1, 3]));
        assert_eq!(compute_set(&[1, 2, 3], SetStat::DesSet), set(&[]));
        assert_eq!(compute_set(&[3, 2, 1], SetStat::DesSet), set(&[1, 2]));
    }

    #[test]
    fn test_ascent_set() {
        assert_eq!(compute_set(&[1, 2, 3], SetStat::AscSet), set(&[1, 2]));
        assert_eq!(compute_set(&[3, 2, 1], SetStat::AscSet), set(&[]));
    }

    #[test]
    fn test_descent_bottom_set() {
        // [3,1,4,2]: descent at pos 1 gives bottom w[1]=1, descent at pos 3 gives bottom w[3]=2
        assert_eq!(
            compute_set(&[3, 1, 4, 2], SetStat::DesBottomSet),
            set(&[1, 2])
        );
        // [4,2,3,1]: descents at 1 (4>2) and 3 (3>1), bottoms are 2 and 1
        assert_eq!(
            compute_set(&[4, 2, 3, 1], SetStat::DesBottomSet),
            set(&[1, 2])
        );
    }

    #[test]
    fn test_descent_top_set() {
        // [3,1,4,2]: descent at pos 1 gives top w[0]=3, descent at pos 3 gives top w[2]=4
        assert_eq!(compute_set(&[3, 1, 4, 2], SetStat::DesTopSet), set(&[3, 4]));
    }

    #[test]
    fn test_excedance_set() {
        // [2,3,1]: pos 1 (2>1), pos 2 (3>2), pos 3 (1<3)
        assert_eq!(compute_set(&[2, 3, 1], SetStat::ExcSet), set(&[1, 2]));
        assert_eq!(compute_set(&[1, 2, 3], SetStat::ExcSet), set(&[]));
    }

    #[test]
    fn test_fixed_point_set() {
        assert_eq!(compute_set(&[1, 2, 3], SetStat::FixSet), set(&[1, 2, 3]));
        assert_eq!(compute_set(&[2, 1, 3], SetStat::FixSet), set(&[3]));
        assert_eq!(compute_set(&[2, 3, 1], SetStat::FixSet), set(&[]));
    }

    #[test]
    fn test_peak_set() {
        // [1,3,2]: peak at position 2 (1 < 3 > 2)
        assert_eq!(compute_set(&[1, 3, 2], SetStat::PeakSet), set(&[2]));
        assert_eq!(compute_set(&[1, 2, 3], SetStat::PeakSet), set(&[]));
        // [2,4,1,3,5]: peak at position 2 (2<4>1)
        assert_eq!(compute_set(&[2, 4, 1, 3, 5], SetStat::PeakSet), set(&[2]));
        assert_eq!(compute_set(&[1, 2, 2, 1], SetStat::PeakSet), set(&[]));
        assert_eq!(compute_set(&[1, 2, 2, 1], SetStat::WeakPeakSet), set(&[3]));
    }

    #[test]
    fn test_valley_set() {
        // [3,1,2]: valley at position 2 (3 > 1 < 2)
        assert_eq!(compute_set(&[3, 1, 2], SetStat::ValleySet), set(&[2]));
        assert_eq!(compute_set(&[1, 2, 3], SetStat::ValleySet), set(&[]));
    }

    #[test]
    fn test_lrmin_lrmax_set() {
        // [3,1,4,2]: LR-min at pos 1 (3), pos 2 (1)
        assert_eq!(compute_set(&[3, 1, 4, 2], SetStat::LrminSet), set(&[1, 2]));
        // [1,2,3]: LR-max at all positions
        assert_eq!(compute_set(&[1, 2, 3], SetStat::LrmaxSet), set(&[1, 2, 3]));
    }

    #[test]
    fn test_rlmin_rlmax_set() {
        // [3,2,1]: scanning right-to-left: val 1 (pos 3) is min, then 2>1, 3>1 — only pos 3
        assert_eq!(compute_set(&[3, 2, 1], SetStat::RlminSet), set(&[3]));
        // [1,2,3]: scanning right-to-left: 3(pos3), 2(pos2) new min, 1(pos1) new min
        assert_eq!(compute_set(&[1, 2, 3], SetStat::RlminSet), set(&[1, 2, 3]));
        // [1,2,3]: RL-max at pos 3 only (scanning right-to-left: 3 is first and biggest)
        assert_eq!(compute_set(&[1, 2, 3], SetStat::RlmaxSet), set(&[3]));
    }

    #[test]
    fn test_long_swaps() {
        // [1,2,3]: 1 left of 2 adjacent, 2 left of 3 adjacent → 0 swaps
        assert_eq!(compute(&[1, 2, 3], Stat::Swaps), 0);
        // [3,2,1]: 1 right of 2, 2 right of 3 → 0 swaps
        assert_eq!(compute(&[3, 2, 1], Stat::Swaps), 0);
        // [1,3,2]: 1 left of 2 but not adjacent (gap=2) → swap; 2 right of 3 → no
        assert_eq!(compute(&[1, 3, 2], Stat::Swaps), 1);
        // [2,1,3]: inv: 1→pos1, 2→pos0, 3→pos2. i=1: pos1>pos0→no. i=2: pos0<pos2, gap=2→swap
        assert_eq!(compute(&[2, 1, 3], Stat::Swaps), 1);
        // [3,1,4,2]: inv: 1→pos1, 2→pos3, 3→pos0, 4→pos2
        // i=1: pos1 < pos3, gap=2 → swap
        // i=2: pos3 > pos0 → no
        // i=3: pos0 < pos2, gap=2 → swap
        assert_eq!(compute(&[3, 1, 4, 2], Stat::Swaps), 2);
    }

    #[test]
    fn test_long_swaps_set() {
        assert_eq!(compute_set(&[1, 2, 3], SetStat::SwapsSet), set(&[]));
        assert_eq!(compute_set(&[1, 3, 2], SetStat::SwapsSet), set(&[1]));
        assert_eq!(compute_set(&[3, 1, 4, 2], SetStat::SwapsSet), set(&[1, 3]));
    }

    #[test]
    fn test_swaps_set_consistent_with_scalar() {
        let perms: Vec<Vec<u8>> = {
            let mut p = vec![vec![1u8, 2, 3, 4, 5]];
            let mut current = vec![1u8, 2, 3, 4, 5];
            while next_perm(&mut current) {
                p.push(current.clone());
            }
            p
        };
        for pi in &perms {
            assert_eq!(
                compute_set(pi, SetStat::SwapsSet).len(),
                compute(pi, Stat::Swaps),
                "mismatch on {:?}",
                pi
            );
        }
    }

    #[test]
    fn test_descent_set_bitmask() {
        // [1,2,3]: no descents
        assert_eq!(descent_set_bitmask(&[1, 2, 3]), 0);
        // [3,2,1]: descents at positions 1,2 (0-indexed: 0,1)
        assert_eq!(descent_set_bitmask(&[3, 2, 1]), 0b11);
        // [1,3,2,4]: descent at position 2 (0-indexed: 1)
        assert_eq!(descent_set_bitmask(&[1, 3, 2, 4]), 0b10);
        // consistency: bitmask popcount == Des count
        let perms: Vec<Vec<u8>> = {
            let mut p = vec![vec![1u8, 2, 3, 4]];
            let mut current = vec![1u8, 2, 3, 4];
            while next_perm(&mut current) {
                p.push(current.clone());
            }
            p
        };
        for pi in &perms {
            assert_eq!(
                descent_set_bitmask(pi).count_ones() as usize,
                compute(pi, Stat::Des),
                "mismatch on {:?}",
                pi
            );
        }
    }

    #[test]
    fn test_checked_descent_set_bitmask_rejects_unrepresentable_descents() {
        let mut high_descent = vec![1; 66];
        high_descent[64] = 2;
        assert_eq!(checked_descent_set_bitmask(&high_descent), None);

        let mut long_but_representable = vec![1; 66];
        long_but_representable[0] = 2;
        assert_eq!(
            checked_descent_set_bitmask(&long_but_representable),
            Some(1)
        );
    }

    #[test]
    #[should_panic(expected = "descent at position 65 or later")]
    fn test_descent_set_bitmask_reports_unrepresentable_descent() {
        let mut word = vec![1; 66];
        word[64] = 2;
        descent_set_bitmask(&word);
    }

    #[test]
    fn test_swaps_known_alternating_polys() {
        // From the paper: H_n(t) coefficient vectors for alternating perms by swaps
        // H_4 = [1, 3, 1], H_5 = [1, 7, 7, 1]
        // Verify by brute force on small n
        let expected: Vec<Vec<i64>> = vec![
            vec![1, 3, 1],    // n=4
            vec![1, 7, 7, 1], // n=5
        ];
        for (idx, coeffs) in expected.iter().enumerate() {
            let n = (idx + 4) as u8;
            let mut poly = vec![0i64; coeffs.len()];
            let mut perm = (1..=n).collect::<Vec<u8>>();
            loop {
                let is_alt = (0..n as usize - 1).all(|i| {
                    if i % 2 == 0 {
                        perm[i] < perm[i + 1]
                    } else {
                        perm[i] > perm[i + 1]
                    }
                });
                if is_alt {
                    let s = compute(&perm, Stat::Swaps);
                    poly[s] += 1;
                }
                if !next_perm(&mut perm) {
                    break;
                }
            }
            assert_eq!(&poly, coeffs, "H_{} mismatch", n);
        }
    }

    #[test]
    fn test_set_stat_consistent_with_scalar() {
        // The cardinality of the set stat should match the scalar stat
        let perms: Vec<Vec<u8>> = {
            let mut p = vec![vec![1u8, 2, 3, 4]];
            let mut current = vec![1u8, 2, 3, 4];
            while next_perm(&mut current) {
                p.push(current.clone());
            }
            p
        };
        for pi in &perms {
            assert_eq!(
                compute_set(pi, SetStat::DesSet).len(),
                compute(pi, Stat::Des)
            );
            assert_eq!(
                compute_set(pi, SetStat::TieSet).len(),
                compute(pi, Stat::Tie)
            );
            assert_eq!(
                compute_set(pi, SetStat::PeakSet).len(),
                compute(pi, Stat::Peak)
            );
            assert_eq!(
                compute_set(pi, SetStat::WeakPeakSet).len(),
                compute(pi, Stat::WeakPeak)
            );
            assert_eq!(
                compute_set(pi, SetStat::ExcSet).len(),
                compute(pi, Stat::Exc)
            );
            assert_eq!(
                compute_set(pi, SetStat::FixSet).len(),
                compute(pi, Stat::Fix)
            );
            assert_eq!(
                compute_set(pi, SetStat::ValleySet).len(),
                compute(pi, Stat::Valley)
            );
            assert_eq!(
                compute_set(pi, SetStat::LrminSet).len(),
                compute(pi, Stat::Lrmin)
            );
            assert_eq!(
                compute_set(pi, SetStat::LrmaxSet).len(),
                compute(pi, Stat::Lrmax)
            );
            assert_eq!(
                compute_set(pi, SetStat::RlminSet).len(),
                compute(pi, Stat::Rlmin)
            );
            assert_eq!(
                compute_set(pi, SetStat::RlmaxSet).len(),
                compute(pi, Stat::Rlmax)
            );
            assert_eq!(
                compute_set(pi, SetStat::SwapsSet).len(),
                compute(pi, Stat::Swaps)
            );
        }
    }

    #[test]
    fn test_mahonian_equidistribution() {
        // inv and maj are equidistributed over S_n
        // For S_4, check that the generating polynomials match
        let perms: Vec<Vec<u8>> = {
            let mut p = vec![vec![1u8, 2, 3, 4]];
            // Generate all permutations of [1,2,3,4]
            let mut current = vec![1u8, 2, 3, 4];
            while next_perm(&mut current) {
                p.push(current.clone());
            }
            p
        };
        let mut inv_dist = [0u32; 7]; // max inv for S_4 is 6
        let mut maj_dist = [0u32; 7];
        for pi in &perms {
            inv_dist[compute(pi, Stat::Inv)] += 1;
            maj_dist[compute(pi, Stat::Maj)] += 1;
        }
        assert_eq!(inv_dist, maj_dist);
    }
}

// Helper for test
#[cfg(test)]
fn next_perm(perm: &mut [u8]) -> bool {
    let n = perm.len();
    if n <= 1 {
        return false;
    }
    let mut i = n - 2;
    loop {
        if perm[i] < perm[i + 1] {
            break;
        }
        if i == 0 {
            return false;
        }
        i -= 1;
    }
    let mut j = n - 1;
    while perm[j] <= perm[i] {
        j -= 1;
    }
    perm.swap(i, j);
    perm[i + 1..].reverse();
    true
}
