//! Parking functions: enumeration and validation.

use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};

/// Generate all parking functions of size n.
///
/// A parking function is a sequence (a_1, ..., a_n) with 1 <= a_i <= n
/// such that the sorted sequence b satisfies b_i <= i for all i.
/// The count is (n+1)^{n-1}.
pub fn all_parking_functions(n: u8) -> Vec<Vec<u8>> {
    let m = n as usize;
    let mut result = Vec::new();
    let mut current = Vec::with_capacity(m);
    let mut freq = vec![0usize; m];
    gen_pf(n, m, &mut current, &mut freq, &mut result);
    result
}

/// Descent polynomials of Smirnov words, refined by the final letter.
///
/// The output is indexed by the final letter in increasing order.  If
/// `S[r][j]` denotes the polynomial for words of length `r` ending in `j`,
/// appending `j` gives
///
/// `S[r + 1][j] = sum_{i < j} S[r][i] + t * sum_{i > j} S[r][i]`.
///
/// # Panics
///
/// Panics unless both parameters are positive.
pub fn smirnov_last_letter_descent_polynomials_bigint(
    length: usize,
    alphabet_size: usize,
) -> Vec<Vec<BigInt>> {
    assert!(length > 0, "Smirnov word length must be positive");
    assert!(alphabet_size > 0, "alphabet size must be positive");

    let mut refined = vec![vec![BigInt::from(1)]; alphabet_size];
    for _ in 1..length {
        let mut next = Vec::with_capacity(alphabet_size);
        for final_letter in 0..alphabet_size {
            let mut polynomial = vec![BigInt::zero()];
            for previous_letter in 0..final_letter {
                add_polynomial_assign(&mut polynomial, &refined[previous_letter], 0);
            }
            for previous_letter in final_letter + 1..alphabet_size {
                add_polynomial_assign(&mut polynomial, &refined[previous_letter], 1);
            }
            next.push(polynomial);
        }
        refined = next;
    }
    refined
}

/// The Chow polynomial of the noncrossing partition lattice `NC_{n+1}`.
///
/// Equivalently, this is the descent polynomial of tieless parking functions
/// of length `n`.  It is computed as `1/(n+1)` times the descent polynomial of
/// Smirnov words of length `n` on an `(n+1)`-letter alphabet.
pub fn noncrossing_chow_polynomial_bigint(n: usize) -> Vec<BigInt> {
    if n == 0 {
        return vec![BigInt::from(1)];
    }

    let refined = smirnov_last_letter_descent_polynomials_bigint(n, n + 1);
    let mut polynomial = vec![BigInt::zero()];
    for final_letter_polynomial in refined {
        add_polynomial_assign(&mut polynomial, &final_letter_polynomial, 0);
    }

    let divisor = BigInt::from(n + 1);
    for coefficient in &mut polynomial {
        assert!(
            (&*coefficient % &divisor).is_zero(),
            "normalized Smirnov coefficient is not integral"
        );
        *coefficient /= &divisor;
    }
    trim_trailing_zeros(&mut polynomial);
    polynomial
}

/// Checked `i64` wrapper for [`noncrossing_chow_polynomial_bigint`].
///
/// # Panics
///
/// Panics if a coefficient does not fit in `i64`.
pub fn noncrossing_chow_polynomial(n: usize) -> Vec<i64> {
    noncrossing_chow_polynomial_bigint(n)
        .into_iter()
        .map(|coefficient| {
            coefficient
                .to_i64()
                .expect("noncrossing Chow coefficient does not fit in i64")
        })
        .collect()
}

fn add_polynomial_assign(target: &mut Vec<BigInt>, source: &[BigInt], shift: usize) {
    let required_length = source.len() + shift;
    if target.len() < required_length {
        target.resize(required_length, BigInt::zero());
    }
    for (exponent, coefficient) in source.iter().enumerate() {
        target[exponent + shift] += coefficient;
    }
}

fn trim_trailing_zeros(polynomial: &mut Vec<BigInt>) {
    while polynomial.len() > 1 && polynomial.last().is_some_and(Zero::is_zero) {
        polynomial.pop();
    }
}

fn gen_pf(
    n: u8,
    m: usize,
    current: &mut Vec<u8>,
    freq: &mut Vec<usize>,
    result: &mut Vec<Vec<u8>>,
) {
    let placed = current.len();

    if placed == m {
        result.push(current.clone());
        return;
    }

    let remaining = m - placed - 1;

    for v in 1..=n {
        freq[v as usize - 1] += 1;

        // Pruning: for each k, we need cumsum(freq[0..=k]) + remaining >= k+1
        let mut feasible = true;
        let mut cumul = 0usize;
        for (k, &freq_k) in freq.iter().enumerate() {
            cumul += freq_k;
            if cumul + remaining < k + 1 {
                feasible = false;
                break;
            }
        }

        if feasible {
            current.push(v);
            gen_pf(n, m, current, freq, result);
            current.pop();
        }

        freq[v as usize - 1] -= 1;
    }
}

// ---------------------------------------------------------------------------
// Run-sorted parking function generator (fused pruning, zero allocation)
// ---------------------------------------------------------------------------

/// How to break a word into maximal ascending runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunBreak {
    /// Strictly ascending: break when w[i] >= w[i+1] (ties break the run).
    StrictAsc,
    /// Non-decreasing: break when w[i] > w[i+1] (ties continue the run).
    NonDecr,
}

/// How to order consecutive runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RunSort {
    /// min(R_1) < min(R_2) < ... < min(R_k).
    StrictMin,
    /// min(R_1) <= min(R_2) <= ... <= min(R_k).
    WeakMin,
    /// R_1 <_lex R_2 <_lex ... <_lex R_k.
    Lex,
}

/// Iterate over all run-sorted parking functions of size `n`, invoking
/// `callback` on each one.  Nothing is stored in memory beyond the
/// current partial word and O(n) bookkeeping.
///
/// The run-sorted condition and PF feasibility are checked incrementally
/// so that the vast majority of non-qualifying branches are pruned early.
pub fn for_each_runsorted_pf(
    n: u8,
    run_break: RunBreak,
    run_sort: RunSort,
    callback: &mut impl FnMut(&[u8]),
) {
    if n == 0 {
        return;
    }
    let m = n as usize;
    let mut word = Vec::with_capacity(m);
    let mut freq = vec![0usize; m];
    gen_rs_pf(
        n, m, &mut word, &mut freq, run_break, run_sort, 0, // cur_run_start
        0, // prev_run_start
        0, // prev_run_end
        0, // num_completed_runs
        callback,
    );
}

fn gen_rs_pf(
    n: u8,
    m: usize,
    word: &mut Vec<u8>,
    freq: &mut [usize],
    run_break: RunBreak,
    run_sort: RunSort,
    cur_run_start: usize,
    prev_run_start: usize,
    prev_run_end: usize,
    num_completed: usize,
    callback: &mut impl FnMut(&[u8]),
) {
    let pos = word.len();

    if pos == m {
        // Leaf: for lex, verify the final (still-open) run against its predecessor.
        if run_sort == RunSort::Lex && num_completed > 0 {
            if word[cur_run_start..] <= word[prev_run_start..prev_run_end] {
                return;
            }
        }
        callback(word);
        return;
    }

    let remaining = m - pos - 1;

    for v in 1..=n {
        freq[v as usize - 1] += 1;

        // PF feasibility: sorted prefix condition.
        let mut feasible = true;
        let mut cumul = 0usize;
        for (k, &fk) in freq.iter().enumerate() {
            cumul += fk;
            if cumul + remaining < k + 1 {
                feasible = false;
                break;
            }
        }

        if feasible {
            if pos == 0 {
                // First position: start first run.
                word.push(v);
                gen_rs_pf(n, m, word, freq, run_break, run_sort, 0, 0, 0, 0, callback);
                word.pop();
            } else {
                let is_descent = match run_break {
                    RunBreak::StrictAsc => v <= word[pos - 1],
                    RunBreak::NonDecr => v < word[pos - 1],
                };

                if !is_descent {
                    // Continue current run.
                    word.push(v);
                    gen_rs_pf(
                        n,
                        m,
                        word,
                        freq,
                        run_break,
                        run_sort,
                        cur_run_start,
                        prev_run_start,
                        prev_run_end,
                        num_completed,
                        callback,
                    );
                    word.pop();
                } else {
                    // Descent: current run word[cur_run_start..pos] just completed;
                    // new run starts at `pos` with value `v`.
                    let ok = match run_sort {
                        RunSort::StrictMin => {
                            // v = min(new run); word[cur_run_start] = min(just-completed run).
                            v > word[cur_run_start]
                        }
                        RunSort::WeakMin => v >= word[cur_run_start],
                        RunSort::Lex => {
                            // 1. If there is an earlier run, verify the just-completed
                            //    run against it.
                            let prev_ok = num_completed == 0
                                || word[cur_run_start..pos] > word[prev_run_start..prev_run_end];
                            // 2. Early prune: if v < first element of just-completed
                            //    run, the new run will be lex-smaller for sure.
                            prev_ok && v >= word[cur_run_start]
                        }
                    };

                    if ok {
                        word.push(v);
                        gen_rs_pf(
                            n,
                            m,
                            word,
                            freq,
                            run_break,
                            run_sort,
                            pos,           // new cur_run_start
                            cur_run_start, // new prev_run_start
                            pos,           // new prev_run_end
                            num_completed + 1,
                            callback,
                        );
                        word.pop();
                    }
                }
            }
        }

        freq[v as usize - 1] -= 1;
    }
}

/// Check if a sequence is a parking function.
pub fn is_parking_function(a: &[u8]) -> bool {
    let mut sorted: Vec<u8> = a.to_vec();
    sorted.sort();
    sorted
        .iter()
        .enumerate()
        .all(|(i, &v)| v != 0 && (v as usize) <= a.len() && (v as usize) <= i + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parking_function_counts() {
        // |PF(n)| = (n+1)^{n-1}
        assert_eq!(all_parking_functions(1).len(), 1); // 2^0 = 1
        assert_eq!(all_parking_functions(2).len(), 3); // 3^1 = 3
        assert_eq!(all_parking_functions(3).len(), 16); // 4^2 = 16
        assert_eq!(all_parking_functions(4).len(), 125); // 5^3 = 125
        assert_eq!(all_parking_functions(5).len(), 1296); // 6^4 = 1296
    }

    #[test]
    fn test_is_parking_function() {
        assert!(is_parking_function(&[1, 1, 2]));
        assert!(is_parking_function(&[1, 2, 1]));
        assert!(is_parking_function(&[2, 1, 1]));
        assert!(is_parking_function(&[1, 1, 1]));
        assert!(!is_parking_function(&[1, 3, 3]));
        assert!(!is_parking_function(&[2, 2, 2]));
        assert!(is_parking_function(&[1]));
        assert!(is_parking_function(&[]));
        assert!(!is_parking_function(&[0]));
        assert!(!is_parking_function(&[0, 1]));
        assert!(!is_parking_function(&[2]));
        assert!(!is_parking_function(&[1, 1, 4]));
    }

    #[test]
    fn test_all_are_valid() {
        for pf in all_parking_functions(4) {
            assert!(is_parking_function(&pf), "Invalid PF: {:?}", pf);
        }
    }

    #[test]
    fn test_smirnov_last_letter_descent_polynomials() {
        let refined = smirnov_last_letter_descent_polynomials_bigint(3, 4);
        let expected = vec![vec![0, 6, 3], vec![0, 8, 1], vec![1, 8], vec![3, 6]];
        let expected: Vec<Vec<BigInt>> = expected
            .into_iter()
            .map(|row| row.into_iter().map(BigInt::from).collect())
            .collect();
        assert_eq!(refined, expected);
    }

    #[test]
    fn test_noncrossing_chow_rows() {
        let expected = [
            vec![1],
            vec![1, 1],
            vec![1, 7, 1],
            vec![1, 31, 31, 1],
            vec![1, 116, 391, 116, 1],
            vec![1, 407, 3480, 3480, 407, 1],
        ];
        for (n, row) in expected.into_iter().enumerate() {
            assert_eq!(noncrossing_chow_polynomial(n + 1), row);
        }
    }

    #[test]
    fn test_chow_rows_match_tieless_parking_functions() {
        for n in 1..=5u8 {
            let mut brute = vec![0i64; n as usize];
            for word in all_parking_functions(n) {
                if word.windows(2).any(|pair| pair[0] == pair[1]) {
                    continue;
                }
                let descents = word.windows(2).filter(|pair| pair[0] > pair[1]).count();
                brute[descents] += 1;
            }
            while brute.len() > 1 && brute.last() == Some(&0) {
                brute.pop();
            }
            assert_eq!(noncrossing_chow_polynomial(n as usize), brute);
        }
    }

    #[test]
    fn test_chow_row_sums() {
        for n in 1..=12usize {
            let total: BigInt = noncrossing_chow_polynomial_bigint(n).into_iter().sum();
            assert_eq!(total, BigInt::from(n).pow((n - 1) as u32));
        }
    }

    /// Verify the fused run-sorted generator matches naive filter for all 6 variants.
    #[test]
    fn test_runsorted_matches_naive() {
        fn ascending_runs(w: &[u8]) -> Vec<&[u8]> {
            let mut runs = Vec::new();
            let mut start = 0;
            for i in 1..w.len() {
                if w[i] <= w[i - 1] {
                    runs.push(&w[start..i]);
                    start = i;
                }
            }
            runs.push(&w[start..]);
            runs
        }
        fn nondecr_runs(w: &[u8]) -> Vec<&[u8]> {
            let mut runs = Vec::new();
            let mut start = 0;
            for i in 1..w.len() {
                if w[i] < w[i - 1] {
                    runs.push(&w[start..i]);
                    start = i;
                }
            }
            runs.push(&w[start..]);
            runs
        }
        fn mins_strict(runs: &[&[u8]]) -> bool {
            runs.len() <= 1
                || (1..runs.len()).all(|i| runs[i].iter().min() > runs[i - 1].iter().min())
        }
        fn mins_weak(runs: &[&[u8]]) -> bool {
            runs.len() <= 1
                || (1..runs.len()).all(|i| runs[i].iter().min() >= runs[i - 1].iter().min())
        }
        fn lex(runs: &[&[u8]]) -> bool {
            runs.len() <= 1 || (1..runs.len()).all(|i| runs[i] > runs[i - 1])
        }

        let variants: Vec<(RunBreak, RunSort, Box<dyn Fn(&[u8]) -> bool>)> = vec![
            (
                RunBreak::StrictAsc,
                RunSort::StrictMin,
                Box::new(|w: &[u8]| mins_strict(&ascending_runs(w))),
            ),
            (
                RunBreak::StrictAsc,
                RunSort::WeakMin,
                Box::new(|w: &[u8]| mins_weak(&ascending_runs(w))),
            ),
            (
                RunBreak::StrictAsc,
                RunSort::Lex,
                Box::new(|w: &[u8]| lex(&ascending_runs(w))),
            ),
            (
                RunBreak::NonDecr,
                RunSort::StrictMin,
                Box::new(|w: &[u8]| mins_strict(&nondecr_runs(w))),
            ),
            (
                RunBreak::NonDecr,
                RunSort::WeakMin,
                Box::new(|w: &[u8]| mins_weak(&nondecr_runs(w))),
            ),
            (
                RunBreak::NonDecr,
                RunSort::Lex,
                Box::new(|w: &[u8]| lex(&nondecr_runs(w))),
            ),
        ];

        for n in 1..=6u8 {
            let all = all_parking_functions(n);
            for (rb, rs, filter) in &variants {
                let naive_count = all.iter().filter(|w| filter(w)).count();
                let mut fused_count = 0usize;
                for_each_runsorted_pf(n, *rb, *rs, &mut |_| fused_count += 1);
                assert_eq!(
                    fused_count, naive_count,
                    "Mismatch at n={}, rb={:?}, rs={:?}: fused={}, naive={}",
                    n, rb, rs, fused_count, naive_count
                );
            }
        }
    }
}
