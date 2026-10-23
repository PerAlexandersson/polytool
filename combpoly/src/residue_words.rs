//! Rotation CSPs for words with parity constraints on residue subsequences.
//!
//! A length-`m*a` word is split into its `m` subsequences of positions modulo
//! `m`. Letters `0..zero_letters` have bit zero; the remaining `one_letters`
//! letters have bit one. We require exactly one subsequence to have odd bit
//! sum. One-step word rotation preserves this family.
//!
//! If `E_a(q), O_a(q)` enumerate even/odd bit-sum length-`a` words by major
//! index, a CSP polynomial for this rotation is
//! `[m]_(q^a) O_a(q) E_a(q)^(m-1)`. It enumerates the statistic
//! `a*r + sum_j maj(row_j)`, where `r` is the zero-based odd residue class.
//! This follows by induction from the diagonal `C_a` rotation of the rows.
//! Ordinary major index on the interleaved word need not give a CSP.
use num_bigint::BigInt;
use num_traits::{One, Zero};
use polytool::Polynomial;

fn trim(mut coefficients: Vec<BigInt>) -> Vec<BigInt> {
    if coefficients.is_empty() {
        return vec![BigInt::zero()];
    }
    while coefficients.len() > 1 && coefficients.last().is_some_and(Zero::is_zero) {
        coefficients.pop();
    }
    coefficients
}

/// Major-index polynomials for even and odd bit-sum words, respectively.
///
/// The ordered alphabet is `0..zero_letters+one_letters`, with bit one on
/// its final `one_letters` letters. Major index sums the one-based positions
/// of strict descents. Coefficients are arbitrary-precision integers in
/// ascending degree order; zero is `[0]`. For length zero the answer is
/// `([1], [0])`. Either bit class may be empty; an empty alphabet is allowed
/// and has no positive-length words. Panics on `usize` overflow in alphabet
/// size or maximum degree.
///
/// A last-letter/parity DP uses `O(k^2 * length^3)` arithmetic operations
/// and `O(k * length^2)` coefficient storage, where `k` is alphabet size.
pub fn parity_word_major_polynomials_bigint(
    length: usize,
    zero_letters: usize,
    one_letters: usize,
) -> (Vec<BigInt>, Vec<BigInt>) {
    let alphabet = zero_letters
        .checked_add(one_letters)
        .expect("alphabet size overflow");
    if length == 0 {
        return (vec![BigInt::one()], vec![BigInt::zero()]);
    }
    if alphabet == 0 {
        return (vec![BigInt::zero()], vec![BigInt::zero()]);
    }
    let maximum_degree = length
        .checked_mul(length - 1)
        .expect("major-index degree overflow")
        / 2;
    let mut current = vec![vec![vec![BigInt::zero(); maximum_degree + 1]; alphabet]; 2];
    for letter in 0..alphabet {
        current[usize::from(letter >= zero_letters)][letter][0] += 1;
    }
    for prefix_length in 1..length {
        let mut next = vec![vec![vec![BigInt::zero(); maximum_degree + 1]; alphabet]; 2];
        let degree_bound = prefix_length * (prefix_length - 1) / 2;
        for (new_parity, by_last) in next.iter_mut().enumerate() {
            for (letter, target) in by_last.iter_mut().enumerate() {
                let old_parity = new_parity ^ usize::from(letter >= zero_letters);
                for (last, polynomial) in current[old_parity].iter().enumerate() {
                    let shift = if last > letter { prefix_length } else { 0 };
                    for degree in 0..=degree_bound {
                        target[degree + shift] += &polynomial[degree];
                    }
                }
            }
        }
        current = next;
    }
    let mut result = vec![vec![BigInt::zero(); maximum_degree + 1]; 2];
    for parity in 0..2 {
        for polynomial in &current[parity] {
            for (coefficient, value) in result[parity].iter_mut().zip(polynomial) {
                *coefficient += value;
            }
        }
    }
    (trim(result.remove(0)), trim(result.remove(0)))
}

fn spaced_q_integer(modulus: usize, width: usize) -> Polynomial<BigInt> {
    let degree = (modulus - 1)
        .checked_mul(width)
        .expect("q-integer degree overflow");
    let mut coefficients = vec![BigInt::zero(); degree + 1];
    for residue in 0..modulus {
        coefficients[residue * width] = BigInt::one();
    }
    Polynomial::new(coefficients)
}

/// CSP polynomial for rotation of length-`modulus*width` words with exactly
/// one odd bit-sum residue subsequence, for arbitrary loop multiplicities.
///
/// Returns `[modulus]_(q^width) O_width(q) E_width(q)^(modulus-1)` in ascending
/// degree order. This is the generating polynomial of
/// `width*r + sum_j maj(row_j)`, with `r` the unique odd residue class.
/// Reuses [`parity_word_major_polynomials_bigint`] and dense exact polynomial
/// multiplication. Empty alphabets or no odd letters give the zero polynomial.
/// Requires positive modulus and width; panics on degree overflow.
pub fn single_odd_residue_csp_polynomial_bigint(
    modulus: usize,
    width: usize,
    zero_letters: usize,
    one_letters: usize,
) -> Vec<BigInt> {
    assert!(
        modulus > 0 && width > 0,
        "positive modulus and width required"
    );
    let (even, odd) = parity_word_major_polynomials_bigint(width, zero_letters, one_letters);
    let even = Polynomial::new(even);
    let mut result = spaced_q_integer(modulus, width) * Polynomial::new(odd);
    for _ in 1..modulus {
        result = result * even.clone();
    }
    trim(result.coeffs().to_vec())
}

/// Compact binary CSP polynomial for the same rotation family:
/// `[modulus]_(q^width) * product_(j=1..width-1) (1+q^j)^modulus`.
///
/// This gives the same root values as
/// [`single_odd_residue_csp_polynomial_bigint`] with one letter of each bit.
/// It generally differs as a polynomial from that statistic's distribution.
/// Positive modulus and width are required; coefficients are ascending BigInts.
/// Uses dense polynomial products; output degree is
/// `width*(modulus-1) + modulus*width*(width-1)/2`.
pub fn binary_single_odd_residue_product_bigint(modulus: usize, width: usize) -> Vec<BigInt> {
    assert!(
        modulus > 0 && width > 0,
        "positive modulus and width required"
    );
    let mut result = spaced_q_integer(modulus, width);
    for j in 1..width {
        let mut factor = vec![BigInt::zero(); j + 1];
        factor[0] = BigInt::one();
        factor[j] = BigInt::one();
        let factor = Polynomial::new(factor);
        for _ in 0..modulus {
            result = result * factor.clone();
        }
    }
    trim(result.coeffs().to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial_builder::build_generating_polynomial;
    use crate::statistics::{compute, Stat};
    use polytool::cyclic_sieving::root_of_unity_evaluation_bigint;

    fn words(length: usize, alphabet: usize) -> Vec<Vec<u8>> {
        (0..alphabet.pow(length as u32))
            .map(|mut code| {
                (0..length)
                    .map(|_| {
                        let letter = (code % alphabet) as u8;
                        code /= alphabet;
                        letter
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn parity_major_matches_independent_enumeration() {
        for zero_letters in 0..=2 {
            for one_letters in 0..=2 {
                let alphabet = zero_letters + one_letters;
                for length in 0..=5 {
                    let (even, odd) =
                        parity_word_major_polynomials_bigint(length, zero_letters, one_letters);
                    for (parity, expected) in [even, odd].into_iter().enumerate() {
                        let subset: Vec<Vec<u8>> = words(length, alphabet)
                            .into_iter()
                            .filter(|word| {
                                word.iter().filter(|&&v| v as usize >= zero_letters).count() % 2
                                    == parity
                            })
                            .collect();
                        let actual: Vec<BigInt> = build_generating_polynomial(&subset, Stat::Maj)
                            .into_iter()
                            .map(BigInt::from)
                            .collect();
                        assert_eq!(expected, actual);
                    }
                }
            }
        }
    }

    #[test]
    fn residue_statistic_and_rotation_match() {
        for modulus in 1..=4 {
            for width in 1..=3 {
                let n = modulus * width;
                for (zeros, ones) in [(1usize, 1usize), (1, 2), (2, 1), (0, 1)] {
                    if (zeros + ones).pow(n as u32) > 100_000usize {
                        continue;
                    }
                    let polynomial =
                        single_odd_residue_csp_polynomial_bigint(modulus, width, zeros, ones);
                    let subset: Vec<(Vec<u8>, usize)> = words(n, zeros + ones)
                        .into_iter()
                        .filter_map(|word| {
                            let mut parity = vec![false; modulus];
                            for (i, &letter) in word.iter().enumerate() {
                                parity[i % modulus] ^= letter as usize >= zeros;
                            }
                            if parity.iter().filter(|&&odd| odd).count() != 1 {
                                return None;
                            }
                            let r = parity.iter().position(|&odd| odd).unwrap();
                            let major_sum: usize = (0..modulus)
                                .map(|j| {
                                    let row: Vec<u8> =
                                        word.iter().skip(j).step_by(modulus).copied().collect();
                                    compute(&row, Stat::Maj)
                                })
                                .sum();
                            Some((word, width * r + major_sum))
                        })
                        .collect();
                    let mut histogram = vec![BigInt::zero(); polynomial.len()];
                    for (_, value) in &subset {
                        histogram[*value] += 1;
                    }
                    assert_eq!(trim(histogram), polynomial);
                    for shift in 0..n {
                        let fixed = subset
                            .iter()
                            .filter(|(word, _)| (0..n).all(|i| word[i] == word[(i + shift) % n]))
                            .count();
                        assert_eq!(
                            root_of_unity_evaluation_bigint(&polynomial, n, shift).integer_value,
                            Some(BigInt::from(fixed)),
                            "modulus={modulus} width={width} loops=({zeros},{ones}) shift={shift}"
                        );
                        if zeros == 1 && ones == 1 {
                            let compact = binary_single_odd_residue_product_bigint(modulus, width);
                            assert_eq!(
                                root_of_unity_evaluation_bigint(&compact, n, shift).integer_value,
                                Some(BigInt::from(fixed))
                            );
                        }
                    }
                }
            }
        }
    }
}
