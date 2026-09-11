//! Standard polynomial sequences for combinatorial research and testing.
//!
//! Each function returns a vector of coefficient vectors (ascending degree order),
//! indexed starting from n=0 (or n=1 where noted).
//!
//! These sequences exercise many properties in [`crate::real_rootedness`]:
//!
//! | Sequence | Palindromic | Gamma-positive | Real-rooted | Consecutive interlace |
//! |---|---|---|---|---|
//! | Eulerian A_n | yes | yes | yes | yes |
//! | Narayana N_n | yes | yes | yes | yes |
//! | Type B Eulerian B_n | yes | yes | yes | yes |
//! | Chebyshev T_n, U_n | no | no | yes | yes |
//! | Hermite He_n | no | no | yes | yes |

use num_bigint::BigInt;

pub mod literature;

/// Compute Eulerian polynomials A_1(t), A_2(t), ..., A_n(t).
///
/// The Eulerian polynomial A_n(t) = sum_{k=0}^{n-1} A(n,k) t^k where A(n,k) is
/// the number of permutations of {1,...,n} with exactly k descents.
///
/// Recurrence: A_n(t) = (1 + (n-1)t) A_{n-1}(t) + t(1-t) A'_{n-1}(t).
///
/// Properties: palindromic, gamma-positive, real-rooted (all roots negative),
/// and consecutive Eulerian polynomials interlace.
///
/// Returns `polys[i]` = A_{i+1}(t), so `polys[0]` = A_1 = `[1]`.
pub fn eulerian_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(eulerian_polynomials_bigint(max_n))
}

/// Compute Eulerian polynomials with arbitrary-size integer coefficients.
pub fn eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    if max_n == 0 {
        return vec![];
    }
    let mut polys = Vec::with_capacity(max_n);
    polys.push(vec![BigInt::from(1)]); // A_1 = 1

    for n in 2..=max_n {
        let prev = &polys[n - 2];
        let d = prev.len();

        // Derivative of prev
        let mut dp = vec![BigInt::from(0); d.saturating_sub(1)];
        for k in 1..d {
            dp[k - 1] = &prev[k] * usize_to_bigint(k);
        }

        // (1 + (n-1)t) * prev
        let mut term1 = vec![BigInt::from(0); d + 1];
        let n_minus_one = usize_to_bigint(n - 1);
        for k in 0..d {
            term1[k] += &prev[k];
            term1[k + 1] += &prev[k] * &n_minus_one;
        }

        // t(1-t) * dp = t*dp - t^2*dp
        let dp_len = dp.len();
        let mut term2 = vec![BigInt::from(0); dp_len + 2];
        for k in 0..dp_len {
            term2[k + 1] += &dp[k];
            term2[k + 2] -= &dp[k];
        }

        // Sum
        let len = term1.len().max(term2.len());
        let mut result = vec![BigInt::from(0); len];
        for k in 0..term1.len() {
            result[k] += &term1[k];
        }
        for k in 0..term2.len() {
            result[k] += &term2[k];
        }
        trim_trailing_zeros_bigint(&mut result);
        polys.push(result);
    }
    polys
}

fn trim_trailing_zeros_bigint(coeffs: &mut Vec<BigInt>) {
    while coeffs.last() == Some(&BigInt::from(0)) {
        coeffs.pop();
    }
}

fn bigint_polys_to_i64(polys: Vec<Vec<BigInt>>) -> Vec<Vec<i64>> {
    polys
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|c| i64::try_from(&c).expect("sequence coefficient overflow"))
                .collect()
        })
        .collect()
}

fn usize_to_bigint(value: usize) -> BigInt {
    BigInt::from(value)
}

fn bigint_two() -> BigInt {
    BigInt::from(2)
}

fn bigint_zero_vec(len: usize) -> Vec<BigInt> {
    vec![BigInt::from(0); len]
}

fn add_polys_bigint(lhs: &[BigInt], rhs: &[BigInt]) -> Vec<BigInt> {
    let len = lhs.len().max(rhs.len());
    let mut result = bigint_zero_vec(len);
    for (k, c) in lhs.iter().enumerate() {
        result[k] += c;
    }
    for (k, c) in rhs.iter().enumerate() {
        result[k] += c;
    }
    trim_trailing_zeros_bigint(&mut result);
    result
}

fn sub_polys_bigint(lhs: &[BigInt], rhs: &[BigInt]) -> Vec<BigInt> {
    let len = lhs.len().max(rhs.len());
    let mut result = bigint_zero_vec(len);
    for (k, c) in lhs.iter().enumerate() {
        result[k] += c;
    }
    for (k, c) in rhs.iter().enumerate() {
        result[k] -= c;
    }
    trim_trailing_zeros_bigint(&mut result);
    result
}

fn shift_scale_bigint(poly: &[BigInt], scale: &BigInt, shift: usize) -> Vec<BigInt> {
    let mut result = bigint_zero_vec(poly.len() + shift);
    for (k, c) in poly.iter().enumerate() {
        result[k + shift] = c * scale;
    }
    trim_trailing_zeros_bigint(&mut result);
    result
}

fn scale_bigint(poly: &[BigInt], scale: &BigInt) -> Vec<BigInt> {
    let mut result = poly.iter().map(|c| c * scale).collect::<Vec<_>>();
    trim_trailing_zeros_bigint(&mut result);
    result
}

fn binomial_row_bigint(n: usize) -> Vec<BigInt> {
    let mut row = vec![BigInt::from(1); n + 1];
    for j in 1..=n {
        row[j] = &row[j - 1] * usize_to_bigint(n - j + 1) / usize_to_bigint(j);
    }
    row
}

fn exact_div_bigint(numerator: BigInt, denominator: BigInt) -> BigInt {
    debug_assert!(denominator != BigInt::from(0));
    debug_assert_eq!(&numerator % &denominator, BigInt::from(0));
    numerator / denominator
}

/// Compute Narayana polynomials with arbitrary-size integer coefficients.
pub fn narayana_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    if max_n == 0 {
        return vec![];
    }
    let mut polys = Vec::with_capacity(max_n);

    for n in 1..=max_n {
        let binom_n = binomial_row_bigint(n);
        let n_big = usize_to_bigint(n);
        let mut coeffs = Vec::with_capacity(n);
        for k in 0..n {
            let j = k + 1;
            coeffs.push(exact_div_bigint(
                &binom_n[j] * &binom_n[j - 1],
                n_big.clone(),
            ));
        }
        trim_trailing_zeros_bigint(&mut coeffs);
        polys.push(coeffs);
    }
    polys
}

/// Compute type-B Eulerian polynomials with arbitrary-size integer coefficients.
pub fn type_b_eulerian_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]); // B_0 = 1

    for n in 1..=max_n {
        let prev = &polys[n - 1];
        let d = prev.len();

        let mut dp = bigint_zero_vec(d.saturating_sub(1));
        for k in 1..d {
            dp[k - 1] = &prev[k] * usize_to_bigint(k);
        }

        let mut term1 = bigint_zero_vec(d + 1);
        let two_n_minus_one = usize_to_bigint(2 * n - 1);
        for k in 0..d {
            term1[k] += &prev[k];
            term1[k + 1] += &prev[k] * &two_n_minus_one;
        }

        let mut term2 = bigint_zero_vec(dp.len() + 2);
        for k in 0..dp.len() {
            let twice = &dp[k] * bigint_two();
            term2[k + 1] += &twice;
            term2[k + 2] -= &twice;
        }

        polys.push(add_polys_bigint(&term1, &term2));
    }
    polys
}

/// Compute Chebyshev polynomials of the first kind with arbitrary-size integer coefficients.
pub fn chebyshev_polynomials_t_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(1)]);

    for n in 2..=max_n {
        let term1 = shift_scale_bigint(&polys[n - 1], &bigint_two(), 1);
        polys.push(sub_polys_bigint(&term1, &polys[n - 2]));
    }
    polys
}

/// Compute Chebyshev polynomials of the second kind with arbitrary-size integer coefficients.
pub fn chebyshev_polynomials_u_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(2)]);

    for n in 2..=max_n {
        let term1 = shift_scale_bigint(&polys[n - 1], &bigint_two(), 1);
        polys.push(sub_polys_bigint(&term1, &polys[n - 2]));
    }
    polys
}

/// Compute probabilist's Hermite polynomials with arbitrary-size integer coefficients.
pub fn hermite_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut polys = Vec::with_capacity(max_n + 1);
    polys.push(vec![BigInt::from(1)]);
    if max_n == 0 {
        return polys;
    }
    polys.push(vec![BigInt::from(0), BigInt::from(1)]);

    for n in 2..=max_n {
        let term1 = shift_scale_bigint(&polys[n - 1], &BigInt::from(1), 1);
        let term2 = scale_bigint(&polys[n - 2], &usize_to_bigint(n - 1));
        polys.push(sub_polys_bigint(&term1, &term2));
    }
    polys
}

/// Compute Narayana polynomials N_1(t), N_2(t), ..., N_n(t).
///
/// The Narayana polynomial N_n(t) = sum_{k=0}^{n-1} N(n,k+1) t^k where
/// N(n,k) = (1/n) C(n,k) C(n,k-1) is the Narayana number.
///
/// These are the h-vectors of the associahedra and of the type-A Coxeter complex.
///
/// Properties: palindromic, gamma-positive, real-rooted (all roots negative),
/// and row sums are the Catalan numbers.
///
/// Returns `polys[i]` = N_{i+1}(t), so `polys[0]` = N_1 = `[1]`.
pub fn narayana_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(narayana_polynomials_bigint(max_n))
}

/// Compute type-B Eulerian polynomials B_0(t), B_1(t), ..., B_n(t).
///
/// The type-B Eulerian polynomial B_n(t) counts signed permutations of
/// {±1, ..., ±n} by number of descents (in type B sense).
///
/// Recurrence: B_n(t) = (1 + (2n-1)t) B_{n-1}(t) + 2t(1-t) B'_{n-1}(t).
///
/// Properties: palindromic, gamma-positive, real-rooted.
///
/// Returns `polys[i]` = B_i(t), so `polys[0]` = B_0 = `[1]`.
pub fn type_b_eulerian_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(type_b_eulerian_polynomials_bigint(max_n))
}

/// Compute Chebyshev polynomials of the first kind T_0(t), T_1(t), ..., T_n(t).
///
/// Recurrence: T_0 = 1, T_1 = t, T_n = 2t T_{n-1} - T_{n-2}.
///
/// Properties: real-rooted (roots are cos(kπ/n)), consecutive T_n interlace.
///
/// Returns `polys[i]` = T_i(t).
pub fn chebyshev_polynomials_t(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(chebyshev_polynomials_t_bigint(max_n))
}

/// Compute Chebyshev polynomials of the second kind U_0(t), U_1(t), ..., U_n(t).
///
/// Recurrence: U_0 = 1, U_1 = 2t, U_n = 2t U_{n-1} - U_{n-2}.
///
/// Properties: real-rooted, consecutive U_n interlace,
/// and T_n and U_{n-1} interlace (T'_n = n U_{n-1}).
///
/// Returns `polys[i]` = U_i(t).
pub fn chebyshev_polynomials_u(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(chebyshev_polynomials_u_bigint(max_n))
}

/// Compute probabilist's Hermite polynomials He_0(t), He_1(t), ..., He_n(t).
///
/// Recurrence: He_0 = 1, He_1 = t, He_n = t He_{n-1} - (n-1) He_{n-2}.
///
/// Properties: real-rooted, consecutive He_n interlace.
/// These are the orthogonal polynomials for the standard Gaussian measure.
///
/// Returns `polys[i]` = He_i(t).
pub fn hermite_polynomials(max_n: usize) -> Vec<Vec<i64>> {
    bigint_polys_to_i64(hermite_polynomials_bigint(max_n))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::real_rootedness::*;

    #[test]
    fn test_eulerian_values() {
        let polys = eulerian_polynomials(5);
        assert_eq!(polys[0], vec![1]); // A_1 = 1
        assert_eq!(polys[1], vec![1, 1]); // A_2 = 1 + t
        assert_eq!(polys[2], vec![1, 4, 1]); // A_3 = 1 + 4t + t^2
        assert_eq!(polys[3], vec![1, 11, 11, 1]); // A_4
        assert_eq!(polys[4], vec![1, 26, 66, 26, 1]); // A_5
    }

    #[test]
    fn test_eulerian_properties() {
        let polys = eulerian_polynomials(10);
        for p in &polys {
            assert!(is_palindromic(p), "Eulerian not palindromic: {:?}", p);
            assert!(is_real_rooted(p), "Eulerian not real-rooted: {:?}", p);
            assert!(is_gamma_positive(p), "Eulerian not gamma-positive: {:?}", p);
        }
        // Consecutive interlacing (degree diff = 1)
        for i in 1..polys.len() {
            let result = check_interlacing(&polys[i - 1], &polys[i]);
            assert_eq!(
                result,
                Some(true),
                "A_{} and A_{} don't interlace",
                i + 1,
                i
            );
        }
    }

    #[test]
    fn test_narayana_values() {
        let polys = narayana_polynomials(5);
        assert_eq!(polys[0], vec![1]); // N_1 = 1
        assert_eq!(polys[1], vec![1, 1]); // N_2 = 1 + t
        assert_eq!(polys[2], vec![1, 3, 1]); // N_3 = 1 + 3t + t^2
        assert_eq!(polys[3], vec![1, 6, 6, 1]); // N_4
        assert_eq!(polys[4], vec![1, 10, 20, 10, 1]); // N_5
    }

    #[test]
    fn test_narayana_properties() {
        let polys = narayana_polynomials(10);
        for p in &polys {
            assert!(is_palindromic(p), "Narayana not palindromic: {:?}", p);
            assert!(is_real_rooted(p), "Narayana not real-rooted: {:?}", p);
            assert!(is_gamma_positive(p), "Narayana not gamma-positive: {:?}", p);
        }
        // Row sums are Catalan numbers: 1, 2, 5, 14, 42, ...
        let catalan = [1i64, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796];
        for (i, p) in polys.iter().enumerate() {
            let sum: i64 = p.iter().sum();
            assert_eq!(sum, catalan[i], "N_{} row sum wrong", i + 1);
        }
    }

    #[test]
    #[should_panic(expected = "sequence coefficient overflow")]
    fn test_narayana_rejects_i64_overflow() {
        let _ = narayana_polynomials(67);
    }

    #[test]
    fn test_type_b_eulerian_values() {
        let polys = type_b_eulerian_polynomials(4);
        assert_eq!(polys[0], vec![1]); // B_0 = 1
        assert_eq!(polys[1], vec![1, 1]); // B_1 = 1 + t
        assert_eq!(polys[2], vec![1, 6, 1]); // B_2 = 1 + 6t + t^2
        assert_eq!(polys[3], vec![1, 23, 23, 1]); // B_3
                                                  // Row sums: (2n)! / 2^n = 1, 1, 3, 15, 105, ...
                                                  // Actually B_n(1) = number of elements in hyperoctahedral group with...
                                                  // Actually sum of type B Eulerian numbers for B_n is 2^n n! / ... let me just verify
                                                  // B_0(1)=1, B_1(1)=2, B_2(1)=8, B_3(1)=48
        assert_eq!(polys[0].iter().sum::<i64>(), 1);
        assert_eq!(polys[1].iter().sum::<i64>(), 2);
        assert_eq!(polys[2].iter().sum::<i64>(), 8);
        assert_eq!(polys[3].iter().sum::<i64>(), 48);
    }

    #[test]
    fn test_type_b_eulerian_properties() {
        let polys = type_b_eulerian_polynomials(8);
        for (i, p) in polys.iter().enumerate() {
            assert!(is_palindromic(p), "B_{} not palindromic", i);
            assert!(is_real_rooted(p), "B_{} not real-rooted", i);
            assert!(is_gamma_positive(p), "B_{} not gamma-positive", i);
        }
    }

    #[test]
    fn test_chebyshev_t_values() {
        let polys = chebyshev_polynomials_t(5);
        assert_eq!(polys[0], vec![1]); // T_0 = 1
        assert_eq!(polys[1], vec![0, 1]); // T_1 = t
        assert_eq!(polys[2], vec![-1, 0, 2]); // T_2 = 2t^2 - 1
        assert_eq!(polys[3], vec![0, -3, 0, 4]); // T_3 = 4t^3 - 3t
        assert_eq!(polys[4], vec![1, 0, -8, 0, 8]); // T_4 = 8t^4 - 8t^2 + 1
    }

    #[test]
    fn test_chebyshev_real_rooted() {
        let t_polys = chebyshev_polynomials_t(12);
        let u_polys = chebyshev_polynomials_u(12);
        for (i, p) in t_polys.iter().enumerate() {
            assert!(is_real_rooted(p), "T_{} not real-rooted", i);
        }
        for (i, p) in u_polys.iter().enumerate() {
            assert!(is_real_rooted(p), "U_{} not real-rooted", i);
        }
    }

    #[test]
    #[should_panic(expected = "sequence coefficient overflow")]
    fn test_chebyshev_t_rejects_i64_overflow() {
        let _ = chebyshev_polynomials_t(64);
    }

    #[test]
    fn test_hermite_values() {
        let polys = hermite_polynomials(5);
        assert_eq!(polys[0], vec![1]); // He_0 = 1
        assert_eq!(polys[1], vec![0, 1]); // He_1 = t
        assert_eq!(polys[2], vec![-1, 0, 1]); // He_2 = t^2 - 1
        assert_eq!(polys[3], vec![0, -3, 0, 1]); // He_3 = t^3 - 3t
        assert_eq!(polys[4], vec![3, 0, -6, 0, 1]); // He_4 = t^4 - 6t^2 + 3
    }

    #[test]
    fn test_hermite_real_rooted() {
        let polys = hermite_polynomials(12);
        for (i, p) in polys.iter().enumerate() {
            assert!(is_real_rooted(p), "He_{} not real-rooted", i);
        }
    }

    #[test]
    fn test_chebyshev_interlacing() {
        // T'_n = n * U_{n-1}, so T_n and U_{n-1} interlace.
        let t_polys = chebyshev_polynomials_t(10);
        let u_polys = chebyshev_polynomials_u(10);
        for n in 2..t_polys.len() {
            let result = check_interlacing(&u_polys[n - 1], &t_polys[n]);
            assert_eq!(
                result,
                Some(true),
                "T_{} and U_{} don't interlace",
                n,
                n - 1
            );
        }
    }

    #[test]
    fn test_hermite_interlacing() {
        // Consecutive Hermite polynomials interlace.
        let polys = hermite_polynomials(10);
        for i in 2..polys.len() {
            let result = check_interlacing(&polys[i - 1], &polys[i]);
            assert_eq!(
                result,
                Some(true),
                "He_{} and He_{} don't interlace",
                i,
                i - 1
            );
        }
    }
}
