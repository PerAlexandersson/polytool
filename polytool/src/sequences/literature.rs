//! Exact polynomial families for reproducible literature-conjecture checks.
//!
//! Coefficients are in ascending degree order. No real-rootedness assertion
//! for a conjectural family is part of these generation APIs.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

/// Delannoy triangle rows D_0(x), ..., D_max_n(x).
///
/// D(n,k) counts paths to (n-k,k) with steps (1,0), (0,1), (1,1).
/// D_0=1, D_1=1+x, D_n=(1+x)D_(n-1)+x D_(n-2).
pub fn delannoy_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut rows = vec![vec![BigInt::one()]];
    for n in 1..=max_n {
        let mut row = vec![BigInt::zero(); n + 1];
        add_linear_multiple(&mut row, &rows[n - 1], 1, 1);
        if n >= 2 {
            add_linear_multiple(&mut row, &rows[n - 2], 0, 1);
        }
        rows.push(row);
    }
    rows
}

/// Row-generating polynomials of the square of the Delannoy triangle, n=0..max_n.
///
/// Uses the four-term recurrence from Section 4 of Mao--Wang,
/// "The Narayana transformation", arXiv:2607.01572v1:
/// G_n=(x+3)G_(n-1)+(x-1)G_(n-2)+(x-1)G_(n-3)+xG_(n-4).
/// Initial rows: 1; 2+x; 5+6x+x^2; 12+25x+10x^2+x^3.
pub fn delannoy_square_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let seeds: [&[i64]; 4] = [&[1], &[2, 1], &[5, 6, 1], &[12, 25, 10, 1]];
    let mut rows: Vec<Vec<BigInt>> = seeds
        .iter()
        .take(max_n.saturating_add(1))
        .map(|r| r.iter().map(|&x| BigInt::from(x)).collect())
        .collect();
    for n in 4..=max_n {
        let mut row = vec![BigInt::zero(); n + 1];
        add_linear_multiple(&mut row, &rows[n - 1], 3, 1);
        add_linear_multiple(&mut row, &rows[n - 2], -1, 1);
        add_linear_multiple(&mut row, &rows[n - 3], -1, 1);
        add_linear_multiple(&mut row, &rows[n - 4], 0, 1);
        rows.push(row);
    }
    rows
}

/// Eulerian matrix-square rows, n=0..max_n, in Mao--Wang's normalization.
///
/// A(0,0)=1; for n>=1, A(n,k) counts permutations of n with k-1 descents.
/// Thus the nonzero-row basis is x times the usual Eulerian polynomial.
/// This is triangular matrix multiplication, NOT coefficientwise squaring.
/// The output retains its structural zero constant term for n>=1.
pub fn eulerian_square_polynomials_bigint(max_n: usize) -> Vec<Vec<BigInt>> {
    let mut triangle = vec![vec![BigInt::one()]];
    for row in super::eulerian_polynomials_bigint(max_n) {
        let mut shifted = vec![BigInt::zero()];
        shifted.extend(row);
        triangle.push(shifted);
    }
    square_triangle(&triangle)
}

/// Hoggatt polynomials H_1^[m](q,t), ..., H_max_n^[m](q,t), exact integer q>=1.
///
/// The coefficient of t^k in row n (1<=n<=max_n, 0<=k<n) is
/// q^(m*k*(k+1)/2) times the volume generating function of plane partitions
/// in a k by (n-1-k) by m box, evaluated at q. This is the normalization of
/// Fang--Zhang--Zhao, arXiv:2505.05873v1, Section 1 and Conjecture 4.1.
///
/// MacMahon's product is evaluated as a quotient of Gaussian binomial
/// products, also valid at q=1. Every division is checked to be exact.
/// m=1 gives product_(i=1)^(n-1)(1+q^i*t); m=2,q=1 gives Narayana rows;
/// m=3,q=1 gives Baxter descent rows. No root property is assumed.
///
/// # Panics
/// Panics for m=0 or q=0, or index/exponent arithmetic exceeding usize.
pub fn hoggatt_polynomials_bigint(max_n: usize, m: usize, q: u32) -> Vec<Vec<BigInt>> {
    assert!(m >= 1, "Hoggatt box height m must be positive");
    assert!(q >= 1, "Hoggatt specialization q must be positive");
    if max_n == 0 {
        return Vec::new();
    }
    let max_top = max_n.checked_add(m).expect("Hoggatt index overflow") - 2;
    let q = BigInt::from(q);
    let mut powers = vec![BigInt::one()];
    for i in 1..=max_top {
        powers.push(&powers[i - 1] * &q);
    }
    let mut gaussian = vec![vec![BigInt::one()]];
    for n in 1..=max_top {
        let mut row = vec![BigInt::one(); n + 1];
        for k in 1..n {
            row[k] = &gaussian[n - 1][k - 1] + &powers[k] * &gaussian[n - 1][k];
        }
        gaussian.push(row);
    }
    (1..=max_n)
        .map(|n| {
            let top = n + m - 2;
            let denominator: BigInt = (1..m).map(|j| &gaussian[top][j]).product();
            (0..n)
                .map(|k| {
                    let numerator: BigInt = (0..m).map(|i| &gaussian[top][k + i]).product();
                    let (count, remainder) = numerator.div_rem(&denominator);
                    assert!(remainder.is_zero(), "nonintegral MacMahon evaluation");
                    let twice_triangular = k.checked_mul(k + 1).expect("Hoggatt exponent overflow");
                    let exponent = (twice_triangular / 2)
                        .checked_mul(m)
                        .expect("Hoggatt exponent overflow");
                    count * pow_usize(&q, exponent)
                })
                .collect()
        })
        .collect()
}

fn pow_usize(base: &BigInt, mut exponent: usize) -> BigInt {
    let mut result = BigInt::one();
    let mut factor = base.clone();
    while exponent != 0 {
        if exponent & 1 != 0 {
            result *= &factor;
        }
        exponent >>= 1;
        if exponent != 0 {
            factor = &factor * &factor;
        }
    }
    result
}

fn add_linear_multiple(out: &mut [BigInt], row: &[BigInt], constant: i32, linear: i32) {
    for (k, value) in row.iter().enumerate() {
        out[k] += value * constant;
        out[k + 1] += value * linear;
    }
}

fn square_triangle(triangle: &[Vec<BigInt>]) -> Vec<Vec<BigInt>> {
    triangle
        .iter()
        .enumerate()
        .map(|(n, row)| {
            let mut out = vec![BigInt::zero(); n + 1];
            for (j, weight) in row.iter().enumerate() {
                for (k, value) in triangle[j].iter().enumerate() {
                    out[k] += weight * value;
                }
            }
            out
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn big(row: &[i64]) -> Vec<BigInt> {
        row.iter().map(|&x| BigInt::from(x)).collect()
    }

    #[test]
    fn delannoy_square_recurrence_matches_matrix_product() {
        let d = delannoy_polynomials_bigint(25);
        assert_eq!(d[3], big(&[1, 5, 5, 1]));
        assert_eq!(delannoy_square_polynomials_bigint(25), square_triangle(&d));
    }

    #[test]
    fn eulerian_square_retains_paper_indexing() {
        let a = eulerian_square_polynomials_bigint(4);
        assert_eq!(a[0], big(&[1]));
        assert_eq!(a[1], big(&[0, 1]));
        assert_eq!(a[2], big(&[0, 2, 1]));
        assert_eq!(a[3], big(&[0, 6, 8, 1]));
        assert_eq!(a[4], big(&[0, 24, 66, 22, 1]));
    }

    #[test]
    fn hoggatt_at_one_matches_narayana_and_baxter() {
        assert_eq!(
            hoggatt_polynomials_bigint(8, 2, 1),
            super::super::narayana_polynomials_bigint(8)
        );
        let h = hoggatt_polynomials_bigint(5, 3, 1);
        assert_eq!(h[0], big(&[1]));
        assert_eq!(h[2], big(&[1, 4, 1]));
        assert_eq!(h[3], big(&[1, 10, 10, 1]));
        assert_eq!(h[4], big(&[1, 20, 50, 20, 1]));
    }

    #[test]
    fn hoggatt_height_one_matches_q_binomial_product() {
        for q in [1, 2, 3, 10] {
            let h = hoggatt_polynomials_bigint(9, 1, q);
            let mut product = vec![BigInt::one()];
            for n in 1..=9 {
                assert_eq!(h[n - 1], product);
                let factor = pow_usize(&BigInt::from(q), n);
                let mut next = vec![BigInt::zero(); product.len() + 1];
                for (i, c) in product.iter().enumerate() {
                    next[i] += c;
                    next[i + 1] += c * &factor;
                }
                product = next;
            }
        }
    }

    #[test]
    fn hoggatt_matches_small_plane_partition_enumeration() {
        // Exhaust all 2x2 arrays with weakly decreasing rows/columns, entries 0..m.
        for m in 1..=3 {
            for q in [1, 2, 3] {
                let mut weight = BigInt::zero();
                for a in 0..=m {
                    for b in 0..=a {
                        for c in 0..=a {
                            for d in 0..=b.min(c) {
                                weight += pow_usize(&BigInt::from(q), a + b + c + d);
                            }
                        }
                    }
                }
                assert_eq!(
                    hoggatt_polynomials_bigint(5, m, q)[4][2],
                    weight * pow_usize(&BigInt::from(q), 3 * m)
                );
            }
        }
    }

    #[test]
    fn empty_family_and_zero_degree_boundaries() {
        assert!(hoggatt_polynomials_bigint(0, 2, 2).is_empty());
        assert_eq!(hoggatt_polynomials_bigint(1, 3, 2), vec![big(&[1])]);
        assert_eq!(delannoy_square_polynomials_bigint(0), vec![big(&[1])]);
        assert_eq!(eulerian_square_polynomials_bigint(0), vec![big(&[1])]);
    }

    #[test]
    #[should_panic(expected = "must be positive")]
    fn zero_q_is_rejected() {
        hoggatt_polynomials_bigint(2, 2, 0);
    }
}
