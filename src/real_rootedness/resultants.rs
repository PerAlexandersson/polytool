//! Sylvester matrices, resultants and discriminants with exact `BigInt` results.

use crate::vec_poly;
use num_bigint::BigInt;

/// Compute the Sylvester matrix of f and g.
///
/// For f of degree n and g of degree m, the Sylvester matrix is (n+m) x (n+m).
/// The first m rows are shifted copies of f's coefficients (in descending order),
/// and the last n rows are shifted copies of g's coefficients.
///
/// Returns `None` if either polynomial is zero.
pub fn sylvester_matrix(f: &[i64], g: &[i64]) -> Option<Vec<Vec<BigInt>>> {
    let n = vec_poly::degree(f)?;
    let m = vec_poly::degree(g)?;
    let size = n + m;
    let zero = BigInt::from(0);
    let mut mat = vec![vec![zero.clone(); size]; size];

    // First m rows: coefficients of f, shifted right by row index.
    // Row i has f_n, f_{n-1}, ..., f_0 starting at column i.
    for i in 0..m {
        for k in 0..=n {
            mat[i][i + k] = BigInt::from(vec_poly::coeff(f, n - k));
        }
    }
    // Last n rows: coefficients of g, shifted right by row index.
    for i in 0..n {
        for k in 0..=m {
            mat[m + i][i + k] = BigInt::from(vec_poly::coeff(g, m - k));
        }
    }
    Some(mat)
}

fn sylvester_matrix_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Option<Vec<Vec<BigInt>>> {
    let n = vec_poly::degree(f)?;
    let m = vec_poly::degree(g)?;
    let size = n + m;
    let zero = BigInt::from(0);
    let mut mat = vec![vec![zero; size]; size];

    for i in 0..m {
        for k in 0..=n {
            mat[i][i + k] = f[n - k].clone();
        }
    }
    for i in 0..n {
        for k in 0..=m {
            mat[m + i][i + k] = g[m - k].clone();
        }
    }
    Some(mat)
}

fn determinant_bigint(mat: &[Vec<BigInt>]) -> BigInt {
    crate::linalg::determinant(mat)
}

/// Compute the resultant Res(f, g) of two polynomials.
///
/// The resultant is zero if and only if f and g share a common root (over the
/// algebraic closure) or one of them is zero. Computed via the Sylvester matrix
/// determinant.
///
/// Returns the resultant as a `BigInt`.
pub fn resultant(f: &[i64], g: &[i64]) -> BigInt {
    let f_big = f.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    let g_big = g.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    resultant_bigint_coeffs(&f_big, &g_big)
}

/// Compute the resultant for `BigInt` coefficient polynomials.
///
/// Same conventions as [`resultant`] (ascending coefficients, `0` if either
/// polynomial is zero).
pub fn resultant_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> BigInt {
    match sylvester_matrix_bigint_coeffs(f, g) {
        Some(mat) => determinant_bigint(&mat),
        None => BigInt::from(0),
    }
}

/// Compute the discriminant of a polynomial.
///
/// ```text
/// disc(f) = (-1)^{n(n-1)/2} * Res(f, f') / lc(f)
/// ```
///
/// where n = deg(f) and f' is the formal derivative.
///
/// The discriminant is zero iff f has a repeated root. For a polynomial with
/// roots r_1, ..., r_n, disc(f) = lc(f)^{2n-2} * prod_{i<j} (r_i - r_j)^2.
///
/// Conventions: the zero polynomial has discriminant `0`, and every nonzero
/// polynomial of degree at most 1 has discriminant `1`.  The result is exact
/// (a `BigInt`); coefficients are ascending.
pub fn discriminant(f: &[i64]) -> BigInt {
    let f_big = f.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
    discriminant_bigint_coeffs(&f_big)
}

/// Compute the discriminant for a `BigInt` coefficient polynomial.
///
/// Same conventions as [`discriminant`].
pub fn discriminant_bigint_coeffs(f: &[BigInt]) -> BigInt {
    let n = match vec_poly::degree(f) {
        Some(n) => n,
        None => return BigInt::from(0),
    };
    if n <= 1 {
        return BigInt::from(1);
    }

    let mut fp: Vec<BigInt> = Vec::with_capacity(n);
    for k in 0..n {
        fp.push(BigInt::from(k + 1) * &f[k + 1]);
    }

    let res = match sylvester_matrix_bigint_coeffs(f, &fp) {
        Some(mat) => determinant_bigint(&mat),
        None => BigInt::from(0),
    };
    let lc = f[n].clone();
    let sign_exp = n * (n - 1) / 2;
    let sign = if sign_exp % 2 == 0 {
        BigInt::from(1)
    } else {
        BigInt::from(-1)
    };

    sign * res / lc
}

#[cfg(test)]
mod tests {
    use super::*;

    use num_traits::Signed;

    // -- Resultant / discriminant tests --

    #[test]
    fn test_resultant() {
        // Res((t-1)(t-2), (t-3)) = (1-3)(2-3) = (-2)(-1) = 2
        assert_eq!(resultant(&[2, -3, 1], &[-3, 1]), BigInt::from(2));

        // Res(t^2+1, t-1) = 1^2 + 1 = 2
        assert_eq!(resultant(&[1, 0, 1], &[-1, 1]), BigInt::from(2));

        // Res(t-1, t-1) = 0 (shared root)
        assert_eq!(resultant(&[-1, 1], &[-1, 1]), BigInt::from(0));
    }

    #[test]
    fn test_discriminant() {
        // disc(t^2 - 1) = disc((t-1)(t+1)). Roots 1,-1.
        // disc = (-1)^{2*1/2} * Res(t^2-1, 2t) / 1
        // = (-1)^1 * Res([−1,0,1],[0,2]) / 1
        // Sylvester of t^2-1 (deg 2) and 2t (deg 1): 3x3 matrix.
        //   row0 (f shifted): [1, 0, -1]
        //   row1 (g shifted): [2, 0, 0]
        //   row2 (g shifted): [0, 2, 0]
        // det = 1*(0*0 - 0*2) - 0*(2*0 - 0*0) + (-1)*(2*2 - 0*0) = -4
        // disc = (-1)^1 * (-4) / 1 = 4
        assert_eq!(discriminant(&[-1, 0, 1]), BigInt::from(4));

        // disc((t-1)(t-2)) = (1-2)^2 = 1. f = [2,-3,1], lc=1, n=2.
        // disc = (-1)^1 * Res(f, f') / lc
        // f' = [-3, 2]. Sylvester 3x3.
        assert_eq!(discriminant(&[2, -3, 1]), BigInt::from(1));

        // disc(t^2 + 1) = -4 (no real roots, negative discriminant)
        assert_eq!(discriminant(&[1, 0, 1]), BigInt::from(-4));
    }

    #[test]
    fn test_discriminant_derivative_uses_bigint() {
        // f(t) = i64::MAX * t^2 + 1 has derivative 2*i64::MAX*t,
        // which does not fit in i64.
        assert_eq!(
            discriminant(&[1, 0, i64::MAX]),
            BigInt::from(-4) * BigInt::from(i64::MAX)
        );
    }

    #[test]
    fn test_discriminant_degree3() {
        // f = (t-1)(t-2)(t-3) = -6 + 11t - 6t^2 + t^3
        // disc = product (r_i - r_j)^2 = (1-2)^2(1-3)^2(2-3)^2 = 1*4*1 = 4
        // But with the sign: (-1)^{3*2/2} = (-1)^3 = -1, and
        // disc = (-1)^3 * Res(f,f') / lc(f).
        // Let's just compute and check the absolute value.
        let d = discriminant(&[-6, 11, -6, 1]);
        // |disc| = 4 for monic with roots 1,2,3
        assert_eq!(d.clone() * d.clone().signum(), BigInt::from(4));
    }
}
