//! Bézoutians with orientation `(f(x)g(y)-f(y)g(x))/(x-y)`.
//!
//! A quadratic coefficient recurrence builds the matrix in O(n²) operations.

use num_bigint::BigInt;
use num_traits::Zero;

/// Build the general integer Bézoutian from ascending coefficient vectors.
///
/// Its order is `max(deg f, deg g)`, treating constants and zero as degree zero.
/// The matrix is symmetric; swapping arguments negates it. A zero argument
/// gives a zero matrix of the appropriate order.
pub fn bezoutian_bigint_coeffs(f: &[BigInt], g: &[BigInt]) -> Vec<Vec<BigInt>> {
    let n = f
        .iter()
        .rposition(|x| !x.is_zero())
        .unwrap_or(0)
        .max(g.iter().rposition(|x| !x.is_zero()).unwrap_or(0));
    let mut b = vec![vec![BigInt::zero(); n]; n];
    let zero = BigInt::zero();
    let coefficient = |p: &[BigInt], k: usize| p.get(k).unwrap_or(&zero).clone();
    // N[i,j] = B[i-1,j] - B[i,j-1]. Start from the right edge,
    // B[i-1,j] = N[i,j] + B[i,j-1], with B[n,*] = 0.
    for i in (1..=n).rev() {
        for j in 0..n {
            b[i - 1][j] =
                coefficient(f, i) * coefficient(g, j) - coefficient(g, i) * coefficient(f, j);
            if i < n && j > 0 {
                let next = b[i][j - 1].clone();
                b[i - 1][j] += next;
            }
        }
    }
    b
}
