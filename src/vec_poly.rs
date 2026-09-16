//! Lightweight polynomial arithmetic on `&[i64]` coefficient vectors.
//!
//! These are convenience functions for the common pattern in research experiments
//! where polynomials are represented as `Vec<i64>` in ascending degree order
//! (`coeffs[i]` = coefficient of t^i). For generic coefficient types, use
//! [`Polynomial<C>`](crate::polynomial::Polynomial) instead.

/// Trim trailing zeros from a coefficient vector.
pub fn trim(p: &[i64]) -> Vec<i64> {
    let mut v = p.to_vec();
    while v.len() > 1 && v.last() == Some(&0) {
        v.pop();
    }
    if v.is_empty() {
        vec![0]
    } else {
        v
    }
}

/// Whether the polynomial is zero.
pub fn is_zero(p: &[i64]) -> bool {
    p.iter().all(|&c| c == 0)
}

/// Degree of the polynomial, or `None` if zero.
pub fn degree(p: &[i64]) -> Option<usize> {
    p.iter().rposition(|&c| c != 0)
}

/// Add two polynomials.
pub fn add(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] += v;
    }
    trim(&r)
}

/// Subtract: a - b.
pub fn sub(a: &[i64], b: &[i64]) -> Vec<i64> {
    let len = a.len().max(b.len());
    let mut r = vec![0i64; len];
    for (i, &v) in a.iter().enumerate() {
        r[i] += v;
    }
    for (i, &v) in b.iter().enumerate() {
        r[i] -= v;
    }
    trim(&r)
}

/// Multiply two polynomials.
pub fn mul(a: &[i64], b: &[i64]) -> Vec<i64> {
    if is_zero(a) || is_zero(b) {
        return vec![0];
    }
    let mut r = vec![0i64; a.len() + b.len() - 1];
    for (i, &av) in a.iter().enumerate() {
        if av == 0 {
            continue;
        }
        for (j, &bv) in b.iter().enumerate() {
            r[i + j] += av * bv;
        }
    }
    trim(&r)
}

/// Multiply by t (shift coefficients up by one degree).
pub fn shift(p: &[i64]) -> Vec<i64> {
    if is_zero(p) {
        return vec![0];
    }
    let mut r = vec![0i64; p.len() + 1];
    for (i, &v) in p.iter().enumerate() {
        r[i + 1] = v;
    }
    r
}

/// Get the coefficient of t^k, or 0 if out of range.
pub fn coeff(p: &[i64], k: usize) -> i64 {
    if k < p.len() {
        p[k]
    } else {
        0
    }
}

/// Evaluate polynomial at an integer point.
pub fn evaluate(p: &[i64], x: i64) -> i64 {
    let mut result = 0i64;
    for &c in p.iter().rev() {
        result = result * x + c;
    }
    result
}

/// Multiply all coefficients by a scalar.
pub fn scale(p: &[i64], c: i64) -> Vec<i64> {
    if c == 0 {
        return vec![0];
    }
    trim(&p.iter().map(|&v| v * c).collect::<Vec<_>>())
}

/// Negate a polynomial.
pub fn neg(p: &[i64]) -> Vec<i64> {
    p.iter().map(|&v| -v).collect()
}

/// Permanent of a matrix whose entries are polynomials (stored as `Vec<i64>`).
///
/// The matrix is `mat[row][col]`, where each entry is a polynomial. The
/// permanent is defined for square matrices; a nonempty matrix must therefore
/// have exactly as many columns as rows, and every row must have the same
/// length. A rectangular or ragged matrix causes a panic. The empty matrix has
/// permanent one, represented by the constant polynomial `[1]`.
///
/// The computation uses the standard subset dynamic program: after processing
/// row `r`, the entry for a column subset records the sum of products obtained
/// by assigning rows `0..=r` to those columns. This avoids the factorial work
/// of a Laplace expansion while retaining exact integer arithmetic.
pub fn permanent(mat: &[Vec<Vec<i64>>]) -> Vec<i64> {
    let n = mat.len();
    if n == 0 {
        return vec![1];
    }

    let columns = mat[0].len();
    assert!(
        mat.iter().all(|row| row.len() == columns),
        "permanent requires a rectangular matrix"
    );
    assert_eq!(
        n, columns,
        "permanent requires a square matrix (got {n} rows and {columns} columns)"
    );
    assert!(
        n < usize::BITS as usize,
        "permanent matrix is too large for subset dynamic programming"
    );

    let subset_count = 1usize << n;
    let full_mask = subset_count - 1;
    let mut current: Vec<Option<Vec<i64>>> = vec![None; subset_count];
    current[0] = Some(vec![1]);

    for row in mat {
        let mut next: Vec<Option<Vec<i64>>> = vec![None; subset_count];
        for (used_columns, partial) in current.into_iter().enumerate() {
            let Some(partial) = partial else {
                continue;
            };
            for (column, entry) in row.iter().enumerate() {
                let bit = 1usize << column;
                if used_columns & bit != 0 {
                    continue;
                }
                let mask = used_columns | bit;
                let term = mul(entry, &partial);
                if let Some(existing) = next[mask].as_mut() {
                    *existing = add(existing, &term);
                } else {
                    next[mask] = Some(term);
                }
            }
        }
        current = next;
    }

    current[full_mask].take().unwrap_or_else(|| vec![0])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim() {
        assert_eq!(trim(&[1, 2, 0, 0]), vec![1, 2]);
        assert_eq!(trim(&[0, 0, 0]), vec![0]);
        assert_eq!(trim(&[3]), vec![3]);
    }

    #[test]
    fn test_is_zero() {
        assert!(is_zero(&[0]));
        assert!(is_zero(&[0, 0, 0]));
        assert!(!is_zero(&[1]));
        assert!(!is_zero(&[0, 1]));
    }

    #[test]
    fn test_degree() {
        assert_eq!(degree(&[0]), None);
        assert_eq!(degree(&[3]), Some(0));
        assert_eq!(degree(&[1, 2, 3]), Some(2));
        assert_eq!(degree(&[1, 0, 0]), Some(0));
    }

    #[test]
    fn test_add() {
        assert_eq!(add(&[1, 2], &[3, 4, 5]), vec![4, 6, 5]);
        assert_eq!(add(&[1, -1], &[-1, 1]), vec![0]);
    }

    #[test]
    fn test_sub() {
        assert_eq!(sub(&[3, 4, 5], &[1, 2]), vec![2, 2, 5]);
        assert_eq!(sub(&[1, 2], &[1, 2]), vec![0]);
    }

    #[test]
    fn test_mul() {
        // (1 + t)(1 + t) = 1 + 2t + t^2
        assert_eq!(mul(&[1, 1], &[1, 1]), vec![1, 2, 1]);
        // (1 + t)(1 - t) = 1 - t^2
        assert_eq!(mul(&[1, 1], &[1, -1]), vec![1, 0, -1]);
        assert_eq!(mul(&[0], &[1, 2, 3]), vec![0]);
    }

    #[test]
    fn test_shift() {
        assert_eq!(shift(&[1, 2, 3]), vec![0, 1, 2, 3]);
        assert_eq!(shift(&[0]), vec![0]);
    }

    #[test]
    fn test_coeff() {
        assert_eq!(coeff(&[10, 20, 30], 0), 10);
        assert_eq!(coeff(&[10, 20, 30], 2), 30);
        assert_eq!(coeff(&[10, 20, 30], 5), 0);
    }

    #[test]
    fn test_evaluate() {
        // p(t) = 1 + 2t + 3t^2, p(2) = 1 + 4 + 12 = 17
        assert_eq!(evaluate(&[1, 2, 3], 2), 17);
        assert_eq!(evaluate(&[5], 100), 5);
        assert_eq!(evaluate(&[0, 0, 1], 3), 9);
    }

    #[test]
    fn test_scale() {
        assert_eq!(scale(&[1, 2, 3], 2), vec![2, 4, 6]);
        assert_eq!(scale(&[1, 2, 3], 0), vec![0]);
        assert_eq!(scale(&[1, 2, 3], -1), vec![-1, -2, -3]);
    }

    #[test]
    fn test_permanent_1x1() {
        let mat = vec![vec![vec![0, 1]]]; // entry = t
        assert_eq!(permanent(&mat), vec![0, 1]);
    }

    #[test]
    fn test_permanent_2x2_identity() {
        // [[1, 0], [0, 1]] -> perm = 1*1 + 0*0 = 1
        let mat = vec![vec![vec![1], vec![0]], vec![vec![0], vec![1]]];
        assert_eq!(permanent(&mat), vec![1]);
    }

    #[test]
    fn test_permanent_2x2_all_ones() {
        // [[1, 1], [1, 1]] -> perm = 1*1 + 1*1 = 2
        let mat = vec![vec![vec![1], vec![1]], vec![vec![1], vec![1]]];
        assert_eq!(permanent(&mat), vec![2]);
    }

    #[test]
    fn test_permanent_2x2_with_t() {
        // [[1, t], [t, 1]] -> perm = 1*1 + t*t = 1 + t^2
        let mat = vec![vec![vec![1], vec![0, 1]], vec![vec![0, 1], vec![1]]];
        assert_eq!(permanent(&mat), vec![1, 0, 1]);
    }

    #[test]
    fn test_permanent_empty() {
        let mat: Vec<Vec<Vec<i64>>> = vec![];
        assert_eq!(permanent(&mat), vec![1]);
    }

    fn naive_permanent(mat: &[Vec<Vec<i64>>]) -> Vec<i64> {
        if mat.is_empty() {
            return vec![1];
        }
        let columns = mat[0].len();
        let mut result = vec![0];
        for column in 0..columns {
            let sub: Vec<Vec<Vec<i64>>> = mat[1..]
                .iter()
                .map(|row| {
                    row.iter()
                        .enumerate()
                        .filter_map(|(index, entry)| (index != column).then_some(entry.clone()))
                        .collect()
                })
                .collect();
            result = add(&result, &mul(&mat[0][column], &naive_permanent(&sub)));
        }
        result
    }

    #[test]
    fn test_permanent_matches_small_naive_oracle() {
        for size in 0..=4 {
            let mat: Vec<Vec<Vec<i64>>> = (0..size)
                .map(|row| {
                    (0..size)
                        .map(|column| {
                            let constant = ((3 * row + 5 * column + 1) % 7) as i64 - 3;
                            let linear = ((row + 2 * column) % 3) as i64 - 1;
                            vec![constant, linear]
                        })
                        .collect()
                })
                .collect();
            assert_eq!(permanent(&mat), naive_permanent(&mat), "size {size}");
        }
    }

    #[test]
    #[should_panic(expected = "permanent requires a square matrix")]
    fn test_permanent_rejects_rectangular_matrix() {
        let mat = vec![
            vec![vec![1], vec![2]],
            vec![vec![3], vec![4]],
            vec![vec![5], vec![6]],
        ];
        let _ = permanent(&mat);
    }

    #[test]
    #[should_panic(expected = "permanent requires a rectangular matrix")]
    fn test_permanent_rejects_ragged_matrix() {
        let mat = vec![vec![vec![1]], vec![vec![2], vec![3]]];
        let _ = permanent(&mat);
    }
}
