//! Polynomial arithmetic on plain coefficient vectors `&[C]`.
//!
//! These functions serve the common research pattern where a polynomial is a
//! `Vec<C>` in ascending degree order: `p[i]` is the coefficient of `t^i`.
//! They are generic over [`CoeffRing`], so the same code works for `i64`,
//! `i128`, `BigInt`, `BigRational` (`Ratio<BigInt>`) and `Ratio<i64>`.  For a
//! normalizing wrapper type with the same operations, use
//! [`Polynomial<C>`](crate::polynomial::Polynomial).
//!
//! # Conventions
//!
//! - **Inputs** may carry trailing zeros, and the empty slice is accepted as
//!   the zero polynomial.
//! - **Outputs** of the arithmetic functions are *canonical*: trailing zeros
//!   are removed, and the zero polynomial is returned as the one-element
//!   vector `[0]` (not as an empty vector).  [`trim`] produces this form;
//!   [`strip_trailing_zeros`] instead produces the empty-vector form used by
//!   [`Polynomial::coeffs`](crate::polynomial::Polynomial::coeffs).
//! - [`degree`] of the zero polynomial is `None`.
//! - Several primitive integer types implement [`CoeffRing`], so a call whose
//!   arguments are *only* integer literals needs a type annotation, for
//!   example `vec_poly::mul(&[1i64, 1], &[1, -1])`.
//!
//! # Overflow
//!
//! All arithmetic goes through the borrowing hooks of [`CoeffRing`].  For the
//! fixed-width coefficient types `i64`, `i128` and `Ratio<i64>` these hooks are
//! checked, so every function here **panics on overflow** (in debug and
//! release builds alike) rather than silently wrapping.  Use `BigInt`
//! coefficients when values may exceed the machine range; for example,
//! convert with `p.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>()`.
//!
//! ```
//! use num_bigint::BigInt;
//! use polytool::vec_poly;
//!
//! // (1 + t)^3 = 1 + 3t + 3t^2 + t^3
//! assert_eq!(vec_poly::pow(&[1i64, 1], 3), vec![1, 3, 3, 1]);
//!
//! // The same code works for BigInt coefficients.
//! let p: Vec<BigInt> = vec![BigInt::from(-1), BigInt::from(0), BigInt::from(1)];
//! let dp = vec_poly::derivative(&p);
//! assert_eq!(dp, vec![BigInt::from(0), BigInt::from(2)]);
//!
//! // Zero is canonically `[0]`.
//! assert_eq!(vec_poly::sub(&[1i64, 2], &[1, 2]), vec![0]);
//! ```

use crate::polynomial::CoeffRing;

// ---------------------------------------------------------------------------
// Normalization and inspection
// ---------------------------------------------------------------------------

/// Return the canonical form of `p`: trailing zeros removed, with the zero
/// polynomial (including the empty slice) represented as `[0]`.
pub fn trim<C: CoeffRing>(p: &[C]) -> Vec<C> {
    match degree(p) {
        Some(d) => p[..=d].to_vec(),
        None => vec![C::zero()],
    }
}

/// Remove all trailing zeros in place, so the zero polynomial becomes the
/// **empty** vector.
///
/// This is the representation used by
/// [`Polynomial::coeffs`](crate::polynomial::Polynomial::coeffs) and by code
/// that tests `p.is_empty()` for zero.  Use [`trim`] for the canonical `[0]`
/// form returned by the arithmetic in this module.
pub fn strip_trailing_zeros<C: CoeffRing>(p: &mut Vec<C>) {
    let len = degree(p).map_or(0, |d| d + 1);
    p.truncate(len);
}

/// Canonicalize an owned vector without copying it.
fn canonical<C: CoeffRing>(mut p: Vec<C>) -> Vec<C> {
    strip_trailing_zeros(&mut p);
    if p.is_empty() {
        p.push(C::zero());
    }
    p
}

/// Whether the polynomial is zero (all coefficients zero, or empty).
pub fn is_zero<C: CoeffRing>(p: &[C]) -> bool {
    p.iter().all(C::is_zero)
}

/// Degree of the polynomial, or `None` for the zero polynomial.
pub fn degree<C: CoeffRing>(p: &[C]) -> Option<usize> {
    p.iter().rposition(|c| !c.is_zero())
}

/// The coefficient of `t^k`, or zero if `k` is past the end of `p`.
pub fn coeff<C: CoeffRing>(p: &[C], k: usize) -> C {
    p.get(k).cloned().unwrap_or_else(C::zero)
}

// ---------------------------------------------------------------------------
// Ring operations
// ---------------------------------------------------------------------------

/// The sum `a + b`.
pub fn add<C: CoeffRing>(a: &[C], b: &[C]) -> Vec<C> {
    let (long, short) = if a.len() >= b.len() { (a, b) } else { (b, a) };
    let mut r = long.to_vec();
    for (x, y) in r.iter_mut().zip(short) {
        x.add_assign_ref(y);
    }
    canonical(r)
}

/// The difference `a - b`.
pub fn sub<C: CoeffRing>(a: &[C], b: &[C]) -> Vec<C> {
    let mut r = a.to_vec();
    if r.len() < b.len() {
        r.resize(b.len(), C::zero());
    }
    for (x, y) in r.iter_mut().zip(b) {
        x.sub_assign_ref(y);
    }
    canonical(r)
}

/// The negation `-p`.
pub fn neg<C: CoeffRing>(p: &[C]) -> Vec<C> {
    canonical(p.iter().map(C::neg_ref).collect())
}

/// The product `a * b` (schoolbook multiplication).
pub fn mul<C: CoeffRing>(a: &[C], b: &[C]) -> Vec<C> {
    let (Some(da), Some(db)) = (degree(a), degree(b)) else {
        return vec![C::zero()];
    };
    let (a, b) = (&a[..=da], &b[..=db]);
    let mut r = vec![C::zero(); da + db + 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            if !y.is_zero() {
                r[i + j].add_mul_assign(x, y);
            }
        }
    }
    canonical(r)
}

/// The scalar multiple `c * p`.
pub fn scale<C: CoeffRing>(p: &[C], c: &C) -> Vec<C> {
    if c.is_zero() {
        return vec![C::zero()];
    }
    canonical(p.iter().map(|x| x.mul_ref(c)).collect())
}

/// Multiply by `t`, shifting every coefficient up by one degree.
///
/// This is the coefficient shift `p(t) ↦ t·p(t)`, not the argument
/// translation `p(t) ↦ p(t + a)` of
/// [`Polynomial::shift`](crate::polynomial::Polynomial::shift).
pub fn shift<C: CoeffRing>(p: &[C]) -> Vec<C> {
    shift_by(p, 1)
}

/// Multiply by `t^k`, shifting every coefficient up by `k` degrees.
pub fn shift_by<C: CoeffRing>(p: &[C], k: usize) -> Vec<C> {
    let Some(d) = degree(p) else {
        return vec![C::zero()];
    };
    let mut r = vec![C::zero(); k];
    r.extend_from_slice(&p[..=d]);
    r
}

/// The power `p^e`, by repeated squaring.  `p^0 = 1`, also for `p = 0`.
pub fn pow<C: CoeffRing>(p: &[C], e: usize) -> Vec<C> {
    let mut result = vec![C::one()];
    if e == 0 {
        return result;
    }
    let mut base = trim(p);
    let mut e = e;
    loop {
        if e & 1 == 1 {
            result = mul(&result, &base);
        }
        e >>= 1;
        if e == 0 {
            return result;
        }
        base = mul(&base, &base);
    }
}

/// Evaluate `p(x)` by Horner's rule.
pub fn evaluate<C: CoeffRing>(p: &[C], x: &C) -> C {
    let Some(d) = degree(p) else {
        return C::zero();
    };
    let mut result = p[d].clone();
    for c in p[..d].iter().rev() {
        result = result.mul_ref(x);
        result.add_assign_ref(c);
    }
    result
}

/// The composition `p(q(t))`, by Horner's rule in the polynomial ring.
pub fn compose<C: CoeffRing>(p: &[C], q: &[C]) -> Vec<C> {
    let Some(d) = degree(p) else {
        return vec![C::zero()];
    };
    let mut result = vec![p[d].clone()];
    for c in p[..d].iter().rev() {
        result = mul(&result, q);
        result[0].add_assign_ref(c);
        result = canonical(result);
    }
    result
}

// ---------------------------------------------------------------------------
// Derivatives and reversal
// ---------------------------------------------------------------------------

/// The formal derivative `p'(t)`.
pub fn derivative<C: CoeffRing>(p: &[C]) -> Vec<C> {
    let Some(d) = degree(p) else {
        return vec![C::zero()];
    };
    if d == 0 {
        return vec![C::zero()];
    }
    canonical(
        (1..=d)
            .map(|k| p[k].mul_ref(&C::from_i64(index_to_i64(k))))
            .collect(),
    )
}

/// The Wronskian `W(f, g) = f'·g - f·g'`.
///
/// Sign convention: if `f` and `g` have positive leading coefficients and `f`
/// interlaces `g` in the directed sense of [`crate::check_weak_interlacing`]
/// (`f ≪ g`), then `W(f, g) ≤ 0` on all of ℝ.
pub fn wronskian<C: CoeffRing>(f: &[C], g: &[C]) -> Vec<C> {
    sub(&mul(&derivative(f), g), &mul(f, &derivative(g)))
}

/// The reversal of `p` with respect to the nominal degree `n`:
/// `t^n · p(1/t)`.
///
/// Coefficient `i` of the result is `p[n - i]` for `0 ≤ i ≤ n`, so the input
/// is read as the length-`n + 1` vector `(p_0, ..., p_n)` and reversed.  The
/// result is canonical: when `p_0 = 0` its degree is below `n`.  Trailing
/// zeros of the input beyond index `n` are allowed.
///
/// Returns `None` if `deg p > n`.  The zero polynomial reverses to zero for
/// every `n`.  See [`reciprocal`] for reversal at the actual degree.
pub fn reverse<C: CoeffRing>(p: &[C], n: usize) -> Option<Vec<C>> {
    match degree(p) {
        None => Some(vec![C::zero()]),
        Some(d) if d > n => None,
        Some(_) => Some(canonical((0..=n).rev().map(|i| coeff(p, i)).collect())),
    }
}

/// The reciprocal polynomial `t^d · p(1/t)` with `d = deg p`.
///
/// Equivalently, the reversed coefficient vector with trailing zeros
/// removed.  The zero polynomial maps to `[0]`.
pub fn reciprocal<C: CoeffRing>(p: &[C]) -> Vec<C> {
    let d = degree(p).unwrap_or(0);
    reverse(p, d).expect("degree is at most the actual degree")
}

fn index_to_i64(k: usize) -> i64 {
    i64::try_from(k).expect("polynomial degree exceeds i64::MAX")
}

// ---------------------------------------------------------------------------
// Permanent
// ---------------------------------------------------------------------------

/// Permanent of a matrix whose entries are polynomials (coefficient vectors).
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
/// of a Laplace expansion while retaining exact arithmetic (fixed-width
/// coefficient types panic on overflow, see the module documentation).
pub fn permanent<C: CoeffRing>(mat: &[Vec<Vec<C>>]) -> Vec<C> {
    let n = mat.len();
    if n == 0 {
        return vec![C::one()];
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
    let mut current: Vec<Option<Vec<C>>> = vec![None; subset_count];
    current[0] = Some(vec![C::one()]);

    for row in mat {
        let mut next: Vec<Option<Vec<C>>> = vec![None; subset_count];
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

    current[full_mask].take().unwrap_or_else(|| vec![C::zero()])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_trim() {
        assert_eq!(trim(&[1i64, 2, 0, 0]), vec![1, 2]);
        assert_eq!(trim(&[0i64, 0, 0]), vec![0]);
        assert_eq!(trim(&[3i64]), vec![3]);
        assert_eq!(trim::<i64>(&[]), vec![0]);
    }

    #[test]
    fn test_strip_trailing_zeros() {
        let mut p = vec![1i64, 0, 2, 0, 0];
        strip_trailing_zeros(&mut p);
        assert_eq!(p, vec![1, 0, 2]);
        let mut z = vec![0i64, 0];
        strip_trailing_zeros(&mut z);
        assert!(z.is_empty());
    }

    #[test]
    fn test_is_zero() {
        assert!(is_zero(&[0i64]));
        assert!(is_zero(&[0i64, 0, 0]));
        assert!(!is_zero(&[1i64]));
        assert!(!is_zero(&[0i64, 1]));
    }

    #[test]
    fn test_degree() {
        assert_eq!(degree(&[0i64]), None);
        assert_eq!(degree(&[3i64]), Some(0));
        assert_eq!(degree(&[1i64, 2, 3]), Some(2));
        assert_eq!(degree(&[1i64, 0, 0]), Some(0));
    }

    #[test]
    fn test_add() {
        assert_eq!(add(&[1i64, 2], &[3, 4, 5]), vec![4, 6, 5]);
        assert_eq!(add(&[1i64, -1], &[-1, 1]), vec![0]);
    }

    #[test]
    fn test_sub() {
        assert_eq!(sub(&[3i64, 4, 5], &[1, 2]), vec![2, 2, 5]);
        assert_eq!(sub(&[1i64, 2], &[1, 2]), vec![0]);
    }

    #[test]
    fn test_mul() {
        // (1 + t)(1 + t) = 1 + 2t + t^2
        assert_eq!(mul(&[1i64, 1], &[1, 1]), vec![1, 2, 1]);
        // (1 + t)(1 - t) = 1 - t^2
        assert_eq!(mul(&[1i64, 1], &[1, -1]), vec![1, 0, -1]);
        assert_eq!(mul(&[0i64], &[1, 2, 3]), vec![0]);
    }

    #[test]
    fn test_shift() {
        assert_eq!(shift(&[1i64, 2, 3]), vec![0, 1, 2, 3]);
        assert_eq!(shift(&[0i64]), vec![0]);
        assert_eq!(shift(&[1i64, 2, 0]), vec![0, 1, 2]);
        assert_eq!(shift_by(&[1i64, 2], 3), vec![0, 0, 0, 1, 2]);
        assert_eq!(shift_by(&[1i64, 2], 0), vec![1, 2]);
        assert_eq!(shift_by::<i64>(&[], 2), vec![0]);
    }

    #[test]
    fn test_coeff() {
        assert_eq!(coeff(&[10i64, 20, 30], 0), 10);
        assert_eq!(coeff(&[10i64, 20, 30], 2), 30);
        assert_eq!(coeff(&[10i64, 20, 30], 5), 0);
    }

    #[test]
    fn test_evaluate() {
        // p(t) = 1 + 2t + 3t^2, p(2) = 1 + 4 + 12 = 17
        assert_eq!(evaluate(&[1i64, 2, 3], &2), 17);
        assert_eq!(evaluate(&[5i64], &100), 5);
        assert_eq!(evaluate(&[0i64, 0, 1], &3), 9);
    }

    #[test]
    fn test_scale() {
        assert_eq!(scale(&[1i64, 2, 3], &2), vec![2, 4, 6]);
        assert_eq!(scale(&[1i64, 2, 3], &0), vec![0]);
        assert_eq!(scale(&[1i64, 2, 3], &-1), vec![-1, -2, -3]);
    }

    #[test]
    fn test_permanent_1x1() {
        let mat = vec![vec![vec![0i64, 1]]]; // entry = t
        assert_eq!(permanent(&mat), vec![0, 1]);
    }

    #[test]
    fn test_permanent_2x2_identity() {
        // [[1, 0], [0, 1]] -> perm = 1*1 + 0*0 = 1
        let mat = vec![vec![vec![1i64], vec![0]], vec![vec![0], vec![1]]];
        assert_eq!(permanent(&mat), vec![1]);
    }

    #[test]
    fn test_permanent_2x2_all_ones() {
        // [[1, 1], [1, 1]] -> perm = 1*1 + 1*1 = 2
        let mat = vec![vec![vec![1i64], vec![1]], vec![vec![1], vec![1]]];
        assert_eq!(permanent(&mat), vec![2]);
    }

    #[test]
    fn test_permanent_2x2_with_t() {
        // [[1, t], [t, 1]] -> perm = 1*1 + t*t = 1 + t^2
        let mat = vec![vec![vec![1i64], vec![0, 1]], vec![vec![0, 1], vec![1]]];
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
            vec![vec![1i64], vec![2]],
            vec![vec![3], vec![4]],
            vec![vec![5], vec![6]],
        ];
        let _ = permanent(&mat);
    }

    #[test]
    #[should_panic(expected = "permanent requires a rectangular matrix")]
    fn test_permanent_rejects_ragged_matrix() {
        let mat = vec![vec![vec![1i64]], vec![vec![2], vec![3]]];
        let _ = permanent(&mat);
    }

    // -- Generic API, edge cases, and differential tests --

    use num_bigint::BigInt;
    use num_rational::Ratio;

    type Q = Ratio<BigInt>;

    /// Deterministic xorshift generator so the randomized tests are reproducible.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }

        fn poly(&mut self, max_len: usize, bound: i64) -> Vec<i64> {
            let len = (self.next() % (max_len as u64 + 1)) as usize;
            (0..len)
                .map(|_| (self.next() % (2 * bound as u64 + 1)) as i64 - bound)
                .collect()
        }
    }

    fn big(p: &[i64]) -> Vec<BigInt> {
        p.iter().map(|&c| BigInt::from(c)).collect()
    }

    fn rat(p: &[i64]) -> Vec<Q> {
        p.iter()
            .map(|&c| Q::from_integer(BigInt::from(c)))
            .collect()
    }

    /// The pre-0.3 `i64` multiplication, kept as an independent oracle.
    fn old_mul_i64(a: &[i64], b: &[i64]) -> Vec<i64> {
        if a.iter().all(|&c| c == 0) || b.iter().all(|&c| c == 0) {
            return vec![0];
        }
        let mut r = vec![0i64; a.len() + b.len() - 1];
        for (i, &av) in a.iter().enumerate() {
            for (j, &bv) in b.iter().enumerate() {
                r[i + j] += av * bv;
            }
        }
        trim(&r)
    }

    #[test]
    fn empty_slices_are_zero() {
        let e: &[i64] = &[];
        assert!(is_zero(e));
        assert_eq!(degree(e), None);
        assert_eq!(coeff(e, 3), 0);
        assert_eq!(add(e, e), vec![0]);
        assert_eq!(add(e, &[1, 2]), vec![1, 2]);
        assert_eq!(sub(e, &[1, 2]), vec![-1, -2]);
        assert_eq!(neg(e), vec![0]);
        assert_eq!(mul(e, &[1, 2]), vec![0]);
        assert_eq!(scale(e, &5), vec![0]);
        assert_eq!(evaluate(e, &7), 0);
        assert_eq!(derivative(e), vec![0]);
        assert_eq!(compose(e, &[1, 1]), vec![0]);
        assert_eq!(compose(&[1i64, 1], e), vec![1]);
        assert_eq!(pow(e, 0), vec![1]);
        assert_eq!(pow(e, 3), vec![0]);
        assert_eq!(reverse(e, 4), Some(vec![0]));
        assert_eq!(reciprocal(e), vec![0]);
        assert_eq!(wronskian(e, &[1, 2, 3]), vec![0]);
        assert_eq!(permanent::<i64>(&[vec![vec![]]]), vec![0]);
    }

    #[test]
    fn outputs_are_canonical() {
        assert_eq!(neg(&[1i64, -2, 0, 0]), vec![-1, 2]);
        assert_eq!(mul(&[1i64, 0, 0], &[2, 0]), vec![2]);
        assert_eq!(scale(&[1i64, 2, 0], &3), vec![3, 6]);
        assert_eq!(sub(&[1i64, 2, 3], &[0, 0, 3]), vec![1, 2]);
    }

    #[test]
    fn derivative_examples() {
        assert_eq!(derivative(&[5i64]), vec![0]);
        assert_eq!(derivative(&[1i64, 2, 3, 4]), vec![2, 6, 12]);
        assert_eq!(derivative(&[1i64, 2, 3, 0]), vec![2, 6]);
        assert_eq!(
            derivative(&rat(&[1, 1, 1])),
            rat(&[1, 2]),
            "BigRational derivative"
        );
    }

    #[test]
    fn reverse_and_reciprocal_conventions() {
        // t^4 p(1/t) for p = 1 + 2t + 3t^2.
        assert_eq!(reverse(&[1i64, 2, 3], 4), Some(vec![0, 0, 3, 2, 1]));
        assert_eq!(reverse(&[1i64, 2, 3], 2), Some(vec![3, 2, 1]));
        // Zero constant term lowers the degree of the reversal.
        assert_eq!(reverse(&[0i64, 1, 2], 2), Some(vec![2, 1]));
        // Degree above the nominal degree.
        assert_eq!(reverse(&[1i64, 2, 3], 1), None);
        // Trailing zeros past n are fine.
        assert_eq!(reverse(&[1i64, 2, 0, 0], 1), Some(vec![2, 1]));
        assert_eq!(reciprocal(&[1i64, 2, 3]), vec![3, 2, 1]);
        assert_eq!(reciprocal(&[0i64, 1, 2, 0]), vec![2, 1]);
        assert_eq!(reciprocal(&[7i64]), vec![7]);
    }

    #[test]
    fn pow_examples() {
        assert_eq!(pow(&[1i64, 1], 0), vec![1]);
        assert_eq!(pow(&[1i64, 1], 1), vec![1, 1]);
        assert_eq!(pow(&[1i64, 1], 4), vec![1, 4, 6, 4, 1]);
        assert_eq!(pow(&[0i64, 2], 3), vec![0, 0, 0, 8]);
        assert_eq!(pow(&[0i64], 0), vec![1]);
    }

    #[test]
    fn compose_examples() {
        // p(t) = t^2 + 1, q(t) = t + 1 -> (t + 1)^2 + 1 = 2 + 2t + t^2.
        assert_eq!(compose(&[1i64, 0, 1], &[1, 1]), vec![2, 2, 1]);
        // Constant p.
        assert_eq!(compose(&[3i64], &[1, 2, 3]), vec![3]);
        // Constant q.
        assert_eq!(compose(&[1i64, 1, 1], &[2]), vec![7]);
    }

    #[test]
    fn generic_mul_matches_old_i64_implementation() {
        let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
        for _ in 0..500 {
            let a = rng.poly(7, 20);
            let b = rng.poly(7, 20);
            let expected = old_mul_i64(&a, &b);
            assert_eq!(mul(&a, &b), expected, "{a:?} * {b:?}");
            assert_eq!(mul(&big(&a), &big(&b)), big(&expected));
            assert_eq!(mul(&rat(&a), &rat(&b)), rat(&expected));
            let a128: Vec<i128> = a.iter().map(|&c| c.into()).collect();
            let b128: Vec<i128> = b.iter().map(|&c| c.into()).collect();
            let e128: Vec<i128> = expected.iter().map(|&c| c.into()).collect();
            assert_eq!(mul(&a128, &b128), e128);
        }
    }

    #[test]
    fn ring_identities_on_random_inputs() {
        let mut rng = Rng(12345);
        for _ in 0..300 {
            let a = rng.poly(6, 9);
            let b = rng.poly(6, 9);
            let c = rng.poly(4, 9);
            assert_eq!(sub(&add(&a, &b), &b), trim(&a));
            assert_eq!(add(&a, &neg(&a)), vec![0]);
            assert_eq!(
                mul(&add(&a, &b), &c),
                add(&mul(&a, &c), &mul(&b, &c)),
                "distributivity"
            );
            // Leibniz rule.
            assert_eq!(
                derivative(&mul(&a, &b)),
                add(&mul(&derivative(&a), &b), &mul(&a, &derivative(&b)))
            );
            let n = degree(&a).unwrap_or(0) + 2;
            let r = reverse(&a, n).unwrap();
            for x in -3i64..=3 {
                // t^n a(1/t) evaluated at x != 0 equals x^n a(1/x); compare
                // over the rationals.
                if x != 0 {
                    let xq = Q::from_integer(BigInt::from(x));
                    let lhs = evaluate(&rat(&r), &xq);
                    let rhs = evaluate(&pow(&rat(&[0, 1]), n), &xq)
                        * evaluate(&rat(&a), &(Q::from_integer(BigInt::from(1)) / xq));
                    assert_eq!(lhs, rhs);
                }
            }
        }
    }

    #[test]
    fn compose_matches_evaluation_at_points() {
        let mut rng = Rng(777);
        for _ in 0..200 {
            let p = rng.poly(5, 6);
            let q = rng.poly(4, 6);
            let composed = compose(&big(&p), &big(&q));
            for x in -4i64..=4 {
                let x = BigInt::from(x);
                assert_eq!(
                    evaluate(&composed, &x),
                    evaluate(&big(&p), &evaluate(&big(&q), &x)),
                    "p = {p:?}, q = {q:?}"
                );
            }
            // Rational evaluation point.
            let x = Q::new(BigInt::from(-3), BigInt::from(2));
            assert_eq!(
                evaluate(&compose(&rat(&p), &rat(&q)), &x),
                evaluate(&rat(&p), &evaluate(&rat(&q), &x))
            );
        }
    }

    #[test]
    fn pow_matches_repeated_multiplication() {
        let mut rng = Rng(4242);
        for _ in 0..100 {
            let p = big(&rng.poly(4, 5));
            let mut expected = vec![BigInt::from(1)];
            for e in 0..7 {
                assert_eq!(pow(&p, e), expected);
                expected = mul(&expected, &p);
            }
        }
    }

    #[test]
    fn wronskian_matches_direct_formula() {
        // W(f, g) = sum_{i, j} (i - j) f_i g_j t^{i + j - 1}.
        let mut rng = Rng(99);
        for _ in 0..300 {
            let f = rng.poly(6, 9);
            let g = rng.poly(6, 9);
            let mut direct = vec![0i64; f.len() + g.len()];
            for (i, &fi) in f.iter().enumerate() {
                for (j, &gj) in g.iter().enumerate() {
                    if i + j >= 1 {
                        direct[i + j - 1] += (i as i64 - j as i64) * fi * gj;
                    }
                }
            }
            assert_eq!(wronskian(&f, &g), trim(&direct), "f = {f:?}, g = {g:?}");
            assert_eq!(wronskian(&f, &f), vec![0]);
            assert_eq!(wronskian(&g, &f), neg(&wronskian(&f, &g)));
            assert_eq!(wronskian(&big(&f), &big(&g)), big(&trim(&direct)));
        }
    }

    #[test]
    fn wronskian_sign_for_interlacing_pairs() {
        // Products of linear factors with interlacing integer roots.
        let from_roots =
            |roots: &[i64]| roots.iter().fold(vec![1i64], |acc, &r| mul(&acc, &[-r, 1]));
        let pairs: [(&[i64], &[i64]); 4] = [
            (&[0], &[-1, 1]),
            (&[-2, 1], &[-3, 0, 2]),
            (&[-1, 0], &[0, 1]), // same degree, f's roots to the left
            (&[-2, -2], &[-3, -2, 0]),
        ];
        for (f_roots, g_roots) in pairs {
            let f = from_roots(f_roots);
            let g = from_roots(g_roots);
            assert_eq!(crate::check_weak_interlacing(&f, &g), Some(true));
            let w = wronskian(&f, &g);
            for x in -40i64..=40 {
                assert!(evaluate(&w, &x) <= 0, "f = {f:?}, g = {g:?}, x = {x}");
            }
        }
    }

    #[test]
    fn i128_exceeds_i64_range() {
        let p = [i128::from(i64::MAX), 1];
        let sq = mul(&p, &p);
        assert_eq!(sq[0], i128::from(i64::MAX) * i128::from(i64::MAX));
        assert_eq!(sq.len(), 3);
    }

    #[test]
    fn permanent_bigint_matches_i64() {
        let mat: Vec<Vec<Vec<i64>>> = (0..3)
            .map(|r| {
                (0..3)
                    .map(|c| vec![r as i64 - c as i64, 1, r as i64])
                    .collect()
            })
            .collect();
        let mat_big: Vec<Vec<Vec<BigInt>>> = mat
            .iter()
            .map(|row| row.iter().map(|e| big(e)).collect())
            .collect();
        assert_eq!(permanent(&mat_big), big(&permanent(&mat)));
    }

    #[test]
    #[should_panic(expected = "i64 overflow in polynomial coefficient multiplication")]
    fn i64_mul_overflow_panics() {
        let _ = mul(&[i64::MAX, 1], &[2]);
    }

    #[test]
    #[should_panic(expected = "i64 overflow in polynomial coefficient addition")]
    fn i64_add_overflow_panics() {
        let _ = add(&[i64::MAX], &[1]);
    }

    #[test]
    #[should_panic(expected = "i64 overflow in polynomial coefficient negation")]
    fn i64_neg_overflow_panics() {
        let _ = neg(&[i64::MIN]);
    }

    #[test]
    #[should_panic(expected = "i64 overflow")]
    fn i64_evaluate_overflow_panics() {
        let _ = evaluate(&[0i64, 0, 1], &(1i64 << 32));
    }

    #[test]
    #[should_panic(expected = "i64 overflow")]
    fn i64_pow_overflow_panics() {
        let _ = pow(&[2i64, 1], 64);
    }

    #[test]
    #[should_panic(expected = "i128 overflow in polynomial coefficient multiplication")]
    fn i128_mul_overflow_panics() {
        let _ = mul(&[i128::MAX], &[2]);
    }

    #[test]
    #[should_panic(expected = "Ratio<i64> overflow")]
    fn ratio_i64_overflow_panics() {
        let big = Ratio::<i64>::from_integer(i64::MAX);
        let _ = mul(&[big], &[Ratio::from_integer(2)]);
    }

    #[test]
    fn i64_near_limits_without_overflow() {
        assert_eq!(add(&[i64::MAX], &[i64::MIN]), vec![-1]);
        assert_eq!(mul(&[i64::MIN], &[1]), vec![i64::MIN]);
        assert_eq!(evaluate(&[0i64, 0, 1], &(1i64 << 31)), 1i64 << 62);
    }
}
