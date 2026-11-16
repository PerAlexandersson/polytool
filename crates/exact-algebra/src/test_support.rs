//! Small matrix constructors shared by the `linalg` unit tests.

use num_bigint::BigInt;

pub(super) fn bi(v: i64) -> BigInt {
    BigInt::from(v)
}

pub(super) fn bi_mat(rows: &[&[i64]]) -> Vec<Vec<BigInt>> {
    rows.iter()
        .map(|row| row.iter().map(|&v| bi(v)).collect())
        .collect()
}

use num_traits::{One, Zero};
pub(super) fn bareiss_leading_principal_minors_bigint(mat: &[Vec<BigInt>]) -> Option<Vec<BigInt>> {
    let n = mat.len();
    if mat.iter().any(|row| row.len() != n) {
        return None;
    }
    if n == 0 {
        return Some(Vec::new());
    }

    let mut a = mat.to_vec();
    let mut denom = BigInt::one();
    let mut minors = Vec::with_capacity(n);

    for k in 0..n {
        let pivot = a[k][k].clone();
        if pivot.is_zero() {
            if k == n - 1 {
                minors.push(pivot);
                break;
            }
            return None;
        }
        minors.push(pivot.clone());
        if k == n - 1 {
            break;
        }

        let mut next_entries = Vec::with_capacity((n - k - 1) * (n - k - 1));
        for i in (k + 1)..n {
            for j in (k + 1)..n {
                let numerator = &a[i][j] * &pivot - &a[i][k] * &a[k][j];
                if &numerator % &denom != BigInt::zero() {
                    return None;
                }
                next_entries.push((i, j, numerator / &denom));
            }
        }
        for (i, j, value) in next_entries {
            a[i][j] = value;
        }
        denom = pivot;
    }

    Some(minors)
}
