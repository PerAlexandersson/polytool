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
