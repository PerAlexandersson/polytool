//! Row-wise clearing of rational denominators with exact integer scales.

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::One;

/// Return `(integer_row, scale)` with `integer_row = scale * row`.
/// Inputs must be canonical rationals (with positive denominators).
/// The scale is the positive least common multiple of the denominators; an
/// empty row has scale one. No primitive-content normalization is performed.
pub fn clear_row_denominators(row: &[BigRational]) -> (Vec<BigInt>, BigInt) {
    let mut scale = BigInt::one();
    for entry in row {
        if !entry.denom().is_one() {
            scale = scale.lcm(entry.denom());
        }
    }
    let entries = if scale.is_one() {
        row.iter().map(|entry| entry.numer().clone()).collect()
    } else {
        row.iter()
            .map(|entry| entry.numer() * (&scale / entry.denom()))
            .collect()
    };
    (entries, scale)
}

/// Clear each row independently, returning integer rows and positive row scales.
/// Rows may have different lengths; this operation does not validate matrix shape.
pub fn clear_matrix_row_denominators(
    matrix: &[Vec<BigRational>],
) -> (Vec<Vec<BigInt>>, Vec<BigInt>) {
    matrix.iter().map(|row| clear_row_denominators(row)).unzip()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn exact_scales_and_empty_rows() {
        let row = vec![
            BigRational::new((-1).into(), 6.into()),
            BigRational::new(1.into(), 4.into()),
        ];
        let (integers, scale) = clear_row_denominators(&row);
        assert_eq!(scale, BigInt::from(12));
        assert_eq!(integers, vec![BigInt::from(-2), BigInt::from(3)]);
        for (x, y) in row.iter().zip(&integers) {
            assert_eq!(
                x * &BigRational::from_integer(scale.clone()),
                BigRational::from_integer(y.clone())
            );
        }
        assert_eq!(clear_row_denominators(&[]), (vec![], BigInt::one()));
        let (_, scales) = clear_matrix_row_denominators(&[row, vec![]]);
        assert_eq!(scales, vec![BigInt::from(12), BigInt::from(1)]);
    }
}
