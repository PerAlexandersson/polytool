//! Compatibility checks relating integer-row clearing to polynomial coefficients.
use super::Q;
use exact_algebra::denominators::clear_row_denominators as integer_row;
use num_bigint::BigInt;
use num_traits::One;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integer_row_keeps_its_raw_negative_denominator_fast_path() {
        let row = [Q::new_raw(BigInt::one(), -BigInt::one())];
        assert_eq!(integer_row(&row), (vec![BigInt::one()], BigInt::one()));
        assert_eq!(
            crate::vec_poly::clear_denominators(&row),
            vec![-BigInt::one()]
        );
    }
}
