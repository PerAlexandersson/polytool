//! Multiply signed integer matrices with proved nonnegative output bounds.

use exact_algebra::nonnegative_matrix_product::multiply_integer_matrices_nonnegative_bounded;
use num_bigint::BigInt;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let left = vec![vec![BigInt::from(2), BigInt::from(-1)]];
    let right = vec![
        vec![BigInt::from(3), BigInt::from(1)],
        vec![BigInt::from(1), BigInt::from(2)],
    ];
    // Directly, the two entries are 2*3-1=5 and 2-2=0.
    let product =
        multiply_integer_matrices_nonnegative_bounded(&left, &right, &[5.into(), 0.into()])?;
    assert_eq!(product, vec![vec![BigInt::from(5), BigInt::from(0)]]);
    println!("{} {}", product[0][0], product[0][1]);
    Ok(())
}
