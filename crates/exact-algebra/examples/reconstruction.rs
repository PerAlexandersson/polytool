//! Bounded integer and rational reconstruction with exact denominator clearing.
use exact_algebra::chinese_remainder::ChineseRemainderAccumulator;
use exact_algebra::denominators::clear_row_denominators;
use exact_algebra::primes::ValidatedPrime;
use exact_algebra::rational_reconstruction::rational_reconstruction;
use num_rational::BigRational;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut reconstruction = ChineseRemainderAccumulator::new(1);
    reconstruction.combine_prime(&[3], ValidatedPrime::new(5)?)?;
    reconstruction.combine_prime(&[5], ValidatedPrime::new(7)?)?;
    let signed = reconstruction.reconstruct_symmetric(&[2.into()])?;
    assert_eq!(signed, vec![num_bigint::BigInt::from(-2)]);
    println!("signed = {signed:?}");
    let rational = rational_reconstruction(&44.into(), &101.into())?;
    assert_eq!(rational, BigRational::new(5.into(), 7.into()));
    println!("rational = {rational}");
    let row = [
        BigRational::new(1.into(), 2.into()),
        BigRational::new((-1).into(), 3.into()),
    ];
    let (integers, scale) = clear_row_denominators(&row);
    assert_eq!(
        integers,
        vec![num_bigint::BigInt::from(3), num_bigint::BigInt::from(-2)]
    );
    assert_eq!(scale, num_bigint::BigInt::from(6));
    println!("integer row = {integers:?}, scale = {scale}");
    Ok(())
}
