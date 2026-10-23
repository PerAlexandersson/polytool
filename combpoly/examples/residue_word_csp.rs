//! One odd residue class, arbitrary bit-loop multiplicities.
use combpoly::residue_words::{
    binary_single_odd_residue_product_bigint, single_odd_residue_csp_polynomial_bigint,
};
use polytool::cyclic_sieving::root_of_unity_evaluation_bigint;

fn main() {
    let (modulus, width) = (3, 3);
    let binary = binary_single_odd_residue_product_bigint(modulus, width);
    let row_statistic = single_odd_residue_csp_polynomial_bigint(modulus, width, 1, 1);
    assert_ne!(binary, row_statistic);
    for shift in 0..modulus * width {
        let first = root_of_unity_evaluation_bigint(&binary, modulus * width, shift);
        let second = root_of_unity_evaluation_bigint(&row_statistic, modulus * width, shift);
        assert_eq!(first.integer_value, second.integer_value);
        println!("shift={shift}: fixed={}", first.integer_value.unwrap());
    }
}
