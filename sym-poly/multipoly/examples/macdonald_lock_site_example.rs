use num_rational::Ratio;
use sym_poly_core::Ring;
use sym_poly_multipoly::{
    lock_polynomial, nonsymmetric_macdonald_filling_formula, nonsymmetric_macdonald_q0,
    permuted_basement_macdonald_filling_formula, MultiPoly,
};

type Q = Ratio<i64>;

fn main() {
    // Assaf--Searles lock with one nonzero row (row 2): the key polynomial kappa_(0,2).
    let lock: MultiPoly<i64> = lock_polynomial(&[0, 2]);
    assert_eq!(lock.coefficient(&[2, 0]), 1);
    assert_eq!(lock.coefficient(&[1, 1]), 1);
    assert_eq!(lock.coefficient(&[0, 2]), 1);

    let permuted = permuted_basement_macdonald_filling_formula::<Q>(&[1, 1, 0, 1], &[2, 4, 1, 3]);
    assert_eq!(permuted.terms().len(), 3);
    assert!(!permuted.coefficient(&[1, 0, 1, 1]).is_zero());
    assert!(!permuted.coefficient(&[1, 1, 1, 0]).is_zero());
    assert!(!permuted.coefficient(&[0, 1, 1, 1]).is_zero());

    let identity_basement = nonsymmetric_macdonald_filling_formula::<Q>(&[1, 0, 2]);
    let q0_operator_side: MultiPoly<i64> = nonsymmetric_macdonald_q0(&[1, 0, 2], &2);
    assert!(!q0_operator_side.is_zero());
    assert!(identity_basement.terms().keys().all(|weight| {
        weight.iter().sum::<u32>() == 3 && weight.len() == q0_operator_side.num_vars()
    }));

    println!("lock_(0,2) = {lock}");
    println!(
        "permuted-basement Macdonald example has {} terms",
        permuted.terms().len()
    );
    println!(
        "identity-basement Macdonald E_(1,0,2) filling formula has {} terms",
        identity_basement.terms().len()
    );
}
