//! Exact C_n and C_n^{-1} on a Schur function and a partition-graph polynomial.

use std::error::Error;

use num_bigint::BigInt;
use num_rational::BigRational;
use sym_poly_core::Partition;
use sym_poly_sym::{c_n, c_n_inverse, partition_graph_symmetric, SymmetricFunction};

fn main() -> Result<(), Box<dyn Error>> {
    let lambda = Partition::new(vec![2, 1]);
    let schur = SymmetricFunction::<BigInt>::schur_symmetric(lambda.clone());
    assert_eq!(c_n(&schur, 3).coefficient(&lambda), BigInt::from(24));
    println!("C_3(s_(2,1)) = {}", c_n(&schur, 3));
    println!("C_3^(-1)(s_(2,1)) = {}", c_n_inverse(&schur, 3)?);

    let graph = partition_graph_symmetric(&Partition::new(vec![2]), 4)?;
    let inverse = c_n_inverse(&graph, 4)?;
    assert_eq!(
        inverse.coefficient(&Partition::new(vec![2])),
        BigRational::new(3.into(), 10.into())
    );
    assert_eq!(
        inverse.coefficient(&Partition::new(vec![1, 1])),
        BigRational::new((-5).into(), 6.into())
    );
    println!("b_(2)^(4) = {}", graph.to_schur_basis());
    println!("C_4^(-1)(b_(2)^(4)) = {inverse}");
    println!("In forgotten basis: {}", inverse.to_forgotten_basis());
    Ok(())
}
