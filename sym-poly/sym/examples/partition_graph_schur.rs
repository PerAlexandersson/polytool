//! Schur expansions of Alexandersson--Shapiro's all-degree partition-graph basis.
//!
//! Run without arguments for degrees 2 through 6 at n = 2*d, or supply
//! `degree [num_vars]` for one degree. Coefficients use the full S_n sum.

use std::error::Error;

use num_bigint::BigInt;
use sym_poly_core::Partition;
use sym_poly_sym::partition_graph_symmetric;

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err("usage: partition_graph_schur [degree [num_vars]]".into());
    }
    let degrees = match args.first() {
        Some(value) => {
            let degree = value.parse::<u32>()?;
            degree..=degree
        }
        None => 2..=6,
    };
    let variables = args
        .get(1)
        .map(|value| value.parse::<usize>())
        .transpose()?;
    let b2 = partition_graph_symmetric(&Partition::new(vec![2]), 4)?.to_schur_basis();
    assert_eq!(
        b2.coefficient(&Partition::new(vec![1, 1])),
        BigInt::from(-10)
    );
    let b3 = partition_graph_symmetric(&Partition::new(vec![3]), 6)?.to_schur_basis();
    assert_eq!(
        b3.coefficient(&Partition::new(vec![2, 1])),
        BigInt::from(-192)
    );

    println!("Full permutation sums b_lambda^(n); stable basis range: n >= 2*d.");
    for degree in degrees {
        let n = variables.unwrap_or(
            usize::try_from(degree)?
                .checked_mul(2)
                .ok_or("variable count overflow")?,
        );
        for index in Partition::all_of_size(degree)
            .into_iter()
            .filter(|p| p.parts().iter().all(|&a| a >= 2))
        {
            let schur = partition_graph_symmetric(&index, n)?.to_schur_basis();
            println!("d={degree}, n={n}, lambda={index}: {schur}");
        }
    }
    Ok(())
}
