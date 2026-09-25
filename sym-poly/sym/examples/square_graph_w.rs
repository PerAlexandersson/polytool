//! W coefficients for Section 2.3's square graphs: W_mu = C_n^{-1}(f_mu).
//! Default: even degrees 2 through 8 at n=2*d; otherwise `degree [num_vars]`.
//! Printed coefficients use the permutation average h_lambda/n!.
//! Each row also prints the direct forgotten expansion, without applying C_n.

use std::error::Error;

use num_bigint::BigInt;
use num_rational::BigRational;
use sym_poly_core::Partition;
use sym_poly_sym::{c_n, square_partition_graph_symmetric};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 2 {
        return Err("usage: square_graph_w [even_degree [num_vars]]".into());
    }
    let degrees: Vec<u32> = if let Some(arg) = args.first() {
        let degree = arg.parse()?;
        if degree % 2 != 0 {
            return Err("degree must be even".into());
        }
        vec![degree]
    } else {
        vec![2, 4, 6, 8]
    };
    let variables = args.get(1).map(|arg| arg.parse::<usize>()).transpose()?;
    for degree in degrees {
        let n = variables.unwrap_or(
            usize::try_from(degree)?
                .checked_mul(2)
                .ok_or("size overflow")?,
        );
        if n < usize::try_from(degree)? {
            return Err("this W expansion example requires n >= degree".into());
        }
        let normalization: BigInt = (1..=n).map(BigInt::from).product();
        for lambda in Partition::all_of_size(degree)
            .into_iter()
            .filter(|p| !p.parts().contains(&1))
        {
            let graph = square_partition_graph_symmetric(&lambda, n)?;
            let coefficients = c_n(&graph, n).to_forgotten_basis();
            assert!(coefficients
                .terms()
                .keys()
                .all(|mu| !mu.parts().contains(&1)));
            let terms: Vec<_> = coefficients
                .terms()
                .iter()
                .map(|(mu, c)| {
                    format!(
                        "({}) W[{}]",
                        BigRational::new(c.clone(), normalization.clone()),
                        mu
                    )
                })
                .collect();
            println!("n={n} lambda={lambda}: h/n! = {}", terms.join(" + "));
            let direct = graph.to_forgotten_basis();
            let terms: Vec<_> = direct
                .terms()
                .iter()
                .map(|(mu, c)| {
                    format!(
                        "({}) f[{}]",
                        BigRational::new(c.clone(), normalization.clone()),
                        mu
                    )
                })
                .collect();
            println!("n={n} lambda={lambda}: h/n! = {}", terms.join(" + "));
        }
    }
    Ok(())
}
