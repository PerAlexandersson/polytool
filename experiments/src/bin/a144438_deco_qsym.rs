//! Test the literal Gessel-fundamental refinement of the A144438
//! decorated-permutation descent enumerator.
//!
//! For `rho in S_n`, let `b_r` be the zero-based position of `n-r+1` in
//! the restriction of `rho` to `{n-r+1, ..., n}`.  We compute
//!
//!     Q_n = sum_rho 2^{#{r : 2 <= r <= n-1, (b_r,b_{r+1})=(0,1)}} F_Des(rho).
//!
//! The specialization `F_Des(rho) -> t^{des(rho)}` is the A144438
//! polynomial `P_{n-1}(t)` from the associated research project.

use std::collections::BTreeMap;

use combpoly::permutation::all_permutations;
use num_bigint::BigInt;
use sym_poly_core::Composition;
use sym_poly_qsym::{descent_set_to_composition, symmetric_qsym_to_sym, QSymBasis, QSymFunction};

fn chronological_minimum_insertion_code(permutation: &[u8]) -> Vec<usize> {
    let n = permutation.len();
    (1..=n)
        .map(|r| {
            let value = (n - r + 1) as u8;
            let position = permutation
                .iter()
                .position(|&entry| entry == value)
                .expect("the input must be a permutation of 1,...,n");
            permutation[..position]
                .iter()
                .filter(|&&entry| entry >= value)
                .count()
        })
        .collect()
}

fn eligible_01_count(permutation: &[u8]) -> usize {
    let code = chronological_minimum_insertion_code(permutation);
    // Vector indices r-1 and r represent b_r and b_{r+1}.  The allowed
    // range 2 <= r <= n-1 becomes 1 <= r-1 < n-1.
    (1..code.len().saturating_sub(1))
        .filter(|&index| code[index] == 0 && code[index + 1] == 1)
        .count()
}

fn descent_composition(permutation: &[u8]) -> Composition {
    let descents = permutation
        .windows(2)
        .enumerate()
        .filter_map(|(index, pair)| (pair[0] > pair[1]).then_some((index + 1) as u32))
        .collect::<Vec<_>>();
    descent_set_to_composition(&descents, permutation.len() as u32)
}

fn deco_descent_qsym(n: u8) -> QSymFunction<BigInt> {
    let mut terms = BTreeMap::new();
    for permutation in all_permutations(n) {
        let composition = descent_composition(&permutation);
        let weight = BigInt::from(1u8) << eligible_01_count(&permutation);
        *terms
            .entry(composition)
            .or_insert_with(|| BigInt::from(0u8)) += weight;
    }
    QSymFunction::from_terms(QSymBasis::Fundamental, terms)
}

fn first_monomial_symmetry_mismatch(
    function: &QSymFunction<BigInt>,
) -> Option<(Composition, BigInt, Composition, BigInt)> {
    let monomial = function.to_monomial_basis();
    for composition in monomial.terms().keys() {
        let reversed = Composition::new(composition.parts().iter().copied().rev().collect());
        let coefficient = monomial.coefficient(composition);
        let reversed_coefficient = monomial.coefficient(&reversed);
        if coefficient != reversed_coefficient {
            return Some((
                composition.clone(),
                coefficient,
                reversed,
                reversed_coefficient,
            ));
        }
    }
    None
}

fn main() {
    for n in 1..=8 {
        let function = deco_descent_qsym(n);
        let symmetric = symmetric_qsym_to_sym(&function).is_some();
        let qsym_schur = function.to_quasisymmetric_schur_basis();
        let first_negative_qs = qsym_schur
            .terms()
            .iter()
            .find(|(_, coefficient)| **coefficient < BigInt::from(0u8));
        println!(
            "n={n}: symmetric={symmetric}, QS-positive={}",
            first_negative_qs.is_none()
        );
        if let Some((alpha, coefficient)) = first_negative_qs {
            println!("  QS obstruction: [QS_{alpha}]={coefficient}");
        }
        if let Some((alpha, a, beta, b)) = first_monomial_symmetry_mismatch(&function) {
            println!("  certificate: [M_{alpha}]={a}, [M_{beta}]={b}");
        }
        if n == 3 {
            println!("  F expansion: {function}");
            println!("  M expansion: {}", function.to_monomial_basis());
            println!("  QS expansion: {qsym_schur}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insertion_code_detects_the_unique_degree_three_decoration() {
        let decorated = all_permutations(3)
            .into_iter()
            .filter(|permutation| eligible_01_count(permutation) == 1)
            .collect::<Vec<_>>();
        assert_eq!(decorated, vec![vec![2, 1, 3]]);
    }

    #[test]
    fn degree_three_is_the_first_nonsymmetric_case() {
        assert!(symmetric_qsym_to_sym(&deco_descent_qsym(1)).is_some());
        assert!(symmetric_qsym_to_sym(&deco_descent_qsym(2)).is_some());
        assert!(symmetric_qsym_to_sym(&deco_descent_qsym(3)).is_none());

        let degree_three = deco_descent_qsym(3).to_monomial_basis();
        assert_eq!(
            degree_three.coefficient(&Composition::new(vec![1, 2])),
            BigInt::from(4u8)
        );
        assert_eq!(
            degree_three.coefficient(&Composition::new(vec![2, 1])),
            BigInt::from(3u8)
        );
    }

    #[test]
    fn quasisymmetric_schur_positivity_first_fails_in_degree_six() {
        for n in 1..=5 {
            assert!(
                deco_descent_qsym(n)
                    .to_quasisymmetric_schur_basis()
                    .terms()
                    .values()
                    .all(|coefficient| coefficient >= &BigInt::from(0u8)),
                "expected QS-positivity in degree {n}"
            );
        }
        let degree_six = deco_descent_qsym(6).to_quasisymmetric_schur_basis();
        assert_eq!(
            degree_six.coefficient(&Composition::new(vec![1, 1, 2, 1, 1])),
            BigInt::from(-8)
        );
    }
}
