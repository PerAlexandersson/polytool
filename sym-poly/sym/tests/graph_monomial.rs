use std::collections::BTreeMap;

use combinatoric_core::next_permutation;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use sym_poly_core::Partition;
use sym_poly_multipoly::MultiPoly;
use sym_poly_sym::{
    partition_graph_symmetric, symmetrized_graph_monomial, Basis, GraphMonomialError,
    SymmetricFunction,
};

fn factorial(n: usize) -> BigInt {
    (1..=n).map(BigInt::from).product()
}

fn part(parts: &[u32]) -> Partition {
    Partition::new(parts.to_vec())
}

fn expected(basis: Basis, terms: &[(&[u32], i64)]) -> SymmetricFunction<BigInt> {
    SymmetricFunction::from_terms(
        basis,
        terms
            .iter()
            .map(|(p, c)| (part(p), BigInt::from(*c)))
            .collect(),
    )
}

// Independent literal S_n expansion for small graphs. No stabilizer formula.
fn literal_permutation_sum(n: usize, edges: &[(usize, usize)]) -> MultiPoly<BigInt> {
    let mut permutation: Vec<_> = (0..n).collect();
    let mut result = MultiPoly::zero(n);
    loop {
        let mut term = MultiPoly::constant(n, BigInt::one());
        for &(u, v) in edges {
            term = term * (MultiPoly::var(n, permutation[u]) - MultiPoly::var(n, permutation[v]));
        }
        result = result + term;
        if !next_permutation(&mut permutation) {
            return result;
        }
    }
}

fn expand_monomials(sf: &SymmetricFunction<BigInt>, n: usize) -> MultiPoly<BigInt> {
    let mut terms = BTreeMap::new();
    for (partition, coefficient) in sf.terms() {
        let mut exponent = partition.parts().to_vec();
        exponent.resize(n, 0);
        exponent.sort_unstable();
        loop {
            terms.insert(exponent.clone(), coefficient.clone());
            if !next_permutation(&mut exponent) {
                break;
            }
        }
    }
    MultiPoly::from_terms(n, terms)
}

fn evaluate(polynomial: &MultiPoly<BigInt>, point: &[i64]) -> BigInt {
    polynomial
        .terms()
        .iter()
        .map(|(exponents, coefficient)| {
            exponents
                .iter()
                .zip(point)
                .fold(coefficient.clone(), |value, (exponent, x)| {
                    value * BigInt::from(*x).pow(*exponent)
                })
        })
        .sum()
}

#[test]
fn literal_symmetrization_includes_isolated_vertices_and_multiplicities() {
    for (n, edges) in [
        (0, vec![]),
        (4, vec![]),
        (4, vec![(0, 1), (0, 2)]),
        (4, vec![(0, 1), (0, 2), (0, 3)]),
        (4, vec![(0, 1), (0, 1)]),
        (4, vec![(1, 3), (3, 1), (3, 2)]),
        (6, vec![(0, 1), (0, 2), (3, 4), (3, 5)]),
    ] {
        let sf = symmetrized_graph_monomial(n, &edges).unwrap();
        assert_eq!(
            expand_monomials(&sf, n).terms(),
            literal_permutation_sum(n, &edges).terms()
        );
    }
}

#[test]
fn degree_two_and_odd_degree_three_schur_expansions() {
    let b2 = partition_graph_symmetric(&part(&[2]), 4).unwrap();
    assert_eq!(b2, expected(Basis::Monomial, &[(&[2], 6), (&[1, 1], -4)]));
    assert_eq!(
        b2.to_schur_basis(),
        expected(Basis::Schur, &[(&[2], 6), (&[1, 1], -10)])
    );
    let b3 = partition_graph_symmetric(&part(&[3]), 6).unwrap();
    assert_eq!(
        b3,
        expected(
            Basis::Monomial,
            &[(&[3], 120), (&[2, 1], -72), (&[1, 1, 1], 72)]
        )
    );
    assert_eq!(
        b3.to_schur_basis(),
        expected(
            Basis::Schur,
            &[(&[3], 120), (&[2, 1], -192), (&[1, 1, 1], 336)]
        )
    );
}

#[test]
fn stars_are_disjoint_and_point_outwards() {
    for (parts, edges) in [
        (vec![3], vec![(0, 1), (0, 2), (0, 3)]),
        (vec![2, 2], vec![(0, 1), (0, 2), (3, 4), (3, 5)]),
        (vec![3, 2], vec![(0, 1), (0, 2), (0, 3), (4, 5), (4, 6)]),
    ] {
        let index = Partition::new(parts);
        let n = 2 * index.size() as usize;
        assert_eq!(
            partition_graph_symmetric(&index, n).unwrap(),
            symmetrized_graph_monomial(n, &edges).unwrap()
        );
    }
}

#[test]
fn two_partition_submatrix_is_diagonal_and_schur_roundtrips() {
    // Definition 1 adds the zero-exponent factorial to Lemma 6's coloring count.
    for degree in 2..=8 {
        let n = 2 * degree as usize;
        let indices: Vec<_> = Partition::all_of_size(degree)
            .into_iter()
            .filter(|p| p.parts().iter().all(|&a| a >= 2))
            .collect();
        for alpha in &indices {
            let sf = partition_graph_symmetric(alpha, n).unwrap();
            let mut multiplicities = BTreeMap::<u32, usize>::new();
            for &a in alpha.parts() {
                *multiplicities.entry(a).or_default() += 1;
            }
            let diagonal = factorial(n - alpha.num_parts())
                * multiplicities
                    .values()
                    .map(|&m| factorial(m))
                    .product::<BigInt>();
            for beta in &indices {
                assert_eq!(
                    sf.coefficient(beta),
                    if alpha == beta {
                        diagonal.clone()
                    } else {
                        BigInt::zero()
                    },
                    "alpha={alpha}, beta={beta}"
                );
            }
            assert_eq!(sf.to_schur_basis().to_monomial_basis(), sf);
        }
    }
}

#[test]
fn translation_invariance_and_homogeneity() {
    for parts in [vec![2], vec![3], vec![2, 2], vec![3, 2]] {
        let index = Partition::new(parts);
        let n = index.size() as usize + index.num_parts();
        let sf = partition_graph_symmetric(&index, n).unwrap();
        assert_eq!(sf.degree(), Some(index.size()));
        let polynomial = expand_monomials(&sf, n);
        let point: Vec<_> = (0..n).map(|i| (i * i + 2 * i) as i64 - 3).collect();
        let translated: Vec<_> = point.iter().map(|x| x + 7).collect();
        assert!(!evaluate(&polynomial, &point).is_zero());
        assert_eq!(
            evaluate(&polynomial, &point),
            evaluate(&polynomial, &translated)
        );
        assert_eq!(evaluate(&polynomial, &vec![2; n]), BigInt::zero());
    }
}

#[test]
fn variable_count_and_bigint_normalization_are_explicit() {
    for n in [3, 4, 6, 30] {
        let sf = partition_graph_symmetric(&part(&[2]), n).unwrap();
        assert_eq!(sf.coefficient(&part(&[2])), factorial(n - 1));
        assert_eq!(
            sf.coefficient(&part(&[1, 1])),
            -BigInt::from(2) * factorial(n - 2)
        );
        assert_eq!(
            sf.to_schur_basis().coefficient(&part(&[1, 1])),
            -BigInt::from(n + 1) * factorial(n - 2)
        );
    }
    assert_eq!(
        partition_graph_symmetric(&Partition::empty(), 0).unwrap(),
        expected(Basis::Monomial, &[(&[], 1)])
    );
    assert_eq!(
        partition_graph_symmetric(&Partition::empty(), 4).unwrap(),
        expected(Basis::Monomial, &[(&[], 24)])
    );
    assert!(Partition::all_of_size(1)
        .iter()
        .all(|p| partition_graph_symmetric(p, 2).is_err()));
}

#[test]
fn invalid_inputs_report_errors_before_graph_allocation() {
    assert_eq!(
        partition_graph_symmetric(&part(&[3, 1]), 8),
        Err(GraphMonomialError::NotTwoPartition)
    );
    assert_eq!(
        partition_graph_symmetric(&part(&[2, 2]), 5),
        Err(GraphMonomialError::TooFewVariables {
            required: 6,
            provided: 5
        })
    );
    assert_eq!(
        partition_graph_symmetric(&part(&[u32::MAX, 2]), 0),
        Err(GraphMonomialError::SizeOverflow)
    );
    assert_eq!(
        symmetrized_graph_monomial(3, &[(0, 3)]),
        Err(GraphMonomialError::VertexOutOfRange {
            vertex: 3,
            num_vars: 3
        })
    );
    assert_eq!(
        symmetrized_graph_monomial(3, &[(2, 2)]),
        Err(GraphMonomialError::Loop { vertex: 2 })
    );
}

#[test]
fn edge_reversal_and_vertex_relabeling() {
    let sf = symmetrized_graph_monomial(6, &[(0, 1), (0, 2), (0, 3)]).unwrap();
    let reversed = symmetrized_graph_monomial(6, &[(1, 0), (0, 2), (0, 3)]).unwrap();
    let relabeled = symmetrized_graph_monomial(6, &[(5, 2), (5, 0), (5, 4)]).unwrap();
    assert_eq!(reversed, -sf.clone());
    assert_eq!(relabeled, sf);
    assert!(symmetrized_graph_monomial(3, &[(0, 1)]).unwrap().is_zero());
}
