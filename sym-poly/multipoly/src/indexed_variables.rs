//! Indexed variable systems for diagonal `S_n` actions.
//!
//! Variables are flattened by `variable = alphabet * n + index`. A permutation
//! `sigma in S_n` acts simultaneously in every alphabet:
//! `x_{a,i} -> x_{a,sigma(i)}`.

use std::collections::BTreeMap;

use sym_poly_core::linear_algebra::Matrix;
use sym_poly_core::sn_action::{
    assert_permutation, conjugacy_class_representatives, simple_transposition,
};
use sym_poly_core::{Field, Partition};

use crate::groebner::GroebnerBasis;
use crate::quotient::{
    permute_variables, quotient_action_matrix_by_permutation, restrict_matrix_to_indices,
    QuotientBasis,
};
use crate::MultiPoly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IndexedVariables {
    num_alphabets: usize,
    num_indices: usize,
}

impl IndexedVariables {
    pub fn new(num_alphabets: usize, num_indices: usize) -> Self {
        assert!(num_alphabets > 0, "there must be at least one alphabet");
        Self {
            num_alphabets,
            num_indices,
        }
    }

    pub fn num_alphabets(&self) -> usize {
        self.num_alphabets
    }

    pub fn num_indices(&self) -> usize {
        self.num_indices
    }

    pub fn num_vars(&self) -> usize {
        self.num_alphabets
            .checked_mul(self.num_indices)
            .expect("indexed variable count overflow")
    }

    pub fn variable_index(&self, alphabet: usize, index: usize) -> usize {
        assert!(alphabet < self.num_alphabets, "alphabet out of range");
        assert!(index < self.num_indices, "index out of range");
        alphabet
            .checked_mul(self.num_indices)
            .and_then(|offset| offset.checked_add(index))
            .expect("indexed variable index overflow")
    }

    pub fn alphabet_and_index(&self, variable: usize) -> (usize, usize) {
        assert!(variable < self.num_vars(), "variable out of range");
        (variable / self.num_indices, variable % self.num_indices)
    }

    pub fn variable_permutation(&self, index_permutation: &[usize]) -> Vec<usize> {
        assert_permutation(index_permutation);
        assert_eq!(
            index_permutation.len(),
            self.num_indices,
            "index permutation has wrong size"
        );

        let mut variable_permutation = vec![0usize; self.num_vars()];
        for alphabet in 0..self.num_alphabets {
            for (source_index, &target_index) in index_permutation.iter().enumerate() {
                let source = self.variable_index(alphabet, source_index);
                let target = self.variable_index(alphabet, target_index);
                variable_permutation[source] = target;
            }
        }
        variable_permutation
    }

    pub fn simple_transposition_variable_permutation(&self, i: usize) -> Vec<usize> {
        self.variable_permutation(&simple_transposition(self.num_indices, i))
    }

    pub fn monomial_multidegree(&self, exponents: &[u32]) -> Vec<u32> {
        assert_eq!(
            exponents.len(),
            self.num_vars(),
            "monomial has wrong number of variables"
        );
        let mut degrees = vec![0u32; self.num_alphabets];
        for (variable, &exponent) in exponents.iter().enumerate() {
            let (alphabet, _) = self.alphabet_and_index(variable);
            degrees[alphabet] = degrees[alphabet]
                .checked_add(exponent)
                .expect("monomial multidegree overflow");
        }
        degrees
    }

    pub fn permute_polynomial<C: Field>(
        &self,
        polynomial: &MultiPoly<C>,
        index_permutation: &[usize],
    ) -> MultiPoly<C> {
        assert_eq!(
            polynomial.num_vars(),
            self.num_vars(),
            "polynomial has wrong number of variables"
        );
        permute_variables(polynomial, &self.variable_permutation(index_permutation))
    }
}

pub fn quotient_basis_multidegrees(
    variables: &IndexedVariables,
    basis: &QuotientBasis,
) -> BTreeMap<Vec<u32>, Vec<usize>> {
    assert_eq!(
        basis.num_vars,
        variables.num_vars(),
        "basis and indexed variables have incompatible sizes"
    );

    let mut by_multidegree = BTreeMap::new();
    for (i, monomial) in basis.monomials.iter().enumerate() {
        by_multidegree
            .entry(variables.monomial_multidegree(monomial))
            .or_insert_with(Vec::new)
            .push(i);
    }
    by_multidegree
}

pub fn quotient_action_matrix_multidegree_blocks<C: Field>(
    variables: &IndexedVariables,
    basis: &QuotientBasis,
    action_matrix: &[Vec<C>],
) -> BTreeMap<Vec<u32>, Matrix<C>> {
    quotient_basis_multidegrees(variables, basis)
        .into_iter()
        .map(|(degree, indices)| (degree, restrict_matrix_to_indices(action_matrix, &indices)))
        .collect()
}

pub fn quotient_action_matrices_by_index_permutation_and_multidegree<C: Field>(
    variables: &IndexedVariables,
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
    index_permutation: &[usize],
) -> Option<BTreeMap<Vec<u32>, Matrix<C>>> {
    assert_eq!(
        gb.num_vars,
        variables.num_vars(),
        "Groebner basis and indexed variables have incompatible sizes"
    );
    let variable_permutation = variables.variable_permutation(index_permutation);
    let action = quotient_action_matrix_by_permutation(gb, basis, &variable_permutation)?;
    if !is_multidegree_preserving_action_matrix(variables, basis, &action) {
        return None;
    }
    Some(quotient_action_matrix_multidegree_blocks(
        variables, basis, &action,
    ))
}

pub fn quotient_action_matrices_by_multidegree_and_cycle_type<C: Field>(
    variables: &IndexedVariables,
    gb: &GroebnerBasis<C>,
    basis: &QuotientBasis,
) -> Option<BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>>> {
    let mut by_degree: BTreeMap<Vec<u32>, BTreeMap<Partition, Matrix<C>>> = BTreeMap::new();

    for (cycle_type, representative) in conjugacy_class_representatives(variables.num_indices()) {
        let blocks = quotient_action_matrices_by_index_permutation_and_multidegree(
            variables,
            gb,
            basis,
            &representative,
        )?;
        for (degree, matrix) in blocks {
            by_degree
                .entry(degree)
                .or_default()
                .insert(cycle_type.clone(), matrix);
        }
    }

    Some(by_degree)
}

pub fn is_multidegree_preserving_action_matrix<C: Field>(
    variables: &IndexedVariables,
    basis: &QuotientBasis,
    action_matrix: &[Vec<C>],
) -> bool {
    if action_matrix.len() != basis.dimension()
        || action_matrix
            .iter()
            .any(|row| row.len() != basis.dimension())
    {
        return false;
    }

    let degrees: Vec<_> = basis
        .monomials
        .iter()
        .map(|monomial| variables.monomial_multidegree(monomial))
        .collect();
    for row in 0..basis.dimension() {
        for col in 0..basis.dimension() {
            if degrees[row] != degrees[col] && !action_matrix[row][col].is_zero() {
                return false;
            }
        }
    }
    true
}

pub fn ideal_generators_are_invariant_under_index_permutation<C: Field>(
    variables: &IndexedVariables,
    generators: &[MultiPoly<C>],
    gb: &GroebnerBasis<C>,
    index_permutation: &[usize],
) -> bool {
    generators.iter().all(|generator| {
        let image = variables.permute_polynomial(generator, index_permutation);
        gb.normal_form(&image).is_zero()
    })
}

pub fn ideal_generators_are_sn_invariant<C: Field>(
    variables: &IndexedVariables,
    generators: &[MultiPoly<C>],
    gb: &GroebnerBasis<C>,
) -> bool {
    if variables.num_indices() < 2 {
        return true;
    }

    (0..variables.num_indices() - 1).all(|i| {
        let transposition = simple_transposition(variables.num_indices(), i);
        ideal_generators_are_invariant_under_index_permutation(
            variables,
            generators,
            gb,
            &transposition,
        )
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{elementary_symmetric_generators, quotient_basis, MonomialOrder, MultiPoly};
    use num_rational::Ratio;
    use sym_poly_core::linear_algebra::matrix_trace;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn test_indexed_variable_conventions() {
        let variables = IndexedVariables::new(2, 3);

        assert_eq!(variables.num_vars(), 6);
        assert_eq!(variables.variable_index(1, 2), 5);
        assert_eq!(variables.alphabet_and_index(4), (1, 1));
        assert_eq!(
            variables.variable_permutation(&[1, 0, 2]),
            vec![1, 0, 2, 4, 3, 5]
        );
        assert_eq!(
            variables.monomial_multidegree(&[2, 1, 0, 0, 3, 1]),
            vec![3, 4]
        );
    }

    #[test]
    fn test_indexed_variable_polynomial_permutation() {
        let variables = IndexedVariables::new(2, 2);
        let f = MultiPoly::monomial(4, vec![1, 0, 0, 2], q(1));

        assert_eq!(
            variables.permute_polynomial(&f, &[1, 0]),
            MultiPoly::monomial(4, vec![0, 1, 2, 0], q(1))
        );
    }

    #[test]
    #[should_panic(expected = "monomial multidegree overflow")]
    fn test_monomial_multidegree_rejects_overflow() {
        let variables = IndexedVariables::new(1, 2);

        let _ = variables.monomial_multidegree(&[u32::MAX, 1]);
    }

    #[test]
    fn test_quotient_multidegree_blocks_for_artin_s2() {
        let variables = IndexedVariables::new(1, 2);
        let generators = elementary_symmetric_generators::<Q>(2);
        let gb = GroebnerBasis::new(generators.clone(), MonomialOrder::Lex);
        let basis = quotient_basis(&gb).unwrap();
        let action = quotient_action_matrix_by_permutation(
            &gb,
            &basis,
            &variables.variable_permutation(&[1, 0]),
        )
        .unwrap();

        assert!(ideal_generators_are_sn_invariant(
            &variables,
            &generators,
            &gb
        ));
        assert!(is_multidegree_preserving_action_matrix(
            &variables, &basis, &action
        ));
        assert_eq!(
            quotient_basis_multidegrees(&variables, &basis),
            BTreeMap::from([(vec![0], vec![0]), (vec![1], vec![1])])
        );
        let blocks = quotient_action_matrices_by_index_permutation_and_multidegree(
            &variables,
            &gb,
            &basis,
            &[1, 0],
        )
        .unwrap();
        assert_eq!(matrix_trace(&blocks[&vec![0]]), q(1));
        assert_eq!(matrix_trace(&blocks[&vec![1]]), q(-1));

        let by_cycle_type =
            quotient_action_matrices_by_multidegree_and_cycle_type(&variables, &gb, &basis)
                .unwrap();
        assert_eq!(matrix_trace(&by_cycle_type[&vec![0]][&p(&[1, 1])]), q(1));
        assert_eq!(matrix_trace(&by_cycle_type[&vec![0]][&p(&[2])]), q(1));
        assert_eq!(matrix_trace(&by_cycle_type[&vec![1]][&p(&[1, 1])]), q(1));
        assert_eq!(matrix_trace(&by_cycle_type[&vec![1]][&p(&[2])]), q(-1));
    }

    #[test]
    fn test_quotient_multidegree_blocks_reject_non_preserving_action() {
        let variables = IndexedVariables::new(1, 2);
        let x1 = MultiPoly::<Q>::var(2, 0);
        let x2 = MultiPoly::<Q>::var(2, 1);
        let generators = vec![
            x1.clone() + x2.clone() + MultiPoly::constant(2, q(-1)),
            x1 * x2,
        ];
        let gb = GroebnerBasis::new(generators, MonomialOrder::Lex);
        let basis = quotient_basis(&gb).unwrap();
        let action = quotient_action_matrix_by_permutation(
            &gb,
            &basis,
            &variables.variable_permutation(&[1, 0]),
        )
        .unwrap();

        assert!(!is_multidegree_preserving_action_matrix(
            &variables, &basis, &action
        ));
        assert!(
            quotient_action_matrices_by_index_permutation_and_multidegree(
                &variables,
                &gb,
                &basis,
                &[1, 0],
            )
            .is_none()
        );
        assert!(
            quotient_action_matrices_by_multidegree_and_cycle_type(&variables, &gb, &basis)
                .is_none()
        );
    }
}
