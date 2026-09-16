//! Degree-by-degree quotients of homogeneous polynomial ideals.
//!
//! This is complementary to the Groebner-based quotient code. For a fixed
//! multidegree it builds the ideal component directly by multiplying
//! homogeneous generators by all monomials of the complementary multidegree and
//! then takes a linear quotient. This is useful for modular representation
//! computations where only traces in bounded degrees are needed.

use std::collections::BTreeMap;

use sym_poly_core::linear_algebra::{quotient_action_matrix, QuotientSpace};
use sym_poly_core::{Field, Ring};

use crate::indexed_variables::IndexedVariables;
use crate::MultiPoly;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GradedQuotientComponent<C: Field> {
    pub multidegree: Vec<u32>,
    pub ambient_monomials: Vec<Vec<u32>>,
    pub quotient: QuotientSpace<C>,
}

impl<C: Field> GradedQuotientComponent<C> {
    pub fn dimension(&self) -> usize {
        self.quotient.dimension()
    }

    pub fn action_matrix_by_index_permutation(
        &self,
        variables: &IndexedVariables,
        index_permutation: &[usize],
    ) -> Vec<Vec<C>> {
        let ambient_index: BTreeMap<_, _> = self
            .ambient_monomials
            .iter()
            .cloned()
            .enumerate()
            .map(|(i, monomial)| (monomial, i))
            .collect();
        let mut ambient_action = sym_poly_core::linear_algebra::zero_matrix::<C>(
            self.ambient_monomials.len(),
            self.ambient_monomials.len(),
        );

        for (col, monomial) in self.ambient_monomials.iter().enumerate() {
            let basis_element: MultiPoly<C> =
                MultiPoly::x_power(variables.num_vars(), monomial.clone());
            let image = variables.permute_polynomial(&basis_element, index_permutation);
            let (image_monomial, coeff) = image
                .terms()
                .iter()
                .next()
                .expect("permuted monomial should be nonzero");
            debug_assert_eq!(*coeff, C::one());
            let row = ambient_index[image_monomial];
            ambient_action[row][col] = C::one();
        }

        quotient_action_matrix(&self.quotient, &ambient_action)
    }
}

pub fn graded_quotient_component<C: Field>(
    variables: &IndexedVariables,
    generators: &[MultiPoly<C>],
    multidegree: &[u32],
) -> GradedQuotientComponent<C> {
    assert_eq!(
        multidegree.len(),
        variables.num_alphabets(),
        "multidegree has wrong number of alphabets"
    );
    assert!(
        generators
            .iter()
            .all(|generator| generator.num_vars() == variables.num_vars()),
        "generators must use the indexed variable count"
    );

    let ambient_monomials = monomials_with_multidegree(variables, multidegree);
    let ambient_index: BTreeMap<_, _> = ambient_monomials
        .iter()
        .cloned()
        .enumerate()
        .map(|(i, monomial)| (monomial, i))
        .collect();
    let mut relations = Vec::new();

    for generator in generators {
        if generator.is_zero() {
            continue;
        };
        let generator_degree = polynomial_multidegree(variables, generator)
            .expect("nonzero generators must be homogeneous in the indexed multigrading");
        let Some(complement_degree) = multidegree_difference(multidegree, &generator_degree) else {
            continue;
        };
        for multiplier in monomials_with_multidegree(variables, &complement_degree) {
            let shifted = multiply_polynomial_by_monomial(generator, &multiplier);
            relations.push(polynomial_coordinate_vector(
                &shifted,
                &ambient_index,
                ambient_monomials.len(),
            ));
        }
    }

    let quotient = QuotientSpace::from_relations(ambient_monomials.len(), &relations);
    GradedQuotientComponent {
        multidegree: multidegree.to_vec(),
        ambient_monomials,
        quotient,
    }
}

pub fn monomials_with_multidegree(
    variables: &IndexedVariables,
    multidegree: &[u32],
) -> Vec<Vec<u32>> {
    assert_eq!(
        multidegree.len(),
        variables.num_alphabets(),
        "multidegree has wrong number of alphabets"
    );
    let mut per_alphabet = Vec::new();
    for &degree in multidegree {
        per_alphabet.push(weak_compositions(degree, variables.num_indices()));
    }

    let mut monomials = Vec::new();
    let mut current = vec![0u32; variables.num_vars()];
    assemble_alphabet_monomials(variables, &per_alphabet, 0, &mut current, &mut monomials);
    monomials.sort();
    monomials
}

pub fn polynomial_multidegree<C: Ring>(
    variables: &IndexedVariables,
    polynomial: &MultiPoly<C>,
) -> Option<Vec<u32>> {
    if polynomial.is_zero() {
        return None;
    }
    let mut degree = None;
    for monomial in polynomial.terms().keys() {
        let monomial_degree = variables.monomial_multidegree(monomial);
        match &degree {
            None => degree = Some(monomial_degree),
            Some(old) if old == &monomial_degree => {}
            Some(_) => return None,
        }
    }
    degree
}

fn assemble_alphabet_monomials(
    variables: &IndexedVariables,
    per_alphabet: &[Vec<Vec<u32>>],
    alphabet: usize,
    current: &mut [u32],
    monomials: &mut Vec<Vec<u32>>,
) {
    if alphabet == variables.num_alphabets() {
        monomials.push(current.to_vec());
        return;
    }

    for composition in &per_alphabet[alphabet] {
        for (index, &exponent) in composition.iter().enumerate() {
            current[variables.variable_index(alphabet, index)] = exponent;
        }
        assemble_alphabet_monomials(variables, per_alphabet, alphabet + 1, current, monomials);
        for index in 0..variables.num_indices() {
            current[variables.variable_index(alphabet, index)] = 0;
        }
    }
}

fn weak_compositions(total: u32, length: usize) -> Vec<Vec<u32>> {
    if length == 0 {
        return if total == 0 { vec![vec![]] } else { vec![] };
    }
    let mut result = Vec::new();
    let mut current = vec![0u32; length];
    weak_compositions_rec(total, 0, &mut current, &mut result);
    result
}

fn weak_compositions_rec(
    remaining: u32,
    index: usize,
    current: &mut [u32],
    result: &mut Vec<Vec<u32>>,
) {
    if index + 1 == current.len() {
        current[index] = remaining;
        result.push(current.to_vec());
        current[index] = 0;
        return;
    }
    for value in 0..=remaining {
        current[index] = value;
        weak_compositions_rec(remaining - value, index + 1, current, result);
    }
    current[index] = 0;
}

fn multidegree_difference(target: &[u32], subtrahend: &[u32]) -> Option<Vec<u32>> {
    assert_eq!(
        target.len(),
        subtrahend.len(),
        "multidegrees have different lengths"
    );
    target
        .iter()
        .zip(subtrahend.iter())
        .map(|(&a, &b)| if a >= b { Some(a - b) } else { None })
        .collect()
}

fn multiply_polynomial_by_monomial<C: Ring>(
    polynomial: &MultiPoly<C>,
    monomial: &[u32],
) -> MultiPoly<C> {
    assert_eq!(
        polynomial.num_vars(),
        monomial.len(),
        "monomial has wrong number of variables"
    );
    let terms = polynomial
        .terms()
        .iter()
        .map(|(exponents, coeff)| {
            (
                exponents
                    .iter()
                    .zip(monomial.iter())
                    .map(|(&a, &b)| a.checked_add(b).expect("monomial exponent overflow"))
                    .collect(),
                coeff.clone(),
            )
        })
        .collect();
    MultiPoly::from_terms(polynomial.num_vars(), terms)
}

fn polynomial_coordinate_vector<C: Ring>(
    polynomial: &MultiPoly<C>,
    ambient_index: &BTreeMap<Vec<u32>, usize>,
    ambient_dimension: usize,
) -> Vec<C> {
    let mut vector = vec![C::zero(); ambient_dimension];
    for (monomial, coeff) in polynomial.terms() {
        let &index = ambient_index
            .get(monomial)
            .expect("relation term should lie in the target multidegree");
        vector[index] = coeff.clone();
    }
    vector
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_rational::Ratio;

    type Q = Ratio<i64>;

    fn q(n: i64) -> Q {
        Q::from_integer(n)
    }

    fn mono(exponents: &[u32], coefficient: i64) -> MultiPoly<Q> {
        MultiPoly::monomial(exponents.len(), exponents.to_vec(), q(coefficient))
    }

    #[test]
    fn test_monomials_with_multidegree_two_alphabets() {
        let variables = IndexedVariables::new(2, 2);
        let monomials = monomials_with_multidegree(&variables, &[1, 1]);

        assert_eq!(
            monomials,
            vec![
                vec![0, 1, 0, 1],
                vec![0, 1, 1, 0],
                vec![1, 0, 0, 1],
                vec![1, 0, 1, 0],
            ]
        );
    }

    #[test]
    fn test_polynomial_multidegree() {
        let variables = IndexedVariables::new(2, 2);
        let homogeneous = mono(&[1, 0, 1, 0], 1) + mono(&[0, 1, 0, 1], 1);
        let inhomogeneous = homogeneous.clone() + mono(&[1, 0, 0, 0], 1);

        assert_eq!(
            polynomial_multidegree(&variables, &homogeneous),
            Some(vec![1, 1])
        );
        assert_eq!(polynomial_multidegree(&variables, &inhomogeneous), None);
    }

    #[test]
    fn test_graded_quotient_component_one_relation() {
        let variables = IndexedVariables::new(1, 2);
        let generator = mono(&[1, 0], 1) + mono(&[0, 1], 1);
        let component = graded_quotient_component(&variables, &[generator], &[1]);

        assert_eq!(component.ambient_monomials.len(), 2);
        assert_eq!(component.dimension(), 1);
    }

    #[test]
    fn test_public_zero_index_graded_quotient_paths() {
        let variables = crate::IndexedVariables::new(2, 0);

        assert_eq!(
            crate::monomials_with_multidegree(&variables, &[0, 0]),
            vec![vec![]]
        );
        assert!(crate::monomials_with_multidegree(&variables, &[1, 0]).is_empty());

        let component = crate::graded_quotient_component::<Q>(&variables, &[], &[0, 0]);
        assert_eq!(component.ambient_monomials, vec![vec![]]);
        assert_eq!(component.dimension(), 1);
        assert_eq!(
            component.action_matrix_by_index_permutation(&variables, &[]),
            vec![vec![q(1)]]
        );

        let positive_degree = crate::graded_quotient_component::<Q>(&variables, &[], &[1, 0]);
        assert!(positive_degree.ambient_monomials.is_empty());
        assert_eq!(positive_degree.dimension(), 0);
    }

    #[test]
    #[should_panic(expected = "monomial exponent overflow")]
    fn test_relation_monomial_shift_rejects_overflow() {
        let polynomial = mono(&[u32::MAX], 1);

        let _ = multiply_polynomial_by_monomial(&polynomial, &[1]);
    }

    #[test]
    #[should_panic(expected = "nonzero generators must be homogeneous")]
    fn test_graded_quotient_component_rejects_inhomogeneous_generator() {
        let variables = IndexedVariables::new(1, 2);
        let inhomogeneous = mono(&[1, 0], 1) + mono(&[0, 0], -1);

        let _ = graded_quotient_component(&variables, &[inhomogeneous], &[1]);
    }
}
