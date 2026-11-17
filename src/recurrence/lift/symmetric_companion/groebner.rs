//! A small exact Gröbner-basis solver over `Q` for the candidate search.
//!
//! Polynomials are maps from exponent vectors to nonzero rationals, ordered
//! lexicographically with variable 0 largest (`BTreeMap` order on the
//! exponent vectors), so the leading term is the last entry.  The solver
//! finds the rational points of a zero-dimensional system by computing a lex
//! basis with Buchberger's algorithm, taking the rational roots of its
//! univariate element in the smallest variable, substituting, and recursing.
//! Every step is exact and bounded by a step budget; the points are only
//! candidates for the caller, which certifies them separately.

use crate::factorization::irreducible_factorization_rational_coeffs_with_budget;
use crate::recurrence::BigRational as Q;
use exact_algebra::groebner::{
    reduced_lexicographic_groebner_basis_with_budget, GroebnerBudget, GroebnerError, Polynomial,
    RationalField, SparsePolynomial,
};
use num_traits::Zero;

type Monomial = Vec<u32>;
const FACTORIZATION_SUBSETS: usize = 1 << 12;

/// Construction and point-extraction adapter over the shared sparse representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct MultivariatePolynomial(Polynomial<RationalField>);
impl MultivariatePolynomial {
    pub(super) fn zero(variables: usize) -> Self {
        Self(Polynomial::zero(variables))
    }
    pub(super) fn add_term(&mut self, monomial: Monomial, coefficient: &Q) {
        self.0.add_term(monomial, coefficient);
    }
    pub(super) fn is_zero(&self) -> bool {
        self.0.is_zero()
    }
    fn is_constant(&self) -> bool {
        self.0
            .terms()
            .last_key_value()
            .is_some_and(|(m, _)| m.iter().all(|&e| e == 0))
    }
    /// Substitute `value` for variable `index`.
    fn substitute(&self, index: usize, value: &Q) -> Self {
        let mut out = Self::zero(self.0.num_vars());
        for (monomial, coefficient) in self.0.terms() {
            let mut reduced = monomial.clone();
            let power = std::mem::take(&mut reduced[index]);
            let mut factor = coefficient.clone();
            for _ in 0..power {
                factor *= value;
            }
            out.add_term(reduced, &factor);
        }
        out
    }

    /// The ascending coefficients in variable `index`, when no other variable
    /// occurs.
    fn univariate_in(&self, index: usize) -> Option<Vec<Q>> {
        let mut coefficients = Vec::new();
        for (monomial, coefficient) in self.0.terms() {
            if monomial
                .iter()
                .enumerate()
                .any(|(i, &e)| i != index && e != 0)
            {
                return None;
            }
            let power = monomial[index] as usize;
            if coefficients.len() <= power {
                coefficients.resize(power + 1, Q::zero());
            }
            coefficients[power] = coefficient.clone();
        }
        Some(coefficients)
    }
}
/// Why the solver stopped without an answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SolverStop {
    /// The step, size or coefficient budget was reached.
    Budget(String),
    /// The system has infinitely many complex solutions.
    PositiveDimensional,
}

fn groebner_basis(
    polynomials: &[MultivariatePolynomial],
    budget: &mut GroebnerBudget,
) -> Result<Vec<MultivariatePolynomial>, SolverStop> {
    let generators: Vec<_> = polynomials.iter().map(|p| p.0.clone()).collect();
    reduced_lexicographic_groebner_basis_with_budget(&generators, budget)
        .map(|basis| basis.into_iter().map(MultivariatePolynomial).collect())
        .map_err(|error| {
            SolverStop::Budget(match error {
                GroebnerError::BasisSizeBudgetExceeded => {
                    "the Gröbner basis grew beyond 400 elements".into()
                }
                other => other.to_string(),
            })
        })
}

/// Rational roots of a nonzero univariate polynomial.
fn rational_roots(coefficients: &[Q]) -> Result<Vec<Q>, SolverStop> {
    let factorization =
        irreducible_factorization_rational_coeffs_with_budget(coefficients, FACTORIZATION_SUBSETS)
            .map_err(|e| SolverStop::Budget(format!("univariate factorization failed: {e}")))?;
    Ok(factorization
        .factors
        .iter()
        .filter(|(factor, _)| factor.len() == 2)
        .map(|(factor, _)| Q::new(-factor[0].clone(), factor[1].clone()))
        .collect())
}

/// All rational points of the system `polynomials = 0` in `variables`
/// unknowns, when its complex solution set is finite.
pub(super) fn rational_points(
    polynomials: &[MultivariatePolynomial],
    variables: usize,
    max_steps: usize,
) -> Result<Vec<Vec<Q>>, SolverStop> {
    let mut budget = GroebnerBudget {
        remaining_steps: max_steps,
        max_coefficient_bits: 4096,
        max_basis_elements: 400,
    };
    let mut points = Vec::new();
    solve_from(
        polynomials.to_vec(),
        variables,
        vec![Q::zero(); variables],
        &mut budget,
        &mut points,
    )?;
    Ok(points)
}

/// Assign variables `index - 1, ..., 0` (variables `>= index` are assigned).
fn solve_from(
    polynomials: Vec<MultivariatePolynomial>,
    index: usize,
    assignment: Vec<Q>,
    budget: &mut GroebnerBudget,
    points: &mut Vec<Vec<Q>>,
) -> Result<(), SolverStop> {
    let basis = groebner_basis(&polynomials, budget)?;
    if basis.iter().any(MultivariatePolynomial::is_constant) {
        return Ok(());
    }
    if index == 0 {
        points.push(assignment);
        return Ok(());
    }
    let variable = index - 1;
    let Some(univariate) = basis
        .iter()
        .find_map(|p| p.univariate_in(variable).filter(|c| c.len() > 1))
    else {
        return Err(SolverStop::PositiveDimensional);
    };
    for root in rational_roots(&univariate)? {
        let substituted: Vec<MultivariatePolynomial> = basis
            .iter()
            .map(|p| p.substitute(variable, &root))
            .filter(|p| !p.is_zero())
            .collect();
        let mut next = assignment.clone();
        next[variable] = root;
        solve_from(substituted, variable, next, budget, points)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    fn poly(variables: usize, terms: &[(&[u32], i64)]) -> MultivariatePolynomial {
        let mut p = MultivariatePolynomial::zero(variables);
        for (m, c) in terms {
            p.add_term(m.to_vec(), &q(*c));
        }
        p
    }

    #[test]
    fn preserves_budget_diagnostics_and_recursive_sharing() {
        let x = poly(1, &[(&[1], 1)]);
        assert_eq!(
            rational_points(&[x.clone(), x.clone()], 1, 0),
            Err(SolverStop::Budget(
                "the Gröbner step budget was exhausted".into()
            ))
        );
        assert_eq!(rational_points(&[x.clone(), x], 1, 1), Ok(vec![vec![q(0)]]));
        let system = vec![
            poly(2, &[(&[2, 0], 1), (&[0, 0], -1)]),
            poly(2, &[(&[1, 1], 1), (&[0, 0], -2)]),
        ];
        let full = rational_points(&system, 2, 1000).unwrap();
        assert_eq!(full.len(), 2);
        let first_success = (0..1000)
            .find(|&steps| rational_points(&system, 2, steps).is_ok())
            .unwrap();
        assert!(first_success > 0);
        assert_eq!(rational_points(&system, 2, first_success).unwrap(), full);
        assert_eq!(
            rational_points(&system, 2, first_success - 1),
            Err(SolverStop::Budget(
                "the Gröbner step budget was exhausted".into()
            ))
        );
    }

    #[test]
    fn solves_small_systems() {
        // x^2 - 1 = 0, x y - 2 = 0: (1, 2), (-1, -2).
        let system = vec![
            poly(2, &[(&[2, 0], 1), (&[0, 0], -1)]),
            poly(2, &[(&[1, 1], 1), (&[0, 0], -2)]),
        ];
        let mut points = rational_points(&system, 2, 10_000).unwrap();
        points.sort();
        assert_eq!(points, vec![vec![q(-1), q(-2)], vec![q(1), q(2)]]);
        // x y = 0 alone is positive dimensional.
        let system = vec![poly(2, &[(&[1, 1], 1)])];
        assert_eq!(
            rational_points(&system, 2, 10_000),
            Err(SolverStop::PositiveDimensional)
        );
        // x^2 - 2 = 0 has no rational point; x = 1, x = 2 is inconsistent.
        let system = vec![poly(1, &[(&[2], 1), (&[0], -2)])];
        assert_eq!(rational_points(&system, 1, 10_000), Ok(vec![]));
        let system = vec![
            poly(1, &[(&[1], 1), (&[0], -1)]),
            poly(1, &[(&[1], 1), (&[0], -2)]),
        ];
        assert_eq!(rational_points(&system, 1, 10_000), Ok(vec![]));
    }
}
