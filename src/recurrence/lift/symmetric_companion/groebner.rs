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
use num_traits::{One, Zero};
use std::collections::BTreeMap;

type Monomial = Vec<u32>;

/// Recombination budget of the univariate factorizations.
const FACTORIZATION_SUBSETS: usize = 1 << 12;

/// A polynomial in `variables` unknowns with rational coefficients.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct MultivariatePolynomial {
    variables: usize,
    terms: BTreeMap<Monomial, Q>,
}

impl MultivariatePolynomial {
    pub(super) fn zero(variables: usize) -> Self {
        Self {
            variables,
            terms: BTreeMap::new(),
        }
    }

    /// Add `coefficient * x^monomial`.
    pub(super) fn add_term(&mut self, monomial: Monomial, coefficient: &Q) {
        debug_assert_eq!(monomial.len(), self.variables);
        if coefficient.is_zero() {
            return;
        }
        let entry = self.terms.entry(monomial.clone()).or_insert_with(Q::zero);
        *entry += coefficient;
        if entry.is_zero() {
            self.terms.remove(&monomial);
        }
    }

    pub(super) fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    fn leading(&self) -> Option<(&Monomial, &Q)> {
        self.terms.last_key_value()
    }

    fn is_constant(&self) -> bool {
        self.leading()
            .is_some_and(|(m, _)| m.iter().all(|&e| e == 0))
    }

    /// `self += coefficient * x^shift * other`.
    fn add_scaled_shifted(&mut self, other: &Self, coefficient: &Q, shift: &[u32]) {
        for (monomial, value) in &other.terms {
            let shifted = monomial.iter().zip(shift).map(|(a, b)| a + b).collect();
            self.add_term(shifted, &(value * coefficient));
        }
    }

    fn monic(mut self) -> Self {
        if let Some((_, lead)) = self.leading() {
            let inverse = lead.recip();
            for value in self.terms.values_mut() {
                *value *= &inverse;
            }
        }
        self
    }

    fn max_bits(&self) -> u64 {
        self.terms
            .values()
            .map(|c| c.numer().bits().max(c.denom().bits()))
            .max()
            .unwrap_or(0)
    }

    /// Substitute `value` for variable `index`.
    fn substitute(&self, index: usize, value: &Q) -> Self {
        let mut out = Self::zero(self.variables);
        for (monomial, coefficient) in &self.terms {
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
        for (monomial, coefficient) in &self.terms {
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

fn divides(a: &[u32], b: &[u32]) -> bool {
    a.iter().zip(b).all(|(x, y)| x <= y)
}

fn quotient_monomial(b: &[u32], a: &[u32]) -> Monomial {
    b.iter().zip(a).map(|(x, y)| x - y).collect()
}

fn lcm(a: &[u32], b: &[u32]) -> Monomial {
    a.iter().zip(b).map(|(x, y)| *x.max(y)).collect()
}

/// Why the solver stopped without an answer.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum SolverStop {
    /// The step, size or coefficient budget was reached.
    Budget(String),
    /// The system has infinitely many complex solutions.
    PositiveDimensional,
}

struct Budget {
    steps: usize,
    max_bits: u64,
}

impl Budget {
    fn spend(&mut self, polynomial: &MultivariatePolynomial) -> Result<(), SolverStop> {
        self.steps = self
            .steps
            .checked_sub(1)
            .ok_or_else(|| SolverStop::Budget("the Gröbner step budget was exhausted".into()))?;
        if polynomial.max_bits() > self.max_bits {
            return Err(SolverStop::Budget(
                "a Gröbner coefficient exceeded the bit budget".into(),
            ));
        }
        Ok(())
    }
}

/// Full reduction of `p` modulo `basis` (monic leading coefficients).
fn reduce(
    p: &MultivariatePolynomial,
    basis: &[MultivariatePolynomial],
    budget: &mut Budget,
) -> Result<MultivariatePolynomial, SolverStop> {
    let mut rest = p.clone();
    let mut remainder = MultivariatePolynomial::zero(p.variables);
    while let Some((monomial, coefficient)) = rest.leading().map(|(m, c)| (m.clone(), c.clone())) {
        match basis
            .iter()
            .find(|g| divides(g.leading().unwrap().0, &monomial))
        {
            Some(g) => {
                let (g_monomial, _) = g.leading().unwrap();
                let shift = quotient_monomial(&monomial, g_monomial);
                rest.add_scaled_shifted(g, &-coefficient, &shift);
                budget.spend(&rest)?;
            }
            None => {
                remainder.add_term(monomial.clone(), &coefficient);
                rest.terms.remove(&monomial);
            }
        }
    }
    Ok(remainder)
}

/// The reduced lex Gröbner basis of `polynomials` (Buchberger, with the
/// coprime-leading-monomial criterion).
fn groebner_basis(
    polynomials: &[MultivariatePolynomial],
    budget: &mut Budget,
) -> Result<Vec<MultivariatePolynomial>, SolverStop> {
    let mut basis: Vec<MultivariatePolynomial> = Vec::new();
    for p in polynomials {
        let r = reduce(p, &basis, budget)?;
        if !r.is_zero() {
            basis.push(r.monic());
        }
    }
    let mut pairs: Vec<(usize, usize)> = (0..basis.len())
        .flat_map(|j| (0..j).map(move |i| (i, j)))
        .collect();
    while !pairs.is_empty() {
        if basis.iter().any(MultivariatePolynomial::is_constant) {
            break;
        }
        // Normal strategy: the pair with the smallest lcm.
        let position = (0..pairs.len())
            .min_by(|&x, &y| {
                let key = |(i, j): (usize, usize)| {
                    let l = lcm(basis[i].leading().unwrap().0, basis[j].leading().unwrap().0);
                    (l.iter().sum::<u32>(), l)
                };
                key(pairs[x]).cmp(&key(pairs[y]))
            })
            .unwrap();
        let (i, j) = pairs.swap_remove(position);
        let (mi, _) = basis[i].leading().unwrap();
        let (mj, _) = basis[j].leading().unwrap();
        if mi.iter().zip(mj).all(|(a, b)| *a == 0 || *b == 0) {
            continue;
        }
        let l = lcm(mi, mj);
        let mut s = MultivariatePolynomial::zero(basis[i].variables);
        s.add_scaled_shifted(&basis[i], &Q::one(), &quotient_monomial(&l, mi));
        let shift_j = quotient_monomial(&l, mj);
        s.add_scaled_shifted(&basis[j], &-Q::one(), &shift_j);
        let r = reduce(&s, &basis, budget)?;
        if !r.is_zero() {
            let index = basis.len();
            basis.push(r.monic());
            pairs.extend((0..index).map(|i| (i, index)));
            if basis.len() > 400 {
                return Err(SolverStop::Budget(
                    "the Gröbner basis grew beyond 400 elements".into(),
                ));
            }
        }
    }
    if let Some(constant) = basis.iter().find(|p| p.is_constant()) {
        return Ok(vec![constant.clone().monic()]);
    }
    // Minimal, then reduced.
    let mut minimal: Vec<MultivariatePolynomial> = Vec::new();
    for (index, p) in basis.iter().enumerate() {
        let lead = p.leading().unwrap().0;
        let redundant = basis.iter().enumerate().any(|(other, q)| {
            let q_lead = q.leading().unwrap().0;
            other != index && divides(q_lead, lead) && (q_lead != lead || other < index)
        });
        if !redundant {
            minimal.push(p.clone());
        }
    }
    let mut reduced = Vec::new();
    for index in 0..minimal.len() {
        let others: Vec<MultivariatePolynomial> = minimal
            .iter()
            .enumerate()
            .filter(|(other, _)| *other != index)
            .map(|(_, p)| p.clone())
            .collect();
        let lead = minimal[index]
            .leading()
            .map(|(m, c)| (m.clone(), c.clone()))
            .unwrap();
        let mut tail = minimal[index].clone();
        tail.terms.remove(&lead.0);
        let mut p = reduce(&tail, &others, budget)?;
        p.add_term(lead.0, &lead.1);
        reduced.push(p.monic());
    }
    Ok(reduced)
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
    let mut budget = Budget {
        steps: max_steps,
        max_bits: 4096,
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
    budget: &mut Budget,
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
