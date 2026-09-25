//! Translation-invariant symmetric polynomials from graph monomials.
//!
//! Implements the partition-graph basis of Alexandersson--Shapiro,
//! *Discriminants, symmetrized graph monomials, and sums of squares*,
//! <https://arxiv.org/abs/1104.0589>, Section 2.2, in every degree.
//!
//! This is a basis of the translation-invariant subspace in **n variables**,
//! not an additional basis of all of `Sym`. The parameter `n` is essential:
//! its coefficients do not stabilize on adding zero variables. For degree
//! `d`, the paper works in the stable range `n >= 2*d`. Constructors also accept
//! smaller `n` when the graph fits, without asserting the basis property there.
//!
//! We use Definition 1 literally: sum over all `n!` permutations, with no
//! division by graph automorphisms or by `n!`. Results are monomial expansions
//! over `BigInt`; call `.to_schur_basis()` for Schur expansions. This conversion
//! reuses the crate's existing i64 transition matrices, so its documented
//! degree/coefficient limits still apply. Divide the resulting coefficients
//! by `n!` over `BigRational` if a permutation *average* is desired.
//!
//! The implementation expands the graph monomial only on incident vertices,
//! groups exponent vectors by partition, and multiplies each group by its
//! stabilizer `(n - length(mu))! * product_r multiplicity_mu(r)!`. It never
//! enumerates `S_n`. Sparse polynomial expansion can nevertheless be exponential
//! in the number of edges; this API is intended for small degrees.

use std::collections::{BTreeMap, BTreeSet};
use std::error::Error;
use std::fmt;

use num_bigint::BigInt;
use num_traits::{One, Zero};
use sym_poly_core::Partition;
use sym_poly_multipoly::MultiPoly;

use crate::{Basis, SymmetricFunction};

/// Invalid graph or partition input to a symmetrized graph monomial.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GraphMonomialError {
    /// A basis index contains a part smaller than two.
    NotTwoPartition,
    /// There are insufficient variables to label every graph vertex.
    TooFewVariables {
        /// Number of graph vertices, excluding added isolated vertices.
        required: usize,
        /// Requested number of variables.
        provided: usize,
    },
    /// An edge endpoint is not in `0..num_vars`.
    VertexOutOfRange {
        /// Invalid endpoint.
        vertex: usize,
        /// Requested number of variables.
        num_vars: usize,
    },
    /// The paper's graphs are loopless.
    Loop {
        /// The vertex at which a loop was supplied.
        vertex: usize,
    },
    /// The degree or vertex count exceeds the indexing type's range.
    SizeOverflow,
}

impl fmt::Display for GraphMonomialError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotTwoPartition => write!(f, "all parts of the index must be at least two"),
            Self::TooFewVariables { required, provided } => {
                write!(
                    f,
                    "the graph needs {required} variables, but only {provided} were supplied"
                )
            }
            Self::VertexOutOfRange { vertex, num_vars } => {
                write!(f, "vertex {vertex} is outside 0..{num_vars}")
            }
            Self::Loop { vertex } => write!(f, "a loop at vertex {vertex} is not allowed"),
            Self::SizeOverflow => write!(f, "graph degree or vertex count exceeds indexing range"),
        }
    }
}

impl Error for GraphMonomialError {}

/// Compute `sum_{sigma in S_n} product_{(u,v) in edges} (x_sigma(u) - x_sigma(v))`.
///
/// Vertices are zero-indexed in `0..num_vars`. Repeat an edge to give it
/// multiplicity; `(u,v)` contributes `x_u - x_v`. Isolated vertices are
/// implicit and contribute to the full permutation-sum normalization.
/// The empty edge list returns `num_vars!`, including `0! = 1`.
///
/// Returns a monomial-basis symmetric function whose specialization to
/// `num_vars` variables is the indicated polynomial. When `num_vars` is smaller
/// than the degree, this is the lift with no monomial terms of length exceeding
/// `num_vars`; discard Schur terms longer than `num_vars` for a finite-variable
/// Schur expansion.
///
/// Errors on loops, out-of-range endpoints, or a degree larger than `u32::MAX`.
pub fn symmetrized_graph_monomial(
    num_vars: usize,
    edges: &[(usize, usize)],
) -> Result<SymmetricFunction<BigInt>, GraphMonomialError> {
    u32::try_from(edges.len()).map_err(|_| GraphMonomialError::SizeOverflow)?;
    let mut vertices = BTreeSet::new();
    for &(u, v) in edges {
        for vertex in [u, v] {
            if vertex >= num_vars {
                return Err(GraphMonomialError::VertexOutOfRange { vertex, num_vars });
            }
            vertices.insert(vertex);
        }
        if u == v {
            return Err(GraphMonomialError::Loop { vertex: u });
        }
    }
    let indices: BTreeMap<_, _> = vertices
        .into_iter()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    let active = indices.len();
    let mut polynomial = MultiPoly::constant(active, BigInt::one());
    for &(u, v) in edges {
        let difference =
            MultiPoly::<BigInt>::var(active, indices[&u]) - MultiPoly::var(active, indices[&v]);
        polynomial = polynomial * difference;
    }

    let mut grouped: BTreeMap<Partition, BigInt> = BTreeMap::new();
    for (exponents, coefficient) in polynomial.terms() {
        *grouped
            .entry(Partition::new(exponents.clone()))
            .or_default() += coefficient;
    }
    grouped.retain(|_, coefficient| !coefficient.is_zero());
    for (partition, coefficient) in &mut grouped {
        *coefficient *= factorial(num_vars - partition.num_parts()) * part_stabilizer(partition);
    }
    Ok(SymmetricFunction::from_terms(Basis::Monomial, grouped))
}

/// The partition-graph basis element `b_lambda` of Section 2.2, in `num_vars` variables.
///
/// Each part `a >= 2` supplies a disjoint outward-oriented star with `a` leaves.
/// Thus the graph has `|lambda|` edges and `|lambda| + length(lambda)` vertices.
/// Uses the full permutation sum; for example in four variables,
/// `b_(2) = 6*m_(2) - 4*m_(1,1) = 6*s_(2) - 10*s_(1,1)`.
///
/// All parts must be at least two, and the graph must fit in `num_vars`.
/// The empty partition is accepted and gives `num_vars!`.
/// Together these elements form the paper's degree-`d` basis when `num_vars >= 2*d`.
#[doc(alias = "b_lambda")]
pub fn partition_graph_symmetric(
    partition: &Partition,
    num_vars: usize,
) -> Result<SymmetricFunction<BigInt>, GraphMonomialError> {
    let degree = two_partition_degree(partition)?;
    let required = degree
        .checked_add(partition.num_parts())
        .ok_or(GraphMonomialError::SizeOverflow)?;
    check_variables(required, num_vars)?;
    let mut edges = Vec::with_capacity(degree);
    let mut next = 0;
    for &part in partition.parts() {
        let center = next;
        next += 1;
        for _ in 0..part {
            edges.push((center, next));
            next += 1;
        }
    }
    symmetrized_graph_monomial(num_vars, &edges)
}

fn two_partition_degree(partition: &Partition) -> Result<usize, GraphMonomialError> {
    if partition.parts().iter().any(|&part| part < 2) {
        return Err(GraphMonomialError::NotTwoPartition);
    }
    let degree = partition
        .parts()
        .iter()
        .try_fold(0u32, |sum, part| sum.checked_add(*part))
        .ok_or(GraphMonomialError::SizeOverflow)?;
    usize::try_from(degree).map_err(|_| GraphMonomialError::SizeOverflow)
}

fn check_variables(required: usize, provided: usize) -> Result<(), GraphMonomialError> {
    if required > provided {
        Err(GraphMonomialError::TooFewVariables { required, provided })
    } else {
        Ok(())
    }
}

fn factorial(n: usize) -> BigInt {
    (1..=n).map(BigInt::from).product()
}

fn part_stabilizer(partition: &Partition) -> BigInt {
    let mut counts = BTreeMap::<u32, usize>::new();
    for &part in partition.parts() {
        *counts.entry(part).or_default() += 1;
    }
    counts.values().map(|&count| factorial(count)).product()
}
