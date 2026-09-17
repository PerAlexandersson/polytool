//! Exact Pieri channels for chordal graphs with a marked clique boundary.
//!
//! Let `R` be a clique in a chordal graph `G`, put `I = V(G) \ R`, and let
//! `a` be Pawlowski's relative perfect-elimination-order operator, normalized
//! so the full operator is `A_R a`.  For `mu` partitioning `|I|`, this module
//! computes the cleared channel
//!
//! ```text
//! F_hat(G,R,mu) = Ch(A_R a z_mu),
//! z_mu = sum_{sigma in S_I} chi^mu(sigma^{-1}) sigma.
//! ```
//!
//! Thus `F_hat = |I|! / f^mu * F`, where `F` is the primitive Pieri-channel
//! function.  If `t_(lambda,mu)` is the eigenvalue of the relative compression
//! `P_R a P_R` on the corresponding Pieri summand, then
//!
//! ```text
//! [s_lambda] F_hat(G,R,mu) = |R|! |I|! t_(lambda,mu).
//! ```
//!
//! Some group-algebra checkers instead report the eigenvalue of the full
//! operator `A_R a = |R|! P_R a`; for that convention the right side is
//! `|I|!` times the reported eigenvalue.
//!
//! In particular, the cleared channels have exactly the signs needed for
//! clique-boundary positivity while using integral arithmetic throughout.

use std::collections::{BTreeMap, HashMap};
use std::fmt;

use combinatoric_core::{all_permutations_zero_indexed, Graph, Partition};
use num_bigint::BigInt;
use sym_poly_core::sn_action::cycle_type;

use crate::kostka::sn_character;
use crate::{Basis, SymmetricFunction};

/// Largest graph order accepted by the exact factorial-time channel routine.
pub const MAX_EXACT_CLIQUE_BOUNDARY_ORDER: usize = 9;

/// Failure modes for exact clique-boundary channel computation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CliqueBoundaryError {
    /// The factorial-time computation was requested above its explicit cap.
    OrderTooLarge { order: usize, maximum: usize },
    /// No perfect elimination order ends in the supplied boundary.
    NoPerfectEliminationOrder,
}

impl fmt::Display for CliqueBoundaryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::OrderTooLarge { order, maximum } => write!(
                formatter,
                "clique-boundary computation has order {order}, above the exact cap {maximum}"
            ),
            Self::NoPerfectEliminationOrder => {
                write!(
                    formatter,
                    "no perfect elimination order ends in the boundary"
                )
            }
        }
    }
}

impl std::error::Error for CliqueBoundaryError {}

type GroupElement = HashMap<Vec<usize>, i64>;

/// Compute all cleared clique-boundary Pieri channels in the Schur basis.
///
/// The map key is the interior partition `mu`; the value is
/// `F_hat(G,R,mu)` in the Schur basis.  The boundary order affects the chosen
/// elimination certificate but not the resulting channels.  Returns an error
/// unless the boundary is a clique admitting a perfect elimination order at
/// the end.  The algorithm is exact and factorial-time, and is deliberately
/// capped by [`MAX_EXACT_CLIQUE_BOUNDARY_ORDER`].
pub fn chordal_clique_boundary_schur_channels(
    graph: &Graph,
    boundary: &[usize],
) -> Result<BTreeMap<Partition, SymmetricFunction<BigInt>>, CliqueBoundaryError> {
    let order = graph.num_vertices();
    if order > MAX_EXACT_CLIQUE_BOUNDARY_ORDER {
        return Err(CliqueBoundaryError::OrderTooLarge {
            order,
            maximum: MAX_EXACT_CLIQUE_BOUNDARY_ORDER,
        });
    }
    let peo = graph
        .perfect_elimination_order_ending_with(boundary)
        .ok_or(CliqueBoundaryError::NoPerfectEliminationOrder)?;
    let boundary_set = boundary
        .iter()
        .copied()
        .collect::<std::collections::BTreeSet<_>>();
    let interior = (0..order)
        .filter(|vertex| !boundary_set.contains(vertex))
        .collect::<Vec<_>>();

    let relative = relative_chordal_operator(graph, boundary.len(), &peo);
    let full = left_antisymmetrize(relative, boundary, order);
    let interior_shapes = Partition::all_of_size(interior.len() as u32);
    let outer_shapes = Partition::all_of_size(order as u32);
    let interior_cycle_indices = interior_shapes
        .iter()
        .enumerate()
        .map(|(index, shape)| (cycle_type_code(shape), index))
        .collect::<HashMap<_, _>>();
    let outer_cycle_indices = outer_shapes
        .iter()
        .enumerate()
        .map(|(index, shape)| (cycle_type_code(shape), index))
        .collect::<HashMap<_, _>>();
    let local_permutations = all_permutations_zero_indexed(interior.len())
        .into_iter()
        .map(|permutation| {
            let cycle_index = interior_cycle_indices[&permutation_cycle_type_code(&permutation)];
            (permutation, cycle_index)
        })
        .collect::<Vec<_>>();
    let mut interior_position = vec![usize::MAX; order];
    for (index, &vertex) in interior.iter().enumerate() {
        interior_position[vertex] = index;
    }
    let mut cycle_pair_weights = vec![vec![0i64; outer_shapes.len()]; interior_shapes.len()];
    for (permutation, coefficient) in &full {
        for (local, interior_cycle_index) in &local_permutations {
            let outer_code = right_composition_cycle_type_code(
                permutation,
                local,
                &interior,
                &interior_position,
            );
            let outer_cycle_index = outer_cycle_indices[&outer_code];
            let entry = &mut cycle_pair_weights[*interior_cycle_index][outer_cycle_index];
            *entry = entry
                .checked_add(*coefficient)
                .expect("order-nine cycle-pair weight fits in i64");
        }
    }

    let inner_characters = interior_shapes
        .iter()
        .map(|shape| {
            interior_shapes
                .iter()
                .map(|cycle| sn_character(shape, cycle))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let outer_characters = outer_shapes
        .iter()
        .map(|shape| {
            outer_shapes
                .iter()
                .map(|cycle| sn_character(shape, cycle))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut channels = BTreeMap::new();
    for (interior_shape_index, interior_shape) in interior_shapes.iter().enumerate() {
        let mut schur_terms = BTreeMap::new();
        for (outer_shape_index, outer_shape) in outer_shapes.iter().enumerate() {
            if !is_vertical_strip(outer_shape, interior_shape, boundary.len()) {
                continue;
            }
            let mut coefficient = 0i128;
            for (interior_cycle_index, row) in cycle_pair_weights.iter().enumerate() {
                let inner_character = inner_characters[interior_shape_index][interior_cycle_index];
                if inner_character == 0 {
                    continue;
                }
                for (outer_cycle_index, &weight) in row.iter().enumerate() {
                    let outer_character = outer_characters[outer_shape_index][outer_cycle_index];
                    if weight != 0 && outer_character != 0 {
                        let contribution = i128::from(weight)
                            * i128::from(inner_character)
                            * i128::from(outer_character);
                        coefficient = coefficient
                            .checked_add(contribution)
                            .expect("order-nine channel coefficient fits in i128");
                    }
                }
            }
            if coefficient != 0 {
                schur_terms.insert(outer_shape.clone(), BigInt::from(coefficient));
            }
        }
        channels.insert(
            interior_shape.clone(),
            SymmetricFunction::from_terms(Basis::Schur, schur_terms),
        );
    }
    Ok(channels)
}

fn cycle_type_code(shape: &Partition) -> u64 {
    let mut code = 0u64;
    for &part in shape.parts() {
        code += 1u64 << (4 * part);
    }
    code
}

fn permutation_cycle_type_code(permutation: &[usize]) -> u64 {
    let mut seen = 0u16;
    let mut code = 0u64;
    for start in 0..permutation.len() {
        if seen & (1u16 << start) != 0 {
            continue;
        }
        let mut current = start;
        let mut length = 0usize;
        while seen & (1u16 << current) == 0 {
            seen |= 1u16 << current;
            length += 1;
            current = permutation[current];
        }
        code += 1u64 << (4 * length);
    }
    code
}

fn right_composition_cycle_type_code(
    permutation: &[usize],
    local: &[usize],
    interior: &[usize],
    interior_position: &[usize],
) -> u64 {
    let mut seen = 0u16;
    let mut code = 0u64;
    for start in 0..permutation.len() {
        if seen & (1u16 << start) != 0 {
            continue;
        }
        let mut current = start;
        let mut length = 0usize;
        while seen & (1u16 << current) == 0 {
            seen |= 1u16 << current;
            length += 1;
            let local_position = interior_position[current];
            let right_image = if local_position == usize::MAX {
                current
            } else {
                interior[local[local_position]]
            };
            current = permutation[right_image];
        }
        code += 1u64 << (4 * length);
    }
    code
}

fn is_vertical_strip(outer: &Partition, inner: &Partition, size: usize) -> bool {
    outer.size() == inner.size() + size as u32
        && (0..outer.num_parts().max(inner.num_parts())).all(|row| {
            let difference = outer.part(row).checked_sub(inner.part(row));
            matches!(difference, Some(0 | 1))
        })
}

fn relative_chordal_operator(graph: &Graph, boundary_size: usize, peo: &[usize]) -> GroupElement {
    let order = graph.num_vertices();
    let identity = (0..order).collect::<Vec<_>>();
    let mut positions = vec![0usize; order];
    for (position, &vertex) in peo.iter().enumerate() {
        positions[vertex] = position;
    }
    let mut result = HashMap::from([(identity, 1i64)]);
    for &vertex in peo[..order - boundary_size].iter().rev() {
        let later_neighbors = graph
            .neighbors(vertex)
            .iter()
            .copied()
            .filter(|&neighbor| positions[neighbor] > positions[vertex])
            .collect::<Vec<_>>();
        let mut next = HashMap::with_capacity(result.len() * (later_neighbors.len() + 1));
        for (permutation, coefficient) in result {
            add_term(&mut next, permutation.clone(), coefficient);
            for &neighbor in &later_neighbors {
                let mut product = permutation.clone();
                product.swap(vertex, neighbor);
                add_term(&mut next, product, -coefficient);
            }
        }
        next.retain(|_, coefficient| *coefficient != 0);
        result = next;
    }
    result
}

fn left_antisymmetrize(element: GroupElement, boundary: &[usize], order: usize) -> GroupElement {
    let mut result = HashMap::new();
    for local in all_permutations_zero_indexed(boundary.len()) {
        let sign = if ((boundary.len() - cycle_type(&local).num_parts()) & 1) == 0 {
            1i64
        } else {
            -1i64
        };
        let mut embedded = (0..order).collect::<Vec<_>>();
        for (index, &vertex) in boundary.iter().enumerate() {
            embedded[vertex] = boundary[local[index]];
        }
        for (permutation, coefficient) in &element {
            let product = permutation
                .iter()
                .map(|&value| embedded[value])
                .collect::<Vec<_>>();
            add_term(&mut result, product, sign * coefficient);
        }
    }
    result.retain(|_, coefficient| *coefficient != 0);
    result
}

fn add_term(element: &mut GroupElement, permutation: Vec<usize>, coefficient: i64) {
    let entry = element.entry(permutation).or_insert(0);
    *entry = entry
        .checked_add(coefficient)
        .expect("order-nine group-algebra coefficient fits in i64");
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::Signed;

    fn p(parts: &[u32]) -> Partition {
        Partition::new(parts.to_vec())
    }

    #[test]
    fn edge_with_singleton_boundary_has_one_channel() {
        let channels = chordal_clique_boundary_schur_channels(&Graph::complete(2), &[1]).unwrap();
        let channel = &channels[&p(&[1])];
        assert_eq!(channel.coefficient(&p(&[1, 1])), BigInt::from(2));
        assert_eq!(channel.terms().len(), 1);
    }

    #[test]
    fn terminal_edge_of_path_four_matches_cleared_channel_identity() {
        let graph = Graph::path(4);
        let channels = chordal_clique_boundary_schur_channels(&graph, &[2, 3]).unwrap();
        assert!(channels[&p(&[2])].is_zero());

        let expected = crate::chromatic_symmetric::<BigInt>(&graph).to_schur_basis();
        let sign_channel = &channels[&p(&[1, 1])];
        for (shape, coefficient) in expected.terms() {
            assert_eq!(sign_channel.coefficient(shape), coefficient * 2);
        }
        assert_eq!(sign_channel.terms().len(), expected.terms().len());
    }

    #[test]
    fn outer_vertex_of_complete_three_sun_has_nonnegative_channels() {
        let graph = Graph::complete_three_sun_clique_blowup([1, 1, 1], [1, 1, 1]);
        let channels = chordal_clique_boundary_schur_channels(&graph, &[5]).unwrap();
        assert!(channels.values().all(|channel| channel
            .terms()
            .values()
            .all(|coefficient| !coefficient.is_negative())));
    }

    #[test]
    fn rejects_nonchordal_boundary_problem_and_large_order() {
        assert_eq!(
            chordal_clique_boundary_schur_channels(&Graph::cycle(4), &[0]),
            Err(CliqueBoundaryError::NoPerfectEliminationOrder)
        );
        assert_eq!(
            chordal_clique_boundary_schur_channels(&Graph::empty(10), &[]),
            Err(CliqueBoundaryError::OrderTooLarge {
                order: 10,
                maximum: MAX_EXACT_CLIQUE_BOUNDARY_ORDER,
            })
        );
    }
}
