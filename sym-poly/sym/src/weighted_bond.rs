//! Möbius symmetric functions from the bond lattice of a graph.
//!
//! The implementation follows the connected-partition recurrence used for
//! the weighted bond symmetric functions of González D'León--Wachs. Set
//! partitions are enumerated exactly, so these routines are intended for
//! small graphs; their cost grows at least with the relevant Bell number.

use std::collections::BTreeMap;
use std::error::Error;
use std::fmt;

use combinatoric_core::{set_partitions, Graph, SetPartition};
use sym_poly_core::{Partition, Ring};

use crate::{Basis, SymmetricFunction};

type GraphKey = (usize, Vec<(usize, usize)>);

/// An invalid input to a weighted bond symmetric-function computation.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum WeightedBondError {
    /// The connected-graph recurrence for `M_G` was requested on a
    /// disconnected graph.
    DisconnectedGraph,
}

impl fmt::Display for WeightedBondError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DisconnectedGraph => {
                write!(
                    formatter,
                    "the graph Möbius symmetric function requires a connected graph"
                )
            }
        }
    }
}

impl Error for WeightedBondError {}

/// Enumerate set partitions whose blocks induce connected subgraphs.
///
/// These are the elements of the graph's bond lattice, represented as
/// canonical set partitions of `{0, ..., graph.num_vertices() - 1}`.
pub fn bond_set_partitions(graph: &Graph) -> impl Iterator<Item = SetPartition> + '_ {
    set_partitions(graph.num_vertices()).filter(|partition| {
        partition
            .blocks()
            .iter()
            .all(|block| graph.induced_subgraph(block).is_connected())
    })
}

fn one<C: Ring>() -> SymmetricFunction<C> {
    SymmetricFunction::complete_h_symmetric(Partition::empty())
}

fn graph_key(graph: &Graph) -> GraphKey {
    let mut edges = graph.edges().to_vec();
    edges.sort_unstable();
    (graph.num_vertices(), edges)
}

fn mobius_connected<C: Ring>(
    graph: &Graph,
    cache: &mut BTreeMap<GraphKey, SymmetricFunction<C>>,
) -> SymmetricFunction<C> {
    let key = graph_key(graph);
    if let Some(function) = cache.get(&key) {
        return function.clone();
    }
    if graph.num_vertices() == 1 {
        return one();
    }

    let mut sum = SymmetricFunction::zero(Basis::CompleteH);
    for partition in bond_set_partitions(graph).filter(|partition| partition.num_blocks() > 1) {
        let block_count = u32::try_from(partition.num_blocks() - 1)
            .expect("a materialized graph has fewer than u32::MAX vertices");
        let complete = SymmetricFunction::complete_h_symmetric(Partition::new(vec![block_count]));
        let product = partition.blocks().iter().fold(one(), |product, block| {
            product.multiply(&mobius_connected(&graph.induced_subgraph(block), cache))
        });
        sum = sum + complete.multiply(&product);
    }
    let result = -sum;
    cache.insert(key, result.clone());
    result
}

/// Compute the graph Möbius symmetric function `M_G` of a connected graph.
///
/// The result is returned in the complete homogeneous basis over coefficient
/// ring `C`. Prefer `BigInt` when coefficient growth is not independently
/// bounded. The single-vertex graph maps to `1`; disconnected and empty graphs
/// return [`WeightedBondError::DisconnectedGraph`].
#[doc(alias = "M_G")]
pub fn graph_mobius_symmetric_function<C: Ring>(
    graph: &Graph,
) -> Result<SymmetricFunction<C>, WeightedBondError> {
    if graph.num_vertices() == 0 || !graph.is_connected() {
        return Err(WeightedBondError::DisconnectedGraph);
    }
    Ok(mobius_connected(graph, &mut BTreeMap::new()))
}

/// Compute the chromatic Möbius symmetric function `Psi_G`.
///
/// This sums, over connected set partitions of `G`, the product of the graph
/// Möbius symmetric functions of the induced blocks. The result is in the
/// complete homogeneous basis. The empty graph maps to `1`.
#[doc(alias = "Psi_G")]
pub fn chromatic_mobius_symmetric_function<C: Ring>(graph: &Graph) -> SymmetricFunction<C> {
    let mut cache = BTreeMap::new();
    bond_set_partitions(graph).fold(
        SymmetricFunction::zero(Basis::CompleteH),
        |sum, partition| {
            let product = partition.blocks().iter().fold(one(), |product, block| {
                product.multiply(&mobius_connected(
                    &graph.induced_subgraph(block),
                    &mut cache,
                ))
            });
            sum + product
        },
    )
}

#[cfg(test)]
mod tests {
    use num_bigint::BigInt;

    use super::*;

    fn coefficient(function: &SymmetricFunction<BigInt>, parts: &[u32]) -> BigInt {
        function.coefficient(&Partition::new(parts.to_vec()))
    }

    #[test]
    fn path_and_complete_graph_match_published_small_examples() {
        let path = Graph::path(3);
        let complete = Graph::complete(3);

        let m_path = graph_mobius_symmetric_function::<BigInt>(&path)
            .unwrap()
            .to_monomial_basis();
        let m_complete = graph_mobius_symmetric_function::<BigInt>(&complete)
            .unwrap()
            .to_monomial_basis();
        assert_eq!(coefficient(&m_path, &[2]), BigInt::from(1));
        assert_eq!(coefficient(&m_path, &[1, 1]), BigInt::from(3));
        assert_eq!(coefficient(&m_complete, &[2]), BigInt::from(2));
        assert_eq!(coefficient(&m_complete, &[1, 1]), BigInt::from(5));

        let psi_path = chromatic_mobius_symmetric_function::<BigInt>(&path).to_elementary_basis();
        let psi_complete =
            chromatic_mobius_symmetric_function::<BigInt>(&complete).to_elementary_basis();
        for (parts, expected) in [
            (&[][..], 1),
            (&[1][..], -2),
            (&[2][..], 1),
            (&[1, 1][..], 1),
        ] {
            assert_eq!(coefficient(&psi_path, parts), BigInt::from(expected));
        }
        for (parts, expected) in [
            (&[][..], 1),
            (&[1][..], -3),
            (&[2][..], 1),
            (&[1, 1][..], 2),
        ] {
            assert_eq!(coefficient(&psi_complete, parts), BigInt::from(expected));
        }
    }

    #[test]
    fn validates_connected_input_and_handles_empty_graph() {
        let disconnected = Graph::empty(2);
        assert_eq!(
            graph_mobius_symmetric_function::<BigInt>(&disconnected),
            Err(WeightedBondError::DisconnectedGraph)
        );
        assert_eq!(
            graph_mobius_symmetric_function::<BigInt>(&Graph::empty(0)),
            Err(WeightedBondError::DisconnectedGraph)
        );
        assert_eq!(
            chromatic_mobius_symmetric_function::<BigInt>(&Graph::empty(0)),
            one()
        );
    }

    #[test]
    fn bond_partitions_cover_boundary_and_small_graphs() {
        assert_eq!(bond_set_partitions(&Graph::empty(0)).count(), 1);
        assert_eq!(bond_set_partitions(&Graph::empty(1)).count(), 1);
        assert_eq!(bond_set_partitions(&Graph::path(3)).count(), 4);
        assert_eq!(bond_set_partitions(&Graph::complete(3)).count(), 5);
    }
}
