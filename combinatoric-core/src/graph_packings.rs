//! Vertex-disjoint packings of connected, not necessarily induced edge-subgraphs.

use crate::Graph;
use num_bigint::BigInt;
use std::collections::{BTreeMap, HashMap};

/// Counts packings by their number of pieces, in ascending coefficient order.
/// Each piece has between one and `max_edges` edges. Distinct edge subsets
/// count separately even when their vertex supports agree. Unused vertices
/// are allowed, isolated vertices are not pieces, and the empty packing has
/// weight one. `max_edges == 0` therefore returns `[1]`.
///
/// Uses exact `BigInt` coefficients, edge-subset enumeration, and memoized
/// vertex-subset deletion. Exponential: at most `2^n` states, each considering
/// all distinct piece supports. For bounded exploratory use, returns an error
/// above 20 vertices or after 1,000,000 candidate edge subsets; no partial
/// polynomial is returned. The edge bound is capped at the graph's edge count.
///
/// ```
/// use combinatoric_core::{Graph, graph_packings::connected_subgraph_packing_polynomial};
/// use num_bigint::BigInt;
/// let p = connected_subgraph_packing_polynomial(&Graph::path(4), 2).unwrap();
/// assert_eq!(p, vec![BigInt::from(1), BigInt::from(5), BigInt::from(1)]);
/// ```
pub fn connected_subgraph_packing_polynomial(
    graph: &Graph,
    max_edges: usize,
) -> Result<Vec<BigInt>, String> {
    if graph.num_vertices() > 20 {
        return Err("connected-subgraph packings: limit is 20 vertices".into());
    }
    let mut supports = BTreeMap::new();
    enumerate_pieces(
        graph.edges(),
        max_edges.min(graph.num_edges()),
        0,
        &mut Vec::new(),
        &mut 0,
        &mut supports,
    )?;
    let pieces: Vec<_> = supports.into_iter().collect();
    Ok(packing_dp(
        (1u64 << graph.num_vertices()) - 1,
        &pieces,
        &mut HashMap::new(),
    ))
}

fn enumerate_pieces(
    edges: &[(usize, usize)],
    bound: usize,
    start: usize,
    chosen: &mut Vec<(usize, usize)>,
    examined: &mut usize,
    supports: &mut BTreeMap<u64, BigInt>,
) -> Result<(), String> {
    if !chosen.is_empty() {
        *examined += 1;
        if *examined > 1_000_000 {
            return Err("connected-subgraph packings: candidate limit exceeded".into());
        }
        let support = chosen
            .iter()
            .fold(0u64, |s, &(a, b)| s | (1 << a) | (1 << b));
        let mut reached = 1u64 << chosen[0].0;
        loop {
            let old = reached;
            for &(a, b) in chosen.iter() {
                let pair = (1u64 << a) | (1u64 << b);
                if reached & pair != 0 {
                    reached |= pair;
                }
            }
            if reached == old {
                break;
            }
        }
        if reached == support {
            *supports.entry(support).or_insert_with(|| BigInt::from(0)) += 1;
        }
    }
    if chosen.len() == bound {
        return Ok(());
    }
    for i in start..edges.len() {
        chosen.push(edges[i]);
        enumerate_pieces(edges, bound, i + 1, chosen, examined, supports)?;
        chosen.pop();
    }
    Ok(())
}

fn packing_dp(
    vertices: u64,
    pieces: &[(u64, BigInt)],
    memo: &mut HashMap<u64, Vec<BigInt>>,
) -> Vec<BigInt> {
    if vertices == 0 {
        return vec![BigInt::from(1)];
    }
    if let Some(p) = memo.get(&vertices) {
        return p.clone();
    }
    let pivot = 1u64 << vertices.trailing_zeros();
    let mut out = packing_dp(vertices ^ pivot, pieces, memo);
    for (support, multiplicity) in pieces {
        if support & pivot != 0 && support & vertices == *support {
            let tail = packing_dp(vertices ^ support, pieces, memo);
            out.resize(out.len().max(tail.len() + 1), BigInt::from(0));
            for (i, c) in tail.iter().enumerate() {
                out[i + 1] += c * multiplicity;
            }
        }
    }
    memo.insert(vertices, out.clone());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // Independent oracle: enumerate selected edges, then inspect their
    // connected components. Each permissible selection is exactly one packing.
    fn edge_subset_oracle(g: &Graph, k: usize) -> Vec<BigInt> {
        let mut out = vec![BigInt::from(0); g.num_vertices() / 2 + 1];
        for mask in 0usize..(1usize << g.num_edges()) {
            let mut parent: Vec<_> = (0..g.num_vertices()).collect();
            fn root(p: &[usize], mut v: usize) -> usize {
                while p[v] != v {
                    v = p[v];
                }
                v
            }
            for (i, &(a, b)) in g.edges().iter().enumerate() {
                if mask & (1 << i) != 0 {
                    let ra = root(&parent, a);
                    let rb = root(&parent, b);
                    parent[ra] = rb;
                }
            }
            let mut sizes = vec![0; g.num_vertices()];
            for (i, &(a, _)) in g.edges().iter().enumerate() {
                if mask & (1 << i) != 0 {
                    sizes[root(&parent, a)] += 1;
                }
            }
            if sizes.iter().all(|&s| s <= k) {
                out[sizes.iter().filter(|&&s| s > 0).count()] += 1;
            }
        }
        while out.len() > 1 && out.last() == Some(&BigInt::from(0)) {
            out.pop();
        }
        out
    }

    #[test]
    fn conventions_and_multiplicity() {
        for (g, k, want) in [
            (Graph::path(0), 2, vec![1]),
            (Graph::path(4), 0, vec![1]),
            (Graph::path(3), 2, vec![1, 3]),
            (Graph::path(4), 2, vec![1, 5, 1]),
            (Graph::complete(3), 2, vec![1, 6]),
            (Graph::complete(3), 3, vec![1, 7]),
        ] {
            assert_eq!(
                connected_subgraph_packing_polynomial(&g, k).unwrap(),
                want.into_iter().map(BigInt::from).collect::<Vec<_>>()
            );
        }
        assert!(connected_subgraph_packing_polynomial(&Graph::path(21), 2).is_err());
    }

    #[test]
    fn exhaustive_small_edge_subset_oracle() {
        for n in 0..=4 {
            let all = Graph::complete(n);
            for mask in 0usize..(1 << all.num_edges()) {
                let edges: Vec<_> = all
                    .edges()
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) != 0)
                    .map(|(_, e)| *e)
                    .collect();
                let g = Graph::new(n, &edges);
                for k in 0..=g.num_edges() + 1 {
                    assert_eq!(
                        connected_subgraph_packing_polynomial(&g, k).unwrap(),
                        edge_subset_oracle(&g, k)
                    );
                }
                assert_eq!(
                    connected_subgraph_packing_polynomial(&g, 1).unwrap(),
                    g.matching_polynomial()
                        .into_iter()
                        .map(BigInt::from)
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[test]
    fn twelve_vertex_tree_independent_oracle() {
        let g = Graph::new(
            12,
            &[
                (0, 6),
                (1, 7),
                (2, 8),
                (3, 9),
                (0, 10),
                (1, 10),
                (4, 10),
                (2, 11),
                (3, 11),
                (5, 11),
                (10, 11),
            ],
        );
        for (k, expected) in [
            (2, vec![1, 27, 142, 212, 117, 21, 1]),
            (3, vec![1, 56, 323, 410, 183, 26, 1]),
        ] {
            let expected: Vec<_> = expected.into_iter().map(BigInt::from).collect();
            assert_eq!(
                connected_subgraph_packing_polynomial(&g, k).unwrap(),
                expected
            );
            assert_eq!(edge_subset_oracle(&g, k), expected);
        }
    }
}
