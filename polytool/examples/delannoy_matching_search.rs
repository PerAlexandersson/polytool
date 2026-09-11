//! Stream a graph6 search for matching-polynomial models of reversed
//! Delannoy-square rows.
//!
//! Isolated vertices do not change a matching polynomial, so exhaustive calls
//! should ask nauty for minimum degree one and scan every feasible order.

use combinatoric_core::Graph;
use num_bigint::BigInt;
use polytool::sequences::literature::delannoy_square_polynomials_bigint;
use serde_json::json;
use std::io::{self, BufRead};

fn target(rank: usize) -> Vec<BigInt> {
    delannoy_square_polynomials_bigint(rank)
        .pop()
        .expect("the requested rank must have a Delannoy-square row")
        .into_iter()
        .rev()
        .collect()
}

fn candidate_polynomial(graph: &Graph) -> Vec<BigInt> {
    graph
        .matching_polynomial()
        .into_iter()
        .map(BigInt::from)
        .collect()
}

fn two_matching_count(graph: &Graph) -> usize {
    let edge_pairs = graph.num_edges() * (graph.num_edges() - 1) / 2;
    let adjacent_pairs = (0..graph.num_vertices())
        .map(|vertex| {
            let degree = graph.degree(vertex);
            degree * degree.saturating_sub(1) / 2
        })
        .sum::<usize>();
    edge_pairs - adjacent_pairs
}

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args[0] == "--help" || args[0] == "-h" {
        eprintln!("Usage: delannoy_matching_search RANK [--witness-limit N] < graphs.g6");
        eprintln!("Input is streamed; blank lines and graph6 headers are skipped.");
        std::process::exit(if args.len() == 1 { 0 } else { 2 });
    }
    let rank = args[0]
        .parse::<usize>()
        .expect("RANK must be a positive integer");
    assert!(rank >= 1, "RANK must be positive");
    let witness_limit = match args.as_slice() {
        [_] => usize::MAX,
        [_, flag, limit] if flag == "--witness-limit" => limit
            .parse::<usize>()
            .expect("witness limit must be a nonnegative integer"),
        _ => {
            eprintln!("Unknown arguments; use --help");
            std::process::exit(2);
        }
    };
    let wanted = target(rank);
    let expected_edges = usize::try_from(&wanted[1]).expect("target edge count does not fit usize");
    assert!(
        expected_edges <= 62,
        "Graph::matching_polynomial uses i64; this search is guarded to at most 62 edges"
    );

    let mut input_graphs = 0_u64;
    let mut correct_edges = 0_u64;
    let mut without_isolates = 0_u64;
    let mut correct_two_matchings = 0_u64;
    let mut exact_polynomials = 0_u64;
    let mut connected_exact_polynomials = 0_u64;
    let mut witnesses = Vec::new();

    for (line_number, line) in io::stdin().lock().lines().enumerate() {
        let line = line.unwrap_or_else(|error| panic!("failed to read stdin: {error}"));
        let graph6 = line.trim();
        if graph6.is_empty() || graph6.starts_with('>') {
            continue;
        }
        input_graphs += 1;
        let graph = Graph::from_graph6(graph6).unwrap_or_else(|error| {
            panic!("invalid graph6 on input line {}: {error}", line_number + 1)
        });
        if graph.num_edges() != expected_edges {
            continue;
        }
        correct_edges += 1;
        if (0..graph.num_vertices()).any(|vertex| graph.degree(vertex) == 0) {
            continue;
        }
        without_isolates += 1;
        if let Some(two_matchings) = wanted.get(2) {
            if BigInt::from(two_matching_count(&graph)) != *two_matchings {
                continue;
            }
        }
        correct_two_matchings += 1;
        let polynomial = candidate_polynomial(&graph);
        if polynomial != wanted {
            continue;
        }
        exact_polynomials += 1;
        let connected = graph.is_connected();
        connected_exact_polynomials += u64::from(connected);
        if witnesses.len() < witness_limit {
            let mut degree_sequence = (0..graph.num_vertices())
                .map(|vertex| graph.degree(vertex))
                .collect::<Vec<_>>();
            degree_sequence.sort_unstable();
            witnesses.push(json!({
                "graph6": graph6,
                "vertices": graph.num_vertices(),
                "connected": connected,
                "degree_sequence": degree_sequence,
                "edges": graph.edges(),
                "matching_polynomial": polynomial.iter().map(ToString::to_string).collect::<Vec<_>>(),
            }));
        }
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "rank": rank,
            "target": wanted.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "expected_edges": expected_edges,
            "counts": {
                "input_graphs": input_graphs,
                "correct_edges": correct_edges,
                "without_isolates": without_isolates,
                "correct_two_matchings": correct_two_matchings,
                "exact_polynomials": exact_polynomials,
                "connected_exact_polynomials": connected_exact_polynomials,
            },
            "witnesses_omitted": exact_polynomials.saturating_sub(witnesses.len() as u64),
            "witnesses": witnesses,
        }))
        .expect("JSON serialization failed")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversed_targets_start_with_the_stated_rows() {
        assert_eq!(target(1), vec![BigInt::from(1), BigInt::from(2)]);
        assert_eq!(
            target(4),
            vec![1, 14, 61, 88, 29]
                .into_iter()
                .map(BigInt::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn two_disjoint_edges_calibrate_rank_one() {
        let path_three = Graph::path(3);
        assert_eq!(candidate_polynomial(&path_three), target(1));
    }
}
