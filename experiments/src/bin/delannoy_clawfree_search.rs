//! Stream an isomorph-free graph6 search for reversed Delannoy-square rows.
//!
//! Typical use with nauty is
//! `nauty-geng -q -l -F 10 20:20 | cargo run -q -p experiments --bin
//! delannoy_clawfree_search -- 3`.

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
        .independence_polynomial()
        .into_iter()
        .map(BigInt::from)
        .collect()
}

fn expected_order(rank: usize) -> usize {
    rank.checked_mul(4)
        .and_then(|value| value.checked_sub(2))
        .expect("rank is too large")
}

fn main() {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args[0] == "--help" || args[0] == "-h" {
        eprintln!("Usage: delannoy_clawfree_search RANK [--witness-limit N] < graphs.g6");
        eprintln!("Input is streamed; blank lines and graph6 headers are skipped.");
        std::process::exit(if args.len() == 1 { 0 } else { 2 });
    }
    let rank = args[0]
        .parse::<usize>()
        .expect("RANK must be a positive integer");
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
    assert!(rank >= 1, "RANK must be positive");
    let order = expected_order(rank);
    assert!(
        order <= 62,
        "Graph::independence_polynomial uses i64; this search is guarded to at most 62 vertices"
    );
    let wanted = target(rank);
    let independent_pairs = wanted.get(2).cloned().unwrap_or_default();
    let possible_edges = order * (order - 1) / 2;
    let expected_edges = possible_edges
        .checked_sub(
            usize::try_from(&independent_pairs)
                .expect("target independent-pair count does not fit usize"),
        )
        .expect("target has too many independent pairs");

    let mut input_graphs = 0_u64;
    let mut correct_order = 0_u64;
    let mut correct_edges = 0_u64;
    let mut claw_free = 0_u64;
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
        if graph.num_vertices() != order {
            continue;
        }
        correct_order += 1;
        if graph.num_edges() != expected_edges {
            continue;
        }
        correct_edges += 1;
        if !graph.is_claw_free() {
            continue;
        }
        claw_free += 1;
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
                "connected": connected,
                "degree_sequence": degree_sequence,
                "edges": graph.edges(),
                "independence_polynomial": polynomial.iter().map(ToString::to_string).collect::<Vec<_>>(),
            }));
        }
    }

    println!(
        "{}",
        serde_json::to_string_pretty(&json!({
            "rank": rank,
            "target": wanted.iter().map(ToString::to_string).collect::<Vec<_>>(),
            "expected_vertices": order,
            "expected_edges": expected_edges,
            "counts": {
                "input_graphs": input_graphs,
                "correct_order": correct_order,
                "correct_edges": correct_edges,
                "claw_free": claw_free,
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
            target(3),
            vec![1, 10, 25, 12]
                .into_iter()
                .map(BigInt::from)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn small_witnesses_have_the_target_polynomials() {
        assert_eq!(candidate_polynomial(&Graph::complete(2)), target(1));

        let wheel_six = Graph::new(
            6,
            &[
                (0, 1),
                (0, 2),
                (0, 3),
                (0, 4),
                (0, 5),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (1, 5),
            ],
        );
        assert!(wheel_six.is_claw_free());
        assert_eq!(candidate_polynomial(&wheel_six), target(2));
    }
}
