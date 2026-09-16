# combinatoric-core

`combinatoric-core` is the foundational Rust library for exact and finite
combinatorics. It is independent of the `sym-poly` crates; `sym-poly-core`
builds on it, followed by `sym-poly-multipoly`, `sym-poly-sym`, and
`sym-poly-qsym`.

## Capability map

| Area | Main API |
| --- | --- |
| Shapes and words | `Partition`, `Composition`, `WeakComposition`, permutations, set partitions |
| Graphs | `Graph` generators, predicates, matchings, independence, graph6 |
| Posets | `Poset`, linear extensions, map counts, order-polytope Ehrhart and h* |
| Exact algebra | `SparseMatrix<BigInt>`, Smith normal form, chain complexes and homology |
| Polynomial models | GT/Kogan-face key polynomials and Ehrhart interpolation |
| Meanders | Noncrossing matchings and rooted meandric permutations |

Partitions are stored in weakly decreasing order; compositions preserve order;
weak compositions preserve their length. Labels are normally zero-based, while
the meander permutation output is explicitly one-based. See the module
rustdoc for algorithm-specific conventions.

## Small checked uses

The graph module documents and tests this matching-polynomial example:

```rust
use combinatoric_core::graph::Graph;

let k4 = Graph::complete(4);
assert_eq!(k4.num_vertices(), 4);
assert_eq!(k4.num_edges(), 6);
assert_eq!(k4.matching_polynomial(), vec![1, 6, 3]);
```

The poset module tests the diamond throughout its map-count and order-polytope
APIs:

```rust
use combinatoric_core::poset::Poset;

let diamond = Poset::new(4, &[(0, 1), (0, 2), (1, 3), (2, 3)]);
assert_eq!(diamond.num_linear_extensions(), 2);
assert_eq!(diamond.count_weak_order_preserving_dp(5), 105);
assert_eq!(diamond.order_polytope_hstar(), vec![1, 1]);
```

Run the crate's focused checks with `cargo test -p combinatoric-core` (the
workspace guide supplies the timeout/nice wrapper for resource-sensitive runs).

## Exactness and scope

`count_syt` uses `BigUint`, order-polytope Ehrhart calculations use exact
`BigRational`, and chain-complex/Smith workflows use `BigInt` with explicit
resource limits and replayable certificates. Other public counts and
coefficient vectors intentionally use `i64`, `u64`, `u128`, or `usize`; they
remain bounded. In particular, choosing `BigInt` downstream does not make
current i64 transition/Kostka data arbitrary precision. Prefer one exact
implementation plus a checked narrow wrapper when adding a new counting API.

Before adding code, search the workspace with `rg` for an existing generator,
counter, or normalization helper. Keep one-use experiments outside the crate;
promote only a stable, documented algorithm with tests and a source-verified
example.
