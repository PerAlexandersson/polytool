# combinatoric-core

`combinatoric-core` is the foundational Rust library for exact and finite
combinatorics. It is independent of the `sym-poly` crates; `sym-poly-core`
builds on it, followed by `sym-poly-multipoly`, `sym-poly-sym`, and
`sym-poly-qsym`.

## Capability map

| Area | Main API |
| --- | --- |
| Shapes and words | `Partition`, `Composition`, `WeakComposition`, permutations, set partitions |
| Graphs | `Graph` generators, union/join operations, predicates, matchings, independence, graph6 |
| Connected-subgraph packings | `graph_packings::connected_subgraph_packing_polynomial`: exact BigInt coefficients by piece count, non-induced pieces bounded by edge count; at most 20 vertices and 1,000,000 candidate pieces |
| Posets | `Poset`, `Poset::from_permutation`, Young shapes, linear extensions, map counts, order-polytope Ehrhart and h* |
| Exact algebra | `SparseMatrix<BigInt>`, Smith normal form, chain complexes and homology |
| Polynomial models | GT/Kogan-face key polynomials and Ehrhart interpolation |
| Meanders | Noncrossing matchings and rooted meandric permutations |
| Colored exchange | `colored_exchange`: fixed-content families, cycle witnesses, exact weighted marginals, and checked decomposition certificates |

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

let counterexample = Graph::claw_free_schur_counterexample();
assert_eq!(counterexample.num_vertices(), 12);
assert!(counterexample.is_claw_free());
```

`Graph::claw_free_schur_counterexample` is the graph6 graph
``K?`CRAWWUXIM`` from J. Prajapati's counterexample to Schur positivity of
chromatic symmetric functions for claw-free graphs.  Its structural fixture
and exact Schur coefficient `-64` are regression-tested in the core and
`sym-poly-sym` crates.

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

`Poset::from_permutation(&[2, 1, 4, 3])` constructs the dominance order
`i <_P j` iff `i < j` and `perm[i] < perm[j]`, with redundant edges
removed. Permutation **values are one-indexed** (`1..=n`); poset vertices
are zero-indexed positions. Invalid permutations panic; empty input is valid.
Run `cargo run -p combinatoric-core --example permutation_poset` for a
checked P-Eulerian/height example. The constructor reuses the shared
permutation validator and Hasse reduction, without enumerating extensions.

The `colored_exchange` module validates finite families of words with fixed
color quotas. Its exchange checker returns a cycle witness or a specific
obstruction; its marginals use arbitrary-precision nonnegative coefficients.
Product, binary-port, and private-block certificates verify exact support
reconstruction and factor exchange. They do not assert a general Lorentzian
theorem. Run `cargo run -p combinatoric-core --example colored_exchange` for
a small checked use; broad support searches belong in downstream consumers.

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
