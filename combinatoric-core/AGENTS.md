# combinatoric-core guide

## Role and boundaries

`combinatoric-core` is the workspace's foundational combinatorics crate. It
has no dependency on the `sym-poly/*` crates. The current dependency direction
is:

```text
combinatoric-core -> sym-poly-core -> sym-poly-multipoly -> sym-poly-sym -> sym-poly-qsym
```

The arrows mean “is used by”; `sym-poly-sym` also uses `sym-poly-multipoly`.
`combpoly` and other workspace crates are separate projects with their own
dependency edges. Do not put
symmetric-function algorithms here merely because they use partitions.

The stable boundary is the public API in `src/lib.rs` and its modules:

- `partition`, `composition`, `permutation`, and `set_partition` provide the
  canonical indexing and enumeration types;
- `graph` provides simple graphs, named graph generators, graph predicates,
  matchings, independence data, and graph6 parsing;
- `poset` provides Hasse-diagram posets, linear extensions, order-preserving
  maps, order polytopes, h*-vectors, and P-Eulerian polynomials;
- `chain_complex`, `sparse_matrix`, and `integer_linear_algebra` provide
  exact integral chain-complex, sparse-matrix, Smith-form, and certificate
  workflows;
- `key_polynomial` provides GT-pattern/Kogan-face key-polynomial and Ehrhart
  computations; `meander` provides noncrossing matchings and rooted meanders;
- `ring` defines the coefficient-ring contract used by downstream polynomial
  crates.

## Conventions and numerical limits

- `Partition::new` sorts parts weakly decreasing and removes zeros;
  `from_sorted` trusts the caller. `Composition::new` removes trailing zeros;
  `WeakComposition` preserves its length because it can encode variables.
- Partition/composition parts and degrees are `u32`; graph, poset, and
  permutation labels are zero-based unless a function explicitly documents a
  one-based convention. Meandric output is intentionally one-based.
- `Partition::count_syt` returns `BigUint`; sparse matrices, Smith forms,
  chain complexes, and order-polytope rational calculations use `BigInt` or
  `BigRational` where the implementation promises exactness.
- Several established combinatorial polynomial/counting APIs still return
  `i64`, `u64`, `u128`, or `usize`. Those are bounded APIs and can panic or
  overflow at their documented machine-type limits. A `BigInt` input or a
  downstream `Ring<BigInt>` does not upgrade current i64 transition/Kostka or
  other machine-integer data to arbitrary precision.
- The Smith reducer is deliberately dense and budgeted; use sparse cancellation
  first and pass explicit limits for residual computations. Exact division in
  integer coefficient rings rejects a non-divisible result.

## Change and verification rules

Search all Rust crates with `rg` before adding an algorithm or public name.
Extend one implementation and, where a narrow wrapper is useful, have it call
the exact implementation; do not maintain duplicate combinatorial algorithms.
Experiments are deliberately disposable and are not part of this crate's API.
Do not add one-use experiment binaries here in place of tests or examples.

Public functions and types should have rustdoc describing indexing, ordering,
normalization, complexity/limits, and failure behavior. Add focused unit tests
for invariants and boundary cases, plus a small `examples/` program when an
API is intended for a paper or website. Keep examples source-verified and
small. Run focused tests before a workspace-wide check, and keep generated
build output outside the Dropbox checkout according to the workspace guide.
