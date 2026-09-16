# Experiment promotion registry

This is a small queue of reusable mathematical kernels found during research.
It is not an inventory of `experiments/src/bin`: most local experiments are
deliberately ignored and disposable. A local-only source path may be absent in
a fresh clone; its entry is a pointer for this workstation until the kernel is
reviewed and promoted into tracked library code.

Statuses are `audit`, `ready`, `in progress`, `promoted`, or `rejected`.
Promotion work must follow `experiments/AGENTS.md` and update this table in the
same checkpoint.

## Current queue

| Priority | Status | Candidate | Proposed owner | Required before promotion |
|---|---|---|---|---|
| P1 | audit | Tracked `experiments/src/matroids.rs`: `BasisMatroid`, `PrefixLpm`, transversal and contingency routines | `combinatoric-core::matroid`, with a possible `matroid::contingency` submodule | Decide whether construction validates basis exchange or exposes an explicitly unchecked type; document labels and exhaustive costs; retain independent tests. |
| P1 | audit, local-only | `derangement_fixed_point_slices.rs`: fixed-point/descent bivariate recurrence | `polytool::sequences::permutation_statistics` | Fix and test the index-zero convention; use `BigInt`; preserve a small brute-force permutation oracle. |
| P1 | audit, local-only | `nn_rook_qdeform.rs`: nesting-refined rook polynomial | `combpoly::rook_placements` | Handle the empty board without underflow; use `BigInt`; stream placements; verify the `q=0` and `q=1` specializations against canonical generators. |
| P2 | promoted | Tracked `sym-poly/sym/examples/weighted_bond_symmetric_site_example.rs`: connected bond partitions and Möbius symmetric functions | `sym_poly_sym::weighted_bond` | Promoted with the foundation set-partition iterator, generic exact coefficients, graph-key memoization, connected-input and Bell-number cost documentation, focused tests, and the retained example as a consumer. |
| P2 | audit, local-only | `matroid_contingency_scan.rs`: orbit-compressed Catalan core-walk coefficient | future matroid contingency module | Replace unchecked `usize` multiplicities, explain the orbit argument, and cross-check the general basis enumerator. |
| P2 | audit, local-only | `edge_triangle_word_sink_polynomials.rs`: frontier transfer for edge/triangle words | specialized `combinatoric-core` graph-family module | Demonstrate a useful advantage over existing general and chordal sink APIs; document the word convention and compare against `Graph` on small cases. |

## Known non-candidates

The following functionality is already canonical and should be reused rather
than promoted again:

- non-nesting rook generators and packet utilities in
  `combpoly::rook_placements`;
- general, chordal, frozen-edge, and tree sink polynomials in
  `combinatoric_core::Graph`;
- sparse matrices, bounded Smith reduction, and finite chain complexes in
  `combinatoric-core`;
- deco sequences in `polytool::sequences::deco`;
- literature Eulerian/Delannoy/Hoggatt sequences in
  `polytool::sequences::literature`;
- linear-extension promotion and orbits in `combinatoric_core::Poset`.

`experiments/src/nn_rook_utils.rs` is compatibility plumbing over the
maintained rook library, not another implementation to promote. Principal-
stratum cell generation remains research-specific; only its generic sparse
algebra and homology machinery belongs in the foundation library, where it is
already maintained.

## Promotion record template

For a new row, record:

```text
candidate:
source availability: tracked | local-only
mathematical contract:
proposed crate/module/API:
exact arithmetic and indexing:
independent oracle:
durable consumers:
research provenance:
status and owner:
```

After promotion, replace local copies in durable consumers, link the source
commit from the associated research note, and keep only genuinely exploratory
drivers in `experiments`.
