# sym-poly-sym

`sym-poly-sym` implements generic symmetric functions in the six classical
bases—monomial, elementary, complete homogeneous, power sum, Schur, and
forgotten. It depends on `sym-poly-core` for rings, partitions, tableaux, and
matrix machinery, and on `sym-poly-multipoly` for selected nonsymmetric/GKM
support. `sym-poly-qsym` is downstream of this crate.

## Capability map

| Area | Main API or module |
| --- | --- |
| Core algebra | `SymmetricFunction<C>`, `Basis`, all basis conversions, products, omega, Hall product |
| Tableaux and characters | `hook_schur`, Kostka/character/Frobenius, shifted LR |
| Graph/chromatic families | `chromatic`, universal-vertex extension, `weighted_bond`, q-chromatic and Hessenberg/GKM routines |
| Translation-invariant graph basis | `partition_graph_symmetric`, `symmetrized_graph_monomial`; `examples/partition_graph_schur.rs` |
| Chordal clique boundaries | `chordal_clique_boundary_schur_channels`; `examples/clique_boundary_channels.rs` |
| q,t families | `llt`, `macdonald`, Petrie and shifted/Grothendieck routines |
| Additional constructions | `lah`, twin GKM, circular/affine helpers |

Functions are indexed by normalized `Partition` values. Use
`SymmetricFunction::<C>::schur_symmetric(partition)` (or another named
constructor), then call `to_*_basis`; multiplication requires both operands
to use the same basis.

## Source-verified uses

The checked `examples/boolean_product_site_example.rs` verifies the Schur
conversion for `(2,1)`:

```rust
use sym_poly_core::Partition;
use sym_poly_sym::SymmetricFunction;

let schur = SymmetricFunction::<i64>::schur_symmetric(Partition::new(vec![2, 1]));
let monomial = schur.to_monomial_basis();
assert_eq!(monomial.coefficient(&Partition::new(vec![2, 1])), 1);
assert_eq!(monomial.coefficient(&Partition::new(vec![1, 1, 1])), 2);
```

For a q-polynomial-valued family, run the checked LLT example with
`cargo run -p sym-poly-sym --example llt_site_example`.
The checked `chromatic_universal_vertex` example verifies the monomial-basis
universal-vertex operator against `Graph::join`.

## Partition-graph basis

`partition_graph_symmetric(&lambda, n)` implements the all-degree partition-graph
basis from Alexandersson--Shapiro, [arXiv:1104.0589, Section 2.2](https://arxiv.org/abs/1104.0589).
Indices have every part at least two. Each part supplies a disjoint outward
star. These span the translation-invariant degree-`d` subspace in `n` variables
in the paper's stable range `n >= 2*d`; they are not a basis of all symmetric
functions. The explicit `n` affects the coefficients. The result uses `BigInt`
and the full permutation sum of Definition 1, including isolated vertices.
Call `.to_schur_basis()` to expand it. For example, `b_(2)` at `n=4` is
`6*s_(2) - 10*s_(1,1)`, and `b_(3)` at `n=6` is
`120*s_(3) - 192*s_(2,1) + 336*s_(1,1,1)`.

```bash
cargo run -p sym-poly-sym --example partition_graph_schur
cargo run -p sym-poly-sym --example partition_graph_schur -- 5 10
```

The first command lists degrees 2--6 at `n=2*d`. Symmetrization groups sparse
graph-monomial terms by exponent partition and includes the stabilizer
factor `(n-length(mu))! * product_r multiplicity_mu(r)!`, avoiding `n!`
permutation enumeration. This factor is necessary for literal Definition 1
normalization (the coloring count in Lemma 6 omits the zero-exponent factorial).
Polynomial expansion can still be exponential in the edge count.

## Exactness and limits

The generic coefficient parameter can be `i64`, `BigInt`, `Ratio<i64>`,
`Ratio<BigInt>`, or compatible `Ring` layers such as `UnivariatePolynomial`.
Power-sum conversions use exact division: integer coefficients panic when a
required quotient is not integral, while rational coefficients retain it.
Current transition matrices and several Kostka/character helpers are stored or
computed through i64 data, so selecting `BigInt` does not make those paths
arbitrary precision. Partition sizes/degrees are `u32`; basis conversion,
tableau enumeration, and q,t families can grow rapidly and may hit machine
integer or allocation limits.

Run focused tests with `cargo test -p sym-poly-sym`, and keep paper/site output
in small checked examples. New public families need module rustdoc, source or
paper examples, and tests for normalization, basis conversion, and known small
values. Search with `rg` before adding another transition, Kostka, tableau, or
character algorithm; one-use experiments belong outside the repository.
