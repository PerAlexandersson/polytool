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
| Schur content operators | `c_n`, `c_n_inverse`, `schur_content_product`; `examples/schur_content_operator.rs` |
| Tableaux and characters | `hook_schur`, Kostka/character/Frobenius, shifted LR |
| Graph/chromatic families | `chromatic`, universal-vertex extension, `weighted_bond`, q-chromatic and Hessenberg/GKM routines |
| Translation-invariant graph bases | `partition_graph_symmetric`, `square_partition_graph_symmetric`, `symmetrized_graph_monomial`; `examples/partition_graph_schur.rs`, `examples/square_graph_w.rs` |
| Chordal clique boundaries | `chordal_clique_boundary_schur_channels`; `examples/clique_boundary_channels.rs` |
| q,t families | `llt`, including `LltColoringStandardization`; `macdonald`, Petrie and shifted/Grothendieck routines |
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
`llt_standardize_coloring` records the LLT-specific decreasing-vertex
tie standardization, its rank-color word, and the mandatory strict-rise set;
the certificate reconstructs the original coloring and works over any ordered
color alphabet.
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

The even-degree square basis from Section 2.3 is available as
`square_partition_graph_symmetric(&lambda, n)`. An even part gives a star
with half as many leaves and doubled edges; odd parts in decreasing order are
paired consecutively into two stars joined by a doubled edge. This includes
indices such as `(3,3)`, not only partitions with even parts. The returned
`BigInt` monomial expansion uses the full permutation sum, not the square of
an already symmetrized polynomial. The graph needs
`degree/2 + length(lambda) - number_of_odd_parts/2` vertices; odd degree and
parts equal to one are rejected. The empty partition gives `n!`.

`cargo run -p sym-poly-sym --example square_graph_w -- 4 8` prints both the
direct forgotten expansion and the expansion in `W_mu = C_n^{-1}(f_mu)` for
each degree-four square basis element in eight variables, normalized by `n!`.
For W coordinates use `c_n(&graph, n).to_forgotten_basis()`; for direct
forgotten coordinates use `graph.to_forgotten_basis()`. Neither expansion is
uniformly positive: see [the checked examples and sign evidence](../../docs/square-graph-expansions.md).

## Schur content operators

`c_n(&value, n)` multiplies each Schur coefficient indexed by `lambda` by
`product_(i,j in lambda)(n+j-i)`, using one-based row/column indices.
`schur_content_product(&lambda, n)` returns this eigenvalue as a `BigInt`.
The forward map accepts any supported coefficient ring; `BigInt` and
`BigRational` avoid fixed-width overflow in these multiplications.

`c_n_inverse(&value, n)` accepts `BigInt` or `BigRational` input and returns
`Result<SymmetricFunction<BigRational>, ContentOperatorError>`. Both operators
accept any classical basis and return Schur expansions. For example,
`C_3(s_(2,1)) = 24*s_(2,1)` and `C_3^(-1)(s_(2,1)) = s_(2,1)/24`.
Integer graph polynomials can be supplied directly to the inverse; the checked
example gives `C_4^(-1)(b_(2)^(4)) = 3/10*s_(2) - 5/6*s_(1,1)`:

```bash
cargo run -p sym-poly-sym --example schur_content_operator
```

The inverse's domain is explicit: a nonzero Schur term with more than `n` rows
is an error because its eigenvalue is zero. When working in exactly `n`
variables, first convert to Schur and discard such terms (they specialize to
zero). The forward map annihilates them. Constants are unchanged, including
at `n=0`; on the admissible Schur span both compositions are the identity.
This inverse is over the rationals, not generally over the integers. Conversion
from non-Schur input retains the transition-table limits described below.

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
