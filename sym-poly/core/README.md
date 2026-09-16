# sym-poly-core

`sym-poly-core` is the shared foundation for the symmetric and
quasisymmetric-function crates and for `sym-poly-multipoly`. It depends on
`combinatoric-core` for the canonical partition, composition, ring, and sparse
matrix adapters; it does not depend on the downstream algebra crates.

## Capability map

| Area | Main API or module |
| --- | --- |
| Coefficients | `Ring`, `Field`, `PrimeField`, `UnivariatePolynomial`, rational functions |
| Indices and sums | `Partition`, `Composition`, `WeakComposition`, `BasisIndex`, `FormalSum` |
| Tableaux/fillings | `Tableau`, `SkewTableau`, standard iterators, `Ssaf`, P–RS insertion |
| Exact matrices | `matrix`, `linear_algebra`, sparse/packed modular row reduction |
| Symmetric-function helpers | `TransitionCache`, Hamel–Goulden cutting strips, CRT, finite `S_n` modules |

`UnivariatePolynomial` stores coefficients in ascending degree order. Partition
and composition conventions are inherited from `combinatoric-core`; weak
composition length remains structural for variables and polynomial bases.

## Source-verified example

The checked `examples/hamel_goulden_site_example.rs` constructs an outside
decomposition and its determinant matrix:

```rust
use sym_poly_core::{ContentInterval, CuttingStripSegment, OutsideDecomposition};

let decomposition = OutsideDecomposition::from_intervals([(-3, 2), (-1, 1)]);
let matrix = decomposition.determinant_matrix();
assert_eq!(matrix[0][1],
    CuttingStripSegment::Segment(ContentInterval::new(-1, 2)));
```

The full checked output is available with
`cargo run -p sym-poly-core --example hamel_goulden_site_example`.

## Exactness and limits

`Ring` supports i64, BigInt, integer/rational coefficients, and polynomial
layers, while exact division in integer rings rejects a non-exact quotient.
The general `matrix` and `TransitionCache` APIs store matrices as
`Vec<Vec<i64>>`; sparse integral linear algebra and certificate workflows can
use `BigInt`, but this does not make current i64 transition/Kostka data
arbitrary precision. `PrimeField` requires a prime const modulus. Exponent and
degree fields remain bounded `u32` values, and constructors/operators may panic
on malformed shapes, zero denominators, or checked overflow.

Run `cargo test -p sym-poly-core` for focused checks. Public additions should
carry rustdoc, invariant/boundary tests, and a small checked example where
useful. Search with `rg` before adding another index, tableau, matrix, or
coefficient algorithm; put one-use experiments outside the repository.
