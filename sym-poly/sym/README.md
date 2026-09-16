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
| Graph/chromatic families | `chromatic`, `weighted_bond`, q-chromatic and Hessenberg/GKM routines |
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
