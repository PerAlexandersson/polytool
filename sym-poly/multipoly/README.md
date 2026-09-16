# sym-poly-multipoly

`sym-poly-multipoly` is the sparse multivariate-polynomial layer. It depends
on `sym-poly-core` and `combinatoric-core`, and is used by `sym-poly-sym` for
selected constructions; it does not depend on the symmetric-function crate.

## Capability map

| Area | Main API or module |
| --- | --- |
| Sparse algebra | `MultiPoly<C>`, `MultiPolyFunction`, monomial orders, division and normal forms |
| Operators | Divided differences, Demazure/θ/π and t-deformed words |
| Polynomial families | Key, atom, Schubert, flagged Schur, slide/glide, lock, Kohnert/Lascoux, Grothendieck |
| q,t families | Nonsymmetric Macdonald filling formula and q=0 Hall–Littlewood specialization |
| Algebraic computation | Gröbner bases, modular reconstruction, finite quotients, indexed-variable `S_n` actions, Lorentzian checks |
| Combinatorial models | SSAF fillings, multiline queues, Borodin–Wheeler weights |

`MultiPoly` uses exponent vectors of a fixed variable count; vectors must have
the same length, and exponents and total degrees are checked `u32` values.
Weak compositions preserve row/variable length. Follow each module's stated
basement, reading, permutation, monomial-order, and operator-index convention.

## Source-verified use

The checked `examples/key_site_example.rs` compares a key polynomial with its
SSAF key-fillings weight counts:

```rust
use std::collections::BTreeMap;
use sym_poly_core::Ssaf;
use sym_poly_multipoly::key_polynomial;

let alpha = vec![1, 0, 2];
let polynomial = key_polynomial::<i64>(&alpha);
let fillings = Ssaf::key_fillings(&alpha);
let mut counts = BTreeMap::new();
for filling in &fillings {
    *counts.entry(filling.weight_vector()).or_insert(0) += 1;
}
assert_eq!(polynomial.terms(), &counts);
```

Run the checked flagged-Schur example with
`cargo run -p sym-poly-multipoly --example flagged_schur_site_example`.

## Exactness and limits

Polynomial coefficients are generic over `Ring`, including i64, BigInt, and
rational/polynomial layers. This does not make every algorithm arbitrary
precision: transition matrices remain i64-backed (as does Kostka data in the
sibling `sym` crate), and a BigInt coefficient choice cannot upgrade them. Exponent
addition and total-degree accumulation panic on `u32` overflow; malformed
variable lengths, zero divisors, singular matrices, and unsupported quotient
conditions likewise fail according to their module contracts. Gröbner,
quotient, filling, and basis-enumeration routines are intended for bounded
exact computations, not unrestricted symbolic workloads.

Run focused checks with `cargo test -p sym-poly-multipoly`. Public additions
need rustdoc, tests for operator/basis identities and boundary cases, and a
small source-verified example when a named family is exposed. Search with `rg`
before adding another divided-difference, filling, transition, or polynomial
generator; keep one-use experiments outside the repository and maintain one
canonical algorithm with ergonomic wrappers where needed.
