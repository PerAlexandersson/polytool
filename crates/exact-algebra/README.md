# exact-algebra

Exact integer and rational arithmetic kernels, independent of polynomial types.
The crate requires Rust 1.82 and depends only on the `num-*` crates.
It is available in-tree and is not currently published to crates.io.

| Module | Purpose |
| --- | --- |
| `primes` | Deterministic word-size primality validation and descending prime streams |
| `modular_arithmetic` | Overflow-safe word-size residue operations |
| `chinese_remainder` | Batched BigInt CRT with signed, squared-Hadamard and nonnegative bounds |
| `rational_reconstruction` | Bounded BigInt rational reconstruction with congruence and uniqueness checks |
| `denominators` | Row-wise denominator clearing and positive row scales |
| `integer_elimination` | Internal fraction-free engine; optional `test-hooks` fallback guard |
| `linear_systems` | Exact solves, rank/nullity and certified rational kernels |
| `dense_modular` | CRT reconstruction of all leading principal minors |

Matrices are row-major `Vec<Vec<_>>`; kernel basis vectors are stored separately.
All indices are zero-based. Reconstruction requires a justified bound on the true
integer or rational target; a modulus threshold alone cannot validate modular input.
Sparse elimination and higher-level quotient adapters remain separate consumers.

```rust
use exact_algebra::chinese_remainder::ChineseRemainderAccumulator;
use exact_algebra::primes::ValidatedPrime;

let mut reconstruction = ChineseRemainderAccumulator::new(1);
reconstruction.combine_prime(&[3], ValidatedPrime::new(5)?)?;
reconstruction.combine_prime(&[5], ValidatedPrime::new(7)?)?;
assert_eq!(reconstruction.reconstruct_symmetric(&[2.into()])?, vec![num_bigint::BigInt::from(-2)]);
# Ok::<(), Box<dyn std::error::Error>>(())
```

The optional `test-hooks` feature exposes a thread-local `ForcedFallback` guard
for exercising rational fallback paths. It is absent from normal builds.
