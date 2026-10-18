# exact-real-roots

Exact real-root counting and real-rootedness tests for dense univariate
polynomials with integer coefficients.

The crate accepts coefficient slices in ascending degree order and uses
arbitrary-precision integer arithmetic. It handles repeated roots by counting
the distinct roots of the square-free part. The implementation combines exact
low-degree formulas, signed pseudo-remainder sequences, and an independent
Uspensky/Descartes path.

```rust
use exact_real_roots::{
    count_real_roots_prs_bigint_coeffs, is_real_rooted_fast_bigint_coeffs,
};
use num_bigint::BigInt;

// (x + 1)^3
let polynomial: Vec<BigInt> = [1, 3, 3, 1].into_iter().map(BigInt::from).collect();
assert!(is_real_rooted_fast_bigint_coeffs(&polynomial));
assert_eq!(count_real_roots_prs_bigint_coeffs(&polynomial), 1);
```

`polytool` re-exports this API from `polytool::root_count` and its crate root.
