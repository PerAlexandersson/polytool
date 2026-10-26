# exact-real-roots

Exact real-root counting and real-rootedness tests for dense univariate
polynomials with integer coefficients.

The crate accepts coefficient slices in ascending degree order and uses
arbitrary-precision integer arithmetic. It handles repeated roots by counting
the distinct roots of the square-free part. The implementation combines exact
low-degree formulas, a signed Collins--Brown subresultant pseudo-remainder
sequence, and an independent Uspensky/Descartes path. The subresultant path
divides out predictable exact factors instead of computing a coefficient GCD
after every pseudo-division. Root-counting calls stream the sequence through
two active polynomials and one reusable scratch buffer; they retain only the
endpoint sign variations and final gcd degree needed by Sturm's theorem.

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

`isolate_real_roots_bigint_coeffs` isolates the distinct real roots by Sturm
bisection: each root is returned either exactly (a `DyadicRational`) or in an
open interval with dyadic endpoints that are not roots, sorted increasingly.

`polytool` re-exports this API from `polytool::root_count` (and the most-used
items from its crate root).

## Algorithm lineage

The subresultant scaling is a clean Rust implementation of the classical
Collins--Brown recurrence described by Brown and Traub and by Cohen. FLINT's
public root-counting documentation is used as a design and performance oracle;
this crate does not translate or incorporate FLINT source code.

- W. S. Brown and J. F. Traub, “On Euclid's Algorithm and the Theory of
  Subresultants,” *JACM* 18 (1971), 505–514,
  <https://doi.org/10.1145/321662.321665>.
- FLINT `fmpz_poly` root-counting documentation,
  <https://flintlib.org/doc/fmpz_poly.html#roots>.
