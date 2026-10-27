# Library usage

Add to your `Cargo.toml`:

```toml
[dependencies]
polytool = { path = "../polytool" }
```

See the [README](../README.md#using-as-a-crate) for a short quick start.  This page
collects the longer library examples; the command-line interface is described in
[cli.md](cli.md).

## Standalone exact root counting

The low-level integer root-counting engine is the independent
`exact-real-roots` crate under `crates/exact-real-roots`. Polytool depends on
it and preserves the existing `polytool::root_count` module and crate-root
re-exports. Projects that only need exact root counts or coefficient-slice
real-rootedness can depend directly on the smaller crate:

```toml
[dependencies]
exact-real-roots = { path = "path/to/polytool/crates/exact-real-roots" }
```

Both crates use ascending-degree coefficient order and arbitrary-precision
integer arithmetic. Repeated roots are handled through the square-free part.

## Real-rootedness and interlacing

```rust
use polytool::*;

// Coefficients in ascending degree order: coeffs[i] = coeff of t^i
let eulerian_4 = [1, 11, 11, 1];

// Default exact PRS/root-counting method
assert!(is_real_rooted(&eulerian_4));

// Sturm chain method (slower, gives root locations)
assert!(is_real_rooted_sturm(&eulerian_4));

// Strict interlacing q ≪ p: the smaller-degree polynomial comes first
let p = [-15, 23, -9, 1];  // (t-1)(t-3)(t-5)
let q = [8, -6, 1];         // (t-2)(t-4)
assert_eq!(check_interlacing(&q, &p), Some(true));
assert_eq!(check_interlacing(&p, &q), None); // wrong order for these degrees

// Weak interlacing g ≪ f (allows shared roots)
let f = [2, -3, 1];   // (t-1)(t-2)
let g = [-1, 1];       // (t-1)
assert_eq!(check_weak_interlacing(&g, &f), Some(true));

// Same-degree interlacing: roots alternate on the real line
// (t-1)(t-3) roots {1,3} and (t-2)(t-4) roots {2,4}: 1 < 2 < 3 < 4
assert_eq!(check_weak_interlacing(&[3, -4, 1], &[8, -6, 1]), Some(true));

// Nested roots do NOT interlace: (t-1)(t-4) vs (t-2)(t-3)
assert_eq!(check_weak_interlacing(&[4, -5, 1], &[6, -5, 1]), Some(false));

// Undirected check: true iff the pair interlaces in some order allowed by
// the degrees; the zero polynomial interlaces every real-rooted polynomial.
assert!(interlaces_either_order(&g, &f));
assert!(interlaces_either_order(&f, &g));
assert!(interlaces_either_order(&[0], &f));
assert!(!interlaces_either_order(&[4, -5, 1], &[6, -5, 1]));
```

For large exact coefficients, use the `BigInt` APIs:

```rust
use num_bigint::BigInt;
use polytool::check_interlacing_bigint_coeffs;

let center = BigInt::from(10).pow(20);
let f = vec![-&center, BigInt::from(1)]; // t - center
let g = vec![
    (&center - 1u32) * (&center + 1u32),
    -BigInt::from(2) * &center,
    BigInt::from(1),
];

assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
```

Root-location helpers are named for what they return.  `approximate_real_roots`
(formerly `real_roots`) returns interval midpoints, rational approximations
accurate to below `1e-6`, and `None` when the polynomial is not real-rooted;
`real_rootedness::isolate_real_roots_rational_coeffs` (formerly
`isolate_real_roots_exact`) returns certified enclosures with multiplicities for
rational coefficients.  Total nonnegativity is tested by
`linalg::is_totally_nonnegative_neville` (all minors, Neville elimination) and
`linalg::is_totally_nonnegative_by_minors` (brute force up to a minor size).
Failures of the diagnostic `check_*`, parsing and report functions are typed
errors that implement `Display` and `std::error::Error`: `TotalPositivityError`,
`InterlacingMatrixError`, `CyclicSievingError`, `ParsePolynomialError`,
`NetworkError`, `oeis::OeisCatalogError`, `renewal::RenewalError`, and, for the
recurrence reports, searches and certificate replays, `recurrence::RecurrenceError`
(with `recurrence::ParseRationalError` for coefficient strings).  Public
functions never return `Result<_, String>`; the `Display` text is the
human-readable message.  Functions that cannot fail on well-formed input but
reject a precondition by panicking (zero divisors, non-square matrices, `i64`
overflow of a convenience wrapper) say so in a `# Panics` rustdoc section.

## Athanasiadis--Wagner interlacing matrices

```rust
use polytool::*;

// A column vector of polynomials P, Q, R in ascending coefficient order.
let polys = vec![
    vec![2, 1],       // 2 + t
    vec![8, 6, 1],    // (2+t)(4+t)
    vec![3, 4, 1],    // (1+t)(3+t)
];

// One finite truncation of the infinite Lace(A) matrix.
let lace = lace_matrix_sequence_i64(&polys, 1, 3).unwrap();
assert_eq!(lace, vec![vec![2, 1, 0], vec![8, 6, 1], vec![3, 4, 1]]);

// This is Athanasiadis--Wagner's pairwise-but-not-fully-interlacing example:
// the 3 x 3 determinant is negative.
assert!(check_lace_sequence_totally_nonnegative_by_minors_i64(&polys, 1, 3, 3).is_err());
```

For a `p x q` polynomial matrix `A`, `lace_matrix_i64(A, r, c)` returns the
`p r` by `q c` truncation whose `(p rb+i, q cb+j)` entry is the coefficient
of `x^(cb-rb)` in `A_ij(x)`.  This follows Definition 3.5 of
Athanasiadis--Wagner, *Veronese sections and interlacing matrices of formal
power series*.  A finite TNN check is computational evidence for the infinite
matrix, not a proof of full interlacing unless a separate finite criterion
applies. See [`doc/athanasiadis-wagner-lace.md`](athanasiadis-wagner-lace.md)
for further notes on the construction.

## Planar networks, chain polynomials and strip-digraph recognition

```rust
use num_bigint::BigInt;
use polytool::*;

let to_rows = |rows: Vec<Vec<i64>>| -> Vec<Vec<BigInt>> {
    rows.into_iter()
        .map(|r| r.into_iter().map(BigInt::from).collect())
        .collect()
};

// TNN check with a structured diagnosis: rows {1,2}, columns {0,1} give -2.
let bad = to_rows(vec![vec![1], vec![1, 1], vec![2, 0, 1]]);
match check_lower_unitriangular_totally_nonnegative(&bad) {
    Err(NetworkError::NotTotallyNonnegative { witness: Some(minor), .. }) => {
        assert_eq!(minor.value, BigInt::from(-2));
    }
    other => panic!("{:?}", other),
}

// Chain polynomials p_0 = 1, p_n = t * sum_{k<n} r_{n,k} p_k.
let ones = to_rows(vec![vec![1], vec![1, 1], vec![1, 1, 1]]);
let report = check_chain_polynomials_bigint(&ones).unwrap();
assert!(report.is_totally_nonnegative() && report.chain_conclusion_holds());

// Search for a planar strip digraph reproducing the rows (1+t)^n.
let rows: Vec<Polynomial<BigInt>> = (1..=4)
    .map(|n| Polynomial::new((0..=n).map(|k| numeric::binomial_bigint(n, k)).collect()))
    .collect();
let result = recognize_brenti_strip_digraph(&rows, &BrentiRecognitionOptions::default()).unwrap();
assert!(result.is_certified());
```

- `tnn_network::check_lower_unitriangular_totally_nonnegative` decides total nonnegativity
  of a lower-unitriangular matrix exactly; on failure `NetworkError` names the
  bad entry or shape problem, or gives a `NegativeMinor` (rows, columns, exact
  value).  `lower_triangular_rows_from_square` converts a square matrix with
  diagnostics.
- `chain_polynomials_bigint`, `chain_polynomials_rational`, and
  `check_chain_polynomials_bigint` implement the chain polynomials of
  Brändén--Saud Maia Leite (arXiv:2412.06595, Theorem 3.7): for a TNN
  lower-unitriangular `R` they are real-rooted in `[-1, 0]` and interlace.
  The report keeps the TNN verdict separate from the polynomial checks.  The
  theorem is about chain polynomials, not rows: `[[1,0,0],[1,1,0],[1,1,1]]` is
  TNN but has the row `1 + t + t^2`.
- `brenti_recognition::recognize_brenti_strip_digraph` searches exactly for a
  `BrentiStripDigraph` (height-dependent transitions, monotone targets) whose
  rows from `(0,0)` are the given polynomials, with rational or integer
  weights and an optional bound on the target span.  It returns a verified
  certificate, an exact proof of impossibility within the bounds, or
  "inconclusive" when a budget runs out.  A certificate covers only the
  supplied prefix and proves that the coefficient matrix is TNN (by
  Lindström--Gessel--Viennot); it does not prove that the rows are
  real-rooted (`[1 + t, 1 + t + t^2]` is recognized).
- `pf_network::recognize_pf_network` searches exactly for a
  `HeightInvariantNetwork` (steps `(length, rise, weight)` per column, the
  same at every height) whose rows are the given polynomials.  Such a network
  certifies that **each row is a PF sequence** (real-rooted, nonpositive
  zeros): its planarity makes the Toeplitz matrix of every row a path matrix
  with fully compatible sources and sinks (Brenti 1995).
  `build_pf_rows_certificate` and `PfRowsCertificate::verify_against_bigint_polynomials`
  re-check a given network (exact planarity test, rows, explicitly rebuilt
  Toeplitz windows); `toeplitz_path_matrix` and `toeplitz_window` expose the
  two sides of that comparison.
- `left_production::check_left_production_matrix` computes the left
  production matrix `Q_N` of a triangle with nonzero diagonal and decides
  whether it is TNN; if so (Chen--Fu--Ruan, arXiv:2512.08369), the rows are
  PF and the triangle and its reversal are TNN.  `LeftProductionCertificate::verify`
  re-checks a claimed certificate, `unified_network_toeplitz_window` evaluates
  the paper's Theorem 4 window, and the report includes an explicit planar
  network for `Q_N`.  Failure of TNN proves nothing about the rows.

```rust
use num_bigint::BigInt;
use polytool::*;

// PF certificates for the rows (1+t)^n.
let rows: Vec<Polynomial<BigInt>> = (1..=4)
    .map(|n| Polynomial::new((0..=n).map(|k| numeric::binomial_bigint(n, k)).collect()))
    .collect();
assert!(recognize_pf_network(&rows, &PfNetworkOptions::default()).unwrap().is_certified());
assert!(check_left_production_matrix(&rows, 1).unwrap().rows_certified_pf());

// 1 + t + t^2 is never certified.
let bad = vec![Polynomial::from_i64_coeffs(&[1, 1]), Polynomial::from_i64_coeffs(&[1, 1, 1])];
assert!(!recognize_pf_network(&bad, &PfNetworkOptions::default()).unwrap().is_certified());
assert!(!check_left_production_matrix(&bad, 1).unwrap().rows_certified_pf());
```

## Polynomial arithmetic

```rust
use polytool::Polynomial;

let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]);  // 1 + t
let q = p.clone() * p.clone();                         // 1 + 2t + t^2
assert_eq!(q.evaluate(&3), 16);
assert!(q.is_palindromic());
assert_eq!(q.derivative(), Polynomial::from_i64_coeffs(&[2, 2]));
```

For plain coefficient vectors, `polytool::vec_poly` provides the same
operations on `&[C]` for any `CoeffRing` (`i64`, `i128`, `BigInt`,
`BigRational`, ...). Results are trimmed, and zero is returned as `[0]`.
Fixed-width coefficient types panic on overflow instead of wrapping.

```rust
use num_bigint::BigInt;
use polytool::vec_poly;

let p = [1i64, 1];                                   // 1 + t
assert_eq!(vec_poly::pow(&p, 3), vec![1, 3, 3, 1]);
assert_eq!(vec_poly::compose(&[0, 0, 1], &p), vec![1, 2, 1]); // (1 + t)^2
assert_eq!(vec_poly::reverse(&[1, 2], 3), Some(vec![0, 0, 2, 1])); // t^3 p(1/t)
assert_eq!(vec_poly::wronskian(&[0, 1], &[-1, 0, 1]), vec![-1, 0, -1]);

let big: Vec<BigInt> = vec![BigInt::from(2), BigInt::from(0), BigInt::from(1)];
assert_eq!(vec_poly::derivative(&big), vec![BigInt::from(0), BigInt::from(2)]);
```

The module also covers the small helpers that research code otherwise
re-implements: exact division by `t` and by `t - a`
(`exact_div_by_t`, `exact_div_by_t_minus_one`, `exact_div_by_t_minus`; these
return `None` rather than dropping a nonzero remainder), `sum` and
`prefix_sums` of many polynomials, multiplication by a linear factor
`a + b t` and its powers (`mul_linear`, `mul_linear_pow`),
`has_nonnegative_coefficients`, and division-free evaluation at fractions
(`evaluate_homogeneous` gives `den^n p(num/den)`; `evaluate_at_ratio`
evaluates an integer polynomial at a `BigRational` point).

```rust
use polytool::vec_poly;

assert_eq!(vec_poly::exact_div_by_t_minus_one(&[1i64, 0, 0, -1]), Some(vec![-1, -1, -1]));
assert_eq!(vec_poly::exact_div_by_t(&[1i64, 2]), None);
let rows = vec![vec![1i64], vec![0, 1], vec![0, 0, 1]];
let prefix = vec_poly::prefix_sums(&rows); // prefix[k] = rows[0] + ... + rows[k-1]
assert_eq!(vec_poly::sub(&prefix[3], &prefix[1]), vec_poly::sum(&rows[1..]));
assert_eq!(vec_poly::mul_linear_pow(&[1i64], &1, &-1, 2), vec![1, -2, 1]); // (1 - t)^2
assert!(vec_poly::has_nonnegative_coefficients(&[1i64, 0, 3]));
```

Exact binomials and factorials live in `polytool::numeric`:
`binomial_bigint`, `factorial_bigint`, `falling_factorial_bigint`,
`multinomial_bigint`, and the checked `binomial_u128`.  Its q-analogues return
ascending `Vec<BigInt>` polynomials in `q`: `q_integer(n)`, `q_factorial(n)`,
`q_binomial(n, k)` (Gaussian binomial), `q_catalan(n)` (MacMahon's maj/inv
`[2n choose n]_q / [n+1]_q`) and `q_narayana(n, k)`
(`q^(k(k-1)) [n,k]_q [n,k-1]_q / [n]_q`).

## Standard sequences and indexing

`polytool::sequences` generators take `max_n` and return `rows[n]` for
`n = 0..=max_n` (the n-th member of the family, with `rows[0]` the `n = 0`
member), so Eulerian rows are `A_0 = 1, A_1 = 1, A_2 = 1 + t, ...`.  The
OEIS triangles A008292 and A001263 start at `n = 1`, i.e. `rows[1..]`; the
same holds for `literature::hoggatt_polynomials_bigint`, whose `rows[0]` is the
zero polynomial.  The one exception follows its source's normalization and
says so in its rustdoc: the `deco` histogram layers take `max_size` and return
`rows[N - 1]` for size `N`.  Besides the classical
families the module provides permutation-statistic generators over `BigInt`:
`derangement_descent_polynomials_bigint` (derangements of `[n]` by descents
in one-line notation, OEIS A219836; this is not the excedance triangle
A046739), `fixed_point_descent_polynomials_bigint` (the bivariate
fixed-point/descent array `rows[n][j][k]`),
`binomial_eulerian_polynomials_bigint`, `refined_eulerian_polynomials_bigint`
(descents by first letter) and `s_eulerian_polynomial_bigint`.
`shifted_legendre_polynomials_bigint` gives the shifted Legendre polynomials
`P~_n(x) = P_n(2x - 1) = sum_k (-1)^(n+k) C(n,k) C(n+k,k) x^k` on `[0, 1]`
(`P~_n(1) = 1`, coefficients up to sign in OEIS A063007).

## Characteristic polynomials and closed walks

```rust
use num_bigint::BigInt;
use polytool::linalg::{characteristic_polynomial_bigint, closed_walk_counts_bigint};

// Path on three vertices: det(tI - A) = t^3 - 2t.
let a: Vec<Vec<BigInt>> = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
    .iter()
    .map(|r| r.iter().map(|&x| BigInt::from(x)).collect())
    .collect();
let chi = characteristic_polynomial_bigint(&a); // ascending: [0, -2, 0, 1]
let walks = closed_walk_counts_bigint(&a, 4); // tr(A^k), k = 0..=4: [3, 0, 4, 0, 8]
# let _ = (chi, walks);
```

`characteristic_polynomial_bigint` uses the Faddeev--LeVerrier recurrence over
the integers (exact divisions only).  `is_ultra_log_concave_infinite` tests
ULC(∞) (`k a_k^2 >= (k+1) a_(k-1) a_(k+1)`), the weakest member of the
degree-based ULC(n) family tested by `is_ultra_log_concave_of_degree`.

## Recurrence search

```rust
use polytool::recurrence::*;

// Eulerian polynomials A_1, A_2, ..., A_10
let polys: Vec<Vec<i64>> = vec![
    vec![1],
    vec![1, 1],
    vec![1, 4, 1],
    vec![1, 11, 11, 1],
    vec![1, 26, 66, 26, 1],
    vec![1, 57, 302, 302, 57, 1],
    vec![1, 120, 1191, 2416, 1191, 120, 1],
    vec![1, 247, 4293, 15619, 15619, 4293, 247, 1],
    // ...
];

// Adaptive search: tries small parameter spaces first
let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
if let Some(res) = result {
    println!("{}", res.recurrence);
}

// Deterministic candidate budget with an exact termination reason.
let outcome = find_recurrence_adaptive_with_budget(
    &polys,
    &AdaptiveSearchOptions::default(),
    AdaptiveSearchBudget::limited(100),
);
match outcome {
    AdaptiveSearchOutcome::Found(result) => println!("{}", result.recurrence),
    AdaptiveSearchOutcome::NoRecurrence(_) => println!("search space exhausted"),
    AdaptiveSearchOutcome::BudgetExhausted(summary) => println!(
        "budget exhausted after {} candidates",
        summary.diagnostics.considered_candidates
    ),
}

// Or search with specific parameters
let opts = RecurrenceOptions {
    rec_len: 2,
    var_deg: 1,
    idx_deg: 1,
    diff_deg: 1,
    ..Default::default()
}.with_alternating_sign(true);
if let Some(rec) = find_polynomial_recurrence(&polys, &opts) {
    println!("{}", rec);
}
```

## Shared-denominator vector recurrences

`find_vector_recurrence_with_denominator_rational` (and its integer wrapper)
searches exact polynomial-vector transitions of the form

```text
q(n,x) F_(n+1) = M(n,x,D_x) F_n + G(n,x).
```

Use `VectorRecurrenceDenominatorOptions` to set the matrix, forcing, and
common-factor degree bounds.  The returned
`SharedDenominatorVectorRecurrence` evaluates by exact polynomial division and
returns `SharedDenominatorVectorRecurrenceEvaluationError` for dimension,
zero-denominator, or non-polynomial-quotient failures.  The
fitting functions return `SharedDenominatorVectorRecurrenceFitError` and
joint diagnostics rather than the per-row diagnostics used by the normalized
vector API.  The
`find_companion_vector_recurrence_with_denominator[_rational]` variants support
lagged systems; companion shift rows contain `q(n,x) I` in the numerator.

The denominator is normalized canonically: its first nonzero coefficient in
`n^0 x^0, n^0 x^1, ..., n^1 x^0, ...` order is one.  With
`require_unique = true`, all possible normalized pivots are checked and an
ambiguous projective fit is rejected.  Existing `q = 1` vector APIs and their
options are unchanged.

For Hermite--Biehler or Euclidean-chain experiments, the library also exposes a
signed Sturm continued-fraction certificate checker over exact rationals:

```rust
use num_bigint::BigInt;
use num_rational::BigRational;
use polytool::check_sturm_continued_fraction_certificate;

let q = |n: i64| BigRational::from_integer(BigInt::from(n));
let p0 = vec![q(1), q(3), q(1)]; // x^2 + 3x + 1
let p1 = vec![q(2), q(1)];       // x + 2
let certificate =
    check_sturm_continued_fraction_certificate(&p0, &p1).expect("Sturm-CF certificate");
assert_eq!(certificate.steps, 2);
```

To build the MCP server from the repository root:

```sh
cargo build --release --manifest-path mcp/Cargo.toml --target-dir target
```

The MCP binary is at `target/release/polytool-mcp`.

## Coupled Weyl-matrix recurrences

The library also fits exact first-order systems

```text
F_(n+1) = M(n,x,D_x) F_n + G(n,x),
```

where `F_n` is a vector of polynomials. Every matrix entry is stored in the
normal form `sum_d c_d(n,x) D_x^d`; this is general for the Weyl algebra because
`D_x x = x D_x + 1`. Input has shape
`states[index][component][x_degree]`, and `first_index` explicitly fixes the
meaning of `n`:

```rust
use polytool::recurrence::{find_vector_recurrence, VectorRecurrenceOptions};

// F_n = [n]_x, including [0]_x = 0.  The affine recurrence is
// F_(n+1) = x F_n + 1.
let states = (0..=9)
    .map(|n| vec![if n == 0 { vec![0] } else { vec![1; n] }])
    .collect::<Vec<_>>();
let options = VectorRecurrenceOptions {
    var_deg: 1,
    idx_deg: 0,
    homogeneous: false,
    forcing_var_deg: 0,
    forcing_idx_deg: 0,
    held_out_transitions: 2,
    ..Default::default()
};
let fit = find_vector_recurrence(&states, 0, &options)?;
assert_eq!(fit.rows[0].rank, fit.rows[0].unknowns);
# Ok::<(), Box<dyn std::error::Error>>(())
```

Rows of `M` are solved independently over exact rationals. The result reports
the exact rank and nullity of each row; `require_unique` defaults to `true`, so
cross-component dependencies are rejected instead of being hidden by an
arbitrary choice of free parameters. Final complete transitions are held out
of the fit and then verified exactly.

For a lag-`r` system, `find_companion_vector_recurrence` returns the first-order
recurrence on `(F_(n-r+1), ..., F_n)`. It fits only the unknown final block row
and inserts the shift-identity rows directly. Unit tests recover two published
systems: the up-down-run Eulerian pair of
[Ma--Ma--Yeh--Yeh (2022)](https://doi.org/10.1016/j.disc.2021.112716), and the
lag-two type-B `1/k`-Eulerian system of
[Ma et al. (2020)](https://doi.org/10.37236/9089), specialized to `k=1`.

## Other properties

```rust
use polytool::*;

let coeffs = [1, 11, 11, 1];

assert!(is_palindromic(&coeffs));
assert!(is_unimodal(&coeffs));
assert!(is_log_concave(&coeffs));
assert!(is_ultra_log_concave(&coeffs));
// ULC with respect to a nominal degree n ≥ deg p (e.g. a matroid rank)
assert!(is_ultra_log_concave_of_degree(&coeffs, 5));
assert!(is_gamma_positive(&coeffs));

// Gamma coefficients: p(t) = Σ γ_i t^i (1+t)^{d-2i}
let gamma = gamma_coefficients(&coeffs).unwrap();
assert_eq!(gamma, vec![1, 8]);

// Resultant and discriminant
let disc = discriminant(&[1, 0, -3, 1]);
```

Hurwitz stability (all roots in the open left half-plane) is decided exactly
by `polytool::hurwitz`: `hurwitz_matrix_bigint` (classical layout, built from
the ascending coefficients as given), `hurwitz_leading_minors_bigint` (the
Hurwitz determinants `Δ_1, ..., Δ_n`) and `is_hurwitz_stable_bigint` (which
normalizes a negative leading coefficient first).  A zero constant term is a
root at the origin and makes the polynomial unstable; nonzero constants are
stable and the zero polynomial is not.  `bezoutian_bigint_coeffs(f, g)` gives
the symmetric `max(deg f, deg g)`-square Bézoutian for arbitrary degrees, with
the orientation of `bezout_matrix_bigint_coeffs` (entry `[i][j]` is the
coefficient of `x^i y^j` in `(f(x) g(y) - f(y) g(x)) / (x - y)`) and rank
`max(deg f, deg g) - deg gcd(f, g)`.

```rust
use num_bigint::BigInt;
use polytool::{bezoutian_bigint_coeffs, is_hurwitz_stable_bigint};

let big = |v: &[i64]| v.iter().map(|&c| BigInt::from(c)).collect::<Vec<_>>();
assert!(is_hurwitz_stable_bigint(&big(&[1, 3, 3, 1])));   // (1 + t)^3
assert!(!is_hurwitz_stable_bigint(&big(&[2, 1, 1, 1])));  // Δ_2 < 0
let b = bezoutian_bigint_coeffs(&big(&[-1, 0, 1]), &big(&[0, 1]));
assert_eq!(b, vec![big(&[1, 0]), big(&[0, 1])]);
```

To treat a zero polynomial as vacuously interlacing in either position, write
`vec_poly::is_zero(p) || vec_poly::is_zero(q) || check_weak_interlacing(p, q)
== Some(true)`; `interlaces_either_order` already lets zero interlace every
real-rooted polynomial.

## Enumerating real-rooted integer polynomials

`real_rooted_enumeration::enumerate_real_rooted_integer_polynomials_bigint`
returns all real-rooted `p ∈ Z[x]` of a given degree with at least three
prescribed leading coefficients (given in descending order; results are
ascending vectors in lexicographic order of the free coefficients).  It walks
the scaled derivative tower `p_i = p^(n-i)`, where
`p_i = ∫_0^x p_(i-1) + (n-i)! a_i`, and at each level keeps exactly the
constants allowed by Lemma 2.5 of Greaves--Syatriadi: for real-rooted `q` and
`P = ∫_0^x q`, `P + C` is real-rooted iff `C >= -P(µ)` at the roots with
`q'(µ) < 0`, `C <= -P(µ)` at those with `q'(µ) > 0`, and `C = -P(µ)` at
multiple roots.  The paper's floating-point endpoint step is replaced by exact
root isolation, dyadic enclosures of the critical values refined until they
are narrower than the factorial modulus (or contain none of its multiples),
and an exact real-rootedness test of the at most one undecided candidate at
each end.  Every emitted polynomial is re-checked exactly.

```rust
use num_bigint::BigInt;
use polytool::real_rooted_enumeration::{
    enumerate_real_rooted_integer_polynomials_bigint,
    real_rooted_antiderivative_integer_constants_bigint, RealRootedEnumerationOptions,
    RealRootedEnumerationStatus,
};

let prefix = [1, 0, -6].map(BigInt::from); // x^4 - 6x^2 + a_3 x + a_4
let options = RealRootedEnumerationOptions {
    max_results: Some(10_000),
    ..Default::default()
};
let result = enumerate_real_rooted_integer_polynomials_bigint(4, &prefix, &options).unwrap();
assert_eq!(result.status, RealRootedEnumerationStatus::Complete);
assert_eq!(result.polynomials.len(), 64);

// The kernel on its own: x^3 - 3x + C is real-rooted iff -2 <= C <= 2.
let q = [-3, 0, 3].map(BigInt::from);
let range = real_rooted_antiderivative_integer_constants_bigint(&q).unwrap().unwrap();
assert_eq!((range.lower, range.upper), (BigInt::from(-2), BigInt::from(2)));
```

`RealRootedEnumerationOptions` bounds examined nodes, returned results and the
bit length of tower coefficients; exhausting one gives
`RealRootedEnumerationStatus::ResourceLimit(limit)` with the results found so
far (an initial segment of the complete output), never an empty "complete"
answer.  `real_rooted_antiderivative_constant_multiples_bigint` restricts the
kernel to multiples of a modulus, and `antiderivative_exact_bigint` /
`antiderivative_rational_bigint` give zero-constant antiderivatives (the
integral version reports any non-exact division).  The search statistics are
deterministic counters; no wall-clock times are recorded.  The MCP tool
`enumerate_real_rooted` exposes the operation with conservative ceilings (see
[the MCP notes](mcp.md)).
