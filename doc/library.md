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

## Ore algebra and image recurrences

The exact Ore algebra of shifts `S` and derivatives `D` acting on sequences of
polynomials (`S n = (n + 1) S`, `D x = x D + 1`) is the independent
`ore-algebra` crate under `crates/ore-algebra` (dependencies `num-*` only).
Besides the arithmetic it computes runnable annihilators of image sequences
with identity certificates.  Polytool wraps it in
`recurrence::image_recurrence_of_operator`, which turns a runnable rule of
`P` into a runnable rule, with seeds and a replayable certificate, for
`Q_n = lambda(n) T_n(P_(n+s))`:

```rust
use polytool::recurrence::{
    image_recurrence_of_operator, BivarPolyJson, ImageOperatorSpecification,
    ImageRecurrenceOptions, ImageRowMap, OperatorTermJson, PreserverParentLibrary,
};

let (_, eulerian) = PreserverParentLibrary::named_family_records()
    .into_iter()
    .find(|(name, _)| name == "eulerian")
    .unwrap();
let term = |deriv_order: usize| OperatorTermJson {
    deriv_order,
    coeff: BivarPolyJson { coeffs: vec![vec!["1".into()]] },
};
// T = 1 + D: A_n + A_n'.
let specification = ImageOperatorSpecification::new(ImageRowMap::WeylOperator {
    terms: vec![term(0), term(1)],
    denominator: vec!["1".into()],
});
let image = image_recurrence_of_operator(
    &eulerian,
    eulerian.first_index,
    &specification,
    &ImageRecurrenceOptions::default(),
)
.unwrap();
assert_eq!(image.order, 2);
image.verify().unwrap();
```

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
`literature::one_over_k_eulerian_polynomials_rational(max_n, &k)` gives the
`1/k`-Eulerian polynomials (Savage--Viswanathan) for an exact rational `k`
(`A_0 = 1`, `A_(n+1) = (1 + k n x) A_n + k x (1 - x) A_n'`, `k = 1` the
Eulerian rows), and `one_over_k_eulerian_symmetric_parts_rational` their
symmetric decompositions `(a_n, b_n)` with respect to `d_n = n - 1`.
`literature::peak_polynomials_bigint` gives the peak polynomials
`W_n(z) = sum z^(pk(sigma))` over the permutations of `[n]` (OEIS A008303 from
`n = 1`; by complementation also the valley polynomials), with
`peak_polynomials_iter`.

Families with a recurrence from the literature expose it next to the
generator as `(Recurrence, seeds, first_index)` (`literature::RecurrenceWithSeeds`):
`peak_polynomials_recurrence`
(`W_n = (2 + (n - 2) z) W_(n-1) + 2z (1 - z) W'_(n-1)` from `W_1 = 1`, the
generator's own rule), and `one_over_k_eulerian_symmetric_part_a_recurrence(&k)`
and `one_over_k_eulerian_symmetric_part_b_recurrence(&k)` (second-order rules
with second derivatives from `n = 1`, eliminated from the coupled step of the
decomposition; the derivation is in the rustdoc).  The recurrences feed the
real-rootedness routes directly instead of a fitted rule: the peak polynomials
are proved by a gamma lift from the Eulerian polynomials, and the symmetric
parts `a_n` (`k = 2, 3`, from `n = 2`) by the symmetric-companion lift.

## Lazy rows

Recurrence-based generation can be consumed one polynomial at a time, without
choosing a final index and without recomputing earlier rows.
`Recurrence::generate_rows_rational_iter(&initial_rows, first_index)` returns a
`recurrence::RecurrenceRowIterator`: the seed rows come first (as
`P_{first_index}, ...`), then one exactly evaluated row per `next()`, each
item a `Result<Vec<BigRational>, RecurrenceEvaluationError>`.  Nothing is
evaluated until a row is requested, only the last `max_offset()` rows are kept
(the polynomials themselves may grow), and a failing step yields its error once
and ends the iterator (it is fused), so an earlier valid prefix stays
consumable.  `generate_rows_rational` collects a prefix of this iterator.

The same holds for the OEIS catalog:
`OeisSequenceDefinition::generate_rows_iter()` and `generate_rows_from_iter(first_row)` yield
`Result<Vec<BigInt>, OeisSequenceError>` over the authoritative prefix, the
recurrence's initial rows and the computed rows, with displayed row indices,
checked integer conversion and triangle padding exactly as in `generate_rows`
and `generate_rows_from`, which collect from them.  A new iterator still has to
traverse the rows before its start (there is no random access); continuing an
iterator you keep does not.  Named recurrence families have endless lazy
counterparts that yield member `n` as the `n`-th item:
`sequences::eulerian_polynomials_iter`, `type_b_eulerian_polynomials_iter`,
`chebyshev_polynomials_t_iter`, `chebyshev_polynomials_u_iter`,
`hermite_polynomials_iter`, and in `sequences::literature`
`delannoy_polynomials_iter` and `delannoy_square_polynomials_iter`.  Closed-form
families (Narayana, shifted Legendre, Hoggatt) and the matrix- or
array-valued recurrences keep their batch generators.

```rust
use polytool::oeis::by_id;

let sequence = by_id("A008292").unwrap();
let mut rows = sequence.generate_rows_iter().unwrap();
let first_three: Vec<_> = rows.by_ref().take(3).collect::<Result<_, _>>().unwrap();
let fourth = rows.next().unwrap().unwrap(); // continues; rows 1-3 are not recomputed
assert_eq!(fourth, [1, 11, 11, 1].map(num_bigint::BigInt::from).to_vec());
assert_eq!(first_three.len(), 3);
```

See `examples/lazy_rows.rs`.

### OEIS real-rootedness flags

A catalog record may carry a real-rootedness flag (`OeisCatalogRealRootedness`
in the data file, documented with the CLI in `doc/cli-oeis.md`).
`OeisSequenceDefinition::real_rootedness()` returns it as an
`OeisRealRootedness { status, as_of }` (`None`: unknown), and
`oeis::entries_with_real_rootedness_flag()` lists every flagged entry.
`OeisRealRootednessStatus::LeanProved` means that every row is zero or
real-rooted by a machine-checked Lean proof made outside polytool, for a
definition whose rows were shown equal to the catalog rows for every index;
`KnownNotRealRooted` means that an exact check found a row with a non-real
root.  The flag is never reported as a polytool proof.

```rust
use polytool::oeis::{self, OeisRealRootednessStatus};

let flag = oeis::by_id("A001263").unwrap().real_rootedness().unwrap();
assert_eq!(flag.status, OeisRealRootednessStatus::LeanProved);
assert!(oeis::by_id("A157153").unwrap().real_rootedness().is_none());
```

The Veronese parent provers (next sections) use a `LeanProved` entry as an
attestation: a parent whose rule, first index and initial rows equal the
catalog entry's `recurrence_parts()` at `recurrence_first_index` exactly gets
the verdict `ConditionalOnAttestation` with the label
`oeis_catalog:<A-number>:lean_proved`, provided polytool also proves that every
row has coefficients of one sign (the flag gives real-rootedness only); a
parent that differs in one coefficient or seed does not match.

### OEIS evidence records

A catalog record may also carry finite evidence about its rows
(`OeisCatalogEvidence` in the data file; the definitions and the CLI are in
`doc/cli-oeis.md`).  `OeisSequenceDefinition::evidence()` returns it as an
`OeisEvidence { strongest, rows_checked, first_failure, as_of }`: the strongest
`OeisEvidenceProperty` (real-rooted, ultra-log-concave, log-concave, unimodal,
or `NoProperty`) that holds on every one of the first `rows_checked` rows, and
the first failing displayed row of each stronger property.  This is evidence
about a finite segment, not a proof, and it is independent of the
real-rootedness flag.  `compute_evidence(rows, as_of)` recomputes the record
exactly from the entry's own rows; `OeisEvidenceAccumulator` builds it row by
row, and `row_has_evidence_property` evaluates one row.

```rust
use polytool::oeis::{self, OeisEvidenceProperty};

// The Eulerian triangle: every checked row is real-rooted.
let eulerian = oeis::by_id("A008292").unwrap();
let evidence = eulerian.compute_evidence(30, "2026-11-02").unwrap().unwrap();
assert_eq!(evidence.strongest, OeisEvidenceProperty::RealRooted);

// Row 2 of A122431 is 3 + 3t + t^2: not real-rooted, still log-concave.
let other = oeis::by_id("A122431").unwrap();
let evidence = other.compute_evidence(30, "2026-11-02").unwrap().unwrap();
assert_eq!(evidence.strongest, OeisEvidenceProperty::LogConcave);
assert_eq!(evidence.first_failure.real_rooted, Some(2));
// `other.evidence()` returns the record stored in the data file, if any.
```

### B-files from any rows

`oeis::flatten_rows_to_bfile` turns coefficient rows (from a catalog entry or any
recurrence) into OEIS b-file text.  `BFileOptions` sets the first index (the
OEIS offset), the row orientation, a policy for zeros at the ends of each row
and a term cap that never splits a row; `BFileStoppingRule` chooses between the
longest prefix with at most the cap many terms (the default) and the whole row
nearest to the cap.

```rust
use polytool::oeis::{self, BFileOptions};

let rows = oeis::by_id("A008292").unwrap().generate_rows(4).unwrap();
let options = BFileOptions { first_index: 1, max_terms: Some(8), ..Default::default() };
let file = oeis::flatten_rows_to_bfile(&rows, &options).unwrap();
assert_eq!(file.text, "1 1\n2 1\n3 1\n4 1\n5 4\n6 1\n");
assert_eq!(file.rows_included, 3);
```

### OEIS text from any recurrence

`oeis::OeisRecurrencePresentation` re-indexes a recurrence to the rows of an
OEIS triangle (`OeisTriangleReading`: the OEIS row number of the first row,
ascending or descending powers, zero trimming, regular-triangle padding) and
formats the `%e` block, the `%F` sentence and Mathematica and PARI/GP programs
that list the same rows; the recurrence index `n` in the text is the OEIS row
number.

```rust
use polytool::oeis::{BFileRowOrientation, OeisRecurrencePresentation, OeisTriangleReading};
use polytool::recurrence::RecurrenceJson;

let record: RecurrenceJson = serde_json::from_str(r#"{
  "schema": "polytool.recurrence.v1", "first_index": 0, "initial_polynomials": [["1"]],
  "recurrence": {"terms": [
    {"offset": 1, "deriv_order": 0, "sign": "none", "coeff": {"coeffs": [["1", "-2"], ["0", "2"]]}},
    {"offset": 1, "deriv_order": 1, "sign": "none", "coeff": {"coeffs": [["0", "2", "-2"]]}}]}
}"#).unwrap();
let reading = OeisTriangleReading {
    first_row_index: 0,
    orientation: BFileRowOrientation::Reversed,
    regular_triangle: true,
    ..OeisTriangleReading::default()
};
// A156919 is this triangle reversed from the row n = 1 on.
let presentation =
    OeisRecurrencePresentation::from_recurrence_json(&record, Some(1), reading).unwrap();
let text = presentation.oeis_text(3).unwrap();
assert_eq!(text.triangle_example, "Triangle begins:\n  0: 1;\n  1: 2,  1;\n  2: 4, 10, 1;\n");
assert!(text.recurrence_sentence.contains("P(n,x) = (1+2*n*x)*P(n-1,x)"));
```

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

**Left-hand side variation.**  With `AdaptiveSearchOptions::try_laurent`
(builder `with_laurent`; off by default in the library so that internal
callers such as the lift searches keep their outcomes, on by default in
`polytool recurrence` and the MCP recurrence tools) the adaptive search adds a
third stage after its two
historical passes: the same shapes with a left-hand monomial `x^e`
(`e <= max_laurent_x_power`, default 2; `RecurrenceOptions::leader_x_power`,
so the rule reads `x^e f(n, t) P(n) = ...`) and with the left-hand factor
normalized on a coefficient other than its constant term (`x^j` or `n^i`
fixed to 1).  This reaches factors whose constant term is zero, for example
`n P(n) = (2n - 1)(2t - 1) P(n-1) - (n - 1) P(n-2)` for the shifted Legendre
polynomials.  Results the historical passes find are unchanged; the stage
runs only when they find nothing.  The rule keeps the whole left-hand factor
as its `denominator`, rows are exact quotients (a non-divisible row is
`RecurrenceEvaluationError::NonPolynomialQuotient`), and the JSON record of
the rule carries `x^e` in its denominator.

## Veronese-section lifts

`recurrence::find_veronese_section_lifts` looks for parents `P_n(x)` whose
coefficient section `V_(r,j)(P_n)(y) = sum_k [x^(r k + j)] P_n(x) y^k` is a
given recurrence-defined sequence `Q_n(y)`.  A real-rooted parent with
coefficients of one sign transfers real-rootedness to its sections, so a lift
with status `ProvedProductParent` proves that every `Q_n` is zero or
real-rooted with nonpositive roots.  The input recurrence must have
derivative-free, non-alternating terms with coefficients independent of `n`,
a constant (or no) denominator and no inhomogeneous term.

```rust
use polytool::recurrence::*;

let q = |n: i64| BigRational::from_integer(n.into());
// Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n, Q_0 = Q_1 = 1.
let recurrence = Recurrence {
    terms: vec![
        RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(2)]] }),
        RecurrenceTerm::new(2, 0, BivarPoly { coeffs: vec![vec![q(-1), q(1)]] }),
    ],
    denominator: None,
    inhomogeneous: None,
};
let rows = vec![vec![q(1)], vec![q(1)]];
let report = find_veronese_section_lifts(
    &recurrence,
    &rows,
    0,
    &VeroneseSectionLiftOptions::default(), // max_modulus 4, 12 parent rows, 256 candidates
)?;
let best = &report.lifts[0]; // proved lifts first
assert_eq!((best.modulus, best.residue), (2, 0));
assert_eq!(
    format_characteristic_polynomial(&best.parent_characteristic_polynomial, "x"),
    "L - 1 - x"
); // P_n = (1 + x) P_(n-1), P_0 = 1
assert_eq!(best.status, VeroneseSectionLiftStatus::ProvedProductParent);
// The parent as a polytool.recurrence.v1 record for the other recurrence tools.
let parent_json = best.parent_recurrence_json();
# let _ = parent_json;
# Ok::<(), RecurrenceError>(())
```

The report lists the input characteristic polynomial `chi(L, y)`, the lifts
(modulus, residue, parent characteristic polynomial `G(L, x)`, section
characteristic polynomial `N(L, y)`, cofactor `chi / N`, parent recurrence and
seeds, exact parent-row checks and status), the rejected candidates with a
reason, moduli whose factorization failed, and whether the candidate budget
was exhausted.  The report types implement `Serialize`
(`polytool.recurrence-veronese-lifts.v1`, exact rationals as strings).  The
algorithm and the argument that `d` identifying rows prove the identity for
all `n` are in [algorithms.md](algorithms.md#veronese-section-lifts).

### Certifying one candidate (index-dependent and differential inputs)

`recurrence::certify_veronese_section_lift` checks one candidate
`Q_n = lambda(n) V_(r,j)(P_(n+s))` for an input whose rule may have
index-dependent coefficients, derivatives on earlier rows and a denominator
`f(n, y)`.  The candidate (`VeroneseSectionCandidate`) gives the modulus
(2 to 4), the residue, the parent as a `polytool.recurrence.v1` record whose
denominator is free of `x`, an index shift `|s| <= 4` and a gauge
(`VeroneseSectionGauge::{None, Ratio, Closed}`).  The result
(`VeroneseSectionCertification`, `polytool.recurrence-veronese-certificate.v1`)
records the identity status (`Certified` with a replayable
`VeroneseForwardDivisionCertificate`, `SupportedByRows`, `Refuted`), the
parent proof (periodic product, paper theorem, root-domain claim, the
three-term, Liu--Wang, Liu--Wang derivative-lag or Ma--Wang interval theorem
below, a caller-supplied
`VeroneseParentAttestation`, or a `lean_proved` OEIS catalog flag, both of which
only give a conditional verdict) and the verdict.  `verify()` replays the certificate.

```rust
use polytool::recurrence::*;

let q = |n: i64| BigRational::from_integer(n.into());
// m (m - 1) Q_m = 2 (m - 1) Q_(m-1) + (y - 1) Q_(m-2), Q_0 = Q_1 = 1:
// the even sections of (1 + x)^m divided by m!.
let input = Recurrence {
    terms: vec![
        RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(-2)], vec![q(2)]] }),
        RecurrenceTerm::new(2, 0, BivarPoly { coeffs: vec![vec![q(-1), q(1)]] }),
    ],
    denominator: Some(BivarPoly { coeffs: vec![vec![q(0)], vec![q(-1)], vec![q(1)]] }),
    inhomogeneous: None,
};
let parent = Recurrence {
    terms: vec![RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(1), q(1)]] })],
    denominator: None,
    inhomogeneous: None,
};
let candidate = VeroneseSectionCandidate {
    modulus: 2,
    residue: 0,
    parent: RecurrenceJson::from_recurrence_rational(&parent, 0, &[vec![q(1)]], None),
    index_shift: 0,
    // lambda(n) / lambda(n - 1) = 1 / n, lambda(0) = 1.
    gauge: VeroneseSectionGauge::Ratio {
        numerator: vec!["1".into()],
        denominator: vec!["0".into(), "1".into()],
        initial_value: "1".into(),
    },
};
let result = certify_veronese_section_lift(
    &input,
    &[vec![q(1)], vec![q(1)]],
    0,
    &candidate,
    &VeroneseSectionCertifyOptions::default(),
)?;
assert_eq!(result.verdict, VeroneseSectionVerdict::Proved);
result.verify()?;
# Ok::<(), RecurrenceError>(())
```

The transport conventions, the pseudo-division and the soundness argument are
in [algorithms.md](algorithms.md#certifying-one-candidate-v2).

### Real-rootedness theorems for recurrences

`recurrence::check_three_term_recurrence_theorem`,
`recurrence::check_liu_wang_theorem`,
`recurrence::check_liu_wang_derivative_lag_theorem` and
`recurrence::check_ma_wang_interval_theorem` apply a paper theorem to a
family given by a rule, a first index and seed rows (the rule applies from
`first_index + seeds.len()` on), and check every hypothesis exactly for all
indices.  The result (`RealRootednessTheoremReport`,
`polytool.recurrence-real-rootedness-theorem.v1`) gives the theorem and its
reference, the matched normal form with its coefficients as rational
functions of `n`, the hypotheses as `LiftHypothesis` entries (verified, failed
with the smallest failing index, or unresolved, with the sign checker's
replayable decision where one applies), the status (`proved`,
`hypothesis_failed`, `unresolved`, `not_applicable`) and, when proved, the
conclusion (first index, interlacing of consecutive rows, simple roots, root
region).  `verify()` replays it.

- Three-term: `f(n) P_n = c_1(n, x) P_(n-1) + c_2(n, x) P_(n-2)` with
  `deg_x c_1 <= 1` and `c_2 = gamma(n) w(x)`, `w` constant (Hermite, Laguerre,
  Chebyshev, Legendre) or `kappa (x + v_0)^2` (then `Q_n(y) =
  y^(d_n) P_n(1/y - v_0)` is the three-term family; Narayana).  With
  `RealRootednessTheoremOptions::root_region` it also certifies that every root
  lies in an open interval.
- Liu--Wang: `f(n) P_n = c_0(n, x) P_(n-1) + c_1(n, x) P'_(n-1)` with
  `deg_x c_0 <= 1`, `deg_x c_1 <= 2`; it proves every root `< 0` (Eulerian,
  type-B Eulerian).
- Liu--Wang derivative lag: `f(n) P_n = c_0 P_(n-1) + c_1 P'_(n-1) +
  c_2 P_(n-2)` with `deg_x c_0 <= 1`, `deg_x c_1, deg_x c_2 <= 2`,
  `V_n = c_1 / f` and `W_n = c_2 / f` nonpositive left of `0` (one of them
  negative there), vanishing at `0`, and growth one of the degrees; it proves
  strict interlacing of consecutive rows and every root `< 0` (Delannoy,
  deco, OEIS A046739).  When a power `x^e` divides the two seed rows of the
  induction it proves the rows divided by `x^e` and concludes every root in
  `(-inf, 0]` (derangement excedances, OEIS A271697, where `P_2 = 0`); zero
  seed rows are allowed and named in the conclusion.  The report records the
  mirrored row route of the Lean library
  [RealRooted](https://github.com/PerAlexandersson/RealRooted) in
  `lean_route` (`halfGrowth`, `sign` for three-term rules, or `LiuWang`) and
  its declarations in `lean_declarations` (both omitted from the JSON for the
  other checks).  The binomial Eulerian rule (OEIS A046802) is refused: its
  lag coefficient `(1 - n) x` is positive left of `0` (witness `n = 2`).
- Ma--Wang interval: `f(n) P_n = c_0(n, x) P_(n-1) + gamma(n) (x - L)(x - U)
  P'_(n-1)` with `deg_x c_0 <= 1` and rational `L < U` read off the rule; it
  proves every root in the closed interval `[L, U]` (an end open when no row
  has a root there) and weak interlacing of consecutive rows, by the exact
  two-endpoint criterion of Ma and Wang (EJC 15 (2008) R17, Theorem 2 and
  Corollary 3; RealRooted `affineEndpointStep_closed`): with
  `a_n = [x^1] c_0 / gamma > 0`, the thresholds
  `c_0(n, U) / ((U - L) gamma) + m_U(n) >= 0` and
  `-c_0(n, L) / ((U - L) gamma) + m_L(n) >= 0`, where the endpoint
  multiplicities `m` of the previous row are tracked exactly (each grows by
  one after an index where its threshold vanishes).  Ordered set partitions
  `k! S(n, k)` (OEIS A019538, A131689), the Ward numbers (A134991, A269939)
  and the Stirling--Frobenius rows (A225472, A225473, A225476) are proved
  this way; a failing threshold names the row that leaves the interval.
  `lean_declarations` lists the RealRooted step theorems (no row route).

When a theorem fails at the supplied seeds, the checks retry with the
induction start moved up to `RealRootednessTheoremOptions::maximum_start_delay`
rows later (default 3): the skipped rows are generated by the rule and enter
as seed rows, checked exactly, so the conclusion still covers every
`n >= first_index`.  The report records the delay in `start_delay` (omitted
from the JSON when zero) and `verify()` replays it.  The Eulerian rule from
`A_0 = 1` (where `A_1 = A_0` does not raise the degree) is proved this way
with `start_delay = 1`.

A rule with a left-hand factor `x^e f(n)` (`e >= 1`, `f` free of `x`), as the
adaptive search returns it with `try_laurent`, is checked through the exactly
equivalent rule `f(n) P_n = sum (c / x^e) D^d P_(n-r)` when `x^e` divides every
right-hand coefficient (terms with the same lag, derivative order and sign are
added first): both rules give the same rows, because `x^e` cancels in `Q[x]` at
every index with `f(n) != 0` and both keep a row with `f(n) = 0` as a seed.  The
report keeps the supplied rule in `recurrence` and records the original and the
normalized rule in `left_hand_factor_normalization`
(`LeftHandFactorNormalization`); the normal form, hypotheses and conclusion are
those of the normalized rule, and `verify()` re-derives the normalization.
When `x^e` does not divide the right-hand side, the division by `x^e` is exact
on the rows only and not an identity of the rule, and the check reports
`not_applicable` with that reason.  The parent provers of the lifts, the lift
rule preparation and the Rolle root-domain matcher normalize the same way.

```rust
use polytool::recurrence::*;
use polytool::real_rootedness::RealInterval;

let q = |n: i64| BigRational::from_integer(n.into());
// Legendre: n P_n = (2n - 1) x P_(n-1) - (n - 1) P_(n-2), P_0 = 1, P_1 = x.
let legendre = Recurrence {
    terms: vec![
        RecurrenceTerm::new(1, 0, BivarPoly { coeffs: vec![vec![q(0), q(-1)], vec![q(0), q(2)]] }),
        RecurrenceTerm::new(2, 0, BivarPoly { coeffs: vec![vec![q(1)], vec![q(-1)]] }),
    ],
    denominator: Some(BivarPoly { coeffs: vec![vec![q(0)], vec![q(1)]] }),
    inhomogeneous: None,
};
let options = RealRootednessTheoremOptions {
    root_region: Some(RealInterval::closed(q(-1), q(1))),
    ..Default::default()
};
let report =
    check_three_term_recurrence_theorem(&legendre, 0, &[vec![q(1)], vec![q(0), q(1)]], &options)?;
assert!(report.is_proved());
let conclusion = report.conclusion.as_ref().unwrap();
assert_eq!(conclusion.interlacing, Some(RowInterlacing::Strict));
assert_eq!(conclusion.root_region.as_ref().unwrap().to_string(), "(-1, 1)");
report.verify()?;
# Ok::<(), RecurrenceError>(())
```

The parent provers of `certify_veronese_section_lift` try both theorems with
the region `(-inf, 0)` after the root-domain route.  The theorem statements,
the coordinate change and the root-region certificate are in
[algorithms.md](algorithms.md#real-rootedness-theorems-for-recurrences).

### The common lift summary

Every lift report carries a `lift_summary` (`recurrence::LiftSummary`): the
operation (`veronese_section` with modulus and residue, `symmetric_companion`
with `d_(n0)`, `hadamard_factor` with the factor, `row_factor` with `F` and
the exponent schedule, `gamma_polynomial` with the symmetry schedule and the
scale), the parent record, the identity status (`certified` with its
method, `supported_by_rows`, `unsupported`, `refuted`), the root transfer
(theorem, the parent root region and the input root region as
`RootInterval`s; the Veronese route needs parent roots in `(-inf, 0]`, the
symmetric-companion route an interlacing pair instead), the hypotheses with
their exact sign decisions, the parent proof (`proved` with a route, `open`,
`failed`, `attested`) and the verdict of `lift_verdict`: `proved` exactly when
the identity is certified, the parent is proved and every hypothesis is
verified, `route_failed` when the identity is refuted, the parent fails or a
hypothesis fails (which says nothing against the input).
`VeroneseSectionLift::lift_summary`,
`VeroneseSectionCertification::recompute_lift_summary` and
`SymmetricCompanionLift::recompute_lift_summary` compute it from the detailed
report, which stays authoritative.

The summary also records `input_interlacing` (`LiftInterlacingClaim`):
`proved` with a `LiftInterlacing` (kind `strict` or `weak`, orientation
`previous_interlaces_next` for `Q_(n-1) << Q_n`, `next_interlaces_previous`
or `unordered`, the first index `from_index`, and the theorem with its
premises), or `not_proved` with the reason.  It is an extra proved property:
the verdict does not depend on it.  The claim is transferred from the
interlacing that the parent's proving route concludes, which the parent
provers now record in `VeroneseParentProof::interlacing` (the window theorems
of route 2 and the three-term, Liu--Wang, derivative-lag and Ma--Wang
checks, all `P_(m-1) << P_m`), by one cited theorem per operation, and only
when the identity is certified and every transfer hypothesis is verified:

- operations independent of `n` that preserve real-rootedness on the whole
  line (the preserver-image claim, Hadamard multiplier sequences such as
  `1/k!` and `1/(alpha)_k`, basis maps whose theorem covers the whole line):
  the Hermite--Kakeya--Obreschkoff theorem, weak and `unordered`;
- Veronese sections (`VeroneseSectionCertification`): Athanasiadis and
  Wagner, arXiv:2404.12989, Corollary 5.6(b), weak, in the parent's
  orientation;
- row factors with a constant exponent: a common real-rooted factor keeps
  the cofactor's interlacing and orientation;
- `n`-dependent operators (`binom(n, k)`, `F^(n)`), gamma lifts and
  symmetric companions: `not_proved`.

### Parent provers and prepared rules

The building blocks the lifts share are public for applications that chain
lifts over their own parent records:

- `recurrence::prepare_lift_rule` validates a scalar rule with its supplied
  rows (positive offsets, homogeneous, rows reproduced, leader decided
  nonzero at every later index) and returns a `PreparedLiftRule` with the
  normalized rule, the first index and the identifying seeds;
- `recurrence::prove_parent_real_rootedness` runs the parent provers on a
  `ParentRecurrenceRecord` (rule, seeds, first index) and returns the
  `VeroneseParentProof`: every row zero or real-rooted with coefficients of
  one sign, for every `n`, with the route, or `not_proved` / `refuted`;
- `recurrence::find_bounded_root_region_report` returns a replayable
  `RealRootednessTheoremReport` proving every row real-rooted in a bounded
  interval (the Ma--Wang interval check, then the three-term theorem with a
  chosen bounded region), or `None`;
- `PreserverParentLibrary::add_parent_record` adds one named parent record
  to a preserver-image parent library, with its proved root region.

## Symmetric-companion lifts

`recurrence::find_symmetric_companion_lifts` recognizes a palindromic sequence
`A_n`, given by a scalar recurrence and seeds, as the symmetric part of a
parent `P_n = A_n + x B_n` with a first-order rule
`P_(n+1) = (a_n + b_n x) P_n + c_n x (1 - x) D P_n`, and certifies the
reconstruction for every index.  The parameters are searched as exact
rational polynomials of degree at most one in `n`; the symmetry schedule
`d_n = d_(n0) + (n - n0)` is supplied or inferred.  The input recurrence may
have derivatives; its left factor must be free of `x`, with no `(-1)^n`
factors and no inhomogeneous term.

```rust
use polytool::recurrence::*;

let q = |n: i64| BigRational::from_integer(n.into());
let term = |offset, derivative, rows: &[&[i64]]| {
    RecurrenceTerm::new(
        offset,
        derivative,
        BivarPoly::from_coefficients_normalized(
            rows.iter().map(|r| r.iter().map(|&c| q(c)).collect()).collect(),
        ),
    )
};
// The symmetric parts of the 1/2-Eulerian polynomials (target index m):
// A_m = [1 + (4m - 7) x] A_(m-1) + 4x(1 - x) A'_(m-1)
//       - [(6m - 15) x + (4m^2 - 22m + 30) x^2] A_(m-2)
//       - x(1 - x)[6 + (8m - 26) x] A'_(m-2) - 4x^2 (1 - x)^2 A''_(m-2).
let recurrence = Recurrence {
    terms: vec![
        term(1, 0, &[&[1, -7], &[0, 4]]),
        term(1, 1, &[&[0, 4, -4]]),
        term(2, 0, &[&[0, 15, -30], &[0, -6, 22], &[0, 0, -4]]),
        term(2, 1, &[&[0, -6, 32, -26], &[0, 0, -8, 8]]),
        term(2, 2, &[&[0, 0, -4, 8, -4]]),
    ],
    denominator: None,
    inhomogeneous: None,
};
let seeds = vec![vec![q(1), q(1)], vec![q(1), q(7), q(1)]]; // A_2, A_3
let report = find_symmetric_companion_lifts(
    &recurrence,
    &seeds,
    2,
    &SymmetricCompanionLiftOptions::default(),
)?;
let lift = report.proving_lift().expect("certified, theorem applies");
// P_(n+1) = (1 + 2n x) P_n + 2x(1 - x) P_n', P_2 = 1 + 2x, B_2 = 1.
assert_eq!(lift.certificate.companion_seed, vec!["1".to_string()]);
let replay = verify_symmetric_companion_certificate(&lift.certificate, 10)?;
assert!(replay.is_certified() && replay.theorem_applies());
# Ok::<(), RecurrenceError>(())
```

Each lift carries a replayable `SymmetricCompanionCertificate` (source and
seeds, schedule, the parameter witnesses `a, b, c, delta, epsilon` as
polynomials in `n`, the eliminated operator with its derivation from the
source rule, the companion seed and the parent as a `polytool.recurrence.v1`
record) and the verifier's statuses, kept apart: the identity (`certified`,
`candidate`, `unsupported` for a vanishing coupling, `rejected`), the
real-rootedness route through the preservation theorem (`theorem_applies`
only when `a_n >= 0`, `0 <= delta_n <= c_n`, `d_n >= 1` and the initial pair is
nonnegative with `B <= A`), and the Lean replay (`not_attempted`).  A failed
hypothesis means only that this route fails.  The report serializes as
`polytool.recurrence-symmetric-companion.v1` (exact rationals as strings).
The algorithm and the soundness argument are in
[algorithms.md](algorithms.md#symmetric-companion-lifts).

## Hadamard-factor lifts

`recurrence::find_hadamard_factor_lifts` writes the rows
`Q_n = sum_k q(n, k) x^k` of a recurrence-defined sequence as
`Q_n = lambda(n) sum_k u(n, k) p(n, k) x^k` with a known factor `u` and a
parent `P_n = sum_k p(n, k) x^k` whose rule is smaller than the input's (lower
order, or the same order with smaller degrees).  The menu
(`HadamardFactorLiftOptions`) has `binom(n + a, k + b)` with `|a|, |b| <= 2`,
`k!`, `1/k!`, `(n - k)!`, `1/(n - k)!`, `(alpha)_k` and `1/(alpha)_k` for the
listed `alpha`, and optionally products of two (`HadamardFactor`,
`HadamardFactorTerm`).  Undetermined quotient coefficients (where `u = 0`) are
set to zero or continued by a fitted hypergeometric ratio, and the rows may be
normalised by their leading coefficients `lambda(n)` (a
`VeroneseSectionGauge::Ratio`).

```rust
use polytool::oeis;
use polytool::recurrence::{find_hadamard_factor_lifts, HadamardFactorLiftOptions};

// Narayana (A001263, rows sum_k N(n+1, k+1) x^k): u = binom(n+1, k+1),
// P_n = (1 + x)^(n+1) / (n + 1), proved end to end.
let narayana = oeis::by_id("A001263").unwrap();
let (rule, rows) = narayana.recurrence_parts().unwrap();
let report = find_hadamard_factor_lifts(
    &rule,
    &rows,
    narayana.recurrence_first_index,
    &HadamardFactorLiftOptions::default(),
)?;
let lift = report.proved_lifts().next().unwrap();
assert_eq!(lift.factor_text, "binom(n+1, k+1)");
report.verify()?;
```

A lift (`HadamardFactorLift`) records the factor, gauge, parent
(`polytool.recurrence.v1`), the identity status
(`HadamardFactorIdentityStatus::{Certified, SupportedByRows, Refuted}`, the
certificate `HadamardFactorCertificate` replayable by `verify`), the
multiplier-sequence justification of the factor (`None` when no transfer
theorem is recorded, as for `k!` or `(alpha)_k`), the parent provers' record
and the `lift_summary`.  Rejected factors (`HadamardFactorRejection`: invalid
for the domain, vanishing on the support, no smaller parent, parent budget)
are listed as near misses.  The report is
`polytool.recurrence-hadamard-factor-lifts.v1`.  When the division by the
parent leaves a remainder, it is removed by a row annihilator of the parent:
that of a product parent, or, for any first-order parent, an operator `A(n)`
certified by an exact invariance identity under the parent step (fields
`row_annihilator_step_multiplier`, `row_annihilator_step_quotient`,
`row_annihilator_start_base`, `row_annihilator_restart_bases` of the
certificate), or, for a parent of order two or more, an invariant annihilator
system on the row windows `(P_n, ..., P_(n+k-1))` (field
`row_annihilator_system`: generators `A^(l) = sum_(i < k) A^(l)_i(n) S^i` with
their invariance identities under the companion step, the quotients `Z_l`, a
multiplier `mu_2(n, x)`, the start and restart bases).  This proves the Baxter
triangle A056939 (factor `binom(n+2, k+1)`, hypergeometric parent proved by
Liu--Wang, or the second-order Jacobi parent proved by the three-term theorem)
and certifies the `binom(n+1, k)` identities of the descent polynomials of
321-avoiding permutations and the peak polynomials of 123-avoiding
permutations, whose parents are second-order rules with derivative terms.  The
certificate, its boundary conditions and the transfer theorem are in
[algorithms.md](algorithms.md#hadamard-factor-lifts).

## Row-factor lifts

`recurrence::find_row_factor_lifts` finds a fixed polynomial `F` and a
schedule `e(n) = alpha n + beta` with `Q_n = F^(e(n)) S_n` (common irreducible
factors of the first rows, exact multiplicities), derives the cofactor rule
by substitution, and proves the cofactors with the parent provers.  A lift
(`RowFactorLift`, report `polytool.recurrence-row-factor-lifts.v1`) records
`F`, `e`, whether `F` is real-rooted (with roots `<= 0`), the cofactor rule and
seeds, the identity (`RowFactorCertificate`, replayable), the cofactor proof
and the `lift_summary`; schedules that fail (a row not divisible by
`F^(e(n))`, a cofactor rule with an `x`-dependent leader) are listed in
`rejected`.  See [algorithms.md](algorithms.md#row-factor-lifts).

## Gamma lifts

`recurrence::find_gamma_lifts` recognises the rows as rescaled gamma
polynomials of a palindromic parent, `Q_n = lambda(n) Gamma_(n+s)(c z)` with
`P_m = sum_j Gamma_(m,j) x^j (1 + x)^(d_m - 2j)`.  Parents come from a menu
(`GammaLiftParent::standard_menu`: Eulerian, type-B Eulerian, Narayana; add
`GammaLiftParent::one_over_k_eulerian_symmetric_part(&k)` or your own rule
with its symmetry schedule); detection reads `c` and `lambda` (constant or of
constant ratio) from the rows, and every match is certified by
`certify_gamma_lift`, which callers can also run on a candidate of their own.

```rust
use polytool::recurrence::*;
use polytool::sequences::literature::peak_polynomials_recurrence;

// Peak polynomials: W_n = 2^(n-1) Gamma(A_n)(z/4).
let (rule, seeds, first_index) = peak_polynomials_recurrence();
let report = find_gamma_lifts(&rule, &seeds, first_index, &GammaLiftSearchOptions::default())?;
let lift = report.proved_lifts().next().unwrap();
assert_eq!(lift.candidate.parent.name, "eulerian");
assert_eq!(lift.lift_summary.root_transfer.input_root_region.to_string(), "(-inf, 0]");
# Ok::<(), RecurrenceError>(())
```

A certification (`GammaLiftCertification`, schema
`polytool.recurrence-gamma-lift.v1`; the search report is
`polytool.recurrence-gamma-lifts.v1`) records the parent symmetry status (the
reversal induction), the transported gamma rule, the identity
(`GammaLiftCertificate`, replayable by `verify()`), the parent proof and the
`lift_summary` (new `LiftOperation::GammaPolynomial`).  The transfer is the
gamma criterion: a palindromic `P` is real-rooted with roots `<= 0` if and only
if its gamma polynomial is, so the input roots lie in `(-inf, 0]` for `c > 0`.
See [algorithms.md](algorithms.md#gamma-lifts).

## Preserver-image lifts

`recurrence::find_preserver_image_lifts` recognises the rows as images
`Q_n = lambda(n) T_n(P_(n+s))` of a proved parent `P` under a
real-rootedness preserver `T_n` from a fixed menu (`PreserverOperator`):
Euler-operator polynomials `f(n, theta)` (multiplier sequences `k -> f(n, k)`),
Laguerre--Pólya operators `phi(n, D)` such as `P + c(n) P'`, `D^j o x^m` and
`x^m D^j`, affine substitutions `P(a x + b)`, reversals `x^(d_n) P(1/x)`, Möbius
maps such as the `f`-polynomial, and hypergeometric multipliers (`1/k!`,
`1/(alpha)_k`, `binom(n + a, k + b)`, ...).  Parents come from a
`PreserverParentLibrary` (the named families of `sequences` with their
defining recurrences, optionally the OEIS catalog entries the parent provers
prove; the proofs are computed once per library, so share one library across
calls).  The menu also has the Euler step
`kappa (theta + a)(theta + b) + x (u + v theta)` (RealRooted
`IsPFPolynomial.generalStep_of_natDegree_le`; it maps polynomials with
nonnegative coefficients and real roots to such polynomials when
`kappa, a, b, u > 0`, `v >= 0` and the comparison defect is negative on
`s <= 0`), fitted from five linear coefficients.  A parent may also be proved
real-rooted in a bounded interval `[L, U]` (`PreserverParentRegion::Interval`,
`PreserverParent::interval_proof`): by the Ma--Wang interval check, or by the
three-term theorem with a root region.  Affine, reversal and Möbius maps and
the Rolle operators then give the exact image of `[L, U]` as the input region
(for example `T_n(2x - 1)` for the Chebyshev polynomials `T_n` is rooted in
`(0, 1)`), and the interval report is recorded in the certification
(`interval_parent_proof`) and replayed by `verify()`.  Detection fits the operator parameters exactly from a few rows (row
by row linear algebra, then rational functions of `n`); every match is
certified by `certify_preserver_image_lift`, which callers can also run on a
candidate of their own.

```rust
use polytool::recurrence::*;
use std::sync::Arc;

let library = Arc::new(PreserverParentLibrary::named_families()?);
let options = PreserverImageSearchOptions { library: Some(library), ..Default::default() };
// `rule`, `rows` (rows[0] is Q_(first_index)) and `first_index`: the input.
let report = find_preserver_image_lifts(&rule, &rows, first_index, &options)?;
for lift in report.proved_lifts() {
    println!("{} from {}", lift.candidate.operator.text(), lift.candidate.parent_name);
}
# Ok::<(), RecurrenceError>(())
```

A certification (`PreserverImageCertification`, schema
`polytool.recurrence-preserver-image.v1`; the search report is
`polytool.recurrence-preserver-image-lifts.v1`) records the identity (a
`PreserverImageCertificate`: pseudo-division of the gauged input relation
composed with, or conjugated through, the operator by the shifted parent
relation, with row annihilators for a remainder, or the Hadamard conjugation
for a hypergeometric multiplier; `verify()` regenerates the input and parent
rows from the recorded rules and seeds and re-checks the identifying rows and
exceptional bases), the parent proof (roots `<= 0`, a bounded interval, or
real-rootedness on the whole line by the three-term theorem), the Lean route and RealRooted
declarations of the operator's theorem, and the `lift_summary` (new
`LiftOperation::PreserverImage`) whose hypotheses are the side conditions of
the theorem, decided exactly for every `n` (fixed polynomials exactly,
index-dependent coefficients of degree `<= 2` by the sign checker; anything
else stays open).  It also records a `PreserverInterlacingClaim`: consecutive
input rows interlace (weakly, orientation not determined) only for an operator
independent of `n` that preserves every real-rooted polynomial, from the
parent's interlacing proved by a replayed theorem report
(Hermite--Kakeya--Obreschkoff); otherwise `none`.  See [algorithms.md](algorithms.md#preserver-image-lifts).

A candidate is also an image specification of `image_recurrence_of_operator`
(`PreserverImageCandidate::image_operator_specification`; operators convert
with `ImageRowMap::try_from`, every operator except the multipliers that
depend on `n`), and `PreserverImageCandidate::image_recurrence` computes the
runnable rule of its image from the parent rule.  With
`PreserverImageOptions::image_rule` set, a certified identity also carries
that rule (`PreserverImageCertification::image_rule`), which `verify()`
replays and compares with the input rows.

## Basis maps and basis-map lifts

A `BasisMap` is the linear map `T_P(sum p_k u^k) = sum p_k P_(i0 + k)(t)` of a
family given by a rule, its first index `i0` and seed rows (serialized as
`{ "name", "family": polytool.recurrence.v1 }`).  `apply` evaluates it,
`expand` inverts it exactly on triangular rows (strictly increasing degrees),
`tail(s)` shifts the first index, and `intertwining_relation()` gives the
identity `sum t^i D^d T(phi(theta) u^m p) = 0` derived from the rule (for a
rule free of `n`, `T(u^R p) = sum_j A_j T(u^j p)`).

`check_basis_map_preserves_pf` recognises the family's shape and decides the
hypotheses of the matching theorem exactly; each proved `BasisMapConclusion`
says that every nonzero polynomial with all roots in `input_interval` is sent
to zero or to a real-rooted polynomial with all roots in `output_interval`,
in every degree.  Interval ends are exact `QuadraticNumber`s `a + b sqrt(d)`.

```rust
use polytool::recurrence::*;

let map = BasisMap::from_oeis("A035607")?;
let report = check_basis_map_preserves_pf(&map, &BasisMapTheoremOptions::default())?;
let pf = report.pf_conclusion().unwrap();
assert_eq!(pf.input_interval.to_string(), "[-1 - sqrt(2), 0]");
// T((u + 4)^3) = 125 + 150 t + 32 t^2 + 2 t^3 is not real-rooted: -4 is outside.
# Ok::<(), RecurrenceError>(())
```

`find_basis_map_lifts(rule, rows, first_index, &options)` writes the input as
`Q_n = T_P(R_n)` for a map of a `BasisMapLibrary` (default: the admitted OEIS
catalog maps; build one with `BasisMapLibrary::oeis_catalog()` or
`BasisMapLibrary::build` and share it): it expands the rows, rejects maps whose
coefficient rows have a root outside every source interval (exactly), fits a
derivative-free rule for `R` and certifies the candidate with
`certify_basis_map_lift`.  The certification (`BasisMapLiftCertification`,
schema `polytool.recurrence-basis-map-lift.v1`) records the identity (a
`BasisMapDivisionCertificate`, replayed by `verify()`), the map theorem, the
conclusion used, the root proof of `R` (`BasisMapCoefficientRootProof`: the
parent provers on `R_n(v + U')` and `R_n(L' - v)`) and the `lift_summary`
(`LiftOperation::BasisMap`).  `apply_basis_map_to_family(&map, &family,
&options)` is the forward move: the rule of `T_P(R)` (transported for maps of
order one; for higher order computed from the transported system by
`ore_algebra::find_component_annihilator`, `BasisMapImageConstruction::Eliminated`,
and fitted only when that finds nothing within its budget) with the same
certification.  The summary's `root_transfer.input_root_region` is the map
theorem's output interval as a `RootInterval`: each irrational end records its
exact value (`RootEndpoint::exact`, an `ExactQuadraticEndpoint`) next to a
rational bound strictly outside it, so `format_root_interval_exact` prints
`[-3-2*sqrt(2), -3+2*sqrt(2)]` for the A008288 (Delannoy) map while readers of
the rational values see a slightly larger correct interval; `verify()` also
accepts records written without the exact form.  See
[algorithms.md](algorithms.md#basis-map-lifts).

## Exact sign decisions for an integer index

`decide_index_rational_function_sign(p, q, &range, condition)` decides
`p(n) / q(n) R 0` for every integer `n` of an `IntegerRange`, exactly.  A
vanishing denominator is a failure, a failure names the smallest failing
index, and `undecided` happens only over the `IndexSignBudget`.

```rust
use num_bigint::BigInt;
use polytool::{
    decide_index_rational_function_sign, BigRational, IndexSignCondition, IndexSignDecision,
    IntegerRange,
};

let q = |n: i64| BigRational::from_integer(BigInt::from(n));
// (n^2 - 5n + 7) / (2n - 5) > 0 for n >= 3: the numerator has no real root
// and the denominator changes sign between 2 and 3.
let decision = decide_index_rational_function_sign(
    &[q(7), q(-5), q(1)],
    &[q(-5), q(2)],
    &IntegerRange::starting_at(3),
    IndexSignCondition::Positive,
);
let IndexSignDecision::Holds { certificate } = &decision else { unreachable!() };
certificate.verify().unwrap();
// From n = 0 it fails at n = 0 (value -7/5).
let decision = decide_index_rational_function_sign(
    &[q(7), q(-5), q(1)],
    &[q(-5), q(2)],
    &IntegerRange::starting_at(0),
    IndexSignCondition::Positive,
);
assert_eq!(decision.failing_index(), Some(&BigInt::from(0)));
```

The algorithm is in [algorithms.md](algorithms.md#exact-sign-decisions-for-an-integer-index).

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
