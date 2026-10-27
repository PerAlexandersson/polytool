# Algorithm notes

## Real-rootedness algorithms

The default real-rootedness check uses an adaptive exact path from
`root_count`. Degrees at most four are decided directly by exact discriminant
and quartic-invariant signs, including the repeated-root cases; they do not
construct a remainder sequence or square-free quotient. For higher degrees,
the signed PRS of `p` and `p'` simultaneously counts
the distinct real roots and ends in their gcd, whose degree determines the
square-free degree. Thus the general check avoids constructing the square-free
quotient and running a second Euclidean pass. The polynomial is real-rooted
precisely when the root count equals the square-free degree. Public operations
that do need the quotient use the same fraction-free primitive integer GCD.
Degree five omits content normalization because its short chain does not repay
the scaling bookkeeping. From degree six onward, a signed Collins--Brown
subresultant PRS divides out predictable exact factors, avoiding both
coefficient swell and a coefficient-content GCD after every pseudo-division.
The root-counting path streams this sequence through two active polynomials
and one reusable pseudo-remainder buffer. It records only the requested Sturm
endpoint variations and the final gcd degree rather than retaining the full
chain. Pseudo-division also moves out the leading coefficient and skips the
large-integer multiply/subtract pair that would cancel it exactly.

For one-signed coefficient polynomials, the backend first tries cheap
coefficient filters (Kurtz, sufficient; Newton, necessary) and then counts the
distinct real roots with the subresultant PRS: after removing powers of `t`,
such polynomials can only have negative real roots, so only the signs at
`±∞` are needed. Subresultant PRS is the general backend. An earlier rule
switched large palindromic inputs (Eulerian, type-B Eulerian, Touchard) to
Uspensky; after the Collins--Brown PRS that rule made those families about
10--26 times slower at degree 79, so it was removed. `polytool bench
real-rootedness` tracks these cases.

For algorithm comparisons, `is_real_rooted_uspensky_bigint_coeffs` uses an
independent exact Uspensky/Descartes path.  It applies a strict Fujiwara root
bound, separates roots into dyadic magnitude bands (using the reciprocal
polynomial below `1`), and applies homographic subdivision with Descartes sign
variations.  It does not use finite fields, floating point, or approximate
roots. The explicit PRS and Uspensky entry points remain available for direct
comparison.

Explicit matrix comparison paths remain available as
`is_real_rooted_bezout`, `is_real_rooted_bezout_squarefree`, and
`is_real_rooted_hermite`.

All interlacing checks use the **Bézout matrix** (Fisk, Cor. 9.145). For
polynomials f (degree d) and g (degree d−1), the Bézout matrix B(f,g) is the
d×d symmetric matrix with (i,j) entry equal to the coefficient of x^i y^j in
`(f(x)g(y) - f(y)g(x)) / (x-y)`.

**Theorem:** f and g are both real-rooted and g strictly interlaces f
if and only if B(f,g) is positive definite.

This reduces interlacing to a single exact matrix definiteness check, avoiding
root isolation entirely. It is 100–400× faster than Sturm chains at degree 15+
for interlacing checks.  The definiteness implementation uses fraction-free
BigInt Bareiss elimination below dimension `30`, and switches to a
CRT-over-prime-fields path for larger matrices.  The modular path reconstructs
the leading principal minors exactly using Hadamard bounds.  For speed
comparisons, the explicit entry points are `is_positive_definite_bareiss`,
`is_positive_definite_modular`, and `modular_leading_principal_minors_bigint`.

Run the local comparison example with:

```bash
cargo run --release -p polytool --example bench_bareiss_vs_modular_linalg
```

For **same-degree** polynomials, `check_weak_interlacing` reduces to the
deg+1 case by extending one polynomial with a root far to the right
(using the Cauchy bound for root radius). When all coefficients are
positive (all roots negative), multiplying by t suffices, skipping the
Cauchy bound computation.

See `doc/bezout-interlacing.md` for the Mathematica reference implementation.

## Interlacing matrices

Athanasiadis--Wagner's interlacing matrix construction packages a polynomial
sequence or polynomial matrix into an infinite block Toeplitz matrix `Lace(A)`.
For a column vector, total nonnegativity of this infinite matrix is called
full interlacing and implies pairwise interlacing of the entries.  The converse
fails, so this gives a stronger finite-experiment target than checking
ordinary pairwise interlacing.

The crate implements finite truncations through `lace_matrix_i64`,
`lace_matrix_bigint`, and the sequence wrappers
`lace_matrix_sequence_i64` / `lace_matrix_sequence_bigint`.  Exact finite TNN
checks are available with `check_lace_total_nonnegative_i64` and the faster
Neville-elimination wrappers.

## Recognizing weighted renewal recurrences

`renewal::recognize_weighted_renewal(&Recurrence)` recognizes
`P_n = sum_j r_j t^j P_(n-1-ell*j)`, with positive integer `ell`,
nonnegative rational `r_j`, `r_0=1`, and an exactly real-rooted seed
`R(z)=sum_j r_j z^j`. It normalizes a constant LHS and combines duplicate
terms. Derivatives, alternating signs and index-dependent coefficients are
conservatively rejected. Failure means this criterion did not apply.

The returned `WeightedRenewal` provides the canonical generator `rows(count)`
for `1/(1-x R(t x^ell))`, and checks supplied initial rows with
`matches_initial_rows(first_index, rows)`. At least `1+ell*degree(R)` rows,
starting at zero, must match. Its paper theorem gives oriented weak
interlacing `P_n << P_(n+s)` for `degree(R)*s <= degree(R)*ell+1`.
`max_interlacing_step()` returns this bound, or `None` for the constant
family. Repeated roots are allowed; no numerical root finding is used.

An example accepts the existing recurrence JSON format on stdin:

```sh
cargo run -p polytool --example recognize_renewal < recurrence.json
```

Its JSON distinguishes structural recognition from matching initial data.
Neither proves that a recurrence guessed from finite rows holds for the
intended combinatorial sequence. The all-index identity remains a caller
obligation: this Rust feature is not a formal certificate.

## Variable-predecessor recurrences

`renewal::VariablePredecessorRecurrence` represents
`Q_0 = 1`, `Q_n = Q_(n-1) + w(n) t Q_(p(n))` (`n >= 1`) with the predecessor
rule and the weight kept apart from the fixed-offset `RecurrenceTerm`. If
`w(n) >= 0`, `0 <= p(n) < n` and `p` is nondecreasing, every `Q_n` has
nonnegative coefficients and constant term 1, is real-rooted, and
`Q_i << Q_n` for `p(n) <= i <= n`; for `p(n) = max(0, n-L)` this gives
`Q_A << Q_B` whenever `0 <= A <= B` and `B - A <= L` (a window of `L+1` rows,
not derived from transitivity). The paper proof is in
[RealRooted#2037](https://github.com/PerAlexandersson/RealRooted/issues/2037)
and is not Lean checked.

Supported symbolic forms are `PredecessorRule::FixedDelay { delay: L }` with
`PredecessorWeight::Polynomial` (coefficients in ascending powers of `n`) or
`PredecessorWeight::Periodic` (`w(n) = values[(n-1) mod period]`).
`certify()` returns a `VariablePredecessorCertificate` or an error naming a
violating index:

- the delay must satisfy `L >= 1` (symbolic: `p` is nondecreasing and
  `0 <= p(n) < n`);
- a polynomial weight is proven nonnegative on every integer `n >= 1`, never
  sampled: the certificate stores the least integer `c >= 1` such that
  `w(c+m)` has nonnegative coefficients in `m` (so `w >= 0` on `[c, oo)`) and
  the exact values `w(1), ..., w(c-1)`. A negative leading coefficient yields a
  witness beyond the Cauchy root bound. Certificates needing `c > 65536` are
  rejected as not certified;
- a periodic weight has each residue value checked.

`PredecessorRule::Explicit` and `PredecessorWeight::Explicit` lists certify only
the supplied prefix: the certificate scope is `prefix` with `last_index`, the
conclusion says so, and `certify_all_indices()` rejects them.
`verify_certificate` replays a certificate (any valid shift is accepted).
`rows(count)` is the canonical generator (exact `BigRational` rows), and
`recognize_fixed_delay_predecessor(&Recurrence)` recognizes
`P_n = P_(n-1) + w(n) t P_(n-L)` with `w` polynomial in `n`. Note that the
canonical family uses `Q_(max(0, n-L))`, so for `n < L` it adds `w(n) t Q_0`
where a fixed-offset recurrence would use a zero row.

## Conditional theorem matches

`renewal::match_recurrence_theorems(&Recurrence, first_index, rows)` checks a
found or supplied recurrence against `weighted_renewal_window` and
`fixed_delay_predecessor_window` and returns explicit conditional matches plus
rejection reasons. `RecurrenceTheoremReport` (schema
`polytool.recurrence-theorems.v1`) wraps them with the recurrence record and
its `finite_evidence` (`adaptive_search_fit` or `supplied_recurrence`, status
`finite_evidence_only`). Each match lists the exactly checked hypotheses, the
conclusion, the window (`d*s <= d*ell+1`, or `B - A <= L`), the boundary
contract with a finite comparison of the supplied rows (`canonical`,
`canonical_after_relabel` when the first supplied row plays index 0,
`mismatch` or `insufficient_rows`), `remaining_obligations` (always including
the all-`n` validity of the recurrence) and `proof_status:
"paper_proof_not_lean_checked"`.

```sh
polytool recurrence --theorem-matches < rows.txt          # text report
polytool recurrence --json --theorem-matches < rows.txt   # theorem report JSON
polytool recurrence --json < rows.txt | polytool recurrence-theorems --recurrence -
```

Default `recurrence` output is unchanged: a match is only noted on stderr, and
`--json` still emits the plain `polytool.recurrence.v1` record. With
`--theorem-matches --json` the record is nested as `recurrence_json`.
`recurrence-theorems` exits 0 with a match, 1 without, 2 on invalid input. The
MCP `find_recurrence` response gains the additive `theorem_matches` field, and
the MCP tool `recurrence_theorems` checks a supplied record.
