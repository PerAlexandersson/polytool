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
checks are available with `check_lace_totally_nonnegative_by_minors_i64` and the faster
Neville-elimination wrappers (`check_lace_totally_nonnegative_neville_i64`,
`check_lace_totally_nonnegative_neville_bigint`).  They return `InterlacingMatrixError`;
`InterlacingMatrixError::NotTotallyNonnegative` wraps the `TotalPositivityError`
naming the failing minor or Neville elimination step (the same error type as
`linalg::check_totally_nonnegative_neville` and `linalg::check_total_positivity`).

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

## Veronese-section lifts

`recurrence::find_veronese_section_lifts(&Recurrence, initial_rows,
first_index, &VeroneseSectionLiftOptions)` recognizes a recurrence-defined
sequence `Q_n(y)` as a coefficient section of a parent sequence `P_n(x)`:

```text
V_(r,j)(P)(y) = sum_k [x^(r k + j)] P(x) y^k,   Q_n = V_(r,j)(P_n)  (n >= n0).
```

The transfer theorem: if `P` is real-rooted with coefficients of one sign
(all roots in `(-inf, 0]`), every section is zero or real-rooted with
nonpositive roots.  The coefficient sequence of `P` is a Pólya frequency
sequence (Aissen--Schoenberg--Whitney, Edrei), and the Toeplitz matrix of
`a_(r k + j)` is the submatrix of the Toeplitz matrix of `a` with rows `r i`
and columns `r k + j`, so it is totally nonnegative as well.

Supported input: derivative-free terms without `(-1)^n` factors and with
coefficients independent of `n`, an absent or nonzero constant denominator, no
inhomogeneous term.  With `Q_m = sum_k c_k(y) Q_(m-k)` and order `d`, the
characteristic polynomial is `chi(L, y) = L^d - sum_k c_k(y) L^(d-k)`.

A parent with monic characteristic polynomial `G(L, x)` (`P_(n+e) = -sum_i
g_i(x) P_(n+i)`) has all `r` sections satisfying the recurrence with
characteristic polynomial `N = Norm_r(G) = prod_(w^r=1) G(L, w x)`, written in
`y = x^r`: `N(E, x^r)` is a multiple of `G(E, x)` with coefficients in
`Q[x^r]`, which commute with taking sections.  `N` is computed exactly as the
determinant of multiplication by `G` on `1, x, ..., x^(r-1)` in
`Q[L, y][x]/(x^r - y)` (cofactor expansion, `r <= 6`).

For each modulus `r = 2..max_modulus` (at most `d`):

1. factor `F(L, x) = chi(L, x^r)` over `Q` (`polytool::factorization`); a
   failure records the modulus in `skipped_moduli`;
2. enumerate the monic divisors `G` of `F`, nonconstant with
   `r deg_L G <= d`, by increasing `L`-degree and then `x`-degree, at most
   `max_parent_candidates` per modulus (`exhausted_budget` when reached);
3. keep `G` when `N` divides `chi` in `Q[y][L]` (the quotient is the reported
   `cofactor`);
4. for each residue `j`, solve the `r e` seed sections `s_(i,a)(y)` of
   `P_(n0+i) = sum_a x^a s_(i,a)(x^r)` from `V_j(P_n) = Q_n` on the first `d`
   rows.  The map is linear with polynomial matrix `O(y)` (column `(i,a)` is
   `V_j` of the rows from the basis seed `x^a` at position `i`), solved by
   fraction-free elimination with row pivoting over `Q[y]`; the seeds must be
   polynomials (`seeds_not_polynomial`), rank below `r e` is not supported
   (`seeds_not_determined`), and the rows from the solved seeds must match all
   `d` identifying rows (`row_mismatch`);
5. check the first `parent_rows_to_check` parent rows exactly (nonzero, one
   coefficient sign, real-rooted after clearing denominators);
6. an order-one parent `G = L - c(x)` gives `P_n = c^(n-n0) P_(n0)`; when `c`
   and `P_(n0)` are real-rooted with coefficients of one sign the lift is
   `proved_product_parent`.  Order-two and higher parents with passing rows
   are `parent_real_rootedness_open` (finite evidence only); a failing row
   gives `parent_row_zero`, `parent_row_mixed_signs` or
   `parent_row_not_real_rooted` with its index.

Why `d` matching rows prove the identity for all `n`: `S_n = V_j(P_n) - Q_n`
satisfies the `chi` recurrence at every target `>= n0 + d` (the sections
satisfy `N`, and `N | chi`), `chi` is monic of order `d`, and `S` vanishes at
the `d` identifying indices, so `S = 0`.  A few further rows are compared as a
defensive self-check.

Example: `Q_(n+2) = 2 Q_(n+1) + (y - 1) Q_n`, `Q_0 = Q_1 = 1` has
`chi = (L - 1)^2 - y`, and `chi(L, x^2) = (L - 1 - x)(L - 1 + x)`.  The norm of
`L - 1 - x` is `chi` itself, the seed is `P_0 = 1`, and the parent is
`(1 + x)^n`: proved.  The conjugate `L - 1 + x` (parent `(1 - x)^n`) is listed
with `parent_row_mixed_signs` at index 1.  The trivial lift `Q_n(x^r)` has a
norm of degree `r d > d` and never appears.  For the even sections of the
Fibonacci polynomials (`chi = L^4 - 2L^3 + L^2 - y`) the parent
`G = L^2 - L - x` with seeds `1, 1` is found with an open status.

Limitations of this version: no `n`-dependent or differential recurrences,
modulus at most 6, undetermined seeds are rejected rather than parametrized,
and only order-one parents are proved.  Natural next steps are proved parent
templates of higher order (for example interlacing-preserving matrix
recurrences), parameterized seeds when `O(y)` is rank deficient, and an exact
replay format for the identification certificate.  The CLI command is
`polytool recurrence-veronese-lifts`; the MCP tool is
`find_veronese_section_lifts`.

### Certifying one candidate (v2)

`recurrence::certify_veronese_section_lift(&input, rows, n0, &candidate,
&options)` certifies a candidate `Q_n = lambda(n) V_(r,j)(P_(n+s))` for an
input with index-dependent coefficients, derivatives on earlier rows and a
denominator `f(n, y)` (its leading row carries no derivative), and a parent
whose leader is free of `x`.  Bounds: `r <= 4`, `|s| <= 4`, derivative order
at most 6, at most 10000 operator terms, at most 64 compared rows; hitting one
is a typed error.

**Transport.**  Write `P(x) = sum_(a<r) x^a F_a(x^r)`, `y = x^r`; an
operator matrix acts by `(A F)_b = sum_a A[b][a] F_a`.  Multiplication by `x`
is `X` with `X[0][r-1] = y` and `X[b][b-1] = 1`; differentiation is `DD` with
`DD[b][b+1] = (b+1) + r y D_y` (`b < r-1`) and `DD[r-1][0] = r D_y`; `n` and
`S` act diagonally.  Since `DD X - X DD = I`, the map
`phi(sum c(n, x) D^d S^i) = sum c(n, X) DD^d S^i` is a homomorphism of Ore
algebras, and the sections of `(L P)_n` are `(phi(L) F)_n`.  In Euler form,
`x^s a(n, theta_x)` sends component `b` to `(b + s) mod r` with the factor
`y^floor((b+s)/r) a(n, r theta_y + b)`.  These conventions are tested against
direct section computations on rows.

**Gauge.**  With `lambda(n)/lambda(n-1) = a(n)/b(n)` (checked nonzero for
every `n > n0`, and `lambda(n0) != 0`), the parent relation
`sum c_(i,d)(n, x) D^d P_(n+i) = 0` (leader `g(n) S^k`) becomes
`sum c_(i,d)(n+s, x) prod_(t=i+1..k) a(n+t) prod_(t=1..i) b(n+t) D^d P~_(n+i) = 0`
for `P~_n = lambda(n) P_(n+s)`; its leader stays free of `x`.

**Pseudo-division.**  `e_j^T M` (`M` the input rule as a relation at base `n`)
is divided by `Phi = phi(L~)`, whose shift-`k` part is `g~(n) S^k I`: while
the remainder has a shift `m >= k`, multiply it (and the recorded quotients)
by `g~(n + m - k)` and subtract `c o S^(m-k) o Phi` for its top coefficient
row `c`.  The result is `w(n) e_j^T M = sum_nu H_nu S^nu Phi + R`, `R` of
shift `< k`, rechecked with exact operator arithmetic.

**Soundness.**  The first rows (at least the `K` identifying rows) are
compared exactly; a difference refutes the candidate at that index.  With
`R = 0`, `w(n) (M V_j(P~))_n = sum_nu H_nu (phi(L~) F~)_(n+nu) = 0` at every base
`n` where `L~` holds (`n >= max(n0, p0 + Kp - k - s)`).  `w` depends on `n`
only; its integer roots, and the bases between `n0 + K - R` and that bound,
are evaluated exactly on the compared rows.  So `V_j(P~)` satisfies `M`
wherever `Q` does (bases `n >= n0 + K - R`); the difference `S = Q - V_j(P~)`
then satisfies `f(m, y) S_m = (earlier rows)` for every target `m >= n0 + K`,
`f(m, .)` is never the zero polynomial there (checked from the exact integer
roots of its coefficients), and `S` vanishes at the identifying rows, so
`S = 0`.  With `R != 0` the identity is only `supported_by_rows`.

**Parent.**  The parent statement (every used `P_m` is zero or real-rooted
with coefficients of one sign) is proved by a periodic product
`P_(n+k) = c(x) P_n`, by `match_recurrence_theorems` with a canonical
boundary, or by a root-domain claim of `transform_recurrence` (maps: identity,
`f_polynomial` and `reversal` with a degree schedule fitted to the first rows,
the inverse `f`-map; accepted only for an unconditional all-index claim that
the source is real-rooted in an interval ending at or below 0).  The Eulerian
parent A008292 is proved by the `f_polynomial` route.  A caller attestation
matching the parent record exactly gives a conditional verdict.  A checked
parent row that fails refutes the parent statement.

Examples (tests): the even sections of the Eulerian polynomials satisfy a
monic order-3 successor recurrence with second derivatives on the earlier
rows; its forward certificate and the `f_polynomial` parent proof give
`proved`.  The even Narayana sections satisfy an order-4 recurrence with a
leader depending on `y` (multiplier `w(n) = (n+3)(n+4)(n+5)`); the identity is
certified and the order-two Narayana parent stays open.  Gauged and shifted
binomial sections (`V_(2,0)((1+x)^n)/n!`, `2^n V_(3,1)((1+x)^(n+2))`) are
proved.

Limitations: inputs whose leading row carries a derivative, the reverse
(reconstruction) certificate, candidate search and Lean replay are not part
of this stage.

## Symmetric-companion lifts

`recurrence::find_symmetric_companion_lifts(&Recurrence, initial_rows,
first_index, &SymmetricCompanionLiftOptions)` inverts the symmetric
decomposition: given a scalar recurrence `M` and seeds for a palindromic
sequence `A_n` (symmetric at `d_n = d_(n0) + (n - n0)`, an ambient index that
allows degree deficits and zero rows), it looks for a parent rule

```text
P_(n+1) = (a_n + b_n x) P_n + c_n x (1 - x) D P_n,   P_n = A_n + x B_n.
```

With `delta_n = b_n - a_n - c_n d_n`, `epsilon_n = c_n - delta_n` and
`H(d, a, c) f = [a + (a + c d) x] f + c x (1 - x) f'`, the decomposition
evolves as `A_(n+1) = U_n A_n + epsilon_n x B_n`,
`B_(n+1) = delta_n A_n + V_n B_n` with `U_n = H(d_n, a_n, c_n)` and
`V_n = H(d_n - 1, a_n + delta_n, c_n)`.  Eliminating
`B_n = R_n / (epsilon_n x)`, `R_n = A_(n+1) - U_n A_n`, gives the operator

```text
E = eps_n S^2 - [eps_n U_(n+1) + eps_(n+1) Vt_n] S
    + [eps_(n+1) Vt_n U_n - eps_n eps_(n+1) delta_n x],   Vt_n = V_n - c_n (1 - x).
```

**Search.**  `a_n`, `c_n`, `epsilon_n` are exact rational polynomials of
degree at most `max_parameter_degree_in_n` (0 or 1), tried in increasing
degree on `fit_rows` rows.  `R_n(0) = 0` fixes `a_n` by exact interpolation of
`A_(n+1)(0) / A_n(0)` (an unknown when every constant term vanishes; ratios
that fit no polynomial of the degree are an obstruction).  Substituting `B` in
the second equation gives `G_n = eps_n R_(n+1) - eps_(n+1) Vt_n R_n -
eps_n eps_(n+1) delta_n x A_n = 0`, of degree at most three in the unknown
coefficients with an `epsilon` factor in every monomial.  Its `x`-coefficients
form a homogeneous linear system in those monomials (at most
`max_linearized_monomials` columns).  A one-dimensional exact kernel gives the
point directly (scale from `t k[e]^2 = k[e^2]`, coordinates `k[v e] / k[e]`,
full consistency check).  Otherwise the independent equations, together with
`z epsilon(n0) = 1`, are solved exactly: lex Gröbner basis (Buchberger), the
rational roots of its univariate element (`polytool::factorization`),
substitution and recursion, within `max_groebner_steps`.  A positive
dimensional system or a budget stop is reported as `underdetermined` or
`budget_exceeded` (`exhausted_budget`).  The search finds rational solutions
only, and it stops after the first degree with a certified lift.

**Certificate.**  `verify_symmetric_companion_certificate` replays a recorded
certificate without the search.  It recomputes `E`, the companion seed and the
parent from the parameters, then establishes for every `n >= n0`:

1. `(E A)_n = 0`: an exact Ore-algebra identity `g(n) E = sum_nu H_nu S^nu M`
   (right division when `M` has a constant leader, `g = 1`; for an order-two
   source with an `x`-free leader `f`, `f(n + 2) E = eps_n M`).  `M` holds at
   every base whose leading row the source rule generates; the remaining
   bases and the roots of `g` are checked on exact rows.
2. `eps_n != 0` on the domain (the exact index sign checker); a vanishing
   coupling is `unsupported`, never cancelled.
3. `R_n(0) = 0`: the derivative terms of `E` vanish at `x = 0` and the
   projection factors as `(eps_n S - g_n) o (S - a_n)`, so
   `w_n = A_(n+1)(0) - a_n A_n(0)` satisfies `eps_n w_(n+1) = g_n w_n`; with
   `w_(n0) = 0` it vanishes identically.
4. The second equation: `(E A)_n = eps_n E_(n+1) + eps_n eps_(n+1) x F_n -
   eps_(n+1) Vt_n E_n` (with `E_n`, `F_n` the residuals of the two equations)
   is checked as two operator identities (the `B` component is
   `Vt_n x = x V_n`); with `E_n = 0` by definition, `F_n = 0`.
5. Symmetry by induction on pairs from the two seeds, since `H(d, a, c)` maps
   symmetric-at-`d` to symmetric-at-`(d + 1)` and `x` adds two; `d_(n0) >= 1`.
6. The collapse identities `L_n = U_n + delta_n x`,
   `L_n x = eps_n x + x V_n` give the parent rule.

A row replay of `self_check_rows` indices is a defensive check only.

**Real-rootedness.**  For `d >= 1`, `r >= 0` and `0 <= eta <= c`, the step
`(a, b) -> (H(d, r, c) a + (c - eta) x b, eta a + H(d - 1, r + eta, c) b)` maps
pairs with nonnegative coefficients, `deg a <= d`, `deg b <= d - 1` and
`b <= a` to such pairs (partial fractions of `b / a` give `R' <= 0` and
`R + x R' >= 0`, which control the sign of `b+` at the roots of `a+`; boundary
cases by limits), and `a + x b` is then real-rooted (Brändén--Solus,
arXiv:1808.04141, Theorem 2.6).  The route applies when the identity is
certified, `a_n >= 0` and `0 <= delta_n <= c_n` for all integers `n >= n0`
(decided exactly by the index sign checker, below), `d_(n0) >= 1`, and `A_(n0)`, `B_(n0)` are
nonnegative with `B_(n0) <= A_(n0)`.  Here `f <= g` is directed weak
interlacing with `g` holding the rightmost root, which for nonzero rows is
`check_weak_interlacing(f, g) == Some(true)`; a zero member is accepted when
the other is real-rooted.  Real-rootedness of the parent alone is never used
for `A_n`.

Worked fixture (Ma--Ma--Yeh--Yeh, arXiv:2102.00899): the symmetric parts of
the `1/k`-Eulerian polynomials, given only their eliminated recurrence and
`A_2 = 1 + x`, `A_3 = 1 + (3k + 1) x + x^2`, give `a_n = 1`, `b_n = k n`,
`c_n = k`, `delta_n = k - 1`, `epsilon_n = 1`, `B_2 = k - 1`, `P_2 = 1 + k x`
for `k = 1` (zero companion), `2`, `3/2`, and the theorem applies.  With
`A_3 = 1 + 8x + x^2` and `k = 2` the same rule with `B_2 = 2` is a separately
valid certificate.

Limitations: parameters of degree at most one in `n`; only rational
solutions; a derivation of `E` from orders other than one (with a constant
leader) or two is not attempted, which leaves the identity a `candidate`; no
CLI or MCP adapter and no integration into `search_representations` yet.

## Exact sign decisions for an integer index

`decide_index_rational_function_sign(p, q, range, condition)` decides
`p(n) / q(n) R 0` for every integer `n` in `[first, last]` or `[first, +inf)`,
with `R` among `> 0`, `>= 0`, `< 0`, `<= 0`, `!= 0`, `= 0`, without floating
point.

1. **Taylor shift.**  If the coefficients of `p(first + u)` have one sign, and
   those of `q(first + u)` too with a nonzero constant term, the sign of `f` is
   known at `u = 0` and for every real `u > 0`; when the condition holds for
   both, the shifted coefficients are the certificate.
2. **Sign samples.**  Otherwise each nonzero polynomial gets a Sturm chain
   (`sturm::SturmChain`, on a positive integer multiple).  With `U` the
   ceiling of the Cauchy bound, every real root lies in `(-U, U)`.  The part
   `[max(first, -U), min(last, U)]` of the range is bisected at integer
   midpoints: a window `(l, h)` with no root (one Sturm count) is not split,
   and a window of width one with roots is recorded as a root cell.  The
   samples are `first`, `last`, every window end, and `l + 1` for a root-free
   window both of whose ends are roots.  So between consecutive samples
   `s < t` with `t - s >= 2` the polynomial has no root and does not vanish at
   both ends, and beyond the last sample of an unbounded range it has no root
   and does not vanish at the sample.  Hence its sign at every integer of the
   range is the sign at a sample.  The samples of `p` and `q` are merged (a
   refinement keeps both properties), `p` and `q` are evaluated exactly at
   each, and the smallest failing index is read off (a vanishing denominator
   fails every condition).

Each polynomial needs `O(deg * log U)` Sturm counts, independent of the length
of the range, so the decision is complete; the budget (degree, bisection
steps) is the only source of `undecided`.  The certificate replay recomputes
the signs, recounts the roots in every gap of length at least two and beyond
the last sample (up to the Cauchy bound), checks the root cells and the
condition, and rechecks the Taylor-shift data.

## Common lift summary

Every lift report carries a `LiftSummary`: operation, parent, identity status,
root transfer (theorem, parent and input root regions), hypotheses (sign
hypotheses with the decision above), parent proof and verdict.  The verdict is
computed by one function, `lift_verdict`: `route_failed` if the identity is
refuted, the parent fails or a hypothesis fails; `proved` if the identity is
certified, the parent is proved and every hypothesis is verified;
`conditional_on_attestation` with an attested parent instead; `open`
otherwise.  The mappings preserve each operation's own verdict: a Veronese
v1 lift is proved exactly for `proved_product_parent` (its identity is
certified by the identifying rows), a v2 certification exactly for
`VeroneseSectionVerdict::Proved`, and a symmetric-companion lift exactly when
the preservation theorem applies (its parent proof is the preservation
induction).  A `route_failed` or `open` summary says nothing against
real-rootedness of the input.
