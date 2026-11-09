# Recurrence search

For a complete generation/search/replay pipeline, see [CLI workflows](workflows.md).
`recurrence-generate`, `recurrence-degree`, `recurrence-prove`, `recurrence-theorems` and
`recurrence-veronese-lifts` share a bounded JSON reader:
they accept plain recurrence records, theorem reports with `recurrence_json`,
and recurrence representation menus. Select a zero-based `representations`
entry with `--candidate N`; a menu with multiple entries requires this option.
`--recurrence FILE` names the record; without it (or with `--recurrence -`)
the record is the command input, standard input or the file named by the
universal `--input FILE` option. Extraction does not verify a report's claims. `recurrence-generate --format rows`
names the default coefficient-list output explicitly.

Part of the [CLI reference](cli.md).  See also
[recurrence transformations](cli-recurrence-transforms.md) and
[searching equivalent representations](cli-recurrence-representations.md).

## Search for a recurrence

Given a sequence of polynomials (one per line), search for a polynomial
recurrence `f(n,t) P_n(t) = Σ c_{r,d}(n,t) D^d P_{n-r}(t)`. With
`--try-alternating-sign`, the search also allows terms
`(-1)^n c_{r,d}(n,t) D^d P_{n-r}(t)`:

```sh
polytool recurrence < polys.txt
```

The adaptive search orders candidates by a weighted parameter count. Ordinary
coefficient unknowns are cheapest; derivative, denominator, alternating-sign,
and inhomogeneous unknowns are delayed. A candidate is only solved when the
available fitting prefix has at least `unknowns + min_margin` scalar equations.
By default the last input row is reserved for exact verification, the fit uses
the first solvable prefix plus one extra row, and modular prefiltering rejects
impossible candidates before exact rational solving. Use `--no-held-out-rows` to fit
against all input rows.

Options:

```
--skip-prefix <k>    Ignore the first k input polynomials before searching
--first-index <n>    Index of the first input row (nonnegative integer)
--min-recurrence-length <k>    Minimum recurrence depth to try (default: 1)
--max-recurrence-length <k>    Maximum recurrence depth (default: 5)
--min-variable-degree <d>    Minimum degree in t for coefficients (default: 0)
--max-variable-degree <d>    Maximum degree in t for coefficients (default: 3)
--min-index-degree <d>    Minimum degree in n for coefficients (default: 0)
--max-index-degree <d>    Maximum degree in n for coefficients (default: 3)
--min-derivative-order <d>   Minimum derivative order (default: 0)
--max-derivative-order <d>   Maximum derivative order (default: 2)
--try-inhomogeneous      Also try inhomogeneous recurrences
--min-inhomogeneous-variable-degree Minimum degree in t for the inhomogeneous term
--max-inhomogeneous-variable-degree Maximum degree in t for the inhomogeneous term
--min-inhomogeneous-index-degree Minimum degree in n for the inhomogeneous term
--max-inhomogeneous-index-degree Maximum degree in n for the inhomogeneous term
--try-denominator        Allow a nontrivial LHS factor f(n,t)
--try-alternating-sign   Also allow right-hand-side terms multiplied by (-1)^n
--no-laurent         Skip the left-hand side variation stage (x^e factors and
                      normalization charts; on by default in this command)
--max-laurent-power <e>
                      Largest left-hand power x^e (default: 2)
--max-denominator-variable-degree  Max degree in t for f(n,t) (default: 2, implies --try-denominator)
--max-denominator-index-degree  Max degree in n for f(n,t) (default: 2, implies --try-denominator)
--min-equation-margin <k>     Require equations >= unknowns + k (default: 1)
--fit-extra-rows <k> Extra rows after the first solvable prefix (default: 1)
--no-held-out-rows          Fit all rows instead of reserving held-out verification rows
--no-modular-prefilter
                      Disable default modular candidate rejection
--max-candidates <n> Inspect at most n adaptive candidate configurations
--max-work <n>       Stop after n deterministic work units across the whole search
--json               Emit recurrence JSON with initial conditions
--latex              Emit LaTeX with subscripts, initial conditions and starting index
--python             Emit exact standalone Python code using Fraction arithmetic
--format json        Alias for --json
--format latex       Alias for --latex
--format python      Alias for --python
--verbose            Print each candidate tried
```

The modular prefilter is enabled by default. It is often much faster on false
candidates because it rejects a candidate when every usable fixed large-prime
reduction is inconsistent. Use `--no-modular-prefilter` only when comparing
against the exact-only search path.

`--max-candidates N` is a deterministic iteration budget, not a wall-clock
timeout. A candidate is counted as soon as its parameter configuration is
taken from the adaptive iterator, before fit-row, structural, modular, or
exact-solve filtering. Thus `0` evaluates no candidates, candidate `N` may
still succeed, and budget exhaustion is reported only when an `(N+1)`st
candidate exists. In text mode exhaustion is a distinct diagnostic and exit
status `3`; with `--format json`, the JSON status is `budget_exhausted` and
includes `candidates_considered`, `candidates_tried`, and `max_candidates`.
If all configured candidates have been checked, the JSON status is instead
`not_found`. Omitting the option preserves the unbounded behavior.

`--max-work N` bounds what happens *inside* the candidates. Work is counted in
abstract units computed from the dimensions of the systems that are solved
(a modular or exact solve of `r` equations in `c` unknowns costs
`r * c * (min(r, c) + 1)`; held-out and lift checks, derivative and modular
precomputation, and each candidate and its structural prechecks are charged
too), never from the clock, so the same input, options and `N` stop at the same
step on any machine. A stage is charged before it runs and runs only if its
whole cost fits, so the units used never exceed `N`; the first stage that does
not fit ends the search as `budget_exhausted` (exit status `3`). The JSON status
carries `work_units` and `max_work`, and a found record has `search.work_units`.
The option combines with `--max-candidates` (whichever limit is hit first
ends the search); it is not supported with `--representations`, which has its
own `--max-candidates` bound. A pipeline that fixes one `N` is both reproducible
and bounded in time.

By default, retained rows are numbered from 1, including after
`--skip-prefix`. With explicit `--first-index N`, the original first input row
is `P_N`, so `--skip-prefix K` makes the first retained row `P_(N+K)`.
Every coefficient, denominator, forcing term and alternating sign is shifted
to that numbering. Text, Python, JSON, representation LaTeX and theorem-match
reports use the same indices. Indices must be nonnegative and fit the supported
integer range.

Successful JSON search metadata includes `search.verification_equations`:
the number of dense coefficient slots checked on held-out rows. For each row
we count from degree zero through the largest degree of an evaluated summand
on either side, before cancellation between summands. Derivatives, denominators
and forcing terms contribute their actual degrees; a zero identity contributes
one slot. This counts equations, including zero or dependent ones, **not**
independent constraints or a probability of correctness. It is zero with
`--no-held-out-rows`; older JSON records without the field remain readable (default 0).
Representation entries expose the same count directly as well as in their
nested recurrence JSON.

### Several recurrence representations

Use `--representations` to compare lag length against derivative order instead
of stopping at the first match:

```sh
polytool recurrence --representations --json \
  --max-recurrence-length 3 --max-derivative-order 2 --max-variable-degree 3 --max-index-degree 1 \
  --try-denominator --max-denominator-variable-degree 0 --max-denominator-index-degree 1 \
  --max-candidates 2000 < rows.txt > representations.json
```

Without `--json`, the command prints an aligned comparison table followed by
numbered formulas and their checked index ranges. `lag` and `deriv` are actual
orders; `deg(n,t)` and `lhs(n,t)` are the coefficient/forcing and left-factor
degrees. `coeffs` counts nonzero scalar coefficients, `fit`/`check` count fitting
and held-out rows, and `eqs` counts held-out coefficient equations, including
dependent ones. A `*` marks the frontier. Table/formula indices match the JSON
menu, so the same `--candidate N` selects the displayed formula. Checked
singular indices, when present, are listed next to its formula.

The additive Rust APIs are
`recurrence::find_recurrence_representations_rational` and the `i64` convenience
wrapper `find_recurrence_representations`. Both take rows, existing
`AdaptiveSearchOptions`, and an `AdaptiveSearchBudget`, and return
`Result<RecurrenceRepresentations, RecurrenceError>`. The MCP tool
`find_recurrence_representations` accepts the same inputs and controls as
`find_recurrence`. Existing single-result APIs and schemas are unchanged.

The `polytool.recurrence-representations.v1` report contains:

- `representations`: one preferred candidate per length/order bound pair,
  deduplicated when formula and checked domain agree. Each includes text,
  LaTeX, a reusable `recurrence_json` **object**, actual nonzero coefficient
  complexity, fit/held-out row counts, and checked starting indices.
- `frontier`: indices into that array. These candidates are nondominated in
  actual lag, derivative order, coefficient and denominator degrees, use of
  forcing/alternating terms, number of nonzero coefficients, and checked start.
  The frontier is relative to returned candidates, not a minimality theorem.
- `coverage`: requested length/order pairs, candidate indices where found,
  configuration counts, and `found`, `not_found`, `incomplete`, or
  `insufficient_rows` statuses.
- `search`, `source_rows`, and global budget/coverage counters.
  `status` is `complete` or `budget_exhausted`; CLI exhaustion exits with code 3
  **after** emitting the partial report.

Within each pair the existing weighted-complexity ordering chooses the first
successful fit. The scalar path shares derivative/modular caches, and solved pairs are
skipped without spending more of the global configuration budget.
Denominators are strictly opt-in in this API. Its denominator fallback reuses
the vector fitter's projective normalization (with a per-shape system), so denominators such as `n` or
`t`, whose constant coefficient is zero, are supported too.

Every result has `evidence: "exact_fit"`: finite rows, even held-out rows, do
not prove an identity for all indices. `complete` means this representative
search finished, **not** that all recurrences have been classified.
`not_found` is not a nonexistence certificate: fitting-prefix selection,
equation-count safeguards, and choice of a particular solution can miss
representations. On budget exhaustion, unresolved cells are conservatively
labelled `incomplete`. This API does not return a complete coefficient kernel
or an Ore-ideal basis.

Rows after `skip_prefix` are renumbered from 1 unless `--first-index` is
explicitly supplied, as described above. `first_checked_index` is the
first checked target in the output numbering; `first_checked_input_row` refers to
the original input. JSON retains **every** initial row before the checked
domain, even when the actual lag is smaller. We do not cancel polynomial
common factors or extend a checked domain implicitly. Vanishing denominators
on supplied target rows are listed in `singular_checked_indices`; a formula
with such an index cannot generate that row. Future singularities are not
classified, and the row generator reports them if encountered.

For a checked example using the canonical OEIS generator:

```sh
cargo run -p polytool --example recurrence_representations
```

For A171608 with the displayed affine index bounds and index-only denominator,
the menu includes length/order `(3,0)`, `(2,1)`, and `(1,2)`. Tests fit 40 rows,
then independently replay each exported recurrence through row 60. This is
finite exact verification, not an all-rank certificate.

To generate rows or check conditional theorem matches for one entry, select
its zero-based index directly:

```sh
polytool recurrence-generate --recurrence representations.json --candidate 0 --count 60
polytool recurrence-theorems --recurrence representations.json --candidate 0 --json
```

`--representations --theorem-matches` is rejected with this follow-up command;
the flag is not silently ignored. The separate theorem command checks the
selected recurrence and its exported initial rows. A conditional theorem match
does not establish the recurrence for all indices of the intended sequence.

### Complete recurrence spaces and certificates

Use `--space` for the full rational coefficient space in specified supports.
The existing single-result finder and `--representations` retain their finite-fit
semantics. The new mode reads a structured request rather than polynomial text:

```sh
polytool recurrence --space < request.json > spaces.json
polytool recurrence --verify-space < spaces.json
cargo run -p polytool --example recurrence_spaces > a171608-spaces.json
```

A small request for `P[n] = X^n` is:

```json
{
  "source": {"coefficients": [["0", "1"]], "seeds": [["1"]]},
  "domain_start": 0,
  "shapes": [
    {"kind": "box", "length": 1, "derivative_order": 1,
     "index_degree": 1, "variable_degree": 1, "leading_index_degree": 1}
  ]
}
```

`source.coefficients[i]` is `a[i](X)` in the defining, homogeneous monic rule
`P[n+k] = sum(i<k) a[i](X) P[n+i]`; exactly `k` seeds are required.
Coefficients are ascending-power rational **strings**. `source.origin` defaults
to zero. The source is valid from that origin and its coefficients must be
independent of `n`. Alternatively, `source_recurrence` accepts an existing
`polytool.recurrence.v1` object satisfying these restrictions. Extra initial
rows and optional `rows` are checked for consistency. A guessed source is an
assumption, not proof that an external sequence satisfies it.

Omit the source to analyze finite `rows`, with `row_origin` defaulting to zero.
`domain_start` is the first base index tested, must not precede the origin,
and defaults to zero; set it explicitly for a nonzero origin. All coefficient
polynomials use the **original** index, including on a tail.

A box includes every RHS monomial `n^u X^v D^s(P[n+j])` for
`j < length`, `s <= derivative_order`, `u <= index_degree`, and
`v <= variable_degree`. Its leader multiplies `P[n+length]` and has separate
`leading_index_degree` and `leading_variable_degree` bounds. All omitted degree
bounds default to zero. For sparse supports use
`{"kind":"explicit","target_lag":1,"monomials":[{"lag":1},{"lag":0,"variable_degree":1}]}`.
Each monomial has `lag`, `derivative_order`, `index_degree`, `variable_degree`;
omitted powers are zero. Target-lag derivatives are forbidden. Columns are
sorted and the report records that ordering. Target columns enter the residual
positively, other columns negatively.

The schema `polytool.recurrence-spaces.v1` echoes the request and returns:

- A **complete kernel basis** per completed shape, including relations with zero
  successor coefficient; its `basis` entries are coefficient vectors.
- `successor_projection_rank` and representatives with independent projections
  onto the leader coordinates. Their linear combinations together with the
  zero-leader relations describe every recurrence in that shape.
- Exact rank witnesses (an inverse pivot minor and an identity minor of the
  kernel basis), plus row-space dual witnesses for `absent_in_shape`.
- Evidence `conditional_on_source` for all-rank identities, `exact_fit` for
  finite-row positives, or `rows_only_obstruction` for bounded nonexistence.
  Budget or row shortages carry `evidence: "none"` and never mean absence.
- A frontier among compatible **requested boxes**, holding coefficient bounds
  fixed. It is not an unrestricted minimal-order claim or an Ore ideal basis.

Why finitely many rows suffice: let `chi(S)` be the monic order-`k` source
operator. For support exponent pair `(u,s)`, set `h = 1 + max(u+s)`, including
leader columns. Commutators with `n` and differentiation show that `chi(S)^h`
annihilates every residual in this support, without a simple-root assumption.
Thus zero residuals at `domain_start` through `domain_start+h*k-1` are necessary
and sufficient for an all-rank identity on that domain. We generate rows through
`domain_start+h*k-1+target_lag`. Original source order, including a nonminimal
source or zero oldest coefficient, is safe. See
[RealRooted #1978](https://github.com/PerAlexandersson/RealRooted/issues/1978)
for the formalization proposal; this implementation does not emit a Lean proof.
For an E-by-M residual matrix, elimination uses O(E M min(E,M) + M^3)
rational operations and O(E M + M^2) storage; integer bit growth is separately
bounded. The row-generation cost depends on source order and polynomial degrees.

The Rust entry points are `find_recurrence_spaces(&RecurrenceSpaceRequest)` and
`verify_recurrence_spaces(&expected_request, &report)`. Verification regenerates
residuals and checks the supplied kernel/rank/dual witnesses, reducing only the
small leader projection. It binds to the caller's expected request. CLI replay
checks against the embedded request; it does not authenticate the source.
`ConstantLagSource::generate_rows` is the canonical exact source generator.
`representative.to_recurrence(result)` converts original forward indices to the
existing backward evaluator, but callers must provide the required seeds.

`singular_ranks` lists **all nonnegative integer base ranks** at which the whole
leader vanishes as a polynomial in X. Exact Sturm analysis is limited to degree
at most 8 and coefficient bit size at most 128; outside these diagnostic bounds
it is `null`, meaning unknown. Such ranks need extra seeds for generation.
`requires_x_division` separately records a variable-dependent leader. We never
cancel boundary factors or silently extend the declared domain.

Default limits per shape/phase are 512 generated rows, 128 columns, degree 512,
500,000 matrix entries, 20,000,000 arithmetic-work units, and 4096-bit
coefficients, with at most 16 shapes. `limits` can override individual fields;
MCP and WASM allow only reductions and cap JSON input at 1 MiB. Generation,
elimination, and witness replay have separate work accounting, so replay can independently reject a report on resource grounds. Such a
rejection is not a mathematical counterexample.
CLI exit codes: 0 for completed claims, 3 if any shape is incomplete, and 2 for
invalid input or a rejected certificate. `--verify-space` reports the number of
incomplete shapes and only verifies completed claims. These two flags do not
combine with legacy search controls; all controls are in the JSON request.

### Browser verification

The recurrence browser regression drives the actual page, Web Worker and WASM
at desktop and 375-pixel widths. It checks guided full-search controls, exact
row parsing and origins, example formulas, clipboard copying, backend replay,
partial-result labels, cancellation, tab routing, stale-input invalidation,
renderer-error recovery and the older finder modes. In the Transform tab it
checks the examples, the evidence-label styles, clipboard and download export,
theorem and root sub-reports, replay of a computed report and of imported
original and tampered reports, implicit and vector results, stale-setting
invalidation, cancellation, and the hand-off from the Recurrence tab. It also
drives the representation-search mode: checklist and bounds, the tree, budget
and equal-sequence markers, duplicate and rejected maps, replay, a tampered or
relabelled import, and the phone layout.
It uses a temporary localhost server (serving `index.html`, `app.css`,
`js/*.js` and `favicon.svg` from `web/`, and `pkg/` from `POLYTOOL_WASM_DIR`)
and closes its own browser and server.
No hosting deployment is involved.

With Playwright/Chromium installed and a no-modules WASM build available:

```sh
NODE_PATH="$(npm root -g)" \
POLYTOOL_WASM_DIR=/external/path/to/wasm-bundle \
  node polytool/web/tests/recurrence_browser.cjs
```

Run from the parent workspace, or omit the initial `polytool/` in a standalone
checkout. `POLYTOOL_PLAYWRIGHT` can name the installed Playwright module and
`POLYTOOL_CHROMIUM` an existing Chromium executable. Optional
`POLYTOOL_BROWSER_ARTIFACTS` writes screenshots outside the source checkout.
The ordinary `node --test polytool/web/tests/*.test.cjs` checks remain useful, but
the browser test additionally checks real DOM events and visible rendering.
A viewport test does not substitute for real-device or screen-reader testing.

### Automatic variable degree, collective leaders, and indexed sources

These extensions use the same `--space` / `--verify-space` entry points and
`polytool.recurrence-spaces.v1` report. Existing requests do not enable the new
leader-combination output. The new fields are optional on replay.

To remove a guessed RHS degree in X, use this shape:

```json
{"kind":"automatic_variable_degree","length":1,"derivative_order":1,
 "index_degree":1,"leading_index_degree":1,"leading_variable_degree":0}
```

A defining source is required. Index-degree, derivative-order and leader bounds
remain prescribed. The algorithm computes a sufficient RHS variable degree
from the finite polynomial residual matrix. `automatic_degree_bound` records
the bound and its inputs. The **leader projection** then covers every RHS
X-degree; the returned kernel is the finite space at the computed cap, not the
whole polynomial module. A completed negative result therefore excludes all
RHS X-degrees under those remaining bounds. A large computed cap can exceed
resource limits and yield `incomplete_budget` with the computed degree in the
reason; this never becomes a negative mathematical result.

The bound is elementary. If `H(X)p=b(X)` has a polynomial solution, H has
rank rho, its entry degrees are at most D, and deg b <= B, then a solution
exists with degree at most `max(B+(rho-1)D, rho*D-1, 0)`. Divide nonpivot
coordinates by a nonzero rank-minor determinant and transfer their quotients
to pivot coordinates using its adjugate; Cramer's rule bounds what remains.
We use the safe upper bound `rho <= min(number of residuals, RHS columns)`.
The zero RHS-matrix case is handled separately. This preserves every allowed
leader, not only one positive witness.

Set `"include_leader_coverage": true` to request a combined recurrence. If the
independent projected representatives are `L_i` with leaders `delta_i`, the
report returns `sum_i delta_i L_i`, with leader `sum_i delta_i^2`. Its leader
vanishes identically in X at a real rank iff every original leader does. Lag
and derivative order are preserved, while coefficient degrees can increase.
For index-only leaders this coefficient is positive at every regular real
rank. The output `leader_coverage` contains the polynomial multipliers,
combined support/formula, and common singular ranks (or null if unclassified).
These ranks are unavoidable within the polynomial module generated by the
reported identities, not necessarily among arbitrary larger supports.
The combined recurrence has the same evidence and domain as its input kernel;
finite fits do not become infinite identities. Optional combination output can
be omitted if its own resource budget is exceeded. Rust callers can convert
it with `coverage.to_recurrence(result.target_lag)` for exact row generation.

For a monic source polynomial in n, use `indexed_source`, whose coefficient
arrays are `[lag][n_power][X_power]`, with ascending powers and rational strings:

```json
{
  "indexed_source": {
    "coefficients": [[["-2"], ["-2"]], [["0", "2"]]],
    "seeds": [["1"], ["0", "2"]]
  },
  "include_leader_coverage": true,
  "shapes": [{"kind":"box","length":1,"derivative_order":1,"variable_degree":1}]
}
```

This defines the physicists' Hermite polynomials by
`H[n+2] = 2X H[n+1] - 2(n+1) H[n]`; the result certifies
`H[n+1] = 2X H[n] - H'[n]`. The source holds from its `origin` (default zero).
Exactly k seeds are required; optional consistency rows and original-index
semantics are unchanged. `source`, `source_recurrence` and `indexed_source`
are mutually exclusive. The existing `source_recurrence` adapter still accepts
only constant-index sources; use `indexed_source` for the new coefficient form.
Nonmonic/rational-index source denominators are not yet supported here.

`IndexedLagSource::generate_rows` is the canonical generator. It is tested
against exact Hermite and factorial-times-monomial rows, and returned Hermite
lifts replay in the older independent row evaluator. Checked examples:

```sh
cargo run -p polytool --example recurrence_spaces -- hermite > hermite-space.json
cargo run -p polytool --example recurrence_spaces -- automatic > automatic-space.json
polytool recurrence --verify-space < hermite-space.json
```

For maximum derivative order s, the k(s+1)-dimensional derivative state has
a polynomial transition matrix T(n,X). Each basic residual column is a row
of V(n,X) times that state; X-powers are omitted because they commute with
shift. We find a polynomial dependence among the matrices
`W_0=V`, `W_(j+1)=W_j(n+1)T(n)`. With m distinct basic columns its order is
at most `m*k*(s+1)`. The emitted `propagation.coefficients` annihilate every
residual in the support. The verifier checks that polynomial identity directly,
without rerunning the dependence search.

`propagation.regular_from` is beyond every nonnegative integer zero of the
whole leading polynomial. `cutoff_method` is `exact_integer_zeros` when the
bounded exact diagnostic applies, or `coefficient_root_bound` for a safe
Cauchy bound on a nonzero X-coefficient; in the latter case `exceptional_ranks`
is null. The sufficient prefix ends at `regular_from + order - 1`.
Thus a late exceptional rank increases the checked prefix; the constant-source
h*k shortcut is not applied to indexed sources. The result's `h` still records
the support filtration, while `residual_count` comes from the propagator.
This uniform propagator also justifies `automatic_variable_degree` for indexed
sources. Zero-order propagation certificates are supported.

The initial propagation implementation uses exact polynomial elimination,
with at most cubic-in-`m*k*(s+1)` polynomial row work plus transition products.
Degrees and coefficient sizes can grow substantially; entry/degree/bit/work
limits apply before successful certification. This is not a claim of polynomial
bit complexity or a shortest propagator. No source hypothesis or certificate
is promoted to a trusted Lean axiom.

### Complete spaces after factorial coefficient rescaling

The optional request field `coefficient_rescaling` searches the rescaled rows
`Q[n] = T_m(P[n])`, where `T_m(X^k) = X^k/(k!)^m` divides coefficient `k` by
`(k!)^m` for a fixed integer `m >= 1`. It uses the same `--space` /
`--verify-space` entry points, MCP tools and `polytool.recurrence-spaces.v1`
report; requests without the field behave exactly as before.

```json
{
  "source": {"coefficients": [["1", "1"]], "seeds": [["1"]], "origin": 0},
  "coefficient_rescaling": {"kind": "factorial_power", "power": 1},
  "shapes": [
    {"kind": "box", "length": 1, "derivative_order": 1,
     "index_degree": 1, "variable_degree": 1, "leading_index_degree": 1},
    {"kind": "box", "length": 2, "index_degree": 1, "variable_degree": 1,
     "leading_index_degree": 1}
  ]
}
```

For canonical A007318, `P[n] = (1+X)^n` from origin 0, this gives
`Q[n] = sum_k binom(n,k) X^k/k!` and certifies one-dimensional complete
kernels spanned by

```text
(1 + n) Q[n+1] = (1 + n + X) Q[n+0] + (X) D^1(Q[n+0])                 n >= 0
(2 + n) Q[n+2] = (-1 - n) Q[n+0] + (3 + 2n + X) Q[n+1]                n >= 0
```

The second is `(n+1)Q[n+1] = (2n+1+X)Q[n] - nQ[n-1]` for `n >= 1` in the
original forward-index convention. Formulas name the rescaled rows `Q`.
Each result adds `rescaling_bound`: the rescaling, source order `r`, the
support maxima `U` (index degree), `A` (variable degree) and `J` (derivative
order) over all columns including the leader, the exponent
`H = U + m(A+J) + 1` (also reported as `h`), the rescaled initial rows
`Q[domain_start + j]` for `j < target_lag` (here `Q[0] = 1`, respectively
`Q[0] = 1, Q[1] = 1 + X`), and a completeness statement.

Why `r H` residuals suffice: multiply a residual by `X^J` and pull it back by
`T_m^-1`. With `theta = X D`, a term `n^u X^a D^j S^i` becomes
`n^u X^s theta^(j falling) (prod_(h=1..s) (theta+h))^m S^i` with
`s = a + J - j >= 0`, whose derivative order is at most `m(A+J)`. As for the
untransformed bound, `chi^(u+d+1)` annihilates `n^u X^v D^d P[n+i]`, so
`chi^H` annihilates the pulled-back residual; this monic recurrence has order
`r H`. Multiplication by `X^J` and `T_m^-1` are injective, so the first `r H`
zero residuals from `domain_start` imply that all later residuals vanish. The
returned kernel is therefore the complete space of all-index identities for
`Q` in the requested support, conditional on the supplied source. This is a
sufficient uniform count; sharper counts would need separate proofs.

Validated hypotheses: `m >= 1` (and at most `limits.max_degree`), a
constant-index monic `source` or `source_recurrence`, and `domain_start` not
before the source origin. Rows-only requests, `indexed_source` and
`automatic_variable_degree` shapes are rejected with this field. Optional
consistency `rows` remain rows of the source `P`. Row generation budgets
include the target-lag shift, and a rescaled coefficient above
`max_coefficient_bits` yields an explicit `incomplete_budget` result, never an
absence claim. Verification regenerates the source rows and replays the
rescaling, support, bound and certificates against the expected request.
Rust callers can apply the exact row map and its inverse with
`RecurrenceCoefficientRescaling::{rescale_coefficients, unscale_coefficients}`.

```sh
cargo run -p polytool --example recurrence_spaces -- rescaled > rescaled-space.json
polytool recurrence --verify-space < rescaled-space.json
```

### Machine-readable recurrence JSON and row generation

To save a recurrence with enough initial conditions to regenerate the sequence,
use JSON output:

```sh
polytool recurrence --json < polys.txt > recurrence.json
```

The JSON schema stores exact rational coefficients as strings, the recurrence
terms, the LHS denominator if present, and enough initial polynomial rows
to reach the checked recurrence domain.
Generate more rows from the saved recurrence with:

```sh
polytool recurrence-generate --recurrence recurrence.json --count 100 > rows-100.txt
```

Use `--additional-count n` instead of `--count n` to keep all initial rows and append
`n` newly generated rows. By default generation preserves the record's
`first_index`. To relabel it, use:

```sh
polytool recurrence-generate --recurrence recurrence.json --first-index 0 --count 100
```

This shifts the rule as well as its labels: it generates exactly the same rows,
with a different starting index. It does not skip rows or evaluate the original
rule at a different starting point. Rust callers can use
`RecurrenceJson::reindexed(first_index)` for the same operation.

### LaTeX presentation

```sh
polytool recurrence --first-index 0 --latex < polys.txt
```

The output contains one formula per line, with subscripts, exact fractions,
grouped coefficients and derivative notation. For example:

```latex
% Candidate fitted to finite rows; validity for all n is unproved.
$P_{0} = 1$
$P_{n} = (1 + t) P_{n-1} \quad (n \geq 1)$
```

A recurrence with index-dependent coefficients looks like this (seeds
`1` and `1 + t`):

```latex
$P_{0} = 1$
$P_{1} = 1 + t$
$(n + 2) P_{n} = (2n + 1)(1 + t) P_{n-1} - (n - 1)(1 - t)^{2} P_{n-2} \quad (n \geq 2)$
```

All supplied seed rows in the exported record are printed. The starting index
is immediately after them, including any extra seeds needed before the
recurrence can generate rows. `--first-index` and `--skip-prefix` use the same
conventions as JSON and Python export. The LaTeX comment records the finite-fit
status; the displayed continuation is a candidate recurrence for the input
sequence.

Every coefficient is factored exactly over $\mathbb{Q}$ in `n` and `t`
(`BivarPoly::irreducible_factorization`): polynomials in one variable and
polynomials monic in `n` go straight to the exact factorization routines, any
other coefficient is first sheared by `t = v + c n` (with `c` chosen so the top
homogeneous part is nonzero at `(1, c)`), which makes it monic up to a constant.
The factors are mapped back, written as primitive integer polynomials whose first
displayed term (ascending powers of `t`, then descending powers of `n`) is
positive, repeated factors carry exponents, and the rational constant and overall
sign are written in front.  The product is multiplied out and compared with the
coefficient before it is printed.  If the factorization budget (recombination
subsets) is exhausted, that coefficient stays expanded and a `%` comment line
names it; an unfinished factorization is never presented as complete.

`--json` includes the same text in the optional `latex` field. Generation
uses the exact recurrence data and ignores that presentation field. Older
records without it remain accepted. With `--theorem-matches`, the field is in
the nested `recurrence_json`; theorem notes go to stderr in LaTeX mode.

`--representations --latex` prints each definition with a numbered LaTeX
comment. Representation JSON retains its bare formula in `latex` and includes
the complete definition in `recurrence_json.latex`.

Rust callers can use `Recurrence::to_latex()` for a bare expanded formula,
`Recurrence::to_latex_factored()` for the factored one, or
`RecurrenceJson::to_latex()` (factored) and `RecurrenceJson::to_latex_expanded()`
for the complete presentation. The latter checks
seed availability and index overflow. `RecurrenceJson::with_latex()` refreshes
the optional stored presentation; `reindexed()` refreshes it when present.

### B-files from a recurrence

`recurrence-generate --format bfile` writes the generated rows as an OEIS b-file
(`index value` lines, coefficients in ascending powers within each row, so a
triangle is read by rows).  Every coefficient must be an integer.

```sh
polytool recurrence-generate --recurrence recurrence.json --format bfile \
  --first-term-index 1 --max-terms 1000 > b008292.txt
```

- `--first-term-index K` is the index of the first term (the OEIS offset; default 0).
- `--max-terms N` caps the number of terms and never splits a row: by default
  the output is the longest prefix of whole rows with at most `N` terms, and the
  command fails if even the first row is longer.  `--nearest-row` instead picks
  the whole row nearest to `N` terms, which may exceed it (a tie goes to the
  shorter prefix).  Without `--count` or `--additional-count`, rows are generated until
  the cap is reached.
- `--descending` writes each row from its highest coefficient down to its
  constant term.
- `--trim-zeros none|low|high|both` drops zero coefficients at the constant
  (low) end, the leading (high) end or both ends of each stored row, before
  `--descending` is applied.  The default keeps every coefficient, as
  `oeis generate --format bfile` does.

- `--from-index K` starts at the row with absolute index `K` of the record
  (default its first index); `--regular-triangle` reads the rows as a regular
  triangle (see below).

For a catalog entry, `polytool oeis generate A008292 --count 12 --format bfile
--max-terms N` writes the same lines with the catalog's own offset; both
commands use the library function `oeis::flatten_rows_to_bfile`.

### OEIS text from a recurrence

The formats `oeis-example`, `oeis-formula`, `mathematica`, `pari` and
`oeis-json` write fields for an OEIS submission from the rows of the record:

- `oeis-example`: the `%e` block `Triangle begins:`, one line per row with the
  row number on the left and right-aligned columns;
- `oeis-formula`: the `%F` sentence, the recurrence in plain ASCII with its
  initial rows, e.g. `The row polynomials satisfy P(n,x) = (2+x)*P(n-1,x) +
  (1+x)*P(n-2,x) for n >= 2, with P(0,x) = 1, P(1,x) = 1 + x.` (derivatives as
  `d/dx P(n-1,x)`, a left-hand factor as `(n+1)*P(n,x) = ...`);
- `mathematica`: a `%t` program, a memoized `P[n_, x_] := P[n, x] =
  Expand[...]` (a left-hand factor is divided out with `Together`) and
  `Table[CoefficientList[P[n, x], x], {n, o, N}] // Flatten`;
- `pari`: a `%o` PARI/GP program listing the same rows with `T(N)`;
- `oeis-json`: all four (schema `polytool.oeis-text.v1`).

The rows from `--from-index` on are read as OEIS rows: `--first-term-index` is the OEIS
row number of the first row (default: the first row keeps its own index),
`--descending` lists each row in descending powers, `--trim-zeros` drops zeros
at the ends of each ascending row, and `--regular-triangle` gives row `i`
(counted from the first) exactly `L0 + i` entries, padding a shorter row with
zeros at its high end and stopping before a row that is longer.  The index
`n` of the formula and programs is the OEIS row number (the recurrence is
shifted to it); `--count N` sets the number of rows shown (default 10).

```sh
# A156919: the reversed rows from n = 1 on, offset 0.
polytool recurrence-generate --recurrence rec.json --format oeis-formula \
  --from-index 1 --first-term-index 0 --descending --regular-triangle
```

The web interface exports this compact JSON record and links back to these
generation instructions instead of emitting a large standalone Python program.

Use Python output when a sequence extension script is more convenient than JSON:

```sh
polytool recurrence --python < polys.txt > recurrence.py
```

The generated Python uses exact `fractions.Fraction` arithmetic and dense
coefficient lists in ascending powers of `t`.

## Degree laws and leading coefficients

`polytool recurrence-degree --recurrence recurrence.json [--json]` decides the
degree law `deg P_n = D(n)` and the sign of the leading coefficient of every
row, for all `n`, exactly (library: `certify_leading_coefficient`, see
[Library usage](library.md#degree-laws-and-leading-coefficients-of-recurrence-rows)).
`D(n)` is periodic by residue class (affine, half growth, Veronese sections,
two-step products); it is read off exactly computed rows (`--count N`, default
36; `--max-period P`, default 6) and proved for every later index by one of
the RealRooted routes `cone`, `recurrence`, `subseq`, `twoStep`, `cancel` or
`coefficients`.  The text output lists the conclusion, every hypothesis with
its status, the attempts of a family that is not proved, and the RealRooted
declarations; `--json` emits `polytool.recurrence-leading-coefficient.v1`, and
`polytool recurrence-degree --verify < certificate.json` replays it.  Exit
status 0 proved (or a valid replay), 1 not proved (a failed route is not a
disproof), 2 invalid input or replay.

```sh
polytool recurrence --json < rows.txt | polytool recurrence-degree --recurrence - --json > degree.json
polytool recurrence-degree --verify < degree.json
```

## General proof discovery

`polytool recurrence-prove --recurrence recurrence.json --json` checks eight
general theorem families (the parent cascade, including the iterated-preserver
and hypergeometric-row criteria) and the RealRooted row routes fixed-delay
window (`variable_delay`) and residue classes (`residue_class_subsequences`)
where the rule has their form, then the tail after rows that are zero or not
real-rooted and certified alternative representations
(`--max-representation-candidates N`, 0 disables them), and emits a replayable
report with exact hypotheses, conclusions and unresolved obligations. A
residue-class proof may leave whole residue classes of zero rows: the report
then omits `real_rooted_from_index` and names them in
`real_rooted_outside_zero_residue_classes`. Use
`--candidate N` for a menu and
`polytool recurrence-prove --verify --json < proof.json` to replay. See
[General recurrence proof discovery](proof-discovery.md) for the supported
rules, evidence boundaries, product certificates and Lean correspondence.
