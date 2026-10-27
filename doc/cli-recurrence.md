# Recurrence search

Part of the [CLI reference](cli.md).  See also
[recurrence transformations](cli-recurrence-transforms.md) and
[searching equivalent representations](cli-recurrence-representations.md).

## Search for a recurrence

Given a sequence of polynomials (one per line), search for a polynomial
recurrence `f(n,t) P_n(t) = Σ c_{r,d}(n,t) D^d P_{n-r}(t)`. With
`--alternating-sign`, the search also allows terms
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
impossible candidates before exact rational solving. Use `--no-verify` to fit
against all input rows.

Options:

```
--skip-prefix <k>    Ignore the first k input polynomials before searching
--min-rec-len <k>    Minimum recurrence depth to try (default: 1)
--max-rec-len <k>    Maximum recurrence depth (default: 5)
--min-var-deg <d>    Minimum degree in t for coefficients (default: 0)
--max-var-deg <d>    Maximum degree in t for coefficients (default: 3)
--min-idx-deg <d>    Minimum degree in n for coefficients (default: 0)
--max-idx-deg <d>    Maximum degree in n for coefficients (default: 3)
--min-diff-deg <d>   Minimum derivative order (default: 0)
--max-diff-deg <d>   Maximum derivative order (default: 2)
--inhomogeneous      Also try inhomogeneous recurrences
--min-inhomo-var-deg Minimum degree in t for the inhomogeneous term
--max-inhomo-var-deg Maximum degree in t for the inhomogeneous term
--min-inhomo-idx-deg Minimum degree in n for the inhomogeneous term
--max-inhomo-idx-deg Maximum degree in n for the inhomogeneous term
--denominator        Allow a nontrivial LHS factor f(n,t)
--alternating-sign   Also allow right-hand-side terms multiplied by (-1)^n
--max-denom-var-deg  Max degree in t for f(n,t) (default: 2, implies --denominator)
--max-denom-idx-deg  Max degree in n for f(n,t) (default: 2, implies --denominator)
--min-margin <k>     Require equations >= unknowns + k (default: 1)
--fit-extra-rows <k> Extra rows after the first solvable prefix (default: 1)
--no-verify          Fit all rows instead of reserving held-out verification rows
--no-modular-prefilter
                      Disable default modular candidate rejection
--max-candidates <n> Inspect at most n adaptive candidate configurations
--json               Emit recurrence JSON with initial conditions
--python             Emit exact standalone Python code using Fraction arithmetic
--format json        Alias for --json
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

### Several recurrence representations

Use `--representations` to compare lag length against derivative order instead
of stopping at the first match:

```sh
polytool recurrence --representations --json \
  --max-rec-len 3 --max-diff-deg 2 --max-var-deg 3 --max-idx-deg 1 \
  --denominator --max-denom-var-deg 0 --max-denom-idx-deg 1 \
  --max-candidates 2000 < rows.txt > representations.json
```

The additive Rust APIs are
`recurrence::find_recurrence_representations_rational` and the `i64` convenience
wrapper `find_recurrence_representations`. Both take rows, existing
`AdaptiveSearchOptions`, and an `AdaptiveSearchBudget`, and return
`Result<RecurrenceRepresentations, String>`. The MCP tool
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

Rows after `skip_prefix` are renumbered from 1. `first_checked_index` is the
first checked target in this numbering; `first_checked_input_row` refers to
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

To generate rows from one entry, extract its nested object:

```sh
jq '.representations[0].recurrence_json' representations.json > recurrence.json
polytool recurrence-generate --recurrence recurrence.json --rows 60
```

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
The ordinary `node --test polytool/web/*.test.cjs` checks remain useful, but
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
terms, the LHS denominator if present, and the minimal initial polynomial rows.
Generate more rows from the saved recurrence with:

```sh
polytool recurrence-generate --recurrence recurrence.json --rows 100 > rows-100.txt
```

Use `--additional n` instead of `--rows n` to keep all initial rows and append
`n` newly generated rows.

The web interface exports this compact JSON record and links back to these
generation instructions instead of emitting a large standalone Python program.

Use Python output when a sequence extension script is more convenient than JSON:

```sh
polytool recurrence --python < polys.txt > recurrence.py
```

The generated Python uses exact `fractions.Fraction` arithmetic and dense
coefficient lists in ascending powers of `t`.
