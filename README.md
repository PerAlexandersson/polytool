# polytool

Dense univariate polynomial toolkit for combinatorial research. Provides
real-rootedness testing by adaptive subresultant integer root counting,
Bézout-matrix interlacing checks, exact Hurwitz stability, log-concavity,
gamma-positivity, resultants, Ehrhart/h\*-vector conversion, recurrence search
for polynomial sequences, finite Athanasiadis--Wagner interlacing matrices,
and standard sequences — all with exact arithmetic.

Rename note: this crate was renamed from `polynomial-tools` to `polytool`
before publication. Use `polytool` in Cargo manifests, Rust imports, CLI
documentation, and MCP server configuration.

## Installation

From the repository root:

```sh
cargo build --release -p polytool
```

The CLI binary is at `target/release/polytool`.

Use `polytool --version` (or `polytool -V`) for a single deterministic line:

```text
polytool 0.5.0 (git 0123456789ab)
```

The version comes from the crate manifest and the lowercase 12-hex-digit Git
commit is captured at build time.  If the source is built without trustworthy
Git metadata, the suffix is `(git unavailable)` rather than a fabricated
revision.  Reproducible packaging environments may set `POLYTOOL_GIT_COMMIT`
to a full or at least 12-digit hexadecimal revision; an empty or malformed
value deliberately selects the unavailable form.

On the generated standalone `main` branch, the equivalent command is:

```sh
cargo build --release
```

## Using as a crate

The library crate is usable today from another Rust project by depending on
this directory:

```toml
[dependencies]
polytool = { path = "../polytool" }
```

From GitHub, depend on the repository's standalone `main` branch:

```toml
[dependencies]
polytool = { git = "https://github.com/PerAlexandersson/polytool.git" }
```

After publishing to crates.io, the dependency would be:

```toml
[dependencies]
polytool = "0.3"
```

The public API is re-exported from `polytool`, so a downstream crate
can write:

```rust
use polytool::{check_weak_interlacing, is_real_rooted};

assert!(is_real_rooted(&[1, 11, 11, 1]));
assert_eq!(check_weak_interlacing(&[-1, 1], &[2, -3, 1]), Some(true));
```

Coefficient vectors are always in ascending degree order. Most convenience
functions accept `&[i64]`; use the `*_bigint_coeffs` variants when coefficients
may exceed `i64`. Interlacing functions return `Option<bool>`: `None` means the
directed degree relation is not valid for that test, not that interlacing failed.
For the undirected question "do `f` and `g` interlace in some order?", use
`interlaces_either_order`, which returns a plain `bool`.

### Standalone exact root counting

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

### Shared-denominator vector recurrences

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

To install the MCP server into a user bin directory and print a client
configuration snippet:

```sh
./mcp/install.sh
```

## CLI usage

### Input format

Polynomials are given as comma-separated integer coefficients in
**ascending degree order**: `a_0, a_1, ..., a_d` represents
`a_0 + a_1 t + ... + a_d t^d`.

A text file `polys.txt` might look like:

```
1, 11, 11, 1
1, 26, 66, 26, 1
1, 57, 302, 302, 57, 1
```

The commands `interlacing`, `interlacing-profile`, `properties`,
`gamma-expansion`, `bernstein-expansion`, `family-check`, `sequence`, `coefficient-tests`,
`hstar-to-ehrhart`, `ehrhart-to-hstar`, `hstar-inequalities`,
`cyclic-sieving`, `cyclic-sieving-sequence`, `bicyclic-sieving`, `tnn-network`,
`brenti-recognize`, `pf-network`, and `left-production` also accept `--json` for
machine-readable output.
Large integer coefficients in JSON output are serialized as strings.

The exact property commands `real-rooted`, `interlacing`,
`interlacing-profile`, `properties`, `gamma-expansion`, and the non-recurrence
checks inside `family-check` accept arbitrary-size integer coefficients.
Standard sequence generation, resultants, discriminants, Ehrhart conversion,
Stapledon decomposition, and the standalone `recurrence` command also use
arbitrary-size exact integers in the CLI. The main remaining compatibility
exception is `family-check --recurrence`, which currently requires
coefficients that fit in `i64`.

### Check real-rootedness

```sh
polytool real-rooted < polys.txt
```

`real-rooted` accepts arbitrary-size integer coefficients in dense-list and
expanded-polynomial input.

Output (one line per polynomial):

```
1 + 11t + 11t^2 + t^3: real-rooted
1 + 26t + 66t^2 + 26t^3 + t^4: real-rooted
1 + 57t + 302t^2 + 302t^3 + 57t^4 + t^5: real-rooted
```

### Check interlacing

Test whether consecutive polynomials interlace:

```sh
polytool interlacing < polys.txt
polytool interlacing --json < polys.txt
```

Count how many consecutive previous rows each polynomial interlaces, stopping
the backward scan at the first failure:

```sh
polytool interlacing-profile < polys.txt
polytool interlacing-profile --json < polys.txt
```

The JSON output includes `previous_count`, `checked_previous_count`, and
`interlacing_previous_count`, plus the checked pair reports.
Both commands accept arbitrary-size integer coefficients in dense-list and
expanded-polynomial input.

### Check unimodality, log-concavity, palindromicity, gamma-positivity

```sh
polytool properties < polys.txt
polytool properties --json < polys.txt
```

`properties` accepts arbitrary-size integer coefficients in both dense-list
and expanded-polynomial input.  When all coefficients fit in `i64`, it uses the
same fast exact path as the rest of the CLI; otherwise it falls back to exact
BigInt property checks.

Output includes all properties for each polynomial:

```
1 + 11t + 11t^2 + t^3: real-rooted, palindromic, gamma-positive [1, 8], unimodal, log-concave, ultra-log-concave
```

For a palindromic polynomial, print the full gamma expansion:

```sh
polytool gamma-expansion < polys.txt
polytool gamma-expansion --json < polys.txt
```

`gamma-expansion` accepts arbitrary-size integer coefficients.  In JSON output,
both input coefficients and gamma coefficients are serialized as strings.

Example output:

```text
1 + 11t + 11t^2 + t^3: gamma [1, 8]; expansion: (1+t)^3 + 8 t (1+t)
```

### Expand in the Bernstein basis

The standard degree-$n$ Bernstein basis is

```text
binomial(n,j) t^j (1-t)^(n-j),  0 <= j <= n.
```

The command uses the input polynomial's degree by default.  Pass `--degree`
to elevate every input row to a common ambient degree.  Inputs and outputs may
have arbitrary-size exact rational coefficients.

```sh
echo '1, 2, 3' | polytool bernstein-expansion --degree 3
echo '1, 2, 3' | polytool bernstein --degree 3 --json
```

The text output is the exact Bernstein coordinate row:

```text
1, 5/3, 10/3, 6
```

### Coefficient-only real-rootedness tests

```sh
polytool coefficient-tests < polys.txt
polytool coefficient-tests --json < polys.txt
```

This reports Newton inequalities and Kurtz's sufficient condition.  When the
Kurtz condition holds, the polynomial has distinct real roots.

### Generate standard sequences

```sh
polytool sequence eulerian 5
polytool sequence narayana 5 --json
```

Supported sequence names are `eulerian`, `narayana`, `type-b-eulerian`,
`chebyshev-t`, `chebyshev-u`, and `hermite`.
Generated coefficients use arbitrary-size integers.  The CLI and MCP tool list
`A_1, ..., A_n` and `N_1, ..., N_n` for `eulerian` and `narayana`, and
`B_0, ..., B_n`, `T_0, ..., T_n`, ... for the others; the library functions
`polytool::sequences::eulerian_polynomials_bigint` and
`narayana_polynomials_bigint` additionally include the row `n = 0`.

### Generate recurrence-backed OEIS families

The bundled OEIS catalog exposes one Rust function per A-number, for example
`polytool::oeis::A008292()`, together with dynamic lookup by identifier.  Each
entry is generated from a sparse exact recurrence description and its initial
rows; installed binaries do not need an OEIS database, Python, or SymPy.

```sh
polytool oeis list
polytool oeis info A008292 --json
polytool oeis generate A008292 --rows 8
polytool oeis generate A008292 --rows 8 --format triangle
polytool oeis generate A008292 --rows 40 --format bfile --max-terms 1000
```

Generation formats are `rows`, `triangle`, `polynomial`, `json`, `jsonl`,
`csv`, and `bfile`.  `--start-row` selects a later displayed OEIS row.  Strict
b-file output is enabled only when the bundled complete-row prefix has been
matched against OEIS data; `--max-terms` stops before a row that would cross
the cap.  The JSON returned by `oeis info` keeps the human-readable
`recurrence` string and also includes `recurrence_data`, a directly consumable
`polytool.recurrence.v1` object with exact rational coefficients and initial
polynomials.  Catalog entries do not fabricate recurrence-search diagnostics,
so this embedded object omits the optional `search` member.  The catalog
distinguishes holdout-`verified` recurrences from
`validated` recurrences whose generated rows match the current OEIS prefix but
whose original fitting/holdout provenance is unavailable.  Both are enabled by
default. The original generator rejected imported queue and Lean definitions
that could not be aligned with complete OEIS rows; in particular, the catalog
does not publish a recurrence fitted to a row fragment or to coefficient data
with positional zeros removed. The `experimental` status and `--include-experimental` option
remain available for explicitly curated future entries, but the current
bundled catalog contains none.

#### Catalog data, provenance, and adding sequences

The catalog's single source of truth is the tracked data file
`data/oeis-catalog.jsonl`: one JSON object per sequence, sorted by A-number,
with a fixed key order so that diffs stay small.  The build script
(`build.rs`) validates it and generates the static catalog, the accessor
functions and their rustdoc titles; no external data or tooling is needed to
build.  The record schema is public as `polytool::oeis::OeisCatalogRecord`.
Each record holds the executable recurrence (`layout`, `first_row`,
`flattened_offset`, `bfile_prefix_verified`, `prefix_rows`,
`recurrence_first_index`, `recurrence_first_width`, `initial_rows`,
`recurrence`, with recurrence coefficients as sparse `[n_degree, t_degree,
"value"]` monomials) and a `provenance` object with the `status`, the
`source_rows` / `verification_rows` counts and the `rows_sha256` fingerprint.
The `source` field records where the recurrence came from:

| `source` | Entries | Origin |
| --- | --- | --- |
| `benchmark_fixture` | 73 | `polytool recurrence --json` on `fixtures/recurrence-benchmarks` (held-out verified); also records `fixture_slug` |
| `lean_definition` | 630 | Fitted by an earlier version of the polytool recurrence search and transcribed into generated Lean definitions; validated against the OEIS data prefix (16 rows); records the last compared row as `validation_row` |
| `curated_queue` | 71 | Fitted by recurrence search on database-cached rows; `verified` when checked on held-out rows |
| `polytool_recurrence_search` | new | `polytool oeis add` (below); records `fit_rows`, `checked_terms`, `search_options` and the polytool version |

The first three tiers were produced by an earlier offline generator, which
combined the benchmark fixtures, a curated recurrence queue, and generated
Lean definitions with a local OEIS copy; its output was migrated to the data
file unchanged (a test compares every entry's definition and first 20 rows
with fingerprints captured before the migration).

To add sequences, run from the source directory with a local copy of the OEIS
in its internal format (`seq/A123/A123456.seq` files):

```sh
polytool oeis add A008275 --oeis-data ~/oeis --dry-run   # print the entry only
polytool oeis add A008275 A001497 --oeis-data ~/oeis      # write data/oeis-catalog.jsonl
polytool oeis add --ids-file ids.txt --max-candidates 20000
polytool oeis add A055151 --layout table --row-widths 1+1/2 --replace
polytool oeis add A123456 --bfile b123456.txt --name "..." --layout triangle
polytool oeis check --oeis-data ~/oeis                    # re-validate every entry
```

`oeis add` reads the name, offset, keywords and terms (or an explicit b-file),
cuts the terms into complete rows (regular triangles for `tabl`; any other
sequence needs `--layout triangle` or `--layout table --row-widths W+S[/R]`),
runs the adaptive recurrence search with the options of `polytool recurrence`
on a prefix of the rows (`--fit-rows`, by default all but `max(2, rows/4)`),
and then replays the candidate entry through the catalog's own row generator
and compares it with every available term, held-out rows and a trailing
partial row included.  The status is `verified` when held-out rows were
checked and `validated` otherwise.  Each A-number gets one outcome: `added`,
`updated` (with `--replace`), `already_present`, `no_recurrence`,
`layout_unsupported`, `insufficient_data`, `mismatch`, or `unavailable`.
Rebuild polytool to bundle new entries.  `oeis check` regenerates entries
(bundled, or from `--catalog FILE`) and compares them term by term with the
local OEIS data, listing entries it could not check.  A finite match is
evidence, not a proof that the recurrence holds for every row.

Release packages omit repository-only agent notes, handoffs, CI configuration,
and benchmark output; they include the catalog data file, from which the build
script generates the catalog. The recurrence fixtures stay
in the package deliberately: the installed CLI's `bench recurrence-fixtures`
command uses them as its default self-contained corpus.

### Check Family H PF/Jensen pencils

```sh
polytool pf-pencil --case A036969 --degree 5
polytool pf-pencil --case A036969 --min-degree 1 --max-degree 10 --json
polytool pf-pencil --all-family-h --max-degree 8 --json
```

`pf-pencil` builds the fixed Family H coefficient actions used in
a separate real-rootedness proof project.  For an old-row degree `d`, it forms
the integer-cleared Jensen endpoint kernels

```text
J_alpha,d(t) = sum_k binomial(d,k) alpha(k,d) t^k
t J_beta,d(t) = t sum_k binomial(d,k) beta(k,d) t^k
```

and checks exact real-rootedness of both endpoints and of the finite positive
pencil specializations `J_alpha,d + lambda * t J_beta,d`.  The default
lambda list is `0,1,10,100`; override it with `--lambdas 0,2,50`.

Supported built-in cases are:

```text
A036969 A071951 A080248 A156289 A160562 A269945
A166960 A166961 A166962 A166972 A191935
A371081 A371259 A390433 A198204
```

JSON output uses schema `polytool.pf-pencil.v1`.  Coefficients and
lambda values are serialized as strings.  Each item reports the case id,
degree, alpha/beta labels, common denominator, endpoint coefficient vectors,
endpoint real-rootedness booleans, lambda checks, and `finite_evidence_ok`.
Common-factor/residual extraction is not yet implemented; the full endpoint
polynomials are reported instead.

### Check a family in one pass

`family-check` reports properties, consecutive weak/strict interlacing, and
optionally an adaptive recurrence search.  Requirement flags make the command
exit nonzero at the first failed requested condition.
Property checks and consecutive interlacing accept arbitrary-size integer
coefficients.  Recurrence search currently requires coefficients that fit in
`i64`; oversized rows are preserved in the report and recurrence is marked
unavailable rather than silently dropped.

```sh
polytool family-check \
  --require-real-rooted \
  --require-weak-interlacing \
  --recurrence \
  --json < polys.txt
```

### Planar networks: TNN triangles and strip digraphs

`tnn-network` reads monic rows `P_0, P_1, ...` with `deg P_n = n` (or
`P_1, P_2, ...`, in which case `P_0 = 1` is prepended) and decides exactly
whether the coefficient triangle `r_{n,k} = [t^k] P_n` is lower-unitriangular
and totally nonnegative (TNN).  On success it prints the weights
`lambda[n][k]` of the canonical planar network `Gamma_N(lambda)` (the weight of
the vertical edge `(n+1,k) -> (n,k)`; horizontal edges have weight 1) and
re-evaluates the network's path matrix against the input (reconstruct →
evaluate → compare).  On failure it says why: a shape problem (not monic, wrong
degree, entry above the diagonal), a negative entry, or an explicit negative
minor with zero-based rows, columns and exact value.  With `--matrix` each line
is a matrix row instead (a square lower-triangular matrix, or triangular rows
with `n + 1` entries in row `n`).  `--chain` adds the chain polynomials
`p_0 = 1`, `p_n = t * sum_{k<n} r_{n,k} p_k` and checks that each is
real-rooted with zeros in `[-1, 0]` and that `p_n` interlaces `p_{n+1}`.

```sh
printf '1,1\n1,2,1\n1,3,3,1\n' | polytool tnn-network --chain
```

```text
TNN network check
- Input: 3 polynomial row(s); P_0 = 1 prepended
- Coefficient matrix (triangular rows):
  0: 1
  1: 1, 1
  2: 1, 2, 1
  3: 1, 3, 3, 1
- Verdict: lower-unitriangular and totally nonnegative (TNN)
- Canonical network weights lambda[n][k] (edge (n+1,k) -> (n,k)):
  n=0: 1
  n=1: 1, 1
  n=2: 1, 1, 1
- Certificate (reconstruct -> evaluate -> compare): the network's path matrix equals the input

Chain polynomials (Brändén--Saud Maia Leite, arXiv:2412.06595, Theorem 3.7):
  p_0 = 1: real-rooted, zeros in [-1, 0]
  p_1 = t: real-rooted, zeros in [-1, 0]
  p_2 = t + 2t^2: real-rooted, zeros in [-1, 0]
  p_3 = t + 6t^2 + 6t^3: real-rooted, zeros in [-1, 0]
- Consecutive interlacing p_n << p_(n+1): all hold
- Theorem conclusion (real-rooted in [-1, 0] and interlacing): holds

Note: TNN of the coefficient matrix does not imply that the row polynomials are real-rooted: ...
```

The guardrail example `[[1,0,0],[1,1,0],[1,1,1]]` (rows `1`, `1 + t`,
`1 + t + t^2`) is TNN with weights `lambda = [[1], [1, 0]]`, although its last
row has no real zeros.  The report makes no claim about the rows:

```sh
printf '1,0,0\n1,1,0\n1,1,1\n' | polytool tnn-network --matrix
printf '1,0,0\n1,1,0\n2,0,1\n' | polytool tnn-network --matrix --json
```

The second matrix is rejected with the witness `rows [1, 2], columns [0, 1]`,
minor `-2`.  JSON output (schema `polytool.tnn-network.v1`) has the fields
`input_mode`, `input_rows`, `prepended_p0`, `coefficient_matrix`, `tnn`,
`failure` (`kind`, `message`, `row`, `col`, `value`, `witness` with `rows`,
`cols`, `value`), `lambda`, `certificate_verified`, `chain` (`theorem`, `error`,
`polynomials` with `real_rooted` and `roots_in_minus_one_zero`, `interlacing`,
`conclusion_holds`, `contradicts_theorem`) and the fixed `note`.  Integers and
rationals are strings.  Exit status: 0 if the triangle is TNN and the
certificate re-evaluates, 1 if it is not lower-unitriangular TNN, 2 for invalid
input or options.

`brenti-recognize` searches exactly for a planar strip digraph (edges from
column `x` to `x + 1` that depend only on the height, with monotone targets)
whose path counts from `(0,0)` to `(n,k)` are the coefficients of the input
rows `P_1, P_2, ...` (`P_0 = 1` implicit; `--first-row 0` reads `P_0, P_1, ...`).
Bounds mirror `BrentiRecognitionOptions`: `--max-span <k>`,
`--weights rational|integer`, `--max-search-nodes <n>`,
`--max-integer-nodes <n>`.  The verdict is `certified` (the digraph is
re-verified against every row), `impossible_within_bounds` (exact) or
`inconclusive` (a budget ran out); exit status 0, 1 and 3 respectively, and 2
for invalid input or options.

```sh
printf '0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n' | polytool brenti-recognize --weights integer
printf '1,1\n1,4,1\n1,11,11,1\n' | polytool brenti-recognize --json
```

```text
Strip-digraph recognition
- Input: 4 row(s), digraph rows 1..4; P_0 = 1 implicit
- Bounds: integer weights, target span unbounded, 1000000 search nodes, 10000 integer nodes per column
- Verdict: certified
- Digraph (height -> target: weight):
  0 -> 1: 1
  1 -> 1: 1, 2: 1
  2 -> 2: 2, 3: 1
  3 -> 3: 3, 4: 1
  4 -> (no edges)
- Certificate re-verified against the input rows: yes
```

The Stirling rows above are certified; the Eulerian rows are
`impossible_within_bounds` (no nonnegative weights produce target height 1,
even without the planarity constraint).  JSON
output (schema `polytool.brenti-recognize.v1`) has `options`, `rows`,
`verdict`, `certified`, `digraph` (`num_heights`, `outgoing[y]` as
`{to, weight}`), `row_indices`, `certificate_verified`, `supports_checked`,
`obstruction`, `reason`, `message` and the fixed `note`.

Guardrails for both commands:

- TNN of the coefficient matrix does **not** make the row polynomials
  real-rooted: `[[1,0,0],[1,1,0],[1,1,1]]` is TNN and has the row
  `1 + t + t^2`.  Brändén--Saud Maia Leite (arXiv:2412.06595, Theorem 3.7)
  concerns the chain polynomials: they are real-rooted in `[-1, 0]` and
  interlace when the triangle is lower-unitriangular TNN.
- A strip-digraph certificate proves, by Lindström--Gessel--Viennot, that the
  coefficient matrix of `P_0, ..., P_N` is TNN.  It does not prove that the
  rows are real-rooted or Pólya frequency sequences (`[1 + t, 1 + t + t^2]` is
  certified): the transitions depend on the height, so the digraph is
  invariant under shifts in `x`, not in `y`, and Brenti's PF theorem needs
  height-shift invariance.
- Both commands certify only the supplied finite prefix.

### Planar networks: Pólya frequency rows

Two commands certify the *rows* themselves: each certified row is a Pólya
frequency (PF) sequence, i.e. real-rooted with only nonpositive zeros.

`pf-network` searches exactly for a planar network on `{0..N} x Z` that is
invariant under **height** shifts (Brenti, *Combinatorics and total
positivity*, JCTA 1995).  Into each column `n` it has steps
`(length, rise, weight)`, standing for the edges `(n - length, y) -> (n, y + rise)`
at every height `y`, so
`P_n = sum weight * t^rise * P_(n - length)` with `P_0 = 1`.  When every
positive-weight step is primitive (`gcd(length, rise) = 1`) and no two edges
meet away from a common endpoint, the Toeplitz matrix of every row is the path
matrix of fully compatible sources `(0, i)` and sinks `(n, j)`, so it is
totally nonnegative (Lindström--Gessel--Viennot) and the row is PF
(Aissen--Schoenberg--Whitney).  The search tries, column by column, every
inclusion-minimal set of steps that is planar together with the edges already
chosen and admits nonnegative weights for that row; this is exact for the
bounds `--max-length <L>` and `--max-rise <R>` (default 3 and 3).  Other
options: `--first-row <0|1>`, `--weights rational|integer`,
`--max-search-nodes <n>`, `--max-integer-nodes <n>`.  The certificate is
re-verified exactly: structure, planarity, every row, and the path matrix
rebuilt over the explicit lattice against the `(deg + 1) x (deg + 1)` Toeplitz
window of the row.  The PF conclusion comes from planarity and height
invariance, for Toeplitz sections of every size; a finite Toeplitz window that
happens to be TNN would prove nothing (the `3 x 3` window of `1 + t + t^2` is
TNN).  A row that is not real-rooted rules out every network.  Verdicts and
exit codes are as for `brenti-recognize` (0 certified, 1
`impossible_within_bounds`, 2 invalid input, 3 inconclusive).

```sh
printf '1,1\n1,2,1\n1,3,3,1\n' | polytool pf-network
printf '1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n' | polytool pf-network --json
```

```text
Height-invariant planar network (PF rows)
- Input: 3 row(s), network rows 1..3; P_0 = 1 implicit
- Bounds: rational weights, step length <= 3, rise <= 3, 1000000 search nodes, 10000 integer nodes per row
- Verdict: certified (every supplied row is PF)
- Steps into column n (length, rise): weight
  n=1: (1, 0): 1, (1, 1): 1
  n=2: (1, 0): 1, (1, 1): 1
  n=3: (1, 0): 1, (1, 1): 1
- Certificate re-verified (planarity, rows, Toeplitz windows): yes
```

The Eulerian rows `A_1, A_2, A_3` are certified with Brenti's recursive
pattern (steps `(1,0)`, `(1,1)`, `(2,1)`; `A_3 = (1+t) A_2 + 6t A_1`), but no
network with steps of length and rise at most 8 reproduces `A_1, ..., A_4`
(rational or integer weights); the same holds for the Narayana rows
`N_1, ..., N_4`.  Likewise no such network reproduces the Bell polynomials
`B_1, ..., B_4` (Stirling numbers of the second kind), which
`left-production` does certify.  These are statements about the stated
bounds only.  JSON output (schema
`polytool.pf-network.v1`) has `options`, `rows`, `verdict`, `certified`,
`network` (`last_column`, `steps[n-1]` as `{length, rise, weight}`),
`row_indices`, `certificate_verified`, `supports_checked`, `obstruction`,
`reason`, `message` and the fixed `note`.

`left-production` implements the unified planar network of Chen--Fu--Ruan
(arXiv:2512.08369).  Rows `P_0, ..., P_N` with `deg P_n = n` form a lower
triangular matrix `A` with nonzero diagonal; its left production matrix
`Q = A (1 ⊕ A^{-1})` is determined row by row by
`a_{n,0} = q_{n,0}`, `a_{n,k} = sum_{i=k..n} q_{n,i} a_{i-1,k-1}`, and row `n`
of `A` depends only on the leading block `Q_n`.  If `Q_N` is totally
nonnegative, then (Theorems 3 and 4 of the paper, applied to `Q_N ⊕ I`) the
rows `P_0, ..., P_N` are PF and `A_N` and its reversal `[a_{n,n-k}]` are TNN.
The command computes `Q_N` exactly, decides TNN (with a negative-minor witness
on failure), re-verifies the certificate (regenerates `A_N` from `Q_N`,
re-checks TNN, and compares the unified-network window
`M_{n,n}[n..2n; 0..n] = T_n(P_n)^T` for every row), and prints an explicit
planar network for `Q_N = L D` (canonical network weights of `L`, diagonal
`D`).  It covers Stirling numbers of both kinds, Lah, `r`-Whitney and
idempotent numbers, and Brenti's `n`-recursive triangles (Pascal, Delannoy,
derangements) without naming the family.  The criterion is sufficient, not
necessary: `Q_4` of the Eulerian triangle has `q_{4,2} = -5` although the rows
are real-rooted, and the Narayana triangle fails at `Q_4`.  `--first-row 0`
reads `P_0, P_1, ...`; by default `P_0 = 1` is prepended.  Exit status 0 if
certified, 1 if not, 2 for invalid input.

```sh
printf '0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n' | polytool left-production
printf '1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n' | polytool left-production --json
```

```text
Left production matrix (Chen--Fu--Ruan)
- Input: 4 row(s); P_0 = 1 prepended
- Production matrix Q_N (triangular rows):
  0: 1
  1: 0, 1
  2: 0, 1, 1
  3: 0, 1, 2, 1
  4: 0, 1, 3, 3, 1
- Verdict: Q_4 is TNN; certificate re-verified
- Conclusion: rows P_0..P_4 are PF (real-rooted, nonpositive zeros); the triangle and its reversal are TNN
- Planar network for Q_N = L D (lambda[n][k] of Gamma_N(lambda) for L):
  n=0: 0
  n=1: 0, 1
  n=2: 0, 1, 1
  n=3: 0, 1, 1, 1
  column scales D = diag(1, 1, 1, 1, 1)
```

JSON output (schema `polytool.left-production.v1`) has `first_row_index`,
`rows`, `triangle`, `production_matrix`, `tnn`, `failure_reason`,
`negative_minor` (`rows`, `cols`, `value`), `certificate_verified`,
`rows_certified_pf`, `production_network` (`lambda`, `column_scales`),
`message` and the fixed `note`.

Guardrails for these two commands:

- They certify that the supplied rows are PF.  `pf-network` says nothing
  about the coefficient matrix of the rows (`tnn-network` and
  `brenti-recognize` concern that matrix, and do not certify PF rows).
- `impossible_within_bounds` from `pf-network` is relative to the step
  bounds, unless the obstruction is a row that is not real-rooted.  A
  `left-production` failure says nothing about the rows.
- Both certify only the supplied rows; `1 + t + t^2` is never certified.

### Search for a recurrence

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

#### Several recurrence representations

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

#### Complete recurrence spaces and certificates

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

#### Browser verification

The recurrence browser regression drives the actual page, Web Worker and WASM
at desktop and 375-pixel widths. It checks guided full-search controls, exact
row parsing and origins, example formulas, clipboard copying, backend replay,
partial-result labels, cancellation, tab routing, stale-input invalidation,
renderer-error recovery and the older finder modes.
It uses a temporary localhost server and closes its own browser and server.
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

#### Automatic variable degree, collective leaders, and indexed sources

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

#### Complete spaces after factorial coefficient rescaling

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

#### Machine-readable recurrence JSON and row generation

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

#### Coordinate transformations and recurrence transport

`polytool recurrence-transform` applies an exact change of polynomial
coordinates to a sequence, transports its defining recurrence into the new
coordinates, and can search independently there. It reads one JSON request on
stdin and prints a `polytool.recurrence-transform.v1` report; `--verify`
replays a report exactly from its embedded request (exit status 2 on a
mismatch or invalid input), and `--compact` prints compact JSON.

| `kind` | Row map | Precondition / inverse image |
| --- | --- | --- |
| `multiply_power {power}` | `Q_n = x^r P_n` | inverse needs `x^r` to divide `Q_n` |
| `divide_power {power}` | `Q_n = P_n / x^r` | `x^r` divides `P_n` |
| `reversal {schedule}` | `Q_n = x^(d_n) P_n(1/x)` | `deg P_n <= d_n` |
| `affine {a, b}` | `Q_n = P_n(a x + b)`, `a != 0` | none |
| `mobius {a, b, c, e, schedule}` | `Q_n = (c x + e)^(d_n) P_n((a x + b)/(c x + e))`, `a e - b c != 0` | `deg P_n <= d_n` |
| `f_polynomial {schedule}` | `Q_n = (1 + x)^(d_n) P_n(x/(1 + x))` | `deg P_n <= d_n` |
| `gamma {schedule}` | `P_n = sum_i Q_(n,i) x^i (1 + x)^(d_n - 2i)` | `x^(d_n) P_n(1/x) = P_n`; inverse needs `deg Q_n <= d_n/2` |
| `symmetric_decomposition {schedule}` | `P_n = a_n + x b_n`, state `(a_n, b_n)` | `deg P_n <= d_n`; inverse needs `a_n`, `b_n` symmetric at `d_n`, `d_n - 1` |

Rational parameters are strings such as `"3/2"`. A schedule
`{"first_index": n0, "first_degree": d0, "step": h, "last_index": n1}` declares
`d_n = d0 + h (n - n0)` on the index domain `n0 <= n (<= n1)`, so
`d_(n+1) = d_n + h`; a decreasing schedule needs `last_index`, and negative
degrees, empty domains and rows outside the domain are rejected. The declared
bound is used even when a row has lower degree (`f_polynomial` calls
`decomposition::f_polynomial` with `d_n`, never the observed degree), so
missing top coefficients and endpoint factors are preserved exactly.

The source is a defining `polytool.recurrence.v1` rule with seeds
(`source_recurrence`), a bundled OEIS family (`oeis_id`), or explicit `rows`
with `row_origin`; with a source, `row_count` rows are generated (default 10).
Every row is mapped, the inverse reconstructs it exactly
(`reconstruction_verified`), and a row violating a precondition rejects the
request. Each condition is reported as `certified_all_indices` when the source
rule proves it for the whole domain (the transported rule is denominator-free
and the seeds qualify; for gamma, the reversal-transported rule equals the
source rule) and otherwise as `checked_supplied_rows`, which is finite
evidence.

A rule `f(n,x) P_n = sum s(n) c(n,x) D^k P_(n-r) + g(n,x)` is transported
exactly. For a Möbius map with `phi = (a x + b)/ell`, `ell = c x + e`,
`Delta = a e - b c`, a first-order rule `R P_(n+1) = A P_n + B P_n'` becomes

```text
R Q_(n+1) = ell^h (A(phi) - d_n c ell B(phi)/Delta) Q_n + ell^(h+2) B(phi)/Delta Q_n'
```

Higher derivatives and lags follow from the product and chain rules; the
identity is multiplied by the power of `ell` that clears every denominator and
a common power of `ell` is then divided out exactly, never evaluated at the
pole. An additive term `g` is transported only when `ell^(d_n)` does not
depend on `n` (`h = 0`, or constant `ell` with `ell^h = 1`). The symmetric
decomposition turns a first-order rule without `(-1)^n` factors into the
complete two-component shared-denominator vector system for `(a_n, b_n)`.
Gamma coordinates are not transported (the substitution `y = x/(1+x)^2` is
quadratic); search the gamma rows instead.

The transported rule (`evidence: "transported_identity"`, a replayable
`polytool.recurrence.v1` record with transformed seeds) is an exact
consequence of the supplied source rule and is kept separately from the
optional independent search (`"search": {}` or explicit bounds,
`evidence: "finite_fit"`). Canonical Eulerian rows A008292 (`E_1 = 1`, degree
`n - 1`) under the `f`-polynomial map give `Q_1 = 1`,
`Q_(n+1) = (1 + 2x) Q_n + x (1 + x) Q_n'`:

```sh
echo '{"oeis_id":"A008292","row_count":8,"search":{},
  "transformation":{"kind":"f_polynomial",
    "schedule":{"first_index":1,"first_degree":0,"step":1}}}' \
  | polytool recurrence-transform > report.json
polytool recurrence-transform --verify < report.json
```

In Rust, `RecurrenceTransformation::{map_row, inverse_row, transport,
record}` act on single rows and rules, and
`recurrence::{transform_recurrence, verify_recurrence_transform}` build and
replay reports. The MCP tools `transform_recurrence` and
`verify_recurrence_transform` take `{"request": {...}}` and
`{"report": {...}}`. Root-location information is opt-in (next section).

#### Root-domain information for transformed results

Algebraic discovery never needs an interval or a real-rootedness proof. An
optional `"root_domain": {}` (or explicit fields) in a `recurrence-transform`
request adds a `root_domain` section to the report; requests and reports
without it are unchanged. The section keeps apart

- `implications`: theorem adapters for the map, each with its `direction`
  (`equivalent`, `source_to_target`, `target_to_source`), hypotheses with a
  status (`certified_all_indices`, `checked_supplied_rows`, `failed_on_rows`,
  ...), degree bounds, coefficient-sign hypotheses, endpoint conventions, zero
  cases, index domain and, for one-directional maps, why the converse fails;
- `claims` about the source rows `P_n` or the transformed rows, each with an
  `evidence` label: `finite_exact_observation` (exact root isolation of the
  listed rows only), `supplied_assumption` (from the request, checked on the
  rows), `proposed_interval` (a candidate window), `completed_recurrence_induction`
  (a closed induction over the transported rule for every index) and
  `certified_family_invariant` (transferred to every index by a map theorem
  whose conditions are certified); unproved dependencies are listed in
  `conditional_on`;
- `recurrence_theorem`: the match of the transported rule against Rolle's
  theorem for derivatives of products, `Q_m = s(m) D^j (g Q_(m-1))`, with every
  interval attempt (bounded by `max_candidate_intervals`);
- `observations`: exact rational root enclosures with multiplicities for the
  first `max_observed_rows` rows.

| Adapter | Statement | Direction |
| --- | --- | --- |
| `monomial_factor` | `Q_n = x^r P_n`; only the multiplicity of the root 0 changes | equivalent |
| `real_mobius_equivalence` (`reversal`, `mobius`, `f_polynomial`), `real_affine_equivalence` | real Möbius maps with `a e - b c != 0` and `deg P_n <= d_n`; the pole `-e/c` is a root of multiplicity `d_n - deg P_n`, roots at `a/c` lower `deg Q_n`; intervals and rays map exactly | equivalent |
| `nonnegative_coefficient_window` | nonnegative coefficients put the real roots of `P_n` in `(-inf, 0]`; for the `f`-polynomial every real root of `Q_n` lies in `[-1, 0]` (`-1` records missing top degree, `0` the valuation) | proposes a window only |
| `gamma_quarter_bound` | `P_n` is real-rooted iff every root of the gamma polynomial is real and at most `1/4` (`G = 1 - z` at `d = 2` gives `1 + x + x^2`) | equivalent |
| `hermite_kakeya_obreschkoff_symmetric_pair` | real-rooted interlacing `a_n`, `x b_n` give real-rooted `P_n = a_n + x b_n` | components to `P_n` only |
| `rolle_derivative_of_product` | roots of `g` and `Q_(m-1)` in an interval `I` put the roots of `s(m) D^j(g Q_(m-1))` in `I` (zero allowed) | recurrence step |

Real-rooted components alone do not suffice for the symmetric pair
(`a = (1 + x)^3`, `b = x` give `1 + 3x + 4x^2 + x^3`), and its converse fails
(`(1 + 2x)^2` has `a = 1 + x + x^2`). Invertible maps need not preserve
real-rootedness both ways: for `P = (xQ)'`, real-rooted `Q` gives real-rooted
`P`, but `(1 + x)^2` has the antecedent `1 + x + x^2/3`.

The induction tries supplied `candidate_intervals`, then the default window,
then the hull of the roots of `g` and the seeds, checking `g` and every seed
exactly. Canonical Eulerian rows A008292 under the `f`-polynomial map become
`Q_(n+1) = (x(1 + x) Q_n)'`, `Q_1 = 1`; the induction closes on `[-1, 0]`
(`completed_recurrence_induction`), and since `deg P_n <= d_n` is certified
and `deg Q_n = d_n` follows from the degree formula, the Möbius equivalence
gives real-rooted Eulerian polynomials with roots in `(-inf, 0]`
(`certified_family_invariant`):

```sh
echo '{"oeis_id":"A008292","row_count":8,"root_domain":{},
  "transformation":{"kind":"f_polynomial",
    "schedule":{"first_index":1,"first_degree":0,"step":1}}}' \
  | polytool recurrence-transform > report.json
polytool recurrence-transform --verify < report.json
```

`root_domain` fields: `candidate_intervals`, `assumptions`, `default_window`
(default `true`), `max_candidate_intervals` (8), `max_observed_rows` (32),
`max_observed_degree` (128). An interval is `{"lower": {...}, "upper": {...}}`
with endpoints `{"value": "-1", "per_index": "0", "closed": true}`; an omitted
end is infinite (a ray) and `value + per_index * n` depends on the index (for
`P_n = x + n` the roots `-n` escape while the `f`-polynomial roots stay in
`(-1, 0]`). An assumption is `{"subject": "source" | "transformed",
"component": 0, "property": "real_rooted" | "real_rooted_in" | "real_roots_in",
"interval": ..., "first_index": ..., "last_index": ..., "note": ...}`; it is
checked exactly on the observed rows (`assumption_checks`), transferred through
equivalences with its label in `conditional_on`, and replayed with the report.
Finite exact checks are never promoted to statements about every index.

`polytool root-enclosures` (MCP `exact_root_enclosures`) is the exact interval
interface: for `{"polynomials": [["0", "0", "1", "3", "3", "1"]],
"max_width": "1/1000", "interval": {...}}` it returns enclosures with
multiplicities (`lower == upper` is an exact rational root; otherwise exactly
one root lies strictly inside and neither endpoint is a root), degree,
valuation, `real_rooted` and exact interval verdicts. The MCP `real_roots`
midpoints are representatives only. In Rust use
`real_rootedness::{isolate_real_roots_exact, RealInterval}`; fixed-cap diagonal
multiplier tests use `recurrence::jensen_fixed_cap_check` (Jensen polynomial
`T[(1 + x)^N]`, Schur--Malo composition; evidence `fixed_cap`, never a proof
for larger degrees).

#### Implicit differential relations, diagonal maps and reconstruction

Some coordinate changes turn a successor rule into an implicit relation with
derivatives on the newest row, or into a differential equation within one row.
An implicit relation is `sum c(n,x) D^d Q[n+i] = 0` at base index `n`
(JSON `{"terms": [{"shift": i, "deriv_order": d, "coeff": {"coeffs": ...}}]}`),
and a coefficient condition is `sum c(n) [x^k] Q[n+i] = 0`
(`{"terms": [{"shift": i, "coefficient_index": k, "coeff": [...]}]}`).
`polytool recurrence-implicit` maps a sequence with a defining rule
(`source_recurrence` or `oeis_id`) by one of two registry maps, which
`recurrence-transform` also accepts for row mapping:

| `kind` | Row map | Inverse |
| --- | --- | --- |
| `diagonal {weight}` | `Q_n = T_n P_n`, `T_n(x^k) = lambda(n,k) x^k` | divide coefficient `k` by `lambda(n,k)` |
| `derivative_antecedent {g}` | `Q_n` with `P_n = D^m(g Q_n)/m!`, `m = deg g` | `P_n = D^m(g Q_n)/m!` |

Weights are `factorial_power {power}` (`1/(k!)^m`), `binomial {m}`
(`binom(k+m,m)`, so `Q = W_m P` with `W_m = D^m x^m/m!`), `inverse_binomial {m}`
(`P = W_m Q`) and `product {scale, factorial_exponent, factors}` for
`scale (k!)^p prod f_j(n,k)^(e_j)`. A weight factor with a certified zero
(for example `k - 2`) is rejected; a factor whose zero set depends on `n`
is reported as `checked_supplied_rows` in `map_conditions`, and the
transported components then carry the same finite evidence. The antecedent
map is triangular with diagonal `lc(g) binom(k+m,m)`, so antecedents are found
by exact back substitution.

With `theta = x D`, a relation multiplied by `x^J` is written as
`sum x^s a(n,theta) S^i` and transported by
`x^s a(n,theta) S^i -> x^s a(n,theta) lambda(n,theta+s)/lambda(n+i,theta) S^i`.
One common multiplier `H(n,theta) = lcm b_(i,s)(n, theta - s)` clears the
denominators of every column (`H(theta) x^s = x^s H(theta+s)`). For the
original residual `Z_n`, `Z_n = 0` is equivalent to `H(n,theta) Z_n = 0`
together with `[x^k] Z_n = 0` at every `k >= 0` with `H(n,k) = 0`. These
exceptional indices are certified uniformly in `n` (exact nonnegative integer
roots of an index-free factor, or a one-signed shifted factor); otherwise the
clearing record says `checked_supplied_rows`, lists the factor and the zeros
seen on the supplied rows, and makes no equivalence claim. The derivative
antecedent transports a relation by substitution, `L T_g Q = 0`.

The report keeps evidence apart. `transports` holds each transported relation
(`transported_identity`), its clearing certificate and, when `P_n = T Q_n` with
a differential operator `T` (`W_m`, `T_g`), a reconstruction certificate
`L o T = sum_nu H_nu o S^nu o M` whose multipliers `H_nu = sum_d h_d(n,x) D^d`
record their differentiations. `description` is a generating description:
seeds, relations, `exceptional_index_conditions` and supplied
`boundary_conditions`. `assessment` gives the status (`equivalent_description`,
`certified_consequence`, `finite_evidence`, `inconsistent_seeds`) and a
uniqueness argument that is separate from annihilation. The argument picks a
relation of positive leader shift and compares coefficients of a difference
of two solutions. In the degree grading the top coefficient is multiplied by
`psi(n,k)`. In the valuation grading the lowest coefficient is multiplied by
`phi(n,k)`, and every certified exceptional index needs a condition that pins
`[x^k]` of the new row. Singular leaders and uncovered indices are listed as
reasons. `additional_relations` stores supplied or transported row equations
without treating them as successor rules. `"search": {}` fits a successor rule
to `Q` (`finite_fit`). The fit becomes `certified_by_reconstruction` when the
identity `L o T = sum H_nu S^nu M` holds, the source rule determines `P` from
its seeds, the generated `T Q*` agrees with the source seeds, and `T` is
injective.

Canonical binomial rows A007318 divided coefficientwise by `k!` give
`Q'_(n+1) = Q'_n + Q_n` with the exceptional condition
`Q_(n+1)(0) = Q_n(0)` and `Q_0 = 1`, an equivalent description; without the
constant condition (`--assess` on the edited description) the status drops to
`certified_consequence`. The row equation `x Q'' + (1 + x) Q' - n Q = 0` can be
stored through `target_relations` or transported from `(1 + x) P' - n P = 0`
through `source_relations`. Canonical A065826 divided by `k + 1` gives the
Eulerian rows A008292, and the fitted Eulerian rule is proved by
reconstruction:

```sh
echo '{"oeis_id":"A065826","row_count":9,"search":{},
  "transformation":{"kind":"diagonal","weight":{"kind":"inverse_binomial","m":1}}}' \
  | polytool recurrence-implicit > report.json
polytool recurrence-implicit --verify < report.json
```

`--verify` replays the report exactly and re-checks every reconstruction
certificate symbolically and on the rows. In Rust, see
`recurrence::{ImplicitRelation, CoefficientCondition, GeneratingDescription,
assess_generating_description, DiagonalWeight,
RecurrenceTransformation::transport_implicit, find_reconstruction_certificate,
ReconstructionCertificate, implicit_transform_recurrence,
verify_recurrence_implicit}`. The MCP tools are `transform_recurrence_implicit`,
`verify_recurrence_implicit` and `assess_generating_description`. A
general search over implicit relations is not part of this release: fits are
successor rules, and implicit relations come from transport or input.

### Scout BKW equal-modulus loci

For a polynomial recurrence with characteristic symbol

```text
F(x,z) = a_0(x) + a_1(x) z + ... + a_r(x) z^r,
```

`bkw-scout` numerically scans a complex rectangle for points where the two
dominant characteristic roots have nearly equal modulus.  This is a scout for
Beraha--Kahane--Weiss accumulation obstructions; it does not certify dominance,
amplitudes, or eventual non-real-rootedness.

Give the symbol as the z-coefficient polynomials in ascending z-degree:

```sh
polytool bkw-scout \
  --symbol '1; -x; 1' \
  --box -3 3 -2 2 \
  --grid 61 \
  --top 10
```

The same symbol can be supplied on stdin, one coefficient polynomial per line.
Use `--format json` for machine-readable output, and `--mathematica` to print a
`Reduce` skeleton for exact equal-modulus follow-up.

### Compute resultant/discriminant

```sh
# Resultant of two polynomials (given as two lines)
echo "1, 0, 1
1, -1" | polytool resultant

# Discriminant of a single polynomial
echo "1, 0, -3, 1" | polytool discriminant
```

Inputs for these CLI commands may have arbitrary-size integer coefficients.
The exact resultant or discriminant output is printed as an integer.

### Enumerate real-rooted integer polynomials

`enumerate-real-rooted` lists **every** real-rooted integer polynomial
`p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n` of degree `n` whose first `t`
leading coefficients are prescribed, for `3 <= t <= n + 1` and `a_0 != 0` of
either sign.  Reducible polynomials and polynomials with repeated roots are
included.  The search is the Greaves--Syatriadi derivative-tower algorithm
(arXiv:2504.09241) with certified endpoints (see the library section below);
nothing is read from stdin.

```sh
# All real-rooted x^3 - 3x + c: c = -2, ..., 2
polytool enumerate-real-rooted --degree 3 --leading-desc 1,0,-3

# Monic, trace 0, a_2 = -15 (degree 6), as JSON with exact string coefficients
polytool enumerate-real-rooted --degree 6 --leading-desc 1,0,-15 --format json
```

The prefix `--leading-desc` is in **descending** order (`a_0, a_1, ...`), while
every output polynomial is in the usual **ascending** order `c_0, ..., c_n`
(`c_k = a_(n-k)`), one per line, so the rows can be piped into
`polytool real-rooted`.  Rows are sorted lexicographically by
`(a_t, ..., a_n)`.  Budgets: `--max-results` (default 10000), `--max-nodes`
(default 1000000) and `--max-coefficient-bits` (default none); pass `none` to
lift one.  The status is `complete` (the full set, possibly empty) or
`resource_limit` (an initial segment of the full sorted set; JSON
`limit_reached` names the budget).  Invalid arguments exit with code 2 and,
with JSON output, print `{"status":"invalid_input", ...}`.  Text and JSON
output report the active limits and the search statistics.  With only two
prescribed coefficients the set can be infinite (`x^(n-2) (x^2 - k)` for all
`k >= 0`), which is why `t >= 3` is required.  A complete enumeration is a
certified finite computation for one prefix, not a statement about a family.

### Ehrhart ↔ h\*-vector conversion

```sh
# h*-vector → Ehrhart polynomial
echo "1, 8, 35, 32, 9" | polytool hstar-to-ehrhart

# Ehrhart polynomial → h*-vector (coefficients as rationals: num/den)
echo "1, 2, 2" | polytool ehrhart-to-hstar
```

For `hstar-to-ehrhart`, h\*-vector entries may be arbitrary-size integers.
For `ehrhart-to-hstar`, input coefficients may be exact rationals and the
resulting h\*-vector entries are arbitrary-size integers.

Check named h\*-vector inequalities:

```sh
echo "1, 4, 1, 0" | polytool hstar-inequalities --dimension 3
echo "1, 20, 1, 0" | polytool hstar-inequalities --dimension 3 --json
```

The report includes basic Ehrhart h\*-conditions, Stanley and Hibi inequalities,
the Balletti--Higashitani universal Scott inequality when applicable, and
Stapledon-derived nonnegativity checks.  Failures include the named inequality
and a reference.

### Cyclic sieving checks

For one polynomial and one cyclic group order:

```sh
echo "1, 1" | polytool cyclic-sieving --order 2 --fixed-counts "2,0"
```

Without fixed counts, the command profiles exact root-of-unity evaluations.  For
a sequence, the default candidate orders for row index `n` are
`n-2,n-1,n,n+1,n+2,n+3`:

```sh
polytool cyclic-sieving-sequence --first-index 2 < polys.txt
```

Optional fixed-count files use lines of the form:

```text
5 7: 10, 0, 2, 0, 2, 0, 0
```

#### Orbit decomposition and the Alexandersson--Amini criterion

Every `cyclic-sieving` run also decomposes the polynomial into cyclic orbits,
following P. Alexandersson and N. Amini, *The cone of cyclic sieving
phenomena*, Discrete Math. 342 (2019) 1581--1601 (arXiv:1804.01447).  All
arithmetic is exact and integral.  Reduce `f` modulo `q^n - 1` to residues
`c_0, ..., c_{n-1}`.

* `f(zeta_n^d)` is an integer for all `d` iff `c_j = c_{gcd(j,n)}` for all `j`
  (Desarmenien); otherwise the non-constant gcd classes are listed.
* If integral, `f = sum_{r | n} a_r h_r (mod q^n - 1)` with
  `h_r = 1 + q^r + ... + q^{n-r}` and `a_r = sum_{e | r} mu(r/e) c_e`.  A CSP
  action needs exactly `a_r` orbits of size `n/r`, i.e.
  `S_k = k a_{n/k} = sum_{j | k} mu(k/j) f(zeta_n^j)` elements in orbits of
  size `k`, and `f(zeta_n^d) = sum_{(n/r) | d} (n/r) a_r`.
* For `f` in `N[q]` that is integral, a `C_n`-action making `(X, C_n, f)` a CSP
  exists iff all `a_r >= 0` (`thm:cspCond`).  The natural-coefficient hypothesis
  and the cone condition are reported separately.

The paper's counterexample `f = q^5 + 3q^3 + q + 10`, `n = 6` has natural
coefficients and integral evaluations `15, 8, 12, 5, 12, 8` (`d = 0..5`) but
`S_k = 8, 4, -3, 6` for `k = 1, 2, 3, 6`:

```sh
echo "10, 1, 0, 3, 0, 1" | polytool cyclic-sieving --order 6
#   orbit size 3: -1 orbit(s), -3 element(s)  [a_2 = -1]
#   verdict: no C_n-action makes (X, C_n, f) a CSP
```

To check a concrete action, give its multiset of orbit sizes instead of fixed
counts.  Rotation of 2-subsets of `[4]` has orbits of sizes 4 and 2 and
`[4 choose 2]_q = 1 + q + 2q^2 + q^3 + q^4`:

```sh
echo "1,1,2,1,1" | polytool cyclic-sieving --order 4 --orbit-sizes 4,2
```

With `--json`, each item keeps its `report` object unchanged and gains
`orbit_decomposition` (residues, gcd classes, `integral`, `violating_classes`,
`orbits` with `orbit_count`/`element_count`, `fixed_points`,
`natural_coefficients`, `cone_condition`, `realizable`) and `orbit_sizes_check`
(`null` unless `--orbit-sizes` is given).

Library API (`polytool::`, all `BigInt`): `cyclic_sieving_orbit_decomposition_bigint(coeffs, n)`,
`cyclic_sieving_check_orbit_sizes(coeffs, n, orbit_sizes)`,
`orbit_counts_from_fixed_points(fixed, n)` (Moebius inversion; errors if the
table is not that of an action), `orbit_counts_from_orbit_sizes`,
`fixed_points_from_orbit_sizes`, `fixed_points_from_orbit_counts`,
`reduce_modulo_cyclic_bigint`.

#### Bicyclic sieving

For `G = C_{n1} x C_{n2}` and `f(q, t)` (Barcelo--Reiner--Stanton), the
polynomial is a JSON coefficient matrix on stdin, `f[i][j]` being the
coefficient of `q^i t^j`; entries are integers or decimal strings:

```sh
echo '[[1,0],[0,1]]' | polytool bicyclic-sieving --order1 2 --order2 2 \
  --fixed-table '[[2,0],[0,2]]'
```

The command reduces `f` modulo `(q^{n1} - 1, t^{n2} - 1)`, evaluates it exactly
at all `(zeta_{n1}^{d1}, zeta_{n2}^{d2})` and compares with the table
`fix(g1^{d1} g2^{d2})`.  If `c_{j1,j2} = c_{gcd(j1,n1),gcd(j2,n2)}` the
polynomial lies in the product basis `h_{r1}(q) h_{r2}(t)`; its coordinates
`a_{r1,r2}` (two-dimensional Moebius inversion) are reported together with the
evaluations they imply.  **They are not orbit counts and may be negative**: for
the diagonal action of `C_2 x C_2` on two points, `f = 1 + qt` has coordinates
`a_{1,1}, a_{1,2}, a_{2,1}, a_{2,2} = 1, -1, -1, 2`.  Integrality at all points
is reported separately from product form, because they differ in general when
`gcd(n1, n2) > 1` (for `n1 = n2 = 3`, `sum (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j`
takes only the values 0 and 9 on `mu_3 x mu_3` but is not in the product
span).  Library: `bicyclic_sieving_report_bigint(matrix, n1, n2, fixed_table)`
(limit `n1 n2 <= 2500`).

### Stapledon decomposition

```sh
echo "1, 11, 11, 1" | polytool stapledon 3
```

CLI Stapledon input coefficients may be arbitrary-size integers.

### Benchmark recurrence and interlacing performance

Run the recurrence fixture timing suite:

```sh
polytool bench recurrence-fixtures --repeat 3
polytool bench recurrence-fixtures --only 23_sparse --repeat 5
polytool bench recurrence-fixtures --only oeis --repeat 1
polytool bench recurrence-fixtures --only oeis --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/oeis.md
polytool bench recurrence-fixtures --only oeis --repeat 3 --format json \
  > bench-results/recurrence-fixtures/oeis.json
polytool bench compare old.json new.json --top 10
```

Run consecutive interlacing timings for a generated sequence:

```sh
polytool bench interlacing --sequence eulerian --max-n 30 --repeat 5
```

Run the real-rootedness regression suite on positive-coefficient families of
degree 20--80:

```sh
polytool bench real-rootedness --quick            # small degrees, seconds; used by `cargo test`
polytool bench real-rootedness                    # degrees 20, 40, 60, 80 (a few minutes)
polytool bench real-rootedness --degrees 20,40 --repeat 5
polytool bench real-rootedness --only narayana
polytool bench real-rootedness --format json > rr-baseline.json
polytool bench real-rootedness --baseline rr-baseline.json --max-slowdown 2
```

The suite covers Eulerian, Narayana and type B Eulerian polynomials, h\*-vectors
(the unit cube, products of two simplices, the cross-polytope, each also
round-tripped through the Ehrhart polynomial), and positive-coefficient
counterexamples (a fixed quintic, `1 + t^n`, Eulerian polynomials times
`1 + t + t^2`, `(1+t)^m (1 + 2t + 2t^2)`, and an Eulerian polynomial with a
dented middle coefficient).  Each check has an exact expected verdict known
independently of the implementation: `real_rooted` (exact default test),
`interlacing` (strict `f << g`, including the `none` verdict for degree gaps),
`root_count` (distinct real roots by exact root isolation) and
`hstar_identity` (Ehrhart/h\* conversion round trips).  A wrong verdict, a
verdict that changes between repeats, or a verdict that differs from the
`--baseline` run is reported on stderr and makes the command exit with status 1;
`--max-slowdown <x>` additionally fails checks whose baseline median was at
least 1 ms and became more than `x` times slower.  The output is TSV by
default (`case family kind degree coeff_bits expected verdict status repeat
min_us median_us`, plus `ratio_vs_baseline` with `--baseline`) or JSON
(`polytool.bench.real-rootedness.v1`).  The suite is deterministic; only the
timings vary.  The `--quick` mode is exercised by
`tests/cli_bench_real_rootedness.rs`, so `cargo test` catches verdict
regressions.  Saved timing runs are not checked in; keep them next to the
other results under `bench-results/`.

The recurrence and interlacing benchmark subcommands print tab-separated
output so results can be redirected to a log or pasted into project notes. The
recurrence fixture benchmark can also append fixture/category summaries with
`--summary` and write a Markdown report with `--report <path.md>`. Use
`--format json` for machine-readable per-run records, fixture/category summaries, and adaptive
search diagnostics. The diagnostic records include rejection counters, modular
lift counters, and cumulative timing buckets such as derivative precomputation,
modular prefiltering, modular lifting, exact solving, and held-out
verification. `polytool bench compare` compares two such JSON runs and
highlights fixture/category speedups or regressions. The recurrence fixture
suite contains synthetic stress tests and natural OEIS-derived recurrences from
a curated recurrence queue.

## Library usage

Add to your `Cargo.toml`:

```toml
[dependencies]
polytool = { path = "../polytool" }
```

### Real-rootedness and interlacing

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

### Athanasiadis--Wagner interlacing matrices

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
assert!(check_lace_sequence_total_nonnegative_i64(&polys, 1, 3, 3).is_err());
```

For a `p x q` polynomial matrix `A`, `lace_matrix_i64(A, r, c)` returns the
`p r` by `q c` truncation whose `(p rb+i, q cb+j)` entry is the coefficient
of `x^(cb-rb)` in `A_ij(x)`.  This follows Definition 3.5 of
Athanasiadis--Wagner, *Veronese sections and interlacing matrices of formal
power series*.  A finite TNN check is computational evidence for the infinite
matrix, not a proof of full interlacing unless a separate finite criterion
applies. See [`doc/athanasiadis-wagner-lace.md`](doc/athanasiadis-wagner-lace.md)
for further notes on the construction.

### Planar networks, chain polynomials and strip-digraph recognition

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
match check_lower_unitriangular_tnn(&bad) {
    Err(NetworkError::NotTnn { witness: Some(minor), .. }) => {
        assert_eq!(minor.value, BigInt::from(-2));
    }
    other => panic!("{:?}", other),
}

// Chain polynomials p_0 = 1, p_n = t * sum_{k<n} r_{n,k} p_k.
let ones = to_rows(vec![vec![1], vec![1, 1], vec![1, 1, 1]]);
let report = check_chain_polynomials_bigint(&ones).unwrap();
assert!(report.is_tnn() && report.chain_conclusion_holds());

// Search for a planar strip digraph reproducing the rows (1+t)^n.
let rows: Vec<Polynomial<BigInt>> = (1..=4)
    .map(|n| Polynomial::new((0..=n).map(|k| numeric::binomial_bigint(n, k)).collect()))
    .collect();
let result = recognize_brenti_strip_digraph(&rows, &BrentiRecognitionOptions::default()).unwrap();
assert!(result.is_certified());
```

- `tnn_network::check_lower_unitriangular_tnn` decides total nonnegativity
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

### Polynomial arithmetic

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

### Standard sequences and indexing

`polytool::sequences` generators take `max_n` and return `rows[n]` for
`n = 0..=max_n` (the n-th member of the family, with `rows[0]` the `n = 0`
member), so Eulerian rows are `A_0 = 1, A_1 = 1, A_2 = 1 + t, ...`.  The
OEIS triangles A008292 and A001263 start at `n = 1`, i.e. `rows[1..]`.  The
two exceptions follow their sources' normalization and say so in their
rustdoc: `literature::hoggatt_polynomials_bigint` (`rows[i] = H_{i+1}`) and
the `deco` histogram layers (`rows[N - 1]` for size `N`).  Besides the classical
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

### Characteristic polynomials and closed walks

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

### Recurrence search

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

### Coupled Weyl-matrix recurrences

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

### Other properties

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

### Enumerating real-rooted integer polynomials

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
deterministic counters; no wall-clock times are recorded.  The operation is
not yet exposed through the MCP server.

## MCP server

The repository also contains a local Model Context Protocol server in
`mcp/`. It exposes the exact polynomial routines to MCP clients over stdio.

Build it from the repository root:

```sh
cargo build --release --manifest-path mcp/Cargo.toml --target-dir target
```

Install it locally:

```sh
./mcp/install.sh
```

Check the installed binary:

```sh
polytool-mcp --help
polytool-mcp --version
```

Example MCP client configuration:

```json
{
  "mcpServers": {
    "polytool": {
      "command": "/absolute/path/to/polytool/target/release/polytool-mcp"
    }
  }
}
```

The MCP server is tools-only: no resources, prompts, HTTP server, sampling, or
filesystem access. It returns structured JSON and also includes the same compact
JSON as text content for client compatibility.

Available MCP tools:

- `parse_polynomials`
- `check_polynomial_family`
- `polynomial_properties`
- `check_interlacing_pair`
- `check_interlacing_sequence`
- `check_interlacing_profile`
- `real_roots`
- `find_recurrence`
- `generate_recurrence_rows`
- `resultant`
- `discriminant`
- `ehrhart_hstar`
- `analyze_decomposition`
- `generate_sequence`
- `check_tnn_network`
- `recognize_brenti_strip_digraph`
- `recognize_pf_network`
- `check_left_production_matrix`
- `transform_recurrence`
- `verify_recurrence_transform`
- `transform_recurrence_implicit`
- `verify_recurrence_implicit`
- `assess_generating_description`

`check_tnn_network` and `recognize_brenti_strip_digraph` return the same
fields as the JSON output of `polytool tnn-network` and
`polytool brenti-recognize` (without the `schema` tag), including the fixed
guardrail `note`: TNN of the coefficient matrix, and a strip-digraph
certificate, say nothing about real-rootedness of the row polynomials.
`recognize_pf_network` and `check_left_production_matrix` mirror
`polytool pf-network` and `polytool left-production`; their certificates do
prove PF (real-rooted) rows, and their notes say what they do not prove (at
most 64 rows of degree at most 64 and steps of length and rise at most 8 for
`recognize_pf_network`; at most 32 rows for `check_left_production_matrix`).

The MCP property, interlacing, resultant, discriminant, Ehrhart conversion,
sequence-generation, and recurrence tools accept arbitrary-size integer
coefficients. JSON integers may be provided directly when they fit the client
stack, and larger exact integers should be sent as strings.

`find_recurrence` exposes all adaptive-search controls as top-level MCP
arguments. For example, a compact exact search can use:

```json
{
  "coefficients": [[1], [2], [4], [8], [16]],
  "max_rec_len": 1,
  "max_var_deg": 0,
  "max_idx_deg": 0,
  "max_diff_deg": 0,
  "include_code": false
}
```

The legacy nested `options` object remains accepted, with top-level values
taking precedence. Omitting `include_code` preserves the full generated-code
response.

See [`mcp/README.md`](mcp/README.md) for request schemas, examples, and
development notes.

## Development

Run the core tests:

```sh
cargo test -p polytool
```

Run the MCP server tests:

```sh
cargo test --manifest-path mcp/Cargo.toml
```

Run the web wrapper tests:

```sh
cargo test --manifest-path web/Cargo.toml
node web/tests/string_safety.mjs
```

### Web arbitrary-precision contract and build

The recurrence panel also offers **Compare representations**. It uses the
visible length, coefficient-degree and derivative-order bounds, and a global
configuration limit (default 2000). A table lets the user select a formula,
copy its text/LaTeX/JSON, and extend the sequence with that selected formula.
It displays checked starting indices, held-out row counts, and partial coverage.
The interface preserves the task-oriented UI and row-extension feature deployed
from `878113b`, which had not yet reached the parent integration branch.

The **Full search** mode uses the entered polynomial rows by default. Choose
the maximum number of previous terms and derivative order; coefficient bounds
are in a separate expandable section. The first-row index is explicit. The
browser searches every length/order combination within these bounds, with up
to 16 combinations in the guided form. Examples use defining recurrences for
A171608, Hermite polynomials and automatic polynomial-degree selection.

Results lead with typeset formulas and a copy button, with additional formulas
and search bounds expandable. Short labels distinguish finite-row matches,
identities following from a defining recurrence, bounded absence and unfinished
searches. Certificate exports and linear-algebra diagnostics are kept out of
the webpage; the CLI and library retain their complete reports and replay APIs.

**Advanced: request JSON** supports custom defining recurrences, sparse
supports and resource limits. Editing the JSON switches to custom mode and
disables the guided form; the main input is then unused. Any relevant input or
settings change cancels active work and clears stale results. The existing
WASM export `parse_and_format(input)` supplies exact decimal-string coefficient
rows without computing root properties. The form rejects the entire input if
any row is invalid. Rational rows remain
available through the custom JSON request.
The WASM exports `find_recurrence_spaces(request_json)` and
`verify_recurrence_spaces(report_json)` share the Rust backend.

The browser wrapper parses polynomial input with `parse_polynomials_bigint`.
Every integer coefficient returned across its string-based WASM/JSON boundary
is a canonical decimal string, including input, gamma, decomposition, and
recurrence-initial coefficients. Resultants and discriminants are decimal
strings as well. Browser code must keep these values as strings or convert
them to JavaScript `BigInt` for exact arithmetic; converting them to `Number`
can round values above `2^53`.

Build the no-modules bundle with `wasm-pack`:

```sh
wasm-pack build web --target no-modules --release --out-dir pkg
```

`web/build.sh` runs the same build from inside `web/`. The page `web/index.html`
loads `pkg/polytool_web.js` and `pkg/polytool_web_bg.wasm`; serve the `web/`
directory (for example with `python3 -m http.server`) to try it locally.

## Examples

The tracked programs in `examples/` are run with
`cargo run [--release] --example NAME`:

- `bench_bareiss_vs_modular_linalg`: compares BigInt Bareiss and modular/CRT
  leading principal minors.
- `bench_bezout_vs_sturm`: benchmarks Bézout-matrix and Sturm-chain
  real-rootedness and interlacing tests on Eulerian polynomials.
- `bench_positive_real_rooted`: benchmarks exact real-rootedness checks on
  representative polynomial families.
- `bench_streaming_prs`: reproduces the deterministic corpus used to benchmark
  exact root counting.
- `benchmark_recurrence_fixtures`: benchmarks adaptive recurrence search on the
  fixture suite in `fixtures/recurrence-benchmarks/`.
- `brenti_recognition`: recognizes the Pascal rows as rows of a planar strip
  digraph and checks the chain polynomials of a small TNN matrix.
- `generate_recurrence_benchmarks`: regenerates the recurrence benchmark
  fixtures.
- `literature_sequence_rows`: exports exact literature families (Delannoy,
  Eulerian-square, Hoggatt) as polytool-compatible coefficient rows.
- `recognize_renewal`: reads recurrence JSON on stdin and reports whether it is
  a weighted renewal recurrence.

## Changelog

### Unreleased

Performance:

- The default exact real-rootedness test no longer switches large
  one-signed palindromic inputs to the Uspensky path and counts real roots
  directly (signs at infinity only); degree-79 Eulerian, type-B Eulerian and
  Touchard checks are about 10, 15 and 26 times faster (6.9 s to 0.72 s,
  14.5 s to 0.94 s, 4.9 s to 0.19 s).

Additions (no breaking changes):

- New module `pf_network`: height-shift-invariant planar networks
  (`HeightInvariantNetwork`, `PfStep`, `PfStepType`) with an exact planarity
  test, `build_pf_rows_certificate` / `PfRowsCertificate` (rows plus
  explicitly rebuilt Toeplitz windows), and `recognize_pf_network`, an exact
  column-by-column search over minimal planar supports with step-length,
  rise and budget bounds.  A certificate proves that every matched row is a
  Pólya frequency sequence (Brenti 1995): real-rooted with nonpositive zeros.
  CLI `pf-network`, MCP tool `recognize_pf_network`.
- New module `left_production` (Chen--Fu--Ruan, arXiv:2512.08369):
  `left_production_matrix`, `triangle_from_left_production_matrix`,
  `unified_network_toeplitz_window`, `check_production_tnn`,
  `LeftProductionCertificate::verify` and `check_left_production_matrix`.  A
  TNN left production matrix certifies PF rows and a TNN triangle and
  reversal for the supplied prefix; the report includes a negative-minor
  witness on failure and an explicit planar network on success.  CLI
  `left-production`, MCP tool `check_left_production_matrix`.

- New module `real_rooted_enumeration`:
  `enumerate_real_rooted_integer_polynomials_bigint` enumerates all
  real-rooted integer polynomials of fixed degree with at least three
  prescribed leading coefficients (Greaves--Syatriadi `AllRealRooted`, with
  certified exact endpoints instead of floating point), with node, result and
  coefficient-bit budgets; the kernel
  `real_rooted_antiderivative_integer_constants_bigint` (and
  `..._constant_multiples_bigint`) returns the integer shifts making an
  antiderivative real-rooted.  CLI command `enumerate-real-rooted`.

- `recurrence::find_recurrence_representations_rational`: bounded menus of
  candidate recurrence representations (lags, derivatives, polynomial
  coefficients).
- `recurrence::{find_recurrence_spaces, verify_recurrence_spaces,
  ConstantLagSource, IndexedLagSource}`: complete bounded recurrence spaces
  with exact certificates and replay, automatic right-hand-side degree
  bounds, collective leaders and polynomial-index propagation.
- CLI `--space` / `--verify-space`, matching MCP tools, and a guided
  full-search mode in the web app.
- Recurrence spaces: optional request field `coefficient_rescaling`
  (`recurrence::RecurrenceCoefficientRescaling`, kind `factorial_power`)
  certifies complete spaces for `Q[n] = T_m(P[n])`, `T_m(X^k) = X^k/(k!)^m`,
  from `r(U+m(A+J)+1)` residuals of a constant-index source; results carry a
  replayable `rescaling_bound` (`RecurrenceRescalingBound`).
- `recurrence::{RecurrenceTransformation, DegreeSchedule, transform_recurrence,
  verify_recurrence_transform}`: exact coordinate changes (`x^r`
  multiplication/removal, reversal, affine and Möbius maps, `f`-polynomial,
  gamma coordinates, symmetric decomposition) with inverses, image conditions
  and affine degree schedules on explicit index domains; exact transport of
  derivative recurrences (a vector system for the symmetric decomposition);
  replayable `polytool.recurrence-transform.v1` reports keeping the
  transported rule and an independent fit separately.  CLI
  `recurrence-transform` and MCP tools `transform_recurrence` /
  `verify_recurrence_transform`.  `decomposition::f_polynomial` now uses a
  homogeneous Horner scheme (same results and degree-bound semantics).
- Optional root-domain information for transformed results: request field
  `root_domain` (`recurrence::RootDomainRequest`) adds a `root_domain` report
  section (`RootDomainReport`) with theorem adapters (`root_implication`:
  monomial factors, real affine/Möbius maps, nonnegative-coefficient window,
  gamma quarter bound, Hermite--Kakeya--Obreschkoff symmetric pairs) recording
  direction, hypotheses, degree/sign/endpoint/zero conventions and failing
  converses; evidence-labelled claims; a Rolle recurrence match
  `Q_m = s(m) D^j(g Q_(m-1))` closing inductions (canonical Eulerian rows on
  `[-1, 0]`); and exact root observations.  Replayed by
  `verify_recurrence_transform`; reports without the field are unchanged.
- `real_rootedness::{isolate_real_roots_exact, RealInterval,
  RealRootIsolation, RealRootEnclosure}`: exact rational root enclosures with
  multiplicities, degree and valuation, and exact interval/ray tests;
  `recurrence::root_enclosures` with CLI `root-enclosures` and MCP tool
  `exact_root_enclosures`.  `recurrence::jensen_fixed_cap_check`: fixed-cap
  Jensen-polynomial test of diagonal multipliers (`fixed_cap` evidence).
- Implicit differential relations: `recurrence::{ImplicitRelation,
  CoefficientCondition, GeneratingDescription, assess_generating_description,
  DiagonalWeight, find_reconstruction_certificate, ReconstructionCertificate,
  implicit_transform_recurrence, verify_recurrence_implicit}` and
  `RecurrenceTransformation::{transport_implicit, reconstruction_operator}`.
  New registry maps `diagonal` (`T_n(x^k) = lambda(n,k) x^k`: factorial,
  `W_m`, inverse `W_m`, products of polynomial factors) and
  `derivative_antecedent` (`P_n = D^m(g Q_n)/m!`, exact back substitution).
  Diagonal transport uses one common clearing multiplier with certified
  exceptional-index conditions; uncertified sets stay visible. Generating
  descriptions record seeds, exceptional-index and boundary conditions, and a
  uniqueness argument separate from annihilation. Reconstruction certificates
  `L o T = sum H_nu S^nu M` record their differentiations and are replayed.
  `polytool.recurrence-implicit.v1` reports; CLI `recurrence-implicit`
  (`--verify`, `--assess`); MCP tools `transform_recurrence_implicit`,
  `verify_recurrence_implicit`, `assess_generating_description`.
- `vec_poly::{exact_div_by_t, exact_div_by_t_minus_one, exact_div_by_t_minus,
  sum, prefix_sums, mul_linear, mul_linear_pow, has_nonnegative_coefficients,
  evaluate_homogeneous, evaluate_at_ratio}`: exact division by `t` and
  `t - a` (`None` on a nonzero remainder), sums and exclusive prefix sums,
  products with powers of `a + b t`, a coefficient sign test, and
  division-free evaluation at fractions.
- New module `hurwitz` with `hurwitz_matrix_bigint`,
  `hurwitz_leading_minors_bigint` and `is_hurwitz_stable_bigint` (exact
  Routh--Hurwitz criterion).
- `bezoutian_bigint_coeffs`: Bézoutian for arbitrary degree pairs; agrees
  with `bezout_matrix_bigint_coeffs` when `deg f = deg g + 1`.
- `sequences::shifted_legendre_polynomials_bigint`: shifted Legendre
  polynomials on `[0, 1]` with `P~_n(1) = 1`.
- The `check_weak_interlacing` documentation now spells out the idiom for
  treating a zero polynomial as vacuously interlacing.
- CLI `tnn-network`: TNN verdict for the coefficient triangle of monic rows
  (or, with `--matrix`, of a lower-triangular matrix), with the failing shape
  condition, negative entry or negative-minor witness; canonical network
  weights `lambda` and a reconstruct → evaluate → compare certificate;
  `--chain` adds the Brändén--Saud Maia Leite chain polynomials with exact
  real-rootedness, `[-1, 0]` and interlacing checks.  TNN of the triangle is
  reported separately from, and does not imply, real-rootedness of the rows.
- CLI `brenti-recognize`: exact planar strip-digraph search with
  `--first-row`, `--max-span`, `--weights`, `--max-search-nodes` and
  `--max-integer-nodes`; verdict `certified` (exit 0),
  `impossible_within_bounds` (exit 1) or `inconclusive` (exit 3).  A
  certificate proves TNN of the coefficient matrix, not real-rooted or PF rows.
- MCP tools `check_tnn_network` and `recognize_brenti_strip_digraph` with the
  same structured results (at most 128 rows for `check_tnn_network`; at most
  64 rows of degree at most 64 and bounded search budgets for
  `recognize_brenti_strip_digraph`).
- CLI `bench real-rootedness`: deterministic, verdict-checked regression suite
  (real-rootedness, strict interlacing, distinct-real-root counts and h\*
  conversion round trips) on Eulerian, Narayana, type B Eulerian, h\*-vector
  and positive-coefficient counterexample families of degree 20--80, with a
  `--quick` mode run by `cargo test`, TSV/JSON output and `--baseline`
  comparison.
- Internal: `linalg::sparse_modular` is split into the row type, a private
  `solver` submodule and a `tests` submodule; public paths are unchanged.
- `renewal::{VariablePredecessorRecurrence, PredecessorRule,
  PredecessorWeight, VariablePredecessorCertificate,
  recognize_fixed_delay_predecessor}`: symbolic variable-predecessor
  recurrences `Q_n = Q_(n-1) + w(n) t Q_(p(n))` with exact, replayable
  hypotheses certificates (fixed delay; polynomial weights proven nonnegative
  on all integers `n >= 1` by a Taylor-shift certificate, periodic weights by
  residue) and a canonical row generator; finite lists certify only a prefix.
- `renewal::{match_recurrence_theorems, RecurrenceTheoremReport}`: explicit
  conditional theorem matches (`weighted_renewal_window`,
  `fixed_delay_predecessor_window`) for found or supplied recurrences, schema
  `polytool.recurrence-theorems.v1`, kept separate from finite fit evidence.
- CLI `recurrence --theorem-matches` and `recurrence-theorems --recurrence
  <file|->`; default `recurrence` output is unchanged apart from a stderr note.
  MCP: additive `theorem_matches` field of `find_recurrence` (and the
  `check_polynomial_family` recurrence) and the tool `recurrence_theorems`.

### 0.5.0

Breaking changes:

- `tnn_network::NetworkError::NotTnn` carries an optional explicit negative
  minor (`witness: Option<NegativeMinor>`), and `NetworkError` has the new
  variants `NotSquare` and `NonzeroAboveDiagonal`.
- The `brenti_sequence` documentation no longer claims that a valid strip
  digraph makes its row polynomials real-rooted.  Its transitions depend on the
  height (the digraph is invariant in `x`), so by Lindström--Gessel--Viennot it
  certifies that the coefficient matrix is totally nonnegative; heights
  `0 -> {0, 1}`, `1 -> {2}` give the row `1 + t + t^2`.

Additions:

- `tnn_network::{check_lower_unitriangular_tnn, NegativeMinor,
  lower_triangular_rows_from_square, chain_polynomials_bigint,
  chain_polynomials_rational, check_chain_polynomials_bigint}`: TNN diagnostics
  with explicit failing minors, and the chain polynomials of a
  lower-unitriangular matrix (real-rooted in `[-1, 0]` and interlacing when the
  matrix is TNN; Brändén--Saud Maia Leite, Theorem 3.7).
- `brenti_recognition::recognize_brenti_strip_digraph`: exact search for a
  planar strip digraph reproducing a finite prefix of rows (certificate,
  impossibility within the stated bounds, or inconclusive).
- `BrentiStripDigraph::{from_weight_matrix, weight_matrix}`.
- `exact-real-roots` 0.1.1: `isolate_real_roots_bigint_coeffs` with exact
  dyadic isolating intervals.

### 0.4.0

Breaking changes:

- `sequences::eulerian_polynomials(_bigint)` and
  `sequences::narayana_polynomials(_bigint)` now follow the module-wide
  `rows[n]` convention, include the row `n = 0` (`A_0 = N_0 = [1]`) and return
  `max_n + 1` rows; previously they returned `max_n` rows with
  `rows[i] = A_{i+1}` (resp. `N_{i+1}`).  Migrate with `eulerian_polynomials(n)[1..]`
  to get the old rows.  The CLI `sequence` command and the MCP `generate_sequence`
  inputs still list `A_1..A_n` and `N_1..N_n`.

Additions:

- `sequences::{derangement_descent_polynomials_bigint,
  fixed_point_descent_polynomials_bigint, binomial_eulerian_polynomials_bigint,
  refined_eulerian_polynomials_bigint, s_eulerian_polynomial_bigint}` and an
  "Indexing conventions" section in the `sequences` module documentation.
- `linalg::{characteristic_polynomial_bigint, closed_walk_counts_bigint}`.
- `numeric::{q_integer, q_factorial, q_binomial, q_catalan, q_narayana}`.
- `is_ultra_log_concave_infinite` / `_bigint_coeffs` (ULC(∞)).

### 0.3.0

Breaking changes:

- `vec_poly` is generic over `CoeffRing` instead of fixed to `i64`.
  Scalars are now passed by reference (`scale(&p, &c)`, `evaluate(&p, &x)`),
  every arithmetic result is trimmed (`neg` and `shift` included), and calls
  whose arguments are only integer literals need a type annotation such as
  `&[1i64, 1]`.
- `i64` (and the new `i128`) coefficient arithmetic in `vec_poly` and in
  `Polynomial` arithmetic panics on overflow instead of wrapping in release
  builds; the same holds for `Ratio<i64>`.
- `CoeffRing` gains borrowing arithmetic hooks (`add_ref`, `mul_ref`,
  `add_mul_assign`, ...) with default implementations; a `CoeffRing` impl for
  `i128` is added, which makes integer-literal type inference ambiguous where
  it previously defaulted to `i64`.

Additions:

- `vec_poly::{derivative, reverse, reciprocal, compose, pow, wronskian,
  shift_by, strip_trailing_zeros}` and `Polynomial::pow`.
- `numeric`: `binomial_bigint`, `factorial_bigint`, `falling_factorial_bigint`,
  `multinomial_bigint`, `binomial_u128`.
- `interlaces_either_order` / `interlaces_either_order_bigint_coeffs`
  (undirected, zero-tolerant weak interlacing) and
  `is_ultra_log_concave_of_degree` / `_bigint_coeffs` (ULC at a prescribed
  nominal degree).
- Expanded documentation of the directed and same-degree semantics of
  `check_weak_interlacing`.

## Support

Report issues in the Git repository, or contact Per Alexandersson
(@PerAlexandersson, <per.w.alexandersson@gmail.com>).

## Algorithm notes

### Real-rootedness algorithms

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

### Interlacing matrices

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

### Recognizing weighted renewal recurrences

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

### Variable-predecessor recurrences

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

### Conditional theorem matches

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
