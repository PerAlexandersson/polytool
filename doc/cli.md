# CLI usage

For complete examples see [CLI workflows](workflows.md). `polytool --help`
lists the commands grouped by task with a "Start here" block, and
`polytool help <command> [<subcommand>]` (or `<command> --help`) prints the
full page of a command: description, options, input, examples, exit status
and deprecated spellings. `sequence --format rows` emits coefficient lists
suitable for piping into the checking and recurrence commands. Invalid options
and input return a nonzero exit status; result fields still distinguish
negative mathematical answers from failed commands.

## Help, options and compatibility

- **Topics.** `polytool help <topic>` explains concepts shared by many
  commands: `getting-started`, `input`, `coefficient-order`, `indexing`,
  `output-formats`, `json-schemas`, `exit-codes`, `budgets`,
  `recurrence-json` and `compatibility`.
- **Machine-readable index.** `polytool help --json` (schema
  `polytool.commands.v1`) lists every command and subcommand with its flags
  (value placeholder, default, choices, deprecated aliases), positionals,
  examples, exit statuses and the topics. The documentation generator and the
  shell completion use the same specification.
- **Completion.** `polytool completions bash|zsh|fish` completes commands,
  subcommands, flags, fixed flag values and `polytool help` topics.
- **Uniform options.** `--flag value` and `--flag=value` are equivalent. Every
  command with structured output accepts `--format <...>` and `--json`
  (`--json` = `--format json`); JSON-only commands accept `--format json` and
  `--compact`. Every command that reads input accepts `--input <file|->`
  instead of standard input.
- **Errors.** A usage error is reported on stderr with exit status 2, suggests
  the nearest flag, command or topic (`did you mean`), and ends with
  ``Run `polytool help <command>` for usage.``
- **Renamed flags.** Since 0.6.0, polytool spells out the abbreviated flags and makes
  the index vocabulary consistent (`--first-index` labels the first input row,
  `--from-index` selects where output starts, `--count` is the number of rows
  to produce). The old spellings keep working and print one warning on stderr
  per use; set `POLYTOOL_DEPRECATION_WARNINGS=0` to silence them. Standard
  output is unchanged. `polytool help compatibility` lists every old -> new
  pair; the main ones are:

  | Old | New | Commands |
  |-----|-----|----------|
  | `--min/max-rec-len` | `--min/max-recurrence-length` | recurrence, oeis add |
  | `--min/max-var-deg` | `--min/max-variable-degree` | recurrence, oeis add |
  | `--min/max-idx-deg` | `--min/max-index-degree` | recurrence, oeis add |
  | `--min/max-diff-deg` | `--min/max-derivative-order` | recurrence, oeis add |
  | `--min/max-inhomo-var-deg`, `--min/max-inhomo-idx-deg` | `--min/max-inhomogeneous-variable-degree`, `--min/max-inhomogeneous-index-degree` | recurrence, oeis add |
  | `--max-denom-var-deg`, `--max-denom-idx-deg` | `--max-denominator-variable-degree`, `--max-denominator-index-degree` | recurrence, oeis add |
  | `--inhomogeneous`, `--denominator`, `--alternating-sign` | `--try-inhomogeneous`, `--try-denominator`, `--try-alternating-sign` | recurrence, oeis add |
  | `--min-margin`, `--no-verify`, `--max-laurent-x-power` | `--min-equation-margin`, `--no-held-out-rows`, `--max-laurent-power` | recurrence, oeis add |
  | `--rows` | `--count` | recurrence-generate, recurrence-degree, oeis generate, oeis evidence |
  | `--additional` | `--additional-count` | recurrence-generate |
  | `--first-row K` | `--from-index K` | recurrence-generate |
  | `--start-row K` | `--from-index K` | oeis generate |
  | `--offset` | `--first-term-index` | recurrence-generate, oeis add |
  | `--reversed` | `--descending` | recurrence-generate |
  | `--first-row <0\|1>` | `--first-index <0\|1>` | brenti-recognize, pf-network, left-production |
  | `--max-integer-nodes`, `--max-span`, `--max-length` | `--max-integer-search-nodes`, `--max-target-span`, `--max-step-length` | brenti-recognize, pf-network |
  | `--leading-desc` | `--leading-coefficients-descending` | enumerate-real-rooted |
  | `--recurrence` (switch) | `--find-recurrence` | family-check |
  | `--offsets` | `--order-offsets` | cyclic-sieving-sequence |
  | `--max-iter`, `--tol`, `--min-imag`, `--grid-re`, `--grid-im` | `--max-iterations`, `--tolerance`, `--min-imaginary-part`, `--real-grid-points`, `--imaginary-grid-points` | bkw-scout |
  | `--max-n` | `--max-index` | bench interlacing |

  JSON output is additive: `oeis generate --json` and the OEIS catalog JSON
  carry `first_index` next to `first_row`, and the planar-network reports
  carry `first_index` next to `first_row_index`.

The `polytool` command-line interface reads polynomials from standard input and
writes text or, with `--json`, machine-readable output.  This page covers the
input format and the general-purpose commands.  The commands for specific areas
are documented on separate pages:

- [OEIS catalog](cli-oeis.md): recurrence-backed OEIS families, catalog data
  and provenance, adding sequences.
- [Planar networks](cli-networks.md): TNN triangles, strip digraphs and Pólya
  frequency rows.
- [Recurrence search](cli-recurrence.md): `polytool recurrence`, several
  representations, complete spaces and certificates, browser verification,
  automatic variable degree, rescaling, machine-readable JSON.
- [Recurrence transformations and implicit relations](cli-recurrence-transforms.md):
  coordinate transformations, root-domain information, implicit differential
  relations, the lossy derivative map, Veronese-section lifts.
- [Searching equivalent representations](cli-recurrence-representations.md):
  the bounded queue search over representations.
- [Ehrhart conversion, cyclic sieving and Stapledon decomposition](cli-ehrhart-sieving.md).

For the Rust API see [library.md](library.md); for the MCP server see
[mcp.md](mcp.md).

## Input format

Polynomials are given as comma-separated integer coefficients in
**ascending degree order**: `a_0, a_1, ..., a_d` represents
`a_0 + a_1 t + ... + a_d t^d`.

Parse errors identify the one-based physical input line, including skipped
comments and blank lines. Rational-list readers also identify the coefficient
position. These source locations are separate from zero-based report row indices.

A text file `polys.txt` might look like:

```
1, 11, 11, 1
1, 26, 66, 26, 1
1, 57, 302, 302, 57, 1
```

Every command with structured output accepts `--json` (or `--format json`)
for machine-readable output; `polytool help <command>` lists its formats.
Large integer coefficients in JSON output are serialized as strings.

The exact property commands `real-rooted`, `interlacing`,
`interlacing-profile`, `properties`, `gamma-expansion`, and the non-recurrence
checks inside `family-check` accept arbitrary-size integer coefficients.
Standard sequence generation, resultants, discriminants, Ehrhart conversion,
Stapledon decomposition, the standalone `recurrence` command, and
`family-check --find-recurrence` also use arbitrary-size exact integers in the
CLI.

## Check real-rootedness

```sh
polytool real-rooted < polys.txt
polytool real-rooted --json < polys.txt
```

The root-only JSON response is `{"items":[{"index":0,"degree":2,
"real_rooted":true}]}`. Indices are zero-based positions among parsed
polynomials. Trailing zero coefficients do not affect the degree; the zero
polynomial has degree `null`, and a nonzero constant has degree `0`.
This command performs only the real-rootedness test, without computing the
additional properties reported by `properties --json`.

`real-rooted` accepts arbitrary-size integer coefficients in dense-list and
expanded-polynomial input.

Output (one line per polynomial):

```
1 + 11t + 11t^2 + t^3: real-rooted
1 + 26t + 66t^2 + 26t^3 + t^4: real-rooted
1 + 57t + 302t^2 + 302t^3 + 57t^4 + t^5: real-rooted
```

## Check interlacing

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

## Check unimodality, log-concavity, palindromicity, gamma-positivity

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

## Expand in the Bernstein basis

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

## Coefficient-only real-rootedness tests

```sh
polytool coefficient-tests < polys.txt
polytool coefficient-tests --json < polys.txt
```

This reports Newton inequalities and Kurtz's sufficient condition.  When the
Kurtz condition holds, the polynomial has distinct real roots.

## Generate standard sequences

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

## Check Family H PF/Jensen pencils

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

## Check a family in one pass

`family-check` reports properties, consecutive weak/strict interlacing, and
optionally an adaptive recurrence search.  Requirement flags make the command
exit nonzero at the first failed requested condition.
Property checks, consecutive interlacing and the recurrence search all accept
arbitrary-size integer coefficients; the gamma expansion of a palindromic row
is exact even when its gamma coefficients exceed `i64`.

```sh
polytool family-check \
  --require-real-rooted \
  --require-weak-interlacing \
  --find-recurrence \
  --json < polys.txt
```

## Scout BKW equal-modulus loci

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
Use `--json` (or `--format json`) for machine-readable output, and
`--mathematica` for a `Reduce` skeleton for exact equal-modulus follow-up.
Together these options put the skeleton in a `mathematica` string field in
the same JSON document. Floating-point options must be finite; `--tolerance` and
`--max-iterations` must be positive, and `--min-imaginary-part` nonnegative. Invalid options
are rejected before reading or scouting the symbol. If numerical overflow
prevents valid JSON, the command fails with a diagnostic and no partial result.

## Compute resultant/discriminant

```sh
# Resultant of two polynomials (given as two lines)
echo "1, 0, 1
1, -1" | polytool resultant

# Discriminant of a single polynomial
echo "1, 0, -3, 1" | polytool discriminant
```

Inputs for these CLI commands may have arbitrary-size integer coefficients.
The exact resultant or discriminant output is printed as an integer.

## Enumerate real-rooted integer polynomials

`enumerate-real-rooted` lists **every** real-rooted integer polynomial
`p(x) = a_0 x^n + a_1 x^(n-1) + ... + a_n` of degree `n` whose first `t`
leading coefficients are prescribed, for `3 <= t <= n + 1` and `a_0 != 0` of
either sign.  Reducible polynomials and polynomials with repeated roots are
included.  The search is the Greaves--Syatriadi derivative-tower algorithm
(arXiv:2504.09241) with certified endpoints (see [Enumerating real-rooted integer polynomials](library.md#enumerating-real-rooted-integer-polynomials));
nothing is read from stdin.

```sh
# All real-rooted x^3 - 3x + c: c = -2, ..., 2
polytool enumerate-real-rooted --degree 3 --leading-coefficients-descending 1,0,-3

# Monic, trace 0, a_2 = -15 (degree 6), as JSON with exact string coefficients
polytool enumerate-real-rooted --degree 6 --leading-coefficients-descending 1,0,-15 --format json
```

The prefix `--leading-coefficients-descending` is in **descending** order (`a_0, a_1, ...`), while
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

## Benchmark recurrence and interlacing performance

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
polytool bench interlacing --sequence eulerian --max-index 30 --repeat 5
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
`tests/integration/cli_bench_real_rootedness.rs`, so `cargo test --features full-tests` catches verdict
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
