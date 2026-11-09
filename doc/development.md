# Development

## Documentation site

Source API documentation uses Rust's `//!` and `///` comments. Build it with
`cargo doc --locked --no-deps`. To combine the guides, executable-derived
command reference and searchable rustdoc for all three core libraries, install
Python 3.10 or later and Pandoc, then run from the repository root:

```sh
python3 scripts/build-docs.py --output /tmp/polytool-docs
python3 -m http.server --directory /tmp/polytool-docs 8000
```

The output directory must be new or empty and outside the source workspace.
The builder uses a fresh external Cargo target directory, checks local guide
links and rejects symbolic links in the site. It renders every `doc/*.md`
guide, gathers command names from `polytool help --json`, and gathers details
from each command's `--help`, including OEIS maintenance subcommands. A new
command is one `CommandSpec` in `src/bin/polytool/commands/` (registered in
`commands.rs`): its flags, help text, examples and exit codes drive the clap
parser, the help pages, `polytool help --json` and completion. The CLI tests
check that every page has usage, examples and exit status, that every
documented flag reaches its handler, and compare the pages with the snapshots
in `tests/snapshots/cli/help` (`POLYTOOL_BLESS_CLI_HELP=1` re-blesses them
after a deliberate change).

`build-docs.py` options: `--target-dir DIR` reuses a Cargo target directory
outside the workspace (its `doc/` is rebuilt each run); `--binary PATH` takes
help output from an existing `polytool` instead of building one;
`--command-source scrape` forces the `--help` scraping path.

The command reference has two sources. When `polytool help --json` describes
commands richly (entries with `usage`, `description`, `options`, `examples` or
`exit_codes`, nested `subcommands`, and an optional top-level `topics` list of
`{name, title, body}`), the builder writes one page per command under
`doc/commands/` and per topic under `doc/topics/`, indexed by the command
reference. If the JSON only lists names and summaries (or the binary has no
`help --json`), it falls back to scraping each command's `--help`. Every
`doc/*.md` guide is included automatically (the Python guide also gets a
navigation entry); guides not named on the front page are listed under "More
guides". Rustdoc for `polytool`, `ore-algebra` and
`exact-real-roots` is published under `/api/`.

## Continuous integration

Three workflows in `.github/workflows/` run on pushes to `main` and on pull
requests (newer pushes cancel running pull-request runs); `python.yml` runs
only when the bindings, the library sources or its own files change.

`ci.yml`:

- `core`: formatting, OEIS catalog data tests, clippy, tests, package check,
  `ore-algebra` tests.
- `companions`: MCP formatting, clippy and tests; `cargo check` of the web crate.
- `msrv`: matrix over `polytool`, `ore-algebra`, `exact-real-roots`, `mcp` and
  `web`; each job reads the package's declared `rust-version` (currently 1.82,
  and 1.89 for `mcp`) and runs `cargo check --lib --bins` on that toolchain.
  A package without a committed `Cargo.lock` (`web`, `mcp`) is first resolved
  by stable Cargo with `CARGO_RESOLVER_INCOMPATIBLE_RUST_VERSIONS=fallback`,
  so dependency releases that need a newer Rust are skipped. Raise an MSRV
  only in the manifest.
- `mcp-release`: `cargo build --release` of `polytool-mcp`, the stdio smoke test
  in release mode, and the binary uploaded as artifact `polytool-mcp-linux-x86_64`.
- `web`: `cargo check --target wasm32-unknown-unknown`, `wasm-pack build
  --target no-modules` (as `web/build.sh`), the CLI, then
  `web/tests/string_safety.mjs` and `node --test web/tests/*.test.cjs` with
  `POLYTOOL_WASM_DIR` and `POLYTOOL_CLI` set, so the real WASM exports and
  browser/CLI parity are tested. The bundle is uploaded as `polytool-web-pkg`.

All Rust jobs use `Swatinem/rust-cache`. The same commands run locally; for
example the web job is `wasm-pack build --release --target no-modules
--out-dir DIR` in `web/` followed by the Node commands above.

`docs.yml` builds the site with `scripts/build-docs.py` (which fails on broken
local links, unresolved guide links and symbolic links), checks that the front
page, command reference and the three rustdoc entry points exist, and uploads
the Pages artifact. Pull requests only build. A push to `main`, or a manual run
on `main` with `publish` left enabled, deploys to GitHub Pages through the
`github-pages` environment. One-time setup: in the repository settings set
Pages source to "GitHub Actions". The site is separate from the WebAssembly
application.

## Checks

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
node --test web/tests/*.test.cjs
```

The Node tests read the page as a browser does: `index.html`, then `app.css`
and `js/*.js` (`POLYTOOL_HTML` may name a staged `index.html`, whose sibling
`app.css` and `js/` are used). With `POLYTOOL_WASM_DIR` naming a built
no-modules bundle, the Node tests also
run the real WASM exports. With `POLYTOOL_CLI` naming a `polytool` binary as
well, `recurrence-transforms.test.cjs` checks that the CLI and the browser
produce identical transform and implicit reports and that each replays the
other's reports.

## Adding a real-rootedness theorem route

The theorem routes of `recurrence::check_real_rootedness_theorem`, of the
proof discovery (`recurrence::discover_recurrence_proofs`) and of the parent
cascade (`recurrence::prove_parent_real_rootedness`) are listed once, in the
private registry `src/recurrence/lift/parent_theorems/registry.rs`
(`ROUTES`, in the order of `RealRootednessTheorem::ALL`).  A new route needs:

1. one module under `src/recurrence/lift/parent_theorems/` with the check at
   the supplied seeds, its public wrapper `check_..._theorem` (a call of
   `check_real_rootedness_theorem`) and its registry entry `ROUTE`: the
   theorem, the check, the start policy (`MovesStart` or `OwnStarts`), the
   preference among proofs at different starts, the certificate name, the
   discovery policy (`Always` or `WhenApplicable`), the residue-class and
   index-section menu flags, and the `ParentPolicy` (stage, requested region,
   whether an attempt that does not apply is recorded);
2. one line in `ROUTES` (the position is the order of the discovery attempts,
   the menus and the cascade stage);
3. the enum variants: `RealRootednessTheorem` (with `ALL` and `as_str`),
   `VeroneseParentProofRoute` (with `theorem` and `of_theorem`) and, when its
   reports carry a certificate, `RouteCertificate` (with `theorem`, the
   replay dispatch and an accessor); the compiler points at every match;
4. the exports in `src/recurrence/lift.rs`, the prose lists (the module
   documentation of `parent_theorems` and `parent_proof`, `doc/library.md`,
   `doc/proof-discovery.md`, the CLI and MCP help), a CHANGELOG entry, and a
   re-blessed fingerprint file whose changed entries are those the route now
   proves or records.

The registry's unit tests and `tests/route_registry.rs` check that every
theorem has one route in the frozen order, that the names are the serialized
names, and that every theorem has its parent route.

## Proof-discovery fingerprints

`tests/discovery_fingerprints.rs` is the evidence gate for refactorings of the
theorem routes. For every bundled OEIS catalog family (860 entries) it hashes
the complete JSON of `recurrence::discover_recurrence_proofs` (default
options) and of `recurrence::prove_parent_real_rootedness` (first index, no
attestations, root-domain route enabled) and compares both with
`tests/fixtures/discovery-report-fingerprints.tsv`; a behaviour-preserving
change leaves every line unchanged, and a change of results re-blesses the
file with every changed entry explained. A sample of 48 cheaper entries
runs in the normal test suite; the ignored full sweep takes about 7 minutes in
release mode (each entry has a limit of 120 s, `POLYTOOL_DISCOVERY_ENTRY_SECONDS`):

```text
cargo test --release --test discovery_fingerprints -- --ignored every_
```

Re-bless only for an intentional change of results:
`POLYTOOL_BLESS_DISCOVERY_FINGERPRINTS=1 cargo test --release --test
discovery_fingerprints -- --ignored every_`.

Integration tests share the helpers in `tests/common/mod.rs` (running the
`polytool` binary with standard input, reading its output, exact-rational
constructors); include them with `mod common;`.

## Ore-algebra workloads

`tests/ore_algebra_workloads.rs` runs a fixed set of real workloads of the
`ore-algebra` crate (image recurrences of the preserver menu shapes,
basis-map elimination, annihilator systems, difference and substitution
image annihilators, right pseudo-division of large dividends) and compares
a fingerprint of each complete output with
`tests/fixtures/ore-algebra-workload-fingerprints.tsv`.  The Hadamard-factor
workloads are compared only by the ignored full test.  Timings (median and
minimum of `ORE_WORKLOAD_REPEAT` runs, filtered by `ORE_WORKLOAD_ONLY`):

```text
cargo test --release --test ore_algebra_workloads -- --ignored --nocapture time_
```

Re-bless the fingerprints only for an intentional change of results:
`POLYTOOL_BLESS_ORE_WORKLOADS=1 cargo test --release --test
ore_algebra_workloads -- --ignored every_`.

## Linear-algebra workloads

`tests/linalg_workloads.rs` runs a fixed set of real workloads of the exact
linear algebra (`linalg`) and of the recurrence fit built on it: adaptive
searches on catalog families (Eulerian, Narayana, type-B Eulerian, shifted
Legendre with `try_laurent`, Hermite with denominators, a sequence without a
recurrence), the Hadamard-factor lifts of the 321- and 123-avoiding
permutation families, `solve_linear_system_with_rank` and
`rational_kernel_with_certificate` on integer, sparse and rational matrices
and on Hankel-type matrices of Eulerian rows.  The kernel workloads also
record the smallest `max_operations` and `max_coefficient_bits` at which the
kernel still succeeds, so the budget accounting is pinned as well as the
result.  The test compares a fingerprint of each complete output (every rule,
counter and certificate, no timings) with
`tests/fixtures/linalg-workload-fingerprints.tsv`; the slower items (the
default-option searches and the lifts) are compared only by the ignored full
test.  Timings (median and minimum of `LINALG_WORKLOAD_REPEAT` runs, filtered
by `LINALG_WORKLOAD_ONLY`):

```text
cargo test --release --test linalg_workloads -- --ignored --nocapture time_
```

Re-bless the fingerprints only for an intentional change of results:
`POLYTOOL_BLESS_LINALG_WORKLOADS=1 cargo test --release --test linalg_workloads
-- --ignored every_`.

To profile, `valgrind --tool=callgrind` on the test binary (with
`LINALG_WORKLOAD_ONLY` set to one workload) works where `perf` does not.  The
elimination used to spend about 80% of these workloads in the gcd reductions of
`num_rational`; `src/linalg/integer_elimination.rs` now does it on integer rows
(see its module documentation), and the remaining cost is `num_bigint`
multiplication and exact division.

The integer update now combines delayed target-row scaling with elimination and
sets the eliminated pivot-column entry to zero without computing cancelling
products. The integer square solver uses the same integer back substitution as
the general solver, reducing fractions only for final solution coordinates.
These changes preserve the integer update results and rational certificate budgets.

Rational square systems now clear denominators once per augmented row using
the same helper as general integer elimination. Integer-only rows retain their
direct conversion path, and rational elimination remains an exact fallback.
`benchmark_rational_square --size 48 --repeat 5` measures integer, row-scaled,
mixed-denominator, fractional-RHS and Hilbert systems with known exact solutions.
Setup and solution comparison are outside the timer; denominator clearing is
included. Use separate binaries and alternate their order for comparisons.
The `singular_early` and `singular_late` cases check rejection near the first
and last pivot; their exact expected result is `None`. Proven singularity
returns immediately, while a non-exact division still requests the fallback.
The `bidiagonal` and `block_diagonal` cases measure sparse nonsingular systems.
Square elimination now also keeps a level for each row, skipping rows with a
zero pivot-column entry. A selected pivot row is brought to the current level;
target updates fuse scaling and elimination, dividing by the target's own
level. The existing integer rescaling helper checks every exact division.

Exact solves stop when a row proves `0 = c` with nonzero `c`; this path returns
no rank. Kernel and certificate computations still request complete reduction.
`benchmark_exact_solving` measures six bounded recurrence searches (including
500- and 1,500-candidate cutoffs) and emits JSON Lines with exact outcome keys.
Build separate release executables before and after a change, alternate their
order, and compare outcome keys and search options as well as timings. A cutoff
is an incomplete search, not evidence that no recurrence exists.

For integral source rows and derivatives, exact recurrence fitting constructs
the matrix with `BigInt` entries and feeds the existing integer elimination
directly. Fractional sources retain rational construction and row denominator
clearing. Both paths use the same column layout, normalization charts and
verification; integral input does not require integral solution coefficients.

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
- `benchmark_exact_solving`: benchmarks bounded recurrence searches;
  `--repeat N` and `--only SUBSTRING` select repetitions and cases.
- `benchmark_rational_square`: benchmarks exact square solving with integer and
  fractional inputs; `--size N`, `--repeat N` and `--only SUBSTRING` select cases.
- `brenti_recognition`: recognizes the Pascal rows as rows of a planar strip
  digraph and checks the chain polynomials of a small TNN matrix.
- `generate_recurrence_benchmarks OUT_DIR`: regenerates the synthetic
  recurrence benchmark fixtures into a scratch directory.
- `lazy_rows`: pulls rows lazily from the OEIS catalog, a plain recurrence and a
  named family (`.by_ref().take(k)`, continuing the same iterator, and the
  equivalent batch form).
- `literature_sequence_rows`: exports exact literature families (Delannoy,
  Eulerian-square, Hoggatt) as polytool-compatible coefficient rows.
- `recognize_renewal`: reads recurrence JSON on stdin and reports whether it is
  a weighted renewal recurrence.

See also [web.md](web.md) for the web wrapper and its build.
