# Development

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

Exact solves stop when a row proves `0 = c` with nonzero `c`; this path returns
no rank. Kernel and certificate computations still request complete reduction.
`benchmark_exact_solving` measures six bounded recurrence searches (including
500- and 1,500-candidate cutoffs) and emits JSON Lines with exact outcome keys.
Build separate release executables before and after a change, alternate their
order, and compare outcome keys and search options as well as timings. A cutoff
is an incomplete search, not evidence that no recurrence exists.

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
