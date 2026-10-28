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
