# Rust repository health audit — 2026-09-16

Status: proposal only. The deliberately ignored experiment forest remains
local and disposable; this note does not propose tracking or cataloguing it.

## Immediate correctness fixes

1. Replace the `u32 as i32` arithmetic in the two skew-Schur Jacobi--Trudi
   determinants in `sym-poly/sym/src/symmetric_function.rs` with checked
   `i64` arithmetic. Parts above `i32::MAX` currently produce a wrong term or
   overflow. Add one-row boundary tests for the h- and e-determinants.
2. Give `sym-poly/core/src/crt.rs` an honest overflow contract. Its public
   `i128` CRT and rational-reconstruction arithmetic wraps in release when
   valid moduli or intermediate products exceed `i128`; modular Gröbner
   lifting reaches this code. Prefer a `BigInt` core, or return an explicit
   overflow error from checked arithmetic, with near-limit tests.
3. Use `checked_add` in
   `sym-poly/multipoly/src/indexed_variables.rs::monomial_multidegree`.
   Multiple exponents in one alphabet can silently wrap a `u32` degree in
   release builds.
4. Define or reject the zero-index case for `IndexedVariables`. With zero
   indices, graded-quotient weak-composition enumeration has no base case and
   recurses until stack overflow. Test the public component API at degree zero
   and positive degree.
5. Check degree multiplication in
   `SymmetricFunction::plethysm_power_sum`; `k * part` currently wraps in
   release. Document the representable-degree limit.
6. Reconcile `Composition` and `WeakComposition`. `Composition::new` preserves
   interior zeros, weak-composition enumeration constructs such values, while
   ribbon and descent-set conversions assume positive parts and may underflow.
   Either enforce strong compositions or make preconditions/errors explicit.

## Repository guardrails

1. Separate durable integration checks from the local scratch crate. This
   workstation has 971 auto-discovered `experiments/src/bin` targets, only ten
   of which are tracked. Cargo still sees all ignored files, so
   `cargo test --workspace --all-targets --locked` currently fails on the
   ignored `uig_path_recurrence_probe.rs`, which uses a removed API. A clean
   clone and this workstation therefore have different target graphs.

   Preferred design: a small tracked integration-check crate with explicit
   `[[bin]]` targets, plus an ignored standalone scratch crate outside the root
   workspace. A smaller interim step is safe `default-members` and CI commands
   that explicitly exclude experiments.
2. Add root CI. The only workflow is `polytool/.github/workflows/ci.yml`, which
   becomes active only in the generated standalone Polytool branch. Use a
   package matrix for maintained foundation, Sym/QSym/multipoly, Combpoly,
   Polynomial Lab, and selected integration checks; keep retired/database and
   scratch workloads separate.
3. Resolve the Combpoly dual-history layout. The parent repository tracks all
   `combpoly/` files as ordinary files, while `combpoly/.git` maintains a
   second identical repository with no remote. Choose one source of truth:
   remove the inner metadata after preserving history, or convert it into a
   documented submodule/subtree. `kostka` is already a proper submodule;
   `ehrcalc` is an intentionally ignored independent repository.
4. Mark internal packages `publish = false`. For intended public packages,
   add explicit license/readme/repository metadata and choose/test an MSRV via
   `rust-version`. At present all 15 workspace packages have implicit publish
   eligibility and none declares an MSRV.
5. Remove the stale `combinatoric-core/Cargo.lock`; the root workspace ignores
   it and it has not been updated since the initial import. Retain lockfiles
   only at actual standalone workspace boundaries.
6. Make generated OEIS catalog verification reproducible. The 47,605-line
   generated Rust file has a generator, but the default generator inputs are
   absolute local OEIS, project, and Lean paths, and standalone CI does not run
   `--check`. Check in a normalized source snapshot/checksum mode or provide a
   self-contained drift check.
7. Slim the public Polytool package and remove local target surprises. The
   package contains 255 files and 8.0 MiB uncompressed, including 5.8 MiB of
   recurrence fixtures plus internal handoff/CI material. Also, the ignored
   `polytool/examples/ferrers_nofactor_interlace.rs` is auto-discovered locally
   and causes a packaging warning. Move that probe to the scratch area and use
   `include`/`exclude` metadata for release contents.
8. Archive old chronology from the root and Polytool handoffs. Each is about
   700 lines/39 KiB; keep a short current opening and move completed history to
   an archive file.
9. Track the MySQL dependency's future-incompatibility warning. It enters only
   through `experiments` and legacy `ktt-search` via `proc-macro-error2 2.0.1`.
   Isolating those crates keeps maintained-library checks clean while the
   dependency is upgraded upstream.

## Performance work, after correctness and CI

1. Stream graph matching/independent-set counts instead of materializing every
   object before computing `matching_polynomial`, `independence_polynomial`,
   and chordal sink polynomials. Add pruned traversals for perfect,
   noncrossing, and nonnesting matchings.
2. Construct `Graph::line_graph` from incident-edge lists. This improves
   `O(m^2)` pairwise comparison to `O(sum_v degree(v)^2)`, linear for bounded
   degree graphs.
3. Store transition matrices behind `Arc` in the Sym and multipoly caches;
   warm cache hits currently clone the full dense matrix.
4. Stream tableau weights and P-partition descent compositions rather than
   materializing every tableau or linear extension. Likewise, make
   `Poset::num_linear_extensions` count leaves directly.
5. Stream Combpoly's `poly` and recurrence CLI paths. They currently retain
   the entire permutation family and then a second vector of statistic values.
6. Replace `polytool::vec_poly::permanent` Laplace expansion with subset DP,
   and define rectangular/ragged input behavior.
7. Accumulate QSym quasi-shuffle multiplicities during recursion instead of
   materializing duplicate paths.
8. Add targeted benchmarks before optimizing: line graphs, graph polynomial
   enumeration, cold/warm basis conversion, flagged tableaux, antichain
   P-partitions, Combpoly streaming, and polynomial permanents.

## Audit evidence

- `cargo fmt --all -- --check` passed.
- The standalone Polytool export is byte-identical to `origin/main` and passes
  locked standalone metadata without sibling dependencies.
- Maintained package tests completed successfully during an
  experiments-excluded workspace run, but the long recurrence-fixture suite
  was interrupted; this is not a claim of a complete workspace test pass.
- No `unsafe` blocks, `todo!`, or `unimplemented!` branches were found in the
  audited maintained libraries.
