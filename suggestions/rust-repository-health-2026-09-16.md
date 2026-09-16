# Rust repository health audit — 2026-09-16

Status: implemented on 2026-09-16, except for the explicitly retained follow-up
items below. The deliberately ignored experiment forest remains local and
disposable; it was neither tracked nor catalogued.

## Completed correctness fixes

1. Replaced the `u32 as i32` arithmetic in the two skew-Schur Jacobi--Trudi
   determinants in `sym-poly/sym/src/symmetric_function.rs` with checked
   `i64` arithmetic, with one-row boundary tests for the h- and e-determinants.
2. Gave `sym-poly/core/src/crt.rs` an honest overflow contract. Its public
   `i128` API now uses `BigInt` intermediates and returns an explicit overflow
   error when a result cannot be represented, with near-limit tests.
3. Added `checked_add` in
   `sym-poly/multipoly/src/indexed_variables.rs::monomial_multidegree`.
   Multiple exponents in one alphabet can no longer silently wrap a `u32`
   degree in release builds.
4. Defined the zero-index case for `IndexedVariables`. With zero
   indices, degree zero has the single empty monomial and positive degrees
   yield no monomials; both public component cases are tested.
5. Checked degree multiplication in
   `SymmetricFunction::plethysm_power_sum`, with the representable-degree limit
   enforced rather than wrapping in release.
6. Reconciled `Composition` and `WeakComposition` with checked strong- and
   weak-composition APIs, so ribbon and descent-set conversions no longer
   silently receive invalid zero parts.

These fixes are in pushed commit `4a897bf`; the CI-only CRT test helper repair
is in `f7d136f`.

## Completed repository guardrails

1. Separated durable integration checks from the local scratch crate. This
   workstation has 961 ignored `experiments/src/bin` targets and ten tracked
   targets. The standalone experiments workspace now declares its durable
   targets explicitly, so local scratch files cannot alter the CI target graph.
2. Added root CI. The only prior workflow was
   `polytool/.github/workflows/ci.yml`, which
   becomes active only in the generated standalone Polytool branch. Use a
   package matrix covers maintained foundation, Sym/QSym/multipoly, Combpoly,
   Polynomial Lab, and selected integration checks; retired/database and
   scratch workloads remain separate.
3. Resolved the Combpoly dual-history layout. The parent repository tracks all
   `combpoly/` files as ordinary files. The inner Git history was preserved in
   a verified bundle and metadata backup before its nested `.git` was removed.
4. Marked internal packages `publish = false`; intended public packages now
   declare license/readme/repository metadata and a tested `rust-version`.
5. Removed the stale `combinatoric-core/Cargo.lock`; the root workspace ignores
   it. Lockfiles now remain only at actual standalone workspace boundaries.
6. Made generated OEIS catalog verification reproducible. The 47,605-line
   generated Rust catalog now has a self-contained checksum drift check used by
   CI, without depending on absolute local OEIS, project, or Lean paths.
7. Slimmed the public Polytool package and removed local target surprises. The
   release excludes internal handoff/CI material while retaining the recurrence
   fixtures required by the installed CLI. The ignored example probe moved to
   the scratch workspace and no longer affects package discovery.
8. Archived old chronology from the root and Polytool handoffs. Each was about
   700 lines/39 KiB; each now has a short current opening and an archive file.
9. Isolated and documented the MySQL dependency's future-incompatibility
   warning. It enters only
   through `experiments` and legacy `ktt-search` via `proc-macro-error2 2.0.1`.
   Isolating those crates keeps maintained-library checks clean while the
   dependency is upgraded upstream.

These guardrails are in pushed commit `e777c57`. The root workflow now also
checks the declared Polytool/web and MCP MSRVs. The deliberately local 961-file
scratch forest remains ignored; CI checks only the ten tracked experiment
targets from a clean checkout.

## Completed performance work

1. Streamed graph matching/independent-set counts instead of materializing every
   object before computing `matching_polynomial`, `independence_polynomial`,
   and chordal sink polynomials. Perfect, noncrossing, and nonnesting matchings
   use pruned traversals.
2. Constructed `Graph::line_graph` from incident-edge lists. This improves
   `O(m^2)` pairwise comparison to `O(sum_v degree(v)^2)`, linear for bounded
   degree graphs.
3. Stored transition matrices behind `Arc` in the Sym and multipoly caches;
   warm cache hits no longer clone the full dense matrix.
4. Streamed flagged-tableau weights and P-partition descent compositions rather
   than materializing every tableau or linear extension.
   `Poset::num_linear_extensions` now counts leaves directly.
5. Streamed Combpoly's `poly` and recurrence CLI paths. They previously retained
   the entire permutation family and then a second vector of statistic values.
6. Replaced `polytool::vec_poly::permanent` Laplace expansion with subset DP,
   with defined rectangular/ragged input behavior.
7. Accumulated QSym quasi-shuffle multiplicities during recursion instead of
   materializing duplicate paths.

## Retained follow-up

1. Add targeted benchmarks for line graphs, graph polynomial
   enumeration, cold/warm basis conversion, flagged tableaux, antichain
   P-partitions, Combpoly streaming, and polynomial permanents. The optimized
   paths have exhaustive or oracle-equivalence tests, but durable performance
   baselines still deserve their own focused change.
2. Pay down the 38 pre-existing strict-Clippy findings in
   `combinatoric-core`. Keep that cleanup separate from behavior changes.
3. Upgrade the experiments/KTT MySQL dependency when its upstream
   `proc-macro-error2 2.0.1` future-incompatibility warning is resolved.

## Audit evidence

- `cargo fmt --all -- --check` passed for the completed checkpoints.
- `cargo doc --locked --workspace --no-deps` passes with
  `RUSTDOCFLAGS="-D warnings"`; root CI now keeps rustdoc links warning-free.
- The standalone Polytool export is byte-identical to `origin/main` and passes
  locked standalone metadata without sibling dependencies.
- Every maintained package passed its focused suite; the long Polytool library
  suite passed locally. Root CI separately covers each maintained package,
  tracked experiments, and legacy Kostka/KTT.
- No `unsafe` blocks, `todo!`, or `unimplemented!` branches were found in the
  audited maintained libraries.
