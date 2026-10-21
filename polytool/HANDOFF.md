# Polytool Handoff

## Verified: automatic bounds and polynomial-index propagation — 2026-10-21

The derivative-lift lane completed the three authorized scoped extensions on
`feat/recurrence-representations-20261020`, based on 578654d. Ownership is
released after this checkpoint. Owned files: `polytool/src/recurrence/space.rs`,
new `propagation.rs`, recurrence-space integration/CLI tests and example,
MCP/web adapters and tests, browser HTML, parent/MCP READMEs, library guide
and root/crate handoff openings. Shared Rust, Lean/OEIS and shared research
checkouts remain untouched. No workers, merge, subtree publication,
installed-MCP replacement or live deployment occurred.

`automatic_variable_degree` computes a sufficient RHS X-degree preserving
all allowed leaders; index, derivative and leader bounds stay explicit.
`include_leader_coverage` combines identities with leader multipliers to
obtain a sum-of-squares leader, preserving every regular rank within the
reported space, at possibly higher coefficient degrees. `indexed_source`
adds monic Q[n,X] pure-lag sources, a canonical exact generator, and explicit
common residual propagators. Replay verifies the polynomial identity and
cutoff, using exact exceptional ranks or a conservative Cauchy bound with
roots marked unknown. Old APIs/schema behavior is preserved. Resource
exhaustion is incomplete; automatic caps can exceed the configured budget.
The installed Polytool MCP predates this API; exact checks use the owned
crate/CLI, while polynomial-lab remains the proof-search ledger.

All suites passed: 360 library tests, existing integration/fixture suites,
25 recurrence-space and three CLI tests, 34 MCP library plus binary/stdio,
13 native web and five doctests. After incremental coefficient-storage
limits, the 25+3 focused tests passed again. Strict all-target Clippy,
strict Rustdoc, touched-file formatting, standalone metadata, bundled OEIS
checksum, package boundaries and whitespace passed. Hermite and automatic
example reports both replay through the CLI. All 50 Node checks pass against
optimized WASM, plus string safety. No visual-browser check is claimed.
Bundle: `/cargo-target/polytool-recurrence-continuation-web`; logs:
`/cargo-target/recurrence-continuation-*`. wasm-pack retains its existing
missing license-file packaging notice.

The browser adds Hermite/automatic examples and displays the degree scope,
propagation cutoff and combined identity. RealRooted #1978 and the isolated
research note record the proofs and future Lean lemmas. General index-degree,
spectral/localization and unrestricted-order backends, nonmonic singular
sources and Lean tactic integration remain future work. No survey expansion
or new real-rootedness theorem is claimed.

## Verified: complete recurrence spaces — 2026-10-21

The isolated recurrence-representations lane adds `recurrence --space` and
`--verify-space`, shared library/MCP/WASM APIs, and a **Complete spaces** browser
view. Reports contain full bounded kernels, successor projections, source-
conditional certificates and finite-row absence duals. Original indices,
source consistency, singular ranks and resource shortages remain explicit.
Legacy finite-fit APIs are unchanged. Exact scope and all passing checks are
in the parent handoff; ownership is released. No installed MCP replacement,
shared-checkout edit, merge, subtree publication or live deployment occurred.
RealRooted #1978 tracks the remaining Lean checker/tactic work. Spectral and
unrestricted-order algorithms remain future stages.

## Verified recurrence representation checkpoint — 2026-10-20

The isolated `feat/recurrence-representations-20261020` lane adds the library,
CLI `--representations`, MCP tool and browser comparison view. Source ownership
is released after the checkpoint; exact scope/checks are in the root handoff.
The shared checkout and installed MCP remain untouched. All outputs are finite
exact fits with explicit domains and coverage; complete-kernel/all-rank source
certification remains future work. The optimized browser bundle passes all
48 Node checks. Website staging is separate; no live deployment or merge.

Full chronological history through 2026-09-16 is archived in
[`doc/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md`](doc/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md).

## Current state

`master` in the parent Rust repository is canonical. The public `main` branch
is a generated standalone subtree and must be updated only with
`../scripts/sync-polytool-main.sh` after the parent checkpoint is committed,
pushed, clean for eligible paths, and verified.

The repository-health batch owns Polytool manifest/package metadata, the
self-contained generated-OEIS checksum check, and the completed subset-DP
permanent optimization in `src/vec_poly.rs`. Source ownership is released. No
release or deployment has been authorized.

## Verification expectations

- Keep exact arithmetic and explicit overflow/shape contracts.
- Run focused library/CLI tests before broader checks.
- Run `python3 scripts/build_oeis_catalog.py --check-bundled` without requiring
  external OEIS or Lean sources.
- Verify `cargo package --list` excludes internal agent notes, handoffs, CI,
  benchmark output, and generator tooling while retaining the recurrence
  fixtures used by the installed benchmark CLI.
- Keep targets outside the source checkout and preserve unrelated changes.
