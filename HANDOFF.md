# Handoff

## Active Polytool GitHub issues #1 and #2 (2026-09-06)

The coding worker owns the Polytool issue implementation in the isolated
worktree `/tmp/polytool-issues-1-2-20260906` on branch
`fix/polytool-issues-1-2-20260906`, based on freshly fetched
`origin/master` commit `6486f93`.  The shared `/workspace/rust` checkout and
its four local-only commits remain untouched.

Ownership is limited to `HANDOFF.md`, `polytool/HANDOFF.md`,
`polytool/build.rs`, `polytool/src/version.rs`, `polytool/src/lib.rs`,
`polytool/src/recurrence.rs`, `polytool/src/bin/polytool.rs`, focused new or
existing tests under `polytool/tests/`, `polytool/README.md`,
`polytool/mcp/src/lib.rs`, and `polytool/mcp/README.md`.  No other worker
currently claims these files.  The separately owned derangement experiment
and SymCat weighted-bond example will not be touched.

The task is to implement deterministic CLI build-version reporting and an
exact candidate budget for adaptive recurrence finding, keep library/CLI/MCP
outcomes aligned, verify the canonical branch, publish it without including
the shared checkout's local-only history, project `polytool/` to standalone
`main`, and close GitHub issues #1 and #2 only after green CI.

Implementation checkpoints `92f7612` and `87ef7f8` add the reproducible
12-hex-digit build commit with an honest `git unavailable` fallback, plus an
outer-candidate recurrence budget and distinct found/search-space-exhausted/
budget-exhausted library, CLI, and MCP outcomes.  Focused zero, small,
exact-boundary, success, fallback, and unbounded regressions pass.

Established verification is green with external
`CARGO_TARGET_DIR=/cargo-target/ai-projects`, `timeout 60s`, and `nice -n 10`:
318 non-OEIS library tests; 61 recurrence tests; 19 CLI BigInt tests; 7 CLI
OEIS tests; 6 version/budget CLI tests; 2 overfit fixtures; 5 interlacing API
tests; 22 MCP library tests plus its binary and documentation targets; 5
Polytool doctests; the focused imported-OEIS validation replay; strict
Polytool/MCP Clippy; Cargo metadata; formatting; and `git diff --check`.

Branch publication is at PR #3.  A local standalone subtree preflight at
split commit `a7bcd0a` passes all six version/budget CLI tests and reports
`polytool 0.2.1-rc.5 (git a7bcd0abf054)`.  GitHub reports no checks because
the repository contains no `.github` workflow and `master` has no required
status checks.  Per the user's green-CI condition, PR #3 has not been merged,
standalone `main` has not been updated, and issues #1/#2 remain open pending
authorization to treat the comprehensive local checks as the merge gate or
to add CI outside the task's allowed file scope.

## Final Polytool integration (2026-09-06)

The final integration worker merged clean local `master` commit `29c075e` into
`fix/polytool-review-20260906` at merge commit `a801eb1`. Its parents are
`41d57bd` (all review fixes plus the earlier `origin/master` merge) and
`29c075e` (all monorepo audit checkpoints), so both histories are preserved.
No textual conflict occurred: the Polytool handoff and root handoff are distinct
files, and the other LLT work is under `sym-poly/sym`.

Integrated verification used external
`CARGO_TARGET_DIR=/cargo-target/ai-projects`, `timeout 60s`, and `nice -n 10`.
The Polytool 312-test non-OEIS library suite, 57 recurrence tests, 62 linear
algebra tests, focused OEIS validation replay, CLI/fixture/interlacing tests,
documentation tests, all Polytool MCP targets, and strict Clippy passed. The
changed Combpoly, Multipoly, Sym, Kostka, and experiments packages also passed
their focused suites; all five tracked experiment binaries compiled. Root
Cargo metadata and `git diff --check` passed.

This handoff accompanies the user-authorized publication of the integrated
branch to monorepo `master` and the subsequent standalone Polytool `main`
projection. Exact final remote commit identities are recorded in this root
handoff after publication. Integration ownership is released after the final
verification; Ehrcalc remains untouched.

Publication completed successfully. The integrated code and Polytool handoff
were pushed non-forced to monorepo `origin/master` at `93b3971`, and canonical
`/workspace/rust` was fast-forwarded to the same commit without changing any
ignored experiment or nested Ehrcalc file. Running
`scripts/sync-polytool-main.sh` projected `polytool/` to standalone
`origin/main` commit `fed0702`. This root-only publication record does not alter
the projected subtree; current local and remote refs are verified again after
this handoff checkpoint is pushed. No force push was used.

## Completed monorepo integration audit (2026-09-06)

The user-authorized audit adopted the pre-existing dirty checkout on local
`master` at `d0a1ca3`. Ownership is now released. The focused checkpoints are:

- `cd34ef3` records ownership and provenance;
- `3f755f8` synchronizes the clean nested `combpoly` commits `b42aecf`,
  `c6c13c3`, and `56bc239` into the parent snapshot;
- `33c2fb0` adds connective-K Grothendieck/Lascoux operators and expansions;
- `65a1cfc` preserves three reproducible multipoly source probes;
- `19c5435` adds the LLT/chromatic Dyck recursions and Petrie functions;
- `5297104` advances only the `kostka` gitlink to clean published
  `origin/main` commit `2481c97`;
- `dd47899` restores the `experiments` manifest and stable shared helpers so a
  clean workspace loads and its five tracked binaries compile, while retaining
  the disposable experiment tree as ignored material.

The standalone `ehrcalc/` directory was deliberately excluded and is now
ignored by the parent repository. It has its own `.git`, upstream, guide, dirty
formatting changes, and an untracked `ktt-search/` tree containing generated
JSON reports. Those nested-repository files and results remain untouched and
must be reviewed and committed in Ehrcalc itself. Other generated logs/results,
all disposable experiment binaries, the actively owned
`derangement_index_parity_recurrence.rs`, and the separately owned SymCat
example also remain untouched and untracked.

The supervisor ownership check was attempted, but `supervisor-tool report`
could not run because this container has no Docker executable. No handoff entry
claimed active ownership of the adopted files. The audit used external
`CARGO_TARGET_DIR=/cargo-target/ai-projects`, `timeout 60s`, and `nice -n 10`.
Verification completed as follows:

- `cargo test -q -p combpoly`: 194 passed;
- `cargo test -q -p sym-poly-multipoly --lib`: 198 passed;
- all three multipoly probes ran successfully, and all multipoly examples
  compiled;
- `cargo test -q -p sym-poly-sym`: all unit, integration, example, and doctest
  targets passed;
- `cargo test -q -p kostka`: 34 passed;
- `cargo test -q -p experiments --lib`: 19 passed, strict library Clippy
  passed, and all five tracked experiment binaries compiled;
- a detached clean worktree initialized the Kostka submodule, loaded the root
  manifest, and compiled all five tracked experiment binaries;
- the broad `cargo test --workspace --lib` attempt passed every preceding
  package and reached 307 of 308 Polytool tests before the 60-second guard;
  the focused owned-package suites above completed independently.

Strict Clippy passed on the adopted multipoly, sym, and experiments code after
allowing only known warnings in pre-existing files. The combpoly strict-Clippy
attempt was blocked by unrelated existing warnings in `catalan.rs` and
`lattice_path_matroid.rs`; its full test suite passed. Nothing was pushed.

## Circular LLT highest-reachable-vertex expansion

Commit `39fdd41` adds the ordinary Mathematica-compatible HRV, the lifted
circular HRV, and the circular vertical-strip orientation formula
`sum q^asc e_lambda` to `sym-poly/sym/src/llt.rs`.  The independent integration
test `sym-poly/sym/tests/circular_hrv.rs` compares the formula with direct LLT
coloring enumeration and an independent monomial-to-elementary conversion for
every circular area sequence through rank 4, plus a strict-edge example.  All
four original integration tests and all 21 LLT unit tests pass.  Follow-up
commit `e5fa03c` adds four genuinely circular rank-8 unicellular examples and
four genuinely circular rank-8 vertical-strip examples with 2--3 admissible
strict corner edges.  Commit `9477786` adds the exact rank-8 example used in
the LLT-flip manuscript: it checks the displayed HRV fiber, all five
coefficients in the full elementary expansion, and equality with direct LLT
coloring enumeration.  All seven integration tests pass.  The pre-existing
uncommitted Dyck-recursion changes were preserved, and no active ownership
remains from this task.

Verification:

```text
CARGO_TARGET_DIR=/mnt/2TB-Babel/ai-storage/cargo-target timeout 60s nice -n 10 \
  cargo test -p sym-poly-sym --test circular_hrv -- --nocapture
```

## Matroid contingency-array scan

The completed matroid contingency-array scan used
`experiments/src/matroids.rs` and the disposable binary
`experiments/src/bin/matroid_contingency_scan.rs`.  The reusable increment adds
graphic and arbitrary-transversal constructors and an exact common-column-sum
polynomial method to `BasisMatroid`; the binary uses these to test small
graphic and transversal matroids.  All 17 focused `matroids` library tests
pass.  The completed search ranges and exact commands are recorded in
`../projects/matroid-contingency-arrays/notes/initial-graphic-transversal-scan.md`.
The reusable `experiments/src/matroids.rs` helper is now tracked so the workspace
member builds from a clean checkout; the scan binary remains intentionally
ignored/disposable. No active file ownership remains from this task.

The same scanner now also has a `catalan` mode, and `BasisMatroid::catalan`
constructs the standard Dyck-path transversal matroid.  The completed Catalan
ranges, exact boundary coefficient, and OEIS matches are in
`../projects/matroid-contingency-arrays/notes/catalan-matroids-and-oeis.md`.
The follow-up proof search adds `BasisMatroid::column_sum_offset_polynomial`,
the scanner flags `--check-row-interlacing` and
`--check-imbalance-packet`, a fixed-Catalan-order row scan, and an
orbit-compressed Catalan orthant recurrence.  All 19 focused matroid tests
pass.  The exact recurrence, completed class-level checks, and next proof
lemma are recorded in
`../projects/matroid-contingency-arrays/notes/catalan-core-proof-search.md`.
The shared matroid helper is tracked; the scanner and its generated output
remain intentionally ignored/disposable. No active file ownership remains from
this task.

## Derangement index-parity recurrence search

The main Rust worker owns only the new disposable experiment
`experiments/src/bin/derangement_index_parity_recurrence.rs` for the current
task.  It uses the new exact vector Weyl-recurrence API to test the split
`E_m=P_(2m)`, `O_m=P_(2m+1)` for normalized reciprocal derangement descent
polynomials.  Existing experiment binaries and all unrelated dirty files
remain untouched.

The probe verifies both the affine and homogeneous staggered half-step systems
through `D_25` and performs rank-aware searches with three held-out
transitions.  It finds no strict same-index simultaneous recurrence of
derivative order at most one through `deg_x<=3`, `deg_m<=2`; the preferred
homogeneous `E/O` half-step recurrence uses only first derivatives and no
forcing coordinate.

The same probe checks the recurrence's operator parts with exact `BigInt`
real-rootedness and weak interlacing.  Through `2<=m<=8` it finds the stable
chain `A_m << C_m << B_m << D_m`, the expected placements of both sums, and
`E_(m+1) << O_(m+1)`.  Two attempted all-summands cone decompositions fail,
so the probe also tests the smaller tail state

```text
G_m = (E_m - 1)/x,  H_m = O_m/x.
```

Writing `R_a[f]=(2+ax)f+x(1-x)f'`, the exact recurrence is
`H_m=2m+R_(2m-1)[G_m]`, `G_(m+1)=R_(2m)[H_m]`.  Through `m=8` the focused scan
finds

```text
G_m << H_m << G_(m+1) << H_(m+1),
```

as well as `H_m << E_m=1+xG_m` for `m>=3`.  Both adjacent tail pairs also
pass coefficient LR and coefficientwise Wronskian-sign checks.  The preferred
proof target is now the boundary packet
`G_(m+1) << H_(m+1) << 1+xG_(m+1)`; the middle transition
`H_m << G_(m+1)` is a routine Rolle/Ma--Wang step.  Use
`--interlacing-only` for this focused scan.

The probe also has `--generic-boundary`.  It rules out a black-box one-input
lemma: `h=1+3x+x^2` has distinct negative roots, but with
`g=R_4[h]=2+13x+13x^2+2x^3`, both `1+xg` and `6+R_5[g]` are not real-rooted.
The proof must therefore use the predecessor sandwich or another quantitative
spacing invariant.

The same scan finds a useful homogeneous-core fork.  For
`J_m=H_m-2m=R_(2m-1)[G_m]` and `m>=3`, both
`G_m << J_m << E_m` and `G_m << H_m << E_m` pass, while `J_m` and `H_m` do
not interlace in either direction.  Thus the proof-facing state is
`G_m -> {J_m,H_m} -> E_m`, with the first core edge supplied by Rolle.

Verification:

```text
timeout 60s nice -n 10 cargo run -q -p experiments --bin \
  derangement_index_parity_recurrence

timeout 60s nice -n 10 cargo run -q -p experiments --bin \
  derangement_index_parity_recurrence -- --interlacing-only

timeout 60s nice -n 10 cargo run -q -p experiments --bin \
  derangement_index_parity_recurrence -- --generic-boundary
```

## SymCat weighted-bond example

The main website worker owns only
`sym-poly/sym/examples/weighted_bond_symmetric_site_example.rs`.  It implements
the weighted-bond recurrence locally to verify the (P_3) and (K_3) examples
used on `symmetricfunctions.com`; it does not change the public Rust API.

At the time of this website task, the other modified Rust files remained owned
by their existing workers. The later monorepo audit adopted and committed the
previously unowned `chromatic.rs` and `lib.rs` changes; it did not touch this
separately owned example.

Verification:

```text
timeout 60s nice -n 10 cargo run -q -p sym-poly-sym \
  --example weighted_bond_symmetric_site_example
```
