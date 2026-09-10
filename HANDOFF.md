# Handoff

## Sparse integral homology — completed 2026-09-10

Owner released task files after commits `63e9861` and `5bd050b`; no push was
made. Shared CSR/mutable sparse storage, bounded cancellation, bounded exact
Smith replay/verification, and compatibility re-exports are in
`combinatoric-core`; the principal-stratum automaton, count DP, streamed
input-order modular boundary path and driver are in `experiments`. The driver
records/replays cancellation certificates against the residual and, with its
recorded Smith operations, verifies claimed Smith factors before printing
groups. Resource guards cover cells, boundary terms, mutable shape slots,
NNZ, coefficient bits, dense Smith entries/shape/operations/bits, and invalid
finite-field/UCT moduli.

Verification: `cargo test -p combinatoric-core --lib` (239 passed),
`cargo test -p sym-poly-core --lib` (123 passed before final API-only fixes),
`cargo test -p experiments --lib` (25 passed), and target-only
`cargo check -p experiments --bin principal_stratum_homology`, all with
`CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10`.
Release calibrations, guarded and without larger cases: d18 `(3,1,1,3)` gave
8280 cells, Euler -2, `H7=Z^2`, 41ms cancellation; d24 `(3,1,1,5)` gave
76384, Euler 0, `H7=H8=Z`, 829ms cancellation; d26 gave 212900, Euler 0,
`H8=H9=Z`, 3628ms cancellation. d26 first correctly rejected a 2,000,000
boundary-term guard; its actual initial NNZ was 1,706,568 under a 4,000,000
guard. Claude Code 2.1.267/Opus completed a read-only review of `63e9861`
over 15 minutes; its no-blocker result and resolutions are in
`docs/PRINCIPAL_STRATUM_HOMOLOGY_REVIEW.md`, with log
`/tmp/principal-stratum-claude-63e9861.log`. Unrelated ignored
`polytool/scripts/__pycache__/` remains untouched.

## Sparse integral homology implementation — active 2026-09-10

The assigned implementation worker owns the task's shared library modules,
model driver, exact tests, narrowly required registrations and this opening
entry; the host supervisor performs read-only coordination and verification.
Checkpoints 1, 3, and 4 now have their first verified shared implementation:
`combinatoric-core::sparse_matrix` provides generic shape-preserving CSR and a
bidirectionally indexed mutable `BigInt` form; `chain_complex` provides
graded unit cancellation with independently replayable SHA-256-bound pivot
certificates; and `integer_linear_algebra` provides a bounded exact BigInt
Smith reducer whose optional operation certificate is replayed in tests.
The full `combinatoric-core --lib` suite (235 tests) passes. No task source
files were dirty at launch;
unrelated `polytool/scripts/__pycache__/` is preserved. Claude review is
mandatory, read-only, and starts with a monitored 30-minute allowance after
the implementation is frozen. Its launching worker owns the child process. No
pushes, publications or Abacus jobs are authorized.

The frozen implementation checkpoints are `cc0dbf7` (CSR/mutable sparse
storage), `c21e104` (bounded BigInt Smith forms and integral cancellation),
and `75eb540` (principal-stratum driver, field baseline, compatibility
adapters). Focused verification passed: 235 `combinatoric-core` library tests,
123 `sym-poly-core` library tests, and 21 `experiments` library tests. The
driver exactly reproduces the small d=4 calibration and the d=18 target:
8,280 cells, Euler -2, and residual integral `H_7 = Z^2` after 4,139 unit
pivots in 8.6 seconds. The bounded d=24 field run completed: 76,384 cells,
Euler 0, and exact F_251 dimensions `b_7=b_8=1`; the integral reducer stayed
at roughly 163 MiB RSS but hit its 10-minute cap before producing a result.
The d=26 field-only attempt stayed below roughly 290 MiB RSS but likewise hit
its 10-minute cap during validation without output. These are resource-limited
incomplete runs, not negative results or homology claims. A read-only Claude
review is now required before release; see
`docs/PRINCIPAL_STRATUM_HOMOLOGY_REVIEW.md`.

Claude's completed read-only Opus review found and the implementation worker
fixed a malformed-Smith-certificate acceptance bug, a quadratic sparse
composition validation path, and missing user-visible cancellation-certificate
and universal-coefficient checks. The malicious Bézout certificate now has a
regression test. A focused Claude re-review of this material correction is in
progress; do not treat this paragraph as final review clearance yet.

That focused re-review found same-index Bézout, self-add and rectangular-index
certificate edge cases; all are now guarded with regression tests. A final
guard-only Claude confirmation remains the last review action.

Final Claude confirmation found the guards correct and only asked for two
rectangular rejection tests, which were added and passed. Review is complete:
Claude Code 2.1.267 requested `opus` in read-only plan mode reviewed commits
through `bc70287`; no reviewer edits occurred. The implementation worker has
released ownership of all task files. Remaining limits are intentional: dense
Smith is budgeted, d=24 integral cancellation and d=26 field validation hit
their 10-minute resource caps, and benchmark logs are under `/tmp`, not
Dropbox. Exact successful commands include:

```text
CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10 \
  cargo test -p combinatoric-core --lib       # 235 passed
CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10 \
  cargo test -p sym-poly-core --lib           # 123 passed
CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10 \
  cargo test -p experiments --lib             # 21 passed
CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10 \
  cargo run -q -p experiments --bin principal_stratum_homology -- \
  --omega 3,1,1,3 --d 18 --max-cells 10000 --max-nnz 250000 \
  --max-reduction-nnz 1000000 --max-pivots 10000
```

## Principal-stratum homology plan (2026-09-10)

The user requested a careful plan, with reusable sparse matrices and Smith
normal form in the shared library. See docs/PRINCIPAL_STRATUM_HOMOLOGY_PLAN.md.
Planning only: no code, dependencies, builds, workers or profiles changed.
The plan uses additive combinatoric-core modules and a separate model driver,
reuses existing modular solvers, and preserves standalone Polytool packaging.
It specifies exact unit cancellation, bounded Smith fallback, certificate
replay, torsion tests, and separate timing/memory measurements for each stage.
Original computation scripts were not attached to Boris's forwarded email;
the reported homology remains an independent-reproduction target. The
supervisor's small automaton DP verified cell counts only. Ownership of this
entry and the plan is released after the documentation checkpoint.

## Path-IC frozen-core certificate (2026-09-06)

The optional `--verify-frozen-core` scan in
`experiments/src/bin/path_ic_crystal_wall.rs` tests insertion-word
independence of the positioned subword on vertices `1,2,3` and the
inverse-column event which emits each complemented core vertex.  Through
`P_17`, all `198812` four-color source/repair incidences and `596436` event
signatures agree within each length type and recording tableau.  An
independent all-standard-recording audit through `P_11` covers `9190`
incidences.  The same scan verifies that all `27` forced-debt word graphs are
connected to the explicit canonical word by `528` legal length-two/three
edges of `12` types, with maximum distance six.  Four focused binary tests,
strict target-only Clippy with `--no-deps`, JSON parsing, and the calibrated
release run pass.  No reusable API changed.  Ownership is released after the
checkpoint commit.

## Path-IC crystal-wall tool (2026-09-06)

Reusable checkpoint `6ddb666` adds `sym-poly-core::p_rs`: an abstract
`PInsertionOrder`, one exact traced implementation of inverse column
Algorithm 2, six event kinds, boundary certificates, reverse complement, and
six focused tests.  It also corrects `Tableau::evacuation` to the direct
shape-preserving jeu-de-taquin algorithm; the old RSK shortcut transposed
non-self-conjugate shapes.  Two regression tests cover that correction.  All
`123` `sym-poly-core` tests pass.  Strict whole-crate Clippy is blocked by
pre-existing warnings in unrelated core/dependency modules; ordinary Clippy
reports no warning in the new `p_rs.rs` or changed `tableau.rs` code.

The separate tracked research driver
`experiments/src/bin/path_ic_crystal_wall.rs` supplies only path-IC data:
the forced-debt ternary word grammar, path order/ladders, duplicate carry,
raw recording reconstruction, and the excluded wall census.  Its three unit
tests and strict target-only Clippy pass.  The release scan through `P_17`
finds exactly `27` types, `282` insertion words, `989` distinct rejected
tableaux, and `7852` incidences; all keep the complemented maximum active
through a boundary ladder and have `p_Q(1)>m_1(Q)`, split into `5060` ladder
copies and `2792` ladder moves.  Its per-type TSV agrees with the independent
Python trace on all `27` rows and `13` compared columns.  The architecture and
next finite-state-transducer increment are recorded in
`../projects/Line-graph-chromatics/crystal-wall-rust-tool-plan-2026-09-06.md`.
Follow-up normalization shows that every participating ladder begins with
adjacent active/source letters `AS`, never contains the distinguished maximum,
and has no `AA`.  In the row-two-starts-three obstruction its merge word is
strictly alternating.  In the multiple-row-one-twos obstruction it has at
most one `SS`; only `90/7852` incidences use that exceptional transition.
Ladder copy/move is exactly whether the merge ends in `A`/`S`.  These regular
languages are now asserted by the driver and reduce the transducer to an
alternating state plus one optional source-doubling state.  Ownership is
released after the follow-up checkpoint; no unrelated experiment binary or
library module was changed.

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
