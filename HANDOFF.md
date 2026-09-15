# Handoff

## Completed: linear-extension promotion — 2026-09-15

`combinatoric-core::Poset` now validates linear extensions and computes
Schützenberger promotion, individual promotion orbits, and the complete orbit
decomposition.  Orbit computations reuse one reachability matrix, and tests
cover invalid words, the two-extension V-poset orbit, chains, antichains, and
the empty poset.  All 257 crate unit tests and two documentation tests pass;
rustfmt, `git diff --check`, and strict target Clippy pass when the crate's
documented pre-existing lint classes are allowed.  Unqualified strict Clippy
still reports 38 pre-existing warnings outside this increment.  No manifest,
lockfile, other crate, experiment, or generated file changed.  The unrelated
untracked `polytool/scripts/__pycache__/` remains untouched, and ownership is
released.  The focused commit is `48543b5`; the routine push to
`origin/master` failed because this environment has no accepted GitHub SSH key.

## Completed: connected-PEO flagged extension — 2026-09-15

Worker `/root` reclaims only the ignored experiment
`experiments/src/bin/chordal_clawfree_flagged_atom_tests.rs` and this opening
entry for higher-color stress tests of the corrected connected-PEO flag
direction.  The binary now supports bounded graph-id ranges.  Four exact
chunks exhaust all `7449` ordered graphs at order seven with flags at most
five: all `2458170` graph/flag rows are atom-positive.  Focused tests, rustfmt,
and target Clippy pass.  No shared library source, manifest, lockfile, or
unrelated dirty file was changed.  Ownership is released.

## Completed: flagged chordal claw-free atom scan — 2026-09-15

Worker `/root` owns only the new ignored experiment
`experiments/src/bin/chordal_clawfree_flagged_atom_tests.rs` and this opening
entry.  The binary reuses the existing graph predicates and multivariate
Demazure-atom machinery to test both monotone flag directions along a fixed
reverse perfect-elimination order.  Its three tests, rustfmt, and target
Clippy pass.  Unqualified strict Clippy is blocked by pre-existing warnings in
dependency crates and the experiments library.  Exact results and compact TSV
certificates are recorded in the separate `Line-graph-chromatics` project.
No shared library source, manifest, lockfile, or unrelated dirty file was
changed; the pre-existing untracked `polytool/scripts/__pycache__/` remains
untouched.  Ownership is released.

## Completed: exact Bernstein-basis conversion — 2026-09-15

Worker `/root` added exact standard Bernstein-basis construction and coordinate
conversion to Polytool, including rational coefficients, degree elevation, and
the `bernstein-expansion`/`bernstein` CLI with text and JSON output.  Focused and
full library/CLI tests pass, as does strict Clippy on the library and binary.
The local checkpoint is commit `04d647f`; `git push origin master` was attempted
but the container has no GitHub SSH key (`Permission denied (publickey)`).
Ownership is released.  The pre-existing untracked
`polytool/scripts/__pycache__/` remains untouched.

## Completed: deco subexceedant polynomial generators — 2026-09-15

The recurrence-backed joint block/eligible-site distribution from the A144438
project now has canonical exact generators in `polytool::sequences::deco`.
The module provides `BigInt` APIs and checked `i64` wrappers for the joint
histograms, the uniform `w`-specializations, fixed eligible-site and selected-
join layers, and the proved gamma recurrence.  Tests recover the ordinary
Eulerian family at `w=0`, the first deco rows at `w=1`, the first layer rows,
and verify gamma expansion against independent coefficient extraction through
size 20.  The source commit is `5e07336`, and the research project records that
exact revision.

At the user's request, the formerly divergent local `master` was integrated
with the 24 newer commits from `origin/master` by a non-rewriting merge; local
recovery branch `backup/pre-origin-sync-20260915-55f9546` preserves the
pre-merge tip.  Both chronological handoff conflicts retain both sides, and
the graph conflict retains the local Cartesian/rook/triangular constructors
beside the incoming line-graph recognizer.  On the merged tree, all 253
`combinatoric-core` tests, 350 Polytool library tests, five Polytool
documentation tests, six version/budget CLI tests, 28 MCP library tests, two
MCP binary tests, the MCP stdio smoke test, and seven web tests pass.  Strict
target-only Clippy passes with the documented pre-existing warning classes
allowed; unqualified strict Clippy still reports those existing warnings.
Formatting and `git diff --check` pass.  The pre-existing untracked
`polytool/scripts/__pycache__/` remains untouched.  Merge commit `e76bda9`
was pushed non-forced to canonical `origin/master` over authenticated HTTPS;
the configured SSH remote was not changed.  Ownership is released.

## Completed: standard graph-family acyclic-sink sequences — 2026-09-13

This batch adds a documented geometric triangular-lattice patch constructor,
an exact deletion-contraction acyclic-orientation count, and an exact general
acyclic sink-polynomial generator using the sink-set identity and a shared
cache.  The exhaustive test agrees with direct orientation enumeration for
every labeled graph through five vertices; all 253 library tests and both
documentation tests pass.  Strict whole-crate Clippy remains blocked by 39
unrelated pre-existing warnings and reports no warning in the new code.  The
focused standard-family binary also passes strict Clippy.  The accidental
workspace-formatter touch to `polytool/web/src/lib.rs` was reversed exactly.
The pre-existing untracked `polytool/scripts/__pycache__/` remains untouched.
Ownership is released; no push, PR, or publication was made.

## Completed: edge/triangle block-tree documentation correction — 2026-09-13

This documentation-only follow-up qualifies that the constructor covers the
nontrivial connected graphs in the stated class, while the one-vertex graph is
the sole degenerate exception.  Rustfmt and the documentation tests pass.  No
remote action was taken, and ownership is released.

## Completed: disjoint-union graph operation — 2026-09-13

Commit `f74c2be` adds the missing generic disjoint-union operation required by
the componentwise OEIS families, with exact graph and chordal-sink
multiplicativity tests.  It was fast-forwarded into local `master`.  All 252
library tests and both documentation tests pass with the external target
directory.  No remote action was taken, and ownership is released.

## Completed: chordal graph API and reusable family generators — 2026-09-13

Commit `aedd875` promotes reusable code from the chordal claw-free sink project
to `combinatoric-core/src/graph.rs`: maximum-cardinality-search perfect
elimination certificates and `is_chordal`; exact `BigInt` and checked `i64`
chordal acyclic-sink generators using the proved independent-set expansion;
generic Cartesian products; path powers; rook and triangular graphs; and a
validated edge/triangle block-tree root constructor.  Documentation marks the
rook and triangular families as generally nonchordal.

All 251 library tests and both documentation tests pass with the external
target directory.  New exhaustive tests compare the chordality result with an
independent simplicial-deletion implementation and compare the new sink method
with direct orientation enumeration on every labeled graph through five
vertices.  Closed-form tests cover paths, `K_n-e`, balanced double-star line
graphs, and friendship line graphs.  Strict whole-crate Clippy remains blocked
by 40 unrelated pre-existing warnings; it reported no warning in the newly
added code.  The clean feature branch was fast-forwarded into local `master`.
No push, PR, or publication was made.  Ownership is released.  The unrelated
untracked `polytool/scripts/__pycache__/` remains untouched.

## Completed: frozen chordal claw-free sink scan — 2026-09-12

Worker `/root` owned only the new ignored experiment
`experiments/src/bin/frozen_chordal_clawfree_sinks.rs` and this handoff entry.
The exact scan verifies two order-six frozen counterexamples, checks all
`105398` frozen states on all `64` natural unit interval area sequences through
order five, and checks the unfrozen polynomial on all `1291` connected chordal
claw-free graphs among the `273193` connected unlabelled graphs through order
nine.  The three focused tests pass, and the bounded release scans completed.
The experiment uses the existing
`Graph::acyclic_sink_polynomial_with_frozen_edges` and `polytool` exact
real-rootedness checks.
No shared library source, manifest, lockfile, or unrelated dirty file is
owned; in particular `combinatoric-core/src/graph.rs` remains untouched.  The
pre-existing untracked `polytool/scripts/__pycache__/` remains untouched.
Ownership is released; no push or publication was made.

## Completed: A144438 quasisymmetric-refinement test — 2026-09-11

Worker `/root` owned only new
`experiments/src/bin/a144438_deco_qsym.rs` and this handoff entry for the
user-authorized test of the decorated-permutation descent-set refinement in
the shared QSym library.  No existing Rust source, Cargo manifest, generated
cache, or unrelated experiment is owned.  The pre-existing untracked
`polytool/scripts/__pycache__/` remains untouched.  The exact computation
finds the first symmetry failure in degree three, where the monomial
coefficients at `(1,2)` and `(2,1)` are `4` and `3`.  Quasisymmetric Schur
positivity holds through degree five and fails in degree six, with coefficient
`-8` at `(1,1,2,1,1)`.  All three focused binary tests pass with the external
target directory and reduced priority; the bounded scan stopped at degree
eight.  Results are recorded in the separate A144438 project.  Ownership is
released after the focused checkpoint; no push or publication was made.

## Completed: factored Delannoy rank-six search — 2026-09-11

Worker `agent-research-delann-p-39629cb1` reclaims only
`polytool/examples/delannoy_matching_search.rs` and this handoff entry for the
user-requested continuation of research issue PerAlexandersson/research#24.
The example now has a narrow tested `--divide-one-plus-t` option, allowing the
newly found
`R_6=(1+t)Q_5` decomposition to reuse the existing streamed Graph matching
search. No library module, Cargo manifest, or graph implementation is claimed.
The focused example test suite passes (3 tests), and the release example
builds with `CARGO_TARGET_DIR=/cargo-target/ai-projects`. It exhaustively
checked 5,596,003 trees at quotient-root order 22 and 57,359,412 unicyclic
graphs at order 21, finding no exact root. Research details are recorded under
`/workspace/research/notes/delannoy-clawfree/`. Ownership is released after
the focused checkpoint. The unrelated untracked
`polytool/scripts/__pycache__/` remains untouched; no Rust push was made.

## Completed: Delannoy claw-free graph search — 2026-09-11

Worker `agent-research-delann-p-39629cb1` owns only
`polytool/Cargo.toml`, its one-line generated `Cargo.lock` dependency update, new
`polytool/examples/delannoy_clawfree_search.rs`, and this handoff entry for
research issue PerAlexandersson/research#24. The example will stream graph6
input, obtain reversed Delannoy-square targets from the existing canonical
`polytool` generator, and verify claw-freeness plus exact independence
coefficients through `combinatoric-core::Graph`. The Cargo edit is limited to
the example's dev-dependency. The worker also owns new
`polytool/examples/delannoy_matching_search.rs` for the explicitly authorized,
separate matching-graph search. No `graph.rs` edit is currently planned.

Code checkpoint `389222e` adds the two streaming examples. The direct
claw-free scan is exhaustive through n=4 (315 witnesses at n=4), and the
isolate-free matching scan is exhaustive through n=4 (six matching witnesses
at n=4). A restricted n=5 scan checked all 12,413,039 connected graphs with 18
edges and 17--19 vertices and found one exact matching witness; its line graph
is an independently checked claw-free witness. Results and exact scope are in
research commit `82bce66` under
`/workspace/research/notes/delannoy-clawfree/`.

Both focused example test suites and `git diff --check` pass. Strict Clippy is
blocked by unrelated preexisting warnings in `combinatoric-core` and
`polytool`; whole-workspace rustfmt checking is blocked by unrelated existing
formatting in an experiment and `polytool/web`. The two new examples were
formatted directly. The unrelated untracked
`polytool/scripts/__pycache__/` and all other worker files remain untouched.
Ownership is released. No push was made.

## Completed: explicit literature polynomial generators — 2026-09-11

Host supervisor owns `polytool/src/sequences.rs` (one module declaration),
new `polytool/src/sequences/literature.rs`, new
`polytool/examples/literature_sequence_rows.rs`, and this entry. Implement
exact reusable Eulerian/Delannoy matrix-square and Hoggatt family APIs with
independent small tests for the user-requested research screen. No other
worker owns these files; preserve unrelated `polytool/scripts/__pycache__/`.
No public push or changes to other crates, shared Lean builds or worker
profiles. Research evidence belongs in the separate private research repo.
Initial APIs committed as `8a3443e`. Follow-up exact regression tests check
the Hoggatt rank-shift identity (including the missing leading-term boundary)
and a Delannoy three-lag differential identity through degree120. All 23
sequence tests pass (312 other library tests filtered). The private research
note also verifies the differential identity by rational generating-function
algebra; no real-rootedness conjecture is claimed proved. Ownership released
after this focused checkpoint. No public push.

## Principal-stratum conjecture checks — 2026-09-10

Host read-only computation found a counterexample to the draft's
`eq:H-repeat3`: for `(3,1,1,3,1,1,3)` at d=23, exact integral homology is
`H_8=Z^7`, `H_9=Z^2`, with no other groups or torsion. The formula predicts
`7s^8+3s^9+s^10`, rather than the computed `7s^8+2s^9`.
Integral cancellation/Smith replay and unreduced ranks over F_251 and F_2
agree. The d=21 calibration matches; `(3,3,1,5)` at d=28 is integrally
acyclic as predicted. General resonance and torsion-freeness remain open.
Part (i) of the absorption conjecture follows from later draft theorems.
See `docs/PRINCIPAL_STRATUM_CONJECTURE_CHECK.md` for proof, qualifications,
commands, and next targets. Host owns this note/entry through the focused
documentation checkpoint, then releases them. No source edits, worker
launches, remote jobs, manuscript changes, or outbound mail.

## Sparse integral homology — completed 2026-09-10

Task ownership is released following host review closure `c853ea7`; no push was
made. Shared CSR/mutable sparse storage, bounded cancellation, bounded exact
Smith replay/verification, and compatibility re-exports are in
`combinatoric-core`; the principal-stratum automaton, count DP, streamed
input-order modular boundary path and driver are in `experiments`. The driver
records/replays cancellation certificates against the residual and, with its
recorded Smith operations, verifies claimed Smith factors before printing
groups. Resource guards cover cells, boundary terms, mutable shape slots,
NNZ, coefficient bits, dense Smith entries/shape/operations/bits, and invalid
finite-field/UCT moduli.

Final verification after `b134824`: `cargo test -p combinatoric-core --lib`
(241 passed), `cargo test -p sym-poly-core --lib` (123 passed), and
`cargo test -p experiments --lib` (25 passed), all with
`CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10`. The
same command prefix also ran the release driver smoke:
`cargo run --release -q -p experiments --bin principal_stratum_homology --
--omega 1,1 --d 4 --max-cells 20 --max-nnz 100 --max-reduction-nnz 100
--max-pivots 20 --record-certificate`; it reported 9 cells, Euler -1,
`H_3 = Z`, and a matching residual certificate.
Release calibrations, guarded and without larger cases: d18 `(3,1,1,3)` gave
8280 cells, Euler -2, `H7=Z^2`, 41ms cancellation; d24 `(3,1,1,5)` gave
76384, Euler 0, `H7=H8=Z`, 829ms cancellation; d26 gave 212900, Euler 0,
`H8=H9=Z`, 3628ms cancellation. d26 first correctly rejected a 2,000,000
boundary-term guard; its actual initial NNZ was 1,706,568 under a 4,000,000
guard. Claude Code 2.1.267/Opus completed a read-only review of `63e9861`
over 15 minutes; its focused re-review found and `d29227c` resolved one
u64-primality witness defect. Results and resolutions are in
`docs/PRINCIPAL_STRATUM_HOMOLOGY_REVIEW.md`, with log
`/tmp/principal-stratum-claude-63e9861.log`; the final `d29227c` micro-review
also returned no blocker. Unrelated ignored
`polytool/scripts/__pycache__/` remains untouched.

The final missing-map fix avoids allocating absent shaped-zero differentials in
validation and Smith assembly; `differential_or_zero_with_limits` provides
checked zero retrieval, and metadata-only complexes with `10^12` generators
are covered without generator allocation. Commit `b134824` also checks shape
slots and NNZ before cloning a stored map, with 1-by-1/nonzero and 9-by-0
stored-map budget regressions; the trusted infallible helper is retained.
Focused `chain_complex` tests pass (7 tests). Claude Code 2.1.267/Opus
reviewed `b134824` read-only in plan mode and found no correctness blocker;
log: `/tmp/principal-stratum-claude-stored-retrieval-b134824-retry.log`.

Host direct-executable measurements supersede earlier worker timing claims:
d18/d24/d26 integral-and-replay elapsed 0.131/2.452/10.787 seconds; d26 peak
RSS was 482912 KiB. The host observed no source edits or Cargo builds.

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
## Completed Polytool GitHub issues #1 and #2 (2026-09-07)

The coding worker implemented the Polytool issues in the isolated
worktree `/tmp/polytool-issues-1-2-20260906` on branch
`fix/polytool-issues-1-2-20260906`, based on freshly fetched
`origin/master` commit `6486f93`.  The shared `/workspace/rust` checkout and
its four local-only commits remain untouched.

Implementation ownership was limited to `HANDOFF.md`, `polytool/HANDOFF.md`,
`polytool/build.rs`, `polytool/src/version.rs`, `polytool/src/lib.rs`,
`polytool/src/recurrence.rs`, `polytool/src/bin/polytool.rs`, focused new or
existing tests under `polytool/tests/`, `polytool/README.md`,
`polytool/mcp/src/lib.rs`, and `polytool/mcp/README.md`.  No other worker
claimed these files.  The separately owned derangement experiment and SymCat
weighted-bond example were not touched.  Ownership is released by this final
handoff after canonical and standalone publication verification.

The task is to implement deterministic CLI build-version reporting and an
exact candidate budget for adaptive recurrence finding, keep library/CLI/MCP
outcomes aligned, verify the canonical branch, publish it without including
the shared checkout's local-only history, project `polytool/` to standalone
`main`, and close GitHub issues #1 and #2 after the accepted verification gate.

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
status checks.

On 2026-09-07 the user clarified that the earlier green-CI wording reflected
supervisor caution, not a user-imposed condition, and accepted the complete
local monorepo and standalone-subtree checks as the merge gate.  The earlier
claim that user authorization was still required was therefore incorrect.

PR #3 merged without force as canonical monorepo commit `0675132`, whose
parents are prior `origin/master` `6486f93` and task head `b5c7976`; none of
the shared checkout's four local-only commits entered the merge.  The
documented sync script produced the first published standalone projection
`49fff55`.  The actual monorepo binary reported
`polytool 0.2.1-rc.5 (git 0675132a0408)`, while the standalone binary reported
`polytool 0.2.1-rc.5 (git 49fff5546586)`.  Both returned structured
`budget_exhausted` output with exact zero/one candidate counts and exit status
3; boundary success also passed, as did all six standalone CLI tests.

Release handoff commit `8245e3a` was pushed non-forced to canonical
`origin/master`, and the documented sync script projected it to standalone
`origin/main` commit `7e4ac54`.  GitHub issues #1 and #2 were then closed as
completed with separate notes citing PR #3, implementation checkpoints,
published refs, and verified behavior.  This root-only final-ref record does
not alter `polytool/`, so the standalone split remains `7e4ac54`.  Ownership
is fully released.
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
