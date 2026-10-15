# Handoff

## Active: explicit ordinary bounce lift — 2026-10-15

`/root` owns only `experiments/src/bin/marked_coloring_hrv_bijection.rs`
and this handoff on branch
`exp/verticalstripllt-marked-coloring-bijection`.  The current checkpoint is
`675b791`.  The task is to replace the prototype's exact-bucket bounce
fallback with reversible first-/second-bounce maps while retaining labelled
HRV fibre fillings.  No maintained Rust API, Lean source, paper source, or
circular code is owned.  Build output remains under
`/cargo-target/ai-projects`; the shared Rust checkout is untouched.

The single-bounce NN and DN cases are now explicit in the owned prototype.
They port the paper's coloring involutions and Lean's `bounceMatch` /
`secondMatch` orientation correspondences, with their explicit inverses.
Labelled columns are transported by the proved identity-or-adjacent-swap HRV
conjugacy: each fibre's actual color set is moved and placed increasingly.
Dedicated K3 and `NDNDEE` second-bounce tests pass.  The full rank-five census
still passes all 277,302 certificates: 162,160 certificates traverse an
explicit first bounce, 23,872 traverse an explicit second bounce, and 78,076
still encounter the residual generalized-bounce fallback.  The next change is
to replace only that residual fallback by canonical matching at the smaller
recursion branches, since the published many-bounce derivation uses algebraic
cancellation and does not provide a local object map.

## Completed prototype: ordinary marked-coloring/HRV bijection — 2026-10-15

Checkpoint `bedd1a2` on task-local branch
`exp/verticalstripllt-marked-coloring-bijection` contains the prototype and
releases this handoff, `.gitignore`, and the new experiment artifact.
The shared checkout and its existing branch remain untouched.  No maintained
crate API, Lean source, paper source, circular LLT code, or external action is
owned.  The prototype tests a reversible certificate-producing map from an
ordinary packed marked coloring to an HRV orientation with increasing fibre
columns, preserving packed content and marked-ascent degree, beginning through
rank five.  Build output must use `/cargo-target/ai-projects`.

The isolated prototype is
`experiments/src/bin/marked_coloring_hrv_bijection.rs`, retained as a narrowly
tracked reproducibility artifact.  It enumerates ordinary Schroeder paths and
decodes normal/strict edges by the same zero-based convention as
`VerticalStripLLT/Path.lean`; every natural edge set is independently checked
against `sym_poly_sym::unit_interval_edges`.  Ordinary HRV maps reuse
`sym_poly_sym::llt::highest_reachable_vertices`.

The implemented bijection is hybrid and fully reversible.  It follows the
positive shifted Alexandersson--Sulzgruber reduction objectwise for column
bases, diagonal-return products, and both branches of the unicellular
relation.  At the first east-free prefix ending in a diagonal (the first
genuine bounce seed), it uses sorted rank matching inside the exact bucket
`(ordered packed-color content, q-degree)`.  The certificate labels this as a
canonical fallback; it is not claimed to be the missing local bounce rule.
Both forward and inverse maps replay the same reduction decisions.

Independent source and target enumerators verify validity, exact color
multisets, q-degree, injectivity, surjectivity, and trace reversal for every
ordinary decorated path through rank five: 258 paths, 9,271 exact buckets,
and 277,302 source/target certificates.  Rank counts are
`1, 1, 3, 11, 45, 197` paths and `1, 1, 8, 144, 4,796, 272,352`
certificates.  The complete K3 map is printed by `--k3-map`; a readable trace
is printed by, for example,
`--path NNNEEE --trace-coloring 1,3,2 --trace-mask 3`.

The first exact obstruction to the naive local rule `eta = A` already occurs
for K3.  With three distinct colors and q-degree two, the orientation with
ascending edges `12,13` has two marked-coloring sources but three HRV-column
targets.  Thus a genuine lift must mix orientations; the mismatch is not a
sorting artifact.  The next task is to replace the bounce fallback with the
explicit first-/second-bounce orientation correspondences while transporting
actual labelled HRV fibres, not just their size partition.

Focused commands (all with `CARGO_TARGET_DIR=/cargo-target/ai-projects`) pass:

```bash
cargo check --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection
cargo test --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection
cargo run --quiet --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection -- --max-n 5
```

Strict Clippy for the requested binary is clean when allowing the pre-existing
`clippy::manual_is_multiple_of` finding in the unrelated
`experiments/src/principal_stratum.rs`; unqualified `-D warnings` stops at that
existing library diagnostic before checking the binary.  No push was made,
preserving the Rust workspace's recorded no-push boundary.

## Completed: monotone-pattern descent/peak recurrence — 2026-10-14

Host `/root` completed and releases `combpoly/src/monotone_pattern.rs`,
`combpoly/src/lib.rs`, `combpoly/README.md`, the checked
`combpoly/examples/increasing_pattern_avoidance.rs`, `docs/LIBRARY_GUIDE.md`,
and this handoff.  The new arbitrary-precision APIs generate descent and
ordinary peak polynomials on permutations avoiding `12...m` by an exact RSK
recording-tableau recurrence.  The shape has first-row length at most `m-1`;
the implementation's initial at-most-rows orientation was rejected by the
independent coefficient test and corrected before checkpointing.

For `m=4`, both families agree coefficient-for-coefficient with direct
`Av(1234)` enumeration through `n=8`, and descent-row sums reproduce the known
counts through `n=10`.  Exact polytool checks find real-rootedness and directed
weak adjacent interlacing through `n=30` for both statistics.  A scan through
`n=60` was stopped after one minute because exact root certification was slow;
generation itself was immediate and no failure was found.  The complete
combpoly suite passes (205 library tests, 2 binary tests, and doctests), as do
the checked example, strict Rustdoc, workspace rustfmt, and diff checks.
Research-only scan and `1342` skew-block drivers remain ignored.  No other
maintained source was changed and the existing no-push boundary remains.

## Completed: claw-free Schur counterexample fixture — 2026-10-14

`/root` completed and releases `combinatoric-core/src/graph.rs`,
`combinatoric-core/README.md`, `sym-poly/sym/src/chromatic.rs`, and this
handoff.  `Graph::claw_free_schur_counterexample()` is the 12-vertex,
22-edge graph6 graph ``K?`CRAWWUXIM`` from J. Prajapati,
*A counterexample to the claw-free Schur-positivity conjecture*,
arXiv:2607.26364.  Core regressions check its graph6 identity, connectedness,
claw-freeness, and line-graph property.  The symmetric-function regression
checks `[s_(3,3,3,3)]X_G=-64` exactly.  The public constructor is documented
in the core README.

The shared checkout was clean at `4349e40` before these edits.  Full checks
pass: 288 core unit tests, 12 core integration tests, 4 core doctests; 137
sym unit tests, 100 sym integration tests with 3 benchmark tests ignored, and
2 sym doctests.  Strict Rustdoc for both crates, workspace rustfmt, and diff
checks also pass.  Exact commands used `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo fmt --all --check
cargo test -p combinatoric-core
cargo test -p sym-poly-sym
RUSTDOCFLAGS='-D warnings' cargo doc -p combinatoric-core -p sym-poly-sym --no-deps
git diff --check
```

No dependencies or unrelated APIs changed.  The focused checkpoint is local
to the existing branch; no push was made, preserving the standing no-push
boundary.

## Completed: public permutation-poset constructor — 2026-09-30

The 3-Matroids Dependent-order worker completed and releases only
`combinatoric-core/src/poset.rs`, its README, the new checked
`combinatoric-core/examples/permutation_poset.rs`, library/promotion guide
rows, this handoff, and four existing ignored permutation-poset helper
bodies. User explicitly authorized this narrow shared-library promotion.
The checkout was clean; supervisor report shows no live Rust worker.
Use the existing one-indexed permutation validator and Hasse reduction;
zero-based poset vertices are permutation positions. No push.

`Poset::from_permutation(&[usize])` is public in `combinatoric_core::poset`.
It accepts a permutation of 1..=n (including empty), validates with the
existing helper, and returns the reduced dominance order on positions 0..n.
All 5914 permutations of sizes 0..=7 pass independent cover-interval,
dominance-reachability and LIS/height checks. Empty/monotone, invalid-input,
diamond and paired-block polynomial fixtures also pass. Full core suite:
287 unit tests, 12 integration tests, 4 doctests; checked example and
strict Rustdoc pass, as do targeted rustfmt and diff checks.

Exact commands (with `CARGO_TARGET_DIR=/cargo-target/ai-projects` and `nice -n 10`):

```bash
cargo test -p combinatoric-core permutation_poset
cargo test -p combinatoric-core
cargo run -p combinatoric-core --example permutation_poset
cargo check --manifest-path experiments/Cargo.toml --bin stembridge_hstar --bin stembridge_maxim --bin stembridge_delete --bin permutation_poset_w_eulerian
RUSTDOCFLAGS='-D warnings' cargo doc -p combinatoric-core --no-deps
```

The four ignored consumers now delegate to the public constructor via a
lossless u8-to-usize conversion; all compile. Their ignored files remain
untracked, per policy, while the maintained API/example are committed.
The pre-existing unused `perm` warning in stembridge_hstar is unchanged.
No mathematical search, dependency update, or unrelated source edit.

## Completed bounded fit: ladder recurrences — 2026-09-30

Host `/root` completed and releases this handoff and ignored
`experiments/src/bin/connected_subgraph_packings.rs`. User explicitly rejected
building a transfer model and requested existing polytool recurrence fitting
from the polynomial rows. No transfer implementation or mathematical library
change is authorized. Reuse current ladder rows through n=10; hold out the
last row and distinguish a fitted relation from proof. No public push.
Polytool fits an order-6, coefficient-degree-4 recurrence for k=2 using
P_0..P_9, and P_10 passes the exact held-out identity. For k=3 no fit passes
the same held-out check within order<=7 and coefficient degree<=order.
Training only on P_0..P_8 does not return a held-out-valid k=2 candidate in
that box either; the result is not certified unique or proved for all n.
Research #39 / RECURRENCES.md records the formula, initial rows and scope.
No transfer code or maintained mathematical source was changed. Release
consumer build, targeted formatting and diff checks pass.

## Completed: rectangular-grid packing check — 2026-09-30

Host `/root` completed and releases only this handoff and the ignored consumer
`experiments/src/bin/connected_subgraph_packings.rs` for the requested grid
check. Reuse `Graph::grid` and the canonical packing generator at `455c5ad`;
no maintained mathematical source changes. Scan rectangles with at most
20 vertices for k=2,3. Exact polytool Rust APIs replace unavailable MCP
adapters; private research Markdown records evidence. No public push or
worker launch; previous ownership is released and checkout starts clean.
All 35 rectangles up to transposition with 1<=rows<=cols and rows*cols<=20
pass for each of k=2 and k=3. The default exact PRS and rational Sturm tests
agree on all 70 cases. Bounds: 1x1..20, 2x2..10, 3x3..6, 4x4..5. This is
finite evidence only, not a proof for arbitrary rectangles. The ignored
consumer accepts `<k> grids`; no mathematical library code changed.
Release build, targeted rustfmt, and diff checks pass. No public push.

## Completed: bounded connected-subgraph packing check — 2026-09-30

Host `/root` completed and releases this handoff, `combinatoric-core/src/lib.rs`, new
`combinatoric-core/src/graph_packings.rs`, `combinatoric-core/README.md`,
and ignored `experiments/src/bin/connected_subgraph_packings.rs` only.
Also releases `docs/LIBRARY_GUIDE.md` and `docs/EXPERIMENT_PROMOTION.md` to
register the new canonical generator without publishing research results.
The shared checkout was clean at `24559ad`; no live Rust worker appeared
in the supervisor report. User explicitly requested a small Rust root check.
Reuse Graph and polytool; add a documented BigInt generator with independent
edge-subset tests, then run bounded scans. Preserve the existing no-push
boundary. No other workers, models, profiles, or files are assigned.
The BigInt generator is canonical in `graph_packings`; its independent
edge-subset oracle verifies all graphs through four vertices and the specific
12-vertex tree found by the scan. For k=2, its coefficients are
`[1,27,142,212,117,21,1]`; exact PRS and rational Sturm both reject real-rootedness,
with square-free degree 6 and exactly 4 real roots. The same tree fails for
k=3. Private research issue #39 records the graph and reproducibility details.
Full verification: 283 unit tests, 12 integration tests, 3 doctests pass;
targeted rustfmt and `git diff --check` pass. The scan driver remains ignored.
No public push, broader search, or worker launch was performed.

## Completed: induced P4 predicate and unit-interval enumeration — 2026-09-30

`Graph::induced_p4_witness` returns a deterministic four-vertex certificate in
path order, and `Graph::is_p4_free` recognizes cographs in `O(n^4)` time and
constant auxiliary space.  Tests cover positive and negative graph families,
relabeling, union/join closure, and natural unit interval examples.

An ignored task-local `experiments/src/bin/uig_cograph_enumeration.rs` scanned
all Catalan area sequences through length 12 using the existing streaming
generator.  The P4-free counts for `n=0..12` are

```text
1, 1, 2, 5, 13, 33, 82, 202, 497, 1224, 3017, 7439, 18343
```

and the connected counts are

```text
0, 1, 1, 2, 4, 7, 11, 16, 22, 29, 37, 46, 56.
```

For `n>=1`, the latter are A152947, `1 + binomial(n-1,2)`.  The classification
and proof are in symmetricfunctions.com commit `5cb596c`, which also derives
the all-case generating function `(1-x)^3/(1-4x+5x^2-3x^3)`.

The focused tests and full `combinatoric-core` suite pass: 280 unit tests, 12
integration tests, and 2 doctests.  Formatting and `git diff --check` pass.
Strict package Clippy remains blocked by 38 pre-existing diagnostics; none is
in the new predicate or tests.  The ignored experiment is not part of this
checkpoint.  Ownership is released, and the existing no-push boundary remains
in force.

## Completed: square-graph W/forgotten expansions — 2026-09-25

`square_partition_graph_symmetric(&lambda, n)` implements Section 2.3's
even-degree square basis: doubled stars for even parts; consecutive decreasing
odd parts give paired stars joined by a doubled edge. It returns the full S_n
sum over BigInt, not the square of an already symmetrized graph polynomial.
The new `OddDegree` error complements existing index/vertex validation.
`examples/square_graph_w.rs` prints both direct forgotten and W coefficients,
normalized by n!; W coefficients reuse `c_n(...).to_forgotten_basis()`.

`docs/square-graph-expansions.md` records exact examples and the degree-four
counterexamples to uniform positivity in both bases. Direct forgotten:
H_(4)=(n-2)! [6 f_(2,2)-4 f_(3,1)-(n+1) f_(4)]. This has a direct algebraic
proof, not just interpolation. All 141 direct-forgotten cases through even
degree 12 at n=d,2d,2d+3 roundtrip; 123 have mixed signs, with the remaining
18 exactly the indices (2^k). The note proves (-1)^k H_(2^k) strictly
forgotten-positive for n>=2k. Mixed signs for every other index remain only
a bounded observation. W tests check 78 cases through degree 10, no singleton
indices, and exact inverse recovery.

Six new integration tests include independent literal permutation sums for
small doubled/glued graphs, multiple odd-pair ordering, boundary errors,
translation checks, and direct-forgotten/W regressions. Full crate tests pass
(238 passed, 3 pre-existing ignored); checked example, strict Rustdoc, targeted
formatting, and diff checks pass. Focused Clippy reports no diagnostics in the
owned graph files/example (unrelated warnings remain). Logs are
`/tmp/square-graph-{tests,doc,clippy}.log`, `/tmp/square-graph-clippy.jsonl`;
the larger probe output is `/tmp/square-graph-forgotten-results.txt`.
Paper-cache still had no copy; the source PDF was read through the browser.

All task files are released after the local checkpoint on the existing branch.
No other worker, Lean build, or publication; the no-push boundary still holds.

## Completed: Schur content operators — 2026-09-25

`sym-poly-sym::{c_n, c_n_inverse, schur_content_product}` implements the
user-defined eigenvalue `product_(i,j in lambda)(n+j-i)` with one-based box
coordinates. Content products use BigInt; forward application preserves the
coefficient ring, while the inverse accepts BigInt/BigRational and returns
BigRational. Both return Schur expansions. The inverse reports
`ContentOperatorError` on nonzero terms with more than n rows: finite-variable
users must first discard these vanishing Schur terms. Empty products and
constants are one/unchanged, including n=0. No signed narrowing of n occurs.
Existing non-Schur transition-table limits remain explicitly documented.

Eight new tests cover independent box products, the full usize variable-count
boundary, arbitrary-precision products, zero/constants, the exact inverse
domain after basis conversion, fractional coefficients, all six input bases,
both composition identities, mixed-degree linearity, and graph families through
degree 7. Full crate tests pass (232 passed, 3 pre-existing ignored), including
the new doctest. The checked `schur_content_operator` example, strict Rustdoc,
formatting, and diff checks pass. Focused Clippy reports no findings in new
files; unrelated existing warnings remain. Logs are
`/tmp/content-operator-{tests,doc,clippy}.log` and
`/tmp/content-operator-clippy.jsonl`.

All task files are released after the local checkpoint on
`feat/translation-invariant-graph-bases`. No Lean work, worker launches, or
publication occurred. The prior no-push boundary below still applies.

## Completed: forgotten-basis investigation — 2026-09-25

All-degree partition graphs have strictly positive forgotten coefficients
after multiplying by `(-1)^(|lambda|-length(lambda))`, using `f=omega(m)`.
`docs/partition-graph-forgotten.md` records the star-expansion and set-partition
inclusion-exclusion proof, normalization, examples, and exact evidence.
All 41 two-partitions through degree 10 were checked at three variable counts
(123 cases), with complete support, the predicted all-ones coefficient,
monomial roundtrips, and an independent Schur-conjugation conversion route.
Only this handoff and the mathematical note were changed; existing library
generators and conversions were reused by a bounded temporary Rust driver.
These note files are released. The prior no-push boundary below still applies.

## Completed: all-degree partition-graph basis — 2026-09-25

`sym-poly-sym::partition_graph_symmetric(&lambda, n)` implements the
all-degree basis from Alexandersson--Shapiro, arXiv:1104.0589 Section 2.2.
The user explicitly selected this family; the even-degree square-graph basis
is not included. The generator uses disjoint outward stars, BigInt coefficients,
explicit variable count, and Definition 1's full permutation sum. The reusable
`symmetrized_graph_monomial` expands via MultiPoly and exponent stabilizers,
without enumerating S_n; Schur conversion reuses the existing Sym API.
The checked `partition_graph_schur` example lists degrees 2--6 at n=2*d,
or accepts `degree [num_vars]`. README and library-guide routing are updated.

Eight new tests pass: literal permutation-sum agreement (including isolates,
repeated edges and disconnected stars), known even/odd Schur expansions,
all 2-partition diagonal coefficients and Schur roundtrips through degree 8,
translation invariance, n=30 BigInt coefficients, boundary/errors, orientation,
and relabeling. Full package tests pass (223 passed, 3 pre-existing ignored);
the subsequently strengthened nonzero translation probes also pass. The
example, strict Rustdoc, formatting, and diff checks pass. Focused Clippy has
no findings in the new module/tests/example; existing unrelated warnings remain.
Logs are `/tmp/partition-graph-{sym-tests,doc,clippy}.log` and
`/tmp/partition-graph-clippy.jsonl`.

All task files are released after the local checkpoint on
`feat/translation-invariant-graph-bases`, based on `62d6e6f`. No push was made:
this branch includes 13 pre-existing unpublished master commits, including a
prior handoff's explicitly unpublished work, which this task does not publish.
The i64 limitation in existing Schur transition data is documented and a
separate exactness follow-up is proposed in
`suggestions/sym-transition-exactness-2026-09-25.md`.
Paper-cache lookup found no copy; its arXiv fetch failed with HTTP 406 from
the metadata endpoint. The PDF was read through the browser instead.
No Lean work or additional worker was involved.

## Completed: colored exchange library — 2026-09-22

The user-authorized 3-Matroids integration adds `combinatoric-core::colored_exchange`:
validated fixed-content families, deterministic symmetric cycle witnesses,
explicit search-budget exhaustion, exact BigUint marginals, and verified
product/binary-port/private-block certificates. The API has no 128-word cap.
Its example demonstrates witnesses, obstructions, weighted coefficients,
and decompositions. README and library-guide routing are updated.

Claude implemented the module, export, 11 integration tests, and example.
Root reviewed and integrated it, fixed actual traversal-budget accounting,
removed finite sentinels from unbounded searches, compressed arbitrary raw
color labels, rejected unused private-certificate data, and guarded corrupt
binary count overflow. A twelfth integration test covers these defects.
All eight unit tests and twelve integration tests pass. Independent brute
force agrees on all 4,095 quota-(2,1,1) supports (405 successes); tests also
cover larger-than-128-word families and arbitrary-precision coefficients.
The example, strict Rustdoc, and focused formatting pass. Clippy reports
pre-existing findings outside the owned files (44 in the lib-test target,
38 duplicated from the library); none concern this addition.

All task files are released after the local checkpoint. No nested Claude
process remains active. No push/publication was authorized or performed.
The project consumer and further exact classification replay are recorded in
`/workspace/projects/3-Matroids/RUST_INTEGRATION.md` and its handoff.
No analytic log-concavity conclusion is built into the certificate API.

## Completed: A355089 statistic refinements — 2026-09-18

`combpoly::parity_pattern` now uses one generic maximum-and-gap traversal for
exact counts and five exact polynomial refinements: descents, ascents, major
index, ascent-major index, and left-to-right minima. The public
`parity_preserving_123_statistic_polynomials` API returns ascending-degree
`BigUint` rows, and the checked `a355089_statistics` example emits rows in
Polytool's input format. Independent permutation generation verifies every
refinement through `n=8`; specialization at `t=1` verifies all five through
`n=20`. Reverse-complement and the disjoint record-set decomposition prove
that the left-to-right-minimum polynomial is palindromic at every even length;
this has a focused regression through `n=20`. All 200 Combpoly library tests,
both examples, strict Rustdoc,
formatting, and focused strict Clippy pass.

Project-local Polytool searches found no scalar or even/odd-vector recurrence
in the recorded bounded spaces. The ignored
`experiments/src/bin/a355089_stat_recurrence.rs` remains a use-once search
driver, not maintained library functionality. No generic permutation or
Polytool source was changed. All files from this batch are released.

## Completed: parity-alternating 123-avoidance DP — 2026-09-18

`combpoly::parity_pattern` now provides an exact `BigUint`, count-only memoized
DP for both parity orientations of `321`-avoidance and the canonical A355089
generator obtained by value complement. The checked `a355089` example
reproduces generation. Its independent brute-force test covers `n <= 9`, all
published A355089 terms through `n=30` are regression-tested, and the preserving
orientation agrees with A354208 through `n=20`. All 197 Combpoly library tests,
the example target, strict Rustdoc, formatting, and focused Clippy pass. Strict
whole-crate Clippy initially reported only 14 pre-existing findings outside the
owned files; allowing precisely those existing lint classes exposed no new
finding. The ignored `experiments/src/bin/parity_123.rs` and generic
`combpoly/src/permutation.rs` were profiled/read but not modified. All files
from this task are released.

Full chronology through the start of the current repository-health batch is in
[`docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md`](docs/HANDOFF_ARCHIVE_THROUGH_2026-09-16.md).

## Active: repository health remediation — 2026-09-16

`/root` (2026-09-17) has released this handoff and the ignored project
consumer `experiments/src/bin/chordal_clawfree_boundary_scan.rs` after a
bounded higher-rank transfer run.  No maintained library source was owned or
changed.

The chordal-claw-free Schur project has released all Rust files.  Reusable
work now includes `Graph::join`, the monomial-basis universal-vertex operator
`chromatic_symmetric_add_universal_vertex`, and cached degree-dependent
permutations, cycle indices, and character tables in the exact clique-boundary
routine.  The universal-vertex API documents and tests its exact positive
action on every two-column Schur function.  The family-specific scan remains in the ignored
`experiments/src/bin/chordal_clawfree_boundary_scan.rs` consumer.  Its new
`positive-compression` mode refutes the naive model obtained by leaving the
new universal label unused in the old harmonic module; that model cannot
produce the lower outer two-row shapes present in the true transfer.  The
`algebra-cone` mode implements exact insertion of the universal choice in
each PEO transposition code and shows that applying this map to the naive
orthogonal primitive-idempotent section also gives the wrong refined kernel.
The further `orbit-compression` mode tests a diagonal projector onto
`c`-subsets containing exactly `c-b` boundary labels.  It matches all `29`
endpoint columns through order five and all `18` endpoint columns in the
order-six spanning samples, but fails for mixed boundary content already at
`(n,k,c,b)=(5,2,2,1)`: at `C=1` it gives `(1,9)/10`, rather than the exact
`(7,9)/16`.  This is project-specific finite evidence and a refuted general
model, not shared-library functionality.  During this checkpoint `/root`
temporarily owned only this handoff and the ignored experiment driver; both
are released.  The failed compression models are project-specific refuted
routes.  Its
`cone-transform` mode reconstructs exact rational candidate transfer matrices:
all fourteen matrices through order five have full input rank, and every
negative entry has an outer source shape of first row at least three.  The
same statement now holds for the full-rank order-seven spanning
reconstructions at boundary sizes two, three, and four (`28`, `19`, and `10`
input coordinates, respectively).  The ignored driver permits these bounded
order-seven runs; this project-specific cap change is not maintained-library
functionality.  The
order-six unlabeled universal-clique audit covered 17 relevant graphs and 361
positive boundary states without a failure; the attempted order-seven audit
was interrupted before completion and has no claimed result.

Verification passed: 269 `combinatoric-core` tests, 135 `sym-poly-sym` unit
tests plus all integration and documentation tests, strict Rustdoc, the
checked universal-vertex example, and the project-specific exact transform
sweeps.  Focused strict Clippy reports only the 38 pre-existing
`combinatoric-core` findings and 36 pre-existing `sym-poly-sym` findings; the
changed source files add none.

The chordal-claw-free Schur project has again released its ignored experiment
driver.  `experiments/src/bin/chordal_clawfree_boundary_scan.rs` now has a
`detail` mode for individual channel tables and a `universal` mode for the
universal-vertex preservation probe.  The latter finds the unrestricted claim
false (the edgeless three-vertex graph cones to a claw), but finds no failure
among 3,181 boundary-positive chordal claw-free inputs through order five when
the cone also remains claw-free.  These are project-specific proof-search
checks and remain outside the maintained workspace.

The chordal-claw-free Schur project has released all Rust file ownership after
promoting its reusable core.  `combinatoric-core` now provides the complete
`3`-sun clique-blowup generator and clique-rooted perfect elimination orders;
`sym-poly-sym` now provides exact cleared chordal clique-boundary Schur
channels through order nine.  `chromatic_symmetric` uses one stable-set-
partition traversal shared by all monomial types.  The checked example and
library guide expose the new boundary API.  The family-specific scan remains
an ignored driver at
`experiments/src/bin/chordal_clawfree_boundary_scan.rs`, outside the maintained
workspace.

Verification passed: all 268 `combinatoric-core` tests, all 134
`sym-poly-sym` unit tests plus its integration and documentation tests, strict
Rustdoc for `sym-poly-sym`, the checked boundary example, and the explicit
negative-intermediate/positive-completion cross-check.  Focused Clippy adds no
warning in the new module; its 36 reports are pre-existing elsewhere in
`sym-poly-sym`.  The project scan found no negative outer-private-clique
channel through total order eight (28 blowups, 243 nonzero coefficients),
whereas the overlapping central boundary is already negative for the bare
`3`-sun.  This sharply identifies the outer private clique as the viable
interface for the next proof.
The shared-library checkpoint is committed locally.  Its routine push to the
existing SSH remote failed because this environment has no usable public-key
authentication; no remote or credential configuration was changed.

The host supervisor owns final integration, CI observation, documentation,
suggestions, and this handoff. All bounded source workers have completed and
released their files; there is no active source ownership in this repository.
Commits `4a897bf`, `e777c57`, and `f7d136f` are pushed on `master`.

Kostka retirement is complete in standalone commit `95cc655` and workspace
commits `b4f1174` and `06340f3`; GitHub Actions run `35095705694` passed. The
historical `kostka` submodule and database-backed `KTT-search` application remain in Git
for reproducibility, but are excluded from the maintained root workspace and
have no dedicated CI job. New GT/Ehrhart work routes to Ehrcalc;
LR/Kostka/Schur compatibility work routes to lrcalc-rs. The Kostka repository
README and GitHub description identify both successors. The tracked-experiments
CI job still checks out the frozen submodule because four retained GT
experiment binaries consume it.

The final performance follow-up streams Combpoly permutation statistics,
flagged-tableau weights, QSym quasi-shuffle multiplicities, and both general
and chordal graph sink-polynomial independent sets. Noncrossing and nonnesting
matchings now prune invalid branches directly. Full affected-package tests pass:
266 combinatoric-core, 194 Combpoly library plus 2 binary, 204 multipoly, and
136 QSym tests, together with all affected examples.

The 961 ignored experiment binaries remain local and untracked. Broad builds
must not compete with the active eight-process KTT computation; focused checks
use external targets, reduced priority, and separate CPU affinity when safe.

The root workspace now excludes the standalone `experiments` and retired
`KTT-search` workspaces, each with a committed lockfile. The Combpoly inner Git
metadata was preserved as a verified bundle and recoverable metadata directory
under `/mnt/2TB-Babel/ai-storage/source/repo-metadata-backups/` before removal
from the working tree. The nested Kostka metadata-only checkpoint is
`a4f01fe`, pushed to `origin/main`.

The remaining deliberate debt is narrow: add durable performance benchmarks,
resolve the 38 pre-existing strict-Clippy findings in combinatoric-core, and
upgrade the isolated experiments/KTT MySQL chain when its upstream
future-incompatibility warning is resolved. See
[`suggestions/rust-repository-health-2026-09-16.md`](suggestions/rust-repository-health-2026-09-16.md).

The complete workspace documentation builds with `RUSTDOCFLAGS="-D warnings"`;
the root CI enforces this alongside the declared Polytool MSRVs. Kostka's two
documentation repairs are pushed separately as `ec1383c`.

The first MSRV run found post-1.82 `is_multiple_of` calls in Polytool; all
library, CLI, and example occurrences are now 1.82-compatible. Strict Clippy
also passes for Polytool and its MCP companion. The full 356-test Polytool
library suite passes after those changes, as does an all-target check on the
exact Rust 1.82.0 toolchain for Polytool and its web companion.

Polytool's `oeis info <A-number> --json` now embeds `recurrence_data`, the
standard `polytool.recurrence.v1` object with exact coefficients and initial
polynomials, while preserving the existing human-readable recurrence and
language exports. Focused CLI tests, recurrence JSON fixture replay, formatting,
and strict binary Clippy pass.

The OEIS catalog's former 30 experimental entries are now reconciled. Nineteen
queue recurrences with complete OEIS-row alignment are correctly `validated`;
the corrected A099040 import is also validated. Ten entries that cannot align
to complete OEIS rows are no longer bundled: A103328 has reversed structural
zeros and already has its valid recurrence in OEIS, seven are stale row-boundary
artifacts, and A266178/A266298 lost cellular-automaton zero positions. The
generator now rejects unaligned queue and Lean imports. Regeneration also drops
the corrected-away A062154 entry. The catalog contains 774 entries (125
verified, 649 validated, zero experimental). Generator checks, 9 focused OEIS
library tests, 9 CLI tests, all 31 MCP tests, formatting, and strict all-target
Clippy pass.

## Durable policy

- `master` is the canonical monorepo branch; `polytool/main` is generated only
  by `scripts/sync-polytool-main.sh` after its prerequisites pass.
- Use root Cargo commands for maintained crates and explicit `--manifest-path`
  commands for experiments or legacy KTT.
- Keep ignored experiments disposable. Promote only stable reusable kernels
  into documented and tested library APIs.
- Keep build output in `/cargo-target/ai-projects`, not this source tree.
- Preserve unrelated worker changes and record new ownership here before
  adopting a dirty file.
