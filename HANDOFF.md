# Handoff

## Completed: explicit cycle certificates and normalization obstruction — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No shared library,
Lean/paper source, or circular code is owned.  Implementation checkpoint:
`528410f` (`Certify Hall bump cycles and code obstruction`).

Every inductive bridge step now records a deterministic minimum decomposition
of the matrix bump into zero, one, or two simple alternating cycles, the
intermediate bumped matrix, and the exact maximum-color row deletion.  The
test reconstructs the intermediate matrix by applying the signed cycle
incidences, performs the recorded deletion, and compares the result with the
next certified residual object.  Hence cycle validity is now explicit rather
than re-established by enumerating all radius-two candidates.  Across all
rank-at-most-five certificates (tails counted in every containing
certificate), source cycle counts are `{0:4854, 1:492, 2:63}` and target
counts are `{0:4996, 1:345, 2:68}`.  The first selected two-cycle source step
is `mu=(2,2)`, word `2134`: toggle
`-r1c2,-r2c1,+r1c1,+r2c2`, then
`-r1c3,-r2c4,+r1c4,+r2c3`, then delete color 4 from row 1.

The signed drops cannot be converted into an ordinary bounded Carlitz code by
adding a potential depending only on residual row shape.  Nonnegativity alone
is feasible through rank five, but the labels range from 0 to 6.  The sharper
bounds `0 <= label <= rank-1` hold through rank four and fail first at rank
five.  The exact contradictory two-edge constraint uses the same shape move
`(2,1,1,1) -> (1,1,1,1)`: one certified transition has drop 0, forcing the
potential difference to give a label at most 4, while another has drop -6,
forcing it to be at least 6.  Concrete witnesses are `mu=(4,1)`,
`11112 -> 11112` (drop 0, identity bump) and
`43214 -> 43214` (drop -6, switch
`-r2c4,-r4c2,+r2c2,+r4c4`).  Allowing the potential also to depend on packed
color content still has an explicit negative constraint cycle, so the bump
choice itself carries essential code information.  The next clean-rule search
should therefore classify the recorded zero/one/two-cycle moves rather than
try another shape-only normalization.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 24 passed, 38.79s
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Completed: inductive two-cycle Hall bridge — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No shared library,
Lean/paper source, or circular code is owned.  Implementation checkpoint:
`b0b80d0` (`Make Hall bump bridge inductive`).

The endpoint code bridge below was not yet genuinely recursive.  Its first
paired-trace shape mismatch is `mu=(3,1)`, source `1323`, at depth one:
source residual shape `(1,1,1)` versus target residual shape `(2,1)`.  Adding
the residual partition to every code step makes this obstruction visible.
With at most one alternating cycle per deletion, the shape-enriched graph has
perfect matchings in 106/107 fixed-content buckets through rank five; the
first exact failure is `mu=(3,2)`, squarefree source `12345`.  Allowing two
successive simple alternating cycles restores all 107/107 buckets.

Choosing these shape-enriched matchings independently at every rank is still
not inductive: the first mismatch is `mu=(3,1)`, source `3312`, whose chosen
rank-three residual pair differs from the separately chosen rank-three map.
`HallRecursiveBumpBridge` now fixes this by induction.  A rank-`n` edge is
admissible only when its source and target first-step drops and residual row
shapes agree *and* the residual pair is exactly an edge of the already chosen
rank-`n-1` bridge.  Deterministic perfect matchings still exist for every one
of the 19 Hall diagrams and all 1119 packed strict colorings through rank
five.  Thus the resulting endpoint map factors through its recorded smaller
objects at every depth.

Each certificate contains the target coloring, the common sequence
`(signed degree drop, residual row shape)`, and the complete source and target
matrix traces.  Exhaustive checks verify forward/inverse recovery, exact
packed color content, strictness, valid at-most-two-cycle transitions, one
incidence deleted per step, recursive agreement with the independently built
residual bridge, and
`asc(source)-shift = sum(drops) = charge(target)`.  The CLI prints the
certificate and readable traces.  This is now a genuine finite recursive
bijection on the tested Hall--Littlewood subset, but the per-content perfect
matching and enumeration of all radius-two cycle moves remain global search;
the next research target is a closed `PullUp`-style selector producing the
same inductive edge locally.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 24 passed
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Superseded checkpoint: endpoint one-cycle Hall bump-code bridge — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No shared library,
Lean/paper source, or circular code is owned.

The experiment now constructs an actual reversible Hall statistic-switching
bridge, not just a generating-function comparison.  At each rank it deletes
the current maximum incidence; when that maximum occurs in a longest
non-singleton strict-chain row, it may first apply one simple alternating
incidence cycle.  Recursing this candidate operation on both the source
ascent-minus-shift statistic and the target charge statistic produces sets of
signed drop codes.  In every fixed-content bucket through rank five, the
bipartite graph joining states with a common code has a perfect matching:
107/107 buckets and all 1119 packed strict colorings.  Allowing two cycles at
one deletion adds no matching power in this range.

`HallBumpCodeBridge` chooses a deterministic matching, records the target and
the common signed drop code, and stores the inverse certificate.  The code can
be replayed deterministically as a source and target matrix-deletion trace.
The exhaustive test checks inverse recovery, strictness, exact packed color
content, trace length and rank decrease, and
`asc(source)-shift = sum(code) = charge(target)` for every partition through
rank five.  The report includes the first nontrivial readable witness
`mu=(3,1): 3124 -> 2314`, code `[-2,2,1,0]`, with source trace
`[24|3|1] -> [1|2|3] -> [1|2] -> [1] -> []` and target trace
`[14|2|3] -> [1|2|3] -> [1|2] -> [1] -> []`.

Two qualifications are durable.  Fixing the target to the earlier naive
deletion code still fails: only 32/38 rank-at-most-four buckets match, first at
`mu=(3,1)` with isolated source `4213`; bumping must occur on both sides.
Also, this is a canonical finite matching among recursively local code
choices, not yet a closed Gillespie-style `PullUp` selector.  Negative entries
such as the first `-2` above show that these raw degree drops are not yet a
standard nonnegative Carlitz code.  The next clean-rule search should inspect
the selected paired traces and replace the bucket matching by a deterministic
local choice, or isolate its first obstruction.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 24 passed
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Completed: recursive deletion-code obstruction search — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No shared library,
Lean/paper source, or circular code is owned.  The `--hall-code-report` probe
now compares canonical colored-incidence deletion codes, adaptive degree-drop
standardization, occurrence-level charge labels, and blockwise Foata maps.

Deleting the largest color from the last eligible strict-chain row and
recording the successive statistic drops is an exact source/target code for
every rank-three Hall diagram (all 13 fixed-content buckets).  The first
failure is the complete graph `mu=(4)`, content `(1,1,2)`: the unique bad word
is `2313`, whose naive source code is `[2,1,0,0]`; the proved complete-graph
Foata correction sends it to `2331`, whose target code is `[1,2,0,0]`.  Thus
the first repair is exactly the required weakly increasing carry among two
equal largest letters.  Choosing repeated maximum occurrences by minimum
current degree drop implements that carry.

The adaptive rule then reaches the first genuinely non-Carlitz obstruction at
`mu=(3,1)`, squarefree content: `3124` and `4123` both acquire code
`[0,1,0,0]`.  Their first deletions have equal statistic drop zero but leave
different row-margin partitions, `(3)` and `(2,1)`.  Adding the residual shape
to every code entry restores source injectivity in 100 of 107 buckets through
rank five, but source and target enriched-code sets already disagree for
`mu=(2,1)`.  Therefore the shape is not extra common-code data: the deletion
transition itself must perform a reversible bump to the canonical corner.
This matches the role of Gillespie's `psi`/ZeroBump operation and isolates the
next implementation target: a matrix-row bumping sequence before deletion.

Two tempting alternatives are now ruled out precisely.  Raw charge labels are
not injective even in small complete graphs.  Applying the proved Foata map
independently inside clique blocks ignores essential cross-block data: through
rank five only 1023/1119 outputs remain strict and only 743/1023 of those have
the right degree; the first degree failure is already `mu=(2,1)`, `123`.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 23 passed
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Completed: alternating-cycle correction census — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No shared library,
Lean/paper source, or circular code is owned.  The experiment now measures
fixed-margin corrections intrinsically in the strict-chain/color incidence
matrix, requiring every already degree-correct state to remain fixed.

A single alternating cycle is exactly enough for every `Gamma_(3,2)` state:
among its 101 packed strict colorings, the deterministic matching fixes all 29
raw matches and uses 55 two-chain switches and 17 genuine three-chain cycles.
For the complete at-most-three-color subfamily (23 states), an explicit
reversible rule is recorded and checked: four transpositions
`11223<->21123`, `12132<->31122`, `11233<->31123`,
`12133<->31132`, one directed three-cycle
`12233->31223->32132->12233`, and all other states fixed.  It preserves exact
content and strictness, and raw dual RSK after correction has the original
ascent degree.

One cycle is not universal.  The first forced obstruction is already
`mu=(3,1)`, squarefree content, degree three: the singleton bucket forces
`1234 -> 3412`, whose incidence difference has distance four and decomposes
into two simple alternating cycles.  Exhaustively over all 19 partitions and
all 1119 packed strict colorings through rank five, however, a deterministic
degree-correcting matching that fixes every raw match always exists using at
most **two successive simple alternating cycles**.  This is a useful bounded
recursive target, not yet a uniform local selection rule.  In particular it
rules out both a universal involution and a universal one-cycle recipe; the
next conceptual step remains a Carlitz/Foata-style deletion code that chooses
the first cycle and carries enough state to reconstruct the second.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 22 passed
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Completed: Hall--Littlewood matrix/Foata shortcut — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`; checkpoint `8f46480`
contains the implementation and tests.  No shared library, Lean/paper source,
or circular code is owned.  This checkpoint isolates the exact
statistic-switching problem after finding both a clean complete-graph rule and
the natural ungraded rule for every Hall--Littlewood diagram.

For complete graphs there is now a direct, reversible map on every packed
coloring, not just distinct colors.  For a packed word `w`, reverse it, apply
the inverse generalized Foata transformation, and RSK the resulting positive
word as `(P,Q)`; the output is the color SSYT `P` together with charge tableau
`transpose(evacuation(Q))`.  The inverse applies inverse RSK, generalized
Foata, and reversal.  It preserves exact content and satisfies
`asc_Kn(w)=charge(T)`.  The inverse Foata transform is now algorithmic rather
than a bounded enumeration.  Exhaustive tests through rank six prove forward
and inverse validity, injectivity, and (through rank five) equality with the
independently enumerated full target set.  A readable example is available as
`--foata-trace-coloring 1,1,2`.

For general `Gamma_mu`, the strict-edge graph is a disjoint union of paths of
lengths `mu'`.  A strict coloring is therefore exactly a binary matrix with
rows=strict chains and columns=colors.  Dual RSK (sort by color increasingly,
chain label decreasingly on ties; row-insert chain labels; transpose the
recording tableau) is an exact ungraded bijection to all pairs `(T,P)` with
`content(T)=mu'`, conjugate shapes, and the exact packed color content.  The
implementation independently checks the entire source and target sets for all
19 partitions through rank five: 1119 objects, no collision or omission.

The remaining issue is precisely a statistic switch on fixed-margin binary
matrices.  Raw dual RSK first fails for `mu=(3)`, coloring `112`: source ascent
degree 2, but raw charge 1.  Minimum-distance degree-correcting matchings exist
in every fixed content bucket; their maximum fixed-margin matrix distance grows
to five through rank five, so a fixed local switch is not the rule.  This is the
maximal-slice Macdonald `q,t` switch underlying the brief proof of Proposition
37, not an omitted RSK convention: horizontally reflecting the corresponding
right-justified `Gamma_mu` filling gives shape `mu`, maximal major index, and
`inv=asc-A`; the desired correction sends max-major fillings of `mu` to
max-inversion fillings of `mu'` while preserving content.  Gillespie,
arXiv:1503.02109v2, supplies arbitrary-content hook maps and squarefree maps,
but explicitly leaves the general arbitrary-content switch open (Conjecture
15 and Question 64).  Its repeated-content zero-slice results for two/three
rows do not imply this maximal slice; even shape `(2,2)`, content `(2,2)` is
outside that deduction.  Thus forcing a purported universal local rule would
hide a genuinely unresolved bijection problem.

The first nonhook case nevertheless has a clean rule.  For `mu=(2,2)`, write
the two strict-chain color sets as `S1={c1,c3}` and `S2={c2,c4}`.  If all four
colors are distinct `a<b<c<d`, with `a` in `S1` and `d` in `S2`, exchange the
middle colors `b,c`; otherwise fix the matrix.  Equivalently this toggles
`(a,c,b,d) <-> (a,b,c,d)`.  This is an involution preserving both margins, and
exhaustive packed tests prove that dual RSK after the switch has the source
degree for every `Gamma_(2,2)` coloring.  A simple involution cannot be the
general pattern: the joint table of `(source degree, raw target degree)` is
already asymmetric for `mu=(4)`, content `[1,1,2]`, and for `mu=(3,2)` even in
squarefree content.  A Carlitz/Foata-style recursive code, rather than pairwise
swapping all discrepancies, is therefore the more plausible extension.

The standardized Pieri-switch audit gives the same warning from the HRV side:
only 36 of 90 current Hall survivors through rank five are already aligned,
and reserving all aligned survivors is impossible first at `mu=(3)`.  The
`--hall-switch-report` command now records this obstruction instead of
aborting.  `--hall-matrix-report --max-n 5` prints the matrix audit, degree
errors, and minimal correction radii.

Current verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                         # 18 passed
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

## Completed: Hall--Littlewood charge bridge and charge repair — 2026-10-16

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs`,
`sym-poly/core/src/tableau.rs`, and this handoff on branch
`exp/verticalstripllt-marked-coloring-bijection`.  No Lean/paper source,
circular code, or other maintained API changed; the existing Rust branch was
preserved and was not pushed.  Implementation checkpoint `543e1de` contains
the charge repair, Hall--Littlewood bridge, tests, and initial handoff record.

The exact-current cached source is arXiv `1903.03998v3`, tree SHA-256
`32e82b0272e790383b66d9b0f7b0acc1727b58e4b86746e92be1a9f2cb6001b0`, at
`/papers/math-co/arxiv-extract/arxiv-sources/1903.03998v3`.  Definition 36 has
two literal endpoint errors: its clique intervals overlap and the displayed
triangular bound produces too many edges.  The implementation uses the unique
interpretation agreeing with the stated `binom(mu_j,2)` count and the printed
`Gamma_432` figure: disjoint clique blocks, cross-block offsets
`i+k <= mu_j-2`, and all outer corners strict.  The figure's normal-area
sequence `[0,1,2,3,2,2,2,1,1]`, complete-graph cases, and strict-chain cases
are regression-tested.

`--hall-littlewood-report --max-n 5` audits all 19 partitions through rank
five.  In every case the final positive HRV survivor multiplicities agree
exactly, degree by degree and conjugate shape by conjugate shape, with
`q^A K_{nu,mu'}(q)` computed from charge tableaux.  A reversible
`HallLittlewoodBridge` now pairs survivors with those tableaux inside each
exact bucket while leaving the arbitrary insertion tableau `P` unchanged.
Thus the original coloring polynomial `G_{Gamma_mu}` maps to pairs `(T,P)`
where `T` has content `mu'`, `P` has shape conjugate to `shape(T)`, and the
degree is `A+charge(T)`.  This is a genuine canonical-rank certificate, but is
not yet the desired local HRV-to-charge shortcut.  `--hall-partition 2,2`
prints a human-readable forward/reverse bridge.

The restricted family simplifies but does not eliminate tableau transport.
Exact `Q` fails first for `mu=(2,2)` and also for `(3,2)` and `(2,2,1)` through
rank five.  For `(2,2)`, the forced degree-one cancellation is the one-edge
toggle `eta: {12,23}->{23}`, with fixed `B={23}`, together with the elementary
dual-equivalence move `d_2(Q)`; neither mark transfer nor triangle surgery is
needed.  Charge extraction for content `mu'` has successive subword lengths
`mu_1,mu_2,...`, exactly the clique-block sizes, making Pieri switching between
these rounds the strongest current shortcut candidate.  Fixed strict-chain
RSK already fails on `mu=(2,2)`, coloring `(1,3,2,4)`: the two chain orders
land in degrees three and one rather than the required degree two.

The audit exposed a real shared-library bug: `Tableau::charge()` discarded the
positions selected into every standard subword and computed
`maj(inverse(reverse(w)))` instead of SymCat's
`maj(reverse(inverse(w)))`.  `word_charge` and `word_cocharge` now preserve
subword order, use the correct charge labels, and normalize nonpartition
content via RSK plus the existing Lascoux--Schuetzenberger involutions.  Exact
tests include row/column tableaux, SymCat's four-subword charge-eight example,
cocharge complements, and nonpartition Sage reference values.  The experiment
now uses this shared API; its duplicate workaround was removed.

Verification with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --locked -p sym-poly-core                         # 133 passed
cargo test --release --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection                       # 15 passed
cargo clippy --locked -p sym-poly-core --no-deps -- -D warnings \
  -A clippy::suspicious_arithmetic_impl -A clippy::needless_range_loop \
  -A clippy::while_let_loop -A clippy::type_complexity \
  -A clippy::new_without_default
cargo clippy --locked --manifest-path experiments/Cargo.toml \
  --bin marked_coloring_hrv_bijection --no-deps -- -D warnings \
  -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

Unqualified strict Clippy remains blocked by pre-existing warnings in
`combinatoric-core` and unrelated `sym-poly-core` files; the scoped owned-file
gates above are clean.

## Completed: fixed-HRV Boolean cancellation and residual census — 2026-10-15

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  No maintained API,
Lean/paper source, circular code, or external action changed.  An independent
Astra review derived the fixed-HRV row factorization and checked the K4
triangle obstruction; Astra changed no files.

The actual Schur involution now has a proved first stage.  At fixed HRV map
`r`, retained marks `B`, and tableaux `P,Q`, it scans rows and toggles the first
unmarked redundant normal edge.  Its critical states have every lower-HRV edge
marked; forced rows have every tied-HRV edge marked; and each unforced row has
either every tied edge marked or only its first tied edge present and unmarked.
The implementation checks the partner has the same `r,B,Q,P`, reverses sign,
and toggles back.  Exhaustively through rank five, this cancels 142,348 pairs
and compresses 308,351 skeletons to 13,258 positive plus 10,397 negative
critical states while preserving every signed `(q-degree,shape)` bucket.

The remaining critical states are still matched deterministically within
combined `(eta,B)` radius five.  Residual radius one suffices through rank
three, while radius five is necessary at rank four and sufficient through all
197 rank-five paths.  A toggle-plus-triangle graph is not yet enough: it fails
on 8 rank-four and 76 rank-five paths.  The residual census has 10,397 pairs;
4,857 preserve `Q`, while the rest require dual-equivalence distance 1--4.
Thus the next conceptual step is a Pieri/tableau-switching rule that transports
`Q` across an HRV fibre split/merge or triangle surgery.

The earlier “exact-Q K3 obstruction” was really an exact-`(B,Q)` obstruction.
Allowing `B` to move repairs complete graphs, but literal fixed `Q` is still
impossible in general.  The first counterexample is path `NNDNEEE`, q-degree
one, shape `[2,1,1]`, `Q=[1 3]/[2]/[4]`, where the fixed-Q group has one
positive and two negative critical skeletons.  This is regression-tested.

The two-stage `J` retains the full end-to-end audit: 48,840 packed colorings
bijection with 48,840 surviving SSYT objects through rank five, with exact
content and q-degree and replayable forward/inverse certificates.  Transported
`S_i` pass through rank four (4,023 involutions, 631 distant commutations,
2,163 braid relations).  Thirteen focused tests pass.  Useful commands:

```bash
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --boolean-report --max-n 5
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --reduced-move-report --max-n 5
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --exact-q-report --path NNDNEEE
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --tableau-report --path NNNEEE --trace-coloring 1,2,1
```

## Completed: bounded-radius repair of Schur cancellation — 2026-10-15

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  Checkpoint `42daa4f`
contains the repair.  The actual `J`
cancellation used by the coloring-to-SSYT map no longer rank-matches arbitrary
states in an entire `(q-degree, shape)` bucket.  It constructs the bipartite
graph of opposite-sign skeletons at combined Hamming distance at most five in
`(eta,B)` and uses a deterministic augmenting-path matching.  It may change
`Q` within the same shape but always leaves `P` fixed, so crystal equivariance
is unchanged.  No maintained API, Lean/paper source, or circular code changed.

The data isolate why weaker repairs cannot work.  Exact `(B,Q)` already has a
K3 obstruction: for `B={12}`, `Q=[1 2]/[3]`, the sole state is the negative
`eta={12,13}`.  Allowing one elementary dual-equivalence move fixes K3 but
fails first on K4 at `B={12}`, `eta={12,34}`,
`Q=[1 3]/[2 4]`.  Even arbitrary `Q` of the same shape and every odd
orientation distance cannot match this K4 group while literal `B` is fixed.
Thus a valid rule must move retained q-marks as well as orientation edges.

When only `|B|`, shape, and `P` are fixed, a matching of combined `(eta,B)`
radius one suffices through rank three; radius five suffices for every one of
the 45 rank-four and 197 rank-five paths.  Radius three genuinely fails on a
K4 group.  A simple lexicographic greedy pairing needs radius 13 at rank five,
so the augmenting-path rerouting is currently essential; finding an explicit
local description of that rerouting remains the next conceptual step.

With radius-five `J` installed in the final map, the independent rank-five
audit again gives 46,957 colorings on each side and 2,556 positive skeletons;
all ranks through five give 48,840 objects and 2,861 survivors.  The new
alternating paths total 1,385,562 steps and have maximum length 18,258.
Transported `S_i` again pass through rank four: 4,023 involution/weight checks,
631 distant commutations, and 2,163 braid relations.  Eleven focused tests,
strict Clippy, rustfmt, and diff checks pass.  Reproduction:

```bash
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --local-schur-report --max-n 5
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --schur-report --skip-si --max-n 5
```

## Completed: coloring-to-SSYT involution-principle prototype — 2026-10-15

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  Checkpoint `07ab9ff`
contains the implementation.  No maintained crate
API, Lean/paper source, circular code, or external action was changed.

The new `--tableau-report` experiment expands `q -> q-1`, orders HRV fibres by
increasing sink, reads each increasing fibre bottom-to-top as a column-crystal
word, and applies reversible RSK.  A signed e-object is stored exactly as
`(eta, B_eta, P, Q)`: `Q` and `eta` recover the ordered column boundaries, so
inverse RSK reconstructs the HRV filling.  The subset `B` is transported by
rank between the marked source ascents and the ascending edges of `eta`.

The obvious source sign involution `I` toggles the first ascent outside `B`.
Its transported version is not crystal-compatible: the first obstruction is
already K2, where it pairs the `e_2` column shape with the row contribution in
`e_1 e_1`.  Across ranks through four, it preserves insertion shape in only
5,862 of 11,456 nonfixed signed states and commutes in 23,488 of 36,186
nonzero crystal transitions.  This failure is retained as a regression.

At checkpoint `07ab9ff`, the initial alternate rule `J` was a canonical global
rank matching inside each `(q-degree, insertion shape)` bucket.  It left `P`
unchanged and supported the involution-principle map, but was not local.  The
bounded-radius checkpoint above supersedes that matching while retaining the
same positive multiplicities and reversible coloring-to-SSYT construction.

`--schur-report --skip-si --max-n 5` exhaustively checks all 258 ordinary
decorated paths: 308,351 signed LR skeletons cancel to 2,861 positive
survivors, and 48,840 packed strict colorings biject with 48,840 independently
enumerated surviving SSYT objects.  Rank five alone has 304,692 skeletons,
2,556 survivors, and 46,957 objects on each side.  Exact content and q-degree,
forward/inverse paths, injectivity, and independent surjectivity all pass; the
longest alternating certificate has 15,886 steps.

Crystal reflections/Bender--Knuth operators `S_i` are defined on `P` and
conjugated through the final bijection.  Exhaustively through rank four, 4,023
involutivity/weight/q-degree checks, 631 distant commutations, and 2,163 braid
relations pass.  The all-path rank-five bijection was checked separately;
repeating all Coxeter words there was stopped as redundant after the runtime
became disproportionate.  Ten focused tests, strict owned-binary Clippy,
rustfmt, and diff checks pass.  Useful commands:

```bash
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection --schur-report --skip-si --max-n 5
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection --tableau-report --max-n 4
/cargo-target/ai-projects/release/marked_coloring_hrv_bijection \
  --tableau-report --path NNEE --trace-coloring 2,1
```

## Completed: direct-HRV shortcut data audit — 2026-10-15

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on
branch `exp/verticalstripllt-marked-coloring-bijection`.  Checkpoint `d67dd18`
adds `--shortcut-report`; ordinary traces now print the one-based HRV map and
the actual vertex set of every fibre.  No maintained API, Lean/paper source,
or circular code was changed.

The exhaustive report through rank five covers all 277,302 source states.
Orienting exactly the marked edges gives the exact output HRV map in 175,695
states, the same sink set in 189,405, and the same fibre-size partition in
181,559.  Thus none is a direct rule.  Of 6,954 `(path,A)` groups, 3,924 have
multiple output HRV maps as the coloring varies, so a shortcut reproducing
this recursive map must use the coloring data.  The first failure is K3:
coloring `112`, marks `{13,23}` has marked
reachability `1->3,2->3,3->3`, while the bijection gives eta `{12,13}` and HRV
`1->3,2->2,3->3`; its HRV vertices are 2 and 3.

Keeping the source word fixed is also impossible, independently of the chosen
prototype matching.  For K4 (`NNNNEEEE`), word `1123` and q-degree 3 give 10
marked-coloring sources but only 8 HRV orientations on which that same word is
increasing along every fibre.  Across rank five, 6,244 of 127,284 fixed-word
buckets fail the count comparison.  Hence a direct reversible shortcut cannot
keep the coloring word fixed: it must permit colors to move among
vertices/fibres.  The simple greedy increasing-subsequence partition is much
weaker (exact fibre sets in 32,236 states).

Reproduction:

```bash
/cargo-target/ai-projects/debug/marked_coloring_hrv_bijection --shortcut-report --max-n 5
/cargo-target/ai-projects/debug/marked_coloring_hrv_bijection \
  --path NNNEEE --trace-coloring 1,1,2 --trace-mask 6
```

Seven focused tests, strict owned-binary Clippy, rustfmt, and diff checks pass.
The shortcut census and its minimal K3/K4 obstructions are regression-tested.

## Completed: recursive ordinary marked-coloring/HRV bijection — 2026-10-15

`/root` completed and releases
`experiments/src/bin/marked_coloring_hrv_bijection.rs` and this handoff on the
isolated branch `exp/verticalstripllt-marked-coloring-bijection`.  Checkpoints
`775bad7` and `0355ab0` replace the original whole-path bounce fallback.  No
maintained Rust API, Lean source, paper source, circular code, or external
action is owned.  The shared Rust checkout and its pre-existing branch remain
untouched; no push was made, preserving the recorded no-push boundary.

The ordinary prototype is now recursive at every node and has no global
source-to-target fallback.  Single NN and DN bounces port the paper's coloring
involutions and the explicit `bounceMatch` / `secondMatch` correspondences from
`BounceMatching.lean` and `SecondBounceMatching.lean`, including their
inverses.  Actual labelled column colors are transported using the identity or
adjacent-swap HRV conjugacies proved in `BounceTransport.lean` and
`SecondBounceTransport.lean`: a fibre's color set is carried across and placed
in increasing vertex order.  This preserves the ordered packed content, not
only the fibre-size partition.

For generalized NN/DN/ND bounces, the scanner follows `BounceTrace.lean` and
the three smaller-successor recurrences in `GeneralizedFirstBounce.lean`,
`GeneralizedSecondBounceReduction.lean`, and
`GeneralizedThirdBounceReduction.lean`.  The published many-bounce argument is
an algebraic cancellation proof and supplies no local object map.  The
prototype therefore uses the requested alternate rule: canonical rank
matching inside exact `(ordered content, q-degree)` buckets only between the
parent and its smaller recursion-branch objects, recurses on the chosen child,
then uses an independent target-side branch matching.  Branch kind, q-offset,
rank, and inverse replay are recorded in the trace.  An unsupported bounce is
now an error rather than a silent whole-path match.

Independent enumerators again verify source/target validity, exact colors and
degree, injectivity, surjectivity, and inverse trace replay through rank five:
258 ordinary decorated paths, 9,271 exact buckets, and 277,302 certificates.
Residual unsupported/global-fallback count is zero.  Across certificates,
175,432 traverse an explicit first bounce, 29,476 an explicit second bounce,
and 78,076 a generalized recursion-branch match (categories can overlap).
Dedicated traces/tests cover K3, the second-bounce path `NDNDEE`, and the third
bounce path `NNDDEE`.  The exact K3 obstruction to `eta = A` remains: for
content `[1,1,1]`, degree two, and orientation ascents `{12,13}`, there are two
sources but three targets.

Checks pass with `CARGO_TARGET_DIR=/cargo-target/ai-projects`:

```bash
cargo test --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection
cargo build --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection
/cargo-target/ai-projects/debug/marked_coloring_hrv_bijection --max-n 5
cargo clippy --locked --manifest-path experiments/Cargo.toml --bin marked_coloring_hrv_bijection -- -D warnings -A clippy::manual_is_multiple_of
cargo fmt --all --check
git diff --check
```

The single allowed Clippy lint is the pre-existing unrelated
`principal_stratum.rs` diagnostic documented below; the owned binary is clean.

## Earlier prototype checkpoint (superseded) — 2026-10-15

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
