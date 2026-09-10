# Shared sparse integer homology: implementation plan

Date: 2026-09-10. Status: planned; no implementation started.

## Objective and boundaries

Compute the reduced integral homology of Boris Shapiro's principal-stratum
complexes, with reusable sparse matrices, integer chain cancellation and Smith
normal form in the shared library. Keep the composition-specific model and
research conjectures outside the generic algebra implementation.

The source is the forwarded email UID 32 and its attachment
principal_poset_resonance_v6_stronger.tex. Only TeX and PDF were attached: no
original program, timings or cancellation certificates. Its homology outputs
are reproduction targets, not independently verified facts. Do not email,
publish the manuscript, deploy, run Abacus jobs or launch workers as part of
implementing this plan without the relevant user authorization.

## 1. Existing implementation and placement decision

Inspected the actual manifests and source, not just the older root guide:

- combinatoric-core supplies Composition and Ring, and already depends on
  num-bigint, num-rational, num-integer and num-traits.
- sym-poly-core depends on combinatoric-core and re-exports its foundational
  types. It supplies Field, PrimeField, dense exact linear algebra, generic
  sparse RREF/rank/kernel/quotient routines, and experimental packed rows.
- sym-poly-core::sparse_linear_algebra::sparse_rank currently computes full
  RREF, including elimination above pivots; useful as a reference, potentially
  wasteful for large rank-only calls.
- polytool::linalg already supplies streaming modular elimination, optional
  solution computation, input/NNZ/initial-Markowitz/dynamic-Markowitz ordering,
  and elimination statistics. Dynamic Markowitz currently collects the matrix
  and rescans active rows; its overhead must be measured, not assumed cheap.
- Polytool is also distributed by a standalone subtree projection. Do not add
  a path dependency from it into another workspace crate for this task.
- No generic chain-cancellation or Smith-normal-form implementation was found
  in the searched tracked and ignored source files.

Proposed additive layout (no new foundational crate required):

    combinatoric-core/src/
      sparse_matrix.rs            # generic sparse storage and operations
      integer_linear_algebra.rs   # exact Smith invariant factors / certificates
      chain_complex.rs            # graded complexes, unit cancellation, homology
    experiments/src/
      principal_stratum.rs        # Boris-specific cells, boundary and automaton
      bin/principal_stratum_homology.rs

Add narrow adapters in sym-poly-core as needed, with compatibility re-exports
for new foundational types. Keep its existing SparseVector and public APIs
working. The research driver may depend on both combinatoric-core and Polytool
to use the existing modular solver. No generic elimination algorithm is copied
into the driver. Do not migrate all existing algebra modules in this change.

If future consumers need a common Field trait at the foundational layer,
perform a separate compatible move with re-exports; it is not required for
Ring-based storage or the initial BigInt algorithms.

## 2. Shared sparse-matrix contract

Expose a coefficient-generic SparseMatrix<C>, with explicit row and column
counts even for 0-by-n, n-by-0 and all-zero matrices. Initial immutable format:
CSR (row offsets, sorted column indices, coefficients). A triplet/row builder
aggregates duplicate entries and removes exact zeros before finalization.

Required operations:

- Checked constructors and deterministic row iteration.
- Transpose and matrix-vector application; streaming adapters for existing
  sparse row solvers, preserving shape without dense conversion.
- Exact sparse composition-is-zero checking with an early witness, rather
  than constructing a potentially large full matrix product.
- Explicit, budget-checked conversion to dense matrices for small references
  and reduced Smith problems only.
- Dimension, index, allocation-budget and coefficient-conversion errors.

CSR is not the cancellation data structure. Provide a shared mutable sparse
form with access to both a row and the rows incident to a column. Start with
deterministic row maps and column incidence sets, storing each coefficient
once; measure their allocator overhead. Hide storage behind methods so a
compact edge arena can replace maps later if profiling justifies it. Maintain
bidirectional incidence and NNZ counters on every update and remove zero edges.

Use usize indexing first; packed u32 indices require checked conversions.
BigInt is the canonical integer coefficient type. Do not silently inherit
unchecked i64 arithmetic from the generic Ring implementation. A checked
small-integer fast path is a later measured optimization, not a second
mathematical algorithm.

## 3. Problem-specific state generation

Input: a nonempty composition omega of strictly positive parts and d with
d >= sum(omega) and the same parity. Validate explicitly: Composition::new
currently strips trailing zeros but does not reject internal zeros.

Generate descendants eta by merges of adjacent parts and insertion of 2,
with total weight at most d. Grade each cell by

    k = d - sum(eta) + length(eta).

Keep two independently implemented enumerators:

1. Simple closure/BFS for small reference cases, using a deduplicating set.
2. Production DFS through the deterministic membership automaton. Its state
   is the subset of positions in omega that could have been consumed. Retain
   all possible witnesses; a greedy single witness is incorrect. Prune a
   prefix only using a proved accepting-completion bound and remaining weight.

Intern each cell once, assign deterministic degree-local IDs, and keep the
lookup index separate from matrix coefficients. Do not build an exhaustive
poset/transitive-closure object: only cells and direct signed moves are needed.
Generate one degree at a time where useful; avoid storing all boundaries for
the initial modular rank pass.

For a length-r cell, the boundary has r-1 merge contributions and, when the
weight bound allows them, r+1 insertion contributions, each with sign (-1)^j.
Aggregate repeated targets before reduction. Check target membership and that
each term lowers k by one. Carry a count-only automaton DP as an independent
cell-count and Euler-characteristic check; it does not compute homology.

## 4. Shared chain-complex API and field baseline

A finite based chain complex records generator counts per integer degree and
D_k of shape c_(k-1)-by-c_k. Missing end maps are shaped zero maps. Validate
dimensions and D_(k-1) D_k = 0 exactly; supply a witness on failure. The
generic layer assumes no particular geometric or reduced-homology convention:
the model adapter is responsible for including any augmentation generators.

First compute field Betti numbers using existing rank engines:

    b_k = c_k - rank(D_k) - rank(D_(k+1)).

Generate boundary columns as rows of D_k transpose: the existing streaming
modular solver can consume them with RHS zero and compute_solution=false.
No explicit transposition or kernel basis is required for dimensions. Compare
streaming, NNZ order and existing Markowitz modes on the same fixed cases.
Use sym-poly-core's rational rank implementation only for small oracle cases
or sufficiently reduced matrices, never silently densify a large input.

Result records must distinguish exact F_p Betti numbers, exact rational Betti
numbers, exact integral groups, and incomplete computations. Matching several
primes is evidence, not a torsion-free or rational-rank certificate.

## 5. Shared integral unit cancellation

Own this algorithm in chain_complex, not the research binary. For a current
unit coefficient a = D[v,u] = +/-1, remove u and v and apply

    D'[y,x] = D[y,x] - D[y,u] * a^(-1) * D[v,x].

Maintain the entire graded differential coherently, including deletion of
incident rows/columns in neighboring degrees. Never independently simplify
adjacent matrices and reuse inconsistent bases. Cross-check the optimized
local update against a simple global-differential reference implementation.

Scheduling:

- Prefer zero-fill/low-fill unit pairs using current incidence counts.
- Use a deterministic priority queue with lazy invalidation/version stamps;
  recheck both coefficient and score before applying a queued pivot.
- Limit expensive score rebuilding. Track observed fill-in, not only estimates.
- Use one reducer thread per complex initially. Parallelize distinct input
  cases only within an explicit global memory budget.

Return the reduced complex, degree counts and reduction statistics. Optionally
stream a replayable certificate of stable pivot IDs and exact unit values;
the verifier reconstructs updates from the original complex without sharing
the pivot-selection implementation. Check certificate input digest, schema
and algorithm versions. Full chain maps/cycle representatives are not MVP
requirements and must not be stored implicitly at quadratic cost.

If all remaining differentials vanish, the surviving counts certify integral
free homology. If no unit pivot remains, report a residual complex, not a
torsion-free answer and not an assertion that no other reduction is possible.

## 6. Smith normal form: planned shared fallback

Implement a reusable BigInt Smith API after a small exact reference reducer;
execute it on a research case only when unit cancellation leaves a nonzero
residual. This is not a private one-off routine. Default output is nonzero
invariant factors a_1,...,a_r with a_i > 0 and a_i dividing a_(i+1), plus
shape and rank. Preserve the distinction between a diagonalization and Smith
normal form: a diagonal such as diag(2,3) must normalize to diag(1,6).

Use integer Euclidean/extended-gcd row and column transformations, with only
unimodular operations. Check that each finalized pivot divides the trailing
block; fix a violation rather than prematurely declaring Smith form. Begin
with a deterministic dense BigInt implementation for small residuals, behind
explicit entry and memory budgets. Sparse pre-reduction precedes conversion.
If the residual is too large, return a structured limit result; no unbounded
dense fallback. Large sparse Smith algorithms or external backends are a later
decision driven by residual sizes, not part of the first implementation.

### Important: invariants do not require full basis transformations

For a finite free integer chain complex, torsion(H_k) equals
torsion(coker D_(k+1)). Indeed, the exact sequence

    0 -> H_k -> coker D_(k+1) -> im D_k -> 0

splits abstractly because im D_k is free abelian. Thus the nonunit Smith
factors of D_(k+1), together with

    free_rank(H_k) = c_k - rank(D_k) - rank(D_(k+1)),

determine the integral group. Separate Smith invariants of the residual
boundary matrices suffice for abstract groups. Actual compatible cycle
representatives require extra basis information and are a separate option.

Offer an optional elementary-operation certificate; explicit dense U,V with
U A V = S are small-matrix/debug options with size limits. Verify unimodularity
through allowed operations, positivity/divisibility and reconstruction. Do
not confuse field row reduction or diagonal entries of an arbitrary
triangular matrix with integral invariant factors.

## 7. Tests and acceptance gates

Shared sparse storage:

- Duplicate/sign cancellation; ragged/out-of-range input rejection; empty
  dimensions; transpose twice; sparse/dense small products and applications.
- Randomized operation sequences preserve mutable/immutable equivalence,
  bidirectional incidence, exact zeros and NNZ counters.

Smith and chain algorithms:

- Rectangular, rank-deficient, negative, zero and large-BigInt matrices.
- Known diagonals, including diag(2,3) -> (1,6), and unimodular perturbations
  of prescribed Smith forms using fixed random seeds.
- For tiny matrices, compare products of the first j invariant factors with
  gcds of all j-by-j minors; this is a bounded independent oracle, not a
  production algorithm.
- Optional comparisons with an already available independent CAS; no new
  installation required for the core test suite.
- Complexes with free homology, Z/2 torsion, odd-prime torsion, mixed groups,
  adjacent nonzero homology and contractible pairs.
- Cancellation replay agrees with unreduced exact results on small cases;
  tests deliberately exercise fill-in, duplicate incidences and stale pivots.
- Different valid pivot orders yield identical groups, without requiring
  identical residual matrices. Reject tampered certificates.
- Verify the universal-coefficient dimension relation against F_p results:
  dim H_k(F_p) = free_rank(H_k) + p-torsion factor counts in H_k and H_(k-1).

Model calibration (always verify, do not hardcode homology as fact):

- omega=(1,1), d=4: cell counts by degree 4,3,2,1 are 1,4,3,1;
  integral H_3=Z, other degrees zero.
- omega=(3,1,1,3), d=18: 8,280 cells, Euler=-2; manuscript target H_7=Z^2.
- omega=(3,1,1,5), d=24: 76,384 cells, Euler=0; manuscript target
  H_7=H_8=Z. Euler zero must never short-circuit the homology computation.
- First extension: omega=(3,1,1,5), d=26: 212,900 cells; conjectural target
  H_8=H_9=Z. Stop and investigate any mismatch; do not change the algorithm
  merely to fit the conjecture.

The supervisor independently obtained the preceding cell counts by automaton
DP, not by homology computation. Further planning counts: d=28: 578,786;
d=30: 1,540,449; d=32: 4,025,371; d=34: 10,350,929.

## 8. Implementation checkpoints and benchmarks

1. Shared sparse storage, minimal shaped complex type and invariants tests.
2. Model generator/boundary plus existing modular rank baseline; cross-check
   enumerators, exact D^2 and small rational ranks before performance claims.
3. Shared integral unit reducer and independent certificate replay; reproduce
   both manuscript counterexamples or report the first discrepancy.
4. Shared bounded Smith fallback and integral homology assembler, with genuine
   torsion tests even if the first research examples need no Smith step.
5. Measure and extend to d=26. Only then consider d=28/30, compact storage,
   specialized F_2 bitsets or structural Morse matching.

At each checkpoint run focused tests, formatting/diff checks and relevant
consumer regressions, then make a focused commit. The public API is additive;
do not break standalone Polytool packaging. Before implementation, inspect
current ownership and manifests again. One worker owns generic algebra files;
another worker must not concurrently modify those files or shared manifests.
This plan does not start any workers.

Report separately: enumeration/indexing time, boundary time, validation time,
unit-reduction time, residual rank/Smith time, cell counts, initial/peak/final
NNZ, peak RSS, pivot update count, maximum coefficient bit length, and residual
dimensions. Include selected pivot mode, commit identity and exact input.
Do not compare automaton counting time with full complex generation time.

All runs use external Docker Cargo output. Initial measured executions get
bounded time/memory/cell/NNZ limits; failures preserve enough diagnostics to
resume investigation, never return partial output as a completed homology
answer. Keep checkpoints/certificates outside Dropbox with a small provenance
record in the project. No broad rank sweeps or uncontrolled parallel builds.

## 9. Deferred research and non-goals

No proof of the resonance conjecture is implied by these computations.
Structural finite-state Morse reductions are a later research route, requiring
an acyclicity/correctness proof before simultaneous cancellations. Do not
assume reductions at d remain valid at d+2: truncation changes boundary terms.
Reusable cell/transition data may be cached, but reductions need separate
validation or a proved filtration-compatible construction.

Defer full representatives, cup products/cohomology rings, GPU kernels,
distributed elimination, blanket algebra-module migration and unrestricted
sparse Smith normal form. Reassess only after the baseline identifies the
actual bottleneck. Original scripts and certificates from Boris would improve
comparison, but are not needed to build an independent implementation.

## References checked

- Sage chain-complex homology and Smith-based algorithms:
  https://doc.sagemath.org/html/en/reference/homology/sage/homology/chain_complex.html
- Sage's sparse preprocessing before Smith computation:
  https://doc.sagemath.org/html/en/reference/homology/sage/homology/matrix_utils.html
- Kozlov, Discrete Morse Theory for free chain complexes:
  https://arxiv.org/abs/cs/0504090
