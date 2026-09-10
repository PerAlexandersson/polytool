# Independent supervisor review notes

2026-09-10. Owner: host supervisor; implementation worker should return fixes
in source/tests, not edit this note. These are observations of in-progress
code, not findings against a final reviewed release. Recheck before closing.

Review continuation note (10:32 UTC): the first Claude run DID inspect source
and run independent tests. Its empty final-report file was not evidence of a
read-permission failure. Its session is
`328b0217-1a44-4d3b-ae9b-51909d779e73`; the transcript is under
`/home/dev/.claude/projects/-workspace-rust/` in Docker. It spent considerable
time searching over omega to infer undocumented command parameters, despite
those parameters already being in the plan. Prefer resuming its context and
requesting findings over repeating that search. Host observation only: the
host did not kill, restart, or signal the reviewer or its Python child.

Calibration inputs (also specified in the plan): d=18 uses omega=(3,1,1,3);
d=24 and d=26 use omega=(3,1,1,5). There is no need to search over omega to
recover these parameters. The observed d=24 debug timings were enumeration
0.75 s, assembly plus validation 122 s, and field ranks 5.9 s; an incomplete
ten-minute d=26 validation should prompt fixing the repeated column scans,
not treating generation or modular elimination as the measured bottleneck.

1. Initial compose_is_zero used left.row(middle) rather than its column.
   That orientation is now corrected in the working tree, but column_entries
   scans the entire matrix for each requested column. Use rowwise sparse
   accumulation or transpose/index once before large chain validation. Add a
   rectangular, nonzero cancellation product test, not only A times zero.
2. Smith certificate replay currently permits AddRowMultiple with target ==
   source, and the analogous column operation. Those are not generally
   unimodular. Reject such forged operations (or explicitly verify the rare
   unit scaling cases). Require distinct coordinates for two-row/column
   Bezout transformations as well.
3. The replay Bezout check x*a+y*b=g alone does not imply that truncated integer
   divisions a/g and b/g define a unimodular operation. Require exact
   divisibility and determinant +/-1. Concrete forged certificate: a=b=1,
   x=y=1,g=2 passes the identity but creates a zero second row/column through
   truncation. Include regression tests for this and same-index operations.
4. Initial Smith replay sets the dense budget to the computed entry count,
   giving the caller no effective dense allocation cap. Apply explicit replay
   limits, including coefficient sizes and input preflight. Similarly check
   limits on an already diagonal huge-coefficient input, before any operation.
5. Sparse zero/builder finish currently allocate rows+1 without a fallible
   checked shape/budget API. A zero-column matrix can also bypass an entry-only
   dense budget while allocating arbitrarily many row vectors. Decide and
   document realistic checked shape/allocation limits.
6. The initial unit reducer computes total NNZ and scans every coefficient for
   maximum bit length after each pivot. Profile this global work; use maintained
   NNZ totals and inspect changed entries for coefficient guards. Large
   calibration cases must not be abandoned without checking these avoidable
   repeated scans. Keep exact validation, rather than skipping D^2 to mask cost.

The host did not run concurrent Cargo builds or edit the worker's source.
These notes supplement, and do not replace, the required Claude review.

Further observations at the first model implementation (~12 minutes):

7. Unit cancellation currently rebuilds candidate_queue over the entire
   complex after every pivot. That is not the planned local/lazy update and
   will multiply work by tens of thousands of pivots. Profile and replace;
   merely having a BinaryHeap does not make this an incremental scheduler.
8. The model currently uses BFS production and DFS closure as its second
   enumerator. Different traversal order of the same moves is not the planned
   independent membership-automaton validation. Add the actual subset automaton
   and test accepted words independently; do not describe closure DFS as that
   oracle. Cell-count-only DP should independently confirm the large counts.
9. The modular path currently materializes/transposes every boundary and
   collects all rows before invoking the streaming solver. This is not the
   planned boundary-on-demand baseline. Implement the stream or document and
   justify measured limits, and enforce cell/NNZ limits during generation,
   not only after building the entire complex.

Further observations (~24 minutes):

10. Validate the modulus before converting coefficients: bigint_mod_u64 uses
    mod_floor(prime) before the existing solver checks primality, so prime=0
    can panic on a nonzero boundary. Reject invalid/composite moduli with a
    clear error, including zero, one and four regression tests.
11. Initial reducer/replay NNZ and coefficient budgets must be checked even
    when there are no unit pivots. Check limits before cloning/allocating where
    possible. max_cells=0 must reject the initial BFS cell as well.
12. Canonical digest encoding should be unambiguous: prefix sections, counts
    and BigInt byte lengths rather than separating arbitrary integer bytes
    with a zero byte (which can occur inside their representation).
13. Grading arithmetic degree +/- 1 (including in errors/certificates) must
    not panic or wrap at i32 extrema. Decide a checked supported range at
    construction and reject hostile out-of-range certificate degrees.

Recheck after 19632cb (first review-fix checkpoint):

14. The new replay current-entry checks index matrix[first][first]. RowBezout
    validates first as a row, not as a column, so a tall rectangular input can
    panic there; ColumnBezout has the dual wide-matrix problem. Enforce every
    accessed index, or verify the genuinely general unimodular operation
    without requiring pivot-specific entries. Same-index additive and Bezout
    operations from point 2 still need rejection. Add adversarial rectangular
    and zero-shape tests; valid generated operation logs must still replay.
15. --record-certificate currently prints only an input digest and pivot count,
    then discards the actual pivot sequence. That is not a production replay
    artifact. Provide bounded emission/reloading and an independent replay
    command (or explicit in-run replay plus a documented artifact API). Bind
    claimed residual/Smith factors to replay, not just record an unused log.
16. universal_coefficient_dimension now exposes another public prime=0 panic
    through mod_floor. Validate its input or accept an already validated field
    type; test the previous-degree torsion term, not only current torsion.

After bc70287, the new guards look correct, but its purported rectangular
RowBezout regression uses a 1-by-2 matrix with first=1: that already failed
the old row-index check. The actual regression should be a 2-by-1 matrix,
first=1, second=0, a=b=g=2, x=1,y=0; the old code accepts both row indices
then panics reading column 1. Include the dual ColumnBezout 1-by-2 case and
test both column and row self-add/Bezout rejection.

Completion-pass observations (11:09 UTC; recheck after edits settle):

17. The new UCT prime check trial-divides up to sqrt(n). A valid large u64
    prime can require billions of divisions per homology degree. Use bounded
    fast validation or a documented validated-prime API; do not add this cost
    after the existing modular solver already handled the same input quickly.
    The existing Polytool linalg validator uses modular exponentiation, not
    trial division. Preserve standalone Polytool packaging if reusing code.
    The experiment has a second copy of the same trial-division helper. A
    host smoke check on the release binary, omega=(1,1), d=4,
    prime=18446744073709551557, hit a deliberately short two-second guard
    before field output. This is only nine cells; the check should be cheap.
    The host's normal small run passed integral/rational/field agreement and
    replay. Zero cell/NNZ budgets, moduli 0 and 4, and d=2147483648 all returned
    clean errors. No Cargo build or source edit was run by the host.
18. The new digest has section counts and value byte lengths, but still needs
    each differential's NNZ/entry count before its variable entry sequence.
    A one-byte D marker alone is not a delimiter for arbitrary binary indices.

19. Blocking check on the new Miller--Rabin draft: bases 2,3,5,7,11,13,17
    do NOT suffice for all u64 inputs. Exact independently checked witness:
    341550071728321 = 10670053 * 32010157. It passes all seven chosen strong
    probable-prime tests (also base 19) and is not divisible by the trial
    primes through 37. Base 23 rejects it. The existing Polytool implementation
    tests all twelve prime bases through 37, not only the first seven. Restore
    a genuinely complete u64 witness set and add this composite regression,
    along with a large valid prime, before calling the public check exact.

20. Missing-zero materialization still bypasses resource budgets before the
    bounded Smith call. With counts {0:10^12,1:0}, no stored differentials,
    construction is cheap, but integral_homology calls differential_or_zero(1)
    and tries to allocate 10^12 row offsets before Smith sees its limits.
    validate() likewise materializes a missing adjacent map unnecessarily.
    Treat absent maps as rank/factors zero and skip products with absent maps,
    borrowing stored matrices instead of allocating shaped zeros/clones.
    Add metadata-only large-free-complex tests that allocate no generators.
    Do not execute an allocation-bomb reproduction before fixing this.
    Also check differential_or_zero(i32::MIN): constructor grade validation
    does not protect this independent public method argument from degree-1
    overflow; a checked/bounded retrieval path should cover it.

## Independent release-binary calibration checks

Host verification on 2026-09-10 during the completion pass. These invoked the
worker-built release executable directly; the host did not run Cargo or edit
source. Elapsed times are Bash `time`, exclude compilation, and include field
ranks, integral computation, UCT comparison and cancellation-certificate replay.
The existing release binary passed all three and replay matched each residual.

| omega | d | cells | nonzero integral groups | elapsed seconds |
|---|---:|---:|---|---:|
| (3,1,1,3) | 18 | 8,280 | H7 = Z^2 | 0.131 |
| (3,1,1,5) | 24 | 76,384 | H7 = Z, H8 = Z | 2.452 |
| (3,1,1,5) | 26 | 212,900 | H8 = Z, H9 = Z | 10.787 |

All groups listed were torsion-free; all other returned degrees were zero.
d26 field-only independently completed in 5.981 s. Its first construction
attempt with max-nnz=2000000 returned an orderly upper-bound-budget error;
max-nnz=3000000 allowed completion with canonical NNZ=1706568. This limit
counts a conservative contribution bound, not final aggregated NNZ.

Each run used `docker exec docker-setup-app-1 bash -c` with Bash `time`,
`timeout 60s nice -n 10`, and executable
`/cargo-target/ai-projects/release/principal_stratum_homology`.
Exact arguments for the integral-and-replay runs were:

```text
--omega 3,1,1,3 --d 18 --max-cells 10000 --max-nnz 250000 --max-reduction-nnz 1000000 --max-pivots 10000 --record-certificate
--omega 3,1,1,5 --d 24 --max-cells 80000 --max-nnz 1000000 --max-reduction-nnz 3000000 --max-pivots 100000 --record-certificate
--omega 3,1,1,5 --d 26 --max-cells 220000 --max-nnz 3000000 --max-reduction-nnz 4000000 --max-pivots 110000 --record-certificate
```

d24 stage times (ms): enumeration 18, assembly/validation 625, field 665,
unit reduction 737, Smith below 1. NNZ initial/peak/final=561173/561173/20;
38187 pivots; residual has five generators each in degrees 7 and 8.
d26 stage times (ms): enumeration 52, assembly/validation 2193, field 3577,
unit reduction 3307, Smith 1. NNZ initial/peak/final=1706568/1706568/90;
106440 pivots; residual has ten generators each in degrees 8 and 9.
Smith reduction was therefore exercised on nonzero residuals in these two
research calibrations, not only on synthetic torsion fixtures.

A second d26 integer/replay run sampled its own process VmHWM from
`/proc/<pid>/status` every 50 ms under the same 60-second outer guard.
Observed peak RSS was 482912 KiB (about 472 MiB); exit status was zero and
the groups/replay agreed again. Reduction reported maximum coefficient size
13 bits, 3498762 stale candidates and 134599 score revalidations. This was
an ephemeral read-only monitoring shell, not a new source or artifact file.
