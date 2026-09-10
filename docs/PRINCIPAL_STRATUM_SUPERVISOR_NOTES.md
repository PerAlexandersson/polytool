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
18. The new digest has section counts and value byte lengths, but still needs
    each differential's NNZ/entry count before its variable entry sequence.
    A one-byte D marker alone is not a delimiter for arbitrary binary indices.
