# Independent supervisor review notes

2026-09-10. Owner: host supervisor; implementation worker should return fixes
in source/tests, not edit this note. These are observations of in-progress
code, not findings against a final reviewed release. Recheck before closing.

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
