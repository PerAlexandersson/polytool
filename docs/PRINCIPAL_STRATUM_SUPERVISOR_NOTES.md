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
