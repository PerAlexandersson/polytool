# Principal-stratum homology review record

## Frozen implementation under review

- Initial reviewed commit: `2303a1c` (plus the two preceding shared-library checkpoints).
- Reviewer: Claude Code 2.1.267, requested model `opus`, read-only plan mode.
- Scope: `combinatoric-core` sparse storage, chain cancellation and Smith
  reduction; narrow `sym-poly-core` adapters; and the principal-stratum model
  and driver in `experiments`.

The reviewer is instructed not to edit files.  Its requested focus is sparse
zero/shape handling, Smith termination/divisibility and unimodularity,
coefficient growth, cancellation basis consistency and certificates, abstract
homology assembly, model signs, and calibrated resource bounds.  The final
result and any implemented resolutions are appended here.

## Initial review (completed)

Claude's independent review found that Smith certificate replay accepted a
malformed Bézout step with a non-dividing claimed gcd, hence a non-unimodular
two-by-two transform. It also identified the quadratic CSR composition scan,
the driver hardcoding cancellation-certificate recording off, and lack of a
universal-coefficient consistency check. The review's mathematical assessment
was positive for the model signs, cancellation formula, Smith pivot reduction,
and torsion-degree convention.

The follow-up patch adds divisibility, positivity, and current-entry checks to
Bézout replay (with a malicious certificate regression), uses CSR transposes
on both inputs for composition checking, exposes `--record-certificate`, and
checks field dimensions against the integral groups through universal
coefficients. A focused Claude re-review is required for this patch.

## Focused re-review

Claude confirmed the transpose-based composition check and the
universal-coefficient/certificate-driver changes, then found three further
certificate-replay edge cases: same-index Bézout operations, self-addition
operations, and rectangular transverse indexing. The final guard patch rejects
all three and adds regression tests. A final targeted confirmation has been
requested on that guard-only patch.

## Final targeted confirmation

Claude found the final guards sound and valid recorder-produced certificates
still replay. It noted only that the initial rectangular rejection test did not
exercise either new transverse guard. That test coverage was corrected with a
2-by-1 row-Bézout and 1-by-2 column-Bézout rejection case; the focused Smith
and principal-stratum suites passed afterward. No repository file was edited
by Claude.

## Completion review (2026-09-10)

- Reviewer: Claude Code 2.1.267, requested `opus`, read-only plan mode.
- Reviewed checkpoint: `63e9861`; log retained outside Dropbox at
  `/tmp/principal-stratum-claude-63e9861.log`.
- Result: no blocking correctness finding. Claude independently ran the fixed
  d=18/d=24/d=26 inputs and confirmed the count, Euler, field/UCT and integral
  outputs, plus cancellation replay.

The review found and the follow-up `5bd050b` resolves these material issues:
bounded deterministic Miller--Rabin is now shared in combinatoric-core;
cancellation and Smith replay have explicit shape limits; default Smith replay
does not infer unlimited limits; integral assembly verifies recorded Smith
operations against claimed factors; the digest frames each differential with
its NNZ count; and the membership automaton is cross-checked against move
closure on the nontrivial `(3,1,1,3)` prefixes. Focused suites passed after
the fixes. The remaining implementation limits are the documented explicit
budgets; no larger calibration was run.

## Focused re-review resolution

The focused Opus re-review of `5bd050b` found that its first seven-base
Miller--Rabin set was not deterministic for all `u64`, exhibiting
`341550071728321`. Commit `d29227c` uses the full witness set through 37 and
adds both that composite and a large-prime regression; the focused core suite
then passed (240 tests). This is the final resolved review finding.

Claude's final micro-review of `d29227c` returned no blocker: it independently
confirmed the complete witness bound, the pseudoprime rejection, and the
M61 positive control. Log: `/tmp/principal-stratum-claude-primality-d29227c.log`.

## Missing-map resource review

Claude Opus reviewed `a56ac68` read-only and found no blocker: absent maps are
correctly treated as zero without allocation, checked retrieval handles
`i32::MIN`, and Smith/homology logic is unchanged. Its suggested huge-metadata
validation regression was added in `4d22c7e`. Log:
`/tmp/principal-stratum-claude-missing-map-a56ac68.log`.

## Stored bounded-retrieval review

Claude Code 2.1.267, requested `opus` in read-only plan mode, reviewed
`b134824` without editing files. It found no correctness blocker: the bounded
retrieval API checks `rows + 1` against the CSR shape-slot limit and NNZ
against its limit before cloning a stored matrix, returns the established
structured `ReductionLimit` errors, and leaves the trusted infallible helper
unchanged. The review independently ran the focused chain-complex suite (7
passed). Its minor notes were documentation/coverage polish only, so no source
change was warranted for this deliberately narrow completion pass. Log:
`/tmp/principal-stratum-claude-stored-retrieval-b134824-retry.log`.
