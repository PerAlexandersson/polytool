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
