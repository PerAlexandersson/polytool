# Implementation assignment: sparse integral homology

The user approved docs/PRINCIPAL_STRATUM_HOMOLOGY_PLAN.md for implementation.
Use the assigned Terra model for implementation, and Claude for a subsequent
independent code/mathematical review. Remain on your assigned account/model.
Read /workspace/AGENTS.md, /workspace/rust/AGENTS.md, HANDOFF.md and the complete
plan before acting. Inspect codex mcp list in Docker. The root Rust CLAUDE.md
does not exist; AGENTS.md is the current guide.

You own implementation and integration of all five planned checkpoints:
shared sparse storage, model/field baseline, shared integral cancellation,
shared Smith fallback with torsion tests, and calibrated bounded benchmarks.
Own only task files in combinatoric-core, the needed narrow sym-poly-core
adapters/re-exports, principal-stratum experiment module/binary/tests, their
required manifest/lib registrations, task docs and the opening HANDOFF entry.
Preserve unrelated work, particularly polytool/scripts/__pycache__, nested
Ehrcalc, and all other experiments. No other worker owns these task files.
Do not broadly format or refactor existing modules. Commit verified logical
increments regularly and stage exact files only. No pushes or publications.

The user especially requires reusable sparse matrices and Smith normal form
in the shared library. Do not replace the plan with a private script. Reuse
the existing modular solver; preserve standalone Polytool dependencies. Exact
BigInt algorithms first; no unchecked narrowing or floating-point proofs.

The plan is a specification, not an excuse to stop after a prototype. Complete
the feasible implementation and tests, and honestly document any resource
limits or justified scope changes. Do not claim a calibration or performance
result without running it. No artificial formal goal is requested.

Use Docker Cargo cache outside Dropbox. Serialize builds; run focused checks
first. For a new potentially runaway execution use reduced priority and a
60-second first-run limit, then enlarge only after observing behavior. Monitor
memory and NNZ. Existing compiler jobs may have longer bounded timeouts.
Attempt d=18 and d=24 calibrations with resource guards, then d=26 if safe.
Stop and investigate discrepancies instead of hardcoding expected homology.

## Independent Claude review (required before final handoff)

When implementation and focused tests are ready, checkpoint and freeze files.
Inspect the installed Claude CLI and run a read-only reviewer on the actual
diff and tests, preferably its high-capability model if available. The user
explicitly asks that Claude receive ample time: start with a monitored
30-minute timeout, not a 60-second timeout. Extend if it is making meaningful
progress; do not kill a healthy reviewer merely because it is slow.

You own the nested Claude process, its monitoring/retry/cleanup and logs. Keep
large logs outside Dropbox. Record model, reviewed commit, result, findings
and resolutions in a small task review note. Do not report review as completed
merely because it was launched. Do not read or print credentials. If Claude
authentication/usage blocks the review, report the precise non-secret blocker
and continue other verification without substituting an unrequested model.

Reviewer focus: sparse shape/zero edge cases, termination and divisibility in
Smith reduction, genuinely unimodular operations, integer overflow, graded
basis consistency of cancellations, fill-in limits, certificate independence
and rejection of malformed certificates, invariant-factor-to-homology logic,
field/integer evidence labels, resource guards, and model boundary signs.
Implement fixes as Terra, rerun tests, and seek a focused Claude re-review for
material correctness fixes. Reviewer must not edit task files.

Finish with commits, exact test/benchmark commands and outputs, review status,
remaining limits, and current ownership in HANDOFF.md. The host supervisor
will inspect results independently; it must not restart or kill your Claude
child. No remote Abacus jobs, new worker spawns, emails, deployments or pushes.
