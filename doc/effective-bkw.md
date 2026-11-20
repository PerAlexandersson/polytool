# Effective dominant-pair certificates

`bkw::effective::certify_effective_bkw` proves all-index splitting from a
rational generating function `N(t,z)/D(t,z)`, with `D(t,0)=1`. It returns a
replayable `EffectiveBkwReport`; `EffectiveBkwReport::verify` reconstructs
all arithmetic and compares the complete report. Recurrence inputs reuse
`recurrence::rational_generating_function`, including every supplied seed.
No finite sample, limit-set hypothesis or floating computation is a proof.

We factor `D=(1-lambda*z)(1-conj(lambda)*z)H`, choosing the continuous
upper-half-plane root `lambda`, and put `F=N/H`. The exact pair expansion is
`P_n=2 Re(A lambda^n)+rest_n`, where `A` is the pair amplitude. We use strict
remainder domination to fix signs at successive first crossings of the
continuous phase. The intermediate value theorem gives distinct zeros.
We close the count with the all-index degree bound and exterior signs,
then check every row below the cutoff exactly. First crossings give ordered
samples even when the phase is not monotone.

## Suppliers and supported geometry

- Order two: proper numerator, bounded pair interval with simple collisions
  at negative and positive eigenvalues. The remainder is zero. A real
  nonzero amplitude or a lift in the lower half-plane certifies the winding.
- Order three: one real root `r` and one conjugate pair on a negative ray.
  Exact inequalities check `|r|<|lambda|`, `Re A>0`, `Im A<0` and
  `4|A|^2>=C^2`. They give strict domination for `n>=1`, even when the
  modulus ratio tends to one at negative infinity. A rational positive
  double root and negative simple root of the rescaled characteristic,
  an implicit Taylor expansion, and a positive order of
  `(Im A/Re A)^2` certify the edge limits. The finite end must be a positive
  simple pair collision with phase limit `-pi/2`.
  Positive-infinity rays use the exact reflection
  `Q_n(u)=(-1)^n P_n(-u)`, recorded in `coordinate_transport`.
- Order four: negative characteristic discriminant and negative constant
  coefficient give a unique positive pair-product branch `b=|lambda|^2`.
  The cofactor roots are real and have opposite signs. Two exact Jury
  inequalities check a uniform modulus bound. Supplied `q,m,beta` bound
  the amplitude and the numerator weights on the whole interval.
- Orders through eight: a supplied polynomial pair factorization
  `D=(1-2a z+b z^2)H` is checked identically. Rational-polynomial `a,b` and
  coefficient bounds for `H` supply a compact template even when the
  automatic pair-product parametrization is unavailable.

For the compact supplier, `uniform_constants` has `modulus_ratio=q` in
`(0,1)`, `amplitude_lower_bound=m>0`, and `numerator_weights=beta`, including
`beta_0`. It certifies `|N_j|/|lambda|^j <= beta_j` and `|N/H|>=m` at the
reciprocal pair root. With `d=deg_z D`, its exact tail is

```text
T_n = sum_j beta_j sum_(k>=n+2-j) binomial(k+d-3,d-3) q^k.
```

Positive-series tails decrease. The first checked cutoff at most 4096 with
`T_n<m` works for every later index. `uniform_tail_cutoff` exposes the same
exact rational calculation for arbitrary orders through 32.

For a supplied factorization, `pair_factorization` contains ascending
`real_part=a(t)`, `norm_squared=b(t)` and `cofactor=H(t,z)` arrays. Its
uniform constants additionally contain `cofactor_weights=gamma`, one for
each `H_j`, `j>=1`. The checker proves
`|H_j|/|lambda|^j <= gamma_j` and `sum gamma_j/q^j<1`; the triangle inequality
then puts every cofactor eigenvalue strictly inside `q|lambda|`. This
Cauchy supplier is conservative and permits complex cofactor roots.

## Replay evidence

The certificate records the characteristic and branch equations, branch
selection, expansion formulas, domination supplier, cutoff and exact tail,
phase endpoints/winding, exterior signs, periodic zero-count law, and all
finite rows below the cutoff with exact degree and distinct-root counts.
`expansion_data` stores the rational pair parameters and cofactor coefficients
as structured bivariate numerator/denominator arrays, with the cubic real
remainder root and amplitude when present. `phase_limits` stores rational
multiples of pi for both eigenvalue and amplitude arguments. These fields
can be consumed directly without parsing the displayed explanation strings.
All certificate coordinates and finite coefficients belong to the recorded
chart family. Identity transport uses the input family; reflected transport
records `N(-u,-z)/D(-u,-z)` and preserves root multiplicities exactly.
`root_region` is stated in the input variable. If the pair interval has
nonpositive upper endpoint, all phase roots and left exterior roots are
negative. We extend this conclusion over the finite prefix by checking a
nonzero constant coefficient and coefficients of one sign after the exact
real-root count. Reflection transports a negative ray to a positive ray.
`degree_certificate` records the reduced top-coefficient GF, the periodic
degree cap and the residue classes with positive top coefficient.
Each branch sign includes its polynomial, exact resultant, rational sample,
branch-root enclosure, and boxes around all off-branch resultant zeros.

A zero of a branch condition must be a zero of its resultant. Exact Sturm
isolation separates these parameter roots from algebraic endpoints. A
root-free condition has a constant sign on the selected continuous branch.
We refine each isolated root inside its current interval. Overlap with
another root's enclosure does not identify the roots.
For each interior resultant root, rational interval arithmetic checks
opposite signs on the walls of a branch box and a strict condition sign
throughout that box. Uniqueness of the selected real or positive branch
keeps it in the box, excluding that critical root. A possible on-branch
zero of a strict condition yields `not_proved`. For weak inequalities,
`nonnegative_cell_samples` instead checks the selected branch in every
open resultant cell; continuity extends nonnegativity across critical
points. A negative cell rejects the proposed inequality.

Endpoint collision formulas are checked at exact algebraic roots. Cubic
infinity evidence records `s=1/t`, the scaled characteristic, the simple
root's implicit Taylor coefficients, and the orders used in the phase limit.
The degree proof checks `deg_t D_j<=j`, `deg_t N_j<=j` and the leading
coefficient GF after rescaling. Current top-coefficient templates are a
positive geometric sequence or `1/(1-z^2)`. The latter gives the parity
degree law needed for the 4-by-n tiling example.

The certificate is shaped for the RealRooted `EffectiveBKW` API:
`PairExpansion`, its `Dom` suppliers, phase crossings, `EdgeLimits`,
`ZeroCount`, and `Exterior`. Recorded declaration names are API leads;
Lean replay remains `not_attempted`. This checker is an exact arithmetic
implementation of the analytic sign-grid argument, not a Lean proof emitter.

## Examples and discovery

The checked Rust example [effective_bkw.rs](../examples/effective_bkw.rs)
constructs and replays the tiling certificate and prints its JSON.

The checked examples are Chebyshev `1/(1-2tz+z^2)` on `(-1,1)`,
the 4-by-n tiling GF
`(1-tz)/(1-(1+t)z-(2t+t^2)z^2+(t^2+t^3)z^3)` on its negative pair ray,
and the path-cube transfer GF

```text
(1-2tz-2tz^2+t^2(t-1)z^3) /
(1-(1+3t)z-2t(1-t)z^2-t(2-t-t^2)z^3-t^2(1-t)^2z^4).
```

For the last family, `q=4/5`, `m=13/100` and
`beta=(1,8/5,3/5,1/2)` give cutoff 39. The phase supplies `n-1`
interior zeros and the positive monic coefficient supplies one exterior zero.
The tiling supplier needs no asymptotic cutoff: for `n>=1` it supplies
`n-1` sign changes, one further zero for even `n`, and degree `<=n-1` for
odd `n`, hence `2 floor(n/2)` simple real zeros.

`check_effective_bkw_theorem` is the `effective_bkw` registry route.
`discover_recurrence_proofs` and the parent prover offer it automatically;
`RouteCertificate::EffectiveBkw` embeds the complete report.
Discovery records a successful route only when it strengthens earlier
conclusions. Existing proof-report verification replays it and ties its GF
to the exact recurrence and every seed. Consumers of recurrence discovery
therefore receive this route without a separate adapter or a sequence-name
lookup. Without supplied constants, the quartic supplier searches `q=k/16`
for `1<=k<=15`, `m=2^(-k)` for `0<=k<=16`, and each `beta_j` in
`{0} union {2^k : 0<=k<=16}`. It accepts the first certified bound for
each obligation. Exact sample checks only reject candidates; full-interval
sign certificates establish every accepted bound and are retained for replay.
These grids are independent of the family and numerator degree.
For the path-cube example this search finds `q=13/16`, `m=1/8` and
`beta=(1,2,1,1)`, giving cutoff 45 without supplied constants.
The certified negative ray makes the tiling and path-cube proofs eligible
in the existing parent-prover cascade. Downstream consumers can also replay
theorem reports whose region is the whole real line through their existing
theorem-report ingestion.
Older v1 reports may omit this added route; replay checks their own recorded
conclusions and does not inherit the new proof.

## Limits

The numerator must be nonzero and proper. Polynomial parts, general winding
across multiple half-planes, multiple pair intervals, finite exterior
amplitude-zero counts, triple collisions and
general algebraic higher-order factorization discovery are unsupported.
The degree templates and endpoint orientations are deliberately narrow.
Automatic inputs have order at most four; supplied polynomial factorizations
have order at most eight. Parameter degree is capped at eight, branch box
refinement at 128 bits/steps and cutoff search at 4096. Failure of any
sufficient condition is reported with its obligation. Simplicity is proved;
consecutive interlacing is not claimed.

Recurrence records with more prescribed seeds than their order can produce
an improper numerator and are then unsupported. The route does not silently
discard these seeds or advance its claimed starting index. Supporting a
shifted tail also requires separate prefix evidence and replay binding to
the original complete seed record.

The quartic proof and its two exact negative controls use the existing
`full-tests` tier because elimination is costly in debug builds. Run their
complete focused suite with
`cargo test --release -p polytool --features full-tests --lib bkw::effective::tests`.
