# Effective dominant-pair certificates

`bkw::effective::certify_effective_bkw` proves all-index splitting from a
rational generating function `N(t,z)/D(t,z)`, with `D(t,0)=1`. It returns a
replayable `EffectiveBkwReport`; `EffectiveBkwReport::verify` reconstructs
all arithmetic and compares the complete report. Recurrence inputs reuse
`recurrence::rational_generating_function`, including every supplied seed.
Zero prefixes and polynomial parts of recurrence GFs can be removed by exact
series shifting. `source_recurrence` retains the complete original record,
`skipped_prefix` records every removed row, and replay recomputes the shift.
The nonzero simple-root conclusion starts at the shifted `first_index`.
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
  The original degree-one/parity chart checks `|r|<|lambda|`, `Re A>0`, `Im A<0` and
  `4|A|^2>=C^2`. They give strict domination for `n>=1`, even when the
  modulus ratio tends to one at negative infinity. A rational positive
  double root and negative simple root of the rescaled characteristic,
  an implicit Taylor expansion, and a positive order of
  `(Im A/Re A)^2` certify the edge limits. The finite end must be a positive
  simple pair collision with phase limit `-pi/2`.
  Positive-infinity rays use the exact reflection
  `Q_n(u)=(-1)^n P_n(-u)`, recorded in `coordinate_transport`.
  Newton charts below instead fix the lift by `Re A>0` and calculate both
  endpoint arguments from exact weighted expansions.
- Order four: negative characteristic discriminant and negative constant
  coefficient give a unique positive pair-product branch `b=|lambda|^2`.
  The cofactor roots are real and have opposite signs. Two exact Jury
  inequalities check a uniform modulus bound. Supplied `q,m,beta` bound
  the amplitude and the numerator weights on the whole interval.
- Orders through eight: a supplied polynomial pair factorization
  `D=(1-2a z+b z^2)H` is checked identically. Rational-polynomial `a,b` and
  coefficient bounds for `H` supply a compact template even when the
  automatic pair-product parametrization is unavailable.
- Automatic orders five through eight reuse exact characteristic factorization
  over `Q[t]`. Every quadratic candidate is tried against the compact/Cauchy
  supplier, including when the cofactor contains another nonreal pair.
  If no polynomial pair is certified,
  Newton identities construct the pair-product resultant
  `prod_(i<j)(b-lambda_i*lambda_j)`. A quadratic subresultant recovers the
  pair and cofactor. Nonvanishing resultant discriminant and constant term
  fix the continuous largest positive branch; full-interval pair and Cauchy
  inequalities must still hold. The stored real-part and cofactor
  denominators are checked strictly nonzero before quotient arithmetic can
  cancel any of their factors. At each endpoint, the cofactor obtained by
  deflating the double root must satisfy the Cauchy bounds there as well.
  This binds the collision to the dominant pair rather than to two cofactor
  roots. Branch collisions produce `not_proved`.

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
When constants are omitted for a polynomial higher-order pair, each cofactor
weight is searched in `{0} union {2^k : -16<=k<=16}`, then the ratio grid
is checked against the strict Cauchy budget. The accepted constants are stored
in `automatic_uniform_constants` for replay and inspection.
An identically zero cofactor coefficient receives `gamma=0` with an exact
zero-identity certificate; the Cauchy coefficient inequality is weak.

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
The degree proof retains the positive-geometric and parity certificates and
otherwise uses the top Newton edge of the denominator. For slope `g/p` in
lowest terms, substitute `t=u^(-p), z=z*u^g` and multiply the numerator by
`u^h` to remove its negative powers. Exact rational coefficient layers in
`u` identify the first eventual nonzero layer `e_r` in each index residue.
Earlier layers are zero in that residue or have exactly checked finite
polynomial support, giving
`deg P_n=(g*n+h-e_r)/p` once the selected layer is positive. Positivity
comes from a nonnegative renewal expansion, a checked binomial power of the
Newton face, or nonnegative layer numerators
after cancelling Newton-face powers. Binomial faces additionally use an
exact polynomial formula for each coefficient residue, with a replayable
integer-index sign certificate. A finite support seed in each residue
modulo a positive face step proves positivity for every later index.

`degree_certificate.newton` stores the weights, layers, support seeds, law
start and exact earlier degrees. For Newton certificates,
`leading_numerator/leading_denominator` is the GF of the actual leading
coefficients, including its exactly computed finite correction. Legacy
certificates retain their original diagonal top-coefficient GF.
`zero_count.divisor` defaults to one; the general count is
`(slope*n+offsets[n % period])/divisor`.

Quadratic rays and cubic single-remainder rays support fractional Newton
scaling at infinity. An implicit Taylor expansion at a rational simple
scaled root supplies exact phase limits. For a binomial edge `u^d-c`, put
`k=d*g/p`. When `gcd(d,k)=1`, choose nonnegative `beta<d` and positive
`alpha` with `d*alpha-k*beta=1`. The rational dilations
`r=c^alpha*u*s^(-g)`, `t=-c^beta*s^(-p)` normalize the limit to `u^d-1`.
This algebraic normalizer works for any coprime `d,k`; the ray supplier
currently calls it only for cubics. For cubics with `3|k`, rational
normalization is possible exactly when `c` is a rational cube: any rational
scales would imply `c=(eigen_scale/parameter_scale^(k/3))^3`.
The existing linear-factor search already handles this rational-root case;
otherwise an algebraic coefficient field is required.

Rational squared cosine/tangent limits are recognized by exact cyclotomic
divisibility, without value tables. For `v=cos(theta)^2`, the unit-circle
root `exp(2 i theta)` satisfies `x^2-2*(2v-1)*x+1`. A rational multiple of
pi is equivalent to a cyclotomic factor. Its degree is at most two, and
`phi(m)^2>=m/2` bounds its order by eight. Thus the finite order search is
complete for every rational `v`, rather than an angle-denominator guard.
Tangent squares use `cos(theta)^2=1/(1+tan(theta)^2)`; signs fix the lift.
Algebraic squared trigonometric inputs require a different coefficient type.
The crossing count
uses strictly interior integer phase levels. On a common residue period,
their count and both exterior signs must equal the exact degree, and their
affine increments must agree. Finite rows cover all earlier indices.

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
A118424 now stops at its direct Newton certificate. Its older affine proof
is retained explicitly in the affine-coordinate regression: transport by
`P_n(u+1)`, replay the fixed-delay theorem on `(-infinity,0)`, and transport
the roots back to `(-infinity,1)`. This auxiliary evidence is available even
though the default discovery report no longer runs that later stage.
For a recorded open basis-map parent, replay can restrict the parent routes
to none to retain the original conditional root-bound requirement. It still
recomputes the identity, transfer and hypotheses and never promotes an open
parent to a proved one.

## Limits

Direct GF input requires a nonzero proper numerator. General winding
across multiple half-planes, multiple pair intervals, finite exterior
amplitude-zero counts and triple collisions are unsupported.
Higher-order automatic suppliers require compact geometry; higher-order ray
edges with a degenerating gap remain unsupported.
Compact charts with eigenvalue arguments pi and zero require degree slope
one; a sublinear Newton cap is rejected before elimination. Newton layers whose
positivity cannot be certified are rejected. Newton search has 64 layers,
1024 support coefficients and phase period at most 1024; cubic Taylor limits
use sixteen terms. General algebraic scaled edge roots and algebraic squared
trigonometric inputs are unsupported. All inputs have order at most eight.
Parameter degree is capped at eight, branch box
refinement at 128 bits/steps and cutoff search at 4096. Failure of any
sufficient condition is reported with its obligation. Simplicity is proved;
consecutive interlacing is not claimed.
Before exact GF construction, the registry also caps seed `j` at degree
`8*(j+1)`. This is a deliberately generous resource guard tied to the
coefficient-degree cap, not a necessary real-rootedness condition or the
Newton degree law. An oversized seed gives `not_applicable`; accepted inputs
still require the complete exact degree proof.

Removed prefix rows must be zero or simple real-rooted. Their removal never
asserts that a zero row is nonzero. Nonzero prefix rows with repeated roots
remain unsupported by this simple-root route.

The quartic proof and its two exact negative controls use the existing
`full-tests` tier because elimination is costly in debug builds. Run their
complete focused suite with
`cargo test --release -p polytool --features full-tests --lib bkw::effective::tests`.
