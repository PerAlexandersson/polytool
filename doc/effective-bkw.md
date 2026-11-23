# Effective dominant-pair certificates

`bkw::effective::certify_effective_bkw` proves all-index splitting from a
rational generating function `N(t,z)/D(t,z)`, with `D(t,0)=1`. It returns a
replayable `EffectiveBkwReport`; `EffectiveBkwReport::verify` checks stored
spectral, sign-grid and finite-row data, then regenerates the remaining
phase, domination and degree closure and compares the complete report. Recurrence inputs reuse
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
then check every row below the cutoff exactly. Rows of degree zero are
nonzero constants by the exact positive leading-coefficient certificate and
need no alternating samples. First crossings give ordered samples even when
the phase is not monotone.


## Pointwise remainder supplier

For a cubic, the strict inequality `4|A|^2>C^2` gives domination for every
`n>=0`; the existing weak inequality with a strict modulus gap remains a
fallback for `n>=1`. No remainder gives domination from zero automatically.

For a quartic with two real cofactor roots `x,y`, put
`h(mu)=mu^2+s*mu+c`, `q(mu)=mu^2-2a*mu+b`, and reduce
`Nbar(mu)/q(mu)` modulo `h` to `U*mu+V`. The real residues satisfy

```text
Cx+Cy = U,
(Cx-Cy)^2 = (-s*U+2*V)^2/(s^2-4c).
```

Thus two symmetric rational inequalities certify
`4|A|^2>(Cx+Cy)^2` and `4|A|^2>(Cx-Cy)^2`. Since
`|Cx|+|Cy|=max(|Cx+Cy|,|Cx-Cy|)`, they imply strict domination from zero.
Every denominator and the real-root and modulus conditions are checked on
the whole branch. In particular `b-c>0` is required in addition to
`b+c>0` and `(b+c)^2>b*s^2`; the last two alone do not exclude two real
roots outside the pair modulus on the same side.

We also try cutoffs `N=1,...,16`, replacing residues by `Cx*x^N,Cy*y^N`
and the pair amplitude by `A*lambda^N`. Exact reduction of
`mu^N*(U*mu+V)` modulo `h` supplies the same two inequalities. Strict
individual modulus gaps make the normalized absolute sum decrease for all
later indices. Samples only reject candidates; acceptance requires both
full-interval resultant/sign certificates. The compact uniform supplier
remains a fallback when these sufficient inequalities do not certify.

A nonreal quadratic cofactor contributes `2 Re(B*mu^n)`. Put
`Y_h=c-s^2/4>0`, and reduce `Nbar mod h=u*mu+v` and `q mod h=e*mu+f`.
Then

```text
4|B|^2 = (c*u^2+v^2-s*u*v)/(Y_h*(f^2-s*e*f+c*e^2)).
```

An exact inequality `|A|^2*b^N>|B|^2*c^N`, together with `b>c`, certifies
all later rows. The norm and residue denominators are nonzero because the
pairs have distinct moduli. On a quartic ray the characteristic discriminant
can be positive. The largest positive pair product is still continuous:
it is the maximum of the continuous conjugate-pair norms and positive
real-root products. A positive mixed product lies below the larger pair
norm. Collisions between smaller products therefore need not exclude this
chart. The characteristic roots remain simple, and the selected pair and
all its rational denominators are checked strictly.

A positive real part of `A` can also fix the bounded-interval lift when
`Im F` changes sign. Both endpoint arguments then enter the same exact
residue-class phase closure. This closure requires an exact degree and a
positive actual leading coefficient. The older parity upper-bound shortcut
keeps its original phase argument; its vanishing top coefficient cannot
supply generic exterior signs.

`domination_cutoff` records this independent start. Older certificates
that omit it are rejected explicitly and must be regenerated, including
enclosing theorem reports with moved seeds or rule normalization.
`cutoff` is the joint start after degree/phase closure,
and `finite_rows` contains exactly the earlier local rows. They can differ
when the eventual degree law starts later than domination. In the classical
charts, degree-zero rows allow both starts to be zero.

## Suppliers and supported geometry

### Replay contract

`polytool.effective-bkw.v1` identifies the report layout and exact replay
protocol. A proved report retains its complete request, supplier evidence,
finite prefix and optional original recurrence. Verification reconstructs
and compares every field; editing a witness or its source is rejected.
Certificates without `domination_cutoff` are rejected. Current suppliers
record `Some(N)`, including the uniform fallback. Stored legacy fixtures are
negative migration tests; current pointwise fixtures remain replay targets.

The schema remains v1 because its wire layout is unchanged; accepted proof
certificates now require the previously optional domination bound.
The schema is not an implementation commit: exact replay can reject stale
evidence when an arithmetic algorithm or witness label changes. Retain the
producing source version with archived reports. For migration, generate a
fresh report from the archived request and verify it; never bypass comparison
or edit the old witnesses. In particular, newer supplied-pair ray reports
whose pair-discriminant check was labelled as a full characteristic
discriminant must be regenerated to obtain the corrected label. The stored
pointwise fixtures are maintained replay compatibility targets.
An unsupported `not_proved` result records the current search outcome;
its diagnostic and completeness are allowed to change as suppliers improve.
The verifier still rejects any recorded result that differs from its current
reconstruction. Prefix-shift failures produce source-bound `not_proved`
reports with exit 0, while malformed input remains an error with exit 2.
Order-two no-remainder requests reject supplied `uniform_constants` rather
than ignoring them; legacy reports cannot bypass this rejection.

- Order two: proper numerator, bounded pair interval with simple collisions
  at negative and positive eigenvalues. The remainder is zero. A real
  nonzero amplitude or a lift in the lower half-plane certifies the winding.
  Its Newton ray chart also has zero remainder. Automatic endpoint selection
  includes `t=0` when the characteristic is nonreal on the left component;
  all interval, amplitude and phase checks still apply.
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
  Newton charts fix the lift by `Re A>0` or `Im A<0` and calculate endpoint
  arguments from exact weighted expansions. A generic positive triple point
  is also supported: the endpoint characteristic is `(mu-r0)^3`, `p_t<0`,
  and `Nbar(r0)>0`. Its Puiseux roots give `theta->0` and `arg A->-2pi/3`
  in the certified lower-half-plane lift. No uniform modulus gap is needed
  at this edge or at a cubic equimodular Newton face.
- Order four: negative characteristic discriminant and negative constant
  coefficient give a unique positive pair-product branch `b=|lambda|^2`.
  The cofactor roots are real and have opposite signs. By default we first
  try exact pointwise remainder domination, described below. Supplied
  `q,m,beta` retain the compact uniform supplier. Quartic Newton rays also
  allow same-sign real cofactor roots or a second conjugate pair. The largest
  positive pair product,
  its strict denominator sign, nonreal-pair test and strict cofactor gaps
  identify the dominant conjugates. A simple rational limiting pair-product
  root gives exact fractional angles at infinity; a finite positive dominant
  double collision supplies the second phase limit. Both limits close the
  degree separately in every residue class.
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

For a positive constant numerator, a bounded higher-order chart can also use
real-root barriers. We count exactly `d-2` real characteristic roots in
`(-rho,rho)` at one rational parameter. Nonvanishing characteristic
discriminant and signed barriers at `lambda=+/-rho` keep those roots there
throughout the interval. A pair-product barrier gives `|lambda|>s>rho`.
Vieta's formulas then give the uniform Cauchy weights
`gamma_j=binomial(d-2,j)*(rho/s)^j`. If their Cauchy budget is below one,
`|H(1/lambda)-1|<1`: the constant numerator gives a right-half-plane lift
of `F=N/H` and the exact lower bound `|F|>=N/(1+sum gamma_j)`.
Both endpoint collisions are still checked against the same cofactor
weights, the phase count closes the exact degree, and all earlier rows are
checked exactly. The sample root count becomes an all-interval argument
through the discriminant and barrier certificates. The automatic grid uses
`s=k/16<=1`, `rho/s=j/16` and `q=l/16`; other cases retain the general supplier.

The algebraic supplier reuses resultants of known denominator products by
exact multiplicativity after checking each product identity. Repeated factors may be reduced
for critical-root isolation, while retaining every root and the complete
original resultant in the certificate. Regression checks compare complete
sign certificates, including off-branch boxes, against direct elimination.
The original fraction arithmetic and deflation order are retained for replay.

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

The default pointwise supplier has cutoff 0 for all three examples. The
positive constant row at n=0 follows from the degree certificate. For the
last family, the optional uniform bounds `q=4/5`, `m=13/100` and
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
lookup. If pointwise domination is not certified, the quartic uniform
fallback searches `q=k/16`
for `1<=k<=15`, `m=2^(-k)` for `0<=k<=16`, and each `beta_j` in
`{0} union {2^k : 0<=k<=16}`. It accepts the first certified bound for
each obligation. Exact sample checks only reject candidates; full-interval
sign certificates establish every accepted bound and are retained for replay.
These grids are independent of the family and numerator degree.
For the path-cube example this search finds `q=13/16`, `m=1/8` and
`beta=(1,2,1,1)`, giving uniform cutoff 45. The default exact pointwise
bound instead uses cutoff 0 and two amplitude resultants of degrees 322
and 378.
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
amplitude-zero counts and nongeneric multiple collisions are unsupported.
The positive generic cubic triple point described above is supported.
A constant numerator and a binomial Newton face with `k>4` equal limiting
modes give an exact obstruction to global domination on a negative-infinity
ray. At every local `n=k*m`, the pair has normalized amplitude `2/k` and the
remainder `(k-2)/k`, so no cutoff works. This does not refute real-rootedness;
a section argument, merged pairs or domination only at phase samples may
still apply. At a zero endpoint, a unique minimum valuation/lag face with
`k>4` gives the same obstruction for a constant numerator. At `n=k*m`, the
unique lowest-power coefficient and the limiting residues give the same
pair/remainder ratio. Deflating a zero root rescales both sides of the
inequality equally and cannot remove this obstruction.

Automatic suppliers above order four require compact geometry; their ray
edges and multi-mode pointwise remainder bounds remain unsupported.
Compact charts with eigenvalue arguments pi and zero require degree slope
one. Above order four, a sublinear Newton cap is rejected before expensive
elimination until a higher-order ray supplier is available. Newton layers
whose positivity cannot be certified are rejected. Newton search has 64 layers,
1024 support coefficients and phase period at most 1024; cubic and quartic
Taylor limits use sixteen terms. The pointwise index search stops at 16. General algebraic scaled edge roots and algebraic squared
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
remain unsupported by this simple-root route. Growing zero multiplicities
and zero-eigenvalue finite Newton edges need a valuation-aware conclusion;
the generic nonzero triple-point edge does not cover them.

The quartic proof and its two exact negative controls use the existing
`full-tests` tier because elimination is costly in debug builds. Run their
complete focused suite with
`cargo test --release -p polytool --features full-tests --lib bkw::effective::tests`.

## Stored-witness replay and remaining regeneration

`EffectiveBkwReport::verify` directly checks the stored characteristic and
pair/cofactor factorization in the branch quotient, every recorded sign
resultant, rational sample/root enclosure, critical-root box and wall sign,
nonnegative-cell coverage, and every finite-prefix row identity, degree and
exact distinct-real-root count. These checks consume stored data; they do
not call the sign-grid producer. Invalid intervals and zero spectral
denominators are rejected before quotient arithmetic.

The continuous branch selection, phase endpoint limits, amplitude lift,
Newton degree/leading-sign closure,
domination constants/cutoffs, source-recurrence prefix binding and descriptive
fields still use current deterministic supplier regeneration and complete
comparison. Those supplier-specific obligations do not yet expose a standalone
proof-object interface; separating them requires further certificate protocol
work. Thus this is a partially independent checker, and replay cost can
still include the original supplier search. A supplier bug in these remaining
obligations is not excluded by replay. Stored proofs missing
`domination_cutoff` are rejected explicitly; regenerate them from their request
with the current API or CLI. Omitting the field never enables older guards.
