# Exponential-family rows

`polytool exp-rows` generates the canonical labelled-block family
`F_n(t)=n![z^n]exp(A(z)+t B(z))`, starting with `F_0=1`.
All arithmetic is exact. Coefficient vectors ascend in powers of `t`.

```sh
polytool exp-rows --marked tail:1 --rows 0..=12 --check interlacing
polytool exp-rows --unmarked sizes:1,2 --marked tail:3 --rows 0..=14 --json
polytool exp-rows --unmarked sizes:1,2 --marked sizes:7 --rows 304,305,1000 --check real-rooted --summary --json
polytool exp-rows --unmarked sizes:1,2 --marked sizes:7 --rows 0..=305 --check first-failure --summary --json
```

Weights mean `A=sum alpha_r*z^r/r!`, `B=sum beta_r*z^r/r!`.
For example `weights:0,2,0,0,0,720` describes `z^2+z^6`, not
`2*z^2+720*z^6`. `sizes:2,6` instead describes `z^2/2!+z^6/6!`.
Sizes are sets; repetitions do not add multiplicity. No size-zero blocks
are permitted. Negative rational weights can be generated, but do not satisfy
the nonnegative proof-route hypotheses.

| Specification | Formal series |
| --- | --- |
| `none` or `{"kind":"sizes","sizes":[]}` | zero |
| `sizes:1,2` | `z+z^2/2!` |
| `weights:1,1/2` | `z+z^2/4` |
| `tail:3` | `exp(z)-1-z-z^2/2!` |
| `{"kind":"exponential","scale":"2","rate":"3","minimum":2}` | `2*(exp(3*z)-1-3*z)` |
| `{"kind":"polynomial_substitution","seed":SPEC,"coefficients":["2","3","1"]}` | `3*SPEC+SPEC^2` |
| `{"kind":"ode_substitution","seed":SPEC,"a":"1","b":"0","lambda":"1"}` | `exp(SPEC)-1` |

A JSON specification can be supplied directly to either weight flag. For
larger inputs use `--input request.json` (or `--input -` for stdin):

```json
{
  "family": {
    "unmarked": {"kind":"sizes","sizes":[1]},
    "marked": {
      "kind":"ode_substitution",
      "seed":{"kind":"sizes","sizes":[2]},
      "a":"1", "b":"2", "lambda":"0"
    }
  },
  "rows":[6,12], "check":"real_rooted", "summary":false,
  "finite_rank":[]
}
```

The CLI flag uses `real-rooted`/`first-failure`; JSON enum values use
`real_rooted`/`first_failure`. JSON input is exclusive with the family, rows,
check and summary flags. The output schema is `polytool.exp-rows.v1`; exact
integers and rationals are strings. `--verify --input report.json --json`
recomputes the full report and rejects any mismatch. A stored report from a
negative check is still replayable successfully. Exit codes: 0 for generation,
passing checks or verified replay; 1 for a requested finite failure; 2 for
invalid input or a replay mismatch. Route obstructions do not change the exit
code of generation-only requests.

Dense ranks stream the insertion recurrence
`F_(n+1)=sum_r binomial(n,r-1)*(alpha_r+t*beta_r)*F_(n+1-r)`.
Finite supports retain a bounded history window. Sparse ranks share the
scalar table `b_j=j![z^j]exp(A)` and avoid earlier polynomial rows. For a
single marked size `k` the exact coefficient is
`n!*beta_k^m*b_(n-k*m)/((n-k*m)!*m!*(k!)^m)`.
For several sizes, truncated products extract `n![z^n]exp(A)*B^m/m!`.
Sparse extraction can still require the scalar prefix through the largest
rank; it does not promise constant memory or sublinear work in that rank.

`--check interlacing` checks real-rootedness and both directed weak
interlacings for **adjacent selected ranks**. A gap such as `305,1000` is
not an adjacent-rank pair and produces no pair check. Either orientation
passes; `null` records an inadmissible directed degree relation. Repeated
roots, constants and zero rows follow the existing root/interlacing API
conventions. `first-failure` checks real-rootedness only and stops at the
first selected failure; a sparse selection makes no claim about omitted ranks.
`--summary` omits row coefficients, retaining the first failure's integer
row, square-free part, degree and distinct PRS real-root count. That exact
failure certificate replays with `ExponentialFailureCertificate::verify`;
it is a root-engine replay record, not a standalone signed remainder chain.

## Symbolic routes and evidence levels

The exponential-series route registry is separate from the recurrence-record
registry: the input is a symbolic series definition, not a fitted recurrence.
`exponential_row_routes` never treats equality of finite prefixes as a formal
series identity. Structured ODE composition itself supplies the exact identity;
recognized exponential tails with threshold one or two supply equivalent ODEs.
Finite coefficient lists are also checked for exact quadratic substitutions of
singleton and nonnegative quartic seeds, by full polynomial multiplication.
Other unstructured compositions remain open.

- `theorem_backed`: nonnegative singleton/binomial seeds; `A=z` with nonnegative
  quadratic, cubic or quartic marked series; threshold-two big blocks; and
  exact `Phi(B_seed)` with `Phi(0)=0`, `Phi'=a+b*x+lambda*Phi`, all parameters
  nonnegative. Both unmarked and inner series must have structurally
  nonnegative coefficients, and the unchanged unmarked series must match
  the seed's hypotheses. These are Rust identity/hypothesis checks against
  public theorems, not a new Lean compilation or axiom audit.
- `paper_proof`: for one positive marked size `k` prime or twice a prime,
  a nonzero unmarked weight at a size not dividing `k` is an all-rank
  obstruction. Passing the divisibility filter supplies only a necessary
  condition. The filter recognizes finite unmarked supports, tails and positive exponential
series, using the shared deterministic primality checker for word-size k.
- `under_review`: the analogous `k=8` necessity is explicitly under review;
  even a violation is not issued as a certified obstruction.
- `finite_rank`: multiplier-symbol or pencil certificates, described below.
  They never establish all-rank real-rootedness.

Public theorem references:
[Substitution.lean](https://github.com/PerAlexandersson/RealRooted/blob/main/RealRooted/ExponentialRows/Substitution.lean),
`RealRooted.ExponentialRows.insertionRow_splits_subst` and
`insertionRow_interl_subst`; the seed declarations
`insertionRow_affine_X_X`, `insertionRow_isRealRooted_threshold_two` and
`insertionRow_strictInterl_threshold_two`; and
[QuarticGrace.lean](https://github.com/PerAlexandersson/RealRooted/blob/main/RealRooted/ExponentialRows/QuarticGrace.lean),
`insertionRow_isRealRooted_quartic_of_nonneg`. The binomial identity extends
algebraically to nonnegative singleton weights, including zero; the cited
seed declaration states its unit specialization. Touchard rows follow from
the monomial seed and `Phi=exp(x)-1`. Quartic splitness alone supplies no
consecutive-interlacing claim. Interlacings are transferred only when the
seed route supplies them.

The paper obstruction source is
[RealRooted #2415](https://github.com/PerAlexandersson/RealRooted/issues/2415).
The new exponential registry distinguishes evidence levels, so the optional
`polynomial_substitution` spec accepts `seed` and ordinary `coefficients` of
`P`. The route requires `A=z`, an admissible seed, positive leading coefficient
and only real nonpositive zeros of `P`, checked exactly. It is labelled
`paper_proof`, with status `paper_proof`, and cites the
[unreviewed round-four argument](https://github.com/PerAlexandersson/RealRooted/issues/2415#issuecomment-7141592218).
Later ODE substitutions retain the paper evidence level; they do not promote
the seed to a theorem-backed result. The `rational_ode_substitution` specification takes the ODE fields plus
`d`, and defines `(1-d*x)Phi'=a+b*x+lambda*Phi`, with zero constant term.
For `d>0`, its route also uses the unreviewed round-four paper argument and
requires `A=z` and nonnegative parameters. The block-derivative coefficients
satisfy `phi_(r+1)=(lambda+r*d)*phi_r+b*[r=1]`, starting at `phi_1=a`.
For example, the following marked specification gives `-log(1-z)`:

```json
{"kind":"rational_ode_substitution", "seed":{"kind":"sizes","sizes":[1]},
 "a":"1", "b":"0", "lambda":"0", "d":"1"}
```

The choice `lambda="1"` instead gives `z/(1-z)`. At `d=0` we use the theorem-backed
ODE route, preserving the seed's evidence level. No paper closure is
registered as a Lean theorem.

## Finite-rank reference certificates

A request's `finite_rank` array can contain entries of this form:

```json
{"index":6,"reference":["1","2","1"],"pencil":null}
```

At the selected rank we set `gamma_i=target_i/reference_i`. Every reference
coefficient must be positive; reference and target lengths must agree.
The reference must be real-rooted. The symbol route tests
`sum_i binomial(d,i)*gamma_i*t^i`, with `d=degree(reference)`, for nonnegative
coefficients and real-rootedness. This is the finite Pólya–Schur multiplier
criterion; see J. Borcea and P. Brändén,
[The Lee–Yang and Pólya–Schur Programs I, Theorem 1.1](https://arxiv.org/abs/0809.0401).
The binomial normalization matters; testing the unnormalized ratios would
be a different and generally invalid sufficient condition.

For a one-free-factor pencil supply, for example,
`"pencil":{"a":"1","b":"1","factor":["1","1"]}`.
We define the diagonal map by `T(t^i)=gamma_i*t^i`. The checker verifies `reference=(a+b*t)S` exactly, `a,b>=0`, and nonnegative,
real-rooted images `T(S)` and `T(t*S)`. Either directed weak interlacing of the two images
certifies the whole **nonnegative** pencil `u*T(S)+v*T(t*S)` at this rank.
See the compatibility criterion in Brändén,
[Unimodality, log-concavity, real-rootedness and beyond, §8](https://arxiv.org/abs/1410.6601).
A failed sufficient test is inconclusive. An invalid identity or undefined
ratio rejects the request. Certificates embed multipliers and symbols/images
and replay with `FiniteRankCertificate::verify`.

## Library and MCP

Library entry points live in `polytool::exponential_rows`:
`BlockWeightSpec`, `ExponentialFamily::{rows_iter,rows,row,selected_rows}`,
`check_exponential_rows`, `exponential_row_routes`, `certify_finite_rank`.
Rational row coefficients use `BigRational`; JSON boundary coefficients use
strings. The existing OEIS catalog and its attestations are unchanged.
A355144 is an independent regression source for threshold-three generation;
its catalog attestation is not generalized to every exponential family.

MCP tool `generate_exponential_rows` accepts `{"request":{...}}` with the
shared request above. It returns the same report, with finite failures in
result fields. `tool_help {"tool":"generate_exponential_rows"}` returns
this guide. Bounds are enforced before row generation/root checking:
JSON <=65536 bytes, 1..512 rows, ranks <=1000, generated degree <=256,
checked degree <=64, total selected coefficient bound <=65536, weight lists
and size sets <=1000, rational text <=80 bytes and numerator/denominator
<=256 bits. Composition depth <=4 and rank <=128; polynomial composition degree <=64. At most 16 reference
certificates, each of degree <=64. Use CLI/library for larger computations.
