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
| `exp(z)-1-z-z^2/2` | exact expression, here a threshold-three tail |
| `coeffs:[0,1,1/2]` | finite ordinary series `z+z^2/2` |
| `coeffs:exponential:[0,1,1]` | finite exponential weights, also `z+z^2/2` |
| `{"kind":"series_expression","expression":"log(1/(1-z))"}` | `sum_{r>=1} z^r/r` |
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

## Exact series expressions

Both weight flags accept expressions directly. JSON uses
`{"kind":"series_expression","expression":"exp(z)-1"}`. Expressions
always use ordinary powers; `weights:` retains its size-one exponential
normalization. Inline coefficient lists explicitly include the constant slot.
They describe finite polynomials, with omitted higher coefficients zero.
`coeffs:ordinary:[...]` is an alias for `coeffs:[...]`;
`coeffs:exponential:[w0,w1,...]` means `sum w_r*z^r/r!`.
Rational list entries may also be quoted strings. Request files use the same
JSON specification and grammar as the flags.

```sh
polytool exp-rows --marked 'exp(z)-1' --rows 0..=20 --check interlacing
polytool exp-rows --unmarked 'z+z^2/2' --marked 'exp(z)-1-z-z^2/2' --rows 0..=14 --json
polytool exp-rows --unmarked z --marked 'z^2*(1+z)^3' --rows 8,12 --json
polytool exp-rows --unmarked z --marked 'compose(exp(x)-1-x,z^2)' --rows 0..=12 --json
polytool exp-rows --unmarked z --marked 'log(1/(1-z))' --rows 0..=12 --check real-rooted
polytool exp-rows --unmarked z --marked '1/(1-z)-1' --rows 0..=12 --json
polytool exp-rows --unmarked z --marked 'z^2*(1+z^4)' --rows 0..=8 --check first-failure
```

The last command stops at rank 6 and prints this final line.

```text
First failure n=6: not real-rooted (degree 3, 1 distinct real root, square-free degree 3)
```

For an interlacing failure, the message identifies the comparison with the
previous row. The root counts in parentheses describe the row at the reported rank.

Grammar (whitespace is ignored between tokens):

```text
sum       := product (("+" | "-") product)*
product   := unary (("*" | "/") unary)*
unary     := ("+" | "-") unary | atom ("^" nonnegative_integer)?
atom      := nonnegative_integer | "z" | "x" | "(" sum ")"
           | "exp(" sum ")" | "log(" sum ")" | "inv(" sum ")"
           | "compose(" sum "," sum ")"
           | "trunc(" sum "," nonnegative_integer ")"
           | "derivative(" sum ")" | coefficient_list
```

`x` and `z` denote the same formal variable; in `compose(outer,inner)`
the first argument is evaluated in its own formal variable and then
substituted. Multiplication is explicit. Powers accept a single nonnegative
integer exponent (`-z^2` means `-(z^2)`); `z^1/2` means `z/2`.
Use parentheses for other grouping. `exp(f)` requires `f(0)=0`,
`log(f)` requires `f(0)=1`, every divisor and `inv(f)` require nonzero
constant term, and composition requires zero inner constant term.
`exp(1)`, `log(2+z)`, `1/z`, and `compose(z,1+z)` are errors.
These restrictions keep every computed coefficient rational and every
operation formal at zero. Final block series must have constant term zero.
Syntax errors include byte offsets; domain errors state the failed condition.

`trunc(f,k)` denotes the **finite polynomial** consisting of degrees zero
through k inclusive, so subsequent operations know all discarded coefficients
are zero. `derivative(f)` computes enough extra terms for the requested
output, including when nested. Ordinary expression evaluation requests only
coefficients through the largest selected rank, plus derivative lookahead;
no infinite series or analytic approximation is constructed.

Default library/CLI expression limits are 8192 bytes, 512 tree nodes,
depth 64, exponent/truncation degree and list degree at most 1024,
evaluation degree at most 4096 (including derivative lookahead), and 65536
bits per numerator/denominator in literals and intermediate expression results.
The reusable evaluator accepts custom `SeriesExpressionLimits`; MCP applies
tighter bounds below. These bounds reject excessive inputs with errors.

## Symbolic routes and evidence levels

The exponential-series route registry is separate from the recurrence-record
registry: the input is a symbolic series definition, not a fitted recurrence.
`exponential_row_routes` never treats equality of finite prefixes as a formal
series identity. Structured ODE composition itself supplies the exact identity;
recognized exponential tails with threshold one or two supply equivalent ODEs.
Finite coefficient lists are also checked for exact quadratic substitutions of
singleton and nonnegative quartic seeds, by full polynomial multiplication.
Expression recognition normalizes bounded exponential polynomials and verifies
full rational or exponential-polynomial differential identities. A finite
prefix may propose ODE parameters, but a full symbolic identity is required
to accept them. `compose(Phi,seed)` supports the existing ODE and polynomial
classes; explicit exponential inner subtrees also recognize expressions such
as `exp(z^2)-1-z^2` by exact substitution. Unmatched free variables reject
that identity. Recognized finite expressions enter the same seed and divisibility
routes as finite weight lists. Recognition is conservative: symbolic
polynomials have degree at most 128, exponential normal forms at most 16 rates,
and symbolic coefficients at most 4096 bits. Exceeding recognition bounds
leaves generation available and supplies no proof route.

Reports add optional `coefficient_signs` entries for expression inputs and
actual negative prefix witnesses. Status `nonnegative` means a structural
all-coefficient sign check succeeded; `negative` includes the first negative
ordinary coefficient and its degree in the inspected prefix; `unknown` means
no global sign proof was found, even if the prefix has no negative entries.
`checked_through` is the inclusive requested maximum rank. Negative coefficients
and unknown signs never justify a nonnegative theorem route. These diagnostics
do not change generation or finite-check exit status.

For example, `--unmarked 'z-z^3' --marked z --rows 3` prints these sign lines.

```text
unmarked coefficient signs: negative through rank 3, first negative coefficient at degree 3: -1
marked coefficient signs: nonnegative through rank 3, no negative coefficient
```

- `theorem_backed`: nonnegative singleton/binomial seeds; `A=z` with nonnegative
  quadratic, cubic or quartic marked series; threshold-two big blocks; and
  exact `Phi(B_seed)` with `Phi(0)=0`, `(1-d*x)Phi'=a+b*x+lambda*Phi`, all
  four parameters nonnegative (including the original `d=0` ODE route);
  powers `B_seed^(2^k)`; and `a*B_seed^2+c*B_seed^4` with `a,c>=0`.
  Both unmarked and inner series must have structurally
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
The exponential registry distinguishes evidence levels, so the optional
`polynomial_substitution` spec accepts `seed` and ordinary `coefficients` of
`P`. The paper route requires `A=z`, an admissible seed, and either positive
leading coefficient with only real nonpositive zeros of `P`, or
`P(x)-P(0)=a*x^2+b*x^3+c*x^4` with `a,b,c>=0`, checked exactly.
Outside the proved faces below, it is labelled `paper_proof`, with status
`paper_proof`, and cites the
[unreviewed round-four argument](https://github.com/PerAlexandersson/RealRooted/issues/2415#issuecomment-7141592218).
Arbitrary powers (for example `B_seed^3` and `B_seed^6`) and the full quartic
`a*B_seed^2+b*B_seed^3+c*B_seed^4` remain paper-only closure cases. The general
PF outer-polynomial route does not claim interlacing transfer. A primitive
quartic seed at `A=z` still has its separate splitness theorem; this does not
prove the full quartic substitution closure for arbitrary admissible seeds.

The proved polynomial faces are recognized before the PF outer-polynomial
test: `x^2+x^4` is not PF, but has the even-quartic closure theorem.
The route identifiers are `power_two_substitution` for exact `B_seed^(2^k)`
and `even_quartic_substitution` for nonnegative `a*B_seed^2+c*B_seed^4`.
The constant `P(0)` is subtracted by the structured specification.
For example, `compose(x^4,exp(z)-1)` and `compose(x^2+x^4,exp(z)-1)`
both have theorem-backed routes. Explicit powers such as `(exp(z)-1)^4`
also retain their seed witness. The declarations are
`RealRooted.ExponentialRows.insertionRow_splits_pow_two_pow`,
`insertionRow_interl_pow_two_pow`, `insertionRow_splits_even_quartic_subst`
and `insertionRow_interl_even_quartic_subst` in
[SubstitutionPowers.lean at b09ae7ff6](https://github.com/PerAlexandersson/RealRooted/blob/b09ae7ff6562bcaada3f998622737f7091feb9f3/RealRooted/ExponentialRows/SubstitutionPowers.lean).

The `rational_ode_substitution` specification takes the ODE fields plus
`d`, and defines `(1-d*x)Phi'=a+b*x+lambda*Phi`, with zero constant term.
For `d>0`, route `rational_ode_substitution` cites
`RealRooted.ExponentialRows.insertionRow_splits_subst_of_rational_ode` and,
when transferring directed weak interlacing,
`insertionRow_interl_subst_of_rational_ode` in
[SubstitutionRational.lean at b09ae7ff6](https://github.com/PerAlexandersson/RealRooted/blob/b09ae7ff6562bcaada3f998622737f7091feb9f3/RealRooted/ExponentialRows/SubstitutionRational.lean).
It requires nonnegative parameters and structurally nonnegative unmarked/seed
coefficients, with the unchanged unmarked series matching the seed route.
The block-derivative coefficients
satisfy `phi_(r+1)=(lambda+r*d)*phi_r+b*[r=1]`, starting at `phi_1=a`.
For example, the following marked specification gives `-log(1-z)`:

```json
{"kind":"rational_ode_substitution", "seed":{"kind":"sizes","sizes":[1]},
 "a":"1", "b":"0", "lambda":"0", "d":"1"}
```

The choice `lambda="1"` instead gives `z/(1-z)`. Expressions `-log(1-z)`
and `z/(1-z)` recognize the same proved rational-ODE routes.
At `d=0` we use the original ODE route. Every proved closure preserves the
seed's evidence level: a paper-only premise keeps the resulting route at
`paper_proof`. Directed weak interlacing transfers only with a cited
interlacing theorem and a seed route that supplies it. Quartic splitness
alone does not establish this premise. These reports check symbolic identities
and theorem hypotheses in Rust; they do not replay Lean proofs.

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
and size sets <=1000, structured-spec rational text <=80 bytes and numerator/denominator
<=256 bits. Composition depth <=4 and rank <=128; polynomial composition degree <=64. At most 16 reference
certificates, each of degree <=64. Expression specs additionally require rank
<=128, input <=2048 bytes, tree <=128 nodes, depth <=16, exponent/truncation
and list degree <=64, derivative lookahead degree <=144, and numerator/denominator
<=4096 bits in literals and intermediate expression results (replacing the
structured-weight rational bounds). Use CLI/library
for larger computations.
