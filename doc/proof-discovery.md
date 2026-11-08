# General recurrence proof discovery

`recurrence-prove` discovers proofs by matching general theorem hypotheses.
It uses the recurrence, its starting index and its seeds. There are no
sequence-name tests, catalogue lookups or caller attestations in this search.

```sh
polytool recurrence --json < rows.txt > recurrence.json
polytool recurrence-prove --recurrence recurrence.json --json > proof.json
polytool recurrence-prove --verify --json < proof.json
```

For a representation menu, select a candidate explicitly:

```sh
polytool recurrence-prove --recurrence menu.json --candidate 0
```

The first command fits a recurrence to finite data. A successful proof in the
second command concerns the sequence **defined by that recurrence and its
seeds**. To apply it to an independently defined combinatorial sequence, we
still need an all-index recurrence identity and agreement of the seeds.
Neither fitting nor held-out checks prove that identity.

## Available general rules

| Rule | Input shape and main checks | Conclusion |
| --- | --- | --- |
| Product induction | `B(n) P_n = A(n,x) P_(n-k)`, `k >= 1`; nonvanishing scalar denominator; uniformly nonzero real-rooted multipliers; every seed nonzero and real-rooted | Nonzero real-rooted rows; no automatic interlacing |
| Three-term/Favard | Affine preceding-row coefficient and scalar lag-two coefficient, also supported reciprocal coordinates; positive normalized recurrence coefficients and interlacing seeds | Real-rootedness, simple roots and strict interlacing; optional root region |
| Liu–Wang first order | `P_n = A_n P_(n-1) + V_n P'_(n-1)`; degree, sign, seed and root-location induction | Real-rootedness, negative roots and weak interlacing |
| Derivative-lag | An additional `W_n P_(n-2)` term; degree/sign conditions, including supported half-growth patterns | Real-rootedness and the reported interlacing/root domain |
| Ma–Wang interval | Derivative coefficient `kappa_n (x-L)(x-U)`; endpoint conditions with multiplicities | Real-rootedness in a fixed interval and weak interlacing |
| Liu–Wang three-term | Polynomial lag coefficient, optionally a derivative; degree growth zero or one, appropriate sign and root-location checks | Real-rootedness and weak interlacing; root domain when certified |
| Fixed-delay window (RealRooted `variableDelay`) | `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)`, `L >= 2`; `A/f > 0`, `B/f >= 0` for every index, constant first row, seed rows as exact window steps | Nonnegative coefficients, roots in `(-inf, 0)`, `P_A << P_B` for `A <= B <= A + L` |
| Residue classes (RealRooted `subseq`) | Every lag a multiple of `p >= 2`; each class `P_(first + j + p m)` proved by one of the rules above (zero classes allowed) | Zero (named) or real-rooted rows, hull of the class root regions; no interlacing of consecutive rows |

The dispatcher tries the first six and records every result; the last two
are recorded only when the rule has their form. `proved` means that
at least one attempt checked every hypothesis. `not_proved` means that this
search did not find a proof; failure of a sufficient hypothesis is not a
counterexample. Each attempt records its own starting index, interlacing
range, root domain and unresolved obligations.

`--max-start-delay N` allows a theorem induction to start later after checking
the skipped rows exactly (default 3, maximum 16 in CLI/MCP). The library also
accepts a requested root region through `RealRootednessTheoremOptions`.
The report schema is `polytool.recurrence-proof-discovery.v1`.
Exit codes are 0 for a proof, 1 for no proof, and 2 for invalid input.
With `--verify`, 0 means the report replays, including a valid `not_proved`
report; the replay output preserves the discovery status.

## Degree laws

The Liu--Wang, derivative-lag and Liu--Wang three-term criteria need a degree
law with leading coefficients of one sign.  Each first tries its own exact
degree step; where that is too narrow (a negative lag multiplier of the top
coefficients, top terms that cancel, a degree-step table outside the
supported cases) it takes the hypothesis `degree_law` from the shared
leading-coefficient module (`certify_leading_coefficient`, also the command
`polytool recurrence-degree`), whose certificate is recomputed when the
report is replayed.  The hypergeometric-row criterion does the same for its
first coefficient in the reversal coordinates (the leading coefficient,
hypothesis `first_coefficient`).  Reports proved by the local steps are
unchanged.  The module works without any real-rootedness hypothesis, so it
can also be run first to see the degree schedule, the exceptional initial
rows and the leading-coefficient signs of a fitted recurrence.

## Product certificates

For each recurrence index, the multiplier checker supports:

- `g(n) h(x)`, with arbitrary-degree fixed `h` checked exactly and `g(n)`
  proved nonzero on the entire integer range;
- `a(n)x^2+b(n)x+c(n)`, including linear and constant cases, using
  `b(n)^2-4a(n)c(n) >= 0` and `a(n)^2+b(n)^2+c(n)^2 > 0` throughout the
  range; this allows a quadratic to drop to a nonzero linear or constant;
- products of these factors discovered by exact bivariate factorization
  when the leading coefficient in `x` is a nonzero constant.

All coefficients are rational. Sign statements use exact integer-range
certificates, not a finite sample of indices. For example, the multiplier
`x^2+n-10000` fails at `n=10001`, even though the first 10000 steps pass.
Quadratics such as `x^2-n` need no rational linear factorization.

Automatic bivariate factorization is limited to degree 12 in each variable
and 256 recombination subset trials. An unsupported factor or exhausted
factorization budget leaves an unresolved obligation. Nonconstant leading
coefficients are supported directly for degree at most two and for separable
multipliers; general higher-degree factorization over `Q(n)` is not provided.
The current conclusion requires nonzero rows, so zero seeds and zero
multipliers are not accepted by this rule. Repeated roots are allowed.

Requested open root intervals can be checked for separable factors and
linear/quadratic factors of fixed degree. For a quadratic with real roots,
the vertex and endpoint evaluation decide whether both roots lie beyond an
endpoint. Degree drops with a requested interval may remain unresolved or
fail a sufficient fixed-degree hypothesis even when real-rootedness alone
is proved. In the existing parent-transfer pipeline, the new product rule
is tried directly with negative roots required; a proof of arbitrary real
roots alone is not used for a transfer requiring nonpositive roots.

## Lean correspondence and evidence

The report records `lean_replay: "not_attempted"`. Rust checks and certificate
replay do not generate a Lean kernel certificate. Existing Lean declarations
are reusable building blocks; normalizations, exact arithmetic certificates,
index shifts and inductions still require formal replay when one is wanted.
The inventory below was source-inspected at RealRooted revision
`1fc0e36cc67a28fa0fa0be4a1bab818926e020fd` (the two row routes at
`0d0bf80e278275b57f8a1ba4b12edce102b56747`).

| General building block | Primary Lean source | Additional correspondence work |
| --- | --- | --- |
| `RealRooted.productSequence_splits` | [Product induction](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/Interlacing/Euclid.lean#L118) | Scalar division, residue classes and multiplier certificates |
| `RealRooted.favardInterlacing` | [Favard](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/Favard.lean#L21) | Its normalized constant/linear seeds cover a subcase; general seeds, reciprocal coordinates and region certificates need their step arguments |
| Derivative recurrence interlacing | [Derivative recurrence](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/DerivativeRecurrence/Interlacing.lean) | Degree, leading signs and root-location induction; do not substitute a global-sign theorem for a sign condition only at roots |
| `derivLag_strictInterl_of_nonnegCoeffs` and Liu–Wang steps | [Derivative-lag](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/DerivativeRecurrence/QuadraticLagStrict.lean) | Report-specific `halfGrowth`/`sign` conditions and any power-of-x normalization |
| `affineEndpointStep_closed` and endpoint multiplicities | [Ma–Wang endpoint step](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/MaWang/EndpointTwoSided.lean#L298) | Affine coordinate change and multiplicity-counter induction |
| Three-term sign and window steps | [Three-term recurrence](https://github.com/PerAlexandersson/RealRooted/blob/1fc0e36cc67a28fa0fa0be4a1bab818926e020fd/RealRooted/ThreeTermRecurrence/Interlacing.lean) | Degree and root-region invariants appropriate to the selected route |
| `strictInterl_of_scaled_fixedDelay_recurrence`, `hasNonnegCoeffs_of_scaled_fixedDelay_recurrence` (route `variableDelay`) | [Fixed-delay window](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/Interlacing/VariableDelay.lean#L280) | Division by the constant first row, the seed rows as window steps, the sign certificates of `a` and `b` |
| `exists_mul_add_lt` (route `subseq`) | [Residue reassembly](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/ThreeTermRecurrence/HalfGrowth.lean#L127) | The class rules in the class index and each class report's own correspondence; Lean's `rr_row_subseq_splits` covers periods 2 to 6 and class rules with one or two lags |

The existing theorem reports retain their paper references and applicable
Lean declaration names. No theorem is labeled Lean-replayed merely because
its general step appears in the formal library.

## Library and MCP

The library entry point is
`polytool::recurrence::discover_recurrence_proofs(&record, &options)`.
It returns `RecurrenceProofDiscoveryReport`; `report.verify()` recomputes the
complete search and replays every attempt's exact certificates. Altered
conclusions, missing attempts and mismatched seeds are rejected.
`check_product_recurrence_theorem` is also available independently.

MCP tools:

- `discover_recurrence_proofs`: `recurrence_json` as JSON text, optional
  `maximum_start_delay` (default 3, maximum 16);
- `verify_recurrence_proofs`: `report_json` as JSON text.

Both tools use the same library code as the CLI and retain the MCP input
limits. The existing `recurrence-theorems` command remains the separate
renewal/predecessor-window matcher, with its existing schema and evidence labels.

## Extending coverage

New discovery rules should match a general structural hypothesis and emit
replayable obligations. Named families belong in regression tests only.
Further useful general routes include derivative-of-product induction,
Euler-bidiagonal degree patterns, matrix interlacing with a closed refinement,
and compositions of preserving transformations. Several of these already
have specialized library adapters; composing their certificates into this
direct search and formal Lean replay are further work. This command does not
yet dispatch every theorem or transformation in RealRooted.
