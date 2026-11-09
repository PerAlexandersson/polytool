# General recurrence proof discovery

`recurrence-prove` discovers proofs by matching general theorem hypotheses.
It uses the recurrence, its starting index and its seeds. There are no
sequence-name tests, catalogue lookups or caller attestations in this search.
`recurrence-theorems` is a different command: the conditional
renewal/predecessor-window matcher, whose schema is unchanged.

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
| Liu–Wang three-term | Polynomial lag coefficient, optionally a derivative; degree growth zero or one, appropriate sign and root-location checks; also the degenerate case `A_n = a_n g(x)`, `B_n = b_n g(x)^2` with proportional seeds, where every row is `c_n g^(n-N+2) P_(N-2)` and the check proves `c_n != 0` (a positive sequence, or a hypergeometric closed form `c_n / c_(n-1) = r(n)`, which allows `a_n` to vanish at isolated indices) | Real-rootedness and weak interlacing; root domain when certified |
| Iterated preserver | First-order or single-lag rules `P_n = a_n(x) S_n(P_(n-L))`: a real-rooted multiplier after the identity or a preserver of the menu (Laguerre, Pólya–Schur, Rolle, Euler step, Brändén's Weyl step), also in gauge, translation and reversal coordinates; earlier seed rows checked exactly | Real-rootedness with roots in `(-inf, 0]`, `[0, inf)` or on the line |
| Hypergeometric row | Rows whose coefficient ratio is a product of factors linear in `k` (a row relation proved for every `n`), the image of `(1 + x)^M(n)` under nonnegative multiplier sequences | Real-rootedness with roots in `(-inf, 0]` |
| Fixed-delay window (RealRooted `variableDelay`) | `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)`, `L >= 2`; `A/f > 0`, `B/f >= 0` for every index, constant first row, seed rows as exact window steps | Nonnegative coefficients, roots in `(-inf, 0)`, `P_A << P_B` for `A <= B <= A + L` |
| Residue classes (RealRooted `subseq`) | Every lag a multiple of `p >= 2`; each class `P_(first + j + p m)` proved by one of the rules above (zero classes allowed, but not only zero classes) | Zero or real-rooted rows: the zero classes are named in `conclusion.zero_residue_classes`, every other row is nonzero and real-rooted; hull of the class root regions; no interlacing of consecutive rows |
| Index-section lift (RealRooted `pathPower`) | Multipliers depending on `x` only; `P_(first + s + m) = T(Q_(p m))` for every `m`, with `T` the identity or the reversal `x^(d m + e) P(1/x)`, `p <= 5`, `s` exceptional prefix rows, and a parent `Q` with a constant-coefficient rule found by factoring `chi_T(w^p)`: `Res_w(psi(w), z - w^p) = chi_T(z)` and polynomial parent seeds (an equal-order recurrence argument, nothing fitted); the parent proved by one of the rules above | Nonzero real-rooted rows; the parent region (identity) or its reciprocal image closed at `0` (reversal); no interlacing |

The dispatcher tries the first eight, in the order of the parent cascade of
the lift provers, and records every result; the two RealRooted row routes
after them are recorded only when the rule has their form, and the
index-section lift only when it finds a section identity (its search and the
all-index argument are documented in the module `index_section` of
`recurrence::lift::parent_theorems` and in `library.md`). `proved` means that at least one
attempt checked every hypothesis. `not_proved` means that this search did
not find a proof; failure of a sufficient hypothesis is not a
counterexample. Each attempt records its own starting index, interlacing
range, root domain and unresolved obligations, and the RealRooted
declarations of its route where they are recorded (`lean_route`,
`lean_declarations`). `real_rooted_from_index` is the first index `N` with
every row `P_n`, `n >= N`, proved nonzero and real-rooted.

**Zero residue classes.** A residue class whose seed rows are all zero is
zero at every index (for example `P_n = x^2 P_(n-2)` with `P_0 = 1`,
`P_1 = 0`: every odd row is zero; OEIS A003982, A132918 and A286100 have
such a class). The residue-class conclusion then records
`zero_residue_classes: {period, residues}` (residues of the index `n`
modulo the period, not of its offset from the first index): every row
`P_n`, `n >= conclusion.first_index`, with `n mod period` in `residues` is
zero, and only the other rows from `conclusion.first_index` on are nonzero
and real-rooted. Since no index has every later row nonzero, the
discovery report then omits `real_rooted_from_index` and sets
`real_rooted_outside_zero_residue_classes: {from_index,
zero_residue_classes}`; the status stays `proved` (the attempt proves its
stated conclusion) and the statement names the zero rows. A consumer that
turns a report into the claim "real-rooted from index `N` on" must read
`real_rooted_from_index`, or check
`RealRootednessTheoremConclusion::has_zero_residue_classes` on an attempt.
The field is omitted for every other conclusion, so their JSON is
unchanged. A rule whose classes are all zero is not applicable.

## When no criterion proves the record

**Tail.** When one of the first rows (the seeds and the next
`max-start-delay` rows) is zero or not real-rooted, no criterion can prove
the record from its first index. Discovery then reruns every criterion on
the *tail record*: the same rule, starting after the last such row, with
seeds that are rows of the family (generated exactly where needed), so that
the rule computes the same rows. A proof of the tail is the statement
"every row from index `N` on is nonzero and real-rooted", reported as
`real_rooted_from_index = N` with the tail under `tail` (its `reason` names
the excluded row). Typical inputs are OEIS triangles with prefix rows that
are not real-rooted or are zero. (Zero seed rows that the derivative-lag
and iterated-preserver criteria handle themselves give such a proof
directly.)

**Alternative representations.** When neither the record nor its tail is
proved, the recurrence finder (`recurrence --representations`) runs on
exact generated rows in two bounded passes: first-order rules with
derivatives and a left-hand factor `f(n)` (the Liu–Wang and Ma–Wang
shapes), and rules of lag at most two. A common factor `c(n)` of the
fitted coefficients is cancelled (rows where `c` vanishes stay seeds).
Each fitted rule is only a candidate: it is used once
`certify_recurrence_equality` proves, for every index, that it defines
the family's rows (schema `polytool.recurrence-equality.v1`):

- right pseudo-division in the Ore algebra, `w(n) E = H o L + R`, of one rule
  by the other (`ore-algebra`), so that the rows of `L` satisfy `E` at every
  base where `w(n) != 0`;
- when the remainder `R` is not zero and `L` has order one, an invariant
  row annihilator `A(n)` with `nu A(n+1) o T(n) = B(n) o A(n)` and
  `mu R = Z o A` (as for the Hadamard-factor lifts), so that `R` kills
  every row;
- exact row checks at the integer roots of `w`, at the restart bases of the
  annihilator induction and at the bases before the rule starts, equal seed
  rows, and left-hand factors that do not vanish where a rule computes a
  row.

Either direction may succeed: the family's rows satisfy the alternative
rule, or the alternative rows satisfy the defining rule. A candidate not
certified against the family is compared with the alternatives certified
before it (equality is transitive). The criteria then run on the certified
alternative; the report keeps both certificates (`equality` and
`attempts`). This proves, for example, the convolution triangles A030524,
A049324–A049326 and A383348 (constant-coefficient lag rules) by Ma–Wang
through a derivative rule `n P_n = (a + b n + 2 x) P_(n-1) + x (x + c)
P'_(n-1)`, A321966 and A271704 through first-order and derivative-lag
rules, and A180013 by the product route. The search is bounded:
`--max-representation-candidates N` (default 300; 0 disables the stage)
counts finder configurations per pass deterministically, at most four
fitted rules are certified, and the Ore-algebra arithmetic has a term
budget. Finding nothing proves nothing.

**Not applicable.** A rule with an inhomogeneous term matches none of the
criteria (they all need a homogeneous rule); the status is then
`not_applicable` with `not_applicable_reason`.

**Known limitation: left-hand factors depending on `x`.** A left-hand
factor `x^e f(n)` is normalized exactly when `x^e` divides every
right-hand coefficient. Other left-hand factors that depend on `x` (for
example `prod_j (x + n + 2 + j)`) are not supported: the criteria report
`not_applicable` with that reason, and the equality certificates reject
them. Normalizing by exact division of the whole rule by a factor `h(n, x)`
dividing every coefficient would be the natural extension; no catalog rule
needs it yet.

`--max-start-delay N` allows a theorem induction to start later after checking
the skipped rows exactly (default 3, maximum 16 in CLI/MCP). The library also
accepts a requested root region through `RealRootednessTheoremOptions`.
The report schema is `polytool.recurrence-proof-discovery.v1`.
Exit codes are 0 for a proof, 1 for no proof, and 2 for invalid input.
With `--verify`, 0 means the report replays, including a valid `not_proved`
report; the replay output preserves the discovery status.

`--verify` checks a saved report's claims, not byte equality with a
recomputation. The search is rerun on the recorded record and options
(same recurrence, seeds, start delay bound and root region); the recorded
attempts must appear, equal and in order, among the recomputed attempts
(a subsequence; each criterion at most once), and each recorded attempt
must equal its recomputed attempt or replay on its own (a report of an
earlier version of that check). A recorded attempt that proves nothing
where the current check proves is accepted as superseded; a proved attempt
always replays. The six criteria of the first schema version must be
recorded. The tail record is re-derived from the supplied record, every
alternative's equality certificate and attempts are recomputed, and the
summary fields (`status`, `real_rooted_from_index`,
`real_rooted_outside_zero_residue_classes`, `trivially_real_rooted`,
`statement`) must follow from the recorded
attempts alone, so a report never claims more than its attempts prove. A
recomputation that tries or proves more is fine: reports saved before the
iterated-preserver and hypergeometric-row criteria (or any later one) were
added keep replaying. Extra or repeated attempts, reordered attempts and
attempts flipped to proved are rejected.

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
| `LinearRecurrence.ofMatrix_pow_isSolution_subsequence`, `RealRooted.Tactic.reflectedRecurrence_isSolution`, `RealRooted.Tactic.prependedRecurrence_isSolution` (route `pathPower`) | [Subsequence recurrence](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/Mathlib/Algebra/LinearRecurrence/Subsequence.lean#L121), [path-power route](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/Tactic/Row/PathPower.lean) | Lean's route fixes the parent (path-power clique polynomials, the `CubicParent` instance of Tran--Zumba); here any constant-coefficient parent found by factoring `chi_T(w^p)` and proved by a polytool theorem, whose report is the parent's correspondence |
| `exists_mul_add_lt` (route `subseq`) | [Residue reassembly](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/ThreeTermRecurrence/HalfGrowth.lean#L127) | The class rules in the class index and each class report's own correspondence; Lean's `rr_row_subseq_splits` covers periods 2 to 6 and class rules with one or two lags |

The existing theorem reports retain their paper references and applicable
Lean declaration names. No theorem is labeled Lean-replayed merely because
its general step appears in the formal library.

## Library and MCP

The library entry point is
`polytool::recurrence::discover_recurrence_proofs(&record, &options)`.
It returns `RecurrenceProofDiscoveryReport`; `report.verify()` replays the
recorded routes (above) and every attempt's exact certificates. Altered
conclusions, missing attempts and mismatched seeds are rejected.
`discover_recurrence_proofs_with_options` takes the bounds of the
alternative-representation stage (`AlternativeRepresentationOptions`), and
`certify_recurrence_equality` certifies two definitions of one sequence on
its own. `check_product_recurrence_theorem` is also available independently.

MCP tools:

- `discover_recurrence_proofs`: `recurrence_json` as JSON text, optional
  `max_start_delay` (default 3, maximum 16) and
  `max_representation_candidates` (default 300, maximum 2000, 0 disables the
  alternative-representation stage);
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
