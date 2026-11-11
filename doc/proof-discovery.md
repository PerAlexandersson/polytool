# General recurrence proof discovery

`recurrence-prove` discovers proofs by matching general theorem hypotheses.
It uses the recurrence, its starting index and its seeds. There are no
sequence-name tests or caller attestations in this search; only one stage
(catalog lifts, below) consults the bundled OEIS catalog, and it keeps
a catalog flag apart from a proof.
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
| Iterated preserver | First-order or single-lag rules `P_n = a_n(x) S_n(P_(n-L))`: a real-rooted multiplier after the identity or a preserver of the menu (Laguerre, Pólya–Schur, Rolle, Euler step, Brändén's Weyl step), also in gauge, translation and reversal coordinates; lag-one second-order steps can also be factored into first-order menu or Liu–Wang steps; earlier seed rows checked exactly | Real-rootedness with roots in `(-inf, 0]`, `[0, inf)` or on the line |
| Hypergeometric row | Rows whose coefficient ratio is a product of factors linear in `k` (a row relation proved for every `n`), the image of `(1 + x)^M(n)` under nonnegative multiplier sequences | Real-rootedness with roots in `(-inf, 0]` |
| Fixed-delay window (RealRooted `variableDelay`) | `f(n) P_n = A(n) P_(n-1) + B(n) x P_(n-L)`, `L >= 2`; `A/f > 0`, `B/f >= 0` for every index, constant first row, seed rows as exact window steps | Nonnegative coefficients, roots in `(-inf, 0)`, `P_A << P_B` for `A <= B <= A + L` |
| Residue classes (RealRooted `subseq`) | Every lag a multiple of `p >= 2`; each class `P_(first + j + p m)` proved by one of the rules above (zero classes allowed, but not only zero classes) | Zero or real-rooted rows: the zero classes are named in `conclusion.zero_residue_classes`, every other row is nonzero and real-rooted; hull of the class root regions; no interlacing of consecutive rows |
| Index-section lift (RealRooted `pathPower`) | Multipliers depending on `x` only; `P_(first + s + m) = T(Q_(p m))` for every `m`, with `T` the identity or the reversal `x^(d m + e) P(1/x)`, `p <= 5`, `s` exceptional prefix rows, and a parent `Q` with a constant-coefficient rule found by factoring `chi_T(w^p)`: `Res_w(psi(w), z - w^p) = chi_T(z)` and polynomial parent seeds (an equal-order recurrence argument, nothing fitted); the parent proved by one of the rules above | Nonzero real-rooted rows; the parent region (identity) or its reciprocal image closed at `0` (reversal); no interlacing |
| Factor gauge (RealRooted `powerForm`, `closedForm`, `factorPower`) | `f(n) P_n = sum_j kappa_j(n) q(x)^j P_(n-j)`, order `L >= 2`, one monic `q` (a constant allowed, found exactly as a root of a lag multiplier); `P_n = q^(n-e) g S_n` (`e` the last seed row, `g` the gcd of the gauged window `q^(e-m) P_m`) with `S_n` in the fixed span of the residual window, its coordinates moved by the scalar recurrence. `S_n` real-rooted for every index: degree at most 1; `B^2 - 4AC >= 0` (degree 2); the formal cubic discriminant `>= 0` (degree 3); palindromic quartics through `u = x + 1/x` (outer coefficient of one sign and one of three sign cases). For degree 2 and more the coefficients of `S_n` must be polynomials in `n` satisfying the scalar recurrence identically (otherwise `not_applicable`); `S_n != 0` by a coefficient polynomial, a dominant rational characteristic root, or the shared leading-coefficient certificate; `q`, `g` and the seeds checked exactly | Nonzero real-rooted rows; `(-inf, 0)` or `(-inf, 0]` when every factor and coefficient sign is certified; weak interlacing only for a constant residual and linear `q` |
| Positive lag-two Euler (RealRooted `lagTwoEuler`) | `theta=x D`, `m=n-first_index`: `P_n=[(a+b theta)+x(a+b(m-1)-b theta)]P_(n-1)+c x(theta+1)(m-1-theta)P_(n-2)`; full operator identity, nonzero denominator, `a>0`, `b-a>=0`, `c>=0` decided for every integer; positive constant/linear seeds and exact later seed steps | `n-first_index` simple negative roots, adjacent strict interlacing and auxiliary `(theta+1)(m-theta)P_(n-1) StrictInterl P_n`; replayable normalized coefficients and seed steps |
| Binomial-perturbed Euler orbit (RealRooted issue #2513) | `m=n+offset`, `theta=u D`: `Q_m=[((2m-3)/(m-2))(2u+1)+s(u+1)theta]Q_(m-1)-((m-1)/(m-2))(2u+1)[(2u+1)+s(u+1)theta]Q_(m-2)`; exact rational operator identity, `s>=1`, `kappa>=s+4`, canonical seeds `1`, `1+2u`, `1+kappa*u+kappa*u^2`, rule starts at `m>=3`, all scalar denominators nonzero; also exact Eulerian Mobius transport | `m` simple roots in `(-1,0)` in Q coordinates or `(-inf,0)` in original coordinates; strict directed adjacent interlacing; replayable operator, coordinate and seed certificate |

The binomial-perturbed Euler-orbit route is recorded only when the operators
match, after the existing routes. Its Eulerian Mobius map is exactly
`Q_m(u)=(1+u)^m P_n(u/(1+u))`, with inverse
`P_n(t)=(1-t)^m Q_m(t/(1-t))`. We transport the full rule, check the
inductive degree cap and compare every supplied seed with the canonical
family. The increasing map `t=u/(1+u)` sends `(-1,0)` onto `(-inf,0)`;
degrees, simple roots and directed interlacing are preserved. This proves
A141689 `(s,kappa)=(1,5)`, A141690 `(1,8)`, A141696 `(1,10)`, A141697
`(1,18)` and A147565 `(2,6)`. The first three catalog rules have `m=n+1`;
the last two have `m=n`. Finite exact checks through degree 60 are sanity
checks, separate from the theorem. The reference is
[RealRooted issue #2513](https://github.com/PerAlexandersson/RealRooted/issues/2513);
`lean_declarations` names `RealRooted.strictInterl_of_orbitPerturbation_recurrence`
and, for rows in the Eulerian coordinate, the Cayley-transform interlacing lemma
`RealRooted.strictInterl_cayleyTransform_succ_of_roots_mem_Icc`. Lean replay
remains `not_attempted`.

The positive lag-two Euler route follows the existing routes, preserving
their proof priority. It proves A157147–A157151 (`a=1`, `b=m>=1`, `c=m`)
and A157637 (`c=n/2`). A negative lag correction fails `c>=0`: it does
not prove A157152 (row 8 is not real-rooted) or the open A157153–A157156.

The dispatcher tries the first eight, in the order of the parent cascade of
the lift provers, and records every result; the RealRooted row routes
after them are recorded only when the rule has their form, and the
index-section lift only when it finds a section identity (its search and the
all-index argument are documented in the module `index_section` of
`recurrence::lift::parent_theorems` and in `library.md`), and the factor
gauge only when the rule has its form and the residual shape is decided
(module `factor_gauge`). `proved` means that at least one
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
- otherwise (neither rule divides the other, or a left-hand factor depends
  on `x`), an *invariant module*: the remainder of one rule modulo the
  other, together with relations of the rows fitted exactly to rows (a
  first-order rule and a row annihilator, used only as seeds), is closed
  under the shift in the module of operators of shift below the order of
  the target rule over `Q(n, x)<D>`, kept in an echelon (Hermite) form.
  The generators `G_l` may carry derivatives on any row, including the
  leading one, and the identities `W_i (S o G_i) = Y_i o L + sum_l
  Z_(i,l) o G_l` and `W E = Y o L + sum_l Z_l o G_l` prove by induction
  that every generator, and then the rule, holds for every index once the
  generators kill the rows at a start base and after the bases where a
  multiplier vanishes (exact row checks, recorded in `checked_bases`; the
  generators are in `invariant_module`);
- exact row checks at the integer roots of `w`, at the restart bases of the
  annihilator induction and at the bases before the rule starts, equal seed
  rows, and left-hand factors that are not the zero polynomial where a rule
  computes a row (they may depend on `x`).

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
fitted rules are certified, and the Ore-algebra arithmetic has term,
coefficient-size and work budgets (the invariant-module search stops
early on large intermediate operators). Finding nothing proves nothing.

**Catalog lifts.** When nothing above proves the record (and no root region
is requested), the fourth stage connects it to a family of the bundled OEIS
catalog (919 entries, no external database) by an identity certified for
every index, and transfers that family's statement. Two lifts are tried:

1. *Normalized identification.* A catalog family `S` with
   `T_n = c(n) x^(alpha n + beta) S_(n + s)` for every `n` from some index
   on (`certify_normalized_recurrence_equality`: index shift, dropped
   leading rows, a power of `x` with an affine exponent, a hypergeometric
   scalar proved nonzero for every index). Candidates come from a
   fingerprint index of the catalog, built once per process: the first 16
   rows of every family with the power of `x` removed, denominators cleared,
   the content divided out and a positive leading coefficient, so that the
   scalar and the power of `x` do not change them and the index shift only
   moves them. A family is a candidate when at least 6 of the record's first
   16 rows match at one shift; a family with the record's own rule (up to
   that shift) is skipped. The index is only a hint: each candidate's
   equality is certified, the family is proved from scratch by the stages
   above (never by another lift: depth one), and the proved statement
   "every row from `f` on is nonzero and real-rooted with every root in
   `R`" is carried to the record by
   `transfer_real_rootedness_along_normalization` (a power of `x` adds the
   root `0` to the region; rows the identity does not cover are checked
   exactly). For example A123125 (Eulerian rows with an extra factor `x` and
   a first row `1`) is the catalog family A008292 shifted by a power of `x`,
   which Liu--Wang proves.
2. *Backward basis-map lift.* `find_basis_map_lifts` with the library of
   proved basis maps of the catalog and the published Narayana and Eulerian
   transformations, on the record (on its tail record when discovery has a
   tail): `Q_n = T_P(R_n)` for every `n`, the map theorem sends
   polynomials with roots in a source `S` to ones with roots in a target
   `J`, and the coefficient family `R` is proved to have its roots in `S`
   by the shared parent provers (the same registry theorems). When the
   identity is certified but those provers leave `R` open, `R` is
   identified with a catalog family as in 1, and its transferred region must
   lie in the source of a conclusion of the map theorem. For example A067613
   is `T(R)` for the Hermite basis map of A060821.

The record's statement is "every row from `N` on is nonzero and real-rooted
with every root in the lift's region". Nonzero rows come from the identity
(a nonzero scalar and a power of `x` times a nonzero parent row) for an
identification with a polytool proof; otherwise (an attested parent, which
only states "zero or real-rooted", and the basis-map routes) from the
record's own leading-coefficient certificate (`recurrence-degree`), with the
rows before its `holds_from` checked exactly.

*Polytool proofs and catalog attestations.* A parent is either proved by
polytool (its discovery report or theorem route is replayed) or only
asserted by the catalog flag `lean_proved` of the identified family (label
`oeis_catalog:<A-number>:lean_proved`, a Lean proof made outside polytool),
also through the catalog-flag route of the parent provers for a basis-map
coefficient family. A lift with a polytool parent is preferred. A lift
that rests only on a flag gives the status `conditional_on_attestation`
(exit status 1), never `proved`: `real_rooted_from_index` is omitted and the
conditional claim is `catalog_lifts.conclusion` with its labels in
`conditional_on`.

The stage is bounded: `--max-lift-candidates N` (default 6; 0 disables the
stage) counts the catalog families sent to certification (the record's and
the coefficient families' together, more matching rows first, then catalog
order); the backward search tries at most 512 maps, certifies at most 16
candidates, runs only when the record's first 12 rows have degree at most
24 (exact root location of the expanded rows is its cost), and keeps one
geometric map `g (a + d t)^k` per pair (`g`, sign of `d`), since the lifts
through the others differ from it by an affine change of the variable. The
fingerprint index and the map library are built once per process; on a
record the other stages leave unproved the stage costs about one second.
The report records
every identification (catalog id, the normalized-equality certificate, the
parent's own discovery report, the transfer or the reason none applies),
the basis-map search with its certifications and the coefficient
identifications, and the `conclusion` (`route`: `normalized_identification`,
`basis_map` or `basis_map_with_identified_coefficients`; `catalog_family`
or `basis_map`; `parent_proof`: `polytool` or `catalog_attestation`;
`from_index`, `root_region`, the leading-coefficient certificate and the
exactly checked rows; `statement`).

**Affine coordinates.** When still nothing gives a polytool proof (a
conditional catalog lift does not stop the stage), the criteria run on the
record in exact affine coordinates. For rational `alpha != 0` and `beta`
put `x = alpha u + beta` and `Q_n(u) = P_n(alpha u + beta)`. By the chain
rule (`D_u = alpha D_x`) the rule

```text
f(n, x) P_n = sum_j sign_j(n) c_j(n, x) D_x^(d_j) P_(n - r_j)
```

becomes, for every index,

```text
f(n, alpha u + beta) Q_n = sum_j sign_j(n) alpha^(-d_j) c_j(n, alpha u + beta) D_u^(d_j) Q_(n - r_j)
```

with the seeds `Q_k(u) = P_k(alpha u + beta)`: an identity of the two
records, not a finite check. `P_n` is zero exactly when `Q_n` is, and the
roots of `P_n` are `alpha r + beta` for the roots `r` of `Q_n`, so

- nonzero rows, real-rootedness, zero residue classes and simple roots carry
  over row by row (the summary fields `real_rooted_from_index` and
  `real_rooted_outside_zero_residue_classes` are those of the inner proof);
- a root region `R` of the rows `Q_n` gives `alpha R + beta` (ends swapped for
  `alpha < 0`; an index-dependent end `v + p n` becomes
  `alpha v + beta + alpha p n`, an exact quadratic end is mapped exactly);
- interlacing: for `alpha > 0` the relation of the inner theorem holds for
  `P_(n-1)` and `P_n` (`interlacing_orientation: preserved`). For
  `alpha < 0` the order of the roots is reversed (`mirrored`): in increasing
  order the roots of `P_n` and `P_(n-1)` alternate starting with the
  smallest root of `P_n`. That is `P_(n-1)` interlacing `P_n` when
  `deg P_n = deg P_(n-1) + 1` and `P_n` interlacing `P_(n-1)` when the
  degrees are equal; the stage records the orientation and claims nothing
  more. The inner statement (for example the windows `P_A << P_B` of the
  fixed-delay route) is about the rows `Q_n`.

The candidate coordinates come from the rule alone, in a deterministic
order: the nonzero rational roots `beta` of the factor free of `n` of each
coefficient polynomial (the term coefficients in the order of the rule,
then the left-hand factor), each with `alpha = 1` and then `alpha = -1`,
and last the reflection `x = -u`. A coefficient `x - 1` thus suggests
`u = x - 1`. The iterated-preserver criterion uses the same root helper for
its translation coordinates. At most `--max-affine-coordinates N`
coordinates are tried (default 8; 0 disables the stage); a requested root
region is pulled back to `u`. In each coordinate every criterion runs on
the transformed record, and on the transformed tail record when discovery
has a tail; no other stage runs there (depth one), and the parents of
catalog lifts are proved without this stage. The search stops at the first
coordinate with a proof. On a record that nothing proves, the stage costs
about as much as running the criteria once per coordinate.

The report field `affine_coordinates` holds `options`,
`candidates_available` (the size of the menu), `coordinates` (each with
`coordinate` = `slope`, `intercept`, `reason`; the transformed `recurrence`
and its `attempts`; an optional transformed `tail`; `proved`) and the
`conclusion` (`coordinate`, `theorem`, `on_tail`, `from_index`,
`zero_residue_classes`, `root_region` mapped to `x`, `simple_roots`,
`interlacing_from`, `interlacing`, `interlacing_orientation`, `statement`).
For example A118424 and A118884, `P_n = 2 P_(n-1) + (x - 1) P_(n-L)` with
`L = 3, 4`, are the fixed-delay window `Q_n = 2 Q_(n-1) + u Q_(n-L)` in
`u = x - 1`: every row is real-rooted with every root in `(-inf, 1)`.

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
Exit codes are 0 for a proof, 1 for no proof (including
`conditional_on_attestation`), and 2 for invalid input.
With `--verify`, 0 means the report replays, including a valid `not_proved`
report; the replay output preserves the discovery status.

`--verify` checks a saved report's claims, not byte equality with a
recomputation. The search is rerun on the recorded record and options
(same recurrence, seeds, start delay bound and root region); the recorded
attempts must be of exactly the criteria of the recomputed attempts, in
order, and each recorded attempt must equal its recomputed attempt or
replay on its own (a report of an earlier version of that check). A
recorded attempt that proves nothing where the current check proves is
accepted as superseded; a proved attempt always replays. The tail record is
re-derived from the supplied record, every
alternative's equality certificate and attempts are recomputed, and the
summary fields (`status`, `real_rooted_from_index`,
`real_rooted_outside_zero_residue_classes`, `trivially_real_rooted`,
`statement`) must follow from the recorded
attempts alone, so a report never claims more than its attempts prove.
The catalog-lift stage is replayed from its certificates: every normalized
equality (its source must be the named catalog family's record), every
parent report (a depth-one discovery on that record, replayed itself), every
transfer and its premise (a proved conclusion of the parent report, or the
catalog flag of the named entry), the basis-map search with every
certification, and each coefficient identification against its lift; the
conclusion is then recomputed from the recorded attempts. The
affine-coordinate stage is replayed by recomputing the candidate menu (the
recorded coordinates must be its first ones, in order, up to the bound or
the first proof), every transformed record and the reference attempts on
it; the recorded attempts are replayed like those of the record and the
conclusion is recomputed. A
recomputation that proves more is fine; a report saved before a criterion
was added no longer replays and is recomputed. Missing, extra or repeated
attempts, reordered attempts and attempts flipped to proved are rejected.

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
`0d0bf80e278275b57f8a1ba4b12edce102b56747`; the positive lag-two Euler
route at `5d30bc9fc69d3473c7b7d97afadf04921be45f7f`).

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
| `forall_eq_C_mul_pow_of_rec2`, `forall_eq_C_mul_pow_of_rec3`, `interlaces_of_forall_eq_C_mul_pow` (route `powerForm`), `splits_C_add_C_mul_X_add_C_mul_X_sq`, `splits_palindromic_quartic` (route `closedForm`), `eq_pow_mul_of_factoredThreeTerm` (route `factorPower`) | [Power form](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/ProductSequence/Interlacing.lean#L140), [quadratic and palindromic quartic](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/QuadraticRoot.lean#L160), [factor power](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/Recurrence/FactorPower.lean#L19) | Lean's routes fit a closed form per family (orders two and three, a linear `q`, positive palindromic coefficients); here the gauge, the residual span and the coefficient polynomials come from the rule and the seeds exactly, for any order, with the mirror and straddling palindromic cases, the cubic discriminant (no Lean counterpart) and a dominant-root argument for nonzero rows |
| `exists_mul_add_lt` (route `subseq`) | [Residue reassembly](https://github.com/PerAlexandersson/RealRooted/blob/0d0bf80e278275b57f8a1ba4b12edce102b56747/RealRooted/ThreeTermRecurrence/HalfGrowth.lean#L127) | The class rules in the class index and each class report's own correspondence; Lean's `rr_row_subseq_splits` covers periods 2 to 6 and class rules with one or two lags |
| `RealRooted.lagTwoEuler_recurrence_certificate`, `RealRooted.lagTwoEuler_recurrence_simpleNegRooted` (route `lagTwoEuler`) | [Lag-two Euler recurrence](https://github.com/PerAlexandersson/RealRooted/blob/5d30bc9fc69d3473c7b7d97afadf04921be45f7f/RealRooted/EulerianTriangles/LagTwoRecurrence.lean), [operators](https://github.com/PerAlexandersson/RealRooted/blob/5d30bc9fc69d3473c7b7d97afadf04921be45f7f/RealRooted/EulerianTriangles/LagTwoOperators.lean), [triangle corollary](https://github.com/PerAlexandersson/RealRooted/blob/5d30bc9fc69d3473c7b7d97afadf04921be45f7f/RealRooted/EulerianTriangles/LagTwoTriangle.lean) | Exact theta identity, scalar division and index shift; positive constant/linear seeds and piecewise seed-step scalars; the three rational-function sign certificates. The triangle corollary specializes `a=1`, `b=m>=1`, `c=m lambda_n>=0` |

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
alternative-representation stage (`AlternativeRepresentationOptions`),
`discover_recurrence_proofs_with_stages` also those of the catalog-lift
stage (`CatalogLiftOptions`; the report field `catalog_lifts` is a
`CatalogLiftSearch`), `discover_recurrence_proofs_with_stage_options`
takes the bounds of every later stage (`RecurrenceProofStageOptions`,
including `AffineCoordinateOptions`; the report field `affine_coordinates`
is an `AffineCoordinateSearch`), and
`certify_recurrence_equality` certifies two definitions of one sequence on
its own. `check_product_recurrence_theorem` is also available independently.

MCP tools:

- `discover_recurrence_proofs`: `recurrence_json` as JSON text, optional
  `max_start_delay` (default 3, maximum 16),
  `max_representation_candidates` (default 300, maximum 2000, 0 disables the
  alternative-representation stage) and `max_lift_candidates` (default 6,
  maximum 100, 0 disables the catalog-lift stage) and
  `max_affine_coordinates` (default 8, maximum 64, 0 disables the
  affine-coordinate stage);
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
