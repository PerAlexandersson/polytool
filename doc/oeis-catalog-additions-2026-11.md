# OEIS catalog additions, November 2026

The 2026-11-09 update adds 59 families to the catalog, bringing its size from
860 to 919. We used the [catalog maintenance commands](cli-oeis.md) for every
addition. All pre-existing catalog records are unchanged.

## Sources and validation

The `lean_proved` flags come from a snapshot of an external collection of
Lean proofs of real-rootedness for OEIS sequences, one file per A-number.
We inspected each
definition and its repository-local import closure for unfinished proofs;
there were no executable `sorry` or `admit` declarations in these closures.
Mentions in comments are not proof holes. No Lean build was performed for this
catalog update.

The OEIS snapshot is dated 2026-06-29. Its supporting b-files were Git LFS
pointers, so the maintenance commands used the complete `%S/%T/%U` prefixes
in the internal-format records. Every available term was compared, including
partial final rows. Each entry records the term source, checked term count,
last complete validation row and SHA-256 digest. No network term retrieval was
used.

We added 54 entries with `oeis add` recurrence search. Each has status
`verified`, at least two complete held-out OEIS rows, recorded fitting bounds
and the polytool version. The separate `recurrence` searches used deterministic
work limits, with successive bounded searches up to 3,000,000,000 work units
and 10,000 configurations per pass. We retained exceptional initial rows as
prefixes where necessary.

Five entries use the documented supplied-recurrence workflow. Their provenance
is `supplied_recurrence`, with status `validated` and zero catalog fitting
holdouts; their separate finder checks do not change that provenance:

| Entry | Rule source | Independent finder holdouts | OEIS validation |
| --- | --- | --- | --- |
| A140956 | A non-alternating rule recovered from the explicit Lean product recurrence | 28 of 42 source-derived rows | All 55 terms; 10 complete rows |
| A154229 | The explicit Lean lag-two rule, with squared triangular-number weight | 14 of 22 source-derived rows, after skipping the first row | All 45 terms; 9 complete rows |
| A154230 | The explicit Lean lag-two rule, with fourth-power-sum weight | 13 of 22 source-derived rows | All 45 terms; 9 complete rows |
| A174266 | The composition of the three declared `A174266RawStep` operators, divided by 6 | 14 of 22 source-derived rows | All 51 terms; 6 complete rows |
| A198204 | The OEIS-fitted second-derivative rule, cross-checked with the Lean bidiagonal form | 2 of 9 OEIS rows | All 45 terms; 9 complete rows |

For the four source-derived supplied rules, the row generator was polytool's exact
recurrence evaluator. Its rows agree with every complete OEIS row before the
independent finder check. For A154230 the finder used three extra fitting rows
to avoid a spurious relation on the shortest solvable prefix.
The targeted A140956 finder pass allowed 10,000,000,000 work units and used
17,697,524.

## Proof flags and definition bridges

We dated 42 new `lean_proved` flags 2026-11-09. A finite OEIS match alone does
not justify a flag. We checked the all-index correspondence with the upstream
definition for 50 additions, retaining eight of those proof results in these
notes without adding their catalog flags, as explained below:

- Explicit recurrence definitions, products and their index shifts account
  for the elementary products, Chebyshev and Favard families, the three
  supplied Lean rules, A198204 and the A36407x families. The latter use
  `A36407xNative_eq_generated` in the upstream source.
- A038455 follows by cancellation in the declared consecutive-factor product:
  its row ratio is `(x+2n)(x+2n+1)/(x+n+1)`. The polynomial divisor is nonzero.
- A140956 uses a non-alternating representation so the catalog remains usable
  by the implicit-relation parent library. For its zero-based product rows,
  eliminating the parity sign gives
  `-(n-1) P_n = x P_(n-1) + n((n-1)^2-x^2) P_(n-2)`.
  The retained four seeds precede the nonsingular computed range.
- A067613 uses two successive Hermite differential steps. Its two-step
  operator is `(2x^2-1) - 2x D + D^2/2`, with the declared parity scale.
- Differentiating the rational generating function of `A122765Core` gives
  `(1-z)/(1-xz+z^2)^2` for A122765 and
  `2(1-z)/(1-xz+z^2)^3` for A122766. Substitution verifies their recorded
  recurrence operators identically, including the initial rows.
- For A060693, A090181, A103371, A118963, A131198, A132813, A178301 and
  A342982, exact factorial-ratio simplification verifies the coefficient
  recurrence identically. We checked the lower and upper coefficient boundary
  identities separately, with out-of-range coefficients zero.
- A098173 follows from the characteristic identity for every fourth binomial
  section, after the reversal and monomial shift in its definition. A131427,
  A253180 and A256061 use the Catalan ratio `(4n-2)/(n+1)`, together with the
  monomial, Touchard or ordered-Touchard recurrence, respectively.

Where a source definition begins after the first catalog row, at most three
leading rows precede it. Those rows are constants, linear polynomials or the
displayed elementary quadratic factors and are real-rooted exactly. Scalar
leaders of the flagged rules have no zeros at computed indices. Polynomial
leaders in the product identities are nonzero divisors, so equality and the
initial rows determine the generated sequence uniquely.

Nine additions deliberately have no proof flag:

| Entries | Reason |
| --- | --- |
| A010816, A089799, A106459, A184366, A197870 | The OEIS row grouping gives constant polynomials with trailing zeros; the upstream definitions give shifted monomials. We did not equate those different row definitions. |
| A108263, A108759 | The irregular-table recurrences passed all finite checks, but we did not establish an all-index identity with the upstream closed forms. |
| A335340 | The upstream definition includes a forced factor of `X` omitted from the catalog's OEIS row coordinates. |
| A337991 | The recurrence passed the finite checks; its all-index equality with the upstream affine-transform definition was not certified. |

Eight further entries with a checked definition bridge also carry the flag:
A049310, A053120, A090181, A098173, A137286, A181332, A244419 and A353279.
Recording them changes no proof: it only adds these ids as near-miss
candidates to the diagnostic text of failed attestation attempts in eleven
existing parent reports (A056242, A098172, A098173, A110319, A143395,
A151509, A151511, A162303, A198204, A321331, A390725), whose statuses and
routes are unchanged; the parent-proof fingerprints were re-blessed
accordingly.

We did not add any `known_not_real_rooted` flag. The known incomplete proofs
A157153, A157154, A157155, A157156, A141689 and A141690 were not among the
additions, and their existing records were left untouched.

## Unadded families

Twenty of the original 29 stretch families remain outside the catalog:

```text
A039811 A039812 A048743 A055314 A090657 A090683 A111636 A129062
A130191 A141618 A141686 A155467 A155493 A175124 A245111 A320280
A326600 A383149 A383206 A384968
```

Here D-finite means P-recursive in the row index over `Q(x)`. A mixed
shift–derivative row recurrence is a weaker condition. We classify nineteen
of these families as non-D-finite in this sense:

| Entries | Obstruction |
| --- | --- |
| A048743, A055314, A090657, A090683, A141618, A141686, A155467, A155493, A245111, A384968 | Their fixed-column ordinary generating functions have poles at reciprocals of unbounded integers. Stirling or Eulerian column formulas give the poles; for A384968 the largest exponential base in column `k` of the Bell-row square is `floor(k^2/4)`. We first reverse A055314 in degree `n`. |
| A039811, A039812, A130191, A383206 | A fixed-column exponential generating function is an iterated exponential of infinite order (or its derivative). Entire D-finite functions have finite order. |
| A129062, A320280, A383149 | Fixed-column exponential generating functions contain `log(2-exp(t))` or `log((3-exp(2t))/2)`, giving infinitely many finite singularities. |
| A111636 | At positive rational `x`, a middle term grows as `exp(c n^2-O(n))`, exceeding the `C^n (n!)^s` bound for P-recursive sequences. |
| A326600 | Its leading coefficient is the Bell number `B_n`. Reversing in degree `n` and extracting the constant coefficient gives a sequence whose exponential generating function has infinite order. |

For the first group, an annihilating operator over `Q[x,z]`, after removing
its common power of `x`, would force every fixed-column generating function
to have finite singularities in one fixed finite set. The unbounded sets of
poles contradict this. Fixed-column extraction and degree-`n` reversal
preserve P-recursiveness, which gives the other obstructions.

A175124 is D-finite by finite hypergeometric sums on each parity and the
gamma transform. Its available relation includes derivatives of the current
row through order four and has only finite validation. We retain no catalog
entry for it.

## Padded-table additions, 2026-11-22

A266178 and A266298 are now included, bringing the catalog to 921 entries.
Both use `layout: table` with explicit coordinate lengths
`row_lengths: {first: 1, step: 2, every: 1}`. Row `n`, counted from zero,
therefore retains all `2n+1` OEIS coordinates, including trailing zeros.
Existing catalog records and generated rows are unchanged.

Both entries were added with the supplied-recurrence maintenance command.
Their common seeds are `P_0 = 1` and `P_1 = 1+x`. A266178 has
`P_n = P_(n-2)` for `n >= 2`; A266298 has `P_n = P_(n-1)` for `n >= 2`.
These are exactly the upstream Lean definitions, with identical indices and
unit leaders. Induction gives equality for every index; zero padding changes
no polynomial coefficient. Both carry `lean_proved`, dated 2026-11-22.
Neither upstream theorem file nor its repository-local import closure
contains executable `sorry`, `admit` or a custom `axiom`. This source audit
does not constitute a fresh Lean build.

Every one of the 10,000 locally available consecutive b-file terms matched
for each entry: 100 complete rows, of which 98 were computed by the
recurrence. Each b-file also agrees with the corresponding 1,000-term local
snapshot prefix and the internal-format `%S/%T/%U` prefix. No network
retrieval was used. Provenance is `supplied_recurrence`, status `validated`,
with zero fitting holdouts, `terms_source: bfile`, `checked_terms: 10000`,
the exact term digest, and the full last compared row (`n = 99`). These
finite OEIS checks are separate from the all-index definition bridge and
the recorded upstream proof flag.

## Stretch additions, 2026-11-22

We add A056940, A140982, A142470, A155497, A155516, A155537, A202750,
A278880 and A370232, bringing the catalog from 921 to 930 entries. We use
`oeis add --recurrence-json` for each entry. All 921 earlier catalog records
remain byte-identical. Each new record has source `supplied_recurrence`,
status `validated`, zero fitting holdouts, `terms_source: bfile`, and its
checked term count, term digest and last complete validation row.

The b-files used here reproduce every available term of the local OEIS
snapshot. We also compare the generated rows with internal-format prefixes
made from that same snapshot; these are two input formats for the same
terms, rather than independent data sources. Partial final rows are included.
We use no network retrieval.

| Entry | Terms checked | Complete rows | Partial final terms | Lag / derivative order | Default `recurrence-prove` status |
| --- | ---: | ---: | ---: | --- | --- |
| A056940 | 1000 | 44 | 10 | 1 / 3 | `not_proved` |
| A140982 | 48 | 9 | 3 | 1 / 3 | `not_proved` |
| A142470 | 1000 | 44 | 10 | 1 / 3 | `not_proved` |
| A155497 | 1000 | 44 | 10 | 1 / 3 | `not_proved` |
| A155516 | 1000 | 44 | 10 | 1 / 3 | `not_proved` |
| A155537 | 1000 | 44 | 10 | 4 / 0 | `not_proved` |
| A202750 | 231 | 21 | 0 | 1 / 3 | `not_proved` |
| A278880 | 1000 | 44 | 10 | 1 / 3 | `not_proved` |
| A370232 | 55 | 10 | 0 | 2 / 1 | `not_proved` |

For eight entries, exact factorial-ratio simplification proves the
coefficient identity of the supplied recurrence, with the low and high
coefficient boundaries checked separately. For A155537 we use the Narayana
recurrence
`(n+2)N_n = (2n+1)(1+x)N_(n-1) - (n-1)(1-x)^2 N_(n-2)`
and eliminate the two states `N_n` and `2^(n+1)N_n`. Their sum equals the
upstream definition. Every identically vanishing polynomial leader occurs
within the prescribed initial rows. The nonsingular recurrence and matching
seeds therefore identify all rows with the upstream formulas. Exact replay
through native row 80 provides an additional finite check.

All nine carry the existing `real_rootedness` attestation fields
`status: lean_proved`, `as_of: 2026-11-22`. Their upstream proofs are Lean
files in a separate formalization corpus, at a fixed revision.
We inspected all nine files and their locally available import closure:
none contains executable `sorry`, `admit` or a custom `axiom`. The direct
file hashes agree with the source snapshot used for the recurrence
identities. We performed no new Lean build; these flags record the upstream
proofs and the all-index definition bridges.

## Polytool proof discovery

With the default `recurrence-prove --json` options, 49 of the 59 added
recurrence records have status `proved`. The other ten are A038455, A067613,
A098173, A108263, A108759, A122765, A140956, A174266, A184366 and A198204.
This dispatcher result is separate from the upstream Lean flags: it proves
properties of the supplied recurrence and seeds, without a catalog attestation.
A122766 is proved through a certified alternative representation, using the
Ma--Wang interval and iterated-preserver routes.

The nine stretch additions each return `not_proved` with the default
`recurrence-prove --json` options. We retain their upstream Lean
attestations separately. Adding these flags leaves every pre-existing
discovery report unchanged; 21 parent-proof reports gain near-miss catalog
labels in failed attestation diagnostics. The older fingerprints also omitted
labels for the already catalogued A266178 and A266298. Removing the nine
stretch labels and these two padded-table labels recovers every older parent
fingerprint exactly; statuses, routes and all other report fields are unchanged.
