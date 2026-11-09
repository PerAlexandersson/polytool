# OEIS catalog additions, November 2026

The 2026-11-09 update adds 59 families to the catalog, bringing its size from
860 to 919. We used the [catalog maintenance commands](cli-oeis.md) for every
addition. All pre-existing catalog records are unchanged.

## Sources and validation

The Lean source snapshot is
[`sqrt-of-2/real-rooted-oeis-proofs` at `caca416aac4b9dd76e95512da1ef6d9ffd3a39b3`](https://github.com/sqrt-of-2/real-rooted-oeis-proofs/tree/caca416aac4b9dd76e95512da1ef6d9ffd3a39b3).
The definitions are in `ProofsOeisStretch/<A-number>.lean`. We inspected each
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

The bounded searches did not produce a validated rule for these 29 entries:

```text
A039811 A039812 A048743 A055314 A056940 A090657 A090683 A111636
A129062 A130191 A140982 A141618 A141686 A142470 A155467 A155493
A155497 A155516 A155537 A175124 A202750 A245111 A278880 A320280
A326600 A370232 A383149 A383206 A384968
```

A266178 and A266298 have simple recurrences for their polynomial values, but
their OEIS table rows retain increasing numbers of trailing zeros. The current
table evaluator trims those coordinates, so `oeis add` rejected both with
fitted-row mismatches. A regular-triangle reinterpretation would describe
different rows. Supporting explicit irregular padding requires a separate
schema and generator change.

## Polytool proof discovery

With the default `recurrence-prove --json` options, 49 of the 59 added
recurrence records have status `proved`. The other ten are A038455, A067613,
A098173, A108263, A108759, A122765, A140956, A174266, A184366 and A198204.
This dispatcher result is separate from the upstream Lean flags: it proves
properties of the supplied recurrence and seeds, without a catalog attestation.
A122766 is proved through a certified alternative representation, using the
Ma--Wang interval and iterated-preserver routes.
