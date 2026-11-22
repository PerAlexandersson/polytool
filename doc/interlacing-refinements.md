# Interlacing refinements

We seek coordinates `v_n=T_n u_n` and a transfer `v_(n+1)=M_n v_n`
that support an interlacing induction. The source is a defining scalar
recurrence, a defining vector Weyl system, or finitely many supplied vectors.
The library entry points are `recurrence::discover_interlacing_refinements`,
`verify_interlacing_refinements`, and `check_affine_interlacing_preserver`.
The CLI command is `interlacing-refine`; the MCP tool is
`discover_interlacing_refinements`.

## Inputs and conclusions

A request has exactly one source:

- `recurrence`: a `polytool.recurrence.v1` record with seeds;
- `source`: a `TransformVectorRecurrenceJson` defining a square homogeneous
  system, with a positive index-only common denominator;
- `states`: vectors of ascending rational-string coefficient rows, beginning
  at `first_index`. Optional `target_rows` supply their proposed collapse.

For a scalar recurrence of lag `k`, the source companion vector consists of
its last `k` seed rows. The report records the index of its newest row as
`first_index`. Its conclusions apply from that index onward.

A candidate has `name`, optional `coordinates`, optional `transition`, and
`projection`. Coordinates use the existing vector-system encoding:
`coordinates[i][j][d].coeffs[a][b]` is the coefficient of `n^a*x^b*D^d`
in `T[i][j]`. Rectangular and singular coordinate maps are permitted.
A missing map means the identity. A missing transfer is fitted under the
specified degree limits. A supplied transfer uses the same vector-system
encoding as `source`, including its exact initial vector.

`projection` gives nonnegative constant weights recovering the target; an
empty list means sum. With a defining system, `source_projection` specifies
the target in source coordinates; an empty list selects the last component.
A zero projection is rejected. Polynomial weights require their own
preservation argument and are not accepted as scalar collapse weights.

Each candidate reports separately:

1. the exact intertwining identity, or its finite-fit status;
2. target recovery;
3. the nonnegative, real-rooted, pairwise ordered initial vector;
4. the initial degree conditions, including zero coordinates;
5. preservation of the ordered PF invariant.

`proved` certifies a real-rooted or zero target for every subsequent index.
`finite_fit` describes supplied rows without an infinite identity.
`unresolved` identifies an outstanding obligation. `rejected` identifies a
failed candidate, with an exact preservation witness when one is found.
Search failure does not refute real-rootedness of the target family.

The report schema is `polytool.interlacing-refinement.v1`. Replay recomputes
its candidates, identities, fits and preservation data exactly. The field
`lean_replay` is `not_attempted`.

## Algebra and search bounds

For denominator-free systems, we check
`T(n+1) C(n) = M(n) T(n)` in the Ore algebra. With source denominator `q`
and transfer denominator `r`, we check `r*T(n+1)*C = q*M*T`;
both denominators are independent of `x`, and `r` must be positive on the
whole integer index ray. Operator products use `D*X=X*D+I`.
Target recovery is a second operator identity.

Supplied candidates precede the automatic proposals. We derive companion
coordinates, last-step differences whose sum recovers the newest row,
derivative coordinates retaining a target component, and sparse polynomial
shears. The third-order template below has a separate exact scale solver.
For a fixed map, coefficient comparison gives rational linear equations
for the transfer. A bounded search over active supports seeks nonnegative
solutions. Row-only fitting reuses the vector recurrence fitter, with its
last transition held out when enough rows are available.

Default limits are dimension 5, 64 candidates, transition degrees 1 in
`n` and `x`, derivative order 1, 500,000 equation cells per fit and 32
active-support attempts. Hard limits are dimension 5, 256 candidates,
index degree 4, variable degree 3, derivative order 2, 2,000,000 equation
cells and 64 support attempts. Input coordinates have coefficient degrees
at most 8; seed polynomials have degree at most 64. At most 128 finite
vectors are accepted. A bounded support search does not decide general
nonnegative feasibility. Rejecting a fitted transfer leaves open the
possibility of another transfer for the same coordinates.

## Preservation

We use the zero-aware convention: nonzero coordinates are real-rooted and
have nonnegative coefficients; every earlier coordinate interlaces every
later one. For a polynomial block `[[a,b],[c,d]]`, the required relation is
`L*b+d << L*a+c` for every `L=lambda*x+mu`, `lambda,mu>0`.
We check all weakly ordered row and column pairs, including repeated
indices. This is the zero-aware form of P. Brändén's Theorem 8.5 in
[Unimodality, log-concavity, real-rootedness and beyond](https://arxiv.org/abs/1410.6601).
Related Lean: `RealRooted.PreservesInterlacingSeq0.iff_nonneg_and_2x2`.

For rational affine entries, the existing exact checker reduces the
Wronskian condition to quadratic copositivity. With polynomial dependence
on `n`, we put `n=first_index+t`. We certify the boundary separately and
expand each copositivity matrix in powers of `t`. Copositivity of every
coefficient matrix proves the condition for all `t>0`. This is a sufficient
uniform test; a failed expansion remains unresolved unless the bounded
rational screen finds a failed pencil. Polynomial entries of higher degree
remain outside this proof engine.

The differential template `[[0,I],[c(n)X,a(n)I+b(n)XD]]`, with
`a,b,c>=0`, has a separate proof. For an ordered PF pair `(f,g)`,
`g` interlaces `X*f` and `theta(g)`. We combine these with `a*g` by the
common-interleaver rule. If `g=0`, the output is `(0,c*X*f)`.
One-coordinate operators `a+b*theta` and `X*(a+b*D)` use Euler/Rolle and
multiplication. These are operator proofs; the polynomial-entry pencil
criterion is not applied to a differential symbol.

## A sharp three-coordinate example

For

```text
P_(n+2)=(a+b*x)P_(n+1)+e*x*P_n+d*x^2*P_(n-1),
```

we use `v_n=(P_n, P_n+kappa*x*P_(n-1), P_(n+1))`. Its transfer is

```text
[          0,             0,       1 ]
[    kappa*x,             0,       1 ]
[ (e-d/kappa)*x, (d/kappa)*x, a+b*x ].
```

The scalar companion indexes this displayed state by its newest row,
whose index is `n+1`.

Universal ordered-PF preservation is equivalent to `a,b,d>=0` and
`a*kappa^2-e*kappa+d<=0`, with `kappa>0`.
To see the remaining inequality, put `u=e-d/kappa` and `v=d/kappa`.
The block with rows 2,3 and columns 1,3 compares
`(lambda+b)*x+mu+a` with `x*(kappa*lambda*x+kappa*mu+u)`.
Its required root order is equivalent to

```text
lambda*(u-kappa*a)+b*(kappa*mu+u)>=0  for all lambda,mu>0.
```

Letting `lambda` increase shows `u>=kappa*a`; this condition also suffices.
The other nonrepeated blocks compare nonnegative affine polynomials with
multiples of `x`, or multiples of `x` with `x` times a positive affine
polynomial. Repeated rows compare the entries in reverse column order;
repeated columns give a common factor `x` times an affine polynomial, or
an affine polynomial. These comparisons hold for `u,v,a,b>=0`.
Thus the displayed inequality, together with the nonnegative entries,
settles every weakly ordered block.
 The function
`solve_scaled_refinement_interval` solves this inequality, including all
zero-parameter cases and exact quadratic endpoints. When `a,d>0`, the
interval exists exactly when `e>0` and `e^2>=4*a*d`.
For `a=b=1,d=e=4`, the unique scale is `kappa=2`; scale 1 fails.
`solve_common_scaled_refinement_interval` intersects the exact intervals
for a finite list of parameter tuples. A particular orbit must also satisfy
the initial ordering condition.
For index-dependent coefficients, a supplied common scale is checked on
the whole index ray by the symbolic preservation engine.

## Seeded half-step windows

The theorem route `half_step_lag` recognizes

```text
A_m=1                       for m<k,
A_(q+k)=(1+xD)A_(q+k-1)+c_q*x*A_q,  k>=1, c_q>0.
```

We check the operator identity, positivity on the entire integer rule
range, unit initial seeds and positive lag steps for every extra seed.
The seeded half-step window theorem gives PF rows, degree `floor(m/k)`,
positive leading coefficient and directed interlacing at every distance
at most `k`. This is a reachable-orbit theorem.
Related Lean: `RealRooted.HalfStepRows.strictInterl_add_of_le`,
`RealRooted.HalfStepRows.isPFPolynomial`, and
`RealRooted.HalfStepRows.natDegree_eq_div`.

The route proves A059022 and A059023 and also applies to A008299,
A124324 and A136011. Its certificate records the lag and the exact extra
seed weights. Normal recurrence proof discovery and replay include it.
