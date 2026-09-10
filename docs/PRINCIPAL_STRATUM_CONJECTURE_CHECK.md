# First conjecture checks

2026-09-10. Host supervisor owns this note and the corresponding handoff entry;
ownership is released after the documentation checkpoint. No source edits,
worker launches, remote jobs, manuscript edits, or outbound mail.

We inspect Boris's forwarded email UID 32, attachment
`principal_poset_resonance_v6_stronger.tex`, using its internal labels below.
This is a first pass, not a verification of every theorem in the draft.
The background cellular model is from Katz--Shapiro--Welker,
[arXiv:2112.15205](https://arxiv.org/abs/2112.15205).

## A counterexample to the repeated-3 formula

For `omega = (3,1,1,3,1,1,3)`, equation `eq:H-repeat3` proposes

\[
\mathcal H_\omega(s,t)=
\frac{s^6t^{19}(1+s^2t^2+s^4t^4)(1+st^2)^3}
     {(1-s^2t^4)^4}.
\]

Its coefficient at ambient degree 23 is
`7 s^8 + 3 s^9 + s^10`: the numerator contributes
`3 s^8 + 3 s^9 + s^10`, and the denominator contributes another `4 s^8`.
We instead compute the reduced integral groups

\[
\widetilde H_8=\mathbb Z^7,\qquad
\widetilde H_9=\mathbb Z^2,
\]

with all other groups zero and no torsion. Thus `eq:H-repeat3` is false.
Both polynomials specialize to 5 at `s=-1`; their difference is
`s^9(1+s)`. The signed enumerator cannot detect this difference.
This is a numerical difference, not a construction of a cancelling pair
of classes in the conjectured model.

The earlier degree 21 case gives `H_7 = Z^3`, `H_8 = Z`, agreeing with
the draft and strengthening its reported finite-field evidence to an exact
integral calculation.

| Pattern | Ambient degree | Cells | Nonzero reduced integral groups |
|---|---:|---:|---|
| `(3,1,1,3,1,1,3)` | 21 | 54,234 | `H_7 = Z^3`, `H_8 = Z` |
| `(3,1,1,3,1,1,3)` | 23 | 193,695 | `H_8 = Z^7`, `H_9 = Z^2` |
| `(3,3,1,5)` | 28 | 242,788 | none |

The last row agrees with `eq:H-3315`: neither `28-20=8` nor `28-24=4`
is divisible by 6, so its proposed series predicts zero. This single check
does not prove that formula.

These results do not disprove the general bigraded resonance conjecture,
torsion-freeness, or the pure-factor lift: the repeated-3 Euler series has a
nonconstant numerator and is outside the stated pure-factor subclass.

## Verification and reproduction

We use the existing reviewed release executable from source through
`b134824`; the repository was at `ce0c05b` before these documentation edits.
We do not build or change code. Each run is serial and bounded by
`timeout 60s nice -n 10`, with these common arguments:

```text
--max-cells 250000 --max-nnz 4000000
--max-reduction-nnz 4000000 --max-pivots 125000 --record-certificate
```

Docker executable:
`/cargo-target/ai-projects/release/principal_stratum_homology`.
Example host command:

```bash
docker exec docker-setup-app-1 timeout 60s nice -n 10 \
  /cargo-target/ai-projects/release/principal_stratum_homology \
  --omega 3,1,1,3,1,1,3 --d 23 \
  --max-cells 250000 --max-nnz 4000000 \
  --max-reduction-nnz 4000000 --max-pivots 125000 --record-certificate
```

All three integral runs validate the chain condition, compare unreduced
modular homology over `F_251`, replay unit cancellations, and verify recorded
Smith operations when needed. The degree 23 counterexample is independently
rechecked by the unreduced field path over `F_2`, with the same ranks.
These are independent reduction paths within the same implementation, not
an independent implementation of the cell boundary.

For degree 23, initial NNZ is 1,601,252; 96,830 unit pivots leave chain ranks
2, 20, 13 in degrees 7, 8, 9, with 151 residual nonzeros. Exact Smith
reduction gives the groups above. Certificate replay matches the residual.
The displayed SHA-256 value binds the input; it is not an exported proof file:
`ff8769b385dc66192d4d3c56ff0196926a16893e978e5690f726dc8e2b4e804d`.
The commands reproduce the in-run verification; no certificate export was
added. Degree 21 reduces to four isolated generators; `(3,3,1,5)` at
degree 28 reduces to the empty complex by 121,394 unit pivots.

Related Rust: `experiments/src/principal_stratum.rs`;
`experiments/src/bin/principal_stratum_homology.rs`;
`combinatoric-core/src/chain_complex.rs::integral_homology`.

## Part of absorption already follows from the draft

Part (i) of `conj:absorption` says that, for `d > a+1`, the inclusions
from the singleton `(a+1)` into `(1,a)` and `(a,1)` induce homology
isomorphisms. Assuming the later theorems of the draft, this is a corollary:

- If `a` is odd, the proof of `thm:two-odd-resonance` applies to `(1,a)`.
  It identifies the singleton stratum as the closed collision locus and
  its complement as `R^2` times the space of nonnegative monic polynomials
  of degree `d-a-1 > 0`. That complement is Borel--Moore acyclic, so the
  long exact sequence gives the isomorphism induced by the inclusion.
  Reversal gives `(a,1)`.
- If `a` is even, `thm:extreme-even-acyclic` makes both two-part strata
  acyclic. The singleton has odd weight and is acyclic for `d>a+1` by
  `thm:odd-singleton`. The inclusions therefore induce isomorphisms
  between zero groups.

The strict degree bound matters: the statement fails at `d=a+1`, as the
draft observes. This argument does not settle part (ii), concerning
`(1,a,1)`. In particular, matching abstract homology groups alone would not
establish that its specified inclusion induces an isomorphism.

## Remaining priorities and formulation issues

1. Determine the actual primitive data for the repeated-3 family. We should
   inspect a structured reduced differential or filtration, rather than fit
   another numerator to just three coefficients. More exact degrees can
   reject candidates, but cannot by themselves prove a rational formula.
2. Test the pure-factor families separately. The previously verified
   `(3,1,1,5)` degree 26 result agrees with its formula. Its degree 30
   prediction is a useful later target, subject to a separate cell/memory
   estimate; we did not run it in this first pass.
3. For anchored support, an exact proof that the signed enumerator is not
   identically zero for every anchored pattern would suffice. This is a
   sufficient route, not an equivalence: an identically zero Euler series
   does not force all homology to vanish.
4. The finite automaton proves rationality of cell counts and the Euler
   specialization, not rationality of the homology series. The missing
   ingredient is control of differential ranks or a finite chain-level
   structure preserved by the proposed stabilization operators.
5. In `conj:pure-lift`, `(-1)^h` determines only the parity of `h`, whereas
   the conclusion uses its actual value. The initial homological degree
   needs an independent definition. A normalization of the period list is
   also needed: `Phi_2(t)=1`, so adding a period 2 changes the proposed lift
   without changing its Euler hypothesis.

No claim about the full torsion question, unbounded pole order, or the
remaining local absorption law is made by this first pass.
