# Conditional canonical search and the Assaf restoration gap

## Purpose

The local Rust scheduler in `marked_coloring_hrv_bijection.rs` follows the
frontiers in Assaf's proof and produces a replayable certificate whenever its
deterministic local search succeeds.  Small-rank success does not prove that
this particular local scheduler succeeds in every rank.  This note isolates a
conditional construction: if Assaf reachability is valid for the LLT graphs
in question, a canonical finite search makes the rectifier deterministic and
certificate-producing.  The reachability premise itself is not settled here.

The older source is S. Assaf, *Dual equivalence graphs and a combinatorial proof of
LLT and Macdonald positivity*, arXiv:1005.3759v5, especially Lemmas 60--63 and
Theorems 64--65.  In the source these are the four local restoration lemmas,
the restoration theorem, and the final conversion of a D graph into a dual
equivalence graph.  The later paper *Dual equivalence graphs II:
Transformations on locally Schur positive graphs*, arXiv:1704.07039v1, is
decisive for status: it proves only a highest-color transformation theorem,
introduces additional axioms 4'(a,b), and leaves open whether those conditions
are sufficient for the full iterative transformation.  It also describes the
general LLT Schur-positivity problem as open.  Consequently the older Theorem
65 is not an acceptable black box for an unconditional all-rank construction.

## The finite search graph

Fix a finite D graph

\[
  \mathcal H=(V,\sigma,E_2,\ldots,E_{n-1}).
\]

Give `V` its natural total order.  For the LLT application, vertices are
permutations written in lexicographic order.  Define an **edge state** to be a
choice of the colored partial matchings on this same ordered vertex set, with
the same signature function.  Encode a state by the sorted list

\[
  (i,v,E_i(v))\qquad (2\leq i<n,\ v<E_i(v)).
\]

There are finitely many such states.  In fact, axiom 1 determines which
vertices must have an `i`-neighbor, so only finitely many perfect matchings of
those vertices can occur.

Make a directed search graph whose vertices are edge states reachable from
`H`.  An arc is one valid instance of one of Assaf's transformations
`phi_i`, `psi_i`, `gamma_i`, or `theta_i`.  Its label records

1. the transformation kind and color;
2. the distinguished vertex or component;
3. every rooted package isomorphism used by the transformation; and
4. the before/after colored pairs.

The arc is admitted only after checking the hypotheses of the corresponding
definition and checking the resulting edge state.  We may restrict throughout
to the states satisfying axioms 1, 2, 3 and 5, since the construction in
Theorem 65 maintains these axioms.  We must **not** require local Schur
positivity at every intermediate state: Theorem 65 explicitly warns that its
intermediate graphs need not be D graphs.

Order arc labels lexicographically in the order

\[
  \phi<\psi<\gamma<\theta,
\]

then by color, distinguished vertex, and the finite certificate data.  Run
breadth-first search, expanding arcs in this order, until reaching an edge
state satisfying all dual-equivalence axioms.  Call the first terminal state
`Rect(H)` and retain its shortest lexicographically first path as
`Cert(H)`.

## Conditional termination and correctness

**Conditional proposition.** If a dual equivalence edge state is reachable
from `H` by the four admitted transformations, then `Rect(H)` and `Cert(H)`
are defined.

**Proof.** The search has finitely many edge states, so breadth-first search
terminates after visiting each reachable state at most once.  By the premise,
at least one terminal state is reachable.  Breadth-first search consequently
reaches a terminal state, and the fixed orders make both the state and its
certificate canonical.  Every arc can be checked directly from its stored
local package data, so replaying `Cert(H)` independently verifies the result.
∎

This proves termination of the search once reachability is known.  It does not
prove reachability, and therefore is not the missing all-rank proof.  Nor does
it claim that the current, much smaller depth-first local search always
follows a retrospectively repaired history.

## Why retrospective restoration causes no gap

Lemmas 60--63 have the following common form.  If a move of color `i` would
destroy local Schur positivity in one of the relevant four- or five-letter
windows, another specified move of color `i`, `i+1`, or `i+2` can be applied
first, after which the original move can be replayed safely.  Theorem 64 uses these lemmas
retrospectively: it returns to the state just before an offending lower-color
move, inserts the higher-color repair, and then continues.  Independence from
the lower restriction ensures that the already constructed restriction is
unchanged.

The older paper presents this as an existence proof for a path from the
*original* graph, not a
monotone rewrite system on the current final edge state.  Therefore a greedy
current-state scheduler is not formally justified by Theorem 64.  The finite
search explores histories from the original graph, but it succeeds only if
the retrospective argument really supplies a compatible history.  The later
paper's weaker theorem shows that this implication needs an independent proof
for our vertical-strip LLT subclass.

The repair table used in that existence argument is:

| offending move | failed window | repair inserted before it |
| --- | --- | --- |
| `phi_i` | `E_(i-1) union E_i` | `phi_i` at the forced lower witness (Lemma 60) |
| `phi_i` or `psi_i` | `E_i union E_(i+1)` | `phi_(i+1)` (Lemma 61) |
| `phi_i` or `psi_i` | `E_(i-1) union E_i union E_(i+1)` | `phi_(i+1)` or `psi_(i+1)` (Lemma 62) |
| `phi_i` or `psi_i` | `E_i union E_(i+1) union E_(i+2)` | `phi_(i+1)` or `psi_(i+2)` (Lemma 63) |

The first row is the downward four-letter repair; the remaining rows are the
upper four-letter and the two five-letter repairs.

## Exact LLT grading

For a vertical-strip LLT datum, let `H_d` be the induced signed graph on the
standard colorings of q-degree `d`.  The initial LLT involutions preserve
q-degree, so each `H_d` is a union of connected components of the initial D
graph and hence is itself a D graph.  Conditionally apply the canonical
rectifier separately to every nonempty `H_d`.  Whenever every search succeeds,
it follows formally that

- no rectified edge crosses q-degree;
- every terminal standard-dual-equivalence component has a single exact
  q-degree; and
- the transformation certificate never needs an extra assertion that an
  Assaf package happens to respect the LLT statistic.

This componentwise formulation is cleaner than rectifying the full graph and
testing q-degree preservation afterward.

## From a successful rectification to a crystal

Each connected component of `Rect(H_d)` is a standard dual equivalence graph
of a unique partition shape `lambda`.  Order equal-shape components by their
least LLT vertex and use this order as the copy index.  The full colored graph
and signatures determine the signature-preserving isomorphism with
`SYT(lambda)`; equivalently, map the unique dominant-signature vertex first
and propagate along colored edges.  This gives a canonical bijection

\[
  \{\text{standard LLT colorings}\}
    \longleftrightarrow
  \bigsqcup_{d,\lambda,c}\operatorname{SYT}(\lambda).
\]

For an arbitrary positive-integer coloring `kappa`, use decreasing vertex
order inside equal-color ties to obtain its LLT standardization `pi` and its
weakly increasing rank-color word `a`.  Map `pi` to a standard tableau `T`,
then replace rank `r` in `T` by `a_r`.  The descent compatibility of LLT
standardization makes the result semistandard.  Standardizing the resulting
SSYT and replaying the inverse standard-component isomorphism recovers
`kappa`.  Thus this is a genuine reversible bijection, not merely a Schur
coefficient identity, and it preserves the exact color multiset and q-degree.

Finally, transport the usual Kashiwara operators on each SSYT copy through
this bijection.  Because the target is a disjoint union of ordinary SSYT
crystals, the transported operators satisfy all crystal axioms for every
instance whose rectification certificate exists.  An all-rank claim still
requires reachability for every vertical-strip LLT graph.

## What this does and does not settle

This gives an intrinsic conditional algorithm: it is determined only by the
finite LLT signed graph, the natural order on its vertices, and Assaf's four
transformations.  Every successful run produces a complete, independently
replayable certificate, but success in all ranks is not proved.

It is also not the desired short local formula on a coloring.  Canonical
global search may be exponentially large, and the transported operator can
alter colors outside the two colors visible in the SSYT operator.  Before
implementing global search, the prototype must add the later paper's missing
guards: LSP_6 and axioms 4'(a,b), checked after every move.  The policy remains
fail-fast.  The first failed guard or first graph from which no guarded move is
available is the obstruction to analyze; testing a larger rank before that
would add no proof evidence.

The prototype now implements the first part of this audit.  Every accepted
`phi`, `psi`, `gamma`, and `theta` candidate must preserve LSP_4, LSP_5,
LSP_6, and axiom 4'(a).  The complete rank-five suite and the exhaustive
rank-six stress test pass these stronger guards.  This is evidence about those
finite graphs only.  Axiom 4'(b) is deliberately pending because its definition
uses maximal flat chains; substituting the nearby forbidden diagram without
proving equivalence would recreate the same ambiguity this audit is meant to
remove.
