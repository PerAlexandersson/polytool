# The real-rooted atlas

`polytool::atlas` and `polytool atlas` maintain a proof DAG of integer
polynomial sequences `P_n`, `n >= n0`, each proved real-rooted, with a root
region, by polytool's own automatic routes.  The data lives in a directory of
plain JSONL files; every proof step carries a replayable polytool certificate,
and `polytool atlas audit` replays all of them from the roots.

## Data model

- **Nodes** (`nodes.jsonl`, schema `polytool.atlas-node.v1`) are sequences,
  not recurrences.  A node stores its first 15 rows (ascending exact integer
  coefficients, trailing zeros removed), its first index, a fingerprint, its
  definitions, a description (name, polytool family, OEIS catalog ids), its
  proved properties and its observations.  A catalog id is listed under
  `oeis_ids` when the catalog rule defines the whole entry, and under
  `oeis_tail_ids` when the entry starts with prefix rows that the rule does
  not produce: the node is then the tail after those rows (for example
  A122431, whose prefix row `3 + 3x + x^2` is not real-rooted while every
  row of the tail `x^n (1 + x)^3` is).
  - The **id** `n-<hash>` is the 64-bit FNV-1a hash of the first definition
    in canonical form: the rule after the left-hand factor normalization of
    the lifts, with its first index and exactly the identifying seeds.
  - **Definitions**: the first one defines the node.  A further definition is
    admitted only with an all-index equality certificate (an edge of
    operation `definition_equality`).  Two definitions that agree on the
    stored rows only stay separate nodes with a `candidate_equal` link.
  - **Properties** are claims with the edges that prove them:
    `real_rooted` (every row from `from_index` on is zero or real-rooted with
    every root in an index-independent region), `interlacing` (consecutive
    rows, with kind and orientation), and the derived claims
    `one_signed_coefficients` and `hadamard_kernel` (Malo's theorem: rows
    with real roots of one sign are Hadamard kernels for real-rooted
    polynomials), and the map property `basis_map_preserver` (the conclusions
    of the basis-map theorem for `u^k -> P_k`).
  - **Observations** are facts about the stored rows only: degrees and the
    degree pattern (`deg P_n = a n + b`, half growth `floor((a n + b)/2)`, or
    irregular), palindromicity, gamma-positivity, nonnegative coefficients and
    the row sums.  They are evidence, not proofs.
- **Edges** (`edges.jsonl`, `polytool.atlas-edge.v1`) are proof steps:
  premises (nodes in the role of a sequence or a map, each with the claim the
  step uses), the conclusion node, the claims, the operation, the certificate,
  the cited theorem, the RealRooted route and declarations where known, and
  the verification state (`polytool: replayed`, `lean: not_attempted`).
  Several edges into one node are all kept.
- **Links** (`links.jsonl`, `polytool.atlas-link.v1`) are evidence that is not
  a proof: `candidate_equal`, `fingerprint_related` (equal up to content,
  sign, reversal or an index shift) and `oeis_match` (triangle, reversed
  triangle or row sums, with the alignment).

Files are append-only: a node is updated by appending a new revision with the
same id (the last revision wins); nothing is deleted.  Ids are content hashes,
so a repeated computation adds nothing.

## Soundness

Each claim of an edge is read off its certificate after the certificate has
been replayed (`atlas::analyze_certificate`):

| operation | certificate | replay | claims and premises |
|---|---|---|---|
| `theorem_check` | `RealRootednessTheoremReport` (three-term, Liu--Wang, derivative lag, Ma--Wang) | `verify()`, must be proved | the conclusion's region (whole line when none) and interlacing; no premises |
| `parent_provers` | the record and its `VeroneseParentProof` | exact recomputation of the provers, proved by a polytool route (never an attestation) | roots in `(-inf, 0]` and the route's interlacing; no premises |
| `preserver_image`, `basis_map`, `gamma_polynomial` | the lift's certification | `verify()` | identity certified for every index and every transfer hypothesis verified; claim: the input region; premise: the parent with the transfer theorem's parent region |
| `veronese_section` | the section rule and its certification | recomputation from the rule and the candidate, then `verify()` | as for the lifts |
| `obreschkoff_combination` | a preserver certification of `P_n + P_(n+1) = T_n(P_n)` | `verify()` and recomputation of `T_n` | real-rooted on the line; premises: the parent's real-rootedness and interlacing |
| `basis_map_theorem` | `BasisMapTheoremReport` | `verify()`, must be proved | the map's conclusions; no premises |
| `definition_equality` | a preserver certification of the identity map | `verify()` | the new definition equals the old for every index |

A lift's own parent proof is never used: its parent must be an atlas node
whose established region lies inside the transfer theorem's parent region,
from an index at most the first parent index the lift uses.  So every claim
rests on polytool certificates alone; Lean-flagged OEIS entries are
benchmarks and match targets, not seeds.

The audit replays the edges in file order, checks the claims against the
certificate, the conclusion's definition against the certificate's sequence
(an exact canonical match, or the same rule with the same rows on the seeds,
which gives the same rows for every later index), every premise against the
claims already established, every node property against its cited edges, the
rows of every node against its first definition (and the id, fingerprint and
observations against their recomputation), and every claim exactly on the
stored rows (exact real-root isolation).  It exits with status 1 on any
problem.  `--sample K` replays only the incoming edges of `K` pseudo-random
nodes, of every node with several incoming edges and of every node in a
candidate-equal link, and takes the other edges as recorded; the row, id,
fingerprint and observation checks and the exact check of every property on
the stored rows still run for every node.

## Commands

```text
polytool atlas seed DIR [--jobs N] [--item-timeout S] [--only LIST] [--limit N]
                        [--no-lifts] [--lift-batch N] [--lift-passes N]
polytool atlas grow DIR [--jobs N] [--item-timeout S] [--batch N]
                        [--max-tasks N] [--max-basis-maps N] [--operations LIST]
                        [--per-operation N] [--parents FILE]
polytool atlas oeis-index --oeis-data OEIS_DIR --cache FILE
polytool atlas match-oeis DIR --cache FILE [--oeis-data OEIS_DIR]
polytool atlas audit DIR [--json] [--sample K [--seed S]]
polytool atlas report DIR [--json]
```

- `seed` tries every named family with a defining recurrence and every
  bundled OEIS catalog entry with a runnable rule.  Each source runs the root
  provers in its own child process, killed after `--item-timeout` seconds
  (reported as `timeout`).  Sources that no root route proves are then tried
  by the lifts (preserver images and gamma polynomials with the atlas nodes
  as parents, basis maps of the catalog), in batches that share the parent
  libraries, with the same per-item limit; a basis-map lift's coefficient
  family is tried by the root provers first so that it can become a node.
  Results are integrated in source order.
- `grow` plans one forward round (`atlas::grow::plan_forward_round`) with a
  budget: preserver images (`P'`, `P + P'`, `(theta + 1) P`,
  `(theta + 1)^2 P`, `P(x + 1)`, and for an affine degree pattern the
  reversal and `(x + 1)^(d_n) P(x/(x + 1))`), basis maps of the first catalog
  maps on derivative-free nodes, gamma polynomials of palindromic nodes,
  Veronese sections `r = 2, 3`, and Obreschkoff combinations of nodes with
  interlacing rows and a first-order rule.  Nodes are taken smallest first.
  Image rules are computed exactly (`image_recurrence_of_operator`,
  `apply_basis_map_to_family`) or fitted and then certified (gamma, Veronese);
  an image whose stored rows are not integers is dropped.  `--operations`,
  `--per-operation` and `--parents FILE` (node ids, one per line) restrict
  the plan, so that one round can be run in parts over the same parents;
  timed-out tasks are listed and not rerun.
- `oeis-index` writes the terms of every entry of a local OEIS copy to a
  cache file (keep it outside the data directory); `match-oeis` matches every
  node's flattened rows (ascending and reversed, with and without the zeros
  below the lowest coefficient, from each of the first rows) and its row sums
  against it.  A hit needs at least 10 agreeing terms with at least three
  distinct values; hits are links, never identifications.
- `report` prints the counts: nodes by source, edges by operation, links,
  OEIS hits by target, nodes matched to entries outside the catalog, and
  nodes with several proofs.

## Limitations of this version

- The Hadamard-factor, row-factor and symmetric-companion lifts are not
  recorded: their certificates name the input rule only by its relation, so
  the atlas cannot tie them to a node's definition by replay alone.
- Multiplier sequences and Schur--Szegő kernels beyond Malo's theorem are not
  yet node properties; the degree pattern is observed, not proved.
- Forward images are not normalized by a gauge: an image with rational rows
  is dropped.
