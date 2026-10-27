# Planar networks

Part of the [CLI reference](cli.md).

## Planar networks: TNN triangles and strip digraphs

`tnn-network` reads monic rows `P_0, P_1, ...` with `deg P_n = n` (or
`P_1, P_2, ...`, in which case `P_0 = 1` is prepended) and decides exactly
whether the coefficient triangle `r_{n,k} = [t^k] P_n` is lower-unitriangular
and totally nonnegative (TNN).  On success it prints the weights
`lambda[n][k]` of the canonical planar network `Gamma_N(lambda)` (the weight of
the vertical edge `(n+1,k) -> (n,k)`; horizontal edges have weight 1) and
re-evaluates the network's path matrix against the input (reconstruct →
evaluate → compare).  On failure it says why: a shape problem (not monic, wrong
degree, entry above the diagonal), a negative entry, or an explicit negative
minor with zero-based rows, columns and exact value.  With `--matrix` each line
is a matrix row instead (a square lower-triangular matrix, or triangular rows
with `n + 1` entries in row `n`).  `--chain` adds the chain polynomials
`p_0 = 1`, `p_n = t * sum_{k<n} r_{n,k} p_k` and checks that each is
real-rooted with zeros in `[-1, 0]` and that `p_n` interlaces `p_{n+1}`.

```sh
printf '1,1\n1,2,1\n1,3,3,1\n' | polytool tnn-network --chain
```

```text
TNN network check
- Input: 3 polynomial row(s); P_0 = 1 prepended
- Coefficient matrix (triangular rows):
  0: 1
  1: 1, 1
  2: 1, 2, 1
  3: 1, 3, 3, 1
- Verdict: lower-unitriangular and totally nonnegative (TNN)
- Canonical network weights lambda[n][k] (edge (n+1,k) -> (n,k)):
  n=0: 1
  n=1: 1, 1
  n=2: 1, 1, 1
- Certificate (reconstruct -> evaluate -> compare): the network's path matrix equals the input

Chain polynomials (Brändén--Saud Maia Leite, arXiv:2412.06595, Theorem 3.7):
  p_0 = 1: real-rooted, zeros in [-1, 0]
  p_1 = t: real-rooted, zeros in [-1, 0]
  p_2 = t + 2t^2: real-rooted, zeros in [-1, 0]
  p_3 = t + 6t^2 + 6t^3: real-rooted, zeros in [-1, 0]
- Consecutive interlacing p_n << p_(n+1): all hold
- Theorem conclusion (real-rooted in [-1, 0] and interlacing): holds

Note: TNN of the coefficient matrix does not imply that the row polynomials are real-rooted: ...
```

The guardrail example `[[1,0,0],[1,1,0],[1,1,1]]` (rows `1`, `1 + t`,
`1 + t + t^2`) is TNN with weights `lambda = [[1], [1, 0]]`, although its last
row has no real zeros.  The report makes no claim about the rows:

```sh
printf '1,0,0\n1,1,0\n1,1,1\n' | polytool tnn-network --matrix
printf '1,0,0\n1,1,0\n2,0,1\n' | polytool tnn-network --matrix --json
```

The second matrix is rejected with the witness `rows [1, 2], columns [0, 1]`,
minor `-2`.  JSON output (schema `polytool.tnn-network.v1`) has the fields
`input_mode`, `input_rows`, `prepended_p0`, `coefficient_matrix`, `tnn`,
`failure` (`kind`, `message`, `row`, `col`, `value`, `witness` with `rows`,
`cols`, `value`), `lambda`, `certificate_verified`, `chain` (`theorem`, `error`,
`polynomials` with `real_rooted` and `roots_in_minus_one_zero`, `interlacing`,
`conclusion_holds`, `contradicts_theorem`) and the fixed `note`.  Integers and
rationals are strings.  Exit status: 0 if the triangle is TNN and the
certificate re-evaluates, 1 if it is not lower-unitriangular TNN, 2 for invalid
input or options.

`brenti-recognize` searches exactly for a planar strip digraph (edges from
column `x` to `x + 1` that depend only on the height, with monotone targets)
whose path counts from `(0,0)` to `(n,k)` are the coefficients of the input
rows `P_1, P_2, ...` (`P_0 = 1` implicit; `--first-row 0` reads `P_0, P_1, ...`).
Bounds mirror `BrentiRecognitionOptions`: `--max-span <k>`,
`--weights rational|integer`, `--max-search-nodes <n>`,
`--max-integer-nodes <n>`.  The verdict is `certified` (the digraph is
re-verified against every row), `impossible_within_bounds` (exact) or
`inconclusive` (a budget ran out); exit status 0, 1 and 3 respectively, and 2
for invalid input or options.

```sh
printf '0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n' | polytool brenti-recognize --weights integer
printf '1,1\n1,4,1\n1,11,11,1\n' | polytool brenti-recognize --json
```

```text
Strip-digraph recognition
- Input: 4 row(s), digraph rows 1..4; P_0 = 1 implicit
- Bounds: integer weights, target span unbounded, 1000000 search nodes, 10000 integer nodes per column
- Verdict: certified
- Digraph (height -> target: weight):
  0 -> 1: 1
  1 -> 1: 1, 2: 1
  2 -> 2: 2, 3: 1
  3 -> 3: 3, 4: 1
  4 -> (no edges)
- Certificate re-verified against the input rows: yes
```

The Stirling rows above are certified; the Eulerian rows are
`impossible_within_bounds` (no nonnegative weights produce target height 1,
even without the planarity constraint).  JSON
output (schema `polytool.brenti-recognize.v1`) has `options`, `rows`,
`verdict`, `certified`, `digraph` (`num_heights`, `outgoing[y]` as
`{to, weight}`), `row_indices`, `certificate_verified`, `supports_checked`,
`obstruction`, `reason`, `message` and the fixed `note`.

Guardrails for both commands:

- TNN of the coefficient matrix does **not** make the row polynomials
  real-rooted: `[[1,0,0],[1,1,0],[1,1,1]]` is TNN and has the row
  `1 + t + t^2`.  Brändén--Saud Maia Leite (arXiv:2412.06595, Theorem 3.7)
  concerns the chain polynomials: they are real-rooted in `[-1, 0]` and
  interlace when the triangle is lower-unitriangular TNN.
- A strip-digraph certificate proves, by Lindström--Gessel--Viennot, that the
  coefficient matrix of `P_0, ..., P_N` is TNN.  It does not prove that the
  rows are real-rooted or Pólya frequency sequences (`[1 + t, 1 + t + t^2]` is
  certified): the transitions depend on the height, so the digraph is
  invariant under shifts in `x`, not in `y`, and Brenti's PF theorem needs
  height-shift invariance.
- Both commands certify only the supplied finite prefix.

## Planar networks: Pólya frequency rows

Two commands certify the *rows* themselves: each certified row is a Pólya
frequency (PF) sequence, i.e. real-rooted with only nonpositive zeros.

`pf-network` searches exactly for a planar network on `{0..N} x Z` that is
invariant under **height** shifts (Brenti, *Combinatorics and total
positivity*, JCTA 1995).  Into each column `n` it has steps
`(length, rise, weight)`, standing for the edges `(n - length, y) -> (n, y + rise)`
at every height `y`, so
`P_n = sum weight * t^rise * P_(n - length)` with `P_0 = 1`.  When every
positive-weight step is primitive (`gcd(length, rise) = 1`) and no two edges
meet away from a common endpoint, the Toeplitz matrix of every row is the path
matrix of fully compatible sources `(0, i)` and sinks `(n, j)`, so it is
totally nonnegative (Lindström--Gessel--Viennot) and the row is PF
(Aissen--Schoenberg--Whitney).  The search tries, column by column, every
inclusion-minimal set of steps that is planar together with the edges already
chosen and admits nonnegative weights for that row; this is exact for the
bounds `--max-length <L>` and `--max-rise <R>` (default 3 and 3).  Other
options: `--first-row <0|1>`, `--weights rational|integer`,
`--max-search-nodes <n>`, `--max-integer-nodes <n>`.  The certificate is
re-verified exactly: structure, planarity, every row, and the path matrix
rebuilt over the explicit lattice against the `(deg + 1) x (deg + 1)` Toeplitz
window of the row.  The PF conclusion comes from planarity and height
invariance, for Toeplitz sections of every size; a finite Toeplitz window that
happens to be TNN would prove nothing (the `3 x 3` window of `1 + t + t^2` is
TNN).  A row that is not real-rooted rules out every network.  Verdicts and
exit codes are as for `brenti-recognize` (0 certified, 1
`impossible_within_bounds`, 2 invalid input, 3 inconclusive).

```sh
printf '1,1\n1,2,1\n1,3,3,1\n' | polytool pf-network
printf '1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n' | polytool pf-network --json
```

```text
Height-invariant planar network (PF rows)
- Input: 3 row(s), network rows 1..3; P_0 = 1 implicit
- Bounds: rational weights, step length <= 3, rise <= 3, 1000000 search nodes, 10000 integer nodes per row
- Verdict: certified (every supplied row is PF)
- Steps into column n (length, rise): weight
  n=1: (1, 0): 1, (1, 1): 1
  n=2: (1, 0): 1, (1, 1): 1
  n=3: (1, 0): 1, (1, 1): 1
- Certificate re-verified (planarity, rows, Toeplitz windows): yes
```

The Eulerian rows `A_1, A_2, A_3` are certified with Brenti's recursive
pattern (steps `(1,0)`, `(1,1)`, `(2,1)`; `A_3 = (1+t) A_2 + 6t A_1`), but no
network with steps of length and rise at most 8 reproduces `A_1, ..., A_4`
(rational or integer weights); the same holds for the Narayana rows
`N_1, ..., N_4`.  Likewise no such network reproduces the Bell polynomials
`B_1, ..., B_4` (Stirling numbers of the second kind), which
`left-production` does certify.  These are statements about the stated
bounds only.  JSON output (schema
`polytool.pf-network.v1`) has `options`, `rows`, `verdict`, `certified`,
`network` (`last_column`, `steps[n-1]` as `{length, rise, weight}`),
`row_indices`, `certificate_verified`, `supports_checked`, `obstruction`,
`reason`, `message` and the fixed `note`.

`left-production` implements the unified planar network of Chen--Fu--Ruan
(arXiv:2512.08369).  Rows `P_0, ..., P_N` with `deg P_n = n` form a lower
triangular matrix `A` with nonzero diagonal; its left production matrix
`Q = A (1 ⊕ A^{-1})` is determined row by row by
`a_{n,0} = q_{n,0}`, `a_{n,k} = sum_{i=k..n} q_{n,i} a_{i-1,k-1}`, and row `n`
of `A` depends only on the leading block `Q_n`.  If `Q_N` is totally
nonnegative, then (Theorems 3 and 4 of the paper, applied to `Q_N ⊕ I`) the
rows `P_0, ..., P_N` are PF and `A_N` and its reversal `[a_{n,n-k}]` are TNN.
The command computes `Q_N` exactly, decides TNN (with a negative-minor witness
on failure), re-verifies the certificate (regenerates `A_N` from `Q_N`,
re-checks TNN, and compares the unified-network window
`M_{n,n}[n..2n; 0..n] = T_n(P_n)^T` for every row), and prints an explicit
planar network for `Q_N = L D` (canonical network weights of `L`, diagonal
`D`).  It covers Stirling numbers of both kinds, Lah, `r`-Whitney and
idempotent numbers, and Brenti's `n`-recursive triangles (Pascal, Delannoy,
derangements) without naming the family.  The criterion is sufficient, not
necessary: `Q_4` of the Eulerian triangle has `q_{4,2} = -5` although the rows
are real-rooted, and the Narayana triangle fails at `Q_4`.  `--first-row 0`
reads `P_0, P_1, ...`; by default `P_0 = 1` is prepended.  Exit status 0 if
certified, 1 if not, 2 for invalid input.

```sh
printf '0,1\n0,1,1\n0,1,3,1\n0,1,7,6,1\n' | polytool left-production
printf '1,1\n1,4,1\n1,11,11,1\n1,26,66,26,1\n' | polytool left-production --json
```

```text
Left production matrix (Chen--Fu--Ruan)
- Input: 4 row(s); P_0 = 1 prepended
- Production matrix Q_N (triangular rows):
  0: 1
  1: 0, 1
  2: 0, 1, 1
  3: 0, 1, 2, 1
  4: 0, 1, 3, 3, 1
- Verdict: Q_4 is TNN; certificate re-verified
- Conclusion: rows P_0..P_4 are PF (real-rooted, nonpositive zeros); the triangle and its reversal are TNN
- Planar network for Q_N = L D (lambda[n][k] of Gamma_N(lambda) for L):
  n=0: 0
  n=1: 0, 1
  n=2: 0, 1, 1
  n=3: 0, 1, 1, 1
  column scales D = diag(1, 1, 1, 1, 1)
```

JSON output (schema `polytool.left-production.v1`) has `first_row_index`,
`rows`, `triangle`, `production_matrix`, `tnn`, `failure_reason`,
`negative_minor` (`rows`, `cols`, `value`), `certificate_verified`,
`rows_certified_pf`, `production_network` (`lambda`, `column_scales`),
`message` and the fixed `note`.

Guardrails for these two commands:

- They certify that the supplied rows are PF.  `pf-network` says nothing
  about the coefficient matrix of the rows (`tnn-network` and
  `brenti-recognize` concern that matrix, and do not certify PF rows).
- `impossible_within_bounds` from `pf-network` is relative to the step
  bounds, unless the obstruction is a row that is not real-rooted.  A
  `left-production` failure says nothing about the rows.
- Both certify only the supplied rows; `1 + t + t^2` is never certified.
