# Stationary interlacing-matrix orbits

`recurrence::check_matrix_orbit` proves nonzero real-rootedness of the scalar
sequence defined by a recurrence and its seeds. An optional supplied system
has the form

\[
 v_{n+1}=M(x)v_n,\qquad P_n=\ell v_n.
\]

We check a matrix annihilating identity for the scalar recurrence, every
prescribed scalar seed, and the hypotheses of P. Brändén, *Unimodality,
log-concavity, real-rootedness and beyond*, Theorem 8.5 and Lemma 8.4
([arXiv:1410.6601](https://arxiv.org/abs/1410.6601)). The theorem preserves
ordered interlacing vectors with nonnegative coefficients. Nonnegative
constant projections of such vectors are real-rooted by compatibility.

## Realizations and basis changes

Automatic recognition uses no sequence identifiers or catalog lookup. It
matches the exact characteristic polynomial of threshold and
separated-characteristic shapes of dimensions one through three, and solves
coefficient equations for constant/affine seed coordinates. Thresholds range
from zero to the dimension, in increasing order, with markers in
`{0,1,1+x}`. The separated-characteristic parameters range over
`a in {0,1,2}`, `b in {0,-1,-2}`. A companion realization is also tried for
orders one through four.

Automatic projections are coordinate vectors and the all-ones row. The seed
solver chooses one exact particular solution with free variables set to zero;
it does not solve the full cone-feasibility problem for an underdetermined
seed space. Every prescribed scalar row is checked after reconstruction.

For each realization we search all state orderings and nonnegative rational
unit triangular bases. The finite menu contains elementary shears with entry
`1/2`, `1` or `2`, and triangular 0/1 patterns; in dimension five these
patterns have at most two nonzero off-diagonal entries. At most 4096 basis
candidates are checked per request. The diagonal is normalized to one. Finding
nothing means that the bounded search found no proof.

A supplied system may have dimension one through five. Its polynomial matrix
and seed entries have degree at most 24. The scalar recurrence must be
stationary, homogeneous and derivative-free, of order at most eight, with
enough seeds and at most 32 prescribed rows. The successful transformed
matrix must be affine in the polynomial variable.

For a permutation matrix \(R\) and triangular basis \(T\), we check

\[
 G=TRMR^{-1}T^{-1},\qquad w_0=TRv_0,
 \qquad \ell'=\ell R^{-1}T^{-1}.
\]

The entries of \(G\), the coefficients of \(w_0\), and the constant
weights of \(\ell'\) must be nonnegative. The seed is checked as a fully
ordered interlacing vector, using exact root and interlacing decisions.
Positivity of \(T\) alone does not establish positivity of the recovered
projection.

For a recurrence \(P_{n+r}=\sum_{j=1}^r a_j(x)P_{n+r-j}\), the identity

\[
 M^r=\sum_{j=1}^r a_j(x)M^{r-j}
\]

and every prescribed scalar seed are checked over \(\mathbb Q[x]\).
Recurrence uniqueness then identifies every orbit projection with the
recurrence-defined row. A finite fit to another sequence still requires an
independent all-index identity.

## Universal affine pencil decision

For each weakly ordered row pair and column pair, write the corresponding
entries as \(\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\)
and put \(L=\lambda x+\mu\), where \(\lambda,\mu>0\). We decide

\[
 Lb+d\preceq La+c
\]

for every positive parameter pair. Here \(f\preceq g\) means directed
weak interlacing with positive leading signs. Zero is allowed in either
position when its nonzero partner is real-rooted. Repeated index pairs
retain coordinate real-rootedness in singleton or zero-row cases.

All four entries are affine and coefficientwise nonnegative. Lemma 8.4
reduces real-rootedness of the pencils to \(b\preceq d\) and
\(a\preceq c\), which are rational cross-product and degree comparisons.
Their pencil degrees are constant on the positive parameter domain. For
nonzero pencils we require the directed degree difference to be zero or one.

For these real-rooted polynomials, of degree at most two, proper position is
equivalent to nonpositivity on the real line of the Wronskian
\(W(f,g)=f'g-fg'\). This follows by canceling the common real-rooted factor:
the residues of \(g/f\) must be nonpositive, and its polynomial part has
nonnegative slope, so its zeros alternate with the poles in the required
direction. A remaining multiple root would force a Wronskian sign change.

Write

\[
 A=W(b,a),\quad C=W(d,c),\quad K=W(b,c)+W(d,a),
 \qquad bc-da=D_0+D_1x+D_2x^2.
\]

Then the Wronskian has coefficients

\[
 u=\lambda(A\lambda+D_2),\quad
 v=\lambda(2A\mu+K+D_1),\quad
 w=A\mu^2+K\mu+C+\lambda D_0
\]

as a quadratic in \(x\). Nonpositivity is equivalent to
\(u\le0\), \(w\le0\), \(v^2-4uw\le0\), including degree drops.
For every positive parameter pair, these reduce to
\(A,D_0,D_2\le0\) and copositivity of

\[
 Q_2=\begin{pmatrix}-A&-K/2\\-K/2&-C\end{pmatrix}
\]

and

\[
 Q_3=\begin{pmatrix}
 4AD_0&-2AD_1&H/2\\
 -2AD_1&4AD_2&2D_2K\\
 H/2&2D_2K&4D_2C
 \end{pmatrix},\qquad
 H=4AC+4D_2D_0-(K+D_1)^2.
\]

Copositivity means \(y^TQy\ge0\) for every nonnegative vector \(y\).
Homogenization and continuity include the boundary faces of the positive
parameter domain.

We minimize the rational quadratic on the unit simplex. For every nonempty
support \(S\), we solve

\[
 2Q_{SS}y=\tau\mathbf1,\qquad \mathbf1^Ty=1,
\]

and retain solutions with nonnegative coordinates. The minimum of their
values decides copositivity, using at most seven exact linear systems for
\(Q_3\). Singular stationary systems may be omitted: a minimizer of
minimal support cannot have a singular bordered system. Otherwise a kernel
vector \((h,t)\) satisfies \(\mathbf1^Th=0\); multiplying its stationary
equation by the minimizer gives \(t=0\). Moving along \(h\) preserves the
value until a coordinate vanishes, contradicting minimal support.

The certificate retains both rational matrices and their exact simplex
minima. This is a universal decision for affine entries; finite parameter
sampling is used only as a test of the implementation.

## Nonvanishing, replay and conclusions

Since all transformed coefficients and projection weights are nonnegative,
evaluation at one detects nonzero polynomials exactly. Their supports evolve
by Boolean matrix action. We enumerate these finitely many support masks
until the first repeat and require a nonempty projected support at every
step. This proves nonvanishing for every index, including zero coordinates
and periodic support patterns.

`MatrixOrbitReport` embeds the recurrence request and, on success, a
`MatrixOrbitCertificate` with the original realization, permutation, basis,
transformed system, universal pencil checks, support masks and matching Lean
declarations. `verify_matrix_orbit` replays a positive certificate from its
recorded realization without repeating the search. Negative reports are
checked by rerunning their bounded search.

The standard recurrence route is `matrix_orbit`, with additive payload
`certificate.matrix_orbit`. Its conclusion states nonzero real-rootedness
only. It claims neither simplicity nor consecutive interlacing. For example,
the constant all-ones matrix preserves interlacing, but sends the middle
coordinate of

\[
 ((x+4)(x+2),(x+3)(x+1),x(x+2))^T
\]

to \(3x^2+12x+11\), whose two roots both lie strictly between the roots
\(-3,-1\) of the former middle coordinate. The two rows do not interlace.

Where the original matrix has the corresponding shape, metadata cites
`RealRooted.GustafssonSolus.orbit_certificate`,
`RealRooted.HaglundZhang.orbit_certificate`, or
`RealRooted.separatedCharacteristicRows_isPFPolynomial`. These identify
related preservation theorems; `lean_replay` remains `not_attempted`.
