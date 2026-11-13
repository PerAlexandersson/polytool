The data in this directory describes exact algebraic predicates for real
Hessenberg Toeplitz symbols
\[
b(z)=z^{-1}+b_0+b_1z+b_2z^2+b_3z^3.
\]
We distinguish two predicates: reality of the limit spectrum, defined by the
two smallest roots of \(z(b(z)-t)\), and existence of a radius satisfying
circle nonnegativity. An explicit circle certificate applies the theorem
"circle criterion for Hessenberg Toeplitz symbols (unpublished manuscript)":
for a real symbol b, if some R>0 satisfies Im b(R exp(i theta))>=0 for
0<theta<pi, every standard row has simple real zeros and consecutive rows
strictly interlace. The standard rows are defined by P_0=1, P_m=0 for m<0 and
P_n(t)=(t-b_0)*P_(n-1)(t)-sum_(j=1)^s b_j*P_(n-j-1)(t).

We remove \(b_0\) by a real spectral shift. For \(b_1\ne0\), we choose
\(\rho=|b_1|^{-1/2}\) and replace the symbol by
\(\rho(b(\rho z)-b_0)\). Its first tail is
\(\sigma=\operatorname{sgn}(b_1)\). Scaling roots by a positive real number
preserves their modulus ordering, and the corresponding spectral change is
real and invertible. The reflection b(z) -> -b(-z) changes the sign of the second tail and
negates the spectral parameter. It preserves both predicates. We therefore use the rational coordinates
\[
c=\frac{b_2^2}{|b_1|^3},\qquad d=\frac{b_3}{|b_1|^2}.
\]
These coordinates avoid representing an irrational normalized second tail.
For \(b_1=0\), we use \(\rho=1\), \(c=b_2^2\), and \(d=b_3\).
Trailing zero tails reduce the bandwidth before chart selection.

`regions.wls` starts from a tied quadratic factor of the denominator and
encodes dominance of the two smallest roots. `dominant-pair.wls` uses the
reciprocal characteristic polynomial to reduce the quartic query to two
real unknowns. Both scripts include the boundary case of the quadratic
Schur test. `strict-pair.wls` removes the four-equal-modulus branch on the
actual quartic chart: Vieta forces a real spectral parameter there, so this
branch cannot witness an off-real tie. It retains equality in the Schur gap,
including three-root ties. The first script computes the circle predicate independently by minimizing
its polynomial on the interval [-1,1]. Each invocation has a 250-second
internal limit and is intended to run under a 300-second process limit.

`export-regions.wls` replaces radical comparisons by integer-polynomial
sign conditions. It verifies equivalence on the physical parameter chart
before export. `regions-v1.json` records the variable order, normalization,
Mathematica version, exact Boolean expressions, and input/output SHA-256
hashes. Unresolved queries supply no region verdict.

The geometric definition is due to P. Schmidt and F. Spitzer,
*The Toeplitz Matrices of an Arbitrary Laurent Polynomial*,
Mathematica Scandinavica 8 (1960), 15–38,
[doi:10.7146/math.scand.a-10588](https://doi.org/10.7146/math.scand.a-10588).
For the real-limit-spectrum all-row assumption, we use the explicit
Hessenberg restriction of the unproved implication discussed in
B. Shapiro and F. Štampach,
[Correction (2023)](https://doi.org/10.1007/s00365-022-09614-0),
to *Non-Self-Adjoint Toeplitz Matrices Whose Principal Submatrices Have Real
Spectrum*. A region calculation establishes geometry; the separately named
`hessenberg_real_limit_spectrum` hypothesis supplies its conditional all-row
conclusion.

`sample-regions.wls` obtains rational representatives from every returned
CAD disjunct on both sides of each predicate. `cells-v1.json` retains these
cells, points, and the sampler hashes. `check-instances.wls` checks batches
against the original denominator-root query, or a direct quantified circle
polynomial. Row checks use the library's existing exact recurrence generator
and real-root counting.

To regenerate a raw artifact, we run one of the following under a resource
guard, with at most two Wolfram kernels in total:

```sh
timeout 300 wolframscript -file regions.wls s2plus qe-s2plus.json
timeout 300 wolframscript -file strict-pair.wls s3plus qe-strict-s3plus.json
timeout 300 wolframscript -file dominant-pair.wls s3zeroplus qe-s3zeroplus.json
timeout 300 wolframscript -file regions.wls circle3plus qe-circle3plus.json
timeout 300 wolframscript -file export-regions.wls qe-s2plus.json region-s2plus.json
timeout 300 wolframscript -file sample-regions.wls qe-s2plus.json cells-s2plus.json
```

We repeat the appropriate commands for the `plus`, `minus`, and `zero`
charts and both circle bandwidths. `assemble.py ARTIFACT_DIRECTORY` validates
all 15 completed artifacts and their generator hashes before writing the
bundles. Rust evaluates the bundled integer-polynomial expressions; it does
not invoke Mathematica. No bandwidth-four region is bundled.

Related Rust: [`Recurrence::generate_rows_rational_iter`](../../src/recurrence/base/iterate.rs)
is the canonical exact row generator. The integration tests in
[`bkw_regions.rs`](../../tests/integration/bkw_regions.rs) and
[`bkw_certification.rs`](../../tests/integration/bkw_certification.rs)
exercise the region, row, boundary, degeneracy and replay checks.
