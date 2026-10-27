# Ehrhart conversion, cyclic sieving and Stapledon decomposition

Part of the [CLI reference](cli.md).

## Ehrhart ↔ h\*-vector conversion

```sh
# h*-vector → Ehrhart polynomial
echo "1, 8, 35, 32, 9" | polytool hstar-to-ehrhart

# Ehrhart polynomial → h*-vector (coefficients as rationals: num/den)
echo "1, 2, 2" | polytool ehrhart-to-hstar
```

For `hstar-to-ehrhart`, h\*-vector entries may be arbitrary-size integers.
For `ehrhart-to-hstar`, input coefficients may be exact rationals and the
resulting h\*-vector entries are arbitrary-size integers.

Check named h\*-vector inequalities:

```sh
echo "1, 4, 1, 0" | polytool hstar-inequalities --dimension 3
echo "1, 20, 1, 0" | polytool hstar-inequalities --dimension 3 --json
```

The report includes basic Ehrhart h\*-conditions, Stanley and Hibi inequalities,
the Balletti--Higashitani universal Scott inequality when applicable, and
Stapledon-derived nonnegativity checks.  Failures include the named inequality
and a reference.

## Cyclic sieving checks

For one polynomial and one cyclic group order:

```sh
echo "1, 1" | polytool cyclic-sieving --order 2 --fixed-counts "2,0"
```

Without fixed counts, the command profiles exact root-of-unity evaluations.  For
a sequence, the default candidate orders for row index `n` are
`n-2,n-1,n,n+1,n+2,n+3`:

```sh
polytool cyclic-sieving-sequence --first-index 2 < polys.txt
```

Optional fixed-count files use lines of the form:

```text
5 7: 10, 0, 2, 0, 2, 0, 0
```

### Orbit decomposition and the Alexandersson--Amini criterion

Every `cyclic-sieving` run also decomposes the polynomial into cyclic orbits,
following P. Alexandersson and N. Amini, *The cone of cyclic sieving
phenomena*, Discrete Math. 342 (2019) 1581--1601 (arXiv:1804.01447).  All
arithmetic is exact and integral.  Reduce `f` modulo `q^n - 1` to residues
`c_0, ..., c_{n-1}`.

* `f(zeta_n^d)` is an integer for all `d` iff `c_j = c_{gcd(j,n)}` for all `j`
  (Desarmenien); otherwise the non-constant gcd classes are listed.
* If integral, `f = sum_{r | n} a_r h_r (mod q^n - 1)` with
  `h_r = 1 + q^r + ... + q^{n-r}` and `a_r = sum_{e | r} mu(r/e) c_e`.  A CSP
  action needs exactly `a_r` orbits of size `n/r`, i.e.
  `S_k = k a_{n/k} = sum_{j | k} mu(k/j) f(zeta_n^j)` elements in orbits of
  size `k`, and `f(zeta_n^d) = sum_{(n/r) | d} (n/r) a_r`.
* For `f` in `N[q]` that is integral, a `C_n`-action making `(X, C_n, f)` a CSP
  exists iff all `a_r >= 0` (`thm:cspCond`).  The natural-coefficient hypothesis
  and the cone condition are reported separately.

The paper's counterexample `f = q^5 + 3q^3 + q + 10`, `n = 6` has natural
coefficients and integral evaluations `15, 8, 12, 5, 12, 8` (`d = 0..5`) but
`S_k = 8, 4, -3, 6` for `k = 1, 2, 3, 6`:

```sh
echo "10, 1, 0, 3, 0, 1" | polytool cyclic-sieving --order 6
#   orbit size 3: -1 orbit(s), -3 element(s)  [a_2 = -1]
#   verdict: no C_n-action makes (X, C_n, f) a CSP
```

To check a concrete action, give its multiset of orbit sizes instead of fixed
counts.  Rotation of 2-subsets of `[4]` has orbits of sizes 4 and 2 and
`[4 choose 2]_q = 1 + q + 2q^2 + q^3 + q^4`:

```sh
echo "1,1,2,1,1" | polytool cyclic-sieving --order 4 --orbit-sizes 4,2
```

With `--json`, each item keeps its `report` object unchanged and gains
`orbit_decomposition` (residues, gcd classes, `integral`, `violating_classes`,
`orbits` with `orbit_count`/`element_count`, `fixed_points`,
`natural_coefficients`, `cone_condition`, `realizable`) and `orbit_sizes_check`
(`null` unless `--orbit-sizes` is given).

Library API (`polytool::`, all `BigInt`): `cyclic_sieving_orbit_decomposition_bigint(coeffs, n)`,
`cyclic_sieving_check_orbit_sizes(coeffs, n, orbit_sizes)`,
`orbit_counts_from_fixed_points(fixed, n)` (Moebius inversion; errors if the
table is not that of an action), `orbit_counts_from_orbit_sizes`,
`fixed_points_from_orbit_sizes`, `fixed_points_from_orbit_counts`,
`reduce_modulo_cyclic_bigint`.

### Bicyclic sieving

For `G = C_{n1} x C_{n2}` and `f(q, t)` (Barcelo--Reiner--Stanton), the
polynomial is a JSON coefficient matrix on stdin, `f[i][j]` being the
coefficient of `q^i t^j`; entries are integers or decimal strings:

```sh
echo '[[1,0],[0,1]]' | polytool bicyclic-sieving --order1 2 --order2 2 \
  --fixed-table '[[2,0],[0,2]]'
```

The command reduces `f` modulo `(q^{n1} - 1, t^{n2} - 1)`, evaluates it exactly
at all `(zeta_{n1}^{d1}, zeta_{n2}^{d2})` and compares with the table
`fix(g1^{d1} g2^{d2})`.  If `c_{j1,j2} = c_{gcd(j1,n1),gcd(j2,n2)}` the
polynomial lies in the product basis `h_{r1}(q) h_{r2}(t)`; its coordinates
`a_{r1,r2}` (two-dimensional Moebius inversion) are reported together with the
evaluations they imply.  **They are not orbit counts and may be negative**: for
the diagonal action of `C_2 x C_2` on two points, `f = 1 + qt` has coordinates
`a_{1,1}, a_{1,2}, a_{2,1}, a_{2,2} = 1, -1, -1, 2`.  Integrality at all points
is reported separately from product form, because they differ in general when
`gcd(n1, n2) > 1` (for `n1 = n2 = 3`, `sum (zeta_3^{i+j} + zeta_3^{-(i+j)}) q^i t^j`
takes only the values 0 and 9 on `mu_3 x mu_3` but is not in the product
span).  Library: `bicyclic_sieving_report_bigint(matrix, n1, n2, fixed_table)`
(limit `n1 n2 <= 2500`).

## Stapledon decomposition

```sh
echo "1, 11, 11, 1" | polytool stapledon 3
```

CLI Stapledon input coefficients may be arbitrary-size integers.
