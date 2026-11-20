# Exact univariate stability

The `exact-stability` crate and the `polytool::stability` facade count roots
of a polynomial over the Gaussian rationals, ℚ(i), in circular regions.
Coefficients are in ascending degree order. All decisions are exact;
root approximations and numerical tolerances are unnecessary.

## CLI and MCP

```sh
polytool stability '(1+t)^3' --region left
polytool stability 't^2+(1+2I)t-I' --region upper --json
polytool stability 't-I' --mobius '1,0,0,1' --convention inside
polytool stability 't-I' --disk-center I --radius 1 --certificate report.json
polytool stability verify --input report.json
```

Positional arguments form a batch; otherwise read one polynomial per line
from stdin or `--input FILE`. Blank/comment (`#`) lines are ignored.
Both adapters use `exact_stability::parse`: expressions or ascending lists
such as `[-I,1+2I,1]`, with integer/rational parts, `I`, parentheses, implicit
products and natural powers. Floating/scientific literals, nonconstant
denominators and multivariate expressions fail. Use `--` before negative
positional expressions. Lowercase `i` and `j` are reserved: use `I` for the
imaginary unit and another name for the variable. Expression nesting is bounded by 64.

`--region upper|lower|left|right|disk|exterior` defaults to upper. It is
exclusive with `--mobius A,B,C,D` (image of the upper half-plane under
`(Az+B)/(Cz+D)`) and paired `--disk-center CENTER --radius R` selectors.
Maps must be nonsingular and radii positive exact rationals.
`--convention avoid|inside` overrides the region default.

`--format summary|json|full` defaults to readable summaries. `--json` aliases
`--format json`, emitting `{items:[...]}` with zero-based index, degree,
counts, verdict, meaning and method. Zero explicitly has null degree/counts,
passes avoidance by convention and fails inside. Constants pass both.
`full` and `--certificate FILE` include canonical ascending input/options and
all certificate entries under schema `polytool.stability.v1`. Exact integers
and rational parts are strings. `stability verify` recomputes every field;
exit statuses are 0 verified, 1 mismatch, 2 malformed input.

The schema identifier selects the replay contract. The `version` field is a
nonempty producer-version string for provenance and may differ from the
verifier's package version. Replay still compares every mathematical field and
rejects extra fields. Incompatible changes to the wire structure or canonical
certificate computation require a new schema identifier; unsupported schemas
are rejected rather than interpreted using the current contract.

`--max-degree N` and `--max-coefficient-bits N` can lower the CLI
defaults/maxima 2000 and 1000000. Bits apply to each part's numerator and
denominator, including intermediate parser results. CLI text is capped at
16 MiB. Dense exact arithmetic can
be expensive below these admission bounds.

MCP `stability_univariate` accepts `polynomials` (1–512 expressions or
ascending arrays), optional `region`, `mobius` (four Gaussian values), paired
`disk_center`/`radius`, `convention`, and `detail:summary|full`. Gaussian
values are JSON integers, exact strings or objects such as
`{"real":"1/2","imaginary":"-3/4"}`. JSON floats fail; object parts and
radius are rational strings. Degree is capped at 32, part bits at 256,
total coefficients at 65536, request text at 1 MiB and coefficient literals
at 16384 bytes. A 50-second elapsed-time guard rejects the whole response at
checkpoints within each item: Möbius substitution, gcd reductions, repeated-root
layers and between counting stages, including before and after inertia. It does
not interrupt individual arithmetic, Sturm or inertia operations. Expiry never
returns partial counts or a verdict.

Default items are about 150–200 bytes with no certificate. `detail:"full"`
includes evidence. Summary mode constructs no full JSON report and stores no
evidence. With `detail:"full"`, configured storage adds `report_id` and `report` metadata;
the stored canonical report always contains full evidence and can be replayed
by the CLI. `read_report` selects evidence by JSON pointer; `tool_help` with
`{"tool":"stability_univariate"}` explains the complete contract.

## API and conventions

```rust
use polytool::stability::{
    check_univariate_stability, GaussianRational, RootCounts, StabilityRegion,
};

let polynomial: Vec<_> = [1i64, 3, 3, 1]
    .into_iter().map(GaussianRational::from).collect();
let report = check_univariate_stability(&polynomial, &StabilityRegion::Left).unwrap();
assert_eq!(report.counts, Some(RootCounts { inside: 3, boundary: 0, outside: 0 }));
assert!(report.verdict);
assert!(report.verify(&polynomial));
```

The entry points are:

```text
check_univariate_stability(&[GaussianRational], &StabilityRegion)
    -> Result<StabilityReport, StabilityError>
check_univariate_stability_with_convention(
    &[GaussianRational], &StabilityRegion, StabilityConvention
) -> Result<StabilityReport, StabilityError>
check_univariate_stability_with_checkpoint(
    &[GaussianRational], &StabilityRegion, StabilityConvention, &impl Fn() -> bool
) -> Result<StabilityReport, StabilityError>
StabilityReport::verify(&self, &[GaussianRational]) -> bool
```

`GaussianRational::new(real, imaginary)` takes two `BigRational` values.
Integer and real rational conversions are available through `From`.
`checked_divide` returns `None` for a zero divisor. The public parts are
`real` and `imaginary`; no decimal approximation is needed.

| Region | Open interior | Default verdict |
| --- | --- | --- |
| `Upper` | Im z > 0 | No roots inside |
| `Lower` | Im z < 0 | No roots inside |
| `Left` | Re z < 0 | All roots inside (Hurwitz) |
| `Right` | Re z > 0 | All roots inside |
| `Disk` | absolute value of z < 1 | All roots inside (Schur) |
| `Exterior` | absolute value of z > 1 | All roots inside |
| `CenteredDisk { center, radius }` | absolute value of z − center < radius | All roots inside |
| `Mobius(Box::new(map))` | map of the open upper half-plane | No roots inside |

`Exterior` is the open complement of the closed unit disk; the unit circle
is reported separately as boundary. The radius must be strictly positive.
`StabilityConvention::{Avoid, Inside}` overrides the defaults without changing
any counts. Counts always include root multiplicity. Boundary roots pass
avoidance and fail all-inside. Every nonzero constant passes both conventions.

The zero polynomial vanishes everywhere and has no finite root multiset.
Its report has `degree: None`, `counts: None`, `certificate: None` and method
`ZeroPolynomial`. Its avoidance verdict is true by the polynomial-stability
convention; its all-inside verdict is false, agreeing with the strict Hurwitz
API. This exception is explicit in `StabilityReport::meaning`.

Trailing zeros are ignored. There are no library size limits: callers should
bound degree and coefficient bit size before requesting expensive computations.

## Algorithm

For a nonzero polynomial p of degree n, choose a Möbius map
M(w) = (aw+b)/(cw+d) from the upper half-plane to the requested region.
The named maps are respectively w, −w, iw, −iw, (w−i)/(w+i) and
(w+i)/(w−i). A centered disk uses center + radius (w−i)/(w+i).

1. Homogeneous Horner substitution computes
   q(w) = (cw+d)^n p(M(w)). The degree drop counts the root multiplicity
   at the boundary point M(∞) = a/c; these roots must not be discarded.
2. Multiply by the conjugate leading coefficient and jointly clear real
   and imaginary denominators with a positive multiplier. Write the integer
   result as g+ih, where the imaginary part has smaller degree.
3. A primitive integer pseudo-remainder gcd extracts the common real factor
   d. The existing `exact-real-roots` Sturm/PRS and squarefree kernels count
   its real roots with multiplicity: each successive layer divides by its
   squarefree part and lowers every multiplicity by one. The remaining roots
   of d occur in conjugate pairs, one root in each half-plane per pair.
4. Divide g and h by d and clear denominators jointly. Build the Bézoutian
   with orientation (g(x)h(y)−g(y)h(x))/(x−y). Its negative inertia counts
   upper-half-plane roots and its positive inertia counts lower-half-plane
   roots of the coprime quotient.
5. Add the common factor's conjugate pairs to both open-region counts, and
   its real roots and the Möbius degree drop to the boundary count.

Real input in the upper or lower half-plane needs only Sturm counting:
the imaginary part is zero, so the coprime quotient is constant and the
Bézoutian is empty. For strict left-half-plane verdicts on real integer
input, `polytool::hurwitz::is_hurwitz_stable_bigint` remains an independent
Routh–Hurwitz criterion. Differential tests check their agreement.
Schur–Cohn is used only as a test oracle, not as another production engine.

## Möbius maps and replay

`MobiusTransform::new(a, b, c, d)` rejects zero determinant exactly.
`inverse`, `compose`, `evaluate` and `transform_polynomial` expose the
same coordinate conventions as the engine. `evaluate` returns `None` at a
finite pole. `StabilityRegion::membership_sign` returns +1 inside, 0 on the
boundary and −1 outside, using the inverse map and exact arithmetic.
`StabilityRegion::from(map)` constructs the boxed Möbius variant directly.

Every nonzero report contains the transformed polynomial and degree drop,
the normalized integer parts and gcd, the successive Sturm layers, the
reduced Bézoutian, and its exact block LDLᵀ certificate. The decomposition
lives in `exact_algebra::inertia`: its permutation records original
zero-based indices, and it satisfies A[permutation, permutation] = L D Lᵀ.
The diagonal blocks are 1×1 scalars or 2×2 blocks with zero diagonal and
nonzero off-diagonal entries. Such a 2×2 block has one eigenvalue of each sign.

`SymmetricLdlt::verify` checks shape, unit lower triangularity, permutation
validity and the factorization identity. `StabilityReport::verify` performs
deterministic exact replay from the supplied polynomial and compares every
stored field, so altered counts or proof data are rejected. It has the same
cost as constructing the report; this API does not claim a cheaper replay.

## Cost and scope

Homogeneous substitution and the coefficient recurrence for the Bézoutian
each take O(n²) arithmetic operations. Dense rational LDLᵀ takes O(n³)
operations and O(n²) storage. Bit complexity depends on coefficient growth;
these are arithmetic-operation bounds, not constant-time integer operations.
The repeated-root filtration performs at most n Sturm/squarefree passes on
strictly decreasing degrees. Certificates retain all filtration layers and
the dense matrix/factorization, using O(n²) exact entries overall.

This API covers univariate polynomials only. Command-line and MCP adapters,
their input parsers, size caps and JSON schemas are separate interfaces.
