# polytool

Dense univariate polynomial toolkit for combinatorial research. Provides
real-rootedness testing by primitive integer root counting, Bézout-matrix
interlacing checks, log-concavity, gamma-positivity, resultants,
Ehrhart/h\*-vector conversion, recurrence search for polynomial sequences,
finite Athanasiadis--Wagner interlacing matrices, and standard sequences —
all with exact arithmetic.

Rename note: this crate was renamed from `polynomial-tools` to `polytool`
before publication. Use `polytool` in Cargo manifests, Rust imports, CLI
documentation, and MCP server configuration.

## Installation

From the current workspace root:

```sh
cargo build --release -p polytool
```

The CLI binary is at `target/release/polytool`.

Use `polytool --version` (or `polytool -V`) for a single deterministic line:

```text
polytool 0.2.1-rc.5 (git 0123456789ab)
```

The version comes from the crate manifest and the lowercase 12-hex-digit Git
commit is captured at build time.  If the source is built without trustworthy
Git metadata, the suffix is `(git unavailable)` rather than a fabricated
revision.  Reproducible packaging environments may set `POLYTOOL_GIT_COMMIT`
to a full or at least 12-digit hexadecimal revision; an empty or malformed
value deliberately selects the unavailable form.

On the generated standalone `main` branch, the equivalent command is:

```sh
cargo build --release
```

## Using as a crate

The library crate is usable today from another Rust project by depending on
this directory:

```toml
[dependencies]
polytool = { path = "/home/paxinum/Dropbox/AI-projects/rust/polytool" }
```

From the monorepo `master` branch, Cargo can also depend on the package by
name:

```toml
[dependencies]
polytool = { git = "https://github.com/PerAlexandersson/polytool.git", branch = "master", package = "polytool" }
```

The repository's generated standalone `main` branch supports the shorter form:

```toml
[dependencies]
polytool = { git = "https://github.com/PerAlexandersson/polytool.git" }
```

After publishing to crates.io, the dependency would be:

```toml
[dependencies]
polytool = "0.1"
```

The public API is re-exported from `polytool`, so a downstream crate
can write:

```rust
use polytool::{check_weak_interlacing, is_real_rooted};

assert!(is_real_rooted(&[1, 11, 11, 1]));
assert_eq!(check_weak_interlacing(&[2, -3, 1], &[-1, 1]), Some(true));
```

Coefficient vectors are always in ascending degree order. Most convenience
functions accept `&[i64]`; use the `*_bigint_coeffs` variants when coefficients
may exceed `i64`. Interlacing functions return `Option<bool>`: `None` means the
directed degree relation is not valid for that test, not that interlacing failed.

For Hermite--Biehler or Euclidean-chain experiments, the library also exposes a
signed Sturm continued-fraction certificate checker over exact rationals:

```rust
use num_bigint::BigInt;
use num_rational::BigRational;
use polytool::check_sturm_continued_fraction_certificate;

let q = |n: i64| BigRational::from_integer(BigInt::from(n));
let p0 = vec![q(1), q(3), q(1)]; // x^2 + 3x + 1
let p1 = vec![q(2), q(1)];       // x + 2
let certificate =
    check_sturm_continued_fraction_certificate(&p0, &p1).expect("Sturm-CF certificate");
assert_eq!(certificate.steps, 2);
```

To build the MCP server from the workspace root:

```sh
cargo build --release -p polytool-mcp
```

The MCP binary is at `target/release/polytool-mcp`.

To install the MCP server into a user bin directory and print a client
configuration snippet:

```sh
./polytool/mcp/install.sh
```

On the generated standalone `main` branch, run:

```sh
./mcp/install.sh
```

## CLI usage

### Input format

Polynomials are given as comma-separated integer coefficients in
**ascending degree order**: `a_0, a_1, ..., a_d` represents
`a_0 + a_1 t + ... + a_d t^d`.

A text file `polys.txt` might look like:

```
1, 11, 11, 1
1, 26, 66, 26, 1
1, 57, 302, 302, 57, 1
```

The commands `interlacing`, `interlacing-profile`, `properties`,
`gamma-expansion`, `family-check`, `sequence`, `coefficient-tests`,
`hstar-to-ehrhart`, `ehrhart-to-hstar`, `hstar-inequalities`,
`cyclic-sieving`, and `cyclic-sieving-sequence` also accept `--json` for
machine-readable output.
Large integer coefficients in JSON output are serialized as strings.

The exact property commands `real-rooted`, `interlacing`,
`interlacing-profile`, `properties`, `gamma-expansion`, and the non-recurrence
checks inside `family-check` accept arbitrary-size integer coefficients.
Standard sequence generation, resultants, discriminants, Ehrhart conversion,
Stapledon decomposition, and the standalone `recurrence` command also use
arbitrary-size exact integers in the CLI. The main remaining compatibility
exception is `family-check --recurrence`, which currently requires
coefficients that fit in `i64`.

### Check real-rootedness

```sh
polytool real-rooted < polys.txt
```

`real-rooted` accepts arbitrary-size integer coefficients in dense-list and
expanded-polynomial input.

Output (one line per polynomial):

```
1 + 11t + 11t^2 + t^3: real-rooted
1 + 26t + 66t^2 + 26t^3 + t^4: real-rooted
1 + 57t + 302t^2 + 302t^3 + 57t^4 + t^5: real-rooted
```

### Check interlacing

Test whether consecutive polynomials interlace:

```sh
polytool interlacing < polys.txt
polytool interlacing --json < polys.txt
```

Count how many consecutive previous rows each polynomial interlaces, stopping
the backward scan at the first failure:

```sh
polytool interlacing-profile < polys.txt
polytool interlacing-profile --json < polys.txt
```

The JSON output includes `previous_count`, `checked_previous_count`, and
`interlacing_previous_count`, plus the checked pair reports.
Both commands accept arbitrary-size integer coefficients in dense-list and
expanded-polynomial input.

### Check unimodality, log-concavity, palindromicity, gamma-positivity

```sh
polytool properties < polys.txt
polytool properties --json < polys.txt
```

`properties` accepts arbitrary-size integer coefficients in both dense-list
and expanded-polynomial input.  When all coefficients fit in `i64`, it uses the
same fast exact path as the rest of the CLI; otherwise it falls back to exact
BigInt property checks.

Output includes all properties for each polynomial:

```
1 + 11t + 11t^2 + t^3: real-rooted, palindromic, gamma-positive [1, 8], unimodal, log-concave, ultra-log-concave
```

For a palindromic polynomial, print the full gamma expansion:

```sh
polytool gamma-expansion < polys.txt
polytool gamma-expansion --json < polys.txt
```

`gamma-expansion` accepts arbitrary-size integer coefficients.  In JSON output,
both input coefficients and gamma coefficients are serialized as strings.

Example output:

```text
1 + 11t + 11t^2 + t^3: gamma [1, 8]; expansion: (1+t)^3 + 8 t (1+t)
```

### Coefficient-only real-rootedness tests

```sh
polytool coefficient-tests < polys.txt
polytool coefficient-tests --json < polys.txt
```

This reports Newton inequalities and Kurtz's sufficient condition.  When the
Kurtz condition holds, the polynomial has distinct real roots.

### Generate standard sequences

```sh
polytool sequence eulerian 5
polytool sequence narayana 5 --json
```

Supported sequence names are `eulerian`, `narayana`, `type-b-eulerian`,
`chebyshev-t`, `chebyshev-u`, and `hermite`.
Generated coefficients use arbitrary-size integers.

### Generate recurrence-backed OEIS families

The bundled OEIS catalog exposes one Rust function per A-number, for example
`polytool::oeis::A008292()`, together with dynamic lookup by identifier.  Each
entry is generated from a sparse exact recurrence description and its initial
rows; installed binaries do not need an OEIS database, Python, or SymPy.

```sh
polytool oeis list
polytool oeis info A008292 --json
polytool oeis generate A008292 --rows 8
polytool oeis generate A008292 --rows 8 --format triangle
polytool oeis generate A008292 --rows 40 --format bfile --max-terms 1000
```

Generation formats are `rows`, `triangle`, `polynomial`, `json`, `jsonl`,
`csv`, and `bfile`.  `--start-row` selects a later displayed OEIS row.  Strict
b-file output is enabled only when the bundled complete-row prefix has been
matched against OEIS data; `--max-terms` stops before a row that would cross
the cap.  The catalog distinguishes holdout-`verified` recurrences from
`validated` recurrences whose generated rows match the current OEIS prefix but
whose original fitting/holdout provenance is unavailable.  Both are enabled by
default.  Entries without a safe prefix or row-layout match are `experimental`,
hidden by default, and require `--include-experimental`.

The source generator is `scripts/build_oeis_catalog.py`.  It imports the
machine-readable recurrence benchmark fixtures, supplements them from the
curated `OEIS-polynomials/sequences` queue, and imports the 728 canonical
generated definitions in `real-rooted-oeis-proofs`.  The latter are converted
from their restricted Lean expression grammar into the same sparse exact
representation and replayed against local OEIS data.  The generator rejects
queue recurrences that do not reproduce their cached rows and emits
`src/oeis_catalog_generated.rs`.  Run it with `--check` in verification jobs
to detect drift.

### Check Family H PF/Jensen pencils

```sh
polytool pf-pencil --case A036969 --degree 5
polytool pf-pencil --case A036969 --min-degree 1 --max-degree 10 --json
polytool pf-pencil --all-family-h --max-degree 8 --json
```

`pf-pencil` builds the fixed Family H coefficient actions used in the
`real-rooted-oeis` proof work.  For an old-row degree `d`, it forms the
integer-cleared Jensen endpoint kernels

```text
J_alpha,d(t) = sum_k binomial(d,k) alpha(k,d) t^k
t J_beta,d(t) = t sum_k binomial(d,k) beta(k,d) t^k
```

and checks exact real-rootedness of both endpoints and of the finite positive
pencil specializations `J_alpha,d + lambda * t J_beta,d`.  The default
lambda list is `0,1,10,100`; override it with `--lambdas 0,2,50`.

Supported built-in cases are:

```text
A036969 A071951 A080248 A156289 A160562 A269945
A166960 A166961 A166962 A166972 A191935
A371081 A371259 A390433 A198204
```

JSON output uses schema `polytool.pf-pencil.v1`.  Coefficients and
lambda values are serialized as strings.  Each item reports the case id,
degree, alpha/beta labels, common denominator, endpoint coefficient vectors,
endpoint real-rootedness booleans, lambda checks, and `finite_evidence_ok`.
Common-factor/residual extraction is not yet implemented; the full endpoint
polynomials are reported instead.

### Check a family in one pass

`family-check` reports properties, consecutive weak/strict interlacing, and
optionally an adaptive recurrence search.  Requirement flags make the command
exit nonzero at the first failed requested condition.
Property checks and consecutive interlacing accept arbitrary-size integer
coefficients.  Recurrence search currently requires coefficients that fit in
`i64`; oversized rows are preserved in the report and recurrence is marked
unavailable rather than silently dropped.

```sh
polytool family-check \
  --require-real-rooted \
  --require-weak-interlacing \
  --recurrence \
  --json < polys.txt
```

### Search for a recurrence

Given a sequence of polynomials (one per line), search for a polynomial
recurrence `f(n,t) P_n(t) = Σ c_{r,d}(n,t) D^d P_{n-r}(t)`. With
`--alternating-sign`, the search also allows terms
`(-1)^n c_{r,d}(n,t) D^d P_{n-r}(t)`:

```sh
polytool recurrence < polys.txt
```

The adaptive search orders candidates by a weighted parameter count. Ordinary
coefficient unknowns are cheapest; derivative, denominator, alternating-sign,
and inhomogeneous unknowns are delayed. A candidate is only solved when the
available fitting prefix has at least `unknowns + min_margin` scalar equations.
By default the last input row is reserved for exact verification, the fit uses
the first solvable prefix plus one extra row, and modular prefiltering rejects
impossible candidates before exact rational solving. Use `--no-verify` to fit
against all input rows.

Options:

```
--skip-prefix <k>    Ignore the first k input polynomials before searching
--min-rec-len <k>    Minimum recurrence depth to try (default: 1)
--max-rec-len <k>    Maximum recurrence depth (default: 5)
--min-var-deg <d>    Minimum degree in t for coefficients (default: 0)
--max-var-deg <d>    Maximum degree in t for coefficients (default: 3)
--min-idx-deg <d>    Minimum degree in n for coefficients (default: 0)
--max-idx-deg <d>    Maximum degree in n for coefficients (default: 3)
--min-diff-deg <d>   Minimum derivative order (default: 0)
--max-diff-deg <d>   Maximum derivative order (default: 2)
--inhomogeneous      Also try inhomogeneous recurrences
--min-inhomo-var-deg Minimum degree in t for the inhomogeneous term
--max-inhomo-var-deg Maximum degree in t for the inhomogeneous term
--min-inhomo-idx-deg Minimum degree in n for the inhomogeneous term
--max-inhomo-idx-deg Maximum degree in n for the inhomogeneous term
--denominator        Allow a nontrivial LHS factor f(n,t)
--alternating-sign   Also allow right-hand-side terms multiplied by (-1)^n
--max-denom-var-deg  Max degree in t for f(n,t) (default: 2, implies --denominator)
--max-denom-idx-deg  Max degree in n for f(n,t) (default: 2, implies --denominator)
--min-margin <k>     Require equations >= unknowns + k (default: 1)
--fit-extra-rows <k> Extra rows after the first solvable prefix (default: 1)
--no-verify          Fit all rows instead of reserving held-out verification rows
--no-modular-prefilter
                      Disable default modular candidate rejection
--max-candidates <n> Inspect at most n adaptive candidate configurations
--json               Emit recurrence JSON with initial conditions
--python             Emit exact standalone Python code using Fraction arithmetic
--format json        Alias for --json
--format python      Alias for --python
--verbose            Print each candidate tried
```

The modular prefilter is enabled by default. It is often much faster on false
candidates because it rejects a candidate when every usable fixed large-prime
reduction is inconsistent. Use `--no-modular-prefilter` only when comparing
against the exact-only search path.

`--max-candidates N` is a deterministic iteration budget, not a wall-clock
timeout. A candidate is counted as soon as its parameter configuration is
taken from the adaptive iterator, before fit-row, structural, modular, or
exact-solve filtering. Thus `0` evaluates no candidates, candidate `N` may
still succeed, and budget exhaustion is reported only when an `(N+1)`st
candidate exists. In text mode exhaustion is a distinct diagnostic and exit
status `3`; with `--format json`, the JSON status is `budget_exhausted` and
includes `candidates_considered`, `candidates_tried`, and `max_candidates`.
If all configured candidates have been checked, the JSON status is instead
`not_found`. Omitting the option preserves the unbounded behavior.

#### Machine-readable recurrence JSON and row generation

To save a recurrence with enough initial conditions to regenerate the sequence,
use JSON output:

```sh
polytool recurrence --json < polys.txt > recurrence.json
```

The JSON schema stores exact rational coefficients as strings, the recurrence
terms, the LHS denominator if present, and the minimal initial polynomial rows.
Generate more rows from the saved recurrence with:

```sh
polytool recurrence-generate --recurrence recurrence.json --rows 100 > rows-100.txt
```

Use `--additional n` instead of `--rows n` to keep all initial rows and append
`n` newly generated rows.

The web interface exports this compact JSON record and links back to these
generation instructions instead of emitting a large standalone Python program.

Use Python output when a sequence extension script is more convenient than JSON:

```sh
polytool recurrence --python < polys.txt > recurrence.py
```

The generated Python uses exact `fractions.Fraction` arithmetic and dense
coefficient lists in ascending powers of `t`.

### Scout BKW equal-modulus loci

For a polynomial recurrence with characteristic symbol

```text
F(x,z) = a_0(x) + a_1(x) z + ... + a_r(x) z^r,
```

`bkw-scout` numerically scans a complex rectangle for points where the two
dominant characteristic roots have nearly equal modulus.  This is a scout for
Beraha--Kahane--Weiss accumulation obstructions; it does not certify dominance,
amplitudes, or eventual non-real-rootedness.

Give the symbol as the z-coefficient polynomials in ascending z-degree:

```sh
polytool bkw-scout \
  --symbol '1; -x; 1' \
  --box -3 3 -2 2 \
  --grid 61 \
  --top 10
```

The same symbol can be supplied on stdin, one coefficient polynomial per line.
Use `--format json` for machine-readable output, and `--mathematica` to print a
`Reduce` skeleton for exact equal-modulus follow-up.

### Compute resultant/discriminant

```sh
# Resultant of two polynomials (given as two lines)
echo "1, 0, 1
1, -1" | polytool resultant

# Discriminant of a single polynomial
echo "1, 0, -3, 1" | polytool discriminant
```

Inputs for these CLI commands may have arbitrary-size integer coefficients.
The exact resultant or discriminant output is printed as an integer.

### Ehrhart ↔ h\*-vector conversion

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

### Cyclic sieving checks

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

### Stapledon decomposition

```sh
echo "1, 11, 11, 1" | polytool stapledon 3
```

CLI Stapledon input coefficients may be arbitrary-size integers.

### Benchmark recurrence and interlacing performance

Run the recurrence fixture timing suite:

```sh
polytool bench recurrence-fixtures --repeat 3
polytool bench recurrence-fixtures --only 23_sparse --repeat 5
polytool bench recurrence-fixtures --only oeis --repeat 1
polytool bench recurrence-fixtures --only oeis --repeat 3 \
  --summary --report bench-results/recurrence-fixtures/oeis.md
polytool bench recurrence-fixtures --only oeis --repeat 3 --format json \
  > bench-results/recurrence-fixtures/oeis.json
polytool bench compare old.json new.json --top 10
```

Run consecutive interlacing timings for a generated sequence:

```sh
polytool bench interlacing --sequence eulerian --max-n 30 --repeat 5
```

Both benchmark subcommands print tab-separated output so results can be
redirected to a log or pasted into project notes. The recurrence fixture
benchmark can also append fixture/category summaries with `--summary` and write
a Markdown report with `--report <path.md>`. Use `--format json` for
machine-readable per-run records, fixture/category summaries, and adaptive
search diagnostics. The diagnostic records include rejection counters, modular
lift counters, and cumulative timing buckets such as derivative precomputation,
modular prefiltering, modular lifting, exact solving, and held-out
verification. `polytool bench compare` compares two such JSON runs and
highlights fixture/category speedups or regressions. The recurrence fixture
suite contains synthetic stress tests and natural OEIS-derived recurrences from
the curated `real-rooted-oeis` sequence queue.

## Library usage

Add to your `Cargo.toml`:

```toml
[dependencies]
polytool = { path = "../polytool" }
```

### Real-rootedness and interlacing

```rust
use polytool::*;

// Coefficients in ascending degree order: coeffs[i] = coeff of t^i
let eulerian_4 = [1, 11, 11, 1];

// Default exact PRS/root-counting method
assert!(is_real_rooted(&eulerian_4));

// Sturm chain method (slower, gives root locations)
assert!(is_real_rooted_sturm(&eulerian_4));

// Strict interlacing: deg(p) = deg(q) + 1
let p = [-15, 23, -9, 1];  // (t-1)(t-3)(t-5)
let q = [8, -6, 1];         // (t-2)(t-4)
assert_eq!(check_interlacing(&p, &q), Some(true));

// Weak interlacing (allows shared roots)
let f = [2, -3, 1];   // (t-1)(t-2)
let g = [-1, 1];       // (t-1)
assert_eq!(check_weak_interlacing(&f, &g), Some(true));

// Same-degree interlacing: roots alternate on the real line
// (t-1)(t-3) roots {1,3} and (t-2)(t-4) roots {2,4}: 1 < 2 < 3 < 4
assert_eq!(check_weak_interlacing(&[3, -4, 1], &[8, -6, 1]), Some(true));

// Nested roots do NOT interlace: (t-1)(t-4) vs (t-2)(t-3)
assert_eq!(check_weak_interlacing(&[4, -5, 1], &[6, -5, 1]), Some(false));
```

For large exact coefficients, use the `BigInt` APIs:

```rust
use num_bigint::BigInt;
use polytool::check_interlacing_bigint_coeffs;

let center = BigInt::from(10).pow(20);
let f = vec![-&center, BigInt::from(1)]; // t - center
let g = vec![
    (&center - 1u32) * (&center + 1u32),
    -BigInt::from(2) * &center,
    BigInt::from(1),
];

assert_eq!(check_interlacing_bigint_coeffs(&f, &g), Some(true));
```

### Athanasiadis--Wagner interlacing matrices

```rust
use polytool::*;

// A column vector of polynomials P, Q, R in ascending coefficient order.
let polys = vec![
    vec![2, 1],       // 2 + t
    vec![8, 6, 1],    // (2+t)(4+t)
    vec![3, 4, 1],    // (1+t)(3+t)
];

// One finite truncation of the infinite Lace(A) matrix.
let lace = lace_matrix_sequence_i64(&polys, 1, 3).unwrap();
assert_eq!(lace, vec![vec![2, 1, 0], vec![8, 6, 1], vec![3, 4, 1]]);

// This is Athanasiadis--Wagner's pairwise-but-not-fully-interlacing example:
// the 3 x 3 determinant is negative.
assert!(check_lace_sequence_total_nonnegative_i64(&polys, 1, 3, 3).is_err());
```

For a `p x q` polynomial matrix `A`, `lace_matrix_i64(A, r, c)` returns the
`p r` by `q c` truncation whose `(p rb+i, q cb+j)` entry is the coefficient
of `x^(cb-rb)` in `A_ij(x)`.  This follows Definition 3.5 of
Athanasiadis--Wagner, *Veronese sections and interlacing matrices of formal
power series*.  A finite TNN check is computational evidence for the infinite
matrix, not a proof of full interlacing unless a separate finite criterion
applies.

### Polynomial arithmetic

```rust
use polytool::Polynomial;

let p = Polynomial::<i64>::from_i64_coeffs(&[1, 1]);  // 1 + t
let q = p.clone() * p.clone();                         // 1 + 2t + t^2
assert_eq!(q.evaluate(3), 16);
assert!(q.is_palindromic());

let deriv = q.derivative();  // 2 + 2t
let gcd = Polynomial::<i64>::gcd(&q, &deriv);
```

### Recurrence search

```rust
use polytool::recurrence::*;

// Eulerian polynomials A_1, A_2, ..., A_10
let polys: Vec<Vec<i64>> = vec![
    vec![1],
    vec![1, 1],
    vec![1, 4, 1],
    vec![1, 11, 11, 1],
    vec![1, 26, 66, 26, 1],
    vec![1, 57, 302, 302, 57, 1],
    vec![1, 120, 1191, 2416, 1191, 120, 1],
    vec![1, 247, 4293, 15619, 15619, 4293, 247, 1],
    // ...
];

// Adaptive search: tries small parameter spaces first
let result = find_recurrence_adaptive(&polys, &AdaptiveSearchOptions::default());
if let Some(res) = result {
    println!("{}", res.recurrence);
}

// Deterministic candidate budget with an exact termination reason.
let outcome = find_recurrence_adaptive_with_budget(
    &polys,
    &AdaptiveSearchOptions::default(),
    AdaptiveSearchBudget::limited(100),
);
match outcome {
    AdaptiveSearchOutcome::Found(result) => println!("{}", result.recurrence),
    AdaptiveSearchOutcome::NoRecurrence(_) => println!("search space exhausted"),
    AdaptiveSearchOutcome::BudgetExhausted(summary) => println!(
        "budget exhausted after {} candidates",
        summary.diagnostics.considered_candidates
    ),
}

// Or search with specific parameters
let opts = RecurrenceOptions {
    rec_len: 2,
    var_deg: 1,
    idx_deg: 1,
    diff_deg: 1,
    ..Default::default()
}.with_alternating_sign(true);
if let Some(rec) = find_polynomial_recurrence(&polys, &opts) {
    println!("{}", rec);
}
```

### Coupled Weyl-matrix recurrences

The library also fits exact first-order systems

```text
F_(n+1) = M(n,x,D_x) F_n + G(n,x),
```

where `F_n` is a vector of polynomials. Every matrix entry is stored in the
normal form `sum_d c_d(n,x) D_x^d`; this is general for the Weyl algebra because
`D_x x = x D_x + 1`. Input has shape
`states[index][component][x_degree]`, and `first_index` explicitly fixes the
meaning of `n`:

```rust
use polytool::recurrence::{find_vector_recurrence, VectorRecurrenceOptions};

// F_n = [n]_x, including [0]_x = 0.  The affine recurrence is
// F_(n+1) = x F_n + 1.
let states = (0..=9)
    .map(|n| vec![if n == 0 { vec![0] } else { vec![1; n] }])
    .collect::<Vec<_>>();
let options = VectorRecurrenceOptions {
    var_deg: 1,
    idx_deg: 0,
    homogeneous: false,
    forcing_var_deg: 0,
    forcing_idx_deg: 0,
    held_out_transitions: 2,
    ..Default::default()
};
let fit = find_vector_recurrence(&states, 0, &options)?;
assert_eq!(fit.rows[0].rank, fit.rows[0].unknowns);
# Ok::<(), Box<dyn std::error::Error>>(())
```

Rows of `M` are solved independently over exact rationals. The result reports
the exact rank and nullity of each row; `require_unique` defaults to `true`, so
cross-component dependencies are rejected instead of being hidden by an
arbitrary choice of free parameters. Final complete transitions are held out
of the fit and then verified exactly.

For a lag-`r` system, `find_companion_vector_recurrence` returns the first-order
recurrence on `(F_(n-r+1), ..., F_n)`. It fits only the unknown final block row
and inserts the shift-identity rows directly. Unit tests recover two published
systems: the up-down-run Eulerian pair of
[Ma--Ma--Yeh--Yeh (2022)](https://doi.org/10.1016/j.disc.2021.112716), and the
lag-two type-B `1/k`-Eulerian system of
[Ma et al. (2020)](https://doi.org/10.37236/9089), specialized to `k=1`.

### Other properties

```rust
use polytool::*;

let coeffs = [1, 11, 11, 1];

assert!(is_palindromic(&coeffs));
assert!(is_unimodal(&coeffs));
assert!(is_log_concave(&coeffs));
assert!(is_ultra_log_concave(&coeffs));
assert!(is_gamma_positive(&coeffs));

// Gamma coefficients: p(t) = Σ γ_i t^i (1+t)^{d-2i}
let gamma = gamma_coefficients(&coeffs).unwrap();
assert_eq!(gamma, vec![1, 8]);

// Resultant and discriminant
let disc = discriminant(&[1, 0, -3, 1]);
```

## MCP server

The repository also contains a local Model Context Protocol server in
`mcp/`. It exposes the exact polynomial routines to MCP clients over stdio.

Build it from the workspace root:

```sh
cargo build --release -p polytool-mcp
```

Install it locally:

```sh
./polytool/mcp/install.sh
```

Check the installed binary:

```sh
polytool-mcp --help
polytool-mcp --version
```

Example MCP client configuration:

```json
{
  "mcpServers": {
    "polytool": {
      "command": "/absolute/path/to/rust/target/release/polytool-mcp"
    }
  }
}
```

The MCP server is tools-only: no resources, prompts, HTTP server, sampling, or
filesystem access. It returns structured JSON and also includes the same compact
JSON as text content for client compatibility.

Available MCP tools:

- `parse_polynomials`
- `check_polynomial_family`
- `polynomial_properties`
- `check_interlacing_pair`
- `check_interlacing_sequence`
- `check_interlacing_profile`
- `real_roots`
- `find_recurrence`
- `generate_recurrence_rows`
- `resultant`
- `discriminant`
- `ehrhart_hstar`
- `analyze_decomposition`
- `generate_sequence`

The MCP property, interlacing, resultant, discriminant, Ehrhart conversion,
sequence-generation, and recurrence tools accept arbitrary-size integer
coefficients. JSON integers may be provided directly when they fit the client
stack, and larger exact integers should be sent as strings.

`find_recurrence` exposes all adaptive-search controls as top-level MCP
arguments. For example, a compact exact search can use:

```json
{
  "coefficients": [[1], [2], [4], [8], [16]],
  "max_rec_len": 1,
  "max_var_deg": 0,
  "max_idx_deg": 0,
  "max_diff_deg": 0,
  "include_code": false
}
```

The legacy nested `options` object remains accepted, with top-level values
taking precedence. Omitting `include_code` preserves the full generated-code
response.

See [`mcp/README.md`](mcp/README.md) for request schemas, examples, and
development notes.

## Development

Run the core tests:

```sh
cargo test -p polytool
```

Run the MCP server tests:

```sh
cargo test -p polytool-mcp
```

Run the web wrapper tests:

```sh
cargo test -p polytool-web
node web/tests/string_safety.mjs
```

### Web arbitrary-precision contract and deployment staging

The browser wrapper parses polynomial input with `parse_polynomials_bigint`.
Every integer coefficient returned across its string-based WASM/JSON boundary
is a canonical decimal string, including input, gamma, decomposition, and
recurrence-initial coefficients. Resultants and discriminants are decimal
strings as well. Browser code must keep these values as strings or convert
them to JavaScript `BigInt` for exact arithmetic; converting them to `Number`
can round values above `2^53`.

Build the deployable no-modules bundle from the monorepo root with an external
target directory:

```sh
CARGO_TARGET_DIR=/cargo-target/ai-projects timeout 60s nice -n 10 \
  wasm-pack build polytool/web --target no-modules --release \
  --out-dir pkg
```

After review and explicit deployment authorization, stage exactly the files
referenced by `web/index.html` into the website checkout. The backup option
keeps replaced files recoverable, and these commands deliberately do not use
the website's historical `--delete` deployment target:

```sh
site_stage=/home/paxinum/Dropbox/webpages/poly.symmetricfunctions.com/www
install -d "$site_stage/pkg"
cp --backup=numbered web/index.html "$site_stage/index.html"
cp --backup=numbered web/favicon.svg "$site_stage/favicon.svg"
cp --backup=numbered web/pkg/polytool_web.js "$site_stage/pkg/polytool_web.js"
cp --backup=numbered web/pkg/polytool_web_bg.wasm \
  "$site_stage/pkg/polytool_web_bg.wasm"
cmp web/index.html "$site_stage/index.html"
cmp web/pkg/polytool_web.js "$site_stage/pkg/polytool_web.js"
cmp web/pkg/polytool_web_bg.wasm "$site_stage/pkg/polytool_web_bg.wasm"
```

The legacy website Makefile still names `polynomial_tools_web*`; do not use
that stale assembly rule for the current `polytool_web*` bundle without first
updating and reviewing the website repository itself.

The repository's `main` branch is generated from this directory with
`git subtree split`; the monorepo `master` branch remains canonical. The core
crate, `mcp/`, `web/`, and `Cargo.lock` are therefore published together.

## Support

Report issues in the Git repository, or contact Per Alexandersson
(@PerAlexandersson, <per.w.alexandersson@gmail.com>).

## Algorithm notes

### Real-rootedness algorithms

The default real-rootedness check uses an adaptive exact path from
`root_count`. It replaces the polynomial by its square-free part and counts
distinct real roots exactly over `BigInt`; the polynomial is real-rooted
precisely when this count equals the square-free degree. For one-signed
coefficient polynomials, the backend first tries cheap coefficient filters and
then counts positive roots of `f(-t)`, since such polynomials can only have
non-positive real roots. Primitive PRS remains the general backend. Uspensky is
selected conservatively at degree 35 or above when the endpoint coefficients
are equal and an interior coefficient is at least `4^degree` times as large as
an endpoint. This benchmark-derived, scale-invariant rule captures the large
Eulerian, type-B Eulerian, and Touchard cases tested without penalizing
Narayana polynomials or evenly spaced linear-factor products.

For algorithm comparisons, `is_real_rooted_uspensky_bigint_coeffs` uses an
independent exact Uspensky/Descartes path.  It applies a strict Fujiwara root
bound, separates roots into dyadic magnitude bands (using the reciprocal
polynomial below `1`), and applies homographic subdivision with Descartes sign
variations.  It does not use finite fields, floating point, or approximate
roots. The explicit PRS and Uspensky entry points remain available for direct
comparison.

Explicit matrix comparison paths remain available as
`is_real_rooted_bezout`, `is_real_rooted_bezout_squarefree`, and
`is_real_rooted_hermite`.

All interlacing checks use the **Bézout matrix** (Fisk, Cor. 9.145). For
polynomials f (degree d) and g (degree d−1), the Bézout matrix B(f,g) is the
d×d symmetric matrix with (i,j) entry equal to the coefficient of x^i y^j in
`(f(x)g(y) - f(y)g(x)) / (x-y)`.

**Theorem:** f and g are both real-rooted and g strictly interlaces f
if and only if B(f,g) is positive definite.

This reduces interlacing to a single exact matrix definiteness check, avoiding
root isolation entirely. It is 100–400× faster than Sturm chains at degree 15+
for interlacing checks.  The definiteness implementation uses fraction-free
BigInt Bareiss elimination below dimension `30`, and switches to a
CRT-over-prime-fields path for larger matrices.  The modular path reconstructs
the leading principal minors exactly using Hadamard bounds.  For speed
comparisons, the explicit entry points are `is_positive_definite_bareiss`,
`is_positive_definite_modular`, and `modular_leading_principal_minors_bigint`.

Run the local comparison example with:

```bash
cargo run --release -p polytool --example bench_bareiss_vs_modular_linalg
```

For **same-degree** polynomials, `check_weak_interlacing` reduces to the
deg+1 case by extending one polynomial with a root far to the right
(using the Cauchy bound for root radius). When all coefficients are
positive (all roots negative), multiplying by t suffices, skipping the
Cauchy bound computation.

See `doc/bezout-interlacing.md` for the Mathematica reference implementation.

### Interlacing matrices

Athanasiadis--Wagner's interlacing matrix construction packages a polynomial
sequence or polynomial matrix into an infinite block Toeplitz matrix `Lace(A)`.
For a column vector, total nonnegativity of this infinite matrix is called
full interlacing and implies pairwise interlacing of the entries.  The converse
fails, so this gives a stronger finite-experiment target than checking
ordinary pairwise interlacing.

The crate implements finite truncations through `lace_matrix_i64`,
`lace_matrix_bigint`, and the sequence wrappers
`lace_matrix_sequence_i64` / `lace_matrix_sequence_bigint`.  Exact finite TNN
checks are available with `check_lace_total_nonnegative_i64` and the faster
Neville-elimination wrappers.
