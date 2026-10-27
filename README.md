# polytool

Dense univariate polynomial toolkit for combinatorial research. Provides
real-rootedness testing by adaptive subresultant integer root counting,
Bézout-matrix interlacing checks, exact Hurwitz stability, log-concavity,
gamma-positivity, resultants, Ehrhart/h\*-vector conversion, recurrence search
for polynomial sequences, finite Athanasiadis--Wagner interlacing matrices,
and standard sequences — all with exact arithmetic.

Rename note: this crate was renamed from `polynomial-tools` to `polytool`
before publication. Use `polytool` in Cargo manifests, Rust imports, CLI
documentation, and MCP server configuration.

## Features

- Exact real-rootedness testing by adaptive subresultant integer root
  counting (the engine is the standalone `exact-real-roots` crate),
  Bézout-matrix and Sturm-chain interlacing tests, and exact Hurwitz stability.
- Unimodality, log-concavity, ultra-log-concavity, palindromicity and
  gamma-positivity (with the gamma expansion), Bernstein-basis expansion,
  resultants and discriminants.
- Ehrhart/h\*-vector conversion, h\*-inequalities, cyclic and bicyclic sieving
  checks, Stapledon decomposition.
- Recurrence search for polynomial sequences (with certified complete spaces,
  coordinate transformations and implicit differential relations), renewal
  recognition, and a bundled catalog of recurrence-backed OEIS families.
- Planar-network certificates: TNN triangles, strip digraphs and Pólya
  frequency rows.
- Finite Athanasiadis--Wagner interlacing matrices, standard sequences
  (Eulerian, Narayana, Chebyshev, Hermite, ...), and enumeration of real-rooted
  integer polynomials with a prescribed leading prefix.
- A command-line tool, an MCP server and a browser (WebAssembly) wrapper, all
  over the same exact-arithmetic library.

## Installation

From the repository root:

```sh
cargo build --release -p polytool
```

The CLI binary is at `target/release/polytool`.

Use `polytool --version` (or `polytool -V`) for a single deterministic line:

```text
polytool 0.6.0 (git 0123456789ab)
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
polytool = { path = "../polytool" }
```

From GitHub, depend on the repository's standalone `main` branch:

```toml
[dependencies]
polytool = { git = "https://github.com/PerAlexandersson/polytool.git" }
```

After publishing to crates.io, the dependency would be:

```toml
[dependencies]
polytool = "0.6"
```

The public API is re-exported from `polytool`, so a downstream crate
can write:

```rust
use polytool::{check_weak_interlacing, is_real_rooted};

assert!(is_real_rooted(&[1, 11, 11, 1]));
assert_eq!(check_weak_interlacing(&[-1, 1], &[2, -3, 1]), Some(true));
```

Coefficient vectors are always in ascending degree order. Most convenience
functions accept `&[i64]`; use the `*_bigint_coeffs` variants when coefficients
may exceed `i64`. Interlacing functions return `Option<bool>`: `None` means the
directed degree relation is not valid for that test, not that interlacing failed.
For the undirected question "do `f` and `g` interlace in some order?", use
`interlaces_either_order`, which returns a plain `bool`.

More library examples (real-rootedness and interlacing, Athanasiadis--Wagner
interlacing matrices, planar networks, polynomial arithmetic, recurrence
search, enumeration of real-rooted integer polynomials, ...) are in
[doc/library.md](doc/library.md).  API documentation is on
[docs.rs](https://docs.rs/polytool).

## CLI quick start

Polynomials are read from standard input, one per line, as comma-separated
integer coefficients in ascending degree order (`1, 11, 11, 1` is
`1 + 11t + 11t^2 + t^3`):

```sh
echo "1, 11, 11, 1" | polytool real-rooted
polytool interlacing --json < polys.txt
polytool properties < polys.txt
polytool sequence eulerian 5
polytool oeis generate A008292 --rows 8 --format triangle
polytool recurrence < polys.txt
echo "1, 0, -3, 1" | polytool discriminant
```

Many commands accept `--json` for machine-readable output, and large integer
coefficients are serialized as strings.  See [doc/cli.md](doc/cli.md) for the
input format and the full command reference.

## MCP server and web app

- `mcp/` contains a local Model Context Protocol server (`polytool-mcp`) that
  exposes the exact routines to MCP clients over stdio.  Build, install and
  client configuration are described in [mcp/README.md](mcp/README.md); the
  remaining notes on the tool set are in [doc/mcp.md](doc/mcp.md).
- `web/` contains a browser wrapper over the same Rust code (WebAssembly); see
  [doc/web.md](doc/web.md) for the arbitrary-precision contract and the build.
- `crates/exact-real-roots` is the standalone exact root-counting crate that
  polytool depends on; see
  [crates/exact-real-roots/README.md](crates/exact-real-roots/README.md).

## Documentation

- [doc/cli.md](doc/cli.md): CLI input format and general commands
  (real-rootedness, interlacing, properties, Bernstein basis, sequences,
  family checks, BKW scouting, resultants, enumeration, benchmarks).
  - [doc/cli-oeis.md](doc/cli-oeis.md): bundled OEIS catalog, provenance,
    adding sequences.
  - [doc/cli-networks.md](doc/cli-networks.md): planar networks, TNN triangles,
    strip digraphs, Pólya frequency rows.
  - [doc/cli-recurrence.md](doc/cli-recurrence.md): recurrence search, complete
    spaces and certificates, machine-readable recurrence JSON.
  - [doc/cli-recurrence-transforms.md](doc/cli-recurrence-transforms.md):
    coordinate transformations, root-domain information, implicit relations.
  - [doc/cli-recurrence-representations.md](doc/cli-recurrence-representations.md):
    search for equivalent representations.
  - [doc/cli-ehrhart-sieving.md](doc/cli-ehrhart-sieving.md): Ehrhart and
    h\*-vectors, cyclic and bicyclic sieving, Stapledon decomposition.
- [doc/library.md](doc/library.md): Rust library usage.
- [doc/mcp.md](doc/mcp.md) and [mcp/README.md](mcp/README.md): MCP server.
- [doc/web.md](doc/web.md): web wrapper contract and build.
- [doc/development.md](doc/development.md): tests and checked examples.
- [doc/algorithms.md](doc/algorithms.md): algorithm notes.
- [doc/bezout-interlacing.md](doc/bezout-interlacing.md) and
  [doc/athanasiadis-wagner-lace.md](doc/athanasiadis-wagner-lace.md):
  mathematical notes.
- [CHANGELOG.md](CHANGELOG.md): release history.

## Support

Report issues in the Git repository, or contact Per Alexandersson
(@PerAlexandersson, <per.w.alexandersson@gmail.com>).
