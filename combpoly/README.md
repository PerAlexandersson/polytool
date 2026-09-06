# combpoly

A CLI tool for exploring combinatorial polynomials from permutations, words, and parking functions.

Generate combinatorial objects, filter by pattern avoidance and other criteria,
compute generating polynomials from classical statistics, check real-rootedness
and log-concavity, and search for polynomial recurrences.

## Building

```bash
cargo build --release
```

The binary is at `target/release/combpoly`.

## Subcommands

### `poly` -- Compute a generating polynomial

```bash
# Eulerian polynomial A_4(t) = 1 + 11t + 11t^2 + t^3
combpoly poly --perms 4 --stat des

# Ascent polynomial over 312-avoiding permutations of S_7, check real-rootedness
combpoly poly --perms 7 --avoiding 312 --stat asc --real-rooted

# Descent polynomial over a Bruhat lower ideal
combpoly poly --bruhat-ideal 3,2,1 --stat des --real-rooted

# Multiset rook-Eulerian: words with content (2,2,1) on board 22233
combpoly poly --content 2,2,1 --board 22233 --stat asc --real-rooted

# Parking functions of size 5
combpoly poly --parking 5 --stat des

# Noncrossing Chow polynomial: descents on tieless parking functions
combpoly poly --parking 5 --tieless --stat des

# Narayana polynomial: descents on strict-peakless, tieless parking functions
combpoly poly --parking 4 --strict-peakless --tieless --stat des

# Words on a skew board
combpoly poly --content 2,2,1 --board 33333 --skew 11 --stat asc
```

### `scan` -- Scan permutations with ideal polynomials

For each permutation in S_n (optionally filtered), compute the generating polynomial
over its Bruhat or weak lower ideal and report.

```bash
# For each 312-avoiding perm in S_7, check excedance poly over Bruhat ideal
combpoly scan --size 7 --avoiding 312 --ideal bruhat --stat exc --real-rooted

# Weak order scan, stop on first non-real-rooted example
combpoly scan --size 8 --ideal weak --stat des --real-rooted --halt
```

### `list` -- Enumerate objects

```bash
# List 312-avoiding permutations of S_5
combpoly list --perms 5 --avoiding 312

# List derangements of S_4
combpoly list --perms 4 --derangement

# List words with content (2,2,1) on board 22233
combpoly list --content 2,2,1 --board 22233
```

### `recurrence` -- Find polynomial recurrences

Search for a recurrence P_n(t) = f(n,t) P_{n-1}(t) + ... among a family
of polynomials indexed by n.

```bash
# Auto-search for Eulerian recurrence
combpoly recurrence --perms --max-n 8 --stat des --auto

# Search for recurrence of 132-avoiding descent polynomials
combpoly recurrence --perms --max-n 9 --avoiding 132 --stat des --auto

# Manual search with specific parameters
combpoly recurrence --perms --max-n 9 --avoiding 321 --stat des \
  --rec-len 2 --var-deg 1 --idx-deg 1 --diff-deg 0

# Search with LHS factor f(n)*P(n) = ..., needed for Catalan-type recurrences
combpoly recurrence --perms --max-n 9 --stat des \
  --auto --denom-idx-deg 1
```

### `lpm-hstar-table` -- Lattice-path matroid h*-tables from Dyck paths

Enumerate Dyck area sequences of a fixed semilength, build the lattice-path
matroid from the peak intervals, and compute the exact base-polytope
`h*`-vector by Ehrhart inversion using the alcoved cyclic-interval
inequalities and Ehrhart--Macdonald reciprocity.  Weak and strict lattice-point
counts use interval-sum prefix bounds; the slower all-subset rank-inequality
backend remains in the library as a correctness oracle for small tables.  For
library use, `lpm_hstar_from_area_sequence` computes the `h*`-vector directly
from an area sequence without enumerating bases.  Basis counts in the table are
computed by the LGV determinant for the peak-interval presentation, not by
listing bases.

```bash
# Complete area-sequence <-> h* table for semilength 4
combpoly lpm-hstar-table --semilength 4 --real-rooted

# CSV output for spreadsheet/comparison work
combpoly lpm-hstar-table --semilength 5 --format csv --real-rooted

# Complete table for every semilength up to 5
combpoly lpm-hstar-table --max-semilength 5 --format csv --real-rooted

# Include the full basis list
combpoly lpm-hstar-table --semilength 4 --bases --format json

# Inspect the cyclic-interval hyperplanes for one area sequence
combpoly lpm-hyperplanes --area 0,1,1,2
```

## Object sources

| Flag | Description |
|------|-------------|
| `--perms N` | All permutations of S_n |
| `--parking N` | All parking functions of length n |
| `--bruhat-ideal P` | Bruhat lower ideal of permutation P |
| `--weak-ideal P` | Weak lower ideal of permutation P |
| `--content A [--board L] [--skew M]` | Words with content A on board L (skew M) |

Permutations can be specified as `321` or `3,2,1` (use commas for entries >= 10).

## Filters

| Flag | Description |
|------|-------------|
| `--avoiding PAT` | Avoid pattern PAT (repeatable) |
| `--arrow-avoiding PATTERN` | Avoid an arrow pattern (repeatable) |
| `--alternating` | Only alternating (up-down) permutations |
| `--derangement` | Only derangements (no fixed points) |
| `--starts-with V` | First element is V |
| `--ends-with V` | Last element is V |
| `--tieless` | No equal adjacent letters |
| `--strict-peakless` | No indices with w_{i-1} < w_i > w_{i+1} |
| `--weak-peakless` | No indices with w_{i-1} <= w_i > w_{i+1} |

The first six filters are primarily permutation filters.  The last three are
word filters and are especially useful with `--parking`.

## Statistics

| Name | Description |
|------|-------------|
| `des` | Descents: #{i : w_i > w_{i+1}} |
| `asc` | Ascents: #{i : w_i < w_{i+1}} |
| `tie` | Ties: #{i : w_i = w_{i+1}} |
| `exc` | Excedances: #{i : w_i > i} |
| `peak` | Strict peaks: #{i : w_{i-1} < w_i > w_{i+1}} |
| `weak-peak` | Weak-left peaks: #{i : w_{i-1} <= w_i > w_{i+1}} |
| `valley` | Valleys: #{i : w_{i-1} > w_i < w_{i+1}} |
| `inv` | Inversions: #{(i,j) : i < j, w_i > w_j} |
| `coinv` | Coinversions: #{(i,j) : i < j, w_i < w_j} |
| `maj` | Major index: sum of descent positions |
| `comaj` | Comajor index: sum of (n-i) at descent positions |
| `fix` | Fixed points: #{i : w_i = i} |
| `cyc` | Number of cycles |
| `lrmin` | Left-to-right minima |
| `lrmax` | Left-to-right maxima |
| `rlmin` | Right-to-left minima |
| `rlmax` | Right-to-left maxima |
| `charge` | Charge (Lascoux-Schutzenberger) |
| `cocharge` | Cocharge = C(n,2) - charge |
| `lis` | Length of longest increasing subsequence |
| `lds` | Length of longest decreasing subsequence |

## Polynomial checks

| Flag | Description |
|------|-------------|
| `--real-rooted` | Check if polynomial has only real roots (Sturm chains, exact arithmetic) |
| `--log-concave` | Check if coefficient sequence is log-concave |

## Exploration binaries

The `src/bin/` directory contains research exploration tools:

| Binary | Purpose |
|--------|---------|
| `backtrack_explore` | Backtrack polynomials for all (pattern, stat) pairs |
| `backtrack_len4` | Length-4 pattern systematic study |
| `backtrack_len45` | Length-4 and length-5 pattern study |
| `catalan_backtrack` | Catalan search orders exploration |
| `catalan_verify` | Fast Catalan verification (n <= 11+) |
| `peak_*` | Various peak statistic explorations |
| `exc_*` | Excedance explorations |
| `pf_*` | Parking function explorations |
| `multiset_explore` | Multiset rook-Eulerian polynomials |

Build all binaries: `cargo build --release`

Run a specific binary: `cargo run --release --bin backtrack_explore -- 9 -r`

## Research context

This tool supports computational exploration for:

- **Backtrack permutation polynomials**: Generating polynomials from
  backtracking search over pattern-avoiding permutations.
  See `paper/Backtrack-permutations.tex`.
- **Rook-Eulerian polynomials**: Descent polynomials over Bruhat intervals
  of 312-avoiding permutations and multiset generalizations.
- **Real-rootedness**: Checking whether generating polynomials from
  combinatorial statistics have only real roots.

## License

MIT
