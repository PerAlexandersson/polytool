# polytool for Python

Python bindings of [polytool](https://github.com/PerAlexandersson/polytool), a
Rust library for exact computations with univariate polynomials and polynomial
sequences: real-rootedness, real-root isolation, interlacing, coefficient
properties (log-concavity, unimodality, gamma-positivity), recurrence search,
real-rootedness proofs for recurrence-defined families, and a bundled catalog
of OEIS polynomial families.

All arithmetic is exact: coefficients are Python `int` and
`fractions.Fraction` values in ascending degree order, and floats are refused.

```python
import polytool

p = polytool.Polynomial([1, 11, 11, 1])          # 1 + 11t + 11t^2 + t^3
p.is_real_rooted(), p.gamma_coefficients()       # (True, [1, 8])

rows = polytool.oeis.generate("A008292", 8)       # Eulerian triangle
recurrence = polytool.find_recurrence(rows)
recurrence.prove_real_rooted().status            # 'proved'
```

Build from source with [maturin](https://www.maturin.rs/) (Rust 1.85 or
newer):

```sh
pip install maturin
maturin develop --release        # in this directory, inside a virtualenv
```

The full guide is
[doc/python.md](https://github.com/PerAlexandersson/polytool/blob/main/doc/python.md).
The bindings are a thin adapter over the Rust library: names, conventions and
JSON schemas are those of the library and the `polytool` command-line tool.
