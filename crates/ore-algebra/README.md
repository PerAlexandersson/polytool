# ore-algebra

Exact Ore-algebra arithmetic for sequences of polynomials `F_n(x)`: the
shift `S` (`(S F)_n = F_(n+1)`), the derivative `D = d/dx`, and
multiplication by polynomials in the index `n` and the variable `x`, with
the commutation rules `S n = (n + 1) S` and `D x = x D + 1`.

- `BivariatePolynomial`: exact polynomials in `(n, y)` over the rationals;
- `OreOperator`: normal-ordered operators `sum c_(i,d)(n, x) D^d S^i`, their
  composition and their action on rows of polynomials;
- the Euler form with `theta = x D` (`to_euler_form`, `from_euler_form`);
- `exact_kernel_vectors_by_group`: exact kernel vectors of integer matrices
  by modular elimination, Chinese remaindering, rational reconstruction and
  an exact final check;
- `right_pseudo_division`: right pseudo-division by an operator whose leader
  is free of `x` and of derivatives;
- `RowMapStep`, `pull_back`, `push_forward`: maps between sequences (Weyl
  operators, gauges, affine maps, reversals, coefficient multipliers) and the
  transport of annihilators through them;
- `find_annihilator_system`: additional annihilators of the rows of a
  recurrence-defined sequence (row annihilators, runnable rules of lower
  order), fitted to the rows and returned with joint invariance identities;
- `find_image_annihilator`: a runnable annihilator of an image sequence,
  pushed forward or found as a common left multiple from an increasing
  ansatz, with a certificate `w(n, x) E = H o L + Z o A` that replays on its
  own (or the report that none was found within the budget).

Everything is exact (`num_rational::BigRational`); coefficient vectors are
in ascending degree order.

```rust
use ore_algebra::{BivariatePolynomial, OreOperator, Rational};

let r = |v: i64| Rational::from_integer(v.into());
// S o n = (n + 1) S.
let s = OreOperator::single(1, 0, BivariatePolynomial::one());
let n = BivariatePolynomial::from_index_polynomial(&[r(0), r(1)]);
let n_plus_one = BivariatePolynomial::from_index_polynomial(&[r(1), r(1)]);
let expected = OreOperator::single(1, 0, n_plus_one);
assert_eq!(s.compose(&OreOperator::single(0, 0, n)), expected);
```

`polytool` builds its recurrence transport, lifts and image recurrences on
this crate.
