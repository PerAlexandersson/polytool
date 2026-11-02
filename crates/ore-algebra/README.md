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
  own (or the report that none was found within the budget);
- `DifferentialDifferenceOperator`: operators `sum c(n, x) D^d E^e S^i` with
  the variable shift `E f(x) = f(x + 1)` (`E x = (x + 1) E`, integer powers
  of `E`), with composition, action on rows, pseudo-division by a recurrence
  (`difference_right_pseudo_division`) and in `E`
  (`variable_shift_pseudo_division`);
- the falling-factorial basis: the basis change `(x)_k <-> x^k`
  (`falling_factorial_to_power_basis`, `power_to_falling_factorial_basis`),
  its conjugation `x -> x (1 + D)`, `E -> 1 + D` between difference and
  differential operators (`conjugate_to_power_basis`,
  `conjugate_to_falling_factorial_basis`) and the falling-factorial Euler
  form with `X = x E^(-1)` and `Theta = x (1 - E^(-1))`;
- `DifferenceRowMapStep`, `pull_back_difference`, `push_forward_difference`,
  `find_difference_image_annihilator`: maps with `E` (`f(x) -> f(x + c)`,
  `Delta^j`, basis changes, the existing steps, substitutions) and certified
  runnable image annihilators found by push-forward or as common left
  multiples;
- `LaurentOreOperator`: coefficients in `Q[n, x, x^(-1)]`, clearing to
  polynomial operators, and the reversal, power gauge `x^(w_n)` and dilation
  as conjugations;
- `Substitution`, `pull_back_substitution`, `push_forward_substitution`:
  affine maps, reversals (degree slopes of either sign), weighted Möbius
  maps, the `f`-polynomial map and power gauges, with exact row maps,
  inverses, pullbacks and forward rules.

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
