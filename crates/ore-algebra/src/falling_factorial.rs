//! The falling-factorial basis `(x)_k = x (x - 1) ... (x - k + 1)` and the
//! bridge between difference operators (in `E`) and differential operators
//! (in `D`).
//!
//! # The basis change and its conjugation
//!
//! The coefficients of `(x)_k` are [`crate::falling_factorial`]`(k)`.
//! Let `phi` be the linear bijection of `Q[x]` with `phi((x)_k) = x^k`
//! ([`falling_factorial_to_power_basis`]; its inverse is
//! [`power_to_falling_factorial_basis`]).  The Touchard polynomials are
//! `phi(x^n) = sum_k S(n, k) x^k`, because `x^n = sum_k S(n, k) (x)_k`.
//! Since `(x + 1)_k = (x)_k + k (x)_(k-1)`, `x (x)_k = (x)_(k+1) + k (x)_k`,
//! `x (x - 1)_k = (x)_(k+1)` and `Delta (x)_k = k (x)_(k-1)`:
//!
//! ```text
//! phi o E = (1 + D) o phi,           phi o x = x (1 + D) o phi,
//! phi o x E^(-1) = x o phi,          phi o (E - 1) = D o phi.
//! ```
//!
//! So `psi: x -> x (1 + D), E -> 1 + D` (and `n -> n`, `S -> S`) is an
//! algebra homomorphism from the operators free of `D` with nonnegative
//! powers of `E` to the operators free of `E`, with `phi o P = psi(P) o phi`
//! ([`conjugate_to_power_basis`]); it respects `E x = (x + 1) E` because
//! `(1 + D) x = x (1 + D) + 1`.  In the other direction
//! `psi^(-1): x -> x E^(-1), D -> E - 1` respects `D x = x D + 1`, so it is
//! an algebra homomorphism from the operators free of `E` to the operators
//! free of `D`, with `phi o psi^(-1)(Q) = Q o phi`
//! ([`conjugate_to_falling_factorial_basis`]); both compositions
//! `psi o psi^(-1)` and `psi^(-1) o psi` are the identity on generators.
//!
//! # The falling-factorial Euler form
//!
//! The image of `x` under `psi^(-1)` is the raising operator
//! `X = x E^(-1)` (`X (x)_k = (x)_(k+1)`), and that of `theta = x D` is the
//! difference Euler operator `Theta = x (1 - E^(-1)) = x nabla`, which is
//! diagonal: `Theta (x)_k = k (x)_k`.  Applying `psi^(-1)` to the Euler
//! form `x^J psi(P) = sum x^s a_(i,s)(n, theta) S^i` of `psi(P)` gives
//!
//! ```text
//! X^J P = sum X^s a_(i,s)(n, Theta) S^i,
//! ```
//!
//! which acts on the falling-factorial coefficients `f_k` of a row
//! `F = sum f_k (x)_k` by `(x)_k -> a(n, k) (x)_(k+s)`, just as the Euler
//! form acts on power coefficients
//! ([`to_falling_factorial_euler_form`], [`from_falling_factorial_euler_form`]).

use crate::bivariate::BivariatePolynomial;
use crate::difference_row_map::DifferenceRowMapError;
use crate::differential_difference::DifferentialDifferenceOperator;
use crate::euler::{from_euler_form, to_euler_form, EulerOperator};
use crate::operator::{rational, OreOperator};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Zero};

/// Stirling numbers of the second kind `S(m, k)`, `m, k <= max`.
fn stirling_second_kind(max: usize) -> Vec<Vec<Rational>> {
    let mut table = vec![vec![Rational::zero(); max + 1]; max + 1];
    table[0][0] = Rational::one();
    for m in 1..=max {
        for k in 1..=m {
            table[m][k] = &table[m - 1][k - 1] + rational(k) * &table[m - 1][k];
        }
    }
    table
}

/// The ascending coefficients of `(x)_k` (signed Stirling numbers of the
/// first kind), for `k <= max`.
fn falling_factorials(max: usize) -> Vec<Vec<Rational>> {
    let mut out = vec![vec![Rational::one()]];
    for k in 1..=max {
        let next = poly::mul(&out[k - 1], &[-rational(k - 1), Rational::one()]);
        out.push(next);
    }
    out
}

/// `phi(F)` for `phi((x)_k) = x^k`: the falling-factorial coefficients of
/// `F` (given by its ascending power coefficients), as power coefficients.
/// Canonical (zero is `[0]`).
pub fn falling_factorial_to_power_basis(row: &[Rational]) -> Vec<Rational> {
    let row = poly::trim(row);
    let table = stirling_second_kind(row.len() - 1);
    let mut out = vec![Rational::zero(); row.len()];
    for (j, c) in row.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        for (k, s) in table[j].iter().enumerate().take(j + 1) {
            out[k] += c * s;
        }
    }
    poly::trim(&out)
}

/// `phi^(-1)(F)`, `x^k -> (x)_k`, as ascending power coefficients.
/// Canonical (zero is `[0]`).
pub fn power_to_falling_factorial_basis(row: &[Rational]) -> Vec<Rational> {
    let row = poly::trim(row);
    let factorials = falling_factorials(row.len() - 1);
    let mut out = vec![Rational::zero()];
    for (k, c) in row.iter().enumerate() {
        if !c.is_zero() {
            out = poly::add(&out, &poly::scale(&factorials[k], c));
        }
    }
    poly::trim(&out)
}

/// The raising operator `X = x E^(-1)`, `X (x)_k = (x)_(k+1)`.
pub fn falling_factorial_raising_operator() -> DifferentialDifferenceOperator {
    DifferentialDifferenceOperator::single(
        0,
        0,
        -1,
        BivariatePolynomial::monomial(0, 1, Rational::one()),
    )
}

/// The difference Euler operator `Theta = x (1 - E^(-1))`,
/// `Theta (x)_k = k (x)_k`.
pub fn falling_factorial_euler_operator() -> DifferentialDifferenceOperator {
    DifferentialDifferenceOperator::multiplication(BivariatePolynomial::monomial(
        0,
        1,
        Rational::one(),
    ))
    .subtract(&falling_factorial_raising_operator())
}

/// `sum_b c_b(n) image^b` for `c(n, x) = sum_b c_b(n) x^b`, with `image` an
/// operator that commutes with polynomials in `n` (computed with the
/// composition `compose`).
fn substitute_variable<T: Clone>(
    c: &BivariatePolynomial,
    powers: &mut Vec<T>,
    image: &T,
    compose: impl Fn(&T, &T) -> T,
    left_multiply: impl Fn(&T, &BivariatePolynomial) -> T,
    add: impl Fn(&T, &T) -> T,
    zero: T,
) -> T {
    let degree = c.variable_degree().unwrap_or(0);
    while powers.len() <= degree {
        let next = compose(powers.last().expect("nonempty"), image);
        powers.push(next);
    }
    let mut out = zero;
    for (b, power) in powers.iter().enumerate().take(degree + 1) {
        let column = c.variable_coefficient(b);
        if poly::is_zero(&column) {
            continue;
        }
        out = add(
            &out,
            &left_multiply(power, &BivariatePolynomial::from_index_polynomial(&column)),
        );
    }
    out
}

/// `psi(P)` for an operator `P` free of `D` with nonnegative powers of `E`:
/// `x -> x (1 + D)`, `E -> 1 + D` (module documentation), so that
/// `phi o P = psi(P) o phi` on every row.
///
/// # Errors
///
/// [`DifferenceRowMapError::DerivativeInBasisConjugation`] when `P` involves
/// `D`, [`DifferenceRowMapError::NegativeVariableShift`] when it involves a
/// negative power of `E` (multiply by a power of `E` on the left first).
pub fn conjugate_to_power_basis(
    op: &DifferentialDifferenceOperator,
) -> Result<OreOperator, DifferenceRowMapError> {
    if !op.is_derivative_free() {
        return Err(DifferenceRowMapError::DerivativeInBasisConjugation);
    }
    if op.min_variable_shift() < 0 {
        return Err(DifferenceRowMapError::NegativeVariableShift);
    }
    // x (1 + D) and 1 + D.
    let one_plus_derivative = {
        let mut t = OreOperator::identity();
        t.add_term((0, 1), BivariatePolynomial::one());
        t
    };
    let x_image =
        one_plus_derivative.left_multiply(&BivariatePolynomial::monomial(0, 1, Rational::one()));
    let mut x_powers = vec![OreOperator::identity()];
    let mut e_powers = vec![OreOperator::identity()];
    let mut out = OreOperator::default();
    for ((i, _, e), c) in &op.0 {
        let e = usize::try_from(*e).expect("nonnegative power");
        while e_powers.len() <= e {
            let next = e_powers
                .last()
                .expect("nonempty")
                .compose(&one_plus_derivative);
            e_powers.push(next);
        }
        let coefficient = substitute_variable(
            c,
            &mut x_powers,
            &x_image,
            OreOperator::compose,
            OreOperator::left_multiply,
            OreOperator::add,
            OreOperator::default(),
        );
        let term = coefficient
            .compose(&e_powers[e])
            .compose(&OreOperator::single(*i, 0, BivariatePolynomial::one()));
        out.add_assign(&term);
    }
    Ok(out)
}

/// `psi^(-1)(Q)` for an operator `Q` free of `E`: `x -> x E^(-1)`,
/// `D -> E - 1` (module documentation), so that
/// `phi o psi^(-1)(Q) = Q o phi` on every row.  The result is free of `D`.
pub fn conjugate_to_falling_factorial_basis(op: &OreOperator) -> DifferentialDifferenceOperator {
    let x_image = falling_factorial_raising_operator();
    let difference = DifferentialDifferenceOperator::forward_difference(1);
    let mut x_powers = vec![DifferentialDifferenceOperator::identity()];
    let mut d_powers = vec![DifferentialDifferenceOperator::identity()];
    let mut out = DifferentialDifferenceOperator::default();
    for ((i, d), c) in &op.0 {
        while d_powers.len() <= *d {
            let next = d_powers.last().expect("nonempty").compose(&difference);
            d_powers.push(next);
        }
        let coefficient = substitute_variable(
            c,
            &mut x_powers,
            &x_image,
            DifferentialDifferenceOperator::compose,
            DifferentialDifferenceOperator::left_multiply,
            DifferentialDifferenceOperator::add,
            DifferentialDifferenceOperator::default(),
        );
        let term = coefficient
            .compose(&d_powers[*d])
            .compose(&DifferentialDifferenceOperator::index_shift(*i));
        out.add_assign(&term);
    }
    out
}

/// `(J, X^J P)` in the falling-factorial Euler form
/// `sum X^s a_(i,s)(n, Theta) S^i`, keyed by `(shift i, power s of X)`, with
/// coefficients in `(n, Theta)` (module documentation): the Euler form of
/// `psi(P)` read with `x -> X`, `theta -> Theta`.
///
/// # Errors
///
/// As [`conjugate_to_power_basis`].
pub fn to_falling_factorial_euler_form(
    op: &DifferentialDifferenceOperator,
) -> Result<(usize, EulerOperator), DifferenceRowMapError> {
    Ok(to_euler_form(&conjugate_to_power_basis(op)?))
}

/// The operator `sum X^s a_(i,s)(n, Theta) S^i` of a falling-factorial
/// Euler form, normal-ordered (free of `D`).
pub fn from_falling_factorial_euler_form(op: &EulerOperator) -> DifferentialDifferenceOperator {
    conjugate_to_falling_factorial_basis(&from_euler_form(op))
}
