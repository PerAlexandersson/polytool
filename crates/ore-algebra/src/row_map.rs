//! Maps between sequences of polynomials, and the transport of annihilating
//! operators through them.
//!
//! A [`RowMapStep`] sends a sequence `F_n(x)` to a sequence `G_n(x)`.  For
//! an operator `M` acting on `G`, [`pull_back`] computes an operator `E`
//! acting on `F` and a polynomial `b(n)` with
//!
//! ```text
//! (M G)_n = 0   whenever   (E F)_n = 0  and  b(n) != 0,
//! ```
//!
//! and [`push_forward`] turns an annihilator of `F` into a candidate
//! annihilator of `G` (a candidate only: it is certified by pulling it back,
//! see [`crate::find_image_annihilator`]).
//!
//! # The steps and their pullbacks
//!
//! Write `M = sum m_(i,d)(n, x) D^d S^i` with bounds `R >= ` its largest
//! shift, `J >= ` its largest derivative order and `b >= ` the degree in `x`
//! of its coefficients ([`PullbackBounds`]).  The bounds, not the actual
//! shape, enter the formulas, so for fixed bounds the pullback is linear in
//! `M` (over the rationals and over polynomials in `n`).
//!
//! - **Weyl operator** `G_n = T_n(F_n)`, `T` shift-free with polynomial
//!   coefficients: `E = M o T` (composition in the Ore algebra, which is the
//!   composition of the actions), `b = 1`; exactly `(M G)_n = (E F)_n`.
//! - **Gauge** `G_n = mu(n) F_n` with `mu(n)/mu(n - 1) = alpha(n)/beta(n)`:
//!   `E = sum m_(i,d) prod_(t=1..i) alpha(n + t) prod_(t=i+1..R) beta(n + t) D^d S^i`
//!   and `b(n) = prod_(t=1..R) beta(n + t)`, from
//!   `mu(n + i) = mu(n) prod_(t=1..i) alpha(n + t)/beta(n + t)`:
//!   `(M G)_n b(n) = mu(n) (E F)_n` wherever `mu` is finite and nonzero on
//!   `n, ..., n + R`.
//! - **Affine** `G_n(x) = F_n(a x + c)`, `a != 0`: substitute
//!   `x -> (y - c)/a` in every coefficient and multiply the `D^d` term by
//!   `a^d`; then `(M G)_n(x) = (E F)_n(a x + c)`, `b = 1`.
//! - **Reversal** `G_n(x) = x^(d_n) F_n(1/x)` with `d_n = d_0 + h n`
//!   (`h >= 0` an integer): with `C_n(K)(x) = x^(d_n) K(1/x)`,
//!   `D o C_n = C_n o rho(n)`, `rho(n) = d_n y - y^2 D`,
//!   `c(n, x) C_n(K) = C_n(c(n, 1/y) K)` and `C_(n+i)(K) = C_n(y^(-h i) K)`;
//!   multiplying by `y^N`, `N = b + h R`, to clear negative powers gives
//!   `E = sum y^(N - h i) m_(i,d)(n, 1/y) rho(n + i)^d S^i` and
//!   `(M G)_n = C_n(y^(-N) (E F)_n)`, an injective map of Laurent
//!   polynomials; `b = 1`.  (The identity holds in Laurent polynomials, also
//!   where `deg F_n > d_n`.)
//! - **Coefficient multiplier** `G_n = sum_k u(k) f_(n,k) x^k` with
//!   `u(0) = 1`, `u(k + 1)/u(k) = A(k)/B(k)`, where `A(k)` and `B(k)` are
//!   nonzero at every integer `k >= 0` (so `u` is finite and nonzero).  In
//!   the Euler form `x^J M = sum x^s c_(i,s)(n, theta) S^i` (`theta = x D`,
//!   `s <= S = b + J`) put
//!   `E = sum x^s c_(i,s)(n, theta) prod_(t=0..s-1) B(theta + t) prod_(t=s-S..-1) A(theta + t) S^i`.
//!   The coefficient of `x^m` in `x^s c(n, theta) G_(n+i)` is
//!   `c(n, m - s) u(m - s) f_(n+i, m-s)`, and
//!   `u(m - s) = u(m) prod_(j=m-s..m-1) B(j)/A(j)` for `m - s >= 0`, so the
//!   coefficient of `x^m` in `(E F)_n` is `V(m)/u(m)` times that of
//!   `x^J (M G)_n`, with `V(m) = prod_(j=m-S..m-1) A(j)` (terms with
//!   `m - s < 0` vanish on both sides).  Hence `(M G)_n = 0` when
//!   `(E F)_n = 0`, provided `A` is also nonzero at the integers
//!   `-S, ..., -1` (so `V(m) != 0` for every `m >= 0`); `b = 1`.  The bound
//!   `S` is returned by [`PullbackBounds::multiplier_euler_bound`].
//!
//! The pushforwards: a gauge pushes forward by
//! `m'_(i,d) prod_(t=i+1..R) alpha(n + t) prod_(t=1..i) beta(n + t)` (the
//! annihilator of `F` read on `G = mu F`), an affine map and a coefficient
//! multiplier by the pullback through the inverse map (`x -> (x - c)/a`,
//! the ratio `B/A`), a reversal by the pullback through itself (it is an
//! involution of Laurent polynomials); a Weyl operator has no pushforward.

use crate::bivariate::BivariatePolynomial;
use crate::euler::{add_euler_term, from_euler_form, to_euler_form, EulerOperator};
use crate::operator::OreOperator;
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, Signed, Zero};
use std::fmt;

/// One map `F -> G` between sequences of polynomials (module
/// documentation).  Polynomials are ascending coefficient vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RowMapStep {
    /// `G_n = T_n(F_n)` for a shift-free operator `T`.
    WeylOperator(OreOperator),
    /// `G_n = mu(n) F_n` with `mu(n)/mu(n - 1) = numerator(n)/denominator(n)`.
    Gauge {
        /// `alpha(n)`.
        numerator: Vec<Rational>,
        /// `beta(n)`.
        denominator: Vec<Rational>,
    },
    /// `G_n(x) = F_n(slope x + intercept)`, `slope != 0`.
    Affine {
        /// `a`.
        slope: Rational,
        /// `c`.
        intercept: Rational,
    },
    /// `G_n(x) = x^(d_n) F_n(1/x)` with `d_n = degree[0] + degree[1] n`,
    /// `degree[1]` a nonnegative integer.
    Reversal {
        /// `d_n` as an ascending polynomial of degree at most one in `n`.
        degree: Vec<Rational>,
    },
    /// `x^k -> u(k) x^k` with `u(0) = 1` and
    /// `u(k + 1)/u(k) = numerator(k)/denominator(k)`.
    CoefficientMultiplier {
        /// `A(k)`.
        numerator: Vec<Rational>,
        /// `B(k)`.
        denominator: Vec<Rational>,
    },
}

/// Why a [`RowMapStep`] cannot be used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RowMapError {
    /// A Weyl operator has a shift.
    ShiftedWeylOperator,
    /// A gauge or multiplier polynomial is zero.
    ZeroPolynomial,
    /// An affine slope is zero.
    ZeroSlope,
    /// A reversal degree is not `d_0 + h n` with an integer `h >= 0`.
    InvalidReversalDegree,
}

impl fmt::Display for RowMapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::ShiftedWeylOperator => "a Weyl-operator step must be shift-free",
            Self::ZeroPolynomial => "a gauge or multiplier polynomial is zero",
            Self::ZeroSlope => "an affine slope must be nonzero",
            Self::InvalidReversalDegree => {
                "a reversal degree must be d_0 + h n with an integer h >= 0"
            }
        })
    }
}

impl std::error::Error for RowMapError {}

/// Bounds on an operator `M` that fix the pullback formulas: its largest
/// shift, derivative order and coefficient degree in `x` may be smaller.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PullbackBounds {
    /// `R`, at least the largest shift.
    pub top_shift: usize,
    /// `J`, at least the largest derivative order.
    pub derivative_order: usize,
    /// `b`, at least the degree in `x` of every coefficient.
    pub variable_degree: usize,
}

impl PullbackBounds {
    /// The exact bounds of `op`.
    pub fn of(op: &OreOperator) -> Self {
        Self {
            top_shift: op.max_shift(),
            derivative_order: op.max_derivative(),
            variable_degree: op
                .0
                .values()
                .filter_map(BivariatePolynomial::variable_degree)
                .max()
                .unwrap_or(0),
        }
    }

    /// `S = b + J`, the largest power of `x` in the Euler form `x^J M`.
    pub fn multiplier_euler_bound(&self) -> usize {
        self.variable_degree + self.derivative_order
    }

    /// Bounds on the pullback through `step` (valid upper bounds, used to
    /// chain several steps).
    pub(crate) fn after(&self, step: &RowMapStep) -> Self {
        match step {
            RowMapStep::WeylOperator(t) => {
                let bounds = Self::of(t);
                Self {
                    top_shift: self.top_shift,
                    derivative_order: self.derivative_order + bounds.derivative_order,
                    variable_degree: self.variable_degree + bounds.variable_degree,
                }
            }
            RowMapStep::Gauge { .. } | RowMapStep::Affine { .. } => *self,
            RowMapStep::Reversal { degree } => {
                let step = reversal_step(degree).unwrap_or(0);
                Self {
                    top_shift: self.top_shift,
                    derivative_order: self.derivative_order,
                    variable_degree: self.variable_degree
                        + step * self.top_shift
                        + 2 * self.derivative_order,
                }
            }
            RowMapStep::CoefficientMultiplier {
                numerator,
                denominator,
            } => {
                let s = self.multiplier_euler_bound();
                let growth = poly::degree(numerator)
                    .unwrap_or(0)
                    .max(poly::degree(denominator).unwrap_or(0));
                let theta_degree = self.derivative_order + s * growth;
                Self {
                    top_shift: self.top_shift,
                    derivative_order: theta_degree,
                    variable_degree: s + theta_degree,
                }
            }
        }
    }
}

/// The integer slope `h >= 0` of a reversal degree `d_0 + h n`.
fn reversal_step(degree: &[Rational]) -> Option<usize> {
    let degree = poly::trim(degree);
    if degree.len() > 2 {
        return None;
    }
    let h = degree.get(1).cloned().unwrap_or_else(Rational::zero);
    if !h.is_integer() || h.is_negative() {
        return None;
    }
    usize::try_from(h.to_integer()).ok()
}

impl RowMapStep {
    /// Check the parameters (module documentation); the nonvanishing of
    /// gauge and multiplier polynomials at integers is the caller's to
    /// decide.
    pub fn validate(&self) -> Result<(), RowMapError> {
        match self {
            Self::WeylOperator(t) => {
                if t.max_shift() != 0 {
                    return Err(RowMapError::ShiftedWeylOperator);
                }
            }
            Self::Gauge {
                numerator,
                denominator,
            }
            | Self::CoefficientMultiplier {
                numerator,
                denominator,
            } => {
                if poly::is_zero(numerator) || poly::is_zero(denominator) {
                    return Err(RowMapError::ZeroPolynomial);
                }
            }
            Self::Affine { slope, .. } => {
                if slope.is_zero() {
                    return Err(RowMapError::ZeroSlope);
                }
            }
            Self::Reversal { degree } => {
                reversal_step(degree).ok_or(RowMapError::InvalidReversalDegree)?;
            }
        }
        Ok(())
    }

    /// The image `G_n` of one row `F_n` (exact; `None` where a reversal
    /// degree is not a nonnegative integer at `n` or below the row's degree).
    pub fn apply_row(&self, row: &[Rational], n: usize) -> Option<Vec<Rational>> {
        let index = Rational::from_integer(n.into());
        let row = poly::trim(row);
        Some(match self {
            Self::WeylOperator(t) => t.apply_row(&row, n),
            Self::Gauge { .. } => return None,
            Self::Affine { slope, intercept } => {
                poly::compose(&row, &[intercept.clone(), slope.clone()])
            }
            Self::Reversal { degree } => {
                let d = poly::evaluate(degree, &index);
                if !d.is_integer() || d.is_negative() {
                    return None;
                }
                let d = usize::try_from(d.to_integer()).ok()?;
                if poly::degree(&row).is_some_and(|degree| degree > d) {
                    return None;
                }
                let mut out = vec![Rational::zero(); d + 1];
                for (k, c) in row.iter().enumerate() {
                    out[d - k] = c.clone();
                }
                poly::trim(&out)
            }
            Self::CoefficientMultiplier {
                numerator,
                denominator,
            } => {
                let mut u = Rational::one();
                let mut out = Vec::with_capacity(row.len());
                for (k, c) in row.iter().enumerate() {
                    if k > 0 {
                        let kk = Rational::from_integer((k - 1).into());
                        let b = poly::evaluate(denominator, &kk);
                        if b.is_zero() {
                            return None;
                        }
                        u = u * poly::evaluate(numerator, &kk) / b;
                    }
                    out.push(c * &u);
                }
                poly::trim(&out)
            }
        })
    }
}

/// `p(n)` as a coefficient.
fn index_polynomial(p: &[Rational]) -> BivariatePolynomial {
    BivariatePolynomial::from_index_polynomial(p)
}

/// `prod_(t = from..=to) p(n + t)` (one when `from > to`), for `t` signed.
fn shifted_product(p: &[Rational], from: i64, to: i64) -> BivariatePolynomial {
    let p = index_polynomial(p);
    (from..=to).fold(BivariatePolynomial::one(), |acc, t| {
        acc.multiply(&p.shift_index(t))
    })
}

/// `prod_(t = from..=to) q(theta + t)` as a polynomial in `theta` (the
/// variable of a [`BivariatePolynomial`]).
fn theta_shifted_product(q: &[Rational], from: i64, to: i64) -> BivariatePolynomial {
    let q = BivariatePolynomial::from_variable_polynomial(q);
    (from..=to).fold(BivariatePolynomial::one(), |acc, t| {
        acc.multiply(&q.shift(&Rational::zero(), &Rational::from_integer(t.into())))
    })
}

/// The gauge pullback (module documentation).
fn gauge_pullback(
    op: &OreOperator,
    alpha: &[Rational],
    beta: &[Rational],
    top: usize,
) -> OreOperator {
    let mut out = OreOperator::default();
    for ((i, d), c) in &op.0 {
        let factor = shifted_product(alpha, 1, *i as i64).multiply(&shifted_product(
            beta,
            *i as i64 + 1,
            top as i64,
        ));
        out.add_term((*i, *d), c.multiply(&factor));
    }
    out
}

/// The gauge pushforward (module documentation).
fn gauge_pushforward(op: &OreOperator, alpha: &[Rational], beta: &[Rational]) -> OreOperator {
    let top = op.max_shift();
    let mut out = OreOperator::default();
    for ((i, d), c) in &op.0 {
        let factor = shifted_product(alpha, *i as i64 + 1, top as i64)
            .multiply(&shifted_product(beta, 1, *i as i64));
        out.add_term((*i, *d), c.multiply(&factor));
    }
    out
}

/// Substitute `y -> (y - c)/a` and multiply the `D^d` term by `a^d`.
fn affine_pullback(op: &OreOperator, a: &Rational, c: &Rational) -> OreOperator {
    let inverse = [-(c / a), a.recip()];
    let mut out = OreOperator::default();
    for ((i, d), coefficient) in &op.0 {
        let rows: Vec<Vec<Rational>> = coefficient
            .rows()
            .iter()
            .map(|row| poly::trim(&poly::compose(row, &inverse)))
            .collect();
        let mut power = Rational::one();
        for _ in 0..*d {
            power *= a;
        }
        out.add_term((*i, *d), BivariatePolynomial::new(rows).scale(&power));
    }
    out
}

/// `rho(n) = d_n y - y^2 D`.
fn reversal_derivative(degree: &[Rational]) -> OreOperator {
    let mut rho = OreOperator::single(0, 1, BivariatePolynomial::monomial(0, 2, -Rational::one()));
    rho.add_term(
        (0, 0),
        index_polynomial(degree).multiply(&BivariatePolynomial::monomial(0, 1, Rational::one())),
    );
    rho
}

/// The reversal pullback with `y^N`, `N = top` (module documentation);
/// requires `deg_x m_(i,d) + h i <= top` for every term.
fn reversal_pullback(
    op: &OreOperator,
    degree: &[Rational],
    step: usize,
    top: usize,
) -> OreOperator {
    let rho = reversal_derivative(degree);
    let mut powers = vec![OreOperator::identity()];
    let mut out = OreOperator::default();
    for ((i, d), c) in &op.0 {
        while powers.len() <= *d {
            let next = rho.compose(powers.last().expect("nonempty"));
            powers.push(next);
        }
        let shift = top - step * i;
        let rows: Vec<Vec<Rational>> = c
            .rows()
            .iter()
            .map(|row| {
                let mut reversed = vec![Rational::zero(); shift + 1];
                for (j, value) in row.iter().enumerate() {
                    reversed[shift - j] = value.clone();
                }
                reversed
            })
            .collect();
        let term = powers[*d]
            .shifted(*i)
            .left_multiply(&BivariatePolynomial::new(rows));
        out.add_assign(&term);
    }
    out
}

/// The multiplier pullback with the Euler bounds `(J, S)` (module
/// documentation); requires `J >=` the largest derivative order and
/// `S >=` the largest power of `x` in the Euler form `x^J M`.
fn multiplier_pullback(
    op: &OreOperator,
    numerator: &[Rational],
    denominator: &[Rational],
    derivative_bound: usize,
    euler_bound: usize,
) -> OreOperator {
    let (j_op, euler) = to_euler_form(op);
    let extra = derivative_bound - j_op;
    let mut out = EulerOperator::new();
    for ((i, s0), c) in &euler {
        let s = s0 + extra;
        let factor = theta_shifted_product(denominator, 0, s as i64 - 1).multiply(
            &theta_shifted_product(numerator, s as i64 - euler_bound as i64, -1),
        );
        add_euler_term(&mut out, (*i, s), c.multiply(&factor));
    }
    from_euler_form(&out)
}

/// The pullback of `op` through `steps` (applied to `F` first to last, so
/// pulled back last to first) with bounds `bounds` on `op`: the operator
/// `E` on the source sequence and the base factor `b(n)` (module
/// documentation), or the first invalid step.
///
/// # Panics
///
/// Panics when `bounds` is smaller than the actual shape of `op`.
pub fn pull_back(
    op: &OreOperator,
    steps: &[RowMapStep],
    bounds: PullbackBounds,
) -> Result<(OreOperator, Vec<Rational>), RowMapError> {
    let actual = PullbackBounds::of(op);
    assert!(
        op.is_zero()
            || (actual.top_shift <= bounds.top_shift
                && actual.derivative_order <= bounds.derivative_order
                && actual.variable_degree <= bounds.variable_degree),
        "pullback bounds below the operator's shape"
    );
    let mut current = op.clone();
    let mut current_bounds = bounds;
    let mut base_factor = vec![Rational::one()];
    for step in steps.iter().rev() {
        step.validate()?;
        current = match step {
            RowMapStep::WeylOperator(t) => current.compose(t),
            RowMapStep::Gauge {
                numerator,
                denominator,
            } => {
                let top = current_bounds.top_shift;
                let product = shifted_product(denominator, 1, top as i64);
                base_factor = poly::mul(&base_factor, &product.variable_coefficient(0));
                gauge_pullback(&current, numerator, denominator, top)
            }
            RowMapStep::Affine { slope, intercept } => affine_pullback(&current, slope, intercept),
            RowMapStep::Reversal { degree } => {
                let h = reversal_step(degree).ok_or(RowMapError::InvalidReversalDegree)?;
                let top = current_bounds.variable_degree + h * current_bounds.top_shift;
                reversal_pullback(&current, degree, h, top)
            }
            RowMapStep::CoefficientMultiplier {
                numerator,
                denominator,
            } => multiplier_pullback(
                &current,
                numerator,
                denominator,
                current_bounds.derivative_order,
                current_bounds.multiplier_euler_bound(),
            ),
        };
        current_bounds = current_bounds.after(step);
    }
    Ok((current, poly::trim(&base_factor)))
}

/// A candidate annihilator of `G` from an annihilator `op` of `F` (module
/// documentation): `None` for a Weyl-operator step (use
/// [`crate::find_image_annihilator`]) or an invalid step.
pub fn push_forward(op: &OreOperator, step: &RowMapStep) -> Option<OreOperator> {
    step.validate().ok()?;
    Some(match step {
        RowMapStep::WeylOperator(_) => return None,
        RowMapStep::Gauge {
            numerator,
            denominator,
        } => gauge_pushforward(op, numerator, denominator),
        RowMapStep::Affine { slope, intercept } => {
            affine_pullback(op, &slope.recip(), &-(intercept / slope))
        }
        RowMapStep::Reversal { degree } => {
            let h = reversal_step(degree)?;
            let bounds = PullbackBounds::of(op);
            reversal_pullback(op, degree, h, bounds.variable_degree + h * bounds.top_shift)
        }
        RowMapStep::CoefficientMultiplier {
            numerator,
            denominator,
        } => {
            let bounds = PullbackBounds::of(op);
            multiplier_pullback(
                op,
                denominator,
                numerator,
                bounds.derivative_order,
                bounds.multiplier_euler_bound(),
            )
        }
    })
}
