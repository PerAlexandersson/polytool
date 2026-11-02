//! Maps between sequences of polynomials that involve the variable shift
//! `E`, and the transport of annihilating operators with `E` through them.
//!
//! A [`DifferenceRowMapStep`] sends a sequence `F_n(x)` to a sequence
//! `G_n(x)`.  For an operator `M` acting on `G` (a
//! [`DifferentialDifferenceOperator`]), [`pull_back_difference`] computes an
//! operator `P` acting on `F` and a polynomial `b(n)` with
//!
//! ```text
//! (M G)_n = 0   whenever   (P F)_n = 0  and  b(n) != 0,
//! ```
//!
//! and [`push_forward_difference`] turns an annihilator of `F` into an
//! annihilator of `G` where one exists without a search.
//!
//! # The steps and their pullbacks
//!
//! The bounds ([`DifferencePullbackBounds`]) are upper bounds on the shape
//! of `M`; they, not the actual shape, enter the formulas, so for fixed
//! bounds the pullback is linear in `M` (over the rationals and over
//! polynomials in `n`).  Write `M = sum_e M_e E^e` with `M_e` free of `E`
//! (the operators commute, so `E^e` may stand on the right).
//!
//! - **Operator** `G_n = T_n(F_n)`, `T` shift-free (for example `E^c`,
//!   `f(x) -> f(x + c)`, or `Delta^j`): `P = M o T`, `b = 1`, and exactly
//!   `(M G)_n = (P F)_n`.
//! - **Falling factorial to power basis** `G_n = phi(F_n)`,
//!   `phi((x)_k) = x^k`: for `M` free of `E`, `P = psi^(-1)(M)` (see
//!   [`crate::conjugate_to_falling_factorial_basis`]) and
//!   `(M G)_n = phi((P F)_n)`, `phi` injective; `b = 1`.
//! - **Power to falling-factorial basis** `G_n = phi^(-1)(F_n)`: for `M` free
//!   of `D`, with `c = max(0, -e_min)` from the bounds, `E^c M` has
//!   nonnegative powers of `E` and `P = psi(E^c M)`
//!   ([`crate::conjugate_to_power_basis`]):
//!   `(M G)_n = E^(-c) phi^(-1)((P F)_n)`, both maps injective; `b = 1`.
//! - **An existing [`RowMapStep`]**.  A gauge `G_n = mu(n) F_n` commutes
//!   with `E`, and an affine map `G_n(x) = F_n(a x + c)` with an integer
//!   slope `a` satisfies `E^e G_n = (E^(a e) F_n)(a x + c)`; so `M_e E^e`
//!   pulls back to `P_e E^(a e)` (`a = 1` for a gauge) with `P_e` the
//!   pullback of `M_e` ([`crate::pull_back`], with the common bounds, so
//!   the base factor `b(n)` of a gauge is the same for every `e`).  A Weyl
//!   operator `T` gives `P = M o T`.  A reversal, a coefficient multiplier
//!   or an affine map with a non-integer slope does not commute with `E` in
//!   this way and needs `M` free of `E`; then [`crate::pull_back`] applies
//!   as it is.
//! - **A [`Substitution`]** (affine, reversal, Möbius, `f`-polynomial, power
//!   gauge) needs `M` free of `E`: [`crate::pull_back_substitution`], whose
//!   left factor is injective, so `b = 1`.
//!
//! # Push-forwards
//!
//! An operator `T = c E^e` (`c` a nonzero constant) pushes an annihilator
//! `L` of `F` forward to the conjugate `E^e L E^(-e)` (`x -> x + e` in every
//! coefficient); other operators have none (use
//! [`crate::find_difference_image_annihilator`]).  The basis changes push
//! forward by `psi(E^c L)` (`L` free of `D`, `c` making the powers of `E`
//! nonnegative) and `psi^(-1)(L)` (`L` free of `E`); a gauge by the formula
//! of [`crate::push_forward`] with the common top shift; an affine map with
//! slope `a = +-1` componentwise with `E^e -> E^(a e)`; the other existing
//! steps, for `L` free of `E`, by [`crate::push_forward`]; a substitution,
//! for `L` free of `E`, by [`crate::push_forward_substitution`].

use crate::bivariate::BivariatePolynomial;
use crate::differential_difference::DifferentialDifferenceOperator;
use crate::falling_factorial::{
    conjugate_to_falling_factorial_basis, conjugate_to_power_basis,
    falling_factorial_to_power_basis, power_to_falling_factorial_basis,
};
use crate::operator::OreOperator;
use crate::row_map::{pull_back, push_forward, PullbackBounds, RowMapError, RowMapStep};
use crate::substitution::{
    bounds_after_substitution, pull_back_substitution, push_forward_substitution, Substitution,
    SubstitutionError,
};
use crate::univariate as poly;
use crate::Rational;
use num_traits::{One, ToPrimitive, Zero};
use std::fmt;

/// One map `F -> G` between sequences of polynomials (module
/// documentation).  Polynomials are ascending coefficient vectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DifferenceRowMapStep {
    /// `G_n = T_n(F_n)` for a shift-free operator `T` (with `E` allowed).
    Operator(DifferentialDifferenceOperator),
    /// `G_n = phi(F_n)` with `phi((x)_k) = x^k` (the falling-factorial
    /// coefficients of `F_n` become power coefficients).
    FallingFactorialToPowerBasis,
    /// `G_n = phi^(-1)(F_n)`, `x^k -> (x)_k`.
    PowerToFallingFactorialBasis,
    /// An existing step (module documentation for how it meets `E`).
    Ore(RowMapStep),
    /// A substitution of the variable ([`Substitution`]); it needs operators
    /// free of `E` (for an affine map with `E`, use [`Self::Ore`]).
    Substitution(Substitution),
}

/// Why a [`DifferenceRowMapStep`] or a pullback cannot be used.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DifferenceRowMapError {
    /// An existing step is invalid.
    Ore(RowMapError),
    /// A substitution is invalid.
    Substitution(SubstitutionError),
    /// An operator step has a shift.
    ShiftedOperator,
    /// The conjugation to the power basis needs an operator free of `D`.
    DerivativeInBasisConjugation,
    /// The conjugation to the power basis needs nonnegative powers of `E`.
    NegativeVariableShift,
    /// The conjugation to the falling-factorial basis needs an operator free
    /// of `E`.
    VariableShiftInBasisConjugation,
    /// A reversal, a coefficient multiplier, an affine map with a
    /// non-integer slope or a substitution needs an operator free of `E`.
    VariableShiftThroughOreStep,
    /// The bounds are below the operator's shape.
    BoundsBelowShape,
}

impl fmt::Display for DifferenceRowMapError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ore(error) => write!(f, "invalid step: {error}"),
            Self::Substitution(error) => write!(f, "invalid substitution: {error}"),
            Self::ShiftedOperator => f.write_str("an operator step must be shift-free"),
            Self::DerivativeInBasisConjugation => {
                f.write_str("the conjugation to the power basis needs an operator free of D")
            }
            Self::NegativeVariableShift => f.write_str(
                "the conjugation to the power basis needs nonnegative powers of the shift E",
            ),
            Self::VariableShiftInBasisConjugation => f.write_str(
                "the conjugation to the falling-factorial basis needs an operator free of E",
            ),
            Self::VariableShiftThroughOreStep => f.write_str(
                "a reversal, a coefficient multiplier, an affine map with a non-integer \
                 slope or a substitution needs an operator free of E",
            ),
            Self::BoundsBelowShape => f.write_str("the pullback bounds are below the shape"),
        }
    }
}

impl std::error::Error for DifferenceRowMapError {}

impl From<RowMapError> for DifferenceRowMapError {
    fn from(error: RowMapError) -> Self {
        Self::Ore(error)
    }
}

impl From<SubstitutionError> for DifferenceRowMapError {
    fn from(error: SubstitutionError) -> Self {
        Self::Substitution(error)
    }
}

/// Bounds on an operator `M` that fix the pullback formulas: its largest
/// shift, derivative order, coefficient degree in `x` and range of powers
/// of `E` may be smaller.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct DifferencePullbackBounds {
    /// At least the largest shift.
    pub top_shift: usize,
    /// At least the largest derivative order.
    pub derivative_order: usize,
    /// At most the smallest power of `E`.
    pub lowest_variable_shift: i64,
    /// At least the largest power of `E`.
    pub highest_variable_shift: i64,
    /// At least the degree in `x` of every coefficient.
    pub variable_degree: usize,
}

impl DifferencePullbackBounds {
    /// The exact bounds of `op`.
    pub fn of(op: &DifferentialDifferenceOperator) -> Self {
        Self {
            top_shift: op.max_shift(),
            derivative_order: op.max_derivative(),
            lowest_variable_shift: op.min_variable_shift(),
            highest_variable_shift: op.max_variable_shift(),
            variable_degree: op.max_variable_degree(),
        }
    }

    /// Whether `op` is zero or within the bounds.
    pub fn contains(&self, op: &DifferentialDifferenceOperator) -> bool {
        if op.is_zero() {
            return true;
        }
        let actual = Self::of(op);
        actual.top_shift <= self.top_shift
            && actual.derivative_order <= self.derivative_order
            && actual.lowest_variable_shift >= self.lowest_variable_shift
            && actual.highest_variable_shift <= self.highest_variable_shift
            && actual.variable_degree <= self.variable_degree
    }

    /// The bounds of the components free of `E`.
    pub fn ore_bounds(&self) -> PullbackBounds {
        PullbackBounds {
            top_shift: self.top_shift,
            derivative_order: self.derivative_order,
            variable_degree: self.variable_degree,
        }
    }

    fn with_ore_bounds(&self, bounds: PullbackBounds) -> Self {
        Self {
            top_shift: bounds.top_shift,
            derivative_order: bounds.derivative_order,
            variable_degree: bounds.variable_degree,
            ..*self
        }
    }

    /// Bounds on the pullback through `step` (valid upper bounds, used to
    /// chain several steps).
    fn after(&self, step: &DifferenceRowMapStep) -> Self {
        match step {
            DifferenceRowMapStep::Operator(t) => Self {
                top_shift: self.top_shift,
                derivative_order: self.derivative_order + t.max_derivative(),
                lowest_variable_shift: self.lowest_variable_shift + t.min_variable_shift(),
                highest_variable_shift: self.highest_variable_shift + t.max_variable_shift(),
                variable_degree: self.variable_degree + t.max_variable_degree(),
            },
            // X^b (E - 1)^d with b <= variable degree, d <= derivative order.
            DifferenceRowMapStep::FallingFactorialToPowerBasis => Self {
                top_shift: self.top_shift,
                derivative_order: 0,
                lowest_variable_shift: -(self.variable_degree as i64),
                highest_variable_shift: self.derivative_order as i64,
                variable_degree: self.variable_degree,
            },
            // (x (1 + D))^b (1 + D)^(e + c).
            DifferenceRowMapStep::PowerToFallingFactorialBasis => {
                let c = (-self.lowest_variable_shift).max(0);
                let top = (self.highest_variable_shift + c).max(0) as usize;
                Self {
                    top_shift: self.top_shift,
                    derivative_order: self.variable_degree + top,
                    lowest_variable_shift: 0,
                    highest_variable_shift: 0,
                    variable_degree: self.variable_degree,
                }
            }
            DifferenceRowMapStep::Substitution(substitution) => Self {
                lowest_variable_shift: 0,
                highest_variable_shift: 0,
                ..self.with_ore_bounds(bounds_after_substitution(self.ore_bounds(), substitution))
            },
            DifferenceRowMapStep::Ore(step) => {
                let ore = self.ore_bounds().after(step);
                let (low, high) = match step {
                    RowMapStep::WeylOperator(_) | RowMapStep::Gauge { .. } => {
                        (self.lowest_variable_shift, self.highest_variable_shift)
                    }
                    RowMapStep::Affine { slope, .. } => match integer_slope(slope) {
                        Some(a) => {
                            let (p, q) = (
                                a * self.lowest_variable_shift,
                                a * self.highest_variable_shift,
                            );
                            (p.min(q), p.max(q))
                        }
                        None => (0, 0),
                    },
                    RowMapStep::Reversal { .. } | RowMapStep::CoefficientMultiplier { .. } => {
                        (0, 0)
                    }
                };
                Self {
                    lowest_variable_shift: low,
                    highest_variable_shift: high,
                    ..self.with_ore_bounds(ore)
                }
            }
        }
    }
}

/// The slope as an `i64` when it is an integer.
fn integer_slope(slope: &Rational) -> Option<i64> {
    slope
        .is_integer()
        .then(|| slope.to_integer().to_i64())
        .flatten()
}

impl DifferenceRowMapStep {
    /// `f(x) -> f(x + amount)` (an affine step with slope one).
    pub fn translation(amount: Rational) -> Self {
        Self::Ore(RowMapStep::Affine {
            slope: Rational::one(),
            intercept: amount,
        })
    }

    /// `f -> Delta^order f`, `(Delta f)(x) = f(x + 1) - f(x)`.
    pub fn forward_difference(order: usize) -> Self {
        Self::Operator(DifferentialDifferenceOperator::forward_difference(order))
    }

    /// Check the parameters: an operator step must be shift-free, an
    /// existing step as in [`RowMapStep::validate`].
    pub fn validate(&self) -> Result<(), DifferenceRowMapError> {
        match self {
            Self::Operator(t) => {
                if t.max_shift() != 0 {
                    return Err(DifferenceRowMapError::ShiftedOperator);
                }
            }
            Self::Ore(step) => step.validate()?,
            Self::Substitution(substitution) => substitution.validate()?,
            Self::FallingFactorialToPowerBasis | Self::PowerToFallingFactorialBasis => {}
        }
        Ok(())
    }

    /// The image `G_n` of one row `F_n` (exact; `None` where an existing step
    /// gives none, see [`RowMapStep::apply_row`]).
    pub fn apply_row(&self, row: &[Rational], n: usize) -> Option<Vec<Rational>> {
        Some(match self {
            Self::Operator(t) => t.apply_row(row, n),
            Self::FallingFactorialToPowerBasis => falling_factorial_to_power_basis(row),
            Self::PowerToFallingFactorialBasis => power_to_falling_factorial_basis(row),
            Self::Ore(step) => return step.apply_row(row, n),
            Self::Substitution(substitution) => return substitution.apply_row(row, n),
        })
    }
}

/// Collect `sum_e P_e E^(factor e)` from the components `P_e`.
fn recombine(
    components: impl IntoIterator<Item = (i64, OreOperator)>,
    factor: i64,
) -> DifferentialDifferenceOperator {
    let mut out = DifferentialDifferenceOperator::default();
    for (e, component) in components {
        for ((i, d), c) in component.0 {
            out.add_term((i, d, e * factor), c);
        }
    }
    out
}

/// The powers of `E` occurring in `op`.
fn variable_shifts(op: &DifferentialDifferenceOperator) -> Vec<i64> {
    let mut shifts: Vec<i64> = op.0.keys().map(|(_, _, e)| *e).collect();
    shifts.sort_unstable();
    shifts.dedup();
    shifts
}

/// The pullback of `op` through one existing step (module documentation);
/// multiplies `base_factor` by the step's base factor.
fn pull_back_through_ore_step(
    op: &DifferentialDifferenceOperator,
    step: &RowMapStep,
    bounds: &DifferencePullbackBounds,
    base_factor: &mut Vec<Rational>,
) -> Result<DifferentialDifferenceOperator, DifferenceRowMapError> {
    let ore_bounds = bounds.ore_bounds();
    if let RowMapStep::WeylOperator(t) = step {
        return Ok(op.compose(&DifferentialDifferenceOperator::from_ore_operator(t)));
    }
    // The base factor depends on the bounds only.
    let (_, factor) = pull_back(
        &OreOperator::default(),
        std::slice::from_ref(step),
        ore_bounds,
    )?;
    *base_factor = poly::mul(base_factor, &factor);
    let exponent_factor = match step {
        RowMapStep::Gauge { .. } => Some(1),
        RowMapStep::Affine { slope, .. } => integer_slope(slope),
        _ => None,
    };
    let Some(exponent_factor) = exponent_factor.or(op.is_variable_shift_free().then_some(1)) else {
        return Err(DifferenceRowMapError::VariableShiftThroughOreStep);
    };
    let mut components = Vec::new();
    for e in variable_shifts(op) {
        let component = op.variable_shift_component(e);
        let (pulled, _) = pull_back(&component, std::slice::from_ref(step), ore_bounds)?;
        components.push((e, pulled));
    }
    Ok(recombine(components, exponent_factor))
}

/// The pullback of `op` through `steps` (applied to `F` first to last, so
/// pulled back last to first) with bounds `bounds` on `op`: the operator `P`
/// on the source sequence and the base factor `b(n)` (module
/// documentation).
///
/// # Errors
///
/// [`DifferenceRowMapError::BoundsBelowShape`] when `op` is not within
/// `bounds`, an invalid step, or a step that does not apply to the operator
/// reached at that point (module documentation).
pub fn pull_back_difference(
    op: &DifferentialDifferenceOperator,
    steps: &[DifferenceRowMapStep],
    bounds: DifferencePullbackBounds,
) -> Result<(DifferentialDifferenceOperator, Vec<Rational>), DifferenceRowMapError> {
    if !bounds.contains(op) {
        return Err(DifferenceRowMapError::BoundsBelowShape);
    }
    let mut current = op.clone();
    let mut current_bounds = bounds;
    let mut base_factor = vec![Rational::one()];
    for step in steps.iter().rev() {
        step.validate()?;
        current = match step {
            DifferenceRowMapStep::Operator(t) => current.compose(t),
            DifferenceRowMapStep::FallingFactorialToPowerBasis => {
                let ore = current
                    .to_ore_operator()
                    .ok_or(DifferenceRowMapError::VariableShiftInBasisConjugation)?;
                conjugate_to_falling_factorial_basis(&ore)
            }
            DifferenceRowMapStep::PowerToFallingFactorialBasis => {
                let c = (-current_bounds.lowest_variable_shift).max(0);
                DifferentialDifferenceOperator::from_ore_operator(&conjugate_to_power_basis(
                    &current.variable_shifted(c),
                )?)
            }
            DifferenceRowMapStep::Substitution(substitution) => {
                let ore = current
                    .to_ore_operator()
                    .ok_or(DifferenceRowMapError::VariableShiftThroughOreStep)?;
                DifferentialDifferenceOperator::from_ore_operator(&pull_back_substitution(
                    &ore,
                    substitution,
                    current_bounds.ore_bounds(),
                )?)
            }
            DifferenceRowMapStep::Ore(step) => {
                pull_back_through_ore_step(&current, step, &current_bounds, &mut base_factor)?
            }
        };
        current_bounds = current_bounds.after(step);
    }
    Ok((current, poly::trim(&base_factor)))
}

/// `prod_(t = from..=to) p(n + t)` (one when `from > to`).
fn shifted_product(p: &[Rational], from: i64, to: i64) -> BivariatePolynomial {
    let p = BivariatePolynomial::from_index_polynomial(p);
    (from..=to).fold(BivariatePolynomial::one(), |acc, t| {
        acc.multiply(&p.shift_index(t))
    })
}

/// An annihilator of `G` from an annihilator `op` of `F` (module
/// documentation), or `None` when the step has no push-forward for `op` or
/// is invalid.
pub fn push_forward_difference(
    op: &DifferentialDifferenceOperator,
    step: &DifferenceRowMapStep,
) -> Option<DifferentialDifferenceOperator> {
    step.validate().ok()?;
    match step {
        DifferenceRowMapStep::Operator(t) => {
            let (key, c) = t.0.iter().next().filter(|_| t.term_count() == 1)?;
            let constant = c.as_constant().filter(|c| !c.is_zero());
            (key.0 == 0 && key.1 == 0 && constant.is_some())
                .then(|| op.conjugate_by_variable_shift(key.2))
        }
        DifferenceRowMapStep::FallingFactorialToPowerBasis => {
            let c = (-op.min_variable_shift()).max(0);
            conjugate_to_power_basis(&op.variable_shifted(c))
                .ok()
                .map(|q| DifferentialDifferenceOperator::from_ore_operator(&q))
        }
        DifferenceRowMapStep::PowerToFallingFactorialBasis => {
            Some(conjugate_to_falling_factorial_basis(&op.to_ore_operator()?))
        }
        DifferenceRowMapStep::Substitution(substitution) => {
            push_forward_substitution(&op.to_ore_operator()?, substitution)
                .ok()
                .map(|q| DifferentialDifferenceOperator::from_ore_operator(&q))
        }
        DifferenceRowMapStep::Ore(RowMapStep::WeylOperator(_)) => None,
        DifferenceRowMapStep::Ore(RowMapStep::Gauge {
            numerator,
            denominator,
        }) => {
            let top = op.max_shift() as i64;
            let mut out = DifferentialDifferenceOperator::default();
            for ((i, d, e), c) in &op.0 {
                let i = *i as i64;
                let factor = shifted_product(numerator, i + 1, top).multiply(&shifted_product(
                    denominator,
                    1,
                    i,
                ));
                out.add_term((i as usize, *d, *e), c.multiply(&factor));
            }
            Some(out)
        }
        DifferenceRowMapStep::Ore(step) => {
            if let Some(ore) = op.to_ore_operator() {
                return push_forward(&ore, step)
                    .map(|q| DifferentialDifferenceOperator::from_ore_operator(&q));
            }
            let RowMapStep::Affine { slope, .. } = step else {
                return None;
            };
            let a = integer_slope(slope).filter(|a| a.abs() == 1)?;
            let mut components = Vec::new();
            for e in variable_shifts(op) {
                components.push((e, push_forward(&op.variable_shift_component(e), step)?));
            }
            Some(recombine(components, a))
        }
    }
}
