//! First-order vector recurrences over the Weyl algebra: the public types,
//! error enums and exact evaluation of fitted systems.

use super::*;
use num_traits::One;
use num_traits::Zero;
use std::fmt;

/// A first-order affine vector recurrence
/// `F_{n+1} = matrix(n, x, D_x) F_n + forcing(n, x)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorRecurrence {
    /// Row-major operator matrix: target component, then source component.
    pub matrix: Vec<Vec<WeylOperator>>,
    /// Optional additive polynomial vector. `None` means a homogeneous system.
    pub forcing: Option<Vec<BivarPoly>>,
}

impl VectorRecurrence {
    /// Number of source and target components in the square system.
    pub fn dimension(&self) -> usize {
        self.matrix.len()
    }
}

/// Fixed search bounds for a first-order vector Weyl recurrence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorRecurrenceOptions {
    /// Maximum degree in `x` of every matrix-entry coefficient.
    pub var_deg: usize,
    /// Maximum degree in the source index `n` of every matrix entry.
    pub idx_deg: usize,
    /// Maximum derivative order in every matrix entry.
    pub diff_deg: usize,
    /// If false, fit an additive polynomial forcing vector.
    pub homogeneous: bool,
    /// Maximum degree in `x` of each forcing component.
    pub forcing_var_deg: usize,
    /// Maximum degree in `n` of each forcing component.
    pub forcing_idx_deg: usize,
    /// Number of final complete transitions excluded from fitting and checked
    /// independently.
    pub held_out_transitions: usize,
    /// Reject non-identifiable rows instead of returning an arbitrary
    /// particular solution with free parameters set to zero.
    pub require_unique: bool,
}

/// Search bounds for a vector recurrence with one common polynomial factor on
/// the left hand side.  The fitted identity is
/// `q(n,x) F_(n+1) = M(n,x,D_x) F_n + G(n,x)`.
///
/// This is deliberately a separate options type: adding fields to
/// [`VectorRecurrenceOptions`] would break downstream struct literals using
/// the original, normalized (`q = 1`) API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorRecurrenceDenominatorOptions {
    /// Maximum degree in `x` of every matrix-entry coefficient.
    pub var_deg: usize,
    /// Maximum degree in the source index `n` of every matrix-entry coefficient.
    pub idx_deg: usize,
    /// Maximum derivative order in every matrix entry.
    pub diff_deg: usize,
    /// If false, fit an additive polynomial forcing vector.
    pub homogeneous: bool,
    /// Maximum degree in `x` of each forcing component.
    pub forcing_var_deg: usize,
    /// Maximum degree in `n` of each forcing component.
    pub forcing_idx_deg: usize,
    /// Number of final complete transitions `F_n -> F_{n+1}` excluded from fitting and
    /// checked independently.
    pub held_out_transitions: usize,
    /// Reject non-identifiable rows instead of returning an arbitrary particular
    /// solution with free parameters set to zero.
    pub require_unique: bool,
    /// Maximum degree in `x` of the common factor `q`.
    pub denominator_var_deg: usize,
    /// Maximum degree in the source index `n` of `q`.
    pub denominator_idx_deg: usize,
}

/// Deprecated alias of [`VectorRecurrenceDenominatorOptions`].
#[deprecated(since = "0.6.0", note = "use `VectorRecurrenceDenominatorOptions`")]
pub type SharedDenominatorVectorRecurrenceOptions = VectorRecurrenceDenominatorOptions;

impl Default for VectorRecurrenceDenominatorOptions {
    fn default() -> Self {
        Self {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            homogeneous: true,
            forcing_var_deg: 1,
            forcing_idx_deg: 1,
            held_out_transitions: 2,
            require_unique: true,
            denominator_var_deg: 0,
            denominator_idx_deg: 0,
        }
    }
}

impl Default for VectorRecurrenceOptions {
    fn default() -> Self {
        Self {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            homogeneous: true,
            forcing_var_deg: 1,
            forcing_idx_deg: 1,
            held_out_transitions: 2,
            require_unique: true,
        }
    }
}

/// Rank information for one independently fitted output row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorRecurrenceRowDiagnostics {
    /// Target component fitted by this linear system.
    pub output_component: usize,
    /// Number of coefficient equations used for fitting.
    pub equations: usize,
    /// Number of unknown normal-form coefficients.
    pub unknowns: usize,
    /// Exact rank of the coefficient matrix.
    pub rank: usize,
    /// `unknowns - rank`; zero means the row is identifiable in this model.
    pub nullity: usize,
}

/// A fitted recurrence and the diagnostics used to accept it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VectorRecurrenceFit {
    /// The fitted recurrence. For a companion fit the matrix acts on the stacked state
    /// `(F_{n-lag+1}, ..., F_n)`.
    pub recurrence: VectorRecurrence,
    /// Rank diagnostics of each fitted output row. A companion fit solves only its final
    /// block row, so only that block's components are listed.
    pub rows: Vec<VectorRecurrenceRowDiagnostics>,
    /// Number of transitions `F_n -> F_{n+1}` used for fitting.
    pub fit_transitions: usize,
    /// Number of final transitions reserved for independent exact verification.
    pub held_out_transitions: usize,
}

/// A vector recurrence with a common polynomial denominator on its left.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedDenominatorVectorRecurrence {
    /// Common left factor `q(n, x)`.
    pub denominator: BivarPoly,
    /// Row-major operator matrix: target component, then source component.
    pub matrix: Vec<Vec<WeylOperator>>,
    /// Optional additive polynomial vector, one entry per target component. `None` means a
    /// homogeneous system.
    pub forcing: Option<Vec<BivarPoly>>,
}

/// Deprecated alias of [`SharedDenominatorVectorRecurrence`].
#[deprecated(since = "0.6.0", note = "use `SharedDenominatorVectorRecurrence`")]
pub type VectorRecurrenceWithDenominator = SharedDenominatorVectorRecurrence;

impl SharedDenominatorVectorRecurrence {
    /// Number of target components (rows of `matrix`); for a square system this is also
    /// the number of source components.
    pub fn dimension(&self) -> usize {
        self.matrix.len()
    }

    /// Evaluate the recurrence exactly, dividing every numerator component by
    /// the common evaluated denominator.
    pub fn evaluate_next_rational(
        &self,
        state: &[Vec<BigRational>],
        source_index: usize,
    ) -> Result<Vec<Vec<BigRational>>, SharedDenominatorVectorRecurrenceEvaluationError> {
        let dimension = self.dimension();
        if self.matrix.iter().any(|row| row.len() != dimension) {
            return Err(VectorRecurrenceEvaluationError::NonSquareMatrix.into());
        }
        evaluate_shared_denominator_rectangular_rational(self, state, source_index)
    }

    pub(in crate::recurrence) fn holds_rectangular_identity_rational(
        &self,
        state: &[Vec<BigRational>],
        next_state: &[Vec<BigRational>],
        source_index: usize,
    ) -> bool {
        let Ok(numerator) = evaluate_rectangular_vector_map_rational(
            &self.matrix,
            self.forcing.as_deref(),
            state,
            source_index,
        ) else {
            return false;
        };
        let denominator = bivar_eval_n(&self.denominator, source_index);
        numerator.len() == next_state.len()
            && numerator.iter().zip(next_state).all(|(lhs, target)| {
                trim_poly_rational(poly_mul_rational(&denominator, target))
                    == trim_poly_rational(lhs.clone())
            })
    }

    /// Check one exact square-system transition.
    pub fn holds_transition_rational(
        &self,
        state: &[Vec<BigRational>],
        next_state: &[Vec<BigRational>],
        source_index: usize,
    ) -> bool {
        self.evaluate_next_rational(state, source_index)
            .is_ok_and(|actual| vector_states_equal(&actual, next_state))
    }
}

fn evaluate_shared_denominator_rectangular_rational(
    recurrence: &SharedDenominatorVectorRecurrence,
    state: &[Vec<BigRational>],
    source_index: usize,
) -> Result<Vec<Vec<BigRational>>, SharedDenominatorVectorRecurrenceEvaluationError> {
    let numerator = evaluate_rectangular_vector_map_rational(
        &recurrence.matrix,
        recurrence.forcing.as_deref(),
        state,
        source_index,
    )?;
    let denominator = bivar_eval_n(&recurrence.denominator, source_index);
    if poly_is_zero_rational(&denominator) {
        return Err(SharedDenominatorVectorRecurrenceEvaluationError::ZeroDenominator);
    }
    numerator
        .iter()
        .map(|component| {
            poly_div_exact_rational(component, &denominator).map_err(|error| match error {
                RecurrenceEvaluationError::ZeroDenominator => {
                    SharedDenominatorVectorRecurrenceEvaluationError::ZeroDenominator
                }
                RecurrenceEvaluationError::NonPolynomialQuotient => {
                    SharedDenominatorVectorRecurrenceEvaluationError::NonPolynomialQuotient
                }
                _ => SharedDenominatorVectorRecurrenceEvaluationError::NonPolynomialQuotient,
            })
        })
        .collect()
}

/// Result of a shared-denominator vector recurrence search.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedDenominatorVectorRecurrenceFit {
    /// The fitted recurrence.
    pub recurrence: SharedDenominatorVectorRecurrence,
    /// Rank and equation counts of the joint linear system.
    pub diagnostics: SharedDenominatorVectorRecurrenceDiagnostics,
    /// Number of transitions `F_n -> F_{n+1}` used for fitting.
    pub fit_transitions: usize,
    /// Number of final transitions reserved for independent exact verification.
    pub held_out_transitions: usize,
    /// Rank after fixing the canonical denominator pivot to one.
    pub rank: usize,
    /// Nullity after fixing the canonical denominator pivot.
    pub nullity: usize,
    /// Exponents `(i, j)` of the denominator monomial `n^i x^j` that was normalized to one;
    /// all earlier monomials in the ordering `(n^0 x^0, n^0 x^1, ..., n^1 x^0, ...)` are zero.
    pub denominator_pivot: (usize, usize),
}

/// Diagnostics for the joint shared-denominator solve.  Unlike
/// [`VectorRecurrenceRowDiagnostics`], these values describe one global linear
/// system, not independently fitted output rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SharedDenominatorVectorRecurrenceDiagnostics {
    /// Total number of coefficient equations, the sum of `equations_by_output`.
    pub total_equations: usize,
    /// Number of coefficient equations contributed by each output component, in the order
    /// of `output_components`.
    pub equations_by_output: Vec<usize>,
    /// Number of unknown coefficients after the denominator pivot is fixed: later
    /// denominator monomials plus the matrix and forcing coefficients of every output row.
    pub unknowns: usize,
    /// Exact rank of the coefficient matrix of the joint system.
    pub rank: usize,
    /// `unknowns - rank`; zero means the model is identifiable once the pivot is fixed.
    pub nullity: usize,
    /// Component indices of the fitted output rows. A companion fit reports the final
    /// block row, so these indices start at `(lag - 1)` times the component count.
    pub output_components: Vec<usize>,
}

/// Deprecated alias of [`SharedDenominatorVectorRecurrenceFit`].
#[deprecated(since = "0.6.0", note = "use `SharedDenominatorVectorRecurrenceFit`")]
pub type VectorRecurrenceWithDenominatorFit = SharedDenominatorVectorRecurrenceFit;

/// Why a fixed-bound vector recurrence search could not return a model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VectorRecurrenceFitError {
    /// The state list is empty or its first state has no components.
    EmptyStateVector,
    /// A companion fit was requested with `lag == 0`.
    InvalidLag,
    /// There is no fitting transition left after reserving the held-out ones.
    NotEnoughTransitions {
        /// Number of available transitions (`states.len() - 1`, or `states.len() - lag` for a
        /// companion fit).
        transitions: usize,
        /// Requested number of held-out transitions.
        held_out: usize,
    },
    /// A state has a different number of components from the first state.
    InconsistentStateDimension {
        /// Zero-based position of the offending state in the input list.
        state: usize,
        /// Number of components of the first state.
        expected: usize,
        /// Number of components of the offending state.
        found: usize,
    },
    /// The linear system has no solution in the requested search space.
    NoSolution {
        /// Zero-based index of the output component whose system failed.
        output_component: usize,
    },
    /// Unique identification was required, but a row has free parameters.
    Underdetermined {
        /// Zero-based index of the output component whose system is underdetermined.
        output_component: usize,
        /// Exact rank of that system.
        rank: usize,
        /// Number of unknowns of that system.
        unknowns: usize,
    },
    /// The fitted recurrence does not reproduce a held-out transition.
    HeldOutVerificationFailed {
        /// Source index `n` of the first failing held-out transition `F_n -> F_{n+1}`.
        transition: usize,
    },
    /// A source index `first_index + offset` does not fit in `usize`.
    IndexOverflow {
        /// Index of the first state.
        first_index: usize,
        /// Zero-based position of the state whose index overflowed.
        offset: usize,
    },
}

impl fmt::Display for VectorRecurrenceFitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyStateVector => write!(f, "vector states must have at least one component"),
            Self::InvalidLag => write!(f, "companion lag must be at least one"),
            Self::NotEnoughTransitions {
                transitions,
                held_out,
            } => write!(
                f,
                "need at least one fitting transition in addition to {held_out} held-out \
                 transitions; got {transitions} total"
            ),
            Self::InconsistentStateDimension {
                state,
                expected,
                found,
            } => write!(
                f,
                "state {state} has {found} components; expected {expected}"
            ),
            Self::NoSolution { output_component } => write!(
                f,
                "no recurrence in the requested search space for output component \
                 {output_component}"
            ),
            Self::Underdetermined {
                output_component,
                rank,
                unknowns,
            } => write!(
                f,
                "output component {output_component} is not identifiable: rank {rank} for \
                 {unknowns} unknowns"
            ),
            Self::HeldOutVerificationFailed { transition } => write!(
                f,
                "fitted recurrence failed exact held-out transition {transition}"
            ),
            Self::IndexOverflow {
                first_index,
                offset,
            } => write!(
                f,
                "source index overflow: {first_index} + {offset} does not fit in usize"
            ),
        }
    }
}

impl std::error::Error for VectorRecurrenceFitError {}

/// Why a shared-denominator vector recurrence search could not return a model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SharedDenominatorVectorRecurrenceFitError {
    /// A failure shared with the denominator-free fitter.
    Base(VectorRecurrenceFitError),
    /// No nonzero denominator in the requested degrees gives a usable recurrence, that is,
    /// one consistent with the data and nonzero at every source index.
    NoSharedDenominator,
    /// The normalized system has free parameters although unique identification was required.
    SharedDenominatorUnderdetermined {
        /// Exact rank of the normalized system (denominator pivot fixed to one).
        rank: usize,
        /// Number of unknowns of the normalized system.
        unknowns: usize,
    },
    /// The fitted recurrence does not reproduce a held-out transition.
    HeldOutVerificationFailed {
        /// Source index `n` of the first failing held-out transition `F_n -> F_{n+1}`.
        transition: usize,
    },
}

impl From<VectorRecurrenceFitError> for SharedDenominatorVectorRecurrenceFitError {
    fn from(error: VectorRecurrenceFitError) -> Self {
        Self::Base(error)
    }
}

impl fmt::Display for SharedDenominatorVectorRecurrenceFitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Base(error) => error.fmt(f),
            Self::NoSharedDenominator => write!(
                f,
                "no nonzero shared denominator recurrence in the requested search space"
            ),
            Self::SharedDenominatorUnderdetermined { rank, unknowns } => write!(
                f,
                "shared-denominator recurrence is not identifiable: rank {rank} for \
                 {unknowns} unknowns"
            ),
            Self::HeldOutVerificationFailed { transition } => write!(
                f,
                "shared-denominator recurrence failed exact held-out transition {transition}"
            ),
        }
    }
}

impl std::error::Error for SharedDenominatorVectorRecurrenceFitError {}

/// Errors from applying a vector recurrence to one state.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VectorRecurrenceEvaluationError {
    /// Some matrix row does not have as many operators as there are rows.
    NonSquareMatrix,
    /// The state has the wrong number of components.
    StateDimension {
        /// Number of components the matrix acts on.
        expected: usize,
        /// Number of components supplied.
        found: usize,
    },
    /// The forcing vector has the wrong number of components.
    ForcingDimension {
        /// Number of matrix rows.
        expected: usize,
        /// Number of forcing components supplied.
        found: usize,
    },
}

impl fmt::Display for VectorRecurrenceEvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonSquareMatrix => write!(f, "vector recurrence matrix must be square"),
            Self::StateDimension { expected, found } => {
                write!(f, "state has {found} components; expected {expected}")
            }
            Self::ForcingDimension { expected, found } => {
                write!(f, "forcing has {found} components; expected {expected}")
            }
        }
    }
}

impl std::error::Error for VectorRecurrenceEvaluationError {}

/// Errors from applying a recurrence with a shared polynomial denominator.
/// The dimension and forcing errors retain the unchanged base vector error
/// type, while denominator-specific failures are represented separately.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SharedDenominatorVectorRecurrenceEvaluationError {
    /// Dimension or forcing error shared with the denominator-free evaluator.
    Base(VectorRecurrenceEvaluationError),
    /// The shared denominator evaluates to the zero polynomial at the source index.
    ZeroDenominator,
    /// Some numerator component is not exactly divisible by the evaluated denominator.
    NonPolynomialQuotient,
}

impl From<VectorRecurrenceEvaluationError> for SharedDenominatorVectorRecurrenceEvaluationError {
    fn from(error: VectorRecurrenceEvaluationError) -> Self {
        Self::Base(error)
    }
}

impl fmt::Display for SharedDenominatorVectorRecurrenceEvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Base(error) => error.fmt(f),
            Self::ZeroDenominator => {
                write!(f, "evaluated shared denominator is the zero polynomial")
            }
            Self::NonPolynomialQuotient => write!(
                f,
                "shared-denominator recurrence did not divide exactly to a polynomial"
            ),
        }
    }
}

impl std::error::Error for SharedDenominatorVectorRecurrenceEvaluationError {}

pub(in crate::recurrence) fn validate_vector_states<C>(
    states: &[Vec<Vec<C>>],
) -> Result<usize, VectorRecurrenceFitError> {
    let dimension = states.first().map_or(0, Vec::len);
    if dimension == 0 {
        return Err(VectorRecurrenceFitError::EmptyStateVector);
    }
    for (state, components) in states.iter().enumerate() {
        if components.len() != dimension {
            return Err(VectorRecurrenceFitError::InconsistentStateDimension {
                state,
                expected: dimension,
                found: components.len(),
            });
        }
    }
    Ok(dimension)
}

pub(in crate::recurrence) fn i64_vector_states_to_rational(
    states: &[Vec<Vec<i64>>],
) -> Vec<Vec<Vec<BigRational>>> {
    states
        .iter()
        .map(|state| {
            state
                .iter()
                .map(|poly| i64_poly_to_rational(poly))
                .collect()
        })
        .collect()
}

pub(in crate::recurrence) fn zero_bivar() -> BivarPoly {
    BivarPoly {
        coeffs: vec![vec![BigRational::zero()]],
    }
}

pub(in crate::recurrence) fn extract_vector_bivar(
    solution: &[BigRational],
    start: usize,
    idx_deg: usize,
    var_deg: usize,
) -> BivarPoly {
    let width = var_deg + 1;
    BivarPoly {
        coeffs: (0..=idx_deg)
            .map(|i| {
                (0..=var_deg)
                    .map(|j| solution[start + i * width + j].clone())
                    .collect()
            })
            .collect(),
    }
}

fn apply_weyl_operator_rational(
    operator: &WeylOperator,
    source: &[BigRational],
    source_index: usize,
) -> Vec<BigRational> {
    let mut result = rational_zero_poly();
    for (deriv_order, coeff) in operator.coefficients.iter().enumerate() {
        if coeff.is_zero() {
            continue;
        }
        let derivative = poly_nth_derivative_rational(source, deriv_order);
        let evaluated_coeff = bivar_eval_n(coeff, source_index);
        let term = poly_mul_rational(&evaluated_coeff, &derivative);
        poly_add_scaled_assign(&mut result, &term, &BigRational::one());
    }
    trim_poly_rational(result)
}

pub(in crate::recurrence) fn evaluate_rectangular_vector_map_rational(
    matrix: &[Vec<WeylOperator>],
    forcing: Option<&[BivarPoly]>,
    state: &[Vec<BigRational>],
    source_index: usize,
) -> Result<Vec<Vec<BigRational>>, VectorRecurrenceEvaluationError> {
    if matrix.iter().any(|row| row.len() != state.len()) {
        return Err(VectorRecurrenceEvaluationError::StateDimension {
            expected: matrix.first().map_or(0, Vec::len),
            found: state.len(),
        });
    }
    if let Some(forcing) = forcing {
        if forcing.len() != matrix.len() {
            return Err(VectorRecurrenceEvaluationError::ForcingDimension {
                expected: matrix.len(),
                found: forcing.len(),
            });
        }
    }

    let mut result = Vec::with_capacity(matrix.len());
    for (output, row) in matrix.iter().enumerate() {
        let mut polynomial = rational_zero_poly();
        for (operator, source) in row.iter().zip(state) {
            let term = apply_weyl_operator_rational(operator, source, source_index);
            poly_add_scaled_assign(&mut polynomial, &term, &BigRational::one());
        }
        if let Some(forcing) = forcing {
            let term = bivar_eval_n(&forcing[output], source_index);
            poly_add_scaled_assign(&mut polynomial, &term, &BigRational::one());
        }
        result.push(trim_poly_rational(polynomial));
    }
    Ok(result)
}

impl VectorRecurrence {
    /// Evaluate `F_{n+1}` from `F_n`, with matrix coefficients evaluated at
    /// the supplied source index `n`.
    pub fn evaluate_next_rational(
        &self,
        state: &[Vec<BigRational>],
        source_index: usize,
    ) -> Result<Vec<Vec<BigRational>>, VectorRecurrenceEvaluationError> {
        let dimension = self.dimension();
        if self.matrix.iter().any(|row| row.len() != dimension) {
            return Err(VectorRecurrenceEvaluationError::NonSquareMatrix);
        }
        if state.len() != dimension {
            return Err(VectorRecurrenceEvaluationError::StateDimension {
                expected: dimension,
                found: state.len(),
            });
        }
        evaluate_rectangular_vector_map_rational(
            &self.matrix,
            self.forcing.as_deref(),
            state,
            source_index,
        )
    }

    /// Check one exact transition `F_n -> F_{n+1}`.
    pub fn holds_transition_rational(
        &self,
        state: &[Vec<BigRational>],
        next_state: &[Vec<BigRational>],
        source_index: usize,
    ) -> bool {
        self.evaluate_next_rational(state, source_index)
            .is_ok_and(|actual| vector_states_equal(&actual, next_state))
    }
}

pub(in crate::recurrence) fn vector_states_equal(
    lhs: &[Vec<BigRational>],
    rhs: &[Vec<BigRational>],
) -> bool {
    lhs.len() == rhs.len()
        && lhs.iter().zip(rhs).all(|(left, right)| {
            trim_poly_rational(left.clone()) == trim_poly_rational(right.clone())
        })
}
