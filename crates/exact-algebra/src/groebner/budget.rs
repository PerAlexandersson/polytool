//! Subtraction observers and the bounded rational computation policy.

use super::{RationalField, SparsePolynomial};
use std::convert::Infallible;

/// Failure of a bounded rational Gröbner computation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum GroebnerError {
    /// No divisor subtraction remains in the shared step budget.
    StepBudgetExhausted,
    /// A coefficient numerator or denominator exceeded the bit limit.
    CoefficientBitBudgetExceeded,
    /// Appending an S-pair remainder exceeded the basis size limit.
    BasisSizeBudgetExceeded,
    /// Input generators have different variable counts.
    IncompatibleVariableCounts,
}
impl std::fmt::Display for GroebnerError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            Self::StepBudgetExhausted => "the Gröbner step budget was exhausted",
            Self::CoefficientBitBudgetExceeded => "a Gröbner coefficient exceeded the bit budget",
            Self::BasisSizeBudgetExceeded => "the Gröbner basis size budget was exhausted",
            Self::IncompatibleVariableCounts => "generators have different variable counts",
        })
    }
}
impl std::error::Error for GroebnerError {}

/// Mutable limits for rational lexicographic Gröbner computations.
///
/// A step is spent after each successful divisor subtraction, including input
/// reduction and tail interreduction. Moving an irreducible term is free. The
/// coefficient check follows the step debit and inspects only the current
/// dividend. Initial generators and monic scaling are not size checked. The
/// basis limit is checked only after appending a new S-pair remainder. Reusing
/// this value across calls shares the remaining steps, including after errors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GroebnerBudget {
    /// Remaining successful divisor subtractions; the last available step succeeds.
    pub remaining_steps: usize,
    /// Maximum bit length of a coefficient numerator or denominator after subtraction.
    pub max_coefficient_bits: u64,
    /// Maximum basis length after appending an S-pair remainder.
    pub max_basis_elements: usize,
}

pub(super) trait ReductionObserver<P> {
    type Error;
    fn after_subtraction(&mut self, polynomial: &P) -> Result<(), Self::Error>;
    fn after_basis_append(&mut self, _length: usize) -> Result<(), Self::Error> {
        Ok(())
    }
}

pub(super) struct Unlimited;
impl<P> ReductionObserver<P> for Unlimited {
    type Error = Infallible;
    #[inline]
    fn after_subtraction(&mut self, _: &P) -> Result<(), Infallible> {
        Ok(())
    }
}
impl<P: SparsePolynomial<Field = RationalField>> ReductionObserver<P> for GroebnerBudget {
    type Error = GroebnerError;
    fn after_subtraction(&mut self, polynomial: &P) -> Result<(), GroebnerError> {
        self.remaining_steps = self
            .remaining_steps
            .checked_sub(1)
            .ok_or(GroebnerError::StepBudgetExhausted)?;
        if polynomial
            .terms()
            .values()
            .any(|c| c.numer().bits().max(c.denom().bits()) > self.max_coefficient_bits)
        {
            return Err(GroebnerError::CoefficientBitBudgetExceeded);
        }
        Ok(())
    }
    fn after_basis_append(&mut self, length: usize) -> Result<(), GroebnerError> {
        if length > self.max_basis_elements {
            Err(GroebnerError::BasisSizeBudgetExceeded)
        } else {
            Ok(())
        }
    }
}
