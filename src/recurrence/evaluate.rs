//! Replaying a recurrence: index shifts, initial-row requirements and exact
//! row generation.

use super::*;
use num_traits::One;
use std::fmt;

impl Recurrence {
    /// Largest recurrence offset used by any term.
    pub fn max_offset(&self) -> usize {
        self.terms.iter().map(|term| term.offset).max().unwrap_or(0)
    }

    /// Number of available initial rows needed for robust JSON generation.
    ///
    /// A denominator recurrence cannot be evaluated at an index where the
    /// evaluated LHS factor is the zero polynomial.  If such singular indices
    /// occur inside the known prefix, keep those rows as initial conditions so
    /// `recurrence-generate` can still reproduce and extend the checked rows.
    pub fn generation_initial_count(&self, first_index: usize, available_rows: usize) -> usize {
        let count = self.max_offset().min(available_rows);
        let mut required = count;
        if let Some(denominator) = &self.denominator {
            for offset in count..available_rows {
                let Some(n) = first_index.checked_add(offset) else {
                    return available_rows;
                };
                if poly_is_zero_rational(&bivar_eval_n(denominator, n)) {
                    required = offset + 1;
                }
            }
        }
        required
    }

    /// Substitute `n + shift` for every occurrence of the recurrence index.
    ///
    /// Polynomial references keep the same offsets.  Alternating factors are
    /// adjusted as `(-1)^(n + shift) = (-1)^shift (-1)^n`.
    pub fn shift_index(&self, shift: i64) -> Self {
        let odd_shift = shift.rem_euclid(2) == 1;
        Self {
            terms: self
                .terms
                .iter()
                .map(|term| {
                    let mut coefficient = term.coeff.shift_index(shift);
                    if odd_shift && term.sign == RecurrenceSign::AlternatingN {
                        coefficient = coefficient.negated();
                    }
                    RecurrenceTerm {
                        offset: term.offset,
                        deriv_order: term.deriv_order,
                        sign: term.sign,
                        coeff: coefficient,
                    }
                })
                .collect(),
            denominator: self
                .denominator
                .as_ref()
                .map(|polynomial| polynomial.shift_index(shift)),
            inhomogeneous: self
                .inhomogeneous
                .as_ref()
                .map(|polynomial| polynomial.shift_index(shift)),
        }
    }
}

/// Errors that can occur while replaying a recurrence to generate rows.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RecurrenceEvaluationError {
    /// More initial rows are needed before the recurrence can be evaluated.
    InsufficientInitialRows {
        /// Number of rows needed, the largest term offset.
        required: usize,
        /// Number of rows supplied.
        provided: usize,
    },
    /// An offset points before the first available initial/generated row.
    OffsetBeforeInitialRow {
        /// Index `n` of the row being evaluated.
        n: usize,
        /// Offset `r` of the term that refers to `P_{n-r}`.
        offset: usize,
        /// Index of the first available row.
        first_index: usize,
    },
    /// The evaluated LHS factor is the zero polynomial.
    ZeroDenominator,
    /// A denominator recurrence produced a rational function rather than a polynomial.
    NonPolynomialQuotient,
    /// The requested row index cannot be represented by `usize`.
    IndexOverflow {
        /// Index of the first available row.
        first_index: usize,
        /// Number of rows already available; `first_index + row_offset` overflowed.
        row_offset: usize,
    },
}

impl fmt::Display for RecurrenceEvaluationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InsufficientInitialRows { required, provided } => write!(
                f,
                "insufficient initial rows: need at least {required}, got {provided}"
            ),
            Self::OffsetBeforeInitialRow {
                n,
                offset,
                first_index,
            } => write!(
                f,
                "cannot evaluate P({n}): offset {offset} refers before first index {first_index}"
            ),
            Self::ZeroDenominator => write!(f, "evaluated denominator is the zero polynomial"),
            Self::NonPolynomialQuotient => write!(
                f,
                "recurrence evaluation did not divide exactly to a polynomial"
            ),
            Self::IndexOverflow {
                first_index,
                row_offset,
            } => write!(
                f,
                "row index overflow: {first_index} + {row_offset} does not fit in usize"
            ),
        }
    }
}

impl std::error::Error for RecurrenceEvaluationError {}

impl Recurrence {
    /// Evaluate the next row after `rows`, where `rows[0]` is `P(first_index)`.
    pub fn evaluate_next_rational(
        &self,
        rows: &[Vec<BigRational>],
        first_index: usize,
    ) -> Result<Vec<BigRational>, RecurrenceEvaluationError> {
        let required = self.max_offset();
        if rows.len() < required {
            return Err(RecurrenceEvaluationError::InsufficientInitialRows {
                required,
                provided: rows.len(),
            });
        }

        let n = first_index.checked_add(rows.len()).ok_or(
            RecurrenceEvaluationError::IndexOverflow {
                first_index,
                row_offset: rows.len(),
            },
        )?;
        let mut rhs = rational_zero_poly();

        for term in &self.terms {
            let Some(reference_index) = n.checked_sub(term.offset) else {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            };
            if reference_index < first_index {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            }
            let row_index = reference_index - first_index;
            let Some(ref_poly) = rows.get(row_index) else {
                return Err(RecurrenceEvaluationError::OffsetBeforeInitialRow {
                    n,
                    offset: term.offset,
                    first_index,
                });
            };

            let deriv = poly_nth_derivative_rational(ref_poly, term.deriv_order);
            let coeff = bivar_eval_n(&term.coeff, n);
            let product = poly_mul_rational(&coeff, &deriv);
            let sign = match term.sign {
                RecurrenceSign::None => BigRational::one(),
                RecurrenceSign::AlternatingN if n % 2 == 1 => {
                    BigRational::from_integer(BigInt::from(-1))
                }
                RecurrenceSign::AlternatingN => BigRational::one(),
            };
            poly_add_scaled_assign(&mut rhs, &product, &sign);
        }

        if let Some(inhomogeneous) = &self.inhomogeneous {
            let inh = bivar_eval_n(inhomogeneous, n);
            poly_add_scaled_assign(&mut rhs, &inh, &BigRational::one());
        }

        if let Some(denominator) = &self.denominator {
            poly_div_exact_rational(&rhs, &bivar_eval_n(denominator, n))
        } else {
            Ok(trim_poly_rational(rhs))
        }
    }

    /// Generate `total_rows` rows, including the supplied initial rows.
    pub fn generate_rows_rational(
        &self,
        initial_rows: &[Vec<BigRational>],
        first_index: usize,
        total_rows: usize,
    ) -> Result<Vec<Vec<BigRational>>, RecurrenceEvaluationError> {
        if total_rows <= initial_rows.len() {
            return Ok(initial_rows[..total_rows].to_vec());
        }
        let required = self.max_offset();
        if initial_rows.len() < required {
            return Err(RecurrenceEvaluationError::InsufficientInitialRows {
                required,
                provided: initial_rows.len(),
            });
        }

        let mut rows = initial_rows.to_vec();
        while rows.len() < total_rows {
            let next = self.evaluate_next_rational(&rows, first_index)?;
            rows.push(next);
        }
        Ok(rows)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::test_support::*;

    #[test]
    fn generation_initial_count_keeps_singular_denominator_rows() {
        let recurrence = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly {
                    coeffs: vec![vec![BigRational::one()]],
                },
            )],
            denominator: Some(BivarPoly {
                coeffs: vec![vec![br(20, 1)], vec![br(-9, 1)], vec![br(1, 1)]],
            }),
            inhomogeneous: None,
        };

        assert_eq!(recurrence.generation_initial_count(1, 7), 5);
        assert_eq!(recurrence.generation_initial_count(1, 4), 4);
    }

    #[test]
    fn recurrence_generation_reports_index_overflow() {
        let recurrence = Recurrence {
            terms: vec![RecurrenceTerm::new(
                1,
                0,
                BivarPoly::from_coefficients_normalized(vec![vec![BigRational::one()]]),
            )],
            denominator: None,
            inhomogeneous: None,
        };
        assert_eq!(
            recurrence.evaluate_next_rational(&[vec![BigRational::one()]], usize::MAX),
            Err(RecurrenceEvaluationError::IndexOverflow {
                first_index: usize::MAX,
                row_offset: 1,
            })
        );
    }
}
