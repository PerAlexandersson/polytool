//! Signed Euclidean/Sturm continued-fraction certificates.
//!
//! This module checks a finite certificate of the form
//!
//! ```text
//! P_(i-1) = q_i P_i - P_(i+1),
//! q_i(x) = alpha_i x + beta_i,  alpha_i,beta_i > 0,
//! ```
//!
//! with positive-coefficient signed remainders.  Such a certificate is useful
//! for turning exact Euclidean data for a Hermite--Biehler pair into a Sturm
//! chain with negative roots and hence interlacing.
//!
//! Entry points: [`check_sturm_continued_fraction_certificate`] (default
//! options) and [`check_sturm_continued_fraction_certificate_with_options`].
//! Inputs are exact `BigRational` coefficient slices in ascending degree order;
//! the check runs exact Euclidean division, never panics on well-formed input,
//! and reports the first violated condition as a [`SturmContinuedFractionError`]
//! naming the (1-based) division step.  Here `P_0` is the first and `P_1` the
//! second polynomial passed in.  For general interlacing tests use
//! [`crate::check_weak_interlacing`]; for the Sturm chain itself see
//! [`crate::sturm`].

use crate::vec_poly::{coeff, degree, is_zero, neg, scale, sub, trim};
use num_rational::BigRational;
use num_traits::{Signed, Zero};
use std::error::Error;
use std::fmt;

/// Options for checking a signed Euclidean/Sturm-CF certificate.
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub struct SturmContinuedFractionOptions {
    /// Require every signed remainder to have strictly positive coefficients.
    pub require_positive_remainders: bool,
    /// Require the roots `-beta_i/alpha_i` of the linear quotients to be
    /// strictly increasing with `i`.
    pub require_increasing_quotient_roots: bool,
}

impl Default for SturmContinuedFractionOptions {
    fn default() -> Self {
        Self {
            require_positive_remainders: true,
            require_increasing_quotient_roots: false,
        }
    }
}

/// Summary returned by a successful signed Euclidean/Sturm-CF check.
#[derive(Debug, Clone, Eq, PartialEq)]
pub struct SturmContinuedFractionSummary {
    /// Degree of `P_0`.
    pub degree_left: usize,
    /// Degree of `P_1`.
    pub degree_right: usize,
    /// Number of Euclidean quotients checked.
    pub steps: usize,
    /// Largest bit size of any quotient numerator.
    pub max_quotient_numerator_bits: u64,
    /// Largest bit size of any quotient denominator.
    pub max_quotient_denominator_bits: u64,
    /// Roots `-beta_i/alpha_i` of the linear quotients, if requested.
    pub quotient_roots: Option<Vec<BigRational>>,
}

/// Error from a failed signed Euclidean/Sturm-CF check.
#[derive(Debug, Clone, Eq, PartialEq)]
pub enum SturmContinuedFractionError {
    /// `P_0` or `P_1` has a coefficient that is not strictly positive (only
    /// checked when positive remainders are required).
    InitialPolynomialHasNonpositiveCoefficient {
        /// Which input polynomial: `0` for `P_0`, `1` for `P_1`.
        index: usize,
    },
    /// The Euclidean quotient is not of degree exactly 1.
    QuotientDegree {
        /// 1-based division step.
        step: usize,
        /// Actual quotient degree (`None` if the quotient is zero).
        degree: Option<usize>,
    },
    /// The linear quotient `alpha x + beta` has `alpha <= 0` or `beta <= 0`.
    QuotientNotPositiveLinear {
        /// 1-based division step.
        step: usize,
    },
    /// The roots `-beta_i/alpha_i` of consecutive quotients are not strictly
    /// increasing.
    QuotientRootsNotIncreasing {
        /// 1-based division step of the later quotient.
        step: usize,
        /// Root of the previous quotient.
        previous: Box<BigRational>,
        /// Root of this quotient (`<= previous`).
        current: Box<BigRational>,
    },
    /// The signed remainder `P_(i+1)` has a coefficient that is not strictly
    /// positive.
    SignedRemainderHasNonpositiveCoefficient {
        /// 1-based division step.
        step: usize,
    },
    /// The degree of the signed remainder is not exactly one less than that of
    /// the divisor.
    SignedRemainderDegreeDrop {
        /// 1-based division step.
        step: usize,
        /// Degree of the divisor `P_i`.
        previous_degree: usize,
        /// Degree of the signed remainder `P_(i+1)`.
        next_degree: usize,
    },
    /// The divisor is the zero polynomial.
    ZeroDivisor {
        /// 1-based division step.
        step: usize,
    },
    /// The signed remainder is unexpectedly the zero polynomial.
    ZeroSignedRemainder {
        /// 1-based division step.
        step: usize,
    },
}

impl fmt::Display for SturmContinuedFractionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InitialPolynomialHasNonpositiveCoefficient { index } => {
                write!(f, "initial polynomial P_{index} has a nonpositive coefficient")
            }
            Self::QuotientDegree { step, degree } => {
                write!(f, "step {step}: quotient has degree {degree:?}, expected 1")
            }
            Self::QuotientNotPositiveLinear { step } => {
                write!(f, "step {step}: quotient is not positive linear")
            }
            Self::QuotientRootsNotIncreasing {
                step,
                previous,
                current,
            } => write!(
                f,
                "step {step}: quotient roots are not strictly increasing: {previous} >= {current}"
            ),
            Self::SignedRemainderHasNonpositiveCoefficient { step } => write!(
                f,
                "step {step}: signed remainder has a nonpositive coefficient"
            ),
            Self::SignedRemainderDegreeDrop {
                step,
                previous_degree,
                next_degree,
            } => write!(
                f,
                "step {step}: signed remainder degree {next_degree} does not drop from {previous_degree} by one"
            ),
            Self::ZeroDivisor { step } => write!(f, "step {step}: zero divisor"),
            Self::ZeroSignedRemainder { step } => write!(f, "step {step}: zero signed remainder"),
        }
    }
}

impl Error for SturmContinuedFractionError {}

/// Check the full signed Euclidean/Sturm-CF certificate for `(P_0,P_1)`.
///
/// Coefficients are in ascending degree order.
pub fn check_sturm_continued_fraction_certificate(
    left: &[BigRational],
    right: &[BigRational],
) -> Result<SturmContinuedFractionSummary, SturmContinuedFractionError> {
    check_sturm_continued_fraction_certificate_with_options(
        left,
        right,
        SturmContinuedFractionOptions::default(),
    )
}

/// Check a signed Euclidean/Sturm-CF certificate for `(P_0,P_1)` with options.
///
/// Coefficients are in ascending degree order.  The check uses exact rational
/// Euclidean division and the sign convention
/// `P_(i-1)=q_i P_i-P_(i+1)`.
pub fn check_sturm_continued_fraction_certificate_with_options(
    left: &[BigRational],
    right: &[BigRational],
    options: SturmContinuedFractionOptions,
) -> Result<SturmContinuedFractionSummary, SturmContinuedFractionError> {
    let mut previous = trim(left);
    let mut current = trim(right);
    let mut steps = 0usize;
    let mut max_numerator_bits = 0u64;
    let mut max_denominator_bits = 0u64;
    let mut previous_quotient_root: Option<BigRational> = None;
    let mut quotient_roots = Vec::new();

    if options.require_positive_remainders {
        for (index, polynomial) in [&previous, &current].into_iter().enumerate() {
            if !all_coefficients_positive(polynomial) {
                return Err(
                    SturmContinuedFractionError::InitialPolynomialHasNonpositiveCoefficient {
                        index,
                    },
                );
            }
        }
    }

    while !is_zero(&current) {
        let (quotient, remainder) = div_rem(&previous, &current);
        steps += 1;
        if degree(&quotient) != Some(1) {
            return Err(SturmContinuedFractionError::QuotientDegree {
                step: steps,
                degree: degree(&quotient),
            });
        }

        let intercept = coeff(&quotient, 0);
        let slope = coeff(&quotient, 1);
        if !slope.is_positive() || !intercept.is_positive() {
            return Err(SturmContinuedFractionError::QuotientNotPositiveLinear { step: steps });
        }
        for coeff in [&slope, &intercept] {
            max_numerator_bits = max_numerator_bits.max(coeff.numer().bits());
            max_denominator_bits = max_denominator_bits.max(coeff.denom().bits());
        }

        if options.require_increasing_quotient_roots {
            let root = -intercept / slope;
            if let Some(previous_root) = &previous_quotient_root {
                if previous_root >= &root {
                    return Err(SturmContinuedFractionError::QuotientRootsNotIncreasing {
                        step: steps,
                        previous: Box::new(previous_root.clone()),
                        current: Box::new(root),
                    });
                }
            }
            previous_quotient_root = Some(root.clone());
            quotient_roots.push(root);
        }

        if is_zero(&remainder) {
            break;
        }
        let next = neg(&remainder);
        if options.require_positive_remainders && !all_coefficients_positive(&next) {
            return Err(
                SturmContinuedFractionError::SignedRemainderHasNonpositiveCoefficient {
                    step: steps,
                },
            );
        }

        let current_degree =
            degree(&current).ok_or(SturmContinuedFractionError::ZeroDivisor { step: steps })?;
        let next_degree = degree(&next)
            .ok_or(SturmContinuedFractionError::ZeroSignedRemainder { step: steps })?;
        if next_degree + 1 != current_degree {
            return Err(SturmContinuedFractionError::SignedRemainderDegreeDrop {
                step: steps,
                previous_degree: current_degree,
                next_degree,
            });
        }

        previous = current;
        current = next;
    }

    Ok(SturmContinuedFractionSummary {
        degree_left: degree(left).unwrap_or(0),
        degree_right: degree(right).unwrap_or(0),
        steps,
        max_quotient_numerator_bits: max_numerator_bits,
        max_quotient_denominator_bits: max_denominator_bits,
        quotient_roots: options
            .require_increasing_quotient_roots
            .then_some(quotient_roots),
    })
}

fn div_rem(a: &[BigRational], b: &[BigRational]) -> (Vec<BigRational>, Vec<BigRational>) {
    assert!(!is_zero(b), "division by zero polynomial");
    let mut remainder = trim(a);
    let b = trim(b);
    if remainder.len() < b.len() {
        return (vec![BigRational::zero()], remainder);
    }

    let mut quotient = vec![BigRational::zero(); remainder.len() - b.len() + 1];
    let b_lc = b[b.len() - 1].clone();
    while remainder.len() >= b.len() && !is_zero(&remainder) {
        let degree_diff = remainder.len() - b.len();
        let coeff = remainder[remainder.len() - 1].clone() / &b_lc;
        quotient[degree_diff] += &coeff;
        let mut subtract = vec![BigRational::zero(); degree_diff];
        subtract.extend(scale(&b, &coeff));
        remainder = sub(&remainder, &subtract);
    }
    (trim(&quotient), trim(&remainder))
}

fn all_coefficients_positive(p: &[BigRational]) -> bool {
    p.iter().all(BigRational::is_positive)
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    fn q(n: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(n))
    }

    #[test]
    fn checks_basic_positive_linear_chain() {
        // P_0 = (x+1)(x+2)-1 = x^2+3x+1, P_1=x+2.
        let p0 = vec![q(1), q(3), q(1)];
        let p1 = vec![q(2), q(1)];
        let summary =
            check_sturm_continued_fraction_certificate(&p0, &p1).expect("certificate should hold");
        assert_eq!(summary.degree_left, 2);
        assert_eq!(summary.degree_right, 1);
        assert_eq!(summary.steps, 2);
    }

    #[test]
    fn detects_nonpositive_signed_remainder() {
        // P_0 = (x+1)(x+2)+1 gives Euclidean remainder +1, so the signed
        // remainder is -1.
        let p0 = vec![q(3), q(3), q(1)];
        let p1 = vec![q(2), q(1)];
        let err = check_sturm_continued_fraction_certificate(&p0, &p1)
            .expect_err("certificate should fail");
        assert_eq!(
            err,
            SturmContinuedFractionError::SignedRemainderHasNonpositiveCoefficient { step: 1 }
        );
    }

    #[test]
    fn optionally_checks_quotient_root_order() {
        // q_1=x+1 has root -1, while q_2=x+3 has root -3, so the optional
        // increasing-root condition fails.
        let p0 = vec![q(2), q(4), q(1)];
        let p1 = vec![q(3), q(1)];
        let err = check_sturm_continued_fraction_certificate_with_options(
            &p0,
            &p1,
            SturmContinuedFractionOptions {
                require_positive_remainders: true,
                require_increasing_quotient_roots: true,
            },
        )
        .expect_err("quotient roots should not be increasing");
        assert!(matches!(
            err,
            SturmContinuedFractionError::QuotientRootsNotIncreasing { step: 2, .. }
        ));
    }
}
