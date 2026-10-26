//! Coefficient-only real-rootedness criteria.
//!
//! These checks are exact filters on coefficient sequences.  Kurtz's condition is
//! sufficient for distinct real negative roots when all coefficients are
//! positive; Newton's inequalities are necessary for real-rootedness of a
//! nonnegative coefficient polynomial.

use crate::vec_poly::{degree, strip_trailing_zeros};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

/// One coefficient inequality check.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientInequalityCheck {
    pub index: usize,
    pub lhs: BigInt,
    pub rhs: BigInt,
    pub comparison: &'static str,
    pub holds: bool,
}

/// Detailed report for one coefficient criterion.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientCriterionReport {
    pub name: &'static str,
    pub reference: &'static str,
    pub applicable: bool,
    pub holds: bool,
    pub implies_real_rooted: bool,
    pub reason: Option<String>,
    pub checks: Vec<CoefficientInequalityCheck>,
}

/// Detailed coefficient-test report for one polynomial.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientTestReport {
    pub coefficients: Vec<BigInt>,
    pub degree: Option<usize>,
    pub newton: CoefficientCriterionReport,
    pub kurtz: CoefficientCriterionReport,
}

fn strip_initial_zeros(coeffs: &[BigInt]) -> &[BigInt] {
    match coeffs.iter().position(|c| !c.is_zero()) {
        Some(index) => &coeffs[index..],
        None => &[],
    }
}

/// Detailed Newton-inequality report for a nonnegative coefficient sequence.
pub fn newton_inequality_report_bigint(coeffs: &[BigInt]) -> CoefficientCriterionReport {
    let mut p = strip_initial_zeros(coeffs).to_vec();
    strip_trailing_zeros(&mut p);
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => {
            return CoefficientCriterionReport {
                name: "Newton inequalities",
                reference: "Newton inequalities for real-rooted nonnegative coefficient sequences",
                applicable: false,
                holds: true,
                implies_real_rooted: false,
                reason: Some("degree less than 2 after removing initial zeros".to_string()),
                checks: Vec::new(),
            };
        }
    };

    if p.iter().take(d + 1).any(|c| c.is_negative()) {
        return CoefficientCriterionReport {
            name: "Newton inequalities",
            reference: "Newton inequalities for real-rooted nonnegative coefficient sequences",
            applicable: true,
            holds: false,
            implies_real_rooted: false,
            reason: Some("one or more coefficients are negative".to_string()),
            checks: Vec::new(),
        };
    }

    let mut binom = vec![BigInt::one(); d + 1];
    for k in 1..=d {
        binom[k] = &binom[k - 1] * BigInt::from(d - k + 1) / BigInt::from(k);
    }

    let mut checks = Vec::new();
    for k in 1..d {
        let lhs = p[k].pow(2) * &binom[k - 1] * &binom[k + 1];
        let rhs = &p[k - 1] * &p[k + 1] * &binom[k] * &binom[k];
        checks.push(CoefficientInequalityCheck {
            index: k,
            holds: lhs >= rhs,
            lhs,
            rhs,
            comparison: ">=",
        });
    }
    let holds = checks.iter().all(|check| check.holds);
    CoefficientCriterionReport {
        name: "Newton inequalities",
        reference: "Newton inequalities for real-rooted nonnegative coefficient sequences",
        applicable: true,
        holds,
        implies_real_rooted: false,
        reason: None,
        checks,
    }
}

/// Detailed Kurtz-condition report.
pub fn kurtz_condition_report_bigint(coeffs: &[BigInt]) -> CoefficientCriterionReport {
    let mut p = strip_initial_zeros(coeffs).to_vec();
    strip_trailing_zeros(&mut p);
    let d = match degree(&p) {
        Some(d) if d >= 2 => d,
        _ => {
            return CoefficientCriterionReport {
                name: "Kurtz condition",
                reference:
                    "Kurtz, A sufficient condition for all the roots of a polynomial to be real",
                applicable: false,
                holds: false,
                implies_real_rooted: false,
                reason: Some("degree less than 2 after removing initial zeros".to_string()),
                checks: Vec::new(),
            };
        }
    };

    if p.iter().take(d + 1).any(|c| !c.is_positive()) {
        return CoefficientCriterionReport {
            name: "Kurtz condition",
            reference: "Kurtz, A sufficient condition for all the roots of a polynomial to be real",
            applicable: true,
            holds: false,
            implies_real_rooted: false,
            reason: Some("Kurtz requires all coefficients to be positive".to_string()),
            checks: Vec::new(),
        };
    }

    let mut checks = Vec::new();
    for k in 1..d {
        let lhs = p[k].pow(2);
        let rhs = BigInt::from(4) * &p[k - 1] * &p[k + 1];
        checks.push(CoefficientInequalityCheck {
            index: k,
            holds: lhs > rhs,
            lhs,
            rhs,
            comparison: ">",
        });
    }
    let holds = checks.iter().all(|check| check.holds);
    CoefficientCriterionReport {
        name: "Kurtz condition",
        reference: "Kurtz, A sufficient condition for all the roots of a polynomial to be real",
        applicable: true,
        holds,
        implies_real_rooted: holds,
        reason: None,
        checks,
    }
}

/// Run all coefficient-only criteria currently implemented.
pub fn coefficient_test_report_bigint(coeffs: &[BigInt]) -> CoefficientTestReport {
    let mut coefficients = coeffs.to_vec();
    strip_trailing_zeros(&mut coefficients);
    let degree = degree(&coefficients);
    CoefficientTestReport {
        newton: newton_inequality_report_bigint(&coefficients),
        kurtz: kurtz_condition_report_bigint(&coefficients),
        coefficients,
        degree,
    }
}
