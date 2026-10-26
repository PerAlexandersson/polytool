//! Coefficient-only real-rootedness criteria.
//!
//! These checks are exact filters on coefficient sequences.  Kurtz's condition is
//! sufficient for distinct real negative roots when all coefficients are
//! positive; Newton's inequalities are necessary for real-rootedness of a
//! nonnegative coefficient polynomial.
//!
//! Entry points: [`coefficient_test_report_bigint`] (both criteria),
//! [`newton_inequality_report_bigint`] and [`kurtz_condition_report_bigint`].
//!
//! Conventions: coefficients are ascending `BigInt` (`coeffs[i]` is the
//! coefficient of `t^i`); factors of `t` (initial zeros) are divided out and
//! trailing zeros ignored before testing, and all comparisons are exact integer
//! comparisons.  Neither criterion decides real-rootedness in general: a failed
//! Newton check proves the polynomial is **not** real-rooted only for
//! nonnegative coefficients, and a passed Kurtz check proves real-rootedness
//! (distinct negative roots) only for positive coefficients; use
//! [`crate::is_real_rooted_bigint_coeffs`] for an exact decision.

use crate::vec_poly::{degree, strip_trailing_zeros};
use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};

/// One coefficient inequality check.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientInequalityCheck {
    /// Coefficient index `k` (position after removing initial zeros) around
    /// which the inequality is centered, `1 <= k <= d - 1`.
    pub index: usize,
    /// Left-hand side, an exact integer (scaled so that no division occurs).
    pub lhs: BigInt,
    /// Right-hand side, an exact integer.
    pub rhs: BigInt,
    /// The required relation `lhs <comparison> rhs`: `">="` (Newton) or `">"`
    /// (Kurtz).
    pub comparison: &'static str,
    /// Whether the relation holds.
    pub holds: bool,
}

/// Detailed report for one coefficient criterion.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientCriterionReport {
    /// Name of the criterion (`"Newton inequalities"` or `"Kurtz condition"`).
    pub name: &'static str,
    /// Literature reference for the criterion.
    pub reference: &'static str,
    /// Whether the degree (after removing initial zeros) is at least 2, so the
    /// criterion has something to test.  When `false`, `checks` is empty.
    pub applicable: bool,
    /// Whether the criterion's inequalities (and sign hypothesis) hold.  When
    /// `applicable` is `false` this is `true` for Newton and `false` for Kurtz.
    pub holds: bool,
    /// Whether `holds` proves real-rootedness (only Kurtz, when it holds).
    pub implies_real_rooted: bool,
    /// Explanation when the criterion is not applicable or its sign hypothesis
    /// (nonnegative for Newton, positive for Kurtz) fails; otherwise `None`.
    pub reason: Option<String>,
    /// The individual inequalities `k = 1, ..., d - 1`; empty if none were
    /// evaluated.
    pub checks: Vec<CoefficientInequalityCheck>,
}

/// Detailed coefficient-test report for one polynomial.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoefficientTestReport {
    /// The input coefficients with trailing zeros removed (initial zeros kept).
    pub coefficients: Vec<BigInt>,
    /// Degree of the polynomial, `None` for the zero polynomial.
    pub degree: Option<usize>,
    /// The Newton-inequality report.
    pub newton: CoefficientCriterionReport,
    /// The Kurtz-condition report.
    pub kurtz: CoefficientCriterionReport,
}

fn strip_initial_zeros(coeffs: &[BigInt]) -> &[BigInt] {
    match coeffs.iter().position(|c| !c.is_zero()) {
        Some(index) => &coeffs[index..],
        None => &[],
    }
}

/// Detailed Newton-inequality report for a nonnegative coefficient sequence.
///
/// For degree `d >= 2` (after removing initial zeros) checks
/// `p_k^2 C(d,k-1) C(d,k+1) >= p_{k-1} p_{k+1} C(d,k)^2` for `1 <= k <= d-1`,
/// which is Newton's inequality `p_k^2 >= (1 + 1/k)(1 + 1/(d-k)) p_{k-1} p_{k+1}`
/// cleared of denominators.  A negative coefficient is reported as `holds =
/// false` without evaluating inequalities (the criterion assumes nonnegative
/// coefficients).
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
///
/// For a polynomial of degree `d >= 2` (after removing initial zeros) with all
/// coefficients positive, checks `p_k^2 > 4 p_{k-1} p_{k+1}` for `1 <= k <= d-1`.
/// If all hold, all roots are real, negative and distinct.
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
