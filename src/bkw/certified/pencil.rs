//! Discriminant-cell decision for one rational polynomial pencil.

use super::*;
use algebra::{cells, derivative_in_outer, resultant};

pub(super) fn decide(u: &P, v: &P, upper: &Q) -> Result<PencilCertificate, RecurrenceError> {
    let common = u.gcd(v);
    let a = u.exact_div(&common);
    let b = v.exact_div(&common);
    let degree = a.degree().unwrap_or(0).max(b.degree().unwrap_or(0));
    let outer: Vec<_> = (0..=degree)
        .map(|i| P::new(vec![a.coeff(i), b.coeff(i)]))
        .collect();
    let critical = if degree > 0 {
        resultant(&outer, &derivative_in_outer(&outer)) * outer[degree].clone()
    } else {
        P::one()
    };
    if critical.is_zero() {
        return Err(RecurrenceError::Unsupported(
            "identically singular reduced pencil".into(),
        ));
    }
    let (intervals, mut parameters) = cells(&critical, &Q::zero(), upper);
    // Near a degree drop the zero set can lose roots to infinity. Every finite
    // root of a nonzero boundary specialization is a limit of nearby roots.
    // Therefore successful adjacent cells certify irrational boundaries too.
    // The primitive pencil cannot vanish identically unless it is scalar.
    let mut witness = None;
    let mut samples = Vec::new();
    for s in parameters.drain(..) {
        let p = u.clone() + v.scale(&s);
        let sample = PencilSample {
            parameter: s.to_string(),
            polynomial: encode(&p),
            root_counts: counts(&p),
        };
        if !sample.root_counts.real_rooted() && witness.is_none() && !p.is_zero() {
            witness = Some(sample.clone());
        }
        samples.push(sample);
    }
    for s in [Q::zero(), upper.clone()] {
        let p = u.clone() + v.scale(&s);
        samples.push(PencilSample {
            parameter: s.to_string(),
            polynomial: encode(&p),
            root_counts: counts(&p),
        });
    }
    let all_real_rooted = samples
        .iter()
        .all(|s| s.root_counts.degree.is_none() || s.root_counts.real_rooted());
    Ok(PencilCertificate {
        constant: encode(u),
        slope: encode(v),
        interval: ["0".into(), upper.to_string()],
        critical_polynomial: encode(&critical),
        critical_intervals: intervals
            .iter()
            .map(|(a, b)| [a.to_string(), b.to_string()])
            .collect(),
        samples,
        all_real_rooted,
        witness,
    })
}

/// Decide whether every nonzero member of U+sV is real-rooted for 0<=s<=upper.
///
/// Repeated parameter-independent factors are removed before forming the
/// discriminant, but included in every exact sample root count. Critical
/// algebraic parameters inherit reality from adjacent cells by continuity;
/// rational endpoints are also checked explicitly. An identically zero
/// specialization contributes no finite tie point.
///
/// # Errors
/// Invalid coefficients or bounds, an identically zero pencil, or a scalar
/// pencil that vanishes identically at an interior parameter.
pub fn certify_real_rooted_pencil(
    constant: &[String],
    slope: &[String],
    upper: &str,
) -> Result<PencilCertificate, RecurrenceError> {
    let (u, v) = (decode(constant)?, decode(slope)?);
    let upper = parse_rational_coeff(upper).map_err(RecurrenceError::from)?;
    if upper <= Q::zero() || (u.is_zero() && v.is_zero()) {
        return Err(RecurrenceError::InvalidInput(
            "nonzero pencil and positive upper bound required".into(),
        ));
    }
    let common = u.gcd(&v);
    let (a, b) = (u.exact_div(&common), v.exact_div(&common));
    if a.degree().unwrap_or(0) == 0 && b.degree().unwrap_or(0) == 0 && !b.is_zero() {
        let zero = -a.coeff(0) / b.coeff(0);
        if zero > Q::zero() && zero < upper {
            return Err(RecurrenceError::Unsupported(
                "scalar pencil vanishes identically at an interior parameter".into(),
            ));
        }
    }
    decide(&u, &v, &upper)
}
