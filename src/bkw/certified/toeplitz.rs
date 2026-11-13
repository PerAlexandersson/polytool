//! Exact Toeplitz circle certificates for an unpublished manuscript theorem.
//!
//! Circle criterion: let b(z)=z^-1+b_0+sum_(j=1)^s b_j*z^j have real
//! coefficients, with standard rows P_0=1, P_m=0 for m<0 and
//! P_n(t)=(t-b_0)*P_(n-1)(t)-sum_(j=1)^s b_j*P_(n-j-1)(t).
//! If some R>0 satisfies Im b(R*exp(i*theta))>=0 for
//! 0<theta<pi, every standard Hessenberg Toeplitz row has simple real zeros
//! and consecutive rows strictly interlace.
//! Source: circle criterion for Hessenberg Toeplitz symbols (unpublished manuscript).
//! We check nonnegativity of Q_R exactly and apply this theorem.

use super::*;
use algebra::{chain, nonnegative};

/// Available proved criteria for standard Hessenberg Toeplitz rows.
///
/// This list is the extension point for future theorem criteria.
pub const TOEPLITZ_CRITERIA: &[&str] = &["circle"];

/// Theorem source for circle certificates, distinguished from published sources.
pub const CIRCLE_THEOREM_SOURCE: &str =
    "circle criterion for Hessenberg Toeplitz symbols (unpublished manuscript)";

pub(super) fn hessenberg_tails(d: &[P]) -> Option<Vec<String>> {
    if d.len() < 2
        || d[1].degree() != Some(1)
        || d.iter().skip(2).any(|p| p.degree().is_some_and(|j| j > 0))
    {
        return None;
    }
    Some(d.iter().skip(2).map(|p| p.coeff(0).to_string()).collect())
}

fn circle_polynomial(tails: &[Q], radius: &Q) -> P {
    let mut result =
        P::constant(tails.first().cloned().unwrap_or_else(Q::zero) * radius - radius.recip());
    let mut previous = P::one();
    let mut current = P::new(vec![Q::zero(), q(2)]);
    for (j, coefficient) in tails.iter().enumerate().skip(1) {
        result = result + current.scale(&(coefficient * radius.pow((j + 1) as i32)));
        let next = P::new(vec![Q::zero(), q(2)]) * current.clone() - previous;
        previous = current;
        current = next;
    }
    result
}

pub(super) fn recognize_and_certify(
    d: &[P],
    options: &CertificationOptions,
) -> Result<Option<CircleCertificate>, RecurrenceError> {
    if d.len() < 3 || hessenberg_tails(d).is_none() {
        return Ok(None);
    }
    let tails: Vec<_> = d.iter().skip(2).map(|p| p.coeff(0)).collect();
    let mut radii = Vec::new();
    if let Some(radius) = &options.circle_radius {
        radii.push(parse_rational_coeff(radius).map_err(RecurrenceError::from)?);
    }
    radii.push(Q::one());
    // b'(R)=0 iff -1+sum l*b_l R^(l+1)=0. Isolate exactly;
    // sample around the smallest positive critical radius with no floats.
    // Reflection sends b_j to (-1)^(j+1)b_j and replaces u by -u
    // in Q_R. Searching both critical polynomials preserves the predicate.
    for reflected in [false, true] {
        let mut critical = vec![-Q::one(), Q::zero()];
        critical.extend(tails.iter().enumerate().map(|(j, b)| {
            let coefficient = b * q((j + 1) as i64);
            if reflected && j % 2 == 1 {
                -coefficient
            } else {
                coefficient
            }
        }));
        let roots = chain(&P::new(critical)).isolate_roots(&Q::new(1.into(), 4096.into()));
        if let Some((a, b)) = roots.into_iter().find(|(a, _)| a > &Q::zero()) {
            radii.extend([
                a.clone(),
                (&a + &b) / q(2),
                b.clone(),
                &b + Q::new(1.into(), 4096.into()),
            ]);
        }
    }
    // For a tridiagonal symbol every R>=1/sqrt(b1) works. Doubling is
    // an exact candidate search over a wide range of radii.
    if tails.len() == 1 && tails[0] > Q::zero() {
        let mut r = Q::one();
        for _ in 0..128 {
            if &tails[0] * &r * &r >= Q::one() {
                radii.push(r);
                break;
            }
            r *= q(2);
        }
    }
    for r in radii {
        let polynomial = circle_polynomial(&tails, &r);
        if nonnegative(&polynomial, &q(-1), &Q::one()) {
            return Ok(Some(CircleCertificate {
                affine: encode(&(-d[1].clone())),
                tails: tails.iter().map(ToString::to_string).collect(),
                radius: r.to_string(),
                polynomial: encode(&polynomial),
                root_counts: counts(&polynomial),
                provenance: CIRCLE_THEOREM_SOURCE.into(),
                conclusion:
                    "every row has simple real zeros and consecutive rows strictly interlace".into(),
            }));
        }
    }
    Ok(None)
}
