//! No-remainder ray charts with exact rational multiples of pi at infinity.

use super::*;

pub(super) fn infinity_sign(p: &P) -> i8 {
    match p.degree() {
        None => 0,
        Some(d) => {
            if p.coeff(d).is_positive() == (d % 2 == 0) {
                1
            } else {
                -1
            }
        }
    }
}

/// Test rational squared cosines by cyclotomic divisibility, without a
/// trigonometric lookup table. If c=cos(2 theta), exp(2 i theta) solves
/// x^2-2c*x+1. A rational theta/pi therefore has cyclotomic degree <=2.
/// For every m, phi(m)^2 >= m/2: in the prime product for phi(m)^2/m,
/// each odd-prime factor is >=1, and the factor at 2 is >=1/2.
/// Thus m<=2*degree^2 is a complete order bound, not a search heuristic.
pub(super) fn angle_from_cosine_squared(square: &Q, sign: i8) -> Result<Q, RecurrenceError> {
    if square.is_negative() || square > &q(1) || (sign == 0 && !square.is_zero()) {
        return Err(failure("invalid squared cosine or cosine sign"));
    }
    let cosine_twice = square * q(2) - q(1);
    let equation = P::new(vec![q(1), -q(2) * cosine_twice, q(1)]);
    let degree = equation.degree().unwrap();
    for order in 1..=2 * degree * degree {
        let cyclotomic = P::new(
            crate::cyclotomic_polynomial_bigint(order)
                .into_iter()
                .map(Q::from_integer)
                .collect(),
        );
        if equation.div_rem(&cyclotomic).1.is_zero() {
            // Apart from +1, a cyclotomic polynomial of degree <=2 has
            // its upper-circle root at exp(2*pi*i/order). The two conjugate
            // roots give the same cosine; the supplied cosine sign selects
            // theta in [0,pi], rather than choosing a numerical root.
            let angle = if order == 1 {
                q(0)
            } else {
                Q::new(1.into(), order.into())
            };
            return Ok(if sign < 0 { q(1) - angle } else { angle });
        }
    }
    Err(failure("squared cosine gives no rational multiple of pi"))
}

pub(super) fn angle_from_tangent_squared(square: &Q, sign: i8) -> Result<Q, RecurrenceError> {
    if square.is_negative() || (sign == 0 && !square.is_zero()) {
        return Err(failure("invalid squared tangent or tangent sign"));
    }
    let angle = angle_from_cosine_squared(&(q(1) / (q(1) + square)), 1)?;
    Ok(angle * q(sign.into()))
}

pub(super) fn certify(
    n: &B,
    d: &[P],
    request: &EffectiveBkwRequest,
    degree: &DegreeBudget,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, ChartError> {
    let a = d[1].scale(&(q(-1) / q(2)));
    let b = d[2].clone();
    let y = b.clone() - a.pow(2);
    let disc = y.scale(&q(-4));
    let domain = Domain::select(&disc, 3, &request.interval)?;
    if domain.lower.is_some() || domain.upper.is_none() {
        return Err(failure(
            "quadratic Newton ray needs a negative-infinity edge and a finite collision",
        )
        .into());
    }
    let mut out = empty_certificate(d, &variable(), &domain, degree);
    if !legacy {
        set_domination_cutoff(&mut out, 0);
    }
    if domain.upper_is_nonpositive() {
        out.root_region = EffectiveBkwRootRegion::NegativeRay;
    }
    let n0 = P::new(n.rows().first().cloned().unwrap_or_default());
    let n1 = P::new(n.rows().get(1).cloned().unwrap_or_default());
    let h = n0.clone() * a.clone() + n1;
    for (label, p) in [
        ("quadratic imaginary part squared", &y),
        ("quadratic amplitude real part", &n0),
    ] {
        out.signs.push(check_univariate(label, p, &domain, 1)?);
    }
    if !cubic::endpoint_is_zero(&domain, &disc, false)? || domain.endpoint_sign(&a, false)? != 1 {
        return Err(failure("quadratic finite edge is not a positive pair collision").into());
    }
    let right_sign = -domain.endpoint_sign(&h, false)?;
    let ad = a.degree().unwrap_or(0);
    let bd = b.degree().unwrap_or(0);
    let theta = if a.is_zero() || 2 * ad < bd {
        q(1) / q(2)
    } else if 2 * ad == bd {
        let ratio = a.coeff(ad).pow(2) / (b.coeff(bd) * q(if bd % 2 == 0 { 1 } else { -1 }));
        angle_from_cosine_squared(&ratio, infinity_sign(&a))?
    } else {
        return Err(failure("quadratic infinity pair scaling is inconsistent").into());
    };
    let alpha = if h.is_zero() {
        q(0)
    } else {
        let hd = h.degree().unwrap();
        let nd = n0.degree().unwrap();
        let yd = y.degree().unwrap();
        match (2 * hd).cmp(&(2 * nd + yd)) {
            std::cmp::Ordering::Less => q(0),
            std::cmp::Ordering::Greater => q((-infinity_sign(&h)).into()) / q(2),
            std::cmp::Ordering::Equal => {
                let square = h.coeff(hd).pow(2)
                    / (n0.coeff(nd).pow(2) * y.coeff(yd) * q(if yd % 2 == 0 { 1 } else { -1 }));
                angle_from_tangent_squared(&square, -infinity_sign(&h))?
            }
        }
    };
    out.phase_limits = EffectiveBkwPhaseLimits {
        left_eigenvalue_argument: theta.to_string(),
        right_eigenvalue_argument: "0".into(),
        left_amplitude_argument: alpha.to_string(),
        right_amplitude_argument: (q(right_sign.into()) / q(2)).to_string(),
    };
    out.branch_selection =
        "explicit quadratic pair; positive imaginary part throughout the ray".into();
    out.expansion_data.pair_real_part = R::of(t(&a)).wire();
    out.expansion_data.pair_norm_squared = R::of(t(&b)).wire();
    out.expansion_data.cofactor = vec![R::of(c(1)).wire()];
    out.expansion.push("lambda=a+i*sqrt(b-a^2); A=N_0/2-i*(N_0*a+N_1)/(2*sqrt(b-a^2)); P_n=2 Re(A lambda^n) for every n>=0".into());
    out.phase.push(format!("Re A>0 fixes its continuous lift in (-pi/2,pi/2); polynomial leading terms give infinity limits theta={theta}*pi, arg A={alpha}*pi; finite collision gives arg A={right_sign}*pi/2"));
    out.domination = "no_remainder: proper quadratic numerator, exact partial fractions".into();
    degree::close_phase(&mut out)?;
    Ok(out)
}
