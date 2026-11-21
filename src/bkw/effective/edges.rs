//! Exact Newton-face limits on a simple rational algebraic branch.

use super::*;

/// The branch variable grows as s^-growth and t=-s^-weight. A simple
/// positive limiting root identifies the selected largest positive branch.
/// A rational chart instead uses its exact dummy branch r=0.
/// Taylor coefficients are exact consequences of the implicit function theorem.
#[allow(clippy::too_many_arguments)]
pub(super) fn pair_infinity(
    p: &B,
    a: &R,
    b: &R,
    re: &R,
    im: &R,
    y2: &R,
    growth: usize,
    weight: usize,
) -> Result<(Q, Q, Vec<String>), RecurrenceError> {
    let (scaled, _) = cubic::newton_scaled(p, growth, weight, &q(1), &q(1));
    let limit = P::new(scaled.evaluate_variable(&q(0)));
    let root = if p == &variable() {
        q(0)
    } else {
        let factors =
            crate::factorization::irreducible_factorization_rational_coeffs(limit.coeffs())
                .map_err(|e| failure(format!("pair Newton-face factorization: {e}")))?;
        let root_box = branch_box(&scaled, &q(0), Branch::LargestPositive, &(q(1) / q(256)))?;
        factors
            .factors
            .iter()
            .find_map(|(f, m)| {
                if *m == 1 && f.len() == 2 {
                    let x = -Q::new(f[0].clone(), f[1].clone());
                    if x.is_positive() && root_box.lo <= x && x <= root_box.hi {
                        return Some(x);
                    }
                }
                None
            })
            .ok_or_else(|| {
                failure("largest positive Newton pair-product root is not simple rational")
            })?
    };
    let derivative = limit.derivative().evaluate(&root);
    if derivative.is_zero() {
        return Err(failure("pair Newton-face root is not simple"));
    }
    let mut u = P::new(vec![root.clone()]);
    for j in 1..=16 {
        let coefficient = -cubic::truncated_evaluate(&scaled, &u).coeff(j) / &derivative;
        let mut row = u.coeffs().to_vec();
        row.resize(j + 1, q(0));
        row[j] = coefficient;
        u = P::new(row);
    }
    let leading_poly = |f: &B| -> Result<(i64, Q), RecurrenceError> {
        let f = rem(f, p);
        if f.is_zero() {
            return Ok((0, q(0)));
        }
        let (scaled, m) = cubic::newton_scaled(&f, growth, weight, &q(1), &q(1));
        let series = cubic::truncated_evaluate(&scaled, &u);
        let j = series
            .coeffs()
            .iter()
            .position(|c| !c.is_zero())
            .ok_or_else(|| failure("pair Newton leading order exceeds sixteen terms"))?;
        Ok((j as i64 - m as i64, series.coeff(j)))
    };
    let leading = |r: &R| -> Result<(i64, Q), RecurrenceError> {
        let (on, cn) = leading_poly(&r.n)?;
        let (od, cd) = leading_poly(&r.d)?;
        if cd.is_zero() {
            return Err(failure("pair Newton rational denominator is zero"));
        }
        Ok((on - od, cn / cd))
    };
    let (oa, ca) = leading(a)?;
    let (ob, cb) = leading(b)?;
    let (oy, cy) = leading(y2)?;
    let (or, cr) = leading(re)?;
    let (oi, ci) = leading(im)?;
    if !cb.is_positive() || !cy.is_positive() || (cr.is_zero() && ci.is_zero()) {
        return Err(failure("pair Newton limits have inconsistent signs"));
    }
    let cos2 = if ca.is_zero() || 2 * oa > ob {
        q(0)
    } else if 2 * oa == ob {
        &ca * &ca / &cb
    } else {
        return Err(failure("pair Newton real part exceeds its norm"));
    };
    let theta = quadratic::angle_from_cosine_squared(&cos2, if ca.is_negative() { -1 } else { 1 })?;
    let exponent = 2 * oi - oy - 2 * or;
    let sign_im = if ci.is_negative() { -1 } else { 1 };
    let alpha = if ci.is_zero() || (!cr.is_zero() && exponent > 0) {
        if cr.is_positive() {
            q(0)
        } else {
            q(sign_im)
        }
    } else if cr.is_zero() || exponent < 0 {
        q(sign_im) / q(2)
    } else {
        quadratic::angle_from_cosine_squared(
            &(&cr * &cr * &cy / (&cr * &cr * &cy + &ci * &ci)),
            if cr.is_positive() { 1 } else { -1 },
        )? * q(sign_im)
    };
    Ok((theta.clone(),alpha.clone(),vec![format!("Newton pair face: t=-s^(-{weight}), r=u*s^(-{growth}); scaled branch {:?}, simple selected limit {root}, derivative {derivative}; Taylor {:?}; leading (order,coefficient) a=({oa},{ca}), b=({ob},{cb}), Y=({oy},{cy}), Re numerator=({or},{cr}), Im numerator=({oi},{ci}); theta={theta}*pi, arg A={alpha}*pi",encode_b(&scaled),encode(&u))]))
}
