//! Exact pointwise domination for two real modes or a conjugate cofactor pair.

use super::*;

// Denominators in the quartic formulas are powers of already positive
// b and b*(b-c). Remove only exact products of these factors, keeping the
// sign of the remaining rational constant. This avoids introducing large
// artificial resultant factors by multiplying by the denominator again.
fn positive_pair_difference(left: &[&R], right: &[&R], p: &B, positive_factors: &[B]) -> B {
    use crate::recurrence::generating_function::exact_divide;
    let nonconstant =
        |f: &&B| f.index_degree().unwrap_or(0) > 0 || f.rows().iter().any(|r| r.len() > 1);
    let factors: Vec<_> = positive_factors.iter().filter(nonconstant).collect();
    let strip = |mut f: B| {
        if !f.is_zero() {
            for factor in &factors {
                while let Ok(quotient) = exact_divide(&f, factor) {
                    f = quotient;
                }
            }
        }
        f
    };
    let product = |values: &[&R]| {
        values.iter().fold((B::one(), B::one()), |(n, d), r| {
            (n.multiply(&r.n), d.multiply(&r.d))
        })
    };
    let (ln, ld) = product(left);
    let (rn, rd) = product(right);
    let numerator = ln.multiply(&rd).subtract(&rn.multiply(&ld));
    let raw_denominator = ld.multiply(&rd);
    let denominator = strip(raw_denominator.clone());
    if denominator.index_degree().unwrap_or(0) == 0
        && denominator.rows().first().is_some_and(|r| r.len() == 1)
    {
        let sign = if denominator.rows()[0][0].is_positive() {
            q(1)
        } else {
            q(-1)
        };
        rem(&strip(numerator).scale(&sign), p)
    } else {
        rem(&numerator.multiply(&raw_denominator), p)
    }
}

// A nonreal quadratic cofactor is a second conjugate pair. With
// Nbar mod h=u*mu+v and q mod h=e*mu+f, its residue norm is
// 4|B|^2=(c*u^2+v^2-s*u*v)/((c-s^2/4)*(f^2-s*e*f+c*e^2)).
#[allow(clippy::too_many_arguments)]
fn pointwise_conjugate_remainder(
    n: &B,
    p: &B,
    a: &R,
    b: &R,
    h: &[R],
    nn: &R,
    hn: &R,
    domain: &Domain,
    which: Branch,
) -> Result<(usize, Vec<EffectiveBkwSignCertificate>), RecurrenceError> {
    let s = &h[1];
    let c0 = &h[2];
    let yh = c0.sub(&s.square().scale(&(q(1) / q(4))));
    let mut signs = vec![
        check_sign(
            "nonreal cofactor pair",
            p,
            &yh.signed_polynomial(p),
            domain,
            which,
            1,
        )?,
        check_sign(
            "strict gap between conjugate-pair norms",
            p,
            &b.sub(c0).signed_polynomial(p),
            domain,
            which,
            1,
        )?,
    ];
    let nc = |j| R::of(t(&P::new(n.rows().get(j).cloned().unwrap_or_default())));
    let u = nc(2)
        .sub(&nc(1).mul(s))
        .add(&nc(0).mul(&s.square().sub(c0)));
    let v = nc(3).sub(&nc(1).mul(c0)).add(&nc(0).mul(s).mul(c0));
    let e = s.neg().sub(&a.scale(&q(2)));
    let f = b.sub(c0);
    let qq = f.square().sub(&s.mul(&e).mul(&f)).add(&c0.mul(&e.square()));
    let norm = c0.mul(&u.square()).add(&v.square()).sub(&s.mul(&u).mul(&v));
    let y = b.sub(&a.square());
    let g = b.mul(&b.sub(c0));
    let factors: Vec<_> = [b, &g]
        .into_iter()
        .filter(|r| r.d == B::one())
        .map(|r| r.n.clone())
        .collect();
    let sample = branch_box(p, &domain.sample, which, &(q(1) / q(65536)))?;
    let mut bp = b.clone();
    let mut cp = R::of(c(1));
    for cutoff in 0..=16 {
        let condition =
            positive_pair_difference(&[nn, &bp, &yh, &qq], &[hn, &y, &norm, &cp], p, &factors);
        if eval_b(&condition, &I::point(domain.sample.clone()), &sample).sign() >= 0 {
            match check_sign(
                &format!("pointwise conjugate amplitude |A*lambda^N|^2>|B*mu^N|^2, N={cutoff}"),
                p,
                &condition,
                domain,
                which,
                1,
            ) {
                Ok(sign) => {
                    signs.push(sign);
                    return Ok((cutoff, signs));
                }
                Err(RecurrenceError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }
        bp = bp.mul(b);
        cp = cp.mul(c0);
    }
    Err(failure(
        "pointwise conjugate-remainder cutoff search exhausted (16)",
    ))
}

#[allow(clippy::too_many_arguments)]
pub(super) fn certify(
    n: &B,
    p: &B,
    a: &R,
    b: &R,
    h: &[R],
    nn: &R,
    hn: &R,
    domain: &Domain,
    which: Branch,
) -> Result<(usize, Vec<EffectiveBkwSignCertificate>, &'static str), RecurrenceError> {
    match pointwise_two_real_modes(n, p, a, b, h, nn, hn, domain, which) {
        Ok((cutoff, signs)) => Ok((cutoff, signs, "pointwise_two_real_remainders")),
        Err(RecurrenceError::Unsupported(_)) => {
            pointwise_conjugate_remainder(n, p, a, b, h, nn, hn, domain, which)
                .map(|(cutoff, signs)| (cutoff, signs, "pointwise_conjugate_remainder"))
        }
        Err(error) => Err(error),
    }
}

// If h(mu)=mu^2+s*mu+c has distinct real roots x,y, reduce
// Nbar(mu)/q(mu) modulo h to U*mu+V. Then the real residues obey
// Cx+Cy=U and (Cx-Cy)^2=(-s*U+2*V)^2/(s^2-4c).
// Both inequalities below are symmetric; no individual root is selected.
#[allow(clippy::too_many_arguments)]
fn pointwise_two_real_modes(
    n: &B,
    p: &B,
    a: &R,
    b: &R,
    h: &[R],
    nn: &R,
    hn: &R,
    domain: &Domain,
    which: Branch,
) -> Result<(usize, Vec<EffectiveBkwSignCertificate>), RecurrenceError> {
    let s = &h[1];
    let c0 = &h[2];
    let disc = s.square().sub(&c0.scale(&q(4)));
    let real_roots = match check_sign(
        "negative cofactor product implies distinct real roots",
        p,
        &c0.signed_polynomial(p),
        domain,
        which,
        -1,
    ) {
        Ok(sign) => sign,
        Err(RecurrenceError::Unsupported(_)) => check_sign(
            "real cofactor discriminant",
            p,
            &disc.signed_polynomial(p),
            domain,
            which,
            1,
        )?,
        Err(error) => return Err(error),
    };
    let mut signs = vec![real_roots];
    let y2 = b.sub(&a.square());
    let nc = |j| R::of(t(&P::new(n.rows().get(j).cloned().unwrap_or_default())));
    // Nbar=n0*mu^3+n1*mu^2+n2*mu+n3.
    let u = nc(2)
        .sub(&nc(1).mul(s))
        .add(&nc(0).mul(&s.square().sub(c0)));
    let v = nc(3).sub(&nc(1).mul(c0)).add(&nc(0).mul(s).mul(c0));
    let e = s.neg().sub(&a.scale(&q(2)));
    let f = b.sub(c0);
    let denominator = f.square().sub(&s.mul(&e).mul(&f)).add(&c0.mul(&e.square()));
    let mut up = u.mul(&f).sub(&v.mul(&e));
    let mut vp = v.mul(&f.sub(&e.mul(s))).add(&u.mul(&e).mul(c0));
    // nn/hn=|F(1/lambda)|^2, 4|A|^2=b*nn/(y2*hn).
    let mut pair_power = b.clone();
    for (name, value) in [
        ("cofactor root product below pair modulus", b.sub(c0)),
        (
            "cofactor roots strictly inside pair modulus, first",
            b.add(c0),
        ),
        (
            "cofactor roots strictly inside pair modulus, second",
            b.add(c0).square().sub(&b.mul(&s.square())),
        ),
    ] {
        signs.push(check_sign(
            name,
            p,
            &value.signed_polynomial(p),
            domain,
            which,
            1,
        )?);
    }
    // Since h has two distinct real roots and the pair is nonreal, h and
    // q are coprime. Thus Q=q(x)q(y)>0 and |H(1/lambda)|^2>0 without a
    // second elimination of their squared norms.
    let g = b.mul(&b.sub(c0));
    let positive_factors: Vec<_> = [b, &g]
        .into_iter()
        .filter(|r| r.d == B::one())
        .map(|r| r.n.clone())
        .collect();
    let sample = branch_box(p, &domain.sample, which, &(q(1) / q(65536)))?;
    for cutoff in 0..=16 {
        let minus = s.neg().mul(&up).add(&vp.scale(&q(2)));
        let conditions = [
            (
                "pointwise amplitude plus",
                positive_pair_difference(
                    &[nn, &pair_power, &denominator, &denominator],
                    &[hn, &y2, &up, &up],
                    p,
                    &positive_factors,
                ),
            ),
            (
                "pointwise amplitude minus",
                positive_pair_difference(
                    &[nn, &pair_power, &denominator, &denominator, &disc],
                    &[hn, &y2, &minus, &minus],
                    p,
                    &positive_factors,
                ),
            ),
        ];
        // This sample may reject a candidate; both global sign certificates
        // are required for acceptance. Multiplying residues by mu^N gives
        // the exact remainder at N. Individual strict modulus gaps then
        // make the sum of absolute normalized modes decrease for n>=N.
        if !conditions
            .iter()
            .any(|(_, f)| eval_b(f, &I::point(domain.sample.clone()), &sample).sign() < 0)
        {
            let attempt: Result<Vec<_>, _> = conditions
                .iter()
                .map(|(name, f)| check_sign(&format!("{name}, n={cutoff}"), p, f, domain, which, 1))
                .collect();
            match attempt {
                Ok(evidence) => {
                    signs.extend(evidence);
                    return Ok((cutoff, signs));
                }
                Err(RecurrenceError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }
        let next_u = vp.sub(&s.mul(&up));
        vp = c0.mul(&up).neg();
        up = next_u;
        pair_power = pair_power.mul(b);
    }
    Err(failure(
        "pointwise two-real-mode cutoff search exhausted (16)",
    ))
}
