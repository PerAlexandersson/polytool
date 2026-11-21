//! Fractional quartic ray phases with exact pointwise remainder bounds.

use super::compact::{collision_root, quartic_amplitude_forms};
use super::cubic::{endpoint_rational_sign, substitute};
use super::*;

// Ray closure uses a phase of A itself; the compact path below uses the
// equivalent half-integer levels of (n+1)*theta+arg F.
#[allow(clippy::too_many_arguments)]
pub(super) fn certify(
    n: &B,
    d: &[P],
    degree: &DegreeBudget,
    domain: &Domain,
    characteristic: &B,
    p: &B,
    a: &R,
    b: &R,
    h: &[R],
    which: Branch,
    mut out: EffectiveBkwCertificate,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    let newton = degree
        .data
        .newton
        .as_ref()
        .ok_or_else(|| failure("pair ray requires a Newton degree law"))?;
    let y2 = b.sub(&a.square());
    let nc: Vec<_> = n
        .rows()
        .iter()
        .map(|row| R::of(t(&P::new(row.clone()))))
        .collect();
    let (mut nn, mut hn, mut imf, mut ref_) = quartic_amplitude_forms(n, d);
    if matches!(which, Branch::Zero) {
        let at_b = |value: &R| substitute(&value.n, b).div(&substitute(&value.d, b));
        nn = at_b(&nn);
        hn = at_b(&hn);
        imf = at_b(&imf);
        ref_ = at_b(&ref_);
    }

    let g = b.square().sub(&R::of(t(&d[4])));
    // Re F = b*ref/(2*hn), Im F/Im lambda=b*g*imf/hn.
    // Drop the common positive factor b/hn when fixing arg A.
    let rf = ref_.scale(&(q(1) / q(2)));
    let inf = g.mul(&imf);
    let re = rf.add(&a.mul(&inf));
    let im = y2.mul(&inf).sub(&a.mul(&rf));
    let real_lift = check_sign(
        "ray amplitude real part",
        p,
        &re.signed_polynomial(p),
        domain,
        which,
        1,
    );
    let positive_real = real_lift.is_ok();
    out.signs.push(match real_lift {
        Ok(sign) => sign,
        Err(RecurrenceError::Unsupported(_)) => check_sign(
            "ray amplitude lower half-plane",
            p,
            &im.signed_polynomial(p),
            domain,
            which,
            -1,
        )?,
        Err(error) => return Err(error),
    });
    let (x, identity) = collision_root(characteristic, domain, false)?;
    if endpoint_rational_sign(domain, &x, false)? != 1 {
        return Err(failure("ray finite collision is not positive"));
    }
    let s = R::of(t(&d[1])).add(&x.scale(&q(2)));
    let c0 = R::of(t(&d[2]))
        .add(&x.scale(&q(2)).mul(&s))
        .sub(&x.square());
    // The double root must belong to the dominant pair, independently of
    // which characteristic roots collide at this endpoint.
    let cofactor_y = c0.sub(&s.square().scale(&(q(1) / q(4))));
    let complex_cofactor = endpoint_rational_sign(domain, &cofactor_y, false)? == 1;
    let bounds = if complex_cofactor {
        vec![x.square().sub(&c0)]
    } else {
        vec![
            x.square().sub(&c0),
            x.square().add(&c0),
            x.square()
                .add(&c0)
                .square()
                .sub(&x.square().mul(&s.square())),
        ]
    };
    for bound in bounds {
        if endpoint_rational_sign(domain, &bound, false)? != 1 {
            return Err(failure("ray collision cofactor is not strictly smaller"));
        }
    }
    let inv = R::of(c(1)).div(&x);
    let numerator = nc.iter().rev().fold(R::of(c(0)), |v, c| v.mul(&inv).add(c));
    let denominator = R::of(c(1)).add(&s.mul(&inv)).add(&c0.mul(&inv.square()));
    let sign = endpoint_rational_sign(domain, &numerator.div(&denominator), false)?;
    if sign == 0 || (!positive_real && sign != 1) {
        return Err(failure(
            "ray finite amplitude limit is incompatible with its lift",
        ));
    }
    let (theta, alpha, evidence) = edges::pair_infinity(
        p,
        a,
        b,
        &re,
        &im,
        &y2,
        2 * newton.weight_numerator,
        newton.weight_denominator,
    )?;
    out.phase.extend(evidence);
    out.phase.push(format!("finite positive dominant collision: {identity}; F sign {sign}; theta->0, arg A->{}*pi/2 in the certified half-plane",-sign));
    out.phase_limits.left_eigenvalue_argument = theta.to_string();
    out.phase_limits.left_amplitude_argument = alpha.to_string();
    out.phase_limits.right_amplitude_argument = (-q(sign.into()) / q(2)).to_string();
    out.expansion_data = EffectiveBkwExpansionData {
        pair_real_part: a.wire(),
        pair_norm_squared: b.wire(),
        cofactor: h.iter().map(R::wire).collect(),
        remainder_root: None,
        remainder_amplitude: None,
    };
    out.expansion.push("D=(1-2*a*z+b*z^2)H, lambda=a+i*sqrt(b-a^2), F=N(1/lambda)/H(1/lambda), A=F*lambda/(2*i*Im lambda); proper partial fractions give P_n=2 Re(A lambda^n)+Cx*x^n+Cy*y^n".into());
    let (cutoff, signs, method) = pointwise::certify(n, p, a, b, h, &nn, &hn, domain, which)?;
    out.signs.extend(signs);
    set_domination_cutoff(&mut out, cutoff);
    out.domination = format!("{method}: exact residue norms at N={cutoff} and strict root-modulus gaps imply domination for every n>=N");
    degree::close_phase(&mut out)?;
    Ok(out)
}
