//! Uniform real-cofactor bounds from exact root counts and parameter barriers.

use super::*;

pub(super) struct BarrierData<'a> {
    pub d: &'a [P],
    pub characteristic: &'a B,
    pub p: &'a B,
    pub a: &'a R,
    pub b: &'a R,
    pub cofactor: &'a [R],
    pub domain: &'a Domain,
    pub which: Branch,
    pub numerator: &'a Q,
}

fn at_index(p: &B, x: &Q) -> P {
    p.rows()
        .iter()
        .rev()
        .fold(P::zero(), |v, row| v.scale(x) + P::new(row.clone()))
}

fn barrier(
    name: &str,
    p: &P,
    domain: &Domain,
) -> Result<EffectiveBkwSignCertificate, RecurrenceError> {
    let sign = if p.evaluate(&domain.sample).is_positive() {
        1
    } else {
        -1
    };
    if p.is_zero() || p.evaluate(&domain.sample).is_zero() || !domain.critical_roots(p)?.is_empty()
    {
        return Err(failure("uniform real-root barrier crosses the interval"));
    }
    check_univariate(name, p, domain, sign)
}

pub(super) fn certify(
    data: BarrierData<'_>,
    mut out: EffectiveBkwCertificate,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    let BarrierData {
        d,
        characteristic,
        p,
        a,
        b,
        cofactor,
        domain,
        which,
        numerator,
    } = data;
    let order = d.len() - 1;
    let m = order - 2;
    let sample_polynomial = characteristic.evaluate_variable(&domain.sample);
    let sturm = crate::sturm::SturmChain::from_bigint_coeffs(&vec_poly::clear_denominators(
        &sample_polynomial,
    ));
    if sturm.count_real_roots() != m {
        return Err(failure(
            "real-root barrier supplier requires exactly one nonreal pair",
        ));
    }
    let sample = branch_box(p, &domain.sample, which, &(q(1) / q(1 << 32)))?;
    let mut chosen = None;
    'radius: for k in (1..=16).rev() {
        let s = q(k) / q(16);
        if sample.lo <= &s * &s {
            continue;
        }
        let lower = match barrier(
            "pair-product lower-radius barrier",
            &at_index(p, &(&s * &s)),
            domain,
        ) {
            Ok(c) => c,
            Err(RecurrenceError::Unsupported(_)) => continue,
            Err(e) => return Err(e),
        };
        'relative: for j in 1..16 {
            let relative = q(j) / q(16);
            let gamma: Vec<_> = (1..=m)
                .map(|i| {
                    Q::from_integer(crate::numeric::binomial_bigint(m, i)) * relative.pow(i as i32)
                })
                .collect();
            let Some(ratio) = (1..16).map(|i| q(i) / q(16)).find(|candidate| {
                gamma
                    .iter()
                    .enumerate()
                    .map(|(i, g)| g / candidate.pow((i + 1) as i32))
                    .sum::<Q>()
                    < q(1)
            }) else {
                break;
            };
            let rho = &s * &relative;
            if sturm.count_roots_in(&(-rho.clone()), &rho) != m {
                continue;
            }
            let mut signs = vec![lower.clone()];
            for x in [-rho.clone(), rho.clone()] {
                match barrier(
                    &format!("real-cofactor root barrier lambda={x}"),
                    &at_index(characteristic, &x),
                    domain,
                ) {
                    Ok(c) => signs.push(c),
                    Err(RecurrenceError::Unsupported(_)) => continue 'relative,
                    Err(e) => return Err(e),
                }
            }
            chosen = Some((s, rho, ratio, gamma, signs));
            break 'radius;
        }
    }
    let (s, rho, ratio, gamma, signs) =
        chosen.ok_or_else(|| failure("uniform real-cofactor radius grid exhausted"))?;
    out.signs.extend(signs);
    // Simplicity fixes the real/nonreal pattern on the connected interval.
    // No real root can leave (-rho,rho) without crossing a checked barrier.
    // All real products are <rho^2<s^2<b, so the largest positive product
    // is the unique conjugate-pair norm. Vieta gives the Cauchy weights.
    out.phase.push(format!("at t={}, exact Sturm count is {m} real roots, all in (-{rho},{rho}), polynomial {sample_polynomial:?}; characteristic simplicity and the two nonvanishing real-root barriers preserve this count throughout J. The selected largest positive pair product is >{s}^2; all real-root products have modulus <{rho}^2. Thus |h_j|/|lambda|^j <= binomial({m},j)*({rho}/{s})^j, gamma={gamma:?}",domain.sample));
    let sum: Q = gamma.iter().cloned().sum();
    let amplitude = numerator / (q(1) + &sum);
    out.phase.push(format!("sum gamma={sum}<1 gives |H(1/lambda)-1|<=sum gamma<1, hence Re H>0. N={numerator}>0 gives Re F>0 and |F|>={amplitude}; its continuous argument lies in (-pi/2,pi/2)"));
    for left in [true, false] {
        let (x, identity) = collision_root(characteristic, domain, left)?;
        if endpoint_rational_sign(domain, &x, left)? != if left { -1 } else { 1 } {
            return Err(failure(
                "real-cofactor barrier endpoint is not the required dominant collision",
            ));
        }
        let outside = x.square().sub(&R::of(B::constant(&rho * &rho)));
        if endpoint_rational_sign(domain, &outside, left)? != 1 {
            return Err(failure(
                "endpoint double root is inside the real-cofactor barrier",
            ));
        }
        // All m real cofactor roots have endpoint limits in [-rho,rho].
        // A double root outside that interval must be the conjugate pair's
        // collision. Its norm limit is >=s^2, so the same Vieta/Cauchy bounds
        // extend to this endpoint and still give H(1/x)>0.
        out.phase.push(format!("{} endpoint: {identity}; |x|>{rho} excludes every real-cofactor limit, so this is the dominant pair collision. Its norm limit is >= {s}^2; the same Cauchy weights extend by continuity, H(1/x)>0 and F>0",if left {"left"} else {"right"}));
    }
    let weights = vec![numerator.clone()];
    let (cutoff, tail) = cutoff(order, &ratio, &amplitude, &weights)?;
    set_domination_cutoff(&mut out, cutoff);
    out.tail_at_cutoff = tail.to_string();
    out.automatic_uniform_constants = Some(EffectiveBkwUniformConstants {
        modulus_ratio: ratio.to_string(),
        amplitude_lower_bound: amplitude.to_string(),
        numerator_weights: weights.iter().map(ToString::to_string).collect(),
        cofactor_weights: Some(gamma.iter().map(ToString::to_string).collect()),
    });
    out.domination=format!("uniform_compact_tail: exact real-root barriers give q={ratio}, m={amplitude}, beta={weights:?}; complete-homogeneous Cauchy tail T_N={tail}<m and positive-series tails decrease");
    out.lean_declarations
        .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_tail".into());
    out.phase_limits.left_eigenvalue_argument = "1".into();
    out.phase_limits.right_eigenvalue_argument = "0".into();
    out.phase_limits.left_amplitude_argument = "1/2".into();
    out.phase_limits.right_amplitude_argument = "-1/2".into();
    out.phase.push(
        "arg A=arg F+theta-pi/2 fixes the continuous lift; F is positive at both finite collisions"
            .into(),
    );
    out.expansion_data = EffectiveBkwExpansionData {
        pair_real_part: a.wire(),
        pair_norm_squared: b.wire(),
        cofactor: cofactor.iter().map(R::wire).collect(),
        remainder_root: None,
        remainder_amplitude: None,
    };
    out.expansion=vec!["D=(1-2*a*z+b*z^2)H checked exactly in the pair-product quotient; proper partial fractions give the pair expansion and the complete-homogeneous cofactor remainder".into()];
    super::super::degree::close_phase(&mut out)?;
    Ok(out)
}
