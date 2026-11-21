//! Uniform Cauchy bounds and exact tail cutoffs for the compact supplier.

use super::*;

pub(super) struct UniformData<'a> {
    pub order: usize,
    pub supplied: bool,
    pub which: Branch,
    pub n: &'a B,
    pub p: &'a B,
    pub b: &'a R,
    pub cofactor: &'a [R],
    pub nn: &'a R,
    pub hn: &'a R,
    pub nc: &'a [R],
    pub domain: &'a Domain,
    pub endpoint_cofactors: &'a [(bool, R, Vec<R>)],
}

pub(super) fn certify(
    data: UniformData<'_>,
    request: &EffectiveBkwRequest,
    out: &mut EffectiveBkwCertificate,
    legacy: bool,
) -> Result<(), RecurrenceError> {
    let UniformData {
        order,
        supplied,
        which,
        n,
        p,
        b,
        cofactor,
        nn,
        hn,
        nc,
        domain,
        endpoint_cofactors,
    } = data;
    let mut automatic_cofactor_weights = None;
    let mut checked_cofactor_weights = None;
    let (ratio, amplitude, weights) = if let Some(constants) = &request.uniform_constants {
        let (ratio, amplitude, weights) = parse_uniform_constants(constants, n.rows().len())?;
        if supplied || order > 4 {
            let bounds = constants.cofactor_weights.as_ref().ok_or_else(|| {
                failure("supplied compact factorization requires cofactor weights")
            })?;
            if bounds.len() != order - 2 {
                return Err(RecurrenceError::InvalidInput(
                    "cofactor weight count mismatch".into(),
                ));
            }
            let mut budget = Q::zero();
            let mut gammas = Vec::new();
            for (j, (weight, h)) in bounds.iter().zip(cofactor.iter().skip(1)).enumerate() {
                let j = j + 1;
                let gamma = rational(weight)?;
                if gamma.is_negative() {
                    return Err(RecurrenceError::InvalidInput(
                        "negative cofactor weight".into(),
                    ));
                }
                budget += &gamma / ratio.pow(j as i32);
                gammas.push(gamma.clone());
                let f = b.pow(j).scale(&(&gamma * &gamma)).sub(&h.square());
                out.signs.push(check_sign(
                    &format!("Cauchy cofactor weight {j}"),
                    p,
                    &f.signed_polynomial(p),
                    domain,
                    which,
                    0,
                )?);
            }
            if budget >= Q::one() {
                return Err(failure("cofactor Cauchy budget sum gamma_j/q^j is not < 1"));
            }
            checked_cofactor_weights = Some(gammas);
        } else {
            let gap = b.scale(&(&ratio * &ratio)).add(&cofactor[2]);
            for (name, f) in [
                ("uniform modulus gap, part one", gap.clone()),
                (
                    "uniform modulus gap, part two",
                    gap.square()
                        .sub(&b.mul(&cofactor[1].square()).scale(&(&ratio * &ratio))),
                ),
            ] {
                out.signs.push(check_sign(
                    name,
                    p,
                    &f.signed_polynomial(p),
                    domain,
                    which,
                    0,
                )?);
            }
        }
        out.signs.push(check_sign(
            "uniform amplitude lower bound",
            p,
            &nn.sub(&hn.scale(&(&amplitude * &amplitude)))
                .signed_polynomial(p),
            domain,
            which,
            0,
        )?);
        for (j, (weight, coefficient)) in weights.iter().zip(nc).enumerate() {
            let f = b
                .pow(j)
                .scale(&(weight * weight))
                .sub(&coefficient.square());
            out.signs.push(check_sign(
                &format!("numerator weight {j}"),
                p,
                &f.signed_polynomial(p),
                domain,
                which,
                0,
            )?);
        }
        (ratio, amplitude, weights)
    } else {
        let search = UniformSearch {
            polynomial: p,
            domain,
            branch: which,
            sample_branch: branch_box(p, &domain.sample, which, &(q(1) / q(1 << 32)))?,
        };
        let (ratio, signs) = if supplied || order > 4 {
            let mut gammas = Vec::new();
            for (j, h) in cofactor.iter().enumerate().skip(1) {
                if rem(&h.n, p).is_zero() {
                    // The already checked denominator is nonzero. An
                    // identically zero coefficient has the exact weak
                    // Cauchy bound gamma=0; no strict sign is requested.
                    gammas.push(q(0));
                    out.signs.push(check_sign(
                        &format!("Cauchy cofactor weight {j}"),
                        p,
                        &B::zero(),
                        domain,
                        which,
                        0,
                    )?);
                    continue;
                }
                let (gamma, signs) = search.choose(
                    "cofactor weight",
                    std::iter::once(q(0)).chain((-16..=16).map(|k| {
                        if k < 0 {
                            q(1) / q(1 << -k)
                        } else {
                            q(1 << k)
                        }
                    })),
                    |gamma| {
                        vec![(
                            format!("Cauchy cofactor weight {j}"),
                            b.pow(j)
                                .scale(&(gamma * gamma))
                                .sub(&h.square())
                                .signed_polynomial(p),
                        )]
                    },
                )?;
                gammas.push(gamma);
                out.signs.extend(signs);
            }
            let ratio = (1..16)
                .map(|k| q(k) / q(16))
                .find(|ratio| {
                    gammas
                        .iter()
                        .enumerate()
                        .map(|(j, gamma)| gamma / ratio.pow((j + 1) as i32))
                        .sum::<Q>()
                        < Q::one()
                })
                .ok_or_else(|| {
                    failure(
                        "automatic cofactor Cauchy budget is not < 1 on the rational ratio grid",
                    )
                })?;
            automatic_cofactor_weights = Some(gammas.iter().map(ToString::to_string).collect());
            checked_cofactor_weights = Some(gammas);
            (ratio, Vec::new())
        } else {
            search.choose("modulus", (1..16).map(|k| q(k) / q(16)), |ratio| {
                let ratio_squared = ratio * ratio;
                let gap = b.scale(&ratio_squared).add(&cofactor[2]);
                vec![
                    (
                        "uniform modulus gap, part one".into(),
                        gap.signed_polynomial(p),
                    ),
                    (
                        "uniform modulus gap, part two".into(),
                        gap.square()
                            .sub(&b.mul(&cofactor[1].square()).scale(&ratio_squared))
                            .signed_polynomial(p),
                    ),
                ]
            })?
        };
        out.signs.extend(signs);
        let (amplitude, signs) = search.choose(
            "amplitude",
            (0..=16).map(|k| q(1) / q(1 << k)),
            |amplitude| {
                vec![(
                    "uniform amplitude lower bound".into(),
                    nn.sub(&hn.scale(&(amplitude * amplitude)))
                        .signed_polynomial(p),
                )]
            },
        )?;
        out.signs.extend(signs);
        let mut weights = Vec::new();
        for (j, coefficient) in nc.iter().enumerate() {
            let (weight, signs) = search.choose(
                "numerator weight",
                std::iter::once(q(0)).chain((0..=16).map(|k| q(1 << k))),
                |weight| {
                    vec![(
                        format!("numerator weight {j}"),
                        b.pow(j)
                            .scale(&(weight * weight))
                            .sub(&coefficient.square())
                            .signed_polynomial(p),
                    )]
                },
            )?;
            weights.push(weight);
            out.signs.extend(signs);
        }
        (ratio, amplitude, weights)
    };
    if matches!(which, Branch::LargestPositive) {
        let gamma = checked_cofactor_weights
            .as_ref()
            .ok_or_else(|| failure("missing endpoint cofactor bounds"))?;
        for (left, x, h) in endpoint_cofactors {
            check_collision_cofactor(domain, x, h, gamma, *left)?;
            out.phase.push(format!("{} double root is the dominant pair collision: endpoint cofactor satisfies the same Cauchy weights and strict ratio budget q={ratio}",if *left {"left"} else {"right"}));
        }
    }
    let (cutoff, tail) = cutoff(order, &ratio, &amplitude, &weights)?;
    out.cutoff = cutoff;
    if !legacy {
        out.domination_cutoff = Some(cutoff);
    }
    out.tail_at_cutoff = tail.to_string();
    if let Some(gamma) = automatic_cofactor_weights {
        out.automatic_uniform_constants = Some(EffectiveBkwUniformConstants {
            modulus_ratio: ratio.to_string(),
            amplitude_lower_bound: amplitude.to_string(),
            numerator_weights: weights.iter().map(ToString::to_string).collect(),
            cofactor_weights: Some(gamma),
        });
    }
    let weight_strings: Vec<_> = weights.iter().map(ToString::to_string).collect();
    out.domination=format!("uniform_compact_tail: q={ratio}, m={amplitude}, beta={weight_strings:?}; T_n=sum_j beta_j sum_(k>=n+2-j) binomial(k+d-3,d-3) q^k; T_N<m and positive-series tails decrease");
    out.lean_declarations
        .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_tail".into());
    Ok(())
}
