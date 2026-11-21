//! Pointwise real-mode domination and compact uniform fallback suppliers.

use super::cubic::{endpoint_is_zero, endpoint_rational_sign, substitute};
use super::*;

mod barrier;
mod uniform;

fn verify_pair_factorization(
    d: &[P],
    p: &B,
    a: &R,
    b: &R,
    cofactor: &[R],
) -> Result<(), RecurrenceError> {
    let pair = [R::of(c(1)), a.scale(&q(-2)), b.clone()];
    let mut product = vec![R::of(c(0)); d.len()];
    for (i, u) in pair.iter().enumerate() {
        for (j, v) in cofactor.iter().enumerate() {
            product[i + j] = product[i + j].add(&u.mul(v));
        }
    }
    for (j, value) in product.iter().enumerate() {
        if !rem(&value.sub(&R::of(t(&d[j]))).n, p).is_zero() {
            return Err(failure("exact pair/cofactor factorization identity failed"));
        }
    }
    Ok(())
}

fn select_pair(
    d: &[P],
    supplied: Option<&EffectiveBkwPairFactorization>,
    ray: bool,
) -> Result<(B, R, R, Vec<R>, Branch), RecurrenceError> {
    let order = d.len() - 1;
    Ok(if let Some(pair) = supplied {
        let h = decode_b(&pair.cofactor)?;
        if h.index_degree() != Some(order - 2) || h.rows().first() != Some(&vec![Q::one()]) {
            return Err(RecurrenceError::InvalidInput(
                "pair cofactor must have degree d-2 and constant coefficient one".into(),
            ));
        }
        (
            variable(),
            R::of(t(&decode(&pair.real_part)?)),
            R::of(t(&decode(&pair.norm_squared)?)),
            h.rows()
                .iter()
                .map(|r| R::of(t(&P::new(r.clone()))))
                .collect(),
            Branch::Zero,
        )
    } else if order == 2 {
        (
            variable(),
            R::of(t(&d[1])).scale(&(q(-1) / q(2))),
            R::of(t(&d[2])),
            vec![R::of(c(1))],
            Branch::Zero,
        )
    } else if order > 4 {
        higher::algebraic_pair(d)?
    } else {
        let p = pair_product_polynomial(d);
        let b = R::of(variable());
        let a = b
            .mul(&R::of(t(&d[1])).neg().mul(&b).add(&R::of(t(&d[3]))))
            .div(&b.square().sub(&R::of(t(&d[4]))).scale(&q(2)));
        let c0 = R::of(t(&d[4])).div(&b);
        let c1 = R::of(t(&d[1])).add(&a.scale(&q(2)));
        (
            p,
            a,
            b,
            vec![R::of(c(1)), c1, c0],
            if ray {
                Branch::LargestPositive
            } else {
                Branch::UniquePositive
            },
        )
    })
}

fn pair_product_polynomial(d: &[P]) -> B {
    let (b, c, e, f) = (d[1].clone(), d[2].clone(), d[3].clone(), d[4].clone());
    B::new(
        vec![
            f.pow(3),
            -(c.clone() * f.pow(2)),
            b.clone() * e.clone() * f.clone() - f.pow(2),
            c.clone() * f.clone().scale(&q(2)) - e.pow(2) - b.pow(2) * f.clone(),
            b * e - f,
            -c,
            P::one(),
        ]
        .into_iter()
        .map(|p| p.coeffs().to_vec())
        .collect(),
    )
}

#[derive(Clone)]
struct Z {
    re: R,
    im: R,
}
impl Z {
    fn add(&self, r: &Self) -> Self {
        Self {
            re: self.re.add(&r.re),
            im: self.im.add(&r.im),
        }
    }
    fn mul(&self, r: &Self, y2: &R) -> Self {
        Self {
            re: self.re.mul(&r.re).sub(&self.im.mul(&r.im).mul(y2)),
            im: self.re.mul(&r.im).add(&self.im.mul(&r.re)),
        }
    }
    fn norm_squared(&self, y2: &R) -> R {
        self.re.square().add(&self.im.square().mul(y2))
    }
}
fn eval_complex(coefficients: &[R], z: &Z, y2: &R) -> Z {
    coefficients.iter().rev().fold(
        Z {
            re: R::of(B::zero()),
            im: R::of(B::zero()),
        },
        |v, c| {
            v.mul(z, y2).add(&Z {
                re: c.clone(),
                im: R::of(B::zero()),
            })
        },
    )
}

// Reduce the reciprocal evaluation before taking norms. With g=b^2-D_4,
// 2a=an/g, W=lambda^3*N(1/lambda) and V=lambda*h(lambda),
// F=b*(U*lambda+V0)/(L*lambda+M). All four displayed coefficients
// below are polynomials; clearing the common positive g removes artificial
// denominator factors before the resultant calculation.
pub(super) fn quartic_amplitude_forms(n: &B, d: &[P]) -> (R, R, R, R) {
    let b = variable();
    let g = b.power(2).subtract(&t(&d[4]));
    let an = b.multiply(&t(&d[1]).scale(&q(-1)).multiply(&b).add(&t(&d[3])));
    let n0 = t(&P::new(n.rows().first().cloned().unwrap_or_default()));
    let n1 = t(&P::new(n.rows().get(1).cloned().unwrap_or_default()));
    let n2 = t(&P::new(n.rows().get(2).cloned().unwrap_or_default()));
    let n3 = t(&P::new(n.rows().get(3).cloned().unwrap_or_default()));
    let u = n0
        .multiply(&an.power(2).subtract(&b.multiply(&g.power(2))))
        .add(&n1.multiply(&an).multiply(&g))
        .add(&n2.multiply(&g.power(2)));
    let v = n3
        .multiply(&g)
        .subtract(&b.multiply(&n1).multiply(&g))
        .subtract(&b.multiply(&n0).multiply(&an))
        .multiply(&g);
    let l = b
        .multiply(&an)
        .multiply(&t(&d[1]))
        .multiply(&g)
        .add(&b.multiply(&an.power(2)).scale(&q(2)))
        .subtract(&g.power(3));
    let m = b
        .power(2)
        .multiply(&g)
        .multiply(&t(&d[1]).multiply(&g).add(&an.scale(&q(2))))
        .scale(&q(-1));
    let norm = |u: &B, v: &B| {
        g.multiply(&b.multiply(&u.power(2)).add(&v.power(2)))
            .add(&an.multiply(u).multiply(v))
    };
    let nn = b.power(2).multiply(&norm(&u, &v));
    let hn = norm(&l, &m);
    let im = u.multiply(&m).subtract(&v.multiply(&l));
    let re = g
        .multiply(&b.multiply(&u).multiply(&l).add(&v.multiply(&m)))
        .scale(&q(2))
        .add(&an.multiply(&u.multiply(&m).add(&v.multiply(&l))));
    (R::of(nn), R::of(hn), R::of(im), R::of(re))
}

// A sample can reject a proposal cheaply. Only full-interval certificates
// accept one; retain those certificates so the successful checks are not repeated.
struct UniformSearch<'a> {
    polynomial: &'a B,
    domain: &'a Domain,
    branch: Branch,
    sample_branch: I,
}
impl UniformSearch<'_> {
    fn choose(
        &self,
        label: &str,
        candidates: impl IntoIterator<Item = Q>,
        conditions: impl Fn(&Q) -> Vec<(String, B)>,
    ) -> Result<(Q, Vec<EffectiveBkwSignCertificate>), RecurrenceError> {
        for candidate in candidates {
            let conditions = conditions(&candidate);
            if conditions.iter().any(|(_, f)| {
                eval_b(
                    &rem(f, self.polynomial),
                    &I::point(self.domain.sample.clone()),
                    &self.sample_branch,
                )
                .sign()
                    < 0
            }) {
                continue;
            }
            let mut signs = Vec::new();
            let mut rejected = false;
            for (name, f) in conditions {
                match check_sign(&name, self.polynomial, &f, self.domain, self.branch, 0) {
                    Ok(sign) => signs.push(sign),
                    Err(RecurrenceError::Unsupported(_)) => {
                        rejected = true;
                        break;
                    }
                    Err(error) => return Err(error),
                }
            }
            if !rejected {
                return Ok((candidate, signs));
            }
        }
        Err(failure(format!(
            "automatic compact {label} search exhausted its rational grid"
        )))
    }
}

/// Quotient arithmetic may cancel a zero numerator before a later sign
/// condition sees its denominator. Check the stored pair/cofactor fractions
/// first, strictly throughout the selected branch, independently of b-a^2.
fn norm_root_cover(k: &P, known: &[P]) -> P {
    let mut remaining = k.clone();
    let mut retained = P::one();
    // Strip only exact factors, retaining every critical root. The stored
    // resultant remains complete and unchanged.
    for factor in known {
        if !factor.degree().is_some_and(|d| d > 0 && d <= 128) {
            continue;
        }
        let mut used = false;
        while remaining.degree() >= factor.degree() {
            let (quotient, remainder) = remaining.div_rem(factor);
            if !remainder.is_zero() {
                break;
            }
            remaining = quotient;
            used = true;
        }
        if used {
            retained = retained * factor.clone();
        }
    }
    remaining * retained
}

pub(super) fn check_pair_denominators(
    p: &B,
    a: &R,
    cofactor: &[R],
    domain: &Domain,
    which: Branch,
) -> Result<Vec<EffectiveBkwSignCertificate>, RecurrenceError> {
    let monic = p.rows().last() == Some(&vec![q(1)]);
    let mut signs = Vec::new();
    let mut pair_norm = P::one();
    let mut cofactor_norms: Vec<P> = Vec::new();
    for (j, value) in std::iter::once(a).chain(cofactor).enumerate() {
        if !monic {
            if value.d != B::one() {
                signs.push(check_sign(
                    &format!("pair/cofactor denominator {j} is nonzero"),
                    p,
                    &value.d.power(2),
                    domain,
                    which,
                    1,
                )?);
            }
            continue;
        }
        // For monic p, Res(p,rem(uv,p))=Res(p,u)*Res(p,v).
        // Reuse a previous norm only after checking the exact product identity.
        let norm = if value.d == B::one() {
            P::one()
        } else if j >= 2 {
            let k = j - 1;
            let mut product = a.d.multiply(&cofactor[k - 1].d);
            let mut norm = pair_norm.clone() * cofactor_norms[k - 1].clone();
            if k >= 2 {
                product = product.multiply(&cofactor[k - 2].d);
                norm = norm * cofactor_norms[k - 2].clone();
            }
            if rem(&product, p) == value.d {
                norm
            } else {
                resultant(p, &value.d)
            }
        } else {
            resultant(p, &value.d)
        };
        if value.d != B::one() {
            let k = norm.pow(2);
            let mut known = vec![pair_norm.clone(), norm.clone()];
            known.extend(cofactor_norms.iter().cloned());
            let cover = norm_root_cover(&k, &known);
            signs.push(check_sign_with_resultant_cover(
                &format!("pair/cofactor denominator {j} is nonzero"),
                p,
                &value.d.power(2),
                domain,
                which,
                1,
                |_| (k, cover),
            )?);
        }
        if j == 0 {
            pair_norm = norm;
        } else {
            cofactor_norms.push(norm);
        }
    }
    Ok(signs)
}

pub(super) fn certify(
    n: &B,
    d: &[P],
    request: &EffectiveBkwRequest,
    degree: &DegreeBudget,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    if (d.len() > 5
        || (d.len() == 5
            && !legacy
            && degree.data.newton.is_some()
            && degree.count.slope < degree.count.divisor))
        && request.pair_factorization.is_none()
    {
        // Test every quadratic factor: a smaller complex cofactor need not
        // be the dominant pair, nor share the required collision endpoints.
        for pair in higher::polynomial_pairs(d).unwrap_or_default() {
            match certify_with_pair(n, d, request, degree, Some(&pair), legacy) {
                Ok(certificate) => return Ok(certificate),
                Err(RecurrenceError::Unsupported(_)) => {}
                Err(error) => return Err(error),
            }
        }
    }
    certify_with_pair(n, d, request, degree, None, legacy)
}

fn certify_with_pair(
    n: &B,
    d: &[P],
    request: &EffectiveBkwRequest,
    degree: &DegreeBudget,
    automatic: Option<&EffectiveBkwPairFactorization>,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    let order = d.len() - 1;
    let characteristic = B::new(d.iter().rev().map(|p| p.coeffs().to_vec()).collect());
    let supplied = request.pair_factorization.as_ref().or(automatic);
    let disc = if let Some(pair) = supplied {
        let a = decode(&pair.real_part)?;
        let b = decode(&pair.norm_squared)?;
        (a.pow(2) - b).scale(&q(4))
    } else {
        discriminant(&characteristic)
    };
    let ray_candidate = !legacy
        && order == 4
        && degree.data.newton.is_some()
        && degree.count.slope < degree.count.divisor;
    let ray_discriminant_sign = if ray_candidate {
        quadratic::infinity_sign(&disc)
    } else {
        -1
    };
    let domain_disc = if ray_discriminant_sign == 1 {
        disc.scale(&q(-1))
    } else {
        disc.clone()
    };
    let domain = Domain::select(
        &domain_disc,
        if ray_candidate {
            3
        } else if supplied.is_some() {
            2
        } else {
            order
        },
        &request.interval,
    )
    .map_err(RecurrenceError::from)?;
    let ray = domain.lower.is_none() && domain.upper.is_some() && order == 4;
    if !ray && (domain.lower.is_none() || domain.upper.is_none()) {
        return Err(failure(
            "uniform/no-remainder phase supplier currently needs two finite collision endpoints",
        ));
    }
    let (p, a, b, cofactor, which) = select_pair(d, supplied, ray)?;
    let mut out = empty_certificate(d, &p, &domain, degree);
    if matches!(which, Branch::LargestPositive) {
        let branch_discriminant = discriminant(&p);
        if (!ray
            && (branch_discriminant.is_zero()
                || !domain.critical_roots(&branch_discriminant)?.is_empty()))
            || !domain
                .critical_roots(&P::new(p.rows()[0].clone()))?
                .is_empty()
        {
            return Err(failure(
                "higher-order pair-product branches collide or cross zero in the interval",
            ));
        }
        let sign = if branch_discriminant.evaluate(&domain.sample).is_positive() {
            1
        } else {
            -1
        };
        if !ray {
            out.signs.push(check_univariate(
                "pair-product branch separation",
                &branch_discriminant,
                &domain,
                sign,
            )?);
        } else {
            out.phase.push("the simple characteristic roots have a fixed real/nonreal pattern; the largest positive pair product is the maximum of continuous conjugate-pair norms and positive real-root products; any other positive mixed product is bounded by the larger conjugate-pair norm. Thus this maximum is continuous even when off-branch products collide".into());
        }
        let constant = P::new(p.rows()[0].clone());
        let sign = if constant.evaluate(&domain.sample).is_positive() {
            1
        } else {
            -1
        };
        out.signs.push(check_univariate(
            "pair-product branches avoid zero",
            &constant,
            &domain,
            sign,
        )?);
        // This is exactly the path whose R::combine reduces modulo the
        // relation without a gcd. Polynomial suppliers have unit denominators;
        // legacy quartic fractions retain their earlier gcd/sign evidence.
        out.signs
            .extend(check_pair_denominators(&p, &a, &cofactor, &domain, which)?);
    }
    if domain.upper_is_nonpositive() {
        out.root_region = EffectiveBkwRootRegion::NegativeRay;
    }
    out.branch_selection=if order==2||supplied.is_some(){"r=0 (rational pair parameters, no algebraic branch needed)"} else if order>4 {"b is the largest positive pair-product root; monic pair-product resultant, nonvanishing branch discriminant and constant term fix its continuous selection; quadratic subresultant and nonreal-pair test identify the conjugates"}else{"b is the unique positive real pair-product root: quartic discriminant < 0 and constant coefficient < 0 give a nonreal pair and two opposite-sign real cofactor roots"}.into();
    out.signs.push(check_univariate(
        if ray && supplied.is_some() {
            "nonvanishing supplied pair discriminant on ray"
        } else if ray {
            "nonvanishing characteristic discriminant on ray"
        } else {
            "characteristic discriminant < 0"
        },
        &disc,
        &domain,
        if ray { ray_discriminant_sign } else { -1 },
    )?);
    if order == 4 && supplied.is_none() {
        out.signs.push(check_univariate(
            if ray {
                "nonzero characteristic constant on ray"
            } else {
                "opposite-sign cofactor roots"
            },
            &d[4],
            &domain,
            if ray && d[4].evaluate(&domain.sample).is_positive() {
                1
            } else {
                -1
            },
        )?);
    }
    if ray {
        out.branch_selection = if matches!(which,Branch::Zero) {"r=0; exact polynomial pair parameters"} else {"largest positive quartic pair product; nonreal-pair and strict cofactor-modulus tests identify the dominant conjugates"}.into();
        out.signs.push(check_sign(
            "quartic pair denominator b^2-D4",
            &p,
            &b.square().sub(&R::of(t(&d[4]))).signed_polynomial(&p),
            &domain,
            which,
            1,
        )?);
    }
    let y2 = b.sub(&a.square());
    if !legacy
        && order > 4
        && supplied.is_none()
        && request.uniform_constants.is_none()
        && n.index_degree() == Some(0)
        && n.rows()[0].len() == 1
        && n.rows()[0][0].is_positive()
    {
        verify_pair_factorization(d, &p, &a, &b, &cofactor)?;
        match barrier::certify(
            barrier::BarrierData {
                d,
                characteristic: &characteristic,
                p: &p,
                a: &a,
                b: &b,
                cofactor: &cofactor,
                domain: &domain,
                which,
                numerator: &n.rows()[0][0],
            },
            out.clone(),
        ) {
            Ok(certificate) => return Ok(certificate),
            Err(RecurrenceError::Unsupported(_)) => {}
            Err(error) => return Err(error),
        }
    }
    out.signs.push(check_sign(
        "nonreal pair b-a^2",
        &p,
        &y2.signed_polynomial(&p),
        &domain,
        which,
        1,
    )?);
    verify_pair_factorization(d, &p, &a, &b, &cofactor)?;
    let nc: Vec<R> = n
        .rows()
        .iter()
        .map(|r| R::of(t(&P::new(r.clone()))))
        .collect();
    let (nn, hn, imaginary, real) = if order == 4 && supplied.is_none() {
        quartic_amplitude_forms(n, d)
    } else {
        let reciprocal = Z {
            re: a.div(&b),
            im: R::of(c(-1)).div(&b),
        };
        let nv = eval_complex(&nc, &reciprocal, &y2);
        let hv = eval_complex(&cofactor, &reciprocal, &y2);
        (
            nv.norm_squared(&y2),
            hv.norm_squared(&y2),
            nv.im.mul(&hv.re).sub(&nv.re.mul(&hv.im)),
            nv.re.mul(&hv.re).add(&nv.im.mul(&hv.im).mul(&y2)),
        )
    };
    if ray {
        return quartic_ray::certify(
            n,
            d,
            degree,
            &domain,
            &characteristic,
            &p,
            &a,
            &b,
            &cofactor,
            which,
            out,
        );
    }
    let imag_poly = imaginary.signed_polynomial(&p);
    let imaginary_zero = imag_poly.is_zero();
    let mut positive_real_lift = false;
    if !imaginary_zero {
        match check_sign("Im(N/H)/Im(lambda) < 0", &p, &imag_poly, &domain, which, -1) {
            Ok(sign) => out.signs.push(sign),
            Err(RecurrenceError::Unsupported(_)) if !legacy => {
                let real_a = if order == 4 && supplied.is_none() {
                    let g = b.square().sub(&R::of(t(&d[4])));
                    real.scale(&(q(1) / q(2))).add(&a.mul(&g).mul(&imaginary))
                } else {
                    real.add(&a.mul(&imaginary))
                };
                out.signs.push(check_sign(
                    "positive real part of pair amplitude",
                    &p,
                    &real_a.signed_polynomial(&p),
                    &domain,
                    which,
                    1,
                )?);
                positive_real_lift = true;
            }
            Err(error) => return Err(error),
        }
    }
    let mut endpoint_signs = vec![];
    let mut endpoint_cofactors = Vec::new();
    for left in [true, false] {
        if !endpoint_is_zero(&domain, &disc, left)? {
            return Err(failure("phase endpoint is not a characteristic collision"));
        }
        let (x, identity) = if supplied.is_some() {
            (
                a.clone(),
                "supplied pair discriminant is zero and a!=0 at the endpoint".into(),
            )
        } else {
            collision_root(&characteristic, &domain, left)?
        };
        if endpoint_rational_sign(&domain, &x, left)? != if left { -1 } else { 1 } {
            return Err(failure("pair collision does not have theta endpoints pi,0"));
        }
        let inv = R::of(c(1)).div(&x);
        let nr = nc.iter().rev().fold(R::of(c(0)), |v, c| v.mul(&inv).add(c));
        let cr = if supplied.is_some() {
            cofactor
                .iter()
                .rev()
                .fold(R::of(c(0)), |v, c| v.mul(&inv).add(c))
        } else if order == 2 {
            R::of(c(1))
        } else {
            let mut h = vec![R::of(c(1))];
            for j in 1..=order - 2 {
                let mut coefficient = R::of(t(&d[j])).add(&x.scale(&q(2)).mul(&h[j - 1]));
                if j >= 2 {
                    coefficient = coefficient.sub(&x.square().mul(&h[j - 2]));
                }
                h.push(coefficient);
            }
            if order > 4 {
                endpoint_cofactors.push((left, x.clone(), h.clone()));
            }
            h.iter().rev().fold(R::of(c(0)), |v, c| v.mul(&inv).add(c))
        };
        if cr.n.is_zero() {
            return Err(failure("endpoint pair cofactor vanishes identically"));
        }
        let sign = endpoint_rational_sign(&domain, &nr.div(&cr), left)?;
        endpoint_signs.push(sign);
        out.phase.push(format!(
            "{} endpoint: {}; theta -> {}; F=N/H has sign {}",
            if left { "left" } else { "right" },
            identity,
            if left { "pi" } else { "0" },
            sign
        ));
    }
    let winding = if positive_real_lift {
        0
    } else if imaginary_zero {
        if endpoint_signs[0] != endpoint_signs[1] {
            return Err(failure("real amplitude endpoint signs disagree"));
        }
        out.signs.push(check_sign(
            "nonvanishing real amplitude",
            &p,
            &real.signed_polynomial(&p),
            &domain,
            which,
            endpoint_signs[1],
        )?);
        0
    } else {
        match endpoint_signs.as_slice() {
            [-1, 1] => -1,
            [1, -1] => 1,
            _ => 0,
        }
    };
    if !positive_real_lift
        && degree.data.newton.is_none()
        && (degree.count.slope != 1
            || degree.count.divisor != 1
            || degree.count.offsets != [0]
            || !(winding == 0 || (winding == -1 && endpoint_signs[1] == 1)))
    {
        return Err(failure(
            "phase winding and exterior sign do not close the degree budget",
        ));
    }
    if !positive_real_lift {
        out.phase.push(format!("continuous lift of arg F has winding {winding}; (n+1)*theta+arg F supplies n+{winding} sign changes by first half-integer-pi crossings"));
    }
    if winding == -1 && !positive_real_lift {
        out.lean_declarations.push(
            "RealRooted.EffectiveBKW.splits_nodup_of_strictAnti_alternating_of_leftEscape".into(),
        );
        out.exterior.push("right amplitude positive; leftmost sample sign (-1)^(n-1), positive top coefficient gives sign (-1)^n at -infinity, hence one further real zero".into());
    }
    let delta = |sign: i8| {
        if sign > 0 {
            q(0)
        } else if imaginary_zero {
            q(1)
        } else {
            q(-1)
        }
    };
    out.phase_limits.left_amplitude_argument = (q(1) / q(2) + delta(endpoint_signs[0])).to_string();
    out.phase_limits.right_amplitude_argument =
        (-q(1) / q(2) + delta(endpoint_signs[1])).to_string();
    if positive_real_lift {
        out.phase_limits.left_amplitude_argument = (q(endpoint_signs[0].into()) / q(2)).to_string();
        out.phase_limits.right_amplitude_argument =
            (-q(endpoint_signs[1].into()) / q(2)).to_string();
        out.phase.push("Re A>0 fixes the lift in (-pi/2,pi/2); the finite simple pair collisions give these one-sided amplitude arguments, independently of the winding of F".into());
    }
    out.expansion_data = EffectiveBkwExpansionData {
        pair_real_part: a.wire(),
        pair_norm_squared: b.wire(),
        cofactor: cofactor.iter().map(R::wire).collect(),
        remainder_root: None,
        remainder_amplitude: None,
    };
    out.expansion=vec![
        format!("Re(lambda)=({:?})/({:?}); |lambda|^2=({:?})/({:?})",encode_b(&a.n),encode_b(&a.d),encode_b(&b.n),encode_b(&b.d)),
        format!("H(z) coefficients {:?}; D=(1-2*a*z+b*z^2)*H; equality checked in the branch quotient",cofactor.iter().map(|r|(encode_b(&r.n),encode_b(&r.d))).collect::<Vec<_>>()),
        "sin(theta)*|lambda|^(-n)*P_n = Im(exp(i*(n+1)*theta)*F(1/lambda))-R_n; A=F/(1-conj(lambda)/lambda)".into(),
    ];
    let pointwise = if !legacy && order == 4 && request.uniform_constants.is_none() {
        match pointwise::certify(n, &p, &a, &b, &cofactor, &nn, &hn, &domain, which) {
            Ok(evidence) => Some(evidence),
            Err(RecurrenceError::Unsupported(_)) => None,
            Err(error) => return Err(error),
        }
    } else {
        None
    };
    if order == 2 {
        if !legacy {
            set_domination_cutoff(&mut out, 0);
        }
        out.domination = "no_remainder: proper quadratic numerator, R_n=0 for n>=0".into();
        out.lean_declarations
            .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_rest_eq_zero".into());
    } else if let Some((cutoff, signs, method)) = pointwise {
        out.signs.extend(signs);
        set_domination_cutoff(&mut out, cutoff);
        out.domination = format!("{method}: exact cofactor residue norms at N={cutoff} and strict root-modulus gaps give domination for every n>=N");
        if method == "pointwise_two_real_remainders" {
            out.expansion.push("Nbar/q mod h=U*mu+V; Cx+Cy=U, (Cx-Cy)^2=(-s*U+2*V)^2/(s^2-4c); |Cx|+|Cy|=max(|Cx+Cy|,|Cx-Cy|) for real cofactor roots".into());
        } else {
            out.expansion.push("Nbar mod h=u*mu+v; q mod h=e*mu+f; 4|B|^2=(c*u^2+v^2-s*u*v)/((c-s^2/4)*(f^2-s*e*f+c*e^2)); rest_n=2 Re(B*mu^n)".into());
        }
    } else {
        uniform::certify(
            uniform::UniformData {
                order,
                supplied: supplied.is_some(),
                which,
                n,
                p: &p,
                b: &b,
                cofactor: &cofactor,
                nn: &nn,
                hn: &hn,
                nc: &nc,
                domain: &domain,
                endpoint_cofactors: &endpoint_cofactors,
            },
            request,
            &mut out,
            legacy,
        )?;
    }
    if positive_real_lift && out.domination_cutoff.is_none() {
        out.domination_cutoff = Some(out.cutoff);
    }
    if degree.data.newton.is_some() || positive_real_lift {
        super::degree::close_phase(&mut out)?;
    }
    Ok(out)
}

/// A characteristic collision can belong to the cofactor. After deflating
/// its double root, these bounds prove that the double root is dominant,
/// hence it is the endpoint of the selected conjugate pair. The caller has
/// already checked sum gamma_j/q^j < 1 with q<1.
pub(super) fn check_collision_cofactor(
    domain: &Domain,
    x: &R,
    h: &[R],
    gamma: &[Q],
    left: bool,
) -> Result<(), RecurrenceError> {
    if h.len() != gamma.len() + 1 {
        return Err(failure("endpoint cofactor weight count mismatch"));
    }
    for (j, (coefficient, weight)) in h.iter().skip(1).zip(gamma).enumerate() {
        let bound = x
            .pow(2 * (j + 1))
            .scale(&(weight * weight))
            .sub(&coefficient.square());
        if bound.n.index_degree().unwrap_or(0) > 0 || bound.d.index_degree().unwrap_or(0) > 0 {
            return Err(failure(
                "endpoint cofactor condition retained the branch variable",
            ));
        }
        let n = P::new(bound.n.rows().first().cloned().unwrap_or_default());
        let d = P::new(bound.d.rows().first().cloned().unwrap_or_default());
        let denominator_sign = domain.endpoint_sign(&d, left)?;
        if !n.is_zero()
            && !endpoint_is_zero(domain, &n, left)?
            && domain.endpoint_sign(&n, left)? * denominator_sign < 0
        {
            return Err(failure("characteristic endpoint collision is not the dominant pair: endpoint cofactor Cauchy bound fails"));
        }
    }
    Ok(())
}

pub(super) fn collision_root(
    p: &B,
    domain: &Domain,
    left: bool,
) -> Result<(R, String), RecurrenceError> {
    let mut a = p.clone();
    let mut b = derivative(p);
    while b.index_degree().is_some_and(|d| d > 1) {
        let degree = b.index_degree().unwrap();
        let lead = t(&P::new(b.rows().last().unwrap().clone()));
        let mut r = a;
        while r.index_degree().is_some_and(|d| d >= degree) {
            let shift = r.index_degree().unwrap() - degree;
            let lc = t(&P::new(r.rows().last().unwrap().clone()));
            r = r
                .multiply(&lead)
                .subtract(&b.multiply(&lc).multiply(&B::monomial(shift, 0, Q::one())));
        }
        if r.is_zero() {
            return Err(failure("identically repeated characteristic roots"));
        }
        let content = r
            .rows()
            .iter()
            .fold(P::zero(), |g, row| g.gcd(&P::new(row.clone())));
        r = B::new(
            r.rows()
                .iter()
                .map(|row| P::new(row.clone()).exact_div(&content).coeffs().to_vec())
                .collect(),
        );
        a = b;
        b = r;
    }
    if b.index_degree() != Some(1) {
        return Err(failure("collision has no linear subresultant"));
    }
    let x = R::new(
        t(&(-P::new(b.rows()[0].clone()))),
        t(&P::new(b.rows()[1].clone())),
    );
    for f in [p.clone(), derivative(p)] {
        let value = substitute(&f, &x);
        let np = P::new(value.n.rows().first().cloned().unwrap_or_default());
        if !endpoint_is_zero(domain, &np, left)? {
            return Err(failure("linear subresultant does not give a common characteristic/derivative root at endpoint"));
        }
        endpoint_rational_sign(domain, &R::of(value.d.clone()), left)?;
    }
    let second = substitute(&derivative(&derivative(p)), &x);
    endpoint_rational_sign(domain, &second, left)?;
    Ok((x,format!("linear subresultant {:?}; double root -constant/linear, p=p'=0 and p''!=0 checked at the algebraic endpoint",encode_b(&b))))
}

/// Evaluate the exact positive-series tail used by the uniform supplier.
///
/// `order>=3`, `0<ratio<1` and nonnegative weights are required.
/// Returns an error for invalid data or if no cutoff <=4096 is found.
pub fn uniform_tail_cutoff(
    order: usize,
    ratio: &Q,
    amplitude: &Q,
    weights: &[Q],
) -> Result<(usize, Q), RecurrenceError> {
    if !(3..=32).contains(&order)
        || ratio <= &Q::zero()
        || ratio >= &Q::one()
        || amplitude <= &Q::zero()
        || weights.is_empty()
        || weights.len() > 32
        || weights.iter().any(|w| w.is_negative())
    {
        return Err(RecurrenceError::InvalidInput(
            "invalid uniform tail input".into(),
        ));
    }
    cutoff(order, ratio, amplitude, weights)
}
fn cutoff(
    order: usize,
    ratio: &Q,
    amplitude: &Q,
    weights: &[Q],
) -> Result<(usize, Q), RecurrenceError> {
    let r = order - 2;
    for n in 1.max(weights.len().saturating_sub(2))..=4096 {
        let mut tail = Q::zero();
        for (j, w) in weights.iter().enumerate() {
            let m = n + 2 - j;
            let mut sum = Q::zero();
            for k in 0..r {
                let binomial = crate::numeric::binomial_bigint(m + r - 1, r - 1 - k);
                sum += Q::from_integer(binomial) * ratio.pow(k as i32)
                    / (Q::one() - ratio).pow((k + 1) as i32);
            }
            tail += w * ratio.pow(m as i32) * sum;
        }
        if &tail < amplitude {
            return Ok((n, tail));
        }
    }
    Err(failure(
        "uniform tail cutoff exceeds the exact search guard (4096)",
    ))
}
