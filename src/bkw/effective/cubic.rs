//! Exact single-remainder supplier, including an edge at negative infinity.

use super::*;

pub(super) fn certify(
    n: &B,
    d: &[P],
    request: &EffectiveBkwRequest,
    degree: &DegreeBudget,
    legacy: bool,
) -> Result<EffectiveBkwCertificate, ChartError> {
    let p = B::new(d.iter().rev().map(|c| c.coeffs().to_vec()).collect());
    let disc = discriminant(&p);
    let domain = Domain::select(&disc, 3, &request.interval)?;
    if domain.lower.is_some() || domain.upper.is_none() {
        return Err(failure(
            "single-remainder phase supplier currently requires (-infinity, finite collision)",
        )
        .into());
    }
    let mut out = empty_certificate(d, &p, &domain, degree);
    if domain.upper_is_nonpositive() {
        out.root_region = EffectiveBkwRootRegion::NegativeRay;
    }
    out.branch_selection = "the unique real root r(t) of the negative-discriminant cubic".into();
    out.signs.push(check_univariate(
        "cubic discriminant < 0",
        &disc,
        &domain,
        -1,
    )?);
    let r = variable();
    let r2 = r.power(2);
    let a = t(&(-d[1].clone())).subtract(&r); // 2 Re lambda
    let b = r2.add(&t(&d[1]).multiply(&r)).add(&t(&d[2]));
    let y = b.scale(&q(4)).subtract(&a.power(2)); // 4 (Im lambda)^2
    let pr = derivative(&p);
    let n0 = t(&P::new(n.rows().first().cloned().unwrap_or_default()));
    let n1 = t(&P::new(n.rows().get(1).cloned().unwrap_or_default()));
    let n2 = t(&P::new(n.rows().get(2).cloned().unwrap_or_default()));
    let row1 = n1.subtract(&t(&d[1]).multiply(&n0));
    let cn = n0.multiply(&r2).add(&n1.multiply(&r)).add(&n2);
    let re = n0.multiply(&pr).subtract(&cn); // 2 pr Re A
    let im = re
        .multiply(&a)
        .subtract(&pr.multiply(&row1).scale(&q(2)))
        .add(&cn.multiply(&r).scale(&q(2)));
    // pr^2 Y (4|A|^2-C^2), with Y=4 (Im lambda)^2.
    let k = re
        .power(2)
        .subtract(&cn.power(2))
        .multiply(&y)
        .add(&im.power(2));
    for (name, f, want) in [
        ("pair imaginary part squared", y.clone(), 1),
        ("simple real root derivative", pr.clone(), 1),
        ("strict modulus gap b-r^2", b.subtract(&r2), 1),
        ("amplitude real part", re.clone(), 1),
        ("amplitude imaginary part", im.clone(), -1),
        (
            "single remainder amplitude dominance 4|A|^2-C^2",
            k.clone(),
            0,
        ),
    ] {
        if (name == "amplitude real part" && !legacy)
            || (degree.data.newton.is_some() && name == "amplitude imaginary part")
        {
            continue;
        }
        out.signs
            .push(check_sign(name, &p, &f, &domain, Branch::UniqueReal, want)?);
    }
    let mut positive_real_lift = true;
    if !legacy {
        let real_lift = check_sign(
            "amplitude real part",
            &p,
            &re,
            &domain,
            Branch::UniqueReal,
            1,
        );
        positive_real_lift = real_lift.is_ok();
        out.signs.push(match real_lift {
            Ok(sign) => sign,
            Err(RecurrenceError::Unsupported(_)) if degree.data.newton.is_some() => check_sign(
                "amplitude lower half-plane",
                &p,
                &im,
                &domain,
                Branch::UniqueReal,
                -1,
            )?,
            Err(error) => return Err(error.into()),
        });
    }
    if !legacy {
        out.domination_cutoff = Some(1);
    }
    // Strict coefficient domination also covers n=0. The phase/degree
    // argument can start later, independently of this pointwise bound.
    if let Ok(sign) = if legacy {
        Err(failure("legacy single-mode supplier"))
    } else {
        check_sign(
            "strict single remainder amplitude dominance 4|A|^2-C^2",
            &p,
            &k,
            &domain,
            Branch::UniqueReal,
            1,
        )
    } {
        out.signs.push(sign);
        set_domination_cutoff(&mut out, 0);
    }
    out.expansion_data = EffectiveBkwExpansionData {
        pair_real_part: R::of(a.clone()).scale(&(q(1) / q(2))).wire(),
        pair_norm_squared: R::of(b.clone()).wire(),
        cofactor: vec![R::of(c(1)).wire(), R::of(r.scale(&q(-1))).wire()],
        remainder_root: Some(R::of(r.clone()).wire()),
        remainder_amplitude: Some(R::new(cn.clone(), pr.clone()).wire()),
    };
    out.phase_limits.left_amplitude_argument = "0".into();
    out.phase_limits.right_amplitude_argument = "-1/2".into();
    out.expansion=vec![
        "lambda=( -D_1-r + i*sqrt(Y) )/2; b=r^2+D_1*r+D_2; Y=4*b-(-D_1-r)^2".into(),
        format!("C=({:?})/p'(r)",encode_b(&cn)),
        format!("Re A=({:?})/(2*p'(r)); Im A=({:?})/(2*p'(r)*sqrt(Y))",encode_b(&re),encode_b(&im)),
        "P_n=2*Re(A*lambda^n)+C*r^n for every n>=0 (proper rational partial fractions, including r=0 by continuity)".into(),
    ];
    // A simple collision and a generic triple point have different phase
    // offsets. At a triple point, delta^3 = p_t*(beta-t), so on the
    // dominant side p_t<0 and arg(delta_pair)=pi/3. The residue denominator
    // has argument 2*pi/3, giving arg A=-2*pi/3 when Nbar(r0)>0.
    let delta = d[1].clone().pow(2) - d[2].scale(&q(3));
    let xn = d[3].scale(&q(9)) - d[1].clone() * d[2].clone();
    if !endpoint_is_zero(&domain, &disc, false)? {
        return Err(failure("right endpoint is not a characteristic collision").into());
    }
    let triple = endpoint_is_zero(&domain, &delta, false)?;
    let right_alpha = if triple {
        let root = R::of(t(&d[1].scale(&(q(-1) / q(3)))));
        let cubic_identity = d[3].clone() - d[1].pow(3).scale(&(q(1) / q(27)));
        let pt = B::new(
            p.rows()
                .iter()
                .map(|row| P::new(row.clone()).derivative().coeffs().to_vec())
                .collect(),
        );
        if degree.data.newton.is_none()
            || !endpoint_is_zero(&domain, &cubic_identity, false)?
            || endpoint_rational_sign(&domain, &root, false)? != 1
            || endpoint_rational_sign(&domain, &substitute(&pt, &root), false)? != -1
            || endpoint_rational_sign(&domain, &substitute(&cn, &root), false)? != 1
            || positive_real_lift
        {
            return Err(failure("generic positive triple edge requires p_t<0, Nbar>0 and the lower-half-plane amplitude lift").into());
        }
        out.phase.push("generic positive triple point: p=(mu-r0)^3 at beta, p_t<0 and Nbar(r0)>0; implicit Puiseux roots delta=c*(beta-t)^(1/3)*omega+O((beta-t)^(2/3)), omega=-1,exp(+-i*pi/3); theta->0, arg A->-2*pi/3 in the certified lower-half-plane lift".into());
        -q(2) / q(3)
    } else {
        if domain.endpoint_sign(&delta, false)? != 1 || domain.endpoint_sign(&xn, false)? != 1 {
            return Err(failure("right endpoint is not a positive simple pair collision").into());
        }
        let x = R::new(t(&xn), t(&delta.scale(&q(2))));
        let simple = R::of(t(&(-d[1].clone()))).sub(&x.scale(&q(2)));
        let ie = substitute(&im, &simple);
        let pe = substitute(&pr, &simple);
        let sign = endpoint_rational_sign(&domain, &ie, false)?;
        if sign == 0
            || (!positive_real_lift && sign != -1)
            || (degree.data.newton.is_none() && sign != -1)
            || endpoint_rational_sign(&domain, &pe, false)? != 1
        {
            return Err(
                failure("finite-edge amplitude limit is not certified in its half-plane").into(),
            );
        }
        q(sign.into()) / q(2)
    };
    out.phase_limits.right_amplitude_argument = right_alpha.to_string();
    if let Some(newton) = &degree.data.newton {
        let (theta, alpha, evidence) = newton_infinity_edge(
            &p,
            &a,
            &re,
            &im,
            &y,
            newton.weight_numerator,
            newton.weight_denominator,
        )?;
        out.phase_limits.left_eigenvalue_argument = theta.to_string();
        out.phase_limits.left_amplitude_argument = alpha.to_string();
        out.phase.extend(evidence);
        if legacy {
            let edge_imaginary_sign = right_alpha * q(2);
            out.phase.push(format!("Re A>0 fixes the lift in (-pi/2,pi/2); finite positive simple collision has arg A -> {edge_imaginary_sign}*pi/2"));
        } else {
            out.phase.push(format!("certified amplitude half-plane fixes the continuous lift; finite edge arg A -> {right_alpha}*pi"));
        }
    } else {
        out.phase.push(format!("right collision: x=({:?})/(2*{:?})>0; p'(r)>0; Im A -> -infinity, Re A bounded; phase -> -pi/2",encode(&xn),encode(&delta)));
        let edge = infinity_edge(&p, &a, &re, &im, &y)?;
        out.phase.extend(edge);
        out.phase.push("arg A in (-pi/2,0); theta in (0,pi); first crossings of j*pi, j=0,...,n-1, give n-1 alternating sign changes".into());
        out.exterior.push("at the leftmost sample the sign is (-1)^(n-1). For positive geometric top coefficients, sign P_n(-infinity)=(-1)^n supplies one extra zero. For the parity top GF it supplies one extra exactly when n is even. Odd degree <= n-1 closes without it".into());
    }
    out.domination="single_real_remainder: 4|A|^2>=C^2 and |r|<|lambda| give strict domination for n>=1, also when C=0; no uniform modulus gap is needed at infinity".into();
    if out.domination_cutoff == Some(0) {
        out.domination = "pointwise_single_real_remainder: 4|A|^2>C^2 and |r|<|lambda| imply |C*r^n|<2|A*lambda^n| for every n>=0".into();
    }
    out.lean_declarations
        .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_single_mode".into());
    out.lean_declarations.push(
        "RealRooted.EffectiveBKW.splits_nodup_of_strictAnti_alternating_of_leftEscape".into(),
    );
    out.lean_declarations
        .push("RealRooted.EffectiveBKW.eventually_phase_gt_atBot".into());
    if degree.data.newton.is_some() {
        super::degree::close_phase(&mut out)?;
    }
    Ok(out)
}

pub(super) fn newton_scaled(
    f: &B,
    growth: usize,
    weight: usize,
    eigen_scale: &Q,
    parameter_scale: &Q,
) -> (B, usize) {
    let m = f
        .rows()
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(move |(j, _)| growth * i + weight * j)
        })
        .max()
        .unwrap_or(0);
    let mut out = B::zero();
    for (i, row) in f.rows().iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if !c.is_zero() {
                out = out.add(&B::monomial(
                    i,
                    m - growth * i - weight * j,
                    c * eigen_scale.pow(i as i32)
                        * parameter_scale.pow(j as i32)
                        * q(if j % 2 == 0 { 1 } else { -1 }),
                ));
            }
        }
    }
    (out, m)
}

pub(super) fn truncated_evaluate(f: &B, u: &P) -> P {
    f.rows().iter().rev().fold(P::zero(), |v, row| {
        let value = v * u.clone() + P::new(row.clone());
        P::new(value.coeffs().iter().take(17).cloned().collect())
    })
}

/// Normalize any monic binomial edge u^d-c with c>0 and coprime (d,k),
/// where k=d*growth/weight. Choose 0<=beta<d with 1+k*beta=d*alpha;
/// then r=c^alpha*u*s^-growth, t=-c^beta*s^-weight gives u^d-1.
/// This is Bezout's identity for arbitrary d,k, not a cubic-only formula.
/// If gcd(d,k)>1 there is no such integer-power dilation. Keep the original
/// scaling and let exact factorization find a rational root if one exists.
/// In the current cubic caller this remaining case has 3|k; any rational
/// normalization would make c=(eigen_scale/parameter_scale^(k/3))^3, so it
/// exists iff c is already a rational cube. Higher-degree noncoprime edges
/// would need a separate perfect-power normalization.
pub(super) fn binomial_newton_dilation(p: &B, growth: usize, weight: usize) -> (Q, Q) {
    if weight == 0 {
        return (q(1), q(1));
    }
    let (initial, _) = newton_scaled(p, growth, weight, &q(1), &q(1));
    let limit = P::new(initial.evaluate_variable(&Q::zero()));
    let Some(degree) = limit.degree().filter(|d| *d >= 2) else {
        return (q(1), q(1));
    };
    if !limit.coeff(0).is_negative()
        || !limit.coeff(degree).is_positive()
        || (1..degree).any(|j| !limit.coeff(j).is_zero())
        || degree * growth % weight != 0
    {
        return (q(1), q(1));
    }
    let k = degree * growth / weight;
    let Some(beta) = (0..degree).find(|beta| (1 + k * beta) % degree == 0) else {
        return (q(1), q(1));
    };
    let alpha = (1 + k * beta) / degree;
    let scale = -limit.coeff(0) / limit.coeff(degree);
    (scale.pow(alpha as i32), scale.pow(beta as i32))
}

fn newton_infinity_edge(
    p: &B,
    a: &B,
    re: &B,
    im: &B,
    y: &B,
    growth: usize,
    weight: usize,
) -> Result<(Q, Q, Vec<String>), RecurrenceError> {
    let (eigen_scale, parameter_scale) = binomial_newton_dilation(p, growth, weight);
    let (scaled_p, _) = newton_scaled(p, growth, weight, &eigen_scale, &parameter_scale);
    let limit = P::new(scaled_p.evaluate_variable(&Q::zero()));
    let factors = crate::factorization::irreducible_factorization_rational_coeffs(limit.coeffs())
        .map_err(|e| failure(format!("infinite-edge factorization: {e}")))?;
    let r0 = factors
        .factors
        .iter()
        .filter(|(_, m)| *m == 1)
        .find_map(|(f, _)| {
            if f.len() == 2 {
                Some(-Q::new(f[0].clone(), f[1].clone()))
            } else {
                None
            }
        })
        .ok_or_else(|| failure("Newton infinite-edge simple real root is not rational"))?;
    let derivative = limit.derivative().evaluate(&r0);
    if derivative.is_zero() {
        return Err(failure("Newton infinite-edge root is not simple"));
    }
    let mut u = P::new(vec![r0.clone()]);
    for j in 1..=16 {
        let coefficient = -truncated_evaluate(&scaled_p, &u).coeff(j) / &derivative;
        let mut row = u.coeffs().to_vec();
        row.resize(j + 1, Q::zero());
        row[j] = coefficient;
        u = P::new(row);
    }
    let leading = |f: &B| -> Result<(i64, Q), RecurrenceError> {
        let (scaled, m) = newton_scaled(f, growth, weight, &eigen_scale, &parameter_scale);
        let series = truncated_evaluate(&scaled, &u);
        let order = series
            .coeffs()
            .iter()
            .position(|c| !c.is_zero())
            .ok_or_else(|| failure("Newton infinity leading order exceeds sixteen Taylor terms"))?;
        Ok((order as i64 - m as i64, series.coeff(order)))
    };
    let (oa, ca) = leading(a)?;
    let (oy, cy) = leading(y)?;
    let (or, cr) = leading(re)?;
    let (oi, ci) = leading(im)?;
    if !cy.is_positive() {
        return Err(failure(
            "Newton infinity signs disagree with the positive branch",
        ));
    }
    let cosine_squared = if 2 * oa > oy {
        q(0)
    } else if 2 * oa == oy {
        &ca * &ca / (&ca * &ca + &cy)
    } else {
        q(1)
    };
    let theta = super::quadratic::angle_from_cosine_squared(
        &cosine_squared,
        if ca.is_positive() { 1 } else { -1 },
    )?;
    let exponent = 2 * oi - oy - 2 * or;
    let alpha = if exponent > 0 {
        if cr.is_positive() {
            q(0)
        } else {
            q(if ci.is_positive() { 1 } else { -1 })
        }
    } else if exponent < 0 {
        q(if ci.is_positive() { 1 } else { -1 }) / q(2)
    } else {
        super::quadratic::angle_from_cosine_squared(
            &(&cy * &cr * &cr / (&cy * &cr * &cr + &ci * &ci)),
            if cr.is_positive() { 1 } else { -1 },
        )? * q(if ci.is_positive() { 1 } else { -1 })
    };
    Ok((theta.clone(),alpha.clone(),vec![format!(
        "t=-{parameter_scale}*s^(-{weight}), r={eigen_scale}*u*s^(-{growth}); scaled characteristic {:?}; simple limit root {r0}, derivative {derivative}; implicit Taylor coefficients {:?}; leading orders a={oa},Y={oy},Re numerator={or},Im numerator={oi}; exact limits theta={theta}*pi, arg A={alpha}*pi",encode_b(&scaled_p),encode(&u))]))
}

pub(super) fn substitute(f: &B, x: &R) -> R {
    f.rows().iter().rev().fold(R::of(B::zero()), |v, row| {
        v.mul(x).add(&R::of(t(&P::new(row.clone()))))
    })
}
pub(super) fn endpoint_rational_sign(
    domain: &Domain,
    r: &R,
    left: bool,
) -> Result<i8, RecurrenceError> {
    if r.n.index_degree().unwrap_or(0) > 0 || r.d.index_degree().unwrap_or(0) > 0 {
        return Err(failure("endpoint substitution retained branch variable"));
    }
    let p = |b: &B| P::new(b.rows().first().cloned().unwrap_or_default());
    let denominator = domain.endpoint_sign(&p(&r.d), left)?;
    let numerator = domain.endpoint_sign(&p(&r.n), left)?;
    Ok(denominator * numerator)
}
pub(super) fn endpoint_is_zero(
    domain: &Domain,
    p: &P,
    left: bool,
) -> Result<bool, RecurrenceError> {
    let e = if left {
        domain.lower.as_ref()
    } else {
        domain.upper.as_ref()
    }
    .ok_or_else(|| failure("finite collision endpoint required"))?;
    let g = p.gcd(&e.polynomial);
    if e.bounds.lo == e.bounds.hi {
        Ok(g.evaluate(&e.bounds.lo).is_zero())
    } else {
        let s =
            crate::sturm::SturmChain::from_bigint_coeffs(&vec_poly::clear_denominators(g.coeffs()));
        Ok(s.count_roots_in(&e.bounds.lo, &e.bounds.hi) == 1)
    }
}
fn scaled(f: &B) -> (B, usize) {
    let m = f
        .rows()
        .iter()
        .enumerate()
        .flat_map(|(i, row)| {
            row.iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .map(move |(j, _)| i + j)
        })
        .max()
        .unwrap_or(0);
    let mut out = B::zero();
    for (i, row) in f.rows().iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if !c.is_zero() {
                out = out.add(&B::monomial(i, m - i - j, c.clone()));
            }
        }
    }
    (out, m)
}
fn series_evaluate(f: &B, u: &P) -> P {
    f.rows()
        .iter()
        .rev()
        .fold(P::zero(), |v, r| v * u.clone() + P::new(r.clone()))
}
fn infinity_edge(p: &B, a: &B, re: &B, im: &B, y: &B) -> Result<Vec<String>, RecurrenceError> {
    let (scaled_p, _) = scaled(p);
    let limit = P::new(scaled_p.evaluate_variable(&Q::zero()));
    let double = limit.gcd(&limit.derivative());
    if double.degree() != Some(1) {
        return Err(failure(
            "infinite-edge scaled characteristic has no rational double root",
        ));
    }
    let pair = -double.coeff(0) / double.coeff(1);
    let simple = limit.exact_div(&double.pow(2));
    if simple.degree() != Some(1) {
        return Err(failure("infinite-edge third root is not simple"));
    }
    let r0 = -simple.coeff(0) / simple.coeff(1);
    if !pair.is_positive() || !r0.is_negative() {
        return Err(failure("infinite-edge root orientation is not theta -> pi"));
    }
    let derivative = limit.derivative().evaluate(&r0);
    let mut u = P::new(vec![r0.clone()]);
    for j in 1..=8 {
        let residual = series_evaluate(&scaled_p, &u);
        let coefficient = -residual.coeff(j) / &derivative;
        let mut row = u.coeffs().to_vec();
        row.resize(j + 1, Q::zero());
        row[j] = coefficient;
        u = P::new(row);
    }
    let order = |f: &B| -> Result<(i64, P), RecurrenceError> {
        let (s, m) = scaled(f);
        let e = series_evaluate(&s, &u);
        let v = e
            .coeffs()
            .iter()
            .position(|c| !c.is_zero())
            .ok_or_else(|| failure("edge expansion requires a higher series order"))?;
        if v > 8 {
            return Err(failure(
                "edge leading order exceeds the checked implicit expansion",
            ));
        }
        Ok((
            v as i64 - m as i64,
            P::new(e.coeffs().iter().take(9).cloned().collect()),
        ))
    };
    let (oi, si) = order(im)?;
    let (oy, sy) = order(y)?;
    let (or, sr) = order(re)?;
    let exponent = 2 * oi - oy - 2 * or;
    if exponent <= 0 {
        return Err(failure(
            "infinite-edge arg A -> 0 is not certified by the implicit series",
        ));
    }
    // The pair's real part divided by t tends to the positive double root.
    let (sa, ma) = scaled(a);
    if ma != 1 || series_evaluate(&sa, &u).coeff(0) != q(2) * &pair {
        return Err(failure("infinite-edge pair scaling identity failed"));
    }
    Ok(vec![
        format!("s=1/t, u=r/t: scaled characteristic {:?}; at s=0 pair root {} (double), third root {} (simple), so theta -> pi",encode_b(&scaled_p),pair,r0),
        format!("implicit simple-root Taylor coefficients through s^8: {:?}; derivative {} != 0",encode(&u),derivative),
        format!("scaled Im numerator {:?}, Y {:?}, Re numerator {:?}; ord((Im A/Re A)^2)={} > 0, hence arg A -> 0 within quadrant IV",encode(&si),encode(&sy),encode(&sr),exponent),
    ])
}
