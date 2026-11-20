//! No-remainder quadratic and uniform compact quartic suppliers.

use super::cubic::{endpoint_is_zero, endpoint_rational_sign, substitute};
use super::*;

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
fn quartic_amplitude_forms(n: &B, d: &[P]) -> (R, R, R, R) {
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

pub(super) fn certify(
    n: &B,
    d: &[P],
    request: &EffectiveBkwRequest,
    degree: &DegreeBudget,
) -> Result<EffectiveBkwCertificate, RecurrenceError> {
    let order = d.len() - 1;
    let characteristic = B::new(d.iter().rev().map(|p| p.coeffs().to_vec()).collect());
    let supplied = request.pair_factorization.as_ref();
    let disc = if let Some(pair) = supplied {
        let a = decode(&pair.real_part)?;
        let b = decode(&pair.norm_squared)?;
        (a.pow(2) - b).scale(&q(4))
    } else {
        discriminant(&characteristic)
    };
    let domain = Domain::select(
        &disc,
        if supplied.is_some() { 2 } else { order },
        &request.interval,
    )
    .map_err(RecurrenceError::from)?;
    if domain.lower.is_none() || domain.upper.is_none() {
        return Err(failure(
            "uniform/no-remainder phase supplier currently needs two finite collision endpoints",
        ));
    }
    let (p, a, b, cofactor, which) = if let Some(pair) = supplied {
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
    } else {
        let p = pair_product_polynomial(d);
        let b = R::of(variable());
        let a = b
            .mul(&R::of(t(&d[1])).neg().mul(&b).add(&R::of(t(&d[3]))))
            .div(&b.square().sub(&R::of(t(&d[4]))).scale(&q(2)));
        let c0 = R::of(t(&d[4])).div(&b);
        let c1 = R::of(t(&d[1])).add(&a.scale(&q(2)));
        (p, a, b, vec![R::of(c(1)), c1, c0], Branch::UniquePositive)
    };
    let mut out = empty_certificate(d, &p, &domain, degree);
    if domain.upper_is_nonpositive() {
        out.root_region = EffectiveBkwRootRegion::NegativeRay;
    }
    out.branch_selection=if order==2||supplied.is_some(){"r=0 (rational pair parameters, no algebraic branch needed)"}else{"b is the unique positive real pair-product root: quartic discriminant < 0 and constant coefficient < 0 give a nonreal pair and two opposite-sign real cofactor roots"}.into();
    out.signs.push(check_univariate(
        "characteristic discriminant < 0",
        &disc,
        &domain,
        -1,
    )?);
    if order == 4 && supplied.is_none() {
        out.signs.push(check_univariate(
            "opposite-sign cofactor roots",
            &d[4],
            &domain,
            -1,
        )?);
    }
    let y2 = b.sub(&a.square());
    out.signs.push(check_sign(
        "nonreal pair b-a^2",
        &p,
        &y2.signed_polynomial(&p),
        &domain,
        which,
        1,
    )?);
    // Check the factorization in the algebraic quotient, coefficient by coefficient.
    let pair = [R::of(c(1)), a.scale(&q(-2)), b.clone()];
    let mut product = vec![R::of(c(0)); order + 1];
    for (i, u) in pair.iter().enumerate() {
        for (j, v) in cofactor.iter().enumerate() {
            product[i + j] = product[i + j].add(&u.mul(v));
        }
    }
    for (j, value) in product.iter().enumerate() {
        if !rem(&value.sub(&R::of(t(&d[j]))).n, &p).is_zero() {
            return Err(failure("exact pair/cofactor factorization identity failed"));
        }
    }
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
    let imag_poly = imaginary.signed_polynomial(&p);
    let imaginary_zero = imag_poly.is_zero();
    if !imaginary_zero {
        out.signs.push(check_sign(
            "Im(N/H)/Im(lambda) < 0",
            &p,
            &imag_poly,
            &domain,
            which,
            -1,
        )?);
    }
    let mut endpoint_signs = vec![];
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
            let c0 = R::of(t(&d[4])).div(&x.square());
            let c1 = R::of(t(&d[1])).add(&x.scale(&q(2)));
            R::of(c(1)).add(&c1.mul(&inv)).add(&c0.mul(&inv.square()))
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
    let winding = if imaginary_zero {
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
    if degree.count.offsets != [0] || !(winding == 0 || (winding == -1 && endpoint_signs[1] == 1)) {
        return Err(failure(
            "phase winding and exterior sign do not close the degree budget",
        ));
    }
    out.phase.push(format!("continuous lift of arg F has winding {winding}; (n+1)*theta+arg F supplies n+{winding} sign changes by first half-integer-pi crossings"));
    if winding == -1 {
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
    if order == 2 {
        out.domination = "no_remainder: proper quadratic numerator, R_n=0 for n>=0".into();
        out.lean_declarations
            .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_rest_eq_zero".into());
    } else {
        let (ratio, amplitude, weights) = if let Some(constants) = &request.uniform_constants {
            let (ratio, amplitude, weights) = parse_uniform_constants(constants, n.rows().len())?;
            if supplied.is_some() {
                let bounds = constants.cofactor_weights.as_ref().ok_or_else(|| {
                    failure("supplied compact factorization requires cofactor weights")
                })?;
                if bounds.len() != order - 2 {
                    return Err(RecurrenceError::InvalidInput(
                        "cofactor weight count mismatch".into(),
                    ));
                }
                let mut budget = Q::zero();
                for (j, (weight, h)) in bounds.iter().zip(cofactor.iter().skip(1)).enumerate() {
                    let j = j + 1;
                    let gamma = rational(weight)?;
                    if gamma.is_negative() {
                        return Err(RecurrenceError::InvalidInput(
                            "negative cofactor weight".into(),
                        ));
                    }
                    budget += &gamma / ratio.pow(j as i32);
                    let f = b.pow(j).scale(&(&gamma * &gamma)).sub(&h.square());
                    out.signs.push(check_sign(
                        &format!("Cauchy cofactor weight {j}"),
                        &p,
                        &f.signed_polynomial(&p),
                        &domain,
                        which,
                        0,
                    )?);
                }
                if budget >= Q::one() {
                    return Err(failure("cofactor Cauchy budget sum gamma_j/q^j is not < 1"));
                }
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
                        &p,
                        &f.signed_polynomial(&p),
                        &domain,
                        which,
                        0,
                    )?);
                }
            }
            out.signs.push(check_sign(
                "uniform amplitude lower bound",
                &p,
                &nn.sub(&hn.scale(&(&amplitude * &amplitude)))
                    .signed_polynomial(&p),
                &domain,
                which,
                0,
            )?);
            for (j, (weight, coefficient)) in weights.iter().zip(&nc).enumerate() {
                let f = b
                    .pow(j)
                    .scale(&(weight * weight))
                    .sub(&coefficient.square());
                out.signs.push(check_sign(
                    &format!("numerator weight {j}"),
                    &p,
                    &f.signed_polynomial(&p),
                    &domain,
                    which,
                    0,
                )?);
            }
            (ratio, amplitude, weights)
        } else {
            if supplied.is_some() {
                return Err(failure(
                    "supplied compact factorization requires uniform constants",
                ));
            }
            let search = UniformSearch {
                polynomial: &p,
                domain: &domain,
                branch: which,
                sample_branch: branch_box(&p, &domain.sample, which, &(q(1) / q(1 << 32)))?,
            };
            let (ratio, signs) =
                search.choose("modulus", (1..16).map(|k| q(k) / q(16)), |ratio| {
                    let ratio_squared = ratio * ratio;
                    let gap = b.scale(&ratio_squared).add(&cofactor[2]);
                    vec![
                        (
                            "uniform modulus gap, part one".into(),
                            gap.signed_polynomial(&p),
                        ),
                        (
                            "uniform modulus gap, part two".into(),
                            gap.square()
                                .sub(&b.mul(&cofactor[1].square()).scale(&ratio_squared))
                                .signed_polynomial(&p),
                        ),
                    ]
                })?;
            out.signs.extend(signs);
            let (amplitude, signs) = search.choose(
                "amplitude",
                (0..=16).map(|k| q(1) / q(1 << k)),
                |amplitude| {
                    vec![(
                        "uniform amplitude lower bound".into(),
                        nn.sub(&hn.scale(&(amplitude * amplitude)))
                            .signed_polynomial(&p),
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
                                .signed_polynomial(&p),
                        )]
                    },
                )?;
                weights.push(weight);
                out.signs.extend(signs);
            }
            (ratio, amplitude, weights)
        };
        let (cutoff, tail) = cutoff(order, &ratio, &amplitude, &weights)?;
        out.cutoff = cutoff;
        out.tail_at_cutoff = tail.to_string();
        let weight_strings: Vec<_> = weights.iter().map(ToString::to_string).collect();
        out.domination=format!("uniform_compact_tail: q={ratio}, m={amplitude}, beta={weight_strings:?}; T_n=sum_j beta_j sum_(k>=n+2-j) binomial(k+d-3,d-3) q^k; T_N<m and positive-series tails decrease");
        out.lean_declarations
            .push("RealRooted.EffectiveBKW.PairExpansion.dom_of_tail".into());
    }
    Ok(out)
}

fn collision_root(p: &B, domain: &Domain, left: bool) -> Result<(R, String), RecurrenceError> {
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
