//! Branch signs from exact resultants, Sturm isolation and rational boxes.

use super::*;

/// A box excluding an off-branch resultant zero.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwBox {
    /// Parameter enclosure of the critical root.
    pub parameter: [String; 2],
    /// Enclosure of the selected branch across this parameter box.
    pub branch: [String; 2],
    /// Signs of the branch equation on its lower and upper walls.
    pub wall_signs: [i8; 2],
    /// Sign of the condition on the entire box.
    pub condition_sign: i8,
}

/// A replayable full-interval sign test on one continuous algebraic branch.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwSignCertificate {
    /// Mathematical obligation.
    pub name: String,
    /// Condition polynomial in the branch variable, each row ascending in t.
    pub polynomial: Vec<Vec<String>>,
    /// Exact resultant against the branch equation, ascending in t.
    pub resultant: Vec<String>,
    /// Rational parameter inside the domain.
    pub sample_parameter: String,
    /// Exact isolating enclosure of the branch at the sample.
    pub sample_branch: [String; 2],
    /// The constant strict sign, zero for an identity, or 1 for a
    /// nonnegative condition proved with `nonnegative_cell_samples`.
    pub sign: i8,
    /// Boxes covering every resultant zero strictly inside the domain.
    pub off_branch_boxes: Vec<EffectiveBkwBox>,
    /// For a weak inequality with branch zeros, a positive sample in each
    /// open resultant cell. Continuity extends nonnegativity over the zeros.
    pub nonnegative_cell_samples: Vec<EffectiveBkwCellSample>,
}

/// An exact selected-branch sign in an open resultant cell.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct EffectiveBkwCellSample {
    /// Rational parameter in the cell.
    pub parameter: String,
    /// Isolating branch-root enclosure at that parameter.
    pub branch: [String; 2],
    /// Strictly positive condition sign on this whole root-free cell.
    pub sign: i8,
}

#[derive(Clone)]
pub(super) struct I {
    pub lo: Q,
    pub hi: Q,
}
impl I {
    pub fn point(x: Q) -> Self {
        Self {
            lo: x.clone(),
            hi: x,
        }
    }
    pub fn new(lo: Q, hi: Q) -> Self {
        Self { lo, hi }
    }
    pub fn add(&self, r: &Self) -> Self {
        Self::new(&self.lo + &r.lo, &self.hi + &r.hi)
    }
    pub fn mul(&self, r: &Self) -> Self {
        let v = [
            &self.lo * &r.lo,
            &self.lo * &r.hi,
            &self.hi * &r.lo,
            &self.hi * &r.hi,
        ];
        Self::new(
            v.iter().min().unwrap().clone(),
            v.iter().max().unwrap().clone(),
        )
    }
    pub fn sign(&self) -> i8 {
        if self.lo.is_positive() {
            1
        } else if self.hi.is_negative() {
            -1
        } else {
            0
        }
    }
    pub fn wire(&self) -> [String; 2] {
        [self.lo.to_string(), self.hi.to_string()]
    }
}
pub(super) fn eval_p(p: &P, x: &I) -> I {
    p.coeffs().iter().rev().fold(I::point(Q::zero()), |v, c| {
        v.mul(x).add(&I::point(c.clone()))
    })
}
pub(super) fn eval_b(p: &B, t: &I, r: &I) -> I {
    p.rows().iter().rev().fold(I::point(Q::zero()), |v, c| {
        v.mul(r).add(&eval_p(&P::new(c.clone()), t))
    })
}
fn sturm(p: &P) -> crate::sturm::SturmChain {
    // Elimination produces high-degree repeated factors. Remove them with
    // the shared integer PRS before building the small rational Sturm chain.
    let squarefree = crate::root_count::squarefree_part_primitive_bigint_coeffs(
        &vec_poly::clear_denominators(p.coeffs()),
    );
    crate::sturm::SturmChain::from_bigint_coeffs(&squarefree)
}
fn squarefree(p: &P) -> P {
    P::new(
        crate::root_count::squarefree_part_primitive_bigint_coeffs(&vec_poly::clear_denominators(
            p.coeffs(),
        ))
        .into_iter()
        .map(Q::from_integer)
        .collect(),
    )
}
fn refine_integer_root(integers: &[num_bigint::BigInt], bounds: &I, epsilon: &Q) -> I {
    if &bounds.hi - &bounds.lo <= *epsilon {
        return bounds.clone();
    }
    let sign = |x: &Q| {
        let v = vec_poly::evaluate_homogeneous(integers, x.numer(), x.denom(), integers.len() - 1)
            .expect("directed degree equals polynomial degree");
        if v.is_positive() {
            1
        } else if v.is_negative() {
            -1
        } else {
            0
        }
    };
    let (mut lo, mut hi) = (bounds.lo.clone(), bounds.hi.clone());
    let lower_sign = sign(&lo);
    while &hi - &lo > *epsilon {
        let mid = (&lo + &hi) / q(2);
        let sm = sign(&mid);
        if sm == 0 {
            lo = mid.clone();
            hi = mid;
        } else if sm == lower_sign {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    I::new(lo, hi)
}
// Bisect the original isolating interval. Overlap with a newly isolated root
// is insufficient to identify it: a neighbouring root's enclosure can cross
// the original interval boundary without its root belonging to that interval.
fn refine_root(p: &P, bounds: &I, epsilon: &Q) -> I {
    if &bounds.hi - &bounds.lo <= *epsilon {
        return bounds.clone();
    }
    let integers = crate::root_count::squarefree_part_primitive_bigint_coeffs(
        &vec_poly::clear_denominators(p.coeffs()),
    );
    refine_integer_root(&integers, bounds, epsilon)
}
fn roots(p: &P, epsilon: &Q) -> Vec<I> {
    let integers = crate::root_count::squarefree_part_primitive_bigint_coeffs(
        &vec_poly::clear_denominators(p.coeffs()),
    );
    let isolated = crate::root_count::isolate_real_roots_bigint_coeffs(&integers);
    let rational =
        |x: &crate::root_count::DyadicRational| Q::new(x.numerator().clone(), x.denominator());
    isolated
        .into_iter()
        .map(|root| {
            refine_integer_root(
                &integers,
                &I::new(rational(root.lower()), rational(root.upper())),
                epsilon,
            )
        })
        .collect()
}
pub(super) fn rational(s: &str) -> Result<Q, RecurrenceError> {
    parse_rational_coeff(s).map_err(Into::into)
}

#[derive(Clone)]
pub(super) struct End {
    pub polynomial: P,
    pub bounds: I,
}
impl End {
    fn wire(&self) -> EffectiveBkwEndpoint {
        EffectiveBkwEndpoint::Algebraic {
            polynomial: encode(&self.polynomial),
            enclosure: self.bounds.wire(),
        }
    }
    pub(super) fn refine(&mut self, epsilon: &Q) {
        self.bounds = refine_root(&self.polynomial, &self.bounds, epsilon);
    }
}
#[derive(Clone)]
pub(super) struct Domain {
    pub wire: EffectiveBkwInterval,
    pub lower: Option<End>,
    pub upper: Option<End>,
    pub sample: Q,
}
impl Domain {
    pub fn upper_is_nonpositive(&self) -> bool {
        let Some(mut end) = self.upper.clone() else {
            return false;
        };
        for bits in 8..=128 {
            if end.bounds.hi <= Q::zero() {
                return true;
            }
            if end.bounds.lo >= Q::zero() {
                return false;
            }
            if end.polynomial.evaluate(&Q::zero()).is_zero() {
                return true;
            }
            end.refine(&Q::new(1.into(), num_bigint::BigInt::one() << bits));
        }
        false
    }
    pub fn select(
        disc: &P,
        order: usize,
        request: &Option<EffectiveBkwInterval>,
    ) -> Result<Self, ChartError> {
        let parse_end = |e: &EffectiveBkwEndpoint| -> Result<Option<End>, RecurrenceError> {
            match e {
                EffectiveBkwEndpoint::Rational { value } => {
                    let x = rational(value)?;
                    Ok(Some(End {
                        polynomial: P::new(vec![-&x, Q::one()]),
                        bounds: I::point(x),
                    }))
                }
                EffectiveBkwEndpoint::Algebraic {
                    polynomial,
                    enclosure,
                } => {
                    let p = decode(polynomial)?.squarefree_part();
                    let i = I::new(rational(&enclosure[0])?, rational(&enclosure[1])?);
                    if p.is_zero()
                        || i.lo > i.hi
                        || if i.lo == i.hi {
                            !p.evaluate(&i.lo).is_zero()
                        } else {
                            p.evaluate(&i.lo).is_zero()
                                || p.evaluate(&i.hi).is_zero()
                                || sturm(&p).count_roots_in(&i.lo, &i.hi) != 1
                        }
                    {
                        return Err(RecurrenceError::InvalidInput(
                            "endpoint is not an isolating enclosure".into(),
                        ));
                    }
                    Ok(Some(End {
                        polynomial: p,
                        bounds: i,
                    }))
                }
                _ => Ok(None),
            }
        };
        let (lower, upper) = if let Some(interval) = request {
            if matches!(interval.lower, EffectiveBkwEndpoint::PositiveInfinity)
                || matches!(interval.upper, EffectiveBkwEndpoint::NegativeInfinity)
            {
                return Err(
                    RecurrenceError::InvalidInput("reversed infinite endpoints".into()).into(),
                );
            }
            (parse_end(&interval.lower)?, parse_end(&interval.upper)?)
        } else {
            let p = disc.squarefree_part();
            let rs = roots(&p, &Q::new(1.into(), (1u64 << 32).into()));
            if order == 3 {
                let e = rs
                    .into_iter()
                    .find(|i| i.hi.is_negative())
                    .ok_or(ChartError::NoNegativeCubicCollision)?;
                (
                    None,
                    Some(End {
                        polynomial: p,
                        bounds: e,
                    }),
                )
            } else {
                let selected = rs
                    .windows(2)
                    .find(|w| {
                        w[0].hi < w[1].lo
                            && disc.evaluate(&((&w[0].hi + &w[1].lo) / q(2))).is_negative()
                    })
                    .ok_or_else(|| failure("no bounded nonreal-pair discriminant component"))?;
                (
                    Some(End {
                        polynomial: p.clone(),
                        bounds: selected[0].clone(),
                    }),
                    Some(End {
                        polynomial: p,
                        bounds: selected[1].clone(),
                    }),
                )
            }
        };
        let sample = match (&lower, &upper) {
            (Some(l), Some(u)) => {
                if l.bounds.hi >= u.bounds.lo {
                    return Err(failure("endpoint enclosures are not strictly ordered").into());
                }
                (&l.bounds.hi + &u.bounds.lo) / q(2)
            }
            (None, Some(u)) => u.bounds.lo.floor() - q(1),
            (Some(l), None) => l.bounds.hi.ceil() + q(1),
            _ => Q::zero(),
        };
        let wire = EffectiveBkwInterval {
            lower: lower
                .as_ref()
                .map(End::wire)
                .unwrap_or(EffectiveBkwEndpoint::NegativeInfinity),
            upper: upper
                .as_ref()
                .map(End::wire)
                .unwrap_or(EffectiveBkwEndpoint::PositiveInfinity),
        };
        Ok(Self {
            wire,
            lower,
            upper,
            sample,
        })
    }
    fn location(&self, p: &P, r: &mut I, epsilon: &Q) -> Result<i8, RecurrenceError> {
        // -1 outside, 0 exactly an endpoint, 1 strictly inside. Refine until
        // ordered, or certify equality through a common polynomial root.
        let mut ds = self.clone();
        for k in 0..128 {
            let eps = epsilon / q(2).pow(k);
            let mut touch = false;
            for (end, is_lower) in [(&mut ds.lower, true), (&mut ds.upper, false)] {
                if let Some(e) = end {
                    if (is_lower && r.hi < e.bounds.lo) || (!is_lower && r.lo > e.bounds.hi) {
                        return Ok(-1);
                    }
                    if r.lo <= e.bounds.hi && r.hi >= e.bounds.lo {
                        let g = p.gcd(&e.polynomial);
                        let lo = r.lo.clone().max(e.bounds.lo.clone());
                        let hi = r.hi.clone().min(e.bounds.hi.clone());
                        if (lo == hi && g.evaluate(&lo).is_zero())
                            || (lo < hi && sturm(&g).count_roots_in(&lo, &hi) > 0)
                        {
                            return Ok(0);
                        }
                        e.refine(&eps);
                        touch = true;
                    }
                }
            }
            if !touch {
                return Ok(1);
            }
            *r = refine_root(p, r, &eps);
        }
        Err(failure(
            "critical root/endpoint separation exceeded 128 refinements",
        ))
    }
    pub fn critical_roots(&self, p: &P) -> Result<Vec<I>, RecurrenceError> {
        let eps = Q::new(1.into(), (1u64 << 24).into());
        let mut out = vec![];
        for mut root in roots(p, &eps) {
            if self.location(p, &mut root, &eps)? == 1 {
                out.push(root);
            }
        }
        Ok(out)
    }
    pub fn endpoint_sign(&self, p: &P, left: bool) -> Result<i8, RecurrenceError> {
        let mut e = if left {
            self.lower.clone()
        } else {
            self.upper.clone()
        }
        .ok_or_else(|| failure("finite endpoint required"))?;
        for k in 8..128 {
            let sign = eval_p(p, &e.bounds).sign();
            if sign != 0 {
                return Ok(sign);
            }
            e.refine(&Q::new(1.into(), num_bigint::BigInt::one() << k));
        }
        Err(failure("could not separate endpoint sign from zero"))
    }
}

#[derive(Clone, Copy)]
pub(super) enum Branch {
    UniqueReal,
    UniquePositive,
    LargestPositive,
    Zero,
}
pub(super) fn branch_box(p: &B, t: &Q, which: Branch, epsilon: &Q) -> Result<I, RecurrenceError> {
    if matches!(which, Branch::Zero) {
        return Ok(I::point(Q::zero()));
    }
    let p = P::new(p.evaluate_variable(t));
    let selected: Vec<_> = roots(&p, epsilon)
        .into_iter()
        .filter(|i| {
            !matches!(which, Branch::UniquePositive | Branch::LargestPositive) || i.lo.is_positive()
        })
        .collect();
    if matches!(which, Branch::LargestPositive) {
        return selected
            .last()
            .cloned()
            .ok_or_else(|| failure("pair-product polynomial has no positive root"));
    }
    if selected.len() != 1 {
        return Err(failure(
            "algebraic branch is not uniquely isolated at the rational sample",
        ));
    }
    Ok(selected[0].clone())
}
pub(super) fn check_sign(
    name: &str,
    p: &B,
    f: &B,
    domain: &Domain,
    which: Branch,
    want: i8,
) -> Result<EffectiveBkwSignCertificate, RecurrenceError> {
    #[cfg(test)]
    if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
        eprintln!("effective BKW sign: {name}");
    }
    let f = rem(f, p);
    let mut cert = EffectiveBkwSignCertificate {
        name: name.into(),
        polynomial: encode_b(&f),
        resultant: vec![],
        sample_parameter: domain.sample.to_string(),
        sample_branch: ["0".into(), "0".into()],
        sign: 0,
        off_branch_boxes: vec![],
        nonnegative_cell_samples: vec![],
    };
    if f.is_zero() {
        if want == 0 {
            return Ok(cert);
        }
        return Err(failure(format!("{name}: condition vanishes identically")));
    }
    let mut sample = None;
    for bits in (8..=128).step_by(8) {
        let root = branch_box(
            p,
            &domain.sample,
            which,
            &Q::new(1.into(), num_bigint::BigInt::one() << bits),
        )?;
        let sign = eval_b(&f, &I::point(domain.sample.clone()), &root).sign();
        if sign != 0 {
            sample = Some((root, sign));
            break;
        }
    }
    // An exact sample only rejects; every acceptance still computes and
    // checks the full resultant and all its critical parameter cells.
    if let Some((_, sign)) = &sample {
        if (want != 0 && *sign != want) || (want == 0 && *sign < 0) {
            return Err(failure(format!(
                "{name}: exact branch sign at t={} is {sign}, expected {}",
                domain.sample,
                if want == 0 {
                    "nonnegative"
                } else if want > 0 {
                    "positive"
                } else {
                    "negative"
                }
            )));
        }
    }
    let k = resultant(p, &f);
    cert.resultant = encode(&k);
    #[cfg(test)]
    if std::env::var_os("POLYTOOL_EFFECTIVE_BKW_TRACE").is_some() {
        eprintln!("effective BKW resultant: {name}, degree {:?}", k.degree());
    }
    if k.is_zero() {
        return Err(failure(format!(
            "{name}: resultant vanishes identically; no branch-specific elimination available"
        )));
    }
    let Some((root, sign)) = sample else {
        if want == 0 {
            return certify_nonnegative_cells(cert, p, &f, &k, domain, which);
        }
        return Err(failure(format!("{name}: sample sign unresolved")));
    };
    cert.sample_branch = root.wire();
    cert.sign = sign;
    let critical_polynomial = squarefree(&k);
    for critical in domain.critical_roots(&critical_polynomial)? {
        let mut covered = None;
        for bits in (16..=128).step_by(8) {
            let eps = Q::new(1.into(), num_bigint::BigInt::one() << bits);
            let tc = refine_root(&critical_polynomial, &critical, &eps);
            let tm = (&tc.lo + &tc.hi) / q(2);
            let rb = branch_box(p, &tm, which, &eps)?;
            // A square-root-width branch box eventually covers parameter
            // motion for every finite derivative, without assuming |r'|<=16.
            let padding = Q::new(1.into(), num_bigint::BigInt::one() << (bits / 2));
            let rb = I::new(rb.lo - &padding, rb.hi + &padding);
            if matches!(which, Branch::UniquePositive | Branch::LargestPositive)
                && !rb.lo.is_positive()
            {
                continue;
            }
            let ws = [
                eval_b(p, &tc, &I::point(rb.lo.clone())).sign(),
                eval_b(p, &tc, &I::point(rb.hi.clone())).sign(),
            ];
            let fs = eval_b(&f, &tc, &rb).sign();
            if ws[0] * ws[1] == -1 && fs == sign {
                covered = Some(EffectiveBkwBox {
                    parameter: tc.wire(),
                    branch: rb.wire(),
                    wall_signs: ws,
                    condition_sign: fs,
                });
                break;
            }
        }
        if let Some(covered) = covered {
            cert.off_branch_boxes.push(covered);
        } else if want == 0 {
            return certify_nonnegative_cells(cert, p, &f, &critical_polynomial, domain, which);
        } else {
            return Err(failure(format!(
                "{name}: a resultant root may lie on the selected branch; box exclusion failed"
            )));
        }
    }
    Ok(cert)
}

fn certify_nonnegative_cells(
    mut cert: EffectiveBkwSignCertificate,
    p: &B,
    f: &B,
    k: &P,
    domain: &Domain,
    which: Branch,
) -> Result<EffectiveBkwSignCertificate, RecurrenceError> {
    let critical = domain.critical_roots(&squarefree(k))?;
    let mut samples = Vec::new();
    if critical.is_empty() {
        samples.push(domain.sample.clone());
    } else {
        let first = &critical[0];
        let last = critical.last().unwrap();
        let before = if let Some(mut end) = domain.lower.clone() {
            for bits in 8..=128 {
                if end.bounds.hi < first.lo {
                    break;
                }
                end.refine(&Q::new(1.into(), num_bigint::BigInt::one() << bits));
            }
            if end.bounds.hi >= first.lo {
                return Err(failure(
                    "weak sign: cannot separate first critical root from endpoint",
                ));
            }
            (end.bounds.hi + &first.lo) / q(2)
        } else {
            first.lo.floor() - q(1)
        };
        samples.push(before);
        for pair in critical.windows(2) {
            if pair[0].hi >= pair[1].lo {
                return Err(failure(
                    "weak sign: critical isolating intervals need further separation",
                ));
            }
            samples.push((&pair[0].hi + &pair[1].lo) / q(2));
        }
        let after = if let Some(mut end) = domain.upper.clone() {
            for bits in 8..=128 {
                if end.bounds.lo > last.hi {
                    break;
                }
                end.refine(&Q::new(1.into(), num_bigint::BigInt::one() << bits));
            }
            if end.bounds.lo <= last.hi {
                return Err(failure(
                    "weak sign: cannot separate last critical root from endpoint",
                ));
            }
            (&last.hi + end.bounds.lo) / q(2)
        } else {
            last.hi.ceil() + q(1)
        };
        samples.push(after);
    }
    cert.off_branch_boxes.clear();
    for sample in samples {
        let mut found = None;
        for bits in (8..=128).step_by(8) {
            let root = branch_box(
                p,
                &sample,
                which,
                &Q::new(1.into(), num_bigint::BigInt::one() << bits),
            )?;
            let sign = eval_b(f, &I::point(sample.clone()), &root).sign();
            if sign != 0 {
                found = Some((root, sign));
                break;
            }
        }
        let (root, sign) =
            found.ok_or_else(|| failure(format!("{}: weak cell sample unresolved", cert.name)))?;
        if sign < 0 {
            return Err(failure(format!(
                "{}: negative branch sign in resultant cell at t={sample}",
                cert.name
            )));
        }
        cert.nonnegative_cell_samples.push(EffectiveBkwCellSample {
            parameter: sample.to_string(),
            branch: root.wire(),
            sign,
        });
    }
    if let Some(sample) = cert.nonnegative_cell_samples.first() {
        cert.sample_parameter = sample.parameter.clone();
        cert.sample_branch = sample.branch.clone();
    }
    cert.sign = 1;
    Ok(cert)
}

pub(super) fn check_univariate(
    name: &str,
    f: &P,
    domain: &Domain,
    want: i8,
) -> Result<EffectiveBkwSignCertificate, RecurrenceError> {
    check_sign(name, &variable(), &t(f), domain, Branch::Zero, want)
}
