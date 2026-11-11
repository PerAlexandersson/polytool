//! Bounded exact right factorization over Q(n)[x] for second-order steps.
//!
//! We enumerate at most 64 monic divisors B of the leading coefficient,
//! up to index scalars, and seek A + B D with deg_x A <= 2. The Riccati
//! remainder is quadratic, not generally linear in A. Descending coefficient
//! equations independent of the still unknown lower coefficients are solved
//! over Q(n) by exact bivariate factorization (budget 256). Unconstrained
//! coefficients are set to zero; nontriangular systems remain unsearched.
//! There are at most 64 coefficient-search nodes per divisor and degree
//! bound and 128 returned products. Input x/index degrees are at most 8/4
//! (x-degree at most 4 after removing the common multiplier), rational
//! search-coefficient numerator and denominator degrees at most 24 and
//! coefficient sizes at most 512 bits. Exceeding any bound is incomplete
//! recognition, never a disproof. Every candidate is multiplied out in the
//! Ore algebra before use. Only order two is searched.

use super::*;

/// One first-order factor `a(n,x) + b(n,x) D` of an iterated-preserver step.
/// Coefficients are rational functions of the unchanged step index `n`.
/// Factors are stored in application order, so factor zero acts first.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct IteratedPreserverFactor {
    /// Ascending x coefficients of `a(n,x)`.
    pub multiplication: Vec<IndexRationalFunction>,
    /// Ascending x coefficients of `b(n,x)`.
    pub derivative: Vec<IndexRationalFunction>,
}

/// A polynomial in x over Q(n), represented by a bivariate numerator and
/// one index denominator. This small adapter keeps all Ore coefficients
/// polynomial after clearing index denominators.
#[derive(Clone, Debug)]
struct K {
    p: Biv,
    q: Vec<Q>,
}

fn exact_index_quotient(a: &[Q], b: &[Q]) -> Vec<Q> {
    let (quotient, remainder) = Polynomial::new(a.to_vec()).div_rem(&Polynomial::new(b.to_vec()));
    debug_assert!(remainder.is_zero());
    vec_poly::trim(quotient.coeffs())
}

impl K {
    fn new(p: Biv, q: Vec<Q>) -> Self {
        debug_assert!(!vec_poly::is_zero(&q));
        let mut gcd = Polynomial::new(q.clone());
        for j in 0..=p.variable_degree().unwrap_or(0) {
            gcd = gcd.gcd(&Polynomial::new(p.variable_coefficient(j)));
        }
        let g = vec_poly::trim(gcd.coeffs());
        let q = exact_index_quotient(&q, &g);
        let scale = q.last().unwrap().recip();
        let mut numerator = Biv::zero();
        for j in 0..=p.variable_degree().unwrap_or(0) {
            let column = vec_poly::scale(
                &exact_index_quotient(&p.variable_coefficient(j), &g),
                &scale,
            );
            numerator =
                numerator.add(&Biv::from_index_polynomial(&column).multiply_by_variable_power(j));
        }
        Self {
            p: numerator,
            q: vec_poly::scale(&q, &scale),
        }
    }

    fn polynomial(p: Biv) -> Self {
        Self {
            p,
            q: vec![Q::one()],
        }
    }
    fn zero() -> Self {
        Self::polynomial(Biv::zero())
    }
    fn one() -> Self {
        Self::polynomial(Biv::one())
    }
    fn add(&self, rhs: &Self) -> Self {
        Self::new(
            self.p
                .multiply(&Biv::from_index_polynomial(&rhs.q))
                .add(&rhs.p.multiply(&Biv::from_index_polynomial(&self.q))),
            vec_poly::mul(&self.q, &rhs.q),
        )
    }
    fn neg(&self) -> Self {
        Self {
            p: self.p.scale(&-Q::one()),
            q: self.q.clone(),
        }
    }
    fn sub(&self, rhs: &Self) -> Self {
        self.add(&rhs.neg())
    }
    fn mul(&self, rhs: &Self) -> Self {
        Self::new(self.p.multiply(&rhs.p), vec_poly::mul(&self.q, &rhs.q))
    }
    fn derivative(&self) -> Self {
        Self::new(self.p.variable_derivative(), self.q.clone())
    }
    fn degree(&self) -> Option<usize> {
        self.p.variable_degree()
    }
    fn shift(&self, j: usize) -> Self {
        Self {
            p: self.p.multiply_by_variable_power(j),
            q: self.q.clone(),
        }
    }
    fn column(&self, j: usize) -> Self {
        Self::new(
            Biv::from_index_polynomial(&self.p.variable_coefficient(j)),
            self.q.clone(),
        )
    }
    fn in_bounds(&self) -> bool {
        self.p.index_degree().unwrap_or(0) <= 24
            && self.q.len() <= 25
            && self
                .p
                .rows()
                .iter()
                .flatten()
                .chain(&self.q)
                .all(|c| c.numer().bits() <= 512 && c.denom().bits() <= 512)
    }
    fn div_rem(&self, rhs: &Self) -> Option<(Self, Self)> {
        let d = rhs.degree()?;
        let leading = rhs.p.variable_coefficient(d);
        let mut remainder = self.clone();
        let mut quotient = Self::zero();
        while let Some(e) = remainder.degree().filter(|&e| e >= d) {
            if !remainder.in_bounds() || !quotient.in_bounds() {
                return None;
            }
            let ratio = Self::new(
                Biv::from_index_polynomial(&vec_poly::mul(
                    &remainder.p.variable_coefficient(e),
                    &rhs.q,
                )),
                vec_poly::mul(&remainder.q, &leading),
            )
            .shift(e - d);
            quotient = quotient.add(&ratio);
            remainder = remainder.sub(&ratio.mul(rhs));
            if remainder.degree().is_some_and(|new| new >= e) {
                return None;
            }
        }
        Some((quotient, remainder))
    }
    fn divide(&self, rhs: &Self) -> Option<Self> {
        let (q, r) = self.div_rem(rhs)?;
        r.p.is_zero().then_some(q)
    }
    fn monic(&self) -> Option<Self> {
        Some(Self::new(
            self.p.clone(),
            self.p.variable_coefficient(self.degree()?),
        ))
    }
    fn gcd(&self, other: &Self) -> Option<Self> {
        let (mut a, mut b) = (self.clone(), other.clone());
        while !b.p.is_zero() {
            let remainder = a.div_rem(&b)?.1;
            (a, b) = (b, remainder);
        }
        a.monic().filter(Self::in_bounds)
    }
    fn functions(&self) -> Vec<IndexRationalFunction> {
        (0..=self.degree().unwrap_or(0))
            .map(|j| {
                IndexRationalFunction::from_polynomials(&self.p.variable_coefficient(j), &self.q)
            })
            .collect()
    }
    fn parse(coefficients: &[IndexRationalFunction]) -> Option<Self> {
        if coefficients.len() > 7 {
            return None;
        }
        let mut result = Self::zero();
        for (j, c) in coefficients.iter().enumerate() {
            let (p, q) = c.parse().ok()?;
            result = result.add(&Self::new(Biv::from_index_polynomial(&p), q).shift(j));
        }
        Some(result)
    }
}

impl IteratedPreserverFactor {
    fn of(a: &K, b: &K) -> Self {
        Self {
            multiplication: a.functions(),
            derivative: b.functions(),
        }
    }
    fn parts(&self) -> Option<(DOp, Vec<Q>)> {
        let a = K::parse(&self.multiplication)?;
        let b = K::parse(&self.derivative)?;
        let denominator = vec_poly::mul(&a.q, &b.q);
        let mut op = DOp::single(0, 0, a.p.multiply(&Biv::from_index_polynomial(&b.q)));
        op.add_term((0, 1), b.p.multiply(&Biv::from_index_polynomial(&a.q)));
        Some((op, denominator))
    }
}

pub(super) fn product_parts(factors: &[IteratedPreserverFactor]) -> Option<(DOp, Vec<Q>)> {
    if factors.len() != 2 {
        return None;
    }
    let mut product = DOp::identity();
    let mut denominator = vec![Q::one()];
    for factor in factors {
        let (op, q) = factor.parts()?;
        product = op.compose(&product);
        denominator = vec_poly::mul(&denominator, &q);
    }
    Some((product, denominator))
}

/// All rational-function roots of a nonzero polynomial in z of degree <=2.
/// Factoring over Q[n,z] yields precisely the linear-in-z factors over Q(n).
fn roots(coefficients: &[K; 3]) -> Vec<K> {
    let mut p = K::zero();
    for (j, c) in coefficients.iter().enumerate() {
        p = p.add(&c.shift(j));
    }
    if p.p.is_zero() || !p.in_bounds() {
        return Vec::new();
    }
    let Ok(factorization) = p.p.to_bivar().irreducible_factorization_with_budget(256) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for (factor, _) in factorization.factors {
        let factor = Biv::from_bivar(&factor);
        if factor.variable_degree() == Some(1) {
            let root = K::new(
                Biv::from_index_polynomial(&vec_poly::neg(&factor.variable_coefficient(0))),
                factor.variable_coefficient(1),
            );
            if root.in_bounds() {
                out.push(root);
            }
        }
    }
    out
}

/// The cleared Riccati remainder: B C0 - A C1 + E A B' + E A^2 - B E A'.
/// Here C2 = B E, and the quotient's constant coefficient is
/// (C1 - E(B' + A))/B. Vanishing is exactly right divisibility.
fn remainder(a: &K, b: &K, e: &K, c: &[K; 3]) -> K {
    b.mul(&c[0])
        .sub(&a.mul(&c[1]))
        .add(&e.mul(a).mul(&b.derivative()))
        .add(&e.mul(a).mul(a))
        .sub(&c[2].mul(&a.derivative()))
}

struct Search<'a> {
    b: &'a K,
    e: &'a K,
    c: &'a [K; 3],
    nodes: usize,
    out: Vec<K>,
}
impl Search<'_> {
    fn coefficients(&mut self, a: K, j: usize) {
        if self.nodes >= 64 || !a.in_bounds() {
            return;
        }
        self.nodes += 1;
        let xj = K::one().shift(j);
        let r0 = remainder(&a, self.b, self.e, self.c);
        let rplus = remainder(&a.add(&xj), self.b, self.e, self.c);
        let rminus = remainder(&a.sub(&xj), self.b, self.e, self.c);
        let half = K::polynomial(Biv::constant(Q::new(BigInt::one(), BigInt::from(2))));
        let r1 = rplus.sub(&rminus).mul(&half);
        let r2 = rplus.add(&rminus).mul(&half).sub(&r0);
        // Unknown lower terms h have degree <=j-1. Bound the degree of
        // h(-C1+E B'+2E(A+z*x^j))+E h^2-C2 h'. Above it the equation is independent
        // of h, so descending substitution loses no solutions unless an
        // unconstrained coefficient is fixed to zero.
        let lower_bound =
            if j == 0 {
                -1
            } else {
                let linear = self.c[1].neg().add(&self.e.mul(&self.b.derivative())).add(
                    &self
                        .e
                        .mul(&a)
                        .mul(&K::polynomial(Biv::constant(Q::from_integer(
                            BigInt::from(2),
                        )))),
                );
                [
                    linear.degree().map(|d| d as i64 + j as i64 - 1),
                    self.e.degree().map(|d| d as i64 + 2 * j as i64 - 1),
                    self.c[2].degree().map(|d| d as i64 + j as i64 - 2),
                ]
                .into_iter()
                .flatten()
                .max()
                .unwrap_or(-1)
            };
        let top = [r0.degree(), r1.degree(), r2.degree()]
            .into_iter()
            .flatten()
            .max()
            .unwrap_or(0);
        let equation = (0..=top)
            .rev()
            .filter(|&k| k as i64 > lower_bound)
            .find_map(|k| {
                let coefficients = [r0.column(k), r1.column(k), r2.column(k)];
                coefficients
                    .iter()
                    .any(|c| !c.p.is_zero())
                    .then_some(coefficients)
            });
        let values = equation.as_ref().map_or_else(|| vec![K::zero()], roots);
        for value in values {
            let next = a.add(&value.shift(j));
            if j > 0 {
                self.coefficients(next, j - 1);
            } else if remainder(&next, self.b, self.e, self.c).p.is_zero() {
                self.out.push(next);
            }
        }
    }
}

pub(super) fn decompositions(terms: &StepCoefficients, f: &[Q]) -> Vec<Decomposition> {
    if terms.keys().max() != Some(&2)
        || terms.values().any(|c| {
            c.len() > 5
                || x_degree(c).unwrap_or(0) > 8
                || c.iter()
                    .flatten()
                    .any(|v| v.numer().bits() > 512 || v.denom().bits() > 512)
        })
    {
        return Vec::new();
    }
    let g = fixed_content(terms);
    let mut c: [K; 3] = std::array::from_fn(|d| {
        K::polynomial(Biv::new(divide_rows(
            terms.get(&d).unwrap_or(&Vec::new()),
            &g,
        )))
    });
    let Some(common) = c[2].gcd(&c[1]).and_then(|g| g.gcd(&c[0])) else {
        return Vec::new();
    };
    for part in &mut c {
        let Some(quotient) = part.divide(&common) else {
            return Vec::new();
        };
        if quotient.degree().unwrap_or(0) > 4 || !quotient.in_bounds() {
            return Vec::new();
        }
        *part = quotient;
    }
    let multiplier = K::new(common.p, vec_poly::mul(&common.q, f));
    if !multiplier.in_bounds() {
        return Vec::new();
    }
    let index_factor: Vec<_> = (0..=multiplier.degree().unwrap_or(0))
        .map(|j| (multiplier.p.variable_coefficient(j), multiplier.q.clone()))
        .collect();
    let Ok(leading) = c[2].p.to_bivar().irreducible_factorization_with_budget(256) else {
        return Vec::new();
    };
    let mut divisors = vec![K::one()];
    for (factor, multiplicity) in leading.factors {
        let p = Biv::from_bivar(&factor);
        let Some(degree) = p.variable_degree().filter(|&d| d > 0) else {
            continue;
        };
        let monic = K::new(p.clone(), p.variable_coefficient(degree));
        let old = divisors.clone();
        let mut power = K::one();
        for _ in 0..multiplicity {
            power = power.mul(&monic);
            for d in &old {
                if divisors.len() >= 64 {
                    return Vec::new();
                }
                divisors.push(d.mul(&power));
            }
        }
    }
    let mut out = Vec::new();
    for b in divisors {
        if b.degree().unwrap_or(0) > 2 {
            continue;
        }
        let Some(e) = c[2].divide(&b).filter(|e| e.degree().unwrap_or(0) <= 2) else {
            continue;
        };
        for degree in 0..=2 {
            let mut search = Search {
                b: &b,
                e: &e,
                c: &c,
                nodes: 0,
                out: Vec::new(),
            };
            search.coefficients(K::zero(), degree);
            for a in search.out {
                let Some(left) = c[1].sub(&e.mul(&b.derivative().add(&a))).divide(&b) else {
                    continue;
                };
                let candidate = Decomposition {
                    fixed_factor: g.clone(),
                    index_factor: index_factor.clone(),
                    preserver: None,
                    factors: vec![
                        IteratedPreserverFactor::of(&a, &b),
                        IteratedPreserverFactor::of(&left, &e),
                    ],
                };
                if candidate.identity_holds(terms, f) && !out.contains(&candidate) {
                    out.push(candidate);
                }
                if out.len() >= 128 {
                    return out;
                }
            }
        }
    }
    out
}

pub(super) fn factor_hypotheses(
    factors: &[IteratedPreserverFactor],
    region: ClosedRegion,
    first: usize,
    previous_degree: Option<&[Q]>,
) -> (Vec<LiftHypothesis>, String) {
    let Some(mut degree) = previous_degree.map(<[Q]>::to_vec) else {
        return (
            vec![unresolved_hypothesis(
                "factor_degrees",
                "the input degree is known",
                "no exact nonnegative degree law".into(),
            )],
            String::new(),
        );
    };
    let mut out = Vec::new();
    for (i, factor) in factors.iter().enumerate() {
        let Some((op, q)) = factor.parts() else {
            out.push(unresolved_hypothesis(
                "factor",
                "a valid first-order factor",
                "invalid coefficients".into(),
            ));
            break;
        };
        let terms: StepCoefficients =
            op.0.iter()
                .map(|((_, d), c)| (*d, c.rows().to_vec()))
                .collect();
        let delta = degree_step(&terms).unwrap_or(-1);
        let lambda = leading_lambda(&terms, delta, &[Q::zero(), Q::one()], &degree);
        let mut conditions: Vec<_> = denominator_hypothesis(&q, first).into_iter().collect();
        conditions.push(sign_hypothesis(
            "image_degree",
            "the intermediate output has exactly the recorded degree",
            "leading coefficient",
            &IndexFunction::new(lambda.clone(), q.clone()),
            first,
            IndexSignCondition::Nonzero,
        ));
        let next_degree = vec_poly::add(&degree, &[Q::from_integer(BigInt::from(delta))]);
        conditions.push(sign_hypothesis(
            "image_nonnegative_degree",
            "the intermediate degree is nonnegative",
            "degree",
            &IndexFunction::polynomial(next_degree.clone()),
            first,
            IndexSignCondition::Nonnegative,
        ));
        let menu = super::decompositions(&terms, &q)
            .into_iter()
            .find_map(|candidate| {
                let mut conditions = multiplier_hypotheses(&candidate, region, first);
                if let Some(preserver) = &candidate.preserver {
                    conditions
                        .extend(preserver_hypotheses(preserver, region, first, Some(&degree)).0);
                }
                (status_of(&conditions) == RealRootednessTheoremStatus::Proved)
                    .then_some(conditions)
            });
        if let Some(menu) = menu {
            conditions.extend(menu);
        } else {
            conditions.extend(liu_wang_conditions(&terms, &q, &degree, region, first));
        }
        for condition in &mut conditions {
            condition.name = format!("factor_{i}_{}", condition.name);
        }
        out.extend(conditions);
        degree = next_degree;
    }
    (out,"; each first-order factor is certified by the existing preserver menu or the Liu--Wang derivative step (Theorem 2.1, g=f'); on the closed ray the latter follows by continuity from strictly negative-root inputs of the same degree, with positive leading degree coefficient keeping the limit nonzero; reflection gives the nonnegative ray. No simplicity or interlacing of the composite is asserted".into())
}

fn liu_wang_conditions(
    terms: &StepCoefficients,
    q: &[Q],
    degree: &[Q],
    region: ClosedRegion,
    first: usize,
) -> Vec<LiftHypothesis> {
    if region == ClosedRegion::RealLine
        || terms.get(&0).and_then(x_degree).unwrap_or(0) > 1
        || terms.get(&1).and_then(x_degree).unwrap_or(0) > 2
        || degree_step(terms) != Some(1)
    {
        return vec![unresolved_hypothesis(
            "preserver",
            "the factor preserves the common region",
            "neither an existing menu theorem nor a degree-one Liu--Wang step applies".into(),
        )];
    }
    let mut parts: [Bivariate; 2] =
        std::array::from_fn(|d| terms.get(&d).cloned().unwrap_or_default());
    if region == ClosedRegion::NonnegativeRay {
        for (d, part) in parts.iter_mut().enumerate() {
            for row in part {
                for (j, c) in row.iter_mut().enumerate() {
                    if (j + d) % 2 == 1 {
                        *c = -c.clone();
                    }
                }
            }
        }
    }
    // Multiplication of a factor by -1 has no effect on root locations.
    // Choose the orientation with positive leading degree coefficient.
    let lambda = vec_poly::add(
        &x_column(&parts[0], 1),
        &vec_poly::mul(degree, &x_column(&parts[1], 2)),
    );
    if IndexFunction::new(lambda, q.to_vec())
        .decide(first, IndexSignCondition::Negative)
        .holds()
    {
        for part in &mut parts {
            for row in part {
                for c in row {
                    *c = -c.clone();
                }
            }
        }
    }
    super::super::liu_wang::step_hypotheses(&parts[0], &parts[1], q, degree, first)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_order_menu_factors_share_a_region_and_replay_in_application_order() {
        // (theta+n) o (1+D): Laguerre multiplier after a phi(D) step.
        let derivative = IteratedPreserverFactor::of(&K::one(), &K::one());
        let multiplier = IteratedPreserverFactor::of(
            &K::polynomial(Biv::from_index_polynomial(&[Q::zero(), Q::one()])),
            &K::one().shift(1),
        );
        let factors = vec![derivative, multiplier];
        for region in [ClosedRegion::NonpositiveRay, ClosedRegion::RealLine] {
            let (hypotheses, _) =
                factor_hypotheses(&factors, region, 1, Some(&[Q::zero(), Q::one()]));
            assert_eq!(
                status_of(&hypotheses),
                RealRootednessTheoremStatus::Proved,
                "{hypotheses:#?}"
            );
        }
        let (product, denominator) = product_parts(&factors).unwrap();
        let expected = factors[1]
            .parts()
            .unwrap()
            .0
            .compose(&factors[0].parts().unwrap().0);
        assert_eq!(product, expected);
        assert_eq!(denominator, vec![Q::one()]);
    }

    #[test]
    fn a_rolle_factor_cannot_silently_kill_an_intermediate_constant() {
        let derivative = IteratedPreserverFactor::of(&K::zero(), &K::one());
        let identity = IteratedPreserverFactor::of(&K::one(), &K::zero());
        let (hypotheses, _) = factor_hypotheses(
            &[derivative, identity],
            ClosedRegion::RealLine,
            1,
            Some(&[Q::zero()]),
        );
        assert_ne!(status_of(&hypotheses), RealRootednessTheoremStatus::Proved);
    }
}
