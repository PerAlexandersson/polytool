//! Conservative full-identity recognition, never certification from a prefix.
use crate::exponential_rows::BlockWeightSpec;
use crate::power_series::{Expr, SeriesExpression, SeriesExpressionLimits};
use crate::{numeric::factorial_bigint, vec_poly};
use num_rational::BigRational as Q;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

type ExponentialPolynomial = BTreeMap<Q, Vec<Q>>;
const DEGREE: usize = 128;
fn bounded(p: &[Q]) -> Option<Vec<Q>> {
    (p.len() <= DEGREE + 1
        && p.iter()
            .all(|q| q.numer().bits() <= 4096 && q.denom().bits() <= 4096))
    .then(|| vec_poly::trim(p))
}
fn polynomial(p: Vec<Q>) -> Option<ExponentialPolynomial> {
    Some(BTreeMap::from([(Q::zero(), bounded(&p)?)]))
}
fn add(
    a: &ExponentialPolynomial,
    b: &ExponentialPolynomial,
    sign: &Q,
) -> Option<ExponentialPolynomial> {
    let mut out = a.clone();
    for (r, p) in b {
        let sum = vec_poly::add(
            out.get(r).map_or(&[], Vec::as_slice),
            &vec_poly::scale(p, sign),
        );
        out.insert(r.clone(), bounded(&sum)?);
    }
    out.retain(|_, p| p.iter().any(|c| !c.is_zero()));
    (out.len() <= 16).then_some(out)
}
fn multiply(a: &ExponentialPolynomial, b: &ExponentialPolynomial) -> Option<ExponentialPolynomial> {
    let mut out = ExponentialPolynomial::new();
    for (r, p) in a {
        for (s, q) in b {
            if p.len() + q.len() > DEGREE + 2 {
                return None;
            }
            let rate = r + s;
            let term = vec_poly::mul(p, q);
            let sum = vec_poly::add(out.get(&rate).map_or(&[], Vec::as_slice), &term);
            out.insert(rate, bounded(&sum)?);
            if out.len() > 16 {
                return None;
            }
        }
    }
    out.retain(|_, p| p.iter().any(|c| !c.is_zero()));
    Some(out)
}
fn derivative(a: &ExponentialPolynomial) -> Option<ExponentialPolynomial> {
    let mut out = ExponentialPolynomial::new();
    for (r, p) in a {
        out.insert(
            r.clone(),
            bounded(&vec_poly::add(
                &vec_poly::derivative(p),
                &vec_poly::scale(p, r),
            ))?,
        );
    }
    Some(out)
}
fn finite_polynomial(e: &Expr) -> Option<Vec<Q>> {
    let mut normal = normal(e)?;
    let p = normal.remove(&Q::zero()).unwrap_or_else(|| vec![Q::zero()]);
    normal.is_empty().then_some(p)
}
fn scalar(e: &Expr) -> Option<Q> {
    let p = finite_polynomial(e)?;
    (p.len() <= 1).then(|| p.first().cloned().unwrap_or_default())
}
fn normal(e: &Expr) -> Option<ExponentialPolynomial> {
    match e {
        Expr::Polynomial(p) => polynomial(p.clone()),
        Expr::Variable => polynomial(vec![Q::zero(), Q::one()]),
        Expr::Add(a, b) => add(&normal(a)?, &normal(b)?, &Q::one()),
        Expr::Subtract(a, b) => add(&normal(a)?, &normal(b)?, &-Q::one()),
        Expr::Negate(a) => multiply(&normal(a)?, &polynomial(vec![-Q::one()])?),
        Expr::Multiply(a, b) => multiply(&normal(a)?, &normal(b)?),
        Expr::Divide(a, b) => {
            let c = scalar(b)?;
            if c.is_zero() {
                None
            } else {
                multiply(&normal(a)?, &polynomial(vec![c.recip()])?)
            }
        }
        Expr::Inverse(a) => {
            let c = scalar(a)?;
            (!c.is_zero())
                .then(|| polynomial(vec![c.recip()]))
                .flatten()
        }
        Expr::Power(a, n) => {
            let mut result = polynomial(vec![Q::one()])?;
            let mut base = normal(a)?;
            let mut n = *n;
            while n > 0 {
                if n % 2 == 1 {
                    result = multiply(&result, &base)?;
                }
                n /= 2;
                if n > 0 {
                    base = multiply(&base, &base)?;
                }
            }
            Some(result)
        }
        Expr::Exp(a) => {
            let p = finite_polynomial(a)?;
            if p.len() > 2 || !p[0].is_zero() {
                return None;
            }
            Some(BTreeMap::from([(
                p.get(1).cloned().unwrap_or_default(),
                vec![Q::one()],
            )]))
        }
        Expr::Derivative(a) => derivative(&normal(a)?),
        Expr::Truncate(a, n) if *n <= DEGREE => polynomial(
            a.evaluate(*n, &SeriesExpressionLimits::default())
                .ok()?
                .coefficients()
                .to_vec(),
        ),
        Expr::Compose(a, b) => {
            let p = finite_polynomial(a)?;
            let inner = normal(b)?;
            let mut out = polynomial(vec![Q::zero()])?;
            for c in p.iter().rev() {
                out = add(
                    &multiply(&out, &inner)?,
                    &polynomial(vec![c.clone()])?,
                    &Q::one(),
                )?;
            }
            Some(out)
        }
        _ => None,
    }
}
#[derive(Clone)]
struct Rational {
    n: Vec<Q>,
    d: Vec<Q>,
}
impl Rational {
    fn new(n: Vec<Q>, d: Vec<Q>) -> Option<Self> {
        Some(Self {
            n: bounded(&n)?,
            d: bounded(&d)?,
        })
    }
    fn polynomial(p: Vec<Q>) -> Option<Self> {
        Self::new(p, vec![Q::one()])
    }
    fn multiply(&self, b: &Self) -> Option<Self> {
        if self.n.len() + b.n.len() > DEGREE + 2 || self.d.len() + b.d.len() > DEGREE + 2 {
            return None;
        }
        Self::new(vec_poly::mul(&self.n, &b.n), vec_poly::mul(&self.d, &b.d))
    }
    fn add(&self, b: &Self, sign: &Q) -> Option<Self> {
        if self.n.len() + b.d.len() > DEGREE + 2
            || b.n.len() + self.d.len() > DEGREE + 2
            || self.d.len() + b.d.len() > DEGREE + 2
        {
            return None;
        }
        Self::new(
            vec_poly::add(
                &vec_poly::mul(&self.n, &b.d),
                &vec_poly::scale(&vec_poly::mul(&b.n, &self.d), sign),
            ),
            vec_poly::mul(&self.d, &b.d),
        )
    }
    fn inverse(&self) -> Option<Self> {
        (!self.n.iter().all(Zero::is_zero)).then(|| Self {
            n: self.d.clone(),
            d: self.n.clone(),
        })
    }
    fn derivative(&self) -> Option<Self> {
        Self::new(
            vec_poly::sub(
                &vec_poly::mul(&vec_poly::derivative(&self.n), &self.d),
                &vec_poly::mul(&self.n, &vec_poly::derivative(&self.d)),
            ),
            vec_poly::mul(&self.d, &self.d),
        )
    }
    fn is_zero(&self) -> bool {
        self.n.iter().all(Zero::is_zero)
    }
}
fn rational(e: &Expr) -> Option<Rational> {
    if let Some(p) = finite_polynomial(e) {
        return Rational::polynomial(p);
    }
    match e {
        Expr::Add(a, b) => rational(a)?.add(&rational(b)?, &Q::one()),
        Expr::Subtract(a, b) => rational(a)?.add(&rational(b)?, &-Q::one()),
        Expr::Multiply(a, b) => rational(a)?.multiply(&rational(b)?),
        Expr::Divide(a, b) => rational(a)?.multiply(&rational(b)?.inverse()?),
        Expr::Inverse(a) => rational(a)?.inverse(),
        Expr::Negate(a) => rational(a)?.multiply(&Rational::polynomial(vec![-Q::one()])?),
        Expr::Derivative(a) => rational_derivative(a),
        Expr::Power(a, n) => {
            let mut out = Rational::polynomial(vec![Q::one()])?;
            if *n > DEGREE {
                return None;
            }
            for _ in 0..*n {
                out = out.multiply(&rational(a)?)?;
            }
            Some(out)
        }
        _ => None,
    }
}
fn rational_derivative(e: &Expr) -> Option<Rational> {
    if let Some(r) = rational(e) {
        return r.derivative();
    }
    match e {
        Expr::Log(a) => rational(a)?
            .derivative()?
            .multiply(&rational(a)?.inverse()?),
        Expr::Add(a, b) => rational_derivative(a)?.add(&rational_derivative(b)?, &Q::one()),
        Expr::Subtract(a, b) => rational_derivative(a)?.add(&rational_derivative(b)?, &-Q::one()),
        Expr::Negate(a) => {
            rational_derivative(a)?.multiply(&Rational::polynomial(vec![-Q::one()])?)
        }
        Expr::Multiply(a, b) => {
            if let Some(c) = scalar(a) {
                rational_derivative(b)?.multiply(&Rational::polynomial(vec![c])?)
            } else {
                rational_derivative(a)?.multiply(&Rational::polynomial(vec![scalar(b)?])?)
            }
        }
        Expr::Divide(a, b) => {
            rational_derivative(a)?.multiply(&Rational::polynomial(vec![scalar(b)?.recip()])?)
        }
        _ => None,
    }
}
fn finite_spec(p: &[Q]) -> BlockWeightSpec {
    BlockWeightSpec::Finite {
        weights: p
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, c)| (c * Q::from_integer(factorial_bigint(i))).to_string())
            .collect(),
    }
}
fn ode(e: &Expr, seed: BlockWeightSpec) -> Option<BlockWeightSpec> {
    let values = e.evaluate(8, &SeriesExpressionLimits::default()).ok()?;
    let c = values.coefficients();
    if !c[0].is_zero() {
        return None;
    }
    // Prefix coefficients propose parameters only. Full symbolic identities below
    // certify them, so perturbations beyond the fitting prefix cannot pass.
    let matrix: Vec<_> = (2..8)
        .map(|n| vec![c[n].clone(), Q::from_integer(n.into()) * &c[n]])
        .collect();
    let rhs: Vec<_> = (2..8)
        .map(|n| Q::from_integer((n + 1).into()) * &c[n + 1])
        .collect();
    let solution = exact_algebra::linear_systems::solve_linear_system(&matrix, &rhs)?;
    let (lambda, d) = (solution[0].clone(), solution[1].clone());
    let a = c[1].clone();
    let b = Q::from_integer(2.into()) * &c[2] - (&lambda + &d) * &a;
    // Negative-parameter ODEs cannot supply these closure routes. In
    // particular a PF polynomial may also satisfy such an ODE; leave it
    // available for the separate polynomial-substitution paper route.
    if [&a, &b, &lambda, &d].iter().any(|q| q.is_negative()) {
        return None;
    }
    let factor = vec![Q::one(), -d.clone()];
    let affine = vec![a.clone(), b.clone()];
    let verified = if let Some(f) = normal(e) {
        let residual = add(
            &multiply(&polynomial(factor.clone())?, &derivative(&f)?)?,
            &multiply(&polynomial(vec![lambda.clone()])?, &f)?,
            &-Q::one(),
        )?;
        add(&residual, &polynomial(affine.clone())?, &-Q::one())?.is_empty()
    } else if let Some(df) = rational_derivative(e) {
        let mut residual = df
            .multiply(&Rational::polynomial(factor)?)?
            .add(&Rational::polynomial(affine)?, &-Q::one())?;
        if !lambda.is_zero() {
            residual = residual.add(
                &rational(e)?.multiply(&Rational::polynomial(vec![lambda.clone()])?)?,
                &-Q::one(),
            )?;
        }
        residual.is_zero()
    } else {
        false
    };
    if !verified {
        return None;
    }
    if d.is_zero() {
        Some(BlockWeightSpec::OdeSubstitution {
            seed: Box::new(seed),
            a: a.to_string(),
            b: b.to_string(),
            lambda: lambda.to_string(),
        })
    } else {
        Some(BlockWeightSpec::RationalOdeSubstitution {
            seed: Box::new(seed),
            a: a.to_string(),
            b: b.to_string(),
            lambda: lambda.to_string(),
            d: d.to_string(),
        })
    }
}
// Replace an explicitly occurring inner subtree by the formal outer variable.
// Reject every unmatched variable: it would introduce an independent variable
// and would not prove a univariate composition identity.
fn outer_expression(e: &Expr, inner: &Expr) -> Option<Expr> {
    if e == inner {
        return Some(Expr::Variable);
    }
    let child = |a: &Expr| outer_expression(a, inner).map(Box::new);
    Some(match e {
        Expr::Polynomial(p) => Expr::Polynomial(p.clone()),
        Expr::Variable => return None,
        Expr::Add(a, b) => Expr::Add(child(a)?, child(b)?),
        Expr::Subtract(a, b) => Expr::Subtract(child(a)?, child(b)?),
        Expr::Multiply(a, b) => Expr::Multiply(child(a)?, child(b)?),
        Expr::Divide(a, b) => Expr::Divide(child(a)?, child(b)?),
        Expr::Negate(a) => Expr::Negate(child(a)?),
        Expr::Power(a, n) => Expr::Power(child(a)?, *n),
        Expr::Exp(a) => Expr::Exp(child(a)?),
        Expr::Log(a) => Expr::Log(child(a)?),
        Expr::Inverse(a) => Expr::Inverse(child(a)?),
        // Differentiation/truncation do not commute with substitution.
        _ => return None,
    })
}
fn implicit_substitution(e: &Expr, whole: &Expr) -> Option<BlockWeightSpec> {
    if let Expr::Exp(inner) = e {
        if !matches!(inner.as_ref(), Expr::Variable) {
            if let (Some(outer), Some(seed)) = (outer_expression(whole, inner), recognize(inner)) {
                if let Some(spec) = ode(&outer, seed) {
                    return Some(spec);
                }
            }
        }
    }
    let children: Vec<&Expr> = match e {
        Expr::Add(a, b) | Expr::Subtract(a, b) | Expr::Multiply(a, b) | Expr::Divide(a, b) => {
            vec![a, b]
        }
        Expr::Negate(a) | Expr::Power(a, _) | Expr::Exp(a) | Expr::Log(a) | Expr::Inverse(a) => {
            vec![a]
        }
        _ => vec![],
    };
    children
        .into_iter()
        .find_map(|a| implicit_substitution(a, whole))
}
fn recognize(e: &Expr) -> Option<BlockWeightSpec> {
    let finite = finite_polynomial(e);
    // Prefer the primitive low-degree seed when the whole identity already
    // reduces to it; otherwise retain an explicit substitution witness.
    if let Some(p) = finite.as_ref().filter(|p| p.len() <= 5 && p[0].is_zero()) {
        return Some(finite_spec(p));
    }
    // Preserve the supplied composition before polynomial expansion, which
    // would otherwise discard the polynomial-substitution proof witness.
    if let Expr::Compose(outer, inner) = e {
        let seed = recognize(inner)?;
        if let Some(spec) = ode(outer, seed.clone()) {
            return Some(spec);
        }
        if let Some(p) = finite_polynomial(outer) {
            return Some(BlockWeightSpec::PolynomialSubstitution {
                seed: Box::new(seed),
                coefficients: p.iter().map(ToString::to_string).collect(),
            });
        }
    }
    if let Some(p) = finite {
        if p[0].is_zero() {
            return Some(finite_spec(&p));
        }
        return None;
    }
    if let Some(mut f) = normal(e) {
        let p = f.remove(&Q::zero()).unwrap_or_else(|| vec![Q::zero()]);
        if f.len() == 1 {
            let (rate, scale) = f.first_key_value()?;
            if scale.len() == 1 && !rate.is_zero() {
                let minimum = p.len().max(1);
                let mut expected = vec![Q::zero(); minimum];
                let mut coefficient = scale[0].clone();
                for (i, c) in expected.iter_mut().enumerate() {
                    if i > 0 {
                        coefficient *= rate;
                        coefficient /= Q::from_integer(i.into());
                    }
                    *c = -coefficient.clone();
                }
                if vec_poly::trim(&p) == vec_poly::trim(&expected) {
                    return Some(BlockWeightSpec::Exponential {
                        scale: scale[0].to_string(),
                        rate: rate.to_string(),
                        minimum,
                    });
                }
            }
        }
    }
    ode(e, BlockWeightSpec::Sizes { sizes: vec![1] }).or_else(|| implicit_substitution(e, e))
}
/// Recognize only bounded exact polynomial/exponential/rational identities.
pub(crate) fn recognize_expression(expression: &SeriesExpression) -> Option<BlockWeightSpec> {
    recognize(&expression.root)
}
