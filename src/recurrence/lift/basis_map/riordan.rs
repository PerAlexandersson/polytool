//! Index-side transport through a Riordan coefficient family (module
//! documentation of `basis_map`, "Index-side transport"): a scalar rule for
//! `H_n = T_P(R_n)` when the basis rule of `P` depends on the index (so
//! `T_P` has no companion form) and the coefficient family has an
//! index-free rule whose generating function is `g(z) / (1 - u f(z))`.

use super::map::{apply_with_rows, BasisMap};
use crate::recurrence::implicit::algebra::{Biv, BivarConversion, DOp};
use crate::recurrence::lift::hadamard_factor::prepare::PreparedLiftRule;
use crate::recurrence::RecurrenceError;
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

type Q = BigRational;

/// Rows of the image checked against the transported rule after its
/// construction (a regression guard; the rule is exact by construction).
const SELF_CHECK_BASES: usize = 4;

fn unsupported(reason: &str) -> RecurrenceError {
    RecurrenceError::Unsupported(reason.to_string())
}

// ---------------------------------------------------------------------------
// Univariate polynomials and rational functions in z
// ---------------------------------------------------------------------------

/// Quotient and remainder of `a` by the nonzero `b`.
fn div_rem(a: &[Q], b: &[Q]) -> (Vec<Q>, Vec<Q>) {
    let b = vec_poly::trim(b);
    let db = vec_poly::degree(&b).unwrap_or(0);
    let lead = b[db].clone();
    let mut r = vec_poly::trim(a);
    let mut quotient = vec![Q::zero()];
    while !vec_poly::is_zero(&r) && vec_poly::degree(&r).unwrap_or(0) >= db {
        let dr = vec_poly::degree(&r).unwrap_or(0);
        let c = &r[dr] / &lead;
        let term = vec_poly::shift_by(&[c], dr - db);
        quotient = vec_poly::add(&quotient, &term);
        r = vec_poly::trim(&vec_poly::sub(&r, &vec_poly::mul(&term, &b)));
    }
    (vec_poly::trim(&quotient), r)
}

fn monic(p: &[Q]) -> Vec<Q> {
    let p = vec_poly::trim(p);
    match p.last() {
        Some(lead) if !lead.is_zero() => vec_poly::scale(&p, &lead.recip()),
        _ => p,
    }
}

/// Monic greatest common divisor (`[0]` for two zeros).
fn gcd(a: &[Q], b: &[Q]) -> Vec<Q> {
    let mut x = monic(a);
    let mut y = monic(b);
    while !vec_poly::is_zero(&y) {
        let (_, r) = div_rem(&x, &y);
        x = y;
        y = monic(&r);
    }
    monic(&x)
}

fn exact_div(a: &[Q], b: &[Q]) -> Vec<Q> {
    let (q, r) = div_rem(a, b);
    debug_assert!(vec_poly::is_zero(&r));
    q
}

/// `num / den` in lowest terms with a monic denominator.
#[derive(Clone, Debug, PartialEq)]
struct Rat {
    num: Vec<Q>,
    den: Vec<Q>,
}

impl Rat {
    fn new(num: &[Q], den: &[Q]) -> Self {
        let num = vec_poly::trim(num);
        if vec_poly::is_zero(&num) {
            return Self::zero();
        }
        let g = gcd(&num, den);
        let num = exact_div(&num, &g);
        let den = exact_div(den, &g);
        let lead = vec_poly::trim(&den).last().cloned().unwrap_or_else(Q::one);
        Self {
            num: vec_poly::scale(&num, &lead.recip()),
            den: vec_poly::scale(&den, &lead.recip()),
        }
    }

    fn zero() -> Self {
        Self {
            num: vec![Q::zero()],
            den: vec![Q::one()],
        }
    }

    fn constant(c: Q) -> Self {
        Self::new(&[c], &[Q::one()])
    }

    fn is_zero(&self) -> bool {
        vec_poly::is_zero(&self.num)
    }

    fn add(&self, other: &Self) -> Self {
        Self::new(
            &vec_poly::add(
                &vec_poly::mul(&self.num, &other.den),
                &vec_poly::mul(&other.num, &self.den),
            ),
            &vec_poly::mul(&self.den, &other.den),
        )
    }

    fn mul(&self, other: &Self) -> Self {
        Self::new(
            &vec_poly::mul(&self.num, &other.num),
            &vec_poly::mul(&self.den, &other.den),
        )
    }

    fn scale(&self, c: &Q) -> Self {
        Self::new(&vec_poly::scale(&self.num, c), &self.den)
    }

    /// `theta = z d/dz` applied.
    fn theta(&self) -> Self {
        let numerator = vec_poly::sub(
            &vec_poly::mul(&vec_poly::derivative(&self.num), &self.den),
            &vec_poly::mul(&self.num, &vec_poly::derivative(&self.den)),
        );
        Self::new(
            &vec_poly::shift_by(&numerator, 1),
            &vec_poly::mul(&self.den, &self.den),
        )
    }
}

/// An operator `sum_b r_b(z) theta^b` with rational coefficients.
type ZOp = BTreeMap<usize, Rat>;

fn zop_add(a: &ZOp, b: &ZOp) -> ZOp {
    let mut out = a.clone();
    for (k, v) in b {
        let entry = out.entry(*k).or_insert_with(Rat::zero);
        *entry = entry.add(v);
    }
    out.retain(|_, v| !v.is_zero());
    out
}

/// `(r theta^b) o (s theta^c) = r sum_i binom(b, i) theta^i(s) theta^(b - i + c)`.
fn zop_compose(a: &ZOp, b: &ZOp) -> ZOp {
    let mut out = ZOp::new();
    for (&p, r) in a {
        for (&c, s) in b {
            let mut derivative = s.clone();
            for i in 0..=p {
                let coefficient = r
                    .mul(&derivative)
                    .scale(&Q::from_integer(binomial(p, i).into()));
                let entry = out.entry(p - i + c).or_insert_with(Rat::zero);
                *entry = entry.add(&coefficient);
                derivative = derivative.theta();
            }
        }
    }
    out.retain(|_, v| !v.is_zero());
    out
}

fn binomial(n: usize, k: usize) -> u128 {
    let mut value: u128 = 1;
    for i in 0..k {
        value = value * (n - i) as u128 / (i + 1) as u128;
    }
    value
}

// ---------------------------------------------------------------------------
// The Riordan form of the coefficient family
// ---------------------------------------------------------------------------

/// `sum_(m >= 0) R_(n0 + m)(u) z^m = N(z) / (a(z) - u b(z))`.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct RiordanForm {
    /// `a`, `a(0) != 0`.
    pub(super) a: Vec<Q>,
    /// `b`, `b(0) = 0`, nonzero.
    pub(super) b: Vec<Q>,
    /// `N`, free of `u`, nonzero.
    pub(super) numerator: Vec<Q>,
}

/// The Riordan form of the coefficient family: its relation must be free of
/// `n` and of derivatives, with a constant leader and coefficients of
/// degree at most one in `u`, and the numerator `N` (determined by the
/// identifying rows) must be free of `u`.
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] naming the first condition that fails.
pub(super) fn riordan_form(coefficient: &PreparedLiftRule) -> Result<RiordanForm, RecurrenceError> {
    let relation = &coefficient.relation;
    let k = relation.max_shift();
    if k == 0 {
        return Err(unsupported("the coefficient rule has no positive shift"));
    }
    let mut a = vec![Q::zero(); k + 1];
    let mut b = vec![Q::zero(); k + 1];
    for (&(shift, derivative), e) in &relation.0 {
        if derivative > 0 {
            return Err(unsupported(
                "the coefficient rule has derivatives in u (not a Riordan family)",
            ));
        }
        if !e.is_index_free() {
            return Err(unsupported(
                "the coefficient rule depends on the index (not a Riordan family)",
            ));
        }
        if e.variable_degree().unwrap_or(0) > 1 {
            return Err(unsupported(
                "the coefficient rule has a coefficient of degree above one in u (not a Riordan \
                 family)",
            ));
        }
        let row = e.rows().first().cloned().unwrap_or_default();
        a[k - shift] += vec_poly::coeff(&row, 0);
        b[k - shift] -= vec_poly::coeff(&row, 1);
    }
    if a[0].is_zero() || !b[0].is_zero() {
        return Err(unsupported(
            "the coefficient rule's leader is not a nonzero constant",
        ));
    }
    let a = vec_poly::trim(&a);
    let b = vec_poly::trim(&b);
    if vec_poly::is_zero(&b) {
        return Err(unsupported(
            "the coefficient rule is free of u (constant rows)",
        ));
    }
    // N = (a - u b) R~ (mod z^L); it has degree below the number of seeds.
    // The coefficient of z^m is the relation at base n0 + m - k, which holds
    // from the first valid base on.
    let identifying = coefficient.valid_from_base() - coefficient.first_index + k;
    let count = identifying + k + 3;
    let rows = coefficient.rows(count)?;
    let mut numerator = Vec::with_capacity(identifying);
    for m in 0..count {
        let mut value = vec![Q::zero()];
        for (i, ai) in a.iter().enumerate().take(m + 1) {
            value = vec_poly::add(&value, &vec_poly::scale(&rows[m - i], ai));
        }
        for (i, bi) in b.iter().enumerate().take(m + 1) {
            value = vec_poly::sub(
                &value,
                &vec_poly::shift_by(&vec_poly::scale(&rows[m - i], bi), 1),
            );
        }
        let value = vec_poly::trim(&value);
        if m >= identifying {
            if !vec_poly::is_zero(&value) {
                return Err(RecurrenceError::VerificationFailed(
                    "the coefficient rows do not satisfy their generating-function identity".into(),
                ));
            }
        } else {
            if vec_poly::degree(&value).unwrap_or(0) > 0 {
                return Err(unsupported(
                    "the generating function's numerator depends on u (not a Riordan family)",
                ));
            }
            numerator.push(value[0].clone());
        }
    }
    let numerator = vec_poly::trim(&numerator);
    if vec_poly::is_zero(&numerator) {
        return Err(unsupported("the coefficient family is zero"));
    }
    Ok(RiordanForm { a, b, numerator })
}

// ---------------------------------------------------------------------------
// The transported rule
// ---------------------------------------------------------------------------

/// The scalar rule `Phi` of `H_n = T_P(R_n)` and the first base from which
/// it holds.
#[derive(Clone, Debug)]
pub(super) struct RiordanTransport {
    /// The relation at base `n` (absolute index).
    pub(super) relation: DOp,
    /// The first base at which `(Phi H)_n = 0` is proved.
    pub(super) valid_from_base: usize,
}

/// `p(n + c)^b` as rows in `(n, t)` times `c_t t^j`.
fn shifted_power_term(value: &Q, j: usize, b: usize, c: &Q) -> Biv {
    let mut rows = vec![Vec::new(); b + 1];
    for (i, row) in rows.iter_mut().enumerate() {
        let mut coefficient = value * Q::from_integer(binomial(b, i).into());
        for _ in 0..b - i {
            coefficient *= c;
        }
        let mut t_row = vec![Q::zero(); j + 1];
        t_row[j] = coefficient;
        *row = t_row;
    }
    Biv::new(rows)
}

/// The index-side transport (module documentation of `basis_map`,
/// "Index-side transport").
///
/// # Errors
///
/// [`RecurrenceError::Unsupported`] when the coefficient family is not a
/// Riordan family, and [`RecurrenceError::VerificationFailed`] if the
/// rule fails the regression check on computed rows.
pub(super) fn riordan_transport(
    coefficient: &PreparedLiftRule,
    map: &BasisMap,
) -> Result<RiordanTransport, RecurrenceError> {
    let form = riordan_form(coefficient)?;
    let n0 = coefficient.first_index;
    let RiordanForm { a, b, numerator } = &form;
    let one = vec![Q::one()];
    // rho = f / (z f') = a (b / z) / (b' a - b a'), gamma = z g' / g.
    let b_over_z = vec_poly::exact_div_by_t(b).ok_or_else(|| unsupported("b(0) != 0"))?;
    let wronskian = vec_poly::sub(
        &vec_poly::mul(&vec_poly::derivative(b), a),
        &vec_poly::mul(b, &vec_poly::derivative(a)),
    );
    let rho = Rat::new(&vec_poly::mul(a, &b_over_z), &wronskian);
    let gamma = Rat::new(
        &vec_poly::shift_by(
            &vec_poly::sub(
                &vec_poly::mul(&vec_poly::derivative(numerator), a),
                &vec_poly::mul(numerator, &vec_poly::derivative(a)),
            ),
            1,
        ),
        &vec_poly::mul(numerator, a),
    );
    // K = rho (theta - gamma), the image of theta_w.
    let mut k_hat = ZOp::new();
    k_hat.insert(1, rho.clone());
    let constant = rho.mul(&gamma).scale(&-Q::one());
    if !constant.is_zero() {
        k_hat.insert(0, constant);
    }
    let mut powers: BTreeMap<String, Vec<ZOp>> = BTreeMap::new();
    let mut shifted_power = |shift: &Q, e: usize| -> ZOp {
        let list = powers.entry(shift.to_string()).or_insert_with(|| {
            let mut base = k_hat.clone();
            if !shift.is_zero() {
                base = zop_add(&base, &ZOp::from([(0, Rat::constant(shift.clone()))]));
            }
            vec![ZOp::from([(0, Rat::constant(Q::one()))]), base]
        });
        while list.len() <= e {
            let next = zop_compose(&list[list.len() - 1], &list[1]);
            list.push(next);
        }
        list[e].clone()
    };
    // Lambda-hat = f(i0 + K, t) - sum_terms f^r c(i0 + r + K, t) D^d,
    // collected as (theta power, derivative, t power) -> rational in z.
    let rule = map.rule();
    let i0 = Q::from_integer(map.first_index().into());
    let mut operator: BTreeMap<(usize, usize, usize), Rat> = BTreeMap::new();
    let mut accumulate = |zop: &ZOp, factor: &Rat, t_poly: &[Q], d: usize| {
        for (&power, r) in zop {
            let base = r.mul(factor);
            for (j, c) in t_poly.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                let entry = operator.entry((power, d, j)).or_insert_with(Rat::zero);
                *entry = entry.add(&base.scale(c));
            }
        }
    };
    let leader = rule
        .denominator
        .as_ref()
        .map(Biv::from_bivar)
        .unwrap_or_else(Biv::one);
    for (e, t_poly) in leader.rows().iter().enumerate() {
        let zop = shifted_power(&i0, e);
        accumulate(&zop, &Rat::constant(Q::one()), t_poly, 0);
    }
    for term in &rule.terms {
        let r = term.offset;
        let f_power = Rat::new(&vec_poly::pow(b, r), &vec_poly::pow(a, r)).scale(&-Q::one());
        let shift = &i0 + Q::from_integer(r.into());
        for (e, t_poly) in Biv::from_bivar(&term.coeff).rows().iter().enumerate() {
            let zop = shifted_power(&shift, e);
            accumulate(&zop, &f_power, t_poly, term.deriv_order);
        }
    }
    operator.retain(|_, v| !v.is_zero());
    // pi(w) = sum_(j < seeds) w^j (f P_j - sum c D^d P_(j - r)) at n = i0 + j;
    // the right-hand side g pi(f) = sum_j pi_j N b^j / a^(j + 1).
    let seed_count = map.seeds().len();
    let p_rows = map.basis_rows(seed_count)?;
    let mut rhs: BTreeMap<usize, Rat> = BTreeMap::new();
    for (j, p_j) in p_rows.iter().enumerate() {
        let n = &i0 + Q::from_integer(j.into());
        let mut value = vec_poly::mul(&leader.evaluate_index(&n), p_j);
        for term in &rule.terms {
            if j < term.offset {
                continue;
            }
            let mut row = p_rows[j - term.offset].clone();
            for _ in 0..term.deriv_order {
                row = vec_poly::derivative(&row);
            }
            value = vec_poly::sub(
                &value,
                &vec_poly::mul(&Biv::from_bivar(&term.coeff).evaluate_index(&n), &row),
            );
        }
        let factor = Rat::new(
            &vec_poly::mul(numerator, &vec_poly::pow(b, j)),
            &vec_poly::pow(a, j + 1),
        );
        for (power, c) in vec_poly::trim(&value).iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            let entry = rhs.entry(power).or_insert_with(Rat::zero);
            *entry = entry.add(&factor.scale(c));
        }
    }
    rhs.retain(|_, v| !v.is_zero());
    // Clear denominators (left multiplication by a polynomial in z).
    let mut common = one.clone();
    for value in operator.values().chain(rhs.values()) {
        let g = gcd(&common, &value.den);
        common = vec_poly::mul(&common, &exact_div(&value.den, &g));
    }
    let clear = |value: &Rat| vec_poly::mul(&value.num, &exact_div(&common, &value.den));
    let mut polynomials: BTreeMap<(usize, usize, usize), Vec<Q>> =
        operator.iter().map(|(k, v)| (*k, clear(v))).collect();
    let mut psi: Vec<Vec<Q>> = rhs.values().map(clear).collect();
    // Remove a common polynomial factor of every coefficient and of psi.
    let mut content = vec![Q::zero()];
    for p in polynomials.values().chain(psi.iter()) {
        content = gcd(&content, p);
    }
    if vec_poly::degree(&content).unwrap_or(0) > 0 {
        for p in polynomials.values_mut().chain(psi.iter_mut()) {
            *p = exact_div(p, &content);
        }
    }
    if polynomials.is_empty() {
        return Err(unsupported("the transported operator vanishes"));
    }
    let top = polynomials
        .values()
        .filter_map(|p| vec_poly::degree(p))
        .max()
        .unwrap_or(0);
    let psi_degree = psi.iter().filter_map(|p| vec_poly::degree(p)).max();
    // [z^M] (Omega H~) = sum c(t) (M - a)^b D^d H_(n0 + M - a) = psi_M, zero
    // for M > deg psi; at base n the indices are n + top - a.
    let mut relation = DOp::default();
    let n0q = Q::from_integer(n0.into());
    for (&(power, d, j), p) in &polynomials {
        for (a_power, value) in p.iter().enumerate() {
            if value.is_zero() {
                continue;
            }
            let shift = top - a_power;
            let c = Q::from_integer(shift.into()) - &n0q;
            relation.add_term((shift, d), shifted_power_term(value, j, power, &c));
        }
    }
    let valid_from_base = match psi_degree {
        Some(degree) if degree + 1 > top => n0 + degree + 1 - top,
        _ => n0,
    };
    // Regression check on computed rows.
    let order = relation.max_shift();
    let count = valid_from_base - n0 + order + SELF_CHECK_BASES;
    let r_rows = coefficient.rows(count)?;
    let width = r_rows.iter().map(Vec::len).max().unwrap_or(1);
    let basis = map.basis_rows(width)?;
    let h_rows: Vec<Vec<Q>> = r_rows.iter().map(|r| apply_with_rows(r, &basis)).collect();
    for base in valid_from_base..valid_from_base + SELF_CHECK_BASES {
        let value = relation.apply(&h_rows, n0, base).ok_or_else(|| {
            RecurrenceError::VerificationFailed("missing image rows in the self-check".into())
        })?;
        if !vec_poly::is_zero(&value) {
            return Err(RecurrenceError::VerificationFailed(format!(
                "the index-side transported rule fails on the image rows at base {base}"
            )));
        }
    }
    Ok(RiordanTransport {
        relation,
        valid_from_base,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::lift::hadamard_factor::prepare::prepare_lift_rule;
    use crate::recurrence::{BivarPoly, Recurrence, RecurrenceTerm};

    fn q(v: i64) -> Q {
        Q::from_integer(v.into())
    }

    fn bivar(rows: Vec<Vec<i64>>) -> BivarPoly {
        BivarPoly::from_coefficients_normalized(
            rows.into_iter()
                .map(|r| r.into_iter().map(q).collect())
                .collect(),
        )
    }

    /// `R_n = (1 + u)^n`.
    fn binomial_family() -> PreparedLiftRule {
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(1, 0, bivar(vec![vec![1, 1]]))],
            denominator: None,
            inhomogeneous: None,
        };
        prepare_lift_rule(&rule, &[vec![q(1)]], 0, "coefficient family").unwrap()
    }

    /// `E_0 = 1`, `E_n = n t E_(n-1) + t (1 - t) E'_(n-1)`.
    fn eulerian_map() -> BasisMap {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(1, 0, bivar(vec![vec![0], vec![0, 1]])),
                RecurrenceTerm::new(1, 1, bivar(vec![vec![0, 1, -1]])),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        BasisMap::new("eulerian", rule, 0, vec![vec![q(1)]]).unwrap()
    }

    #[test]
    fn rational_functions_reduce() {
        let r = Rat::new(&[q(-1), q(0), q(1)], &[q(-2), q(2)]);
        assert_eq!(
            r.num,
            vec![Q::new(1.into(), 2.into()), Q::new(1.into(), 2.into())]
        );
        assert_eq!(r.den, vec![q(1)]);
        // theta(1 / (1 - z)) = z / (1 - z)^2.
        let s = Rat::new(&[q(1)], &[q(1), q(-1)]).theta();
        assert_eq!(s, Rat::new(&[q(0), q(1)], &[q(1), q(-2), q(1)]));
    }

    #[test]
    fn binomial_family_has_its_riordan_form() {
        let form = riordan_form(&binomial_family()).unwrap();
        assert_eq!(form.a, vec![q(1), q(-1)]);
        assert_eq!(form.b, vec![q(0), q(1)]);
        assert_eq!(form.numerator, vec![q(1)]);
    }

    #[test]
    fn binomial_eulerian_transport_is_a046802() {
        // A046802: H_n = (1 + n t) H_(n-1) + t (1 - t) H'_(n-1) + (t - n t) H_(n-2).
        let transport = riordan_transport(&binomial_family(), &eulerian_map()).unwrap();
        assert_eq!(transport.valid_from_base, 0);
        assert_eq!(transport.relation.max_shift(), 2);
        assert_eq!(transport.relation.max_derivative(), 1);
    }

    #[test]
    fn index_dependent_coefficient_families_are_not_riordan() {
        // R_n = (n + u) R_(n-1).
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(1, 0, bivar(vec![vec![0, 1], vec![1]]))],
            denominator: None,
            inhomogeneous: None,
        };
        let prepared = prepare_lift_rule(&rule, &[vec![q(1)]], 0, "coefficient family").unwrap();
        assert!(matches!(
            riordan_form(&prepared),
            Err(RecurrenceError::Unsupported(_))
        ));
    }
}
