//! Invariant row annihilators of first-order parents: a shift-free operator
//! `A(n)` that kills every parent row, certified by an exact invariance
//! identity under the parent's one-step operator (module documentation of
//! `hadamard_factor`, section "Row annihilators").
//!
//! The candidate is found as a right greatest common divisor in the Weyl
//! algebra: the division remainder `r(n)` kills every row (on the checked
//! rows), and so does `r(n + 1) o T(n)`; their right gcd is a smaller
//! annihilator, and the loop stops at the first candidate whose step
//! `A(n + 1) o T(n)` is a left multiple of `A(n)`.  Nothing in the proof
//! depends on how the candidate was found: the replay checks the identities
//! and the rows only.

use super::prepare::left_multiply;
use crate::polynomial::Polynomial;
use crate::recurrence::implicit::algebra::{content_and_primitive, roots_in_domain, Biv, DOp};
use crate::recurrence::lift::section_certify::x_free_leader;
use crate::recurrence::lift::section_operators::check_term_budget;
use crate::recurrence::RecurrenceError;
use crate::vec_poly;
use crate::BigRational;
use num_traits::{One, Zero};

type Q = BigRational;

/// `A(n)` with its invariance certificate
/// `nu(n, x) A(n + 1) o T(n) = B(n) o A(n)`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct InvariantAnnihilator {
    /// `A(n)`, shift-free, of positive order.
    pub(super) annihilator: DOp,
    /// `nu(n, x)`.
    pub(super) step_multiplier: Biv,
    /// `B(n)`, shift-free.
    pub(super) step_quotient: DOp,
}

/// `T(n)` of a first-order parent relation `g(n) S + C_0(n)` (leader free of
/// `x` and `D`): `T = -C_0`, so `g(n) P_(n+1) = T(n) P_n` wherever the
/// relation holds.  `None` for parents of other orders.
pub(super) fn step_operator(parent: &DOp) -> Option<DOp> {
    let (k, _) = x_free_leader(parent).ok()?;
    (k == 1).then(|| parent.column(0).scale(&-Q::one()))
}

/// `A(n + by)`: every coefficient shifted in `n`.
pub(super) fn shift_coefficients(op: &DOp, by: i64) -> DOp {
    let mut out = DOp::default();
    for (key, coefficient) in &op.0 {
        out.add_term(*key, coefficient.shift_n(by));
    }
    out
}

/// `nu A(n + 1) o T(n) - B o A(n)` vanishes.
pub(super) fn invariance_holds(invariant: &InvariantAnnihilator, step: &DOp) -> bool {
    let left = left_multiply(
        &invariant.step_multiplier,
        &shift_coefficients(&invariant.annihilator, 1).compose(step),
    );
    left.sub(&invariant.step_quotient.compose(&invariant.annihilator))
        .is_zero()
}

/// The integer bases `n >= from` where `p(n, .)` vanishes identically (the
/// integer roots of the gcd of its `x`-columns); empty for `p` free of `n`.
pub(super) fn vanishing_bases(p: &Biv, from: usize) -> Vec<usize> {
    let (content, _) = content_and_primitive(p);
    if vec_poly::degree(&content).unwrap_or(0) == 0 {
        return Vec::new();
    }
    roots_in_domain(&content, from, None)
}

/// The gcd of the rows of `p` (its coefficients as polynomials in `x`),
/// monic; zero for zero.
fn x_content(p: &Biv) -> Vec<Q> {
    let mut content: Polynomial<Q> = Polynomial::zero();
    for row in p.rows() {
        content = content.gcd(&Polynomial::new(row.clone()));
    }
    vec_poly::trim(content.coeffs())
}

/// `p` with every row divided by the polynomial `divisor` in `x` (exact by
/// construction).
fn divide_rows(p: &Biv, divisor: &[Q]) -> Biv {
    let divisor = Polynomial::new(divisor.to_vec());
    Biv::new(
        p.rows()
            .iter()
            .map(|row| {
                if row.iter().all(Zero::is_zero) {
                    Vec::new()
                } else {
                    vec_poly::trim(Polynomial::new(row.clone()).div_rem(&divisor).0.coeffs())
                }
            })
            .collect(),
    )
}

/// The degree in `n`, `None` for zero.
fn n_degree(p: &Biv) -> Option<usize> {
    p.rows().len().checked_sub(1)
}

/// The pseudo-remainder of `a` by `b` as polynomials in `n` over `Q[x]`.
fn n_pseudo_remainder(a: &Biv, b: &Biv) -> Biv {
    let Some(b_degree) = n_degree(b) else {
        return a.clone();
    };
    let b_lead = Biv::from_y(&b.rows()[b_degree]);
    let mut remainder = a.clone();
    while let Some(degree) = n_degree(&remainder).filter(|d| *d >= b_degree) {
        let lead = Biv::new(vec![remainder.rows()[degree].clone()]);
        let mut shifted_lead = vec![Vec::new(); degree - b_degree];
        shifted_lead.push(lead.rows().first().cloned().unwrap_or_default());
        remainder = b_lead
            .mul(&remainder)
            .add(&Biv::new(shifted_lead).mul(b).scale(&-Q::one()));
    }
    remainder
}

/// A gcd of two bivariate polynomials up to a rational factor (primitive
/// polynomial remainder sequence in `n` over `Q[x]`).
fn bivariate_gcd(a: &Biv, b: &Biv) -> Biv {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    let (a_content, b_content) = (x_content(a), x_content(b));
    let content = Polynomial::new(a_content.clone()).gcd(&Polynomial::new(b_content.clone()));
    let content = Biv::from_y(content.coeffs());
    let (mut left, mut right) = (divide_rows(a, &a_content), divide_rows(b, &b_content));
    if n_degree(&left) < n_degree(&right) {
        std::mem::swap(&mut left, &mut right);
    }
    loop {
        if n_degree(&right) == Some(0) {
            return content;
        }
        let remainder = n_pseudo_remainder(&left, &right);
        if remainder.is_zero() {
            return content.mul(&right);
        }
        left = right;
        right = divide_rows(&remainder, &x_content(&remainder));
    }
}

/// `p / divisor` when the division is exact (lexicographic division, `n`
/// before `x`), `None` otherwise.
fn bivariate_exact_quotient(p: &Biv, divisor: &Biv) -> Option<Biv> {
    let divisor_degree = n_degree(divisor)?;
    let divisor_row = &divisor.rows()[divisor_degree];
    let divisor_x_degree = vec_poly::degree(divisor_row)?;
    let divisor_lead = &divisor_row[divisor_x_degree];
    let mut remainder = p.clone();
    let mut quotient = Biv::zero();
    while let Some(degree) = n_degree(&remainder) {
        let row = &remainder.rows()[degree];
        let x_degree = vec_poly::degree(row)?;
        if degree < divisor_degree || x_degree < divisor_x_degree {
            return None;
        }
        let term = Biv::monomial(
            degree - divisor_degree,
            x_degree - divisor_x_degree,
            &row[x_degree] / divisor_lead,
        );
        remainder = remainder.add(&term.mul(divisor).scale(&-Q::one()));
        quotient = quotient.add(&term);
    }
    Some(quotient)
}

/// The operator divided by the gcd of its coefficients (bivariate
/// polynomials in `n` and `x`) and made monic in its leading rational
/// coefficient.  Only shortens the coefficients; the proof does not depend
/// on it.
pub(super) fn primitive_operator(op: &DOp) -> DOp {
    let common = op.0.values().fold(Biv::zero(), |gcd, coefficient| {
        bivariate_gcd(&gcd, coefficient)
    });
    let mut out = DOp::default();
    for (key, coefficient) in &op.0 {
        match bivariate_exact_quotient(coefficient, &common) {
            Some(quotient) => out.add_term(*key, quotient),
            None => return op.clone(),
        }
    }
    match out.leading() {
        Some(lead) if !lead.is_zero() => out.scale(&lead.recip()),
        _ => out,
    }
}

/// `lead^j f = Z o g + rem` with `ord rem < ord g` (the shared right
/// pseudo-division of the Weyl algebra).
fn right_remainder(f: &DOp, g: &DOp, max_terms: usize) -> Result<DOp, RecurrenceError> {
    Ok(super::certify::right_pseudo_divide(f, g, max_terms)?.2)
}

/// A right gcd of two shift-free operators over `Q(n, x)`, up to a
/// coefficient factor (primitive parts along the way).
fn right_gcd(f: &DOp, g: &DOp, max_terms: usize) -> Result<DOp, RecurrenceError> {
    let (mut a, mut b) = if f.max_derivative() >= g.max_derivative() {
        (primitive_operator(f), primitive_operator(g))
    } else {
        (primitive_operator(g), primitive_operator(f))
    };
    while !b.is_zero() {
        let remainder = primitive_operator(&right_remainder(&a, &b, max_terms)?);
        check_term_budget(remainder.0.len(), max_terms)?;
        a = b;
        b = remainder;
    }
    Ok(a)
}

/// Search an invariant annihilator dividing `remainder` (module
/// documentation): `None` when the gcd loop reaches order zero or an
/// operator without a derivative.
pub(super) fn find_invariant_annihilator(
    remainder: &DOp,
    step: &DOp,
    max_terms: usize,
) -> Result<Option<InvariantAnnihilator>, RecurrenceError> {
    let mut candidate = primitive_operator(remainder);
    for _ in 0..=remainder.max_derivative() {
        if candidate.is_zero() || candidate.max_derivative() == 0 {
            return Ok(None);
        }
        let next = shift_coefficients(&candidate, 1).compose(step);
        let (step_multiplier, step_quotient, rest) =
            super::certify::right_pseudo_divide(&next, &candidate, max_terms)?;
        if rest.is_zero() {
            return Ok(Some(InvariantAnnihilator {
                annihilator: candidate,
                step_multiplier,
                step_quotient,
            }));
        }
        candidate = right_gcd(&candidate, &rest, max_terms)?;
    }
    Ok(None)
}

/// `mu(n, x)` as `c(n) p(x)` when its primitive part (over `Q[n]`) is free
/// of `n`: `(c, p)`.
pub(super) fn separated_multiplier(mu: &Biv) -> Option<(Vec<Q>, Vec<Q>)> {
    let (content, primitive) = content_and_primitive(mu);
    if mu.is_zero() || !primitive.is_n_free() {
        return None;
    }
    Some((
        content,
        primitive.rows().first().cloned().unwrap_or_default(),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::implicit::algebra::q;

    fn biv(rows: &[&[i64]]) -> Biv {
        Biv::new(
            rows.iter()
                .map(|row| row.iter().map(|c| q(*c)).collect())
                .collect(),
        )
    }

    #[test]
    fn bivariate_gcd_and_exact_quotient() {
        // (n + x) (1 - x) and (n + x) (2 + n): gcd n + x up to a constant.
        let common = biv(&[&[0, 1], &[1]]);
        let a = common.mul(&biv(&[&[1, -1]]));
        let b = common.mul(&biv(&[&[2], &[1]]));
        let gcd = bivariate_gcd(&a, &b);
        let quotient = bivariate_exact_quotient(&gcd, &common).expect("a constant multiple");
        assert!(quotient.as_constant().is_some_and(|c| !c.is_zero()));
        assert_eq!(
            bivariate_exact_quotient(&a, &common),
            Some(biv(&[&[1, -1]]))
        );
        assert_eq!(bivariate_exact_quotient(&a, &biv(&[&[1, 1]])), None);
        // Coprime polynomials.
        let gcd = bivariate_gcd(&biv(&[&[1, 1]]), &biv(&[&[0], &[1]]));
        assert!(gcd.as_constant().is_some());
    }

    #[test]
    fn primitive_operator_removes_a_common_bivariate_factor() {
        // (n + x) (x D - n) -> x D - n.
        let common = biv(&[&[0, 1], &[1]]);
        let mut op = DOp::single(0, 1, common.mul(&biv(&[&[0, 1]])));
        op.add_term((0, 0), common.mul(&biv(&[&[], &[-1]])));
        let mut expected = DOp::single(0, 1, biv(&[&[0, 1]]));
        expected.add_term((0, 0), biv(&[&[], &[-1]]));
        let primitive = primitive_operator(&op);
        let lead = primitive.leading().unwrap();
        let normalized = expected.scale(&expected.leading().unwrap().recip());
        assert_eq!(primitive.scale(&lead.recip()), normalized);
    }

    #[test]
    fn invariant_annihilator_of_a_binomial_parent() {
        // P_(n+1) = (1 + x) P_n, P_n = (1 + x)^n: T = 1 + x.  The remainder
        // r = (1 + x)^2 D^2 - n (n - 1) (a left multiple of
        // (1 + x) D - n) yields the invariant annihilator (1 + x) D - n.
        let step = DOp::single(0, 0, biv(&[&[1, 1]]));
        let mut r = DOp::single(0, 2, biv(&[&[1, 2, 1]]));
        r.add_term((0, 0), biv(&[&[], &[1], &[-1]]));
        let invariant = find_invariant_annihilator(&r, &step, 10_000)
            .unwrap()
            .expect("an invariant annihilator");
        assert_eq!(invariant.annihilator.max_derivative(), 1);
        assert!(invariance_holds(&invariant, &step));
        let rows: Vec<Vec<Q>> = (0..6usize)
            .map(|n| {
                (0..=n)
                    .map(|k| crate::recurrence::implicit::algebra::binomial(n, k))
                    .collect()
            })
            .collect();
        for (n, row) in rows.iter().enumerate() {
            assert!(invariant
                .annihilator
                .apply_row(row, n)
                .iter()
                .all(Zero::is_zero));
        }
        // nu is free of n here: no restart bases.
        assert!(vanishing_bases(&invariant.step_multiplier, 0).is_empty());
        assert_eq!(
            vanishing_bases(&biv(&[&[-3, 0], &[1, 1]]), 0),
            Vec::<usize>::new()
        );
        assert_eq!(vanishing_bases(&biv(&[&[-3], &[1]]), 0), vec![3]);
    }
}
