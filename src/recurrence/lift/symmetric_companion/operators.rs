//! The operators of the coupled system as exact Ore-algebra elements.
//!
//! Every parameter is a polynomial in the index `n` (ascending coefficients,
//! canonical `[0]` for zero).  The operators are [`DOp`] values
//! `sum c_(i,d)(n,x) D^d S^i`; composition follows `S n = (n + 1) S` and
//! `D x = x D + 1`, so the identities checked here are identities of
//! operators, valid for every index and every pair of rows.

use crate::recurrence::implicit::algebra::{Biv, DOp};
use crate::recurrence::BigRational as Q;
use crate::vec_poly;
use num_traits::{One, Zero};

/// The parameters of one candidate: `a_n`, `c_n`, `epsilon_n` as
/// polynomials in `n`, and the schedule `d_n = d0 + (n - n0)`.
///
/// `delta_n = c_n - epsilon_n` and `b_n = a_n + c_n d_n + delta_n` are
/// derived, so the defining relations of the issue hold by construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct CompanionParameters {
    /// `a_n`.
    pub(super) a: Vec<Q>,
    /// `c_n`.
    pub(super) c: Vec<Q>,
    /// `epsilon_n`.
    pub(super) epsilon: Vec<Q>,
    /// First index `n0` of the domain.
    pub(super) first_index: usize,
    /// `d_(n0)`.
    pub(super) symmetry_degree_at_first_index: usize,
}

impl CompanionParameters {
    pub(super) fn new(
        a: Vec<Q>,
        c: Vec<Q>,
        epsilon: Vec<Q>,
        first_index: usize,
        symmetry_degree_at_first_index: usize,
    ) -> Self {
        Self {
            a: vec_poly::trim(&a),
            c: vec_poly::trim(&c),
            epsilon: vec_poly::trim(&epsilon),
            first_index,
            symmetry_degree_at_first_index,
        }
    }

    /// `d_n = d0 + (n - n0)` as a polynomial in `n`.
    pub(super) fn schedule(&self) -> Vec<Q> {
        let d0 = Q::from_integer(self.symmetry_degree_at_first_index.into());
        let n0 = Q::from_integer(self.first_index.into());
        vec_poly::trim(&[d0 - n0, Q::one()])
    }

    /// `delta_n = c_n - epsilon_n`.
    pub(super) fn delta(&self) -> Vec<Q> {
        vec_poly::sub(&self.c, &self.epsilon)
    }

    /// `b_n = a_n + c_n d_n + delta_n`.
    pub(super) fn b(&self) -> Vec<Q> {
        let cd = vec_poly::mul(&self.c, &self.schedule());
        vec_poly::add(&vec_poly::add(&self.a, &cd), &self.delta())
    }

    /// The operators of this candidate.
    pub(super) fn operators(&self) -> SystemOperators {
        SystemOperators::new(self)
    }
}

/// `p(n + 1)`.
pub(super) fn shift_index_polynomial(p: &[Q]) -> Vec<Q> {
    vec_poly::trim(&vec_poly::compose(p, &[Q::one(), Q::one()]))
}

/// `p(value)`.
pub(super) fn evaluate_at_index(p: &[Q], value: usize) -> Q {
    vec_poly::evaluate(p, &Q::from_integer(value.into()))
}

/// The polynomial `p(n) + q(n) x` as a bivariate coefficient.
fn index_and_x(constant: &[Q], linear: &[Q]) -> Biv {
    let len = constant.len().max(linear.len());
    Biv::new(
        (0..len)
            .map(|i| {
                vec![
                    constant.get(i).cloned().unwrap_or_else(Q::zero),
                    linear.get(i).cloned().unwrap_or_else(Q::zero),
                ]
            })
            .collect(),
    )
}

/// The polynomial `p(n)` (free of `x`) as a bivariate coefficient.
pub(super) fn index_only(p: &[Q]) -> Biv {
    index_and_x(p, &[])
}

/// Multiplication by `p(n) + q(n) x`.
fn multiplication(constant: &[Q], linear: &[Q]) -> DOp {
    DOp::single(0, 0, index_and_x(constant, linear))
}

/// Multiplication by `p(n)`.
pub(super) fn scalar(p: &[Q]) -> DOp {
    multiplication(p, &[])
}

/// Multiplication by `x`.
pub(super) fn x_operator() -> DOp {
    DOp::single(0, 0, Biv::from_y(&[Q::zero(), Q::one()]))
}

/// The shift `S`.
pub(super) fn shift_operator() -> DOp {
    DOp::single(1, 0, Biv::one())
}

/// `c(n) x (1 - x) D`.
fn derivation(c: &[Q]) -> DOp {
    let rows = c
        .iter()
        .map(|ci| vec![Q::zero(), ci.clone(), -ci.clone()])
        .collect();
    DOp::single(0, 1, Biv::new(rows))
}

/// `H(d, a, c) = [a + (a + c d) x] + c x (1 - x) D`.
pub(super) fn h_operator(d: &[Q], a: &[Q], c: &[Q]) -> DOp {
    let linear = vec_poly::add(a, &vec_poly::mul(c, d));
    multiplication(a, &linear).add(&derivation(c))
}

/// The operators of the coupled system of one candidate.
#[derive(Debug, Clone)]
pub(super) struct SystemOperators {
    /// `U_n = H(d_n, a_n, c_n)`.
    pub(super) u: DOp,
    /// `V_n = H(d_n - 1, a_n + delta_n, c_n)`.
    pub(super) v: DOp,
    /// `Vtilde_n = V_n - c_n (1 - x)`.
    pub(super) v_tilde: DOp,
    /// The eliminated operator `E` (shift 2, base `n`).
    pub(super) eliminated: DOp,
    /// The parent operator `L_n = (a_n + b_n x) + c_n x (1 - x) D`.
    pub(super) parent: DOp,
    /// `epsilon_n`.
    epsilon: Vec<Q>,
    /// `epsilon_(n+1)`.
    epsilon_next: Vec<Q>,
    /// `delta_n`.
    delta: Vec<Q>,
    /// `a_n`.
    a: Vec<Q>,
}

impl SystemOperators {
    fn new(parameters: &CompanionParameters) -> Self {
        let d = parameters.schedule();
        let a = parameters.a.clone();
        let c = parameters.c.clone();
        let epsilon = parameters.epsilon.clone();
        let delta = parameters.delta();
        let b = parameters.b();
        let d_minus_one = vec_poly::sub(&d, &[Q::one()]);
        let u = h_operator(&d, &a, &c);
        let v = h_operator(&d_minus_one, &vec_poly::add(&a, &delta), &c);
        // c (1 - x) as a multiplication operator.
        let c_one_minus_x = multiplication(&c, &vec_poly::neg(&c));
        let v_tilde = v.sub(&c_one_minus_x);
        let epsilon_next = shift_index_polynomial(&epsilon);
        let s = shift_operator();
        let x = x_operator();
        // E = eps S^2 - eps (S o U) - eps' (Vt o S) + eps' (Vt o U) - eps eps' delta x.
        let eps = scalar(&epsilon);
        let eps_next = scalar(&epsilon_next);
        let product = vec_poly::mul(&vec_poly::mul(&epsilon, &epsilon_next), &delta);
        let eliminated = eps
            .compose(&s.compose(&s))
            .sub(&eps.compose(&u.shifted(1)))
            .sub(&eps_next.compose(&v_tilde.compose(&s)))
            .add(&eps_next.compose(&v_tilde.compose(&u)))
            .sub(&scalar(&product).compose(&x));
        let parent = multiplication(&a, &b).add(&derivation(&c));
        Self {
            u,
            v,
            v_tilde,
            eliminated,
            parent,
            epsilon,
            epsilon_next,
            delta,
            a,
        }
    }

    /// The elimination identity for arbitrary pairs `(A, B)`:
    ///
    /// ```text
    /// E A = eps_n E_(n+1) + eps_n eps_(n+1) x F_n - eps_(n+1) Vt_n E_n,
    /// E_n = A_(n+1) - U_n A_n - eps_n x B_n,   F_n = B_(n+1) - delta_n A_n - V_n B_n,
    /// ```
    ///
    /// checked componentwise as two operator identities (the `A` component
    /// and the vanishing `B` component, which is `Vt_n o x = x o V_n`).
    pub(super) fn elimination_identity_holds(&self) -> bool {
        let s = shift_operator();
        let x = x_operator();
        let eps = scalar(&self.epsilon);
        let eps_next = scalar(&self.epsilon_next);
        let eps_eps_next = scalar(&vec_poly::mul(&self.epsilon, &self.epsilon_next));
        // Residual operators on the A and B components.
        let first_on_a = s.sub(&self.u);
        let first_on_b = scalar(&vec_poly::neg(&self.epsilon)).compose(&x);
        let second_on_a = scalar(&vec_poly::neg(&self.delta));
        let second_on_b = s.sub(&self.v);
        let combine = |first: &DOp, second: &DOp| {
            eps.compose(&s.compose(first))
                .add(&eps_eps_next.compose(&x.compose(second)))
                .sub(&eps_next.compose(&self.v_tilde.compose(first)))
        };
        combine(&first_on_a, &second_on_a) == self.eliminated
            && combine(&first_on_b, &second_on_b).is_zero()
    }

    /// The collapse identities `L_n = U_n + delta_n x` and
    /// `L_n o x = eps_n x + x o V_n`, so `P_n = A_n + x B_n` satisfies
    /// `P_(n+1) = L_n P_n` whenever both coupled equations hold.
    pub(super) fn collapse_identity_holds(&self) -> bool {
        let x = x_operator();
        let on_a = self.u.add(&scalar(&self.delta).compose(&x));
        let on_b = scalar(&self.epsilon).compose(&x).add(&x.compose(&self.v));
        self.parent == on_a && self.parent.compose(&x) == on_b
    }

    /// The constant-term projection of `E`: `Some(E0)` when every term with a
    /// derivative has a coefficient vanishing at `x = 0`, so that
    /// `(E A)_n(0) = (E0 z)_n` with `z_n = A_n(0)`.
    pub(super) fn constant_term_projection(&self) -> Option<DOp> {
        let zero = Q::zero();
        let mut projection = DOp::default();
        for ((shift, derivative), coefficient) in &self.eliminated.0 {
            let at_zero = coefficient.eval_y(&zero);
            if *derivative > 0 {
                if !vec_poly::is_zero(&at_zero) {
                    return None;
                }
            } else {
                projection.add_term((*shift, 0), index_only(&at_zero));
            }
        }
        Some(projection)
    }

    /// The factorization `E0 = (eps_n S - g_n) o (S - a_n)` with
    /// `g_n = eps_(n+1) (a_n - eps_n)`, so `w_n = A_(n+1)(0) - a_n A_n(0)`
    /// satisfies `eps_n w_(n+1) = g_n w_n`.
    pub(super) fn constant_term_factorization_holds(&self, projection: &DOp) -> bool {
        let s = shift_operator();
        let g = vec_poly::mul(&self.epsilon_next, &vec_poly::sub(&self.a, &self.epsilon));
        let left = scalar(&self.epsilon).compose(&s).sub(&scalar(&g));
        let right = s.sub(&scalar(&self.a));
        left.compose(&right) == *projection
    }
}

/// Apply a shift-free operator to one row at index `n`.
pub(super) fn apply_to_row(operator: &DOp, row: &[Q], n: usize) -> Vec<Q> {
    vec_poly::trim(&operator.apply_row(row, n))
}

/// Whether `p` is symmetric with respect to `d`: `deg p <= d` and
/// `x^d p(1/x) = p` (the zero polynomial is symmetric at every `d`).
pub(super) fn is_symmetric_at(p: &[Q], d: usize) -> bool {
    vec_poly::reverse(p, d).is_some_and(|r| vec_poly::trim(&r) == vec_poly::trim(p))
}

/// `R / x` when `R(0) = 0`.
pub(super) fn divide_by_x(p: &[Q]) -> Option<Vec<Q>> {
    let p = vec_poly::trim(p);
    if !p[0].is_zero() {
        return None;
    }
    if p.len() == 1 {
        return Some(vec![Q::zero()]);
    }
    Some(vec_poly::trim(&p[1..]))
}
