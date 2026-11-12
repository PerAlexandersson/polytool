//! Bivariate polynomial operations and transformation gauges.

use super::*;

// ---------------------------------------------------------------------------
// Bivariate helpers: `[n power][x power]`, zero is the empty vector
// ---------------------------------------------------------------------------

pub(in crate::recurrence::transform) type Bp = Vec<Vec<Q>>;

pub(in crate::recurrence::transform) fn bp_trim(mut p: Bp) -> Bp {
    for row in &mut p {
        vec_poly::strip_trailing_zeros(row);
    }
    while p.last().is_some_and(Vec::is_empty) {
        p.pop();
    }
    p
}

pub(in crate::recurrence::transform) fn bp_from(b: &BivarPoly) -> Bp {
    bp_trim(b.coeffs.clone())
}

pub(in crate::recurrence::transform) fn bp_to_bivar(p: Bp) -> BivarPoly {
    BivarPoly::from_coefficients_normalized(bp_trim(p))
}

pub(in crate::recurrence::transform) fn bp_is_zero(p: &Bp) -> bool {
    p.iter().all(|row| vec_poly::is_zero(row))
}

pub(in crate::recurrence::transform) fn bp_constant(c: Q) -> Bp {
    bp_trim(vec![vec![c]])
}

fn bp_add_assign(acc: &mut Bp, p: &Bp) {
    if acc.len() < p.len() {
        acc.resize(p.len(), Vec::new());
    }
    for (target, row) in acc.iter_mut().zip(p) {
        *target = vec_poly::add(target, row);
    }
}

pub(in crate::recurrence::transform) fn bp_add(a: &Bp, b: &Bp) -> Bp {
    let mut sum = a.clone();
    bp_add_assign(&mut sum, b);
    bp_trim(sum)
}

pub(in crate::recurrence::transform) fn bp_scale(p: &Bp, s: &Q) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::scale(row, s)).collect())
}

pub(in crate::recurrence::transform) fn bp_mul(a: &Bp, b: &Bp) -> Bp {
    if bp_is_zero(a) || bp_is_zero(b) {
        return Vec::new();
    }
    let mut result: Bp = vec![Vec::new(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate() {
        for (j, right) in b.iter().enumerate() {
            result[i + j] = vec_poly::add(&result[i + j], &vec_poly::mul(left, right));
        }
    }
    bp_trim(result)
}

pub(in crate::recurrence::transform) fn bp_mul_x(p: &Bp, poly: &[Q]) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::mul(row, poly)).collect())
}

pub(in crate::recurrence::transform) fn bp_dx(p: &Bp) -> Bp {
    bp_trim(p.iter().map(|row| vec_poly::derivative(row)).collect())
}

/// `c0 + c1 n`.
pub(in crate::recurrence::transform) fn bp_affine_index(c0: Q, c1: Q) -> Bp {
    bp_trim(vec![vec![c0], vec![c1]])
}

pub(in crate::recurrence::transform) fn bp_x_degree(p: &Bp) -> Option<usize> {
    p.iter().filter_map(|row| vec_poly::degree(row)).max()
}

pub(in crate::recurrence::transform) fn bp_div_linear(p: &Bp, lambda: &[Q]) -> Option<Bp> {
    let root = -(&lambda[0] / &lambda[1]);
    let inverse_lead = lambda[1].recip();
    p.iter()
        .map(|row| {
            vec_poly::exact_div_by_t_minus(row, &root)
                .map(|quotient| vec_poly::scale(&quotient, &inverse_lead))
        })
        .collect::<Option<Bp>>()
        .map(bp_trim)
}

pub(super) fn bp_is_constant(p: &Bp) -> Option<Q> {
    if bp_is_zero(p) {
        return None;
    }
    let tail_zero = p.iter().enumerate().all(|(i, row)| {
        row.iter()
            .enumerate()
            .all(|(j, c)| (i == 0 && j == 0) || c.is_zero())
    });
    tail_zero.then(|| p[0][0].clone())
}

/// `sum c_(i,j) n^i num^j den^(D - j)` with `D` the `x`-degree of `c`.
pub(in crate::recurrence::transform) fn bp_homogenize(
    p: &Bp,
    numerator: &[Q],
    denominator: &[Q],
) -> (Bp, usize) {
    let degree = bp_x_degree(p).unwrap_or(0);
    let numerator_powers = (0..=degree)
        .map(|j| vec_poly::pow(numerator, j))
        .collect::<Vec<_>>();
    let denominator_powers = (0..=degree)
        .map(|j| vec_poly::pow(denominator, j))
        .collect::<Vec<_>>();
    let result = p
        .iter()
        .map(|row| {
            let mut sum = vec![Q::zero()];
            for (j, coefficient) in row.iter().enumerate() {
                if coefficient.is_zero() {
                    continue;
                }
                let basis = vec_poly::mul(&numerator_powers[j], &denominator_powers[degree - j]);
                sum = vec_poly::add(&sum, &vec_poly::scale(&basis, coefficient));
            }
            sum
        })
        .collect();
    (bp_trim(result), degree)
}

// ---------------------------------------------------------------------------
// Scalar transport
// ---------------------------------------------------------------------------

/// A scalar gauge `Q_n(x) = lambda(x)^(w + h n) P_n(phi(x))` with
/// `phi = num/den`: the linear form `lambda`, the homogenization and the
/// exponents used by the scalar transport (additive terms, clearing and
/// normalization), and the chain rule
/// `P^(k)(phi) = kappa^k lambda^(alpha_k) T_k[Q]`,
/// `alpha_0 = -(w + h n)`, `alpha_(k+1) = alpha_k + sigma`,
/// `T_(k+1) = lambda T_k' + alpha_k lambda' T_k` of the componentwise
/// vector transport (module `transform::vector_transport`).
pub(in crate::recurrence::transform) struct Gauge {
    pub(in crate::recurrence::transform) numerator: Vec<Q>,
    pub(in crate::recurrence::transform) denominator: Vec<Q>,
    /// Exponent of `lambda` contributed by `den^D`: 1 for Möbius (`den = lambda`), 0 otherwise.
    pub(in crate::recurrence::transform) tau: i64,
    pub(in crate::recurrence::transform) lambda: [Q; 2],
    pub(in crate::recurrence::transform) weight_constant: i64,
    pub(in crate::recurrence::transform) step: i64,
    pub(in crate::recurrence::transform) sigma: i64,
    pub(in crate::recurrence::transform) kappa: Q,
}

impl Gauge {
    pub(in crate::recurrence::transform) fn monomial(power: i64) -> Self {
        Self {
            numerator: vec![q(0), q(1)],
            denominator: vec![q(1)],
            tau: 0,
            lambda: [q(0), q(1)],
            weight_constant: power,
            step: 0,
            sigma: -1,
            kappa: q(1),
        }
    }

    pub(in crate::recurrence::transform) fn mobius(
        a: Q,
        b: Q,
        c: Q,
        e: Q,
        weight_constant: i64,
        step: i64,
    ) -> Self {
        let delta = &a * &e - &b * &c;
        Self {
            numerator: vec_poly::trim(&[b, a]),
            denominator: vec_poly::trim(&[e.clone(), c.clone()]),
            tau: 1,
            lambda: [e, c],
            weight_constant,
            step,
            sigma: 1,
            kappa: delta.recip(),
        }
    }

    pub(in crate::recurrence::transform) fn lambda_is_constant(&self) -> bool {
        self.lambda[1].is_zero()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::recurrence::transform) enum PieceKey {
    Lhs,
    Term {
        offset: usize,
        derivative: usize,
        alternating: bool,
    },
    Inhomogeneous,
}

pub(in crate::recurrence::transform) fn checked(value: Option<i64>) -> Result<i64, TransformError> {
    value.ok_or_else(|| TransformError::Limit("exponent overflow".into()))
}
