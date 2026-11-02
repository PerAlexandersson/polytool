//! Internal exact algebra behind implicit differential relations: an
//! adapter over the `ore-algebra` crate.
//!
//! - [`Biv`] (`ore_algebra::BivariatePolynomial`): bivariate polynomials
//!   `sum c_(i,j) n^i y^j` over the rationals, where `y` is `x`, the Euler
//!   operator `theta = x D` or a coefficient index `k` depending on context;
//! - [`DOp`] (`ore_algebra::OreOperator`): normal-ordered operators
//!   `sum c_(i,d)(n,x) D^d S^i` with composition in the Ore algebra
//!   (`S n = (n + 1) S`, `D x = x D + 1`);
//! - the Euler form `sum x^s a_(i,s)(n,theta) S^i` and the conversions
//!   `x^J L <-> Euler form`;
//! - exact zero sets of index polynomials `p(n,k)` on `n >= n_0`, `k >= 0`
//!   (these use the Sturm chains of this crate and stay here).
//!
//! The conversions to and from [`BivarPoly`] are the [`BivarConversion`]
//! extension trait (in scope in every recurrence submodule through the
//! facade).  Nothing here is public outside the recurrence module.

use super::{BigRational as Q, BivarPoly};
use crate::polynomial::Polynomial;
use crate::sturm::SturmChain;
use crate::vec_poly;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeSet;

pub(in crate::recurrence) use ore_algebra::{
    add_euler_term as theta_add, falling_factorial, from_euler_form as from_theta,
    to_euler_form as to_theta, BivariatePolynomial as Biv, EulerOperator as ThetaOp,
    OreOperator as DOp,
};

pub(in crate::recurrence) fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

pub(in crate::recurrence) fn q_usize(value: usize) -> Q {
    Q::from_integer(BigInt::from(value))
}

/// `binom(a, t)` as an exact rational.
pub(in crate::recurrence) fn binomial(a: usize, t: usize) -> Q {
    if t > a {
        return Q::zero();
    }
    let mut value = BigInt::one();
    for i in 0..t {
        value *= BigInt::from(a - i);
        value /= BigInt::from(i + 1);
    }
    Q::from_integer(value)
}

/// `t!` as an exact rational.
pub(in crate::recurrence) fn factorial(t: usize) -> Q {
    let mut value = BigInt::one();
    for i in 2..=t {
        value *= BigInt::from(i);
    }
    Q::from_integer(value)
}

/// The conversions of [`Biv`] to and from the public coefficient type
/// [`BivarPoly`] (an extension trait, since `Biv` is defined in the
/// `ore-algebra` crate).
pub(in crate::recurrence) trait BivarConversion: Sized {
    /// From the public coefficient type.
    fn from_bivar(p: &BivarPoly) -> Self;
    /// To the public coefficient type.
    fn to_bivar(&self) -> BivarPoly;
}

impl BivarConversion for Biv {
    fn from_bivar(p: &BivarPoly) -> Self {
        Self::new(p.coeffs.clone())
    }
    fn to_bivar(&self) -> BivarPoly {
        BivarPoly::from_coefficients_normalized(self.rows().to_vec())
    }
}

// ---------------------------------------------------------------------------
// Zero sets of index polynomials
// ---------------------------------------------------------------------------

/// Nonnegative integer roots of a nonzero univariate rational polynomial.
pub(in crate::recurrence) fn nonnegative_integer_roots(p: &[Q]) -> BTreeSet<BigInt> {
    let p = vec_poly::trim(p);
    let mut roots = BTreeSet::new();
    let Some(degree) = vec_poly::degree(&p) else {
        return roots;
    };
    let valuation = p.iter().position(|c| !c.is_zero()).unwrap_or(0);
    if valuation > 0 {
        roots.insert(BigInt::zero());
    }
    if degree == valuation {
        return roots;
    }
    let rest = &p[valuation..=degree];
    let common = rest.iter().fold(BigInt::one(), |acc, c| acc.lcm(c.denom()));
    let integers = rest
        .iter()
        .map(|c| (c * Q::from_integer(common.clone())).to_integer())
        .collect::<Vec<_>>();
    let chain = SturmChain::from_bigint_coeffs(&integers);
    let half = Q::new(BigInt::one(), BigInt::from(2));
    for (lo, hi) in chain.isolate_roots(&half) {
        let first = lo.numer().div_ceil(lo.denom());
        let last = hi.numer().div_floor(hi.denom());
        let mut candidate = first;
        while candidate <= last {
            if !candidate.is_negative()
                && vec_poly::evaluate(rest, &Q::from_integer(candidate.clone())).is_zero()
            {
                roots.insert(candidate.clone());
            }
            candidate += 1;
        }
    }
    roots
}

/// Nonnegative integer roots of `p(n)` inside the index domain.
pub(in crate::recurrence) fn roots_in_domain(
    p: &[Q],
    first: usize,
    last: Option<usize>,
) -> Vec<usize> {
    nonnegative_integer_roots(p)
        .into_iter()
        .filter_map(|r| usize::try_from(r).ok())
        .filter(|r| *r >= first && last.is_none_or(|l| *r <= l))
        .collect()
}

/// Exact description of the zeros of `p(n,k)` with `n` in the domain and `k >= 0`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(in crate::recurrence) struct ZeroAnalysis {
    /// Indices `n` of the domain where `p(n, .)` vanishes identically.
    pub(in crate::recurrence) singular_indices: Vec<usize>,
    /// Indices `k >= 0` with `p(n, k) = 0` for every `n` of the domain.
    pub(in crate::recurrence) uniform_indices: Vec<usize>,
    /// Whether the zero set is exactly `singular x N  union  domain x uniform`.
    pub(in crate::recurrence) certified: bool,
    /// How the description was obtained.
    pub(in crate::recurrence) certificate: String,
}

/// Univariate polynomial gcd over the rationals (monic, zero for all-zero input).
fn gcd_all(polys: &[Vec<Q>]) -> Vec<Q> {
    let mut gcd = Polynomial::<Q>::zero();
    for p in polys {
        gcd = gcd.gcd(&Polynomial::new(p.clone()));
    }
    vec_poly::trim(gcd.coeffs())
}

/// The content `c(n) = gcd_j [k^j] p` and primitive part `p / c(n)`.
pub(in crate::recurrence) fn content_and_primitive(p: &Biv) -> (Vec<Q>, Biv) {
    let width = p.y_degree().map_or(0, |d| d + 1);
    let columns = (0..width).map(|j| p.y_coefficient(j)).collect::<Vec<_>>();
    let content = gcd_all(&columns);
    if vec_poly::degree(&content).unwrap_or(0) == 0 {
        return (vec![Q::one()], p.clone());
    }
    let divisor = Polynomial::new(content.clone());
    let quotient_columns = columns
        .iter()
        .map(|column| {
            let (quotient, _) = Polynomial::new(column.clone()).div_rem(&divisor);
            vec_poly::trim(quotient.coeffs())
        })
        .collect::<Vec<_>>();
    // Reassemble [n power][k power] from the k-columns.
    let n_len = quotient_columns.iter().map(Vec::len).max().unwrap_or(0);
    let rows = (0..n_len)
        .map(|i| {
            quotient_columns
                .iter()
                .map(|column| column.get(i).cloned().unwrap_or_else(Q::zero))
                .collect()
        })
        .collect();
    (content, Biv::new(rows))
}

pub(in crate::recurrence) fn poly_text(p: &[Q], variable: &str) -> String {
    let p = vec_poly::trim(p);
    let mut parts = Vec::new();
    for (i, c) in p.iter().enumerate() {
        if c.is_zero() {
            continue;
        }
        let monomial = match i {
            0 => String::new(),
            1 => variable.to_string(),
            _ => format!("{variable}^{i}"),
        };
        parts.push(match (monomial.is_empty(), c.is_one(), (-c).is_one()) {
            (true, _, _) => c.to_string(),
            (false, true, _) => monomial,
            (false, _, true) => format!("-{monomial}"),
            _ => format!("{c}{monomial}"),
        });
    }
    if parts.is_empty() {
        return "0".into();
    }
    let mut text = parts[0].clone();
    for part in &parts[1..] {
        match part.strip_prefix('-') {
            Some(rest) => {
                text.push_str(" - ");
                text.push_str(rest);
            }
            None => {
                text.push_str(" + ");
                text.push_str(part);
            }
        }
    }
    text
}

/// Text of a polynomial in `(n, y)`.
pub(in crate::recurrence) fn biv_text(p: &Biv, variable: &str) -> String {
    let mut parts = Vec::new();
    let width = p.y_degree().map_or(0, |d| d + 1);
    for j in 0..width {
        for (i, row) in p.rows().iter().enumerate() {
            let Some(c) = row.get(j) else { continue };
            if c.is_zero() {
                continue;
            }
            let mut monomial = String::new();
            if i > 0 {
                monomial.push_str(&if i == 1 { "n".into() } else { format!("n^{i}") });
            }
            if j > 0 {
                if !monomial.is_empty() {
                    monomial.push(' ');
                }
                monomial.push_str(&if j == 1 {
                    variable.to_string()
                } else {
                    format!("{variable}^{j}")
                });
            }
            parts.push(match (monomial.is_empty(), c.is_one(), (-c).is_one()) {
                (true, _, _) => c.to_string(),
                (false, true, _) => monomial,
                (false, _, true) => format!("-{monomial}"),
                _ => format!("{c} {monomial}"),
            });
        }
    }
    if parts.is_empty() {
        return "0".into();
    }
    let mut text = parts[0].clone();
    for part in &parts[1..] {
        match part.strip_prefix('-') {
            Some(rest) => {
                text.push_str(" - ");
                text.push_str(rest);
            }
            None => {
                text.push_str(" + ");
                text.push_str(part);
            }
        }
    }
    text
}

/// Zeros of `p(n, k)` for `n` in `first..(=last)` and integers `k >= 0`.
///
/// The content `c(n)` (gcd of the `k`-coefficients) gives the indices where
/// `p(n, .)` vanishes identically. For the primitive part `psi`, an index-free
/// `psi(k)` has the uniform zero set of its nonnegative integer roots, and
/// `psi(first + u, k)` with coefficients of one strict sign has no zeros.
/// Anything else is reported as uncertified; no prefix testing is used.
pub(in crate::recurrence) fn analyze_zeros(
    p: &Biv,
    first: usize,
    last: Option<usize>,
) -> ZeroAnalysis {
    if p.is_zero() {
        return ZeroAnalysis {
            singular_indices: Vec::new(),
            uniform_indices: Vec::new(),
            certified: false,
            certificate: "the polynomial is identically zero".into(),
        };
    }
    let (content, primitive) = content_and_primitive(p);
    let singular_indices = roots_in_domain(&content, first, last);
    let content_text = if vec_poly::degree(&content).unwrap_or(0) == 0 {
        String::new()
    } else {
        format!(
            "content c(n) = {} vanishes for domain indices {:?}; ",
            poly_text(&content, "n"),
            singular_indices
        )
    };
    if primitive.is_n_free() {
        let psi = primitive.rows().first().cloned().unwrap_or_default();
        let uniform = nonnegative_integer_roots(&psi)
            .into_iter()
            .filter_map(|r| usize::try_from(r).ok())
            .collect::<Vec<_>>();
        return ZeroAnalysis {
            certificate: format!(
                "{content_text}the index-free part psi(k) = {} has nonnegative integer roots \
                 {uniform:?} (exact Sturm isolation)",
                poly_text(&psi, "k")
            ),
            singular_indices,
            uniform_indices: uniform,
            certified: true,
        };
    }
    let shifted = primitive.shift(&q_usize(first), &Q::zero());
    let mut signs = shifted
        .rows()
        .iter()
        .flatten()
        .filter(|c| !c.is_zero())
        .map(|c| c.is_positive());
    let first_sign = signs.next();
    let one_signed = first_sign.is_some() && signs.all(|s| Some(s) == first_sign);
    let constant = shifted
        .rows()
        .first()
        .and_then(|row| row.first())
        .is_some_and(|c| !c.is_zero());
    if one_signed && constant {
        return ZeroAnalysis {
            certificate: format!(
                "{content_text}psi(n_0 + u, k) with n_0 = {first} has coefficients of one strict \
                 sign and a nonzero constant term, so psi has no zeros with n >= {first}, k >= 0"
            ),
            singular_indices,
            uniform_indices: Vec::new(),
            certified: true,
        };
    }
    ZeroAnalysis {
        certificate: format!(
            "{content_text}psi(n, k) = {} depends on n and has no one-signed certificate; \
             its zero set is not certified uniformly",
            biv_text(&primitive, "k")
        ),
        singular_indices,
        uniform_indices: Vec::new(),
        certified: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn theta_round_trip_and_composition() {
        // L = x D^2 + (1 + x) D - n.
        let mut op = DOp::default();
        op.add_term((0, 2), Biv::from_y(&[q(0), q(1)]));
        op.add_term((0, 1), Biv::from_y(&[q(1), q(1)]));
        op.add_term((0, 0), Biv::from_n(&[q(0), q(-1)]));
        let (j, theta) = to_theta(&op);
        assert_eq!(j, 2);
        assert_eq!(from_theta(&theta), op.mul_x_power(2));
        // D o x = x D + 1.
        let d = DOp::single(0, 1, Biv::one());
        let x = DOp::single(0, 0, Biv::from_y(&[q(0), q(1)]));
        let mut expected = DOp::single(0, 1, Biv::from_y(&[q(0), q(1)]));
        expected.add_term((0, 0), Biv::one());
        assert_eq!(d.compose(&x), expected);
        // S o n = (n + 1) S.
        let s = DOp::single(1, 0, Biv::one());
        let n = DOp::single(0, 0, Biv::from_n(&[q(0), q(1)]));
        assert_eq!(s.compose(&n), DOp::single(1, 0, Biv::from_n(&[q(1), q(1)])));
    }

    #[test]
    fn zero_sets() {
        // k (k - 3): uniform zeros {0, 3}.
        let p = Biv::from_y(&[q(0), q(-3), q(1)]);
        let a = analyze_zeros(&p, 0, None);
        assert!(a.certified);
        assert_eq!(a.uniform_indices, vec![0, 3]);
        // (n - 2)(k + 1): singular at n = 2, no uniform zeros.
        let p = Biv::from_n(&[q(-2), q(1)]).mul(&Biv::from_y(&[q(1), q(1)]));
        let a = analyze_zeros(&p, 0, None);
        assert_eq!(a.singular_indices, vec![2]);
        assert!(a.uniform_indices.is_empty() && a.certified);
        // n + k + 1 has no zeros on the quadrant.
        let p = Biv::new(vec![vec![q(1), q(1)], vec![q(1)]]);
        assert!(analyze_zeros(&p, 0, None).certified);
        // k - n is not uniformly certified.
        let p = Biv::new(vec![vec![q(0), q(1)], vec![q(-1)]]);
        assert!(!analyze_zeros(&p, 0, None).certified);
        assert_eq!(
            nonnegative_integer_roots(&[q(-6), q(11), q(-6), q(1)])
                .into_iter()
                .collect::<Vec<_>>(),
            vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]
        );
    }
}
