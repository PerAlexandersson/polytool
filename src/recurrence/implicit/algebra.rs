//! Internal exact algebra behind implicit differential relations.
//!
//! - [`Biv`]: bivariate polynomials `sum c_(i,j) n^i y^j` over the rationals,
//!   where `y` is `x`, the Euler operator `theta = x D` or a coefficient
//!   index `k` depending on context;
//! - [`DOp`]: normal-ordered operators `sum c_(i,d)(n,x) D^d S^i` with
//!   composition in the Ore algebra (`S n = (n + 1) S`, `D x = x D + 1`);
//! - the Euler form `sum x^s a_(i,s)(n,theta) S^i` and the conversions
//!   `x^J L <-> Euler form`;
//! - exact zero sets of index polynomials `p(n,k)` on `n >= n_0`, `k >= 0`.
//!
//! Everything is exact; nothing here is public outside the recurrence module.

use super::{BigRational as Q, BivarPoly};
use crate::polynomial::Polynomial;
use crate::sturm::SturmChain;
use crate::vec_poly;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use std::collections::{BTreeMap, BTreeSet};

pub(super) fn q(value: i64) -> Q {
    Q::from_integer(BigInt::from(value))
}

pub(super) fn q_usize(value: usize) -> Q {
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

// ---------------------------------------------------------------------------
// Bivariate polynomials
// ---------------------------------------------------------------------------

/// A bivariate polynomial `rows[i][j]` = coefficient of `n^i y^j`, kept
/// trimmed (no trailing zeros, zero is the empty vector) so that equality is
/// structural.
#[derive(Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub(in crate::recurrence) struct Biv(Vec<Vec<Q>>);

impl Biv {
    pub(super) fn new(mut rows: Vec<Vec<Q>>) -> Self {
        for row in &mut rows {
            vec_poly::strip_trailing_zeros(row);
        }
        while rows.last().is_some_and(Vec::is_empty) {
            rows.pop();
        }
        Self(rows)
    }

    pub(super) fn rows(&self) -> &[Vec<Q>] {
        &self.0
    }

    pub(super) fn zero() -> Self {
        Self(Vec::new())
    }

    pub(super) fn constant(c: Q) -> Self {
        Self::new(vec![vec![c]])
    }

    pub(super) fn one() -> Self {
        Self::constant(Q::one())
    }

    /// `c n^n_power y^y_power`.
    pub(super) fn monomial(n_power: usize, y_power: usize, c: Q) -> Self {
        let mut rows = vec![Vec::new(); n_power + 1];
        let mut row = vec![Q::zero(); y_power + 1];
        row[y_power] = c;
        rows[n_power] = row;
        Self::new(rows)
    }

    /// A polynomial in `y` only.
    pub(super) fn from_y(p: &[Q]) -> Self {
        Self::new(vec![p.to_vec()])
    }

    /// A polynomial in `n` only.
    #[cfg(test)]
    pub(super) fn from_n(p: &[Q]) -> Self {
        Self::new(p.iter().map(|c| vec![c.clone()]).collect())
    }

    pub(super) fn from_bivar(p: &BivarPoly) -> Self {
        Self::new(p.coeffs.clone())
    }

    pub(super) fn to_bivar(&self) -> BivarPoly {
        BivarPoly::from_coefficients_normalized(self.0.clone())
    }

    pub(super) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// The constant value, when the polynomial is a constant (zero included).
    pub(super) fn as_constant(&self) -> Option<Q> {
        match self.0.as_slice() {
            [] => Some(Q::zero()),
            [row] if row.len() == 1 => Some(row[0].clone()),
            _ => None,
        }
    }

    pub(super) fn is_n_free(&self) -> bool {
        self.0.len() <= 1
    }

    pub(super) fn is_y_free(&self) -> bool {
        self.0.iter().all(|row| row.len() <= 1)
    }

    pub(super) fn add(&self, other: &Self) -> Self {
        let len = self.0.len().max(other.0.len());
        Self::new(
            (0..len)
                .map(|i| {
                    let a = self.0.get(i).map_or(&[][..], Vec::as_slice);
                    let b = other.0.get(i).map_or(&[][..], Vec::as_slice);
                    vec_poly::add(a, b)
                })
                .collect(),
        )
    }

    pub(super) fn scale(&self, c: &Q) -> Self {
        if c.is_zero() {
            return Self::zero();
        }
        Self::new(
            self.0
                .iter()
                .map(|row| row.iter().map(|x| x * c).collect())
                .collect(),
        )
    }

    pub(super) fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut rows = vec![Vec::new(); self.0.len() + other.0.len() - 1];
        for (i, left) in self.0.iter().enumerate() {
            if left.is_empty() {
                continue;
            }
            for (j, right) in other.0.iter().enumerate() {
                if right.is_empty() {
                    continue;
                }
                rows[i + j] = vec_poly::add(&rows[i + j], &vec_poly::mul(left, right));
            }
        }
        Self::new(rows)
    }

    pub(super) fn pow(&self, exponent: usize) -> Self {
        (0..exponent).fold(Self::one(), |acc, _| acc.mul(self))
    }

    /// `p(n + dn, y + dy)`.
    pub(in crate::recurrence) fn shift(&self, dn: &Q, dy: &Q) -> Self {
        let n_linear = [dn.clone(), Q::one()];
        let mut n_power = vec![Q::one()];
        let mut result = Self::zero();
        for (i, row) in self.0.iter().enumerate() {
            if i > 0 {
                n_power = vec_poly::mul(&n_power, &n_linear);
            }
            if row.is_empty() {
                continue;
            }
            let shifted = if dy.is_zero() {
                row.clone()
            } else {
                vec_poly::compose(row, &[dy.clone(), Q::one()])
            };
            let term = Self::new(
                n_power
                    .iter()
                    .map(|c| vec_poly::scale(&shifted, c))
                    .collect(),
            );
            result = result.add(&term);
        }
        result
    }

    pub(super) fn shift_n(&self, dn: i64) -> Self {
        self.shift(&q(dn), &Q::zero())
    }

    /// `d/dy`.
    pub(super) fn dy(&self) -> Self {
        Self::new(self.0.iter().map(|row| vec_poly::derivative(row)).collect())
    }

    /// The polynomial in `y` obtained at `n = value`.
    pub(super) fn eval_n(&self, value: &Q) -> Vec<Q> {
        let mut result = vec![Q::zero()];
        for row in self.0.iter().rev() {
            result = vec_poly::scale(&result, value);
            result = vec_poly::add(&result, row);
        }
        vec_poly::trim(&result)
    }

    /// The polynomial in `n` obtained at `y = value`.
    pub(in crate::recurrence) fn eval_y(&self, value: &Q) -> Vec<Q> {
        vec_poly::trim(
            &self
                .0
                .iter()
                .map(|row| vec_poly::evaluate(row, value))
                .collect::<Vec<_>>(),
        )
    }

    /// Coefficient of `y^j` as a polynomial in `n`.
    pub(super) fn y_coefficient(&self, j: usize) -> Vec<Q> {
        vec_poly::trim(
            &self
                .0
                .iter()
                .map(|row| row.get(j).cloned().unwrap_or_else(Q::zero))
                .collect::<Vec<_>>(),
        )
    }

    pub(super) fn y_degree(&self) -> Option<usize> {
        self.0.iter().filter_map(|row| vec_poly::degree(row)).max()
    }

    /// Smallest power of `y` with a nonzero coefficient.
    pub(super) fn y_valuation(&self) -> Option<usize> {
        self.0
            .iter()
            .filter_map(|row| row.iter().position(|c| !c.is_zero()))
            .min()
    }

    /// Divide by `y^e` (the caller guarantees divisibility).
    pub(super) fn div_y_power(&self, e: usize) -> Self {
        Self::new(
            self.0
                .iter()
                .map(|row| row.get(e..).map(<[Q]>::to_vec).unwrap_or_default())
                .collect(),
        )
    }

    /// Multiply by `y^e`.
    #[cfg(test)]
    pub(super) fn mul_y_power(&self, e: usize) -> Self {
        Self::new(
            self.0
                .iter()
                .map(|row| {
                    if row.is_empty() {
                        Vec::new()
                    } else {
                        vec_poly::shift_by(row, e)
                    }
                })
                .collect(),
        )
    }

    /// Lexicographically leading coefficient: highest power of `n`, then of `y`.
    pub(super) fn leading(&self) -> Option<&Q> {
        self.0.last().and_then(|row| row.last())
    }

    /// `(lead, p / lead)`; the zero polynomial is returned unchanged with lead 1.
    pub(super) fn monic(&self) -> (Q, Self) {
        match self.leading() {
            Some(lead) => {
                let lead = lead.clone();
                (lead.clone(), self.scale(&lead.recip()))
            }
            None => (Q::one(), Self::zero()),
        }
    }
}

// ---------------------------------------------------------------------------
// Normal-ordered operators
// ---------------------------------------------------------------------------

/// `sum c_(i,d)(n,x) D^d S^i`, keyed by `(shift i, derivative order d)`;
/// zero coefficients are never stored.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct DOp(pub(super) BTreeMap<(usize, usize), Biv>);

impl DOp {
    pub(super) fn add_term(&mut self, key: (usize, usize), coefficient: Biv) {
        if coefficient.is_zero() {
            return;
        }
        let entry = self.0.entry(key).or_default();
        *entry = entry.add(&coefficient);
        if entry.is_zero() {
            self.0.remove(&key);
        }
    }

    pub(super) fn single(shift: usize, derivative: usize, coefficient: Biv) -> Self {
        let mut op = Self::default();
        op.add_term((shift, derivative), coefficient);
        op
    }

    pub(super) fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub(super) fn add(&self, other: &Self) -> Self {
        let mut sum = self.clone();
        for (key, c) in &other.0 {
            sum.add_term(*key, c.clone());
        }
        sum
    }

    pub(super) fn scale(&self, c: &Q) -> Self {
        let mut out = Self::default();
        for (key, coefficient) in &self.0 {
            out.add_term(*key, coefficient.scale(c));
        }
        out
    }

    pub(super) fn sub(&self, other: &Self) -> Self {
        self.add(&other.scale(&q(-1)))
    }

    pub(super) fn max_shift(&self) -> usize {
        self.0.keys().map(|(i, _)| *i).max().unwrap_or(0)
    }

    pub(super) fn max_derivative(&self) -> usize {
        self.0.keys().map(|(_, d)| *d).max().unwrap_or(0)
    }

    /// The terms at shift `i`, as a shift-zero operator.
    pub(super) fn column(&self, shift: usize) -> Self {
        let mut out = Self::default();
        for ((i, d), c) in &self.0 {
            if *i == shift {
                out.add_term((0, *d), c.clone());
            }
        }
        out
    }

    /// The composition `self o other`.
    pub(super) fn compose(&self, other: &Self) -> Self {
        let mut out = Self::default();
        for ((i, a), c) in &self.0 {
            let dn = q_usize(*i);
            for ((j, b), e) in &other.0 {
                let mut derivative = e.shift(&dn, &Q::zero());
                for t in 0..=*a {
                    if derivative.is_zero() {
                        break;
                    }
                    let coefficient = c.mul(&derivative).scale(&binomial(*a, t));
                    out.add_term((i + j, a - t + b), coefficient);
                    derivative = derivative.dy();
                }
            }
        }
        out
    }

    /// `S^nu o self`.
    pub(super) fn shifted(&self, nu: usize) -> Self {
        let mut out = Self::default();
        let dn = q_usize(nu);
        for ((i, d), c) in &self.0 {
            out.add_term((i + nu, *d), c.shift(&dn, &Q::zero()));
        }
        out
    }

    /// Smallest power of `x` dividing every coefficient.
    pub(super) fn x_valuation(&self) -> usize {
        self.0
            .values()
            .filter_map(Biv::y_valuation)
            .min()
            .unwrap_or(0)
    }

    pub(super) fn div_x_power(&self, e: usize) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.div_y_power(e));
        }
        out
    }

    #[cfg(test)]
    pub(super) fn mul_x_power(&self, e: usize) -> Self {
        let mut out = Self::default();
        for (key, c) in &self.0 {
            out.add_term(*key, c.mul_y_power(e));
        }
        out
    }

    /// Leading coefficient: largest `(shift, derivative)` key, then the
    /// lexicographically leading coefficient of its polynomial.
    pub(super) fn leading(&self) -> Option<Q> {
        self.0
            .values()
            .next_back()
            .and_then(|c| c.leading().cloned())
    }

    /// `(lead, e, N)` with `self = lead x^e N`, `N` having no common power of
    /// `x` and leading coefficient one.
    pub(super) fn normalize(&self) -> (Q, usize, Self) {
        let e = self.x_valuation();
        let stripped = self.div_x_power(e);
        let lead = stripped.leading().unwrap_or_else(Q::one);
        let normalized = stripped.scale(&lead.recip());
        (lead, e, normalized)
    }

    /// `(self Q)_n` when rows `Q_(origin)..` cover every referenced index.
    pub(super) fn apply(&self, rows: &[Vec<Q>], origin: usize, n: usize) -> Option<Vec<Q>> {
        let value = q_usize(n);
        let mut sum = vec![Q::zero()];
        for ((i, d), c) in &self.0 {
            let index = n.checked_add(*i)?.checked_sub(origin)?;
            let row = rows.get(index)?;
            let mut derivative = row.clone();
            for _ in 0..*d {
                derivative = vec_poly::derivative(&derivative);
            }
            sum = vec_poly::add(&sum, &vec_poly::mul(&c.eval_n(&value), &derivative));
        }
        Some(vec_poly::trim(&sum))
    }

    /// Apply a shift-free operator to one row at index `n`.
    pub(super) fn apply_row(&self, row: &[Q], n: usize) -> Vec<Q> {
        let value = q_usize(n);
        let mut sum = vec![Q::zero()];
        for ((_, d), c) in &self.0 {
            let mut derivative = row.to_vec();
            for _ in 0..*d {
                derivative = vec_poly::derivative(&derivative);
            }
            sum = vec_poly::add(&sum, &vec_poly::mul(&c.eval_n(&value), &derivative));
        }
        vec_poly::trim(&sum)
    }
}

// ---------------------------------------------------------------------------
// Euler form
// ---------------------------------------------------------------------------

/// `sum x^s a_(i,s)(n,theta) S^i`, keyed by `(shift i, x power s)`, with
/// coefficients in `(n, theta)`.
pub(super) type ThetaOp = BTreeMap<(usize, usize), Biv>;

pub(super) fn theta_add(op: &mut ThetaOp, key: (usize, usize), c: Biv) {
    if c.is_zero() {
        return;
    }
    let entry = op.entry(key).or_default();
    *entry = entry.add(&c);
    if entry.is_zero() {
        op.remove(&key);
    }
}

/// `theta (theta - 1) ... (theta - d + 1)`.
pub(super) fn falling_factorial(d: usize) -> Vec<Q> {
    (0..d).fold(vec![Q::one()], |acc, t| {
        vec_poly::mul(&acc, &[-q_usize(t), Q::one()])
    })
}

/// Stirling numbers of the second kind `S(m, d)` for `m <= max`.
fn stirling2(max: usize) -> Vec<Vec<Q>> {
    let mut table = vec![vec![Q::zero(); max + 1]; max + 1];
    table[0][0] = Q::one();
    for m in 1..=max {
        for d in 1..=m {
            table[m][d] = &table[m - 1][d - 1] + q_usize(d) * &table[m - 1][d];
        }
    }
    table
}

/// `(J, x^J L)` in Euler form with `J` the largest derivative order, using
/// `x^d D^d = theta (theta - 1) ... (theta - d + 1)`.
pub(super) fn to_theta(op: &DOp) -> (usize, ThetaOp) {
    let j_max = op.max_derivative();
    let mut out = ThetaOp::new();
    for ((i, d), c) in &op.0 {
        let falling = falling_factorial(*d);
        for (r, row) in c.rows().iter().enumerate() {
            for (j, cj) in row.iter().enumerate() {
                if cj.is_zero() {
                    continue;
                }
                let s = j + j_max - d;
                let mut rows = vec![Vec::new(); r + 1];
                rows[r] = vec_poly::scale(&falling, cj);
                theta_add(&mut out, (*i, s), Biv::new(rows));
            }
        }
    }
    (j_max, out)
}

/// The normal-ordered form of an Euler-form operator, using
/// `theta^m = sum_d S(m, d) x^d D^d`.
pub(super) fn from_theta(op: &ThetaOp) -> DOp {
    let max_m = op.values().filter_map(Biv::y_degree).max().unwrap_or(0);
    let table = stirling2(max_m);
    let mut out = DOp::default();
    for ((i, s), a) in op {
        for (r, row) in a.rows().iter().enumerate() {
            for (m, alpha) in row.iter().enumerate() {
                if alpha.is_zero() {
                    continue;
                }
                for (d, stirling) in table[m].iter().enumerate().take(m + 1) {
                    if stirling.is_zero() {
                        continue;
                    }
                    out.add_term((*i, d), Biv::monomial(r, s + d, alpha * stirling));
                }
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Zero sets of index polynomials
// ---------------------------------------------------------------------------

/// Nonnegative integer roots of a nonzero univariate rational polynomial.
pub(super) fn nonnegative_integer_roots(p: &[Q]) -> BTreeSet<BigInt> {
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
pub(super) struct ZeroAnalysis {
    /// Indices `n` of the domain where `p(n, .)` vanishes identically.
    pub(super) singular_indices: Vec<usize>,
    /// Indices `k >= 0` with `p(n, k) = 0` for every `n` of the domain.
    pub(super) uniform_indices: Vec<usize>,
    /// Whether the zero set is exactly `singular x N  union  domain x uniform`.
    pub(super) certified: bool,
    /// How the description was obtained.
    pub(super) certificate: String,
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
pub(super) fn content_and_primitive(p: &Biv) -> (Vec<Q>, Biv) {
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

pub(super) fn poly_text(p: &[Q], variable: &str) -> String {
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
pub(super) fn biv_text(p: &Biv, variable: &str) -> String {
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
pub(super) fn analyze_zeros(p: &Biv, first: usize, last: Option<usize>) -> ZeroAnalysis {
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
