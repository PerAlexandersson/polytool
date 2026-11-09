//! [`BasisMap`]: the linear map `T_P(sum p_k u^k) = sum p_k P_k(t)` of a
//! recurrence-defined family, its application, its triangular inverse and
//! its intertwining relation with multiplication by `u`.

use super::super::section_operators::OperatorMatrix;
use crate::recurrence::implicit::algebra::{Biv, BivarConversion, DOp};
use crate::recurrence::{
    format_rational_coeff, Recurrence, RecurrenceError, RecurrenceJson, RecurrenceSign,
};
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_traits::{One, Zero};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Largest number of basis rows generated for one application or
/// expansion.
pub const MAX_BASIS_MAP_ROWS: usize = 512;

/// Rows replayed when a basis map is constructed (a check that the rule can
/// be evaluated; not part of any proof).
const CONSTRUCTION_ROWS: usize = 6;

/// The basis map `T_P : u^k -> P_k(t)` of a polynomial family `P`.
///
/// The family is given by a rule `f(n, t) P_n = sum_j c_j(n, t) D^(d_j)
/// P_(n - r_j)` (homogeneous, without `(-1)^n` factors) with its first index
/// `i0` and identifying seed rows, exactly as [`Recurrence::generate_rows_rational`]
/// replays it; the basis element of `u^k` is the row with index `i0 + k`:
///
/// ```text
/// T_P(sum_k p_k u^k) = sum_k p_k P_(i0 + k)(t).
/// ```
///
/// A different first index gives a different map (a *tail* of the family,
/// [`BasisMap::tail`]).  Construction drops trailing seed rows that the rule
/// reproduces from the earlier ones, so [`BasisMap::seeds`] is the minimal
/// identifying prefix.  The JSON form is `{ "name": ..., "family":
/// polytool.recurrence.v1 }`, so that a map is a first-class value that can
/// be stored and replayed.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(try_from = "BasisMapJson", into = "BasisMapJson")]
pub struct BasisMap {
    name: String,
    rule: Recurrence,
    first_index: usize,
    seeds: Vec<Vec<Q>>,
}

/// The JSON form of a [`BasisMap`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapJson {
    /// A short name (for example an OEIS A-number); free text.
    #[serde(default)]
    pub name: String,
    /// The family: rule, first index (the index of the basis element of
    /// `u^0`) and seed rows.
    pub family: RecurrenceJson,
}

impl TryFrom<BasisMapJson> for BasisMap {
    type Error = RecurrenceError;

    fn try_from(json: BasisMapJson) -> Result<Self, Self::Error> {
        let (rule, first_index, seeds) = json.family.to_recurrence_parts()?;
        Self::new(&json.name, rule, first_index, seeds)
    }
}

impl From<BasisMap> for BasisMapJson {
    fn from(map: BasisMap) -> Self {
        Self {
            family: map.recurrence_json(),
            name: map.name,
        }
    }
}

impl PartialEq for BasisMap {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name && self.recurrence_json() == other.recurrence_json()
    }
}

impl Eq for BasisMap {}

impl BasisMap {
    /// The map of the family `rule` with first index `first_index` and seed
    /// rows `seeds` (`seeds[0]` is the basis element of `u^0`).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] or [`RecurrenceError::Unsupported`]
    /// when the rule has a term of offset zero, a `(-1)^n` factor, an
    /// inhomogeneous term or a zero denominator, when there are fewer seeds
    /// than the order, or when the first rows cannot be replayed.
    pub fn new(
        name: &str,
        rule: Recurrence,
        first_index: usize,
        seeds: Vec<Vec<Q>>,
    ) -> Result<Self, RecurrenceError> {
        if rule.terms.iter().any(|t| t.offset == 0) {
            return Err(RecurrenceError::InvalidInput(
                "a basis-map rule term has offset 0; offsets must be positive".into(),
            ));
        }
        if rule.terms.iter().any(|t| t.sign != RecurrenceSign::None) {
            return Err(RecurrenceError::Unsupported(
                "a basis-map rule has a (-1)^n factor".into(),
            ));
        }
        if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
            return Err(RecurrenceError::Unsupported(
                "a basis-map rule has an inhomogeneous term".into(),
            ));
        }
        if rule.denominator.as_ref().is_some_and(|f| f.is_zero()) {
            return Err(RecurrenceError::InvalidInput(
                "the basis-map rule has a zero denominator".into(),
            ));
        }
        let order = rule.max_offset();
        if seeds.is_empty() || seeds.len() < order {
            return Err(RecurrenceError::InvalidInput(format!(
                "the basis-map rule has order {order} but {} seed rows were supplied",
                seeds.len()
            )));
        }
        let mut seeds: Vec<Vec<Q>> = seeds.iter().map(|row| vec_poly::trim(row)).collect();
        // Drop trailing seeds that the rule reproduces.
        while seeds.len() > order.max(1) {
            let last = seeds.len() - 1;
            match rule.evaluate_next_rational(&seeds[..last], first_index) {
                Ok(row) if vec_poly::trim(&row) == seeds[last] => {
                    seeds.pop();
                }
                _ => break,
            }
        }
        let map = Self {
            name: name.to_string(),
            rule,
            first_index,
            seeds,
        };
        map.basis_rows(CONSTRUCTION_ROWS)?;
        Ok(map)
    }

    /// The map of a bundled OEIS catalog entry, `u^k` sent to the row with
    /// displayed index `first_row + k` (the catalog's first row).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::InvalidInput`] for an unknown A-number, a negative
    /// first row or a catalog rule that [`BasisMap::new`] rejects.
    pub fn from_oeis(id: &str) -> Result<Self, RecurrenceError> {
        let definition = crate::oeis::by_id(id).ok_or_else(|| {
            RecurrenceError::InvalidInput(format!("{id} is not in the bundled OEIS catalog"))
        })?;
        let (rule, rows, first_row) = definition.recurrence_export_parts()?;
        let first_index = usize::try_from(first_row).map_err(|_| {
            RecurrenceError::Unsupported(format!("{id} starts at the negative row {first_row}"))
        })?;
        Self::new(id, rule, first_index, rows)
    }

    /// The name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The family's rule.
    pub fn rule(&self) -> &Recurrence {
        &self.rule
    }

    /// The index of the basis element of `u^0`.
    pub fn first_index(&self) -> usize {
        self.first_index
    }

    /// The minimal identifying seed rows.
    pub fn seeds(&self) -> &[Vec<Q>] {
        &self.seeds
    }

    /// The order of the rule (largest offset).
    pub fn order(&self) -> usize {
        self.rule.max_offset()
    }

    /// The family as `polytool.recurrence.v1`.
    pub fn recurrence_json(&self) -> RecurrenceJson {
        RecurrenceJson::from_recurrence_rational(&self.rule, self.first_index, &self.seeds, None)
    }

    /// Whether no coefficient of the rule (denominator included) depends on
    /// the index `n`.
    pub fn is_index_independent(&self) -> bool {
        let free = |p: &crate::recurrence::BivarPoly| {
            p.coeffs
                .iter()
                .skip(1)
                .all(|row| row.iter().all(Zero::is_zero))
        };
        self.rule.terms.iter().all(|t| free(&t.coeff))
            && self.rule.denominator.as_ref().is_none_or(free)
    }

    /// The tail `u^k -> P_(i0 + shift + k)`.
    ///
    /// # Errors
    ///
    /// As [`BasisMap::new`], or when the rows cannot be generated.
    pub fn tail(&self, shift: usize) -> Result<Self, RecurrenceError> {
        if shift == 0 {
            return Ok(self.clone());
        }
        let count = self.seeds.len().max(self.order());
        let rows = self.basis_rows(shift + count)?;
        Self::new(
            &format!(
                "{} (tail from u^0 -> P_({}))",
                self.name,
                self.first_index + shift
            ),
            self.rule.clone(),
            self.first_index + shift,
            rows[shift..].to_vec(),
        )
    }

    /// The basis rows `P_(i0), ..., P_(i0 + count - 1)`, trimmed (the zero
    /// polynomial as `[0]`).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::LimitExceeded`] above [`MAX_BASIS_MAP_ROWS`] rows,
    /// [`RecurrenceError::InvalidInput`] when the rule cannot be replayed.
    pub fn basis_rows(&self, count: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        if count > MAX_BASIS_MAP_ROWS {
            return Err(RecurrenceError::LimitExceeded(format!(
                "a basis map generates at most {MAX_BASIS_MAP_ROWS} rows"
            )));
        }
        let rows = self
            .rule
            .generate_rows_rational(&self.seeds, self.first_index, count.max(self.seeds.len()))
            .map_err(|error| {
                RecurrenceError::InvalidInput(format!(
                    "the basis family cannot be replayed: {error}"
                ))
            })?;
        Ok(rows
            .into_iter()
            .take(count)
            .map(|row| vec_poly::trim(&row))
            .collect())
    }

    /// `T_P(p) = sum_k p_k P_(i0 + k)` for `p` in ascending coefficients.
    ///
    /// # Errors
    ///
    /// As [`BasisMap::basis_rows`].
    pub fn apply(&self, p: &[Q]) -> Result<Vec<Q>, RecurrenceError> {
        let p = vec_poly::trim(p);
        let rows = self.basis_rows(p.len())?;
        Ok(apply_with_rows(&p, &rows))
    }

    /// The triangular inverse: the coefficients `p_k` with
    /// `sum_k p_k P_(i0 + k) = q`, exactly.
    ///
    /// The basis rows are generated until the first row `P_(i0 + K)` with
    /// `deg P_(i0 + K) >= deg q`; the rows `P_(i0)..P_(i0 + K)` must be
    /// nonzero with strictly increasing degrees, so that the expansion in
    /// their span is unique and found by elimination from the top degree.
    /// The result is checked by applying the map again.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::Unsupported`] when the rows used are not triangular
    /// or `q` is not in their span (for example `q` not divisible by the
    /// common factor `g` of a geometric family `g y^k`), and the errors of
    /// [`BasisMap::basis_rows`].
    pub fn expand(&self, q: &[Q]) -> Result<Vec<Q>, RecurrenceError> {
        let q = vec_poly::trim(q);
        let degree = vec_poly::degree(&q).unwrap_or(0);
        let rows = self.triangular_rows_up_to(degree)?;
        let p = expand_with_rows(&q, &rows)?;
        if apply_with_rows(&p, &rows) != q {
            return Err(RecurrenceError::VerificationFailed(
                "the triangular expansion does not reproduce the polynomial".into(),
            ));
        }
        Ok(p)
    }

    /// Expand several polynomials with one set of basis rows.
    ///
    /// # Errors
    ///
    /// As [`BasisMap::expand`].
    pub fn expand_all(&self, polynomials: &[Vec<Q>]) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let degree = polynomials
            .iter()
            .filter_map(|p| vec_poly::degree(&vec_poly::trim(p)))
            .max()
            .unwrap_or(0);
        let rows = self.triangular_rows_up_to(degree)?;
        polynomials
            .iter()
            .map(|q| {
                let q = vec_poly::trim(q);
                let p = expand_with_rows(&q, &rows)?;
                if apply_with_rows(&p, &rows) != q {
                    return Err(RecurrenceError::VerificationFailed(
                        "the triangular expansion does not reproduce the polynomial".into(),
                    ));
                }
                Ok(p)
            })
            .collect()
    }

    /// Rows `P_(i0)..P_(i0 + K)` with `K` the first index of degree
    /// `>= degree`, required triangular.
    fn triangular_rows_up_to(&self, degree: usize) -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let mut count = (degree + 2).max(self.seeds.len() + 1);
        loop {
            let rows = self.basis_rows(count.min(MAX_BASIS_MAP_ROWS))?;
            if let Some(last) = rows
                .iter()
                .position(|row| vec_poly::degree(row).unwrap_or(0) >= degree && !is_zero_row(row))
            {
                let used = &rows[..=last];
                let mut previous: Option<usize> = None;
                for (k, row) in used.iter().enumerate() {
                    if is_zero_row(row) {
                        return Err(RecurrenceError::Unsupported(format!(
                            "the basis row of u^{k} is zero, so the basis map is not triangular"
                        )));
                    }
                    let d = vec_poly::degree(row).unwrap_or(0);
                    if previous.is_some_and(|p| d <= p) {
                        return Err(RecurrenceError::Unsupported(format!(
                            "the basis rows of u^{} and u^{k} do not have strictly increasing \
                             degrees, so the basis map is not triangular",
                            k - 1
                        )));
                    }
                    previous = Some(d);
                }
                return Ok(used.to_vec());
            }
            if count >= MAX_BASIS_MAP_ROWS {
                return Err(RecurrenceError::LimitExceeded(format!(
                    "no basis row among the first {MAX_BASIS_MAP_ROWS} reaches degree {degree}"
                )));
            }
            count *= 2;
        }
    }

    /// The intertwining relation of the map with multiplication by `u`
    /// (module documentation of `basis_map`, "Intertwining").
    pub fn intertwining_relation(&self) -> BasisMapIntertwining {
        intertwining_relation(self)
    }

    /// The operators `A_0, ..., A_(R-1)` (shift-free, in `(t, D_t)`, free
    /// of `n`) with `P_(k+R) = sum_j A_j P_(k+j)` for every `k >= 0`, for an
    /// index-independent rule with a constant nonzero leader and no seed row
    /// beyond the order.
    pub(super) fn companion_operators(&self) -> Result<Vec<DOp>, RecurrenceError> {
        let order = self.order();
        if order == 0 {
            return Err(RecurrenceError::Unsupported(
                "the basis rule has no terms".into(),
            ));
        }
        if !self.is_index_independent() {
            return Err(RecurrenceError::Unsupported(
                "the basis rule has index-dependent coefficients, so multiplication by u is \
                 not intertwined with a fixed operator (only the theta form holds)"
                    .into(),
            ));
        }
        if self.seeds.len() != order {
            return Err(RecurrenceError::Unsupported(format!(
                "the basis family has {} identifying seed rows for a rule of order {order} (an \
                 exceptional prefix), so the intertwining holds only on a subspace",
                self.seeds.len()
            )));
        }
        let leader = match &self.rule.denominator {
            None => Q::one(),
            Some(f) => Biv::from_bivar(f).as_constant().ok_or_else(|| {
                RecurrenceError::Unsupported(
                    "the basis rule's leader depends on t, so the next row is not an operator \
                     image of the previous ones"
                        .into(),
                )
            })?,
        };
        let mut operators = vec![DOp::default(); order];
        for term in &self.rule.terms {
            let coefficient = Biv::from_bivar(&term.coeff).scale(&leader.recip());
            operators[order - term.offset].add_term((0, term.deriv_order), coefficient);
        }
        Ok(operators)
    }

    /// The companion matrix `C` with `V(u p) = C V(p)`, `V(p) = (T(p),
    /// T(u p), ..., T(u^(R-1) p))`.
    pub(super) fn companion_matrix(&self) -> Result<OperatorMatrix, RecurrenceError> {
        let operators = self.companion_operators()?;
        let order = operators.len();
        let mut matrix = OperatorMatrix::zero(order);
        for i in 0..order.saturating_sub(1) {
            matrix.entries[i][i + 1] = DOp::single(0, 0, Biv::one());
        }
        matrix.entries[order - 1] = operators;
        Ok(matrix)
    }
}

/// `sum_k p_k rows[k]`.
pub(super) fn apply_with_rows(p: &[Q], rows: &[Vec<Q>]) -> Vec<Q> {
    let mut sum = vec![Q::zero()];
    for (coefficient, row) in p.iter().zip(rows) {
        if !coefficient.is_zero() {
            sum = vec_poly::add(&sum, &vec_poly::scale(row, coefficient));
        }
    }
    vec_poly::trim(&sum)
}

fn is_zero_row(row: &[Q]) -> bool {
    row.iter().all(Zero::is_zero)
}

/// Top-down elimination in the span of triangular `rows`.
fn expand_with_rows(q: &[Q], rows: &[Vec<Q>]) -> Result<Vec<Q>, RecurrenceError> {
    let mut remainder = vec_poly::trim(q);
    let mut p = vec![Q::zero(); rows.len()];
    for k in (0..rows.len()).rev() {
        if is_zero_row(&remainder) {
            break;
        }
        let r_degree = vec_poly::degree(&remainder).unwrap_or(0);
        let row = &rows[k];
        let b_degree = vec_poly::degree(row).unwrap_or(0);
        if r_degree > b_degree {
            return Err(RecurrenceError::Unsupported(format!(
                "the polynomial is not in the span of the basis rows (degree {r_degree} is \
                 missed)"
            )));
        }
        if r_degree < b_degree {
            continue;
        }
        let factor = &remainder[r_degree] / &row[b_degree];
        remainder = vec_poly::trim(&vec_poly::sub(&remainder, &vec_poly::scale(row, &factor)));
        p[k] = factor;
    }
    if !is_zero_row(&remainder) {
        return Err(RecurrenceError::Unsupported(
            "the polynomial is not in the span of the basis rows".into(),
        ));
    }
    Ok(vec_poly::trim(&p))
}

/// One term `t^i D_t^d T(phi(theta_u) u^m p)` of a
/// [`BasisMapIntertwining`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapIntertwiningTerm {
    /// The power `i` of `t` multiplying the term.
    pub t_power: usize,
    /// The order `d` of the derivative `D_t^d` applied after `T`.
    pub derivative_order: usize,
    /// The power `m` of `u` multiplying the input.
    pub u_power: usize,
    /// `phi(theta)`, ascending exact coefficients, with `theta = u d/du`
    /// applied to `u^m p` before `T`.
    pub theta_polynomial: Vec<String>,
}

/// The intertwining relation of a [`BasisMap`] with multiplication by `u`:
///
/// ```text
/// sum_terms t^i D_t^d T(phi(theta_u) (u^m p)) = 0
/// ```
///
/// for every polynomial `p` divisible by `u^minimal_input_power` (module
/// documentation of `basis_map`, "Intertwining").  For an index-independent
/// rule the `theta` polynomials are constants and the relation reads
/// `f T(u^R p) = sum_j c_j(t) D^(d_j) T(u^(R - r_j) p)`, which for a rule
/// with a constant leader expresses `T o u^R` through `T o u^j`, `j < R`
/// (the companion form used by the lifts).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapIntertwining {
    /// The order `R` of the rule.
    pub order: usize,
    /// The relation holds for every `p` divisible by `u` to this power
    /// (`K - R` for `K` identifying seeds; `0` without an exceptional
    /// prefix).
    pub minimal_input_power: usize,
    /// Whether every `theta` polynomial is constant (the rule's
    /// coefficients are free of `n`).
    pub index_independent: bool,
    /// The terms (merged by `(t_power, derivative_order, u_power)`).
    pub terms: Vec<BasisMapIntertwiningTerm>,
    /// The relation in words.
    pub statement: String,
}

impl BasisMapIntertwining {
    /// Check the relation on the monomials `p = u^j`,
    /// `minimal_input_power <= j < minimal_input_power + count`, against the
    /// basis rows of `map` (a regression check of the derivation; the
    /// relation itself is an algebraic identity, module documentation).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] naming the first failing
    /// monomial, and the errors of [`BasisMap::basis_rows`].
    pub fn check_on_monomials(&self, map: &BasisMap, count: usize) -> Result<(), RecurrenceError> {
        let max_power = self.terms.iter().map(|t| t.u_power).max().unwrap_or(0);
        let rows = map.basis_rows(self.minimal_input_power + count + max_power)?;
        let theta: Vec<Vec<Q>> = self
            .terms
            .iter()
            .map(|t| {
                t.theta_polynomial
                    .iter()
                    .map(|c| crate::recurrence::parse_rational_coeff(c))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<_, _>>()?;
        for j in self.minimal_input_power..self.minimal_input_power + count {
            let mut sum = vec![Q::zero()];
            for (term, phi) in self.terms.iter().zip(&theta) {
                let power = term.u_power + j;
                let value = vec_poly::evaluate(phi, &Q::from_integer(BigInt::from(power)));
                let mut row = vec_poly::scale(&rows[power], &value);
                for _ in 0..term.derivative_order {
                    row = vec_poly::derivative(&row);
                }
                sum = vec_poly::add(&sum, &vec_poly::shift_by(&row, term.t_power));
            }
            if !is_zero_row(&vec_poly::trim(&sum)) {
                return Err(RecurrenceError::VerificationFailed(format!(
                    "the intertwining relation fails on u^{j}"
                )));
            }
        }
        Ok(())
    }
}

/// `sum_e c_e (theta + shift)^e` as a polynomial in `theta`.
fn shifted_in_theta(coefficients_in_n: &[Q], shift: &Q) -> Vec<Q> {
    // Horner in (theta + shift).
    let linear = [shift.clone(), Q::one()];
    let mut result = vec![Q::zero()];
    for c in coefficients_in_n.iter().rev() {
        result = vec_poly::mul(&result, &linear);
        result = vec_poly::add(&result, std::slice::from_ref(c));
    }
    vec_poly::trim(&result)
}

fn intertwining_relation(map: &BasisMap) -> BasisMapIntertwining {
    use std::collections::BTreeMap;
    let rule = &map.rule;
    let order = rule.max_offset();
    let i0 = Q::from_integer(BigInt::from(map.first_index));
    let mut merged: BTreeMap<(usize, usize, usize), Vec<Q>> = BTreeMap::new();
    let mut add = |key: (usize, usize, usize), phi: Vec<Q>| {
        let entry = merged.entry(key).or_insert_with(|| vec![Q::zero()]);
        *entry = vec_poly::trim(&vec_poly::add(entry, &phi));
    };
    // Coefficient of t^i as a polynomial in n.
    let column = |p: &crate::recurrence::BivarPoly, i: usize| -> Vec<Q> {
        p.coeffs
            .iter()
            .map(|row| row.get(i).cloned().unwrap_or_else(Q::zero))
            .collect()
    };
    let width = |p: &crate::recurrence::BivarPoly| p.coeffs.iter().map(Vec::len).max().unwrap_or(0);
    let leader = rule
        .denominator
        .clone()
        .unwrap_or_else(|| crate::recurrence::BivarPoly {
            coeffs: vec![vec![Q::one()]],
        });
    for i in 0..width(&leader) {
        add((i, 0, order), shifted_in_theta(&column(&leader, i), &i0));
    }
    for term in &rule.terms {
        let shift = &i0 + Q::from_integer(BigInt::from(term.offset));
        for i in 0..width(&term.coeff) {
            let phi = shifted_in_theta(&column(&term.coeff, i), &shift);
            add(
                (i, term.deriv_order, order - term.offset),
                vec_poly::neg(&phi),
            );
        }
    }
    let terms: Vec<BasisMapIntertwiningTerm> = merged
        .into_iter()
        .filter(|(_, phi)| !is_zero_row(phi))
        .map(|((i, d, m), phi)| BasisMapIntertwiningTerm {
            t_power: i,
            derivative_order: d,
            u_power: m,
            theta_polynomial: phi.iter().map(format_rational_coeff).collect(),
        })
        .collect();
    let index_independent = terms.iter().all(|t| t.theta_polynomial.len() <= 1);
    let minimal_input_power = map.seeds.len().saturating_sub(order);
    let statement = format!(
        "sum over the terms of t^i D_t^d T(phi(theta) u^m p) = 0 for every polynomial p \
         divisible by u^{minimal_input_power}, theta = u d/du acting before T: the rule \
         f(n, t) P_n = sum_j c_j(n, t) D^(d_j) P_(n - r_j), valid from n = {}, applied to the \
         monomials u^k of u^{order} p (n = k + {} + {}, so n -> theta + {} on u^(k + R - r_j) \
         after the shift r_j)",
        map.first_index + map.seeds.len(),
        order,
        map.first_index,
        map.first_index,
    );
    BasisMapIntertwining {
        order,
        minimal_input_power,
        index_independent,
        terms,
        statement,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::{BivarPoly, RecurrenceTerm};

    fn q(v: i64) -> Q {
        Q::from_integer(BigInt::from(v))
    }

    fn bivar(rows: Vec<Vec<i64>>) -> BivarPoly {
        BivarPoly::from_coefficients_normalized(
            rows.into_iter()
                .map(|r| r.into_iter().map(q).collect())
                .collect(),
        )
    }

    /// A035607: P_0 = 1, P_1 = 1 + 2t, P_(n+2) = (1 + t) P_(n+1) + t P_n.
    fn a035607() -> BasisMap {
        let rule = Recurrence {
            terms: vec![
                RecurrenceTerm::new(1, 0, bivar(vec![vec![1, 1]])),
                RecurrenceTerm::new(2, 0, bivar(vec![vec![0, 1]])),
            ],
            denominator: None,
            inhomogeneous: None,
        };
        BasisMap::new("A035607", rule, 0, vec![vec![q(1)], vec![q(1), q(2)]]).unwrap()
    }

    #[test]
    fn apply_and_expand_are_inverse() {
        let map = a035607();
        // T((u + 4)^3) = 125 + 150 t + 32 t^2 + 2 t^3 (the survey's failure).
        let image = map.apply(&[q(64), q(48), q(12), q(1)]).unwrap();
        assert_eq!(image, vec![q(125), q(150), q(32), q(2)]);
        let back = map.expand(&image).unwrap();
        assert_eq!(back, vec![q(64), q(48), q(12), q(1)]);
        let relation = map.intertwining_relation();
        assert!(relation.index_independent);
        relation.check_on_monomials(&map, 8).unwrap();
    }

    #[test]
    fn an_index_dependent_rule_has_theta_terms() {
        // Rising factorial basis: P_n = (t + n - 1) P_(n-1), P_0 = 1 (first index 0).
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(1, 0, bivar(vec![vec![-1, 1], vec![1]]))],
            denominator: None,
            inhomogeneous: None,
        };
        let map = BasisMap::new("rising", rule, 0, vec![vec![q(1)]]).unwrap();
        let relation = map.intertwining_relation();
        assert!(!relation.index_independent);
        relation.check_on_monomials(&map, 8).unwrap();
        // The expansion of t^3 in the rising-factorial basis: Stirling numbers.
        assert_eq!(
            map.expand(&[q(0), q(0), q(0), q(1)]).unwrap(),
            vec![q(0), q(1), q(-3), q(1)]
        );
    }

    #[test]
    fn redundant_seeds_are_dropped_and_tails_shift() {
        let rule = Recurrence {
            terms: vec![RecurrenceTerm::new(1, 0, bivar(vec![vec![1, 1]]))],
            denominator: None,
            inhomogeneous: None,
        };
        let map = BasisMap::new(
            "binomial",
            rule,
            0,
            vec![vec![q(1)], vec![q(1), q(1)], vec![q(1), q(2), q(1)]],
        )
        .unwrap();
        assert_eq!(map.seeds().len(), 1);
        let tail = map.tail(2).unwrap();
        assert_eq!(tail.seeds().len(), 1);
        assert_eq!(tail.basis_rows(1).unwrap()[0], vec![q(1), q(2), q(1)]);
        assert!(map.expand(&[q(1), q(0), q(1)]).is_ok());
        // Not in the span of the tail rows (1 + t)^(k + 2).
        assert!(tail.expand(&[q(1), q(0), q(1)]).is_err());
        let json = serde_json::to_string(&map).unwrap();
        let back: BasisMap = serde_json::from_str(&json).unwrap();
        assert_eq!(back, map);
    }
}
