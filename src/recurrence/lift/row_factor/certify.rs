//! The substitution `Q_n = F^(e(n)) S_n` into the input relation and the
//! replayable [`RowFactorCertificate`].

use super::super::hadamard_factor::prepare::{left_multiply, prepare_rule};
use crate::polynomial::Polynomial;
use crate::recurrence::implicit::algebra::{q, Biv, DOp};
use crate::recurrence::lift::section_certify::n_polynomial;
use crate::recurrence::{
    decode_poly, ImplicitRelation, ImplicitRelationJson, RecurrenceError, RecurrenceJson,
};
use crate::vec_poly;
use crate::BigRational;
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// `e(n) = slope n + offset` as a polynomial in `n`.
pub(super) fn exponent_polynomial(slope: usize, offset: i64) -> Biv {
    n_polynomial(&[q(offset), q(slope as i64)])
}

/// `T_d = F^d (D + w F'/F)^d = F^d F^(-w) D^d F^w` for `d = 0..=max`, by
/// `T_(d+1) = F D T_d + (w - d) F' T_d`; `w` is a polynomial in `n`.
fn conjugated_derivatives(f: &[Q], w: &Biv, max: usize) -> Vec<DOp> {
    let multiply_f = DOp::single(0, 0, Biv::from_y(f));
    let derivative = DOp::single(0, 1, Biv::one());
    let f_prime = Biv::from_y(&vec_poly::derivative(f));
    let mut powers = vec![DOp::single(0, 0, Biv::one())];
    for d in 0..max {
        let previous = &powers[d];
        let weight = w.add(&Biv::constant(-q(d as i64))).mul(&f_prime);
        let next = multiply_f
            .compose(&derivative.compose(previous))
            .add(&left_multiply(&weight, previous));
        powers.push(next);
    }
    powers
}

/// `M_S = sum c_(i,d) F^(slope i + J - d) T_d^(e(n) + slope i) S^i`, so that
/// `(M F^e S)_n = F^(e(n) - J) (M_S S)_n` as rational functions, with `J` the
/// largest derivative order of `M`.
pub(super) fn substitute(input: &DOp, f: &[Q], slope: usize, offset: i64) -> DOp {
    let j = input.max_derivative();
    let e = exponent_polynomial(slope, offset);
    let mut out = DOp::default();
    for ((i, d), c) in &input.0 {
        let w = e.add(&Biv::constant(q((slope * i) as i64)));
        let powers = conjugated_derivatives(f, &w, *d);
        let power = Biv::from_y(&vec_poly::pow(f, slope * i + j - d));
        let term = left_multiply(&c.mul(&power), &powers[*d]);
        for ((_, order), coefficient) in &term.0 {
            out.add_term((*i, *order), coefficient.clone());
        }
    }
    out
}

/// `p / g` for a polynomial `g` in `x`, row by row, when exact.
fn divide_biv(p: &Biv, g: &[Q]) -> Option<Biv> {
    let divisor = Polynomial::new(g.to_vec());
    let rows = p
        .rows()
        .iter()
        .map(|row| {
            if row.is_empty() {
                return Some(Vec::new());
            }
            let (quotient, remainder) = Polynomial::new(row.clone()).div_rem(&divisor);
            remainder.is_zero().then(|| quotient.coeffs().to_vec())
        })
        .collect::<Option<Vec<_>>>()?;
    Some(Biv::new(rows))
}

/// `op / F^h` coefficientwise, when exact.
pub(super) fn divide_operator(op: &DOp, f: &[Q], h: usize) -> Option<DOp> {
    let divisor = vec_poly::pow(f, h);
    let mut out = DOp::default();
    for (key, c) in &op.0 {
        out.add_term(*key, divide_biv(c, &divisor)?);
    }
    Some(out)
}

/// The exact identity behind a certified row-factor lift
/// `Q_n = F(x)^(e(n)) S_n`, `e(n) = exponent_slope n + exponent_offset`:
/// substituting into the input relation `M` gives `M_S` with
/// `(M F^e S)_n = F^(e(n) - J) (M_S S)_n`, and `M_S = F^h M'_S` with
/// `h = exponent_slope R + J`, so the leader of `M'_S` is the input's
/// `g(n) S^R`, free of `x`.  [`verify`](Self::verify) recomputes `M_S`, the
/// exact division and the parent relation; the all-index argument is in the
/// `row_factor` module documentation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RowFactorCertificate {
    /// `F`, ascending exact coefficients.
    pub factor: Vec<String>,
    /// The slope of `e(n)`.
    pub exponent_slope: usize,
    /// The offset of `e(n)` (`e(n0) >= 0`).
    pub exponent_offset: i64,
    /// `M`: the input rule at base index `n`.
    pub input_relation: ImplicitRelationJson,
    /// The largest offset `R` of the input rule.
    pub input_leader_shift: usize,
    /// `J`, the largest derivative order of `M`.
    pub derivative_order: usize,
    /// `h = exponent_slope R + J`, the power of `F` removed.
    pub removed_factor_power: usize,
    /// `M'_S = M_S / F^h`, the rule of the cofactor rows `S_n`.
    pub cofactor_relation: ImplicitRelationJson,
    /// The cofactor `S`: rule (equal to `M'_S`), first index `n0` and seeds
    /// `S_n = Q_n / F^(e(n))`.
    pub cofactor: RecurrenceJson,
    /// First base index at which `M` holds for the input rows.
    pub input_valid_from_base: usize,
    /// First base index at which `M'_S` holds for the cofactor rows.
    pub cofactor_valid_from_base: usize,
    /// Bases `input_valid_from_base <= n < cofactor_valid_from_base`, where
    /// `M` was evaluated exactly on the candidate rows and vanished.
    pub exceptional_bases: Vec<usize>,
    /// The identifying input indices `n0, ..., n0 + K - 1`.
    pub identifying_indices: Vec<usize>,
    /// What the certificate proves.
    pub statement: String,
}

impl RowFactorCertificate {
    /// Recompute `M_S` from the input relation, `F` and `e`, check
    /// `M_S = F^h M'_S` exactly with `h = exponent_slope R + J`, that the
    /// leader of `M'_S` is free of `x` and `D`, that the cofactor record's
    /// rule is `M'_S` with a leader nonzero at every target index, and that
    /// `e(n0) >= 0`.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        let Some(&first_index) = self.identifying_indices.first() else {
            return failed("the certificate lists no identifying indices");
        };
        if (first_index as i64)
            .checked_mul(self.exponent_slope as i64)
            .and_then(|v| v.checked_add(self.exponent_offset))
            .is_none_or(|v| v < 0)
        {
            return failed("e(n0) is negative");
        }
        let f = decode_poly(&self.factor)?;
        if vec_poly::degree(&f).unwrap_or(0) == 0 {
            return failed("the factor is constant");
        }
        let input = ImplicitRelation::from_json(&self.input_relation)?.op();
        let r = input.max_shift();
        let j = input.max_derivative();
        if r != self.input_leader_shift
            || j != self.derivative_order
            || input.column(r).max_derivative() != 0
        {
            return failed("the input relation does not have the recorded shape");
        }
        if self.removed_factor_power != self.exponent_slope * r + j {
            return failed("the removed power is not slope R + J");
        }
        let substituted = substitute(&input, &f, self.exponent_slope, self.exponent_offset);
        let cofactor = ImplicitRelation::from_json(&self.cofactor_relation)?.op();
        let power = Biv::from_y(&vec_poly::pow(&f, self.removed_factor_power));
        if left_multiply(&power, &cofactor) != substituted {
            return failed("M_S is not F^h times the cofactor relation");
        }
        let leader = cofactor.column(r);
        if leader.0.len() != 1 || !leader.0.get(&(0, 0)).is_some_and(Biv::is_y_free) {
            return failed("the cofactor leader is not free of x and D");
        }
        let (rule, cofactor_first, rows) = self.cofactor.to_recurrence_parts()?;
        if cofactor_first != first_index {
            return failed("the cofactor does not start at the input's first index");
        }
        let prepared = prepare_rule(&rule, &rows, cofactor_first, "cofactor")?;
        if prepared.relation != cofactor {
            return failed("the cofactor record's rule is not the cofactor relation");
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `F^J (M F^e S)_n = F^(e(n)) (M_S S)_n` on arbitrary rows.
    #[test]
    fn substitution_identity_on_arbitrary_rows() {
        let mut m = DOp::default();
        m.add_term((2, 0), Biv::new(vec![vec![q(2)], vec![q(1)]]));
        m.add_term((1, 1), Biv::new(vec![vec![q(0), q(1), q(-1)]]));
        m.add_term((1, 0), Biv::new(vec![vec![q(1), q(3)], vec![q(0), q(1)]]));
        m.add_term((0, 2), Biv::new(vec![vec![q(0), q(0), q(1)]]));
        let f = vec![q(1), q(2), q(1)];
        let (slope, offset) = (1usize, 1i64);
        let substituted = substitute(&m, &f, slope, offset);
        let rows: Vec<Vec<Q>> = (0..7i64).map(|n| vec![q(n), q(1 - n), q(2)]).collect();
        let n0 = 3usize;
        let lifted: Vec<Vec<Q>> = rows
            .iter()
            .enumerate()
            .map(|(k, s)| vec_poly::mul(&vec_poly::pow(&f, slope * (n0 + k) + offset as usize), s))
            .collect();
        let j = m.max_derivative();
        for base in n0..n0 + 4 {
            let left = vec_poly::mul(&vec_poly::pow(&f, j), &m.apply(&lifted, n0, base).unwrap());
            let e = slope * base + offset as usize;
            let right = vec_poly::mul(
                &vec_poly::pow(&f, e),
                &substituted.apply(&rows, n0, base).unwrap(),
            );
            assert_eq!(vec_poly::trim(&left), vec_poly::trim(&right));
        }
    }
}
