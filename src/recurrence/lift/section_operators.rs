//! The transport map `phi` from operators on a polynomial sequence to
//! operator matrices on its Veronese sections.
//!
//! # Convention
//!
//! For a modulus `r >= 2` write a polynomial as
//! `P(x) = sum_(a=0)^(r-1) x^a F_a(x^r)` and put `y = x^r`, so that
//! `F_a = V_(r, a)(P)` and the section vector is `F = (F_0, ..., F_(r-1))`.
//! An operator matrix `A` acts by `(A F)_b = sum_a A[b][a] F_a` (row `b`,
//! column `a`), and matrices compose as `(A B)[b][a] = sum_c A[b][c] o B[c][a]`.
//!
//! - Multiplication by `x` is the matrix `X` with `(X F)_0 = y F_(r-1)` and
//!   `(X F)_b = F_(b-1)` for `b >= 1`: `X[0][r-1] = y`, `X[b][b-1] = 1`.
//! - Differentiation `D_x` is the matrix `DD` with
//!   `(DD F)_b = ((b+1) + r y D_y) F_(b+1)` for `b < r-1` and
//!   `(DD F)_(r-1) = r D_y F_0`: `DD[b][b+1] = (b+1) + r y D_y`,
//!   `DD[r-1][0] = r D_y`.  (From `D_x (x^a G(x^r)) = a x^(a-1) G(x^r) +
//!   r x^(a-1) y G'(y)` for `a >= 1` and `D_x G(x^r) = r x^(r-1) G'(y)`.)
//! - The index `n` and the shift `S` act diagonally: `n -> n I`, `S -> S I`.
//!
//! `phi(sum c_(i,d)(n, x) D^d S^i) = sum c_(i,d)(n, X) DD^d S^i`, where
//! `c(n, X) = sum_e c_e(n) X^e` and `X^e[b][a] = y^((a+e-b)/r)` when
//! `a + e = b (mod r)` (zero otherwise).  `phi` is a homomorphism of Ore
//! algebras because `DD X - X DD = I` and `X`, `DD` commute with `n` and `S`;
//! since the sections of a polynomial are unique,
//! `L P = 0` holds at an index exactly when `phi(L) F = 0` holds there:
//! the sections of `(L P)_n` are `(phi(L) F)_n`.
//!
//! In Euler form the same map sends `x^s a(n, theta_x)` (with
//! `theta_x = x D_x`) to the matrix that maps component `b` to component
//! `(b + s) mod r` with the factor `y^floor((b+s)/r) a(n, r theta_y + b)`;
//! the tests cross-check both descriptions.
//!
//! Related: `implicit::algebra` (the scalar Ore algebra `DOp`),
//! `section_certify` (pseudo-division by a transported parent system).

use super::super::implicit::algebra::{q_usize, Biv, DOp};
use super::*;
use crate::vec_poly;
use num_traits::Zero;

/// A square matrix of operators in `(n, y, D_y, S)`; `entries[b][a]` is the
/// entry in row `b` and column `a`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct OperatorMatrix {
    pub(super) entries: Vec<Vec<DOp>>,
}

impl OperatorMatrix {
    pub(super) fn zero(size: usize) -> Self {
        Self {
            entries: vec![vec![DOp::default(); size]; size],
        }
    }

    pub(super) fn identity(size: usize) -> Self {
        let mut matrix = Self::zero(size);
        for (b, row) in matrix.entries.iter_mut().enumerate() {
            row[b] = DOp::single(0, 0, Biv::one());
        }
        matrix
    }

    pub(super) fn size(&self) -> usize {
        self.entries.len()
    }

    pub(super) fn add(&self, other: &Self) -> Self {
        Self {
            entries: self
                .entries
                .iter()
                .zip(&other.entries)
                .map(|(left, right)| left.iter().zip(right).map(|(a, b)| a.add(b)).collect())
                .collect(),
        }
    }

    #[cfg(test)]
    pub(super) fn sub(&self, other: &Self) -> Self {
        Self {
            entries: self
                .entries
                .iter()
                .zip(&other.entries)
                .map(|(left, right)| left.iter().zip(right).map(|(a, b)| a.sub(b)).collect())
                .collect(),
        }
    }

    /// The product `self o other` (entries composed in the Ore algebra).
    pub(super) fn compose(&self, other: &Self) -> Self {
        let size = self.size();
        let mut out = Self::zero(size);
        for b in 0..size {
            for a in 0..size {
                let mut sum = DOp::default();
                for c in 0..size {
                    if self.entries[b][c].is_zero() || other.entries[c][a].is_zero() {
                        continue;
                    }
                    sum = sum.add(&self.entries[b][c].compose(&other.entries[c][a]));
                }
                out.entries[b][a] = sum;
            }
        }
        out
    }

    /// `self o S^shift`: every entry is multiplied by `S^shift` on the right.
    pub(super) fn times_shift(&self, shift: usize) -> Self {
        Self {
            entries: self
                .entries
                .iter()
                .map(|row| row.iter().map(|op| right_shift(op, shift)).collect())
                .collect(),
        }
    }

    /// Total number of stored operator terms.
    pub(super) fn term_count(&self) -> usize {
        self.entries.iter().flatten().map(|op| op.0.len()).sum()
    }

    /// `(self F)_b` at base index `n`, where `sections[a][i]` is component
    /// `a` of the row with index `origin + i`; `None` when a referenced row
    /// is missing.
    #[cfg(test)]
    pub(super) fn apply(
        &self,
        sections: &[Vec<Vec<BigRational>>],
        origin: usize,
        n: usize,
    ) -> Option<Vec<Vec<BigRational>>> {
        self.entries
            .iter()
            .map(|row| apply_row_operator(row, sections, origin, n))
            .collect()
    }
}

/// `op o S^shift`: the shift key of every term grows by `shift` (the
/// coefficient stands on the left, so nothing is substituted).
pub(super) fn right_shift(op: &DOp, shift: usize) -> DOp {
    let mut out = DOp::default();
    for ((i, d), c) in &op.0 {
        out.add_term((i + shift, *d), c.clone());
    }
    out
}

/// `sum_a row[a] F_a` at base index `n` (see [`OperatorMatrix::apply`]).
#[cfg(test)]
pub(super) fn apply_row_operator(
    row: &[DOp],
    sections: &[Vec<Vec<BigRational>>],
    origin: usize,
    n: usize,
) -> Option<Vec<BigRational>> {
    let mut sum = vec![BigRational::zero()];
    for (op, component) in row.iter().zip(sections) {
        if op.is_zero() {
            continue;
        }
        sum = vec_poly::add(&sum, &op.apply(component, origin, n)?);
    }
    Some(vec_poly::trim(&sum))
}

/// The multiplication matrix `X^e` (a matrix of multiplication operators).
#[cfg(test)]
pub(super) fn power_of_multiplication_matrix(modulus: usize, exponent: usize) -> OperatorMatrix {
    let mut matrix = OperatorMatrix::zero(modulus);
    for b in 0..modulus {
        for a in 0..modulus {
            if (a + exponent) % modulus == b {
                let y_power = (a + exponent - b) / modulus;
                matrix.entries[b][a] = DOp::single(0, 0, Biv::monomial(0, y_power, q_usize(1)));
            }
        }
    }
    matrix
}

/// The derivative matrix `DD` of the module documentation.
pub(super) fn derivative_matrix(modulus: usize) -> OperatorMatrix {
    let r = q_usize(modulus);
    let mut matrix = OperatorMatrix::zero(modulus);
    for b in 0..modulus - 1 {
        let mut entry = DOp::single(0, 0, Biv::constant(q_usize(b + 1)));
        entry.add_term((0, 1), Biv::monomial(0, 1, r.clone()));
        matrix.entries[b][b + 1] = entry;
    }
    matrix.entries[modulus - 1][0] = DOp::single(0, 1, Biv::constant(r));
    matrix
}

/// `c(n, X)` for a coefficient `c(n, x)`: entry `[b][a]` collects
/// `c_e(n) y^((a+e-b)/r)` over `e = b - a (mod r)`.
fn transport_coefficient(coefficient: &Biv, modulus: usize) -> OperatorMatrix {
    let mut matrix = OperatorMatrix::zero(modulus);
    let Some(degree) = coefficient.y_degree() else {
        return matrix;
    };
    for e in 0..=degree {
        let column = coefficient.y_coefficient(e);
        if column.is_empty() {
            continue;
        }
        for b in 0..modulus {
            let a = (b + modulus * (e / modulus + 1) - e) % modulus;
            let y_power = (a + e - b) / modulus;
            let rows: Vec<Vec<BigRational>> = column
                .iter()
                .map(|c| {
                    let mut row = vec![BigRational::zero(); y_power + 1];
                    row[y_power] = c.clone();
                    row
                })
                .collect();
            matrix.entries[b][a].add_term((0, 0), Biv::new(rows));
        }
    }
    matrix
}

/// The transport `phi(op)` at modulus `r`: an `r x r` operator matrix in
/// `(n, y, D_y, S)` (see the module documentation).
///
/// `max_terms` bounds the number of stored terms of every intermediate
/// matrix; exceeding it is a [`RecurrenceError::LimitExceeded`].
pub(super) fn transport_to_sections(
    op: &DOp,
    modulus: usize,
    max_terms: usize,
) -> Result<OperatorMatrix, RecurrenceError> {
    let derivative = derivative_matrix(modulus);
    let mut derivative_powers = vec![OperatorMatrix::identity(modulus)];
    let mut out = OperatorMatrix::zero(modulus);
    for ((shift, order), coefficient) in &op.0 {
        while derivative_powers.len() <= *order {
            let next = derivative_powers
                .last()
                .expect("the identity is present")
                .compose(&derivative);
            check_term_budget(next.term_count(), max_terms)?;
            derivative_powers.push(next);
        }
        let term = transport_coefficient(coefficient, modulus)
            .compose(&derivative_powers[*order])
            .times_shift(*shift);
        out = out.add(&term);
        check_term_budget(out.term_count(), max_terms)?;
    }
    Ok(out)
}

pub(super) fn check_term_budget(count: usize, max_terms: usize) -> Result<(), RecurrenceError> {
    if count > max_terms {
        return Err(RecurrenceError::BudgetExhausted(format!(
            "an operator in the Veronese section transport has {count} terms, more than the \
             budget {max_terms}"
        )));
    }
    Ok(())
}

/// Component `residue` of the modulus-`modulus` sections of a coefficient
/// row: `V_(r, j)(P)(y) = sum_k [x^(r k + j)] P y^k`, trimmed (`[0]` for zero).
pub(super) fn section_row(row: &[BigRational], modulus: usize, residue: usize) -> Vec<BigRational> {
    vec_poly::trim(
        &row.iter()
            .skip(residue)
            .step_by(modulus)
            .cloned()
            .collect::<Vec<_>>(),
    )
}

#[cfg(test)]
mod tests {
    use super::super::super::implicit::algebra::{from_theta, ThetaOp};
    use super::*;
    use num_bigint::BigInt;

    fn q(value: i64) -> BigRational {
        BigRational::from_integer(BigInt::from(value))
    }

    /// A deterministic pseudo-random operator with `n`-dependent coefficients.
    fn sample_operator(seed: u64, max_shift: usize, max_order: usize) -> DOp {
        let mut state = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 33) % 7) as i64 - 3
        };
        let mut op = DOp::default();
        for shift in 0..=max_shift {
            for order in 0..=max_order {
                let rows: Vec<Vec<BigRational>> = (0..2)
                    .map(|_| (0..3).map(|_| q(next())).collect())
                    .collect();
                op.add_term((shift, order), Biv::new(rows));
            }
        }
        op
    }

    /// Rows of a sample sequence `P_n(x)` with index-dependent coefficients.
    fn sample_rows(count: usize) -> Vec<Vec<BigRational>> {
        (0..count)
            .map(|n| {
                let n = n as i64;
                (0..(n as usize + 4))
                    .map(|k| {
                        let k = k as i64;
                        q((k * k + 3 * n + 1) * if (n + k) % 3 == 0 { -1 } else { 1 })
                    })
                    .collect()
            })
            .collect()
    }

    #[test]
    fn derivative_and_multiplication_matrices_commute_to_identity() {
        for modulus in 2..=5 {
            let x = power_of_multiplication_matrix(modulus, 1);
            let derivative = derivative_matrix(modulus);
            let commutator = derivative.compose(&x).sub(&x.compose(&derivative));
            assert_eq!(
                commutator,
                OperatorMatrix::identity(modulus),
                "r = {modulus}"
            );
            // X^r = y I.
            let mut power = OperatorMatrix::identity(modulus);
            for _ in 0..modulus {
                power = power.compose(&x);
            }
            assert_eq!(power, power_of_multiplication_matrix(modulus, modulus));
            for b in 0..modulus {
                for a in 0..modulus {
                    let expected = if a == b {
                        DOp::single(0, 0, Biv::monomial(0, 1, q(1)))
                    } else {
                        DOp::default()
                    };
                    assert_eq!(power.entries[b][a], expected);
                }
            }
        }
    }

    #[test]
    fn transport_is_a_homomorphism() {
        for modulus in 2..=4 {
            for seed in 0..3u64 {
                let left = sample_operator(seed, 1, 2);
                let right = sample_operator(seed + 17, 1, 1);
                let product = transport_to_sections(&left.compose(&right), modulus, 1 << 20)
                    .expect("within budget");
                let left_matrix = transport_to_sections(&left, modulus, 1 << 20).unwrap();
                let right_matrix = transport_to_sections(&right, modulus, 1 << 20).unwrap();
                assert_eq!(product, left_matrix.compose(&right_matrix), "r = {modulus}");
            }
        }
    }

    #[test]
    fn transport_commutes_with_taking_sections() {
        let rows = sample_rows(9);
        for modulus in 2..=4 {
            let sections: Vec<Vec<Vec<BigRational>>> = (0..modulus)
                .map(|a| {
                    rows.iter()
                        .map(|row| section_row(row, modulus, a))
                        .collect()
                })
                .collect();
            for seed in 0..3u64 {
                let op = sample_operator(seed + 5, 2, 2);
                let matrix = transport_to_sections(&op, modulus, 1 << 20).unwrap();
                for n in 0..5 {
                    let direct = op.apply(&rows, 0, n).expect("rows present");
                    let transported = matrix.apply(&sections, 0, n).expect("rows present");
                    for (b, component) in transported.iter().enumerate() {
                        assert_eq!(
                            &section_row(&direct, modulus, b),
                            component,
                            "r = {modulus}, n = {n}, component {b}"
                        );
                    }
                }
            }
        }
    }

    /// The Euler-form description: `x^s a(n, theta_x)` maps component `b`
    /// to `(b + s) mod r` with factor `y^floor((b+s)/r) a(n, r theta_y + b)`.
    #[test]
    fn transport_agrees_with_the_euler_form() {
        for modulus in 2..=4 {
            for s in 0..=4usize {
                // a(n, theta) = (1 + 2n) theta^2 - theta + 3.
                let a = Biv::new(vec![vec![q(3), q(-1), q(1)], vec![q(0), q(0), q(2)]]);
                let mut euler = ThetaOp::new();
                euler.insert((0, s), a.clone());
                let op = from_theta(&euler);
                let matrix = transport_to_sections(&op, modulus, 1 << 20).unwrap();
                for b in 0..modulus {
                    // a(n, r theta_y + b) in Euler form in y, times y^floor((b+s)/r).
                    let scaled_rows: Vec<Vec<BigRational>> = a
                        .rows()
                        .iter()
                        .map(|row| vec_poly::compose(row, &[q_usize(b), q_usize(modulus)]))
                        .collect();
                    let mut component_euler = ThetaOp::new();
                    component_euler.insert((0, (b + s) / modulus), Biv::new(scaled_rows));
                    let expected = from_theta(&component_euler);
                    let target = (b + s) % modulus;
                    for row in 0..modulus {
                        let entry = &matrix.entries[row][b];
                        if row == target {
                            assert_eq!(entry, &expected, "r = {modulus}, s = {s}, b = {b}");
                        } else {
                            assert!(entry.is_zero());
                        }
                    }
                }
            }
        }
    }
}
