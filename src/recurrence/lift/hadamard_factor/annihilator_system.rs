//! Invariant row-annihilator systems of parents of any order `k`: operators
//! `A^(l)(n) = sum_(i < k) A^(l)_i(n) S^i` that kill every parent row
//! window `(P_n, ..., P_(n+k-1))`, certified by exact invariance identities
//! under the parent's companion step (module documentation of
//! `hadamard_factor`, section "Row annihilators").
//!
//! The candidates are guessed from the parent rows: for a growing ansatz
//! (derivative order, degrees in `n` and in `x`) the exact kernel of the
//! linear conditions `(A P)_n = 0` on the rows from the start base on gives
//! the operators of that size that kill the checked rows; with the unknowns
//! ordered by shift, then derivative order, the kernel basis is in echelon
//! form, and for every shift `c` the kernel vector with the smallest leading
//! term among those led by shift `c` is taken as the pivot of column `c`
//! (an echelon basis of the annihilator module over the rational Weyl
//! algebra `Q(n, x)<D>` when the ansatz is large enough).  The step images
//! `N(A^(l))` and the division remainder are then reduced by the pivots
//! (right pseudo-division, column by column from the top); the system is
//! kept when every reduction ends in zero.  Nothing in the proof depends on
//! how the candidates were found: the replay checks the identities and the
//! base rows only.

use super::annihilator::{shift_coefficients, vanishing_bases};
use super::certify::right_pseudo_divide;
use super::modular_kernel::kernel_vectors_by_group;
use super::prepare::{left_multiply, PreparedRule};
use crate::recurrence::implicit::algebra::{Biv, DOp};
use crate::recurrence::lift::section_certify::x_free_leader;
use crate::recurrence::lift::section_operators::check_term_budget;
use crate::recurrence::RecurrenceError;
use crate::vec_poly;
use crate::BigRational;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};

type Q = BigRational;

/// The ansatz sizes `(derivative order, degree in n, degree in x)` tried in
/// turn, smallest first.
const SYSTEM_ANSATZ_SIZES: [(usize, usize, usize); 6] = [
    (1, 1, 2),
    (1, 2, 2),
    (2, 2, 3),
    (2, 3, 4),
    (3, 3, 4),
    (3, 3, 6),
];

/// Largest number of unknowns of one ansatz.
const MAX_SYSTEM_UNKNOWNS: usize = 256;

/// Largest number of parent rows used by one ansatz.
const MAX_SYSTEM_ROWS: usize = 48;

/// The step image
///
/// ```text
/// N(A) = g(n) S o A - A_(k-1)(n + 1) o L
/// ```
///
/// of `A = sum_(i < k) A_i(n) S^i` under the parent relation
/// `L = g(n) S^k + sum_(i < k) C_i(n) S^i` (leader free of `x` and `D`).  It
/// has shifts below `k` again, and `(N(A) P)_n = g(n) (A P)_(n+1)` at every
/// base `n` where `(L P)_n = 0`.
pub(super) fn step_image(a: &DOp, parent: &DOp) -> Result<DOp, RecurrenceError> {
    let (k, leader) = x_free_leader(parent)?;
    if !a.is_zero() && a.max_shift() >= k {
        return Err(RecurrenceError::VerificationFailed(
            "a row annihilator must have shifts below the parent order".into(),
        ));
    }
    let top = shift_coefficients(&a.column(k - 1), 1);
    let image = left_multiply(&leader, &a.shifted(1)).sub(&top.compose(parent));
    if !image.is_zero() && image.max_shift() >= k {
        return Err(RecurrenceError::VerificationFailed(
            "the step image did not cancel the parent leader".into(),
        ));
    }
    Ok(image)
}

/// An invariant system: generators `A^(l)` with
/// `nu_l(n, x) N(A^(l)) = sum_j B_(l,j)(n) o A^(j)` exactly.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct InvariantAnnihilatorSystem {
    /// `A^(l)`, shifts below the parent order.
    pub(super) generators: Vec<DOp>,
    /// `nu_l(n, x)`, nonzero.
    pub(super) step_multipliers: Vec<Biv>,
    /// `B_(l,j)`, shift-free; `step_quotients[l][j]`.
    pub(super) step_quotients: Vec<Vec<DOp>>,
}

impl InvariantAnnihilatorSystem {
    /// The integer bases `m >= from` where some `nu_l(m, .)` is the zero
    /// polynomial, ascending and without repetition.
    pub(super) fn restart_bases(&self, from: usize) -> Vec<usize> {
        let mut bases: Vec<usize> = self
            .step_multipliers
            .iter()
            .flat_map(|nu| vanishing_bases(nu, from))
            .collect();
        bases.sort_unstable();
        bases.dedup();
        bases
    }
}

/// Every invariance identity `nu_l N(A^(l)) = sum_j B_(l,j) o A^(j)` holds
/// exactly, with every `nu_l` nonzero, every `B_(l,j)` shift-free and every
/// generator nonzero with shifts below the parent order.
pub(super) fn system_invariance_holds(
    system: &InvariantAnnihilatorSystem,
    parent: &DOp,
) -> Result<bool, RecurrenceError> {
    let (k, _) = x_free_leader(parent)?;
    let count = system.generators.len();
    if count == 0 || system.step_multipliers.len() != count || system.step_quotients.len() != count
    {
        return Ok(false);
    }
    if system
        .generators
        .iter()
        .any(|a| a.is_zero() || a.max_shift() >= k)
    {
        return Ok(false);
    }
    for ((generator, nu), quotients) in system
        .generators
        .iter()
        .zip(&system.step_multipliers)
        .zip(&system.step_quotients)
    {
        if nu.is_zero()
            || quotients.len() != count
            || quotients.iter().any(|b| !b.is_zero() && b.max_shift() != 0)
        {
            return Ok(false);
        }
        let mut difference = left_multiply(nu, &step_image(generator, parent)?);
        for (b, a) in quotients.iter().zip(&system.generators) {
            difference = difference.sub(&b.compose(a));
        }
        if !difference.is_zero() {
            return Ok(false);
        }
    }
    Ok(true)
}

/// `multiplier element = sum_c quotients[c] o pivot_c + remainder`.
struct Reduction {
    multiplier: Biv,
    /// By pivot column.
    quotients: Vec<DOp>,
    remainder: DOp,
}

/// An echelon basis: `pivots[c]` has zero columns above `c` and a nonzero
/// column `c`.
struct EchelonBasis {
    pivots: Vec<Option<DOp>>,
    max_terms: usize,
}

impl EchelonBasis {
    /// Reduce the columns from the top down by right pseudo-division until a
    /// column without a pivot, or with a remainder, is reached.
    fn reduce(&self, element: &DOp) -> Result<Reduction, RecurrenceError> {
        let order = self.pivots.len();
        let mut remainder = element.clone();
        let mut multiplier = Biv::one();
        let mut quotients = vec![DOp::default(); order];
        for c in (0..order).rev() {
            let column = remainder.column(c);
            if column.is_zero() {
                continue;
            }
            let Some(pivot) = &self.pivots[c] else {
                break;
            };
            let (m, z, _) = right_pseudo_divide(&column, &pivot.column(c), self.max_terms)?;
            remainder = left_multiply(&m, &remainder).sub(&z.compose(pivot));
            for quotient in &mut quotients {
                *quotient = left_multiply(&m, quotient);
            }
            quotients[c] = quotients[c].add(&z);
            multiplier = multiplier.mul(&m);
            check_term_budget(
                remainder.0.len() + quotients.iter().map(|q| q.0.len()).sum::<usize>(),
                self.max_terms,
            )?;
            if !remainder.column(c).is_zero() {
                break;
            }
        }
        Ok(Reduction {
            multiplier,
            quotients,
            remainder,
        })
    }

    /// The pivots, top column first, with their columns.
    fn generators(&self) -> Vec<(usize, DOp)> {
        (0..self.pivots.len())
            .rev()
            .filter_map(|c| self.pivots[c].clone().map(|p| (c, p)))
            .collect()
    }
}

/// `mu r = sum_l Z_l o A^(l)`: an invariant system with the decomposition of
/// the remainder.
pub(super) struct RemainderDecomposition {
    pub(super) system: InvariantAnnihilatorSystem,
    /// `mu(n, x)`.
    pub(super) multiplier: Biv,
    /// `Z_l`, shift-free, one per generator.
    pub(super) quotients: Vec<DOp>,
}

/// The pivots guessed from the parent rows with the ansatz
/// `sum c_(i,d,a,b) n^a x^b D^d S^i` (`i < k`, `d <= order`, `a <= n_degree`,
/// `b <= x_degree`), one per shift where the kernel has a vector led by
/// that shift; `None` when the ansatz has too many unknowns or the rows do
/// not give enough equations.
fn guess_pivots(
    parent: &PreparedRule,
    k: usize,
    from: usize,
    (order, n_degree, x_degree): (usize, usize, usize),
) -> Result<Option<Vec<Option<DOp>>>, RecurrenceError> {
    let mut unknowns = Vec::new();
    for i in 0..k {
        for d in 0..=order {
            for a in 0..=n_degree {
                for b in 0..=x_degree {
                    unknowns.push((i, d, a, b));
                }
            }
        }
    }
    if unknowns.len() > MAX_SYSTEM_UNKNOWNS {
        return Ok(None);
    }
    let first = parent.first_index;
    let available = MAX_SYSTEM_ROWS.min(super::search::MAX_HADAMARD_FACTOR_ROWS);
    let rows = parent.rows(available)?;
    // Derivatives of every row up to `order`.
    let derivatives: Vec<Vec<Vec<Q>>> = rows
        .iter()
        .map(|row| {
            let mut out = vec![row.clone()];
            for _ in 0..order {
                let next = vec_poly::derivative(out.last().map_or(&[][..], Vec::as_slice));
                out.push(vec_poly::trim(&next));
            }
            out
        })
        .collect();
    let target = unknowns.len() + unknowns.len() / 2 + 8;
    let mut equations: Vec<Vec<Q>> = Vec::new();
    let mut base = from.max(first);
    while equations.len() < target {
        if base + k > first + rows.len() {
            return Ok(None);
        }
        let powers: Vec<Q> = (0..=n_degree)
            .map(|a| Q::from_integer(base.into()).pow(a as i32))
            .collect();
        let width = (0..k)
            .map(|i| derivatives[base + i - first][0].len())
            .max()
            .unwrap_or(0)
            + x_degree;
        for j in 0..width {
            let equation: Vec<Q> = unknowns
                .iter()
                .map(|&(i, d, a, b)| {
                    j.checked_sub(b)
                        .and_then(|t| derivatives[base + i - first][d].get(t))
                        .map_or_else(Q::zero, |c| c * &powers[a])
                })
                .collect();
            if equation.iter().any(|c| !c.is_zero()) {
                equations.push(equation);
            }
        }
        base += 1;
    }
    // Integer equations (denominators cleared row by row).
    let integral: Vec<Vec<BigInt>> = equations
        .iter()
        .map(|equation| {
            let denominator = equation
                .iter()
                .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
            equation
                .iter()
                .map(|value| value.numer() * (&denominator / value.denom()))
                .collect()
        })
        .collect();
    // Unknowns ascend by (shift, derivative order, ...): the vector led by
    // the smallest free column of a shift has the smallest leading term
    // among the kernel vectors led by that shift.
    let shifts: Vec<usize> = unknowns.iter().map(|&(i, ..)| i).collect();
    let vectors = kernel_vectors_by_group(&integral, &shifts, k);
    Ok(Some(
        vectors
            .into_iter()
            .map(|vector| {
                vector.map(|vector| {
                    // Integer coefficients without a common factor.
                    let denominator = vector
                        .iter()
                        .fold(BigInt::one(), |acc, value| acc.lcm(value.denom()));
                    let content = vector
                        .iter()
                        .fold(BigInt::zero(), |acc, value| acc.gcd(value.numer()));
                    let scale = Q::new(denominator, content);
                    let mut op = DOp::default();
                    for (value, &(i, d, a, b)) in vector.iter().zip(&unknowns) {
                        if !value.is_zero() {
                            op.add_term((i, d), Biv::monomial(a, b, value * &scale));
                        }
                    }
                    op
                })
            })
            .collect(),
    ))
}

/// Search an invariant annihilator system whose module contains `remainder`
/// (module documentation), with parent rows checked from the base `from`
/// on: `None` when no ansatz size gives one.
pub(super) fn find_invariant_annihilator_system(
    remainder: &DOp,
    parent: &PreparedRule,
    from: usize,
    max_terms: usize,
) -> Result<Option<RemainderDecomposition>, RecurrenceError> {
    let (k, _) = x_free_leader(&parent.relation)?;
    if remainder.is_zero() || remainder.max_shift() >= k {
        return Ok(None);
    }
    for size in SYSTEM_ANSATZ_SIZES {
        let Some(pivots) = guess_pivots(parent, k, from, size)? else {
            continue;
        };
        let basis = EchelonBasis { pivots, max_terms };
        if let Some(found) = decompose(&basis, remainder, &parent.relation)? {
            return Ok(Some(found));
        }
    }
    Ok(None)
}

/// The invariance identities of the pivots and the decomposition of the
/// remainder, when every reduction ends in zero.
fn decompose(
    basis: &EchelonBasis,
    remainder: &DOp,
    parent: &DOp,
) -> Result<Option<RemainderDecomposition>, RecurrenceError> {
    let generators = basis.generators();
    if generators.is_empty() {
        return Ok(None);
    }
    let slot = |column: usize| generators.iter().position(|(c, _)| *c == column);
    let by_generator = |quotients: Vec<DOp>| -> Vec<DOp> {
        let mut out = vec![DOp::default(); generators.len()];
        for (column, quotient) in quotients.into_iter().enumerate() {
            if let Some(index) = slot(column) {
                out[index] = quotient;
            }
        }
        out
    };
    // The remainder first: it is the cheapest rejection.
    let reduction = basis.reduce(remainder)?;
    if !reduction.remainder.is_zero() {
        return Ok(None);
    }
    let (multiplier, quotients) = (reduction.multiplier, by_generator(reduction.quotients));
    let mut step_multipliers = Vec::with_capacity(generators.len());
    let mut step_quotients = Vec::with_capacity(generators.len());
    for (_, generator) in &generators {
        let reduction = basis.reduce(&step_image(generator, parent)?)?;
        if !reduction.remainder.is_zero() {
            return Ok(None);
        }
        step_multipliers.push(reduction.multiplier);
        step_quotients.push(by_generator(reduction.quotients));
    }
    let system = InvariantAnnihilatorSystem {
        generators: generators.into_iter().map(|(_, a)| a).collect(),
        step_multipliers,
        step_quotients,
    };
    if !system_invariance_holds(&system, parent)? {
        return Err(RecurrenceError::VerificationFailed(
            "the annihilator system search produced a non-invariant system".into(),
        ));
    }
    Ok(Some(RemainderDecomposition {
        system,
        multiplier,
        quotients,
    }))
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
    fn an_invariant_system_of_a_second_order_parent() {
        // P_(n+2) = (1 + x)^2 P_n (rows (1 + x)^n): L = S^2 - (1 + x)^2.
        let mut parent = DOp::single(2, 0, Biv::one());
        parent.add_term((0, 0), biv(&[&[-1, -2, -1]]));
        // B = S - (1 + x) and A = (1 + x) D - n.
        let mut top = DOp::single(1, 0, Biv::one());
        top.add_term((0, 0), biv(&[&[-1, -1]]));
        let mut bottom = DOp::single(0, 1, biv(&[&[1, 1]]));
        bottom.add_term((0, 0), biv(&[&[], &[-1]]));
        // N(B) = -(1 + x) B; N(A) = Z B + (1 + x) A with Z = (1 + x) D - (n + 1).
        let image = step_image(&top, &parent).unwrap();
        assert_eq!(image, DOp::single(0, 0, biv(&[&[-1, -1]])).compose(&top));
        let mut z = DOp::single(0, 1, biv(&[&[1, 1]]));
        z.add_term((0, 0), biv(&[&[-1], &[-1]]));
        let system = InvariantAnnihilatorSystem {
            generators: vec![top.clone(), bottom.clone()],
            step_multipliers: vec![Biv::one(), Biv::one()],
            step_quotients: vec![
                vec![DOp::single(0, 0, biv(&[&[-1, -1]])), DOp::default()],
                vec![z.clone(), DOp::single(0, 0, biv(&[&[1, 1]]))],
            ],
        };
        assert!(system_invariance_holds(&system, &parent).unwrap());
        assert!(system.restart_bases(0).is_empty());
        // A wrong quotient, or a generator of shift k, is refused.
        let mut wrong = system.clone();
        wrong.step_quotients[1][0] = DOp::default();
        assert!(!system_invariance_holds(&wrong, &parent).unwrap());
        let mut wrong = system.clone();
        wrong.generators[0] = top.shifted(1);
        assert!(!system_invariance_holds(&wrong, &parent).unwrap());
        // The reduction by the echelon basis recovers the quotients of
        // r = Z B + (1 + x) A = N(A).
        let basis = EchelonBasis {
            pivots: vec![Some(bottom), Some(top)],
            max_terms: 1_000,
        };
        let remainder = step_image(&system.generators[1], &parent).unwrap();
        let decomposition = decompose(&basis, &remainder, &parent)
            .unwrap()
            .expect("a decomposition");
        assert_eq!(decomposition.system.generators.len(), 2);
        let mut sum = DOp::default();
        for (z, a) in decomposition
            .quotients
            .iter()
            .zip(&decomposition.system.generators)
        {
            sum = sum.add(&z.compose(a));
        }
        assert_eq!(left_multiply(&decomposition.multiplier, &remainder), sum);
        // Restart bases: nu = n - 3 vanishes identically at n = 3.
        let mut restarted = system;
        restarted.step_multipliers[0] = biv(&[&[-3], &[1]]);
        assert_eq!(restarted.restart_bases(0), vec![3]);
    }
}
