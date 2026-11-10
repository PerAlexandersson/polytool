//! Additional annihilators of the rows of a recurrence-defined sequence,
//! proved invariant under the recurrence.
//!
//! The annihilating ideal of one particular sequence `F` can be larger than
//! the left ideal of its recurrence `L`: the binomial rows `(1 + x)^n`
//! satisfy `F_(n+1) = (1 + x) F_n` and also `(1 + x) F_n' = n F_n`, which no
//! left multiple of `L` gives (it depends on the initial row), and the rows
//! of a three-term recurrence often satisfy a first-order relation with a
//! derivative (a contiguity relation).  This module handles annihilators `A`
//! of order `k_A` below the order `k` of `L`, with `(A F)_n = 0` (base `n`)
//! from some index `n_A` on: **row annihilators** (`k_A = 0`, shift-free
//! `A = sum_d a_d(n, x) D^d`) and **lower-order annihilators**
//! (`0 < k_A < k`, with a leader `c(n) S^(k_A)` free of `x` and of
//! derivatives, so runnable rules of lower order).  Several of them can
//! form a **system** that is invariant only jointly.
//!
//! # Invariance
//!
//! Let `L` have the leader `g(n) S^k` (free of `x` and of derivatives), and
//! let `C_1, ..., C_r` be members of these two kinds.
//! [`AnnihilatorSystemCertificate`] records, for a step `j >= 1`, one
//! identity per member
//!
//! ```text
//! W_i(n, x) (S^j o C_i) = Y_i o L + sum_l Z_(i,l) o C_l,      Z_(i,l) of shift < j,
//! ```
//!
//! found by reducing `S^j o C_i` modulo `L` (right pseudo-division), then
//! the members of positive order (right pseudo-division), then the
//! shift-free members (derivative pseudo-division, column by column,
//! [`reduce_by_row_annihilator`]).  Applied to `F` at a base `n` from which
//! `L` holds, it gives
//! `W_i(n, x) (C_i F)_(n+j) = sum_l sum_(s<j) z(n, x) D^d (C_l F)_(n+s)`, so
//! if every `(C_l F)_m = 0` at the bases `m = n_C, ..., max(n_C, b_L) + j - 1`
//! (`b_L` the first base of `L`), then every member holds at every later
//! base by a joint induction, except that at a base `n` where `W_i(n, .)`
//! is the zero polynomial `(C_i F)_(n+j) = 0` must be checked directly
//! (finitely many: the integer roots of the content of `W_i` in `n`).
//! These row checks are the caller's; [`AnnihilatorSystemCertificate::verify`]
//! recomputes the identities.  A lone row annihilator uses `j = k`; a
//! runnable member `B` of order `k - 1`, alone or with a row annihilator
//! `A`, uses `j = 1` (for a three-term recurrence the contiguity relation
//! `B` and the differential equation `A` are invariant only together).
//!
//! # Derivative pseudo-division
//!
//! A shift-free `c` is divided by a shift-free `A` with leader
//! `a_r(n, x) D^r` level by level from a fixed top derivative order down to
//! `r`: multiply the remainder and the quotient on the left by `a_r` and
//! subtract `c_top D^(m - r) o A`, whose top term is `c_top a_r D^m`
//! (`D^(m-r) o a_r = a_r D^(m-r) + (lower orders)`).  Every level multiplies,
//! so the remainder is linear in `c` for a fixed top order.

use crate::bivariate::BivariatePolynomial;
use crate::division::right_pseudo_division;
use crate::kernel::{exact_kernel_vectors_by_group, integer_equation_row};
use crate::operator::OreOperator;
use crate::Rational;
use num_bigint::BigInt;
use num_traits::Zero;
use std::collections::BTreeMap;
use std::fmt;

/// The result `W c = Z o A + R` of [`reduce_by_row_annihilator`] for an
/// operator `c` of shift below `columns`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RowAnnihilatorReduction {
    /// `W(n, x)`, nonzero.
    pub multiplier: BivariatePolynomial,
    /// `Z`, of shift below `columns`.
    pub quotient: OreOperator,
    /// `R`, every column of derivative order below that of `A`.
    pub remainder: OreOperator,
}

/// `(multiplier, quotient, remainder)` of the derivative pseudo-division of
/// a shift-free `c` by a shift-free `divisor` with every level from
/// `derivative_top` down (module documentation).
fn derivative_pseudo_division(
    c: &OreOperator,
    divisor: &OreOperator,
    derivative_top: usize,
) -> (BivariatePolynomial, OreOperator, OreOperator) {
    let order = divisor.max_derivative();
    let lead = divisor.0.get(&(0, order)).cloned().unwrap_or_default();
    let mut remainder = c.clone();
    let mut quotient = OreOperator::default();
    let mut multiplier = BivariatePolynomial::one();
    if derivative_top < order {
        return (multiplier, quotient, remainder);
    }
    for level in (order..=derivative_top).rev() {
        let top = remainder.0.get(&(0, level)).cloned();
        remainder = remainder.left_multiply(&lead);
        quotient = quotient.left_multiply(&lead);
        multiplier = multiplier.multiply(&lead);
        let Some(top) = top else {
            continue;
        };
        let step = OreOperator::single(0, level - order, top);
        remainder.subtract_assign(&step.compose(divisor));
        quotient.add_assign(&step);
    }
    (multiplier, quotient, remainder)
}

/// Reduce an operator `c` of shift below `columns` by the row annihilator
/// `A` (module documentation): column `i` is divided by `A(n + i)` with
/// every derivative level from `derivative_top` down, and all columns are
/// brought to the common multiplier `W = prod_i a_r(n + i, x)^(levels)`.
/// `None` when `A` is zero or has a shift, or `c` has a shift `>= columns`
/// or a derivative order above `derivative_top`.
pub fn reduce_by_row_annihilator(
    c: &OreOperator,
    annihilator: &OreOperator,
    columns: usize,
    derivative_top: usize,
) -> Option<RowAnnihilatorReduction> {
    if annihilator.is_zero() || annihilator.max_shift() != 0 {
        return None;
    }
    if !c.is_zero() && (c.max_shift() >= columns || c.max_derivative() > derivative_top) {
        return None;
    }
    let parts: Vec<_> = (0..columns)
        .map(|i| {
            derivative_pseudo_division(
                &c.column(i),
                &annihilator.shift_index(i as i64),
                derivative_top,
            )
        })
        .collect();
    let multiplier = parts
        .iter()
        .fold(BivariatePolynomial::one(), |acc, (m, _, _)| acc.multiply(m));
    let mut quotient = OreOperator::default();
    let mut remainder = OreOperator::default();
    for (i, (_, q, r)) in parts.iter().enumerate() {
        let others = parts
            .iter()
            .enumerate()
            .filter(|(j, _)| *j != i)
            .fold(BivariatePolynomial::one(), |acc, (_, (m, _, _))| {
                acc.multiply(m)
            });
        for ((_, d), coefficient) in &q.0 {
            quotient.add_term((i, *d), coefficient.multiply(&others));
        }
        for ((_, d), coefficient) in &r.0 {
            remainder.add_term((i, *d), coefficient.multiply(&others));
        }
    }
    Some(RowAnnihilatorReduction {
        multiplier,
        quotient,
        remainder,
    })
}

/// Why an annihilator-system certificate does not replay.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AnnihilatorSystemError {
    /// The recurrence leader is not free of `x` and of derivatives.
    RecurrenceLeaderNotIndexOnly,
    /// A member is zero, of order at least the recurrence's, or of positive
    /// order without a leader free of `x` and of derivatives; or the system
    /// is empty, or the step is zero.
    InvalidMember,
    /// The entries do not match the members.
    EntryCount,
    /// A multiplier is zero.
    ZeroMultiplier,
    /// A member quotient has a shift `>=` the step.
    QuotientShiftTooLarge,
    /// An identity does not hold.
    IdentityFails,
}

impl fmt::Display for AnnihilatorSystemError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::RecurrenceLeaderNotIndexOnly => {
                "the recurrence leader must be g(n) S^k, free of x and of derivatives"
            }
            Self::InvalidMember => {
                "the members must be nonzero, of order below the recurrence's, shift-free or \
                 with a leader free of x and of derivatives, and the step positive"
            }
            Self::EntryCount => "there must be one identity per member",
            Self::ZeroMultiplier => "a multiplier is zero",
            Self::QuotientShiftTooLarge => "a member quotient must have shift below the step",
            Self::IdentityFails => {
                "an identity W (S^j o C_i) = Y o L + sum_l Z_l o C_l does not hold"
            }
        })
    }
}

impl std::error::Error for AnnihilatorSystemError {}

/// The invariance identity of one member `C_i`:
/// `W_i (S^j o C_i) = Y_i o L + sum_l Z_(i,l) o C_l`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnihilatorSystemEntry {
    /// `W_i(n, x)`, nonzero.
    pub multiplier: BivariatePolynomial,
    /// `Y_i`.
    pub recurrence_quotient: OreOperator,
    /// `Z_(i,l)`, one per member, each of shift below the step `j`.
    pub member_quotients: Vec<OreOperator>,
}

/// A jointly invariant system of additional annihilators `C_1, ..., C_r`
/// of a recurrence `L` (module documentation): one identity per member.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AnnihilatorSystemCertificate {
    /// `L`, with a leader `g(n) S^k` free of `x` and of derivatives.
    pub recurrence: OreOperator,
    /// The members: shift-free, or of order below `k` with a leader free
    /// of `x` and of derivatives.
    pub members: Vec<OreOperator>,
    /// The step `j >= 1`.
    pub step: usize,
    /// The identities, one per member.
    pub entries: Vec<AnnihilatorSystemEntry>,
}

/// Whether `member` is admissible for a recurrence of order `order`.
fn admissible(member: &OreOperator, order: usize) -> bool {
    if member.is_zero() {
        return false;
    }
    let k_c = member.max_shift();
    k_c < order.max(1) && (k_c == 0 || crate::division::index_only_leader(member).is_some())
}

impl AnnihilatorSystemCertificate {
    /// Recompute every identity exactly.  The row checks (each `(C_i F)_m = 0`
    /// at the bases `m = n_C, ..., max(n_C, b_L) + j - 1` and at `n + j` for
    /// the bases `n` where `W_i(n, .)` vanishes, with `b_L` the first base of
    /// `L`) are the caller's.
    pub fn verify(&self) -> Result<(), AnnihilatorSystemError> {
        let (order, _) = crate::division::index_only_leader(&self.recurrence)
            .ok_or(AnnihilatorSystemError::RecurrenceLeaderNotIndexOnly)?;
        if self.members.is_empty()
            || self.step == 0
            || !self.members.iter().all(|c| admissible(c, order))
        {
            return Err(AnnihilatorSystemError::InvalidMember);
        }
        if self.entries.len() != self.members.len()
            || self
                .entries
                .iter()
                .any(|e| e.member_quotients.len() != self.members.len())
        {
            return Err(AnnihilatorSystemError::EntryCount);
        }
        for (member, entry) in self.members.iter().zip(&self.entries) {
            if entry.multiplier.is_zero() {
                return Err(AnnihilatorSystemError::ZeroMultiplier);
            }
            if entry
                .member_quotients
                .iter()
                .any(|z| !z.is_zero() && z.max_shift() >= self.step)
            {
                return Err(AnnihilatorSystemError::QuotientShiftTooLarge);
            }
            let left = member.shifted(self.step).left_multiply(&entry.multiplier);
            let mut right = entry.recurrence_quotient.compose(&self.recurrence);
            for (z, c) in entry.member_quotients.iter().zip(&self.members) {
                right.add_assign(&z.compose(c));
            }
            if left != right {
                return Err(AnnihilatorSystemError::IdentityFails);
            }
        }
        Ok(())
    }

    /// The member of positive order, if any: a runnable rule of lower order
    /// than the recurrence.
    pub fn lower_order_member(&self) -> Option<&OreOperator> {
        self.members.iter().find(|c| c.max_shift() > 0)
    }

    /// The shift-free member, if any.
    pub fn row_annihilator(&self) -> Option<&OreOperator> {
        self.members.iter().find(|c| c.max_shift() == 0)
    }
}

/// Reduce `op` modulo `L`, then the members of positive order (right
/// pseudo-division), then the shift-free members (derivative pseudo-division
/// column by column): `W op = Y o L + sum_l Z_l o C_l + R`.
#[allow(clippy::type_complexity)]
fn reduce_modulo_system(
    op: &OreOperator,
    recurrence: &OreOperator,
    members: &[OreOperator],
) -> Option<(
    BivariatePolynomial,
    OreOperator,
    Vec<OreOperator>,
    OreOperator,
)> {
    let mut multiplier = BivariatePolynomial::one();
    let mut recurrence_quotient = OreOperator::default();
    let mut quotients = vec![OreOperator::default(); members.len()];
    let mut remainder = op.clone();
    let apply = |factor: &BivariatePolynomial,
                 multiplier: &mut BivariatePolynomial,
                 recurrence_quotient: &mut OreOperator,
                 quotients: &mut Vec<OreOperator>| {
        *multiplier = multiplier.multiply(factor);
        *recurrence_quotient = recurrence_quotient.left_multiply(factor);
        for q in quotients.iter_mut() {
            *q = q.left_multiply(factor);
        }
    };
    let top = remainder.max_shift();
    let division = right_pseudo_division(&remainder, recurrence, top)?;
    apply(
        &division.multiplier,
        &mut multiplier,
        &mut recurrence_quotient,
        &mut quotients,
    );
    recurrence_quotient.add_assign(&division.quotient);
    remainder = division.remainder;
    for (l, member) in members.iter().enumerate() {
        if member.max_shift() == 0 || remainder.is_zero() {
            continue;
        }
        let division = right_pseudo_division(&remainder, member, remainder.max_shift())?;
        apply(
            &division.multiplier,
            &mut multiplier,
            &mut recurrence_quotient,
            &mut quotients,
        );
        quotients[l].add_assign(&division.quotient);
        remainder = division.remainder;
    }
    for (l, member) in members.iter().enumerate() {
        if member.max_shift() != 0 || remainder.is_zero() {
            continue;
        }
        let columns = remainder.max_shift() + 1;
        let reduction =
            reduce_by_row_annihilator(&remainder, member, columns, remainder.max_derivative())?;
        apply(
            &reduction.multiplier,
            &mut multiplier,
            &mut recurrence_quotient,
            &mut quotients,
        );
        quotients[l].add_assign(&reduction.quotient);
        remainder = reduction.remainder;
    }
    Some((multiplier, recurrence_quotient, quotients, remainder))
}

/// Prove the members jointly invariant with step `step` (module
/// documentation): `None` when a member is not admissible or a reduction
/// leaves a remainder or a quotient of too large a shift.
pub fn certify_annihilator_system(
    recurrence: &OreOperator,
    members: &[OreOperator],
    step: usize,
) -> Option<AnnihilatorSystemCertificate> {
    let (order, _) = crate::division::index_only_leader(recurrence)?;
    if step == 0 || members.is_empty() || !members.iter().all(|c| admissible(c, order)) {
        return None;
    }
    let mut entries = Vec::with_capacity(members.len());
    for member in members {
        let (multiplier, recurrence_quotient, member_quotients, remainder) =
            reduce_modulo_system(&member.shifted(step), recurrence, members)?;
        if !remainder.is_zero() {
            return None;
        }
        entries.push(AnnihilatorSystemEntry {
            multiplier,
            recurrence_quotient,
            member_quotients,
        });
    }
    let certificate = AnnihilatorSystemCertificate {
        recurrence: recurrence.clone(),
        members: members.to_vec(),
        step,
        entries,
    };
    certificate.verify().is_ok().then_some(certificate)
}

/// The budget of [`find_annihilator_system`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AnnihilatorSystemBudget {
    /// Largest derivative order of a member.  Default 3.
    pub max_derivative_order: usize,
    /// Largest coefficient degree in `n` and in `x`.  Default 4.
    pub max_degree: usize,
}

impl Default for AnnihilatorSystemBudget {
    fn default() -> Self {
        Self {
            max_derivative_order: 3,
            max_degree: 4,
        }
    }
}

/// Fit an annihilator of order `order` (top column `c(n) S^order` without
/// `x` and derivatives when `order > 0`; for `order = 0` the top derivative
/// column is the pivot) with derivative order `derivative_order` and
/// coefficient degree `degree` to the rows, exactly; `None` when the rows
/// admit none.
fn fit_annihilator(
    rows: &[Vec<Rational>],
    first_index: usize,
    order: usize,
    derivative_order: usize,
    degree: usize,
) -> Option<OreOperator> {
    // (shift, derivative, n power, x power).
    let mut unknowns = Vec::new();
    let mut groups = Vec::new();
    for i in 0..=order {
        for d in 0..=derivative_order {
            if order > 0 && i == order && d > 0 {
                continue;
            }
            for a in 0..=degree {
                for b in 0..=degree {
                    if order > 0 && i == order && b > 0 {
                        continue;
                    }
                    unknowns.push((i, d, a, b));
                    let top = if order > 0 {
                        i == order
                    } else {
                        d == derivative_order
                    };
                    groups.push(usize::from(top));
                }
            }
        }
    }
    let bases = rows.len().checked_sub(order)?;
    let mut derivatives: Vec<Vec<Vec<Rational>>> = Vec::with_capacity(rows.len());
    for row in rows {
        let mut list = vec![row.clone()];
        for _ in 0..derivative_order {
            let next = crate::univariate::derivative(list.last().expect("nonempty"));
            list.push(next);
        }
        derivatives.push(list);
    }
    let mut equations: BTreeMap<(usize, usize), BTreeMap<usize, Rational>> = BTreeMap::new();
    for base in 0..bases {
        // The powers n^a, a = 0..=degree, of the integer base n.
        let n = BigInt::from(first_index + base);
        let weights: Vec<Rational> = (0..=degree)
            .map(|a| Rational::from_integer(num_traits::pow(n.clone(), a)))
            .collect();
        for (column, (i, d, a, b)) in unknowns.iter().enumerate() {
            let weight = &weights[*a];
            if weight.is_zero() {
                continue;
            }
            for (k, value) in derivatives[base + i][*d].iter().enumerate() {
                if value.is_zero() {
                    continue;
                }
                let slot = equations
                    .entry((base, k + b))
                    .or_default()
                    .entry(column)
                    .or_insert_with(Rational::zero);
                crate::univariate::add_product_into(slot, value, weight);
            }
        }
    }
    let matrix: Vec<Vec<BigInt>> = equations
        .values()
        .map(|entries| integer_equation_row(entries, unknowns.len()))
        .collect();
    let vector = exact_kernel_vectors_by_group(&matrix, &groups, 2)[1].clone()?;
    let mut candidate = OreOperator::default();
    for ((i, d, a, b), value) in unknowns.iter().zip(&vector) {
        if !value.is_zero() {
            candidate.add_term(
                (*i, *d),
                BivariatePolynomial::monomial(*a, *b, value.clone()),
            );
        }
    }
    Some(crate::image_annihilator::normalize_annihilator(&candidate))
}

/// The first annihilator of order `order` that the rows admit, over
/// derivative orders and degrees in increasing order (derivative order from
/// one for a shift-free annihilator), fitted exactly to the rows
/// (`rows[j]` is `F_(first_index + j)`).
///
/// For `order > 0` the top column is `c(n) S^order` (free of `x` and of
/// derivatives, a runnable rule of that order); for `order = 0` it is a row
/// annihilator `sum_d a_d(n, x) D^d` with a derivative.  Coefficients have
/// degree at most `budget.max_degree` in `n` and in `x`.  A fit is
/// evidence only: it holds on the supplied rows, and a caller that uses it
/// must prove it (for example by an invariance identity and exact checks of
/// the base rows).
pub fn fit_rows_annihilator(
    rows: &[Vec<Rational>],
    first_index: usize,
    order: usize,
    budget: &AnnihilatorSystemBudget,
) -> Option<OreOperator> {
    first_fit(rows, first_index, order, budget)
}

fn first_fit(
    rows: &[Vec<Rational>],
    first_index: usize,
    order: usize,
    budget: &AnnihilatorSystemBudget,
) -> Option<OreOperator> {
    let lowest = usize::from(order == 0);
    for derivative_order in lowest..=budget.max_derivative_order {
        for degree in 0..=budget.max_degree {
            // Enough distinct indices for the degree in n (spurious fits are
            // rejected by the invariance proof anyway).
            if rows.len() < degree + derivative_order + order + 6 {
                continue;
            }
            if let Some(candidate) =
                fit_annihilator(rows, first_index, order, derivative_order, degree)
            {
                return Some(candidate);
            }
        }
    }
    None
}

/// A jointly invariant system of additional annihilators of the sequence
/// with rows `rows` (`rows[j]` is `F_(first_index + j)`), annihilated by
/// `recurrence` (module documentation).  Candidates are fitted exactly to
/// the rows: a runnable rule `B` of order `k - 1` with a derivative (for
/// `k >= 2`) and a shift-free row annihilator `A`; the systems `{B}`,
/// `{B, A}` (step 1) and `{A}` (step `k`) are tried in this order and
/// returned only with their identities.  The fit uses only the supplied
/// rows; the caller still checks the base rows.  `None` when nothing is
/// found within the budget.
pub fn find_annihilator_system(
    recurrence: &OreOperator,
    rows: &[Vec<Rational>],
    first_index: usize,
    budget: &AnnihilatorSystemBudget,
) -> Option<AnnihilatorSystemCertificate> {
    let (order, _) = crate::division::index_only_leader(recurrence)?;
    let lower = if order >= 2 {
        first_fit(rows, first_index, order - 1, budget)
    } else {
        None
    };
    let row = first_fit(rows, first_index, 0, budget);
    let mut systems: Vec<(Vec<OreOperator>, usize)> = Vec::new();
    if let Some(b) = &lower {
        systems.push((vec![b.clone()], 1));
        if let Some(a) = &row {
            systems.push((vec![b.clone(), a.clone()], 1));
        }
    }
    if let Some(a) = &row {
        systems.push((vec![a.clone()], order));
    }
    systems
        .into_iter()
        .find_map(|(members, step)| certify_annihilator_system(recurrence, &members, step))
}
