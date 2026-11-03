//! Runnable scalar annihilators of one component of a vector sequence
//! annihilated by a square system of operators, with identity
//! certificates.
//!
//! Let `V = (V_0, ..., V_(r-1))` be a vector of sequences with
//! `(Phi V)_n = 0` (every row `b`: `sum_a (Phi_(b,a) V_a)_n = 0`) at every
//! base `n` from some `b_0` on, where the shift-`k` part of `Phi` is
//! `g(n) S^k I` with `g` a nonzero polynomial in `n` only and no entry has
//! a larger shift ([`crate::system_index_only_leader`]).  A scalar operator
//! `M` acting on the component `V_c` is certified by the identity
//!
//! ```text
//! w(n) (M e_c) = sum_b H_b o Phi_b,
//! ```
//!
//! the right pseudo-division of the row `M e_c` by `Phi`
//! ([`crate::right_pseudo_division_of_row`]) with remainder zero.  Applied
//! to `V` at a base `n >= b_0` every term `(H_b Phi_b V)_n` combines rows
//! `(Phi_b V)_(n + nu)`, `nu >= 0`, which vanish, so `w(n) (M V_c)_n = 0`
//! and `(M V_c)_n = 0` wherever `w(n) != 0` (`Q[x]` has no zero divisors).
//! The finitely many integer roots of `w` are the caller's to check on rows.
//! [`ComponentAnnihilatorCertificate::verify`] recomputes the identity.
//!
//! `M` is searched as in [`crate::find_image_annihilator`]: an ansatz
//! `m_R S^R + sum_(i < R, d <= J) m_(i,d) D^d S^i` with coefficients of
//! bounded degree in `n` and `x`; for a fixed top level the remainder of
//! the division is linear in `M`, its vanishing is a homogeneous linear
//! system, and a kernel vector with a nonzero top-shift entry is computed
//! exactly ([`crate::exact_kernel_vectors_by_group`]).  The shapes run by
//! increasing order (from `k`, within the budget's
//! `max_order_increase`) plus derivative order and increasing coefficient
//! degrees.  A found operator is normalised
//! ([`crate::normalize_annihilator`]) and certified; when nothing is found
//! the search says so, nothing is guessed.

use crate::bivariate::BivariatePolynomial;
use crate::division::{combine_rows, right_pseudo_division_of_row, system_index_only_leader};
use crate::image_annihilator::{
    ansatz_unknowns, is_runnable, kernel_candidate, normalize_annihilator, shapes, AnsatzShape,
    ImageAnnihilatorBudget,
};
use crate::operator::OreOperator;
use std::collections::HashMap;
use std::fmt;

/// Why [`find_component_annihilator`] or a certificate check failed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ComponentAnnihilatorError {
    /// The system is not square, or its shift-`k` part is not `g(n) S^k I`.
    SystemLeaderNotIndexOnly,
    /// The component index is not below the system's size.
    ComponentOutOfRange,
    /// The recorded annihilator is not runnable.
    NotRunnable,
    /// The recorded multiplier is zero or involves `x`.
    InvalidMultiplier,
    /// The identity `w (M e_c) = sum_b H_b o Phi_b` does not hold.
    IdentityFails,
}

impl fmt::Display for ComponentAnnihilatorError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::SystemLeaderNotIndexOnly => {
                "the system must be square with shift-k part g(n) S^k I, g free of x"
            }
            Self::ComponentOutOfRange => "the component index exceeds the system's size",
            Self::NotRunnable => {
                "the annihilator is not runnable (its top shift column is not a single term \
                 without derivatives)"
            }
            Self::InvalidMultiplier => "the multiplier must be a nonzero polynomial in n",
            Self::IdentityFails => "the identity w (M e_c) = sum_b H_b o Phi_b does not hold",
        })
    }
}

impl std::error::Error for ComponentAnnihilatorError {}

/// The certificate of a runnable annihilator of one component (module
/// documentation).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ComponentAnnihilatorCertificate {
    /// `M`, runnable, acting on the component `V_c`.
    pub annihilator: OreOperator,
    /// `c`.
    pub component: usize,
    /// `Phi`, by rows.
    pub system: Vec<Vec<OreOperator>>,
    /// `w(n)`, a nonzero polynomial in `n` only.
    pub multiplier: BivariatePolynomial,
    /// `H_b`, one per row of the system.
    pub quotients: Vec<OreOperator>,
}

impl ComponentAnnihilatorCertificate {
    /// Check `w (M e_c) = sum_b H_b o Phi_b` exactly, with `M` runnable, `w`
    /// a nonzero polynomial in `n` and the system's leader `g(n) S^k I`.
    pub fn verify(&self) -> Result<(), ComponentAnnihilatorError> {
        if !is_runnable(&self.annihilator) {
            return Err(ComponentAnnihilatorError::NotRunnable);
        }
        system_index_only_leader(&self.system)
            .ok_or(ComponentAnnihilatorError::SystemLeaderNotIndexOnly)?;
        if self.component >= self.system.len() {
            return Err(ComponentAnnihilatorError::ComponentOutOfRange);
        }
        if self.multiplier.is_zero() || !self.multiplier.is_variable_free() {
            return Err(ComponentAnnihilatorError::InvalidMultiplier);
        }
        if self.quotients.len() != self.system.len() {
            return Err(ComponentAnnihilatorError::IdentityFails);
        }
        let right = combine_rows(&self.quotients, &self.system);
        let left = component_row(
            &self.annihilator.left_multiply(&self.multiplier),
            self.component,
            self.system.len(),
        );
        if left != right {
            return Err(ComponentAnnihilatorError::IdentityFails);
        }
        Ok(())
    }
}

/// The outcome of [`find_component_annihilator`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ComponentAnnihilatorSearch {
    /// A runnable annihilator with a verified certificate.
    Found {
        /// The certificate.
        certificate: Box<ComponentAnnihilatorCertificate>,
        /// The ansatz shape it was found with.
        shape: AnsatzShape,
    },
    /// No runnable annihilator within the budget (no statement about larger
    /// shapes).
    NotFoundWithinBudget {
        /// The number of ansatz shapes solved.
        shapes_tried: usize,
        /// The largest shape solved, by number of unknowns.
        largest_shape: Option<AnsatzShape>,
    },
}

/// The row with `op` in entry `component` and zero elsewhere.
fn component_row(op: &OreOperator, component: usize, size: usize) -> Vec<OreOperator> {
    let mut row = vec![OreOperator::default(); size];
    row[component] = op.clone();
    row
}

/// Certify a candidate `M` for the component `component` of `system`
/// (module documentation): `None` when it is not runnable or the division
/// leaves a remainder.
pub fn certify_component_annihilator(
    candidate: &OreOperator,
    system: &[Vec<OreOperator>],
    component: usize,
) -> Result<Option<ComponentAnnihilatorCertificate>, ComponentAnnihilatorError> {
    system_index_only_leader(system).ok_or(ComponentAnnihilatorError::SystemLeaderNotIndexOnly)?;
    if component >= system.len() {
        return Err(ComponentAnnihilatorError::ComponentOutOfRange);
    }
    if !is_runnable(candidate) {
        return Ok(None);
    }
    let row = component_row(candidate, component, system.len());
    let Some(division) = right_pseudo_division_of_row(&row, system, candidate.max_shift()) else {
        return Ok(None);
    };
    if division.remainder.iter().any(|op| !op.is_zero()) {
        return Ok(None);
    }
    let certificate = ComponentAnnihilatorCertificate {
        annihilator: candidate.clone(),
        component,
        system: system.to_vec(),
        multiplier: division.multiplier,
        quotients: division.quotients,
    };
    // Defensive: the returned certificate always replays.
    Ok(certificate.verify().is_ok().then_some(certificate))
}

/// Solve one ansatz: a candidate with zero remainder, or `None`.
fn solve_shape(
    shape: &AnsatzShape,
    system: &[Vec<OreOperator>],
    component: usize,
) -> Result<Option<OreOperator>, ComponentAnnihilatorError> {
    let (unknowns, groups) = ansatz_unknowns(shape);
    let mut remainders: HashMap<(usize, usize, usize), Vec<OreOperator>> = HashMap::new();
    for unknown in &unknowns {
        let key = (unknown.shift, unknown.derivative, unknown.variable_power);
        if remainders.contains_key(&key) {
            continue;
        }
        let row = component_row(&unknown.index_free_operator(), component, system.len());
        let division = right_pseudo_division_of_row(&row, system, shape.order)
            .ok_or(ComponentAnnihilatorError::SystemLeaderNotIndexOnly)?;
        remainders.insert(key, division.remainder);
    }
    Ok(kernel_candidate(&unknowns, &groups, |unknown| {
        remainders[&(unknown.shift, unknown.derivative, unknown.variable_power)].as_slice()
    }))
}

/// A runnable annihilator of the component `component` of a vector
/// sequence annihilated by `system`, with a verified certificate, or the
/// report that none was found within `budget` (module documentation).
pub fn find_component_annihilator(
    system: &[Vec<OreOperator>],
    component: usize,
    budget: &ImageAnnihilatorBudget,
) -> Result<ComponentAnnihilatorSearch, ComponentAnnihilatorError> {
    let (order, _) = system_index_only_leader(system)
        .ok_or(ComponentAnnihilatorError::SystemLeaderNotIndexOnly)?;
    if component >= system.len() {
        return Err(ComponentAnnihilatorError::ComponentOutOfRange);
    }
    let mut failed: Vec<AnsatzShape> = Vec::new();
    let mut tried = 0;
    let mut largest = None;
    for shape in shapes(order, budget) {
        if tried >= budget.max_shapes {
            break;
        }
        if failed.iter().any(|f| shape.dominated_by(f)) {
            continue;
        }
        tried += 1;
        largest = Some(shape);
        if let Some(candidate) = solve_shape(&shape, system, component)? {
            for operator in [normalize_annihilator(&candidate), candidate] {
                if let Some(certificate) =
                    certify_component_annihilator(&operator, system, component)?
                {
                    return Ok(ComponentAnnihilatorSearch::Found {
                        certificate: Box::new(certificate),
                        shape,
                    });
                }
            }
        }
        failed.push(shape);
    }
    Ok(ComponentAnnihilatorSearch::NotFoundWithinBudget {
        shapes_tried: tried,
        largest_shape: largest,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Rational;

    fn r(value: i64) -> Rational {
        Rational::from_integer(value.into())
    }

    /// `V_0 = F_(n+1)`, `V_1 = F_n` for the Fibonacci-like polynomials
    /// `F_(n+2) = F_(n+1) + x F_n`, as the system `S V_0 = V_0 + x V_1`,
    /// `S V_1 = V_0`: the first component satisfies the scalar rule of
    /// order two.
    #[test]
    fn a_component_of_a_companion_system() {
        let s = OreOperator::single(1, 0, BivariatePolynomial::one());
        let one = OreOperator::identity();
        let x = OreOperator::single(
            0,
            0,
            BivariatePolynomial::from_variable_polynomial(&[r(0), r(1)]),
        );
        let system = vec![
            vec![s.subtract(&one), x.scale(&r(-1))],
            vec![one.scale(&r(-1)), s.clone()],
        ];
        let search =
            find_component_annihilator(&system, 0, &ImageAnnihilatorBudget::default()).unwrap();
        let ComponentAnnihilatorSearch::Found { certificate, .. } = search else {
            panic!("{search:?}")
        };
        certificate.verify().unwrap();
        // M = S^2 - S - x.
        let mut expected = OreOperator::single(2, 0, BivariatePolynomial::one());
        expected = expected.subtract(&s).subtract(&x);
        assert_eq!(certificate.annihilator, expected);
        // The rows: F_0 = 0, F_1 = 1.
        let mut rows = vec![vec![r(0)], vec![r(1)]];
        for n in 0..6 {
            let next = crate::univariate::add(
                &rows[n + 1],
                &crate::univariate::mul(&[r(0), r(1)], &rows[n]),
            );
            rows.push(next);
        }
        for n in 0..5 {
            let value = certificate.annihilator.apply(&rows, 0, n).unwrap();
            assert!(crate::univariate::is_zero(&value));
        }
        // A tampered certificate fails.
        let mut tampered = (*certificate).clone();
        tampered.multiplier = BivariatePolynomial::constant(r(2));
        assert_eq!(
            tampered.verify(),
            Err(ComponentAnnihilatorError::IdentityFails)
        );
        assert_eq!(
            certify_component_annihilator(&s, &system, 0).unwrap(),
            None,
            "S alone does not annihilate V_0"
        );
    }
}
