//! Exact recognition of `P_n = sum_j r_j t^j P_(n-1-ell*j)`.
//!
//! With canonical initial data, `R(z)=sum r_j z^j`, `R(0)=1`, nonnegative
//! coefficients and real zeros imply `P_n << P_(n+s)` for `degree(R)*s <=
//! degree(R)*ell+1`. This follows by periodic weighted gap subsets and the
//! nondecreasing-predecessor interlacing induction. See Fisk, arXiv:math/0612833v2,
//! Example 16.10, for the unweighted fixed-delay antecedent.
//!
//! Recognition certifies the structure of an input recurrence, not that a
//! recurrence fitted to finitely many rows holds for the intended sequence.

use crate::is_real_rooted_bigint_coeffs;
use crate::recurrence::{BivarPoly, Recurrence, RecurrenceSign};
use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};
use std::collections::BTreeMap;

/// A PF seed and positive delay, checked exactly at construction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeightedRenewal {
    ell: usize,
    seed: Vec<BigRational>,
}

fn constant(p: &BivarPoly) -> Option<BigRational> {
    let mut value = BigRational::zero();
    for (i, row) in p.coeffs.iter().enumerate() {
        for (j, c) in row.iter().enumerate() {
            if i == 0 && j == 0 {
                value = c.clone();
            } else if !c.is_zero() {
                return None;
            }
        }
    }
    Some(value)
}

impl WeightedRenewal {
    /// Validate a seed in ascending degree order. Trailing zeros are removed.
    pub fn new(ell: usize, mut seed: Vec<BigRational>) -> Result<Self, String> {
        if ell == 0 {
            return Err("ell must be positive".into());
        }
        while seed.last().is_some_and(Zero::is_zero) {
            seed.pop();
        }
        if seed.first() != Some(&BigRational::one()) || seed.iter().any(Signed::is_negative) {
            return Err("seed must have constant 1 and nonnegative coefficients".into());
        }
        let scale = seed.iter().fold(BigInt::one(), |a, b| a.lcm(b.denom()));
        let integers: Vec<_> = seed
            .iter()
            .map(|c| c.numer() * (&scale / c.denom()))
            .collect();
        if !is_real_rooted_bigint_coeffs(&integers) {
            return Err("seed is not real-rooted".into());
        }
        ell.checked_mul(seed.len() - 1)
            .and_then(|n| n.checked_add(1))
            .ok_or("recurrence offset overflows usize")?;
        Ok(Self { ell, seed })
    }

    /// Delay parameter in the generating function `1/(1-x R(t x^ell))`.
    pub fn ell(&self) -> usize {
        self.ell
    }

    /// The certified PF seed, in ascending degree order.
    pub fn seed(&self) -> &[BigRational] {
        &self.seed
    }

    /// Largest guaranteed index step; `None` means every step (constant seed).
    pub fn max_interlacing_step(&self) -> Option<usize> {
        let d = self.seed.len() - 1;
        if d == 0 {
            None
        } else {
            Some((d * self.ell + 1) / d)
        }
    }

    /// Canonical rows starting at `P_0=1`, with negative-index rows zero.
    /// This is the library generator for the recognized family.
    pub fn rows(&self, total_rows: usize) -> Vec<Vec<BigRational>> {
        let mut rows: Vec<Vec<BigRational>> = Vec::with_capacity(total_rows);
        for n in 0..total_rows {
            if n == 0 {
                rows.push(vec![BigRational::one()]);
                continue;
            }
            let mut row = rows[n - 1].clone();
            for (j, coefficient) in self.seed.iter().enumerate().skip(1) {
                let offset = 1 + self.ell * j;
                if offset > n {
                    break;
                }
                let previous = &rows[n - offset];
                row.resize(row.len().max(previous.len() + j), BigRational::zero());
                for (k, c) in previous.iter().enumerate() {
                    row[k + j] += coefficient * c;
                }
            }
            while row.len() > 1 && row.last().is_some_and(Zero::is_zero) {
                row.pop();
            }
            rows.push(row);
        }
        rows
    }

    /// Check enough initial rows to determine the recurrence, starting at zero.
    /// Extra supplied rows are checked too. This does not prove a fitted identity.
    pub fn matches_initial_rows(&self, first_index: usize, rows: &[Vec<BigRational>]) -> bool {
        let required = 1 + self.ell * (self.seed.len() - 1);
        if first_index != 0 || rows.len() < required {
            return false;
        }
        self.rows(rows.len())
            .iter()
            .zip(rows)
            .all(|(expected, actual)| {
                let mut end = actual.len();
                while end > 0 && actual[end - 1].is_zero() {
                    end -= 1;
                }
                expected.as_slice() == &actual[..end]
            })
    }
}

/// Recognize a stationary homogeneous renewal recurrence and check its seed.
/// Constant nonzero LHS factors and split/duplicate RHS terms are normalized.
/// Derivatives, alternating signs, index dependence and other lag supports are
/// rejected conservatively. For `R=1`, the immaterial delay is returned as 1.
pub fn recognize_weighted_renewal(recurrence: &Recurrence) -> Result<WeightedRenewal, String> {
    let denominator = match &recurrence.denominator {
        None => BigRational::one(),
        Some(p) => constant(p)
            .filter(|c| !c.is_zero())
            .ok_or("LHS must be a nonzero constant")?,
    };
    if recurrence
        .inhomogeneous
        .as_ref()
        .is_some_and(|p| constant(p) != Some(BigRational::zero()))
    {
        return Err("recurrence must be homogeneous".into());
    }
    let mut terms = BTreeMap::<(usize, usize), BigRational>::new();
    for term in &recurrence.terms {
        if term.deriv_order != 0 || term.sign != RecurrenceSign::None {
            return Err("derivatives and alternating signs are unsupported".into());
        }
        for (i, row) in term.coeff.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                if i != 0 {
                    return Err("coefficients must be index-independent".into());
                }
                *terms.entry((term.offset, j)).or_default() += c / &denominator;
            }
        }
    }
    terms.retain(|_, c| !c.is_zero());
    if terms.remove(&(1, 0)) != Some(BigRational::one()) {
        return Err("coefficient of P_(n-1) must be 1 after normalization".into());
    }
    let mut ell = None;
    let mut seed = vec![BigRational::one()];
    for ((offset, j), c) in terms {
        if j == 0 || offset <= 1 || (offset - 1) % j != 0 {
            return Err("support must satisfy offset = 1 + ell * t_degree".into());
        }
        let delay = (offset - 1) / j;
        if ell.is_some_and(|old| old != delay) {
            return Err("inconsistent delays".into());
        }
        ell = Some(delay);
        seed.resize(seed.len().max(j + 1), BigRational::zero());
        seed[j] += c;
    }
    WeightedRenewal::new(ell.unwrap_or(1), seed)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::recurrence::RecurrenceTerm;
    fn q(n: i64) -> BigRational {
        BigRational::from_integer(n.into())
    }
    fn recurrence(ell: usize, seed: &[i64]) -> Recurrence {
        Recurrence {
            denominator: None,
            inhomogeneous: None,
            terms: seed
                .iter()
                .enumerate()
                .map(|(j, c)| {
                    let mut row = vec![q(0); j + 1];
                    row[j] = q(*c);
                    RecurrenceTerm::new(1 + ell * j, 0, BivarPoly { coeffs: vec![row] })
                })
                .collect(),
        }
    }
    #[test]
    fn pf_seed_and_independent_coefficient_formula() {
        let family = recognize_weighted_renewal(&recurrence(1, &[1, 3, 2])).unwrap();
        // [t^j] P_n = [z^j] (1+3z+2z^2)^(n-j).
        assert_eq!(
            family.rows(6),
            vec![
                vec![q(1)],
                vec![q(1)],
                vec![q(1), q(3)],
                vec![q(1), q(6), q(2)],
                vec![q(1), q(9), q(13)],
                vec![q(1), q(12), q(33), q(12)]
            ]
        );
        assert_eq!(family.max_interlacing_step(), Some(1));
        assert!(family.matches_initial_rows(0, &family.rows(5)));
        let mut wrong = family.rows(5);
        wrong[1] = vec![q(1), q(1)];
        assert!(!family.matches_initial_rows(0, &wrong));
        assert!(!family.matches_initial_rows(1, &family.rows(5)));
        assert!(!family.matches_initial_rows(0, &family.rows(2)));
    }
    #[test]
    fn normalize_duplicates_and_rational_seed() {
        let mut rec = recurrence(2, &[1, 2, 1]);
        rec.terms.extend(rec.terms.clone());
        rec.denominator = Some(BivarPoly {
            coeffs: vec![vec![q(2)]],
        });
        assert_eq!(
            recognize_weighted_renewal(&rec)
                .unwrap()
                .max_interlacing_step(),
            Some(2)
        );
        assert!(WeightedRenewal::new(1, vec![q(1), q(1), q(1) / q(4)]).is_ok());
        assert!(WeightedRenewal::new(1, vec![q(1), q(1), q(1) / q(3)]).is_err());
    }
    #[test]
    fn reject_false_positive_structures() {
        assert!(recognize_weighted_renewal(&recurrence(1, &[1, 1, 1])).is_err());
        assert!(recognize_weighted_renewal(&recurrence(1, &[1, -1])).is_err());
        assert!(recognize_weighted_renewal(&recurrence(0, &[1, 1])).is_err());
        let mut rec = recurrence(1, &[1, 2, 1]);
        rec.terms[2].offset += 1;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.terms[1].deriv_order = 1;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec.terms[1].deriv_order = 0;
        rec.terms[1].sign = RecurrenceSign::AlternatingN;
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.terms[1].coeff.coeffs.push(vec![q(1)]);
        assert!(recognize_weighted_renewal(&rec).is_err());
        rec = recurrence(1, &[1, 1]);
        rec.inhomogeneous = Some(BivarPoly {
            coeffs: vec![vec![q(1)]],
        });
        assert!(recognize_weighted_renewal(&rec).is_err());
    }
    #[test]
    fn constant_and_long_delay_boundaries() {
        let family = recognize_weighted_renewal(&recurrence(1, &[1])).unwrap();
        assert_eq!(family.max_interlacing_step(), None);
        assert_eq!(family.rows(3), vec![vec![q(1)]; 3]);
        let family = recognize_weighted_renewal(&recurrence(3, &[1, 2])).unwrap();
        assert_eq!(family.max_interlacing_step(), Some(4));
        assert_eq!(family.rows(5)[4], vec![q(1), q(2)]);
        assert!(WeightedRenewal::new(usize::MAX, vec![q(1), q(1)]).is_err());
    }
}
