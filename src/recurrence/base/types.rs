//! Core recurrence data types: search options, bivariate coefficients,
//! recurrence terms, scalar recurrences and normal-ordered Weyl operators.

use super::*;
use num_traits::One;
use num_traits::Zero;

/// Options controlling the search space for `find_polynomial_recurrence`.
#[derive(Debug, Clone)]
pub struct RecurrenceOptions {
    /// Max degree of coefficient polynomials in the variable t.
    pub var_deg: usize,
    /// Max degree of coefficient polynomials in the index n.
    pub idx_deg: usize,
    /// Max order of differentiation (0 = no derivatives).
    pub diff_deg: usize,
    /// How many previous terms the recurrence may use.
    pub rec_len: usize,
    /// If false, allow an additive inhomogeneous polynomial g(n,t).
    pub homogeneous: bool,
    /// Max degree in t of the inhomogeneous term g(n,t).
    ///
    /// Only used when `homogeneous == false`.
    pub inhomo_var_deg: usize,
    /// Max degree in n of the inhomogeneous term g(n,t).
    ///
    /// Only used when `homogeneous == false`.
    pub inhomo_idx_deg: usize,
    /// Degree in t of the LHS factor f(n,t) beyond the implicit constant 1.
    pub denom_var_deg: usize,
    /// Degree in n of the LHS factor f(n,t) beyond the implicit constant 1.
    pub denom_idx_deg: usize,
    /// Also allow recurrence coefficient terms multiplied by (-1)^n.
    pub alternating_sign: bool,
    /// Use modular consistency checks to reject candidates before exact solving.
    ///
    /// Modular rejection is used only with a full-column-rank certificate;
    /// rank-deficient images fall back to exact rational solving.
    pub modular_prefilter: bool,
}

impl Default for RecurrenceOptions {
    fn default() -> Self {
        Self {
            var_deg: 1,
            idx_deg: 1,
            diff_deg: 0,
            rec_len: 2,
            homogeneous: true,
            inhomo_var_deg: 1,
            inhomo_idx_deg: 1,
            denom_var_deg: 0,
            denom_idx_deg: 0,
            alternating_sign: false,
            modular_prefilter: true,
        }
    }
}

impl RecurrenceOptions {
    /// Return these options with alternating `(-1)^n` recurrence terms enabled
    /// or disabled.
    pub fn with_alternating_sign(mut self, alternating_sign: bool) -> Self {
        self.alternating_sign = alternating_sign;
        self
    }
}

// ---------------------------------------------------------------------------
// Data types
// ---------------------------------------------------------------------------

/// A polynomial in two variables (n, t).
///
/// Stored as `coeffs[i][j]` = coefficient of n^i t^j.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BivarPoly {
    /// Coefficient array: `coeffs[i][j]` is the coefficient of `n^i t^j`. Rows may have
    /// different lengths; missing entries are zero.
    pub coeffs: Vec<Vec<BigRational>>,
}

impl BivarPoly {
    /// Construct a rectangular coefficient matrix, padding ragged rows with
    /// exact zeros. Empty input is normalized to the zero polynomial.
    pub fn from_coefficients_normalized(mut coeffs: Vec<Vec<BigRational>>) -> Self {
        let width = coeffs.iter().map(Vec::len).max().unwrap_or(0).max(1);
        if coeffs.is_empty() {
            coeffs.push(Vec::new());
        }
        for row in &mut coeffs {
            row.resize(width, BigRational::zero());
        }
        Self { coeffs }
    }

    pub(super) fn coefficient(&self, i: usize, j: usize) -> Option<&BigRational> {
        self.coeffs.get(i).and_then(|row| row.get(j))
    }

    pub(super) fn coefficient_width(&self) -> usize {
        self.coeffs.iter().map(Vec::len).max().unwrap_or(0)
    }

    pub(in crate::recurrence) fn normalized_coefficients(&self) -> Vec<Vec<BigRational>> {
        Self::from_coefficients_normalized(self.coeffs.clone()).coeffs
    }

    /// Substitute `n + shift` for the recurrence index `n`.
    pub fn shift_index(&self, shift: i64) -> Self {
        let width = self.coefficient_width().max(1);
        let mut shifted = vec![vec![BigRational::zero(); width]; self.coeffs.len().max(1)];
        let shift = BigRational::from_integer(BigInt::from(shift));
        let mut power = vec![BigRational::one()];

        for (n_degree, row) in self.coeffs.iter().enumerate() {
            debug_assert_eq!(power.len(), n_degree + 1);
            for (variable_degree, coefficient) in row.iter().enumerate() {
                for (new_n_degree, factor) in power.iter().enumerate() {
                    shifted[new_n_degree][variable_degree] += coefficient * factor;
                }
            }

            let mut next_power = vec![BigRational::zero(); power.len() + 1];
            for (degree, coefficient) in power.iter().enumerate() {
                next_power[degree] += coefficient * &shift;
                next_power[degree + 1] += coefficient;
            }
            power = next_power;
        }

        Self::from_coefficients_normalized(shifted)
    }

    pub(super) fn negated(mut self) -> Self {
        for row in &mut self.coeffs {
            for coefficient in row {
                *coefficient = -coefficient.clone();
            }
        }
        self
    }
}

/// Extra sign factor attached to one recurrence term.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RecurrenceSign {
    /// No extra sign factor.
    #[default]
    None,
    /// Multiply the term by (-1)^n.
    AlternatingN,
}

impl RecurrenceSign {
    pub(in crate::recurrence) fn from_family_index(sign_idx: usize) -> Self {
        match sign_idx {
            0 => Self::None,
            1 => Self::AlternatingN,
            _ => unreachable!("invalid recurrence sign family"),
        }
    }
}

/// One term in a recurrence: sign(n) * coeff(n,t) * D^d P_{n-r}(t).
#[derive(Debug, Clone)]
pub struct RecurrenceTerm {
    /// Recurrence offset r (P_{n-r}).
    pub offset: usize,
    /// Derivative order d.
    pub deriv_order: usize,
    /// Optional multiplicative sign factor.
    pub sign: RecurrenceSign,
    /// Coefficient polynomial c(n,t).
    pub coeff: BivarPoly,
}

impl RecurrenceTerm {
    /// Construct an ordinary non-alternating recurrence term.
    pub fn new(offset: usize, deriv_order: usize, coeff: BivarPoly) -> Self {
        Self {
            offset,
            deriv_order,
            sign: RecurrenceSign::None,
            coeff,
        }
    }

    /// Return this term with a different sign factor.
    pub fn with_sign(mut self, sign: RecurrenceSign) -> Self {
        self.sign = sign;
        self
    }

    /// True when this term is multiplied by `(-1)^n`.
    pub fn is_alternating(&self) -> bool {
        self.sign == RecurrenceSign::AlternatingN
    }
}

/// A polynomial recurrence found by `find_polynomial_recurrence`.
#[derive(Debug, Clone)]
pub struct Recurrence {
    /// RHS terms: each is c(n,t) D^d P_{n-r}.
    pub terms: Vec<RecurrenceTerm>,
    /// LHS factor f(n,t) if non-trivial (None means f = 1).
    pub denominator: Option<BivarPoly>,
    /// Inhomogeneous additive term g(n,t), if present.
    pub inhomogeneous: Option<BivarPoly>,
}

/// A Weyl-algebra operator in normal order,
/// `sum_d coefficients[d](n,x) D_x^d`.
///
/// Putting every derivative to the right is no loss of generality because
/// `D_x x = x D_x + 1`.  The representation is therefore directly usable by
/// exact linear recurrence fitting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WeylOperator {
    /// Polynomial coefficient multiplying each derivative order.
    pub coefficients: Vec<BivarPoly>,
}

impl WeylOperator {
    /// Construct a normal-ordered Weyl operator.
    pub fn new(coefficients: Vec<BivarPoly>) -> Self {
        Self { coefficients }
    }

    /// The zero operator.
    pub fn zero() -> Self {
        Self {
            coefficients: Vec::new(),
        }
    }

    /// The identity operator.
    pub fn identity() -> Self {
        Self {
            coefficients: vec![BivarPoly {
                coeffs: vec![vec![BigRational::one()]],
            }],
        }
    }

    /// True when all normal-form coefficients vanish.
    pub fn is_zero(&self) -> bool {
        self.coefficients.iter().all(BivarPoly::is_zero)
    }
}

// ---------------------------------------------------------------------------
// BivarPoly helpers
// ---------------------------------------------------------------------------

impl BivarPoly {
    /// True when all coefficients are zero.
    pub fn is_zero(&self) -> bool {
        self.coeffs
            .iter()
            .all(|row| row.iter().all(|c| c.is_zero()))
    }

    /// True when the polynomial equals the constant 1.
    pub fn is_one(&self) -> bool {
        if !self.coefficient(0, 0).is_some_and(|c| c.is_one()) {
            return false;
        }
        for (i, row) in self.coeffs.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                if i == 0 && j == 0 {
                    if !c.is_one() {
                        return false;
                    }
                } else if !c.is_zero() {
                    return false;
                }
            }
        }
        true
    }

    /// True when the polynomial is a single monomial.
    pub(super) fn is_monomial(&self) -> bool {
        let count = self
            .coeffs
            .iter()
            .flat_map(|row| row.iter())
            .filter(|c| !c.is_zero())
            .count();
        count <= 1
    }

    /// Number of non-zero terms.
    pub(super) fn num_terms(&self) -> usize {
        self.coeffs
            .iter()
            .flat_map(|row| row.iter())
            .filter(|c| !c.is_zero())
            .count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn constant(value: i64) -> BivarPoly {
        BivarPoly {
            coeffs: vec![vec![BigRational::from_integer(BigInt::from(value))]],
        }
    }

    #[test]
    fn weyl_operator_zero_detection() {
        // The canonical zero, explicit zero coefficients (also ragged and nested), and
        // nonzero operators in any derivative slot.
        assert!(WeylOperator::zero().is_zero());
        assert!(WeylOperator::new(vec![constant(0), constant(0)]).is_zero());
        let ragged_zero = BivarPoly {
            coeffs: vec![vec![], vec![BigRational::zero(), BigRational::zero()]],
        };
        assert!(WeylOperator::new(vec![ragged_zero]).is_zero());
        assert!(!WeylOperator::identity().is_zero());
        assert!(!WeylOperator::new(vec![constant(0), constant(0), constant(3)]).is_zero());
        // A nonzero entry hidden at high degree in the second variable.
        let hidden = BivarPoly {
            coeffs: vec![
                vec![BigRational::zero()],
                vec![BigRational::zero(), BigRational::zero(), BigRational::one()],
            ],
        };
        assert!(!WeylOperator::new(vec![hidden]).is_zero());
    }
}
