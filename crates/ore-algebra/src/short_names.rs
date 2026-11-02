//! Short method names kept for existing call sites in `polytool`, which
//! predate this crate.  They are hidden from the documentation, forward to
//! the spelled-out names and will be removed once those call sites use the
//! spelled-out names; new code should not use them.

use crate::{BivariatePolynomial, OreOperator, Rational};

#[doc(hidden)]
#[allow(missing_docs)]
impl BivariatePolynomial {
    pub fn from_y(p: &[Rational]) -> Self {
        Self::from_variable_polynomial(p)
    }
    pub fn from_n(p: &[Rational]) -> Self {
        Self::from_index_polynomial(p)
    }
    pub fn is_n_free(&self) -> bool {
        self.is_index_free()
    }
    pub fn is_y_free(&self) -> bool {
        self.is_variable_free()
    }
    pub fn mul(&self, other: &Self) -> Self {
        self.multiply(other)
    }
    pub fn pow(&self, exponent: usize) -> Self {
        self.power(exponent)
    }
    pub fn shift_n(&self, index_shift: i64) -> Self {
        self.shift_index(index_shift)
    }
    pub fn dy(&self) -> Self {
        self.variable_derivative()
    }
    pub fn eval_n(&self, value: &Rational) -> Vec<Rational> {
        self.evaluate_index(value)
    }
    pub fn eval_y(&self, value: &Rational) -> Vec<Rational> {
        self.evaluate_variable(value)
    }
    pub fn y_coefficient(&self, j: usize) -> Vec<Rational> {
        self.variable_coefficient(j)
    }
    pub fn y_degree(&self) -> Option<usize> {
        self.variable_degree()
    }
    pub fn y_valuation(&self) -> Option<usize> {
        self.variable_valuation()
    }
    pub fn div_y_power(&self, e: usize) -> Self {
        self.divide_by_variable_power(e)
    }
    pub fn mul_y_power(&self, e: usize) -> Self {
        self.multiply_by_variable_power(e)
    }
}

#[doc(hidden)]
#[allow(missing_docs)]
impl OreOperator {
    pub fn sub(&self, other: &Self) -> Self {
        self.subtract(other)
    }
    pub fn x_valuation(&self) -> usize {
        self.variable_valuation()
    }
    pub fn div_x_power(&self, e: usize) -> Self {
        self.divide_by_variable_power(e)
    }
    pub fn mul_x_power(&self, e: usize) -> Self {
        self.multiply_by_variable_power(e)
    }
}
