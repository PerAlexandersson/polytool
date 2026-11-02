//! Exact Ore-algebra arithmetic for sequences of polynomials.
//!
//! A sequence of polynomials `F_n(x)` is acted on by the shift `S`
//! (`(S F)_n = F_(n+1)`), the derivative `D = d/dx` and multiplication by
//! polynomials in the index `n` and the variable `x`.  These generate the
//! Ore algebra with the commutation rules `S n = (n + 1) S` and
//! `D x = x D + 1`; every element has a unique normal-ordered form
//! `sum c_(i,d)(n, x) D^d S^i` ([`OreOperator`]) with coefficients that are
//! bivariate polynomials over the rationals ([`BivariatePolynomial`]).
//!
//! - [`BivariatePolynomial`]: exact polynomials in `(n, y)`;
//! - [`OreOperator`]: normal-ordered operators, their composition and their
//!   action on rows;
//! - the Euler form ([`to_euler_form`], [`from_euler_form`]) with
//!   `theta = x D`;
//! - [`exact_kernel_vectors_by_group`]: exact kernel vectors of integer
//!   matrices by modular elimination with an exact final check;
//! - [`right_pseudo_division`]: right pseudo-division by an operator whose
//!   leader is free of `x` and of derivatives;
//! - [`RowMapStep`], [`pull_back`], [`push_forward`]: maps between
//!   sequences (Weyl operators, gauges, affine maps, reversals, coefficient
//!   multipliers) and the transport of annihilators through them;
//! - [`find_annihilator_system`]: additional annihilators of the rows of a
//!   recurrence-defined sequence (row annihilators, runnable rules of lower
//!   order), with joint invariance certificates;
//! - [`find_image_annihilator`]: a runnable annihilator of an image
//!   sequence, pushed forward or found as a common left multiple, with an
//!   [`ImageAnnihilatorCertificate`] that replays the identity
//!   `w(n) E = H o L` on its own.
//!
//! Everything is exact (`num_rational::BigRational`); nothing uses floating
//! point.  Coefficient vectors are ascending: entry `i` is the coefficient
//! of `t^i`.
//!
//! ```
//! use ore_algebra::{BivariatePolynomial, OreOperator, Rational};
//!
//! let r = |v: i64| Rational::from_integer(v.into());
//! // D o x = x D + 1.
//! let d = OreOperator::single(0, 1, BivariatePolynomial::one());
//! let x = BivariatePolynomial::from_variable_polynomial(&[r(0), r(1)]);
//! let mut expected = OreOperator::single(0, 1, x.clone());
//! expected.add_term((0, 0), BivariatePolynomial::one());
//! assert_eq!(d.compose(&OreOperator::single(0, 0, x)), expected);
//! ```

mod bivariate;
mod difference_division;
mod difference_image_annihilator;
mod difference_row_map;
mod differential_difference;
mod division;
mod euler;
mod falling_factorial;
mod image_annihilator;
mod kernel;
mod laurent;
mod operator;
mod row_annihilator;
mod row_map;
mod short_names;
mod substitution;
mod univariate;

pub use bivariate::BivariatePolynomial;
pub use difference_division::{
    difference_index_only_leader, difference_right_pseudo_division, variable_shift_pseudo_division,
    DifferentialDifferenceDivision,
};
pub use difference_image_annihilator::{
    certify_difference_image_annihilator, find_difference_image_annihilator,
    is_difference_runnable, normalize_difference_annihilator, DifferenceAnsatzShape,
    DifferenceImageAnnihilatorBudget, DifferenceImageAnnihilatorCertificate,
    DifferenceImageAnnihilatorError, DifferenceImageAnnihilatorMethod,
    DifferenceImageAnnihilatorSearch,
};
pub use difference_row_map::{
    pull_back_difference, push_forward_difference, DifferencePullbackBounds, DifferenceRowMapError,
    DifferenceRowMapStep,
};
pub use differential_difference::{DifferentialDifferenceKey, DifferentialDifferenceOperator};
pub use division::{index_only_leader, right_pseudo_division, RightPseudoDivision};
pub use euler::{add_euler_term, falling_factorial, from_euler_form, to_euler_form, EulerOperator};
pub use falling_factorial::{
    conjugate_to_falling_factorial_basis, conjugate_to_power_basis,
    falling_factorial_euler_operator, falling_factorial_raising_operator,
    falling_factorial_to_power_basis, from_falling_factorial_euler_form,
    power_to_falling_factorial_basis, to_falling_factorial_euler_form,
};
pub use image_annihilator::{
    certify_image_annihilator, find_image_annihilator, is_runnable, normalize_annihilator,
    AnsatzShape, ImageAnnihilatorBudget, ImageAnnihilatorCertificate, ImageAnnihilatorError,
    ImageAnnihilatorMethod, ImageAnnihilatorSearch,
};
pub use kernel::exact_kernel_vectors_by_group;
pub use laurent::{LaurentBivariatePolynomial, LaurentError, LaurentOreOperator};
pub use operator::OreOperator;
pub use row_annihilator::{
    certify_annihilator_system, find_annihilator_system, reduce_by_row_annihilator,
    AnnihilatorSystemBudget, AnnihilatorSystemCertificate, AnnihilatorSystemEntry,
    AnnihilatorSystemError, RowAnnihilatorReduction,
};
pub use row_map::{pull_back, push_forward, PullbackBounds, RowMapError, RowMapStep};
pub use substitution::{
    pull_back_substitution, push_forward_substitution, MobiusMap, Substitution, SubstitutionError,
};

/// The exact rational type used for every coefficient.
pub type Rational = num_rational::BigRational;
