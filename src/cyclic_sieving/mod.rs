//! Exact root-of-unity checks for cyclic sieving.
//!
//! The implementation avoids floating point arithmetic.  For a primitive
//! `r`-th root of unity, `P(zeta) = m` if and only if `P(q) - m` is divisible by
//! the cyclotomic polynomial `Phi_r(q)`.
//!
//! The module also decomposes an integral polynomial into cyclic orbits
//! ([`cyclic_sieving_orbit_decomposition_bigint`]), applies the
//! Alexandersson--Amini existence criterion for a cyclic action, converts
//! between orbit sizes and fixed-point counts, and checks bicyclic sieving
//! ([`bicyclic_sieving_report_bigint`]).  Reference: P. Alexandersson and
//! N. Amini, *The cone of cyclic sieving phenomena*, Discrete Math. 342 (2019)
//! 1581-1601, arXiv:1804.01447.
//!
//! Entry points: [`cyclic_sieving_report_bigint`] (check one polynomial and group
//! order against fixed-point counts),
//! [`cyclic_sieving_sequence_reports_bigint`] (scan a sequence of rows and
//! candidate orders), [`cyclic_sieving_orbit_decomposition_bigint`],
//! [`check_cyclic_sieving_orbit_sizes`] and [`bicyclic_sieving_report_bigint`].
//!
//! Conventions: polynomials are `&[BigInt]` in ascending degree order (`q^i` has
//! index `i`); a cyclic group `C_n` of order `n` has generator `g`, and the
//! "power" `d` in `0..n` stands for `g^d` and the root of unity
//! `zeta_n^d`.  Fixed-point tables are indexed by `d` (entry `d` is
//! `fix(g^d)`, entry `0` the size of the set).  All arithmetic is exact `BigInt`.
//! Functions return `Result<_, CyclicSievingError>` (its `Display` is a
//! human-readable message) for invalid input (zero group order, wrong table
//! length, ...) and never panic on such input, except that [`cyclotomic_polynomial_bigint`],
//! [`reduce_modulo_cyclic_bigint`] and [`root_of_unity_evaluation_bigint`]
//! panic on a zero order.  A full report reduces `P` once per divisor of the
//! group order and its size is linear in the order (one evaluation per
//! power), so [`cyclic_sieving_report_bigint`] and
//! [`cyclic_sieving_sequence_reports_bigint`] reject orders above
//! [`CYCLIC_SIEVING_MAX_ORDER`]; [`bicyclic_sieving_report_bigint`] enforces
//! [`BICYCLIC_MAX_GROUP_SIZE`].

use crate::vec_poly::{degree, is_zero, strip_trailing_zeros};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Zero};
use std::collections::BTreeMap;

/// Largest cyclic group order accepted by [`cyclic_sieving_report_bigint`] and
/// [`cyclic_sieving_sequence_reports_bigint`] (one million).
///
/// A report holds one evaluation per power of the generator, so its size is
/// linear in the order; larger orders are rejected with
/// [`CyclicSievingError::OrderTooLarge`] instead of being allocated.
pub const CYCLIC_SIEVING_MAX_ORDER: usize = 1_000_000;

/// Invalid input to a cyclic or bicyclic sieving function.
///
/// The `Display` text is a human-readable message naming the offending value.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CyclicSievingError {
    /// The cyclic group order is zero.
    NonpositiveOrder,
    /// One of the two orders of a bicyclic group is zero.
    NonpositiveBicyclicOrders,
    /// A fixed-point table does not have exactly `order` entries.
    FixedPointCountLength {
        /// The group order `n`, which is also the required table length.
        order: usize,
        /// The length of the table that was supplied.
        found: usize,
    },
    /// An orbit size is zero or does not divide the group order.
    OrbitSizeDoesNotDivideOrder {
        /// The offending orbit size.
        size: usize,
        /// The group order `n`.
        order: usize,
    },
    /// A prescribed number of orbits is negative.
    NegativeOrbitCount {
        /// The orbit size whose count is negative.
        size: usize,
    },
    /// A fixed-point table is not constant on classes of equal `gcd(d, n)`, so
    /// it is not the table of a `C_n`-action.
    FixedPointsNotConstantOnGcdClasses {
        /// The group order `n`.
        order: usize,
        /// The power `d` (zero-based, in `0..n`) at which the table fails.
        power: usize,
        /// `gcd(d, n)`, the class representative the table is compared with.
        gcd: usize,
        /// The table entry `fix(g^d)`.
        value: BigInt,
        /// The table entry `fix(g^gcd)`.
        class_value: BigInt,
    },
    /// Mobius inversion of a fixed-point table gives a negative number of
    /// elements in orbits of size `orbit_size`.
    NegativeOrbitElementCount {
        /// The group order `n`.
        order: usize,
        /// The orbit size `k`.
        orbit_size: usize,
        /// The negative element count `S_k`.
        element_count: BigInt,
    },
    /// Mobius inversion of a fixed-point table gives a number of elements in
    /// orbits of size `orbit_size` that is not a multiple of it.
    OrbitElementCountNotDivisible {
        /// The group order `n`.
        order: usize,
        /// The orbit size `k`.
        orbit_size: usize,
        /// The element count `S_k`, not divisible by `k`.
        element_count: BigInt,
    },
    /// A group order exceeds [`CYCLIC_SIEVING_MAX_ORDER`].
    OrderTooLarge {
        /// The offending order (a sequence candidate order `index + offset` may
        /// exceed `usize`, hence `u128`).
        order: u128,
        /// The supported maximum, [`CYCLIC_SIEVING_MAX_ORDER`].
        maximum: usize,
    },
    /// In a sequence scan, `first_index + row` does not fit in `isize`.
    SequenceIndexOverflow {
        /// The requested index of row 0.
        first_index: isize,
        /// The first row whose index overflows.
        row: usize,
    },
    /// `order1 * order2` exceeds [`BICYCLIC_MAX_GROUP_SIZE`].
    BicyclicGroupTooLarge {
        /// The order of the first cyclic factor.
        order1: usize,
        /// The order of the second cyclic factor.
        order2: usize,
    },
    /// A bicyclic fixed-point table is not an `order1 x order2` matrix.
    BicyclicFixedTableShape {
        /// The order of the first cyclic factor (required number of rows).
        order1: usize,
        /// The order of the second cyclic factor (required row length).
        order2: usize,
    },
}

impl std::fmt::Display for CyclicSievingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonpositiveOrder => write!(f, "cyclic group order must be positive"),
            Self::NonpositiveBicyclicOrders => write!(f, "cyclic group orders must be positive"),
            Self::FixedPointCountLength { order, found } => write!(
                f,
                "expected {order} fixed-point counts for order {order}, got {found}"
            ),
            Self::OrbitSizeDoesNotDivideOrder { size, order } => {
                write!(f, "orbit size {size} does not divide the order {order}")
            }
            Self::NegativeOrbitCount { size } => {
                write!(f, "negative number of orbits of size {size}")
            }
            Self::FixedPointsNotConstantOnGcdClasses {
                order,
                power,
                gcd,
                value,
                class_value,
            } => write!(
                f,
                "fixed-point counts are not those of a C_{order}-action: fix(g^{power}) = {value} \
                 but fix(g^{gcd}) = {class_value} although gcd({power}, {order}) = {gcd}"
            ),
            Self::NegativeOrbitElementCount {
                order,
                orbit_size,
                element_count,
            } => write!(
                f,
                "fixed-point counts are not those of a C_{order}-action: \
                 S_{orbit_size} = {element_count} is negative"
            ),
            Self::OrbitElementCountNotDivisible {
                order,
                orbit_size,
                element_count,
            } => write!(
                f,
                "fixed-point counts are not those of a C_{order}-action: \
                 S_{orbit_size} = {element_count} is not divisible by {orbit_size}"
            ),
            Self::OrderTooLarge { order, maximum } => write!(
                f,
                "cyclic group order {order} exceeds the supported maximum {maximum}"
            ),
            Self::SequenceIndexOverflow { first_index, row } => write!(
                f,
                "sequence index overflows: first index {first_index} plus row {row}"
            ),
            Self::BicyclicGroupTooLarge { order1, order2 } => write!(
                f,
                "group order {order1} x {order2} exceeds the limit {BICYCLIC_MAX_GROUP_SIZE}"
            ),
            Self::BicyclicFixedTableShape { order1, order2 } => {
                write!(f, "fixed-point table must be a {order1} x {order2} matrix")
            }
        }
    }
}

impl std::error::Error for CyclicSievingError {}

/// Exact evaluation data for one root of unity.
#[derive(Clone, Debug, Eq, PartialEq)]
///
/// Evaluates `P(zeta)` for `zeta = zeta_n^d`, a primitive root of unity of order
/// `root_order = n / gcd(n, d)` (order 1 for `d = 0`).
pub struct RootOfUnityEvaluation {
    /// The group order `n`.
    pub group_order: usize,
    /// The exponent `d`, reduced modulo `n` (so `0 <= power < n`).
    pub power: usize,
    /// Multiplicative order of `zeta_n^d`, i.e. `n / gcd(n, d)`.
    pub root_order: usize,
    /// `P(zeta)` as an integer if it is one (remainder zero or constant),
    /// otherwise `None`.
    pub integer_value: Option<BigInt>,
    /// `P(q) mod Phi_{root_order}(q)` in ascending order with trailing zeros
    /// removed (empty when the remainder is zero); for `root_order = 1` it is
    /// the single entry `[P(1)]`.
    pub remainder: Vec<BigInt>,
}

/// One CSP equality check against a supplied fixed-point count.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingPowerCheck {
    /// The exponent `d` of the group element `g^d`.
    pub power: usize,
    /// Order of `zeta_n^d`.
    pub root_order: usize,
    /// The supplied fixed-point count `fix(g^d)`.
    pub expected_fixed_points: BigInt,
    /// The exact evaluation `P(zeta_n^d)`.
    pub evaluation: RootOfUnityEvaluation,
    /// Whether `P(zeta_n^d)` is an integer equal to `expected_fixed_points`.
    pub holds: bool,
}

/// CSP report for one polynomial and one cyclic group order.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingReport {
    /// The group order `n`.
    pub order: usize,
    /// The polynomial's coefficients with trailing zeros removed.
    pub coefficients: Vec<BigInt>,
    /// The supplied fixed-point counts `fix(g^d)`, `d = 0..n-1`, if any.
    pub fixed_counts: Option<Vec<BigInt>>,
    /// Evaluations `P(zeta_n^d)` for `d = 0..n-1`.
    pub evaluations: Vec<RootOfUnityEvaluation>,
    /// One check per `d` when `fixed_counts` was supplied, otherwise empty.
    pub checks: Vec<CyclicSievingPowerCheck>,
    /// `Some(true)` iff counts were supplied and every check holds; `None` when
    /// no counts were supplied (profile only).
    pub holds: Option<bool>,
}

/// Sequence-level CSP report for one row and several candidate group orders.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CyclicSievingSequenceItem {
    /// Zero-based position of the polynomial in the input list.
    pub row: usize,
    /// Sequence index `first_index + row` of the polynomial.
    pub index: isize,
    /// The polynomial's coefficients with trailing zeros removed.
    pub coefficients: Vec<BigInt>,
    /// One report per candidate order `index + offset` that is positive.
    pub candidate_orders: Vec<CyclicSievingReport>,
}

mod evaluation;
use evaluation::*;
pub use evaluation::{
    cyclic_sieving_report_bigint, cyclic_sieving_sequence_reports_bigint,
    cyclotomic_polynomial_bigint, root_of_unity_evaluation_bigint,
};
mod orbits;
use orbits::*;
#[allow(deprecated)]
pub use orbits::{
    check_cyclic_sieving_orbit_sizes, cyclic_sieving_check_orbit_sizes,
    cyclic_sieving_orbit_decomposition_bigint, fixed_points_from_orbit_counts,
    fixed_points_from_orbit_sizes, orbit_counts_from_fixed_points, orbit_counts_from_orbit_sizes,
    reduce_modulo_cyclic_bigint, OrbitDecomposition, OrbitSizeEntry, OrbitSizesCheck, ResidueClass,
};
mod bicyclic;
pub use bicyclic::{
    bicyclic_sieving_report_bigint, BicyclicSievingReport, BICYCLIC_MAX_GROUP_SIZE,
};
#[cfg(test)]
mod tests;
