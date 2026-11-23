//! Discover and check interlacing refinements of recurrence systems.
//!
//! Supplied coordinates are checked before bounded synthesis. Polynomial and
//! Weyl coordinates use exact Ore composition, including the product rule.
//! Row fits remain finite evidence. A symbolic intertwiner, target recovery,
//! ordered seed and uniform preservation certificate are separate obligations.
//! Entries of preserving transitions are affine in x; a sufficient symbolic
//! copositivity expansion handles polynomial dependence on the source index.

use super::*;
use crate::recurrence::lift::affine_pencils;
use crate::vec_poly as vp;
use crate::Polynomial;
use num_traits::{One, Signed, Zero};
use ore_algebra::{BivariatePolynomial as Biv, OreOperator as Op};
use serde::{Deserialize, Serialize};

type Q = BigRational;
type OperatorMatrix = Vec<Vec<Op>>;
/// Serialized operator matrix: target, source, derivative, then bivariate coefficients.
pub type RefinementCoordinates = Vec<Vec<Vec<BivarPolyJson>>>;
/// Protocol tag of discovery and replay reports.
pub const INTERLACING_REFINEMENT_SCHEMA: &str = "polytool.interlacing-refinement.v1";

mod algebra;
mod preservation;
mod scale;
mod search;
use algebra::*;
pub use preservation::{
    check_affine_interlacing_preserver, AffinePreserverReport, OperatorRejection, PencilRejection,
    SymbolicPencilCertificate,
};
pub use scale::{
    solve_common_scaled_refinement_interval, solve_scaled_refinement_interval,
    ScaledRefinementInterval,
};
use search::*;

/// Explicit search limits; failures within them do not refute real-rootedness.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RefinementSearchLimits {
    /// Maximum state dimension (1..=5).
    pub maximum_dimension: usize,
    /// Maximum candidates tried (1..=256), including supplied candidates.
    pub maximum_candidates: usize,
    /// Maximum transition degree in x (0..=3; proof engine supports affine entries).
    pub variable_degree: usize,
    /// Maximum transition degree in n (0..=4).
    pub index_degree: usize,
    /// Maximum derivative order for coordinate discovery (0..=2).
    pub derivative_order: usize,
    /// Maximum scalar equations times unknowns for each exact fit (at most 2,000,000).
    pub maximum_fit_cells: usize,
    /// Maximum support restrictions tried for a nonnegative symbolic fit (1..=64).
    pub maximum_supports: usize,
}
impl Default for RefinementSearchLimits {
    fn default() -> Self {
        Self {
            maximum_dimension: 5,
            maximum_candidates: 64,
            variable_degree: 1,
            index_degree: 1,
            derivative_order: 1,
            maximum_fit_cells: 500_000,
            maximum_supports: 32,
        }
    }
}

/// A supplied refinement, checked before automatic candidates.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InterlacingRefinementCandidate {
    /// Human-readable candidate name.
    pub name: String,
    /// v_n=T_n u_n; rectangular and singular maps are allowed.
    #[serde(default)]
    pub coordinates: Option<RefinementCoordinates>,
    /// Proposed defining transfer. If absent, fit it exactly under the limits.
    #[serde(default)]
    pub transition: Option<TransformVectorRecurrenceJson>,
    /// Nonnegative constant weights recovering the target; empty means sum.
    #[serde(default)]
    pub projection: Vec<String>,
}

/// Discovery from a defining recurrence, a defining vector system, or finite rows.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InterlacingRefinementRequest {
    /// Scalar rule and seeds, converted to an exact companion Weyl system.
    #[serde(default)]
    pub recurrence: Option<RecurrenceJson>,
    /// Alternative defining vector system with an index-only denominator.
    #[serde(default)]
    pub source: Option<TransformVectorRecurrenceJson>,
    /// Target weights in source coordinates; empty selects the last coordinate.
    #[serde(default)]
    pub source_projection: Vec<String>,
    /// Finite supplied vectors, indexed from first_index (finite evidence only).
    #[serde(default)]
    pub states: Vec<Vec<Vec<String>>>,
    /// Optional target rows for exact finite collapse checks.
    #[serde(default)]
    pub target_rows: Vec<Vec<String>>,
    /// First source index when only finite vectors are supplied.
    #[serde(default)]
    pub first_index: usize,
    /// Supplied candidates, tried in this order before automatic proposals.
    #[serde(default)]
    pub candidates: Vec<InterlacingRefinementCandidate>,
    /// Generate companion, last-step difference, derivative and scaled coordinates.
    #[serde(default = "yes")]
    pub derive_candidates: bool,
    /// Explicit caps on synthesis and arithmetic work.
    #[serde(default)]
    pub limits: RefinementSearchLimits,
}
fn yes() -> bool {
    true
}

/// One candidate with distinct algebraic and preservation obligations.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InterlacingRefinementCheck {
    /// Candidate (including discovered coordinates), sufficient to replay it.
    pub candidate: InterlacingRefinementCandidate,
    /// Exact intertwining identity, finite_fit, failed, or unresolved.
    pub identity: String,
    /// Exact source projection identity, finite_fit, failed, or unresolved.
    pub target_recovery: String,
    /// Whether the initial vector is ordered, nonnegative and coordinatewise real-rooted.
    pub initial_invariant: bool,
    /// Initial coordinate degrees; None denotes zero.
    pub initial_degrees: Vec<Option<usize>>,
    /// Pairwise directed degree compatibility of the initial vector.
    pub degree_conditions: bool,
    /// Fitted or supplied defining transfer, including exact initial vector.
    pub transition: Option<TransformVectorRecurrenceJson>,
    /// Uniform symbolic preservation report; absent for a failed fit.
    pub preservation: Option<AffinePreserverReport>,
    /// proved, finite_fit, unresolved, or rejected. Failure rejects this candidate only.
    pub status: String,
    /// Exact obstruction or description of what remains unproved.
    pub reason: String,
}

/// Replayable bounded search, with exact coefficients serialized as strings.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InterlacingRefinementReport {
    /// [`INTERLACING_REFINEMENT_SCHEMA`].
    pub schema: String,
    /// Complete input and limits.
    pub request: InterlacingRefinementRequest,
    /// proved, finite_fit, or not_proved.
    pub status: String,
    /// Source index from which the all-index conclusion applies.
    pub first_index: usize,
    /// Candidate checks in supplied-then-generated order.
    pub candidates: Vec<InterlacingRefinementCheck>,
    /// Index of the smallest successful refinement, ranked by dimension, degree and sparsity.
    pub selected: Option<usize>,
    /// True when proposal enumeration hit its explicit candidate cap.
    pub candidate_limit_reached: bool,
    /// Exact positive scale interval when the third-order scalar template matches.
    pub scaled_interval: Option<ScaledRefinementInterval>,
    /// Theorem and scope, including zero-aware and repeated-index conventions.
    pub reference: String,
    /// Rust certificate replay does not run Lean; always not_attempted.
    pub lean_replay: String,
}

fn invalid(message: impl Into<String>) -> RecurrenceError {
    RecurrenceError::InvalidInput(message.into())
}
fn numbers(row: &[String]) -> Result<Vec<Q>, RecurrenceError> {
    row.iter()
        .map(|s| parse_rational_coeff(s).map_err(Into::into))
        .collect()
}
fn strings(row: &[Q]) -> Vec<String> {
    vp::trim(row).iter().map(ToString::to_string).collect()
}

/// Discover bounded interlacing refinements. A proved result requires a
/// defining system, an exact Ore intertwiner and target recovery, an ordered
/// seed, and uniform Brändén preservation. Finite vectors alone never prove
/// an infinite family. The collapse is a nonnegative scalar combination;
/// no claim of consecutive interlacing of collapsed rows is made.
///
/// # Errors
/// Malformed inputs, conflicting sources or limits beyond supported caps.
pub fn discover_interlacing_refinements(
    request: &InterlacingRefinementRequest,
) -> Result<InterlacingRefinementReport, RecurrenceError> {
    run_search(request)
}

/// Check every stored transition, Ore identity, seed, degree and preservation
/// pencil without coordinate synthesis or fitting. Unsuccessful candidates
/// without a transition still regenerate their failed bounded fit. Candidate
/// enumeration order/completeness is discovery metadata, not a proof obligation.
///
/// # Errors
/// Invalid mathematical evidence, inconsistent verdict metadata, or an invalid request.
pub fn verify_interlacing_refinements(
    report: &InterlacingRefinementReport,
) -> Result<(), RecurrenceError> {
    search::verify(report)
}

#[cfg(test)]
mod tests;
