//! Report types of [`find_veronese_section_lifts`] and their
//! `polytool.recurrence-veronese-lifts.v1` JSON form.
//!
//! The report types hold exact [`BigRational`] values and a [`Recurrence`];
//! their [`Serialize`] implementations write every rational as an exact
//! integer or `num/den` string and every parent recurrence as a
//! `polytool.recurrence.v1` record with its seeds.

use super::veronese_algebra::format_shift_polynomial;
use super::*;
use serde::ser::{SerializeStruct, Serializer};
use serde::Serialize;

/// JSON schema tag of a serialized [`VeroneseSectionLiftReport`].
pub const VERONESE_SECTION_LIFT_SCHEMA: &str = "polytool.recurrence-veronese-lifts.v1";

/// Status of one Veronese-section lift.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VeroneseSectionLiftStatus {
    /// The parent is `P_(n+k) = c(x) P_n` (characteristic polynomial
    /// `L^k - c(x)`, for example any order-one parent), so
    /// `P_(n0 + k m + i) = c^m P_(n0 + i)`, and `c` and the `k` seeds are
    /// nonzero, real-rooted and have coefficients of one sign.
    /// Every parent row is then real-rooted with nonpositive roots, so every
    /// input row `Q_n`, `n >= n0`, is zero or real-rooted with nonpositive
    /// roots.
    ProvedProductParent,
    /// Every checked parent row is nonzero, real-rooted and has coefficients
    /// of one sign, but real-rootedness of all parent rows is not proved.
    ParentRealRootednessOpen,
    /// The checked parent row with this absolute index `n` is zero.
    ParentRowZero {
        /// Absolute row index `n` of the failing parent row `P_n`.
        index: usize,
    },
    /// The checked parent row with this absolute index has coefficients of
    /// both signs.
    ParentRowMixedSigns {
        /// Absolute row index `n` of the failing parent row `P_n`.
        index: usize,
    },
    /// The checked parent row with this absolute index has non-real roots.
    ParentRowNotRealRooted {
        /// Absolute row index `n` of the failing parent row `P_n`.
        index: usize,
    },
}

impl VeroneseSectionLiftStatus {
    /// Stable snake-case label used in the JSON report.
    pub fn label(&self) -> &'static str {
        match self {
            Self::ProvedProductParent => "proved_product_parent",
            Self::ParentRealRootednessOpen => "parent_real_rootedness_open",
            Self::ParentRowZero { .. } => "parent_row_zero",
            Self::ParentRowMixedSigns { .. } => "parent_row_mixed_signs",
            Self::ParentRowNotRealRooted { .. } => "parent_row_not_real_rooted",
        }
    }

    /// Whether the lift proves that every input row is zero or real-rooted.
    pub fn is_proved(&self) -> bool {
        matches!(self, Self::ProvedProductParent)
    }

    /// Whether a checked parent row fails, so this lift cannot give a proof.
    pub fn is_failing(&self) -> bool {
        self.failing_row_index().is_some()
    }

    /// Absolute index of the failing parent row, for a failing status.
    pub fn failing_row_index(&self) -> Option<usize> {
        match self {
            Self::ParentRowZero { index }
            | Self::ParentRowMixedSigns { index }
            | Self::ParentRowNotRealRooted { index } => Some(*index),
            Self::ProvedProductParent | Self::ParentRealRootednessOpen => None,
        }
    }

    /// Sorting group: proved, open, failing.
    pub(super) fn group(&self) -> usize {
        match self {
            Self::ProvedProductParent => 0,
            Self::ParentRealRootednessOpen => 1,
            _ => 2,
        }
    }

    /// One-sentence explanation of the status.
    pub fn statement(&self) -> String {
        match self {
            Self::ProvedProductParent => "parent P_(n+k) = c P_n, so P_(n0+km+i) = c^m P_(n0+i), \
                 with c and the k seeds real-rooted with coefficients of one sign: every parent \
                 row has only nonpositive real roots, so every input row Q_n (n >= n0) is zero or \
                 real-rooted with nonpositive roots"
                .into(),
            Self::ParentRealRootednessOpen => "every checked parent row is real-rooted with \
                 coefficients of one sign; real-rootedness of all parent rows is open"
                .into(),
            Self::ParentRowZero { index } => {
                format!("parent row P_{index} is zero; this lift gives no proof")
            }
            Self::ParentRowMixedSigns { index } => format!(
                "parent row P_{index} has coefficients of both signs; this lift gives no proof"
            ),
            Self::ParentRowNotRealRooted { index } => {
                format!("parent row P_{index} is not real-rooted; this lift gives no proof")
            }
        }
    }
}

/// Exact checks of the parent rows `P_(n0), ..., P_(n0 + rows_checked - 1)`.
///
/// Indices are absolute row indices `n`; a field is `None` when no checked
/// row fails in that way.  Real-rootedness is checked for every nonzero row,
/// also one with coefficients of both signs.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct VeroneseSectionParentRowCheck {
    /// Number of parent rows checked, starting at the parent's first index.
    pub rows_checked: usize,
    /// First checked row that is the zero polynomial.
    pub first_zero_row_index: Option<usize>,
    /// First checked nonzero row with coefficients of both signs.
    pub first_mixed_sign_row_index: Option<usize>,
    /// First checked nonzero row that is not real-rooted.
    pub first_not_real_rooted_row_index: Option<usize>,
}

/// One parent family `P_n(x)` whose coefficient section is the input sequence:
/// `V_(modulus, residue)(P_n) = Q_n` for every `n >= first_index`.
///
/// Coefficient layouts: a characteristic polynomial `[i][j]` is the
/// coefficient of `L^i x^j` (parent) or `L^i y^j` (section and input), with
/// `y = x^modulus`; each row is trimmed of trailing zeros, so a zero row is
/// empty.  Parent rows are ascending coefficient vectors in `x` (the
/// variable `t` of [`Recurrence`]); the zero row is `[0]`.
#[derive(Debug, Clone)]
pub struct VeroneseSectionLift {
    /// The modulus `r >= 2` of the section.
    pub modulus: usize,
    /// The residue `0 <= j < r` of the section.
    pub residue: usize,
    /// `G(L, x)`, monic in `L`: the parent satisfies
    /// `G(E, x) P_n = 0` for the shift `E P_n = P_(n+1)` and every `n >= n0`.
    pub parent_characteristic_polynomial: Vec<Vec<BigRational>>,
    /// `N(L, y) = prod_{w^r = 1} G(L, w x)` in `y = x^r`, monic of `L`-degree
    /// `r deg_L G`; every section of the parent satisfies it.
    pub section_characteristic_polynomial: Vec<Vec<BigRational>>,
    /// `chi / N`, where `chi` is the input characteristic polynomial.
    pub cofactor: Vec<Vec<BigRational>>,
    /// Parent recurrence `P_n = sum_k c_k(x) P_(n-k)` with `c_k = -g_(e-k)`
    /// for `G = L^e + sum_(i<e) g_i(x) L^i`; derivative-free, independent of
    /// `n`, no denominator.  It always has a term of offset `e`.
    pub parent_recurrence: Recurrence,
    /// Index `n0` of the first parent seed (the input's first index).
    pub parent_first_index: usize,
    /// Parent seeds `P_(n0), ..., P_(n0 + e - 1)`.
    pub parent_initial_rows: Vec<Vec<BigRational>>,
    /// Exact checks of the first parent rows.
    pub parent_row_check: VeroneseSectionParentRowCheck,
    /// Whether the lift proves real-rootedness of the input rows.
    pub status: VeroneseSectionLiftStatus,
}

impl VeroneseSectionLift {
    /// The parent recurrence with its seeds as a `polytool.recurrence.v1`
    /// record, ready for the other recurrence tools.
    pub fn parent_recurrence_json(&self) -> RecurrenceJson {
        RecurrenceJson::from_recurrence_rational(
            &self.parent_recurrence,
            self.parent_first_index,
            &self.parent_initial_rows,
            None,
        )
    }

    /// Degree of the parent recurrence (`deg_L G`).
    pub fn parent_order(&self) -> usize {
        self.parent_characteristic_polynomial
            .len()
            .saturating_sub(1)
    }

    /// Largest power of `x` in `G`.
    pub fn parent_degree_in_x(&self) -> usize {
        self.parent_characteristic_polynomial
            .iter()
            .map(|row| row.len().saturating_sub(1))
            .max()
            .unwrap_or(0)
    }
}

/// Why a candidate parent characteristic polynomial was not accepted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VeroneseSectionLiftRejectionReason {
    /// The section characteristic polynomial `N` does not divide the input
    /// characteristic polynomial `chi` in `Q[y][L]`.
    NormDoesNotDivideCharacteristicPolynomial,
    /// The linear map from the parent seed sections to the first input rows
    /// has rank below `r e`, so the seeds are not determined (not supported).
    SeedsNotDetermined,
    /// The unique seed sections are rational functions, not polynomials.
    SeedsNotPolynomial,
    /// The section of the parent row with this absolute index differs from
    /// the input row (the overdetermined rows disagree).
    RowMismatch {
        /// Absolute row index `n` of the first disagreeing row.
        index: usize,
    },
}

impl VeroneseSectionLiftRejectionReason {
    /// Stable snake-case label used in the JSON report.
    pub fn label(&self) -> &'static str {
        match self {
            Self::NormDoesNotDivideCharacteristicPolynomial => "norm_does_not_divide",
            Self::SeedsNotDetermined => "seeds_not_determined",
            Self::SeedsNotPolynomial => "seeds_not_polynomial",
            Self::RowMismatch { .. } => "row_mismatch",
        }
    }
}

/// A candidate parent that was tried and rejected.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VeroneseSectionLiftRejection {
    /// The modulus `r`.
    pub modulus: usize,
    /// The residue, or `None` when the candidate failed for every residue
    /// (the norm test does not depend on it).
    pub residue: Option<usize>,
    /// `G(L, x)`, `[i][j]` = coefficient of `L^i x^j`, rows trimmed.
    pub parent_characteristic_polynomial: Vec<Vec<BigRational>>,
    /// Why it was rejected.
    pub reason: VeroneseSectionLiftRejectionReason,
}

/// A modulus whose search could not run.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VeroneseSectionLiftSkippedModulus {
    /// The modulus `r`.
    pub modulus: usize,
    /// Why it was skipped (the factorization error text).
    pub reason: String,
}

/// Result of [`find_veronese_section_lifts`].
///
/// A bounded search: an empty `lifts` list means that no lift was found
/// within the options, not that none exists.
#[derive(Debug, Clone)]
pub struct VeroneseSectionLiftReport {
    /// `chi(L, y) = L^d - sum_k c_k(y) L^(d-k)` of the input recurrence
    /// `Q_m = sum_k c_k(y) Q_(m-k)`, monic in `L`; `[i][j]` = coefficient of
    /// `L^i y^j`, rows trimmed.
    pub input_characteristic_polynomial: Vec<Vec<BigRational>>,
    /// Accepted lifts: proved first, then open, then failing; within a group
    /// by modulus, `L`-degree and `x`-degree of `G`, then residue.
    pub lifts: Vec<VeroneseSectionLift>,
    /// Candidates that were tried and rejected, in search order.
    pub rejected_candidates: Vec<VeroneseSectionLiftRejection>,
    /// Moduli whose characteristic polynomial could not be factored.
    pub skipped_moduli: Vec<VeroneseSectionLiftSkippedModulus>,
    /// Whether the candidate budget stopped the search at some modulus.
    pub exhausted_budget: bool,
}

/// Plain-text form of a characteristic polynomial with `[i][j]` = coefficient
/// of `L^i v^j`, highest power of `L` first (`"L^2 - L - x"`).
pub fn format_characteristic_polynomial(
    coefficients: &[Vec<BigRational>],
    variable: &str,
) -> String {
    format_shift_polynomial(coefficients, variable)
}

fn rational_matrix_strings(matrix: &[Vec<BigRational>]) -> Vec<Vec<String>> {
    matrix
        .iter()
        .map(|row| row.iter().map(format_rational_coeff).collect())
        .collect()
}

impl Serialize for VeroneseSectionParentRowCheck {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("VeroneseSectionParentRowCheck", 4)?;
        state.serialize_field("rows_checked", &self.rows_checked)?;
        state.serialize_field("first_zero_row_index", &self.first_zero_row_index)?;
        state.serialize_field(
            "first_mixed_sign_row_index",
            &self.first_mixed_sign_row_index,
        )?;
        state.serialize_field(
            "first_not_real_rooted_row_index",
            &self.first_not_real_rooted_row_index,
        )?;
        state.end()
    }
}

impl Serialize for VeroneseSectionLift {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("VeroneseSectionLift", 13)?;
        state.serialize_field("modulus", &self.modulus)?;
        state.serialize_field("residue", &self.residue)?;
        state.serialize_field("status", self.status.label())?;
        state.serialize_field("failing_row_index", &self.status.failing_row_index())?;
        state.serialize_field("statement", &self.status.statement())?;
        state.serialize_field(
            "parent_characteristic_polynomial",
            &rational_matrix_strings(&self.parent_characteristic_polynomial),
        )?;
        state.serialize_field(
            "parent_characteristic_polynomial_text",
            &format_shift_polynomial(&self.parent_characteristic_polynomial, "x"),
        )?;
        state.serialize_field(
            "section_characteristic_polynomial",
            &rational_matrix_strings(&self.section_characteristic_polynomial),
        )?;
        state.serialize_field(
            "section_characteristic_polynomial_text",
            &format_shift_polynomial(&self.section_characteristic_polynomial, "y"),
        )?;
        state.serialize_field("cofactor", &rational_matrix_strings(&self.cofactor))?;
        state.serialize_field(
            "cofactor_text",
            &format_shift_polynomial(&self.cofactor, "y"),
        )?;
        state.serialize_field("parent_recurrence", &self.parent_recurrence_json())?;
        state.serialize_field("parent_row_check", &self.parent_row_check)?;
        state.end()
    }
}

impl Serialize for VeroneseSectionLiftRejection {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("VeroneseSectionLiftRejection", 6)?;
        state.serialize_field("modulus", &self.modulus)?;
        state.serialize_field("residue", &self.residue)?;
        state.serialize_field(
            "parent_characteristic_polynomial",
            &rational_matrix_strings(&self.parent_characteristic_polynomial),
        )?;
        state.serialize_field(
            "parent_characteristic_polynomial_text",
            &format_shift_polynomial(&self.parent_characteristic_polynomial, "x"),
        )?;
        state.serialize_field("reason", self.reason.label())?;
        let row_index = match self.reason {
            VeroneseSectionLiftRejectionReason::RowMismatch { index } => Some(index),
            _ => None,
        };
        state.serialize_field("row_index", &row_index)?;
        state.end()
    }
}

impl Serialize for VeroneseSectionLiftSkippedModulus {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("VeroneseSectionLiftSkippedModulus", 2)?;
        state.serialize_field("modulus", &self.modulus)?;
        state.serialize_field("reason", &self.reason)?;
        state.end()
    }
}

impl Serialize for VeroneseSectionLiftReport {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        let mut state = serializer.serialize_struct("VeroneseSectionLiftReport", 7)?;
        state.serialize_field("schema", VERONESE_SECTION_LIFT_SCHEMA)?;
        state.serialize_field(
            "input_characteristic_polynomial",
            &rational_matrix_strings(&self.input_characteristic_polynomial),
        )?;
        state.serialize_field(
            "input_characteristic_polynomial_text",
            &format_shift_polynomial(&self.input_characteristic_polynomial, "y"),
        )?;
        state.serialize_field("lifts", &self.lifts)?;
        state.serialize_field("rejected_candidates", &self.rejected_candidates)?;
        state.serialize_field("skipped_moduli", &self.skipped_moduli)?;
        state.serialize_field("exhausted_budget", &self.exhausted_budget)?;
        state.end()
    }
}
