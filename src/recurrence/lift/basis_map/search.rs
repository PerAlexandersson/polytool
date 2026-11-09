//! [`find_basis_map_lifts`]: the backward lift through a library of basis
//! maps (module documentation of `basis_map`, "The backward lift").

use super::super::hadamard_factor::prepare::{prepare_lift_rule, PreparedLiftRule};
use super::certify::{
    certify_prepared, BasisMapLiftCandidate, BasisMapLiftCertification, BasisMapLiftOptions,
    MAX_BASIS_MAP_LIFT_CHECKED_ROWS,
};
use super::library::BasisMapLibrary;
use crate::recurrence::{
    find_recurrence_adaptive_rational_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome, Recurrence, RecurrenceError, RecurrenceJson,
};
use crate::vec_poly;
use crate::BigRational;
use serde::{Deserialize, Serialize};
use std::sync::Arc;

type Q = BigRational;

/// JSON schema tag of [`BasisMapLiftReport`].
pub const BASIS_MAP_LIFTS_SCHEMA: &str = "polytool.recurrence-basis-map-lifts.v1";

/// Options of [`find_basis_map_lifts`].
#[derive(Debug, Clone)]
pub struct BasisMapLiftSearchOptions {
    /// The maps; `None` builds [`BasisMapLibrary::oeis_catalog`] for this
    /// call (share a built library across calls).
    pub library: Option<Arc<BasisMapLibrary>>,
    /// Input rows expanded in each basis and used to fit the coefficient
    /// rule, between 6 and [`MAX_BASIS_MAP_LIFT_CHECKED_ROWS`].  Default 12.
    pub detection_rows: usize,
    /// Candidate budget of each coefficient-rule search
    /// ([`find_recurrence_adaptive_rational_with_budget`]).  Default 400.
    pub coefficient_search_candidates: usize,
    /// Largest order of a coefficient rule.  Default 3.
    pub max_coefficient_order: usize,
    /// Largest degree in `u` of a coefficient rule's coefficients.
    /// Default 2.
    pub max_coefficient_degree: usize,
    /// Largest degree in `n` of a coefficient rule's coefficients.
    /// Default 2.
    pub max_coefficient_index_degree: usize,
    /// Largest number of maps tried; reaching it sets
    /// [`BasisMapLiftReport::exhausted_budget`].  Default 512.
    pub max_maps: usize,
    /// Largest number of candidates sent to certification (each runs the
    /// identity certificate and the parent provers); reaching it sets
    /// [`BasisMapLiftReport::exhausted_budget`].  Default `None` (no bound).
    pub max_certifications: Option<usize>,
    /// Stop at the first proved lift.  Default true.
    pub stop_at_first_proof: bool,
    /// The options of each certification.
    pub certification: BasisMapLiftOptions,
}

impl Default for BasisMapLiftSearchOptions {
    fn default() -> Self {
        Self {
            library: None,
            detection_rows: 12,
            coefficient_search_candidates: 400,
            max_coefficient_order: 3,
            max_coefficient_degree: 2,
            max_coefficient_index_degree: 2,
            max_maps: 512,
            max_certifications: None,
            stop_at_first_proof: true,
            certification: BasisMapLiftOptions::default(),
        }
    }
}

/// Why a map whose basis spans the input rows gave no certified candidate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum BasisMapLiftRejection {
    /// A coefficient row `R_n` has a root outside the source interval of
    /// every conclusion of the map theorem (decided exactly): this map
    /// cannot prove the input (a near miss).
    RootOutsideSource {
        /// The map.
        basis: String,
        /// The absolute index of the first such row.
        index: usize,
    },
    /// No coefficient rule was found within the budget.
    NoCoefficientRule {
        /// The map.
        basis: String,
        /// Detail.
        detail: String,
    },
    /// A coefficient rule was found but cannot be used.
    CoefficientRuleUnusable {
        /// The map.
        basis: String,
        /// Why.
        detail: String,
    },
}

/// The result of [`find_basis_map_lifts`] (schema
/// `polytool.recurrence-basis-map-lifts.v1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapLiftReport {
    /// [`BASIS_MAP_LIFTS_SCHEMA`].
    pub schema: String,
    /// The input: rule, first index and identifying seeds.
    pub input: RecurrenceJson,
    /// Maps tried.
    pub maps_tried: usize,
    /// Maps whose basis does not span the input rows (or whose rows are
    /// not triangular up to the needed degree).
    pub maps_not_spanning: usize,
    /// The certified candidates (every verdict, proved first).
    pub lifts: Vec<BasisMapLiftCertification>,
    /// Maps that span the input rows but gave no candidate.
    pub rejections: Vec<BasisMapLiftRejection>,
    /// Whether [`BasisMapLiftSearchOptions::max_maps`] was reached.
    pub exhausted_budget: bool,
}

impl BasisMapLiftReport {
    /// The first proved lift, if any.
    pub fn proved(&self) -> Option<&BasisMapLiftCertification> {
        self.lifts.iter().find(|lift| lift.is_proved())
    }

    /// Replay every recorded certification
    /// ([`BasisMapLiftCertification::verify`]) and check that each is about
    /// the recorded input.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] for a mismatch, and the errors
    /// of the replays.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        if self.schema != BASIS_MAP_LIFTS_SCHEMA {
            return Err(RecurrenceError::VerificationFailed(
                "unsupported basis-map lifts schema".into(),
            ));
        }
        for lift in &self.lifts {
            if lift.input != self.input {
                return Err(RecurrenceError::VerificationFailed(
                    "a certification is about another input".into(),
                ));
            }
            lift.verify()?;
        }
        Ok(())
    }
}

/// Fit a derivative-free rule to the coefficient rows.
fn fit_coefficient_rule(
    rows: &[Vec<Q>],
    n0: usize,
    options: &BasisMapLiftSearchOptions,
) -> Result<Result<PreparedLiftRule, (bool, String)>, RecurrenceError> {
    let search = AdaptiveSearchOptions {
        max_rec_len: options.max_coefficient_order,
        max_var_deg: options.max_coefficient_degree,
        max_idx_deg: options.max_coefficient_index_degree,
        max_diff_deg: 0,
        max_denom_var_deg: 0,
        max_denom_idx_deg: options.max_coefficient_index_degree.max(1),
        ..AdaptiveSearchOptions::default()
    };
    let outcome = find_recurrence_adaptive_rational_with_budget(
        rows,
        &search,
        AdaptiveSearchBudget::limited(options.coefficient_search_candidates),
    );
    let found = match outcome {
        AdaptiveSearchOutcome::Found(found) => found,
        AdaptiveSearchOutcome::NoRecurrence(_) => {
            return Ok(Err((
                false,
                "no derivative-free rule in the search space".into(),
            )))
        }
        AdaptiveSearchOutcome::BudgetExhausted(_) => {
            return Ok(Err((false, "the rule search ran out of budget".into())))
        }
    };
    // The search indexes the first row as 1.
    let rule: Recurrence = found.recurrence.shift_index(1 - n0 as i64);
    match prepare_lift_rule(&rule, rows, n0, "coefficient family") {
        Ok(prepared) => Ok(Ok(prepared)),
        Err(error) => Ok(Err((true, error.to_string()))),
    }
}

/// Find basis-map lifts of the sequence defined by `recurrence` and `rows`
/// (`rows[0]` is `Q_(first_index)`, ascending coefficients): `Q_n =
/// T_P(R_n)` for a map `T_P` of the library and a coefficient family `R`
/// (module documentation of `basis_map`, "The backward lift").
///
/// For each map (library order): expand the detection rows in its basis
/// (exact triangular inverse; maps whose basis does not span the rows are
/// skipped), reject the map when an expanded row has a root outside every
/// source interval of the map theorem (exact; a near miss), fit a
/// derivative-free rule to the expanded rows (exact adaptive search), and
/// certify the candidate with [`certify_basis_map_lift`](super::certify_basis_map_lift).
/// Bounded and deterministic.  A missing lift says nothing against
/// real-rootedness of the input.
///
/// # Errors
///
/// [`RecurrenceError`] for an input rule that does not reproduce its rows,
/// is not homogeneous or has a leader not proved nonvanishing, and for
/// options beyond their bounds.
pub fn find_basis_map_lifts(
    recurrence: &Recurrence,
    rows: &[Vec<Q>],
    first_index: usize,
    options: &BasisMapLiftSearchOptions,
) -> Result<BasisMapLiftReport, RecurrenceError> {
    if !(6..=MAX_BASIS_MAP_LIFT_CHECKED_ROWS).contains(&options.detection_rows) {
        return Err(RecurrenceError::LimitExceeded(format!(
            "detection_rows must be between 6 and {MAX_BASIS_MAP_LIFT_CHECKED_ROWS}"
        )));
    }
    let input = prepare_lift_rule(recurrence, rows, first_index, "input")?;
    let library = match &options.library {
        Some(library) => library.clone(),
        None => Arc::new(BasisMapLibrary::oeis_catalog()?),
    };
    let detection = input.rows(options.detection_rows)?;
    let mut report = BasisMapLiftReport {
        schema: BASIS_MAP_LIFTS_SCHEMA.into(),
        input: RecurrenceJson::from_recurrence_rational(
            &input.rule,
            first_index,
            &input.seeds,
            None,
        ),
        maps_tried: 0,
        maps_not_spanning: 0,
        lifts: Vec::new(),
        rejections: Vec::new(),
        exhausted_budget: false,
    };
    for entry in library.entries() {
        if report.maps_tried >= options.max_maps
            || options
                .max_certifications
                .is_some_and(|bound| report.lifts.len() >= bound)
        {
            report.exhausted_budget = true;
            break;
        }
        report.maps_tried += 1;
        let name = entry.map.name().to_string();
        let expanded = match entry.map.expand_all(&detection) {
            Ok(expanded) => expanded,
            Err(RecurrenceError::Unsupported(_)) | Err(RecurrenceError::LimitExceeded(_)) => {
                report.maps_not_spanning += 1;
                continue;
            }
            Err(error) => return Err(error),
        };
        // Near-miss filter: some conclusion's source must contain the roots
        // of every nonzero expanded row.
        let fits = |row: &Vec<Q>, source: &super::quadratic::AlgebraicInterval| {
            vec_poly::is_zero(row) || source.contains_all_roots_of(row)
        };
        let admissible = entry
            .theorem
            .conclusions
            .iter()
            .any(|c| expanded.iter().all(|row| fits(row, &c.input_interval)));
        if !admissible {
            let index = expanded
                .iter()
                .position(|row| {
                    !entry
                        .theorem
                        .conclusions
                        .iter()
                        .any(|c| fits(row, &c.input_interval))
                })
                .unwrap_or(0);
            report
                .rejections
                .push(BasisMapLiftRejection::RootOutsideSource {
                    basis: name,
                    index: first_index + index,
                });
            continue;
        }
        let coefficient = match fit_coefficient_rule(&expanded, first_index, options)? {
            Ok(coefficient) => coefficient,
            Err((found, detail)) => {
                report.rejections.push(if found {
                    BasisMapLiftRejection::CoefficientRuleUnusable {
                        basis: name,
                        detail,
                    }
                } else {
                    BasisMapLiftRejection::NoCoefficientRule {
                        basis: name,
                        detail,
                    }
                });
                continue;
            }
        };
        let candidate = BasisMapLiftCandidate {
            basis: entry.map.clone(),
            coefficient_family: RecurrenceJson::from_recurrence_rational(
                &coefficient.rule,
                first_index,
                &coefficient.seeds,
                None,
            ),
        };
        let certification =
            certify_prepared(&input, &coefficient, &candidate, &options.certification)?;
        let proved = certification.is_proved();
        report.lifts.push(certification);
        if proved && options.stop_at_first_proof {
            break;
        }
    }
    report.lifts.sort_by_key(|lift| !lift.is_proved());
    Ok(report)
}
