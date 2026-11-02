//! [`apply_basis_map_to_family`]: the forward move `R -> T_P(R)` with a
//! rule for the image and its certificate (module documentation of
//! `basis_map`, "The forward move").

use super::super::hadamard_factor::prepare::{prepare_rule, PreparedRule};
use super::super::section_operators::OperatorMatrix;
use super::certify::{
    certify_prepared, BasisMapLiftCandidate, BasisMapLiftCertification, BasisMapLiftOptions,
};
use super::map::{apply_with_rows, BasisMap};
use super::transport::{transported_scalar_rule, transported_system};
use crate::recurrence::{
    find_recurrence_adaptive_rational_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome, ImplicitRelation, Recurrence, RecurrenceError, RecurrenceJson,
    MAX_IMAGE_RECURRENCE_SHAPES, MAX_IMAGE_RECURRENCE_UNKNOWNS,
};
use crate::BigRational;
use ore_algebra::{find_component_annihilator, ComponentAnnihilatorSearch, ImageAnnihilatorBudget};
use serde::{Deserialize, Serialize};

type Q = BigRational;

/// Largest number of ansatz shapes of the elimination.
const MAX_ELIMINATION_SHAPES: usize = MAX_IMAGE_RECURRENCE_SHAPES;

/// Largest number of unknowns of one ansatz of the elimination.
const MAX_ELIMINATION_UNKNOWNS: usize = MAX_IMAGE_RECURRENCE_UNKNOWNS;

/// JSON schema tag of [`BasisMapImage`].
pub const BASIS_MAP_IMAGE_SCHEMA: &str = "polytool.recurrence-basis-map-image.v1";

/// Options of [`apply_basis_map_to_family`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BasisMapImageOptions {
    /// Largest number of ansatz shapes solved when the rule of the first
    /// component of the transported system is computed (a map of order two
    /// or more).  Default 60.
    pub elimination_shapes: usize,
    /// Largest number of unknowns of one such ansatz.  Default 400.
    pub elimination_unknowns: usize,
    /// Rows of the image used to fit its rule when the computation finds
    /// nothing within its budget.  Default 20.
    pub fit_rows: usize,
    /// Candidate budget of that fit.  Default 600.
    pub fit_candidates: usize,
    /// Options of the certification of the image.
    pub certification: BasisMapLiftOptions,
}

impl Default for BasisMapImageOptions {
    fn default() -> Self {
        Self {
            elimination_shapes: 60,
            elimination_unknowns: 400,
            fit_rows: 20,
            fit_candidates: 600,
            certification: BasisMapLiftOptions::default(),
        }
    }
}

/// How the image's rule was obtained.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BasisMapImageConstruction {
    /// The transported rule `sum_i e_i(n, A) H_(n+i) = 0` (a map of order
    /// one, `T o u = A o T`): exact by construction.
    Transported,
    /// A map of order two or more: a runnable scalar rule of the first
    /// component `H_n = T(R_n)` of the transported system, computed as a
    /// common left multiple (`ore_algebra::find_component_annihilator`),
    /// then certified by the division of the module documentation.
    Eliminated,
    /// Fitted to the image rows (a map of order two or more, when the
    /// computation finds nothing within its budget), then certified by the
    /// division of the module documentation.
    Fitted,
}

/// The image `H_n = T_P(R_n)` of a coefficient family (schema
/// `polytool.recurrence-basis-map-image.v1`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BasisMapImage {
    /// [`BASIS_MAP_IMAGE_SCHEMA`].
    pub schema: String,
    /// The image: rule, first index `n0` and identifying seeds.
    pub image: RecurrenceJson,
    /// How the rule was obtained.
    pub construction: BasisMapImageConstruction,
    /// The certification of `H_n = T_P(R_n)` for every `n >= n0` and of the
    /// root transfer (with the coefficient family's root proof).
    pub certification: BasisMapLiftCertification,
}

impl BasisMapImage {
    /// Whether the identity is certified and the root transfer proved.
    pub fn is_proved(&self) -> bool {
        self.certification.is_proved()
    }
}

/// The rule of the first component of the transported system `Phi`
/// (order `k`, size `r`): a runnable scalar annihilator `M` with
/// `w(n) (M e_0) = sum_b H_b o Phi_b`, searched up to order `r k`, as the
/// successor rule at target index `m` (base `m - R`).  `None` when nothing
/// is found within the budget or `M` has no derivative-free leader.
fn eliminated_rule(
    system: &OperatorMatrix,
    order: usize,
    options: &BasisMapImageOptions,
) -> Result<Option<Recurrence>, RecurrenceError> {
    let budget = ImageAnnihilatorBudget {
        max_order_increase: (system.size() - 1) * order,
        max_shapes: options.elimination_shapes.min(MAX_ELIMINATION_SHAPES),
        max_unknowns: options.elimination_unknowns.min(MAX_ELIMINATION_UNKNOWNS),
        ..ImageAnnihilatorBudget::default()
    };
    let search = find_component_annihilator(&system.entries, 0, &budget)
        .map_err(|error| RecurrenceError::Unsupported(error.to_string()))?;
    let ComponentAnnihilatorSearch::Found { certificate, .. } = search else {
        return Ok(None);
    };
    Ok(ImplicitRelation::from_op(&certificate.annihilator)
        .to_successor_rule()
        .map(|(rule, _)| rule))
}

/// Apply the basis map to a coefficient family `R` (rule in `u` without
/// derivatives, first index `n0`, seeds): a rule with seeds for `H_n =
/// T_P(R_n)`, `n >= n0`, with its certificate, so that the image can enter
/// a proof database as a new family (module documentation of `basis_map`,
/// "The forward move").  When `R`'s roots are proved to lie in the source
/// interval of a conclusion of the map theorem, the certification's summary
/// is proved and every `H_n` is zero or real-rooted in its target.
///
/// # Errors
///
/// [`RecurrenceError`] when `R` cannot be replayed or validated, has
/// derivatives in `u`, when the map has no companion form, or when no rule
/// for the image is found or certified.
pub fn apply_basis_map_to_family(
    map: &BasisMap,
    family: &RecurrenceJson,
    options: &BasisMapImageOptions,
) -> Result<BasisMapImage, RecurrenceError> {
    let (rule, n0, seeds) = family.to_recurrence_parts()?;
    let rows = rule
        .generate_rows_rational(&seeds, n0, seeds.len() + 2)
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!("the family cannot be replayed: {error}"))
        })?;
    let coefficient = prepare_rule(&rule, &rows, n0, "coefficient family")?;
    let system = transported_system(
        &coefficient.relation,
        map,
        options.certification.max_operator_terms,
    )?;
    let image_rows = |count: usize| -> Result<Vec<Vec<Q>>, RecurrenceError> {
        let r_rows = coefficient.rows(count)?;
        let width = r_rows.iter().map(Vec::len).max().unwrap_or(1);
        let basis = map.basis_rows(width)?;
        Ok(r_rows.iter().map(|r| apply_with_rows(r, &basis)).collect())
    };
    let (image, construction): (PreparedRule, BasisMapImageConstruction) = if system.size() == 1 {
        let (image_rule, _) = transported_scalar_rule(&system)?;
        let rows = image_rows(coefficient.seeds.len() + 2)?;
        (
            prepare_rule(&image_rule, &rows, n0, "image")?,
            BasisMapImageConstruction::Transported,
        )
    } else if let Some(image_rule) = eliminated_rule(&system, coefficient.leader_shift, options)? {
        let rows =
            image_rows(coefficient.seeds.len() + 2 * system.size() * coefficient.leader_shift)?;
        (
            prepare_rule(&image_rule, &rows, n0, "image")?,
            BasisMapImageConstruction::Eliminated,
        )
    } else {
        let rows = image_rows(options.fit_rows.max(8))?;
        let search = AdaptiveSearchOptions {
            max_rec_len: 2 * system.size() * coefficient.leader_shift,
            max_var_deg: 4,
            max_idx_deg: 3,
            max_diff_deg: 0,
            max_denom_var_deg: 0,
            max_denom_idx_deg: 3,
            ..AdaptiveSearchOptions::default()
        };
        let found = match find_recurrence_adaptive_rational_with_budget(
            &rows,
            &search,
            AdaptiveSearchBudget::limited(options.fit_candidates),
        ) {
            AdaptiveSearchOutcome::Found(found) => found,
            _ => {
                return Err(RecurrenceError::Unsupported(
                    "no rule for the image rows was found within the budget".into(),
                ))
            }
        };
        let image_rule = found.recurrence.shift_index(1 - n0 as i64);
        (
            prepare_rule(&image_rule, &rows, n0, "image")?,
            BasisMapImageConstruction::Fitted,
        )
    };
    let candidate = BasisMapLiftCandidate {
        basis: map.clone(),
        coefficient_family: RecurrenceJson::from_recurrence_rational(
            &coefficient.rule,
            n0,
            &coefficient.seeds,
            None,
        ),
    };
    let certification = certify_prepared(&image, &coefficient, &candidate, &options.certification)?;
    if !matches!(
        certification.identity,
        super::certify::BasisMapIdentityStatus::Certified { .. }
    ) {
        return Err(RecurrenceError::VerificationFailed(format!(
            "the image rule is not certified: {:?}",
            certification.lift_summary.identity
        )));
    }
    Ok(BasisMapImage {
        schema: BASIS_MAP_IMAGE_SCHEMA.into(),
        image: RecurrenceJson::from_recurrence_rational(&image.rule, n0, &image.seeds, None),
        construction,
        certification,
    })
}
