//! The catalog-lift stage of proof discovery
//! ([`crate::recurrence::discover_recurrence_proofs`], step 4): when no
//! general theorem proves a record (directly, on its tail or on a certified
//! alternative representation), connect it to a family of the bundled OEIS
//! catalog ([`crate::oeis::catalog`]) whose real-rootedness is known, by an
//! identity certified for every index, and transfer the statement.  No
//! external database is used.
//!
//! # The lifts
//!
//! 1. **Normalized identification.**  A catalog family `S` whose rows equal
//!    the record's up to the normalization menu of
//!    [`certify_normalized_recurrence_equality`]: `T_n = c(n) x^(alpha n +
//!    beta) S_(n + s)` from some index on, with dropped leading rows.
//!    Candidates come from a fingerprint index (below); each candidate's
//!    equality is certified, the family `S` is proved by the general
//!    discovery stages (theorems, tail, alternative representations, never
//!    by another lift), and its statement "every row from `f` on is nonzero
//!    and real-rooted with every root in `R`" is carried to the record by
//!    [`transfer_real_rootedness_along_normalization`].  Rows of the record
//!    the identity or the premise does not cover are checked exactly.
//! 2. **Backward basis-map lift.**  [`find_basis_map_lifts`] with the
//!    library [`BasisMapLibrary::oeis_catalog`] (built once per process), on
//!    the record or, when discovery has a tail, on the tail record (whose
//!    rows are the family's rows from its first index on):
//!    `Q_n = T_P(R_n)` for every `n`, the map theorem sends polynomials with
//!    roots in a source `S` to ones with roots in a target `J`, and the
//!    coefficient family `R` is proved to have its roots in `S` by the
//!    shared parent provers (the registry theorems that discovery also
//!    uses).  When the identity is certified but the parent provers leave
//!    `R` open, `R` is identified with a catalog family as in 1; the
//!    transferred region of `R` must lie in the source of a conclusion of
//!    the map theorem, whose target is then the region of the record.
//!
//! # Polytool proofs and catalog attestations
//!
//! A parent statement is either **polytool's own** (the parent report is a
//! replayable discovery report, or the basis-map certification's parent
//! route is a theorem route) or a **catalog attestation**: the catalog flag
//! [`OeisRealRootednessStatus::LeanProved`] of the identified family (label
//! `oeis_catalog:<A-number>:lean_proved`; a proof made outside polytool),
//! or a basis-map certification whose coefficient family is settled by that
//! flag (parent route 10).  A lift that rests on an attestation is reported
//! with [`CatalogParentProofKind::CatalogAttestation`] and gives the
//! discovery status `conditional_on_attestation`, never `proved`; a lift
//! with a polytool parent is preferred whenever both exist.
//!
//! # Nonzero rows
//!
//! The discovery statement is "every row from `N` on is nonzero and
//! real-rooted".  For a normalized identification with a polytool parent
//! every covered row is a nonzero scalar times a power of `x` times a
//! nonzero parent row, and the uncovered rows are checked exactly.
//! Otherwise (an attested parent, which only states "zero or real-rooted",
//! and the basis-map routes) the record's own leading-coefficient
//! certificate ([`certify_leading_coefficient`]) proves every row from its
//! `holds_from` on nonzero, and the rows before are checked exactly; without
//! such a certificate the lift gives no conclusion.
//!
//! # The fingerprint index
//!
//! For every homogeneous catalog family (its displayed record: prefix rows
//! and initial rows as seeds, [`crate::oeis::OeisSequenceDefinition::recurrence_export_parts`]) the first
//! [`INDEX_ROWS`] rows are reduced to a *normalized fingerprint*: the power
//! of `x` dividing the row is removed, denominators are cleared and the
//! content divided out, and the leading coefficient made positive (rows of
//! degree zero carry no information and are skipped).  These fingerprints
//! are invariant under the scalar and the power of `x` of the normalization,
//! and an index shift only moves them.  The index is built once per process,
//! on first use.  A record's candidates are the families with at least
//! [`MINIMUM_INDEX_MATCHES`] of its first rows matching at one index shift;
//! families with the record's own rule (up to that shift) are skipped, since
//! their discovery is the record's own.  The index is only a hint: every
//! candidate is certified and proved from scratch.
//!
//! # Budgets
//!
//! [`CatalogLiftOptions::max_candidates`] bounds the number of catalog
//! families sent to certification (the record's and the coefficient
//! families' together; `0` disables the stage), and
//! [`CatalogLiftOptions::max_basis_maps`] the maps of the backward search.
//! The backward search runs only when the first 12 rows have degree at
//! most 24, certifies at most 16 candidates, and uses
//! [`BasisMapLibrary::oeis_catalog`] with one
//! geometric map `g (a + d t)^k` per pair (`g`, sign of `d`): the lifts
//! through the other geometric maps have coefficient families that differ
//! from the kept map's by an affine change of the variable, so they would
//! repeat its work (the catalog has dozens of such triangles).
//! Candidates are tried in a deterministic order (more matching rows first,
//! then catalog order).
//!
//! # Replay
//!
//! [`CatalogLiftSearch::verify`] replays every certificate: each normalized
//! equality, each parent discovery report (which must have no catalog-lift
//! stage), each transfer and its premise (a proved conclusion of the parent
//! report, or the catalog flag of the named entry, whose record must be the
//! certificate's source), the basis-map search with every certification,
//! and the coefficient identifications against their lifts; then the
//! conclusion is recomputed from the recorded attempts and must equal the
//! recorded one.

use super::normalized_equality::algebraic_region;
use super::*;
use crate::oeis::{catalog, OeisRealRootednessStatus};
use crate::recurrence::{
    certify_leading_coefficient, LeadingCoefficientCertificate, LeadingCoefficientOptions,
};
use crate::vec_poly;
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Signed, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::sync::{Arc, OnceLock};

type Q = BigRational;

/// Rows of each catalog family (and of a record) reduced to fingerprints.
pub const INDEX_ROWS: usize = 16;

/// Least number of a record's first rows that must match a catalog family
/// at one index shift for it to be a candidate.
pub const MINIMUM_INDEX_MATCHES: usize = 6;

/// Largest number of rows a conclusion checks exactly for being nonzero.
const MAX_NONZERO_CHECKED_ROWS: usize = 64;

/// Largest number of basis-map candidates certified by the stage.
const MAX_BASIS_MAP_CERTIFICATIONS: usize = 16;

/// Largest degree of a detection row for which the stage runs the
/// backward basis-map search (exact root location of the expanded rows is
/// the cost; rows of the catalog triangles have degree at most 12).
const MAX_BASIS_MAP_ROW_DEGREE: usize = 24;

/// Rows generated for the backward basis-map search.
const BASIS_MAP_INPUT_ROWS: usize = 16;

/// Bounds of the catalog-lift stage of
/// [`crate::recurrence::discover_recurrence_proofs_with_stages`].
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogLiftOptions {
    /// Largest number of catalog families sent to the normalized-equality
    /// certificate (candidates for the record and for open basis-map
    /// coefficient families together, in a deterministic order).  `0`
    /// disables the whole stage.  Default 6.
    pub max_candidates: usize,
    /// Largest number of library maps tried by the backward basis-map lift
    /// (library order; `0` skips the basis-map lift).  Default 512.
    pub max_basis_maps: usize,
}

impl Default for CatalogLiftOptions {
    fn default() -> Self {
        Self {
            max_candidates: 6,
            max_basis_maps: 512,
        }
    }
}

impl CatalogLiftOptions {
    /// The options that disable the stage (used for the parents of a lift,
    /// which are proved without lifts).
    pub fn disabled() -> Self {
        Self {
            max_candidates: 0,
            max_basis_maps: 0,
        }
    }

    /// Whether the stage runs (`max_candidates > 0`).
    pub fn is_enabled(&self) -> bool {
        self.max_candidates > 0
    }
}

/// How the parent statement of a lift is known.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CatalogParentProofKind {
    /// polytool proves it (a replayed discovery report or theorem route).
    Polytool,
    /// The catalog's `lean_proved` flag asserts it (a proof made outside
    /// polytool): a lift resting on it is conditional.
    CatalogAttestation,
}

/// The route of a [`CatalogLiftConclusion`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CatalogLiftRoute {
    /// The record equals a catalog family up to a certified normalization.
    NormalizedIdentification,
    /// A backward basis-map lift whose certification settles the
    /// coefficient family itself.
    BasisMap,
    /// A backward basis-map lift whose coefficient family is identified
    /// with a catalog family (the certification's parent provers leave it
    /// open).
    BasisMapWithIdentifiedCoefficients,
}

/// One catalog family tried as a normalized identification (module
/// documentation, lift 1).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogIdentification {
    /// The catalog A-number.
    pub catalog_id: String,
    /// Rows of the identified record whose fingerprints matched the family
    /// in the index (a search hint, not evidence).
    pub matched_rows: usize,
    /// The normalized equality: source the catalog family (its displayed
    /// record), target the identified record.
    pub equality: NormalizedRecurrenceEqualityCertificate,
    /// Proof discovery on the catalog family with the catalog-lift stage
    /// disabled (depth one); present when the equality is certified and the
    /// discovery ran.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent: Option<Box<RecurrenceProofDiscoveryReport>>,
    /// How the transferred parent statement is known; omitted when nothing
    /// was transferred.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent_proof: Option<CatalogParentProofKind>,
    /// The attestation label `oeis_catalog:<A-number>:lean_proved`, for a
    /// [`CatalogParentProofKind::CatalogAttestation`] parent.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub attestation: Option<String>,
    /// The parent statement carried to the identified record.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transfer: Option<NormalizedRealRootednessTransfer>,
    /// Why nothing was transferred (when `transfer` is omitted).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

impl CatalogIdentification {
    /// Whether a statement proved by polytool was transferred.
    pub fn is_proved(&self) -> bool {
        self.transfer.is_some() && self.parent_proof == Some(CatalogParentProofKind::Polytool)
    }
}

/// A catalog identification of the coefficient family of one basis-map
/// lift (module documentation, lift 2).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogCoefficientIdentification {
    /// The position of the lift in [`CatalogBasisMapStage::search`]'s
    /// `lifts`.
    pub lift: usize,
    /// The identification; its target is the lift's coefficient family.
    pub identification: CatalogIdentification,
    /// The first conclusion of the map theorem whose source contains the
    /// transferred region of the coefficient family; omitted when none does.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub map_conclusion: Option<BasisMapConclusion>,
}

/// The backward basis-map part of the stage.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogBasisMapStage {
    /// The backward search on the record; omitted when the record is
    /// outside its input class (see `error`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub search: Option<BasisMapLiftReport>,
    /// Why the search did not run.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    /// Catalog identifications of open coefficient families, in the order
    /// tried.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub coefficient_identifications: Vec<CatalogCoefficientIdentification>,
}

/// What the catalog-lift stage establishes for the record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogLiftConclusion {
    /// The lift.
    pub route: CatalogLiftRoute,
    /// The identified catalog family (identification routes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub catalog_family: Option<String>,
    /// The basis map's name (basis-map routes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub basis_map: Option<String>,
    /// How the parent statement is known.
    pub parent_proof: CatalogParentProofKind,
    /// The attestation labels the conclusion is conditional on (empty for
    /// a polytool parent).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub conditional_on: Vec<String>,
    /// Every row `P_n`, `n >= from_index`, is nonzero and real-rooted with
    /// every root in `root_region` (conditional on `conditional_on`).
    pub from_index: usize,
    /// The root region.
    pub root_region: RootInterval,
    /// The record's leading-coefficient certificate behind the nonzero rows
    /// (attested parents and basis-map routes).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub leading_coefficient: Option<Box<LeadingCoefficientCertificate>>,
    /// Rows checked exactly to be nonzero (not covered by the identity or
    /// the leading-coefficient certificate), ascending.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub checked_rows: Vec<usize>,
    /// The conclusion and the parent, in words.
    pub statement: String,
}

/// The catalog-lift stage (module documentation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogLiftSearch {
    /// The bounds used.
    pub options: CatalogLiftOptions,
    /// The normalized identifications of the record, in the order tried;
    /// the search stops at the first one with a polytool parent.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub identifications: Vec<CatalogIdentification>,
    /// The backward basis-map lift; omitted when not run (an identification
    /// proved the record, or `max_basis_maps` is 0).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub basis_maps: Option<CatalogBasisMapStage>,
    /// The conclusion: the first lift with a polytool parent, else the
    /// first conditional one; omitted when no lift concludes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub conclusion: Option<CatalogLiftConclusion>,
}

impl CatalogLiftSearch {
    /// Whether the conclusion is proved by polytool.
    pub fn is_proved(&self) -> bool {
        self.conclusion
            .as_ref()
            .is_some_and(|c| c.parent_proof == CatalogParentProofKind::Polytool)
    }

    /// Number of catalog families sent to certification.
    pub fn candidates_tried(&self) -> usize {
        self.identifications.len()
            + self
                .basis_maps
                .as_ref()
                .map_or(0, |stage| stage.coefficient_identifications.len())
    }

    /// Replay the stage for the record `record` (module documentation,
    /// "Replay"); `basis_record` is the record of the backward basis-map
    /// search (the tail record when discovery has a tail, else `record`).
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] on any mismatch, and the
    /// errors of the replayed certificates.
    pub fn verify(
        &self,
        record: &RecurrenceJson,
        basis_record: &RecurrenceJson,
        options: &RealRootednessTheoremOptions,
    ) -> Result<(), RecurrenceError> {
        let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
        if !self.options.is_enabled() {
            return failed("a catalog-lift stage is recorded with the stage disabled");
        }
        if self.candidates_tried() > self.options.max_candidates {
            return failed("more catalog candidates than the recorded bound");
        }
        for (position, identification) in self.identifications.iter().enumerate() {
            verify_identification(identification, record, options)?;
            if identification.is_proved() && position + 1 != self.identifications.len() {
                return failed("the identification search continues after a proof");
            }
        }
        if let Some(stage) = &self.basis_maps {
            if self.options.max_basis_maps == 0 {
                return failed("a basis-map lift is recorded with max_basis_maps 0");
            }
            if self
                .identifications
                .iter()
                .any(CatalogIdentification::is_proved)
            {
                return failed("a basis-map lift is recorded after a proved identification");
            }
            verify_basis_map_stage(stage, basis_record, &self.options, options)?;
        }
        let expected = conclude(record, basis_record, self)?;
        if expected != self.conclusion {
            return failed("the catalog-lift conclusion does not follow from its attempts");
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------------

/// Run the stage on `record` (a plain record: rule, first index, seeds);
/// the backward basis-map search runs on `basis_record` (the tail record
/// when discovery has a tail, whose rows are the family's rows from its
/// first index on, else `record`).
pub(in crate::recurrence) fn catalog_lift_search(
    record: &RecurrenceJson,
    basis_record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
    alternative: &AlternativeRepresentationOptions,
    lift: &CatalogLiftOptions,
) -> Result<CatalogLiftSearch, RecurrenceError> {
    let mut search = CatalogLiftSearch {
        options: lift.clone(),
        identifications: Vec::new(),
        basis_maps: None,
        conclusion: None,
    };
    let mut budget = lift.max_candidates;
    let index = catalog_index();
    for candidate in index.candidates(record) {
        if budget == 0 {
            break;
        }
        budget -= 1;
        // A candidate whose certification fails with an error is skipped
        // (the stage never makes discovery fail).
        let Ok(identification) = identify(record, index, &candidate, options, alternative, None)
        else {
            continue;
        };
        let proved = identification.is_proved();
        search.identifications.push(identification);
        if proved {
            break;
        }
    }
    if !search
        .identifications
        .iter()
        .any(CatalogIdentification::is_proved)
        && lift.max_basis_maps > 0
    {
        search.basis_maps = Some(basis_map_stage(
            basis_record,
            index,
            options,
            alternative,
            lift,
            &mut budget,
        )?);
    }
    search.conclusion = conclude(record, basis_record, &search)?;
    Ok(search)
}

/// The backward basis-map lift and the identification of open coefficient
/// families.
fn basis_map_stage(
    record: &RecurrenceJson,
    index: &CatalogIndex,
    options: &RealRootednessTheoremOptions,
    alternative: &AlternativeRepresentationOptions,
    lift: &CatalogLiftOptions,
    budget: &mut usize,
) -> Result<CatalogBasisMapStage, RecurrenceError> {
    let mut stage = CatalogBasisMapStage {
        search: None,
        error: None,
        coefficient_identifications: Vec::new(),
    };
    let search = match backward_search(record, lift) {
        Ok(search) => search,
        Err(error) => {
            stage.error = Some(error.to_string());
            return Ok(stage);
        }
    };
    if search.proved().is_none() {
        'lifts: for (position, certification) in search.lifts.iter().enumerate() {
            if !has_open_coefficient_family(certification) {
                continue;
            }
            let family = &certification.candidate.coefficient_family;
            let sources: Vec<&AlgebraicInterval> = certification
                .map_theorem
                .conclusions
                .iter()
                .map(|c| &c.input_interval)
                .collect();
            for candidate in index.candidates(family) {
                if *budget == 0 {
                    break 'lifts;
                }
                *budget -= 1;
                let Ok(identification) = identify(
                    family,
                    index,
                    &candidate,
                    options,
                    alternative,
                    Some(&sources),
                ) else {
                    continue;
                };
                let map_conclusion =
                    fitting_conclusion(certification, &identification).unwrap_or(None);
                let proved = identification.is_proved() && map_conclusion.is_some();
                stage
                    .coefficient_identifications
                    .push(CatalogCoefficientIdentification {
                        lift: position,
                        identification,
                        map_conclusion,
                    });
                if proved {
                    break 'lifts;
                }
            }
        }
    }
    stage.search = Some(search);
    Ok(stage)
}

/// [`find_basis_map_lifts`] on the record with the shared library.
fn backward_search(
    record: &RecurrenceJson,
    lift: &CatalogLiftOptions,
) -> Result<BasisMapLiftReport, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let rows = basis_map_input_rows(&rule, first, &seeds)?;
    let options = BasisMapLiftSearchOptions::default();
    let largest = rows
        .iter()
        .take(options.detection_rows)
        .filter_map(|row| vec_poly::degree(&vec_poly::trim(row)))
        .max()
        .unwrap_or(0);
    if largest > MAX_BASIS_MAP_ROW_DEGREE {
        return Err(RecurrenceError::LimitExceeded(format!(
            "a detection row has degree {largest}, above the stage's bound \
             {MAX_BASIS_MAP_ROW_DEGREE} for the backward basis-map search"
        )));
    }
    find_basis_map_lifts(
        &rule,
        &rows,
        first,
        &BasisMapLiftSearchOptions {
            library: Some(shared_basis_map_library()?),
            max_maps: lift.max_basis_maps,
            max_certifications: Some(MAX_BASIS_MAP_CERTIFICATIONS),
            ..options
        },
    )
}

fn basis_map_input_rows(
    rule: &Recurrence,
    first: usize,
    seeds: &[Vec<Q>],
) -> Result<Vec<Vec<Q>>, RecurrenceError> {
    rule.generate_rows_rational(seeds, first, BASIS_MAP_INPUT_ROWS.max(seeds.len()))
        .map_err(|error| {
            RecurrenceError::InvalidInput(format!("the rows cannot be generated: {error}"))
        })
}

/// The discovery library (module documentation, "Budgets"), built once per
/// process.
fn shared_basis_map_library() -> Result<Arc<BasisMapLibrary>, RecurrenceError> {
    static LIBRARY: OnceLock<Result<Arc<BasisMapLibrary>, String>> = OnceLock::new();
    LIBRARY
        .get_or_init(|| {
            discovery_basis_map_library()
                .map(Arc::new)
                .map_err(|error| error.to_string())
        })
        .clone()
        .map_err(RecurrenceError::Unsupported)
}

/// [`BasisMapLibrary::oeis_catalog`] with one geometric map `g (a + d t)^k`
/// per pair (`g`, sign of `d`): the backward lifts through the others give
/// the coefficient families `Q_n((u - a) / d) / g`, an affine change of the
/// variable away from the kept map's, so they repeat its work (a search
/// economy; nothing is claimed about the maps left out).
fn discovery_basis_map_library() -> Result<BasisMapLibrary, RecurrenceError> {
    let mut library = BasisMapLibrary::oeis_catalog()?;
    let mut seen = Vec::new();
    library.retain_entries(|entry| match &entry.theorem.shape {
        BasisMapShape::Geometric {
            initial_row, slope, ..
        } => {
            let key = (initial_row.clone(), slope.starts_with('-'));
            if seen.contains(&key) {
                false
            } else {
                seen.push(key);
                true
            }
        }
        _ => true,
    });
    Ok(library)
}

/// A certified basis-map lift whose only missing piece is the coefficient
/// family's root location (open, or settled only by an attestation).
fn has_open_coefficient_family(certification: &BasisMapLiftCertification) -> bool {
    let summary = &certification.lift_summary;
    matches!(
        certification.identity,
        BasisMapIdentityStatus::Certified { .. }
            | BasisMapIdentityStatus::RecurrenceEquality { .. }
    ) && summary
        .hypotheses
        .iter()
        .all(|h| h.status == LiftHypothesisStatus::Verified)
        && matches!(
            summary.parent_proof,
            LiftParentProof::Open { .. } | LiftParentProof::Attested { .. }
        )
        && !certification.map_theorem.conclusions.is_empty()
}

/// The first conclusion of the map theorem whose source contains the
/// region transferred to the coefficient family.
fn fitting_conclusion(
    certification: &BasisMapLiftCertification,
    identification: &CatalogIdentification,
) -> Result<Option<BasisMapConclusion>, RecurrenceError> {
    let Some(transfer) = &identification.transfer else {
        return Ok(None);
    };
    let region = algebraic_region(&transfer.conclusion_region)?;
    Ok(certification
        .map_theorem
        .conclusions
        .iter()
        .find(|c| region.is_subset_of(&c.input_interval))
        .cloned())
}

/// The whole real line as a region.
fn whole_line() -> RootInterval {
    RootInterval {
        lower: None,
        upper: None,
    }
}

/// The attestation label of a flagged catalog entry.
fn attestation_label(id: &str) -> String {
    format!("oeis_catalog:{id}:lean_proved")
}

/// A statement "every row from `from_index` on is nonzero and real-rooted
/// with every root in `region`" given by a proved conclusion of a discovery
/// report.
#[derive(Debug, Clone, PartialEq, Eq)]
struct ProvedClaim {
    from_index: usize,
    region: RootInterval,
    theorem: &'static str,
}

/// The statements given by the proved conclusions of a discovery report
/// without zero residue classes (on the record, its tail and certified
/// alternatives): each with its certified region and with the whole line,
/// the smallest first index first (stable).
fn proved_claims(report: &RecurrenceProofDiscoveryReport) -> Vec<ProvedClaim> {
    let mut attempts: Vec<&RealRootednessTheoremReport> = report.attempts.iter().collect();
    if let Some(tail) = &report.tail {
        attempts.extend(&tail.attempts);
    }
    if let Some(search) = &report.alternative_representations {
        for alternative in search.alternatives.iter().filter(|a| a.proved) {
            attempts.extend(&alternative.attempts);
        }
    }
    let mut claims: Vec<ProvedClaim> = Vec::new();
    for attempt in attempts.into_iter().filter(|attempt| attempt.is_proved()) {
        let Some(conclusion) = attempt.conclusion.as_ref() else {
            continue;
        };
        if conclusion.has_zero_residue_classes() {
            continue;
        }
        // The certified region, then the whole line (weaker, and usable when
        // the region depends on the index).
        let mut regions = vec![whole_line()];
        if let Some(region) = &conclusion.root_region {
            if *region != whole_line() {
                regions.insert(0, region.clone());
            }
        }
        for region in regions {
            let claim = ProvedClaim {
                from_index: conclusion.first_index,
                region,
                theorem: attempt.theorem.as_str(),
            };
            if !claims
                .iter()
                .any(|c| c.from_index == claim.from_index && c.region == claim.region)
            {
                claims.push(claim);
            }
        }
    }
    claims.sort_by_key(|claim| claim.from_index);
    claims
}

/// Certify one candidate, prove its family and transfer a statement
/// (module documentation, lift 1).  With `wanted`, the transferred region
/// must lie in one of those intervals.
fn identify(
    target: &RecurrenceJson,
    index: &CatalogIndex,
    candidate: &IndexCandidate,
    options: &RealRootednessTheoremOptions,
    alternative: &AlternativeRepresentationOptions,
    wanted: Option<&[&AlgebraicInterval]>,
) -> Result<CatalogIdentification, RecurrenceError> {
    let entry = &index.entries[candidate.entry];
    let normalization_options = NormalizedEqualityOptions {
        max_index_shift: (candidate.shift.unsigned_abs() as usize).clamp(4, 32),
        max_dropped_rows: candidate.first_matched_offset.clamp(3, 12),
        ..NormalizedEqualityOptions::default()
    };
    let equality =
        certify_normalized_recurrence_equality(&entry.record, target, &normalization_options)?;
    let mut identification = CatalogIdentification {
        catalog_id: entry.id.to_string(),
        matched_rows: candidate.matched_rows,
        equality,
        parent: None,
        parent_proof: None,
        attestation: None,
        transfer: None,
        reason: None,
    };
    if !identification.equality.is_certified() {
        identification.reason = Some("the normalized equality is not certified".into());
        return Ok(identification);
    }
    let fits = |transfer: &NormalizedRealRootednessTransfer| -> bool {
        wanted.is_none_or(|intervals| {
            algebraic_region(&transfer.conclusion_region)
                .is_ok_and(|region| intervals.iter().any(|w| region.is_subset_of(w)))
        })
    };
    match crate::recurrence::discover_recurrence_proofs_with_stage_options(
        &entry.record,
        options,
        &crate::recurrence::RecurrenceProofStageOptions {
            alternative_representations: alternative.clone(),
            catalog_lifts: CatalogLiftOptions::disabled(),
            affine_coordinates: crate::recurrence::AffineCoordinateOptions::disabled(),
        },
    ) {
        Ok(parent) => {
            for claim in proved_claims(&parent) {
                let Ok(transfer) = transfer_real_rootedness_along_normalization(
                    &identification.equality,
                    NormalizationTransferDirection::SourceToTarget,
                    claim.from_index,
                    &claim.region,
                ) else {
                    continue;
                };
                if fits(&transfer) {
                    identification.transfer = Some(transfer);
                    identification.parent_proof = Some(CatalogParentProofKind::Polytool);
                    break;
                }
            }
            identification.parent = Some(Box::new(parent));
        }
        Err(error) => {
            identification.reason = Some(format!(
                "proof discovery on the catalog family failed: {error}"
            ));
        }
    }
    if identification.transfer.is_none() && entry.flagged {
        if let Ok(transfer) = transfer_real_rootedness_along_normalization(
            &identification.equality,
            NormalizationTransferDirection::SourceToTarget,
            entry.record.first_index,
            &whole_line(),
        ) {
            if fits(&transfer) {
                identification.transfer = Some(transfer);
                identification.parent_proof = Some(CatalogParentProofKind::CatalogAttestation);
                identification.attestation = Some(attestation_label(entry.id));
            }
        }
    }
    if identification.transfer.is_none() && identification.reason.is_none() {
        identification.reason = Some(match wanted {
            None => "polytool does not prove the catalog family and it carries no lean_proved \
                     flag"
                .into(),
            Some(_) => "no proved or attested statement of the catalog family transfers to a \
                        region inside a source of the map theorem"
                .into(),
        });
    }
    Ok(identification)
}

// ---------------------------------------------------------------------------
// Replay
// ---------------------------------------------------------------------------

/// Replay one identification whose target is `target`.
fn verify_identification(
    identification: &CatalogIdentification,
    target: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
) -> Result<(), RecurrenceError> {
    let failed = |what: String| Err(RecurrenceError::VerificationFailed(what));
    let id = &identification.catalog_id;
    let index = catalog_index();
    let Some(entry) = index.entries.iter().find(|e| e.id == id) else {
        return failed(format!("{id} is not an indexed catalog family"));
    };
    if identification.equality.source != entry.record {
        return failed(format!(
            "the identification's source is not the record of {id}"
        ));
    }
    if identification.equality.target != *target {
        return failed(format!(
            "the identification with {id} is about another record"
        ));
    }
    identification.equality.verify()?;
    if let Some(parent) = &identification.parent {
        if !identification.equality.is_certified() {
            return failed(format!(
                "a parent of {id} is recorded without a certified identity"
            ));
        }
        if parent.recurrence != entry.record || parent.catalog_lifts.is_some() {
            return failed(format!(
                "the parent report of {id} is not a depth-one discovery on its record"
            ));
        }
        if parent.maximum_start_delay != options.maximum_start_delay
            || parent.requested_root_region.is_some()
        {
            return failed(format!("the parent report of {id} has other options"));
        }
        parent.verify()?;
    }
    match (&identification.transfer, identification.parent_proof) {
        (None, None) => {
            if identification.attestation.is_some() {
                return failed(format!("an attestation of {id} without a transfer"));
            }
        }
        (Some(transfer), Some(kind)) => {
            transfer.verify()?;
            if transfer.equality != identification.equality
                || transfer.direction != NormalizationTransferDirection::SourceToTarget
            {
                return failed(format!("the transfer from {id} uses another identity"));
            }
            match kind {
                CatalogParentProofKind::Polytool => {
                    let supported = identification.attestation.is_none()
                        && identification.parent.as_ref().is_some_and(|parent| {
                            proved_claims(parent).iter().any(|claim| {
                                claim.from_index == transfer.premise_from_index
                                    && claim.region == transfer.premise_region
                            })
                        });
                    if !supported {
                        return failed(format!(
                            "the premise transferred from {id} is not a proved conclusion of \
                             its parent report"
                        ));
                    }
                }
                CatalogParentProofKind::CatalogAttestation => {
                    if !entry.flagged
                        || identification.attestation.as_deref()
                            != Some(attestation_label(entry.id).as_str())
                        || transfer.premise_from_index != entry.record.first_index
                        || transfer.premise_region != whole_line()
                    {
                        return failed(format!(
                            "the attested premise of {id} is not its catalog flag"
                        ));
                    }
                }
            }
        }
        _ => {
            return failed(format!(
                "the transfer and the parent proof of {id} disagree"
            ))
        }
    }
    Ok(())
}

/// Replay the basis-map part.
fn verify_basis_map_stage(
    stage: &CatalogBasisMapStage,
    record: &RecurrenceJson,
    lift: &CatalogLiftOptions,
    options: &RealRootednessTheoremOptions,
) -> Result<(), RecurrenceError> {
    let failed = |what: &str| Err(RecurrenceError::VerificationFailed(what.to_string()));
    let Some(search) = &stage.search else {
        if !stage.coefficient_identifications.is_empty() {
            return failed("coefficient identifications without a basis-map search");
        }
        return Ok(());
    };
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let rows = basis_map_input_rows(&rule, first, &seeds)?;
    let input = prepare_lift_rule(&rule, &rows, first, "input")?;
    let expected = RecurrenceJson::from_recurrence_rational(&input.rule, first, &input.seeds, None);
    if search.input != expected {
        return failed("the basis-map search is about another record");
    }
    if search.maps_tried > lift.max_basis_maps {
        return failed("the basis-map search tried more maps than the recorded bound");
    }
    search.verify()?;
    if !stage.coefficient_identifications.is_empty() && search.proved().is_some() {
        return failed("coefficient identifications after a proved basis-map lift");
    }
    for (position, coefficient) in stage.coefficient_identifications.iter().enumerate() {
        let Some(certification) = search.lifts.get(coefficient.lift) else {
            return failed("a coefficient identification refers to no recorded lift");
        };
        if !has_open_coefficient_family(certification) {
            return failed("a coefficient identification of a lift that is not open");
        }
        verify_identification(
            &coefficient.identification,
            &certification.candidate.coefficient_family,
            options,
        )?;
        if fitting_conclusion(certification, &coefficient.identification).unwrap_or(None)
            != coefficient.map_conclusion
        {
            return failed(
                "the map conclusion of a coefficient identification is not the first \
                           fitting one",
            );
        }
        let proved = coefficient.identification.is_proved() && coefficient.map_conclusion.is_some();
        if proved && position + 1 != stage.coefficient_identifications.len() {
            return failed("the coefficient identifications continue after a proof");
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// The conclusion
// ---------------------------------------------------------------------------

/// A lift's statement before the nonzero rows are settled.
struct Claim {
    route: CatalogLiftRoute,
    catalog_family: Option<String>,
    basis_map: Option<String>,
    parent_proof: CatalogParentProofKind,
    conditional_on: Vec<String>,
    /// Every row from here on is zero or real-rooted in `region`.
    from_index: usize,
    region: RootInterval,
    /// Rows from here on are nonzero by the identity (normalized
    /// identification with a polytool parent); otherwise the record's
    /// leading-coefficient certificate is used.
    nonzero_from: Option<usize>,
    /// Whether the claim is about the basis-map record (else the record).
    on_basis_record: bool,
    description: String,
}

/// The statement of an identification, when one was transferred.
fn identification_claim(identification: &CatalogIdentification) -> Option<Claim> {
    let transfer = identification.transfer.as_ref()?;
    let kind = identification.parent_proof?;
    let id = &identification.catalog_id;
    let normalization = identification
        .equality
        .normalization
        .as_ref()
        .map(RecurrenceNormalization::describe)
        .unwrap_or_default();
    let (parent, conditional_on, nonzero_from) = match kind {
        CatalogParentProofKind::Polytool => {
            let theorem = identification
                .parent
                .as_ref()
                .and_then(|parent| {
                    proved_claims(parent).into_iter().find(|claim| {
                        claim.from_index == transfer.premise_from_index
                            && claim.region == transfer.premise_region
                    })
                })
                .map_or("", |claim| claim.theorem);
            // Rows from the first row past the exactly checked ones are
            // covered by the identity and the nonzero parent rows.
            let covered = transfer
                .checked_rows
                .last()
                .map_or(transfer.conclusion_from_index, |last| last + 1);
            (
                format!(
                    "polytool proves every row of {id} from index {} on nonzero and real-rooted \
                     ({theorem})",
                    transfer.premise_from_index
                ),
                Vec::new(),
                Some(covered),
            )
        }
        CatalogParentProofKind::CatalogAttestation => (
            format!(
                "the catalog flags every row of {id} as zero or real-rooted (lean_proved: a \
                 proof made outside polytool)"
            ),
            identification.attestation.iter().cloned().collect(),
            None,
        ),
    };
    Some(Claim {
        route: CatalogLiftRoute::NormalizedIdentification,
        catalog_family: Some(id.clone()),
        basis_map: None,
        parent_proof: kind,
        conditional_on,
        from_index: transfer.conclusion_from_index,
        region: transfer.conclusion_region.clone(),
        nonzero_from,
        on_basis_record: false,
        description: format!(
            "the rows equal those of the catalog family {id} = S up to the certified \
             normalization {normalization}, and {parent}"
        ),
    })
}

/// The statement of a basis-map certification that is proved or
/// conditional on an attestation.
fn basis_map_claim(certification: &BasisMapLiftCertification) -> Option<Claim> {
    let summary = &certification.lift_summary;
    let (kind, conditional_on) = match (&summary.verdict, &summary.parent_proof) {
        (LiftVerdict::Proved, _) => (CatalogParentProofKind::Polytool, Vec::new()),
        (LiftVerdict::ConditionalOnAttestation, LiftParentProof::Attested { labels, .. }) => {
            (CatalogParentProofKind::CatalogAttestation, labels.clone())
        }
        _ => return None,
    };
    let name = certification.candidate.basis.name().to_string();
    Some(Claim {
        route: CatalogLiftRoute::BasisMap,
        catalog_family: None,
        basis_map: Some(name.clone()),
        parent_proof: kind,
        conditional_on,
        from_index: summary.input_first_index,
        region: summary.root_transfer.input_root_region.clone(),
        nonzero_from: None,
        on_basis_record: true,
        description: format!(
            "Q_n = T(R_n) for the basis map T of {name}, certified for every index; {}; the \
             coefficient family R: {}",
            summary.root_transfer.theorem,
            summary.parent_proof.statement()
        ),
    })
}

/// The statement of a coefficient identification whose region fits.
fn coefficient_claim(
    certification: &BasisMapLiftCertification,
    coefficient: &CatalogCoefficientIdentification,
) -> Option<Claim> {
    let conclusion = coefficient.map_conclusion.as_ref()?;
    let inner = identification_claim(&coefficient.identification)?;
    let name = certification.candidate.basis.name().to_string();
    Some(Claim {
        route: CatalogLiftRoute::BasisMapWithIdentifiedCoefficients,
        catalog_family: inner.catalog_family,
        basis_map: Some(name.clone()),
        parent_proof: inner.parent_proof,
        conditional_on: inner.conditional_on,
        from_index: certification
            .lift_summary
            .input_first_index
            .max(inner.from_index),
        region: conclusion.output_interval.to_root_interval(),
        nonzero_from: None,
        on_basis_record: true,
        description: format!(
            "Q_n = T(R_n) for the basis map T of {name}, certified for every index; the \
             coefficient family R is zero or real-rooted with every root in {} from index {} on, \
             since {}; and {}",
            format_root_interval_exact(&inner.region),
            inner.from_index,
            inner.description,
            conclusion.statement
        ),
    })
}

/// The conclusion of the recorded attempts: the first claim with a
/// polytool parent whose nonzero rows are settled, else the first such
/// conditional claim.
fn conclude(
    record: &RecurrenceJson,
    basis_record: &RecurrenceJson,
    search: &CatalogLiftSearch,
) -> Result<Option<CatalogLiftConclusion>, RecurrenceError> {
    let mut claims: Vec<Claim> = search
        .identifications
        .iter()
        .filter_map(identification_claim)
        .collect();
    if let Some(stage) = &search.basis_maps {
        if let Some(lifts) = &stage.search {
            claims.extend(lifts.lifts.iter().filter_map(basis_map_claim));
            for coefficient in &stage.coefficient_identifications {
                if let Some(certification) = lifts.lifts.get(coefficient.lift) {
                    claims.extend(coefficient_claim(certification, coefficient));
                }
            }
        }
    }
    let mut conditional = None;
    for claim in claims {
        let kind = claim.parent_proof;
        let subject = if claim.on_basis_record {
            basis_record
        } else {
            record
        };
        // A flag of the record's own family (reached through the parent
        // provers' catalog-flag route on a transported coefficient family)
        // says nothing new: such a claim is circular.
        if kind == CatalogParentProofKind::CatalogAttestation
            && claim
                .conditional_on
                .iter()
                .any(|label| attests_own_family(label, subject))
        {
            continue;
        }
        if kind == CatalogParentProofKind::CatalogAttestation && conditional.is_some() {
            continue;
        }
        let Some(conclusion) = settle_nonzero_rows(subject, claim)? else {
            continue;
        };
        match kind {
            CatalogParentProofKind::Polytool => return Ok(Some(conclusion)),
            CatalogParentProofKind::CatalogAttestation => conditional = Some(conclusion),
        }
    }
    Ok(conditional)
}

/// Whether the attestation `label` (`oeis_catalog:<A-number>:lean_proved`)
/// is the flag of a catalog family with the rule of `record` up to an index
/// shift of at most 16.
fn attests_own_family(label: &str, record: &RecurrenceJson) -> bool {
    let Some(id) = label
        .strip_prefix("oeis_catalog:")
        .and_then(|rest| rest.strip_suffix(":lean_proved"))
    else {
        return false;
    };
    let Ok((rule, _, _)) = record.to_recurrence_parts() else {
        return false;
    };
    catalog_index()
        .entries
        .iter()
        .filter(|entry| entry.id == id)
        .any(|entry| (-16..=16).any(|shift| same_rule_up_to_shift(&entry.rule, &rule, shift)))
}

/// Settle the nonzero rows of a claim (module documentation, "Nonzero
/// rows"); `None` when they cannot be.
fn settle_nonzero_rows(
    record: &RecurrenceJson,
    claim: Claim,
) -> Result<Option<CatalogLiftConclusion>, RecurrenceError> {
    let (rule, first, seeds) = record.to_recurrence_parts()?;
    let (nonzero_from, certificate) = match claim.nonzero_from {
        Some(covered) => (covered, None),
        None => {
            let Ok(certificate) = certify_leading_coefficient(
                &rule,
                first,
                &seeds,
                &LeadingCoefficientOptions::default(),
            ) else {
                return Ok(None);
            };
            match certificate.holds_from.filter(|_| certificate.is_proved()) {
                Some(holds_from) => (holds_from, Some(Box::new(certificate))),
                None => return Ok(None),
            }
        }
    };
    let mut from = claim.from_index.max(first);
    let mut checked_rows = Vec::new();
    if nonzero_from > from {
        if nonzero_from - from > MAX_NONZERO_CHECKED_ROWS {
            return Ok(None);
        }
        let Ok(rows) = rule.generate_rows_rational(&seeds, first, nonzero_from - first) else {
            return Ok(None);
        };
        let start = from;
        for n in start..nonzero_from {
            if vec_poly::is_zero(&vec_poly::trim(&rows[n - first])) {
                checked_rows.clear();
                from = n + 1;
            } else {
                checked_rows.push(n);
            }
        }
    }
    let nonzero = match &certificate {
        Some(certificate) => format!(
            "rows from {} on are nonzero by the record's leading-coefficient certificate",
            certificate.holds_from.unwrap_or(first)
        ),
        None => "covered rows are nonzero (a nonzero scalar and a power of x times a nonzero \
                 parent row)"
            .to_string(),
    };
    let checked = if checked_rows.is_empty() {
        String::new()
    } else {
        format!(", rows {checked_rows:?} checked exactly")
    };
    let statement = format!(
        "every row P_n with n >= {from} is nonzero and real-rooted with every root in {}: {}; \
         {nonzero}{checked}.",
        format_root_interval_exact(&claim.region),
        claim.description
    );
    Ok(Some(CatalogLiftConclusion {
        route: claim.route,
        catalog_family: claim.catalog_family,
        basis_map: claim.basis_map,
        parent_proof: claim.parent_proof,
        conditional_on: claim.conditional_on,
        from_index: from,
        root_region: claim.region,
        leading_coefficient: certificate,
        checked_rows,
        statement,
    }))
}

// ---------------------------------------------------------------------------
// The fingerprint index
// ---------------------------------------------------------------------------

/// One indexed catalog family.
struct IndexedEntry {
    id: &'static str,
    /// The displayed record (prefix rows and initial rows as seeds).
    record: RecurrenceJson,
    /// The rule (for the own-family test).
    rule: Recurrence,
    /// Whether the entry is flagged `lean_proved`.
    flagged: bool,
}

/// The normalized-fingerprint index of the catalog (module documentation).
struct CatalogIndex {
    entries: Vec<IndexedEntry>,
    /// Fingerprint -> (entry, absolute row index).
    rows: HashMap<Vec<BigInt>, Vec<(usize, usize)>>,
}

/// A candidate family for a record.
#[derive(Debug, Clone)]
struct IndexCandidate {
    entry: usize,
    /// Source row index minus record row index of the matches.
    shift: i64,
    matched_rows: usize,
    /// The first record row (offset from its first index) that matched.
    first_matched_offset: usize,
}

/// The index, built on first use.
fn catalog_index() -> &'static CatalogIndex {
    static INDEX: OnceLock<CatalogIndex> = OnceLock::new();
    INDEX.get_or_init(CatalogIndex::build)
}

/// The normalized fingerprint of a row: the power of `x` removed,
/// denominators cleared, content divided out, leading coefficient
/// positive; `None` for rows of degree zero (and the zero row).
fn row_fingerprint(row: &[Q]) -> Option<Vec<BigInt>> {
    let row = vec_poly::trim(row);
    let start = row.iter().position(|c| !c.is_zero())?;
    let row = &row[start..];
    if row.len() < 2 {
        return None;
    }
    let integers = vec_poly::clear_denominators(row);
    let content = integers.iter().fold(BigInt::zero(), |g, c| g.gcd(c));
    let sign = if integers.last().is_some_and(|c| c.is_negative()) {
        -BigInt::one()
    } else {
        BigInt::one()
    };
    let divisor = content * sign;
    Some(integers.into_iter().map(|c| c / &divisor).collect())
}

/// The first [`INDEX_ROWS`] rows of a record.
fn index_rows(rule: &Recurrence, first: usize, seeds: &[Vec<Q>]) -> Option<Vec<Vec<Q>>> {
    rule.generate_rows_rational(seeds, first, INDEX_ROWS.max(seeds.len()))
        .ok()
}

/// Whether the rule data of `a` shifted by `shift` (`a(n + shift)`) is `b`.
fn same_rule_up_to_shift(a: &Recurrence, b: &Recurrence, shift: i64) -> bool {
    let data = |rule: &Recurrence| RecurrenceJson::from_recurrence_rational(rule, 0, &[], None);
    data(&a.shift_index(shift)).recurrence == data(b).recurrence
}

impl CatalogIndex {
    fn build() -> Self {
        let mut index = Self {
            entries: Vec::new(),
            rows: HashMap::new(),
        };
        for definition in catalog() {
            let Ok((rule, seeds, first)) = definition.recurrence_export_parts() else {
                continue;
            };
            let Ok(first) = usize::try_from(first) else {
                continue;
            };
            if rule.inhomogeneous.as_ref().is_some_and(|g| !g.is_zero()) {
                continue;
            }
            let Some(rows) = index_rows(&rule, first, &seeds) else {
                continue;
            };
            let position = index.entries.len();
            for (offset, row) in rows.iter().enumerate() {
                if let Some(key) = row_fingerprint(row) {
                    index
                        .rows
                        .entry(key)
                        .or_default()
                        .push((position, first + offset));
                }
            }
            index.entries.push(IndexedEntry {
                id: definition.id,
                record: RecurrenceJson::from_recurrence_rational(&rule, first, &seeds, None),
                rule,
                flagged: definition
                    .real_rootedness()
                    .is_some_and(|flag| flag.status == OeisRealRootednessStatus::LeanProved),
            });
        }
        index
    }

    /// The candidate families of `record`, best first (module
    /// documentation, "The fingerprint index").
    fn candidates(&self, record: &RecurrenceJson) -> Vec<IndexCandidate> {
        let Ok((rule, first, seeds)) = record.to_recurrence_parts() else {
            return Vec::new();
        };
        let Some(rows) = index_rows(&rule, first, &seeds) else {
            return Vec::new();
        };
        // (entry, shift) -> (matching rows, first matching offset)
        let mut matches: BTreeMap<(usize, i64), (usize, usize)> = BTreeMap::new();
        for (offset, row) in rows.iter().enumerate() {
            let Some(key) = row_fingerprint(row) else {
                continue;
            };
            for &(entry, source_index) in self.rows.get(&key).into_iter().flatten() {
                let shift = source_index as i64 - (first + offset) as i64;
                let slot = matches.entry((entry, shift)).or_insert((0, offset));
                slot.0 += 1;
            }
        }
        let mut best: BTreeMap<usize, IndexCandidate> = BTreeMap::new();
        for ((entry, shift), (count, first_offset)) in matches {
            if count < MINIMUM_INDEX_MATCHES {
                continue;
            }
            let indexed = &self.entries[entry];
            if indexed.record == *record || same_rule_up_to_shift(&indexed.rule, &rule, shift) {
                continue;
            }
            let candidate = IndexCandidate {
                entry,
                shift,
                matched_rows: count,
                first_matched_offset: first_offset,
            };
            let better = best.get(&entry).is_none_or(|known| {
                (count, std::cmp::Reverse(shift.unsigned_abs()))
                    > (
                        known.matched_rows,
                        std::cmp::Reverse(known.shift.unsigned_abs()),
                    )
            });
            if better {
                best.insert(entry, candidate);
            }
        }
        let mut candidates: Vec<IndexCandidate> = best.into_values().collect();
        candidates.sort_by_key(|c| (std::cmp::Reverse(c.matched_rows), c.entry));
        candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(v: i64) -> Q {
        Q::from_integer(v.into())
    }

    #[test]
    fn fingerprints_ignore_scalars_and_powers_of_x() {
        let row = vec![q(0), q(0), q(2), q(6), q(4)];
        let scaled = vec![
            Q::new(1.into(), 3.into()),
            Q::one(),
            Q::new(2.into(), 3.into()),
        ];
        let negated = vec![q(0), q(-1), q(-3), q(-2)];
        let expected = Some(vec![BigInt::from(1), BigInt::from(3), BigInt::from(2)]);
        assert_eq!(row_fingerprint(&row), expected);
        assert_eq!(row_fingerprint(&scaled), expected);
        assert_eq!(row_fingerprint(&negated), expected);
        assert_eq!(row_fingerprint(&[q(0), q(5)]), None);
        assert_eq!(row_fingerprint(&[q(0)]), None);
    }

    #[test]
    fn a_record_is_not_its_own_candidate() {
        let index = catalog_index();
        let entry = index
            .entries
            .iter()
            .find(|e| e.id == "A008292")
            .expect("indexed");
        let candidates = index.candidates(&entry.record);
        assert!(candidates
            .iter()
            .all(|c| index.entries[c.entry].id != "A008292"));
    }
}
