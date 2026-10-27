//! Request, limits and report types of the bounded representation search
//! ([`search_representations`]): the `polytool.recurrence-representation-search.v1`
//! schema, its labels and ranking data.
//!
//! Related modules: `representation_search` (the queue and replay),
//! `representation_map` (certified map normal forms).

use super::representation_map::RepresentationMapRecord;
use super::*;
use serde::{Deserialize, Serialize};

/// Schema tag of [`RepresentationSearchReport`].
pub const RECURRENCE_REPRESENTATION_SEARCH_SCHEMA: &str =
    "polytool.recurrence-representation-search.v1";

/// Label: fitted to finitely many rows, or resting on a condition checked on rows only.
pub const FINITE_FIT_LABEL: &str = "finite_fit";

/// Label: annihilation is proved for every index; uniqueness is not certified.
pub const CERTIFIED_RELATION_LABEL: &str = "certified_annihilating_relation";

/// Label: certified annihilation, sufficient initial/boundary data and a
/// certified uniqueness argument.
pub const EQUIVALENT_DESCRIPTION_LABEL: &str = "equivalent_generating_description";

/// Number of source rows generated when `row_count` is omitted.
pub const DEFAULT_REPRESENTATION_ROW_COUNT: usize = 12;

/// Largest number of transformations in one request.
pub const MAX_REPRESENTATION_TRANSFORMATIONS: usize = 32;

/// Largest number of unknowns of the first-order vector fit of a vector
/// representation (matrix entries and left factor at the search bounds); a
/// larger fit is skipped with the budget event `vector_fit_unknowns`.
pub const MAX_VECTOR_FIT_UNKNOWNS: usize = 512;

/// The search label of a description with the given assessment status
/// (`equivalent_description`, `certified_consequence`, ...).
pub fn representation_label(assessment_status: &str) -> &'static str {
    match assessment_status {
        EQUIVALENT_DESCRIPTION => EQUIVALENT_DESCRIPTION_LABEL,
        CERTIFIED_CONSEQUENCE => CERTIFIED_RELATION_LABEL,
        _ => FINITE_FIT_LABEL,
    }
}

pub(super) fn label_rank(label: &str) -> usize {
    match label {
        EQUIVALENT_DESCRIPTION_LABEL => 0,
        CERTIFIED_RELATION_LABEL => 1,
        _ => 2,
    }
}

fn default_true() -> bool {
    true
}

/// The default bounds of [`RepresentationSearchRequest::search`]: smaller
/// than [`TransformSearchBounds::default`] because they are used in every
/// representation.
pub fn default_representation_search_bounds() -> TransformSearchBounds {
    TransformSearchBounds {
        max_rec_len: 2,
        max_var_deg: 2,
        max_idx_deg: 1,
        max_diff_deg: 1,
        try_denominator: true,
        max_candidates: 2_000,
        ..TransformSearchBounds::default()
    }
}

fn is_true(value: &bool) -> bool {
    *value
}

// ---------------------------------------------------------------------------
// Request
// ---------------------------------------------------------------------------

/// Limits of [`search_representations`]; reaching one is reported, never
/// treated as absence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default, deny_unknown_fields)]
pub struct RepresentationSearchLimits {
    /// Largest number of maps composed from the source (queue depth).
    pub max_depth: usize,
    /// Largest number of representations, the source included.
    pub max_representations: usize,
    /// Largest number of backend calls (transform, implicit, space shape, fit menu).
    pub max_steps: usize,
    /// Largest number of certified parent descriptions transported along one edge.
    pub max_sources_per_edge: usize,
    /// Largest number of complete-space supports tried for one representation.
    pub max_space_shapes: usize,
    /// Largest number of rows.
    pub max_rows: usize,
    /// Largest row degree or degree bound.
    pub max_degree: usize,
    /// Largest candidate budget of one fit.
    pub max_search_candidates: usize,
    /// Largest degree or order bound of a fit.
    pub max_search_degree: usize,
    /// Elementary-operation budget of each complete-space elimination phase.
    pub max_operations: usize,
    /// Largest bit length of a rational numerator or denominator in a complete space.
    pub max_coefficient_bits: u64,
}

impl Default for RepresentationSearchLimits {
    fn default() -> Self {
        Self {
            max_depth: 2,
            max_representations: 16,
            max_steps: 160,
            max_sources_per_edge: 2,
            max_space_shapes: 16,
            max_rows: 64,
            max_degree: 256,
            max_search_candidates: 20_000,
            max_search_degree: 8,
            max_operations: 20_000_000,
            max_coefficient_bits: 4096,
        }
    }
}

impl RepresentationSearchLimits {
    /// The largest limits accepted by service adapters such as the MCP server.
    pub fn service_ceiling() -> Self {
        Self {
            max_depth: 4,
            max_representations: 64,
            max_steps: 640,
            max_sources_per_edge: 4,
            max_space_shapes: 32,
            max_rows: 128,
            max_degree: 512,
            max_search_candidates: 50_000,
            max_search_degree: 8,
            max_operations: 50_000_000,
            max_coefficient_bits: 8192,
        }
    }

    /// Reject limits above a service ceiling.
    pub fn check_ceiling(&self, ceiling: &Self) -> Result<(), String> {
        let pairs = [
            ("max_depth", self.max_depth, ceiling.max_depth),
            (
                "max_representations",
                self.max_representations,
                ceiling.max_representations,
            ),
            ("max_steps", self.max_steps, ceiling.max_steps),
            (
                "max_sources_per_edge",
                self.max_sources_per_edge,
                ceiling.max_sources_per_edge,
            ),
            (
                "max_space_shapes",
                self.max_space_shapes,
                ceiling.max_space_shapes,
            ),
            ("max_rows", self.max_rows, ceiling.max_rows),
            ("max_degree", self.max_degree, ceiling.max_degree),
            (
                "max_search_candidates",
                self.max_search_candidates,
                ceiling.max_search_candidates,
            ),
            (
                "max_search_degree",
                self.max_search_degree,
                ceiling.max_search_degree,
            ),
            (
                "max_operations",
                self.max_operations,
                ceiling.max_operations,
            ),
        ];
        for (name, value, limit) in pairs {
            if value > limit {
                return Err(format!(
                    "{name} = {value} exceeds the service ceiling {limit}"
                ));
            }
        }
        if self.max_coefficient_bits > ceiling.max_coefficient_bits {
            return Err(format!(
                "max_coefficient_bits = {} exceeds the service ceiling {}",
                self.max_coefficient_bits, ceiling.max_coefficient_bits
            ));
        }
        Ok(())
    }

    pub(super) fn transform_limits(&self) -> TransformLimits {
        TransformLimits {
            max_rows: self.max_rows,
            max_degree: self.max_degree,
            max_rule_derivative_order: TransformLimits::default().max_rule_derivative_order,
            max_search_candidates: self.max_search_candidates,
            max_search_degree: self.max_search_degree,
        }
    }

    pub(super) fn space_limits(&self) -> RecurrenceSpaceLimits {
        RecurrenceSpaceLimits {
            max_operations: self.max_operations,
            max_coefficient_bits: self.max_coefficient_bits,
            ..RecurrenceSpaceLimits::default()
        }
    }
}

/// Request for [`search_representations`].
///
/// Give a defining rule (`source_recurrence`), a bundled OEIS identifier
/// (`oeis_id`) or rows only (`rows`, exploratory fits). With a source,
/// optional `rows` must agree with the generated rows.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationSearchRequest {
    /// Defining rule with seeds; its `first_index` is the index origin.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_recurrence: Option<RecurrenceJson>,
    /// Bundled OEIS family (displayed row indices), used as the defining rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oeis_id: Option<String>,
    /// Source rows `P_(row_origin + i)`, ascending rational strings.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rows: Option<Vec<Vec<String>>>,
    /// Original index of `rows[0]`; defaults to the source's first index, else 0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_origin: Option<usize>,
    /// Rows generated from a source when `rows` is absent (default 12).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub row_count: Option<usize>,
    /// The allowed maps, applied at every expanded representation in this order.
    #[serde(default)]
    pub transformations: Vec<RecurrenceTransformation>,
    /// Bounds of the fit menu, of the successor fit of implicit steps and of
    /// the default complete-space menu. Omitted means
    /// [`default_representation_search_bounds`] (lags up to 2, `x`-degree 2,
    /// index degree 1, derivative order 1, left factors allowed, 2000
    /// candidates).
    #[serde(default = "default_representation_search_bounds")]
    pub search: TransformSearchBounds,
    /// Run the fit menu in every scalar representation (default `true`).
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub fits: bool,
    /// Certify complete spaces after factorial diagonal maps, Möbius-class
    /// maps, polynomial diagonal weights and differentiation of a certified
    /// constant-index monic rule (default `true`).
    #[serde(default = "default_true", skip_serializing_if = "is_true")]
    pub certified_spaces: bool,
    /// Explicit complete-space supports; empty means the default menu
    /// (lags `1..=max_rec_len`, derivative orders `0..=max_diff_deg`, uniform
    /// coefficient degrees increasing to the search bounds, stopping at the
    /// first successful degree and skipping supports dominated by a success).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub space_shapes: Vec<RecurrenceSpaceShape>,
    /// Optional root-domain analysis attached to every coordinate-change
    /// step ([`RootDomainRequest`]); omitted means none, and discovery never
    /// depends on it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub root_domain: Option<RootDomainRequest>,
    /// Limits.
    #[serde(default)]
    pub limits: RepresentationSearchLimits,
}

impl RepresentationSearchRequest {
    /// A request for the bundled OEIS family `id` with the given maps.
    pub fn from_oeis(id: &str, transformations: Vec<RecurrenceTransformation>) -> Self {
        Self {
            source_recurrence: None,
            oeis_id: Some(id.to_string()),
            rows: None,
            row_origin: None,
            row_count: None,
            transformations,
            search: default_representation_search_bounds(),
            fits: true,
            certified_spaces: true,
            space_shapes: Vec::new(),
            root_domain: None,
            limits: RepresentationSearchLimits::default(),
        }
    }
}

// ---------------------------------------------------------------------------
// Report
// ---------------------------------------------------------------------------

/// The resolved source.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationSource {
    /// `defining_rule`, `oeis` or `rows_only`.
    pub kind: String,
    /// The defining rule with seeds (original indices), when given.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rule: Option<RecurrenceJson>,
    /// Original index of the first row.
    pub row_origin: usize,
    /// Number of source rows (stored in representation 0).
    pub row_count: usize,
}

/// Ranking data of a description; smaller is simpler.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct DescriptionCost {
    /// Number of rows spanned minus one (largest lag or shift).
    pub shift_length: usize,
    /// Largest derivative order on any row.
    pub derivative_order: usize,
    /// Largest power of `x` in a coefficient.
    pub coefficient_degree: usize,
    /// Largest power of the index `n` in a coefficient (0 means index independent).
    pub index_degree: usize,
    /// Number of components of the state (the number of coupled relations
    /// for the two-relation gamma system).
    pub component_count: usize,
    /// Maps to invert back to the source plus coefficient conditions.
    pub reconstruction_cost: usize,
}

/// One description of a representation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationDescription {
    /// `finite_fit`, `certified_annihilating_relation` or
    /// `equivalent_generating_description`.
    pub label: String,
    /// `defining_rule`, `transported_rule`, `transported_vector_system`,
    /// `transported_gamma_system` (the even part and odd companion of a
    /// rank-two gamma transport), `implicit_transport`,
    /// `reconstruction_certified_fit`, `certified_space`, `finite_fit` or
    /// `vector_fit` (a first-order vector system fitted to a vector
    /// representation).
    pub kind: String,
    /// Evidence of annihilation (`defining`, `transported_identity`,
    /// `certified_by_reconstruction`, `certified_complete_space`,
    /// `checked_supplied_rows`, `finite_fit`).
    pub evidence: String,
    /// Human-readable form.
    pub text: String,
    /// The successor rule with seeds (`polytool.recurrence.v1`, original indices).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub scalar: Option<RecurrenceJson>,
    /// A vector system with its initial state.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector: Option<TransformVectorRecurrenceJson>,
    /// Seeds, relations and conditions as a generating description.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<GeneratingDescription>,
    /// Assessment of `description`: annihilation evidence and uniqueness.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assessment: Option<DescriptionAssessment>,
    /// Ranking data.
    pub cost: DescriptionCost,
    /// Where it comes from: `source`, `steps[i]`, `spaces[i]...` or `fits`.
    pub provenance: String,
    /// Further provenances of the same relation (merged duplicates).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub also_found_by: Vec<String>,
    /// Further explanation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    /// Uniqueness of `vector` with its initial state, for vector systems.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector_assessment: Option<VectorSystemAssessment>,
}

/// Whether a first-order vector system `q(n,x) F_(n+1) = M F_n (+ G)` with
/// its initial state determines the states: it does when `q(n, x)` is a
/// nonzero polynomial in `x` at every index of the domain, certified by a
/// coefficient `[x^j] q` with no integer root in the domain.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct VectorSystemAssessment {
    /// `certified` or `not_certified`.
    pub status: String,
    /// The power `j` of `x` whose coefficient of `q` certifies `q(n, .) != 0`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub denominator_coefficient: Option<usize>,
    /// The argument or why it fails.
    pub statement: String,
}

/// A certified vector system of the parent transported through a
/// componentwise Möbius-class map ([`transport_vector_componentwise`]).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct VectorTransportStep {
    /// The map applied to each component, with the schedule shifted by the
    /// component's degree offset.
    pub component_maps: Vec<RecurrenceTransformation>,
    /// The transported system with the mapped initial state.
    pub transported: TransformVectorRecurrenceJson,
    /// Human-readable form.
    pub text: String,
    /// `transported_identity` when the parent system is certified, every
    /// map condition is certified and the system holds on the rows; else
    /// `checked_supplied_rows`.
    pub evidence: String,
    /// Whether every transition of the new rows from the system's first
    /// index satisfies the transported system.
    pub holds_on_rows: bool,
    /// Number of transitions checked.
    pub checked_transitions: usize,
}

/// One replayable backend call of an edge.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationStep {
    /// Index into the parent's `descriptions` of the certified description
    /// used as the step's defining rule; `None` maps rows only.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_description: Option<usize>,
    /// A `polytool.recurrence-transform.v1` report.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transform: Option<RecurrenceTransformReport>,
    /// A `polytool.recurrence-implicit.v1` report.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub implicit: Option<ImplicitTransformReport>,
    /// A certified relation of the parent (a scalar rule in implicit form or
    /// an implicit-only relation) transported by the implicit machinery
    /// (diagonal maps, derivative antecedents and the lossy `derivative`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub implicit_relation: Option<ImplicitRelationStep>,
    /// A certified vector system of a vector parent transported through a
    /// componentwise Möbius-class map.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub vector_transport: Option<VectorTransportStep>,
    /// Why the step failed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
}

/// A parent relation transported into the new coordinates by
/// [`RecurrenceTransformation::transport_implicit`], with the parent's
/// coefficient conditions carried along where the map allows it.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ImplicitRelationStep {
    /// The parent's relation (parent coordinates) with its domain and evidence.
    pub source_relation: DescribedRelation,
    /// The transported relation (`transported_identity` when the source and
    /// the map hypotheses are certified, else `checked_supplied_rows`).
    pub relation: DescribedRelation,
    /// Clearing data of a diagonal map.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clearing: Option<ClearingCertificate>,
    /// The recorded operator identity of the transport (derivative
    /// antecedents, the lossy derivative map, polynomial inverse weights).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub certificate: Option<ReconstructionCertificate>,
    /// Parent conditions that could not be carried into the new coordinates
    /// (they refer to lost coefficients), as `list[i]` references.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub dropped_conditions: Vec<String>,
    /// Seeds, the transported relation, the exceptional-index conditions of
    /// the clearing and the carried parent conditions.
    pub description: GeneratingDescription,
    /// Assessment of `description`.
    pub assessment: DescriptionAssessment,
    /// For the lossy `derivative` map: the description of the parent through
    /// the new coordinates, with and without the lost data.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lossy: Option<LossyReconstruction>,
    /// Whether the relation and conditions hold at every checked base of the
    /// new rows.
    pub holds_on_rows: bool,
    /// Number of base indices checked.
    pub checked_bases: usize,
}

/// The data a lossy map of the chain forgets, recorded so that the source
/// rows can be rebuilt exactly.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct LostDataRecord {
    /// Position in `chain` of the lossy map.
    pub chain_position: usize,
    /// The derivative order `m`: the coefficients `[x^0..x^(m-1)]` of the
    /// rows entering the map are lost.
    pub order: usize,
    /// `rows[i]`: the lost coefficients of the row with index
    /// `row_origin + i` entering the map, ascending rational strings.
    pub rows: Vec<Vec<String>>,
}

/// Outcome of the fit menu in one representation.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationFitSummary {
    /// `complete`, `budget_exhausted`, `insufficient_rows`, `skipped_budget` or `failed`.
    pub status: String,
    /// Configurations inspected.
    pub configurations_considered: usize,
    /// Fits returned by the menu (before merging with certified relations).
    pub fits_found: usize,
    /// Explanation for `failed`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
}

/// Certified equality of two representations' sequences.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationEquality {
    /// The earlier representation.
    pub representation: usize,
    /// Its equivalent description.
    pub other_description: usize,
    /// This representation's equivalent description.
    pub description: usize,
    /// Why the sequences are equal.
    pub certificate: String,
}

/// One representation: the source or the image of a chain of maps.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationNode {
    /// Index in `representations`.
    pub id: usize,
    /// Number of maps from the source.
    pub depth: usize,
    /// Parent representation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub parent: Option<usize>,
    /// Index into `request.transformations` of the edge map.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transformation: Option<usize>,
    /// Indices into `request.transformations`, in application order from the source.
    pub chain: Vec<usize>,
    /// Normal form of the composite map.
    pub map: RepresentationMapRecord,
    /// Component names (`P` for the source, `Q`, or `a`, `b`; maps after a
    /// vector map act componentwise and keep the names).
    pub components: Vec<String>,
    /// Original index of the first row.
    pub row_origin: usize,
    /// Rows: `[row][component][coefficient]`, ascending rational strings.
    pub states: Vec<Vec<Vec<String>>>,
    /// Whether inverting the whole chain reproduces the source rows exactly
    /// (with the recorded `lost_data` for lossy maps).
    pub reconstruction_verified: bool,
    /// The data forgotten by every lossy map of the chain; empty for chains
    /// of bijective maps.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lost_data: Vec<LostDataRecord>,
    /// Preconditions and image conditions of the edge map (from its first step).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub map_conditions: Vec<TransformConditionStatus>,
    /// Descriptions, ranked.
    pub descriptions: Vec<RepresentationDescription>,
    /// Backend calls of the edge.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub steps: Vec<RepresentationStep>,
    /// Certified complete-space reports, one per tried support.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub spaces: Vec<RecurrenceSpaces>,
    /// Complete-space searches after Möbius-class maps, polynomial diagonal
    /// weights and differentiation (`polytool.recurrence-implicit-search.v1`
    /// reports with one successor support each), from the parent's certified
    /// constant-index monic rule; see [`search_implicit_relations`] for the
    /// finite verification theorems.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub implicit_spaces: Vec<ImplicitRelationSearchReport>,
    /// The fit menu.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fits: Option<RepresentationFitSummary>,
    /// Certified equality with an earlier representation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equal_sequence_as: Option<RepresentationEquality>,
    /// Whether the allowed maps were applied to this representation.
    pub expanded: bool,
    /// Further explanations.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

/// One attempted application of a map.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationEdge {
    /// The representation the map was applied to.
    pub parent: usize,
    /// Index into `request.transformations`.
    pub transformation: usize,
    /// `new_representation`, `duplicate_map`, `rejected` or `skipped_budget`.
    pub outcome: String,
    /// The new or duplicated representation.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representation: Option<usize>,
    /// Certificate of a duplicate, or why the map was rejected or skipped.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

/// A limit that cut work.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationBudgetEvent {
    /// The limit (`max_representations`, `max_steps`, `max_space_shapes`,
    /// `fit_candidates`, `space_operations`, `implicit_fit_candidates`,
    /// `vector_fit_unknowns`).
    pub limit: String,
    /// The representation concerned.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub representation: Option<usize>,
    /// What was not done.
    pub detail: String,
}

/// A description in the global ranking.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RankedDescription {
    /// Representation index.
    pub representation: usize,
    /// Index into that representation's `descriptions`.
    pub description: usize,
    /// Its label.
    pub label: String,
    /// Its text.
    pub text: String,
    /// Its ranking data.
    pub cost: DescriptionCost,
}

/// Versioned report of [`search_representations`]
/// (`polytool.recurrence-representation-search.v1`).
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct RepresentationSearchReport {
    /// `polytool.recurrence-representation-search.v1`.
    pub schema: String,
    /// The request, retained for replay.
    pub request: RepresentationSearchRequest,
    /// `complete` or `budget_exhausted`.
    pub status: String,
    /// The resolved source (the reference for replay).
    pub source: RepresentationSource,
    /// Representations; 0 is the source.
    pub representations: Vec<RepresentationNode>,
    /// Every attempted map application, in queue order.
    pub edges: Vec<RepresentationEdge>,
    /// Every description of every representation, simplest first.
    pub ranked: Vec<RankedDescription>,
    /// Representations at `max_depth` that were not expanded.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub depth_frontier: Vec<usize>,
    /// Limits that cut work.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub budget: Vec<RepresentationBudgetEvent>,
    /// Backend calls made.
    pub steps_used: usize,
    /// Conventions of the report.
    pub conventions: Vec<String>,
}

pub(super) fn conventions() -> Vec<String> {
    [
        "rows are ascending coefficient vectors indexed by the original index n; P is the \
         source and Q (or a, b) the image of a chain of maps",
        "labels: finite_fit (finite evidence), certified_annihilating_relation (holds for \
         every index of the domain), equivalent_generating_description (certified relation \
         plus initial/boundary data and a certified uniqueness argument)",
        "certification is relative to the defining source, the maps actually enumerated and \
         the finite supports searched; finite fits are never promoted",
        "duplicate maps are detected by certified equality of normal forms, never by rows",
        "ranking: shift length, derivative order, coefficient degree, index degree, component \
         count, reconstruction cost, then label",
        "lossy maps (derivative) record the forgotten coefficients in lost_data; the parent is \
         described through them only together with conditions pinning the lost data",
        "vector representations (a, b) are mapped componentwise by scalar bijective maps, a \
         schedule shifted by the component's degree offset (b at d_n - 1); their systems are \
         transported through Möbius-class maps and fitted as first-order vector systems",
        "a rank-two gamma transport is the two-relation system of its even part and odd \
         companion, certified only where the symmetry of the rows is certified",
    ]
    .into_iter()
    .map(String::from)
    .collect()
}
