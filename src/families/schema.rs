//! Strict, versioned family sidecar schema. Exact coefficients are strings.
use crate::exponential_rows::ExponentialFamily;
use crate::recurrence::RecurrenceJson;
use serde::{Deserialize, Serialize};

/// A canonical generator; no coefficients are duplicated for OEIS entries.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum FamilyGenerator {
    /// Reference to the bundled OEIS generator.
    Oeis {
        /// Existing uppercase A-number.
        id: String,
    },
    /// Exact labelled-block exponential family, with weights starting at size one.
    ExpRows {
        /// The series A and B in n![z^n] exp(A+tB).
        family: ExponentialFamily,
    },
    /// A complete polytool.recurrence.v1 definition.
    Recurrence {
        /// Recurrence with exact initial rows and index convention.
        record: Box<RecurrenceJson>,
    },
}

/// Proof vocabulary extends the OEIS flags with explicit unknown and conjectural states.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum FamilyProofStatus {
    /// No assertion has been curated.
    #[default]
    Unknown,
    /// An open conjecture, possibly with finite evidence.
    Conjectured,
    /// Only a finite range is recorded; no uniform claim.
    FiniteEvidence,
    /// An external Lean proof and all-index identification are attested.
    LeanProved,
    /// A cited theorem and all-index identification are attested.
    TheoremProved,
    /// An exact nonreal-root counterexample is known.
    KnownNotRealRooted,
}

/// Public commit-pinned Lean source. Manual attestations until the proof owner exports tags.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyLeanProof {
    /// Public repository URL, currently the public RealRooted repository only.
    pub repository: String,
    /// Fully qualified declaration name.
    pub declaration: String,
    /// Full 40-character lowercase commit hash.
    pub commit: String,
    /// Source permalink under repository/blob/commit/.
    pub permalink: String,
}

/// An inclusive range of displayed row indices; evidence is not a proof.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyFiniteRange {
    /// First displayed row (or left index of a consecutive pair).
    pub first: i64,
    /// Last displayed row (or left index of a consecutive pair).
    pub last: i64,
    /// Method and provenance; reported evidence need not have been replayed locally.
    pub evidence: String,
}

/// A proof attestation or explicitly finite evidence for one property.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyProperty {
    /// Status, separate from the generator's validation status.
    pub status: FamilyProofStatus,
    /// Snapshot date, YYYY-MM-DD.
    pub as_of: Option<String>,
    /// Reference key in this family's references, when a theorem is cited.
    pub theorem: Option<String>,
    /// All-index identification, hypotheses, or other qualification.
    pub statement: String,
    /// Optional manual public Lean attestation, never a private repository link.
    pub lean: Option<FamilyLeanProof>,
    /// Optional paper proof description.
    pub paper_proof: Option<String>,
    /// Optional finite evidence with displayed index range.
    pub finite_range: Option<FamilyFiniteRange>,
    /// Positive lag for interlacing; None for a row property.
    pub lag: Option<usize>,
}

/// A bibliographic item. Locators are preserved from the source audit.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyReference {
    /// Unique local key used by theorem fields.
    pub key: String,
    /// Citation text.
    pub citation: String,
    /// DOI when available.
    pub doi: Option<String>,
    /// arXiv identifier, including version if relevant.
    pub arxiv: Option<String>,
    /// Theorem, equation, or source section.
    pub locator: String,
    /// Public source URL.
    pub url: String,
}

/// Optional model identification; this catalog does not perform recognition.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyModel {
    /// Model family name.
    pub family: String,
    /// Model parameters as a JSON object.
    pub parameters: serde_json::Map<String, serde_json::Value>,
    /// pending, passed, or failed, describing an independent replay.
    pub replay_status: String,
}

/// A directed relation to another family in this catalog.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RelatedFamily {
    /// Existing family ID.
    pub id: String,
    /// Precise relation, e.g. shift, reversal, or a named transform with parameters.
    pub relation: String,
}

/// OEIS identity, or a dated record of searches with no match.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyOeis {
    /// A-number, or None (JSON null), never the literal string "none".
    pub id: Option<String>,
    /// Last search date when available.
    pub checked: Option<String>,
    /// Exact queries used, including whether these were flattened terms or row sums.
    pub queries: Vec<String>,
}

/// Submission bookkeeping, independent of mathematical proof status.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilySubmission {
    /// not_started, draft, ready, submitted, accepted, or withdrawn.
    pub state: String,
    /// Submission date (required once submitted or accepted).
    pub date: Option<String>,
    /// Assigned A-number, required when accepted.
    pub a_number: Option<String>,
}

/// One version-1 family record; empty curated text is deliberately left empty.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyRecord {
    /// Schema version; currently 1.
    pub schema_version: u32,
    /// A###### or PT-F####; PT IDs remain stable after OEIS assignment.
    pub id: String,
    /// Human-readable title.
    pub title: String,
    /// Combinatorial meaning and conventions.
    pub definition: String,
    /// Formula text (TeX permitted).
    pub formula: String,
    /// Generating-function text (TeX permitted); executable series live in generator.
    pub generating_function: String,
    /// Canonical generator.
    pub generator: FamilyGenerator,
    /// Displayed index of the first generated row.
    pub first_row: i64,
    /// Index of the first coefficient in the flattened sequence.
    pub flattened_offset: i64,
    /// OEIS layout keyword: tabl or tabf.
    pub layout: String,
    /// Description of row degrees, empty when not curated.
    pub degree_pattern: String,
    /// First complete rows in ascending powers of t, exact integer strings.
    pub initial_rows: Vec<Vec<String>>,
    /// Bibliographic sources.
    pub references: Vec<FamilyReference>,
    /// Real-rootedness attestation, separate from finite generator checks.
    pub real_rootedness: FamilyProperty,
    /// Interlacing attestation; lag states which pairs are meant.
    pub interlacing: FamilyProperty,
    /// Independent model identifications.
    pub models: Vec<FamilyModel>,
    /// Relations to other catalog IDs.
    pub related: Vec<RelatedFamily>,
    /// OEIS identity or novelty search log.
    pub oeis: FamilyOeis,
    /// Submission status.
    pub submission: FamilySubmission,
    /// One-way public external links, including submission proposals.
    pub links: Vec<String>,
    /// Generator validation status inherited from the OEIS catalog; None otherwise.
    pub generator_status: Option<String>,
}
