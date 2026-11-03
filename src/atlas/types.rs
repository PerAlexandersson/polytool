//! Record types of the atlas: nodes, edges, links, claims and certificates
//! (the JSONL schemas `polytool.atlas-node.v1`, `polytool.atlas-edge.v1`
//! and `polytool.atlas-link.v1`).

use crate::recurrence::{
    BasisMapLiftCertification, BasisMapTheoremReport, GammaLiftCertification,
    InterlacingOrientation, PreserverImageCertification, RealRootednessTheoremReport,
    RecurrenceJson, RootInterval, RowInterlacing, VeroneseParentProof,
    VeroneseSectionCertification,
};
use serde::{Deserialize, Serialize};

/// Schema tag of [`AtlasNode`].
pub const ATLAS_NODE_SCHEMA: &str = "polytool.atlas-node.v1";
/// Schema tag of [`AtlasEdge`].
pub const ATLAS_EDGE_SCHEMA: &str = "polytool.atlas-edge.v1";
/// Schema tag of [`AtlasLink`].
pub const ATLAS_LINK_SCHEMA: &str = "polytool.atlas-link.v1";

/// Number of rows stored per node (`P_(n0), ..., P_(n0 + 14)`).
pub const ATLAS_STORED_ROWS: usize = 15;

/// Lean replay state of every edge in this version.
pub const ATLAS_LEAN_NOT_ATTEMPTED: &str = "not_attempted";

/// One definition of a node's sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasDefinition {
    /// The runnable rule with its first index and identifying seeds,
    /// normalized ([`crate::atlas::canonical_definition`]).
    pub recurrence: RecurrenceJson,
    /// Where the definition comes from: `named_family:<name>`,
    /// `oeis_catalog:<A-number>`, or `derived:<operation>` (the rule
    /// computed by an operation of the atlas).
    pub source: String,
    /// The id of the [`AtlasOperation::DefinitionEquality`] edge that proves
    /// this definition equal to the node's first definition for every
    /// index; `None` for the first definition.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub equality_edge: Option<String>,
}

/// Descriptive fields of a node (the labels of a catalog page).
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasDescription {
    /// A short name (family name, OEIS title, or the operation that created
    /// the node).
    pub name: String,
    /// The polytool named family, when the node is one.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polytool_family: Option<String>,
    /// OEIS catalog ids whose bundled rule is a definition of the node and
    /// whose entry has no prefix rows before the rule's first row
    /// (identifications of the whole entry for every index, not matches).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub oeis_ids: Vec<String>,
    /// OEIS catalog ids whose bundled rule is a definition of the node but
    /// whose entry starts with prefix rows the rule does not produce: the
    /// node is the entry's tail after those rows (the prefix rows are not
    /// part of the node and may fail real-rootedness).
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub oeis_tail_ids: Vec<String>,
    /// A combinatorial interpretation, when recorded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub interpretation: Option<String>,
}

/// A proved statement about a node.  Index conventions: `from_index` is an
/// absolute index of the sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "claim", rename_all = "snake_case", deny_unknown_fields)]
pub enum AtlasClaim {
    /// Every row `P_n`, `n >= from_index`, is zero or real-rooted with every
    /// root in `region` (an index-independent interval).
    RealRooted {
        /// First index covered.
        from_index: usize,
        /// The root region.
        region: RootInterval,
    },
    /// Consecutive rows `P_(m-1)`, `P_m` interlace for every
    /// `m >= from_index`.
    Interlacing {
        /// First index `m`.
        from_index: usize,
        /// Strict or weak.
        kind: RowInterlacing,
        /// Which row carries the largest zero.
        orientation: InterlacingOrientation,
    },
    /// Every nonzero row `P_n`, `n >= from_index`, has coefficients of one
    /// sign (derived: a real-rooted polynomial with every root in
    /// `(-inf, 0]` or every root in `[0, inf)` has coefficients of one sign
    /// up to the alternation `x -> -x`; `alternating` records the second
    /// case).
    OneSignedCoefficients {
        /// First index covered.
        from_index: usize,
        /// `false`: all coefficients of a row have one sign (roots `<= 0`);
        /// `true`: the coefficients of `P_n(-x)` do (roots `>= 0`).
        alternating: bool,
    },
    /// The node's rows are Hadamard kernels (derived; Malo's theorem): for
    /// every `n >= from_index` and every real-rooted `f = sum a_k x^k`, the
    /// Hadamard product `sum a_k b_k x^k` with `P_n = sum b_k x^k` is zero
    /// or real-rooted, since `P_n` has only real roots of one sign.
    HadamardKernel {
        /// First index covered.
        from_index: usize,
    },
    /// The basis map `T_P: u^k -> P_k` sends every nonzero polynomial with
    /// roots in `input_interval` to zero or to a real-rooted polynomial with
    /// roots in `output_interval` (one conclusion of
    /// [`crate::recurrence::check_basis_map_preserves_pf`]).
    BasisMapPreserver {
        /// The source interval, as text.
        input_interval: String,
        /// The target interval, as text.
        output_interval: String,
        /// The statement in words.
        statement: String,
    },
}

impl AtlasClaim {
    /// A short label of the claim kind.
    pub fn kind(&self) -> &'static str {
        match self {
            Self::RealRooted { .. } => "real_rooted",
            Self::Interlacing { .. } => "interlacing",
            Self::OneSignedCoefficients { .. } => "one_signed_coefficients",
            Self::HadamardKernel { .. } => "hadamard_kernel",
            Self::BasisMapPreserver { .. } => "basis_map_preserver",
        }
    }
}

/// How a property of a node is backed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AtlasProvenance {
    /// The claim is a conclusion of every cited edge.
    Edges,
    /// The claim follows from the cited edges' real-rootedness claims by a
    /// stated theorem (one-signed coefficients, Hadamard kernels).
    DerivedFromEdges,
}

/// A proved property: the claim and the edges behind it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasProperty {
    /// The claim.
    pub claim: AtlasClaim,
    /// How the cited edges back it.
    pub provenance: AtlasProvenance,
    /// The ids of the edges (at least one).
    pub edges: Vec<String>,
}

/// The degree pattern of the stored rows (an observation).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "pattern", rename_all = "snake_case", deny_unknown_fields)]
pub enum DegreePattern {
    /// `deg P_n = slope n + intercept` on every nonzero stored row.
    Affine {
        /// The slope.
        slope: i64,
        /// The intercept.
        intercept: i64,
    },
    /// `deg P_n = floor((slope n + intercept) / 2)` on every nonzero stored
    /// row (half growth), and not affine.
    HalfAffine {
        /// The slope of the numerator.
        slope: i64,
        /// The intercept of the numerator.
        intercept: i64,
    },
    /// Neither (the degrees are listed in the observations).
    Irregular,
}

/// Facts observed on the stored rows only (evidence, not proofs).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasObservations {
    /// Number of stored rows observed.
    pub rows_observed: usize,
    /// `deg P_n` of each stored row (`-1` for the zero row).
    pub degrees: Vec<i64>,
    /// The degree pattern.
    pub degree_pattern: DegreePattern,
    /// Whether every nonzero stored row is palindromic (about the center of
    /// its lowest and highest nonzero coefficient).
    pub palindromic: bool,
    /// Whether every nonzero stored row is palindromic and gamma-positive.
    pub gamma_positive: bool,
    /// Whether every coefficient of every stored row is nonnegative.
    pub nonnegative_coefficients: bool,
    /// The row sums `P_n(1)` (exact integers).
    pub row_sums: Vec<String>,
}

/// A node: one integer polynomial sequence (schema
/// `polytool.atlas-node.v1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasNode {
    /// [`ATLAS_NODE_SCHEMA`].
    pub schema: String,
    /// Stable id `n-<hash>` of the first definition.
    pub id: String,
    /// Revision (0 for the first record of the id; the last revision wins).
    pub revision: u32,
    /// Index of the first stored row.
    pub first_index: usize,
    /// The stored rows `P_(first_index), ...` (at most
    /// [`ATLAS_STORED_ROWS`]), ascending exact integer coefficients, trailing
    /// zeros removed (the zero row is empty).
    pub rows: Vec<Vec<String>>,
    /// The fingerprint of the rows ([`crate::atlas::fingerprint_rows`]).
    pub fingerprint: String,
    /// The definitions (the first one defines the node).
    pub definitions: Vec<AtlasDefinition>,
    /// Descriptive fields.
    pub description: AtlasDescription,
    /// Proved properties.
    pub properties: Vec<AtlasProperty>,
    /// Observations on the stored rows.
    pub observations: AtlasObservations,
}

/// The role in which an edge uses a premise node.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AtlasPremiseRole {
    /// The node's rows are transformed.
    Sequence,
    /// The node defines the map that is applied.
    Map,
}

/// One premise of an edge.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasPremise {
    /// The premise node.
    pub node: String,
    /// Its role.
    pub role: AtlasPremiseRole,
    /// The property of the premise the step uses (a requirement: the node
    /// must have a proved claim at least as strong).
    pub uses: AtlasClaim,
}

/// The operation of an edge.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "operation", rename_all = "snake_case", deny_unknown_fields)]
pub enum AtlasOperation {
    /// A real-rootedness theorem check of the node's own rule (a root of the
    /// DAG).
    TheoremCheck {
        /// The theorem (`three_term_recurrence`, `liu_wang`, ...).
        theorem: String,
    },
    /// The shared parent provers (periodic product, paper theorem, root
    /// domain) on the node's own rule (a root of the DAG).
    ParentProvers {
        /// The route that proved it.
        route: String,
    },
    /// `Q_n = lambda(n) T_n(P_(n+s))` for a real-rootedness preserver `T_n`.
    PreserverImage {
        /// The operator, as JSON text of
        /// [`crate::recurrence::PreserverOperator`].
        operator: String,
        /// The index shift `s`.
        index_shift: i64,
    },
    /// `Q_n = T_P(R_n)` for a basis map `T_P: u^k -> P_k`.
    BasisMap {
        /// The map's name.
        map: String,
    },
    /// `Q_n = lambda(n) Gamma_(n+s)(c z)`, the gamma polynomials of a
    /// palindromic parent.
    GammaPolynomial {
        /// The scale `c`.
        scale: String,
        /// The index shift `s`.
        index_shift: i64,
    },
    /// The coefficient section `V_(r,j)(P_(n+s))`.
    VeroneseSection {
        /// The modulus `r`.
        modulus: usize,
        /// The residue `j`.
        residue: usize,
        /// The index shift `s`.
        index_shift: i64,
    },
    /// The Obreschkoff combination `Q_n = P_n + P_(n+1)` of interlacing
    /// consecutive rows.
    ObreschkoffCombination,
    /// The basis map theorem for the node's own rule (a property of the map
    /// the node defines; a root of the DAG).
    BasisMapTheorem,
    /// An all-`n` equality of a new definition with the node's first
    /// definition.
    DefinitionEquality,
}

impl AtlasOperation {
    /// A short label of the operation kind.
    pub fn kind(&self) -> &'static str {
        match self {
            Self::TheoremCheck { .. } => "theorem_check",
            Self::ParentProvers { .. } => "parent_provers",
            Self::PreserverImage { .. } => "preserver_image",
            Self::BasisMap { .. } => "basis_map",
            Self::GammaPolynomial { .. } => "gamma_polynomial",
            Self::VeroneseSection { .. } => "veronese_section",
            Self::ObreschkoffCombination => "obreschkoff_combination",
            Self::BasisMapTheorem => "basis_map_theorem",
            Self::DefinitionEquality => "definition_equality",
        }
    }
}

/// The replayable certificate of an edge.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum AtlasCertificate {
    /// A theorem report ([`RealRootednessTheoremReport::verify`]).
    TheoremCheck {
        /// The report.
        report: Box<RealRootednessTheoremReport>,
    },
    /// The parent provers on a record, replayed by an exact recomputation.
    ParentProvers {
        /// The record (rule, first index, seeds).
        record: RecurrenceJson,
        /// Whether the root-domain route was enabled.
        try_root_domain: bool,
        /// The proof.
        proof: Box<VeroneseParentProof>,
    },
    /// A preserver-image certification
    /// ([`PreserverImageCertification::verify`]).
    PreserverImage {
        /// The certification.
        certification: Box<PreserverImageCertification>,
    },
    /// A basis-map certification ([`BasisMapLiftCertification::verify`]).
    BasisMap {
        /// The certification.
        certification: Box<BasisMapLiftCertification>,
    },
    /// A gamma-lift certification ([`GammaLiftCertification::verify`]).
    GammaPolynomial {
        /// The certification.
        certification: Box<GammaLiftCertification>,
    },
    /// A Veronese-section certification, replayed by recomputing it from
    /// the recorded section rule and candidate (default options without the
    /// root-domain route) and then by
    /// [`VeroneseSectionCertification::verify`].
    VeroneseSection {
        /// The section's rule with its first index and seeds.
        input: RecurrenceJson,
        /// The certification.
        certification: Box<VeroneseSectionCertification>,
    },
    /// The Obreschkoff combination ([`crate::atlas::obreschkoff`]).
    ObreschkoffCombination {
        /// The certificate.
        certificate: Box<crate::atlas::obreschkoff::ObreschkoffCertificate>,
    },
    /// A basis-map theorem report ([`BasisMapTheoremReport::verify`]).
    BasisMapTheorem {
        /// The report.
        report: Box<BasisMapTheoremReport>,
    },
    /// An identity certification `Q_n = P_n` (a preserver image under the
    /// identity substitution `x -> x`, index shift 0, no gauge).
    DefinitionEquality {
        /// The certification.
        certification: Box<PreserverImageCertification>,
    },
}

/// Verification state of an edge.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasVerification {
    /// `replayed` once the certificate replayed when the edge was written.
    pub polytool: String,
    /// Lean replay state ([`ATLAS_LEAN_NOT_ATTEMPTED`]).
    pub lean: String,
}

/// An edge: one proof step (schema `polytool.atlas-edge.v1`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasEdge {
    /// [`ATLAS_EDGE_SCHEMA`].
    pub schema: String,
    /// Stable id `e-<hash>` of the operation, premises, conclusion and
    /// certificate.
    pub id: String,
    /// The operation.
    pub operation: AtlasOperation,
    /// The premises (empty for the roots of the DAG).
    pub premises: Vec<AtlasPremise>,
    /// The node the claims are about.
    pub conclusion: String,
    /// The claims (computed from the certificate).
    pub claims: Vec<AtlasClaim>,
    /// The cited theorem with its reference.
    pub theorem: String,
    /// The RealRooted row route, when one is recorded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lean_route: Option<String>,
    /// RealRooted declarations behind the step, when recorded.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub lean_declarations: Vec<String>,
    /// The certificate.
    pub certificate: AtlasCertificate,
    /// Verification state.
    pub verification: AtlasVerification,
}

/// The target of an OEIS match.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum OeisMatchTarget {
    /// The flattened rows (coefficients ascending).
    Triangle,
    /// The flattened reversed rows (coefficients descending).
    ReversedTriangle,
    /// The row sums `P_n(1)`.
    RowSums,
}

/// How a node's data aligns with an OEIS entry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OeisAlignment {
    /// The matched data.
    pub target: OeisMatchTarget,
    /// The absolute node index of the first row used.
    pub first_row: usize,
    /// Whether zero coefficients below the lowest nonzero one were dropped
    /// from each row (triangles only).
    pub low_zeros_dropped: bool,
    /// Position (0-based) in the OEIS term list of the first matched term.
    pub oeis_term_position: usize,
    /// The OEIS offset of the entry.
    pub oeis_offset: i64,
    /// Number of terms compared (all agree).
    pub terms_compared: usize,
}

/// The kind of a link.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(tag = "link", rename_all = "snake_case", deny_unknown_fields)]
pub enum AtlasLinkKind {
    /// Two nodes whose definitions agree on the stored rows (same first
    /// index), without an all-`n` equality certificate: a conjectured
    /// identity.
    CandidateEqual {
        /// Rows compared.
        rows_compared: usize,
        /// Why no certificate was found.
        reason: String,
    },
    /// Two nodes with the same fingerprint whose rows differ: equal up to
    /// content, sign, reversal and a small index shift.
    FingerprintRelated {
        /// The fingerprint.
        fingerprint: String,
    },
    /// An OEIS entry matched by a node's triangle, reversed triangle or row
    /// sums (finite agreement; a row-sum match is never an identification).
    OeisMatch {
        /// The A-number.
        oeis_id: String,
        /// Whether the entry is in polytool's bundled OEIS catalog.
        in_catalog: bool,
        /// The alignment.
        alignment: OeisAlignment,
    },
}

impl AtlasLinkKind {
    /// A short label of the link kind.
    pub fn kind(&self) -> &'static str {
        match self {
            Self::CandidateEqual { .. } => "candidate_equal",
            Self::FingerprintRelated { .. } => "fingerprint_related",
            Self::OeisMatch { .. } => "oeis_match",
        }
    }
}

/// A link: evidence that is not a proof (schema `polytool.atlas-link.v1`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasLink {
    /// [`ATLAS_LINK_SCHEMA`].
    pub schema: String,
    /// Stable id `l-<hash>` of the kind and the nodes.
    pub id: String,
    /// The nodes (one for an OEIS match, two otherwise).
    pub nodes: Vec<String>,
    /// The kind and its evidence.
    pub kind: AtlasLinkKind,
}

/// A sequence produced by a route, before it is matched to a node.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasSequenceSpec {
    /// The definition (canonical).
    pub definition: RecurrenceJson,
    /// The definition's source ([`AtlasDefinition::source`]).
    pub source: String,
    /// A name for a new node.
    pub name: String,
    /// The polytool family, when the sequence is a named family.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub polytool_family: Option<String>,
    /// The OEIS catalog id, when the definition is a catalog rule.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub oeis_id: Option<String>,
    /// Number of prefix rows of that catalog entry before the rule's first
    /// row (`0`: the definition is the whole entry).
    #[serde(default, skip_serializing_if = "is_zero")]
    pub oeis_prefix_rows: usize,
}

fn is_zero(value: &usize) -> bool {
    *value == 0
}

/// A proof step produced by a route: the certificate and the sequences it
/// concerns, before the premises and the conclusion are matched to nodes
/// ([`crate::atlas::integrate_derivations`]).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtlasDerivation {
    /// The conclusion sequence.
    pub conclusion: AtlasSequenceSpec,
    /// The operation.
    pub operation: AtlasOperation,
    /// The certificate.
    pub certificate: AtlasCertificate,
}
