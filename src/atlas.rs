//! The real-rooted atlas: a proof DAG of integer polynomial sequences, each
//! proved real-rooted (with a root region) by polytool's automatic routes.
//!
//! # Data model
//!
//! - A **node** ([`AtlasNode`]) is a sequence `P_n`, `n >= n0`, of
//!   polynomials in `Z[x]` (not a recurrence).  It stores its first
//!   [`ATLAS_STORED_ROWS`] rows (ascending coefficients, exact integers), its
//!   **definitions** (runnable rules with seeds; a second definition is
//!   admitted only with an all-`n` equality certificate, an
//!   [`AtlasOperation::DefinitionEquality`] edge), a description (name,
//!   polytool family, OEIS catalog ids), its proved **properties**
//!   ([`AtlasProperty`]: a claim and the edges that prove it) and
//!   **observations** made on the stored rows only ([`AtlasObservations`]:
//!   degree pattern, palindromicity, gamma-positivity, coefficient signs).
//!   Its id is a stable hash of its first (normalized) definition.
//! - An **edge** ([`AtlasEdge`]) is a proof step: premises (nodes in the
//!   role of a sequence or of a map, each with the property it uses) give
//!   claims about one node.  It records the operation, the replayable
//!   polytool certificate ([`AtlasCertificate`]), the cited theorem, the
//!   RealRooted route and declarations where known, and the verification
//!   state.  A node may have several incoming edges; all are kept.
//! - A **link** ([`AtlasLink`]) is evidence that is not a proof: two
//!   definitions that agree on the stored rows only (*candidate-equal*),
//!   sequences with the same fingerprint (equal up to content, sign,
//!   reversal or an index shift), and OEIS matches of the triangle, the
//!   reversed triangle or the row sums.
//!
//! # Soundness
//!
//! Every claim of an edge is computed from its certificate by
//! [`analyze_certificate`], which first replays the certificate with the
//! report type's own `verify()` (or, for the parent provers, by an exact
//! recomputation), and then reads the conclusion off the verified record:
//!
//! - theorem checks ([`crate::recurrence::RealRootednessTheoremReport`]) and
//!   the parent provers (periodic product, paper theorem, root domain) have
//!   no premises: they are the roots of the DAG;
//! - a lift (preserver image, basis map, gamma polynomial, Veronese
//!   section) is accepted only when its identity is certified for every
//!   index and every hypothesis of its transfer theorem is verified; its
//!   parent statement is **not** taken from the lift's own parent proof but
//!   is a premise: the parent must be an atlas node whose proved root region
//!   lies inside the region the transfer theorem requires
//!   ([`crate::recurrence::LiftRootTransfer::parent_root_region`]), from an
//!   index at most the first parent index the lift uses;
//! - an Obreschkoff combination `P_n + P_(n+1)` rests on the interlacing of
//!   consecutive parent rows (a proved property of the parent node) and the
//!   Hermite--Kakeya--Obreschkoff theorem; its identity with the image rule
//!   is certified by the preserver machinery ([`obreschkoff`]).
//!
//! [`audit_atlas`] replays every edge in file order from the roots,
//! re-derives every claim, checks every premise against the claims already
//! established, checks that every node property is backed by its cited
//! edges, regenerates every node's rows from its definition and checks every
//! real-rootedness claim exactly on the stored rows.
//!
//! # Files
//!
//! A directory with `nodes.jsonl`, `edges.jsonl` and `links.jsonl`
//! ([`AtlasStore`]): one JSON record per line, append-only.  A node is
//! updated by appending a new revision with the same id (the last revision
//! wins); nothing is deleted.  Records are appended in a deterministic
//! order.  Search indexes (fingerprints, row sums, OEIS terms) are built
//! from the files and never stored with them.
//!
//! # Growth
//!
//! - [`seed`]: every named family of polytool and every OEIS catalog entry
//!   is tried; a node is added only when a theorem check, a parent prover or
//!   a lift from an existing node proves it.  Catalog entries that no route
//!   proves stay outside the atlas (they are benchmarks and match targets).
//! - [`grow`]: one forward round of budgeted operations (preserver images,
//!   basis maps, gamma maps, Veronese sections, Obreschkoff combinations) on
//!   existing nodes, with exact image rules, integer rows and
//!   de-duplication by fingerprint.
//! - [`oeis_match`]: matching of every node's triangle, reversed triangle
//!   and row sums against a local OEIS snapshot (an index built once and
//!   cached outside the data directory).
//!
//! Related: [`crate::recurrence`] (every certificate type), [`crate::oeis`]
//! (the bundled catalog).

mod analysis;
mod audit;
mod error;
pub mod grow;
mod integrate;
pub mod obreschkoff;
pub mod oeis_match;
mod report;
mod rows;
pub mod seed;
mod store;
mod types;

pub use analysis::{
    analyze_certificate, canonical_definition, definition_rows, same_sequence_from,
    CertificateAnalysis, PremiseRequirement,
};
pub use audit::{audit_atlas, audit_atlas_with, AtlasAuditOptions, AtlasAuditReport};
pub use error::AtlasError;
pub use integrate::{integrate_derivations, IntegrationOptions, IntegrationSummary};
pub use report::{atlas_report, AtlasReport};
pub use rows::{
    atlas_hash, degree_pattern, fingerprint_rows, integer_rows, observe_rows, rows_satisfy_claim,
};
pub use store::AtlasStore;
pub use types::*;
