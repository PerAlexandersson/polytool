//! General recurrence proof discovery. No catalogue, sequence IDs or attestations.

use super::*;

/// Schema for a replayable general recurrence proof search.
pub const RECURRENCE_PROOF_DISCOVERY_SCHEMA: &str = "polytool.recurrence-proof-discovery.v1";

/// Complete results of the general recurrence theorem checks.
///
/// Every successful attempt proves its stated conclusion for the sequence
/// defined by `recurrence`. A fitted recurrence still needs an all-index
/// identity for any independently defined sequence. Reports are Rust exact
/// certificates; `lean_replay` never claims they were elaborated in Lean.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RecurrenceProofDiscoveryReport {
    /// [`RECURRENCE_PROOF_DISCOVERY_SCHEMA`].
    pub schema: String,
    /// Defining recurrence and seeds, including untrusted search metadata.
    pub recurrence: RecurrenceJson,
    /// Largest induction-start delay attempted.
    pub maximum_start_delay: usize,
    /// Optional requested root interval (finite ends interpreted as open).
    pub requested_root_region: Option<RootInterval>,
    /// `proved` if at least one general theorem succeeded; `not_proved` otherwise.
    /// The latter does not disprove real-rootedness.
    pub status: String,
    /// What is, and is not, established by a successful search.
    pub statement: String,
    /// Every theorem attempt, including failed and unresolved hypotheses.
    pub attempts: Vec<RealRootednessTheoremReport>,
    /// Always `not_attempted`: exact Rust replay is separate from Lean replay.
    pub lean_replay: String,
}

/// Try all supported general theorem families, independently of sequence names.
///
/// Tries product induction, Favard/three-term, Liu--Wang first order,
/// derivative-lag, Ma--Wang interval and Liu--Wang three-term criteria.
/// This is theorem discovery for a supplied rule, not recurrence fitting.
/// Failed sufficient hypotheses do not disprove real-rootedness. There is no
/// catalogue fallback or external attestation route.
///
/// # Errors
///
/// Invalid recurrence/seeds, index overflow, or errors from a theorem check.
pub fn discover_recurrence_proofs(
    record: &RecurrenceJson,
    options: &RealRootednessTheoremOptions,
) -> Result<RecurrenceProofDiscoveryReport, RecurrenceError> {
    let (rule, first_index, seeds) = record.to_recurrence_parts()?;
    first_index
        .checked_add(seeds.len())
        .and_then(|n| n.checked_add(options.maximum_start_delay))
        .and_then(|n| n.checked_add(GENERATED_ROWS))
        .filter(|n| *n <= i64::MAX as usize)
        .ok_or_else(|| {
            RecurrenceError::InvalidInput("proof-search index range overflows".into())
        })?;
    if seeds.is_empty()
        || rule
            .terms
            .iter()
            .any(|t| t.offset == 0 || t.offset > seeds.len())
    {
        return Err(RecurrenceError::InvalidInput(
            "proof discovery needs seeds for every positive recurrence lag".into(),
        ));
    }
    let mut attempts = Vec::new();
    for theorem in [
        RealRootednessTheorem::ProductRecurrence,
        RealRootednessTheorem::ThreeTermRecurrence,
        RealRootednessTheorem::LiuWang,
        RealRootednessTheorem::LiuWangDerivativeLag,
        RealRootednessTheorem::MaWangInterval,
        RealRootednessTheorem::LiuWangThreeTerm,
    ] {
        attempts.push(check_with_first_proving_start(
            theorem,
            &rule,
            first_index,
            &seeds,
            options,
        )?);
    }
    let proved = attempts.iter().any(RealRootednessTheoremReport::is_proved);
    Ok(RecurrenceProofDiscoveryReport {
        schema: RECURRENCE_PROOF_DISCOVERY_SCHEMA.into(),
        recurrence: record.clone(),
        maximum_start_delay: options.maximum_start_delay,
        requested_root_region: options.root_region.as_ref().map(open_root_interval),
        status: if proved { "proved" } else { "not_proved" }.into(),
        statement: "Each proved attempt establishes its conclusion for the recurrence-defined sequence. Applying it to an independently defined sequence requires an all-index recurrence identity and the stated seeds; finite fitting or held-out checks do not establish that identity. No Lean replay was attempted.".into(),
        attempts,
        lean_replay: LIFT_LEAN_REPLAY_NOT_ATTEMPTED.into(),
    })
}

impl RecurrenceProofDiscoveryReport {
    /// Recompute the complete search and replay each theorem certificate.
    ///
    /// Rejects changed conclusions, dropped attempts, altered evidence labels,
    /// or certificates referring to a different rule or seeds.
    ///
    /// # Errors
    ///
    /// [`RecurrenceError::VerificationFailed`] on a mismatch, or an underlying
    /// recurrence/checker error. Verifying `not_proved` is not a proof.
    pub fn verify(&self) -> Result<(), RecurrenceError> {
        let root_region = self
            .requested_root_region
            .as_ref()
            .map(|r| {
                if !r.is_index_independent() {
                    return Err(RecurrenceError::VerificationFailed(
                        "root region must be index-independent".into(),
                    ));
                }
                r.at(0)
            })
            .transpose()?;
        let recomputed = discover_recurrence_proofs(
            &self.recurrence,
            &RealRootednessTheoremOptions {
                root_region,
                maximum_start_delay: self.maximum_start_delay,
            },
        )?;
        if self != &recomputed {
            return Err(RecurrenceError::VerificationFailed(
                "proof discovery report differs from recomputation".into(),
            ));
        }
        for attempt in &self.attempts {
            attempt.verify()?;
        }
        Ok(())
    }
}
