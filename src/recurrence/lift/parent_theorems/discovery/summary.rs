//! Compact views of discovery evidence. Selection is shared with the full
//! report's status/statement computation; no fields are added to its wire schema.

use super::*;

/// Schema of [`RecurrenceProofSummary`].
pub const RECURRENCE_PROOF_SUMMARY_SCHEMA: &str = "polytool.recurrence-proof-summary.v1";

/// One failed or unresolved hypothesis, for a compact diagnostic.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProofNearMiss {
    /// Stable theorem name.
    pub theorem: String,
    /// Stable hypothesis name.
    pub hypothesis: String,
    /// Evidence excerpt, at most 200 Unicode characters.
    pub evidence: String,
}

/// A compact view of a discovery report, without changing its proof semantics.
/// A missing simplicity or interlacing claim is never promoted to one.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct RecurrenceProofSummary {
    /// [`RECURRENCE_PROOF_SUMMARY_SCHEMA`].
    pub schema: String,
    /// The complete report's status, including conditional and refuted states.
    pub status: String,
    /// Winning stage, or `None` when no route concludes anything.
    pub stage: Option<String>,
    /// Stable theorem name, prefixed for catalog, affine and closure routes.
    pub route: Option<String>,
    /// Concise statement; indices and residue-class scope are separate fields.
    pub conclusion: String,
    /// Every row at or after this index is nonzero and real-rooted.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted_from_index: Option<usize>,
    /// Whether the selected proof establishes simple roots.
    pub simple_roots: bool,
    /// Consecutive interlacing kind, when proved.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub interlacing: Option<RowInterlacing>,
    /// First index of the asserted consecutive interlacing.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub interlacing_from: Option<usize>,
    /// Orientation of affine interlacing, when applicable.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub interlacing_orientation: Option<String>,
    /// Lean declarations referenced by the selected route; no replay is claimed.
    pub lean_declarations: Vec<String>,
    /// Number of recorded theorem attempts of each status across stages.
    pub attempt_counts: BTreeMap<String, usize>,
    /// Interpretation boundary of the summary.
    pub scope: String,
    /// Certified root region of the selected proof.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub root_region: Option<RootInterval>,
    /// Infinitely many zero rows, with the same scope as the full report.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub real_rooted_outside_zero_residue_classes: Option<RealRootedOutsideZeroResidueClasses>,
    /// Attestation labels; present only for a conditional result.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub conditional_on: Vec<String>,
    /// Why no theorem applies.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub not_applicable_reason: Option<String>,
    /// True only for an exact refutation in the separate opt-in stage.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub refuted: Option<bool>,
    /// At most three failed/unresolved hypotheses when not proved.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub near_misses: Vec<ProofNearMiss>,
}

pub(super) enum ProofSource<'a> {
    Attempts(&'static str, &'a [RealRootednessTheoremReport]),
    Catalog(&'a crate::recurrence::CatalogLiftConclusion),
    Affine(&'a AffineCoordinateConclusion),
    Closure(&'a crate::recurrence::StructuralClosureConclusion),
}

impl ProofSource<'_> {
    pub(super) fn narrative(&self) -> &'static str {
        match self {
            Self::Attempts("direct", _) => "the supplied record",
            Self::Attempts("tail", _) => "the tail record",
            Self::Attempts(_, _) => "a certified alternative representation",
            Self::Catalog(_) => "a catalog lift",
            Self::Affine(_) => "an affine coordinate",
            Self::Closure(_) => "a general closure",
        }
    }

    pub(super) fn rows(&self) -> Option<ProvedRows> {
        match self {
            Self::Attempts(_, attempts) => proved_rows(attempts),
            Self::Catalog(c) => (c.parent_proof == CatalogParentProofKind::Polytool)
                .then_some(ProvedRows::EveryRow(c.from_index)),
            Self::Affine(c) => Some(match &c.zero_residue_classes {
                None => ProvedRows::EveryRow(c.from_index),
                Some(zero) => {
                    ProvedRows::OutsideZeroResidueClasses(RealRootedOutsideZeroResidueClasses {
                        from_index: c.from_index,
                        zero_residue_classes: zero.clone(),
                    })
                }
            }),
            Self::Closure(c) => Some(ProvedRows::EveryRow(c.from_index)),
        }
    }
}

impl RecurrenceProofDiscoveryReport {
    // One selection policy for the full report and the compact adapter view.
    pub(super) fn proof_source(&self) -> Option<ProofSource<'_>> {
        if proved_rows(&self.attempts).is_some() {
            return Some(ProofSource::Attempts("direct", &self.attempts));
        }
        if let Some(t) = &self.tail {
            if proved_rows(&t.attempts).is_some() {
                return Some(ProofSource::Attempts("tail", &t.attempts));
            }
        }
        if let Some(s) = &self.alternative_representations {
            for a in s.alternatives.iter().filter(|a| a.proved) {
                if proved_rows(&a.attempts).is_some() {
                    return Some(ProofSource::Attempts(
                        "alternative_representation",
                        &a.attempts,
                    ));
                }
            }
        }
        let catalog = self
            .catalog_lifts
            .as_ref()
            .and_then(|s| s.conclusion.as_ref());
        if let Some(c) = catalog {
            if c.parent_proof == CatalogParentProofKind::Polytool {
                return Some(ProofSource::Catalog(c));
            }
        }
        if let Some(c) = self
            .affine_coordinates
            .as_ref()
            .and_then(|s| s.conclusion.as_ref())
        {
            return Some(ProofSource::Affine(c));
        }
        if let Some(c) = self
            .closure_lifts
            .as_ref()
            .and_then(|s| s.conclusion.as_ref())
        {
            return Some(ProofSource::Closure(c));
        }
        catalog.map(ProofSource::Catalog)
    }

    /// Summarize the recorded evidence without running the search or replay.
    /// The full report remains authoritative and untrusted reports still need
    /// [`Self::verify`]. Serialized summaries are bounded to 2 KiB by shortening
    /// diagnostic excerpts only; proof claims and their scopes are retained.
    pub fn summary(&self) -> RecurrenceProofSummary {
        let mut summary = RecurrenceProofSummary {
            schema: RECURRENCE_PROOF_SUMMARY_SCHEMA.into(),
            status: self.status.clone(),
            stage: None,
            route: None,
            conclusion: match self.status.as_str() {
                "proved" if self.real_rooted_outside_zero_residue_classes.is_some() => {
                    "Nonzero and real-rooted outside the specified zero residue classes."
                }
                "proved" => "Every row from the stated index is nonzero and real-rooted.",
                "conditional_on_attestation" => {
                    "Conditional on external attestation; not proved by polytool."
                }
                "refuted" => {
                    "The all-row real-rootedness claim is refuted by an exact certificate."
                }
                _ => "No all-index real-rootedness proof established.",
            }
            .into(),
            real_rooted_from_index: self.real_rooted_from_index,
            simple_roots: false,
            interlacing: None,
            interlacing_from: None,
            interlacing_orientation: None,
            lean_declarations: vec![],
            attempt_counts: BTreeMap::new(),
            scope: "recurrence-defined sequence; Lean replay not attempted".into(),
            root_region: None,
            real_rooted_outside_zero_residue_classes: self
                .real_rooted_outside_zero_residue_classes
                .clone(),
            conditional_on: vec![],
            not_applicable_reason: self.not_applicable_reason.clone(),
            refuted: (self.status == "refuted").then_some(true),
            near_misses: vec![],
        };
        if self.status != "refuted" {
            match self.proof_source() {
                Some(ProofSource::Attempts(stage, attempts)) => {
                    // Match proved_rows: prefer an all-row conclusion, then the
                    // earliest start; this prevents strengthening a zero-class proof.
                    let a = attempts
                        .iter()
                        .filter(|a| a.is_proved() && a.conclusion.is_some())
                        .min_by_key(|a| {
                            let c = a.conclusion.as_ref().unwrap();
                            (c.has_zero_residue_classes(), c.first_index)
                        })
                        .unwrap();
                    let c = a.conclusion.as_ref().unwrap();
                    summary.stage = Some(stage.into());
                    summary.route = Some(a.theorem.as_str().into());
                    summary.simple_roots = c.simple_roots;
                    summary.interlacing = c.interlacing;
                    summary.interlacing_from = c.interlacing_from;
                    summary.root_region = c.root_region.clone();
                    summary.lean_declarations = a.lean_declarations.clone();
                }
                Some(ProofSource::Catalog(c)) => {
                    summary.stage = Some("catalog_lift".into());
                    let route = serde_json::to_value(c.route).unwrap();
                    summary.route = Some(format!("catalog_lift:{}", route.as_str().unwrap()));
                    if self.status == "conditional_on_attestation" {
                        summary.conditional_on = c.conditional_on.clone();
                    } else {
                        summary.root_region = Some(c.root_region.clone());
                        if let Some(i) = &c.interlacing {
                            summary.interlacing = Some(i.kind);
                            summary.interlacing_from = Some(i.from_index);
                            summary.interlacing_orientation = Some(
                                serde_json::to_value(i.orientation)
                                    .unwrap()
                                    .as_str()
                                    .unwrap()
                                    .into(),
                            );
                        }
                    }
                }
                Some(ProofSource::Affine(c)) => {
                    summary.stage = Some("affine_coordinate".into());
                    summary.route = Some(format!("affine_coordinate:{}", c.theorem.as_str()));
                    summary.simple_roots = c.simple_roots;
                    summary.interlacing = c.interlacing;
                    summary.interlacing_from = c.interlacing_from;
                    summary.interlacing_orientation = c
                        .interlacing_orientation
                        .map(|o| serde_json::to_value(o).unwrap().as_str().unwrap().into());
                    summary.root_region = c.root_region.clone();
                }
                Some(ProofSource::Closure(c)) => {
                    summary.stage = Some("closure".into());
                    let stage = self
                        .closure_lifts
                        .as_ref()
                        .and_then(|s| s.attempts.first())
                        .map(|a| a.stage.as_str())
                        .unwrap_or("structural");
                    summary.route = Some(format!("closure:{stage}"));
                    summary.root_region = Some(c.root_region.clone());
                    if let Some(i) = &c.interlacing {
                        summary.interlacing = Some(i.kind);
                        summary.interlacing_from = Some(i.from_index);
                        summary.interlacing_orientation = Some(
                            serde_json::to_value(i.orientation)
                                .unwrap()
                                .as_str()
                                .unwrap()
                                .into(),
                        );
                    }
                }
                None => {}
            }
        }
        let mut attempts: Vec<_> = self.attempts.iter().collect();
        if let Some(t) = &self.tail {
            attempts.extend(&t.attempts);
        }
        if let Some(s) = &self.alternative_representations {
            for a in &s.alternatives {
                attempts.extend(&a.attempts);
            }
        }
        if let Some(s) = &self.affine_coordinates {
            for a in &s.coordinates {
                attempts.extend(&a.attempts);
                if let Some(t) = &a.tail {
                    attempts.extend(&t.attempts);
                }
            }
        }
        for a in attempts {
            let status = serde_json::to_value(a.status).unwrap();
            *summary
                .attempt_counts
                .entry(status.as_str().unwrap().into())
                .or_default() += 1;
            if self.status != "proved" && summary.near_misses.len() < 3 {
                for h in &a.hypotheses {
                    if h.status != crate::recurrence::LiftHypothesisStatus::Verified
                        && summary.near_misses.len() < 3
                    {
                        summary.near_misses.push(ProofNearMiss {
                            theorem: a.theorem.as_str().into(),
                            hypothesis: h.name.clone(),
                            evidence: h.evidence.chars().take(200).collect(),
                        });
                    }
                }
            }
        }
        // Diagnostic excerpts can contain arbitrary supplied text. Drop them
        // rather than truncating a mathematical claim or a root-region value.
        while serde_json::to_vec(&summary).unwrap().len() > 2048 && !summary.near_misses.is_empty()
        {
            summary.near_misses.pop();
        }
        summary
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summaries_follow_the_frozen_theorem_wire_corpus() {
        let entry = crate::oeis::by_id("A008292").unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let record = RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            None,
        );
        let template = discover_recurrence_proofs(&record, &Default::default()).unwrap();
        let corpus = include_str!("../../../../../tests/fixtures/theorem-report-wire.jsonl");
        let mut checked = 0;
        for line in corpus.lines() {
            let case: serde_json::Value = serde_json::from_str(line).unwrap();
            let Ok(attempt) = serde_json::from_str::<RealRootednessTheoremReport>(
                case["input"].as_str().unwrap(),
            ) else {
                continue;
            };
            let mut report = template.clone();
            report.attempts = vec![attempt.clone()];
            report.finish_summary(None);
            let summary = report.summary();
            assert_eq!(
                serde_json::from_slice::<RecurrenceProofSummary>(
                    &serde_json::to_vec(&summary).unwrap()
                )
                .unwrap(),
                summary
            );
            assert_eq!(summary.status, report.status);
            assert_eq!(
                summary.real_rooted_from_index,
                report.real_rooted_from_index
            );
            assert_eq!(
                summary.real_rooted_outside_zero_residue_classes,
                report.real_rooted_outside_zero_residue_classes
            );
            if let (true, Some(conclusion)) = (attempt.is_proved(), &attempt.conclusion) {
                assert_eq!(summary.stage.as_deref(), Some("direct"));
                assert_eq!(summary.route.as_deref(), Some(attempt.theorem.as_str()));
                assert_eq!(summary.simple_roots, conclusion.simple_roots);
            } else {
                assert!(summary.route.is_none());
            }
            assert!(serde_json::to_vec(&summary).unwrap().len() <= 2048);
            checked += 1;
        }
        assert!(checked >= 20);
    }

    #[test]
    fn stage_selection_preserves_tail_conditional_and_closure_scope() {
        use crate::recurrence::{
            CatalogLiftConclusion, CatalogLiftRoute, StructuralClosureConclusion,
        };
        let entry = crate::oeis::by_id("A008292").unwrap();
        let (rule, seeds) = entry.recurrence_parts().unwrap();
        let record = RecurrenceJson::from_recurrence_rational(
            &rule,
            entry.recurrence_first_index,
            &seeds,
            None,
        );
        let mut report = discover_recurrence_proofs(&record, &Default::default()).unwrap();
        let direct = report.summary();
        report.tail = Some(RecurrenceProofTail {
            first_index: record.first_index,
            reason: "test tail".into(),
            recurrence: record.clone(),
            attempts: std::mem::take(&mut report.attempts),
        });
        report.finish_summary(None);
        let tail = report.summary();
        assert_eq!(tail.stage.as_deref(), Some("tail"));
        assert_eq!(tail.route, direct.route);
        assert_eq!(tail.real_rooted_from_index, direct.real_rooted_from_index);
        report.tail = None;
        let root_region = direct.root_region.unwrap();
        report.catalog_lifts = Some(CatalogLiftSearch {
            options: Default::default(),
            identifications: vec![],
            basis_maps: None,
            conclusion: Some(CatalogLiftConclusion {
                route: CatalogLiftRoute::NormalizedIdentification,
                catalog_family: Some("A008292".into()),
                basis_map: None,
                parent_proof: CatalogParentProofKind::CatalogAttestation,
                conditional_on: vec!["external-attestation".into()],
                from_index: 1,
                root_region: root_region.clone(),
                leading_coefficient: None,
                checked_rows: vec![],
                statement: "conditional test".into(),
                interlacing: None,
            }),
        });
        report.finish_summary(None);
        let conditional = report.summary();
        assert_eq!(conditional.status, "conditional_on_attestation");
        assert_eq!(conditional.conditional_on, vec!["external-attestation"]);
        assert!(conditional.real_rooted_from_index.is_none());
        assert!(!conditional.simple_roots);
        report.closure_lifts = Some(DiscoveryClosureSearch {
            recurrence: record,
            max_candidates: 1,
            depth: 1,
            candidates_tried: 1,
            attempts: vec![],
            conclusion: Some(StructuralClosureConclusion {
                from_index: 7,
                root_region,
                statement: "closure test".into(),
                interlacing: None,
            }),
        });
        report.finish_summary(None);
        let closure = report.summary();
        assert_eq!(closure.status, "proved");
        assert_eq!(closure.stage.as_deref(), Some("closure"));
        assert!(closure.route.unwrap().starts_with("closure:"));
        assert_eq!(closure.real_rooted_from_index, Some(7));
        assert!(closure.conditional_on.is_empty());
    }
}
