//! Curated, submission-bound sequence families over the canonical generators.
//!
//! [`catalog`] derives OEIS entries without duplicating their data and merges
//! strict JSON sidecars. [`validate_catalog`] resolves generators, replays the
//! recorded initial rows and checks cross-references. [`FamilyRecord::generate_rows`]
//! returns exact integer strings in ascending powers, retaining OEIS row widths.
//! Proof statuses are external attestations; finite evidence is not a uniform proof.
mod export;
mod schema;
#[cfg(test)]
mod tests;
mod validation;
pub use export::*;
pub use schema::*;
pub use validation::*;

use crate::oeis::{self, OeisLayout, OeisRealRootednessStatus};
use num_rational::BigRational;
use std::collections::BTreeMap;

/// An invalid family, generator failure, or inconsistent catalog attestation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FamilyCatalogError(pub String);
impl std::fmt::Display for FamilyCatalogError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
impl std::error::Error for FamilyCatalogError {}

pub(super) fn failure(message: impl ToString) -> FamilyCatalogError {
    FamilyCatalogError(message.to_string())
}

include!(concat!(env!("OUT_DIR"), "/family_sidecars.rs"));

/// Parse one strict JSON record; call [`validate_catalog`] to resolve references and rows.
pub fn parse_family(text: &str) -> Result<FamilyRecord, FamilyCatalogError> {
    serde_json::from_str(text).map_err(failure)
}

/// Derive every OEIS record, then replace curated entries with bundled sidecars.
/// OEIS status, audit and date always come from the canonical OEIS catalog.
pub fn catalog() -> Result<Vec<FamilyRecord>, FamilyCatalogError> {
    let mut records: BTreeMap<String, FamilyRecord> = oeis::catalog()
        .iter()
        .map(|entry| {
            let property = inherited_property(entry);
            let record = FamilyRecord {
                schema_version: 1,
                id: entry.id.into(),
                title: entry.name.into(),
                definition: String::new(),
                formula: String::new(),
                generating_function: String::new(),
                generator: FamilyGenerator::Oeis {
                    id: entry.id.into(),
                },
                first_row: entry.first_row,
                flattened_offset: entry.flattened_offset,
                layout: if entry.layout == OeisLayout::RegularTriangle {
                    "tabl"
                } else {
                    "tabf"
                }
                .into(),
                degree_pattern: String::new(),
                initial_rows: vec![],
                references: vec![],
                real_rootedness: property,
                interlacing: FamilyProperty::default(),
                models: vec![],
                related: vec![],
                oeis: FamilyOeis {
                    id: Some(entry.id.into()),
                    checked: None,
                    queries: vec![],
                },
                submission: FamilySubmission {
                    state: "not_started".into(),
                    date: None,
                    a_number: None,
                },
                links: vec![],
                generator_status: Some(entry.status.as_str().into()),
            };
            (record.id.clone(), record)
        })
        .collect();
    let mut seen = std::collections::BTreeSet::new();
    for (filename, text) in FAMILY_SIDECARS {
        let mut record = parse_family(text)?;
        if *filename != format!("{}.json", record.id) || !seen.insert(record.id.clone()) {
            return Err(failure("sidecar filename mismatch or duplicate ID"));
        }
        if let FamilyGenerator::Oeis { id } = &record.generator {
            let entry =
                oeis::by_id(id).ok_or_else(|| failure(format!("unknown generator {id}")))?;
            // Curated fields supply citations, but never override the canonical attestation.
            let inherited = inherited_property(entry);
            record.real_rootedness.status = inherited.status;
            record.real_rootedness.as_of = inherited.as_of;
            record.real_rootedness.statement = inherited.statement;
            if record.real_rootedness.finite_range.is_none() {
                record.real_rootedness.finite_range = inherited.finite_range;
            }
            record.generator_status = Some(entry.status.as_str().into());
        }
        records.insert(record.id.clone(), record);
    }
    Ok(records.into_values().collect())
}

fn inherited_property(entry: &oeis::OeisSequenceDefinition) -> FamilyProperty {
    let mut property = FamilyProperty::default();
    if let Some(flag) = entry.real_rootedness() {
        property.status = match flag.status {
            OeisRealRootednessStatus::LeanProved => FamilyProofStatus::LeanProved,
            OeisRealRootednessStatus::TheoremProved => FamilyProofStatus::TheoremProved,
            OeisRealRootednessStatus::KnownNotRealRooted => FamilyProofStatus::KnownNotRealRooted,
        };
        property.as_of = flag.as_of.map(str::to_owned);
        property.statement = flag.audit.unwrap_or_default().into();
    }
    if let Some(evidence) = entry.evidence() {
        if evidence.strongest == oeis::OeisEvidenceProperty::RealRooted {
            if property.status == FamilyProofStatus::Unknown {
                property.status = FamilyProofStatus::FiniteEvidence;
                property.as_of = Some(evidence.as_of.into());
            }
            property.finite_range = i64::try_from(evidence.rows_checked)
                .ok()
                .and_then(|n| entry.first_row.checked_add(n - 1))
                .map(|last| FamilyFiniteRange {
                    first: entry.first_row,
                    last,
                    evidence: "Canonical OEIS catalog exact finite evidence".into(),
                });
        }
    }
    property
}

impl FamilyRecord {
    /// Generate the first `count` complete rows, starting at `first_row`.
    /// Counts are bounded to 512 (largest exp-rows rank 1000). Non-integral
    /// coefficients are errors, never rounded or silently truncated.
    pub fn generate_rows(&self, count: usize) -> Result<Vec<Vec<String>>, FamilyCatalogError> {
        if count > 512 {
            return Err(failure("row count exceeds 512"));
        }
        if count == 0 {
            return Ok(vec![]);
        }
        let integral =
            |rows: Vec<Vec<BigRational>>| -> Result<Vec<Vec<String>>, FamilyCatalogError> {
                rows.into_iter()
                    .map(|row| {
                        row.into_iter()
                            .map(|q| {
                                if q.is_integer() {
                                    Ok(q.to_integer().to_string())
                                } else {
                                    Err(failure("family coefficient is not integral"))
                                }
                            })
                            .collect()
                    })
                    .collect()
            };
        match &self.generator {
            FamilyGenerator::Oeis { id } => {
                let entry = oeis::by_id(id)
                    .ok_or_else(|| failure(format!("unknown OEIS generator {id}")))?;
                let rows = entry.generate_rows(count).map_err(failure)?;
                Ok(rows
                    .into_iter()
                    .map(|r| r.into_iter().map(|v| v.to_string()).collect())
                    .collect())
            }
            FamilyGenerator::ExpRows { family } => {
                let first = usize::try_from(self.first_row).map_err(failure)?;
                let last = first
                    .checked_add(count - 1)
                    .filter(|n| *n <= 1000)
                    .ok_or_else(|| failure("exp-rows rank exceeds 1000"))?;
                integral(
                    family
                        .rows_iter(last)
                        .map_err(failure)?
                        .skip(first)
                        .take(count)
                        .collect(),
                )
            }
            FamilyGenerator::Recurrence { record } => {
                let (rule, index, seeds) = record.to_recurrence_parts().map_err(failure)?;
                integral(
                    rule.generate_rows_rational(&seeds, index, count)
                        .map_err(failure)?,
                )
            }
        }
    }
}
