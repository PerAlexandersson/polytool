//! Structural validation and optional proof-owner declaration export checks.
use super::*;
use std::collections::BTreeSet;

fn a_number(id: &str) -> bool {
    id.len() == 7 && id.starts_with('A') && id[1..].bytes().all(|b| b.is_ascii_digit())
}
fn family_id(id: &str) -> bool {
    a_number(id)
        || (id.len() == 8 && id.starts_with("PT-F") && id[4..].bytes().all(|b| b.is_ascii_digit()))
}
fn date(value: &str) -> bool {
    if value.len() != 10 || value.as_bytes()[4] != b'-' || value.as_bytes()[7] != b'-' {
        return false;
    }
    let number = |s: &str| {
        s.bytes()
            .all(|b| b.is_ascii_digit())
            .then(|| s.parse::<u32>().ok())
            .flatten()
    };
    let (Some(year), Some(month), Some(day)) = (
        value.get(..4).and_then(number),
        value.get(5..7).and_then(number),
        value.get(8..).and_then(number),
    ) else {
        return false;
    };
    let leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 => {
            if leap {
                29
            } else {
                28
            }
        }
        _ => 0,
    };
    year > 0 && day > 0 && day <= days
}
fn public_url(url: &str) -> bool {
    url.starts_with("https://")
        && url.len() > 8
        && !url.chars().any(char::is_whitespace)
        && !url.contains('@')
        && !url.contains('<')
        && !url.contains('>')
}

/// Check schema constraints, generator identity, row replay and local references.
/// Every generator is exercised for at least six rows; this is a consistency
/// check, not a proof of its combinatorial interpretation or root properties.
pub fn validate_catalog(records: &[FamilyRecord]) -> Result<(), FamilyCatalogError> {
    let mut ids = BTreeSet::new();
    for record in records {
        if !family_id(&record.id) || !ids.insert(record.id.as_str()) {
            return Err(failure("invalid or duplicate family ID"));
        }
    }
    for r in records {
        let ensure = |condition: bool, message: &str| {
            if condition {
                Ok(())
            } else {
                Err(failure(format!("{}: {message}", r.id)))
            }
        };
        ensure(
            r.schema_version == 1 && !r.title.trim().is_empty(),
            "schema_version must be 1 and title nonempty",
        )?;
        ensure(
            ["tabl", "tabf"].contains(&r.layout.as_str()),
            "layout must be tabl or tabf",
        )?;
        ensure(r.initial_rows.len() <= 512, "too many initial rows")?;
        if let FamilyGenerator::Oeis { id } = &r.generator {
            let e = oeis::by_id(id)
                .ok_or_else(|| failure(format!("{}: unresolved OEIS generator", r.id)))?;
            ensure(
                r.id == *id
                    && r.oeis.id.as_ref() == Some(id)
                    && r.first_row == e.first_row
                    && r.flattened_offset == e.flattened_offset
                    && r.layout
                        == if e.layout == OeisLayout::RegularTriangle {
                            "tabl"
                        } else {
                            "tabf"
                        },
                "OEIS generator identity, layout or indices disagree",
            )?;
        } else {
            ensure(
                r.id.starts_with("PT-F"),
                "A-number records must use their canonical OEIS generator",
            )?;
        }
        ensure(
            r.real_rootedness.lag.is_none(),
            "row property cannot have an interlacing lag",
        )?;
        if r.interlacing.status != FamilyProofStatus::Unknown {
            ensure(
                r.interlacing.lag.is_some(),
                "interlacing claim needs a positive lag",
            )?;
        }
        let mut keys = BTreeSet::new();
        for reference in &r.references {
            ensure(
                !reference.key.trim().is_empty()
                    && keys.insert(reference.key.as_str())
                    && !reference.citation.trim().is_empty()
                    && !reference.locator.trim().is_empty()
                    && public_url(&reference.url),
                "invalid or duplicate reference",
            )?;
        }
        for property in [&r.real_rootedness, &r.interlacing] {
            ensure(
                property.as_of.as_deref().is_none_or(date),
                "invalid property date",
            )?;
            if let Some(key) = &property.theorem {
                ensure(keys.contains(key.as_str()), "dangling theorem reference")?;
            }
            if let Some(range) = &property.finite_range {
                ensure(
                    range.first >= r.first_row
                        && range.first <= range.last
                        && !range.evidence.trim().is_empty(),
                    "invalid finite range",
                )?;
            }
            ensure(property.lag != Some(0), "lag must be positive")?;
            if let Some(lean) = &property.lean {
                // An explicit public allowlist prevents accidental private repository links.
                ensure(
                    lean.repository == "https://github.com/PerAlexandersson/RealRooted"
                        && lean.commit.len() == 40
                        && lean
                            .commit
                            .bytes()
                            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                        && !lean.declaration.is_empty()
                        && lean
                            .declaration
                            .chars()
                            .all(|c| c.is_alphanumeric() || "_.'".contains(c))
                        && lean
                            .permalink
                            .starts_with(&format!("{}/blob/{}/", lean.repository, lean.commit))
                        && public_url(&lean.permalink),
                    "invalid public pinned Lean attestation",
                )?;
            }
            if matches!(
                property.status,
                FamilyProofStatus::TheoremProved | FamilyProofStatus::LeanProved
            ) && !matches!(r.generator, FamilyGenerator::Oeis { .. })
            {
                ensure(
                    property.theorem.is_some(),
                    "new proof claim needs a theorem reference",
                )?;
            }
            if matches!(
                property.status,
                FamilyProofStatus::Conjectured | FamilyProofStatus::FiniteEvidence
            ) {
                ensure(
                    !property.statement.trim().is_empty() || property.finite_range.is_some(),
                    "conjecture/evidence needs a statement or range",
                )?;
            }
        }
        for model in &r.models {
            ensure(
                !model.family.trim().is_empty()
                    && ["pending", "passed", "failed"].contains(&model.replay_status.as_str()),
                "invalid model replay status",
            )?;
        }
        for relation in &r.related {
            ensure(
                ids.contains(relation.id.as_str())
                    && relation.id != r.id
                    && !relation.relation.trim().is_empty(),
                "dangling or empty related family",
            )?;
        }
        ensure(
            r.oeis.id.as_deref().is_none_or(a_number) && r.oeis.checked.as_deref().is_none_or(date),
            "invalid OEIS ID/date",
        )?;
        if r.oeis.id.is_none() {
            ensure(
                r.oeis.checked.is_some()
                    && !r.oeis.queries.is_empty()
                    && r.oeis.queries.iter().all(|q| !q.trim().is_empty()),
                "non-OEIS family needs dated queries",
            )?;
        }
        ensure(
            [
                "not_started",
                "draft",
                "ready",
                "submitted",
                "accepted",
                "withdrawn",
            ]
            .contains(&r.submission.state.as_str())
                && r.submission.date.as_deref().is_none_or(date)
                && r.submission.a_number.as_deref().is_none_or(a_number),
            "invalid submission state/date/A-number",
        )?;
        if ["submitted", "accepted"].contains(&r.submission.state.as_str()) {
            ensure(r.submission.date.is_some(), "submitted family needs date")?;
        }
        if r.submission.state == "accepted" {
            ensure(
                r.submission.a_number.is_some() && r.submission.a_number == r.oeis.id,
                "accepted family needs consistent A-number",
            )?;
        }
        ensure(
            r.links.iter().all(|url| public_url(url)),
            "invalid external link",
        )?;
        let rows = r.generate_rows(6.max(r.initial_rows.len()))?;
        ensure(
            rows.starts_with(&r.initial_rows),
            "initial rows do not reproduce",
        )?;
    }
    Ok(())
}

/// Optional pinned export record supplied by the repository that checks the proof.
/// This is a future-facing hook; no network lookup or private repository access occurs.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FamilyLeanDeclaration {
    /// Family tag owned by the proving repository.
    pub family_id: String,
    /// Fully qualified public declaration name.
    pub declaration: String,
    /// real_rootedness or interlacing.
    pub property: String,
    /// Repository URL.
    pub repository: String,
    /// Pinned commit of the export and source.
    pub commit: String,
}

/// Reject manual links absent from a supplied proof-owner declaration list.
/// The list must correspond to the pinned commits; fetching it is the caller's job.
pub fn validate_lean_declarations(
    records: &[FamilyRecord],
    declarations: &[FamilyLeanDeclaration],
) -> Result<(), FamilyCatalogError> {
    for r in records {
        for (name, property) in [
            ("real_rootedness", &r.real_rootedness),
            ("interlacing", &r.interlacing),
        ] {
            if let Some(lean) = &property.lean {
                if !declarations.iter().any(|d| {
                    d.family_id == r.id
                        && d.property == name
                        && d.declaration == lean.declaration
                        && d.repository == lean.repository
                        && d.commit == lean.commit
                }) {
                    return Err(failure(format!(
                        "{}: stale Lean declaration {}",
                        r.id, lean.declaration
                    )));
                }
            }
        }
    }
    Ok(())
}
