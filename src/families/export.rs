//! Reproducible OEIS draft blocks, bounded b-files and proposal TOML skeletons.
use super::*;
use serde::Serialize;

/// Reproducible export; reviewing and submitting the draft remains a separate action.
#[derive(Clone, Debug, Serialize)]
pub struct FamilyOeisExport {
    /// Stable output schema.
    pub schema: String,
    /// Catalog ID, retained even after an A-number is assigned.
    pub family_id: String,
    /// Copyable internal-format blocks. A000000 is a placeholder for a new sequence.
    pub blocks: String,
    /// oeis-contributions proposal skeleton; fill authors, dates, and generator adapter.
    pub proposal_toml: String,
    /// B-file content, flattened indices and exact integers, one pair per line.
    pub bfile: String,
    /// Number of terms actually written, capped at 10,000 and 1 MiB.
    pub bfile_terms: usize,
    /// True when a b-file size or term limit omitted generated coefficients.
    pub bfile_truncated: bool,
}

fn plain(s: &str) -> String {
    s.split_whitespace().collect::<Vec<_>>().join(" ")
}
fn quoted(s: &str) -> String {
    // JSON basic strings are also TOML basic strings, apart from rare JSON escapes.
    serde_json::to_string(s).expect("serializable string")
}
fn array(items: &[String]) -> String {
    format!(
        "[{}]",
        items
            .iter()
            .map(|s| quoted(s))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

/// Export complete generated rows (1..=512) as an OEIS draft and bounded b-file.
/// `%S/%T/%U` occupy at most 80 bytes each (including the tag and ID) and
/// together at most 240 payload bytes; integers are never split. This conservative
/// subset fits OEIS's 260-character data field. The b-file has at most 10,000
/// terms and 1 MiB. It contains the requested rows, not an unbounded extension.
/// Keywords are inferred from the generated prefix only and need human review.
pub fn export_oeis(
    record: &FamilyRecord,
    row_count: usize,
) -> Result<FamilyOeisExport, FamilyCatalogError> {
    if row_count == 0 {
        return Err(failure("export needs at least one row"));
    }
    let rows = record.generate_rows(row_count)?;
    record
        .first_row
        .checked_add(i64::try_from(row_count - 1).map_err(failure)?)
        .ok_or_else(|| failure("displayed row index overflow"))?;
    let id = record.oeis.id.as_deref().unwrap_or("A000000");
    let mut blocks = format!("%N {id} {}\n", plain(&record.title));
    let mut comments = vec![];
    if !record.definition.is_empty() {
        comments.push(plain(&record.definition));
    }
    let p = &record.real_rootedness;
    comments.push(format!(
        "Real-rootedness status: {}. {}",
        serde_json::to_value(p.status).unwrap().as_str().unwrap(),
        plain(&p.statement)
    ));
    if let Some(range) = &p.finite_range {
        comments.push(format!(
            "Finite evidence for rows {}..{}: {}. Finite checks are not an all-index proof.",
            range.first,
            range.last,
            plain(&range.evidence).trim_end_matches('.')
        ));
    }
    for c in &comments {
        blocks.push_str(&format!("%C {id} {c}\n"));
    }
    let formulas: Vec<_> = [&record.formula, &record.generating_function]
        .into_iter()
        .filter(|s| !s.is_empty())
        .map(|s| plain(s))
        .collect();
    for formula in &formulas {
        blocks.push_str(&format!("%F {id} {formula}\n"));
    }
    for reference in &record.references {
        blocks.push_str(&format!(
            "%H {id} {} ({}) {}\n",
            plain(&reference.citation),
            plain(&reference.locator),
            reference.url
        ));
    }
    for link in &record.links {
        blocks.push_str(&format!("%H {id} {link}\n"));
    }
    let mut lines: Vec<String> = vec![String::new()];
    let width = 80 - format!("%S {id} ").len();
    let mut payload = 0;
    for term in rows.iter().flatten() {
        // Reserve a comma for a continued line; all selected terms form one prefix.
        let needed = term.len() + usize::from(!lines.last().unwrap().is_empty());
        if payload + needed > 240 || term.len() + 1 > width {
            break;
        }
        if lines.last().unwrap().len() + needed + 1 > width {
            if lines.len() == 3 {
                break;
            }
            lines.last_mut().unwrap().push(',');
            payload += 1;
            lines.push(String::new());
        }
        let line = lines.last_mut().unwrap();
        if !line.is_empty() {
            line.push(',');
            payload += 1;
        }
        line.push_str(term);
        payload += term.len();
    }
    if lines[0].is_empty() {
        return Err(failure("first term cannot fit an OEIS data line"));
    }
    for (tag, line) in ["S", "T", "U"].iter().zip(&lines) {
        blocks.push_str(&format!("%{tag} {id} {line}\n"));
    }
    let examples: Vec<_> = rows
        .iter()
        .take(7)
        .enumerate()
        .map(|(i, row)| format!("{}: {};", record.first_row + i as i64, row.join(", ")))
        .collect();
    for example in &examples {
        blocks.push_str(&format!("%e {id} {example}\n"));
    }
    let nonnegative = rows.iter().flatten().all(|s| !s.starts_with('-'));
    let keywords = vec![
        if nonnegative { "nonn" } else { "sign" }.to_string(),
        record.layout.clone(),
    ];
    let command = format!(
        "polytool families generate {} --count {row_count} --format terms",
        record.id
    );
    blocks.push_str(&format!("%o {id} (Polytool) {command}\n"));
    for r in &record.related {
        let tag = if r.id.starts_with('A') { "Y" } else { "C" };
        blocks.push_str(&format!(
            "%{tag} {id} Catalog relation {}: {}\n",
            r.id,
            plain(&r.relation)
        ));
    }
    blocks.push_str(&format!(
        "%K {id} {}\n%O {id} {},{}\n",
        keywords.join(","),
        record.flattened_offset,
        rows.iter()
            .flatten()
            .position(|s| {
                s.parse::<num_bigint::BigInt>().is_ok_and(|n| {
                    n > num_bigint::BigInt::from(1) || n < num_bigint::BigInt::from(-1)
                })
            })
            .map(|i| i + 1)
            .unwrap_or(1)
    ));
    let mut bfile = String::new();
    let mut terms = 0;
    let mut truncated = false;
    for term in rows.iter().flatten() {
        let index = record
            .flattened_offset
            .checked_add(i64::try_from(terms).map_err(failure)?)
            .ok_or_else(|| failure("b-file index overflow"))?;
        let line = format!("{index} {term}\n");
        if terms == 10_000 || bfile.len() + line.len() > 1024 * 1024 {
            truncated = true;
            break;
        }
        bfile.push_str(&line);
        terms += 1;
    }
    let target = record.oeis.id.as_deref().unwrap_or("new");
    let slug = format!("{}-{}", target, record.id.to_lowercase());
    let evidence = match p.status {
        FamilyProofStatus::LeanProved | FamilyProofStatus::TheoremProved => "proved",
        FamilyProofStatus::Conjectured => "conjectural",
        _ => "experimental",
    };
    let mut proposal = format!("# Save as proposals/{slug}.toml. Review every field before submission.\n# Add a repository-local generator adapter for: {command}\nid = {}\nkind = {}\ntarget = {}\nstatus = \"draft\"\nevidence = {}\nsummary = {}\nauthors = [\"TODO: author name\"]\n",quoted(&slug),quoted(if target=="new" {"new-sequence"} else {"comment"}),quoted(target),quoted(evidence),quoted(&record.title));
    proposal.push_str("\n[oeis_search]\n");
    proposal.push_str(&format!(
        "date = {}\nqueries = {}\nresult = {}\n",
        quoted(record.oeis.checked.as_deref().unwrap_or("TODO: YYYY-MM-DD")),
        array(&record.oeis.queries),
        quoted(if target == "new" {
            "no match (recheck before submission)"
        } else {
            "existing entry"
        })
    ));
    // Keep the single data field below the proposal validator's ~260-character limit.
    let data = lines.concat();
    proposal.push_str(&format!("\n[data]\ngenerator = \"generators/TODO-polytool-adapter.py\"\nargs = []\noffset = {}\nlayout = {}\nterms = {}\n",record.flattened_offset,quoted(&record.layout),quoted(&data)));
    proposal.push_str(&format!("\n[oeis]\nname = {}\ncomments = {}\nformulas = {}\nexample = {}\nprograms = {}\nlinks = {}\ncrossrefs = {}\nkeywords = {}\n", quoted(&record.title),array(&comments),array(&formulas),array(&examples),array(&[format!("(Polytool) {command}")]),array(&record.links),array(&record.related.iter().map(|r| format!("{}: {}",r.id,r.relation)).collect::<Vec<_>>()),array(&keywords)));
    Ok(FamilyOeisExport {
        schema: "polytool.family-oeis-export.v1".into(),
        family_id: record.id.clone(),
        blocks,
        proposal_toml: proposal,
        bfile,
        bfile_terms: terms,
        bfile_truncated: truncated,
    })
}
