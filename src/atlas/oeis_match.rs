//! OEIS matching of the nodes against a local snapshot of the OEIS.
//!
//! # The terms cache
//!
//! [`build_oeis_terms_cache`] reads a local copy of the OEIS in its internal
//! format (one `A123456.seq` file per entry under `DIR/seq/A123/`, or
//! directly under `DIR/A123/`), keeps the `%S/%T/%U` terms, the offset
//! (`%O`) and the keywords (`%K`) of every entry, and writes one
//! tab-separated line per entry (`A-number`, offset, keywords, terms) to a
//! cache file outside the data directory.  [`OeisTermsIndex::load`] reads the
//! cache and indexes every entry by a hash of the first
//! [`OeisMatchOptions::key_terms`] terms starting at each of its first
//! positions.
//!
//! # Matching
//!
//! For every node, [`match_node`] flattens the stored rows from each of the
//! first few rows on, in four forms (ascending coefficients and the reversed
//! rows, each with and without the zero coefficients below the lowest
//! nonzero one), and the row sums `P_n(1)`, looks the first terms up and
//! compares every overlapping term exactly.  A hit needs at least
//! [`OeisMatchOptions::min_terms`] agreeing terms with at least three
//! distinct values.  Hits are links ([`AtlasLinkKind::OeisMatch`]) with their
//! alignment: finite agreement only.  A row-sum hit is never an
//! identification, and neither is a triangle hit (the identification of a
//! node with a catalog entry is a definition of the node, made only by an
//! equality certificate).

use super::error::AtlasError;
use super::integrate::{add_link, IntegrationSummary};
use super::rows::{atlas_hash, rows_from_strings, strip_low_zeros};
use super::store::AtlasStore;
use super::types::{AtlasLinkKind, AtlasNode, OeisAlignment, OeisMatchTarget};
use num_bigint::BigInt;
use num_traits::Zero;
use std::collections::{BTreeSet, HashMap};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Version tag written as the first line of the cache.
pub const OEIS_TERMS_CACHE_VERSION: &str = "#polytool-oeis-terms-v1";

/// One OEIS entry of the cache.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OeisTermsEntry {
    /// The A-number.
    pub id: String,
    /// The offset (first value of `%O`).
    pub offset: i64,
    /// The keywords.
    pub keywords: Vec<String>,
    /// The terms, comma separated.
    pub terms: String,
}

impl OeisTermsEntry {
    fn terms(&self) -> Vec<BigInt> {
        self.terms
            .split(',')
            .filter_map(|t| t.trim().parse().ok())
            .collect()
    }
}

/// Parse one internal-format file: `(offset, keywords, terms)`.
fn parse_internal(text: &str) -> Option<(i64, Vec<String>, Vec<String>)> {
    let mut offset = None;
    let mut keywords = Vec::new();
    let mut terms = Vec::new();
    for line in text.lines() {
        let mut parts = line.splitn(3, ' ');
        let (Some(tag), Some(_), Some(rest)) = (parts.next(), parts.next(), parts.next()) else {
            continue;
        };
        match tag {
            "%S" | "%T" | "%U" => terms.extend(
                rest.split(',')
                    .map(str::trim)
                    .filter(|t| !t.is_empty())
                    .map(str::to_string),
            ),
            "%O" => offset = rest.split(',').next().and_then(|v| v.trim().parse().ok()),
            "%K" => keywords = rest.split(',').map(|k| k.trim().to_string()).collect(),
            _ => {}
        }
    }
    let valid = terms.iter().all(|t| {
        let digits = t.strip_prefix('-').unwrap_or(t);
        !digits.is_empty() && digits.bytes().all(|b| b.is_ascii_digit())
    });
    (valid && !terms.is_empty()).then(|| (offset.unwrap_or(0), keywords, terms))
}

fn sorted_entries(directory: &Path) -> Result<Vec<std::path::PathBuf>, AtlasError> {
    let mut entries: Vec<_> = std::fs::read_dir(directory)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .collect();
    entries.sort();
    Ok(entries)
}

/// Build the terms cache from a local OEIS copy (module documentation);
/// returns the number of entries written.
///
/// # Errors
///
/// [`AtlasError::Io`] when the snapshot cannot be read or the cache cannot
/// be written.
pub fn build_oeis_terms_cache(oeis_data: &Path, cache: &Path) -> Result<usize, AtlasError> {
    let root = if oeis_data.join("seq").is_dir() {
        oeis_data.join("seq")
    } else {
        oeis_data.to_path_buf()
    };
    if let Some(parent) = cache.parent() {
        std::fs::create_dir_all(parent)?;
    }
    let temporary = cache.with_extension("partial");
    let mut out = BufWriter::new(std::fs::File::create(&temporary)?);
    writeln!(out, "{OEIS_TERMS_CACHE_VERSION}")?;
    let mut count = 0;
    for directory in sorted_entries(&root)? {
        if !directory.is_dir() {
            continue;
        }
        for file in sorted_entries(&directory)? {
            let Some(name) = file.file_name().and_then(|n| n.to_str()) else {
                continue;
            };
            let Some(id) = name.strip_suffix(".seq") else {
                continue;
            };
            let Ok(text) = std::fs::read_to_string(&file) else {
                continue;
            };
            if let Some((offset, keywords, terms)) = parse_internal(&text) {
                writeln!(
                    out,
                    "{id}\t{offset}\t{}\t{}",
                    keywords.join(","),
                    terms.join(",")
                )?;
                count += 1;
            }
        }
    }
    out.flush()?;
    drop(out);
    std::fs::rename(&temporary, cache)?;
    Ok(count)
}

/// Options of the matching.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OeisMatchOptions {
    /// Terms hashed per key.  Default 8.
    pub key_terms: usize,
    /// Least number of agreeing terms of a hit.  Default 10.
    pub min_terms: usize,
    /// OEIS positions indexed per entry (`0..=max_oeis_start`).  Default 3.
    pub max_oeis_start: usize,
    /// Node rows tried as the first row (`0..=max_row_start`).  Default 3.
    pub max_row_start: usize,
    /// Largest number of hits kept per node and target (lowest A-numbers
    /// first).  Default 8.
    pub max_hits_per_target: usize,
}

impl Default for OeisMatchOptions {
    fn default() -> Self {
        Self {
            key_terms: 8,
            min_terms: 10,
            max_oeis_start: 3,
            max_row_start: 3,
            max_hits_per_target: 8,
        }
    }
}

/// The loaded cache with its key index.
#[derive(Debug, Clone, Default)]
pub struct OeisTermsIndex {
    entries: Vec<OeisTermsEntry>,
    keys: HashMap<String, Vec<(u32, u8)>>,
}

fn key_of(terms: &[String]) -> String {
    atlas_hash(&terms.join(","))
}

impl OeisTermsIndex {
    /// Load a cache written by [`build_oeis_terms_cache`] and index it.
    ///
    /// # Errors
    ///
    /// [`AtlasError::Io`] or [`AtlasError::InvalidRecord`] for an unreadable
    /// or foreign file.
    pub fn load(cache: &Path, options: &OeisMatchOptions) -> Result<Self, AtlasError> {
        let reader = BufReader::new(std::fs::File::open(cache)?);
        let mut index = Self::default();
        for (number, line) in reader.lines().enumerate() {
            let line = line?;
            if number == 0 {
                if line != OEIS_TERMS_CACHE_VERSION {
                    return Err(AtlasError::InvalidRecord(
                        "the OEIS terms cache has an unknown version".into(),
                    ));
                }
                continue;
            }
            let mut parts = line.splitn(4, '\t');
            let (Some(id), Some(offset), Some(keywords), Some(terms)) =
                (parts.next(), parts.next(), parts.next(), parts.next())
            else {
                continue;
            };
            let keywords: Vec<String> = keywords
                .split(',')
                .filter(|k| !k.is_empty())
                .map(str::to_string)
                .collect();
            if keywords.iter().any(|k| k == "dead") {
                continue;
            }
            let split: Vec<String> = terms.split(',').map(str::to_string).collect();
            let position = index.entries.len() as u32;
            for start in 0..=options.max_oeis_start {
                if split.len() >= start + options.key_terms {
                    index
                        .keys
                        .entry(key_of(&split[start..start + options.key_terms]))
                        .or_default()
                        .push((position, start as u8));
                }
            }
            index.entries.push(OeisTermsEntry {
                id: id.to_string(),
                offset: offset.parse().unwrap_or(0),
                keywords,
                terms: terms.to_string(),
            });
        }
        Ok(index)
    }

    /// Number of entries.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// One hit of [`match_node`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OeisHit {
    /// The A-number.
    pub id: String,
    /// The alignment.
    pub alignment: OeisAlignment,
}

fn flattened(rows: &[Vec<BigInt>], reversed: bool, drop_low_zeros: bool) -> Vec<BigInt> {
    let mut out = Vec::new();
    for row in rows {
        let row: &[BigInt] = if drop_low_zeros {
            strip_low_zeros(row)
        } else {
            row
        };
        if reversed {
            out.extend(row.iter().rev().cloned());
        } else {
            out.extend(row.iter().cloned());
        }
    }
    out
}

/// The OEIS hits of one node (module documentation), deterministic.
pub fn match_node(
    node: &AtlasNode,
    index: &OeisTermsIndex,
    options: &OeisMatchOptions,
) -> Vec<OeisHit> {
    let Ok(rows) = rows_from_strings(&node.rows) else {
        return Vec::new();
    };
    let mut hits = Vec::new();
    let mut seen: BTreeSet<(String, String)> = BTreeSet::new();
    let targets = [
        (OeisMatchTarget::Triangle, false, false),
        (OeisMatchTarget::Triangle, false, true),
        (OeisMatchTarget::ReversedTriangle, true, false),
        (OeisMatchTarget::ReversedTriangle, true, true),
    ];
    let mut per_target: HashMap<String, usize> = HashMap::new();
    let mut consider = |target: OeisMatchTarget,
                        first_row: usize,
                        low_zeros_dropped: bool,
                        sequence: &[BigInt],
                        hits: &mut Vec<OeisHit>| {
        if sequence.len() < options.key_terms {
            return;
        }
        let key_text: Vec<String> = sequence[..options.key_terms]
            .iter()
            .map(ToString::to_string)
            .collect();
        let Some(candidates) = index.keys.get(&key_of(&key_text)) else {
            return;
        };
        let mut candidates = candidates.clone();
        candidates.sort_by(|a, b| {
            index.entries[a.0 as usize]
                .id
                .cmp(&index.entries[b.0 as usize].id)
        });
        for (position, start) in candidates {
            let entry = &index.entries[position as usize];
            let terms = entry.terms();
            let start = start as usize;
            let available = terms.len().saturating_sub(start).min(sequence.len());
            if available < options.min_terms {
                continue;
            }
            let compared = &sequence[..available];
            if terms[start..start + available] != *compared {
                continue;
            }
            let distinct: BTreeSet<&BigInt> = compared.iter().collect();
            if distinct.len() < 3 {
                continue;
            }
            let target_key = format!("{target:?}");
            if !seen.insert((entry.id.clone(), target_key.clone())) {
                continue;
            }
            let count = per_target.entry(target_key).or_default();
            if *count >= options.max_hits_per_target {
                continue;
            }
            *count += 1;
            hits.push(OeisHit {
                id: entry.id.clone(),
                alignment: OeisAlignment {
                    target,
                    first_row: node.first_index + first_row,
                    low_zeros_dropped,
                    oeis_term_position: start,
                    oeis_offset: entry.offset,
                    terms_compared: available,
                },
            });
        }
    };
    for first_row in 0..=options.max_row_start.min(rows.len().saturating_sub(1)) {
        for (target, reversed, drop) in targets {
            let sequence = flattened(&rows[first_row..], reversed, drop);
            consider(target, first_row, drop, &sequence, &mut hits);
        }
        let sums: Vec<BigInt> = rows[first_row..]
            .iter()
            .map(|row| row.iter().fold(BigInt::zero(), |s, c| s + c))
            .collect();
        consider(OeisMatchTarget::RowSums, first_row, false, &sums, &mut hits);
    }
    hits
}

/// Counts of [`match_atlas`].
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct OeisMatchSummary {
    /// Nodes matched.
    pub nodes: usize,
    /// Hits by target.
    pub hits_by_target: std::collections::BTreeMap<String, usize>,
    /// Links written.
    pub links_added: usize,
    /// Nodes with a triangle or reversed-triangle hit of an entry that is not
    /// in the bundled catalog, with the A-numbers.
    pub new_triangle_matches: Vec<(String, Vec<String>)>,
}

/// Match every node and record the hits as links (written on
/// [`AtlasStore::commit`]).
pub fn match_atlas(
    store: &mut AtlasStore,
    index: &OeisTermsIndex,
    options: &OeisMatchOptions,
) -> OeisMatchSummary {
    let mut summary = OeisMatchSummary::default();
    let mut integration = IntegrationSummary::default();
    let nodes: Vec<AtlasNode> = store.nodes().to_vec();
    for node in &nodes {
        summary.nodes += 1;
        let hits = match_node(node, index, options);
        let mut new_ids = BTreeSet::new();
        for hit in hits {
            let in_catalog = crate::oeis::by_id(&hit.id).is_some();
            *summary
                .hits_by_target
                .entry(format!("{:?}", hit.alignment.target))
                .or_default() += 1;
            if !in_catalog && hit.alignment.target != OeisMatchTarget::RowSums {
                new_ids.insert(hit.id.clone());
            }
            add_link(
                store,
                vec![node.id.clone()],
                AtlasLinkKind::OeisMatch {
                    oeis_id: hit.id,
                    in_catalog,
                    alignment: hit.alignment,
                },
                &mut integration,
            );
        }
        if !new_ids.is_empty() {
            summary
                .new_triangle_matches
                .push((node.id.clone(), new_ids.into_iter().collect()));
        }
    }
    summary.links_added = integration.links_added;
    summary
}
