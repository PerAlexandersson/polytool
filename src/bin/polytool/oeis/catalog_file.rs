//! The catalog data file (`data/oeis-catalog.jsonl`) as seen by the CLI:
//! loading, writing, and turning records into evaluable definitions.

use num_traits::Zero;
use polytool::oeis::{
    OeisCatalogRecord, OeisSequenceDefinition, SparseCoefficient, SparseRecurrence,
    SparseRecurrenceTerm,
};
use polytool::recurrence::{format_rational_coeff, BivarPoly, Recurrence, RecurrenceSign};
use std::path::{Path, PathBuf};

/// Location of the data file relative to the polytool source directory.
pub(crate) const DEFAULT_CATALOG_PATH: &str = "data/oeis-catalog.jsonl";

/// `--catalog FILE`, or `data/oeis-catalog.jsonl` when it exists in the
/// current directory.
pub(crate) fn resolve_catalog_path(explicit: Option<&str>) -> Option<PathBuf> {
    match explicit {
        Some(path) => Some(PathBuf::from(path)),
        None => {
            let path = PathBuf::from(DEFAULT_CATALOG_PATH);
            path.is_file().then_some(path)
        }
    }
}

pub(crate) fn load_catalog(path: &Path) -> Result<Vec<OeisCatalogRecord>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|error| format!("failed to read {}: {error}", path.display()))?;
    OeisCatalogRecord::parse_jsonl(&text).map_err(|error| format!("{}: {error}", path.display()))
}

/// Validate, sort and write the records, one canonical line each.  The file
/// is replaced atomically (written next to it, then renamed).
pub(crate) fn write_catalog(path: &Path, records: &mut [OeisCatalogRecord]) -> Result<(), String> {
    records.sort_by(|left, right| left.id.cmp(&right.id));
    OeisCatalogRecord::validate_all(records)?;
    let mut text = String::new();
    for record in records.iter() {
        text.push_str(&record.to_jsonl_line());
        text.push('\n');
    }
    let temporary = path.with_extension("jsonl.tmp");
    std::fs::write(&temporary, text)
        .map_err(|error| format!("failed to write {}: {error}", temporary.display()))?;
    std::fs::rename(&temporary, path)
        .map_err(|error| format!("failed to replace {}: {error}", path.display()))
}

fn leak_str(value: &str) -> &'static str {
    Box::leak(value.to_owned().into_boxed_str())
}

fn leak_rows(rows: &[Vec<String>]) -> &'static [&'static [&'static str]] {
    let rows = rows
        .iter()
        .map(|row| {
            let row: &'static [&'static str] =
                Box::leak(row.iter().map(|value| leak_str(value)).collect());
            row
        })
        .collect::<Vec<_>>();
    Box::leak(rows.into_boxed_slice())
}

fn leak_monomials(terms: &[(usize, usize, String)]) -> &'static [SparseCoefficient] {
    Box::leak(
        terms
            .iter()
            .map(|(n_degree, t_degree, value)| {
                SparseCoefficient::new(*n_degree, *t_degree, leak_str(value))
            })
            .collect(),
    )
}

/// Build the same definition the build script would generate for `record`.
///
/// The definition borrows `'static` data, so its strings are leaked; this is
/// meant for the short-lived CLI only.
pub(crate) fn leak_definition(record: &OeisCatalogRecord) -> &'static OeisSequenceDefinition {
    let provenance = &record.provenance;
    let terms = record
        .recurrence
        .terms
        .iter()
        .map(|term| SparseRecurrenceTerm {
            offset: term.offset,
            derivative_order: term.derivative_order,
            alternating_sign: term.alternating_sign,
            coefficient: leak_monomials(&term.coefficient),
        })
        .collect::<Vec<_>>();
    Box::leak(Box::new(OeisSequenceDefinition {
        id: leak_str(&record.id),
        name: leak_str(&record.name),
        status: provenance.status,
        layout: record.layout,
        first_row: record.first_row,
        flattened_offset: record.flattened_offset,
        bfile_prefix_verified: record.bfile_prefix_verified,
        fixture_slug: leak_str(provenance.fixture_slug.as_deref().unwrap_or("")),
        fixture_row_offset: record.fixture_row_offset(),
        prefix_rows: leak_rows(&record.prefix_rows),
        recurrence_first_index: record.recurrence_first_index,
        recurrence_first_width: record.recurrence_first_width,
        initial_rows: leak_rows(&record.initial_rows),
        recurrence: SparseRecurrence {
            terms: Box::leak(terms.into_boxed_slice()),
            denominator: record.recurrence.denominator.as_deref().map(leak_monomials),
            inhomogeneous: record
                .recurrence
                .inhomogeneous
                .as_deref()
                .map(leak_monomials),
        },
        source_rows: provenance.source_rows,
        verification_rows: provenance.verification_rows,
        rows_sha256: leak_str(&provenance.rows_sha256),
    }))
}

/// Nonzero monomials `(n_degree, t_degree, value)` of `poly`, ordered by
/// `n` degree and then `t` degree.
pub(crate) fn sparse_monomials(poly: &BivarPoly) -> Vec<(usize, usize, String)> {
    let mut terms = Vec::new();
    for (n_degree, row) in poly.coeffs.iter().enumerate() {
        for (t_degree, value) in row.iter().enumerate() {
            if !value.is_zero() {
                terms.push((n_degree, t_degree, format_rational_coeff(value)));
            }
        }
    }
    terms
}

/// Sparse catalog form of a recurrence found by the search.
pub(crate) fn catalog_recurrence(recurrence: &Recurrence) -> polytool::oeis::OeisCatalogRecurrence {
    polytool::oeis::OeisCatalogRecurrence {
        terms: recurrence
            .terms
            .iter()
            .map(|term| polytool::oeis::OeisCatalogTerm {
                offset: term.offset,
                derivative_order: term.deriv_order,
                alternating_sign: term.sign == RecurrenceSign::AlternatingN,
                coefficient: sparse_monomials(&term.coeff),
            })
            .collect(),
        denominator: recurrence.denominator.as_ref().map(sparse_monomials),
        inhomogeneous: recurrence.inhomogeneous.as_ref().map(sparse_monomials),
    }
}

/// Generate catalog rows of `definition` until at least `term_count`
/// flattened coefficients are available (or `max_rows` rows were produced).
pub(crate) fn generate_flattened(
    definition: &OeisSequenceDefinition,
    term_count: usize,
    max_rows: usize,
) -> Result<Vec<Vec<num_bigint::BigInt>>, String> {
    let mut row_count = 8usize.min(max_rows.max(1));
    loop {
        let rows = definition
            .generate_rows(row_count)
            .map_err(|error| error.to_string())?;
        let available = rows.iter().map(Vec::len).sum::<usize>();
        if available >= term_count || row_count >= max_rows {
            return Ok(rows);
        }
        row_count = (row_count * 2).min(max_rows);
    }
}
