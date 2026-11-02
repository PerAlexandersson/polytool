//! `polytool oeis add --recurrence-json`: record a supplied recurrence after
//! checking it exactly against every available OEIS term.
//!
//! No search runs.  The supplied `polytool.recurrence.v1` object gives the
//! rule, the recurrence index of its first initial row and the initial rows;
//! its first initial row is matched with the first OEIS row after the
//! `--prefix-rows` authoritative prefix.  The candidate record is replayed
//! through the catalog's own row generator, flattened and compared term by
//! term with the OEIS data (the same comparison as `polytool oeis check`); the
//! entry is written only when every term matches.

use super::add::{source_record, string_rows, AddOptions, AddResult};
use super::catalog_file::{catalog_recurrence, generate_flattened, leak_definition};
use super::sha256::sha256_hex;
use super::source::RowWidths;
use polytool::oeis::{
    OeisCatalogProvenance, OeisCatalogRecord, OeisLayout, OeisProvenanceSource, OeisSequenceStatus,
    OeisValidationRow,
};
use polytool::recurrence::{format_rational_coeff, RecurrenceJson};

/// Fewest complete OEIS rows the comparison must cover.
const MIN_COMPARED_ROWS: usize = 5;

/// Fewest rows computed by the recurrence (after the prefix and the initial
/// rows) that the comparison must cover.
const MIN_COMPUTED_ROWS: usize = 2;

/// Upper bound on generated rows while covering the OEIS terms.
const MAX_GENERATED_ROWS: usize = 4096;

fn read_supplied(path: &str) -> Result<RecurrenceJson, String> {
    let text =
        std::fs::read_to_string(path).map_err(|error| format!("failed to read {path}: {error}"))?;
    serde_json::from_str(&text)
        .map_err(|error| format!("{path} is not a polytool.recurrence.v1 object: {error}"))
}

/// Validate the supplied recurrence for `id` and build its record.
pub(crate) fn add_supplied_recurrence(
    options: &AddOptions,
    id: &str,
    path: &str,
    present: bool,
) -> AddResult {
    if present && !options.replace {
        return AddResult::new(id, "already_present", "pass --replace to update it");
    }
    let supplied = match read_supplied(path) {
        Ok(supplied) => supplied,
        Err(error) => return AddResult::new(id, "unavailable", error),
    };
    let (rule, first_index, initial_rows) = match supplied.to_recurrence_parts() {
        Ok(parts) => parts,
        Err(error) => {
            return AddResult::new(
                id,
                "mismatch",
                format!("invalid supplied recurrence: {error}"),
            )
        }
    };
    if initial_rows.is_empty() {
        return AddResult::new(
            id,
            "mismatch",
            "the supplied recurrence has no initial rows",
        );
    }
    let source = match source_record(options, id) {
        Ok(source) => source,
        Err(error) => return AddResult::new(id, "unavailable", error),
    };
    let layout = match options.layout {
        Some(layout) => layout,
        None if source.has_keyword("tabl") => OeisLayout::RegularTriangle,
        None => {
            return AddResult::new(
                id,
                "layout_unsupported",
                format!(
                    "keywords `{}` lack `tabl`; pass --layout triangle or --layout table",
                    source.keywords.join(",")
                ),
            )
        }
    };
    let prefix_count = options.prefix_rows.unwrap_or(0);
    let widths = match (layout, options.row_widths) {
        (OeisLayout::RegularTriangle, _) => Some(RowWidths::TRIANGLE),
        (OeisLayout::Table, widths) => widths,
    };
    let prefix_rows = if prefix_count == 0 {
        Vec::new()
    } else {
        let Some(widths) = widths else {
            return AddResult::new(
                id,
                "layout_unsupported",
                "--prefix-rows with a table layout needs --row-widths W+S[/R]",
            );
        };
        let (rows, _) = widths.split(&source.terms);
        if rows.len() < prefix_count {
            return AddResult::new(
                id,
                "insufficient_data",
                format!(
                    "only {} complete rows for {prefix_count} prefix rows",
                    rows.len()
                ),
            );
        }
        string_rows(&rows[..prefix_count])
    };
    let recurrence_first_width = match layout {
        OeisLayout::RegularTriangle => prefix_count + 1,
        OeisLayout::Table => initial_rows[0].len().max(1),
    };
    let terms_text = source
        .terms
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let mut record = OeisCatalogRecord {
        id: id.to_string(),
        name: source.name.clone(),
        layout,
        first_row: source.offset,
        flattened_offset: source.offset,
        bfile_prefix_verified: true,
        prefix_rows,
        recurrence_first_index: first_index,
        recurrence_first_width,
        initial_rows: initial_rows
            .iter()
            .map(|row| row.iter().map(format_rational_coeff).collect())
            .collect(),
        recurrence: catalog_recurrence(&rule),
        provenance: OeisCatalogProvenance {
            source: OeisProvenanceSource::SuppliedRecurrence,
            status: OeisSequenceStatus::Validated,
            source_rows: 0,
            verification_rows: 0,
            rows_sha256: sha256_hex(terms_text.as_bytes()),
            fixture_slug: None,
            fixture_row_offset: None,
            validation_row: None,
            fit_rows: None,
            checked_terms: Some(source.terms.len()),
            terms_source: Some(source.terms_source.to_string()),
            search_options: None,
            polytool_version: Some(env!("CARGO_PKG_VERSION").to_string()),
        },
        real_rootedness: None,
        evidence: None,
    };
    if let Err(error) = record.validate() {
        return AddResult::new(id, "mismatch", format!("invalid catalog record: {error}"));
    }

    // Replay through the catalog's own generator and compare every term.
    let definition = leak_definition(&record);
    let rows = match generate_flattened(definition, source.terms.len(), MAX_GENERATED_ROWS) {
        Ok(rows) => rows,
        Err(error) => return AddResult::new(id, "mismatch", error),
    };
    let generated = rows.iter().flatten().collect::<Vec<_>>();
    if let Some(position) = (0..source.terms.len().min(generated.len()))
        .find(|&index| *generated[index] != source.terms[index])
    {
        return AddResult::new(
            id,
            "mismatch",
            format!(
                "OEIS term a({}) = {} but the supplied recurrence generates {}",
                source.offset + position as i64,
                source.terms[position],
                generated[position]
            ),
        );
    }
    if generated.len() < source.terms.len() {
        return AddResult::new(
            id,
            "insufficient_data",
            format!(
                "only {} of {} terms generated within {MAX_GENERATED_ROWS} rows",
                generated.len(),
                source.terms.len()
            ),
        );
    }
    let mut covered = 0;
    let mut complete_rows = 0usize;
    for row in &rows {
        if covered + row.len() > source.terms.len() {
            break;
        }
        covered += row.len();
        complete_rows += 1;
    }
    let computed_rows = complete_rows.saturating_sub(prefix_count + initial_rows.len());
    if complete_rows < MIN_COMPARED_ROWS || computed_rows < MIN_COMPUTED_ROWS {
        return AddResult::new(
            id,
            "insufficient_data",
            format!(
                "the OEIS terms cover {complete_rows} complete rows, {computed_rows} of them \
                 computed by the recurrence (need {MIN_COMPARED_ROWS} and {MIN_COMPUTED_ROWS})"
            ),
        );
    }
    let last_row = complete_rows - 1;
    record.provenance.source_rows = complete_rows;
    record.provenance.validation_row = Some(OeisValidationRow {
        index: last_row,
        terms: rows[last_row].iter().map(ToString::to_string).collect(),
    });
    let outcome = if present { "updated" } else { "added" };
    let mut result = AddResult::new(
        id,
        outcome,
        format!(
            "{rule} [supplied; {complete_rows} complete rows ({computed_rows} computed), {} terms \
             checked, validated]",
            source.terms.len()
        ),
    );
    result.record = Some(record);
    result
}
