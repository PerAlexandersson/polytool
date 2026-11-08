//! The bundled OEIS catalog of recurrence-backed polynomial families
//! (`polytool::oeis`), exposed as the Python module `polytool.oeis`.

use polytool::oeis::{self, OeisSequenceDefinition, OeisSequenceStatus};
use polytool::recurrence::RecurrenceJson;
use pyo3::exceptions::PyKeyError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::convert::{polytool_error, serialize_to_python};
use crate::recurrence::Recurrence;

fn entry(id: &str) -> PyResult<&'static OeisSequenceDefinition> {
    oeis::by_id(&id.to_ascii_uppercase())
        .ok_or_else(|| PyKeyError::new_err(format!("unknown bundled OEIS sequence: {id}")))
}

/// The record of an entry, as a `polytool.recurrence.v1` record indexed by
/// the displayed OEIS row numbers (`polytool oeis info --json`,
/// `recurrence_data`).
fn entry_record(entry: &OeisSequenceDefinition) -> PyResult<RecurrenceJson> {
    let (recurrence, initial_rows, first_index) =
        entry.recurrence_export_parts().map_err(polytool_error)?;
    let first_index = usize::try_from(first_index).map_err(|_| {
        polytool_error(format!(
            "{} has a negative displayed row index, which a recurrence record cannot represent",
            entry.id
        ))
    })?;
    Ok(RecurrenceJson::from_recurrence_rational(
        &recurrence,
        first_index,
        &initial_rows,
        None,
    ))
}

/// The A-numbers of the bundled catalog, in catalog order.
///
/// Experimental entries (recurrences without a held-out check) are left out
/// unless ``include_experimental=True``.
#[pyfunction]
#[pyo3(name = "ids", signature = (*, include_experimental = false))]
fn oeis_ids(include_experimental: bool) -> Vec<&'static str> {
    oeis::catalog()
        .iter()
        .filter(|entry| include_experimental || entry.status != OeisSequenceStatus::Experimental)
        .map(|entry| entry.id)
        .collect()
}

/// Information about one bundled sequence, as a ``dict``.
///
/// The keys follow ``polytool oeis info --json``: ``id``, ``name``,
/// ``status`` (``verified``, ``validated`` or ``experimental``), ``layout``,
/// ``first_index`` (the displayed OEIS row number of the first row; the CLI
/// key is ``first_row``), ``flattened_offset``, ``bfile_available``,
/// ``source_rows``, ``verification_rows``, ``rows_sha256``,
/// ``real_rootedness`` (the catalog flag, or ``None``), ``evidence``
/// (finite-segment evidence, or ``None``), ``recurrence`` (text) and
/// ``recurrence_data`` (the ``polytool.recurrence.v1`` record; see
/// :func:`recurrence`).
///
/// Raises ``KeyError`` for an A-number that is not bundled.
#[pyfunction]
#[pyo3(name = "info")]
fn oeis_info<'py>(py: Python<'py>, id: &str) -> PyResult<Bound<'py, PyDict>> {
    let entry = entry(id)?;
    let dict = PyDict::new(py);
    dict.set_item("id", entry.id)?;
    dict.set_item("name", entry.name)?;
    dict.set_item("status", entry.status.as_str())?;
    dict.set_item("layout", entry.layout.as_str())?;
    dict.set_item("first_index", entry.first_row)?;
    dict.set_item("flattened_offset", entry.flattened_offset)?;
    dict.set_item("bfile_available", entry.bfile_prefix_verified)?;
    dict.set_item("source_rows", entry.source_rows)?;
    dict.set_item("verification_rows", entry.verification_rows)?;
    dict.set_item("rows_sha256", entry.rows_sha256)?;
    match entry.real_rootedness() {
        Some(flag) => {
            let value = PyDict::new(py);
            value.set_item("status", flag.status.as_str())?;
            value.set_item("as_of", flag.as_of)?;
            dict.set_item("real_rootedness", value)?;
        }
        None => dict.set_item("real_rootedness", py.None())?,
    }
    match entry.evidence() {
        Some(evidence) => {
            let value = PyDict::new(py);
            value.set_item("strongest", evidence.strongest.as_str())?;
            value.set_item("rows_checked", evidence.rows_checked)?;
            value.set_item(
                "first_failure",
                serialize_to_python(py, &evidence.first_failure)?,
            )?;
            value.set_item("as_of", evidence.as_of)?;
            dict.set_item("evidence", value)?;
        }
        None => dict.set_item("evidence", py.None())?,
    }
    let (recurrence, _, _) = entry.recurrence_export_parts().map_err(polytool_error)?;
    dict.set_item("recurrence", recurrence.to_string())?;
    dict.set_item(
        "recurrence_data",
        serialize_to_python(py, &entry_record(entry)?)?,
    )?;
    Ok(dict)
}

/// Rows of a bundled sequence, exactly as the OEIS lists them.
///
/// Returns ``count`` rows (lists of ``int``, ascending powers of ``t``, with
/// the triangle's padding zeros), starting at the displayed OEIS row number
/// ``from_index`` (default: the first row).  Rows come from the
/// authoritative prefix and then from the bundled recurrence, as
/// ``polytool oeis generate``.  Wrap a row in :class:`polytool.Polynomial`
/// for arithmetic.  Experimental entries need ``include_experimental=True``.
///
/// >>> polytool.oeis.generate("A008292", 4)
/// [[1], [1, 1], [1, 4, 1], [1, 11, 11, 1]]
#[pyfunction]
#[pyo3(name = "generate", signature = (id, count, *, from_index = None, include_experimental = false))]
fn oeis_generate(
    py: Python<'_>,
    id: &str,
    count: usize,
    from_index: Option<i64>,
    include_experimental: bool,
) -> PyResult<Vec<Vec<num_bigint::BigInt>>> {
    let entry = entry(id)?;
    if entry.status == OeisSequenceStatus::Experimental && !include_experimental {
        return Err(polytool_error(format!(
            "{} is experimental; pass include_experimental=True to generate it",
            entry.id
        )));
    }
    let from_index = from_index.unwrap_or(entry.first_row);
    py.detach(|| entry.generate_rows_from(from_index, count))
        .map_err(polytool_error)
}

/// The bundled recurrence of a sequence as a :class:`polytool.Recurrence`,
/// indexed by the displayed OEIS row numbers: its initial polynomials are the
/// authoritative prefix and the recurrence's seed rows.
///
/// Rows of the recurrence are polynomials (no triangle padding); compare with
/// :func:`generate` after dropping trailing zeros.
#[pyfunction]
#[pyo3(name = "recurrence")]
fn oeis_recurrence(id: &str) -> PyResult<Recurrence> {
    Recurrence::from_record(entry_record(entry(id)?)?)
}

/// Build the `oeis` submodule.
pub(crate) fn module(py: Python<'_>) -> PyResult<Bound<'_, PyModule>> {
    let module = PyModule::new(py, "oeis")?;
    module.add(
        "__doc__",
        "The bundled catalog of recurrence-backed OEIS polynomial families.",
    )?;
    module.add_function(wrap_pyfunction!(oeis_ids, &module)?)?;
    module.add_function(wrap_pyfunction!(oeis_info, &module)?)?;
    module.add_function(wrap_pyfunction!(oeis_generate, &module)?)?;
    module.add_function(wrap_pyfunction!(oeis_recurrence, &module)?)?;
    Ok(module)
}
