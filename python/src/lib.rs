//! Python bindings for polytool (the extension module `polytool._polytool`).
//!
//! A thin adapter over the library, like the CLI, the MCP server and the web
//! app: it converts Python values to the exact library types, calls the
//! library, and converts the results back.  No mathematics lives here.
//!
//! - [`convert`]: `int`/`fractions.Fraction` <-> `BigInt`/`BigRational`,
//!   coefficient lists (ascending degree order), and reports to `dict`
//!   through their stable JSON schemas;
//! - [`polynomial`]: the `Polynomial` class;
//! - [`properties`]: real roots, interlacing and coefficient properties;
//! - [`recurrence`]: the `Recurrence` class, `find_recurrence`,
//!   `prove_real_rooted` and `certify_leading_coefficient`;
//! - [`oeis`]: the bundled OEIS catalog (`polytool.oeis`);
//! - [`errors`]: `PolytoolError` and its subclasses.
//!
//! The pure-Python package in `python/polytool/` re-exports this module; the
//! user documentation is `doc/python.md`.

// The `///` comments of the Python-facing items are the Python docstrings
// (reStructuredText with `>>>` examples, checked by the pytest suite), so
// rustdoc's Markdown link resolution does not apply to them.
#![allow(rustdoc::broken_intra_doc_links)]

use pyo3::prelude::*;

pub mod convert;
pub mod errors;
pub mod oeis;
pub mod polynomial;
pub mod properties;
pub mod recurrence;

/// The extension module.  The public API is the package `polytool`.
#[pymodule]
fn _polytool(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("PolytoolError", py.get_type::<errors::PolytoolError>())?;
    m.add(
        "VerificationError",
        py.get_type::<errors::VerificationError>(),
    )?;
    m.add(
        "BudgetExhaustedError",
        py.get_type::<errors::BudgetExhaustedError>(),
    )?;

    m.add_class::<polynomial::Polynomial>()?;
    m.add_class::<properties::RealRootEnclosure>()?;
    m.add_class::<recurrence::Recurrence>()?;
    m.add_class::<recurrence::RecurrenceRows>()?;
    m.add_class::<recurrence::ProofReport>()?;
    m.add_class::<recurrence::LeadingCoefficientCertificate>()?;

    m.add_function(wrap_pyfunction!(properties::is_real_rooted, m)?)?;
    m.add_function(wrap_pyfunction!(properties::count_real_roots, m)?)?;
    m.add_function(wrap_pyfunction!(properties::isolate_real_roots, m)?)?;
    m.add_function(wrap_pyfunction!(properties::interlaces, m)?)?;
    m.add_function(wrap_pyfunction!(properties::interlaces_either_order, m)?)?;
    m.add_function(wrap_pyfunction!(properties::is_log_concave, m)?)?;
    m.add_function(wrap_pyfunction!(properties::is_ultra_log_concave, m)?)?;
    m.add_function(wrap_pyfunction!(properties::is_unimodal, m)?)?;
    m.add_function(wrap_pyfunction!(properties::is_palindromic, m)?)?;
    m.add_function(wrap_pyfunction!(properties::is_gamma_positive, m)?)?;
    m.add_function(wrap_pyfunction!(properties::gamma_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(properties::has_simple_roots, m)?)?;
    m.add_function(wrap_pyfunction!(properties::discriminant, m)?)?;
    m.add_function(wrap_pyfunction!(properties::resultant, m)?)?;
    m.add_function(wrap_pyfunction!(properties::properties, m)?)?;

    m.add_function(wrap_pyfunction!(recurrence::find_recurrence, m)?)?;
    m.add_function(wrap_pyfunction!(recurrence::prove_real_rooted, m)?)?;
    m.add_function(wrap_pyfunction!(
        recurrence::certify_leading_coefficient,
        m
    )?)?;

    m.add_submodule(&oeis::module(py)?)?;
    Ok(())
}
