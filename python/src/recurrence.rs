//! Recurrences: the `polytool.Recurrence` class (a `polytool.recurrence.v1`
//! record), adaptive recurrence search, the general real-rootedness proof
//! search and the degree-law / leading-coefficient certificates.

use polytool::real_rootedness::RealInterval;
use polytool::recurrence::{
    certify_leading_coefficient as library_certify_leading_coefficient, discover_recurrence_proofs,
    find_recurrence_adaptive_rational_with_budget, AdaptiveSearchBudget, AdaptiveSearchOptions,
    AdaptiveSearchOutcome, AdaptiveSearchResult, BigRational,
    LeadingCoefficientCertificate as LibraryLeadingCoefficientCertificate,
    LeadingCoefficientOptions, LeadingCoefficientSign, RealRootednessTheoremOptions,
    RealRootednessTheoremReport, RealRootednessTheoremStatus, Recurrence as LibraryRecurrence,
    RecurrenceError, RecurrenceJson, RecurrenceJsonSearch, RecurrenceOptionsJson,
    RecurrenceProofDiscoveryReport, RecurrenceRowIterator,
};
use polytool::renewal::RecurrenceTheoremReport;
use pyo3::exceptions::{PyStopIteration, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::convert::{
    json_text_from_python, polynomial_rows_from_python, polytool_error, rational_from_python,
    rational_to_python, serialize_to_python, to_json_text,
};
use crate::errors::{BudgetExhaustedError, VerificationError};
use crate::polynomial::Polynomial;

/// A recurrence record parsed into the library's rule, first index and seeds.
struct ParsedRecord {
    record: RecurrenceJson,
    rule: LibraryRecurrence,
    first_index: usize,
    seeds: Vec<Vec<BigRational>>,
}

impl ParsedRecord {
    fn new(record: RecurrenceJson) -> PyResult<Self> {
        let (rule, first_index, seeds) = record.to_recurrence_parts().map_err(polytool_error)?;
        Ok(Self {
            record,
            rule,
            first_index,
            seeds,
        })
    }
}

/// A polynomial sequence defined by a linear recurrence and initial rows.
///
/// The rule reads ``f(n, t) P_n(t) = sum_j s_j(n) c_j(n, t) D^(d_j)
/// P_(n - r_j)(t) [+ g(n, t)]``; ``P_first_index, ...`` are the initial
/// polynomials and the rule applies from the next index on.  A
/// ``Recurrence`` is an immutable wrapper of one record of the stable JSON
/// schema ``polytool.recurrence.v1`` (the output of ``polytool recurrence
/// --json``): build it with :func:`find_recurrence`, :meth:`from_json`, or
/// :func:`polytool.oeis.recurrence`, and convert it back with
/// :meth:`to_json` / :meth:`to_dict`.
///
/// A recurrence found by :func:`find_recurrence` is a *fit* to finitely many
/// rows (it carries ``search`` metadata): proofs about it concern the
/// sequence the recurrence defines, and applying them to an independently
/// defined sequence needs an all-index recurrence identity.
///
/// >>> r = polytool.Recurrence.from_json({
/// ...     "schema": "polytool.recurrence.v1", "first_index": 0,
/// ...     "initial_polynomials": [["1"]],
/// ...     "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
/// ...                               "coeff": {"coeffs": [["1", "1"]]}}]}})
/// >>> r.generate(4)
/// [Polynomial([1]), Polynomial([1, 1]), Polynomial([1, 2, 1]), Polynomial([1, 3, 3, 1])]
#[pyclass(module = "polytool", frozen)]
pub struct Recurrence {
    parsed: ParsedRecord,
}

impl Recurrence {
    pub(crate) fn from_record(record: RecurrenceJson) -> PyResult<Self> {
        Ok(Self {
            parsed: ParsedRecord::new(record)?,
        })
    }

    fn record(&self) -> &RecurrenceJson {
        &self.parsed.record
    }

    /// A `Recurrence`, or a record given as JSON text or a `dict`.
    fn record_from_python(value: &Bound<'_, PyAny>) -> PyResult<RecurrenceJson> {
        if let Ok(recurrence) = value.cast::<Recurrence>() {
            return Ok(recurrence.get().record().clone());
        }
        let text = json_text_from_python(value)?;
        serde_json::from_str(&text)
            .map_err(|error| polytool_error(format!("invalid recurrence JSON: {error}")))
    }

    fn rows_from(&self, from_index: Option<usize>) -> PyResult<RecurrenceRows> {
        let first_index = self.parsed.first_index;
        let from_index = from_index.unwrap_or(first_index);
        let skip = from_index.checked_sub(first_index).ok_or_else(|| {
            PyValueError::new_err(format!(
                "from_index {from_index} is before the first index {first_index}"
            ))
        })?;
        Ok(RecurrenceRows {
            rows: self
                .parsed
                .rule
                .generate_rows_rational_iter(&self.parsed.seeds, first_index),
            pending_skip: skip,
        })
    }
}

#[pymethods]
impl Recurrence {
    /// Read a ``polytool.recurrence.v1`` record from JSON text or a ``dict``.
    ///
    /// Raises :class:`PolytoolError` for an invalid record.
    #[staticmethod]
    fn from_json(source: &Bound<'_, PyAny>) -> PyResult<Self> {
        Self::from_record(Self::record_from_python(source)?)
    }

    /// The record as JSON text: compact, or indented by ``indent`` spaces
    /// (``indent=2`` reproduces ``polytool recurrence --json``).
    #[pyo3(signature = (indent = None))]
    fn to_json(&self, indent: Option<usize>) -> PyResult<String> {
        to_json_text(self.record(), indent)
    }

    /// The record as a ``dict`` (the JSON schema ``polytool.recurrence.v1``;
    /// exact coefficients are decimal strings such as ``"-3/2"``).
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        serialize_to_python(py, self.record())
    }

    /// Index ``n`` of the first initial polynomial ``P_n``.
    #[getter]
    fn first_index(&self) -> usize {
        self.parsed.first_index
    }

    /// The initial polynomials ``P_first_index, ...`` (list of
    /// :class:`Polynomial`).
    #[getter]
    fn initial_polynomials(&self) -> Vec<Polynomial> {
        self.parsed
            .seeds
            .iter()
            .map(|row| Polynomial::from_coefficients(row.clone()))
            .collect()
    }

    /// The largest lag ``r`` of the rule (its order).
    #[getter]
    fn order(&self) -> usize {
        self.parsed.rule.max_offset()
    }

    /// Search metadata of a fitted recurrence (a ``dict``), or ``None`` for a
    /// supplied one.
    #[getter]
    fn search<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.record()
            .search
            .as_ref()
            .map(|search| serialize_to_python(py, search))
            .transpose()
    }

    /// Whether the record comes from a search over finitely many rows.
    #[getter]
    fn is_fitted(&self) -> bool {
        self.record().search.is_some()
    }

    /// The rule as text, for example ``P(n) = (1 + t) P(n-1)``.
    #[getter]
    fn text(&self) -> String {
        self.parsed.rule.to_string()
    }

    /// The LaTeX presentation: one ``$...$`` formula per line, the initial
    /// polynomials first, then the rule with its range ``n >= ...``.
    ///
    /// Coefficients are factored over the rationals (verified by
    /// re-expansion) unless ``expanded=True``.  A fitted record starts with a
    /// ``%`` comment saying that validity for all ``n`` is unproved.
    #[pyo3(signature = (*, expanded = false))]
    fn latex(&self, py: Python<'_>, expanded: bool) -> PyResult<String> {
        let record = self.record();
        py.detach(|| {
            if expanded {
                record.to_latex_expanded()
            } else {
                record.to_latex()
            }
        })
        .map_err(polytool_error)
    }

    /// Notebook display: the formulas of :meth:`latex` in one display block.
    fn _repr_latex_(&self, py: Python<'_>) -> Option<String> {
        let latex = self.latex(py, false).ok()?;
        let lines: Vec<&str> = latex
            .lines()
            .filter(|line| !line.starts_with('%'))
            .map(|line| line.trim().trim_matches('$'))
            .collect();
        Some(format!(
            "$$\\begin{{gathered}}{}\\end{{gathered}}$$",
            lines.join(" \\\\ ")
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "<Recurrence {} with {} initial polynomial{} from index {}{}>",
            self.parsed.rule,
            self.parsed.seeds.len(),
            if self.parsed.seeds.len() == 1 {
                ""
            } else {
                "s"
            },
            self.parsed.first_index,
            if self.is_fitted() { ", fitted" } else { "" }
        )
    }

    fn __str__(&self) -> String {
        self.parsed.rule.to_string()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<Recurrence>()
            .is_ok_and(|other| other.get().record() == self.record())
    }

    fn __hash__(&self) -> PyResult<u64> {
        use std::hash::{Hash, Hasher};
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        to_json_text(self.record(), None)?.hash(&mut hasher);
        Ok(hasher.finish())
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let constructor = py.get_type::<Recurrence>().getattr("from_json")?;
        let arguments = PyTuple::new(py, [to_json_text(self.record(), None)?])?;
        PyTuple::new(py, [constructor, arguments.into_any()])
    }

    /// The polynomials ``P_n`` for ``n = from_index, ..., from_index + count
    /// - 1`` (list of :class:`Polynomial`); ``from_index`` defaults to
    /// :attr:`first_index`.
    ///
    /// Rows are computed exactly, one after the other; the rows before
    /// ``from_index`` are computed but not returned.  Raises
    /// :class:`PolytoolError` when a row cannot be evaluated (for example a
    /// left-hand factor that vanishes or does not divide).
    #[pyo3(signature = (count, *, from_index = None))]
    fn generate(
        &self,
        py: Python<'_>,
        count: usize,
        from_index: Option<usize>,
    ) -> PyResult<Vec<Polynomial>> {
        let mut rows = self.rows_from(from_index)?;
        py.detach(|| {
            rows.discard_pending()?;
            rows.rows
                .by_ref()
                .take(count)
                .map(|row| row.map(Polynomial::from_coefficients))
                .collect::<Result<Vec<_>, _>>()
                .map_err(|error| error.to_string())
        })
        .map_err(polytool_error)
    }

    /// A lazy, endless iterator over ``P_from_index, P_(from_index + 1),
    /// ...`` (:class:`Polynomial` items; ``from_index`` defaults to
    /// :attr:`first_index`).  It keeps only the rows the rule needs.
    #[pyo3(signature = (*, from_index = None))]
    fn iterate(&self, from_index: Option<usize>) -> PyResult<RecurrenceRows> {
        self.rows_from(from_index)
    }

    /// The same sequence with the first initial polynomial relabelled
    /// ``P_first_index``; the index ``n`` in the coefficients is shifted so
    /// that the rows do not change.
    fn reindexed(&self, first_index: usize) -> PyResult<Self> {
        Self::from_record(
            self.record()
                .reindexed(first_index)
                .map_err(polytool_error)?,
        )
    }

    /// Try the general real-rootedness theorems on this recurrence; see
    /// :func:`polytool.prove_real_rooted`.
    #[pyo3(signature = (*, maximum_start_delay = 3, root_region = None))]
    fn prove_real_rooted(
        &self,
        py: Python<'_>,
        maximum_start_delay: usize,
        root_region: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<ProofReport> {
        prove_record(py, self.record(), maximum_start_delay, root_region)
    }

    /// Certify the degree law and leading-coefficient signs of every row; see
    /// :func:`polytool.certify_leading_coefficient`.
    #[pyo3(signature = (*, sample_rows = 36, maximum_period = 6))]
    fn certify_leading_coefficient(
        &self,
        py: Python<'_>,
        sample_rows: usize,
        maximum_period: usize,
    ) -> PyResult<LeadingCoefficientCertificate> {
        certify_record(py, &self.parsed, sample_rows, maximum_period)
    }

    /// Conditional theorem matches of this recurrence (weighted renewal and
    /// fixed-delay predecessor windows), as the ``dict`` of schema
    /// ``polytool.recurrence-theorems.v1`` (``polytool recurrence-theorems
    /// --json``).  A match states paper hypotheses that hold; it is separate
    /// from the finite-fit evidence.
    fn match_theorems<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let report =
            RecurrenceTheoremReport::for_supplied(self.record()).map_err(polytool_error)?;
        serialize_to_python(py, &report)
    }
}

/// A lazy iterator over the rows of a :class:`Recurrence` (from
/// :meth:`Recurrence.iterate`).  Raises :class:`PolytoolError` at a row that
/// cannot be evaluated and stops after it.
#[pyclass(module = "polytool")]
pub struct RecurrenceRows {
    rows: RecurrenceRowIterator,
    pending_skip: usize,
}

impl RecurrenceRows {
    fn discard_pending(&mut self) -> Result<(), String> {
        while self.pending_skip > 0 {
            self.pending_skip -= 1;
            match self.rows.next() {
                Some(Ok(_)) => {}
                Some(Err(error)) => return Err(error.to_string()),
                None => return Err("the recurrence ended before the requested index".into()),
            }
        }
        Ok(())
    }
}

#[pymethods]
impl RecurrenceRows {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self, py: Python<'_>) -> PyResult<Polynomial> {
        let rows = self;
        py.detach(|| {
            rows.discard_pending()?;
            match rows.rows.next() {
                Some(Ok(row)) => Ok(Some(Polynomial::from_coefficients(row))),
                Some(Err(error)) => Err(error.to_string()),
                None => Ok(None),
            }
        })
        .map_err(polytool_error)?
        .ok_or_else(|| PyStopIteration::new_err(()))
    }
}

/// The `polytool.recurrence.v1` record of a search result, as
/// `polytool recurrence --json` prints it: rows after `skip_prefix` labelled
/// from `output_first_index`, with search metadata and the LaTeX
/// presentation.
fn search_record(
    result: &AdaptiveSearchResult,
    rows: &[Vec<BigRational>],
    search: &AdaptiveSearchOptions,
    output_first_index: usize,
    work_units: Option<u64>,
) -> Result<RecurrenceJson, RecurrenceError> {
    let searched = rows.get(search.skip_prefix..).unwrap_or(&[]);
    let initial_count = result.export_initial_count(searched.len());
    let mut record = RecurrenceJson::from_recurrence_rational(
        &result.recurrence,
        1,
        &searched[..initial_count],
        Some(RecurrenceJsonSearch {
            recurrence_text: result.recurrence.to_string(),
            source_rows: rows.len(),
            skip_prefix: search.skip_prefix,
            unknowns: result.num_unknowns,
            weighted_unknowns: result.weighted_unknowns,
            equations: result.num_equations,
            fit_polynomials: result.fit_polynomials,
            verification_polynomials: result.verification_polynomials,
            verification_equations: result.verification_equations,
            candidates_tried: result.candidates_tried,
            work_units: None,
            options: RecurrenceOptionsJson::from(&result.opts),
        }),
    )
    .reindexed(output_first_index)?;
    if let Some(search) = record.search.as_mut() {
        search.work_units = work_units;
    }
    record.with_latex()
}

/// Search for a linear recurrence satisfied by a sequence of polynomials.
///
/// ``rows`` are the polynomials ``P_first_index, P_(first_index + 1), ...``
/// (each a :class:`Polynomial`, a coefficient list or polynomial text).  The
/// adaptive search tries small parameter spaces first and returns the first
/// recurrence that fits the rows and reproduces the held-out rows exactly, as
/// a :class:`Recurrence` with ``search`` metadata; ``None`` when the
/// configured search space contains no recurrence.  The result is a fit to
/// finitely many rows, not a proof for every ``n``.
///
/// The keyword options mirror ``polytool recurrence`` (same defaults; in
/// particular the left-hand side is varied, ``laurent=True``):
///
/// - ``first_index`` (default 1): index of ``rows[0]``; with
///   ``skip_prefix`` the skipped rows still count;
/// - ``skip_prefix``: ignore the first rows;
/// - ``min_/max_recurrence_length`` (1, 5): number of previous rows used;
/// - ``min_/max_variable_degree`` (0, 3): degree in ``t`` of the coefficients;
/// - ``min_/max_index_degree`` (0, 3): degree in ``n`` of the coefficients;
/// - ``min_/max_derivative_order`` (0, 2);
/// - ``inhomogeneous`` with ``min_/max_inhomogeneous_variable_degree`` and
///   ``min_/max_inhomogeneous_index_degree``: allow an additive term
///   (setting a bound turns it on);
/// - ``denominator`` with ``max_denominator_variable_degree`` and
///   ``max_denominator_index_degree``: allow a left-hand factor;
/// - ``alternating_sign``: also try terms with a factor ``(-1)^n``;
/// - ``laurent`` and ``max_laurent_power``: vary the left-hand side
///   (``t^e f(n, t) P_n``, ``e <= max_laurent_power``);
/// - ``min_margin``, ``fit_extra_rows``, ``verify`` (held-out rows; default
///   ``True``) and ``modular_prefilter`` tune the fitting;
/// - ``max_candidates`` and ``max_work``: deterministic budgets.  When one
///   is reached, :class:`BudgetExhaustedError` is raised (attributes
///   ``candidates_tried``, ``candidates_considered``, ``work_units``).
///
/// The search releases the GIL; it cannot be interrupted with Ctrl-C, so use
/// ``max_work`` for a bounded run.
///
/// >>> rows = [[1], [1, 1], [1, 4, 1], [1, 11, 11, 1], [1, 26, 66, 26, 1],
/// ...         [1, 57, 302, 302, 57, 1], [1, 120, 1191, 2416, 1191, 120, 1]]
/// >>> r = polytool.find_recurrence(rows)
/// >>> r.text
/// "P(n) = (1 - t + nt) P(n-1) + (t - t^2) P'(n-1)"
/// >>> r.prove_real_rooted().route
/// 'liu_wang'
#[pyfunction]
#[pyo3(signature = (
    rows,
    *,
    first_index = None,
    skip_prefix = 0,
    min_recurrence_length = None,
    max_recurrence_length = None,
    min_variable_degree = None,
    max_variable_degree = None,
    min_index_degree = None,
    max_index_degree = None,
    min_derivative_order = None,
    max_derivative_order = None,
    inhomogeneous = false,
    min_inhomogeneous_variable_degree = None,
    max_inhomogeneous_variable_degree = None,
    min_inhomogeneous_index_degree = None,
    max_inhomogeneous_index_degree = None,
    denominator = false,
    max_denominator_variable_degree = None,
    max_denominator_index_degree = None,
    alternating_sign = false,
    laurent = true,
    max_laurent_power = None,
    min_margin = None,
    fit_extra_rows = None,
    verify = true,
    modular_prefilter = None,
    max_candidates = None,
    max_work = None,
))]
#[allow(clippy::too_many_arguments)]
pub(crate) fn find_recurrence(
    py: Python<'_>,
    rows: &Bound<'_, PyAny>,
    first_index: Option<usize>,
    skip_prefix: usize,
    min_recurrence_length: Option<usize>,
    max_recurrence_length: Option<usize>,
    min_variable_degree: Option<usize>,
    max_variable_degree: Option<usize>,
    min_index_degree: Option<usize>,
    max_index_degree: Option<usize>,
    min_derivative_order: Option<usize>,
    max_derivative_order: Option<usize>,
    inhomogeneous: bool,
    min_inhomogeneous_variable_degree: Option<usize>,
    max_inhomogeneous_variable_degree: Option<usize>,
    min_inhomogeneous_index_degree: Option<usize>,
    max_inhomogeneous_index_degree: Option<usize>,
    denominator: bool,
    max_denominator_variable_degree: Option<usize>,
    max_denominator_index_degree: Option<usize>,
    alternating_sign: bool,
    laurent: bool,
    max_laurent_power: Option<usize>,
    min_margin: Option<usize>,
    fit_extra_rows: Option<usize>,
    verify: bool,
    modular_prefilter: Option<bool>,
    max_candidates: Option<usize>,
    max_work: Option<u64>,
) -> PyResult<Option<Recurrence>> {
    let rows = polynomial_rows_from_python(rows)?;
    let mut search = AdaptiveSearchOptions::default().with_laurent(laurent);
    search.skip_prefix = skip_prefix;
    let set = |target: &mut usize, value: Option<usize>| {
        if let Some(value) = value {
            *target = value;
        }
    };
    set(&mut search.min_rec_len, min_recurrence_length);
    set(&mut search.max_rec_len, max_recurrence_length);
    set(&mut search.min_var_deg, min_variable_degree);
    set(&mut search.max_var_deg, max_variable_degree);
    set(&mut search.min_idx_deg, min_index_degree);
    set(&mut search.max_idx_deg, max_index_degree);
    set(&mut search.min_diff_deg, min_derivative_order);
    set(&mut search.max_diff_deg, max_derivative_order);
    let inhomogeneous_bounds = [
        min_inhomogeneous_variable_degree,
        max_inhomogeneous_variable_degree,
        min_inhomogeneous_index_degree,
        max_inhomogeneous_index_degree,
    ];
    search.try_inhomogeneous = inhomogeneous || inhomogeneous_bounds.iter().any(Option::is_some);
    set(
        &mut search.min_inhomo_var_deg,
        min_inhomogeneous_variable_degree,
    );
    set(
        &mut search.max_inhomo_var_deg,
        max_inhomogeneous_variable_degree,
    );
    set(
        &mut search.min_inhomo_idx_deg,
        min_inhomogeneous_index_degree,
    );
    set(
        &mut search.max_inhomo_idx_deg,
        max_inhomogeneous_index_degree,
    );
    search.try_denominator = denominator
        || max_denominator_variable_degree.is_some()
        || max_denominator_index_degree.is_some();
    set(
        &mut search.max_denom_var_deg,
        max_denominator_variable_degree,
    );
    set(&mut search.max_denom_idx_deg, max_denominator_index_degree);
    search.try_alternating_sign = alternating_sign;
    set(&mut search.max_laurent_x_power, max_laurent_power);
    set(&mut search.min_margin, min_margin);
    set(&mut search.fit_extra_rows, fit_extra_rows);
    search.no_verify = !verify;
    if let Some(prefilter) = modular_prefilter {
        search.modular_prefilter = prefilter;
    }
    let budget = AdaptiveSearchBudget {
        max_candidates,
        max_work,
    };
    let output_first_index = first_index
        .unwrap_or(1)
        .checked_add(if first_index.is_some() {
            skip_prefix
        } else {
            0
        })
        .filter(|index| *index <= i64::MAX as usize && index.checked_add(rows.len()).is_some())
        .ok_or_else(|| PyValueError::new_err("first index or final row index is too large"))?;
    let outcome =
        py.detach(|| find_recurrence_adaptive_rational_with_budget(&rows, &search, budget));
    match outcome {
        AdaptiveSearchOutcome::Found(result) => {
            let work_units = max_work.map(|_| result.diagnostics.work_units);
            let record = py
                .detach(|| search_record(&result, &rows, &search, output_first_index, work_units))
                .map_err(polytool_error)?;
            Ok(Some(Recurrence::from_record(record)?))
        }
        AdaptiveSearchOutcome::NoRecurrence(_) => Ok(None),
        AdaptiveSearchOutcome::BudgetExhausted(summary) => {
            let error = BudgetExhaustedError::new_err(format!(
                "recurrence search budget exhausted after {} candidates ({} work units)",
                summary.diagnostics.considered_candidates, summary.diagnostics.work_units
            ));
            let value = error.value(py);
            value.setattr("candidates_tried", summary.candidates_tried)?;
            value.setattr(
                "candidates_considered",
                summary.diagnostics.considered_candidates,
            )?;
            value.setattr("work_units", summary.diagnostics.work_units)?;
            Err(error)
        }
    }
}

/// An open root region from `(lower, upper)` with `None` for an infinite end.
fn root_region_from_python(value: &Bound<'_, PyAny>) -> PyResult<RealInterval> {
    let (lower, upper): (Bound<'_, PyAny>, Bound<'_, PyAny>) = value.extract().map_err(|_| {
        PyValueError::new_err("root_region must be a pair (lower, upper); None is infinite")
    })?;
    let end = |value: &Bound<'_, PyAny>| -> PyResult<Option<BigRational>> {
        if value.is_none() {
            Ok(None)
        } else {
            rational_from_python(value).map(Some)
        }
    };
    Ok(RealInterval {
        lower: end(&lower)?,
        lower_closed: false,
        upper: end(&upper)?,
        upper_closed: false,
    })
}

fn prove_record(
    py: Python<'_>,
    record: &RecurrenceJson,
    maximum_start_delay: usize,
    root_region: Option<&Bound<'_, PyAny>>,
) -> PyResult<ProofReport> {
    let options = RealRootednessTheoremOptions {
        root_region: root_region.map(root_region_from_python).transpose()?,
        maximum_start_delay,
    };
    let report = py
        .detach(|| discover_recurrence_proofs(record, &options))
        .map_err(polytool_error)?;
    Ok(ProofReport { report })
}

fn certify_record(
    py: Python<'_>,
    parsed: &ParsedRecord,
    sample_rows: usize,
    maximum_period: usize,
) -> PyResult<LeadingCoefficientCertificate> {
    if maximum_period == 0 {
        return Err(PyValueError::new_err("maximum_period must be at least 1"));
    }
    let options = LeadingCoefficientOptions {
        sample_rows,
        maximum_period,
    };
    let certificate = py
        .detach(|| {
            library_certify_leading_coefficient(
                &parsed.rule,
                parsed.first_index,
                &parsed.seeds,
                &options,
            )
        })
        .map_err(polytool_error)?;
    Ok(LeadingCoefficientCertificate { certificate })
}

/// Try the general real-rootedness theorems on a recurrence.
///
/// ``recurrence`` is a :class:`Recurrence`, or a ``polytool.recurrence.v1``
/// record as JSON text or ``dict``.  The search tries, in order, product
/// induction, the three-term (Favard) theorem, Liu--Wang first order,
/// Liu--Wang derivative lag, the Ma--Wang interval criterion, the
/// Liu--Wang three-term theorem, the iterated-preserver and the
/// hypergeometric-row criteria, and checks every hypothesis exactly for all
/// ``n`` (``polytool recurrence-prove``); then the tail after rows that are
/// zero or not real-rooted, and alternative representations with an
/// all-index equality certificate.  When a theorem fails at the
/// supplied initial rows, its induction may start up to
/// ``maximum_start_delay`` rows later (the skipped rows are checked
/// exactly).  ``root_region=(lower, upper)`` (``None`` for an infinite end)
/// asks for every root in that open interval.
///
/// Returns a :class:`ProofReport`.  ``status == "not_proved"`` only means
/// that no sufficient criterion applied; it is not a disproof.  A proof
/// concerns the sequence the recurrence defines.
#[pyfunction]
#[pyo3(signature = (recurrence, *, maximum_start_delay = 3, root_region = None))]
pub(crate) fn prove_real_rooted(
    py: Python<'_>,
    recurrence: &Bound<'_, PyAny>,
    maximum_start_delay: usize,
    root_region: Option<&Bound<'_, PyAny>>,
) -> PyResult<ProofReport> {
    let record = Recurrence::record_from_python(recurrence)?;
    prove_record(py, &record, maximum_start_delay, root_region)
}

/// Certify the degree law ``deg P_n = D(n)`` and the sign of every leading
/// coefficient of a recurrence family, for all ``n``.
///
/// ``recurrence`` is a :class:`Recurrence` or a record (JSON text or
/// ``dict``).  ``D(n)`` is periodic by residue class (period at most
/// ``maximum_period``); the law is read off the first ``sample_rows`` exact
/// rows and every later index is covered by one of the routes of the
/// RealRooted tactic ``rr_row_leadingCoeff_pos`` (``cone``, ``recurrence``,
/// ``subseq``, ``twoStep``, ``cancel``, ``coefficients``), as
/// ``polytool recurrence-degree``.  Returns a
/// :class:`LeadingCoefficientCertificate`; a certificate that is not proved
/// is not a disproof.
#[pyfunction]
#[pyo3(signature = (recurrence, *, sample_rows = 36, maximum_period = 6))]
pub(crate) fn certify_leading_coefficient(
    py: Python<'_>,
    recurrence: &Bound<'_, PyAny>,
    sample_rows: usize,
    maximum_period: usize,
) -> PyResult<LeadingCoefficientCertificate> {
    let parsed = ParsedRecord::new(Recurrence::record_from_python(recurrence)?)?;
    certify_record(py, &parsed, sample_rows, maximum_period)
}

/// The result of :func:`prove_real_rooted` (schema
/// ``polytool.recurrence-proof-discovery.v1``, the output of
/// ``polytool recurrence-prove --json``).
///
/// ``status`` is ``"proved"`` when at least one theorem checked every
/// hypothesis, else ``"not_proved"`` or ``"not_applicable"`` (no disproof).  ``attempts`` lists every
/// theorem tried with its hypotheses; ``route`` and ``certificate`` name the
/// first proved one.  :meth:`verify` replays the whole search exactly.
#[pyclass(module = "polytool", frozen)]
pub struct ProofReport {
    report: RecurrenceProofDiscoveryReport,
}

impl ProofReport {
    /// The first proved attempt: on the supplied record, else on its tail,
    /// else on a certified alternative representation.
    fn first_proved_attempt(&self) -> Option<&RealRootednessTheoremReport> {
        let direct = self
            .report
            .attempts
            .iter()
            .find(|attempt| attempt.status == RealRootednessTheoremStatus::Proved);
        let tail = || {
            self.report.tail.as_ref().and_then(|tail| {
                tail.attempts
                    .iter()
                    .find(|attempt| attempt.status == RealRootednessTheoremStatus::Proved)
            })
        };
        let alternative = || {
            self.report
                .alternative_representations
                .as_ref()
                .and_then(|search| search.alternatives.iter().find(|a| a.proved))
                .and_then(|alternative| {
                    alternative
                        .attempts
                        .iter()
                        .find(|attempt| attempt.status == RealRootednessTheoremStatus::Proved)
                })
        };
        direct.or_else(tail).or_else(alternative)
    }
}

#[pymethods]
impl ProofReport {
    /// Read a report from JSON text or a ``dict``.
    #[staticmethod]
    fn from_json(source: &Bound<'_, PyAny>) -> PyResult<Self> {
        let text = json_text_from_python(source)?;
        let report = serde_json::from_str(&text)
            .map_err(|error| polytool_error(format!("invalid proof report JSON: {error}")))?;
        Ok(Self { report })
    }

    /// The report as JSON text (compact, or indented by ``indent`` spaces).
    #[pyo3(signature = (indent = None))]
    fn to_json(&self, indent: Option<usize>) -> PyResult<String> {
        to_json_text(&self.report, indent)
    }

    /// The report as a ``dict`` (its JSON schema).
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        serialize_to_python(py, &self.report)
    }

    /// ``"proved"``, ``"not_proved"`` or ``"not_applicable"`` (a rule no
    /// theorem matches, such as an inhomogeneous rule).
    #[getter]
    fn status(&self) -> &str {
        &self.report.status
    }

    /// Whether some theorem proved the conclusion (on the record, its tail
    /// or a certified alternative representation).
    #[getter]
    fn is_proved(&self) -> bool {
        self.first_proved_attempt().is_some()
    }

    /// The first index ``N`` with every row ``P_n``, ``n >= N``, proved
    /// nonzero and real-rooted, or ``None``.
    #[getter]
    fn real_rooted_from_index(&self) -> Option<usize> {
        self.report.real_rooted_from_index
    }

    /// The theorem of the first proved attempt (for example
    /// ``"three_term_recurrence"``), or ``None``; attempts on the record come
    /// first, then those on its tail and on a certified alternative.
    #[getter]
    fn route(&self) -> PyResult<Option<String>> {
        self.first_proved_attempt()
            .map(|attempt| {
                serde_json::to_value(attempt.theorem)
                    .map_err(polytool_error)
                    .map(|value| value.as_str().unwrap_or_default().to_owned())
            })
            .transpose()
    }

    /// The first proved attempt (a ``dict`` of schema
    /// ``polytool.recurrence-real-rootedness-theorem.v1`` with the
    /// hypotheses and the conclusion), or ``None``.
    #[getter]
    fn certificate<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.first_proved_attempt()
            .map(|attempt| serialize_to_python(py, attempt))
            .transpose()
    }

    /// Every theorem attempt, proved or not (list of ``dict``).
    #[getter]
    fn attempts<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        serialize_to_python(py, &self.report.attempts)
    }

    /// Whether the supplied (not fitted) rule has the product form
    /// ``B(n) P_n = A(x, n) P_(n-k)`` and product induction proves it.
    #[getter]
    fn trivially_real_rooted(&self) -> bool {
        self.report.trivially_real_rooted
    }

    /// What a successful search establishes, in words.
    #[getter]
    fn statement(&self) -> &str {
        &self.report.statement
    }

    /// The recurrence the report is about.
    #[getter]
    fn recurrence(&self) -> PyResult<Recurrence> {
        Recurrence::from_record(self.report.recurrence.clone())
    }

    /// Replay the complete search exactly and compare it with the report.
    ///
    /// Returns ``True``; raises :class:`VerificationError` on any mismatch
    /// (a changed conclusion, a dropped attempt, another rule).  Verifying a
    /// ``not_proved`` report is not a proof.
    fn verify(&self, py: Python<'_>) -> PyResult<bool> {
        py.detach(|| self.report.verify())
            .map_err(|error| VerificationError::new_err(error.to_string()))?;
        Ok(true)
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(match self.route()? {
            Some(route) => format!("<ProofReport proved by {route}>"),
            None => format!(
                "<ProofReport {} after {} attempts>",
                self.report.status,
                self.report.attempts.len()
            ),
        })
    }
}

/// The result of :func:`certify_leading_coefficient` (schema
/// ``polytool.recurrence-leading-coefficient.v1``, the output of
/// ``polytool recurrence-degree --json``).
///
/// When proved, every ``P_n`` with ``n >= holds_from`` is nonzero with the
/// degree and leading-coefficient sign of :attr:`degree_law` (exact degree
/// unless ``degree_is_exact`` is false, then only an upper bound).
#[pyclass(module = "polytool", frozen)]
pub struct LeadingCoefficientCertificate {
    certificate: LibraryLeadingCoefficientCertificate,
}

#[pymethods]
impl LeadingCoefficientCertificate {
    /// Read a certificate from JSON text or a ``dict``.
    #[staticmethod]
    fn from_json(source: &Bound<'_, PyAny>) -> PyResult<Self> {
        let text = json_text_from_python(source)?;
        let certificate = serde_json::from_str(&text).map_err(|error| {
            polytool_error(format!(
                "invalid leading-coefficient certificate JSON: {error}"
            ))
        })?;
        Ok(Self { certificate })
    }

    /// The certificate as JSON text (compact, or indented by ``indent``).
    #[pyo3(signature = (indent = None))]
    fn to_json(&self, indent: Option<usize>) -> PyResult<String> {
        to_json_text(&self.certificate, indent)
    }

    /// The certificate as a ``dict`` (its JSON schema).
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        serialize_to_python(py, &self.certificate)
    }

    /// ``"proved"``, ``"not_proved"`` or ``"not_applicable"``.
    #[getter]
    fn status(&self) -> PyResult<String> {
        serde_json::to_value(self.certificate.status)
            .map_err(polytool_error)
            .map(|value| value.as_str().unwrap_or_default().to_owned())
    }

    /// Whether the degree law is proved for every index from ``holds_from``.
    #[getter]
    fn is_proved(&self) -> bool {
        self.certificate.is_proved()
    }

    /// The route (``"cone"``, ``"recurrence"``, ``"subseq"``, ``"twoStep"``,
    /// ``"cancel"`` or ``"coefficients"``), or ``None``.
    #[getter]
    fn route(&self) -> Option<&'static str> {
        self.certificate.route.map(|route| route.lean_route_name())
    }

    /// The degree law as a ``dict`` (``start``, ``period``,
    /// ``base_degrees``, ``growth_per_period``,
    /// ``leading_coefficient_signs``), or ``None``.
    #[getter]
    fn degree_law<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.certificate
            .degree_law
            .as_ref()
            .map(|law| serialize_to_python(py, law))
            .transpose()
    }

    /// ``True`` when the law gives ``deg P_n`` exactly, ``False`` when it is
    /// only an upper bound (route ``coefficients``).
    #[getter]
    fn degree_is_exact(&self) -> bool {
        self.certificate.degree_is_exact
    }

    /// The first index from which the law holds, or ``None``.
    #[getter]
    fn holds_from(&self) -> Option<usize> {
        self.certificate.holds_from
    }

    /// The conclusion in words (present exactly when proved).
    #[getter]
    fn conclusion(&self) -> Option<&str> {
        self.certificate.conclusion.as_deref()
    }

    /// The constant leading coefficient of every row from the induction
    /// window on (for example ``1`` for monic rows), or ``None``.
    #[getter]
    fn exact_leading_coefficient<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.certificate
            .exact_leading_coefficient
            .as_ref()
            .map(|value| rational_to_python(py, value))
            .transpose()
    }

    /// ``deg P_n`` for a proved exact law and ``n >= holds_from``, else
    /// ``None``.
    fn degree_at(&self, n: usize) -> Option<usize> {
        self.certificate.degree_at(n)
    }

    /// The sign of the leading coefficient of ``P_n`` (``"positive"`` or
    /// ``"negative"``) for a proved certificate and ``n >= holds_from``,
    /// else ``None``.
    fn sign_at(&self, n: usize) -> Option<&'static str> {
        self.certificate.sign_at(n).map(|sign| match sign {
            LeadingCoefficientSign::Positive => "positive",
            LeadingCoefficientSign::Negative => "negative",
        })
    }

    /// The recurrence the certificate is about.
    #[getter]
    fn recurrence(&self) -> PyResult<Recurrence> {
        Recurrence::from_record(self.certificate.recurrence.clone())
    }

    /// Replay the recorded route and every sign certificate exactly.
    ///
    /// Returns ``True``; raises :class:`VerificationError` on a mismatch.
    fn verify(&self, py: Python<'_>) -> PyResult<bool> {
        py.detach(|| self.certificate.verify())
            .map_err(|error| VerificationError::new_err(error.to_string()))?;
        Ok(true)
    }

    fn __repr__(&self) -> PyResult<String> {
        Ok(match &self.certificate.conclusion {
            Some(conclusion) => format!("<LeadingCoefficientCertificate proved: {conclusion}>"),
            None => format!("<LeadingCoefficientCertificate {}>", self.status()?),
        })
    }
}
