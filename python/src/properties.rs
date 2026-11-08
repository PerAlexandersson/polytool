//! Real roots, interlacing and coefficient properties of one polynomial.
//!
//! Every function here takes exact coefficients in ascending degree order
//! (a `polytool.Polynomial`, a list of `int`/`Fraction`, or polynomial text)
//! and calls the exact library routine.  Rational input is decided on its
//! smallest positive integer multiple (`convert::integer_multiple`); the
//! predicates bound here do not change under a positive scale.

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use polytool::real_rootedness::{
    isolate_real_roots_rational_coeffs, RealRootEnclosure as LibraryEnclosure,
};
use polytool::BigRational;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyDict;

use crate::convert::{
    coefficients_from_python, integer_coefficients, integer_multiple, rational_from_python,
    rational_to_python, rationals_to_python,
};

/// One distinct real root of a polynomial, enclosed exactly, with its
/// multiplicity.
///
/// When ``lower == upper`` the root is exactly this rational number
/// (``is_exact``).  Otherwise the root is the only root in the open interval
/// ``(lower, upper)`` and neither endpoint is a root.  The endpoints are
/// ``int`` or ``fractions.Fraction``; enclosures of different roots never
/// overlap and are listed from left to right.
#[pyclass(module = "polytool", frozen, eq)]
#[derive(Clone, PartialEq)]
pub struct RealRootEnclosure {
    enclosure: LibraryEnclosure,
}

#[pymethods]
impl RealRootEnclosure {
    /// Lower endpoint (``int`` or ``Fraction``).
    #[getter]
    fn lower<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        rational_to_python(py, &self.enclosure.lower)
    }

    /// Upper endpoint (``int`` or ``Fraction``).
    #[getter]
    fn upper<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        rational_to_python(py, &self.enclosure.upper)
    }

    /// Multiplicity of the root (at least 1).
    #[getter]
    fn multiplicity(&self) -> usize {
        self.enclosure.multiplicity
    }

    /// Whether the root is the rational number ``lower == upper``.
    #[getter]
    fn is_exact(&self) -> bool {
        self.enclosure.is_exact()
    }

    fn __repr__(&self) -> String {
        format!(
            "RealRootEnclosure(lower={}, upper={}, multiplicity={})",
            python_number_text(&self.enclosure.lower),
            python_number_text(&self.enclosure.upper),
            self.enclosure.multiplicity
        )
    }
}

/// `repr` of an exact number as Python prints it (`3`, `Fraction(1, 2)`).
pub(crate) fn python_number_text(value: &BigRational) -> String {
    if value.is_integer() {
        value.numer().to_string()
    } else {
        format!("Fraction({}, {})", value.numer(), value.denom())
    }
}

pub(crate) fn is_real_rooted_coefficients(coefficients: &[BigRational]) -> bool {
    polytool::is_real_rooted_bigint_coeffs(&integer_multiple(coefficients).0)
}

pub(crate) fn count_real_roots_coefficients(coefficients: &[BigRational], distinct: bool) -> usize {
    let isolation = isolate_real_roots_rational_coeffs(coefficients);
    if distinct {
        isolation.roots.len()
    } else {
        isolation.real_root_count
    }
}

pub(crate) fn isolate_real_roots_coefficients(
    coefficients: &[BigRational],
    max_width: Option<BigRational>,
) -> PyResult<Vec<RealRootEnclosure>> {
    let mut isolation = isolate_real_roots_rational_coeffs(coefficients);
    if let Some(width) = max_width {
        if !width.is_positive() {
            return Err(PyValueError::new_err("max_width must be positive"));
        }
        isolation.refine(&width);
    }
    Ok(isolation
        .roots
        .into_iter()
        .map(|enclosure| RealRootEnclosure { enclosure })
        .collect())
}

pub(crate) fn interlaces_coefficients(
    f: &[BigRational],
    g: &[BigRational],
    strict: bool,
) -> Option<bool> {
    let (f, _) = integer_multiple(f);
    let (g, _) = integer_multiple(g);
    if strict {
        polytool::check_interlacing_bigint_coeffs(&f, &g)
    } else {
        polytool::check_weak_interlacing_bigint_coeffs(&f, &g)
    }
}

pub(crate) fn gamma_coefficients_of(
    coefficients: &[BigRational],
    ignore_initial_zeros: bool,
) -> Option<Vec<BigRational>> {
    let (integers, scale) = integer_multiple(coefficients);
    let gamma = if ignore_initial_zeros {
        polytool::gamma_coefficients_ignoring_initial_zeros_bigint_coeffs(&integers)
    } else {
        polytool::gamma_coefficients_bigint_coeffs(&integers)
    }?;
    // gamma(c p) = c gamma(p): undo the positive scale exactly.
    Some(
        gamma
            .into_iter()
            .map(|g| BigRational::new(g, scale.clone()))
            .collect(),
    )
}

/// Which coefficient predicate to apply (they share the scaling step).
#[derive(Clone, Copy)]
pub(crate) enum CoefficientProperty {
    LogConcave,
    UltraLogConcave,
    Unimodal,
    Palindromic,
    PalindromicIgnoringInitialZeros,
    GammaPositive,
    GammaPositiveIgnoringInitialZeros,
    SimpleRoots,
}

pub(crate) fn has_property(coefficients: &[BigRational], property: CoefficientProperty) -> bool {
    use CoefficientProperty::*;
    let (c, _) = integer_multiple(coefficients);
    match property {
        LogConcave => polytool::is_log_concave_bigint_coeffs(&c),
        UltraLogConcave => polytool::is_ultra_log_concave_bigint_coeffs(&c),
        Unimodal => polytool::is_unimodal_bigint_coeffs(&c),
        Palindromic => polytool::is_palindromic_bigint_coeffs(&c),
        PalindromicIgnoringInitialZeros => {
            polytool::is_palindromic_ignoring_initial_zeros_bigint_coeffs(&c)
        }
        GammaPositive => polytool::is_gamma_positive_bigint_coeffs(&c),
        GammaPositiveIgnoringInitialZeros => {
            polytool::is_gamma_positive_ignoring_initial_zeros_bigint_coeffs(&c)
        }
        SimpleRoots => polytool::has_simple_roots_bigint_coeffs(&c),
    }
}

/// The text form of a polynomial: the CLI's format for integer coefficients
/// (`1 + 11t + 11t^2 + t^3`); non-integer coefficients are parenthesized
/// (`1/2 - (3/4)t^2`).
pub(crate) fn polynomial_text(coefficients: &[BigRational]) -> String {
    if let Ok(integers) = integer_coefficients(coefficients, "") {
        return polytool::format_poly_bigint_coeffs(&integers);
    }
    let mut text = String::new();
    for (power, coefficient) in coefficients.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        let magnitude = coefficient.abs();
        let monomial = match power {
            0 => String::new(),
            1 => "t".to_owned(),
            _ => format!("t^{power}"),
        };
        let factor = if power > 0 && magnitude.is_one() {
            String::new()
        } else if power > 0 && !magnitude.is_integer() {
            format!("({magnitude})")
        } else {
            magnitude.to_string()
        };
        let sign = coefficient.is_negative();
        if text.is_empty() {
            text.push_str(if sign { "-" } else { "" });
        } else {
            text.push_str(if sign { " - " } else { " + " });
        }
        text.push_str(&factor);
        text.push_str(&monomial);
    }
    text
}

/// The `polytool properties --json` item of one polynomial, without `index`.
pub(crate) fn properties_dict<'py>(
    py: Python<'py>,
    coefficients: &[BigRational],
) -> PyResult<Bound<'py, PyDict>> {
    use CoefficientProperty::*;
    let gamma = gamma_coefficients_of(coefficients, true);
    let gamma_positive = gamma
        .as_ref()
        .is_some_and(|gamma| gamma.iter().all(|g| !g.is_negative()));
    let dict = PyDict::new(py);
    dict.set_item("polynomial", polynomial_text(coefficients))?;
    dict.set_item("coefficients", rationals_to_python(py, coefficients)?)?;
    dict.set_item("degree", coefficients.len().saturating_sub(1))?;
    dict.set_item("real_rooted", is_real_rooted_coefficients(coefficients))?;
    dict.set_item("simple_roots", has_property(coefficients, SimpleRoots))?;
    dict.set_item(
        "palindromic",
        has_property(coefficients, PalindromicIgnoringInitialZeros),
    )?;
    dict.set_item("gamma_positive", gamma_positive)?;
    match &gamma {
        Some(gamma) => dict.set_item("gamma_coefficients", rationals_to_python(py, gamma)?)?,
        None => dict.set_item("gamma_coefficients", py.None())?,
    }
    dict.set_item("unimodal", has_property(coefficients, Unimodal))?;
    dict.set_item("log_concave", has_property(coefficients, LogConcave))?;
    dict.set_item(
        "ultra_log_concave",
        has_property(coefficients, UltraLogConcave),
    )?;
    Ok(dict)
}

/// Whether every root of ``polynomial`` is real.
///
/// Exact: decided by integer root counting (the library's
/// ``is_real_rooted_bigint_coeffs``); rational coefficients are scaled to
/// integers first.  The zero polynomial and nonzero constants count as
/// real-rooted.
///
/// ``polynomial`` is a :class:`Polynomial`, a list of ``int``/``Fraction``
/// coefficients in ascending degree order, or polynomial text with integer
/// coefficients.
///
/// >>> polytool.is_real_rooted([1, 11, 11, 1])
/// True
#[pyfunction]
pub(crate) fn is_real_rooted(polynomial: &Bound<'_, PyAny>) -> PyResult<bool> {
    let coefficients = coefficients_from_python(polynomial)?;
    Ok(polynomial
        .py()
        .detach(|| is_real_rooted_coefficients(&coefficients)))
}

/// Number of real roots of ``polynomial``, exactly.
///
/// By default roots are counted with multiplicity, so a nonzero polynomial is
/// real-rooted exactly when the count equals its degree.  With
/// ``distinct=True`` each distinct root counts once.  The zero polynomial has
/// no roots (count 0).
///
/// >>> polytool.count_real_roots([1, 2, 1])            # (1 + t)^2
/// 2
/// >>> polytool.count_real_roots([1, 2, 1], distinct=True)
/// 1
#[pyfunction]
#[pyo3(signature = (polynomial, *, distinct = false))]
pub(crate) fn count_real_roots(polynomial: &Bound<'_, PyAny>, distinct: bool) -> PyResult<usize> {
    let coefficients = coefficients_from_python(polynomial)?;
    Ok(polynomial
        .py()
        .detach(|| count_real_roots_coefficients(&coefficients, distinct)))
}

/// Isolate the distinct real roots of ``polynomial`` by exact rational
/// intervals.
///
/// Returns a list of :class:`RealRootEnclosure` from left to right, each
/// with its exact multiplicity.  Every enclosure has width at most ``1/1024``,
/// or at most ``max_width`` (a positive ``int``/``Fraction``) when given.
/// Integer roots and other simple rationals are usually reported exactly
/// (``lower == upper``).  The zero polynomial gives an empty list.
///
/// >>> roots = polytool.isolate_real_roots([0, -2, -4, -1, 2, 1])  # t (t + 1)^2 (t^2 - 2)
/// >>> [(r.lower, r.upper, r.multiplicity) for r in roots[1:3]]
/// [(-1, -1, 2), (0, 0, 1)]
/// >>> roots[3].lower < Fraction(14142, 10000) < roots[3].upper
/// True
#[pyfunction]
#[pyo3(signature = (polynomial, *, max_width = None))]
pub(crate) fn isolate_real_roots(
    polynomial: &Bound<'_, PyAny>,
    max_width: Option<&Bound<'_, PyAny>>,
) -> PyResult<Vec<RealRootEnclosure>> {
    let coefficients = coefficients_from_python(polynomial)?;
    let max_width = max_width.map(rational_from_python).transpose()?;
    polynomial
        .py()
        .detach(|| isolate_real_roots_coefficients(&coefficients, max_width))
}

/// Whether ``f`` interlaces ``g`` (written ``f << g``), exactly.
///
/// With roots ``a_1 >= a_2 >= ...`` of ``f`` and ``b_1 >= b_2 >= ...`` of
/// ``g`` (both real-rooted):
///
/// - ``deg g = deg f + 1``: ``b_1 >= a_1 >= b_2 >= a_2 >= ... >= a_m >= b_(m+1)``;
/// - ``deg g = deg f``: ``b_1 >= a_1 >= b_2 >= a_2 >= ... >= b_m >= a_m``.
///
/// The default test is *weak* interlacing (shared roots allowed).  With
/// ``strict=True`` the inequalities are strict and ``deg g = deg f + 1`` is
/// required (Bezout-matrix criterion).  The answer depends on the roots only,
/// not on constant factors.  The zero polynomial interlaces every real-rooted
/// polynomial.
///
/// Returns ``True`` or ``False`` when the degrees admit the directed test and
/// ``None`` when they do not (for example ``deg f = deg g + 1``: try
/// ``interlaces(g, f)`` or :func:`interlaces_either_order`).  ``None`` does
/// not mean "does not interlace".
///
/// >>> polytool.interlaces([8, -6, 1], [-15, 23, -9, 1])   # (t-2)(t-4) << (t-1)(t-3)(t-5)
/// True
/// >>> polytool.interlaces([-15, 23, -9, 1], [8, -6, 1]) is None
/// True
#[pyfunction]
#[pyo3(signature = (f, g, *, strict = false))]
pub(crate) fn interlaces(
    f: &Bound<'_, PyAny>,
    g: &Bound<'_, PyAny>,
    strict: bool,
) -> PyResult<Option<bool>> {
    let py = f.py();
    let f = coefficients_from_python(f)?;
    let g = coefficients_from_python(g)?;
    Ok(py.detach(|| interlaces_coefficients(&f, &g, strict)))
}

/// Whether ``f`` and ``g`` interlace in some order allowed by their degrees.
///
/// The undirected question: ``True`` when ``f << g`` or ``g << f`` (weak
/// interlacing).  The zero polynomial interlaces every real-rooted
/// polynomial.
#[pyfunction]
pub(crate) fn interlaces_either_order(
    f: &Bound<'_, PyAny>,
    g: &Bound<'_, PyAny>,
) -> PyResult<bool> {
    let py = f.py();
    let (f, _) = integer_multiple(&coefficients_from_python(f)?);
    let (g, _) = integer_multiple(&coefficients_from_python(g)?);
    Ok(py.detach(|| polytool::interlaces_either_order_bigint_coeffs(&f, &g)))
}

/// Whether the coefficients are log-concave: ``a_k^2 >= a_(k-1) a_(k+1)``.
///
/// Exact; the library's ``is_log_concave`` convention (see its rustdoc for
/// internal zeros and signs).
#[pyfunction]
pub(crate) fn is_log_concave(polynomial: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        CoefficientProperty::LogConcave,
    ))
}

/// Whether the coefficients are ultra-log-concave of order the degree ``d``:
/// ``a_k / binom(d, k)`` is log-concave.
#[pyfunction]
pub(crate) fn is_ultra_log_concave(polynomial: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        CoefficientProperty::UltraLogConcave,
    ))
}

/// Whether the coefficient sequence is unimodal (weakly increasing, then
/// weakly decreasing).
#[pyfunction]
pub(crate) fn is_unimodal(polynomial: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        CoefficientProperty::Unimodal,
    ))
}

/// Whether the coefficient list reads the same backwards (``a_k = a_(d-k)``).
///
/// With ``ignore_initial_zeros=True`` the leading zero coefficients (a factor
/// ``t^m``) are removed first, as in ``polytool properties``.
#[pyfunction]
#[pyo3(signature = (polynomial, *, ignore_initial_zeros = false))]
pub(crate) fn is_palindromic(
    polynomial: &Bound<'_, PyAny>,
    ignore_initial_zeros: bool,
) -> PyResult<bool> {
    let property = if ignore_initial_zeros {
        CoefficientProperty::PalindromicIgnoringInitialZeros
    } else {
        CoefficientProperty::Palindromic
    };
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        property,
    ))
}

/// Whether the polynomial is palindromic with a nonnegative gamma expansion
/// ``sum_i gamma_i t^i (1 + t)^(d - 2i)``.
///
/// With ``ignore_initial_zeros=True`` a factor ``t^m`` is removed first.
#[pyfunction]
#[pyo3(signature = (polynomial, *, ignore_initial_zeros = false))]
pub(crate) fn is_gamma_positive(
    polynomial: &Bound<'_, PyAny>,
    ignore_initial_zeros: bool,
) -> PyResult<bool> {
    let property = if ignore_initial_zeros {
        CoefficientProperty::GammaPositiveIgnoringInitialZeros
    } else {
        CoefficientProperty::GammaPositive
    };
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        property,
    ))
}

/// The gamma vector ``[gamma_0, gamma_1, ...]`` of a palindromic polynomial,
/// ``p = sum_i gamma_i t^i (1 + t)^(d - 2i)``, or ``None`` when ``p`` is not
/// palindromic.
///
/// Exact (``int``/``Fraction`` entries).  With ``ignore_initial_zeros=True``
/// a factor ``t^m`` is removed first.
///
/// >>> polytool.gamma_coefficients([1, 11, 11, 1])
/// [1, 8]
#[pyfunction]
#[pyo3(signature = (polynomial, *, ignore_initial_zeros = false))]
pub(crate) fn gamma_coefficients<'py>(
    polynomial: &Bound<'py, PyAny>,
    ignore_initial_zeros: bool,
) -> PyResult<Option<Bound<'py, PyAny>>> {
    let coefficients = coefficients_from_python(polynomial)?;
    gamma_coefficients_of(&coefficients, ignore_initial_zeros)
        .map(|gamma| Ok(rationals_to_python(polynomial.py(), &gamma)?.into_any()))
        .transpose()
}

/// Whether the polynomial has no repeated root (over the complex numbers).
#[pyfunction]
pub(crate) fn has_simple_roots(polynomial: &Bound<'_, PyAny>) -> PyResult<bool> {
    Ok(has_property(
        &coefficients_from_python(polynomial)?,
        CoefficientProperty::SimpleRoots,
    ))
}

/// The discriminant of an integer polynomial (an ``int``).
///
/// Raises ``ValueError`` for a non-integer coefficient.
#[pyfunction]
pub(crate) fn discriminant(polynomial: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    let coefficients =
        integer_coefficients(&coefficients_from_python(polynomial)?, "discriminant")?;
    Ok(polytool::discriminant_bigint_coeffs(&coefficients))
}

/// The resultant ``Res(f, g)`` of two integer polynomials (an ``int``).
///
/// Raises ``ValueError`` for a non-integer coefficient.
#[pyfunction]
pub(crate) fn resultant(f: &Bound<'_, PyAny>, g: &Bound<'_, PyAny>) -> PyResult<BigInt> {
    let f = integer_coefficients(&coefficients_from_python(f)?, "resultant")?;
    let g = integer_coefficients(&coefficients_from_python(g)?, "resultant")?;
    Ok(polytool::resultant_bigint_coeffs(&f, &g))
}

/// The coefficient properties of one polynomial, as a ``dict``.
///
/// The keys are those of an item of ``polytool properties --json``:
/// ``polynomial`` (text), ``coefficients``, ``degree``, ``real_rooted``,
/// ``simple_roots``, ``palindromic`` and ``gamma_positive`` (both ignoring
/// initial zeros), ``gamma_coefficients`` (``None`` unless palindromic),
/// ``unimodal``, ``log_concave`` and ``ultra_log_concave``.  Numbers are exact
/// Python ``int``/``Fraction`` values where the CLI writes decimal strings.
///
/// >>> polytool.properties([1, 4, 1])["gamma_coefficients"]
/// [1, 2]
#[pyfunction]
pub(crate) fn properties<'py>(polynomial: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyDict>> {
    let coefficients = coefficients_from_python(polynomial)?;
    properties_dict(polynomial.py(), &coefficients)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(numerator.into(), denominator.into())
    }

    #[test]
    fn rational_input_is_decided_on_an_integer_multiple() {
        // (t - 1/2)(t - 1/3) = 1/6 - 5/6 t + t^2.
        let p = [q(1, 6), q(-5, 6), q(1, 1)];
        assert!(is_real_rooted_coefficients(&p));
        assert_eq!(count_real_roots_coefficients(&p, false), 2);
        // 1/2 (1 + t)^2: gamma = [1/2, 0].
        let gamma = gamma_coefficients_of(&[q(1, 2), q(1, 1), q(1, 2)], false).unwrap();
        assert_eq!(gamma, [q(1, 2), q(0, 1)]);
    }

    #[test]
    fn integer_text_matches_the_cli_format() {
        let p = [q(1, 1), q(11, 1), q(11, 1), q(1, 1)];
        assert_eq!(polynomial_text(&p), "1 + 11t + 11t^2 + t^3");
        assert_eq!(
            polynomial_text(&[q(1, 2), q(0, 1), q(-3, 4), q(-1, 1)]),
            "1/2 - (3/4)t^2 - t^3"
        );
        assert_eq!(polynomial_text(&[q(0, 1), q(-1, 2)]), "-(1/2)t");
    }

    #[test]
    fn directed_interlacing_returns_none_for_the_wrong_order() {
        let f = [q(8, 1), q(-6, 1), q(1, 1)];
        let g = [q(-15, 1), q(23, 1), q(-9, 1), q(1, 1)];
        assert_eq!(interlaces_coefficients(&f, &g, true), Some(true));
        assert_eq!(interlaces_coefficients(&g, &f, false), None);
    }
}
