//! Conversions between Python objects and the exact library types.
//!
//! Coefficients are exact: Python `int` maps to [`BigInt`], and
//! `fractions.Fraction` (or any `numbers.Rational` with integer `numerator`
//! and `denominator`) maps to [`BigRational`].  Floats are refused, never
//! rounded.  Results come back as `int` when the value is an integer and as
//! `fractions.Fraction` otherwise.  Reports are converted to `dict` through
//! their stable JSON schemas (`json.loads` of the library's serialization).

use num_bigint::BigInt;
use num_traits::Zero;
use polytool::BigRational;
use pyo3::exceptions::{PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyBool, PyFloat, PyInt, PyList, PyString, PyType};

use crate::errors::PolytoolError;
use crate::polynomial::Polynomial;

static FRACTION_TYPE: PyOnceLock<Py<PyType>> = PyOnceLock::new();

fn fraction_type(py: Python<'_>) -> PyResult<&Bound<'_, PyType>> {
    FRACTION_TYPE.import(py, "fractions", "Fraction")
}

/// Whether `value` is an exact scalar that [`rational_from_python`] reads: an
/// `int`, an integer type with `__index__` (NumPy, gmpy2, SageMath), or a
/// rational with `numerator` and `denominator` attributes or methods
/// (`fractions.Fraction`, gmpy2, SageMath) that is not a sequence.
pub(crate) fn is_exact_scalar(value: &Bound<'_, PyAny>) -> PyResult<bool> {
    if value.is_instance_of::<PyInt>() {
        return Ok(true);
    }
    let number_like = value.hasattr("__index__")?
        || (value.hasattr("numerator")? && value.hasattr("denominator")?);
    Ok(number_like && value.try_iter().is_err())
}

/// Read one exact number: `int`, an integer type with `__index__`, or a
/// rational with integer `numerator` and `denominator` (attributes, as in
/// `fractions.Fraction`, or methods, as in SageMath).  Floats and strings are
/// refused.
pub(crate) fn rational_from_python(value: &Bound<'_, PyAny>) -> PyResult<BigRational> {
    if value.is_instance_of::<PyFloat>() {
        return Err(PyTypeError::new_err(
            "polytool needs exact numbers (int or fractions.Fraction), not float; \
             use fractions.Fraction(\"1/3\") for a rational",
        ));
    }
    if value.is_instance_of::<PyString>() {
        return Err(PyTypeError::new_err(
            "expected an exact number (int or fractions.Fraction), not str; \
             use polytool.Polynomial(\"1 + 2t\") to parse polynomial text",
        ));
    }
    if value.is_instance_of::<PyInt>() || value.is_instance_of::<PyBool>() {
        return Ok(BigRational::from_integer(value.extract::<BigInt>()?));
    }
    if value.hasattr("__index__")? {
        if let Ok(integer) = value.extract::<BigInt>() {
            return Ok(BigRational::from_integer(integer));
        }
    }
    if value.hasattr("numerator")? && value.hasattr("denominator")? {
        let part = |name: &str| -> PyResult<BigInt> {
            let mut part = value.getattr(name)?;
            if part.is_callable() {
                part = part.call0()?;
            }
            part.extract().map_err(|_| {
                PyTypeError::new_err(format!(
                    "the {name} of a rational coefficient must be an integer"
                ))
            })
        };
        let (numerator, denominator) = (part("numerator")?, part("denominator")?);
        if denominator.is_zero() {
            return Err(PyZeroDivisionError::new_err(
                "rational coefficient with zero denominator",
            ));
        }
        return Ok(BigRational::new(numerator, denominator));
    }
    Err(PyTypeError::new_err(format!(
        "expected an exact number (int or fractions.Fraction), got {}",
        value.get_type().name()?
    )))
}

/// An exact rational as `int` (integer values) or `fractions.Fraction`.
pub(crate) fn rational_to_python<'py>(
    py: Python<'py>,
    value: &BigRational,
) -> PyResult<Bound<'py, PyAny>> {
    if value.is_integer() {
        return Ok(value.numer().clone().into_pyobject(py)?.into_any());
    }
    fraction_type(py)?.call1((value.numer().clone(), value.denom().clone()))
}

/// A list of exact rationals (see [`rational_to_python`]).
pub(crate) fn rationals_to_python<'py>(
    py: Python<'py>,
    values: &[BigRational],
) -> PyResult<Bound<'py, PyList>> {
    let items = values
        .iter()
        .map(|value| rational_to_python(py, value))
        .collect::<PyResult<Vec<_>>>()?;
    PyList::new(py, items)
}

/// Read a polynomial given as a `Polynomial`, a sequence of exact
/// coefficients in ascending degree order, or polynomial text with integer
/// coefficients (parsed by `polytool::parse_polynomial_bigint`).
///
/// Trailing zero coefficients are dropped; the zero polynomial is empty.
pub(crate) fn coefficients_from_python(value: &Bound<'_, PyAny>) -> PyResult<Vec<BigRational>> {
    untrimmed_coefficients_from_python(value).map(trimmed)
}

/// Read coefficients while retaining supplied padding for recurrence records.
fn untrimmed_coefficients_from_python(value: &Bound<'_, PyAny>) -> PyResult<Vec<BigRational>> {
    if let Ok(polynomial) = value.cast::<Polynomial>() {
        return Ok(polynomial.get().coefficients().to_vec());
    }
    if let Ok(text) = value.cast::<PyString>() {
        let coefficients =
            polytool::parse_polynomial_bigint(&text.to_cow()?).map_err(polytool_error)?;
        return Ok(coefficients
            .into_iter()
            .map(BigRational::from_integer)
            .collect());
    }
    if value.is_instance_of::<PyFloat>() || is_exact_scalar(value)? {
        // A bare number is the constant polynomial.
        return Ok(vec![rational_from_python(value)?]);
    }
    let iterator = value.try_iter().map_err(|_| {
        PyTypeError::new_err(format!(
            "expected a polytool.Polynomial, a list of exact coefficients or polynomial text, \
             got {}",
            value
                .get_type()
                .name()
                .map(|n| n.to_string())
                .unwrap_or_default()
        ))
    })?;
    iterator.map(|item| rational_from_python(&item?)).collect()
}

/// Read a sequence of polynomials (see [`coefficients_from_python`]).
pub(crate) fn polynomial_rows_from_python(
    value: &Bound<'_, PyAny>,
) -> PyResult<Vec<Vec<BigRational>>> {
    coefficient_rows_from_python(value).map(|rows| rows.into_iter().map(trimmed).collect())
}

/// Read coefficient rows without discarding padding, as CLI recurrence input does.
pub(crate) fn coefficient_rows_from_python(
    value: &Bound<'_, PyAny>,
) -> PyResult<Vec<Vec<BigRational>>> {
    if value.is_instance_of::<PyString>() {
        return Err(PyTypeError::new_err(
            "expected a sequence of polynomials, not a single string",
        ));
    }
    value
        .try_iter()
        .map_err(|_| PyTypeError::new_err("expected a sequence of polynomials"))?
        .map(|row| untrimmed_coefficients_from_python(&row?))
        .collect()
}

/// Drop trailing zero coefficients (the zero polynomial becomes empty).
pub(crate) fn trimmed(mut coefficients: Vec<BigRational>) -> Vec<BigRational> {
    while coefficients.last().is_some_and(Zero::is_zero) {
        coefficients.pop();
    }
    coefficients
}

/// The integer polynomial `scale * p` with the smallest positive `scale`
/// clearing every denominator, and that scale.
///
/// The predicates bound here (real-rootedness, interlacing, log-concavity,
/// unimodality, palindromicity, gamma-positivity) do not change under
/// multiplication by a positive constant, so rational input is decided on
/// this integer multiple.
pub(crate) fn integer_multiple(coefficients: &[BigRational]) -> (Vec<BigInt>, BigInt) {
    let (scale, integers) = polytool::vec_poly::clear_denominators_with_scale(coefficients);
    (integers, scale)
}

/// Integer coefficients, or a `ValueError` naming the operation when a
/// coefficient is not an integer.
pub(crate) fn integer_coefficients(
    coefficients: &[BigRational],
    operation: &str,
) -> PyResult<Vec<BigInt>> {
    if let Some(c) = coefficients.iter().find(|c| !c.is_integer()) {
        return Err(PyValueError::new_err(format!(
            "{operation} needs integer coefficients; got {c}"
        )));
    }
    Ok(coefficients.iter().map(BigRational::to_integer).collect())
}

/// A library error as `polytool.PolytoolError` (a `ValueError`).
pub(crate) fn polytool_error(error: impl std::fmt::Display) -> PyErr {
    PolytoolError::new_err(error.to_string())
}

/// Parse JSON text with Python's `json` module (exact `int`s, `str` keys).
pub(crate) fn json_text_to_python<'py>(py: Python<'py>, text: &str) -> PyResult<Bound<'py, PyAny>> {
    py.import("json")?.call_method1("loads", (text,))
}

/// A serializable library value as a Python `dict`/`list` (its JSON schema).
pub(crate) fn serialize_to_python<'py, T: serde::Serialize>(
    py: Python<'py>,
    value: &T,
) -> PyResult<Bound<'py, PyAny>> {
    let text = serde_json::to_string(value).map_err(polytool_error)?;
    json_text_to_python(py, &text)
}

/// JSON text from a `str` (taken as it is) or a JSON-compatible Python value
/// (serialized with `json.dumps`).
pub(crate) fn json_text_from_python(value: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(text) = value.cast::<PyString>() {
        return Ok(text.to_cow()?.into_owned());
    }
    value
        .py()
        .import("json")?
        .call_method1("dumps", (value,))?
        .extract()
}

/// Serialize with serde, compactly or with `indent` spaces per level.
pub(crate) fn to_json_text<T: serde::Serialize>(
    value: &T,
    indent: Option<usize>,
) -> PyResult<String> {
    match indent {
        None => serde_json::to_string(value).map_err(polytool_error),
        Some(width) => {
            let indent = vec![b' '; width];
            let formatter = serde_json::ser::PrettyFormatter::with_indent(&indent);
            let mut buffer = Vec::new();
            let mut serializer = serde_json::Serializer::with_formatter(&mut buffer, formatter);
            value.serialize(&mut serializer).map_err(polytool_error)?;
            String::from_utf8(buffer).map_err(polytool_error)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::One;

    fn q(numerator: i64, denominator: i64) -> BigRational {
        BigRational::new(numerator.into(), denominator.into())
    }

    #[test]
    fn integer_multiple_clears_denominators_with_a_positive_scale() {
        let (integers, scale) = integer_multiple(&[q(1, 2), q(-1, 3), q(2, 1)]);
        assert_eq!(scale, BigInt::from(6));
        assert_eq!(integers, [3, -2, 12].map(BigInt::from).to_vec());
        let (integers, scale) = integer_multiple(&[]);
        assert!(integers.is_empty());
        assert_eq!(scale, BigInt::one());
    }

    #[test]
    fn trimmed_drops_trailing_zeros_only() {
        assert_eq!(trimmed(vec![q(0, 1), q(1, 1), q(0, 1)]), [q(0, 1), q(1, 1)]);
        assert!(trimmed(vec![q(0, 1)]).is_empty());
    }

    #[test]
    fn indented_json_matches_serde_pretty_printing() {
        let value = serde_json::json!({"a": [1, 2]});
        assert_eq!(
            to_json_text(&value, Some(2)).unwrap(),
            serde_json::to_string_pretty(&value).unwrap()
        );
        assert_eq!(to_json_text(&value, None).unwrap(), r#"{"a":[1,2]}"#);
    }
}
