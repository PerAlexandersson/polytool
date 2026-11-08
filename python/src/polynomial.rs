//! The `polytool.Polynomial` class: an immutable univariate polynomial in `t`
//! with exact rational coefficients, over the library's
//! `polytool::Polynomial<BigRational>`.

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use num_traits::Zero;
use polytool::recurrence::BivarPoly;
use polytool::BigRational;
use pyo3::exceptions::{PyTypeError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyFloat, PyTuple};

use crate::convert::{
    coefficients_from_python, integer_coefficients, polytool_error, rational_from_python,
    rational_to_python, rationals_to_python,
};
use crate::properties::{self, CoefficientProperty, RealRootEnclosure};

type LibraryPolynomial = polytool::Polynomial<BigRational>;

/// A polynomial in one variable ``t`` with exact rational coefficients.
///
/// ``Polynomial(coefficients)`` takes the coefficients in ascending degree
/// order (entry ``i`` is the coefficient of ``t^i``) as ``int`` or
/// ``fractions.Fraction`` values; floats are refused.  It also accepts
/// polynomial text with integer coefficients (``"1 + 11t + 11t^2 + t^3"`` or
/// ``"1, 11, 11, 1"``, parsed by the library) or another ``Polynomial``.
/// Trailing zero coefficients are dropped, so ``Polynomial([])`` and
/// ``Polynomial([0])`` are the zero polynomial.
///
/// Polynomials are immutable and hashable.  They support ``+``, ``-``,
/// ``*`` (with polynomials and exact numbers), ``**`` (nonnegative integer
/// exponents), division by a nonzero number with ``/``, and polynomial
/// division with ``//``, ``%`` and ``divmod`` over the rationals.  Calling a
/// polynomial evaluates it at an exact number, or composes it with another
/// polynomial.  The property methods (``is_real_rooted`` and so on) are the
/// module functions of the same name applied to ``self``.
///
/// >>> p = polytool.Polynomial([1, 11, 11, 1])
/// >>> p
/// Polynomial([1, 11, 11, 1])
/// >>> str(p)
/// '1 + 11t + 11t^2 + t^3'
/// >>> p.is_real_rooted(), p.gamma_coefficients()
/// (True, [1, 8])
/// >>> p(Fraction(1, 2))
/// Fraction(75, 8)
#[pyclass(module = "polytool", frozen)]
#[derive(Clone)]
pub struct Polynomial {
    polynomial: LibraryPolynomial,
}

impl Polynomial {
    pub(crate) fn from_coefficients(coefficients: Vec<BigRational>) -> Self {
        Self {
            polynomial: LibraryPolynomial::new(coefficients),
        }
    }

    fn from_library(polynomial: LibraryPolynomial) -> Self {
        Self { polynomial }
    }

    /// The trimmed coefficients (empty for the zero polynomial).
    pub(crate) fn coefficients(&self) -> &[BigRational] {
        self.polynomial.coeffs()
    }
}

/// The right operand of an arithmetic operator: a polynomial, an exact
/// number, or something the operator does not handle.
enum Operand {
    Polynomial(LibraryPolynomial),
    NotSupported,
}

fn operand(value: &Bound<'_, PyAny>) -> PyResult<Operand> {
    if let Ok(polynomial) = value.cast::<Polynomial>() {
        return Ok(Operand::Polynomial(polynomial.get().polynomial.clone()));
    }
    if value.is_instance_of::<PyFloat>() {
        // Refuse loudly rather than with a generic "unsupported operand".
        rational_from_python(value)?;
    }
    if crate::convert::is_exact_scalar(value)? {
        return Ok(Operand::Polynomial(LibraryPolynomial::constant(
            rational_from_python(value)?,
        )));
    }
    Ok(Operand::NotSupported)
}

fn not_implemented(py: Python<'_>) -> Py<PyAny> {
    py.NotImplemented()
}

fn wrap(py: Python<'_>, polynomial: LibraryPolynomial) -> PyResult<Py<PyAny>> {
    Ok(Py::new(py, Polynomial::from_library(polynomial))?.into_any())
}

fn nonzero_divisor(divisor: &LibraryPolynomial) -> PyResult<()> {
    if divisor.is_zero() {
        return Err(PyZeroDivisionError::new_err("polynomial division by zero"));
    }
    Ok(())
}

#[pymethods]
impl Polynomial {
    #[new]
    fn new(coefficients: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self::from_coefficients(coefficients_from_python(
            coefficients,
        )?))
    }

    /// The coefficients in ascending degree order (``int`` or ``Fraction``);
    /// empty for the zero polynomial.
    #[getter(coefficients)]
    fn coefficients_list<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Ok(rationals_to_python(py, self.coefficients())?.into_any())
    }

    /// The degree, or ``None`` for the zero polynomial.
    #[getter]
    fn degree(&self) -> Option<usize> {
        self.polynomial.degree()
    }

    /// The coefficient of the highest power, or ``None`` for the zero
    /// polynomial.
    #[getter]
    fn leading_coefficient<'py>(&self, py: Python<'py>) -> PyResult<Option<Bound<'py, PyAny>>> {
        self.polynomial
            .leading_coefficient()
            .map(|c| rational_to_python(py, &c))
            .transpose()
    }

    /// Whether this is the zero polynomial.
    fn is_zero(&self) -> bool {
        self.polynomial.is_zero()
    }

    /// Whether every coefficient is an integer.
    fn has_integer_coefficients(&self) -> bool {
        self.coefficients().iter().all(BigRational::is_integer)
    }

    /// The coefficient of ``t^k`` (``0`` beyond the degree).
    fn __getitem__<'py>(&self, py: Python<'py>, k: usize) -> PyResult<Bound<'py, PyAny>> {
        rational_to_python(py, &self.polynomial.coeff(k))
    }

    /// Iterate over the coefficients (ascending degree order).
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        Ok(rationals_to_python(py, self.coefficients())?
            .into_any()
            .try_iter()?
            .into_any())
    }

    fn __bool__(&self) -> bool {
        !self.polynomial.is_zero()
    }

    fn __eq__(&self, other: &Bound<'_, PyAny>) -> bool {
        other
            .cast::<Polynomial>()
            .is_ok_and(|other| other.get().coefficients() == self.coefficients())
    }

    fn __hash__(&self) -> u64 {
        let mut hasher = DefaultHasher::new();
        self.coefficients().hash(&mut hasher);
        hasher.finish()
    }

    fn __repr__(&self) -> String {
        let coefficients: Vec<String> = self
            .coefficients()
            .iter()
            .map(properties::python_number_text)
            .collect();
        format!("Polynomial([{}])", coefficients.join(", "))
    }

    fn __str__(&self) -> String {
        properties::polynomial_text(self.coefficients())
    }

    /// LaTeX for Jupyter notebooks, for example ``$1 + 11t + 11t^{2} + t^{3}$``.
    fn _repr_latex_(&self) -> String {
        let latex = BivarPoly {
            coeffs: vec![self.coefficients().to_vec()],
        }
        .to_latex();
        format!("${latex}$")
    }

    fn __reduce__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyTuple>> {
        let class = py.get_type::<Polynomial>();
        let arguments = PyTuple::new(py, [rationals_to_python(py, self.coefficients())?])?;
        PyTuple::new(py, [class.into_any(), arguments.into_any()])
    }

    fn __neg__(&self) -> Self {
        Self::from_library(-self.polynomial.clone())
    }

    fn __pos__(&self) -> Self {
        self.clone()
    }

    fn __add__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => wrap(other.py(), self.polynomial.clone() + q),
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __radd__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__add__(other)
    }

    fn __sub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => wrap(other.py(), self.polynomial.clone() - q),
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __rsub__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => wrap(other.py(), q - self.polynomial.clone()),
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __mul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => wrap(other.py(), self.polynomial.clone() * q),
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __rmul__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        self.__mul__(other)
    }

    fn __pow__(&self, exponent: &Bound<'_, PyAny>, modulo: &Bound<'_, PyAny>) -> PyResult<Self> {
        if !modulo.is_none() {
            return Err(PyTypeError::new_err(
                "pow() with a modulus is not supported for polynomials",
            ));
        }
        let exponent: i64 = exponent.extract()?;
        let exponent = usize::try_from(exponent)
            .map_err(|_| PyValueError::new_err("the exponent must be a nonnegative integer"))?;
        Ok(Self::from_library(self.polynomial.pow(exponent)))
    }

    /// Division by a nonzero exact number (use ``//``, ``%`` or ``divmod``
    /// for polynomial division).
    fn __truediv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        if other.cast::<Polynomial>().is_ok() {
            return Err(PyTypeError::new_err(
                "use // and % (or divmod) to divide by a polynomial",
            ));
        }
        match operand(other)? {
            Operand::Polynomial(q) => {
                let divisor = q.leading_coefficient().unwrap_or_else(BigRational::zero);
                if divisor.is_zero() {
                    return Err(PyZeroDivisionError::new_err("division by zero"));
                }
                let inverse = BigRational::from_integer(1.into()) / divisor;
                wrap(other.py(), self.polynomial.scale(&inverse))
            }
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __floordiv__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => {
                nonzero_divisor(&q)?;
                wrap(other.py(), self.polynomial.div_rem(&q).0)
            }
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __mod__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match operand(other)? {
            Operand::Polynomial(q) => {
                nonzero_divisor(&q)?;
                wrap(other.py(), self.polynomial.div_rem(&q).1)
            }
            Operand::NotSupported => Ok(not_implemented(other.py())),
        }
    }

    fn __divmod__(&self, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        let py = other.py();
        match operand(other)? {
            Operand::Polynomial(q) => {
                nonzero_divisor(&q)?;
                let (quotient, remainder) = self.polynomial.div_rem(&q);
                let pair = PyTuple::new(
                    py,
                    [
                        Py::new(py, Self::from_library(quotient))?,
                        Py::new(py, Self::from_library(remainder))?,
                    ],
                )?;
                Ok(pair.into_any().unbind())
            }
            Operand::NotSupported => Ok(not_implemented(py)),
        }
    }

    /// Evaluate at an exact number (``int``/``Fraction``), or compose with a
    /// polynomial: ``p(q)`` is ``p(q(t))``.
    fn __call__<'py>(&self, x: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let py = x.py();
        if let Ok(inner) = x.cast::<Polynomial>() {
            let composed = self.compose_with(inner.get());
            return Ok(Bound::new(py, composed)?.into_any());
        }
        let x = rational_from_python(x)?;
        rational_to_python(py, &self.polynomial.evaluate(&x))
    }

    /// The composition ``p(q(t))``.
    fn compose(&self, inner: &Bound<'_, PyAny>) -> PyResult<Self> {
        let inner = Self::from_coefficients(coefficients_from_python(inner)?);
        Ok(self.compose_with(&inner))
    }

    /// The ``order``-th derivative (default the first).
    #[pyo3(signature = (order = 1))]
    fn derivative(&self, order: usize) -> Self {
        Self::from_library(self.polynomial.nth_derivative(order))
    }

    /// The shifted polynomial ``p(t + a)`` for an exact number ``a``.
    fn shift(&self, a: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self::from_library(
            self.polynomial.shift(&rational_from_python(a)?),
        ))
    }

    /// The dilated polynomial ``p(c t)`` for an exact number ``c``.
    fn dilate(&self, c: &Bound<'_, PyAny>) -> PyResult<Self> {
        Ok(Self::from_library(
            self.polynomial.dilate(&rational_from_python(c)?),
        ))
    }

    /// The reversal ``t^d p(1/t)``.
    ///
    /// ``d`` is ``degree`` when given (it must be at least the degree of
    /// ``p``) and the degree of ``p`` otherwise.  Reversing at the degree of
    /// a polynomial divisible by ``t`` lowers the degree.
    #[pyo3(signature = (degree = None))]
    fn reverse(&self, degree: Option<usize>) -> PyResult<Self> {
        match degree {
            None => Ok(Self::from_library(self.polynomial.reverse())),
            Some(degree) => self
                .polynomial
                .reverse_with_degree(degree)
                .map(Self::from_library)
                .ok_or_else(|| {
                    PyValueError::new_err(format!(
                        "cannot reverse at degree {degree}: the polynomial has degree {}",
                        self.polynomial.degree().unwrap_or(0)
                    ))
                }),
        }
    }

    /// The monic greatest common divisor with ``other`` over the rationals
    /// (zero only when both are zero).
    fn gcd(&self, other: &Bound<'_, PyAny>) -> PyResult<Self> {
        let other = LibraryPolynomial::new(coefficients_from_python(other)?);
        Ok(Self::from_library(self.polynomial.gcd(&other)))
    }

    /// ``p / gcd(p, p')``: the same roots, each with multiplicity one.
    fn squarefree_part(&self) -> Self {
        Self::from_library(self.polynomial.squarefree_part())
    }

    /// See :func:`polytool.is_real_rooted`.
    fn is_real_rooted(&self, py: Python<'_>) -> bool {
        py.detach(|| properties::is_real_rooted_coefficients(self.coefficients()))
    }

    /// See :func:`polytool.count_real_roots`.
    #[pyo3(signature = (*, distinct = false))]
    fn count_real_roots(&self, py: Python<'_>, distinct: bool) -> usize {
        py.detach(|| properties::count_real_roots_coefficients(self.coefficients(), distinct))
    }

    /// See :func:`polytool.isolate_real_roots`.
    #[pyo3(signature = (*, max_width = None))]
    fn isolate_real_roots(
        &self,
        max_width: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Vec<RealRootEnclosure>> {
        let max_width = max_width.map(rational_from_python).transpose()?;
        properties::isolate_real_roots_coefficients(self.coefficients(), max_width)
    }

    /// Whether ``self << other``; see :func:`polytool.interlaces`.
    #[pyo3(signature = (other, *, strict = false))]
    fn interlaces(&self, other: &Bound<'_, PyAny>, strict: bool) -> PyResult<Option<bool>> {
        let other = coefficients_from_python(other)?;
        Ok(properties::interlaces_coefficients(
            self.coefficients(),
            &other,
            strict,
        ))
    }

    /// See :func:`polytool.is_log_concave`.
    fn is_log_concave(&self) -> bool {
        properties::has_property(self.coefficients(), CoefficientProperty::LogConcave)
    }

    /// See :func:`polytool.is_ultra_log_concave`.
    fn is_ultra_log_concave(&self) -> bool {
        properties::has_property(self.coefficients(), CoefficientProperty::UltraLogConcave)
    }

    /// See :func:`polytool.is_unimodal`.
    fn is_unimodal(&self) -> bool {
        properties::has_property(self.coefficients(), CoefficientProperty::Unimodal)
    }

    /// See :func:`polytool.is_palindromic`.
    #[pyo3(signature = (*, ignore_initial_zeros = false))]
    fn is_palindromic(&self, ignore_initial_zeros: bool) -> bool {
        let property = if ignore_initial_zeros {
            CoefficientProperty::PalindromicIgnoringInitialZeros
        } else {
            CoefficientProperty::Palindromic
        };
        properties::has_property(self.coefficients(), property)
    }

    /// See :func:`polytool.is_gamma_positive`.
    #[pyo3(signature = (*, ignore_initial_zeros = false))]
    fn is_gamma_positive(&self, ignore_initial_zeros: bool) -> bool {
        let property = if ignore_initial_zeros {
            CoefficientProperty::GammaPositiveIgnoringInitialZeros
        } else {
            CoefficientProperty::GammaPositive
        };
        properties::has_property(self.coefficients(), property)
    }

    /// See :func:`polytool.gamma_coefficients`.
    #[pyo3(signature = (*, ignore_initial_zeros = false))]
    fn gamma_coefficients<'py>(
        &self,
        py: Python<'py>,
        ignore_initial_zeros: bool,
    ) -> PyResult<Option<Bound<'py, PyAny>>> {
        properties::gamma_coefficients_of(self.coefficients(), ignore_initial_zeros)
            .map(|gamma| Ok(rationals_to_python(py, &gamma)?.into_any()))
            .transpose()
    }

    /// See :func:`polytool.has_simple_roots`.
    fn has_simple_roots(&self) -> bool {
        properties::has_property(self.coefficients(), CoefficientProperty::SimpleRoots)
    }

    /// See :func:`polytool.discriminant` (integer coefficients only).
    fn discriminant(&self) -> PyResult<num_bigint::BigInt> {
        let coefficients = integer_coefficients(self.coefficients(), "discriminant")?;
        Ok(polytool::discriminant_bigint_coeffs(&coefficients))
    }

    /// See :func:`polytool.properties`.
    fn properties<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        properties::properties_dict(py, self.coefficients())
    }

    /// Coordinates in the Bernstein basis ``binom(d, j) t^j (1 - t)^(d - j)``,
    /// ``j = 0, ..., d``, of degree ``d = degree`` (default the degree of
    /// ``p``; a larger ``d`` is degree elevation).
    #[pyo3(signature = (degree = None))]
    fn bernstein_coordinates<'py>(
        &self,
        py: Python<'py>,
        degree: Option<usize>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let degree = degree.unwrap_or_else(|| self.polynomial.degree().unwrap_or(0));
        let coordinates = polytool::bernstein_basis_coordinates(&self.polynomial, degree)
            .map_err(polytool_error)?;
        Ok(rationals_to_python(py, &coordinates)?.into_any())
    }

    /// Coordinates in the magic basis ``t^i (1 + t)^(d - i)``, ``i = 0, ...,
    /// d``, of degree ``d = degree`` (default the degree of ``p``).
    #[pyo3(signature = (degree = None))]
    fn magic_coordinates<'py>(
        &self,
        py: Python<'py>,
        degree: Option<usize>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let degree = degree.unwrap_or_else(|| self.polynomial.degree().unwrap_or(0));
        let basis = polytool::magic_basis::<BigRational>(degree);
        let coordinates =
            polytool::coordinates_in_basis(&self.polynomial, &basis).map_err(polytool_error)?;
        Ok(rationals_to_python(py, &coordinates)?.into_any())
    }

    /// Coordinates ``c`` with ``p = sum_i c[i] basis[i]``, for a basis of the
    /// polynomials of degree ``< len(basis)`` (a list of polynomials).
    ///
    /// Raises :class:`PolytoolError` when the list is not a basis or ``p`` has
    /// too high a degree.
    fn coordinates_in_basis<'py>(&self, basis: &Bound<'py, PyAny>) -> PyResult<Bound<'py, PyAny>> {
        let py = basis.py();
        let basis = crate::convert::polynomial_rows_from_python(basis)?
            .into_iter()
            .map(LibraryPolynomial::new)
            .collect::<Vec<_>>();
        let coordinates =
            polytool::coordinates_in_basis(&self.polynomial, &basis).map_err(polytool_error)?;
        Ok(rationals_to_python(py, &coordinates)?.into_any())
    }
}

impl Polynomial {
    fn compose_with(&self, inner: &Polynomial) -> Polynomial {
        Self::from_coefficients(polytool::vec_poly::compose(
            self.coefficients(),
            inner.coefficients(),
        ))
    }
}
