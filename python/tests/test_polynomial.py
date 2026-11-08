"""The Polynomial class: exact construction, arithmetic and transforms."""

import pickle
from fractions import Fraction

import pytest

import polytool
from polytool import Polynomial

BIG = 10**100 + 7


def test_construction_trims_and_accepts_exact_numbers():
    assert Polynomial([1, 2, 0, 0]).coefficients == [1, 2]
    assert Polynomial([]).coefficients == []
    assert Polynomial([0]).is_zero()
    assert Polynomial([0]).degree is None
    assert Polynomial([Fraction(1, 2), True]).coefficients == [Fraction(1, 2), 1]
    assert Polynomial((c for c in [3, 4])).coefficients == [3, 4]
    assert Polynomial(Polynomial([5])).coefficients == [5]
    assert Polynomial(7).coefficients == [7]


def test_text_input_uses_the_library_parser():
    assert Polynomial("1 + 11t + 11t^2 + t^3").coefficients == [1, 11, 11, 1]
    assert Polynomial("1, 4, 1").coefficients == [1, 4, 1]
    with pytest.raises(polytool.PolytoolError):
        Polynomial("1 +* t")


def test_floats_and_other_inexact_values_are_refused():
    with pytest.raises(TypeError, match="float"):
        Polynomial([1.5])
    with pytest.raises(TypeError, match="float"):
        Polynomial([1]) + 0.5
    with pytest.raises(TypeError):
        Polynomial([object()])
    with pytest.raises(TypeError):
        Polynomial(["1"])


class MethodRational:
    """A rational with numerator() and denominator() methods, as in SageMath."""

    def __init__(self, numerator, denominator):
        self._numerator, self._denominator = numerator, denominator

    def numerator(self):
        return self._numerator

    def denominator(self):
        return self._denominator


class IndexInteger:
    """An integer type with __index__ only (NumPy, gmpy2 and Sage integers)."""

    def __init__(self, value):
        self.value = value

    def __index__(self):
        return self.value


def test_other_exact_number_types():
    p = Polynomial([MethodRational(1, 3), IndexInteger(10**40)])
    assert p.coefficients == [Fraction(1, 3), 10**40]
    assert (p + IndexInteger(1)).coefficients == [Fraction(4, 3), 10**40]
    assert (MethodRational(2, 3) * Polynomial([3])).coefficients == [2]
    assert p(IndexInteger(0)) == Fraction(1, 3)
    with pytest.raises(ZeroDivisionError):
        Polynomial([MethodRational(1, 0)])


def test_big_integers_and_fractions_round_trip_exactly():
    p = Polynomial([BIG, -BIG, Fraction(BIG, 3)])
    assert p.coefficients == [BIG, -BIG, Fraction(BIG, 3)]
    assert type(p.coefficients[0]) is int
    assert type(p.coefficients[2]) is Fraction
    assert p[2] == Fraction(BIG, 3)
    assert p[10] == 0
    assert list(p) == [BIG, -BIG, Fraction(BIG, 3)]
    assert p(1) == Fraction(BIG, 3)
    assert (p * 3).coefficients == [3 * BIG, -3 * BIG, BIG]
    assert p.leading_coefficient == Fraction(BIG, 3)


def test_arithmetic_with_polynomials_and_numbers():
    p = Polynomial([1, 1])
    q = Polynomial([1, 2])
    assert (p + q).coefficients == [2, 3]
    assert (p - q).coefficients == [0, -1]
    assert (p * q).coefficients == [1, 3, 2]
    assert (p**3).coefficients == [1, 3, 3, 1]
    assert (2 + p).coefficients == [3, 1]
    assert (2 - p).coefficients == [1, -1]
    assert (Fraction(1, 2) * p).coefficients == [Fraction(1, 2), Fraction(1, 2)]
    assert (p / 2).coefficients == [Fraction(1, 2), Fraction(1, 2)]
    assert (-p).coefficients == [-1, -1]
    with pytest.raises(ValueError):
        p ** -1
    with pytest.raises(ZeroDivisionError):
        p / 0
    with pytest.raises(TypeError):
        p / q


def test_polynomial_division_over_the_rationals():
    p = Polynomial([1, 11, 11, 1])
    d = Polynomial([1, 2])
    quotient, remainder = divmod(p, d)
    assert quotient * d + remainder == p
    assert p // d == quotient and p % d == remainder
    assert remainder.degree in (None, 0)
    with pytest.raises(ZeroDivisionError):
        divmod(p, Polynomial([]))
    assert Polynomial([-1, 0, 1]).gcd([1, 2, 1]).coefficients == [1, 1]
    assert Polynomial([1, 2, 1]).squarefree_part().degree == 1


def test_evaluation_composition_and_transforms():
    p = Polynomial([1, 11, 11, 1])
    assert p(Fraction(1, 2)) == Fraction(75, 8)
    assert p(Polynomial([0, 1])) == p
    assert p.compose([1, 1]) == p(Polynomial([1, 1]))
    assert p.shift(1) == p.compose([1, 1])
    assert p.dilate(2).coefficients == [1, 22, 44, 8]
    assert p.derivative().coefficients == [11, 22, 3]
    assert p.derivative(2).coefficients == [22, 6]
    assert Polynomial([0, 1, 2]).reverse().coefficients == [2, 1]
    assert Polynomial([0, 1, 2]).reverse(degree=3).coefficients == [0, 2, 1]
    with pytest.raises(ValueError):
        Polynomial([0, 1, 2]).reverse(degree=1)
    with pytest.raises(TypeError, match="float"):
        p(0.5)


def test_basis_coordinates():
    p = Polynomial([1, 3, 3, 1])  # (1 + t)^3
    assert p.magic_coordinates() == [1, 0, 0, 0]
    assert Polynomial([1, 4, 1]).bernstein_coordinates() == [1, 3, 6]
    basis = [Polynomial([1]), Polynomial([1, 1]), Polynomial([1, 2, 1])]
    assert Polynomial([1, 2, 1]).coordinates_in_basis(basis) == [0, 0, 1]
    with pytest.raises(polytool.PolytoolError):
        Polynomial([0, 0, 0, 1]).coordinates_in_basis(basis)


def test_equality_hash_pickle_and_display():
    p = Polynomial([1, Fraction(-3, 4), 2])
    assert p == Polynomial([1, Fraction(-3, 4), 2, 0])
    assert p != [1, Fraction(-3, 4), 2]
    assert len({p, Polynomial([1, Fraction(-3, 4), 2])}) == 1
    assert pickle.loads(pickle.dumps(p)) == p
    assert eval(repr(p), {"Polynomial": Polynomial, "Fraction": Fraction}) == p
    assert str(p) == "1 - (3/4)t + 2t^2"
    assert str(Polynomial([1, 11, 11, 1])) == "1 + 11t + 11t^2 + t^3"
    assert Polynomial([1, 2])._repr_latex_() == "$1 + 2t$"
    assert not Polynomial([])
    assert Polynomial([1]).has_integer_coefficients()
    assert not p.has_integer_coefficients()
