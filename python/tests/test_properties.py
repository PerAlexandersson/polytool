"""Real roots, interlacing and coefficient properties."""

from fractions import Fraction

import pytest

import polytool
from polytool import Polynomial


def test_real_rootedness():
    assert polytool.is_real_rooted([1, 11, 11, 1])
    assert not polytool.is_real_rooted([1, 0, 1])
    assert polytool.is_real_rooted([])  # zero polynomial
    assert polytool.is_real_rooted([5])
    # (t - 1/2)(t - 1/3): rational input is scaled to integers.
    assert polytool.is_real_rooted([Fraction(1, 6), Fraction(-5, 6), 1])
    assert Polynomial("1 + 4t + t^2").is_real_rooted()


def test_counting_and_isolating_roots():
    # t (t + 1)^2 (t^2 - 2): roots -sqrt 2, -1 (double), 0, sqrt 2.
    p = [0, -2, -4, -1, 2, 1]
    assert polytool.count_real_roots(p) == 5
    assert polytool.count_real_roots(p, distinct=True) == 4
    roots = polytool.isolate_real_roots(p)
    assert [r.multiplicity for r in roots] == [1, 2, 1, 1]
    assert roots[1].is_exact and roots[1].lower == roots[1].upper == -1
    assert roots[2].lower == 0
    for root in roots:
        assert root.lower <= root.upper
        assert root.upper - root.lower <= Fraction(1, 1024)
    narrow = polytool.isolate_real_roots([-2, 0, 1], max_width=Fraction(1, 10**12))
    assert all(r.upper - r.lower <= Fraction(1, 10**12) for r in narrow)
    assert narrow[1].lower ** 2 < 2 < narrow[1].upper ** 2
    assert polytool.isolate_real_roots([]) == []
    with pytest.raises(ValueError):
        polytool.isolate_real_roots([1, 1], max_width=0)
    assert "multiplicity=2" in repr(roots[1])


def test_interlacing_conventions():
    f = [8, -6, 1]  # (t-2)(t-4)
    g = [-15, 23, -9, 1]  # (t-1)(t-3)(t-5)
    assert polytool.interlaces(f, g) is True
    assert polytool.interlaces(f, g, strict=True) is True
    assert polytool.interlaces(g, f) is None
    assert polytool.interlaces_either_order(g, f)
    assert polytool.interlaces([4, -5, 1], [6, -5, 1]) is False
    assert not polytool.interlaces_either_order([4, -5, 1], [6, -5, 1])
    assert Polynomial(f).interlaces(g)


def test_interlacing_with_big_integers():
    center = 10**20
    f = [-center, 1]
    g = [(center - 1) * (center + 1), -2 * center, 1]
    assert polytool.interlaces(f, g, strict=True) is True


def test_coefficient_properties():
    eulerian = [1, 11, 11, 1]
    assert polytool.is_log_concave(eulerian)
    assert polytool.is_ultra_log_concave(eulerian)
    assert polytool.is_unimodal(eulerian)
    assert polytool.is_palindromic(eulerian)
    assert polytool.is_gamma_positive(eulerian)
    assert polytool.gamma_coefficients(eulerian) == [1, 8]
    assert polytool.gamma_coefficients([1, 2, 3]) is None
    assert not polytool.is_unimodal([1, 0, 1])
    # t (1 + t): palindromic after dropping the factor t.
    assert not polytool.is_palindromic([0, 1, 1])
    assert polytool.is_palindromic([0, 1, 1], ignore_initial_zeros=True)
    assert polytool.gamma_coefficients([0, 1, 1], ignore_initial_zeros=True) == [1]
    # Positive scaling: gamma(p / 2) = gamma(p) / 2.
    assert polytool.gamma_coefficients([Fraction(1, 2), Fraction(11, 2), Fraction(11, 2),
                                        Fraction(1, 2)]) == [Fraction(1, 2), 4]
    assert polytool.has_simple_roots([-1, 0, 1])
    assert not polytool.has_simple_roots([1, 2, 1])


def test_discriminant_and_resultant_are_integer_only():
    assert polytool.discriminant([-2, 0, 1]) == 8
    assert polytool.resultant([-1, 1], [-2, 1]) in (1, -1)
    assert Polynomial([1, 2, 1]).discriminant() == 0
    with pytest.raises(ValueError, match="integer"):
        polytool.discriminant([Fraction(1, 2), 1])


def test_properties_dict():
    report = polytool.properties([1, 4, 1])
    assert report == {
        "polynomial": "1 + 4t + t^2",
        "coefficients": [1, 4, 1],
        "degree": 2,
        "real_rooted": True,
        "simple_roots": True,
        "palindromic": True,
        "gamma_positive": True,
        "gamma_coefficients": [1, 2],
        "unimodal": True,
        "log_concave": True,
        "ultra_log_concave": True,
    }
    assert Polynomial([1, 4, 1]).properties() == report
    assert polytool.properties([1, 0, 1])["real_rooted"] is False
