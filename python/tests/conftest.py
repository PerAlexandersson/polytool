"""Shared fixtures: small recurrence records from the library documentation."""

import pytest

EULERIAN_ROWS = [
    [1],
    [1, 1],
    [1, 4, 1],
    [1, 11, 11, 1],
    [1, 26, 66, 26, 1],
    [1, 57, 302, 302, 57, 1],
    [1, 120, 1191, 2416, 1191, 120, 1],
    [1, 247, 4293, 15619, 15619, 4293, 247, 1],
]


def term(offset, coefficients, derivative_order=0):
    return {
        "offset": offset,
        "deriv_order": derivative_order,
        "sign": "none",
        "coeff": {"coeffs": [[str(c) for c in row] for row in coefficients]},
    }


@pytest.fixture
def eulerian_rows():
    return [list(row) for row in EULERIAN_ROWS]


@pytest.fixture
def legendre_record():
    """n P_n = (2n - 1) x P_(n-1) - (n - 1) P_(n-2), P_0 = 1, P_1 = x."""
    return {
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"], ["0", "1"]],
        "recurrence": {
            "terms": [term(1, [[0, -1], [0, 2]]), term(2, [[1], [-1]])],
            "denominator": {"coeffs": [["0"], ["1"]]},
        },
    }


@pytest.fixture
def floor_third_record():
    """P_n = P_(n-1) + x P_(n-3), P_0 = P_1 = P_2 = 1 (OEIS A102547)."""
    return {
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"], ["1"], ["1"]],
        "recurrence": {"terms": [term(1, [[1]]), term(3, [[0, 1]])]},
    }
