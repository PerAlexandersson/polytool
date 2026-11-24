"""The bundled OEIS catalog."""

import pytest

import polytool
from polytool import Polynomial, oeis


def test_ids_and_info():
    ids = oeis.ids()
    assert "A008292" in ids
    assert len(oeis.ids(include_experimental=True)) >= len(ids)
    info = oeis.info("A008292")
    assert info["id"] == "A008292"
    assert info["first_index"] == 1
    assert info["real_rootedness"]["status"] == "lean_proved"
    assert info["recurrence_data"]["schema"] == "polytool.recurrence.v1"
    assert oeis.info("a008292")["id"] == "A008292"
    with pytest.raises(KeyError):
        oeis.info("A000000")


def test_generate_rows_as_listed():
    assert oeis.generate("A008292", 4) == [[1], [1, 1], [1, 4, 1], [1, 11, 11, 1]]
    assert oeis.generate("A008292", 2, from_index=3) == [[1, 4, 1], [1, 11, 11, 1]]
    with pytest.raises(polytool.PolytoolError):
        oeis.generate("A008292", 1, from_index=0)
    big = oeis.generate("A008292", 1, from_index=60)[0]
    assert max(big) > 2**64


@pytest.mark.parametrize(
    "sequence,first_index,source",
    [
        ("A035469", 1, "0803.1659v1"),
        ("A049029", 1, "0803.1659v1"),
        ("A089503", 1, "0803.1659v1"),
        ("A143494", 2, "poissonPolynomial_isPFPolynomial"),
        ("A143495", 3, "poissonPolynomial_isPFPolynomial"),
    ],
)
def test_structural_theorem_attestations(sequence, first_index, source):
    info = oeis.info(sequence)
    assert info["first_index"] == first_index
    flag = info["real_rootedness"]
    assert flag["status"] == "theorem_proved"
    assert flag["as_of"] == "2026-11-24"
    assert source in flag["audit"]


def test_bundled_recurrence_matches_the_rows():
    for sequence in ["A008292", "A001263"]:
        recurrence = oeis.recurrence(sequence)
        info = oeis.info(sequence)
        assert recurrence.first_index == info["first_index"]
        rows = oeis.generate(sequence, 12)
        assert recurrence.generate(12) == [Polynomial(row) for row in rows]


@pytest.mark.parametrize("sequence", ["A266178", "A266298"])
def test_explicit_padded_table_coordinates(sequence):
    rows = oeis.generate(sequence, 8)
    assert [len(row) for row in rows] == list(range(1, 16, 2))
    assert rows[0] == [1]
    for n, row in enumerate(rows[1:], 1):
        assert row[0] == 1
        assert row[1] == (n % 2 if sequence == "A266178" else 1)
        assert row[2:] == [0] * (2 * n - 1)
    assert oeis.generate(sequence, 3, from_index=4) == rows[4:7]
    assert oeis.info(sequence)["row_lengths"] == {"first": 1, "step": 2, "every": 1}
    assert oeis.recurrence(sequence).generate(8) == [Polynomial(row) for row in rows]
