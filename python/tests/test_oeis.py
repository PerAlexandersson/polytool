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


def test_bundled_recurrence_matches_the_rows():
    for sequence in ["A008292", "A001263"]:
        recurrence = oeis.recurrence(sequence)
        info = oeis.info(sequence)
        assert recurrence.first_index == info["first_index"]
        rows = oeis.generate(sequence, 12)
        assert recurrence.generate(12) == [Polynomial(row) for row in rows]
