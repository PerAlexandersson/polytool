"""Parity with the JSON output of the polytool command-line tool.

Set POLYTOOL_CLI to the path of a `polytool` binary built from the same
source to run these tests; they are skipped otherwise.
"""

import json
import os
import subprocess

import pytest

import polytool

CLI = os.environ.get("POLYTOOL_CLI")
pytestmark = pytest.mark.skipif(not CLI, reason="POLYTOOL_CLI is not set")


def run(*args, stdin=""):
    completed = subprocess.run(
        [CLI, *args], input=stdin, capture_output=True, text=True, timeout=300
    )
    assert completed.returncode in (0, 1), completed.stderr
    return completed.stdout


def rows_text(rows):
    return "".join((",".join(str(c) for c in row) or "[]") + "\n" for row in rows)


def as_ints(values):
    return None if values is None else [int(v) for v in values]


def test_properties(eulerian_rows):
    rows = eulerian_rows + [[1, 0, 1], [0, 1, 1], [3, -2, 0, 5], [10**30, 1, 10**30]]
    items = json.loads(run("properties", "--json", stdin=rows_text(rows)))["items"]
    for row, item in zip(rows, items):
        item.pop("index")
        item["coefficients"] = as_ints(item["coefficients"])
        item["gamma_coefficients"] = as_ints(item["gamma_coefficients"])
        assert polytool.properties(row) == item


def test_real_rooted(eulerian_rows):
    rows = eulerian_rows + [[1, 0, 1]]
    items = json.loads(run("real-rooted", "--json", stdin=rows_text(rows)))["items"]
    assert [item["real_rooted"] for item in items] == [
        polytool.is_real_rooted(row) for row in rows
    ]


@pytest.mark.parametrize(
    "cli_args, options",
    [((), {}), (("--first-index", "0", "--skip-prefix", "1"),
                {"first_index": 0, "skip_prefix": 1})],
)
@pytest.mark.parametrize("lean_name", [None, "Rows"])
def test_recurrence_record(eulerian_rows, cli_args, options, lean_name):
    if lean_name is not None:
        cli_args = (*cli_args, "--lean-name", lean_name)
    expected = run("recurrence", "--json", *cli_args, stdin=rows_text(eulerian_rows))
    found = polytool.find_recurrence(eulerian_rows, **options, lean_name=lean_name)
    assert found.to_json(indent=2) == expected.rstrip("\n")
    record = found.to_dict()
    assert "presentation" in record and "normalized_recurrence" in record
    assert ("lean" in record) == (lean_name is not None)
    assert "lean_reason" not in record
    assert polytool.Recurrence.from_json(found.to_json()).to_dict() == record


@pytest.mark.parametrize("sequence", ["A065826", "A216118"])
@pytest.mark.parametrize("first_index", [0, 7])
@pytest.mark.parametrize("lean_name", [None, "Rows"])
def test_recurrence_checked_seeds(sequence, first_index, lean_name):
    rows = polytool.oeis.generate(sequence, 12)
    args = ["--first-index", str(first_index)]
    if lean_name is not None:
        args += ["--lean-name", lean_name]
    expected = run("recurrence", "--json", *args, stdin=rows_text(rows))
    found = polytool.find_recurrence(rows, first_index=first_index, lean_name=lean_name)
    assert found.to_json(indent=2) == expected.rstrip("\n")
    assert len(found.initial_polynomials) == 2
    assert found.to_dict()["presentation"]["start_index"] == first_index + 2
    assert found.generate(12) == [polytool.Polynomial(row) for row in rows]


@pytest.mark.parametrize("lean_name", [None, "Rows"])
def test_recurrence_padded_coefficients_and_work_metadata(eulerian_rows, lean_name):
    rows = [row + [0, 0] for row in eulerian_rows]
    args = ["--first-index", "0", "--skip-prefix", "1", "--max-work", "1000000000"]
    if lean_name is not None:
        args += ["--lean-name", lean_name]
    expected = run("recurrence", "--json", *args, stdin=rows_text(rows))
    found = polytool.find_recurrence(
        rows, first_index=0, skip_prefix=1, max_work=10**9, lean_name=lean_name
    )
    assert found.to_json(indent=2) == expected.rstrip("\n")
    assert found.to_dict()["initial_polynomials"][0] == [str(c) for c in rows[1]]
    assert found.search["work_units"] > 0


@pytest.mark.parametrize("lean_name", [None, "Rows"])
def test_recurrence_lean_refusal(lean_name):
    source = {
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["0", "1"]],
        "recurrence": {
            "denominator": {"coeffs": [["-1", "1"], ["1"]]},
            "terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                       "coeff": {"coeffs": [["0", "1"], ["1"]]}}]},
    }
    rows = [row.coefficients for row in polytool.Recurrence.from_json(source).generate(12)]
    options = {
        "first_index": 0, "laurent": False, "denominator": True,
        "max_recurrence_length": 1, "max_variable_degree": 1,
        "max_index_degree": 1, "max_derivative_order": 0,
        "max_denominator_variable_degree": 1, "max_denominator_index_degree": 1,
    }
    args = ["--first-index", "0", "--no-laurent", "--try-denominator"]
    for name in options:
        if name.startswith("max_"):
            args += ["--" + name.replace("_", "-"), str(options[name])]
    if lean_name is not None:
        args += ["--lean-name", lean_name]
    expected = run("recurrence", "--json", *args, stdin=rows_text(rows))
    found = polytool.find_recurrence(rows, **options, lean_name=lean_name)
    assert found.to_json(indent=2) == expected.rstrip("\n")
    record = found.to_dict()
    assert "lean" not in record
    assert ("lean_reason" in record) == (lean_name is not None)
    if lean_name is not None:
        assert "depends on t" in record["lean_reason"]


@pytest.mark.parametrize("lean_name", [None, "Rows"])
def test_recurrence_future_singularity_preserves_all_available_seeds(lean_name):
    source = {
        "schema": "polytool.recurrence.v1", "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {
            "denominator": {"coeffs": [["-20"], ["1"]]},
            "terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                       "coeff": {"coeffs": [["1"], ["1"]]}}]},
    }
    rows = [row.coefficients for row in polytool.Recurrence.from_json(source).generate(12)]
    options = {
        "first_index": 0, "laurent": False, "denominator": True,
        "max_recurrence_length": 1, "max_variable_degree": 0,
        "max_index_degree": 1, "max_derivative_order": 0,
        "max_denominator_variable_degree": 0, "max_denominator_index_degree": 1,
    }
    args = ["--first-index", "0", "--no-laurent", "--try-denominator"]
    for name in options:
        if name.startswith("max_"):
            args += ["--" + name.replace("_", "-"), str(options[name])]
    if lean_name is not None:
        args += ["--lean-name", lean_name]
    expected = run("recurrence", "--json", *args, stdin=rows_text(rows))
    found = polytool.find_recurrence(rows, **options, lean_name=lean_name)
    assert found.to_json(indent=2) == expected.rstrip("\n")
    record = found.to_dict()
    assert len(found.initial_polynomials) == len(rows)
    assert record["presentation"]["status"] == "insufficient_initial_rows"
    assert record["presentation"]["left_hand_zeros"] == ["20"]
    assert record["presentation"]["start_index"] is None
    assert "lean" not in record
    assert ("lean_reason" in record) == (lean_name is not None)


def test_proof_and_degree_reports(eulerian_rows, legendre_record, tmp_path):
    fitted = polytool.find_recurrence(eulerian_rows)
    for recurrence in [fitted, polytool.Recurrence.from_json(legendre_record)]:
        path = tmp_path / "recurrence.json"
        path.write_text(recurrence.to_json(indent=2))
        prove = json.loads(run("recurrence-prove", "--recurrence", str(path), "--json"))
        assert recurrence.prove_real_rooted().to_dict() == prove
        degree = json.loads(run("recurrence-degree", "--recurrence", str(path), "--json"))
        assert recurrence.certify_leading_coefficient().to_dict() == degree


@pytest.mark.parametrize("sequence", ["A008292", "A001263", "A102547"])
def test_oeis(sequence):
    generated = json.loads(run("oeis", "generate", sequence, "--count", "9", "--format", "json"))
    assert polytool.oeis.generate(sequence, 9) == [as_ints(row) for row in generated["rows"]]
    info = json.loads(run("oeis", "info", sequence, "--json"))
    ours = polytool.oeis.info(sequence)
    assert ours.pop("first_index") == info.pop("first_row")
    for key, value in ours.items():
        assert info[key] == value, key
