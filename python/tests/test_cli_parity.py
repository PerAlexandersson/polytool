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
    return "".join(",".join(str(c) for c in row) + "\n" for row in rows)


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
def test_recurrence_record(eulerian_rows, cli_args, options):
    expected = run("recurrence", "--json", *cli_args, stdin=rows_text(eulerian_rows))
    found = polytool.find_recurrence(eulerian_rows, **options)
    assert found.to_json(indent=2) == expected.rstrip("\n")


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
