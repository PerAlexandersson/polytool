"""Recurrence search, row generation, JSON round trips, proofs and degree laws."""

import json
import pickle
from fractions import Fraction

import pytest

import polytool
from polytool import Polynomial, Recurrence


def test_find_recurrence_reproduces_the_rows(eulerian_rows):
    recurrence = polytool.find_recurrence(eulerian_rows)
    assert isinstance(recurrence, Recurrence)
    assert recurrence.is_fitted
    assert recurrence.first_index == 1
    assert recurrence.order == 1
    assert recurrence.search["source_rows"] == len(eulerian_rows)
    rows = recurrence.generate(len(eulerian_rows))
    assert [row.coefficients for row in rows] == eulerian_rows
    tail = recurrence.generate(2, from_index=4)
    assert [row.coefficients for row in tail] == eulerian_rows[3:5]
    with pytest.raises(ValueError):
        recurrence.generate(1, from_index=0)


def test_first_index_and_skip_prefix_follow_the_cli(eulerian_rows):
    recurrence = polytool.find_recurrence(eulerian_rows, first_index=0, skip_prefix=1)
    assert recurrence.first_index == 1  # the skipped row keeps index 0
    assert recurrence.search["skip_prefix"] == 1
    shifted = recurrence.reindexed(10)
    assert shifted.first_index == 10
    assert shifted.generate(3) == recurrence.generate(3)


def test_search_options_and_outcomes(eulerian_rows):
    assert polytool.find_recurrence([[1], [2], [5], [1, 7]], max_recurrence_length=1,
                                    max_variable_degree=0, max_index_degree=0,
                                    max_derivative_order=0) is None
    with pytest.raises(polytool.BudgetExhaustedError) as caught:
        polytool.find_recurrence(eulerian_rows, max_candidates=0)
    assert caught.value.candidates_considered == 0
    assert caught.value.work_units == 0
    assert isinstance(caught.value, polytool.PolytoolError)
    limited = polytool.find_recurrence(eulerian_rows, max_work=10**9)
    assert limited.search["work_units"] > 0
    with pytest.raises(TypeError):
        polytool.find_recurrence(eulerian_rows, max_rec_len=2)


def test_rational_rows_are_supported():
    # P_n = (t + 1/2) P_(n-1): rows (t + 1/2)^n.
    base = Polynomial([Fraction(1, 2), 1])
    rows = [base**n for n in range(8)]
    recurrence = polytool.find_recurrence(rows, first_index=0)
    assert recurrence.generate(8) == rows
    assert recurrence.generate(1, from_index=9)[0] == base**9


def test_lazy_iteration(eulerian_rows):
    recurrence = polytool.find_recurrence(eulerian_rows)
    rows = recurrence.iterate(from_index=2)
    assert next(rows).coefficients == eulerian_rows[1]
    assert next(rows).coefficients == eulerian_rows[2]
    for _ in range(20):
        last = next(rows)
    assert last.degree == 22


def test_json_round_trips(legendre_record):
    recurrence = Recurrence.from_json(legendre_record)
    assert Recurrence.from_json(json.dumps(legendre_record)) == recurrence
    assert Recurrence.from_json(recurrence.to_json()) == recurrence
    assert Recurrence.from_json(recurrence.to_dict()) == recurrence
    assert json.loads(recurrence.to_json(indent=2)) == recurrence.to_dict()
    assert pickle.loads(pickle.dumps(recurrence)) == recurrence
    assert hash(recurrence) == hash(Recurrence.from_json(legendre_record))
    assert not recurrence.is_fitted and recurrence.search is None
    # Legendre: P_2 = (3x^2 - 1)/2.
    assert recurrence.generate(3)[2].coefficients == [Fraction(-1, 2), 0, Fraction(3, 2)]
    with pytest.raises(polytool.PolytoolError):
        Recurrence.from_json({"schema": "polytool.recurrence.v0"})
    with pytest.raises(polytool.PolytoolError):
        Recurrence.from_json("not json")


def test_latex_presentation(eulerian_rows, legendre_record):
    fitted = polytool.find_recurrence(eulerian_rows)
    assert fitted.latex().startswith("% Candidate fitted to finite rows")
    assert fitted.to_dict()["latex"] == fitted.latex()
    supplied = Recurrence.from_json(legendre_record)
    latex = supplied.latex()
    assert latex.splitlines()[0] == "$P_{0} = 1$"
    assert "n \\geq 2" in latex
    assert supplied.latex(expanded=True).count("$") == latex.count("$")
    display = supplied._repr_latex_()
    assert display.startswith("$$\\begin{gathered}") and "%" not in display
    assert "P(n)" in str(supplied) and "Recurrence" in repr(supplied)


def test_prove_real_rooted_and_replay(legendre_record):
    report = polytool.prove_real_rooted(legendre_record, root_region=(-1, 1))
    assert report.status == "proved" and report.is_proved
    assert report.route == "three_term_recurrence"
    certificate = report.certificate
    assert certificate["status"] == "proved"
    assert certificate["conclusion"]["interlacing"] == "strict"
    assert len(report.attempts) == 6
    assert report.verify() is True
    assert report.recurrence == Recurrence.from_json(legendre_record)
    assert polytool.ProofReport.from_json(report.to_json()).verify()
    data = report.to_dict()
    data["status"] = "not_proved"
    with pytest.raises(polytool.VerificationError):
        polytool.ProofReport.from_json(data).verify()
    assert "three_term_recurrence" in repr(report)


def test_not_proved_is_not_a_disproof():
    # P_n = (1 + t^2) P_(n-1): never real-rooted, so no theorem applies.
    record = {
        "schema": "polytool.recurrence.v1",
        "first_index": 0,
        "initial_polynomials": [["1"]],
        "recurrence": {"terms": [{"offset": 1, "deriv_order": 0, "sign": "none",
                                  "coeff": {"coeffs": [["1", "0", "1"]]}}]},
    }
    report = Recurrence.from_json(record).prove_real_rooted(maximum_start_delay=0)
    assert report.status == "not_proved"
    assert report.route is None and report.certificate is None
    assert report.verify()


def test_leading_coefficient_certificate(floor_third_record):
    certificate = polytool.certify_leading_coefficient(floor_third_record)
    assert certificate.is_proved and certificate.status == "proved"
    assert certificate.degree_at(31) == 10
    assert certificate.sign_at(31) == "positive"
    assert certificate.route in {"cone", "recurrence", "subseq", "twoStep", "cancel",
                                 "coefficients"}
    assert certificate.degree_law["period"] == 3
    assert certificate.verify()
    replayed = polytool.LeadingCoefficientCertificate.from_json(certificate.to_dict())
    assert replayed.verify()
    tampered = certificate.to_dict()
    tampered["degree_law"]["growth_per_period"] += 1
    with pytest.raises(polytool.VerificationError):
        polytool.LeadingCoefficientCertificate.from_json(tampered).verify()
    with pytest.raises(ValueError):
        polytool.certify_leading_coefficient(floor_third_record, maximum_period=0)


def test_theorem_matches(eulerian_rows):
    report = polytool.find_recurrence(eulerian_rows).match_theorems()
    assert report["schema"] == "polytool.recurrence-theorems.v1"
    assert "theorem_matches" in report
