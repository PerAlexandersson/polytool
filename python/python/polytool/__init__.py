"""Exact polynomial, recurrence and real-rootedness tools.

Python bindings of the Rust library polytool.  Every computation is exact:
coefficients are Python ``int`` or ``fractions.Fraction`` values in ascending
degree order (entry ``i`` is the coefficient of ``t^i``), and floats are
refused.  Reports convert to ``dict`` through the same stable JSON schemas as
the ``polytool`` command-line tool.

Main entry points:

- :class:`Polynomial` and the polynomial functions :func:`is_real_rooted`,
  :func:`count_real_roots`, :func:`isolate_real_roots`, :func:`interlaces`,
  :func:`interlaces_either_order`, :func:`is_log_concave`,
  :func:`is_ultra_log_concave`, :func:`is_unimodal`, :func:`is_palindromic`,
  :func:`is_gamma_positive`, :func:`gamma_coefficients`,
  :func:`has_simple_roots`, :func:`discriminant`, :func:`resultant` and
  :func:`properties`;
- :func:`find_recurrence` and the :class:`Recurrence` class;
- :func:`prove_real_rooted` (:class:`ProofReport`) and
  :func:`certify_leading_coefficient` (:class:`LeadingCoefficientCertificate`);
- :mod:`polytool.oeis`, the bundled catalog of OEIS polynomial families.

See ``doc/python.md`` in the polytool repository for a tutorial.
"""

from . import _polytool
from ._polytool import (
    BudgetExhaustedError,
    LeadingCoefficientCertificate,
    Polynomial,
    PolytoolError,
    ProofReport,
    RealRootEnclosure,
    Recurrence,
    RecurrenceRows,
    VerificationError,
    __version__,
    certify_leading_coefficient,
    count_real_roots,
    discriminant,
    find_recurrence,
    gamma_coefficients,
    has_simple_roots,
    interlaces,
    interlaces_either_order,
    is_gamma_positive,
    is_log_concave,
    is_palindromic,
    is_real_rooted,
    is_ultra_log_concave,
    is_unimodal,
    isolate_real_roots,
    properties,
    prove_real_rooted,
    resultant,
)
from . import oeis

__all__ = [
    "BudgetExhaustedError",
    "LeadingCoefficientCertificate",
    "Polynomial",
    "PolytoolError",
    "ProofReport",
    "RealRootEnclosure",
    "Recurrence",
    "RecurrenceRows",
    "VerificationError",
    "__version__",
    "certify_leading_coefficient",
    "count_real_roots",
    "discriminant",
    "find_recurrence",
    "gamma_coefficients",
    "has_simple_roots",
    "interlaces",
    "interlaces_either_order",
    "is_gamma_positive",
    "is_log_concave",
    "is_palindromic",
    "is_real_rooted",
    "is_ultra_log_concave",
    "is_unimodal",
    "isolate_real_roots",
    "oeis",
    "properties",
    "prove_real_rooted",
    "resultant",
]

