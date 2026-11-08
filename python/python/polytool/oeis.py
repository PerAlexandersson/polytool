"""The bundled catalog of recurrence-backed OEIS polynomial families.

- :func:`ids` lists the A-numbers;
- :func:`info` describes one entry (``polytool oeis info --json``);
- :func:`generate` returns rows exactly as the OEIS lists them;
- :func:`recurrence` returns the bundled recurrence as a
  :class:`polytool.Recurrence`.
"""

from ._polytool import oeis as _oeis

ids = _oeis.ids
info = _oeis.info
generate = _oeis.generate
recurrence = _oeis.recurrence

__all__ = ["generate", "ids", "info", "recurrence"]
