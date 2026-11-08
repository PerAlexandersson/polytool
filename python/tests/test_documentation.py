"""The `pycon` examples of doc/python.md run as one doctest session."""

import doctest
import re
from pathlib import Path

import pytest

GUIDE = Path(__file__).resolve().parents[2] / "doc" / "python.md"
BLOCK = re.compile(r"^```pycon\n(.*?)^```$", re.MULTILINE | re.DOTALL)


@pytest.mark.skipif(not GUIDE.exists(), reason="doc/python.md is not next to the package")
def test_guide_examples():
    blocks = BLOCK.findall(GUIDE.read_text(encoding="utf-8"))
    assert blocks, "no pycon examples found"
    parser = doctest.DocTestParser()
    test = parser.get_doctest("\n".join(blocks), {}, "python.md", str(GUIDE), 0)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    runner.run(test)
    assert runner.summarize(verbose=False).failed == 0
