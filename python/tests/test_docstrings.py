"""Every public name is documented, and the docstring examples run."""

import doctest
import inspect
from fractions import Fraction

import pytest

import polytool
from polytool import oeis

PUBLIC = [getattr(polytool, name) for name in polytool.__all__ if name != "__version__"]
PUBLIC += [getattr(oeis, name) for name in oeis.__all__]


def _members(obj):
    yield obj
    if inspect.isclass(obj) and obj.__module__ == "polytool":
        for name, member in vars(obj).items():
            if not name.startswith("_") or name == "_repr_latex_":
                yield member


@pytest.mark.parametrize("obj", PUBLIC, ids=lambda obj: getattr(obj, "__name__", repr(obj)))
def test_public_names_have_docstrings(obj):
    for member in _members(obj):
        doc = getattr(member, "__doc__", None)
        assert doc and doc.strip(), f"{obj!r}: {member!r} has no docstring"


@pytest.mark.parametrize("obj", PUBLIC, ids=lambda obj: getattr(obj, "__name__", repr(obj)))
def test_docstring_examples(obj):
    finder = doctest.DocTestFinder(exclude_empty=True)
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    globs = {"polytool": polytool, "Fraction": Fraction}
    name = getattr(obj, "__name__", "object")
    for test in finder.find(obj, name, globs=globs):
        runner.run(test)
    result = runner.summarize(verbose=False)
    assert result.failed == 0
