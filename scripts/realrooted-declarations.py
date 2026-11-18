#!/usr/bin/env python3
"""Regenerate the standalone RealRooted citation inventory from src/.

Run from any directory, with --check to compare without writing. Includes
comments as well as string literals, so every public citation is reviewed.
This records source spelling, not Lean replay or independent theorem validation.
"""

import argparse
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    names = set()
    for source in (root / "src").rglob("*.rs"):
        names.update(name.rstrip(".") for name in re.findall(
            r"RealRooted\.[A-Za-z_][A-Za-z_0-9.]*", source.read_text()))
    expected = "".join(name + "\n" for name in sorted(names))
    path = root / "tests/fixtures/realrooted-declarations.txt"
    if args.check:
        if path.read_text() != expected:
            raise SystemExit("RealRooted inventory differs; review citations and regenerate")
    else:
        path.write_text(expected)


if __name__ == "__main__":
    main()
