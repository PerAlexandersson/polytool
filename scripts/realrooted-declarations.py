#!/usr/bin/env python3
"""Check citations against an independently generated RealRooted inventory.

Offline check: scripts/realrooted-declarations.py --check
Regenerate from a RealRooted source declaration index (never from Rust alone):
  python3 "$REALROOTED_ROOT/scripts/agent_index.py" --repo "$REALROOTED_ROOT" \
      --output /tmp/realrooted-index.json
  python3 scripts/realrooted-declarations.py --inventory /tmp/realrooted-index.json

The snapshot records source declarations, namespaces and module references at
an upstream commit. It validates names, not Lean proofs or axiom assumptions.
"""

import argparse
import json
from pathlib import Path
import re


def source_citations(root):
    return {name.rstrip(".") for source in (root / "src").rglob("*.rs")
            for name in re.findall(r"RealRooted\.[A-Za-z_][A-Za-z_0-9.]*", source.read_text())}


def from_inventory(inventory, names):
    revision = inventory["revision"]
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("inventory must identify an upstream commit")
    declarations = {d["name"] for d in inventory["declarations"]}
    modules = {d["module"] for d in inventory["declarations"]}
    namespaces = {".".join(n.split(".")[:i]) for n in declarations
                  for i in range(1, len(n.split(".")))}
    rows = []
    for name in sorted(names):
        kind = next((kind for kind, choices in (("declaration", declarations),
                    ("namespace", namespaces), ("module", modules)) if name in choices), None)
        if kind is None:
            raise ValueError(f"citation absent from independent inventory: {name}")
        rows.append(f"{kind}\t{name}\n")
    return (f"# RealRooted commit {revision}\n"
            "# Source declaration inventory; not a proof or axiom audit.\n" + "".join(rows))


def check_snapshot(text, names):
    lines = text.splitlines()
    if not lines or not re.fullmatch(r"# RealRooted commit [0-9a-f]{40}", lines[0]):
        raise ValueError("missing pinned RealRooted revision")
    rows = [line.split("\t") for line in lines if line and not line.startswith("#")]
    if any(len(row) != 2 or row[0] not in {"declaration", "namespace", "module"} for row in rows):
        raise ValueError("invalid declaration inventory snapshot")
    if {row[1] for row in rows} != names:
        raise ValueError("changed citations: regenerate against an independent --inventory input")


def self_test():
    inventory = {"revision": "a" * 40,
                 "declarations": [{"name": "RealRooted.actual", "module": "RealRooted.Module"}]}
    snapshot = from_inventory(inventory, {"RealRooted.actual", "RealRooted.Module"})
    check_snapshot(snapshot, {"RealRooted.actual", "RealRooted.Module"})
    for action in (
        lambda: from_inventory(inventory, {"RealRooted.invented"}),
        lambda: check_snapshot(snapshot, {"RealRooted.invented"}),
        lambda: from_inventory({**inventory, "revision": "unknown"}, {"RealRooted.actual"}),
    ):
        try:
            action()
        except ValueError:
            continue
        raise AssertionError("invalid independent inventory was accepted")
    print("independent inventory self-tests passed (valid input and three negative controls)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--inventory", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    root = Path(__file__).resolve().parent.parent
    names = source_citations(root)
    path = root / "tests/fixtures/realrooted-declarations.txt"
    if args.inventory:
        expected = from_inventory(json.loads(args.inventory.read_text()), names)
        if args.check:
            if path.read_text() != expected:
                raise ValueError("snapshot differs from independent declaration inventory")
        else:
            path.write_text(expected)
    elif args.check:
        check_snapshot(path.read_text(), names)
    else:
        parser.error("regeneration requires --inventory; Rust citations are not a declaration inventory")


if __name__ == "__main__":
    main()
