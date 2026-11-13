"""Validate and assemble completed exact region artifacts; no mathematics is recomputed.

Usage: python3 assemble.py ARTIFACT_DIRECTORY [DESTINATION_DIRECTORY]
The inputs are the raw QE JSON, exported polynomial JSON and rational-cell JSON
written by the adjacent Wolfram scripts. All files are regular files.
"""
import hashlib
import argparse
import json
from pathlib import Path
import re

MODES = ["s2plus", "s2minus", "s2zero", "s3plus", "s3minus", "s3zero",
         "s3zeroplus", "s3zerominus", "s3zerozero", "circle2plus", "circle2minus",
         "circle2zero", "circle3plus", "circle3minus", "circle3zero"]

def digest(value):
    return hashlib.sha256(value).hexdigest()

def validate(expression, dimension):
    operation = expression["op"]
    if operation in ("and", "or"):
        for child in expression["args"]:
            validate(child, dimension)
    elif operation == "not":
        validate(expression["arg"], dimension)
    elif operation == "sign":
        assert expression["relation"] in {"Equal", "Unequal", "Less", "LessEqual", "Greater", "GreaterEqual"}
        for monomial in expression["polynomial"]:
            assert re.fullmatch(r"-?\d+", monomial["coefficient"])
            assert len(monomial["powers"]) == dimension
            assert all(isinstance(power, int) and 0 <= power < 1000 for power in monomial["powers"])
    else:
        assert operation in ("true", "false")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", type=Path)
    parser.add_argument("destination", type=Path, nargs="?", default=Path(__file__).resolve().parent)
    parser.add_argument("--regions-only", action="store_true", help="assemble validated geometry before the separate cell audit")
    options = parser.parse_args()
    inputs, destination = options.artifacts, options.destination
    destination.mkdir(parents=True, exist_ok=True)
    sources = {digest(path.read_bytes()): path.name for path in Path(__file__).resolve().parent.glob("*.wls")}
    entries, provenance, samples, cell_provenance = {}, {}, [], {}
    for mode in MODES:
        raw_path = inputs / (("qe-strict-" if mode in {"s3plus", "s3minus", "s3zero"} else "qe-") + mode + ".json")
        raw_bytes = raw_path.read_bytes()
        raw = json.loads(raw_bytes)
        exported = json.loads((inputs / ("region-" + mode + ".json")).read_text())
        assert raw["resolved"] and exported["equivalence_verified"], mode
        cells = None if options.regions_only else json.loads((inputs / ("cells-" + mode + ".json")).read_text())
        if cells is not None:
            assert cells["complete"], mode
        assert exported["generator_sha256"] == raw["script_sha256"]
        assert raw["script_sha256"] in sources and exported["exporter_sha256"] in sources
        if cells is not None:
            assert cells["script_sha256"] in sources
            cell_provenance[mode] = dict(sampler=sources[cells["script_sha256"]],
                sampler_sha256=cells["script_sha256"], version=cells["version"],
                input_sha256=digest(raw_bytes),
                output_sha256=digest((inputs / ("cells-" + mode + ".json")).read_bytes()))
        assert exported["input_sha256"] == digest(raw["input"].encode())
        assert exported["output_sha256"] == digest(raw_bytes)
        assert exported["variables"] == (["c"] if mode.startswith(("s2", "circle2")) else
                                         ["d"] if mode in {"s3zeroplus", "s3zerominus", "s3zerozero"} else ["c", "d"])
        validate(exported["formula"], len(exported["variables"]))
        entries[mode] = exported
        provenance[mode] = dict(raw_output="qe/" + mode + ".json", generator=sources[raw["script_sha256"]],
                                exporter=sources[exported["exporter_sha256"]], qe_seconds=raw["seconds"],
                                generator_sha256=raw["script_sha256"], exporter_sha256=exported["exporter_sha256"])
        for cell in ([] if cells is None else cells["samples"]):
            assert cell["rational"] and cell["sample_verified"], (mode, cell)
            assert all(re.fullmatch(r"-?\d+(?:/\d+)?", value) for value in cell["point"])
            samples.append(dict(mode=mode, **cell))
        (destination / "qe").mkdir(exist_ok=True)
        (destination / "qe" / (mode + ".json")).write_bytes(raw_bytes)
    result = dict(schema="polytool.bkw-regions.v1", version=1,
                  normalization="b0=0; sigma=sign(b1); c=b2^2/abs(b1)^3,d=b3/abs(b1)^2; b1=0 uses c=b2^2,d=b3",
                  maximum_bandwidth=3, arithmetic="integer polynomial signs at exact rational invariants",
                  entries=entries, provenance=provenance)
    (destination / "regions-v1.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    if not options.regions_only:
        (destination / "cells-v1.json").write_text(json.dumps(dict(schema="polytool.bkw-region-cells.v1", provenance=cell_provenance, samples=samples), indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(dict(regions=len(entries), rational_cells=None if options.regions_only else len(samples), destination=str(destination))))

if __name__ == "__main__":
    main()
