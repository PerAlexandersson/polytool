#!/usr/bin/env python3
"""Generate the Rust OEIS catalog from verified recurrence fixtures.

The benchmark recurrence JSON is the executable source.  Local OEIS internal
format files provide names, offsets, keywords, and a prefix check.  Entries
whose flattened prefix cannot be aligned are still generated, but strict
b-file output is disabled for them.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import subprocess
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURES = PROJECT_ROOT / "fixtures" / "recurrence-benchmarks"
DEFAULT_OEIS_DATA = Path("/home/paxinum/OEIS-data/seq")
DEFAULT_OUTPUT = PROJECT_ROOT / "src" / "oeis_catalog_generated.rs"
DEFAULT_SEQUENCE_LIBRARY = Path(
    "/home/paxinum/Dropbox/AI-projects/projects/OEIS-polynomials/sequences"
)
DEFAULT_LEAN_SEQUENCE_REPOSITORY = Path(
    "/home/paxinum/Dropbox/AI-projects/projects/real-rooted-oeis-proofs/ProofsOeis"
)

# The source queue lives in a Dropbox-synced tree where generated cache files
# are intentionally forbidden.
sys.dont_write_bytecode = True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixtures", type=Path, default=DEFAULT_FIXTURES)
    parser.add_argument("--oeis-data", type=Path, default=DEFAULT_OEIS_DATA)
    parser.add_argument("--sequence-library", type=Path, default=DEFAULT_SEQUENCE_LIBRARY)
    parser.add_argument(
        "--lean-sequences", type=Path, default=DEFAULT_LEAN_SEQUENCE_REPOSITORY
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    return parser.parse_args()


def rust_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=False)


def parse_integer_row(line: str) -> list[str]:
    values = [value.strip() for value in line.split(",")]
    if not values or any(re.fullmatch(r"[+-]?\d+", value) is None for value in values):
        raise ValueError(f"fixture row is not integral: {line[:100]}")
    return values


def oeis_record(root: Path, oeis_id: str) -> dict[str, Any]:
    path = root / oeis_id[:4] / f"{oeis_id}.seq"
    name = ""
    offset: int | None = None
    keywords: list[str] = []
    data: list[int] = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("%N "):
            name = line.split(" ", 2)[2]
        elif line.startswith("%O "):
            offset = int(line.split()[2].split(",", 1)[0])
        elif line.startswith("%K "):
            keywords.extend(line.split(" ", 2)[2].split(","))
        elif line.startswith(("%S ", "%T ", "%U ")):
            for token in line.split(" ", 2)[2].split(","):
                token = token.strip()
                if re.fullmatch(r"[+-]?\d+", token):
                    data.append(int(token))
    if not name or offset is None:
        raise ValueError(f"missing OEIS name or offset for {oeis_id}")
    return {"name": name, "offset": offset, "keywords": keywords, "data": data}


def sparse_coefficients(poly: dict[str, Any] | None) -> list[tuple[int, int, str]] | None:
    if poly is None:
        return None
    terms: list[tuple[int, int, str]] = []
    for n_degree, row in enumerate(poly["coeffs"]):
        for variable_degree, value in enumerate(row):
            if value != "0":
                terms.append((n_degree, variable_degree, value))
    return terms


def rational_text(value: Any) -> str:
    numerator, denominator = value.as_numer_denom()
    if denominator == 1:
        return str(numerator)
    return f"{numerator}/{denominator}"


def sympy_sparse(poly: Any, n_symbol: Any, t_symbol: Any) -> list[tuple[int, int, str]]:
    import sympy as sp

    expanded = sp.Poly(sp.expand(poly), n_symbol, t_symbol, domain=sp.QQ)
    terms = [
        (int(powers[0]), int(powers[1]), rational_text(coefficient))
        for powers, coefficient in expanded.terms()
        if coefficient != 0
    ]
    return sorted(terms)


def lean_generated_expression(text: str, oeis_id: str, sympy: Any) -> tuple[Any, set[str]]:
    """Parse the deliberately small expression grammar emitted by the Lean generator."""
    n_symbol, t_symbol = sympy.symbols("N T")
    text = re.sub(
        rf"{oeis_id} \(n \+ (\d+)\)",
        lambda match: f"Y{match.group(1)}D0",
        text,
    )
    text = re.sub(rf"{oeis_id} n\b", "Y0D0", text)
    derivative = re.compile(r"\((Y\d+D\d+)\)\.derivative")
    while derivative.search(text):
        text = derivative.sub(
            lambda match: re.sub(
                r"D(\d+)$",
                lambda order: f"D{int(order.group(1)) + 1}",
                match.group(1),
            ),
            text,
        )
    text = text.replace("(n : ℝ)", "N").replace("X", "T").replace("^", "**")
    placeholders = set(re.findall(r"Y\d+D\d+", text))
    local_symbols = {
        "N": n_symbol,
        "T": t_symbol,
        "C": lambda value: value,
        **{name: sympy.Symbol(name) for name in placeholders},
    }
    unknown = set(re.findall(r"[A-Za-z_][A-Za-z_0-9]*", text)) - set(local_symbols)
    if unknown:
        raise ValueError(f"{oeis_id}: unsupported Lean identifiers {sorted(unknown)}")
    return sympy.sympify(text, locals=local_symbols), placeholders


def sympy_polynomial_row(poly: Any, t_symbol: Any, sympy: Any) -> list[str]:
    expanded = sympy.Poly(sympy.expand(poly), t_symbol, domain=sympy.QQ)
    if expanded.is_zero:
        return ["0"]
    return [rational_text(expanded.nth(index)) for index in range(expanded.degree() + 1)]


def trim_fraction_row(row: list[Fraction]) -> list[Fraction]:
    while len(row) > 1 and row[-1] == 0:
        row.pop()
    return row or [Fraction(0)]


def add_fraction_rows(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * max(len(left), len(right))
    for index, value in enumerate(left):
        out[index] += value
    for index, value in enumerate(right):
        out[index] += value
    return trim_fraction_row(out)


def multiply_fraction_rows(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim_fraction_row(out)


def derivative_fraction_row(row: list[Fraction], order: int) -> list[Fraction]:
    out = list(row)
    for _ in range(order):
        if len(out) <= 1:
            return [Fraction(0)]
        out = [Fraction(index + 1) * out[index + 1] for index in range(len(out) - 1)]
    return trim_fraction_row(out)


def evaluate_sparse(terms: list[tuple[int, int, str]], n: int) -> list[Fraction]:
    width = max((t_degree for _, t_degree, _ in terms), default=0) + 1
    out = [Fraction(0)] * width
    for n_degree, t_degree, value in terms:
        out[t_degree] += Fraction(value) * n**n_degree
    return trim_fraction_row(out)


def divide_fraction_rows_exact(
    numerator: list[Fraction], denominator: list[Fraction]
) -> list[Fraction]:
    remainder = trim_fraction_row(list(numerator))
    denominator = trim_fraction_row(list(denominator))
    if denominator == [0]:
        raise ValueError("zero recurrence denominator")
    if remainder == [0]:
        return [Fraction(0)]
    if len(remainder) < len(denominator):
        raise ValueError("non-polynomial recurrence quotient")
    quotient = [Fraction(0)] * (len(remainder) - len(denominator) + 1)
    while remainder != [0] and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        factor = remainder[-1] / denominator[-1]
        quotient[shift] += factor
        for index, value in enumerate(denominator):
            remainder[index + shift] -= factor * value
        remainder = trim_fraction_row(remainder)
    if remainder != [0]:
        raise ValueError("non-polynomial recurrence quotient")
    return trim_fraction_row(quotient)


def recurrence_matches_rows(
    rows: list[list[str]],
    initial_count: int,
    first_index: int,
    terms: list[dict[str, Any]],
    denominator: list[tuple[int, int, str]] | None,
    inhomogeneous: list[tuple[int, int, str]] | None,
) -> bool:
    generated = [trim_fraction_row([Fraction(value) for value in row]) for row in rows[:initial_count]]
    try:
        while len(generated) < len(rows):
            n = first_index + len(generated)
            value = [Fraction(0)]
            for term in terms:
                source = generated[len(generated) - int(term["offset"])]
                derivative = derivative_fraction_row(source, int(term["derivative"]))
                coefficient = evaluate_sparse(term["coefficient"], n)
                value = add_fraction_rows(value, multiply_fraction_rows(coefficient, derivative))
            if inhomogeneous is not None:
                value = add_fraction_rows(value, evaluate_sparse(inhomogeneous, n))
            if denominator is not None:
                value = divide_fraction_rows_exact(value, evaluate_sparse(denominator, n))
            generated.append(value)
    except (IndexError, ValueError, ZeroDivisionError):
        return False
    expected = [trim_fraction_row([Fraction(value) for value in row]) for row in rows]
    return generated == expected


def generate_sparse_rows(
    initial_rows: list[list[str]],
    first_index: int,
    total_rows: int,
    terms: list[dict[str, Any]],
    denominator: list[tuple[int, int, str]] | None,
    inhomogeneous: list[tuple[int, int, str]] | None,
) -> list[list[Fraction]]:
    generated = [trim_fraction_row([Fraction(value) for value in row]) for row in initial_rows]
    while len(generated) < total_rows:
        n = first_index + len(generated)
        value = [Fraction(0)]
        for term in terms:
            source = generated[len(generated) - int(term["offset"])]
            derivative = derivative_fraction_row(source, int(term["derivative"]))
            coefficient = evaluate_sparse(term["coefficient"], n)
            value = add_fraction_rows(value, multiply_fraction_rows(coefficient, derivative))
        if inhomogeneous is not None:
            value = add_fraction_rows(value, evaluate_sparse(inhomogeneous, n))
        if denominator is not None:
            value = divide_fraction_rows_exact(value, evaluate_sparse(denominator, n))
        generated.append(value)
    return generated[:total_rows]


def matched_prefix_length(left: list[int], right: list[int], position: int) -> int:
    checked = min(len(left), len(right) - position)
    if checked < 20 or left[:checked] != right[position : position + checked]:
        return 0
    return checked


def complete_rows_in_prefix(rows: list[list[int]], term_count: int) -> int:
    used = 0
    count = 0
    for row in rows:
        if used + len(row) > term_count:
            break
        used += len(row)
        count += 1
    return count


def lean_sequence_entries(
    lean_sequences: Path,
    oeis_data: Path,
    existing_ids: set[str],
) -> list[dict[str, Any]]:
    """Import canonical generated recurrences from real-rooted-oeis-proofs."""
    import sympy as sp

    n_symbol, t_symbol = sp.symbols("N T")
    entries: list[dict[str, Any]] = []
    paths = [
        path
        for path in sorted(lean_sequences.glob("A*.lean"))
        if re.fullmatch(r"A\d{6}\.lean", path.name)
    ]
    for path in paths:
        oeis_id = path.stem
        if oeis_id in existing_ids:
            continue
        source = path.read_text(encoding="utf-8")
        definition = re.search(
            rf"^def {oeis_id} : ℕ → ℝ\[X\]\n(.*?)(?=\n\nexample :)",
            source,
            re.MULTILINE | re.DOTALL,
        )
        if definition is None:
            raise ValueError(f"{oeis_id}: generated definition not found")

        clauses: list[tuple[str, Any, set[str]]] = []
        for line in definition.group(1).splitlines():
            left, right = line.strip().removeprefix("|").split("=>", 1)
            expression, placeholders = lean_generated_expression(right.strip(), oeis_id, sp)
            clauses.append((left.strip(), expression, placeholders))
        recurrence_match = re.fullmatch(r"n \+ (\d+)", clauses[-1][0])
        if recurrence_match is None:
            raise ValueError(f"{oeis_id}: unsupported recurrence clause {clauses[-1][0]}")
        initial_count = int(recurrence_match.group(1))
        if initial_count != len(clauses) - 1:
            raise ValueError(f"{oeis_id}: noncontiguous generated base clauses")

        initial_rows: list[list[str]] = []
        for index, (left, expression, placeholders) in enumerate(clauses[:-1]):
            if left != str(index) or placeholders or n_symbol in expression.free_symbols:
                raise ValueError(f"{oeis_id}: unsupported base clause {left}")
            initial_rows.append(sympy_polynomial_row(expression, t_symbol, sp))

        expression = clauses[-1][1]
        placeholders = sorted(clauses[-1][2])
        placeholder_symbols = [sp.Symbol(name) for name in placeholders]
        constant = sp.expand(
            expression.subs({symbol: 0 for symbol in placeholder_symbols})
        )
        coefficients: dict[str, Any] = {}
        reconstructed = constant
        denominators: list[Any] = []
        for name, symbol in zip(placeholders, placeholder_symbols):
            coefficient = sp.cancel(sp.diff(expression, symbol))
            if any(other in coefficient.free_symbols for other in placeholder_symbols):
                raise ValueError(f"{oeis_id}: nonlinear generated recurrence")
            coefficients[name] = coefficient
            reconstructed += coefficient * symbol
            denominators.append(sp.denom(coefficient))
        if sp.simplify(expression - reconstructed) != 0:
            raise ValueError(f"{oeis_id}: recurrence is not linear in earlier rows")
        if constant != 0:
            denominators.append(sp.denom(sp.cancel(constant)))

        common_denominator = sp.Integer(1)
        for denominator in denominators:
            common_denominator = sp.lcm(common_denominator, denominator)
        output_substitution = {n_symbol: n_symbol - initial_count}
        terms: list[dict[str, Any]] = []
        for name, coefficient in coefficients.items():
            parsed_name = re.fullmatch(r"Y(\d+)D(\d+)", name)
            assert parsed_name is not None
            source_shift = int(parsed_name.group(1))
            offset = initial_count - source_shift
            if offset <= 0:
                raise ValueError(f"{oeis_id}: recurrence refers to its output row")
            numerator = sp.cancel(coefficient * common_denominator).subs(output_substitution)
            terms.append(
                {
                    "offset": offset,
                    "derivative": int(parsed_name.group(2)),
                    "alternating": False,
                    "coefficient": sympy_sparse(numerator, n_symbol, t_symbol),
                }
            )
        terms.sort(key=lambda term: (term["offset"], term["derivative"]))
        denominator = sp.cancel(common_denominator).subs(output_substitution)
        sparse_denominator = sympy_sparse(denominator, n_symbol, t_symbol)
        if sparse_denominator == [(0, 0, "1")]:
            sparse_denominator = None
        inhomogeneous = None
        if constant != 0:
            inhomogeneous = sympy_sparse(
                sp.cancel(constant * common_denominator).subs(output_substitution),
                n_symbol,
                t_symbol,
            )

        generated = generate_sparse_rows(
            initial_rows,
            0,
            16,
            terms,
            sparse_denominator,
            inhomogeneous,
        )
        integer_rows: list[list[int]] = []
        for row in generated:
            if any(value.denominator != 1 for value in row):
                raise ValueError(f"{oeis_id}: nonintegral validation row")
            integer_rows.append([value.numerator for value in row])

        record = oeis_record(oeis_data, oeis_id)
        layout = "RegularTriangle" if "tabl" in record["keywords"] else "Table"
        prefix_rows: list[list[str]] = []
        bfile_verified = False
        matched_terms = 0
        recurrence_first_width = len(integer_rows[0])
        row_start = 0
        flattened_offset = record["offset"]
        compared_rows = integer_rows

        if layout == "RegularTriangle":
            prefix_count = 0
            while prefix_count * (prefix_count + 1) // 2 < len(record["data"]):
                padded_rows: list[list[int]] = []
                for index, row in enumerate(integer_rows):
                    width = prefix_count + index + 1
                    if len(row) > width:
                        padded_rows = []
                        break
                    padded_rows.append(row + [0] * (width - len(row)))
                position = prefix_count * (prefix_count + 1) // 2
                flattened = [value for row in padded_rows for value in row]
                checked = matched_prefix_length(flattened, record["data"], position)
                if checked:
                    bfile_verified = True
                    matched_terms = checked
                    recurrence_first_width = prefix_count + 1
                    row_start = record["offset"]
                    compared_rows = padded_rows
                    cursor = 0
                    for width in range(1, prefix_count + 1):
                        prefix_rows.append(
                            [str(value) for value in record["data"][cursor : cursor + width]]
                        )
                        cursor += width
                    break
                prefix_count += 1
        else:
            flattened = [value for row in integer_rows for value in row]
            for position in range(len(record["data"])):
                checked = matched_prefix_length(flattened, record["data"], position)
                if checked:
                    bfile_verified = True
                    matched_terms = checked
                    flattened_offset = record["offset"] + position
                    break

        matched_rows = complete_rows_in_prefix(compared_rows, matched_terms)
        validation_row_index = None
        validation_row = None
        if matched_rows:
            validation_row_index = len(prefix_rows) + matched_rows - 1
            validation_row = [str(value) for value in compared_rows[matched_rows - 1]]
        rows_hash = hashlib.sha256(
            ",".join(str(value) for value in record["data"]).encode("utf-8")
        ).hexdigest()
        entries.append(
            {
                "id": oeis_id,
                "name": record["name"],
                "status": "Validated" if bfile_verified else "Experimental",
                "layout": layout,
                "row_start": row_start,
                "flattened_offset": flattened_offset,
                "bfile_verified": bfile_verified,
                "fixture_slug": "",
                "fixture_row_offset": len(prefix_rows),
                "prefix_rows": prefix_rows,
                "recurrence_first_index": 0,
                "recurrence_first_width": recurrence_first_width,
                "initial_rows": initial_rows,
                "terms": terms,
                "denominator": sparse_denominator,
                "inhomogeneous": inhomogeneous,
                "source_rows": matched_rows,
                "verification_rows": 0,
                "rows_sha256": rows_hash,
                "validation_row_index": validation_row_index,
                "validation_row": validation_row,
            }
        )
    return entries


def load_rows_jsonl(path: Path) -> list[tuple[int, list[str]]]:
    rows: list[tuple[int, list[str]]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        row = json.loads(line)
        rows.append((int(row["n"]), [str(int(value)) for value in row["coefficients"]]))
    rows.sort()
    if not rows or any(rows[index + 1][0] != rows[index][0] + 1 for index in range(len(rows) - 1)):
        raise ValueError(f"noncontiguous row data in {path}")
    return rows


def queue_entries(
    sequence_library: Path,
    oeis_data: Path,
    existing_ids: set[str],
) -> list[dict[str, Any]]:
    import sympy as sp

    tools = sequence_library / "_tools"
    sys.path.insert(0, str(tools))
    try:
        from recurrence_generator import N, T, _parse_recurrence
    finally:
        sys.path.pop(0)

    entries: list[dict[str, Any]] = []
    for sequence_dir in sorted(sequence_library.glob("A[0-9][0-9][0-9][0-9][0-9][0-9]")):
        oeis_id = sequence_dir.name
        if oeis_id in existing_ids or (sequence_dir / "tags" / "invalid_recurrence").exists():
            continue
        manifest = json.loads((sequence_dir / "manifest.json").read_text(encoding="utf-8"))
        recurrence_data = json.loads(
            (sequence_dir / "recurrence.json").read_text(encoding="utf-8")
        )
        raw_rows_text = (sequence_dir / "rows.jsonl").read_text(encoding="utf-8")
        indexed_rows = load_rows_jsonl(sequence_dir / "rows.jsonl")
        rows = [row for _, row in indexed_rows]
        expression, symbols, max_lag = _parse_recurrence(recurrence_data["text"])

        current = symbols.get("Y_0_0")
        if current is None:
            raise ValueError(f"{oeis_id}: recurrence has no current-row term")
        denominator_expression = sp.expand(expression.coeff(current))
        if denominator_expression == 0:
            raise ValueError(f"{oeis_id}: zero current-row coefficient")

        terms: list[dict[str, Any]] = []
        symbol_values = list(symbols.values())
        for name, symbol in symbols.items():
            _, lag_text, derivative_text = name.split("_")
            lag = int(lag_text)
            derivative = int(derivative_text)
            if lag == 0:
                continue
            coefficient = sp.expand(-expression.coeff(symbol))
            if coefficient != 0:
                terms.append(
                    {
                        "offset": lag,
                        "derivative": derivative,
                        "alternating": False,
                        "coefficient": sympy_sparse(coefficient, N, T),
                    }
                )
        constant = sp.expand(-expression.subs({symbol: 0 for symbol in symbol_values}))
        denominator = sympy_sparse(denominator_expression, N, T)
        if denominator == [(0, 0, "1")]:
            denominator = None
        inhomogeneous = sympy_sparse(constant, N, T) if constant != 0 else None

        index_offset = int(recurrence_data.get("recurrence_index_offset", 1))
        initial_count = max_lag
        for row_index in range(max_lag, len(rows)):
            evaluated = sp.expand(denominator_expression.subs(N, row_index + index_offset))
            if evaluated == 0:
                initial_count = row_index + 1
        if initial_count <= 0 or initial_count > len(rows):
            raise ValueError(f"{oeis_id}: invalid initial row count {initial_count}")
        if not recurrence_matches_rows(
            rows,
            initial_count,
            index_offset,
            terms,
            denominator,
            inhomogeneous,
        ):
            print(
                f"skipping {oeis_id}: stored recurrence does not reproduce cached rows",
                file=sys.stderr,
            )
            continue

        record = oeis_record(oeis_data, oeis_id)
        missing_prefix, bfile_verified = find_prefix_alignment(rows, record)
        if bfile_verified:
            first_row = record["offset"]
        else:
            first_row = int(manifest["row_convention"]["start_n"])
        holdout_tag = (sequence_dir / "tags" / "recurrence_verified_holdout").exists()
        rows_used = int(recurrence_data.get("rows_used_to_find") or 0)
        verification_rows = max(0, len(rows) - rows_used) if holdout_tag else 0
        entries.append(
            {
                "id": oeis_id,
                "name": record["name"],
                "status": "Verified" if holdout_tag else "Experimental",
                "layout": "RegularTriangle" if "tabl" in record["keywords"] else "Table",
                "row_start": first_row,
                "flattened_offset": record["offset"],
                "bfile_verified": bfile_verified,
                "fixture_slug": "",
                "fixture_row_offset": len(missing_prefix),
                "prefix_rows": missing_prefix,
                "recurrence_first_index": index_offset,
                "recurrence_first_width": len(rows[0]),
                "initial_rows": rows[:initial_count],
                "terms": terms,
                "denominator": denominator,
                "inhomogeneous": inhomogeneous,
                "source_rows": len(rows),
                "verification_rows": verification_rows,
                "rows_sha256": hashlib.sha256(raw_rows_text.encode("utf-8")).hexdigest(),
            }
        )
    return entries


def find_prefix_alignment(rows: list[list[str]], record: dict[str, Any]) -> tuple[list[list[str]], bool]:
    flat = [int(value) for row in rows for value in row]
    data = record["data"]
    positions: list[tuple[int, int]] = []
    for position in range(len(data)):
        checked = min(len(flat), len(data) - position)
        if checked >= 20 and flat[:checked] == data[position : position + checked]:
            positions.append((position, checked))
    if not positions:
        return [], False
    position, _ = min(positions)

    if "tabl" in record["keywords"]:
        discriminant = 1 + 8 * position
        root = math.isqrt(discriminant)
        prefix_rows = (root - 1) // 2
        if root * root != discriminant or prefix_rows * (prefix_rows + 1) // 2 != position:
            return [], False
        if any(len(row) != prefix_rows + index + 1 for index, row in enumerate(rows[:10])):
            return [], False
        prefix: list[list[str]] = []
        cursor = 0
        for width in range(1, prefix_rows + 1):
            prefix.append([str(value) for value in data[cursor : cursor + width]])
            cursor += width
        return prefix, True

    # A046739 has the exceptional initial data 0; 1; 1,1 before the regular
    # rows 1,7,1; 1,21,21,1; ... .  Keep this explicit and prefix-checked.
    if record["name"].startswith("Triangle read by rows, related to number of permutations"):
        if position == 4 and data[:4] == [0, 1, 1, 1]:
            return [["0"], ["1"], ["1", "1"]], True
    return [], False


def render_sparse(terms: list[tuple[int, int, str]] | None, indent: str) -> str:
    if terms is None:
        return "None"
    if not terms:
        return "Some(&[])"
    rendered = ",\n".join(
        f'{indent}    SparseCoefficient::new({n}, {t}, {rust_string(value)})'
        for n, t, value in terms
    )
    return f"Some(&[\n{rendered},\n{indent}])"


def render_rows(rows: list[list[str]], indent: str) -> str:
    if not rows:
        return "&[]"
    rendered_rows = []
    for row in rows:
        values = ", ".join(rust_string(value) for value in row)
        rendered_rows.append(f"{indent}    &[{values}]")
    return "&[\n" + ",\n".join(rendered_rows) + f",\n{indent}]"


def render_catalog(
    fixtures: Path,
    oeis_data: Path,
    sequence_library: Path,
    lean_sequences: Path,
) -> str:
    if not sequence_library.is_dir():
        raise SystemExit(f"OEIS sequence library not found: {sequence_library}")
    if not oeis_data.is_dir():
        raise SystemExit(f"OEIS data directory not found: {oeis_data}")
    if not lean_sequences.is_dir():
        raise SystemExit(f"Lean OEIS sequence directory not found: {lean_sequences}")
    with (fixtures / "manifest.tsv").open(encoding="utf-8", newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))

    entries: list[dict[str, Any]] = []
    for item in manifest:
        match = re.search(r"oeis_(a\d{6})", item["slug"], re.IGNORECASE)
        if match is None:
            continue
        oeis_id = match.group(1).upper()
        recurrence_json = json.loads((fixtures / item["json_file"]).read_text(encoding="utf-8"))
        rows_text = (fixtures / item["rows_file"]).read_text(encoding="utf-8")
        rows = [parse_integer_row(line) for line in rows_text.splitlines() if line.strip()]
        record = oeis_record(oeis_data, oeis_id)
        missing_prefix, bfile_verified = find_prefix_alignment(rows, record)
        skipped = int((recurrence_json.get("search") or {}).get("skip_prefix") or 0)
        prefix_rows = missing_prefix + rows[:skipped]
        recurrence = recurrence_json["recurrence"]
        entries.append(
            {
                "id": oeis_id,
                "name": record["name"],
                "status": "Verified",
                "layout": "RegularTriangle" if "tabl" in record["keywords"] else "Table",
                "row_start": record["offset"],
                "flattened_offset": record["offset"],
                "bfile_verified": bfile_verified,
                "fixture_slug": item["slug"],
                "fixture_row_offset": len(missing_prefix),
                "prefix_rows": prefix_rows,
                "recurrence_first_index": recurrence_json["first_index"],
                "recurrence_first_width": len(rows[skipped]),
                "initial_rows": recurrence_json["initial_polynomials"],
                "terms": [
                    {
                        "offset": term["offset"],
                        "derivative": term["deriv_order"],
                        "alternating": term["sign"] == "alternating_n",
                        "coefficient": sparse_coefficients(term["coeff"]),
                    }
                    for term in recurrence["terms"]
                ],
                "denominator": sparse_coefficients(recurrence.get("denominator")),
                "inhomogeneous": sparse_coefficients(recurrence.get("inhomogeneous")),
                "source_rows": int(recurrence_json["search"]["source_rows"]),
                "verification_rows": int(recurrence_json["search"]["verification_polynomials"]),
                "rows_sha256": hashlib.sha256(rows_text.encode("utf-8")).hexdigest(),
            }
        )

    entries.extend(queue_entries(sequence_library, oeis_data, {entry["id"] for entry in entries}))
    entries.extend(
        lean_sequence_entries(
            lean_sequences,
            oeis_data,
            {entry["id"] for entry in entries},
        )
    )
    entries.sort(key=lambda entry: entry["id"])
    lines = [
        "// @generated by scripts/build_oeis_catalog.py; do not edit by hand.",
        "",
        "pub static OEIS_CATALOG: &[OeisSequenceDefinition] = &[",
    ]
    for entry in entries:
        lines.extend(
            [
                "    OeisSequenceDefinition {",
                f'        id: {rust_string(entry["id"])},',
                f'        name: {rust_string(entry["name"])},',
                f'        status: OeisSequenceStatus::{entry["status"]},',
                f'        layout: OeisLayout::{entry["layout"]},',
                f'        first_row: {entry["row_start"]},',
                f'        flattened_offset: {entry["flattened_offset"]},',
                f'        bfile_prefix_verified: {str(entry["bfile_verified"]).lower()},',
                f'        fixture_slug: {rust_string(entry["fixture_slug"])},',
                f'        fixture_row_offset: {entry["fixture_row_offset"]},',
                f'        prefix_rows: {render_rows(entry["prefix_rows"], "        ")},',
                f'        recurrence_first_index: {entry["recurrence_first_index"]},',
                f'        recurrence_first_width: {entry["recurrence_first_width"]},',
                f'        initial_rows: {render_rows(entry["initial_rows"], "        ")},',
                "        recurrence: SparseRecurrence {",
                "            terms: &[",
            ]
        )
        for term in entry["terms"]:
            lines.extend(
                [
                    "                SparseRecurrenceTerm {",
                    f'                    offset: {term["offset"]},',
                    f'                    derivative_order: {term["derivative"]},',
                    f'                    alternating_sign: {str(term["alternating"]).lower()},',
                    f'                    coefficient: {render_sparse(term["coefficient"], "                    ").replace("Some(", "").removesuffix(")")},',
                    "                },",
                ]
            )
        lines.extend(
            [
                "            ],",
                f'            denominator: {render_sparse(entry["denominator"], "            ")},',
                f'            inhomogeneous: {render_sparse(entry["inhomogeneous"], "            ")},',
                "        },",
                f'        source_rows: {entry["source_rows"]},',
                f'        verification_rows: {entry["verification_rows"]},',
                f'        rows_sha256: {rust_string(entry["rows_sha256"])},',
                "    },",
            ]
        )
    lines.append("];\n")
    lines.extend(
        [
            "#[cfg(test)]",
            "pub static OEIS_IMPORTED_VALIDATION_ROWS: &[(&str, usize, &[&str])] = &[",
        ]
    )
    for entry in entries:
        validation_row = entry.get("validation_row")
        if validation_row is None:
            continue
        values = ", ".join(rust_string(value) for value in validation_row)
        lines.append(
            f'    ({rust_string(entry["id"])}, {entry["validation_row_index"]}, &[{values}]),'
        )
    lines.extend(["];", ""])
    for index, entry in enumerate(entries):
        lines.extend(
            [
                "#[allow(non_snake_case)]",
                f'pub fn {entry["id"]}() -> &\'static OeisSequenceDefinition {{',
                f"    &OEIS_CATALOG[{index}]",
                "}",
                "",
            ]
        )
    return "\n".join(lines)


def rustfmt(source: str) -> str:
    completed = subprocess.run(
        ["rustfmt", "--edition", "2021"],
        input=source,
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise SystemExit(f"rustfmt failed:\n{completed.stderr}")
    return completed.stdout


def main() -> int:
    args = parse_args()
    rendered = rustfmt(
        render_catalog(
            args.fixtures,
            args.oeis_data,
            args.sequence_library,
            args.lean_sequences,
        )
    )
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"generated catalog is stale: {args.output}")
        return 0
    args.output.write_text(rendered, encoding="utf-8")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
