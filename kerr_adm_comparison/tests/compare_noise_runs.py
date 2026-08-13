#!/usr/bin/env python3
"""Compare deterministic geometry and error tables from two test runs."""

from __future__ import annotations

import argparse
from pathlib import Path
import re


SECTION = re.compile(r"^(HERMITE|LAGRANGE\+FD|ANALYTIC-FD) / (VALUE|D/DX|D/DY|D/DZ)$")
FIELDS = {
    "alpha",
    "betax",
    "betay",
    "betaz",
    "gammaxx",
    "gammaxy",
    "gammaxz",
    "gammayy",
    "gammayz",
    "gammazz",
    "ALL",
}


def parse(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    query = next(line.strip() for line in lines if "random exterior centers =" in line)

    grid = []
    in_grid = False
    errors = {}
    section = None
    for line in lines:
        if line == "Grid diagnostics":
            in_grid = True
            continue
        if line == "Serial interpolation latency":
            in_grid = False
        if in_grid:
            tokens = line.split()
            if len(tokens) == 8 and tokens[0].isdigit():
                # The last column is the OpenMP team size and may differ.
                grid.append(tuple(tokens[:-1]))

        match = SECTION.match(line)
        if match:
            section = match.groups()
            continue
        tokens = line.split()
        if section is not None and len(tokens) == 9 and tokens[0] in FIELDS:
            errors[(section, tokens[0], tokens[1])] = tuple(tokens[2:])

    if not grid or not errors:
        raise ValueError(f"could not parse comparison data from {path}")
    return query, tuple(grid), errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--errors", choices=("equal", "different"), required=True)
    arguments = parser.parse_args()

    left_query, left_grid, left_errors = parse(arguments.left)
    right_query, right_grid, right_errors = parse(arguments.right)
    if left_query != right_query:
        raise SystemExit("query count or query seed changed")
    if left_grid != right_grid:
        raise SystemExit("grid geometry, query radii, or crossing counts changed")
    if left_errors.keys() != right_errors.keys():
        raise SystemExit("reported error channels differ")

    equal = left_errors == right_errors
    if arguments.errors == "equal" and not equal:
        raise SystemExit("deterministic error tables differ")
    if arguments.errors == "different" and equal:
        raise SystemExit("changing the noise seed did not change any error")
    print(
        f"run comparison passed: geometry identical, errors {arguments.errors}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
