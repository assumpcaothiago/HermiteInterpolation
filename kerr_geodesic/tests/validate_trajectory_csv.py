#!/usr/bin/env python3
"""Validate deterministic structural properties of a geodesic trajectory CSV."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


METHODS = {"EXACT", "HERMITE", "LAGRANGE+FD"}
REQUIRED = {
    "resolution",
    "method",
    "t",
    "x",
    "y",
    "z",
    "px",
    "py",
    "pz",
    "position_error",
    "momentum_error",
    "mass_shell_residual",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--resolutions", required=True)
    arguments = parser.parse_args()
    expected = {int(value) for value in arguments.resolutions.split(",")}
    grouped: dict[tuple[int, str], list[float]] = {}
    with arguments.csv.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not REQUIRED <= set(reader.fieldnames):
            raise ValueError("missing required trajectory columns")
        for line_number, row in enumerate(reader, start=2):
            try:
                resolution = int(row["resolution"])
                method = row["method"]
                numeric = [
                    float(value)
                    for name, value in row.items()
                    if name not in {"resolution", "method"}
                ]
                time = float(row["t"])
                position_error = float(row["position_error"])
                momentum_error = float(row["momentum_error"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid row {line_number}") from error
            if resolution not in expected or method not in METHODS:
                raise ValueError(f"unexpected selector on row {line_number}")
            if not all(math.isfinite(value) for value in numeric):
                raise ValueError(f"nonfinite value on row {line_number}")
            if method == "EXACT" and (position_error != 0.0 or momentum_error != 0.0):
                raise ValueError("exact rows must coincide with the numerical reference")
            grouped.setdefault((resolution, method), []).append(time)
    if set(grouped) != {(resolution, method) for resolution in expected for method in METHODS}:
        raise ValueError("missing resolution/method group")
    lengths = {len(times) for times in grouped.values()}
    if len(lengths) != 1 or next(iter(lengths)) < 2:
        raise ValueError("trajectory groups do not have equal nontrivial lengths")
    final_times = set()
    for times in grouped.values():
        if times[0] != 0.0 or any(right <= left for left, right in zip(times, times[1:])):
            raise ValueError("trajectory times must start at zero and increase")
        final_times.add(times[-1])
    if len(final_times) != 1 or not (next(iter(final_times)) > 0.0):
        raise ValueError("trajectory groups have inconsistent final times")
    print(
        f"trajectory CSV: {len(grouped)} groups, {next(iter(lengths))} rows each"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
