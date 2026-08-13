#!/usr/bin/env python3
"""Exercise trajectory plot-data loading without importing Matplotlib."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE))

from plot_trajectory import METHODS, read_rows  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    arguments = parser.parse_args()
    grouped = read_rows(arguments.csv)
    if not grouped:
        raise ValueError("no resolution groups")
    for resolution, methods in grouped.items():
        if set(methods) != set(METHODS):
            raise ValueError(f"resolution {resolution} has incomplete methods")
        lengths = {len(rows) for rows in methods.values()}
        if len(lengths) != 1 or next(iter(lengths)) < 2:
            raise ValueError(f"resolution {resolution} has inconsistent rows")
    print(f"plot data: {len(grouped)} resolution groups passed")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
