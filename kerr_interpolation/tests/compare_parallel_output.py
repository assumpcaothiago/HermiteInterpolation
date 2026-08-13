#!/usr/bin/env python3
"""Compare deterministic one-thread and multi-thread convergence output."""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path


GRID = re.compile(
    r"^\s*(\d+)\s+(\d+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+"
    r"([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$"
)
NORM = re.compile(
    r"^(tt|tx|ty|tz|xx|xy|xz|yy|yz|zz|ALL)\s+(\d+)\s+"
    r"([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)"
)
SECTIONS = {"VALUE", "D/DX", "D/DY", "D/DZ"}


def parse(
    path: Path,
) -> tuple[
    list[str],
    dict[int, tuple[str, ...]],
    dict[tuple[str, str, int], tuple[float, float]],
]:
    metadata: list[str] = []
    grids: dict[int, tuple[str, ...]] = {}
    norms: dict[tuple[str, str, int], tuple[float, float]] = {}
    section = ""

    for line in path.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped in SECTIONS:
            section = stripped
            continue
        if stripped.startswith(
            (
                "M =",
                "active cube =",
                "side length =",
                "random points =",
                "random cloud is",
            )
        ):
            metadata.append(stripped)
            continue
        match = GRID.match(line)
        if match:
            resolution = int(match.group(1))
            grids[resolution] = match.groups()[1:]
            continue
        match = NORM.match(line)
        if match and section:
            name, resolution_text, _spacing, rms, maximum = match.groups()
            norms[(section, name, int(resolution_text))] = (
                float(rms),
                float(maximum),
            )

    if not metadata or not grids or not norms:
        raise ValueError(f"could not parse complete convergence output from {path}")
    return metadata, grids, norms


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("serial", type=Path)
    parser.add_argument("parallel", type=Path)
    parser.add_argument("--parallel-threads", type=int, required=True)
    arguments = parser.parse_args()

    serial_metadata, serial_grids, serial_norms = parse(arguments.serial)
    parallel_metadata, parallel_grids, parallel_norms = parse(
        arguments.parallel
    )
    if serial_metadata != parallel_metadata:
        raise ValueError("serial and parallel run metadata differ")
    if serial_grids.keys() != parallel_grids.keys():
        raise ValueError("serial and parallel grid resolutions differ")

    for resolution in serial_grids:
        serial = serial_grids[resolution]
        parallel = parallel_grids[resolution]
        if serial[:5] != parallel[:5]:
            raise ValueError(f"grid geometry differs at N={resolution}")
        if serial[5:] != ("1", "1"):
            raise ValueError(f"serial run did not use one thread at N={resolution}")
        expected = str(arguments.parallel_threads)
        if parallel[5:] != (expected, expected):
            raise ValueError(
                f"parallel run did not use {expected} threads at N={resolution}"
            )

    if serial_norms.keys() != parallel_norms.keys():
        raise ValueError("serial and parallel norm channels differ")
    for channel, serial in serial_norms.items():
        parallel = parallel_norms[channel]
        for label, left, right in zip(("L2", "Linf"), serial, parallel):
            if not math.isclose(left, right, rel_tol=5.0e-12, abs_tol=1.0e-30):
                raise ValueError(
                    f"{channel} {label} differs: serial={left}, parallel={right}"
                )
    for section in SECTIONS:
        for resolution in serial_grids:
            if serial_norms[(section, "tz", resolution)] != (0.0, 0.0):
                raise ValueError(
                    f"exact-zero channel {section}/tz is nonzero at N={resolution}"
                )

    print("Serial and OpenMP convergence outputs agree.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
