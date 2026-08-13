#!/usr/bin/env python3
"""Plot trajectories and invariant drift from kerr_geodesic CSV output."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile


REQUIRED_COLUMNS = {
    "resolution",
    "method",
    "t",
    "x",
    "y",
    "radial_error",
    "phase_error",
    "exact_energy_drift",
    "lz_drift",
    "closed_x",
    "closed_y",
}
METHODS = ("EXACT", "HERMITE", "LAGRANGE+FD")


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV written by kerr_geodesic")
    parser.add_argument(
        "--resolution",
        type=int,
        help="resolution to plot (default: finest present in the CSV)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="image path (default: replace the CSV suffix with .png)",
    )
    parser.add_argument("--show", action="store_true", help="show the plot")
    return parser.parse_args()


def read_rows(path: Path) -> dict[int, dict[str, list[dict[str, float]]]]:
    grouped: dict[int, dict[str, list[dict[str, float]]]] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not REQUIRED_COLUMNS <= set(reader.fieldnames):
            raise ValueError(
                "trajectory CSV is missing columns: "
                + ", ".join(sorted(REQUIRED_COLUMNS - set(reader.fieldnames or ())))
            )
        for line_number, raw in enumerate(reader, start=2):
            try:
                resolution = int(raw["resolution"])
                method = raw["method"]
                row = {
                    name: float(raw[name])
                    for name in REQUIRED_COLUMNS
                    if name not in {"resolution", "method"}
                }
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid data on CSV line {line_number}") from error
            if method not in METHODS:
                raise ValueError(f"unknown method on CSV line {line_number}: {method}")
            if not all(math.isfinite(value) for value in row.values()):
                raise ValueError(f"nonfinite data on CSV line {line_number}")
            grouped.setdefault(resolution, {}).setdefault(method, []).append(row)
    if not grouped:
        raise ValueError("trajectory CSV contains no rows")
    for methods in grouped.values():
        if set(methods) != set(METHODS):
            raise ValueError("every resolution must contain all three methods")
        for rows in methods.values():
            rows.sort(key=lambda row: row["t"])
    return grouped


def main() -> int:
    arguments = parse_arguments()
    grouped = read_rows(arguments.csv)
    resolution = arguments.resolution or max(grouped)
    if resolution not in grouped:
        raise ValueError(f"resolution {resolution} is not present in the CSV")
    output = arguments.output or arguments.csv.with_suffix(".png")

    # Keep Matplotlib's cache outside the repository and the user's home tree.
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "kerr-geodesic-mpl")
    )
    import matplotlib

    if not arguments.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"EXACT": "black", "HERMITE": "tab:blue", "LAGRANGE+FD": "tab:orange"}
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))
    orbit_axis, radius_axis, phase_axis, invariant_axis = axes.flat
    exact_rows = grouped[resolution]["EXACT"]
    orbit_axis.plot(
        [row["closed_x"] for row in exact_rows],
        [row["closed_y"] for row in exact_rows],
        color="0.65",
        linestyle="--",
        linewidth=1.5,
        label="closed circular orbit",
    )
    for method in METHODS:
        rows = grouped[resolution][method]
        color = colors[method]
        times = [row["t"] for row in rows]
        orbit_axis.plot(
            [row["x"] for row in rows],
            [row["y"] for row in rows],
            color=color,
            linewidth=1.1,
            label=method,
        )
        radius_axis.plot(
            times,
            [row["radial_error"] for row in rows],
            color=color,
            linewidth=1.0,
            label=method,
        )
        phase_axis.plot(
            times,
            [row["phase_error"] for row in rows],
            color=color,
            linewidth=1.0,
            label=method,
        )
        invariant_axis.plot(
            times,
            [abs(row["exact_energy_drift"]) for row in rows],
            color=color,
            linewidth=1.0,
            label=f"{method} |dE|",
        )
        invariant_axis.plot(
            times,
            [abs(row["lz_drift"]) for row in rows],
            color=color,
            linestyle=":",
            linewidth=1.0,
            label=f"{method} |dLz|",
        )

    orbit_axis.set_aspect("equal", adjustable="box")
    orbit_axis.set_xlabel("x")
    orbit_axis.set_ylabel("y")
    orbit_axis.set_title(f"Equatorial trajectory, N={resolution}")
    radius_axis.set_xlabel("t")
    radius_axis.set_ylabel("r - r0")
    radius_axis.set_title("Radial drift")
    phase_axis.set_xlabel("t")
    phase_axis.set_ylabel("phi - Omega t")
    phase_axis.set_title("Unwrapped phase error")
    invariant_axis.set_xlabel("t")
    invariant_axis.set_ylabel("absolute drift")
    invariant_axis.set_title("Exact Kerr invariants")
    if any(
        row[quantity] != 0.0
        for method in METHODS
        for row in grouped[resolution][method]
        for quantity in ("exact_energy_drift", "lz_drift")
    ):
        invariant_axis.set_yscale("log")
    for axis in axes.flat:
        axis.grid(True, alpha=0.25)
        axis.legend(loc="best", fontsize="small")
    figure.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=160)
    print(f"wrote {output}")
    if arguments.show:
        plt.show()
    plt.close(figure)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        raise SystemExit(f"error: {error}") from error
