#!/usr/bin/env python3
"""Plot an analytic and Hermite-interpolated Kerr z-axis profile from CSV."""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import tempfile


REQUIRED_COLUMNS = {
    "component",
    "quantity",
    "resolution",
    "z",
    "analytic",
    "interpolated",
    "error",
}
THROAT = (2.0 + math.sqrt(3.0)) / 8.0


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="profile CSV from kerr_convergence")
    parser.add_argument(
        "--output",
        type=Path,
        help="image path (default: replace the CSV suffix with .png)",
    )
    parser.add_argument(
        "--show", action="store_true", help="also open an interactive plot window"
    )
    return parser.parse_args()


def read_profile(
    path: Path,
) -> tuple[str, str, dict[int, list[tuple[float, float, float, float]]]]:
    """Read one component/quantity profile grouped by grid resolution."""

    grouped: dict[int, list[tuple[float, float, float, float]]] = {}
    selectors: set[tuple[str, str]] = set()
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not REQUIRED_COLUMNS <= set(reader.fieldnames):
            raise ValueError(
                "profile CSV does not contain the expected columns: "
                + ", ".join(sorted(REQUIRED_COLUMNS))
            )
        for line_number, row in enumerate(reader, start=2):
            try:
                component = row["component"]
                quantity = row["quantity"]
                resolution = int(row["resolution"])
                values = (
                    float(row["z"]),
                    float(row["analytic"]),
                    float(row["interpolated"]),
                    float(row["error"]),
                )
            except (TypeError, ValueError) as error:
                raise ValueError(f"invalid data on CSV line {line_number}") from error
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"nonfinite data on CSV line {line_number}")
            if abs(values[0]) < THROAT:
                raise ValueError("the z-axis profile must not sample the region |z|<s")
            selectors.add((component, quantity))
            grouped.setdefault(resolution, []).append(values)

    if not grouped:
        raise ValueError("profile CSV contains no data rows")
    if len(selectors) != 1:
        raise ValueError("profile CSV must contain exactly one component and quantity")
    for rows in grouped.values():
        rows.sort(key=lambda row: row[0])
    component, quantity = selectors.pop()
    return component, quantity, grouped


def branches(
    rows: list[tuple[float, float, float, float]],
) -> tuple[
    list[tuple[float, float, float, float]],
    list[tuple[float, float, float, float]],
]:
    """Split the two exterior branches across the excluded region |z|<s."""

    return ([row for row in rows if row[0] < 0.0], [row for row in rows if row[0] > 0.0])


def variable_label(component: str, quantity: str) -> str:
    metric = rf"g_{{{component}}}"
    if quantity == "value":
        return rf"${metric}$"
    return rf"$\partial_{{{quantity[1:]}}}{metric}$"


def main() -> int:
    arguments = parse_arguments()
    output = arguments.output or arguments.csv.with_suffix(".png")
    component, quantity, grouped = read_profile(arguments.csv)

    # Matplotlib otherwise tries to create a cache in the user's home folder.
    # A temporary cache keeps this optional plotting tool self-contained.
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "kerr-profile-matplotlib")
    )
    import matplotlib

    if not arguments.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    resolutions = sorted(grouped)
    figure, (value_axis, error_axis) = plt.subplots(
        2,
        1,
        figsize=(8.0, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0)},
    )

    # The analytic values are identical at every resolution, so plot one copy.
    first_rows = grouped[resolutions[0]]
    for branch_index, branch in enumerate(branches(first_rows)):
        value_axis.plot(
            [row[0] for row in branch],
            [row[1] for row in branch],
            color="black",
            linewidth=2.0,
            label="analytic" if branch_index == 0 else None,
        )

    color_map = plt.get_cmap("viridis")
    any_positive_error = False
    for color_index, resolution in enumerate(resolutions):
        color_position = (
            0.5 if len(resolutions) == 1 else color_index / (len(resolutions) - 1)
        )
        color = color_map(color_position)
        for branch_index, branch in enumerate(branches(grouped[resolution])):
            value_axis.plot(
                [row[0] for row in branch],
                [row[2] for row in branch],
                color=color,
                linewidth=1.25,
                label=f"Hermite N={resolution}" if branch_index == 0 else None,
            )
            errors = [abs(row[3]) for row in branch]
            any_positive_error = any_positive_error or any(error > 0.0 for error in errors)
            error_axis.plot(
                [row[0] for row in branch],
                errors,
                color=color,
                linewidth=1.25,
                label=f"N={resolution}" if branch_index == 0 else None,
            )

    for axis in (value_axis, error_axis):
        axis.axvline(-THROAT, color="0.55", linestyle=":", linewidth=1.0)
        axis.axvline(THROAT, color="0.55", linestyle=":", linewidth=1.0)
        axis.axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
        axis.grid(True, alpha=0.25)

    label = variable_label(component, quantity)
    value_axis.set_title(f"Kerr z-axis profile: {label}")
    value_axis.set_ylabel(label)
    value_axis.legend(loc="best")
    error_axis.set_xlabel("z")
    error_axis.set_ylabel(r"$|\mathrm{interpolated}-\mathrm{analytic}|$")
    if any_positive_error:
        error_axis.set_yscale("log")
    else:
        error_axis.text(
            0.5,
            0.5,
            "zero interpolation error",
            ha="center",
            va="center",
            transform=error_axis.transAxes,
        )

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
