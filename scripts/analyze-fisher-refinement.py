#!/usr/bin/env python3

"""Summarize spatial refinement of the 83-region metric-graph FEM test."""

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt


FIELDS = ("global", "temporal", "frontal", "parietal", "occipital")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_table(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {key: [float(row[key]) for row in rows]
            for key in ("time", "min", "max", *FIELDS)}


def crossing(time, values, threshold=50.0):
    for index in range(1, len(time)):
        if values[index - 1] < threshold <= values[index]:
            fraction = ((threshold - values[index - 1])
                        / (values[index] - values[index - 1]))
            return time[index - 1] + fraction * (time[index] - time[index - 1])
    return float("nan")


def main():
    args = arguments()
    pattern = re.compile(r"cells_(\d+)_dt_([0-9]+)p([0-9]+)$")
    cases = []
    for directory in args.root.iterdir():
        match = pattern.match(directory.name)
        table = directory / "fem_biomarkers.csv"
        if match and table.exists():
            cells = int(match.group(1))
            dt = float(f"{match.group(2)}.{match.group(3)}")
            cases.append((cells, dt, read_table(table)))
    if not cases:
        raise RuntimeError("No refinement output was found")
    cases.sort(key=lambda item: (item[1], item[0]))
    spatial = [case for case in cases if abs(case[1] - 0.4) < 1.0e-12]
    reference = max(spatial, key=lambda item: item[0])[2]

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "space_refinement.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ["cells_per_edge", "dofs", "dt", "transient_min",
                  "transient_max", "global_t50", "max_biomarker_difference"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for cells, dt, data in spatial:
            difference = max(
                abs(value - reference[key][index])
                for key in FIELDS
                for index, value in enumerate(data[key])
            )
            writer.writerow({
                "cells_per_edge": cells,
                "dofs": 83 + (cells - 1) * 1130,
                "dt": dt,
                "transient_min": f"{min(data['min']):.12g}",
                "transient_max": f"{max(data['max']):.12g}",
                "global_t50": f"{crossing(data['time'], data['global']):.12g}",
                "max_biomarker_difference": f"{difference:.12g}",
            })

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for cells, _, data in spatial:
        axes[0].plot(data["time"], data["global"], linewidth=1.8,
                     label=f"{cells} cell{'s' if cells != 1 else ''}/edge")
        axes[1].plot(data["time"], data["min"], linewidth=1.8,
                     label=f"{cells} cell{'s' if cells != 1 else ''}/edge")
    axes[0].set_xlabel("time (years)")
    axes[0].set_ylabel("global biomarker (%)")
    axes[0].set_ylim(0.0, 102.0)
    axes[1].set_xlabel("time (years)")
    axes[1].set_ylabel("minimum nodal concentration")
    for axis in axes:
        axis.grid(True, linewidth=0.5, alpha=0.35)
        axis.legend()
    figure.tight_layout()
    figure.savefig(args.output_dir / "space_refinement.png", dpi=220)
    figure.savefig(args.output_dir / "space_refinement.pdf")
    plt.close(figure)

    temporal = [case for case in cases if case[0] == 8]
    temporal.sort(key=lambda item: item[1], reverse=True)
    def difference_at_common_times(coarse, fine):
        fine_by_time = {
            round(time, 10): index for index, time in enumerate(fine["time"])
        }
        return max(
            abs(coarse[key][index] - fine[key][fine_by_time[round(time, 10)]])
            for key in FIELDS
            for index, time in enumerate(coarse["time"])
        )

    time_reference = temporal[-1][2]
    differences_to_reference = [
        difference_at_common_times(data, time_reference)
        for _, _, data in temporal
    ]
    successive = [
        difference_at_common_times(temporal[index][2], temporal[index + 1][2])
        if index + 1 < len(temporal) else None
        for index in range(len(temporal))
    ]
    rates = [None] * len(temporal)
    for index in range(1, len(temporal) - 1):
        rates[index] = math.log2(successive[index - 1] / successive[index])

    with (args.output_dir / "time_refinement.csv").open(
            "w", newline="", encoding="utf-8") as stream:
        fields = ["cells_per_edge", "dofs", "dt", "transient_min",
                  "transient_max", "global_t50", "difference_to_finest",
                  "successive_difference", "rate"]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for index, (cells, dt, data) in enumerate(temporal):
            writer.writerow({
                "cells_per_edge": cells,
                "dofs": 83 + (cells - 1) * 1130,
                "dt": dt,
                "transient_min": f"{min(data['min']):.12g}",
                "transient_max": f"{max(data['max']):.12g}",
                "global_t50": f"{crossing(data['time'], data['global']):.12g}",
                "difference_to_finest": f"{differences_to_reference[index]:.12g}",
                "successive_difference": (
                    "" if successive[index] is None
                    else f"{successive[index]:.12g}"
                ),
                "rate": "" if rates[index] is None else f"{rates[index]:.12g}",
            })

    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for _, dt, data in temporal:
        axes[0].plot(data["time"], data["global"], linewidth=1.8,
                     label=rf"$\Delta t={dt:g}$")
        axes[1].plot(data["time"], data["min"], linewidth=1.8,
                     label=rf"$\Delta t={dt:g}$")
    axes[0].set_xlabel("time (years)")
    axes[0].set_ylabel("global biomarker (%)")
    axes[0].set_ylim(0.0, 102.0)
    axes[1].set_xlabel("time (years)")
    axes[1].set_ylabel("minimum nodal concentration")
    for axis in axes:
        axis.grid(True, linewidth=0.5, alpha=0.35)
        axis.legend()
    figure.tight_layout()
    figure.savefig(args.output_dir / "time_refinement.png", dpi=220)
    figure.savefig(args.output_dir / "time_refinement.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
