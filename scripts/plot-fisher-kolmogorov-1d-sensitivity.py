#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def activation_map(times, concentrations, levels):
    result = np.full((len(levels), concentrations.shape[1]), times[-1])
    for level_index, level in enumerate(levels):
        reached = concentrations >= level
        for x_index in range(concentrations.shape[1]):
            indices = np.flatnonzero(reached[:, x_index])
            if indices.size == 0:
                continue
            first = int(indices[0])
            if first == 0:
                result[level_index, x_index] = times[0]
                continue
            c0 = concentrations[first - 1, x_index]
            c1 = concentrations[first, x_index]
            t0 = times[first - 1]
            t1 = times[first]
            fraction = 0.0 if c1 == c0 else (level - c0) / (c1 - c0)
            result[level_index, x_index] = t0 + fraction * (t1 - t0)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Plot the one-dimensional Fisher-Kolmogorov sensitivity test."
    )
    parser.add_argument("csv", type=Path)
    parser.add_argument(
        "-o", "--output", type=Path,
        default=Path("fisher_kolmogorov_1d_sensitivity.pdf"),
    )
    parser.add_argument(
        "--diagnostics", type=Path,
        help="Optional CSV file for measured and theoretical front speeds.",
    )
    args = parser.parse_args()

    profiles = defaultdict(list)
    with args.csv.open(newline="") as stream:
        for row in csv.DictReader(stream):
            key = (float(row["diffusion"]), float(row["alpha"]), float(row["time"]))
            profiles[key].append((float(row["x"]), float(row["c"])))

    diffusions = sorted({key[0] for key in profiles})
    alphas = sorted({key[1] for key in profiles})
    times = np.array(sorted({key[2] for key in profiles}))
    diffusion_labels = {
        1.0e-4: r"$d=10^{-4}$",
        2.0e-4: r"$d=2\cdot10^{-4}$",
        4.0e-4: r"$d=4\cdot10^{-4}$",
    }
    concentration_levels = np.linspace(0.0, 1.0, 201)
    time_levels = np.linspace(0.0, 20.0, 11)

    plt.rcParams.update({"font.size": 10})
    figure, axes = plt.subplots(
        len(diffusions), len(alphas), figsize=(9.2, 7.0),
        sharex=True, sharey=True,
    )
    figure.subplots_adjust(
        left=0.17, right=0.88, bottom=0.16, top=0.91,
        wspace=0.20, hspace=0.12,
    )

    image = None
    for row, diffusion in enumerate(diffusions):
        for column, alpha in enumerate(alphas):
            axis = axes[row, column]
            first_profile = sorted(profiles[(diffusion, alpha, float(times[0]))])
            x_values = np.array([point[0] for point in first_profile])
            concentrations = np.array([
                [value for _, value in sorted(profiles[(diffusion, alpha, float(time))])]
                for time in times
            ])
            activation = activation_map(times, concentrations, concentration_levels)
            image = axis.contourf(
                x_values, concentration_levels, activation,
                levels=time_levels, cmap="jet_r",
            )
            axis.contour(
                x_values, concentration_levels, activation,
                levels=time_levels, colors="black", linewidths=0.45,
            )
            if row == 0:
                axis.set_title(
                    rf"$\alpha={alpha:g}$", fontsize=11, pad=9
                )
            if column == 0:
                axis.text(
                    -0.40, 0.5, diffusion_labels[diffusion],
                    transform=axis.transAxes, rotation=90,
                    va="center", ha="center", fontsize=10,
                )
                axis.set_ylabel(r"$c$")
            if row == len(diffusions) - 1:
                axis.set_xlabel(r"$x$")
            axis.set(xlim=(-1.0, 1.0), ylim=(0.0, 1.0))

            axis.set_xticks([-1.0, -0.5, 0.0, 0.5, 1.0])
            if row == len(diffusions) - 1:
                axis.set_xticklabels(
                    [r"$-1$", r"$-0.5$", r"$0$", r"$0.5$", r"$1$"]
                )

    colorbar_axis = figure.add_axes([0.31, 0.085, 0.38, 0.022])
    colorbar = figure.colorbar(
        image, cax=colorbar_axis, orientation="horizontal",
        ticks=[0.0, 5.0, 10.0, 15.0, 20.0],
    )
    colorbar.set_label(r"activation time $t$")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=300)
    figure.savefig(args.output.with_suffix(".png"), dpi=300)
    print(f"Saved {args.output}")

    if args.diagnostics is not None:
        rows = []
        for diffusion in diffusions:
            alpha = min(alphas)
            front_positions = []
            for time in times:
                if time < 10.0:
                    continue
                profile = sorted(
                    point for point in profiles[(diffusion, alpha, float(time))]
                    if point[0] >= 0.0
                )
                crossing = None
                for left, right in zip(profile, profile[1:]):
                    x0, c0 = left
                    x1, c1 = right
                    if c0 >= 0.5 and c1 < 0.5:
                        crossing = x0 + (0.5 - c0) * (x1 - x0) / (c1 - c0)
                        break
                if crossing is not None and crossing < 0.9:
                    front_positions.append((float(time), crossing))
            if len(front_positions) < 2:
                continue
            mean_time = sum(point[0] for point in front_positions) / len(front_positions)
            mean_position = sum(point[1] for point in front_positions) / len(front_positions)
            numerator = sum(
                (time - mean_time) * (position - mean_position)
                for time, position in front_positions
            )
            denominator = sum(
                (time - mean_time) ** 2 for time, _ in front_positions
            )
            measured_speed = numerator / denominator
            theoretical_speed = 2.0 * math.sqrt(diffusion * alpha)
            rows.append({
                "diffusion": diffusion,
                "alpha": alpha,
                "measured_speed": measured_speed,
                "theoretical_speed": theoretical_speed,
                "relative_error": abs(measured_speed - theoretical_speed) / theoretical_speed,
                "samples": len(front_positions),
            })

        args.diagnostics.parent.mkdir(parents=True, exist_ok=True)
        with args.diagnostics.open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    "diffusion", "alpha", "measured_speed",
                    "theoretical_speed", "relative_error", "samples",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Saved {args.diagnostics}")


if __name__ == "__main__":
    main()
