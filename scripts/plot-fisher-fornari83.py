#!/usr/bin/env python3

"""Plot the Fornari biomarker curves and compute activation times."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


REGIONAL_CURVES = {
    "temporal": ("Temporal", "#2ca02c", "o"),
    "frontal": ("Frontal", "#d62728", "s"),
    "parietal": ("Parietal", "#ff8c00", "^"),
    "occipital": ("Occipital", "#1f77b4", "D"),
}
CURVE_KEYS = (*REGIONAL_CURVES, "global")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodal", type=Path, required=True)
    parser.add_argument("--fem", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=50.0)
    return parser.parse_args()


def read_curves(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        key: [float(row[key]) for row in rows]
        for key in ("time", *CURVE_KEYS)
    }


def crossing_time(time, values, threshold):
    for index in range(1, len(time)):
        if values[index] >= threshold and values[index - 1] < threshold:
            fraction = ((threshold - values[index - 1])
                        / (values[index] - values[index - 1]))
            return time[index - 1] + fraction * (time[index] - time[index - 1])
    return None


def write_activation_times(path, datasets, threshold):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["method", "region", "threshold_percent", "time_years"])
        for method, data in datasets.items():
            for region in ("temporal", "frontal", "parietal", "occipital", "global"):
                time = crossing_time(data["time"], data[region], threshold)
                writer.writerow([
                    method, region, threshold,
                    "" if time is None else f"{time:.12g}",
                ])


def main():
    args = arguments()
    datasets = {"Nodal reference": read_curves(args.nodal),
                "Metric-graph FEM": read_curves(args.fem)}
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update({"font.size": 11})
    figure, axes = plt.subplots(2, 2, figsize=(11.0, 7.2), sharex="col")
    for column, (title, data) in enumerate(datasets.items()):
        absolute = axes[0, column]
        deviation = axes[1, column]

        # Draw the network average first so it cannot hide regional curves.
        absolute.plot(data["time"], data["global"], color="#555555",
                      linewidth=2.4, linestyle="--", label="Network", zorder=1)
        for key, (label, color, marker) in REGIONAL_CURVES.items():
            absolute.plot(data["time"], data[key], color=color,
                          linewidth=1.8, marker=marker, markersize=3.2,
                          markevery=5, label=label, zorder=2)
            difference = [value - average for value, average
                          in zip(data[key], data["global"])]
            deviation.plot(data["time"], difference, color=color,
                           linewidth=1.8, marker=marker, markersize=3.2,
                           markevery=5, label=label)

        absolute.axhline(args.threshold, color="#999999", linewidth=0.8,
                         linestyle=":")
        deviation.axhline(0.0, color="#777777", linewidth=0.8,
                          linestyle=":")
        # The four regional curves coincide to within the line width. Saying so
        # on the figure stops the reader from hunting for a plotting fault; the
        # reason is quantified in benchmark 23.
        crossings = [crossing_time(data["time"], data[key], args.threshold)
                     for key in REGIONAL_CURVES]
        spread = max(crossings) - min(crossings)
        absolute.annotate(
            f"all four lobes cross {args.threshold:g}%\n"
            f"within {spread:.1e} years",
            xy=(0.04, 0.72), xycoords="axes fraction", fontsize=8.5,
            color="0.3")
        absolute.set_title(title)
        absolute.set_xlim(0.0, 25.0)
        absolute.set_ylim(0.0, 102.0)
        deviation.set_xlim(0.0, 25.0)
        deviation.set_xlabel("time (years)")
        for axis in (absolute, deviation):
            axis.grid(True, which="both", linewidth=0.5, alpha=0.35)

    axes[0, 0].set_ylabel("biomarker abnormality (%)")
    axes[1, 0].set_ylabel("deviation from network\n(percentage points)")
    axes[0, 0].legend(loc="lower right", frameon=True, ncol=2)
    axes[1, 0].legend(loc="best", frameon=True, ncol=2)
    figure.suptitle(
        "Fornari et al. report a temporal-to-occipital separation of about "
        "5.5 years at this $\\alpha$;\nwith the unscaled Laplacian the "
        "separation is smaller by seven decades (see benchmark 23)",
        fontsize=9.5, color="0.3", y=0.985)
    figure.tight_layout(rect=(0, 0, 1, 0.925))
    figure.savefig(args.output_dir / "biomarker_comparison.png", dpi=220)
    figure.savefig(args.output_dir / "biomarker_comparison.pdf")
    plt.close(figure)

    write_activation_times(
        args.output_dir / "activation_times.csv", datasets, args.threshold)


if __name__ == "__main__":
    main()
