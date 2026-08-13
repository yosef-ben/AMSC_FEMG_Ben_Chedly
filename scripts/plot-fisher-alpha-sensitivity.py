#!/usr/bin/env python3

"""Plot conversion-rate sensitivity on the 83-region connectome."""

import argparse
import csv
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style


ALPHAS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)
COLORS = ("#555555", "#1f77b4", "#17a589", "#e3a008", "#e66b2e", "#b6242a")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("output/fisher_kolmogorov/fornari83_alpha"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "benchmarks/20_fisher_kolmogorov_alpha_sensitivity/results"
        ),
    )
    parser.add_argument("--threshold", type=float, default=50.0)
    return parser.parse_args()


def tag(alpha):
    return f"{alpha:.1f}".replace(".", "p")


def read_table(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        key: [float(row[key]) for row in rows]
        for key in ("time", "global", "min", "max")
    }


def read_metric_mass(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "time": [float(row["time"]) for row in rows],
        "global": [float(row["metric_global"]) for row in rows],
    }


def crossing_time(time, values, threshold):
    for index in range(1, len(time)):
        if values[index] >= threshold and values[index - 1] < threshold:
            fraction = ((threshold - values[index - 1])
                        / (values[index] - values[index - 1]))
            return time[index - 1] + fraction * (
                time[index] - time[index - 1]
            )
    return None


def verify(datasets):
    for method, cases in datasets.items():
        reference_time = cases[ALPHAS[0]]["time"]
        for alpha in ALPHAS:
            data = cases[alpha]
            if data["time"] != reference_time:
                raise ValueError(f"Inconsistent time grid for {method}")
            if min(data["min"]) < -2.0e-4 or max(data["max"]) > 1.0 + 1.0e-8:
                raise ValueError(f"Concentration bounds failed for {method}")
        for lower, upper in zip(ALPHAS, ALPHAS[1:]):
            if any(
                high + 1.0e-9 < low
                for low, high in zip(
                    cases[lower]["global"], cases[upper]["global"]
                )
            ):
                raise ValueError(f"Alpha monotonicity failed for {method}")
        drift = max(cases[0.0]["global"]) - min(cases[0.0]["global"])
        if drift > 1.0e-10:
            raise ValueError(f"The alpha=0 mass is not conserved for {method}")


def write_summary(path, datasets, threshold):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "method", "alpha", "time_50_percent", "final_global_percent",
            "global_drift_alpha_zero", "transient_min", "transient_max",
        ])
        for method, cases in datasets.items():
            zero_drift = max(cases[0.0]["global"]) - min(
                cases[0.0]["global"]
            )
            for alpha in ALPHAS:
                data = cases[alpha]
                crossing = crossing_time(
                    data["time"], data["global"], threshold
                )
                writer.writerow([
                    method,
                    alpha,
                    "" if crossing is None else f"{crossing:.12g}",
                    f"{data['global'][-1]:.12g}",
                    f"{zero_drift:.12g}" if alpha == 0.0 else "",
                    f"{min(data['min']):.12g}",
                    f"{max(data['max']):.12g}",
                ])


def main():
    args = arguments()
    datasets = {
        "Nodal reference": {
            alpha: read_table(
                args.root / f"alpha_{tag(alpha)}" / "nodal_biomarkers.csv"
            )
            for alpha in ALPHAS
        },
        "Metric-graph FEM": {
            alpha: read_table(
                args.root / f"alpha_{tag(alpha)}" / "fem_biomarkers.csv"
            )
            for alpha in ALPHAS
        },
    }
    for alpha in ALPHAS:
        metric = read_metric_mass(
            args.root / f"alpha_{tag(alpha)}" / "fem_metric_mass.csv"
        )
        fem = datasets["Metric-graph FEM"][alpha]
        if metric["time"] != fem["time"]:
            raise ValueError("Inconsistent FEM mass time grid")
        fem["global"] = metric["global"]

    verify(datasets)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_summary(
        args.output_dir / "alpha_sensitivity.csv",
        datasets,
        args.threshold,
    )

    figure_style.apply()
    figure, axes = plt.subplots(
        1, 2, figsize=(9.6, 3.9), sharex=True, sharey=True
    )
    for axis, (method, cases) in zip(axes, datasets.items()):
        axis.axhline(args.threshold, color="0.6", linewidth=0.9,
                     linestyle=(0, (4, 3)), zorder=1)
        for alpha, color in zip(ALPHAS, COLORS):
            data = cases[alpha]
            axis.plot(data["time"], data["global"], color=color,
                      linewidth=1.9, solid_capstyle="round", zorder=3)
        # Colour-to-rate column in the empty lower-right region, as in the
        # one-dimensional figures: the numbers live in the caption, the
        # column only ties each colour to its rate.
        for slot, (alpha, color) in enumerate(zip(reversed(ALPHAS),
                                                  reversed(COLORS))):
            figure_style.label_series(
                axis, 64.0, 5.5 + 6.5 * slot,
                fr"$\alpha = {alpha:.1f}$".replace("0.0", "0"),
                color, fontsize=9)
        axis.text(0.03, 0.965, method.lower(), transform=axis.transAxes,
                  fontsize=9.5, fontweight="bold", color="0.35", va="top")
        axis.set_xlim(0.0, 80.0)
        axis.set_ylim(0.0, 102.0)
        axis.set_xticks([0, 20, 40, 60, 80])
        axis.set_yticks([0, args.threshold, 100])
        axis.set_yticklabels(["0", "50", "100"])
        figure_style.xname(axis, "time [years]", y=-0.12)
    axes[0].set_ylabel("biomarker abnormality [%]", labelpad=2)
    for letter, axis in zip("ab", axes):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    figure.tight_layout(w_pad=2.6)
    figure.savefig(args.output_dir / "alpha_sensitivity.png", dpi=220)
    figure.savefig(args.output_dir / "alpha_sensitivity.pdf")
    plt.close(figure)


if __name__ == "__main__":
    main()
