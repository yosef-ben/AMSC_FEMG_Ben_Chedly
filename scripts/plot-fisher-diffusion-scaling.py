#!/usr/bin/env python3

"""Report figure for the Fisher-Kolmogorov diffusion-scaling study."""

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Hue assignment of Fornari et al. figure 7. Checked for colour-vision
# deficiency: every pair separates by at least 10.8 in OKLab (x100) under
# protan, deutan and tritan simulation. Line styles repeat the identity so it
# is never carried by colour alone.
LOBE_STYLE = {
    "temporal": ("#2CA02C", "-"),
    "frontal": ("#D62728", "--"),
    "parietal": ("#FF7F0E", "-."),
    "occipital": ("#1F77B4", ":"),
}
# Lobe separation read from Fornari et al. figure 7 (temporal ~10 yr to
# occipital ~15.5 yr at alpha = 0.5).
FORNARI_SPREAD_YEARS = 5.5


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--panels", type=float, nargs=3,
                        default=[1.0, 0.05, 0.005])
    parser.add_argument("--benchmark-21-coefficients", type=Path,
                        default=Path("benchmarks/21_fisher_kolmogorov_corti83"
                                     "/results/reaction_coefficients.csv"))
    return parser.parse_args()


def mean_reaction_coefficient(path):
    """Mean of the seven Corti et al. regional rates over the 83 vertices."""
    with open(path) as stream:
        values = [float(row["alpha"]) for row in csv.DictReader(stream)]
    return sum(values) / len(values)


def main():
    args = arguments()
    with open(args.results_dir / "diffusion_scaling.csv") as stream:
        rows = list(csv.DictReader(stream))
    with open(args.results_dir / "diffusion_scaling_summary.json") as stream:
        summary = json.load(stream)
    curves = np.load(args.results_dir / "curves.npz")

    scalings = np.array([float(r["diffusion_scaling"]) for r in rows])
    damkohler = np.array([float(r["damkohler"]) for r in rows])
    spread = np.array([float(r["lobe_spread_years"]) for r in rows])
    bounded = np.array([r["fem_corti_crank_nicolson"] == "bounded"
                        for r in rows])

    figure = plt.figure(figsize=(11.0, 7.6))
    grid = figure.add_gridspec(2, 3, height_ratios=[1.0, 1.15], hspace=0.42,
                               wspace=0.28)

    for column, scaling in enumerate(args.panels):
        axes = figure.add_subplot(grid[0, column])
        key = f"rho_{scaling:g}"
        time = curves[f"{key}_time"]
        for lobe, (colour, style) in LOBE_STYLE.items():
            axes.plot(time, curves[f"{key}_{lobe}"], color=colour,
                      linestyle=style, linewidth=1.8, label=lobe)
        axes.plot(time, curves[f"{key}_global"], color="0.35",
                  linestyle=(0, (6, 3)), linewidth=1.2, label="network")
        index = int(np.argmin(np.abs(scalings - scaling)))
        axes.set_title(rf"$\rho = {scaling:g}$,  Da $= {damkohler[index]:.1f}$",
                       fontsize=10)
        axes.set_xlim(0, 40)
        axes.set_ylim(0, 100)
        axes.set_xlabel("time (years)")
        if column == 0:
            axes.set_ylabel("biomarker abnormality (%)")
        axes.grid(alpha=0.25, linewidth=0.6)
        axes.axhline(50, color="0.75", linewidth=0.7, zorder=0)
        if column == 0:
            axes.legend(fontsize=8, loc="lower right", framealpha=0.9)

    axes = figure.add_subplot(grid[1, :])
    unbounded_max = damkohler[~bounded].min() if (~bounded).any() else None
    if unbounded_max is not None:
        axes.axvspan(unbounded_max, damkohler.max() * 1.6, color="#D62728",
                     alpha=0.07, zorder=0)
        axes.text(unbounded_max * 1.2, 2.0e-7,
                  "metric-graph FEM leaves [0,1]\n"
                  "(consistent mass, both schemes)",
                  fontsize=8.5, color="#8B1A1A", va="bottom")

    axes.axhspan(FORNARI_SPREAD_YEARS * 0.85, FORNARI_SPREAD_YEARS * 1.15,
                 color="0.55", alpha=0.16, zorder=0)
    axes.text(0.46, FORNARI_SPREAD_YEARS * 0.55,
              "separation reported by Fornari et al., figure 7",
              fontsize=8.5, color="0.25", va="top", ha="left")

    axes.plot(damkohler, spread, color="#333333", linewidth=1.6,
              marker="o", markersize=5, zorder=3)

    alpha_21 = mean_reaction_coefficient(args.benchmark_21_coefficients)
    for label, scaling, alpha, colour in (
            ("benchmark 19\n" + rf"$\rho=1$, $\alpha={summary['alpha']:g}$",
             1.0, summary["alpha"], "#1F77B4"),
            ("benchmark 21\n" + rf"$\rho=1/\max(w)$, $\bar\alpha={alpha_21:.3f}$",
             summary["benchmark_21_scaling"], alpha_21, "#FF7F0E")):
        value = alpha / (scaling * summary["fiedler_value"])
        axes.axvline(value, color=colour, linewidth=1.4, linestyle="--",
                     zorder=2)
        axes.text(value * 1.08, 30.0, label, fontsize=8.5, color=colour,
                  va="top")

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(0.4, damkohler.max() * 1.6)
    axes.set_ylim(1e-7, 60)
    axes.set_xlabel(r"Damkohler number  Da $= \alpha / (\rho\,\lambda_2)$"
                    "     (reaction rate / graph homogenisation rate)")
    axes.set_ylabel("lobe activation spread (years)")
    axes.grid(alpha=0.25, linewidth=0.6, which="both")
    axes.set_title("Lobe separation is set by the diffusion scaling, "
                   "not by the connectome topology", fontsize=10)

    figure.suptitle("Benchmark 23: diffusion scaling of the 83-region "
                    "Fisher-Kolmogorov connectome", fontsize=12)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, bbox_inches="tight")
    stem = args.output.with_suffix(".pdf")
    figure.savefig(stem, bbox_inches="tight")
    print(f"Written {args.output} and {stem}")


if __name__ == "__main__":
    main()
