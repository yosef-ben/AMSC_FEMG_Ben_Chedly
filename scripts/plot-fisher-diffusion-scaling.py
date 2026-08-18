#!/usr/bin/env python3

"""Report figure for the Fisher-Kolmogorov diffusion-scaling study.

Two layouts are written from the same stored sweep. ``--layout curves`` is
the report figure: the lobe biomarkers at three transport scalings. The
default ``--layout full`` adds the spread-against-Damkohler panel with the
Fornari setting, the published separation and the scalings at which the
finite element solution leaves [0,1]; it stays in the benchmark record.
"""

import argparse
import csv
import json
from pathlib import Path

import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style

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
    parser.add_argument("--layout", choices=("full", "curves"),
                        default="full")
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

    figure_style.apply()
    if args.layout == "curves":
        figure = plt.figure(figsize=(9.8, 3.4))
        grid = figure.add_gridspec(1, 3, wspace=0.24)
    else:
        figure = plt.figure(figsize=(9.8, 7.2))
        grid = figure.add_gridspec(2, 3, height_ratios=[1.0, 1.2],
                                   hspace=0.46, wspace=0.24)

    for column, scaling in enumerate(args.panels):
        axes = figure.add_subplot(grid[0, column])
        key = f"rho_{scaling:g}"
        time = curves[f"{key}_time"]
        axes.axhline(50, color="0.6", linewidth=0.9,
                     linestyle=(0, (4, 3)), zorder=1)
        for lobe, (colour, style) in LOBE_STYLE.items():
            axes.plot(time, curves[f"{key}_{lobe}"], color=colour,
                      linestyle=style, linewidth=1.8, zorder=3)
        axes.plot(time, curves[f"{key}_global"], color="0.35",
                  linestyle=(0, (6, 3)), linewidth=1.2, zorder=2)
        index = int(np.argmin(np.abs(scalings - scaling)))
        axes.text(0.05, 0.955,
                  rf"$\rho = {scaling:g}$" "\n"
                  rf"Da $= {damkohler[index]:.1f}$",
                  transform=axes.transAxes, fontsize=9,
                  fontweight="bold", color="0.35", va="top")
        # The first 25 of the simulated years, the window of the biomarker
        # comparison figure; the curves have saturated by then.
        axes.set_xlim(0, 25)
        axes.set_ylim(0, 102)
        axes.set_xticks([0, 5, 10, 15, 20, 25])
        axes.set_yticks([0, 50, 100])
        figure_style.xname(axes, "t [yr]", y=-0.135, fontsize=9)
        if column == 0:
            axes.set_ylabel("biomarker abnormality [%]", labelpad=2)
            # Colour and line style to lobe, once for the whole figure.
            for slot, (lobe, (colour, _)) in enumerate(LOBE_STYLE.items()):
                figure_style.label_series(axes, 16.5, 42.0 - 10.5 * slot,
                                          lobe, colour, fontsize=8.5)
        else:
            axes.set_yticklabels([])
        axes.text(0.0, 1.05, f"({chr(97 + column)})",
                  transform=axes.transAxes, fontsize=10.5,
                  fontweight="bold", style="italic", va="bottom")

    if args.layout == "curves":
        args.output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.output, dpi=200, bbox_inches="tight")
        figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight")
        print(f"Written {args.output} and its PDF (curves layout)")
        return

    axes = figure.add_subplot(grid[1, :])
    unbounded_max = damkohler[~bounded].min() if (~bounded).any() else None
    if unbounded_max is not None:
        # The shading starts at the first scaling of the sweep at which the
        # finite element solution leaves [0,1]; what it means is stated in
        # the caption of the report, not inside the panel.
        axes.axvspan(unbounded_max, damkohler.max() * 1.6, color="#D62728",
                     alpha=0.07, zorder=0)

    # Grey band: the separation of figure 7 of Fornari et al., described in
    # the caption.
    axes.axhspan(FORNARI_SPREAD_YEARS * 0.85, FORNARI_SPREAD_YEARS * 1.15,
                 color="0.55", alpha=0.16, zorder=0)

    axes.plot(damkohler, spread, color="#333333", linewidth=1.6,
              marker="o", markersize=5, zorder=3)

    # One dashed line marks the setting of the comparison with Fornari et
    # al. (rho = 1); the deterministic model of Corti et al. is placed on
    # this curve in its own section, not here.
    axes.axvline(summary["alpha"] / summary["fiedler_value"], color="#1F77B4",
                 linewidth=1.4, linestyle="--", zorder=2)

    axes.set_xscale("log")
    axes.set_yscale("log")
    axes.set_xlim(0.4, damkohler.max() * 1.6)
    axes.set_ylim(1e-7, 60)
    axes.set_xticks([1, 10, 100])
    axes.set_xticklabels(["1", "10", "100"])
    axes.set_yticks([1e-6, 1e-4, 1e-2, 1])
    axes.set_yticklabels([r"$10^{-6}$", r"$10^{-4}$", r"$10^{-2}$", "1"])
    axes.minorticks_off()
    axes.set_ylabel("lobe activation spread [years]", labelpad=2)
    figure_style.xname(axes, r"Da $= \alpha/(\rho\lambda_2)$", y=-0.10)
    axes.text(0.0, 1.03, "(d)", transform=axes.transAxes, fontsize=10.5,
              fontweight="bold", style="italic", va="bottom")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, bbox_inches="tight")
    stem = args.output.with_suffix(".pdf")
    figure.savefig(stem, bbox_inches="tight")
    print(f"Written {args.output} and {stem}")


if __name__ == "__main__":
    main()
