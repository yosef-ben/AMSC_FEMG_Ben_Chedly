#!/usr/bin/env python3

"""Activation order with and without the regional biology.

In the manner of figure 9 of Fornari et al.: every vertex is coloured by its
activation time, the first stored year at which its concentration reaches
0.5. Panel (a) shows the 83 anatomical-vertex histories of the calibrated
run, each coloured by that time; panels (b) and (c) show the anatomy for the
calibrated regional rates and for the control variant in which the seven
rates are replaced by their vertex mean, so that the reaction budget and the
Damkohler number are unchanged and the only anatomy left in the model is
the connectivity.

The demonstration the figure carries: with the calibrated rates the
occipital group activates last, four clear years after every other group,
which is the clinical tail of the sequence; with the uniform rate the
spread collapses to about three years and the order becomes the
connectivity's, the occipital group moving up and the frontal group,
dragged by its weakly connected pole, falling last. The per-vertex and
per-group activation times of both variants are written next to the figure,
so every number quoted in the report is anchored to a stored file.
"""

import argparse
import csv
import glob
import re
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from scipy.spatial import cKDTree
from vtk.util import numpy_support

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style
from connectome_style import load_nodes, minimal_colourbar
import render_connectome as rc

LEVEL = 0.5
COLOURMAP = plt.cm.YlOrRd_r
GROUPS = ("frontal", "temporal", "parietal", "insular", "limbic",
          "occipital", "subcortical")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regional", type=Path,
                        default=Path("output/fisher_kolmogorov/corti83_refined"))
    parser.add_argument("--uniform", type=Path,
                        default=Path("output/fisher_kolmogorov/corti83_uniform"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--time-step", type=float, default=0.2)
    return parser.parse_args()


def snapshots(directory):
    found = {}
    for path in sorted(glob.glob(str(directory / "solution_*.vtp"))):
        step = int(re.search(r"solution_(\d+)\.vtp$", path).group(1))
        found[step] = Path(path)
    if not found:
        raise SystemExit(f"no solution_*.vtp under {directory}")
    return found


def vertex_histories(directory, coords, time_step):
    """Times, the 83 anatomical-vertex histories, and per-vertex activation."""
    series = snapshots(directory)
    index = None
    times, histories = [], []
    activation = None
    for step in sorted(series):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(series[step]))
        reader.Update()
        polydata = reader.GetOutput()
        values = numpy_support.vtk_to_numpy(
            polydata.GetPointData().GetArray("c"))
        if index is None:
            points = numpy_support.vtk_to_numpy(
                polydata.GetPoints().GetData())
            _, index = cKDTree(points).query(coords)
            activation = np.full(len(values), np.nan)
        reached = (values >= LEVEL) & np.isnan(activation)
        activation[reached] = step * time_step
        times.append(step * time_step)
        histories.append(values[index])
    return (np.array(times), np.array(histories),
            np.asarray(activation)[index])


def main():
    args = arguments()
    nodes = load_nodes()
    coords = np.array([node["coords"] for node in nodes])

    variants = {}
    for label, directory in (("regional", args.regional),
                             ("uniform", args.uniform)):
        variants[label] = vertex_histories(directory, coords, args.time_step)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with open(args.output_dir / "activation_order.csv", "w",
              newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node_id", "name", "region",
                         "regional_years", "uniform_years"])
        for node in nodes:
            k = node["node_id"]
            writer.writerow([k, node["name"], node["region"],
                             f"{variants['regional'][2][k]:.12g}",
                             f"{variants['uniform'][2][k]:.12g}"])
    with open(args.output_dir / "activation_order_groups.csv", "w",
              newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["variant", "region", "mean_activation_years"])
        for label in ("regional", "uniform"):
            act = variants[label][2]
            for group in GROUPS:
                members = [n["node_id"] for n in nodes
                           if n["region"] == group]
                writer.writerow([label, group,
                                 f"{np.mean(act[members]):.12g}"])

    both = np.concatenate([variants["regional"][2], variants["uniform"][2]])
    low, high = float(both.min()), float(both.max())
    table = rc.lookup_table(COLOURMAP, low, high)
    norm = matplotlib.colors.Normalize(low, high)
    scale = rc.common_scale()

    figure_style.apply()
    figure = plt.figure(figsize=(11.4, 3.6))
    grid = figure.add_gridspec(1, 3, width_ratios=[1.12, 1.0, 1.0],
                               left=0.055, right=0.985, top=0.90,
                               bottom=0.16, wspace=0.16)
    curves = figure.add_subplot(grid[0, 0])
    anatomy = [figure.add_subplot(grid[0, 1]),
               figure.add_subplot(grid[0, 2])]

    times, histories, activation = variants["regional"]
    # Late vertices are drawn last so the clinically meaningful tail of the
    # sequence stays visible on top of the bundle.
    for k in np.argsort(activation):
        curves.plot(times, histories[:, k],
                    color=COLOURMAP(norm(activation[k])),
                    linewidth=1.0, alpha=0.9, zorder=2)
    curves.axhline(LEVEL, color="0.6", linewidth=0.9,
                   linestyle=(0, (4, 3)), zorder=1)
    curves.set_xlim(0.0, 60.0)
    curves.set_ylim(0.0, 1.02)
    curves.set_xticks([0, 15, 30, 45, 60])
    curves.set_yticks([0.0, 0.5, 1.0])
    curves.set_ylabel("concentration", labelpad=2)
    figure_style.xname(curves, "time [years]", y=-0.115)
    curves.text(0.03, 0.955, "regional rates", transform=curves.transAxes,
                fontsize=9.5, fontweight="bold", color="0.35", va="top")

    with tempfile.TemporaryDirectory() as scratch:
        panels = []
        for axis, (label, title) in zip(
                anatomy, (("regional", "regional rates"),
                          ("uniform", "uniform mean rate"))):
            panels.append((axis, title, rc.render(
                Path(scratch) / f"{label}.png", "sagittal", coords,
                variants[label][2], table, node_radius=3.4,
                size=(2000, 1650), scale=scale,
                surface_opacity=0.075)[0]))
        box = rc.common_box([panel for _, _, panel in panels])
        for axis, title, panel in panels:
            rc.show_render(axis, panel, box)
            axis.text(0.03, 0.985, title, transform=axis.transAxes,
                      fontsize=9.5, fontweight="bold", color="0.35",
                      va="top")

    for letter, axis in zip("abc", (curves, *anatomy)):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")

    mappable = plt.cm.ScalarMappable(cmap=COLOURMAP, norm=norm)
    bar = minimal_colourbar(figure, mappable, anatomy,
                            rf"activation time [years]: first stored "
                            rf"crossing of $c = {LEVEL:g}$",
                            low=f"{low:g}", high=f"{high:g}",
                            fraction=0.045, pad=0.055, aspect=44,
                            shrink=0.62)
    bar.set_ticks(np.linspace(low, high, 5))
    bar.set_ticklabels([f"{value:.0f}"
                        for value in np.linspace(low, high, 5)])

    figure.savefig(args.output_dir / "activation_order.pdf",
                   facecolor="white")
    figure.savefig(args.output_dir / "activation_order.png", dpi=300,
                   facecolor="white")
    plt.close(figure)

    print(f"Saved activation_order to {args.output_dir}")
    for label in ("regional", "uniform"):
        act = variants[label][2]
        means = {group: np.mean(act[[n["node_id"] for n in nodes
                                     if n["region"] == group]])
                 for group in GROUPS}
        order = sorted(means, key=means.get)
        print(f"  {label:9s} " + "  ".join(
            f"{group}={means[group]:.2f}" for group in order))


if __name__ == "__main__":
    main()
