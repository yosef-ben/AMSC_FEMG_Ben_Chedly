#!/usr/bin/env python3

"""Sagittal section of the spreading on the metric-graph discretisation.

Counterpart of the whole-brain progression figures of Weickenmeier et al. Their
images are sections of a continuum field on a tetrahedral mesh of about 80000
unknowns; the substrate here is a metric graph, so what is drawn is exactly the
degrees of freedom the simulation carried.

The construction follows the reference in three respects. The brain is an
opaque object cut at the mid-sagittal plane, not a translucent volume: the pial
surface is clipped to the far half-space and the graph to the near one, so a
single depth buffer occludes correctly and the network sits inside the anatomy.
The discretisation stays visible as a pale neutral network where the
concentration is low, because the zero of the colour ramp is that same neutral
grey. And the stages run left to right on one fixed colour scale.

Everything is read from the `solution_*.vtp` files written by
`test_fisher_kolmogorov_corti83`. Nothing is resampled, smoothed, thresholded
away or interpolated into the brain volume.
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

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import minimal_colourbar
from render_connectome import (common_box, common_scale, field_colourmap,
                               lookup_table, render_section, show_render)

NODES = Path("data/connectome/fornari83/nodes.csv")
LEVEL = 0.5


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solution-dir", type=Path,
        default=Path("output/fisher_kolmogorov/corti83_refined"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-step", type=float, default=0.2)
    parser.add_argument("--stages", type=float, nargs="+",
                        default=[0.2, 0.4, 0.6],
                        help="network mean concentrations selecting the stages")
    return parser.parse_args()


def series(directory):
    found = {}
    for path in sorted(glob.glob(str(directory / "solution_*.vtp"))):
        step = int(re.search(r"solution_(\d+)\.vtp$", path).group(1))
        found[step] = Path(path)
    if not found:
        raise SystemExit(f"no solution_*.vtp under {directory}")
    return found


def read_field(path, scalar="c"):
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    output = reader.GetOutput()
    array = output.GetPointData().GetArray(scalar)
    values = np.array([array.GetValue(i)
                       for i in range(array.GetNumberOfTuples())])
    points = np.array([output.GetPoint(i)
                       for i in range(output.GetNumberOfPoints())])
    return points, values


def supra_threshold_length(path, level=LEVEL, scalar="c"):
    """Fraction of the network length carrying c above `level`.

    Exact for a P1 field: on a segment with endpoint values a and b the
    supra-threshold fraction is 0, 1, or (max-level)/(max-min).
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    output = reader.GetOutput()
    values = output.GetPointData().GetArray(scalar)
    lines = output.GetLines()
    lines.InitTraversal()
    ids = vtk.vtkIdList()
    above = total = 0.0
    while lines.GetNextCell(ids):
        for k in range(ids.GetNumberOfIds() - 1):
            first, second = ids.GetId(k), ids.GetId(k + 1)
            length = np.linalg.norm(np.array(output.GetPoint(first))
                                    - np.array(output.GetPoint(second)))
            a, b = values.GetValue(first), values.GetValue(second)
            low, high = min(a, b), max(a, b)
            if high <= level:
                fraction = 0.0
            elif low >= level:
                fraction = 1.0
            else:
                fraction = (high - level) / (high - low)
            above += fraction * length
            total += length
    return above / total


def main():
    args = arguments()
    snapshots = series(args.solution_dir)
    steps = sorted(snapshots)

    means = {}
    for step in steps:
        _, values = read_field(snapshots[step])
        means[step] = values.mean()
    chosen = [min(steps, key=lambda step: abs(means[step] - target))
              for target in args.stages]

    # Region vertices, taken as a lookup and not as an interpolation: the
    # nearest degree of freedom to each anatomical coordinate is that vertex.
    with open(NODES, newline="") as stream:
        nodes = list(csv.DictReader(stream))
    coords = np.array([[float(row["x"]), float(row["y"]), float(row["z"])]
                       for row in nodes])
    reference, _ = read_field(snapshots[steps[0]])
    tree = cKDTree(reference)
    gap, index = tree.query(coords)
    if gap.max() > 1e-3:
        raise SystemExit(f"region vertices are not degrees of freedom "
                         f"(largest gap {gap.max():.3g} mm)")
    inside = coords[:, 0] <= 0.0

    colourmap = field_colourmap()
    table = lookup_table(colourmap, 0.0, 1.0)
    scale = common_scale()

    plt.rcParams.update({"font.size": 9.5})
    with tempfile.TemporaryDirectory() as scratch:
        panels, marks, fractions = [], [], []
        for step in chosen:
            _, values = read_field(snapshots[step])
            path, count = render_section(
                Path(scratch) / f"stage_{step}.png", snapshots[step], table,
                node_coords=coords[inside], node_values=values[index][inside],
                level=LEVEL, scale=scale)
            panels.append(path)
            marks.append(count)
            fractions.append(supra_threshold_length(snapshots[step]))

        box = common_box(panels)
        figure, axes = plt.subplots(1, len(chosen), figsize=(4.1 * len(chosen),
                                                             3.6))
        axes = np.atleast_1d(axes)
        for column, (axis, step, fraction) in enumerate(
                zip(axes, chosen, fractions)):
            show_render(axis, panels[column], box)
            axis.set_title(f"$t = {step * args.time_step:g}$ years",
                           fontsize=11, pad=16)
            axis.text(0.5, 1.005,
                      rf"$\bar c = {means[step]:.2f}$,   "
                      rf"{100 * fraction:.0f}% of the network length above "
                      rf"$c = {LEVEL:g}$",
                      transform=axis.transAxes, fontsize=8.5, color="0.35",
                      ha="center", va="bottom")

    mappable = plt.cm.ScalarMappable(
        cmap=colourmap, norm=matplotlib.colors.Normalize(0.0, 1.0))
    bar = minimal_colourbar(figure, mappable, axes.tolist(),
                            "misfolded protein concentration $c$",
                            low="0", high="1", fraction=0.05, pad=0.03,
                            aspect=45, shrink=0.55)
    bar.set_ticks([0.0, 0.25, 0.5, 0.75, 1.0])
    bar.set_ticklabels(["0", "0.25", "0.5", "0.75", "1"])
    figure.subplots_adjust(wspace=0.02)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=260, bbox_inches="tight",
                   facecolor="white")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight",
                   facecolor="white")
    print(f"Saved {args.output}")
    for step, fraction, count in zip(chosen, fractions, marks):
        print(f"  t = {step * args.time_step:5g} y   mean c = "
              f"{means[step]:.4f}   length above {LEVEL:g}: "
              f"{100 * fraction:5.1f}%   level-set marks in section: {count}")


if __name__ == "__main__":
    main()
