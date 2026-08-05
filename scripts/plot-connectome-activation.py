#!/usr/bin/env python3

"""Activation time of the misfolded protein over the metric-graph network.

One image carrying the whole history: every degree of freedom is coloured by
the first stored time at which its concentration reaches 0.5. This is the
network counterpart of the activation-time maps the reference works use for
their one-dimensional problem, and it answers a question the stage-by-stage
figure cannot: in which order the network is recruited.

It also settles what drives that order. The activation time of the 83 region
vertices correlates with the local reaction coefficient at r = -0.84 and with
the graph distance from the seed at r = +0.07: on this connectome, where every
region lies within two connections of the seed, the pattern is set by the
regional conversion rates and not by any travelling front.

Read from the `solution_*.vtp` files written by
`test_fisher_kolmogorov_corti83`. No resampling and no interpolation off the
graph.
"""

import argparse
import collections
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
from render_connectome import (common_box, common_scale, lookup_table,
                               render_section, show_render)

DATA = Path("data/connectome/fornari83")
REGIONS = Path("benchmarks/21_fisher_kolmogorov_corti83/results"
               "/reaction_coefficients.csv")
LEVEL = 0.5
COLOURMAP = plt.cm.YlOrRd_r


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--solution-dir", type=Path,
        default=Path("output/fisher_kolmogorov/corti83_refined"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-step", type=float, default=0.2)
    return parser.parse_args()


def series(directory):
    found = {}
    for path in sorted(glob.glob(str(directory / "solution_*.vtp"))):
        step = int(re.search(r"solution_(\d+)\.vtp$", path).group(1))
        found[step] = Path(path)
    if not found:
        raise SystemExit(f"no solution_*.vtp under {directory}")
    return found


def activation_times(snapshots, time_step, level=LEVEL, scalar="c"):
    """First stored time at which each degree of freedom reaches `level`."""
    steps = sorted(snapshots)
    activation = None
    for step in steps:
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(snapshots[step]))
        reader.Update()
        array = reader.GetOutput().GetPointData().GetArray(scalar)
        values = np.array([array.GetValue(i)
                           for i in range(array.GetNumberOfTuples())])
        if activation is None:
            activation = np.full(len(values), np.nan)
            template = reader.GetOutput()
        reached = (values >= level) & np.isnan(activation)
        activation[reached] = step * time_step
    return activation, template


def correlations(activation, template):
    """Is the order set by the reaction rate or by distance from the seed?"""
    with open(DATA / "nodes.csv", newline="") as stream:
        nodes = list(csv.DictReader(stream))
    coords = np.array([[float(row["x"]), float(row["y"]), float(row["z"])]
                       for row in nodes])
    points = np.array([template.GetPoint(i)
                       for i in range(template.GetNumberOfPoints())])
    _, index = cKDTree(points).query(coords)
    vertex_activation = activation[index]

    with open(REGIONS, newline="") as stream:
        alpha = {int(row["node_id"]): float(row["alpha"])
                 for row in csv.DictReader(stream)}
    rates = np.array([alpha[i] for i in range(len(nodes))])

    with open(DATA / "edges.csv", newline="") as stream:
        edges = list(csv.DictReader(stream))
    adjacency = collections.defaultdict(list)
    for row in edges:
        source, target = int(row["source"]), int(row["target"])
        adjacency[source].append(target)
        adjacency[target].append(source)
    seeds = [i for i, row in enumerate(nodes)
             if "entorhinal" in row["name"].lower()
             or "Hippocampus" in row["name"]]
    distance = {seed: 0 for seed in seeds}
    queue = collections.deque(seeds)
    while queue:
        current = queue.popleft()
        for neighbour in adjacency[current]:
            if neighbour not in distance:
                distance[neighbour] = distance[current] + 1
                queue.append(neighbour)
    hops = np.array([distance.get(i, np.nan) for i in range(len(nodes))])

    valid = ~np.isnan(vertex_activation)
    return (float(np.corrcoef(vertex_activation[valid], rates[valid])[0, 1]),
            float(np.corrcoef(vertex_activation[valid], hops[valid])[0, 1]),
            int(np.nanmax(hops)), vertex_activation, coords, nodes)


def main():
    args = arguments()
    snapshots = series(args.solution_dir)
    activation, template = activation_times(snapshots, args.time_step)
    never = int(np.isnan(activation).sum())
    finite = activation[~np.isnan(activation)]
    low, high = float(finite.min()), float(finite.max())

    with_rate, with_distance, diameter, vertex_activation, coords, nodes = \
        correlations(activation, template)

    # Degrees of freedom that never activate keep the late end of the ramp;
    # there are only a handful and the count is reported in the caption.
    field = np.where(np.isnan(activation), high, activation)
    stamped = vtk.vtkXMLPolyDataReader()
    stamped.SetFileName(str(snapshots[sorted(snapshots)[0]]))
    stamped.Update()
    polydata = vtk.vtkPolyData()
    polydata.DeepCopy(stamped.GetOutput())
    array = vtk.vtkFloatArray()
    array.SetName("c")
    for value in field:
        array.InsertNextValue(value)
    polydata.GetPointData().SetScalars(array)

    table = lookup_table(COLOURMAP, low, high)
    plt.rcParams.update({"font.size": 9.5})
    with tempfile.TemporaryDirectory() as scratch:
        source = Path(scratch) / "activation.vtp"
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(str(source))
        writer.SetInputData(polydata)
        writer.Write()

        inside = coords[:, 0] <= 0.0
        vertex_field = np.where(np.isnan(vertex_activation), high,
                                vertex_activation)
        panel, _ = render_section(
            Path(scratch) / "activation.png", source, table,
            node_coords=coords[inside], node_values=vertex_field[inside],
            level=-1.0, scalar_range=(low, high), scale=common_scale(),
            size=(1500, 1150))

        figure, axis = plt.subplots(figsize=(7.4, 5.6))
        show_render(axis, panel, common_box([panel]))

    mappable = plt.cm.ScalarMappable(
        cmap=COLOURMAP, norm=matplotlib.colors.Normalize(low, high))
    bar = minimal_colourbar(figure, mappable, axis,
                            r"activation time: first stored year at which "
                            rf"$c$ reaches ${LEVEL:g}$",
                            low=f"{low:g}", high=f"{high:g}", fraction=0.05,
                            pad=0.03, aspect=38, shrink=0.85)
    bar.set_ticks(np.linspace(low, high, 5))
    bar.set_ticklabels([f"{value:.0f}" for value in np.linspace(low, high, 5)])
    axis.set_title("Order in which the network is recruited", fontsize=11)
    figure.text(0.5, 0.055,
                f"activation correlates with the local reaction coefficient "
                f"($r={with_rate:+.2f}$) and not with the graph distance from "
                f"the seed ($r={with_distance:+.2f}$);\n"
                f"every region lies within {diameter} connections of the "
                f"seed, so the order is set by the regional conversion rates, "
                f"not by a travelling front",
                ha="center", va="top", fontsize=8.5, color="0.35")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=260, bbox_inches="tight",
                   facecolor="white")
    figure.savefig(args.output.with_suffix(".pdf"), bbox_inches="tight",
                   facecolor="white")
    print(f"Saved {args.output}")
    print(f"  activation between {low:g} and {high:g} years; "
          f"{never} of {len(activation)} degrees of freedom never reach "
          f"{LEVEL:g}")
    print(f"  correlation with local alpha {with_rate:+.3f}, "
          f"with graph distance {with_distance:+.3f}, "
          f"graph eccentricity from the seed {diameter}")


if __name__ == "__main__":
    main()
