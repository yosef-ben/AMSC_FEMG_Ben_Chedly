#!/usr/bin/env python3

"""Plot a stored graph solution in the style of the eigenmode figures.

Reads one solution_*.vtp written by the library, draws the metric graph in
blue at height zero with its mesh nodes and lifts the numerical solution in
red, with the same colours, view and typography as plot-eigenmodes.py, so
that the linear benchmarks of the report share one visual language. Nothing
is recomputed: the curve is the stored point data of the chosen step.
"""

import argparse
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def read_vtp(path):
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    values = [float(value) for value
              in piece.find("./Points/DataArray").text.split()]
    points = [tuple(values[i:i + 3]) for i in range(0, len(values), 3)]
    point_data = {
        array.attrib["Name"]: [float(value) for value
                               in array.text.split()]
        for array in piece.findall("./PointData/DataArray")
    }
    line_arrays = {
        array.attrib["Name"]: [int(value) for value in array.text.split()]
        for array in piece.findall("./Lines/DataArray")
    }
    lines = []
    start = 0
    for offset in line_arrays["offsets"]:
        lines.append(line_arrays["connectivity"][start:offset])
        start = offset
    return points, lines, point_data


def unique_vertices(points, ndigits=10):
    vertices = {}
    for x, y, _ in points:
        vertices[(round(x, ndigits), round(y, ndigits))] = (x, y)
    return list(vertices.values())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solution", type=Path,
                        help="one solution_*.vtp of a stored run")
    parser.add_argument("--scalar", default="u")
    parser.add_argument("--title", default=None)
    parser.add_argument("--zlim", type=float, nargs=2, default=None)
    parser.add_argument("--ztick", type=float, default=1.0)
    parser.add_argument("--elev", type=float, default=22.0)
    parser.add_argument("--azim", type=float, default=-55.0,
                        help="view azimuth; the eigenmode default hides "
                             "slopes aligned with the screen iso-heights, "
                             "rotate when the solution needs it")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    points, lines, data = read_vtp(args.solution)
    solution = data[args.scalar]

    fig = plt.figure(figsize=(4.2, 3.25))
    ax = fig.add_subplot(1, 1, 1, projection="3d", computed_zorder=False)
    for line in lines:
        xs = [points[i][0] for i in line]
        ys = [points[i][1] for i in line]
        ax.plot(xs, ys, [0.0 for _ in line], color="royalblue",
                linewidth=0.16, alpha=0.35, zorder=1)
        ax.plot(xs, ys, [solution[i] for i in line], color="red",
                linewidth=1.55, alpha=1.0, zorder=5)
    for x, y in unique_vertices(points):
        ax.scatter([x], [y], [0.0], s=8, facecolors="white",
                   edgecolors="royalblue", linewidths=0.45, alpha=0.55)

    if args.title:
        ax.set_title(args.title, fontsize=8, pad=2)
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_mid = 0.5 * (min(xs) + max(xs))
    y_mid = 0.5 * (min(ys) + max(ys))
    radius = 0.55 * max(max(xs) - min(xs), max(ys) - min(ys))
    ax.set_xlim(x_mid - radius, x_mid + radius)
    ax.set_ylim(y_mid - radius, y_mid + radius)
    if args.zlim:
        ax.set_zlim(*args.zlim)
    ax.set_xlabel("x", fontsize=7, labelpad=-6)
    ax.set_ylabel("y", fontsize=7, labelpad=-6)
    ax.set_zlabel(args.scalar, fontsize=7, labelpad=-10)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.zaxis.set_major_locator(MultipleLocator(args.ztick))
    ax.tick_params(axis="both", labelsize=6, pad=-2)
    ax.tick_params(axis="z", labelsize=6, pad=-2)
    ax.grid(True, linewidth=0.22, alpha=0.22)
    ax.view_init(elev=args.elev, azim=args.azim)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
