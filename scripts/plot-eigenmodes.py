#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


def read_vtp(path):
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")

    points_text = piece.find("./Points/DataArray").text
    point_values = [float(value) for value in points_text.split()]
    points = [
        tuple(point_values[i:i + 3])
        for i in range(0, len(point_values), 3)
    ]

    point_data = {}
    for array in piece.findall("./PointData/DataArray"):
        name = array.attrib["Name"]
        point_data[name] = [float(value) for value in array.text.split()]

    line_arrays = {
        array.attrib["Name"]: [int(value) for value in array.text.split()]
        for array in piece.findall("./Lines/DataArray")
    }
    connectivity = line_arrays["connectivity"]
    offsets = line_arrays["offsets"]

    lines = []
    start = 0
    for offset in offsets:
        lines.append(connectivity[start:offset])
        start = offset

    return points, lines, point_data


def read_eigenvalues(path):
    eigenvalues = {}
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            eigenvalues[int(row["mode"])] = float(row["eigenvalue"])
    return eigenvalues


def unique_vertices(points, ndigits=10):
    vertices = {}
    for x, y, _ in points:
        vertices[(round(x, ndigits), round(y, ndigits))] = (x, y)
    return list(vertices.values())


def set_equal_xy(ax, points):
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    x_mid = 0.5 * (min(xs) + max(xs))
    y_mid = 0.5 * (min(ys) + max(ys))
    radius = 0.55 * max(max(xs) - min(xs), max(ys) - min(ys))
    ax.set_xlim(x_mid - radius, x_mid + radius)
    ax.set_ylim(y_mid - radius, y_mid + radius)


def format_eigenvalue(value):
    if abs(value) < 1.0e-10:
        return "0.0000"
    return f"{value:.4f}"


def plot_mode(ax, points, lines, phi, eigenvalue, plot_scale):
    max_abs_phi = max(abs(value) for value in phi)
    if max_abs_phi == 0.0:
        scaled_phi = phi
    else:
        scaled_phi = [plot_scale * value / max_abs_phi for value in phi]

    for line in lines:
        xs = [points[i][0] for i in line]
        ys = [points[i][1] for i in line]
        zeros = [0.0 for _ in line]
        zs = [scaled_phi[i] for i in line]
        ax.plot(xs, ys, zeros, color="royalblue", linewidth=0.16, alpha=0.35, zorder=1)
        ax.plot(xs, ys, zs, color="red", linewidth=1.55, alpha=1.0, zorder=5)

    for x, y in unique_vertices(points):
        ax.scatter([x], [y], [0.0], s=8, facecolors="white",
                   edgecolors="royalblue", linewidths=0.45, alpha=0.55)

    ax.set_title(
        "Eigenfunction associated to $\\lambda$ = " + format_eigenvalue(eigenvalue),
        fontsize=8,
        pad=2,
    )
    set_equal_xy(ax, points)
    ax.set_zlim(-1.1, 1.1)
    ax.set_xlabel("x", fontsize=7, labelpad=-6)
    ax.set_ylabel("y", fontsize=7, labelpad=-6)
    ax.set_zlabel("u", fontsize=7, labelpad=-10)
    ax.xaxis.set_major_locator(MultipleLocator(1.0))
    ax.yaxis.set_major_locator(MultipleLocator(1.0))
    ax.zaxis.set_major_locator(MultipleLocator(0.5))
    ax.tick_params(axis="both", labelsize=6, pad=-2)
    ax.tick_params(axis="z", labelsize=6, pad=-2)
    ax.grid(True, linewidth=0.22, alpha=0.22)
    ax.view_init(elev=22, azim=-55)


def main():
    parser = argparse.ArgumentParser(
        description="Plot lifted eigenfunctions from FEMG VTP files."
    )
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--modes", type=int, default=6)
    parser.add_argument("--start-mode", type=int, default=0)
    parser.add_argument("--cols", type=int, default=2)
    parser.add_argument("--scale", type=float, default=1.0)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--diagnostics", type=Path, default=None)
    args = parser.parse_args()

    eigenvalues = read_eigenvalues(args.results_dir / "eigenvalues.csv")
    mode_ids = list(range(args.start_mode, args.start_mode + args.modes))

    rows = math.ceil(len(mode_ids) / args.cols)
    fig = plt.figure(figsize=(4.2 * args.cols, 3.25 * rows))

    diagnostics_rows = []
    for plot_index, mode in enumerate(mode_ids):
        points, lines, data = read_vtp(
            args.results_dir / f"eigenmode_{mode:02d}.vtp"
        )
        phi = data["phi"]
        diagnostics_rows.append({
            "mode": mode,
            "eigenvalue": eigenvalues[mode],
            "min_phi": min(phi),
            "max_phi": max(phi),
            "max_abs_phi": max(abs(value) for value in phi),
        })

        ax = fig.add_subplot(rows, args.cols, plot_index + 1, projection="3d", computed_zorder=False)
        plot_mode(ax, points, lines, phi, eigenvalues[mode], args.scale)

    fig.tight_layout()

    output = args.output
    if output is None:
        output = args.results_dir / "eigenmodes.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Wrote {output}")

    diagnostics = args.diagnostics
    if diagnostics is None:
        diagnostics = args.results_dir / "eigenmode_diagnostics.csv"
    with diagnostics.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=[
                "mode",
                "eigenvalue",
                "min_phi",
                "max_phi",
                "max_abs_phi",
            ],
        )
        writer.writeheader()
        writer.writerows(diagnostics_rows)
    print(f"Wrote {diagnostics}")


if __name__ == "__main__":
    main()