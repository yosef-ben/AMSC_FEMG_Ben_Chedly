#!/usr/bin/env python3

import argparse
import csv
import math
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def log_bounds(values):
    logs = [math.log10(v) for v in values if v > 0.0]
    lo = math.floor(min(logs))
    hi = math.ceil(max(logs))
    if lo == hi:
        hi += 1
    return lo, hi


def polyline(points):
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def main():
    parser = argparse.ArgumentParser(
        description="Dependency-free SVG plot for linear solver timings."
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_rows(args.input)
    methods = []
    for row in rows:
        method = row["method"]
        if method not in methods:
            methods.append(method)

    width = 820
    height = 560
    left = 92
    right = 28
    top = 45
    bottom = 82
    plot_w = width - left - right
    plot_h = height - top - bottom

    dofs_all = [float(row["n_dofs"]) for row in rows]
    times_all = [float(row["seconds"]) for row in rows]
    cg_rows_for_bounds = [row for row in rows if row["method"] == "CG"]
    reference_dofs_for_bounds = [float(row["n_dofs"]) for row in cg_rows_for_bounds]
    reference_value_for_bounds = float(cg_rows_for_bounds[0]["seconds"])
    n0_for_bounds = reference_dofs_for_bounds[0]
    for power in (1, 2, 3):
        times_all.extend(
            reference_value_for_bounds * (n / n0_for_bounds) ** power
            for n in reference_dofs_for_bounds
        )
    x_lo, x_hi = log_bounds(dofs_all)
    y_lo, y_hi = log_bounds(times_all)

    def xmap(value):
        return left + (math.log10(value) - x_lo) / (x_hi - x_lo) * plot_w

    def ymap(value):
        return top + (y_hi - math.log10(value)) / (y_hi - y_lo) * plot_h

    colors = {
        "CG": "#1f77b4",
        "GMRES": "#ff7f0e",
        "BiCGSTAB": "#9467bd",
        "LU": "#d62728",
    }

    svg = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append(f'<text x="{width/2:.1f}" y="24" text-anchor="middle" font-family="Arial" font-size="17">Linear solvers for H+M</text>')

    for exponent in range(x_lo, x_hi + 1):
        x = xmap(10.0 ** exponent)
        svg.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#dddddd" stroke-width="1"/>')
        svg.append(f'<text x="{x:.2f}" y="{top + plot_h + 28}" text-anchor="middle" font-family="Arial" font-size="12">10^{exponent}</text>')

    for exponent in range(y_lo, y_hi + 1):
        y = ymap(10.0 ** exponent)
        svg.append(f'<line x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}" stroke="#dddddd" stroke-width="1"/>')
        svg.append(f'<text x="{left - 14}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="12">10^{exponent}</text>')

    svg.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="black" stroke-width="1.2"/>')

    cg_rows = [row for row in rows if row["method"] == "CG"]
    ref_dofs = [float(row["n_dofs"]) for row in cg_rows]
    ref_value = float(cg_rows[0]["seconds"])
    n0 = ref_dofs[0]
    refs = [
        (1, "#999999", "6,5", "N_dof"),
        (2, "#666666", "12,5,3,5", "N_dof^2"),
        (3, "#333333", "2,5", "N_dof^3"),
    ]
    for power, color, dash, label in refs:
        pts = [(xmap(n), ymap(ref_value * (n / n0) ** power)) for n in ref_dofs]
        svg.append(f'<polyline points="{polyline(pts)}" fill="none" stroke="{color}" stroke-width="1.7" stroke-dasharray="{dash}"/>')

    for method in methods:
        method_rows = [row for row in rows if row["method"] == method]
        pts = [(xmap(float(row["n_dofs"])), ymap(float(row["seconds"]))) for row in method_rows]
        color = colors.get(method, "black")
        svg.append(f'<polyline points="{polyline(pts)}" fill="none" stroke="{color}" stroke-width="2.4"/>')
        for x, y in pts:
            svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.2" fill="{color}" stroke="white" stroke-width="1"/>')

    svg.append(f'<text x="{left + plot_w/2:.1f}" y="{height - 20}" text-anchor="middle" font-family="Arial" font-size="14">number of degrees of freedom</text>')
    svg.append(f'<text x="22" y="{top + plot_h/2:.1f}" text-anchor="middle" transform="rotate(-90 22 {top + plot_h/2:.1f})" font-family="Arial" font-size="14">time [s]</text>')

    legend_x = left + plot_w - 150
    legend_y = top + 18
    svg.append(f'<rect x="{legend_x - 12}" y="{legend_y - 16}" width="150" height="168" fill="white" stroke="#bbbbbb"/>')
    for i, method in enumerate(methods):
        y = legend_y + i * 22
        color = colors.get(method, "black")
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 28}" y2="{y}" stroke="{color}" stroke-width="2.4"/>')
        svg.append(f'<circle cx="{legend_x + 14}" cy="{y}" r="4" fill="{color}" stroke="white" stroke-width="1"/>')
        svg.append(f'<text x="{legend_x + 38}" y="{y + 4}" font-family="Arial" font-size="12">{method}</text>')
    for j, (_, color, dash, label) in enumerate(refs):
        y = legend_y + (len(methods) + j) * 22
        svg.append(f'<line x1="{legend_x}" y1="{y}" x2="{legend_x + 28}" y2="{y}" stroke="{color}" stroke-width="1.7" stroke-dasharray="{dash}"/>')
        svg.append(f'<text x="{legend_x + 38}" y="{y + 4}" font-family="Arial" font-size="12">{label}</text>')

    svg.append("</svg>")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(svg))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
