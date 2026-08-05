#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def read_rows(path):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        return [
            {
                "index": int(row["index"]),
                "combinatorial": float(row["combinatorial_laplacian"]),
                "metric": float(row["metric_laplacian"]),
            }
            for row in reader
        ]


def nice_max(value):
    if value <= 0.0:
        return 1.0
    magnitude = 10.0 ** (len(str(int(value))) - 1)
    return magnitude * (int(value / magnitude) + 1)


def line_points(rows, key, x0, y0, w, h, ymax):
    n = max(1, len(rows) - 1)
    pts = []
    for row in rows:
        x = x0 + w * row["index"] / n
        y = y0 + h * (1.0 - max(0.0, row[key]) / ymax)
        pts.append(f"{x:.2f},{y:.2f}")
    return " ".join(pts)


def draw_panel(svg, x0, y0, w, h, title, series):
    ymax = nice_max(max(max(max(row[key], 0.0) for row in rows) for rows, key, *_ in series))
    svg.append(f'<text x="{x0 + w/2:.1f}" y="{y0 - 12}" text-anchor="middle" font-family="Arial" font-size="15">{title}</text>')
    svg.append(f'<rect x="{x0}" y="{y0}" width="{w}" height="{h}" fill="white" stroke="black" stroke-width="1"/>')

    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = y0 + h * (1.0 - frac)
        value = ymax * frac
        svg.append(f'<line x1="{x0}" y1="{y:.2f}" x2="{x0 + w}" y2="{y:.2f}" stroke="#dddddd" stroke-width="0.8"/>')
        svg.append(f'<text x="{x0 - 8}" y="{y + 4:.2f}" text-anchor="end" font-family="Arial" font-size="10">{value:.2g}</text>')

    max_index = max(max(row["index"] for row in rows) for rows, *_ in series)
    for tick in range(0, max_index + 1, max(1, max_index // 4)):
        x = x0 + w * tick / max(1, max_index)
        svg.append(f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y0 + h}" stroke="#eeeeee" stroke-width="0.8"/>')
        svg.append(f'<text x="{x:.2f}" y="{y0 + h + 18}" text-anchor="middle" font-family="Arial" font-size="10">{tick}</text>')

    for rows, key, color, label, dash in series:
        pts = line_points(rows, key, x0, y0, w, h, ymax)
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="2"{dash_attr}/>')
        for row in rows:
            n = max(1, len(rows) - 1)
            cx = x0 + w * row["index"] / n
            cy = y0 + h * (1.0 - max(0.0, row[key]) / ymax)
            svg.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="2.8" fill="{color}" stroke="white" stroke-width="0.7"/>')

    svg.append(f'<text x="{x0 + w/2:.1f}" y="{y0 + h + 42}" text-anchor="middle" font-family="Arial" font-size="12">eigenvalue index</text>')
    svg.append(f'<text x="{x0 - 52}" y="{y0 + h/2:.1f}" transform="rotate(-90 {x0 - 52} {y0 + h/2:.1f})" text-anchor="middle" font-family="Arial" font-size="12">eigenvalue</text>')


def add_legend(svg, x, y, items):
    svg.append(f'<rect x="{x - 12}" y="{y - 18}" width="270" height="{24 * len(items) + 18}" fill="white" stroke="#bbbbbb"/>')
    for i, (label, color, dash) in enumerate(items):
        yy = y + 24 * i
        dash_attr = f' stroke-dasharray="{dash}"' if dash else ""
        svg.append(f'<line x1="{x}" y1="{yy}" x2="{x + 34}" y2="{yy}" stroke="{color}" stroke-width="2"{dash_attr}/>')
        svg.append(f'<text x="{x + 44}" y="{yy + 4}" font-family="Arial" font-size="12">{label}</text>')


def main():
    parser = argparse.ArgumentParser(description="Dependency-free SVG spectral comparison plot.")
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--graphene", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    star = read_rows(args.star)
    graphene = read_rows(args.graphene)
    tree_comb = read_rows(args.tree_dir / "tree_fixed_length_varying_angle" / "spectral_comparison.csv")
    tree_fixed = read_rows(args.tree_dir / "tree_fixed_length_varying_angle" / "spectral_comparison.csv")
    tree_inv = read_rows(args.tree_dir / "tree_angle_pi4_length_inv" / "spectral_comparison.csv")
    tree_inv2 = read_rows(args.tree_dir / "tree_angle_pi4_length_inv2" / "spectral_comparison.csv")

    width, height = 1120, 820
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">']
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="560" y="28" text-anchor="middle" font-family="Arial" font-size="19">Combinatorial and metric Laplacian spectra</text>')

    draw_panel(svg, 90, 80, 410, 250, "Four-pointed star", [
        (star, "combinatorial", "black", "combinatorial", ""),
        (star, "metric", "#d62728", "metric FEM", ""),
    ])
    draw_panel(svg, 620, 80, 410, 250, "Graphene-like graph", [
        (graphene, "combinatorial", "black", "combinatorial", ""),
        (graphene, "metric", "#d62728", "metric FEM", ""),
    ])
    draw_panel(svg, 185, 470, 620, 250, "Binary tree families", [
        (tree_comb, "combinatorial", "black", "combinatorial", ""),
        (tree_fixed, "metric", "#1f77b4", "metric: length 1, angle pi/(4n)", ""),
        (tree_inv, "metric", "#ff7f0e", "metric: length 1/n", ""),
        (tree_inv2, "metric", "#d62728", "metric: length 1/n^2", ""),
    ])

    add_legend(svg, 835, 505, [
        ("combinatorial", "black", ""),
        ("metric FEM", "#d62728", ""),
        ("tree: length 1", "#1f77b4", ""),
        ("tree: length 1/n", "#ff7f0e", ""),
        ("tree: length 1/n^2", "#d62728", ""),
    ])
    svg.append('<text x="835" y="675" font-family="Arial" font-size="12">Tree curves with the same edge lengths</text>')
    svg.append('<text x="835" y="692" font-family="Arial" font-size="12">overlap: planar angles affect the drawing,</text>')
    svg.append('<text x="835" y="709" font-family="Arial" font-size="12">not the metric Laplacian spectrum.</text>')
    svg.append("</svg>")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(svg))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
