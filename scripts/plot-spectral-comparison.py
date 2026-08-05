#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path):
    rows = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            rows.append({
                "index": int(row["index"]),
                "combinatorial": float(row["combinatorial_laplacian"]),
                "metric": float(row["metric_laplacian"]),
            })
    return rows


def plot_basic_case(ax, path, title):
    rows = read_rows(path)
    index = [row["index"] + 1 for row in rows]
    combinatorial = [row["combinatorial"] for row in rows]
    metric = [row["metric"] for row in rows]

    ax.plot(index, combinatorial, marker="o", markersize=4.0,
            linewidth=2.0, color="black", label="combinatorial")
    ax.plot(index, metric, marker="s", markersize=4.0,
            linewidth=2.2, color="tab:red", label="extended")
    ax.set_title(title)
    ax.set_xlabel("eigenvalue index")
    ax.set_ylabel("eigenvalue")
    ax.set_xlim(index[0] - 0.2, index[-1] + 0.2)
    ax.set_xticks(index)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(loc="upper left", frameon=True, fontsize=8)


def plot_tree_cases(ax, tree_dir):
    cases = [
        ("tree_fixed_length_varying_angle", "length 1, decreasing angle", "tab:blue", "o"),
        ("tree_angle_pi4_length_inv", "length 1/n, angle variants overlap", "tab:orange", "s"),
        ("tree_angle_pi4_length_inv2", "length 1/n^2, angle variants overlap", "tab:red", "d"),
    ]

    first_rows = read_rows(tree_dir / cases[0][0] / "spectral_comparison.csv")
    index = [row["index"] + 1 for row in first_rows]
    combinatorial = [row["combinatorial"] for row in first_rows]
    ax.plot(index, combinatorial, marker="o", markersize=3.8,
            linewidth=2.0, color="black", label="combinatorial")

    for directory, label, color, marker in cases:
        rows = read_rows(tree_dir / directory / "spectral_comparison.csv")
        index = [row["index"] + 1 for row in rows]
        metric = [row["metric"] for row in rows]
        ax.plot(index, metric, marker=marker, markersize=4.0,
                linewidth=2.0, color=color, label=label)

    ax.set_title("Binary tree families")
    ax.set_xlabel("eigenvalue index")
    ax.set_ylabel("eigenvalue")
    ax.set_xlim(index[0] - 0.5, index[-1] + 0.5)
    ax.set_xticks(index)
    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(loc="upper left", frameon=True, fontsize=7)


def main():
    parser = argparse.ArgumentParser(
        description="Plot combinatorial vs metric Laplacian spectra."
    )
    parser.add_argument("--star", type=Path, required=True)
    parser.add_argument("--graphene", type=Path, required=True)
    parser.add_argument("--tree-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    plt.rcParams.update({"font.size": 10})
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.6))
    plot_basic_case(axes[0], args.star, "Four-pointed star")
    plot_basic_case(axes[1], args.graphene, "Graphene-like graph")
    plot_tree_cases(axes[2], args.tree_dir)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
