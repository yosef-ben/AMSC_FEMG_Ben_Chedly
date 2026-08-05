#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_decay(path):
    rows = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            rows.append({
                "time": float(row["time"]),
                "numerical_l2": float(row["numerical_l2"]),
                "exact_l2": float(row["exact_l2"]),
            })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Plot numerical and exact L2 decay for a heat eigenmode test."
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    rows = read_decay(args.csv_file)
    time = [row["time"] for row in rows]
    numerical = [row["numerical_l2"] for row in rows]
    exact = [row["exact_l2"] for row in rows]

    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.plot(time, numerical, color="red", linewidth=1.8, label="numerical $L^2$ norm")
    ax.plot(time, exact, color="black", linestyle="--", linewidth=1.4, label="exact decay")
    ax.scatter([time[0]], [numerical[0]], color="royalblue", s=28, zorder=5,
               label="initial normalization")

    ax.set_xlabel("time")
    ax.set_ylabel(r"$\|u_h(t)\|_{L^2(\Gamma)}$")
    if args.title is not None:
        ax.set_title(args.title)
    ax.grid(True, linewidth=0.35, alpha=0.35)
    ax.legend(frameon=True)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()