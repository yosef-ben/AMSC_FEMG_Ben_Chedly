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
                "time": float(row["time"]),
                "l2_squared": float(row["l2_squared"]),
                "energy": float(row["energy"]),
                "exact_l2_squared": float(row["exact_l2_squared"]),
                "exact_energy": float(row["exact_energy"]),
            })
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Plot mass norm and energy decay for a heat problem on a graph."
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()

    rows = read_rows(args.csv_file)
    time = [row["time"] for row in rows]

    plt.rcParams.update({"font.size": 11})
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))

    axes[0].plot(
        time,
        [row["l2_squared"] for row in rows],
        color="red",
        linewidth=1.8,
        label="numerical",
    )
    axes[0].plot(
        time,
        [row["exact_l2_squared"] for row in rows],
        color="black",
        linestyle="--",
        linewidth=1.3,
        label="modal reference",
    )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$u_h^T M u_h$")
    axes[0].set_title(r"$L^2(\Gamma)$ norm squared")

    axes[1].semilogy(
        time,
        [row["energy"] for row in rows],
        color="red",
        linewidth=1.8,
        label="numerical",
    )
    axes[1].semilogy(
        time,
        [row["exact_energy"] for row in rows],
        color="black",
        linestyle="--",
        linewidth=1.3,
        label="modal reference",
    )
    axes[1].set_xlabel("time")
    axes[1].set_ylabel(r"$u_h^T H u_h$")
    axes[1].set_title("Discrete energy")

    for ax in axes:
        ax.grid(True, linewidth=0.35, alpha=0.35)
        ax.legend(frameon=True)

    if args.title is not None:
        fig.suptitle(args.title, y=1.02)

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
