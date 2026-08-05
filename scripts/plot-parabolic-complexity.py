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
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Plot timing data for the parabolic solver."
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = read_rows(args.csv_file)
    dofs = [row["n_dofs"] for row in rows]
    total = [row["total_seconds"] for row in rows]
    assembly = [row["assembly_seconds"] for row in rows]
    factorization = [row["factorization_seconds"] for row in rows]
    stepping = [row["time_stepping_seconds"] for row in rows]

    reference = [total[0] * (n / dofs[0]) for n in dofs]

    plt.rcParams.update({"font.size": 11})
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(dofs, total, marker="o", linewidth=1.8, color="black",
            label="total")
    ax.plot(dofs, assembly, marker="^", linewidth=1.4, color="tab:green",
            label="assembly")
    ax.plot(dofs, factorization, marker="s", linewidth=1.4, color="tab:red",
            label="factorization")
    ax.plot(dofs, stepping, marker="d", linewidth=1.4, color="tab:blue",
            label="time stepping")
    ax.plot(dofs, reference, "--", linewidth=1.4, color="gray",
            label=r"$N_{dof}$")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"number of degrees of freedom")
    ax.set_ylabel("time [s]")
    ax.set_title("Parabolic solver timing on the graphene-like graph")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(frameon=True)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()