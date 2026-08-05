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


def read_table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def values(rows, key):
    return [row[key] for row in rows]


def plot_spectral(ax, csv_file):
    rows = read_rows(csv_file)
    dofs = values(rows, "n_dofs")
    solve = values(rows, "solve_seconds")
    assembly = values(rows, "assembly_seconds")
    total = values(rows, "total_seconds")
    reference = [solve[-1] * (n / dofs[-1]) ** 3 for n in dofs]

    ax.plot(dofs, total, marker="o", linewidth=1.8, color="black", label="total")
    ax.plot(dofs, solve, marker="s", linewidth=1.6, color="tab:red", label="eigen solve")
    ax.plot(dofs, assembly, marker="^", linewidth=1.4, color="tab:green", label="assembly")
    ax.plot(dofs, reference, "--", linewidth=1.4, color="gray", label=r"$N_{dof}^3$")
    ax.set_title("Generalized eigenvalue solver")
    ax.set_xlabel(r"number of degrees of freedom")
    ax.set_ylabel("time [s]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(frameon=True)


def plot_parabolic(ax, csv_file):
    rows = read_rows(csv_file)
    dofs = values(rows, "n_dofs")
    total = values(rows, "total_seconds")
    assembly = values(rows, "assembly_seconds")
    factorization = values(rows, "factorization_seconds")
    stepping = values(rows, "time_stepping_seconds")
    reference = [total[0] * (n / dofs[0]) for n in dofs]

    ax.plot(dofs, total, marker="o", linewidth=1.8, color="black", label="total")
    ax.plot(dofs, assembly, marker="^", linewidth=1.4, color="tab:green", label="assembly")
    ax.plot(dofs, factorization, marker="s", linewidth=1.4, color="tab:red", label="factorization")
    ax.plot(dofs, stepping, marker="d", linewidth=1.4, color="tab:blue", label="time stepping")
    ax.plot(dofs, reference, "--", linewidth=1.4, color="gray", label=r"$N_{dof}$")
    ax.set_title(r"Parabolic solver, $\theta=1/2$")
    ax.set_xlabel(r"number of degrees of freedom")
    ax.set_ylabel("time [s]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(frameon=True)


def main():
    parser = argparse.ArgumentParser(
        description="Plot spectral and parabolic computational complexity."
    )
    parser.add_argument("--spectral", type=Path, required=True)
    parser.add_argument("--parabolic", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    plt.rcParams.update({"font.size": 11})
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4))
    plot_spectral(axes[0], args.spectral)
    plot_parabolic(axes[1], args.parabolic)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()