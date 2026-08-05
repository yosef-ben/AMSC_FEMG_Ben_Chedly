#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(
        description="Plot linear solver timings on the graphene graph."
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

    colors = {
        "CG": "tab:blue",
        "GMRES": "tab:orange",
        "BiCGSTAB": "tab:purple",
        "LU": "tab:red",
    }
    markers = {
        "CG": "o",
        "GMRES": "s",
        "BiCGSTAB": "D",
        "LU": "^",
    }

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(6.2, 4.8))

    reference_dofs = None
    reference_value = None
    for method in methods:
        method_rows = [row for row in rows if row["method"] == method]
        dofs = [float(row["n_dofs"]) for row in method_rows]
        seconds = [float(row["seconds"]) for row in method_rows]
        ax.plot(
            dofs,
            seconds,
            marker=markers.get(method, "o"),
            linewidth=1.7,
            color=colors.get(method),
            label=method,
        )
        if reference_dofs is None and method == "CG":
            reference_dofs = dofs
            reference_value = seconds[0]

    if reference_dofs is None:
        method_rows = [row for row in rows if row["method"] == methods[0]]
        reference_dofs = [float(row["n_dofs"]) for row in method_rows]
        reference_value = float(method_rows[0]["seconds"])

    n0 = reference_dofs[0]
    references = [
        (1, "--", "0.65", r"$N_{dof}$"),
        (2, "-.", "0.45", r"$N_{dof}^2$"),
        (3, ":", "0.25", r"$N_{dof}^3$"),
    ]
    for power, style, color, label in references:
        ax.plot(
            reference_dofs,
            [reference_value * (n / n0) ** power for n in reference_dofs],
            style,
            linewidth=1.2,
            color=color,
            label=label,
        )

    ax.set_title(r"Linear solvers for $H+M$")
    ax.set_xlabel("number of degrees of freedom")
    ax.set_ylabel("time [s]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", linestyle=":", linewidth=0.4, alpha=0.7)
    ax.legend(frameon=True)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
