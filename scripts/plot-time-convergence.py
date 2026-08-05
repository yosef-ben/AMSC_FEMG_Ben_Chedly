#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter


def read_rows(path):
    rows = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            rows.append({
                "method": row["method"],
                "theta": float(row["theta"]),
                "dt": float(row["dt"]),
                "L2_error": float(row["L2_error"]),
                "rate": None if row["rate"] == "" else float(row["rate"]),
            })
    return rows


def select_method(rows, method):
    return sorted(
        [row for row in rows if row["method"] == method],
        key=lambda row: row["dt"],
    )


def values(rows, key):
    return [row[key] for row in rows]


def main():
    parser = argparse.ArgumentParser(
        description="Plot time convergence for Backward Euler and Crank--Nicolson."
    )
    parser.add_argument("csv_file", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--title", default="Temporal convergence")
    args = parser.parse_args()

    rows = read_rows(args.csv_file)
    be = select_method(rows, "BE")
    cn = select_method(rows, "CN")

    be_dt = values(be, "dt")
    be_error = values(be, "L2_error")
    cn_dt = values(cn, "dt")
    cn_error = values(cn, "L2_error")

    plt.rcParams.update({"font.size": 12})
    fig, ax = plt.subplots(figsize=(6.2, 4.0))

    ax.plot(be_dt, be_error, marker="o", linewidth=1.8,
            color="tab:red", label="Backward Euler")
    ax.plot(cn_dt, cn_error, marker="s", linewidth=1.8,
            color="tab:blue", label="Crank--Nicolson")

    be_reference = [be_error[-1] * (dt / be_dt[-1]) for dt in be_dt]
    cn_reference = [cn_error[-1] * (dt / cn_dt[-1]) ** 2 for dt in cn_dt]

    ax.plot(be_dt, be_reference, "--", color="black", linewidth=1.6, alpha=0.85,
            label=r"$\Delta t$")
    ax.plot(cn_dt, cn_reference, ":", color="black", linewidth=1.8, alpha=0.85,
            label=r"$\Delta t^2$")

    ticks = sorted(set(values(rows, "dt")))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{value:g}" for value in ticks])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlabel(r"$\Delta t$")
    ax.set_ylabel(r"$L^2$ error")
    ax.set_title(args.title)
    ax.grid(True, which="both", linestyle=":", linewidth=0.45, alpha=0.75)
    ax.legend(frameon=True)
    fig.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()