#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


def arguments():
    parser = argparse.ArgumentParser(
        description="Plot regional Fisher-Kolmogorov concentrations."
    )
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("regional_averages.png"))
    return parser.parse_args()


def main():
    args = arguments()
    with args.input.open(newline="") as stream:
        rows = list(csv.DictReader(stream))

    time = [float(row["time"]) for row in rows]
    names = [
        "frontal", "temporal", "parietal", "insular", "limbic",
        "occipital", "subcortical",
    ]
    labels = {
        "frontal": "Frontal lobe",
        "temporal": "Temporal lobe",
        "parietal": "Parietal lobe",
        "insular": "Insular lobe",
        "limbic": "Limbic lobe",
        "occipital": "Occipital lobe",
        "subcortical": "Subcortical nuclei",
    }
    colors = ["#d62728", "#1f77b4", "#2ca02c", "#9467bd",
              "#ff7f0e", "#17becf", "#8c564b"]

    fig, axis = plt.subplots(figsize=(9.0, 5.4))
    for name, color in zip(names, colors):
        axis.plot(time, [float(row[name]) for row in rows],
                  color=color, linewidth=1.8, label=labels[name])
    axis.plot(time, [float(row["global"]) for row in rows],
              color="black", linewidth=2.2, linestyle="--", label="Global mean")
    axis.set_xlabel("Time [years]")
    axis.set_ylabel("Mean relative concentration")
    axis.set_xlim(time[0], time[-1])
    axis.set_ylim(bottom=0.0)
    axis.grid(True, which="both", linestyle=":", linewidth=0.7, alpha=0.65)
    axis.legend(ncol=2, fontsize=9, frameon=True)
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=220)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
