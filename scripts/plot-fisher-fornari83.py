#!/usr/bin/env python3

"""Plot the Fornari biomarker curves and compute activation times.

Figure 7 of Fornari et al. carries a small anatomical inset, the connectome
with the vertices of the four lobes in the lobe colours and every other
vertex in grey; panel (a) reproduces it. The lobe membership mirrors
``classify()`` of ``test_fisher_kolmogorov_fornari83.cpp``, the function
that computes the plotted curves, and the script asserts the same four-lobe
vertex count the test enforces, so the inset and the curves cannot drift
apart.
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figure_style
from connectome_style import load_nodes
import render_connectome as rc

# The four lobes in the colours of figure 7 of Fornari et al.
REGIONAL_CURVES = {
    "temporal": ("temporal", figure_style.LOBE_COLOUR["temporal"]),
    "frontal": ("frontal", figure_style.LOBE_COLOUR["frontal"]),
    "parietal": ("parietal", figure_style.LOBE_COLOUR["parietal"]),
    "occipital": ("occipital", figure_style.LOBE_COLOUR["occipital"]),
}
CURVE_KEYS = (*REGIONAL_CURVES, "global")

# Mirror of classify() in test_fisher_kolmogorov_fornari83.cpp, in the same
# order (the parietal check must see "precuneus" before the occipital check
# sees "cuneus"). The temporal lobe follows Fornari et al. and includes the
# entorhinal, fusiform and parahippocampal regions.
LOBE_KEYWORDS = (
    ("temporal", ("temporal", "bankssts", "entorhinal", "fusiform",
                  "parahippocampal")),
    ("frontal", ("frontal", "orbitofrontal", "parsopercularis",
                 "parsorbitalis", "parstriangularis", "precentral")),
    ("parietal", ("parietal", "postcentral", "precuneus", "supramarginal",
                  "paracentral")),
    ("occipital", ("cuneus", "occipital", "lingual", "pericalcarine")),
)
LOBE_COLOUR_ORDER = ("temporal", "frontal", "parietal", "occipital")
NON_LOBE_GREY = "#B8BCC2"


def lobe_of(name):
    lowered = name.lower()
    for lobe, keywords in LOBE_KEYWORDS:
        if any(keyword in lowered for keyword in keywords):
            return lobe
    return None


def render_inset(path):
    """The connectome with the four lobes coloured, after the inset of
    figure 7 of Fornari et al."""
    nodes = load_nodes()
    lobes = [lobe_of(node["name"]) for node in nodes]
    assert sum(lobe is not None for lobe in lobes) == 58, \
        "four-lobe vertex count disagrees with the C++ test"
    order = [None, *LOBE_COLOUR_ORDER]
    values = np.array([order.index(lobe) for lobe in lobes], dtype=float)
    palette = ListedColormap(
        [NON_LOBE_GREY,
         *(figure_style.LOBE_COLOUR[lobe] for lobe in LOBE_COLOUR_ORDER)])
    table = rc.lookup_table(palette, -0.5, len(order) - 0.5,
                            samples=len(order))
    coords = np.array([node["coords"] for node in nodes])
    rc.render(path, "sagittal_right", coords, values, table,
              node_radius=3.6, scale=rc.common_scale(),
              size=(2200, 1800), surface_opacity=0.075)
    return path


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodal", type=Path, required=True)
    parser.add_argument("--fem", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=50.0)
    return parser.parse_args()


def read_curves(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        key: [float(row[key]) for row in rows]
        for key in ("time", *CURVE_KEYS)
    }


def crossing_time(time, values, threshold):
    for index in range(1, len(time)):
        if values[index] >= threshold and values[index - 1] < threshold:
            fraction = ((threshold - values[index - 1])
                        / (values[index] - values[index - 1]))
            return time[index - 1] + fraction * (time[index] - time[index - 1])
    return None


def write_activation_times(path, datasets, threshold):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["method", "region", "threshold_percent", "time_years"])
        for method, data in datasets.items():
            for region in ("temporal", "frontal", "parietal", "occipital", "global"):
                time = crossing_time(data["time"], data[region], threshold)
                writer.writerow([
                    method, region, threshold,
                    "" if time is None else f"{time:.12g}",
                ])


def main():
    args = arguments()
    datasets = {"Nodal reference": read_curves(args.nodal),
                "Metric-graph FEM": read_curves(args.fem)}
    args.output_dir.mkdir(parents=True, exist_ok=True)

    figure_style.apply()
    figure, axes = plt.subplots(2, 2, figsize=(9.6, 6.6), sharex="col")
    with tempfile.TemporaryDirectory() as scratch:
        inset = axes[0, 0].inset_axes([0.015, 0.44, 0.40, 0.46])
        rc.show_render(inset, render_inset(Path(scratch) / "lobes.png"))
    for column, (title, data) in enumerate(datasets.items()):
        absolute = axes[0, column]
        deviation = axes[1, column]

        absolute.axhline(args.threshold, color="0.6", linewidth=0.9,
                         linestyle=(0, (4, 3)), zorder=1)
        # Draw the network average first so it cannot hide regional curves.
        absolute.plot(data["time"], data["global"], color="0.25",
                      linewidth=2.0, linestyle=(0, (5, 3)), zorder=2)
        for key, (label, color) in REGIONAL_CURVES.items():
            absolute.plot(data["time"], data[key], color=color,
                          linewidth=1.9, solid_capstyle="round", zorder=3)
            difference = [value - average for value, average
                          in zip(data[key], data["global"])]
            deviation.plot(data["time"], difference, color=color,
                           linewidth=1.9, solid_capstyle="round", zorder=3)

        deviation.axhline(0.0, color="0.6", linewidth=0.9,
                          linestyle=(0, (4, 3)), zorder=1)
        # No text is drawn inside the panels: which panel is which model,
        # and the fact that the four regional curves coincide to within the
        # line width, are stated in the caption of the report; the spread
        # itself is written to activation_times.csv and checked there.
        absolute.set_xlim(0.0, 25.0)
        absolute.set_ylim(0.0, 102.0)
        absolute.set_yticks([0, args.threshold, 100])
        absolute.set_yticklabels(["0", "50", "100"])
        # The columns share their time axis, but the upper panels keep
        # their own tick labels so no panel is left without times.
        absolute.tick_params(labelbottom=True)
        deviation.set_xlim(0.0, 25.0)
        # One common range for the two deviation panels, so the nodal and
        # the finite element deviations can be compared directly.
        deviation.set_ylim(-0.32, 0.95)
        deviation.set_yticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
        deviation.set_xticks([0, 5, 10, 15, 20, 25])
        figure_style.xname(deviation, "time [years]", y=-0.13)

    # One colour column serves the whole figure, in the empty lower-right
    # corner of panel (a): the four lobes and the dashed network average.
    entries = [*REGIONAL_CURVES.values(), ("network average", "0.25")]
    # Placed at 14 years so that the longest label, "network average",
    # ends well inside the frame of the panel.
    for slot, (label, colour) in enumerate(entries):
        figure_style.label_series(axes[0, 0], 14.0, 38.0 - 7.6 * slot,
                                  label, colour, fontsize=9)

    axes[0, 0].set_ylabel("biomarker abnormality [%]", labelpad=2)
    axes[1, 0].set_ylabel("deviation from network [pp]", labelpad=2)
    for letter, axis in zip("abcd", (axes[0, 0], axes[0, 1],
                                     axes[1, 0], axes[1, 1])):
        axis.text(0.0, 1.045, f"({letter})", transform=axis.transAxes,
                  fontsize=10.5, fontweight="bold", style="italic",
                  va="bottom")
    figure.tight_layout(w_pad=2.6, h_pad=2.2)
    figure.savefig(args.output_dir / "biomarker_comparison.png", dpi=300)
    figure.savefig(args.output_dir / "biomarker_comparison.pdf")
    plt.close(figure)

    write_activation_times(
        args.output_dir / "activation_times.csv", datasets, args.threshold)


if __name__ == "__main__":
    main()
