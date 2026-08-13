"""House conventions for the 2D line figures.

The style follows the line figures of the reference works, as established on
the one-dimensional figures of benchmark 18: full frame, no grid, sparse
inward ticks, bold sans labels, the axis name written after the last tick,
series labelled directly in their own colour rather than in a legend box,
and no caption drawn into the image.

Every plot script applies `apply()` once and uses `xname` for the horizontal
axis name; nothing here touches the data.
"""

import matplotlib.pyplot as plt

# The four lobes of Fornari et al., in the colours of their figure 7.
LOBE_COLOUR = {
    "temporal": "#2CA02C",
    "frontal": "#D62728",
    "parietal": "#FF7F0E",
    "occipital": "#1F77B4",
}


def apply():
    """Set the shared rcParams; call once before creating figures."""
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 10, "font.weight": "bold",
        "axes.labelweight": "bold", "axes.linewidth": 1.5,
        "axes.grid": False,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.major.size": 5, "ytick.major.size": 5,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.frameon": False,
        "text.color": "black", "axes.labelcolor": "black",
        "xtick.color": "black", "ytick.color": "black",
    })


def xname(axis, text, y=-0.11, fontsize=10):
    """Write the horizontal axis name after the last tick, house style.

    Coordinates are axes fractions, so the label lands just beyond the
    lower-right corner whatever the data range is.
    """
    axis.text(1.01, y, text, fontsize=fontsize, fontweight="bold",
              ha="left", va="center", transform=axis.transAxes,
              clip_on=False)


def label_series(axis, x, y, text, colour, fontsize=9, ha="left",
                 va="center"):
    """Direct coloured label at the end of a curve, instead of a legend."""
    axis.text(x, y, text, fontsize=fontsize, fontweight="bold",
              color=colour, ha=ha, va=va, clip_on=False)
