#!/usr/bin/env python3

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import NullFormatter


convergence_data = pd.read_csv(sys.argv[1], sep=",")
output_path = Path("output/convergence/star_sine_space_convergence.pdf")
output_path.parent.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({"font.size": 14})

plt.plot(
    convergence_data.h,
    convergence_data.L2,
    marker="o",
    label="L2",
)

reference = convergence_data.L2.iloc[0] * (
    convergence_data.h / convergence_data.h.iloc[0]
) ** 2

plt.plot(
    convergence_data.h,
    reference,
    "--",
    label=r"$h^2$",
)

plt.xscale("log")
plt.yscale("log")
plt.xticks(
    convergence_data.h,
    [f"{value:g}" for value in convergence_data.h],
)
plt.gca().xaxis.set_minor_formatter(NullFormatter())
plt.xlabel("h")
plt.ylabel("error")
plt.legend()
plt.grid(True, which="both", linestyle=":")
plt.tight_layout(pad=1.2)

plt.savefig(output_path)
