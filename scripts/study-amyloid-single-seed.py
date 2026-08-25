#!/usr/bin/env python3

"""Does the regional reaction help the amyloid seeding when the target is not
already seeded?

The amyloid experiment of benchmark 27 seeds the 54 neocortical vertices at
once, so the neocortex, which is the first phase of the expected progression,
is filled from the start and only the descent to the deeper structures is
left for the model to produce. With a uniform rate the seeded mantle then
rises as a block, which is exactly the reference picture, and any regional
rate necessarily breaks that. This study removes the coincidence: each of the
54 neocortical vertices is seeded alone, so the model has to spread through
the neocortex before reaching the allocortex and the deep structures, and the
run is repeated with the uniform rate and with the regional field. For every
seed the four phases of the expected progression, neocortex, allocortex with
the insula, deep nuclei and brainstem, are ordered by the mean of the first
crossing of c = 0.5 over their vertices, and the number of the six pairwise
orderings that come out right is recorded, together with the margin between
the last cortical phase and the first one after it.
"""

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path

NODES = Path("data/connectome/fornari83/nodes.csv")
REGIONS = Path("benchmarks/21_fisher_kolmogorov_corti83/results"
               "/reaction_coefficients.csv")
LOBE_WORDS = ("frontal", "temporal", "parietal", "occipital")
PHASES = ("neocortex", "allocortex", "deep", "brainstem")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, default=Path(
        "build-release/test_fisher_kolmogorov_fornari83"))
    parser.add_argument("--rates", type=Path, default=REGIONS)
    parser.add_argument("--alpha", type=float, default=0.5)
    parser.add_argument("--dt", type=float, default=0.4)
    parser.add_argument("--final-time", type=float, default=120.0)
    parser.add_argument("--diffusion-scaling", type=float, default=0.005)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_csv(path):
    with open(path, newline="") as stream:
        return list(csv.DictReader(stream))


def classify(name, group):
    """The four phases of the expected amyloid progression."""
    lowered = name.lower()
    if "brain-stem" in lowered:
        return "brainstem"
    allocortex = ("entorhinal" in lowered or "parahippocampal" in lowered
                  or group == "insular")
    if allocortex:
        return "allocortex"
    if group == "subcortical":
        return "deep"
    if group in LOBE_WORDS:
        return "neocortex"
    return "allocortex"          # the remaining limbic belt


def crossings(path, level=0.5):
    rows = read_csv(path)
    times = [float(row["time"]) for row in rows]
    out = {}
    for k in range(83):
        series = [float(row[f"node_{k}"]) for row in rows]
        out[k] = next((t for t, v in zip(times, series) if v >= level), None)
    return out


def main():
    args = arguments()
    names = {int(r["node_id"]): r["name"] for r in read_csv(NODES)}
    group = {int(r["node_id"]): r["region"] for r in read_csv(REGIONS)}
    phase = {k: classify(names[k], group[k]) for k in names}
    seeds = [k for k in names if phase[k] == "neocortex"]
    members = {p: [k for k in names if phase[k] == p] for p in PHASES}

    rows = []
    with tempfile.TemporaryDirectory() as scratch:
        for field in ("uniform", "regional"):
            for seed in seeds:
                target = Path(scratch) / f"{field}_{seed}"
                command = [str(args.executable), "1", str(args.dt),
                           str(target), str(args.alpha),
                           str(args.final_time),
                           str(args.diffusion_scaling), "be_lumped",
                           str(seed)]
                if field == "regional":
                    command.append(str(args.rates))
                subprocess.run(command, check=True,
                               stdout=subprocess.DEVNULL)
                first = crossings(target / "fem_profiles.csv")
                if any(first[k] is None for k in names):
                    rows.append({"field": field, "seed": names[seed],
                                 "correct_pairs": "", "margin_years": "",
                                 **{p: "" for p in PHASES}})
                    continue
                mean = {p: sum(first[k] for k in members[p]) / len(members[p])
                        for p in PHASES}
                correct = sum(1 for i, a in enumerate(PHASES)
                              for b in PHASES[i + 1:] if mean[a] < mean[b])
                rows.append({
                    "field": field, "seed": names[seed],
                    "correct_pairs": correct,
                    "margin_years": f"{mean['allocortex'] - mean['neocortex']:.4f}",
                    **{p: f"{mean[p]:.4f}" for p in PHASES}})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as stream:
        writer = csv.DictWriter(stream, ["field", "seed", "correct_pairs",
                                         "margin_years", *PHASES])
        writer.writeheader()
        writer.writerows(rows)

    for field in ("uniform", "regional"):
        chosen = [r for r in rows if r["field"] == field and r["correct_pairs"] != ""]
        pairs = [int(r["correct_pairs"]) for r in chosen]
        margins = [float(r["margin_years"]) for r in chosen]
        perfect = sum(1 for p in pairs if p == 6)
        print(f"{field:9s}: {len(chosen)} seeds, {perfect} with all six "
              f"orderings right, mean {sum(pairs) / len(pairs):.2f} of 6, "
              f"mean margin {sum(margins) / len(margins):+.2f} yr")
    print(f"Written {args.output}")


if __name__ == "__main__":
    main()
