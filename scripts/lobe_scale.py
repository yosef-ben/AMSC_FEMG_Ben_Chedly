"""The four biomarker lobes of the connectome and their transport rates.

Shared by the figure and record scripts that reason at the scale of the
lobes: the partition of the 83 regions into the four cortical lobes of
Fornari et al. (by FreeSurfer name, the rule of
test_fisher_kolmogorov_fornari83.cpp) plus the remaining regions, the
connectivity summed over that partition, and the relaxation rates of the
patterns that are constant on the groups, for the nodal model (mass matrix
the identity) and for the metric-graph FEM at one element per connection
with unit lengths (consistent or lumped mass). The slowest of these rates is
the lobe-scale eigenvalue of the report; everything is computed from the
stored node and edge lists and nothing is fitted.
"""

import csv
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

LOBE_KEYS = {
    "temporal": ("temporal", "bankssts", "entorhinal", "fusiform",
                 "parahippocampal"),
    "frontal": ("frontal", "orbitofrontal", "parsopercularis",
                "parsorbitalis", "parstriangularis", "precentral"),
    "parietal": ("parietal", "postcentral", "precuneus", "supramarginal",
                 "paracentral"),
    "occipital": ("cuneus", "occipital", "lingual", "pericalcarine"),
}
GROUPS = ("temporal", "frontal", "parietal", "occipital", "other")
MODELS = ("nodal", "consistent", "lumped")
NODES = Path("data/connectome/fornari83/nodes.csv")
EDGES = Path("data/connectome/fornari83/edges.csv")


def classify(name):
    lowered = name.lower()
    for lobe, keys in LOBE_KEYS.items():
        if any(key in lowered for key in keys):
            return lobe
    return "other"


class LobeGraph:
    """Laplacian, masses and group structure of the 83-region connectome."""

    def __init__(self, nodes=NODES, edges=EDGES):
        with open(nodes, newline="") as stream:
            rows = list(csv.DictReader(stream))
        self.names = {int(r["node_id"]): r["name"] for r in rows}
        self.group = {i: classify(name) for i, name in self.names.items()}
        self.size = 1 + max(self.names)
        n = self.size
        self.laplacian = np.zeros((n, n))
        self.lumped = np.zeros(n)
        self.consistent = np.zeros((n, n))
        self.coupling = {}
        with open(edges, newline="") as stream:
            for edge in csv.DictReader(stream):
                i, j = int(edge["source"]), int(edge["target"])
                w = float(edge["connectivity_weight"])
                self.laplacian[i, i] += w
                self.laplacian[j, j] += w
                self.laplacian[i, j] -= w
                self.laplacian[j, i] -= w
                # One P1 element of unit length per connection.
                self.lumped[i] += 0.5
                self.lumped[j] += 0.5
                self.consistent[i, i] += 1.0 / 3.0
                self.consistent[j, j] += 1.0 / 3.0
                self.consistent[i, j] += 1.0 / 6.0
                self.consistent[j, i] += 1.0 / 6.0
                a, b = self.group[i], self.group[j]
                if a != b:
                    key = tuple(sorted((a, b)))
                    self.coupling[key] = self.coupling.get(key, 0.0) + w
        self.counts = {g: sum(1 for v in self.group.values() if v == g)
                       for g in GROUPS}
        self.projector = np.zeros((n, len(GROUPS)))
        for i, g in self.group.items():
            self.projector[i, GROUPS.index(g)] = 1.0

    def mass(self, model):
        if model == "nodal":
            return np.eye(self.size)
        if model == "consistent":
            return self.consistent
        if model == "lumped":
            return np.diag(self.lumped)
        raise ValueError(model)

    def global_rate(self, model):
        """Smallest non-zero eigenvalue of L relative to the mass matrix."""
        if model == "nodal":
            return float(np.sort(np.linalg.eigvalsh(self.laplacian))[1])
        inverse = np.linalg.inv(np.linalg.cholesky(self.mass(model)))
        reduced = inverse @ self.laplacian @ inverse.T
        return float(np.sort(np.linalg.eigvalsh(reduced))[1])

    def lobe_rates(self, model):
        """Non-zero relaxation rates of the patterns constant on the five
        groups, the compressed graph of the report, for the given mass."""
        coarse_l = self.projector.T @ self.laplacian @ self.projector
        coarse_m = self.projector.T @ self.mass(model) @ self.projector
        values = np.sort(eigh(coarse_l, coarse_m, eigvals_only=True))
        return [float(v) for v in values[1:]]

    def lobe_rate(self, model):
        return self.lobe_rates(model)[0]

    def contrast_rate(self, one, two, model="nodal"):
        """Relaxation rate of the difference between the means of two
        groups, the Rayleigh quotient of that pattern."""
        n = self.size
        x = np.zeros(n)
        for i, g in self.group.items():
            if g == one:
                x[i] = 1.0 / self.counts[one]
            elif g == two:
                x[i] = -1.0 / self.counts[two]
        mass = self.mass(model)
        ones = np.ones(n)
        x -= ones * (ones @ mass @ x) / (ones @ mass @ ones)
        return float(x @ self.laplacian @ x / (x @ mass @ x))

    def fiedler_support(self, count=5):
        """The largest components of the eigenvector of the global rate."""
        values, vectors = np.linalg.eigh(self.laplacian)
        vector = vectors[:, 1]
        order = np.argsort(-np.abs(vector))[:count]
        return [(self.names[int(i)], float(vector[i])) for i in order]


def damkohler_lobe(rho, alpha, lobe_rate):
    return alpha / (rho * lobe_rate)
