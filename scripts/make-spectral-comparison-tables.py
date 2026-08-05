#!/usr/bin/env python3

import csv
from pathlib import Path


ROOT = Path("benchmarks/13_spectral_comparison/results")


def read_case(name):
    path = ROOT / name / "spectral_comparison.csv"
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def fmt(value):
    x = float(value)
    if abs(x) < 5.0e-5:
        return "0"
    return f"{x:.4f}"


def two_column_table(rows, caption, label):
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\begin{tabular}{c c}",
        r"\hline",
        r"Combinatorial Laplace matrix & Extended Laplace matrix \\",
        r"\hline",
    ]
    for row in rows:
        lines.append(
            f"{fmt(row['combinatorial_laplacian'])} & "
            f"{fmt(row['metric_laplacian'])} \\\\"
        )
    lines.extend([
        r"\hline",
        r"\end{tabular}",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        r"\end{table}",
        "",
    ])
    return lines


def tree_main_table(tree_fixed, tree_inv, tree_inv2):
    lines = [
        r"\begin{table}[htbp]",
        r"\centering",
        r"\scriptsize",
        r"\begin{tabular}{c c c c}",
        r"\hline",
        r"Comb. & $\ell_n=1$ & $\alpha=\pi/4$, $\ell_n=1/n$ & "
        r"$\alpha=\pi/4$, $\ell_n=1/n^2$ \\",
        r"\hline",
    ]
    for i in range(len(tree_fixed)):
        lines.append(
            f"{fmt(tree_fixed[i]['combinatorial_laplacian'])} & "
            f"{fmt(tree_fixed[i]['metric_laplacian'])} & "
            f"{fmt(tree_inv[i]['metric_laplacian'])} & "
            f"{fmt(tree_inv2[i]['metric_laplacian'])} \\\\"
        )
    lines.extend([
        r"\hline",
        r"\end{tabular}",
        r"\caption{Spectral comparison for the tree graph. The omitted angle "
        r"variants give the same metric spectra for fixed edge-length laws.}",
        r"\label{tab:spectral_tree_comparison_compact}",
        r"\end{table}",
        "",
    ])
    return lines



def main():
    lines = [
        "% Automatically generated from benchmarks/13_spectral_comparison/results/*.csv",
        "",
    ]
    lines += two_column_table(
        read_case("star"),
        "Spectral comparison for the four-pointed star graph.",
        "tab:spectral_star_comparison",
    )
    lines += two_column_table(
        read_case("graphene"),
        "Spectral comparison for the graphene-like graph.",
        "tab:spectral_graphene_comparison",
    )
    lines += tree_main_table(
        read_case("tree_fixed_length_varying_angle"),
        read_case("tree_angle_pi4_length_inv"),
        read_case("tree_angle_pi4_length_inv2"),
    )

    output = ROOT / "spectral_comparison_tables.tex"
    output.write_text("\n".join(lines))
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
