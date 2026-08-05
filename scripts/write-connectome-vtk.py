#!/usr/bin/env python3

"""Write unique-node VTK time series for the 83-region connectome."""

import argparse
import csv
from pathlib import Path


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--stride", type=int, default=5)
    return parser.parse_args()


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_vtp(path, nodes, edges, concentration, points_only):
    lines = [] if points_only else edges
    with path.open("w", encoding="ascii") as output:
        output.write('<?xml version="1.0"?>\n')
        output.write('<VTKFile type="PolyData" version="0.1" '
                     'byte_order="LittleEndian">\n<PolyData>\n')
        number_of_vertices = len(nodes) if points_only else 0
        output.write(f'<Piece NumberOfPoints="{len(nodes)}" '
                     f'NumberOfVerts="{number_of_vertices}" '
                     f'NumberOfLines="{len(lines)}">\n')
        output.write('<Points><DataArray type="Float64" '
                     'NumberOfComponents="3" format="ascii">\n')
        for node in nodes:
            output.write(f'{node["x"]} {node["y"]} {node["z"]}\n')
        output.write('</DataArray></Points>\n<Verts>\n')
        output.write('<DataArray type="Int64" Name="connectivity" '
                     'format="ascii">\n')
        if points_only:
            output.write(" ".join(str(index) for index in range(len(nodes))))
        output.write("\n")
        output.write('</DataArray><DataArray type="Int64" Name="offsets" '
                     'format="ascii">\n')
        if points_only:
            output.write(" ".join(
                str(index) for index in range(1, len(nodes) + 1)))
        output.write("\n</DataArray></Verts>\n")
        output.write('<Lines><DataArray type="Int64" Name="connectivity" '
                     'format="ascii">\n')
        for edge in lines:
            output.write(f'{edge["source"]} {edge["target"]}\n')
        output.write('</DataArray><DataArray type="Int64" Name="offsets" '
                     'format="ascii">\n')
        output.write(" ".join(str(2 * index) for index in range(1, len(lines) + 1))
                     + "\n</DataArray></Lines>\n")
        output.write('<PointData Scalars="c"><DataArray type="Float64" '
                     'Name="c" format="ascii">\n')
        output.write(" ".join(f"{value:.16g}" for value in concentration))
        output.write('\n</DataArray><DataArray type="Int32" Name="lobe_id" '
                     'format="ascii">\n')
        output.write(" ".join(node["lobe_id"] for node in nodes))
        output.write("\n</DataArray></PointData>\n")
        if lines:
            output.write('<CellData Scalars="connectivity_weight"><DataArray '
                         'type="Float64" Name="connectivity_weight" '
                         'format="ascii">\n')
            output.write(" ".join(edge["connectivity_weight"] for edge in lines))
            output.write("\n</DataArray></CellData>\n")
        output.write('</Piece></PolyData></VTKFile>\n')


def write_pvd(path, prefix, kind, records):
    with path.open("w", encoding="ascii") as output:
        output.write('<?xml version="1.0"?>\n<VTKFile type="Collection" '
                     'version="0.1" byte_order="LittleEndian">\n<Collection>\n')
        for index, time in records:
            output.write(f'<DataSet timestep="{time:.16g}" group="" part="0" '
                         f'file="{prefix}_{kind}_{index:04d}.vtp"/>\n')
        output.write('</Collection></VTKFile>\n')


def classify(name):
    name = name.lower()
    if any(token in name for token in
           ("temporal", "bankssts", "entorhinal", "fusiform", "parahippocampal")):
        return "0"
    if any(token in name for token in
           ("frontal", "orbitofrontal", "parsopercularis", "parsorbitalis",
            "parstriangularis", "precentral")):
        return "1"
    if any(token in name for token in
           ("parietal", "postcentral", "precuneus", "supramarginal", "paracentral")):
        return "2"
    if any(token in name for token in
           ("cuneus", "occipital", "lingual", "pericalcarine")):
        return "3"
    return "-1"


def main():
    args = arguments()
    nodes = read_csv(args.nodes)
    edges = read_csv(args.edges)
    profiles = read_csv(args.profiles)
    for node in nodes:
        node["lobe_id"] = classify(node["name"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for index, profile in enumerate(profiles):
        if index % args.stride != 0 and index != len(profiles) - 1:
            continue
        concentration = [float(profile[f"node_{node}"])
                         for node in range(len(nodes))]
        write_vtp(args.output_dir / f"{args.prefix}_nodes_{index:04d}.vtp",
                  nodes, edges, concentration, True)
        write_vtp(args.output_dir / f"{args.prefix}_edges_{index:04d}.vtp",
                  nodes, edges, concentration, False)
        records.append((index, float(profile["time"])))
    write_pvd(args.output_dir / f"{args.prefix}_nodes.pvd",
              args.prefix, "nodes", records)
    write_pvd(args.output_dir / f"{args.prefix}_edges.pvd",
              args.prefix, "edges", records)


if __name__ == "__main__":
    main()
