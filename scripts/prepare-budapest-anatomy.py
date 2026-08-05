#!/usr/bin/env python3

"""Prepare Budapest anatomical coordinates and a ParaView brain surface."""

import argparse
import csv
import json
import struct
import urllib.request
from pathlib import Path


BASE_URL = "https://pitgroup.org/apps/connectome/render/data"
SOURCE_FILES = ("nodes.json", "lh.pial", "rh.pial")


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--node-metadata", type=Path,
                        default=Path("data/connectome/budapest_lcc_nodes.csv"))
    parser.add_argument("--source-dir", type=Path,
                        default=Path("data/connectome/anatomy"))
    parser.add_argument("--download", action="store_true",
                        help="download the public viewer assets when absent")
    return parser.parse_args()


def download_assets(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for name in SOURCE_FILES:
        path = directory / name
        if not path.exists():
            print(f"Downloading {name} ...")
            urllib.request.urlretrieve(f"{BASE_URL}/{name}", path)


def read_surface(path):
    data = path.read_bytes()
    if len(data) < 16 or data[:3] != b"\xff\xff\xfe":
        raise ValueError(f"Unsupported FreeSurfer surface: {path}")
    offset = data.index(b"\n", 3) + 2
    n_vertices, n_faces = struct.unpack_from(">II", data, offset)
    offset += 8
    coordinate_size = 12 * n_vertices
    face_size = 12 * n_faces
    if offset + coordinate_size + face_size > len(data):
        raise ValueError(f"Truncated FreeSurfer surface: {path}")
    coordinates = data[offset:offset + coordinate_size]
    offset += coordinate_size
    faces = data[offset:offset + face_size]
    return n_vertices, n_faces, coordinates, faces


def write_brain_surface(source_dir):
    left = read_surface(source_dir / "lh.pial")
    right = read_surface(source_dir / "rh.pial")
    output = source_dir / "brain_surface.vtk"
    n_vertices = left[0] + right[0]
    n_faces = left[1] + right[1]

    with output.open("wb") as stream:
        stream.write(b"# vtk DataFile Version 3.0\n")
        stream.write(b"Budapest Reference Connectome pial surface\n")
        stream.write(b"BINARY\nDATASET POLYDATA\n")
        stream.write(f"POINTS {n_vertices} float\n".encode("ascii"))
        stream.write(left[2])
        stream.write(right[2])
        stream.write(f"\nPOLYGONS {n_faces} {4 * n_faces}\n".encode("ascii"))
        for hemisphere, vertex_offset in ((left, 0), (right, left[0])):
            faces = hemisphere[3]
            for face in range(hemisphere[1]):
                indices = struct.unpack_from(">III", faces, 12 * face)
                stream.write(struct.pack(
                    ">IIII", 3, *(index + vertex_offset for index in indices)))
        stream.write(b"\n")
    return output


def write_anatomical_coordinates(node_metadata, source_dir):
    with (source_dir / "nodes.json").open(encoding="utf-8") as stream:
        site_nodes = json.load(stream)
    coordinates = {
        item["dn_fsname"]: (
            float(item["pial_x"]), float(item["pial_y"]), float(item["pial_z"]))
        for item in site_nodes.values()
    }

    with node_metadata.open(newline="", encoding="utf-8") as stream:
        nodes = list(csv.DictReader(stream))
    missing = sorted({node["parent_name"] for node in nodes} - coordinates.keys())
    if missing:
        raise ValueError(f"Missing anatomical coordinates for: {missing}")

    output = node_metadata.parent / "budapest_lcc_anatomical_coordinates.csv"
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("local_id", "parent_name", "x", "y", "z"))
        for node in nodes:
            x, y, z = coordinates[node["parent_name"]]
            writer.writerow((node["local_id"], node["parent_name"],
                             f"{x:.12g}", f"{y:.12g}", f"{z:.12g}"))
    return output, len(nodes), len({node["parent_name"] for node in nodes})


def main():
    args = arguments()
    if args.download:
        download_assets(args.source_dir)
    missing = [name for name in SOURCE_FILES
               if not (args.source_dir / name).exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing {missing}; rerun with --download or provide the assets")

    coordinates, n_nodes, n_locations = write_anatomical_coordinates(
        args.node_metadata, args.source_dir)
    surface = write_brain_surface(args.source_dir)
    print(f"Mapped {n_nodes} graph nodes to {n_locations} anatomical locations.")
    print(f"Coordinates: {coordinates}")
    print(f"Brain surface: {surface}")


if __name__ == "__main__":
    main()
