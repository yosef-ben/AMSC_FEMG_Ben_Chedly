#!/usr/bin/env python3

"""Expected against obtained: staged spreading for two clinical seedings.

The construction of figure 1, panel (c), of Fornari et al. and of figure 7 of
Weickenmeier et al.: the network solution on a translucent brain, spheres
whose radius and blue intensity both grow with the local concentration, three
stages from left to right. The top row seeds the entorhinal cortex, where tau
inclusions first appear; the bottom row seeds the four cortical lobes, where
amyloid-beta deposits first appear (Weickenmeier et al., section 2.7).

The solution shown is the nodal network model at the scaling rho = 0.005
that benchmark 23 identifies as reproducing the published lobe separation
(Da of order one hundred); at this scaling the consistent-mass metric-graph
FEM leaves the physical range and cannot be shown, which is the boundary
documented in the report. The nodal solution stays in [0,1] throughout and
the script asserts as much. The three stages are selected by a rule, not by
eye: the first instants at which the network mean reaches 10, 40 and 80
percent. Nothing is smoothed, clipped or rescaled.
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
import vtk
from vtk.util import numpy_support

sys.path.insert(0, str(Path(__file__).resolve().parent))
from connectome_style import load_nodes
import render_connectome as rc

CASES = (
    ("tau", "tau seeding", "entorhinal cortex"),
    ("amyloid", r"amyloid-$\beta$ seeding", "neocortex"),
)
STAGES = (0.10, 0.40, 0.80)
RADIUS_MIN, RADIUS_MAX = 0.55, 5.2     # mm; radius grows with concentration
CONCENTRATION_MAP = plt.cm.Blues


def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tau", type=Path, required=True,
                        help="nodal_profiles.csv of the entorhinal seeding")
    parser.add_argument("--amyloid", type=Path, required=True,
                        help="nodal_profiles.csv of the neocortical seeding")
    parser.add_argument("--reference-dir", type=Path, default=None,
                        help="directory with the weickenmeier_fig1_*.png "
                             "strips cut by extract-weickenmeier-staging.py; "
                             "when given, the expected-against-obtained "
                             "composite is written as well")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_profiles(path):
    with open(path, newline="") as stream:
        rows = list(csv.DictReader(stream))
    times = np.array([float(row["time"]) for row in rows])
    values = np.array([[float(row[f"node_{k}"]) for k in range(83)]
                       for row in rows])
    return times, values


def _staged_actor(coords, concentration, table):
    """Spheres whose radius and colour both carry the concentration."""
    points = vtk.vtkPoints()
    for point in coords:
        points.InsertNextPoint(*point)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    radii = RADIUS_MIN + (RADIUS_MAX - RADIUS_MIN) * concentration
    scalars = numpy_support.numpy_to_vtk(np.asarray(radii, float), deep=True)
    scalars.SetName("radius")
    polydata.GetPointData().SetScalars(scalars)
    colours = numpy_support.numpy_to_vtk(
        np.asarray(concentration, float), deep=True)
    colours.SetName("concentration")
    polydata.GetPointData().AddArray(colours)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1.0)
    sphere.SetThetaResolution(30)
    sphere.SetPhiResolution(30)
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToScaleByScalar()
    glyph.SetInputArrayToProcess(0, 0, 0, 0, "radius")
    glyph.SetColorModeToColorByScalar()
    glyph.SetInputArrayToProcess(3, 0, 0, 0, "concentration")
    glyph.SetScaleFactor(1.0)
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SelectColorArray("concentration")
    mapper.SetScalarModeToUsePointFieldData()
    mapper.SetScalarRange(0.0, 1.0)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.35)
    actor.GetProperty().SetSpecularPower(24)
    actor.GetProperty().SetDiffuse(0.80)
    actor.GetProperty().SetAmbient(0.30)
    return actor


def render_stage(output, coords, concentration, table, scale):
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetUseDepthPeeling(1)
    renderer.SetMaximumNumberOfPeels(16)
    renderer.SetOcclusionRatio(0.0)
    renderer.AddActor(rc._surface_actor(opacity=0.075))
    renderer.AddActor(_staged_actor(coords, concentration, table))

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.AddRenderer(renderer)
    window.SetSize(1250, 950)

    direction, up = rc.CAMERA["sagittal_right"]
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    focal = np.array(camera.GetFocalPoint())
    camera.SetPosition(*(focal + camera.GetDistance() * np.array(direction)))
    camera.SetViewUp(*up)
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    camera.SetParallelScale(scale)

    light = vtk.vtkLight()
    light.SetLightTypeToCameraLight()
    light.SetPosition(-0.35, 0.45, 1.0)
    light.SetIntensity(0.85)
    renderer.AddLight(light)

    window.Render()
    grabber = vtk.vtkWindowToImageFilter()
    grabber.SetInput(window)
    grabber.SetInputBufferTypeToRGB()
    grabber.ReadFrontBufferOff()
    grabber.Update()
    writer = vtk.vtkPNGWriter()
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    writer.SetFileName(str(output))
    writer.SetInputConnection(grabber.GetOutputPort())
    writer.Write()
    window.Finalize()
    return Path(output)


def composite(reference_dir, panels, box, stage_times, output_dir):
    """Expected above obtained: the two Alzheimer strips of figure 1 of
    Weickenmeier et al. (each cut from the published page by
    extract-weickenmeier-staging.py, drawings adopted there from Jucker and
    Walker), each above the corresponding row of our computed stages."""
    strips = {}
    for case in ("tau", "amyloid"):
        path = reference_dir / f"weickenmeier_fig1_{case}.png"
        strips[case] = plt.imread(path)

    top, bottom, left, right = box
    panel_aspect = (right - left) / (bottom - top)
    # Row heights in units of the figure width: a full-width strip of aspect
    # a is 1/a tall, a row of three renders is (1/3)/panel_aspect tall.
    heights = []
    for case, _, _ in CASES:
        strip_aspect = strips[case].shape[1] / strips[case].shape[0]
        heights.extend([1.0 / strip_aspect, (1.0 / 3.0) / panel_aspect])
    width = 9.6
    figure = plt.figure(figsize=(width, 1.06 * width * sum(heights)))
    grid = figure.add_gridspec(4, 3, height_ratios=heights,
                               hspace=0.10, wspace=0.02,
                               left=0.05, right=0.995,
                               top=0.995, bottom=0.035)

    lookup = {position: path for position, path in panels}
    for pair, (case, label, seed) in enumerate(CASES):
        short = "tau" if case == "tau" else r"amyloid-$\beta$"
        strip_axis = figure.add_subplot(grid[2 * pair, :])
        strip_axis.imshow(strips[case], interpolation="lanczos")
        strip_axis.set_axis_off()
        strip_axis.text(-0.013, 0.5, f"{short}\nexpected",
                        transform=strip_axis.transAxes, fontsize=10,
                        ha="right", va="center", rotation=90,
                        linespacing=1.35)
        for column, target in enumerate(STAGES):
            axis = figure.add_subplot(grid[2 * pair + 1, column])
            rc.show_render(axis, lookup[pair, column], box)
            axis.text(0.5, -0.02,
                      f"t = {stage_times[case, target]:.1f} yr   "
                      f"(mean {target:.0%})",
                      transform=axis.transAxes, fontsize=9,
                      ha="center", va="top", color="0.25")
            if column == 0:
                axis.text(-0.04, 0.5, f"{short}\nthis work",
                          transform=axis.transAxes, fontsize=10,
                          ha="right", va="center", rotation=90,
                          linespacing=1.35)

    figure.savefig(output_dir / "seeding_patterns_expected.pdf",
                   facecolor="white")
    figure.savefig(output_dir / "seeding_patterns_expected.png", dpi=300,
                   facecolor="white")
    plt.close(figure)


def main():
    args = arguments()
    nodes = load_nodes()
    coords = np.array([node["coords"] for node in nodes])
    table = rc.lookup_table(CONCENTRATION_MAP, 0.0, 1.0)
    scale = rc.common_scale()
    sources = {"tau": args.tau, "amyloid": args.amyloid}

    plt.rcParams.update({"font.family": "serif",
                         "font.serif": ["DejaVu Serif"]})
    figure, axes = plt.subplots(2, 3, figsize=(9.6, 5.4))
    stage_times = {}
    with tempfile.TemporaryDirectory() as scratch:
        panels = []
        for row, (case, label, seed) in enumerate(CASES):
            times, values = read_profiles(sources[case])
            means = values.mean(axis=1)
            assert values.min() >= 0.0 and values.max() <= 1.0 + 1e-10, \
                "the nodal solution left [0,1]"
            for column, target in enumerate(STAGES):
                index = int(np.argmax(means >= target))
                stage_times[case, target] = times[index]
                panels.append(((row, column), render_stage(
                    Path(scratch) / f"{case}_{column}.png", coords,
                    values[index], table, scale)))
        box = rc.common_box([path for _, path in panels])
        for (row, column), path in panels:
            rc.show_render(axes[row, column], path, box)
        args.output_dir.mkdir(parents=True, exist_ok=True)
        if args.reference_dir is not None:
            composite(args.reference_dir, panels, box, stage_times,
                      args.output_dir)

    for row, (case, label, seed) in enumerate(CASES):
        axes[row, 0].text(-0.02, 0.5, label, transform=axes[row, 0].transAxes,
                          fontsize=11, ha="right", va="center", rotation=90)
        for column, target in enumerate(STAGES):
            axes[row, column].text(
                0.5, -0.02,
                f"t = {stage_times[case, target]:.1f} yr   "
                f"(mean {target:.0%})",
                transform=axes[row, column].transAxes, fontsize=9,
                ha="center", va="top", color="0.25")

    figure.subplots_adjust(left=0.045, right=0.995, top=0.995, bottom=0.06,
                           hspace=0.14, wspace=0.02)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output_dir / "seeding_patterns.pdf",
                   facecolor="white")
    figure.savefig(args.output_dir / "seeding_patterns.png", dpi=260,
                   facecolor="white")
    plt.close(figure)

    print(f"Saved seeding_patterns to {args.output_dir}")
    for case, _, seed in CASES:
        stages = ", ".join(f"{stage_times[case, t]:.1f}" for t in STAGES)
        print(f"  {case:8s} seed: {seed:18s} stages at {stages} years")


if __name__ == "__main__":
    main()
