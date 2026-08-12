#!/usr/bin/env python3
"""Connectome topology, after figure 6 of Fornari et al.

(a) and (b) are stacked on the left as VTK glass brains in which sphere size
AND colour both encode the degree; (c) is a large square adjacency matrix on
the right with anatomical-group strips along all four edges; two minimal
colour bars close the figure.

Every mark is a one-to-one consequence of the stored CSVs.  The drawing
conventions are listed in CONVENTIONS and printed when the script runs.
"""

import argparse
import csv
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import vtk
from matplotlib.colors import LogNorm, Normalize
from matplotlib.patches import Rectangle
from vtk.util import numpy_support

from connectome_style import REGION_COLOUR, REGION_ORDER, load_edges, load_nodes
import render_connectome as rc


# Both ramps are the blue-to-red rainbow of the reference's own figure, so
# the two degree brains and the adjacency matrix read directly against the
# published panels. The rainbow low end is a dark navy, far from the flat
# grey of an absent connection, so the weak/absent distinction survives.
DEGREE_MAP = plt.cm.jet
ADJACENCY_MAP = plt.cm.jet
ZERO_COLOUR = "#E6E8EC"

RADIUS_MIN, RADIUS_MAX = 1.25, 5.60     # mm, sphere radius in the render
LINE_WIDTH, LINE_GREY, LINE_ALPHA = 1.3, 0.18, 0.26
STRIP_GAP = 0.16                        # cells trimmed at each end of a tick
STRIP_OFFSET, STRIP_THICK = 1.1, 1.5    # cells

CONVENTIONS = """
conventions (all explicit, none of them alters a stored number):
 1. panels (a) and (b) share one colour map but are each normalised to their
    own range, printed at the two ends of the degree bar; sphere radius is
    linear in the same value, from {rmin} mm at the panel minimum to {rmax} mm
    at its maximum, so size and colour carry one and the same quantity;
 2. all {edges} connections are drawn, as flat grey lines of constant width at
    {alpha:.0%} opacity; line width and opacity encode nothing;
 3. the pial surface is the stored brain_surface.vtk, drawn at 7.5% opacity
    with the 12 smoothing passes of the house renderer, as decoration only;
 4. panel (c) is A_IJ on a base-10 logarithmic colour scale over the full
    non-zero range; the {absent} structurally absent pairs and the identically
    zero diagonal are drawn in flat grey, outside the scale;
 5. rows and columns keep the solver's own node order; nothing is permuted,
    clustered or thresholded, and the strips report the group of each index;
 6. a group tick is trimmed by {gap} cell at each end so that neighbouring
    groups separate visually; the trim carries no meaning.
"""


# ----------------------------------------------------------------- rendering

def _node_actor(coords, radii, table):
    """Spheres scaled and coloured by one and the same per-node scalar.

    The scalar stored on the points is the radius, an affine function of the
    degree, so colouring linearly over the radius range is identical to
    colouring linearly over the degree range.
    """
    points = vtk.vtkPoints()
    for point in coords:
        points.InsertNextPoint(*point)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    scalars = numpy_support.numpy_to_vtk(np.asarray(radii, float), deep=True)
    scalars.SetName("radius")
    polydata.GetPointData().SetScalars(scalars)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(1.0)
    sphere.SetThetaResolution(36)
    sphere.SetPhiResolution(36)
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToScaleByScalar()
    glyph.SetScaleFactor(1.0)
    glyph.SetColorModeToColorByScalar()
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SetScalarRange(table.GetRange())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.40)
    actor.GetProperty().SetSpecularPower(30)
    actor.GetProperty().SetDiffuse(0.78)
    actor.GetProperty().SetAmbient(0.30)
    return actor


def _line_actor(coords, edges):
    """All connections, as unshaded translucent lines of constant width.

    Shaded tubes of the house renderer read as a solid mesh at this node
    count and hide the spheres; flat translucent lines keep every connection
    visible while letting the regions read through, as in the reference.
    """
    points = vtk.vtkPoints()
    for point in coords:
        points.InsertNextPoint(*point)
    lines = vtk.vtkCellArray()
    for source, target, _ in edges:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, source)
        line.GetPointIds().SetId(1, target)
        lines.InsertNextCell(line)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(polydata)
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(LINE_GREY, LINE_GREY, LINE_GREY)
    actor.GetProperty().SetLineWidth(LINE_WIDTH)
    actor.GetProperty().SetOpacity(LINE_ALPHA)
    actor.GetProperty().SetLighting(False)
    return actor


def render_brain(output, coords, radii, table, edges, scale,
                 size=(2200, 1720)):
    """Sagittal glass brain, frontal pole to the left, as in the reference."""
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetUseDepthPeeling(1)
    renderer.SetMaximumNumberOfPeels(20)
    renderer.SetOcclusionRatio(0.0)
    renderer.AddActor(rc._surface_actor(opacity=0.075))
    renderer.AddActor(_line_actor(coords, edges))
    renderer.AddActor(_node_actor(coords, radii, table))

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.AddRenderer(renderer)
    window.SetSize(*size)

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


# ------------------------------------------------------------------- strips

def region_runs(regions):
    """Maximal runs of consecutive indices sharing an anatomical group."""
    runs, start = [], 0
    for index in range(1, len(regions) + 1):
        if index == len(regions) or regions[index] != regions[start]:
            runs.append((regions[start], start, index))
            start = index
    return runs


def draw_strips(axis, runs, n):
    """Coloured group ticks along all four edges of the matrix."""
    off, thick = STRIP_OFFSET, STRIP_THICK
    for name, first, last in runs:
        colour = REGION_COLOUR[name]
        lo, hi = first + STRIP_GAP, last - STRIP_GAP
        for x, y, w, h in ((lo, -off - thick, hi - lo, thick),
                           (lo, n + off, hi - lo, thick),
                           (-off - thick, lo, thick, hi - lo),
                           (n + off, lo, thick, hi - lo)):
            axis.add_patch(Rectangle((x, y), w, h, facecolor=colour,
                                     edgecolor="none", clip_on=False,
                                     zorder=5))


# --------------------------------------------------------------------- main

def arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path,
                        default=Path("data/connectome/fornari83/nodes.csv"))
    parser.add_argument("--edges", type=Path,
                        default=Path("data/connectome/fornari83/edges.csv"))
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = arguments()
    output = args.output_dir / "connectome_topology.png"
    nodes = load_nodes()
    edges = load_edges()
    n = len(nodes)

    adjacency = np.zeros((n, n))
    degree = np.zeros(n)
    for source, target, weight in edges:
        adjacency[source, target] += weight
        adjacency[target, source] += weight
        degree[source] += 1
        degree[target] += 1
    weighted = adjacency.sum(axis=1)

    # The per-region degrees behind panels (a) and (b), stored so the numbers
    # on the figure stay traceable to a file.
    with (args.output_dir / "connectivity_diagnostics.csv").open(
            "w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["node_id", "name", "unweighted_degree",
                         "weighted_degree"])
        for index, node in enumerate(nodes):
            writer.writerow([index, node["name"], int(degree[index]),
                             f"{weighted[index]:.6f}"])

    coords = np.array([node["coords"] for node in nodes])
    runs = region_runs([node["region"] for node in nodes])
    nonzero = adjacency[adjacency > 0]

    print(f"nodes {n}   connections {len(edges)}")
    print(f"degree              {degree.min():.6g} .. {degree.max():.6g}")
    print(f"weighted degree     {weighted.min():.6g} .. {weighted.max():.6g}")
    print(f"adjacency non-zero  {nonzero.min():.6g} .. {nonzero.max():.6g}")
    absent = n * n - n - nonzero.size
    print(f"non-zero off-diagonal entries {nonzero.size} of {n * n - n} "
          f"({nonzero.size / (n * n - n):.1%}); "
          f"absent pairs {absent}; diagonal identically 0")
    print(CONVENTIONS.format(rmin=RADIUS_MIN, rmax=RADIUS_MAX,
                             edges=len(edges), alpha=LINE_ALPHA,
                             absent=absent, gap=STRIP_GAP))

    panels = [("(a)", "non-weighted", degree),
              ("(b)", "connectivity-weighted", weighted)]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 11,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })

    # Explicit geometry: the matrix is square and as tall as the panel band,
    # and the two brains are stacked beside it at the aspect of their render.
    width, height = 13.0, 9.45
    top, floor, edge = 0.962, 0.200, 0.012
    band = (top - floor) * height                       # inches
    side = band / width                                 # matrix width, figure
    matrix_rect = [1.0 - edge - side, floor, side, top - floor]

    figure = plt.figure(figsize=(width, height))

    scale = rc.common_scale()
    with tempfile.TemporaryDirectory() as scratch:
        images = []
        for index, (_, _, values) in enumerate(panels):
            span = values.max() - values.min()
            radii = RADIUS_MIN + (RADIUS_MAX - RADIUS_MIN) * \
                (values - values.min()) / span
            table = rc.lookup_table(DEGREE_MAP, RADIUS_MIN, RADIUS_MAX)
            images.append(render_brain(Path(scratch) / f"brain_{index}.png",
                                       coords, radii, table, edges, scale))
        box = rc.common_box(images, pad=6)
        aspect = (box[3] - box[2]) / (box[1] - box[0])

        gap = 0.42 / height                             # between the brains
        brain_h = (top - floor - gap) / 2.0
        brain_w = brain_h * height / width * aspect
        left_span = (matrix_rect[0] - 0.016) - edge
        if brain_w > left_span:                         # width-limited instead
            brain_h *= left_span / brain_w
            brain_w = left_span
        brain_x = edge + 0.5 * (left_span - brain_w)
        for index, (label, caption, _) in enumerate(panels):
            bottom = top - brain_h if index == 0 else floor
            axis = figure.add_axes([brain_x, bottom, brain_w, brain_h])
            rc.show_render(axis, images[index], box)
            axis.text(0.02, 0.99, label, transform=axis.transAxes,
                      fontsize=14, style="italic", va="top", ha="left")
            axis.text(0.5, -0.035, caption, transform=axis.transAxes,
                      fontsize=12.5, va="top", ha="center")
    brain_centre = brain_x + 0.5 * brain_w

    # ------------------------------------------------------- adjacency panel
    axis = figure.add_axes(matrix_rect)
    masked = np.ma.masked_where(adjacency <= 0.0, adjacency)
    cmap = ADJACENCY_MAP.copy()
    cmap.set_bad(ZERO_COLOUR)
    norm = LogNorm(nonzero.min(), nonzero.max())
    # Origin at the lower left, as printed: the intra-hemisphere blocks sit
    # in the lower-left and upper-right quadrants, as the reference describes.
    # The extent must agree with the origin, or it silently re-inverts it.
    axis.imshow(masked, cmap=cmap, norm=norm, origin="lower",
                interpolation="nearest", extent=(0, n, 0, n))
    axis.set_aspect("equal")
    draw_strips(axis, runs, n)

    # Thin dashed boxes bound the eight intra-lobe clusters along the
    # diagonal, four cortical lobes per hemisphere, computed from the region
    # assignment and never typed in. Each box is the longest consecutive run
    # of its lobe in the solver order, so it contains cells of that lobe and
    # of nothing else; the one parietal and one temporal vertex per
    # hemisphere that the atlas enumeration separates from their lobe, by
    # interleaving the limbic belt, stay outside and are located by the
    # strips. The benchmark record names them. The boxes mark structure,
    # not single connections.
    for side in ("right", "left"):
        for lobe in ("frontal", "temporal", "parietal", "occipital"):
            members = sorted(k for k, node in enumerate(nodes)
                             if node["hemisphere"] == side
                             and node["region"] == lobe)
            runs_of_lobe, start = [], members[0]
            for previous, current in zip(members, members[1:]):
                if current != previous + 1:
                    runs_of_lobe.append((start, previous))
                    start = current
            runs_of_lobe.append((start, members[-1]))
            low, high = max(runs_of_lobe, key=lambda run: run[1] - run[0])
            assert all(nodes[k]["region"] == lobe
                       for k in range(low, high + 1)), lobe
            axis.add_patch(Rectangle((low, low), high + 1 - low,
                                     high + 1 - low, facecolor="none",
                                     edgecolor="0.15", linewidth=0.9,
                                     linestyle=(0, (4, 3)), zorder=6))
    margin = 3.6
    axis.set_xlim(-margin, n + margin)
    axis.set_ylim(-margin, n + margin)
    axis.set_axis_off()
    axis.text(-0.005, 1.004, "(c)", transform=axis.transAxes, fontsize=14,
              style="italic", va="bottom", ha="left")
    axis.text(0.5, -0.016,
              "rows from the bottom, columns from the left, in solver order:"
              "   right hemisphere 0-40,"
              "   left 41-81,   brainstem 82", transform=axis.transAxes,
              ha="center", va="top", fontsize=10, color="0.42")

    # ---------------------------------------------------------- colour bars
    matrix_centre = matrix_rect[0] + 0.5 * matrix_rect[2]
    bar_y, bar_h = 0.100, 0.020
    half = 0.128

    cax = figure.add_axes([brain_centre - half, bar_y, 2 * half, bar_h])
    bar = figure.colorbar(
        plt.cm.ScalarMappable(cmap=DEGREE_MAP, norm=Normalize(0, 1)),
        cax=cax, orientation="horizontal")
    bar.outline.set_visible(False)
    cax.set_xticks([])
    cax.set_yticks([])
    figure.text(brain_centre, bar_y + bar_h + 0.012, r"degree  $D_{II}$",
                ha="center", va="bottom", fontsize=12.5)
    for row, (label, _, values) in enumerate(panels):
        y = bar_y - 0.027 - 0.030 * row
        figure.text(brain_centre - half - 0.008, y,
                    f"{label}  {values.min():.5g}",
                    ha="right", va="center", fontsize=10.5, color="0.25")
        figure.text(brain_centre + half + 0.008, y, f"{values.max():.5g}",
                    ha="left", va="center", fontsize=10.5, color="0.25")

    cax = figure.add_axes([matrix_centre - half, bar_y, 2 * half, bar_h])
    bar = figure.colorbar(plt.cm.ScalarMappable(cmap=ADJACENCY_MAP, norm=norm),
                          cax=cax, orientation="horizontal")
    bar.outline.set_visible(False)
    cax.set_xticks([])
    cax.set_yticks([])
    cax.minorticks_off()
    figure.text(matrix_centre, bar_y + bar_h + 0.012,
                r"adjacency  $A_{IJ}$,  logarithmic", ha="center",
                va="bottom", fontsize=12.5)
    numbers_y = bar_y - 0.027
    figure.text(matrix_centre - half - 0.008, numbers_y,
                f"{nonzero.min():.5g}", ha="right", va="center",
                fontsize=10.5, color="0.25")
    figure.text(matrix_centre + half + 0.008, numbers_y,
                f"{nonzero.max():.5g}", ha="left", va="center",
                fontsize=10.5, color="0.25")
    figure.patches.append(Rectangle((matrix_centre - 0.062, numbers_y - 0.010),
                                    0.018, 0.020,
                                    transform=figure.transFigure,
                                    facecolor=ZERO_COLOUR, edgecolor="none"))
    figure.text(matrix_centre - 0.038, numbers_y, "no connection",
                ha="left", va="center", fontsize=10.5, color="0.25")

    handles = [Rectangle((0, 0), 1, 1, facecolor=REGION_COLOUR[name])
               for name in REGION_ORDER]
    figure.legend(handles, list(REGION_ORDER), loc="center",
                  bbox_to_anchor=(matrix_centre, 0.024), ncol=7, frameon=False,
                  fontsize=10.5, handlelength=1.1, handleheight=0.9,
                  handletextpad=0.4, columnspacing=1.2)

    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=200, facecolor="white", bbox_inches="tight")
    figure.savefig(output.with_suffix(".pdf"), facecolor="white",
                   bbox_inches="tight")
    plt.close(figure)
    print("wrote", output)


if __name__ == "__main__":
    main()
