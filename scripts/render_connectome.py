"""Three-dimensional rendering of the connectome, through VTK.

The reference works do not state which tool produced their network figures, but
their visual signature -- a translucent pial surface, shaded spheres at the
regions and shaded tubes along the connections -- is that of a rendering engine
with lighting and depth, not of a two-dimensional plot. VTK is the engine
ParaView is built on, so rendering here and displaying the same VTK output
interactively in ParaView produce the same picture.

Nothing about the data is altered for the rendering: the surface is the file
already written for ParaView, and the vertices are at the anatomical
coordinates of `nodes.csv`, in the same frame.
"""

from pathlib import Path

import numpy as np
import vtk
from vtk.util import numpy_support

SURFACE = Path("data/connectome/anatomy/brain_surface.vtk")

# Camera position and up direction for each anatomical projection, in the
# right/anterior/superior frame of the data. The last three are the
# projections of the brain-network figure of Fornari et al.: the sagittal
# view taken from the right so the frontal pole points left, the top-down
# view with the anterior direction to the left, and the unlabelled oblique.
CAMERA = {
    "sagittal": ((-1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    "coronal": ((0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    "axial": ((0.0, 0.0, 1.0), (0.0, 1.0, 0.0)),
    "sagittal_right": ((1.0, 0.0, 0.0), (0.0, 0.0, 1.0)),
    "longitudinal": ((0.0, 0.0, 1.0), (1.0, 0.0, 0.0)),
    "oblique": ((0.32, 0.60, 0.73), (0.0, 0.0, 1.0)),
}


def lookup_table(colourmap, low, high, samples=256):
    """Convert a matplotlib colormap into a VTK lookup table."""
    table = vtk.vtkLookupTable()
    table.SetNumberOfTableValues(samples)
    table.SetRange(low, high)
    for index in range(samples):
        red, green, blue, _ = colourmap(index / (samples - 1))
        table.SetTableValue(index, red, green, blue, 1.0)
    table.Build()
    return table


def _surface_actor(opacity=0.10, colour=(0.60, 0.65, 0.72)):
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(SURFACE))
    reader.Update()
    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputConnection(reader.GetOutputPort())
    smoother.SetNumberOfIterations(12)
    smoother.NonManifoldSmoothingOn()
    smoother.NormalizeCoordinatesOn()
    smoother.Update()
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(smoother.GetOutputPort())
    normals.SetFeatureAngle(60.0)
    normals.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    mapper.ScalarVisibilityOff()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(*colour)
    actor.GetProperty().SetOpacity(opacity)
    actor.GetProperty().SetSpecular(0.15)
    actor.GetProperty().SetDiffuse(0.85)
    actor.GetProperty().SetAmbient(0.25)
    return actor


def _node_actor(coords, values, table, radius):
    points = vtk.vtkPoints()
    for point in coords:
        points.InsertNextPoint(*point)
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    scalars = numpy_support.numpy_to_vtk(np.asarray(values, dtype=float),
                                         deep=True)
    scalars.SetName("value")
    polydata.GetPointData().SetScalars(scalars)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(28)
    sphere.SetPhiResolution(28)
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToDataScalingOff()
    glyph.SetColorModeToColorByScalar()
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SetScalarRange(table.GetRange())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.45)
    actor.GetProperty().SetSpecularPower(28)
    actor.GetProperty().SetDiffuse(0.75)
    actor.GetProperty().SetAmbient(0.30)
    return actor


def _edge_actor(coords, edges, values, table, radii):
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

    # The tube filter takes a single radius, so the connections are grouped by
    # rounded radius and one tube filter is run per group.
    tubes = []
    unique = {}
    for index, radius in enumerate(radii):
        unique.setdefault(round(radius, 3), []).append(index)

    for radius, members in unique.items():
        subset = vtk.vtkPolyData()
        subset_points = vtk.vtkPoints()
        subset_points.DeepCopy(points)
        subset.SetPoints(subset_points)
        subset_lines = vtk.vtkCellArray()
        subset_values = vtk.vtkFloatArray()
        subset_values.SetName("value")
        for index in members:
            source, target, _ = edges[index]
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, source)
            line.GetPointIds().SetId(1, target)
            subset_lines.InsertNextCell(line)
            subset_values.InsertNextValue(values[index])
        subset.SetLines(subset_lines)
        subset.GetCellData().SetScalars(subset_values)

        tube = vtk.vtkTubeFilter()
        tube.SetInputData(subset)
        tube.SetRadius(max(radius, 1e-3))
        tube.SetNumberOfSides(10)
        tube.CappingOn()
        tube.Update()
        tubes.append(tube.GetOutput())

    append = vtk.vtkAppendPolyData()
    for tube in tubes:
        append.AddInputData(tube)
    append.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(append.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SetScalarRange(table.GetRange())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetSpecular(0.25)
    actor.GetProperty().SetDiffuse(0.85)
    actor.GetProperty().SetAmbient(0.30)
    return actor


def render(output, view, coords, node_values, node_table, node_radius=3.4,
           edges=None, edge_values=None, edge_table=None, edge_radii=None,
           size=(1100, 900), scale=None, surface_opacity=0.10):
    """Render one anatomical view to a PNG and return its path."""
    direction, up = CAMERA[view]
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetUseDepthPeeling(1)
    renderer.SetMaximumNumberOfPeels(12)
    renderer.SetOcclusionRatio(0.0)

    renderer.AddActor(_surface_actor(opacity=surface_opacity))
    if edges:
        renderer.AddActor(_edge_actor(coords, edges, edge_values, edge_table,
                                      edge_radii))
    renderer.AddActor(_node_actor(coords, node_values, node_table,
                                  node_radius))

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.AddRenderer(renderer)
    window.SetSize(*size)

    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    focal = np.array(camera.GetFocalPoint())
    distance = camera.GetDistance()
    camera.SetPosition(*(focal + distance * np.array(direction)))
    camera.SetViewUp(*up)
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    if scale is not None:
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
    return Path(output), camera.GetParallelScale()


def common_scale(margin=1.04):
    """Parallel scale that frames the whole surface in any projection."""
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(SURFACE))
    reader.Update()
    bounds = reader.GetOutput().GetBounds()
    extents = [bounds[1] - bounds[0], bounds[3] - bounds[2],
               bounds[5] - bounds[4]]
    return margin * 0.5 * max(extents)


def _read(path):
    from PIL import Image

    return np.asarray(Image.open(path).convert("RGB"))


def ink_box(path, threshold=250):
    """Bounding box of the drawn content of a rendered PNG."""
    image = _read(path)
    ink = np.any(image < threshold, axis=2)
    if not ink.any():
        return 0, image.shape[0], 0, image.shape[1]
    rows = np.where(ink.any(axis=1))[0]
    columns = np.where(ink.any(axis=0))[0]
    return rows[0], rows[-1] + 1, columns[0], columns[-1] + 1


def common_box(paths, pad=8):
    """Union of the content boxes, so panels of one figure stay aligned.

    Every view is rendered with the same camera scale into the same window, so
    cropping them all to one box keeps their scale and their alignment.
    """
    boxes = [ink_box(path) for path in paths]
    top = max(min(box[0] for box in boxes) - pad, 0)
    bottom = min(box[1] for box in boxes)
    bottom = max(max(box[1] for box in boxes) + pad, bottom)
    left = max(min(box[2] for box in boxes) - pad, 0)
    right = max(box[3] for box in boxes) + pad
    return top, bottom, left, right


def show_render(axis, path, box=None):
    """Place a rendered view on a matplotlib axis, without decorations."""
    image = _read(path)
    if box is None:
        box = ink_box(path)
    top, bottom, left, right = box
    axis.imshow(image[top:min(bottom, image.shape[0]),
                      left:min(right, image.shape[1])],
                interpolation="lanczos")
    axis.set_xticks([])
    axis.set_yticks([])
    axis.set_axis_off()


def render_field(output, view, solution, scalar_range, table, radius=1.15,
                 size=(1200, 950), scale=None, surface_opacity=0.09,
                 scalar="c"):
    """Render the finite element field on the vertices of the extended graph.

    `solution` is one of the `solution_*.vtp` files written by the solver, so
    what is drawn is exactly the degrees of freedom the simulation carried:
    the region vertices and the nodes interpolated along every connection. No
    resampling and no interpolation into the brain volume take place.
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(solution))
    reader.Update()
    polydata = reader.GetOutput()
    polydata.GetPointData().SetActiveScalars(scalar)

    sphere = vtk.vtkSphereSource()
    sphere.SetRadius(radius)
    sphere.SetThetaResolution(10)
    sphere.SetPhiResolution(10)
    glyph = vtk.vtkGlyph3D()
    glyph.SetInputData(polydata)
    glyph.SetSourceConnection(sphere.GetOutputPort())
    glyph.SetScaleModeToDataScalingOff()
    glyph.SetColorModeToColorByScalar()
    glyph.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(glyph.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SetScalarRange(*scalar_range)
    field = vtk.vtkActor()
    field.SetMapper(mapper)
    field.GetProperty().SetSpecular(0.30)
    field.GetProperty().SetSpecularPower(24)
    field.GetProperty().SetDiffuse(0.80)
    field.GetProperty().SetAmbient(0.32)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    renderer.SetUseDepthPeeling(1)
    renderer.SetMaximumNumberOfPeels(12)
    renderer.SetOcclusionRatio(0.0)
    renderer.AddActor(_surface_actor(opacity=surface_opacity))
    renderer.AddActor(field)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetAlphaBitPlanes(1)
    window.SetMultiSamples(0)
    window.AddRenderer(renderer)
    window.SetSize(*size)

    direction, up = CAMERA[view]
    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    focal = np.array(camera.GetFocalPoint())
    distance = camera.GetDistance()
    camera.SetPosition(*(focal + distance * np.array(direction)))
    camera.SetViewUp(*up)
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    if scale is not None:
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


# Ramp anchored so that its zero is the neutral grey of the discretisation:
# healthy graph reads as a wireframe and the lesion emerges from it, which is
# the device that makes the reference figures legible. Approximately linear in
# CIE lightness, so equal steps in c are equal steps in perceived darkness.
FIELD_RAMP = (
    (0.00, (0.815, 0.826, 0.840)),
    (0.25, (0.596, 0.686, 0.812)),
    (0.50, (0.361, 0.533, 0.729)),
    (0.75, (0.161, 0.353, 0.573)),
    (1.00, (0.043, 0.145, 0.318)),
)


def field_colourmap():
    """The single-hue pale-to-saturated ramp used by the section figures."""
    from matplotlib.colors import LinearSegmentedColormap

    return LinearSegmentedColormap.from_list(
        "connectome_field", [(stop, colour) for stop, colour in FIELD_RAMP])


def _clipped_surface(normal=(1.0, 0.0, 0.0), origin=(0.0, 0.0, 0.0),
                     base=(0.855, 0.855, 0.865), cap=(0.790, 0.790, 0.805)):
    """Opaque half of the pial surface, capped at the cut plane.

    The surface is watertight, so the cut produces a genuine section rather
    than a hollow shell, and the cap is drawn slightly darker so that it reads
    as a cut face. The surface is not smoothed: every smoothing pass moves the
    rendered anatomy and buys a caveat this figure does not need.
    """
    reader = vtk.vtkPolyDataReader()
    reader.SetFileName(str(SURFACE))
    reader.Update()

    plane = vtk.vtkPlane()
    plane.SetOrigin(*origin)
    plane.SetNormal(*normal)
    planes = vtk.vtkPlaneCollection()
    planes.AddItem(plane)

    clipper = vtk.vtkClipClosedSurface()
    clipper.SetInputConnection(reader.GetOutputPort())
    clipper.SetClippingPlanes(planes)
    clipper.GenerateFacesOn()
    clipper.SetScalarModeToColors()
    clipper.SetBaseColor(*base)
    clipper.SetClipColor(*cap)
    clipper.Update()

    normals = vtk.vtkPolyDataNormals()
    normals.SetInputConnection(clipper.GetOutputPort())
    normals.SetFeatureAngle(60.0)
    normals.Update()

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(normals.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetAmbient(0.28)
    actor.GetProperty().SetDiffuse(0.72)
    actor.GetProperty().SetSpecular(0.04)
    return actor


def _clip_polydata(polydata, normal=(1.0, 0.0, 0.0), origin=(0.0, 0.0, 0.0)):
    plane = vtk.vtkPlane()
    plane.SetOrigin(*origin)
    plane.SetNormal(*normal)
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(polydata)
    clipper.SetClipFunction(plane)
    clipper.InsideOutOn()
    clipper.Update()
    return clipper.GetOutput()


def level_set_points(polydata, level=0.5, scalar="c"):
    """Exact crossings of `level` along the P1 field, segment by segment.

    For a segment carrying endpoint values a and b with (a-level)(b-level) < 0,
    the P1 interpolant crosses `level` once, at the point returned here. No
    element is truncated and nothing is hidden: this only marks where the level
    set falls.
    """
    values = polydata.GetPointData().GetArray(scalar)
    crossings = []
    lines = polydata.GetLines()
    lines.InitTraversal()
    ids = vtk.vtkIdList()
    while lines.GetNextCell(ids):
        for k in range(ids.GetNumberOfIds() - 1):
            first, second = ids.GetId(k), ids.GetId(k + 1)
            a, b = values.GetValue(first), values.GetValue(second)
            if (a - level) * (b - level) < 0.0:
                weight = (level - a) / (b - a)
                start = np.array(polydata.GetPoint(first))
                end = np.array(polydata.GetPoint(second))
                crossings.append(start + weight * (end - start))
    return np.array(crossings) if crossings else np.zeros((0, 3))


def render_section(output, solution, table, node_coords=None,
                   node_values=None, level=0.5, tube_radius=0.55,
                   node_radius=2.6, scalar_range=(0.0, 1.0), size=(1250, 1000),
                   scale=None, normal=(1.0, 0.0, 0.0), scalar="c"):
    """Sagittal section: the near half of the graph against the far half-brain.

    The graph is clipped to one side of the plane and the pial surface to the
    other, so a single depth buffer gives correct occlusion and the network
    sits inside the anatomy instead of floating over it. Tube and sphere radii
    are constant: thickness encodes nothing, only colour carries a value.
    """
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(solution))
    reader.Update()
    graph = _clip_polydata(reader.GetOutput(), normal=normal)
    graph.GetPointData().SetActiveScalars(scalar)

    renderer = vtk.vtkRenderer()
    renderer.SetBackground(1.0, 1.0, 1.0)
    # The graph is clipped to the near half-space and the surface to the far
    # one, so one depth buffer occludes correctly without any layer trick.
    renderer.AddActor(_clipped_surface(normal=normal))

    tubes = vtk.vtkTubeFilter()
    tubes.SetInputData(graph)
    tubes.SetRadius(tube_radius)
    tubes.SetNumberOfSides(8)
    tubes.CappingOn()
    tubes.Update()
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tubes.GetOutputPort())
    mapper.SetLookupTable(table)
    mapper.SetScalarRange(*scalar_range)
    network = vtk.vtkActor()
    network.SetMapper(mapper)
    network.GetProperty().SetAmbient(0.55)
    network.GetProperty().SetDiffuse(0.45)
    renderer.AddActor(network)

    if node_coords is not None and len(node_coords):
        renderer.AddActor(_node_actor(node_coords, node_values, table,
                                      node_radius))

    marks = level_set_points(graph, level=level, scalar=scalar)
    if len(marks):
        points = vtk.vtkPoints()
        for point in marks:
            points.InsertNextPoint(*point)
        cloud = vtk.vtkPolyData()
        cloud.SetPoints(points)
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1.0)
        sphere.SetThetaResolution(12)
        sphere.SetPhiResolution(12)
        glyph = vtk.vtkGlyph3D()
        glyph.SetInputData(cloud)
        glyph.SetSourceConnection(sphere.GetOutputPort())
        glyph.Update()
        mark_mapper = vtk.vtkPolyDataMapper()
        mark_mapper.SetInputConnection(glyph.GetOutputPort())
        mark_mapper.ScalarVisibilityOff()
        mark_actor = vtk.vtkActor()
        mark_actor.SetMapper(mark_mapper)
        mark_actor.GetProperty().SetColor(1.0, 1.0, 1.0)
        mark_actor.GetProperty().SetAmbient(1.0)
        mark_actor.GetProperty().SetDiffuse(0.0)
        renderer.AddActor(mark_actor)

    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetMultiSamples(8)
    window.AddRenderer(renderer)
    window.SetSize(*size)

    renderer.ResetCamera()
    camera = renderer.GetActiveCamera()
    focal = np.array(camera.GetFocalPoint())
    distance = camera.GetDistance()
    camera.SetPosition(*(focal - distance * np.array(normal)))
    camera.SetViewUp(0.0, 0.0, 1.0)
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    if scale is not None:
        camera.SetParallelScale(scale)

    light = vtk.vtkLight()
    light.SetLightTypeToCameraLight()
    light.SetPosition(-0.3, 0.4, 1.0)
    light.SetIntensity(0.8)
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
    return Path(output), len(marks)
