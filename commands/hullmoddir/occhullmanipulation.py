from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeFace
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.Geom import Geom_Plane
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Ax1
from OCC.Core.GProp import GProp_GProps
from OCC.Core.TopAbs import TopAbs_EDGE, TopAbs_VERTEX
from OCC.Core.BRep import BRep_Tool, BRep_Tool_Curve, BRep_Tool_Pnt
from OCC.Core.GeomConvert import (geomconvert_SurfaceToBSplineSurface,
                                  geomconvert_CurveToBSplineCurve,
                                  GeomConvert_CompCurveToBSplineCurve)
from OCC.Core.TopoDS import TopoDS_Edge, topods_Edge, TopoDS_Wire
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopoDS import topods_Face, TopoDS_Compound
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_NurbsConvert
import numpy as np
import scipy as sci

draft = 8500
precision = 35
props_orig = GProp_GProps(gp_Pnt(0, 0, 0))


def nurbs_surface(shape):
    if type(shape) == TopoDS_Compound:
        face_explorer = TopExp_Explorer(BRepBuilderAPI_NurbsConvert(shape).Shape(), TopAbs_FACE)
        nurbs_faces = []
        while face_explorer.More():
            nurbs_faces.append(face_explorer.Current())
            face_explorer.Next()
        surfaces = [BRep_Tool.Surface(nurbs_face) for nurbs_face in nurbs_faces]
        nurbs_surfs = [geomconvert_SurfaceToBSplineSurface(surface) for surface in surfaces]
        nurbs_surf = nurbs_surfs[0]
    else:
        # TopoDS_Face converted to Nurbs
        nurbs_face = topods_Face(BRepBuilderAPI_NurbsConvert(shape).Shape())
        # GeomSurface obtained from Nurbs face
        surface = BRep_Tool.Surface(nurbs_face)
        # surface is now further converted to a bspline surface
        nurbs_surf = geomconvert_SurfaceToBSplineSurface(surface)
    return nurbs_surf


def bounding_box(shape, tol=1e-6, use_mesh=False):
    """Return the bounding box of the TopoDS_Shape `shape`
        Parameters
        ----------
        shape:
        tol: float
            tolerance of the computed boundingbox
        use_mesh : bool
            a flag that tells whether the shape has first to be meshed before the bbox
            computation. This produces more accurate results
        Original source: https://github.com/tpaviot/pythonocc-demos/blob/master/examples/core_geometry_bounding_box.py#L2
        """
    bbox = Bnd_Box()
    bbox.SetGap(tol)
    if use_mesh:
        mesh = BRepMesh_IncrementalMesh()
        mesh.SetParallelDefault(True)
        mesh.SetShape(shape)
        mesh.Perform()
        if not mesh.IsDone():
            raise AssertionError("Mesh not done.")
    brepbndlib_Add(shape, bbox, use_mesh)

    xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
    return xmin, xmax, ymin, ymax, zmin, zmax, xmax - xmin, ymax - ymin, zmax - zmin


def water_plane_particulars(shape):
    shape_sur = shape
    water_plane = Geom_Plane(gp_Pnt(0, 0, draft), gp_Dir(0, 0, 1))
    water_plane_section = BRepAlgoAPI_Section(shape_sur, water_plane, True).Shape()
    water_plane_explorer = TopExp_Explorer(water_plane_section, TopAbs_EDGE)
    wire = BRepBuilderAPI_MakeWire()
    while water_plane_explorer.More():
        wire.Add(water_plane_explorer.Current())
        if not wire.IsDone():
            edge_explorer = TopExp_Explorer(wire.Edge(), TopAbs_VERTEX)
            while edge_explorer.More():
                last_vert = edge_explorer.Current()
                edge_explorer.Next()
            vertex_explorer = TopExp_Explorer(water_plane_explorer.Current(), TopAbs_VERTEX)
            new_edge = BRepBuilderAPI_MakeEdge(last_vert, vertex_explorer.Current())
            wire.Add(new_edge.Edge())
            wire.Add(water_plane_explorer.Current())
        water_plane_explorer.Next()
    water_plane_curve = wire.Wire()
    bbox_lwl = Bnd_Box()
    bbox_lwl.SetGap(1e-8)
    brepbndlib_Add(water_plane_curve, bbox_lwl, False)
    xmin, ymin, zmin, xmax, ymax, zmax = bbox_lwl.Get()
    return xmin, xmax, ymin, ymax, zmin, zmax, xmax - xmin, ymax - ymin, zmax - zmin


def get_section_curves(shape):
    bbox = bounding_box(shape)
    dim = water_plane_particulars(shape)
    lwl = np.around(dim[6], 2)
    xmin = int(bbox[0])
    xmax = int(bbox[1])
    zmin = bbox[4]
    ymax = bbox[3]
    delta_x = int(lwl / precision)
    water_plane = draft
    y = [0, 1.1 * ymax]
    z = [1.1 * zmin, water_plane]
    subm_stations_planes = []
    for x in list(range(int(xmin), int(xmax), delta_x)):
        edge_pl = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[0], z[0]), gp_Pnt(x, y[1], z[0]))
        edge_pl2 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[1], z[0]), gp_Pnt(x, y[1], z[1]))
        edge_pl3 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[1], z[1]), gp_Pnt(x, y[0], z[1]))
        edge_pl4 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[0], z[1]), gp_Pnt(x, y[0], z[0]))
        wire_pl = BRepBuilderAPI_MakeWire()
        wire_pl.Add(edge_pl.Edge())
        wire_pl.Add(edge_pl2.Edge())
        wire_pl.Add(edge_pl3.Edge())
        wire_pl.Add(edge_pl4.Edge())
        subm_stations_planes += [BRepBuilderAPI_MakeFace(wire_pl.Wire()).Face()]
    subm_station_edges = []
    for sec in subm_stations_planes:
        subm_station_edges.append(BRepAlgoAPI_Section(shape, sec, True).Shape())
    subm_section_curves = []
    for section in subm_station_edges:
        curve = BRepAlgoAPI_Section(shape, section, True)
        explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
        while not isinstance(explorer.Current(), TopoDS_Edge):
            if explorer.More():
                xmin += 0.0001 * delta_x
                edge_pl = BRepBuilderAPI_MakeEdge(gp_Pnt(xmin, y[0], z[0]), gp_Pnt(xmin, y[1], z[0]))
                edge_pl2 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmin, y[1], z[0]), gp_Pnt(xmin, y[1], z[1]))
                edge_pl3 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmin, y[1], z[1]), gp_Pnt(xmin, y[0], z[1]))
                edge_pl4 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmin, y[0], z[1]), gp_Pnt(xmin, y[0], z[0]))
                wire_pl = BRepBuilderAPI_MakeWire()
                wire_pl.Add(edge_pl.Edge())
                wire_pl.Add(edge_pl2.Edge())
                wire_pl.Add(edge_pl3.Edge())
                wire_pl.Add(edge_pl4.Edge())
                face = BRepBuilderAPI_MakeFace(wire_pl.Wire()).Face()
                curve = BRepAlgoAPI_Section(shape, face, True)
                explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
            if not explorer.More():
                xmax -= 0.0001 * delta_x
                edge_pl = BRepBuilderAPI_MakeEdge(gp_Pnt(xmax, y[0], z[0]), gp_Pnt(xmax, y[1], z[0]))
                edge_pl2 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmax, y[1], z[0]), gp_Pnt(xmax, y[1], z[1]))
                edge_pl3 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmax, y[1], z[1]), gp_Pnt(xmax, y[0], z[1]))
                edge_pl4 = BRepBuilderAPI_MakeEdge(gp_Pnt(xmax, y[0], z[1]), gp_Pnt(xmax, y[0], z[0]))
                wire_pl = BRepBuilderAPI_MakeWire()
                wire_pl.Add(edge_pl.Edge())
                wire_pl.Add(edge_pl2.Edge())
                wire_pl.Add(edge_pl3.Edge())
                wire_pl.Add(edge_pl4.Edge())
                face = BRepBuilderAPI_MakeFace(wire_pl.Wire()).Face()
                curve = BRepAlgoAPI_Section(shape, face, True)
                explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
        wire = BRepBuilderAPI_MakeWire(explorer.Current())
        while explorer.More():
            wire.Add(explorer.Current())
            if not wire.IsDone():
                edge_explorer = TopExp_Explorer(wire.Edge(), TopAbs_VERTEX)
                while edge_explorer.More():
                    last_vert = edge_explorer.Current()
                    edge_explorer.Next()
                vertex_explorer = TopExp_Explorer(explorer.Current(), TopAbs_VERTEX)
                new_edge = BRepBuilderAPI_MakeEdge(last_vert, vertex_explorer.Current())
                wire.Add(new_edge.Edge())
                wire.Add(explorer.Current())
            explorer.Next()
        subm_section_curves += [wire.Wire()]
    return subm_section_curves


def get_sections(shape):
    bbox = bounding_box(shape)
    dim = water_plane_particulars(shape)
    xmin = np.around(dim[0], 2)
    xmax = np.around(dim[1], 2)
    zmin = bbox[4]
    ymax = bbox[3]
    water_plane = draft
    y = [0, 1.1 * ymax]
    z = [1.1 * zmin, water_plane]
    subm_section_curves = []

    # Define the number of points in the distribution
    n = precision
    # Generate a sine wave from -pi/2 to pi/2 with n points
    x_sin = np.linspace(-np.pi / 2, np.pi / 2, n)
    y_sin = np.sin(x_sin)
    # Scale and shift the values to fit the range of 0 to 10
    a, b = 0.99 * xmin, xmax * 0.95
    y_sin = (y_sin - np.min(y_sin)) / (np.max(y_sin) - np.min(y_sin)) * (b - a) + a

    for x in list(y_sin):
        edge_pl = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[0], z[0]), gp_Pnt(x, y[1], z[0]))
        edge_pl2 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[1], z[0]), gp_Pnt(x, y[1], z[1]))
        edge_pl3 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[1], z[1]), gp_Pnt(x, y[0], z[1]))
        edge_pl4 = BRepBuilderAPI_MakeEdge(gp_Pnt(x, y[0], z[1]), gp_Pnt(x, y[0], z[0]))
        wire_pl = BRepBuilderAPI_MakeWire()
        wire_pl.Add(edge_pl.Edge())
        wire_pl.Add(edge_pl2.Edge())
        wire_pl.Add(edge_pl3.Edge())
        wire_pl.Add(edge_pl4.Edge())
        subm_stations_planes = BRepBuilderAPI_MakeFace(wire_pl.Wire()).Face()
        curve = BRepAlgoAPI_Section(shape, subm_stations_planes, True)
        explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
        wire = BRepBuilderAPI_MakeWire()
        vertx = []
        while explorer.More():
            wire.Add(explorer.Current())
            if not wire.IsDone():
                edge_explorer = TopExp_Explorer(wire.Edge(), TopAbs_VERTEX)
                while edge_explorer.More():
                    last_vert = edge_explorer.Current()
                    edge_explorer.Next()
                vertex_explorer = TopExp_Explorer(explorer.Current(), TopAbs_VERTEX)
                new_edge = BRepBuilderAPI_MakeEdge(last_vert, vertex_explorer.Current())
                wire.Add(new_edge.Edge())
                wire.Add(explorer.Current())
            explorer.Next()
        close_vert_explorer = TopExp_Explorer(wire.Wire(), TopAbs_VERTEX)
        while close_vert_explorer.More():
            vertx.append(close_vert_explorer.Current())
            close_vert_explorer.Next()
        vertx.sort(key=lambda coord: BRep_Tool_Pnt(coord).Z())
        initial_vert = vertx[0]
        final_vert = vertx[-1]
        closing_edge = BRepBuilderAPI_MakeEdge(BRep_Tool_Pnt(final_vert),
                                               gp_Pnt(BRep_Tool_Pnt(final_vert).X(), 0, BRep_Tool_Pnt(final_vert).Z()))
        closing_edge2 = BRepBuilderAPI_MakeEdge(gp_Pnt(BRep_Tool_Pnt(final_vert).X(), 0, BRep_Tool_Pnt(final_vert).Z()),
                                                gp_Pnt(BRep_Tool_Pnt(final_vert).X(), 0,
                                                       BRep_Tool_Pnt(initial_vert).Z()))
        if closing_edge.IsDone() and closing_edge2.IsDone():
            wire.Add(closing_edge.Edge())
            wire.Add(closing_edge2.Edge())
        else:
            continue
        if wire.IsDone():
            subm_section_curves += [wire.Wire()]
        else:
            pass
    sub_faces = []
    for station in subm_section_curves:
        sub_faces.append(BRepBuilderAPI_MakeFace(station).Face())
    return sub_faces


def get_station_curve_gp_points(shape):
    """ Function to get points relying on station curves, returning them as gp_Pnt object"""
    curves = get_section_curves(shape)
    if type(curves[0]) == TopoDS_Wire:
        section_curves = []
        for curve in curves:
            composite_curve_builder = GeomConvert_CompCurveToBSplineCurve()
            # iterator to edges in the TopoDS_Wire
            edge_explorer = TopExp_Explorer(curve, TopAbs_EDGE)
            while edge_explorer.More():
                # getting the edge from the iterator
                edge = topods_Edge(edge_explorer.Current())
                # edge can be joined only if it is not degenerated (zero length)
                if BRep_Tool.Degenerated(edge):
                    edge_explorer.Next()
                    continue
                # the edge must be converted to Nurbs edge
                nurbs_converter = BRepBuilderAPI_NurbsConvert(edge)
                nurbs_converter.Perform(edge)
                nurbs_edge = topods_Edge(nurbs_converter.Shape())

                # here we extract the underlying curve from the Nurbs edge
                nurbs_curve = BRep_Tool_Curve(nurbs_edge)[0]

                # we convert the Nurbs curve to Bspline curve
                bspline_curve = geomconvert_CurveToBSplineCurve(nurbs_curve)

                # we can now add the Bspline curve to the composite wire curve
                composite_curve_builder.Add(bspline_curve, 1e-3)
                edge_explorer.Next()
            # GeomCurve obtained by the builder after edges are joined
            nurbs = composite_curve_builder.BSplineCurve()
            section_curves += [nurbs]
        curves = section_curves
    pnt_list = []
    pnts = []
    stations = {}
    for curve in curves:
        for points in curve.Poles():
            pnt_list.append(points)
        pnts.append(pnt_list)
        pnt_list = []
    for num in range(len(pnts)):
        stations["station " + str(num)] = pnts[num]
    proj_curv = []
    counter = 0
    pnts_on_station = {}
    p_curve = []
    for station, pnt in stations.items():
        c = curves[counter]
        for i in pnt:
            proj_curv.append(GeomAPI_ProjectPointOnCurve(i, c))
        p_curve.append(proj_curv)
        proj_curv = []
        counter += 1
    for num in range(len(p_curve)):
        pnts_on_station["station " + str(num)] = p_curve[num]
    pn_list = []
    pnts_on_station_coord = {}
    for key, value in pnts_on_station.items():
        for val in value:
            pn_list += [val.Point(1)]
        pnts_on_station_coord[key] = pn_list
        pn_list = []
    n_pnts = 27 - 2
    pnt_list = []
    for lists in pnts_on_station_coord.values():
        step = round(len(lists) / n_pnts)
        if step < 1:
            step = 1
        lists.sort(key=lambda coord: coord.Z())
        li = [lists[0]] + [lists[i] for i in list(range(1, len(lists), step))] + [lists[-1]]
        pnt_list += [li]
    pnts_on_station = {}
    for num in range(len(pnt_list)):
        pnts_on_station["station " + str(num)] = pnt_list[num]
    return pnts_on_station


def station_areas(shape):
    areas = []
    sub_faces = get_sections(shape)
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        areas.append(2 * props_orig.Mass())
    return areas


def coord_x(shape):
    coords_x = []
    sub_faces = get_sections(shape)
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        coords_x.append(props_orig.CentreOfMass().X())
    return coords_x


def volume(shape):
    areas = station_areas(shape)
    coords_x = coord_x(shape)
    vol = np.abs(np.around(sci.integrate.simps(areas, coords_x), 3))
    return vol


def prism_coeff(shape):
    bbox_lwl = water_plane_particulars(shape)
    length = np.around(bbox_lwl[6], 2)
    areas = station_areas(shape)
    max_area = max(areas)
    cp = volume(shape) / (max_area * length)
    return cp


def block_coeff(shape):
    bbox_lwl = water_plane_particulars(shape)
    lwl = np.around(bbox_lwl[6], 2)
    bwl = 2 * np.around(bbox_lwl[7], 2)
    cb = volume(shape) / (lwl * bwl * draft)
    return cb


def long_center_buoyancy(shape):
    areas = station_areas(shape)
    coords_x = coord_x(shape)
    mom_vol = [area * x for area, x in zip(areas, coords_x)]
    int_mom_vol = sci.integrate.simps(mom_vol, coords_x)
    xcb = int_mom_vol / volume(shape)
    return xcb


def forebody_norm_areas(shape):
    bbox_lwl = water_plane_particulars(shape)
    length = np.around(bbox_lwl[6], 2)
    areas = station_areas(shape)
    max_area = max(areas)
    norm_areas = [area / max_area for area in areas]  # Normalizing curve of station areas
    coords_x = coord_x(shape)
    midship = min(coords_x,
                  key=lambda x: abs(x - length / 2))  # Getting midship position/index on x coordinates list
    index_mid = coords_x.index(midship)
    forebody_norm_areas_list = norm_areas[index_mid:len(norm_areas)]
    return forebody_norm_areas_list


def parallel_body_norm(shape):
    bbox_lwl = water_plane_particulars(shape)
    length = np.around(bbox_lwl[6], 2)
    coords_x = coord_x(shape)
    midship = min(coords_x,
                  key=lambda x: abs(x - length / 2))  # Getting midship position/index on x coordinates list
    index_mid = coords_x.index(midship)
    forebody_norm_areas_list = forebody_norm_areas(shape)
    forebody_coord_x = coords_x[index_mid:len(coords_x)]  # X coordinates of forebody part
    forebody_coord_x = [midship - coord for coord in forebody_coord_x]
    norm_forebody_coord_x = [x / (length / 2) for x in forebody_coord_x]
    diff_area = []
    index_parallel = []
    for ind in range(len(forebody_norm_areas_list)):
        tol = 0.01 * max(forebody_norm_areas_list)
        if ind > 0:
            diff_area += [np.abs(np.around(forebody_norm_areas_list[ind] - forebody_norm_areas_list[ind - 1], 3))]
            if diff_area[-1] < tol:
                index_parallel += [ind - 1, ind]
    parallel_body_stations = [index_parallel[0], index_parallel[-1]]
    p = np.abs(np.around(norm_forebody_coord_x[parallel_body_stations[0]] -
                         norm_forebody_coord_x[parallel_body_stations[1]], 3))
    return p


def yb_norm(shape):
    norm_ycg = 0
    forebody_norm_areas_list = forebody_norm_areas(shape)
    for y in range(len(forebody_norm_areas_list)):
        if y > 0:
            norm_ycg += (forebody_norm_areas_list[y] + forebody_norm_areas_list[y - 1]) / 4
    norm_yb = norm_ycg / (len(forebody_norm_areas_list) - 1)
    return norm_yb


def midship_position(shape):
    coords_x = coord_x(shape)
    areas = station_areas(shape)
    # Getting midship position/index on x coordinates list
    index_mid = areas.index(max(areas))
    midship = coords_x[index_mid]
    return midship


def lackenby(shape, delta_cp=0.05, delta_p=0.2, station_starts='aft'):
    """Function to apply Lackenby general case transformation on the curve of station areas"""
    bbox_lwl = water_plane_particulars(shape)
    length = np.around(bbox_lwl[6], 2)
    cp = prism_coeff(shape)
    yb = yb_norm(shape)
    p = parallel_body_norm(shape)
    a = cp * (1 - 2 * yb) - p * (1 - cp)
    c = (1 / a) * (delta_cp - delta_p * ((1 - cp) / (1 - p)))
    d = delta_p / (c * (1.0 - p)) - p
    coords_x = coord_x(shape)
    areas = station_areas(shape)
    # Getting midship position/index on x coordinates list
    index_mid = areas.index(max(areas))
    midship = coords_x[index_mid]
    # X coordinates of forebody part
    forebody_coord_x = coords_x[index_mid:len(coords_x)]
    forebody_coord_x = [coord - midship for coord in forebody_coord_x]
    norm_forebody_coord_x = [x / (length / 2) for x in forebody_coord_x]
    # Applying the method and calculating new coordinates in x-axis
    dx = [c * (1 - x) * (x + d) for x in norm_forebody_coord_x]
    new_norm_x = [x + d_x for x, d_x in zip(norm_forebody_coord_x, dx)]
    aftbody_coord_x = coords_x[0:index_mid]
    new_fore_x = [x * length / 2 for x in new_norm_x]
    new_fore_x = [midship + coord for coord in new_fore_x]
    new_coord_x_lack = aftbody_coord_x + new_fore_x
    # if station_starts == 'fore':
    #     new_coord_x_lack.sort(reverse=False)
    return new_coord_x_lack


def set_new_coordx(shape):
    new_coords_x = lackenby(shape)
    new_stations = get_station_curve_gp_points(shape)
    for station, coords in new_stations.items():
        station_list = list(new_stations.keys())
        ind = station_list.index(station)
        new_x = new_coords_x[ind]
        for x in coords:
            x.SetX(new_x)
    return new_stations


def water_plane_area(shape):
    bwl_list = stations_bwl(shape)
    coords_x = coord_x(shape)
    return sci.integrate.simps(bwl_list, coords_x)


def stations_bwl(shape):
    bwl_list = []
    sub_faces = get_sections(shape)
    for face in sub_faces:
        bbox_face = Bnd_Box()
        bbox_face.SetGap(1e-5)
        brepbndlib_Add(face, bbox_face, False)
        bwl_list.append(2 * bbox_face.Get()[4])
    return bwl_list


def metacentric_radius(shape):
    x_axis = gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0))
    i_xx = []
    sub_faces = get_sections(shape)
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        i_xx.append(props_orig.MomentOfInertia(x_axis))
    return sum(i_xx) / volume(shape)


def vertical_center_buoyancy(shape):
    areas = station_areas(shape)
    zb_list = stations_vertical_center_buoyancy(shape)
    static_mom = [area * z for area, z in zip(areas, zb_list)]
    coords_x = coord_x(shape)
    return sci.integrate.simps(static_mom, coords_x) / volume(shape)


def stations_vertical_center_buoyancy(shape):
    station_zb = []
    sub_faces = get_sections(shape)
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        station_zb.append(props_orig.CentreOfMass().Z())
    return station_zb


def metacentre(shape):
    bbox = bounding_box(shape)
    keel = bbox[4]
    kb = vertical_center_buoyancy(shape) - keel
    bm = metacentric_radius(shape)
    km = bm + kb
    return km


############################################################################
def get_section_curves_opt(shape):
    surf = nurbs_surface(shape)
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    xmin = int(bbox[0])
    xmax = int(bbox[1])
    delta_x = int(length / precision)
    planes = []
    section_curves = []
    for point in list(range(xmin, xmax + delta_x, delta_x)):
        plane = Geom_Plane(gp_Pnt(point, 0, 0), gp_Dir(1, 0, 0))
        planes += [plane]
    for section in planes:
        curve = BRepAlgoAPI_Section(surf, section, True)
        explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
        while not isinstance(explorer.Value(), TopoDS_Edge):
            xmin += 0.1 * delta_x
            curve = BRepAlgoAPI_Section(surf, Geom_Plane(gp_Pnt(xmin, 0, 0),
                                                         gp_Dir(1, 0, 0)), True)
            explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
            if not explorer.More():
                xmax -= 0.1 * delta_x
                curve = BRepAlgoAPI_Section(surf, Geom_Plane(gp_Pnt(xmax, 0, 0),
                                                             gp_Dir(1, 0, 0)), True)
                explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
        wire = BRepBuilderAPI_MakeWire(explorer.Value())

        # joining all the wire edges in a single curve here
        # composite curve builder (can only join Bspline curves)
        composite_curve_builder = GeomConvert_CompCurveToBSplineCurve()

        # iterator to edges in the TopoDS_Wire
        edge_explorer = TopExp_Explorer(wire.Wire(), TopAbs_EDGE)
        while edge_explorer.More():
            # getting the edge from the iterator
            edge = topods_Edge(edge_explorer.Current())

            # edge can be joined only if it is not degenerated (zero length)
            if BRep_Tool.Degenerated(edge):
                edge_explorer.Next()
                continue

            # the edge must be converted to Nurbs edge
            nurbs_converter = BRepBuilderAPI_NurbsConvert(edge)
            nurbs_converter.Perform(edge)
            nurbs_edge = topods_Edge(nurbs_converter.Shape())

            # here we extract the underlying curve from the Nurbs edge
            nurbs_curve = BRep_Tool_Curve(nurbs_edge)[0]

            # we convert the Nurbs curve to Bspline curve
            bspline_curve = geomconvert_CurveToBSplineCurve(nurbs_curve)

            # we can now add the Bspline curve to the composite wire curve
            composite_curve_builder.Add(bspline_curve, 1e-5)
            edge_explorer.Next()

        # GeomCurve obtained by the builder after edges are joined
        nurbs = composite_curve_builder.BSplineCurve()
        section_curves += [nurbs]
    return section_curves


def get_station_curve_gp_points_opt(shape):
    """ Function to get points relying on station curves, returning them as gp_Pnt object"""
    curves = get_section_curves_opt(shape)
    pnt_list = []
    pnts = []
    stations = {}
    for curve in curves:
        for points in curve.Poles():
            pnt_list.append(points)
        pnts.append(pnt_list)
        pnt_list = []
    for num in range(len(pnts)):
        stations["station " + str(num)] = pnts[num]
    proj_curv = []
    counter = 0
    pnts_on_station = {}
    p_curve = []
    for station, pnt in stations.items():
        c = curves[counter]
        for i in pnt:
            proj_curv.append(GeomAPI_ProjectPointOnCurve(i, c))
        p_curve.append(proj_curv)
        proj_curv = []
        counter += 1
    for num in range(len(p_curve)):
        pnts_on_station["station " + str(num)] = p_curve[num]
    pn_list = []
    pnts_on_station_coord = {}
    for key, value in pnts_on_station.items():
        for val in value:
            pn_list += [val.Point(1)]
        pnts_on_station_coord[key] = pn_list
        pn_list = []
    return pnts_on_station_coord
