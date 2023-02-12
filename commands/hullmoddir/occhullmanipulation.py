from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.Geom import Geom_Plane
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.Core.gp import gp_Pnt, gp_Dir
from OCC.Core.TopAbs import TopAbs_EDGE
from OCC.Core.BRep import BRep_Tool, BRep_Tool_Curve
from OCC.Core.GeomConvert import (geomconvert_SurfaceToBSplineSurface,
                                  geomconvert_CurveToBSplineCurve,
                                  GeomConvert_CompCurveToBSplineCurve)
from OCC.Core.TopoDS import TopoDS_Edge, topods_Edge
from OCC.Core.TopExp import TopExp_Explorer

from OCC.Core.TopoDS import topods_Face
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire, BRepBuilderAPI_NurbsConvert
import numpy as np
import scipy as sci


def nurbs_surface(shape):
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


def get_section_curves(shape):
    surf = nurbs_surface(shape)
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    xmin = int(bbox[0])
    xmax = int(bbox[1])
    delta_x = int(length/20)
    planes = []
    section_curves = []
    for point in list(range(xmin, xmax+delta_x, delta_x)):
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


def get_station_curve_gp_points(shape):
    """ Function to get points relying on station curves, returning them as gp_Pnt object"""
    curves = get_section_curves(shape)
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


def station_areas(shape):
    stations = get_station_curve_gp_points(shape)
    station_coord_y = {}
    station_coord_z = {}
    areas = []
    for station, coords in stations.items():
        station_coord_y[station] = []
        station_coord_z[station] = []
        for point in coords:
            station_coord_y[station] += [point.Y()]
            station_coord_z[station] += [point.Z()]
    for y, z in zip(station_coord_y.values(), station_coord_z.values()):
        areas += [np.around(2 * sci.integrate.simps(y, z), 3)]
    return areas


def coord_x(shape, station_starts='aft'):
    """Return a list with x coordinates of the stations"""
    stations = get_station_curve_gp_points(shape)
    coords_x = []
    for station, coords in stations.items():
        for point in coords:
            coords_x += [np.around(point.X(), 2)]
            coords_x = list(dict.fromkeys(coords_x))
    if station_starts == 'fore':
        coords_x.sort(reverse=True)
    return coords_x


def volume(shape):
    coords_x = coord_x(shape)
    areas = station_areas(shape)
    coords_x.sort()
    areas.sort()
    vol = np.abs(np.around(sci.integrate.simps(areas, coords_x), 3))
    return vol


def prism_coeff(shape):
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    areas = station_areas(shape)
    max_area = max(areas)
    cp = volume(shape) / (max_area * length)
    return cp


def block_coeff(shape):
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    breath = np.around(bbox[7], 2)
    depth = np.around(bbox[8], 2)
    cb = volume(shape) / (length * breath * depth)
    return cb


def long_center_buoyancy(shape):
    areas = station_areas(shape)
    coords_x = coord_x(shape)
    mom_vol = [area * x for area, x in zip(areas, coords_x)]
    mom_vol.sort()
    coords_x.sort()
    int_mom_vol = sci.integrate.simps(mom_vol, coords_x)
    xcb = int_mom_vol / volume(shape)
    return xcb


def forebody_norm_areas(shape):
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    areas = station_areas(shape)
    max_area = max(areas)
    norm_areas = [area/max_area for area in areas]  # Normalizing curve of station areas
    coords_x = coord_x(shape)
    midship = min(coords_x,
                  key=lambda x: abs(x - length / 2))  # Getting midship position/index on x coordinates list
    index_mid = coords_x.index(midship)
    forebody_norm_areas_list = norm_areas[index_mid:len(norm_areas)]
    return forebody_norm_areas_list


def parallel_body_norm(shape):
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    coords_x = coord_x(shape)
    midship = min(coords_x,
                  key=lambda x: abs(x - length/2))  # Getting midship position/index on x coordinates list
    index_mid = coords_x.index(midship)
    forebody_norm_areas_list = forebody_norm_areas(shape)
    forebody_coord_x = coords_x[index_mid:len(coords_x)]  # X coordinates of forebody part
    forebody_coord_x = [midship - coord for coord in forebody_coord_x]
    norm_forebody_coord_x = [x / (length/2) for x in forebody_coord_x]
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
        if y > 1:
            norm_ycg += (forebody_norm_areas_list[y] + forebody_norm_areas_list[y - 1]) / 4
    norm_yb = norm_ycg / (len(forebody_norm_areas_list) - 1)
    return norm_yb


def lackenby(shape, delta_cp=0.01, delta_p=0.05, station_starts='aft'):
    """Function to apply Lackenby general case transformation on the curve of station areas"""
    bbox = bounding_box(shape)
    length = np.around(bbox[6], 2)
    cp = prism_coeff(shape)
    yb = yb_norm(shape)
    p = parallel_body_norm(shape)
    a = cp * (1 - 2 * yb) - p * (1 - cp)
    c = (1 / a) * (delta_cp - delta_p * ((1 - cp) / (1 - p)))
    d = delta_p / (c * (1.0 - p)) - p
    coords_x = coord_x(shape)
    # Getting midship position/index on x coordinates list
    midship = min(coords_x,
                  key=lambda x: abs(x - length/2))
    index_mid = coords_x.index(midship)
    # X coordinates of forebody part
    forebody_coord_x = coords_x[index_mid:len(coords_x)]
    forebody_coord_x = [midship - coord for coord in forebody_coord_x]
    norm_forebody_coord_x = [x / (length/2) for x in forebody_coord_x]
    # Applying the method and calculating new coordinates in x-axis
    dx = [c * (1 - x) * (x + d) for x in norm_forebody_coord_x]
    new_norm_x = [x + d_x for x, d_x in zip(norm_forebody_coord_x, dx)]
    aftbody_coord_x = coords_x[0:index_mid]
    new_fore_x = [x * length/2 for x in new_norm_x]
    new_fore_x = [midship - coord for coord in new_fore_x]
    new_coord_x_lack = aftbody_coord_x + new_fore_x
    if station_starts == 'fore':
        new_coord_x_lack.sort(reverse=False)
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
