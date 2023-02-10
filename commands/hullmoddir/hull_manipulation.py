from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepBndLib import brepbndlib_Add
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeWire
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.Geom import Geom_Plane
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnCurve
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.Core.TopAbs import TopAbs_WIRE, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX
from OCC.Core.TopoDS import TopoDS_Edge
from OCC.Core.TopExp import TopExp_Explorer
from pygem.cad.cad_deformation import CADDeformation
import numpy as np
import scipy as sci
CAD = CADDeformation()


class HullManipulation:
    def __init__(self, single_surface, station_starts='aft'):
        self.shape = single_surface
        self.length = np.around(self.bounding_box()[6], 2)
        self.breath = np.around(2*self.bounding_box()[7], 2)
        self.depth = np.around(self.bounding_box()[8], 2)
        self.station_starts = station_starts

    def nurbs_surface(self):
        nurbs_surface = CAD._bspline_surface_from_face(self.shape)  # The single NURBS surface
        return nurbs_surface

    def bounding_box(self, tol=1e-6, use_mesh=False):
        """Return the bounding box of the TopoDS_Shape `shape`
            Parameters
            ----------
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
            mesh.SetShape(self.shape)
            mesh.Perform()
            if not mesh.IsDone():
                raise AssertionError("Mesh not done.")
        brepbndlib_Add(self.shape, bbox, use_mesh)

        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        return xmin, xmax, ymin, ymax, zmin, zmax, xmax - xmin, ymax - ymin, zmax - zmin

    def get_section_curves(self):
        surf = self.nurbs_surface()
        length = np.around(self.bounding_box()[6], 2)
        xmin = int(self.bounding_box()[0])
        xmax = int(self.bounding_box()[1])
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
            nurbs = CAD._bspline_curve_from_wire(wire.Wire())
            section_curves += [nurbs]
        return section_curves

    def get_station_curve_gp_points(self):
        """ Function to get points relying on station curves, returning them as gp_Pnt object"""
        curves = self.get_section_curves()
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

    def station_areas(self):
        stations = self.get_station_curve_gp_points()
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

    def volume(self):
        coord_x = self.coord_x()
        areas = self.station_areas()
        coord_x.sort()
        areas.sort()
        vol = np.abs(np.around(sci.integrate.simps(areas, coord_x), 3))
        return vol

    def prism_coeff(self):
        areas = self.station_areas()
        max_area = max(areas)
        cp = self.volume() / (max_area * self.length)
        return cp

    def block_coeff(self):
        cb = self.volume() / (self.length * self.breath * self.depth)
        return cb

    def coord_x(self):
        """Return a list with x coordinates of the stations"""
        stations = self.get_station_curve_gp_points()
        coord_x = []
        for station, coords in stations.items():
            for point in coords:
                coord_x += [np.around(point.X(), 2)]
                coord_x = list(dict.fromkeys(coord_x))
        if self.station_starts == 'fore':
            coord_x.sort(reverse=True)
        return coord_x

    def long_center_buoyancy(self):
        areas = self.station_areas()
        coord_x = self.coord_x()
        mom_vol = [area * x for area, x in zip(areas, coord_x)]
        mom_vol.sort()
        coord_x.sort()
        int_mom_vol = sci.integrate.simps(mom_vol, coord_x)
        xcb = int_mom_vol / self.volume()
        return xcb

    def forebody_norm_areas(self):
        areas = self.station_areas()
        max_area = max(areas)
        norm_areas = [area/max_area for area in areas]  # Normalizing curve of station areas
        coord_x = self.coord_x()
        midship = min(coord_x,
                      key=lambda x: abs(x - self.length / 2))  # Getting midship position/index on x coordinates list
        index_mid = coord_x.index(midship)
        forebody_norm_areas = norm_areas[index_mid:len(norm_areas)]
        return forebody_norm_areas

    def parallel_body_norm(self):
        coord_x = self.coord_x()
        midship = min(coord_x,
                      key=lambda x: abs(x - self.length/2))  # Getting midship position/index on x coordinates list
        index_mid = coord_x.index(midship)
        forebody_norm_areas = self.forebody_norm_areas()
        forebody_coord_x = coord_x[index_mid:len(coord_x)]  # X coordinates of forebody part
        forebody_coord_x = [midship - coord for coord in forebody_coord_x]
        norm_forebody_coord_x = [x / (self.length/2) for x in forebody_coord_x]
        diff_area = []
        index_parallel = []
        for ind in range(len(forebody_norm_areas)):
            tol = 0.01 * max(forebody_norm_areas)
            if ind > 0:
                diff_area += [np.abs(np.around(forebody_norm_areas[ind] - forebody_norm_areas[ind - 1], 3))]
                if diff_area[-1] < tol:
                    index_parallel += [ind - 1, ind]
        parallel_body_stations = [index_parallel[0], index_parallel[-1]]
        p = np.abs(np.around(norm_forebody_coord_x[parallel_body_stations[0]] -
                             norm_forebody_coord_x[parallel_body_stations[1]], 3))
        return p

    def yb_norm(self):
        norm_ycg = 0
        forebody_norm_areas = self.forebody_norm_areas()
        for y in range(len(forebody_norm_areas)):
            if y > 1:
                norm_ycg += (forebody_norm_areas[y] + forebody_norm_areas[y - 1]) / 4
        norm_yb = norm_ycg / (len(forebody_norm_areas) - 1)
        return norm_yb

    def lackenby(self, delta_cp=0.01, delta_p=0.05):
        """Function to apply Lackenby general case transformation on the curve of station areas"""
        cp = self.prism_coeff()
        yb = self.yb_norm()
        p = self.parallel_body_norm()
        a = cp * (1 - 2 * yb) - p * (1 - cp)
        c = (1 / a) * (delta_cp - delta_p * ((1 - cp) / (1 - p)))
        d = delta_p / (c * (1.0 - p)) - p
        coord_x = self.coord_x()
        # Getting midship position/index on x coordinates list
        midship = min(coord_x,
                      key=lambda x: abs(x - self.length/2))
        index_mid = coord_x.index(midship)
        # X coordinates of forebody part
        forebody_coord_x = coord_x[index_mid:len(coord_x)]
        forebody_coord_x = [midship - coord for coord in forebody_coord_x]
        norm_forebody_coord_x = [x / (self.length/2) for x in forebody_coord_x]
        # Applying the method and calculating new coordinates in x-axis
        dx = [c * (1 - x) * (x + d) for x in norm_forebody_coord_x]
        new_norm_x = [x + d_x for x, d_x in zip(norm_forebody_coord_x, dx)]
        aftbody_coord_x = coord_x[0:index_mid]
        new_fore_x = [x * self.length/2 for x in new_norm_x]
        new_fore_x = [midship - coord for coord in new_fore_x]
        new_coord_x_lack = aftbody_coord_x + new_fore_x
        if self.station_starts == 'fore':
            new_coord_x_lack.sort(reverse=False)
        return new_coord_x_lack

    def set_new_coordx(self):
        new_coords_x = self.lackenby()
        new_stations = self.get_station_curve_gp_points()
        for station, coords in new_stations.items():
            station_list = list(new_stations.keys())
            ind = station_list.index(station)
            new_x = new_coords_x[ind]
            for x in coords:
                x.SetX(new_x)
        return new_stations
