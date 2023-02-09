from hullmoddir.cad_deformation import *
from hullformdir.hullform import *
from hullformdir.shipstability import *
from typing import Dict
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.TopAbs import TopAbs_WIRE, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX
from itertools import product
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec
from OCC.Core.Geom import Geom_Plane
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import TopoDS_Edge
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section, BRepAlgoAPI_Fuse
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRep import BRep_Tool_Pnt
import csv
from hullmoddir.occhelperfun import get_open_mesh_from_TopoDS_using_shape_tesselator
from OCC.Display.SimpleGui import init_display
import openmesh as om
F = False
T = True


def get_relative_path(full_path,referent_file):
    relative_path = ''
    return relative_path


def get_full_path(file_name,referent_file):
    full_path=''
    return full_path


def read_hoc_file(file_path: str):
    shipdata = {}
    with open(file_path, newline='') as csvfile:
        f = csv.DictReader(csvfile)
        shipset = 0
        for row in f:  # there is only one row after header row!!!!!
            shipset = row
        list_str_keys = {"surface_file_relative_path"}
        list_floats_keys = {}
        for key, value in shipset.items():
            if key in list_floats_keys:
                shipdata[key] = float(shipset[key])
            elif key in list_floats_keys:
                shipdata[key] = str(shipset[key])
    return shipdata


def read_igs_file(file_path):
    curves = []
    surfaces = []
    igsreader = IGESControl_Reader()
    igsreader.ReadFile(file_path)
    igsreader.TransferRoots()
    nbr = igsreader.NbShapes()
    # load shape from IGES
    for i in range(1, nbr + 1):
        shp = igsreader.Shape(i)
        if not shp.IsNull():
            if shp.ShapeType() == TopAbs_WIRE:
                curves.append(shp)
            elif shp.ShapeType() == TopAbs_FACE:
                surfaces.append(shp)
    return curves, surfaces


class OCCHullform(HullForm, CADDeformation, ShipStability):
    def __init__(self, fileName, name="", translation=np.zeros(3)):
        super().__init__(fileName, name, translation)
        self.tolerance = 1e-6
        self.surface = None
        self._shipdata: Dict[str, object] = {}
        self._surfaces = []
        self._curves = []
        occ_entities = None
        filename_no_ext, file_extension = os.path.splitext(fileName)
        if file_extension == ".hoc":
            self._shipdata = read_hoc_file()
            surface_path = self._shipdata.get('surface_file_relative_path')
            if surface_path is not None:
                full_path = get_full_path(surface_path,self.filename)
                occ_entities = read_hoc_file(full_path)
        elif file_extension == ".igs":
            occ_entities = read_igs_file(self.filename)
        if occ_entities is not None:
            self._surfaces = occ_entities[1]
            self._curves = occ_entities[0]
        self.regenerateHullHorm()

    def regenerateHullHorm(self):
        if len(self._surfaces) > 0:
            self.mesh = get_open_mesh_from_TopoDS_using_shape_tesselator(self._surfaces[0], mesh_quality=0.1)
            self.mirror_mesh_around_symmetry_plan()

    def modify_form(self):
        # do some change with poles
        bsp_surf = []
        new_surfs = []
        for face in self._surfaces:
            bsp_surf += [self.bspline_surface_from_face(face)]
        for surf in bsp_surf:
            n_poles_u = surf.NbUPoles()
            n_poles_v = surf.NbVPoles()
            # Array which will contain the coordinates of the poles
            poles_coordinates = np.zeros(shape=(n_poles_u * n_poles_v, 3))
            pole_ids = product(range(n_poles_u), range(n_poles_v))
            # Cycle over the poles to get their coordinates
            for pole_id, (u, v) in enumerate(pole_ids):
                pole = surf.Pole(u + 1, v + 1)
                poles_coordinates[pole_id, :] = self.pole_get_components(pole)
            new_coord = np.zeros(shape=(n_poles_u * n_poles_v, 3))
            for pnts in range(len(poles_coordinates)):
                new_coord[pnts, :] = [poles_coordinates[pnts][0] * 1.2,
                                      poles_coordinates[pnts][1],
                                      poles_coordinates[pnts][2]]
            pole_ids = product(range(n_poles_u), range(n_poles_v))
            for pole_id_, (u_, v_) in enumerate(pole_ids):
                new_pole = self.pole_set_components(new_coord[pole_id_, :])
                surf.SetPole(u_ + 1, v_ + 1, new_pole)
            new_surfs += [BRepBuilderAPI_MakeFace(surf, 1e-3).Face()]
        self._surfaces = new_surfs
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()
        print('Hull form changed')

    def cross_section_curves(self):
        if self.mesh is not None:
            bb = self.bbox
            xmin = int(1.025*bb.minCoord[0])
            xmax = int(0.975*bb.maxCoord[0])
            length = abs(xmin-xmax)
            delta_x = int(length/35)
            planes = []
            section_curves = []
            for point in list(range(xmin, xmax, delta_x)):
                plane = Geom_Plane(gp_Pnt(point, 0, 0), gp_Dir(1, 0, 0))
                planes += [plane]
            for section in planes:
                curve = BRepAlgoAPI_Section(self._surfaces[0], section, True)
                explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
                while not isinstance(explorer.Value(), TopoDS_Edge):
                    xmin += 0.01*delta_x
                    curve = BRepAlgoAPI_Section(self._surfaces[0], Geom_Plane(gp_Pnt(xmin, 0, 0),
                                                                              gp_Dir(1, 0, 0)), True)
                    explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
                wire = BRepBuilderAPI_MakeWire(explorer.Value())
                nurbs = self.bspline_curve_from_wire(wire.Wire())
                section_curves += [nurbs]
            if T:
                display, start_display, add_menu, add_function_to_menu = init_display()
                for curv in section_curves:
                    display.DisplayShape(curv, update=True, color='Blue1')
                start_display()

    def extrude_side_shell(self):
        cover_high = 10000
        edges = TopExp_Explorer(self._surfaces[0], TopAbs_EDGE)
        z_ed = []
        while edges.More():
            vertices = TopExp_Explorer(edges.Current(), TopAbs_VERTEX)
            zmin = BRep_Tool_Pnt(vertices.Current()).Z()
            vertices.Next()
            zmax = BRep_Tool_Pnt(vertices.Current()).Z()
            z_ed.append([zmin, zmax])
            if z_ed[-1][0] <= zmin and z_ed[-1][1] <= zmax:
                border = edges.Current()
            edges.Next()
        pris = BRepPrimAPI_MakePrism(border, gp_Vec(0, 0, cover_high)).Shape()
        self._surfaces[0] = BRepAlgoAPI_Fuse(self._surfaces[0], pris).Shape()
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()
        print('Hull shell changed')

    def close_transom(self):
        bb = self.bbox
        xmin = bb.minCoord[0]
        y = 0
        plane_point = np.array([xmin, 0, 0])
        plane_normal = np.array([0, 0, -1])
        fvs = self.mesh.fv_indices()
        fvs = fvs.tolist()
        points = self.mesh.points().tolist()
        new_fvs, new_pts = self.get_mesh_below_inclined_waterline(fvs, points, plane_point, plane_normal)
        self.mesh = om.TriMesh(new_pts, new_fvs)

    def close_cover(self):
        hwaterline = 13000
        plane_point = np.array([0, 0, hwaterline])
        plane_normal = np.array([0, 0, -1])
        fvs = self.mesh.fv_indices()
        fvs = fvs.tolist()
        points = self.mesh.points().tolist()
        new_fvs, new_pts = self.get_mesh_below_inclined_waterline(fvs, points, plane_point, plane_normal)
        self.mesh = om.TriMesh(new_pts, new_fvs)

