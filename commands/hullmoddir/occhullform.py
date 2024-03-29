# from hullmoddir.cad_deformation import *
import os
from copy import copy
from hullformdir.hullform import *
from hullformdir.shipstability import *
from hullmoddir.pygemhullform import *
from typing import Dict
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.Interface import Interface_Static
from OCC.Core.TopAbs import TopAbs_WIRE, TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX
from itertools import product
from OCC.Core.BRepBuilderAPI import (BRepBuilderAPI_MakeFace, BRepBuilderAPI_MakeWire, BRepBuilderAPI_GTransform,
                                     BRepBuilderAPI_Transform, BRepBuilderAPI_NurbsConvert, BRepBuilderAPI_MakeEdge)
from OCC.Core.gp import gp_Pnt, gp_Dir, gp_Vec, gp_Ax2, gp_Trsf, gp_GTrsf
from OCC.Core.Geom import Geom_Plane
from OCC.Core.GeomConvert import GeomConvert_CompCurveToBSplineCurve, geomconvert_CurveToBSplineCurve
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import TopoDS_Edge, TopoDS_Compound, topods_Edge
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section, BRepAlgoAPI_Fuse
from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakePrism
from OCC.Core.BRep import BRep_Tool, BRep_Tool_Pnt, BRep_Tool_Curve
import csv
from hullmoddir.occhelperfun import get_open_mesh_from_TopoDS_using_shape_tesselator
from OCC.Display.SimpleGui import init_display
import hullmoddir.occhullmanipulation as occhm

F = False
T = True


def get_relative_path(full_path, referent_file):
    relative_path = ''
    return relative_path


def get_full_path(file_name, referent_file):
    full_path = ''
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
    Interface_Static.SetCVal("xstep.cascade.unit", "M")
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


class OCCHullform(HullForm, PyGemHullform):
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
                full_path = get_full_path(surface_path, self.filename)
                occ_entities = read_hoc_file(full_path)
        elif file_extension == ".igs":
            occ_entities = read_igs_file(self.filename)
        if occ_entities is not None:
            self._surfaces = occ_entities[1]
            self._curves = occ_entities[0]
            self._original_surfaces = copy(self._surfaces)
        self.regenerateHullHorm()
        self.calc_ship_dims()
        # self.position_form(x_offset = 0.1)
        self.position_form(x_offset = 0)

        # self.calc_ship_dims()

        # self.make_ffd_volume(box_center = np.array([50000,6000,3000]), box_length = np.array([5000,5000,5000]), n_control_points = [3,3,3])
        # self.make_ffd_box_mesh()
        # self.move_ffd_pole(0, [1,1,1], [1,0.8,1.1])
        # self.ffd_deform_surfaces()
        # self.regenerateHullHorm()
        # self.visualise_surface()

    def calc_ship_dims(self):   #uses the mesh to find min and max points on x,y and z, used for testing
        points = self.mesh.points()
        x = points[:,0]
        y = points[:,1]
        z = points[:,2]
        # print(f"min x point :{np.min(x)}, min z point: {np.min(z)}")
        self.L = (np.max(x) - np.min(x)) #/1000      #convert to meters
        self.B = (np.max(y) - np.min(y)) #/1000
        self.H = (np.max(z) - np.min(z)) #/ 1000
        self.T = self.H * 0.75 #/ 1000   #0.75 je procijena gaza


    def new_coords_x(self):
        return occhm.lackenby(self._surfaces[0], delta_cp=0.15, delta_p=0.2, station_starts='aft')

    def new_stations(self):
        return occhm.get_station_curve_gp_points(self._surfaces[0])

    def midship(self):
        return occhm.midship_position(self._surfaces[0])

    def coords_x(self):
        return occhm.coord_x(self._surfaces[0])

    def new_stations_opt(self):
        if type(self._surfaces[0]) == TopoDS_Compound:
            surf_explorer = TopExp_Explorer(self._surfaces[0], TopAbs_FACE)
            shape = surf_explorer.Current()
            points = occhm.get_station_curve_gp_points_opt(shape)
        else:
            points = occhm.get_station_curve_gp_points_opt(self._surfaces[0])
        return points

    def bounding_box(self):
        return occhm.bounding_box(self._surfaces[0])

    def water_plane_dim(self):
        return occhm.water_plane_particulars(self._surfaces[0])

    def nurbs_surface(self):
        return occhm.nurbs_surface(self._surfaces[0])

    def printing_values(self):
        print(f'LWL: {occhm.water_plane_particulars(self._surfaces[0])[6]}')
        print(f'BWL: {2*occhm.water_plane_particulars(self._surfaces[0])[7]}')
        print(f'Draft: {occhm.draft}')
        print(f'Volume: {occhm.volume(self._surfaces[0])}')
        print(f'CB: {occhm.block_coeff(self._surfaces[0])}')
        print(f'CP: {occhm.prism_coeff(self._surfaces[0])}')
        print(f'LCB: {occhm.long_center_buoyancy(self._surfaces[0])}')
        print(f'Water plane area: {occhm.water_plane_area(self._surfaces[0])}')
        print(f'BM: {occhm.metacentric_radius(self._surfaces[0])}')
        print(f'VCB: {occhm.vertical_center_buoyancy(self._surfaces[0])}')
        print(f'KM: {occhm.metacentre(self._surfaces[0])}')

    def regenerateHullHorm(self):
        n_surfaces = len(self._surfaces)
        if n_surfaces > 0:
            meshes = []
            for surf in self._surfaces:
                meshes.append(get_open_mesh_from_TopoDS_using_shape_tesselator(surf, mesh_quality=0.2))
            if n_surfaces > 1:  #n of meshes and n of surfaces is the same
                self.mesh = soft_merge_meshes(meshes)
            self.miror_mesh_around_simetry_plan()

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
            xmin = int(1.025 * bb.minCoord[0])
            xmax = int(0.975 * bb.maxCoord[0])
            length = abs(xmin - xmax)
            delta_x = int(length / 35)
            planes = []
            section_curves = []
            for point in list(range(xmin, xmax, delta_x)):
                plane = Geom_Plane(gp_Pnt(point, 0, 0), gp_Dir(1, 0, 0))
                planes += [plane]
            for section in planes:
                curve = BRepAlgoAPI_Section(self._surfaces[0], section, True)
                explorer = TopExp_Explorer(curve.Shape(), TopAbs_EDGE)
                while not isinstance(explorer.Value(), TopoDS_Edge):
                    xmin += 0.01 * delta_x
                    curve = BRepAlgoAPI_Section(self._surfaces[0], Geom_Plane(gp_Pnt(xmin, 0, 0),
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
            if T:
                display, start_display, add_menu, add_function_to_menu = init_display()
                display.DisplayShape(self._surfaces[0], update=True, color='white')
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
        xmin = float(np.around(bb.minCoord[0], 2))
        transom_plane = BRepBuilderAPI_MakeFace(Geom_Plane(gp_Pnt(xmin, 0, 0), gp_Dir(-1, 0, 0)), 1e-6).Shape()
        transom_plane_section = BRepAlgoAPI_Section(self._surfaces[0], transom_plane, True)
        tp_expl = TopExp_Explorer(transom_plane_section.Shape(), TopAbs_VERTEX)
        tp_expl_edge = TopExp_Explorer(transom_plane_section.Shape(), TopAbs_EDGE)
        transom_plane_edge = tp_expl_edge.Value()
        z_vert = []
        vertices_transom = []
        while tp_expl.More():
            z = BRep_Tool_Pnt(tp_expl.Current()).Z()
            z_vert += [z]
            vertices_transom += [BRep_Tool_Pnt(tp_expl.Current())]
            tp_expl.Next()
        cover_high = max(z_vert)
        transom_point = gp_Pnt(xmin, 0, cover_high)
        edge1 = BRepBuilderAPI_MakeEdge(transom_point, vertices_transom[0]).Shape()
        edge2 = BRepBuilderAPI_MakeEdge(vertices_transom[-1], transom_point).Shape()
        transom_wire = BRepBuilderAPI_MakeWire()
        transom_wire.Add(edge1)
        transom_wire.Add(edge2)
        transom_wire.Add(transom_plane_edge)
        transom_wire_profile = transom_wire.Wire()
        transom_plane_face = BRepBuilderAPI_MakeFace(transom_wire_profile).Face()
        new_surf = BRepAlgoAPI_Fuse(self._surfaces[0], transom_plane_face).Shape()
        self._surfaces[0] = new_surf
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()

    def close_cover(self):
        bb = self.bbox
        xmin = float(np.around(bb.minCoord[0], 2))
        water_draft = 8500
        water_plane = BRepBuilderAPI_MakeFace(Geom_Plane(gp_Pnt(0, 0, water_draft), gp_Dir(0, 0, 1)), 1e-6).Shape()
        water_plane_section = BRepAlgoAPI_Section(self._surfaces[0], water_plane, True)
        wp_expl = TopExp_Explorer(water_plane_section.Shape(), TopAbs_VERTEX)
        wp_expl_edge = TopExp_Explorer(water_plane_section.Shape(), TopAbs_EDGE)
        vertices = []
        while wp_expl.More():
            vertices += [BRep_Tool_Pnt(wp_expl.Current())]
            wp_expl.Next()
        water_plane_wire = BRepBuilderAPI_MakeWire()
        while wp_expl_edge.More():
            water_plane_wire.Add(wp_expl_edge.Current())
            wp_expl_edge.Next()
        water_draft_point = gp_Pnt(xmin, 0, water_draft)
        edge1 = BRepBuilderAPI_MakeEdge(water_draft_point, vertices[0]).Shape()
        water_plane_wire.Add(edge1)
        water_plane_profile = water_plane_wire.Wire()
        water_plane_face = BRepBuilderAPI_MakeFace(water_plane_profile).Face()
        new_surf = BRepAlgoAPI_Fuse(self._surfaces[0], water_plane_face).Shape()
        self._surfaces[0] = new_surf
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()

    def scale_1D(self):
        dir_sc = 'x'
        desired_len = 105000
        bb = self.bbox
        length = bb.maxCoord[0] - bb.minCoord[0]
        factor = desired_len / length
        if dir_sc == 'x':
            direction = gp_Dir(1, 0, 0)
        elif dir_sc == 'y':
            direction = gp_Dir(0, 1, 0)
        elif dir_sc == 'z':
            direction = gp_Dir(0, 0, 1)
        else:
            direction = gp_Dir(1, 0, 0)
        center = gp_Pnt()
        z_dir = gp_Ax2(center, direction)
        z_sc = gp_GTrsf()
        z_sc.SetAffinity(z_dir, factor)
        scale_app = BRepBuilderAPI_GTransform(self._surfaces[0], z_sc, False)
        scale_surf = scale_app.Shape()
        self._surfaces[0] = scale_surf
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()

    def full_scale(self):
        factor = 2
        scaling_transform = gp_Trsf()
        scaling_transform.SetScale(gp_Pnt(), factor)
        scale_app = BRepBuilderAPI_Transform(self._surfaces[0], scaling_transform)
        scale_surf = scale_app.Shape()
        self._surfaces[0] = scale_surf
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()
