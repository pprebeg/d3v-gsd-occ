import OCC
import numpy as np
import openmesh as om
from pygem.cad import *
from copy import copy


from OCC.Core.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCC.Core.Tesselator import *
# from OCC.Extend.TopologyUtils import TopologyExplorer
from OCC.Extend.TopologyUtils import *
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.GeomAbs import GeomAbs_BSplineSurface
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge

from OCC.Core.gp import gp_Pnt, gp_Vec
from OCC.Core.GeomFill import (
    GeomFill_BSplineCurves,
    GeomFill_StretchStyle,
    GeomFill_CoonsStyle,
    GeomFill_CurvedStyle,
)
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.Geom import Geom_BSplineCurve

from OCC.Display.SimpleGui import init_display

from OCC.Extend.ShapeFactory import point_list_to_TColgp_Array1OfPnt, make_face


def make_block(block_dims=np.array([20, 6, 3]), move_vector=np.array([0, 0, 0])):
    mesh = om.TriMesh()
    axes = []
    # stvara 2 tocke na svakoj osi
    for dim in block_dims:
        axes.append(np.linspace(0, dim, 2))

    block_corners = np.asarray(np.meshgrid(*axes)).T.reshape(8, 3)
    block_corners += move_vector

    # shrink block to try to avoid errors in boolean operations
    # block_corners[[4,5,6,7]] += np.array([0,0,-0.00005])			#z+ face is lowered
    # block_corners[[0,1,2,3]] += np.array([0,0,+0.00005])			#z- face is increased

    corner_vertices = []
    for corner in block_corners:
        corner_vertices.append(mesh.add_vertex(corner))

    # x+face
    mesh.add_face(corner_vertices[2], corner_vertices[3], corner_vertices[6])
    mesh.add_face(corner_vertices[3], corner_vertices[7], corner_vertices[6])

    # x-face
    mesh.add_face(corner_vertices[0], corner_vertices[4], corner_vertices[1])
    mesh.add_face(corner_vertices[1], corner_vertices[4], corner_vertices[5])

    # y+face
    mesh.add_face(corner_vertices[3], corner_vertices[1], corner_vertices[5])
    mesh.add_face(corner_vertices[3], corner_vertices[5], corner_vertices[7])

    # y-face
    mesh.add_face(corner_vertices[0], corner_vertices[2], corner_vertices[4])
    mesh.add_face(corner_vertices[2], corner_vertices[6], corner_vertices[4])

    # z+face
    mesh.add_face(corner_vertices[4], corner_vertices[6], corner_vertices[5])
    mesh.add_face(corner_vertices[6], corner_vertices[7], corner_vertices[5])

    # z-face
    mesh.add_face(corner_vertices[2], corner_vertices[0], corner_vertices[1])
    mesh.add_face(corner_vertices[2], corner_vertices[1], corner_vertices[3])

    return mesh
def soft_merge_meshes(meshes, vh_idx_to_sync_list=None):  # meshes je lista sa meshevima, vh_idx_to_sync_list sa lista isog lena ko meshes, svaka sadržava array sa vh_idx koji zelimo syncat
    points = np.empty((0, 3))
    merged_fvi = np.empty((0, 3))

    if vh_idx_to_sync_list is None:
        for mesh in meshes:
            mesh_fvi = mesh.face_vertex_indices()
            if mesh_fvi.size == 0:
                continue
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            points = np.append(points, mesh.points(), axis=0)

        return om.TriMesh(points, merged_fvi)

    else:
        synced_vh_idx = []
        for i in range(len(meshes)):
            mesh = meshes[i]
            mesh_fvi = mesh.face_vertex_indices()
            merged_fvi = np.append(merged_fvi, mesh_fvi + points.shape[0],
                                   axis=0)  # +points.shape[0] je tu da poreda face_vertex_indices sa njihovim indexom u novom arrayu
            synced_vh_idx.append(vh_idx_to_sync_list[i] + points.shape[0])
            points = np.append(points, mesh.points(), axis=0)

        return (om.TriMesh(points, merged_fvi), synced_vh_idx)
def move_mesh(mesh, move_vector):
    return om.TriMesh(mesh.points() + move_vector, mesh.face_vertex_indices())

class ffd_maker():
    def __init__(self):
        self._ffd_volumes = {}
        self._current_volume_id = 0

    def make_ffd_volume(self, box_center, box_length, n_control_points = [2,2,2], box_rotation_angles = [0,0,0], knots_to_add = [0,0,0]):
        ffd_volume = FFD(n_control_points, knots_to_add[0],knots_to_add[1],knots_to_add[2])
        ffd_volume.box_length = box_length
        ffd_volume.box_origin = (box_center - box_length/2)
        ffd_volume.rot_angle = box_rotation_angles
        # _ffd_volume.control_points()
        self._ffd_volumes[self._current_volume_id] = ffd_volume
        self._current_volume_id += 1

    def move_ffd_pole(self, ffd_id, pole_id:list, move_vector):
        current_ffd_volume = self._ffd_volumes[ffd_id]
        # current_ffd_volume.array_mu_x[pole_id[0], pole_id[1], pole_id[2]] = move_vector[0]
        # current_ffd_volume.array_mu_y[pole_id[0], pole_id[1], pole_id[2]] = move_vector[1]
        # current_ffd_volume.array_mu_z[pole_id[0], pole_id[1], pole_id[2]] = move_vector[2]
        current_ffd_volume.array_mu_x[[*pole_id,]] = move_vector[0]
        current_ffd_volume.array_mu_y[[*pole_id,]] = move_vector[1]
        current_ffd_volume.array_mu_z[[*pole_id,]] = move_vector[2]

    def make_ffd_box_mesh(self, pole_size = 500):    #moza ih ne radi na tocno pravom mjestu pogledaj odakle je origin make_boxa
        pole = make_block(np.array([pole_size,]*3), np.array([0, 0, 0]))
        ffd_meshes = {}
        test_visualisation = []
        for ffd_id, ffd_volume in self._ffd_volumes.items():
            ffd_control_points = ffd_volume.control_points().tolist()
            ffd_meshes = []

            for cpoint in ffd_control_points:
                pole_mesh = copy(pole)
                pole_mesh = move_mesh(pole_mesh, cpoint)
                ffd_meshes.append(pole_mesh)
                test_visualisation.append(pole_mesh)
            ffd_meshes[ffd_id] = ffd_meshes

        self.mesh = soft_merge_meshes(test_visualisation+[self.mesh,])


class PyGemHullform(ffd_maker):
    def __init__(self):
        super().__init__()
        self._deformed_surfaces = []
    def visualise_surface(self):
        display, start_display, add_menu, add_function_to_menu = init_display()
        for i in self._surfaces:
            display.DisplayShape(i, update=True)

        # boxes = []
        # for ffd_volume in self._ffd_volumes.values():
        #     box_origin = ffd_volume.box_origin
        #     box_length = ffd_volume.box_length
        #     point = gp_Pnt(*box_origin)
        #     box_visualisation = BRepPrimAPI_MakeBox(*([point, ] + box_length.tolist())).Shape()
        #     # box_visualisation = BRepPrimAPI_MakeBox(point, box_length[0], box_length[1], box_length[2]).Shape()
        #     boxes.append(box_visualisation)
        #
        # for box in boxes:
        #     display.DisplayShape(box, update=True)

        start_display()

    def ffd_deform_surfaces(self):
        for surface in self._surfaces:
            for ffd in self._ffd_volumes.values():
                deformed_surface =  ffd(surface)       #vraca compound?
                deformed_surface = self.regenerate_surface(deformed_surface)
                self._deformed_surfaces.append(deformed_surface)
        self._surfaces = self._deformed_surfaces
        self.regenerateHullHorm()


    def regenerate_surface(self, surface):
        if isinstance(surface, OCC.Core.TopoDS.TopoDS_Compound):      #nekada vraca compound a nekada surface...
            _compound = surface
            _expl = TopologyExplorer(_compound)
            for _face in _expl.faces():   #uvijek bi trebao biti samo jedan surface u compoundu
                _face_to_regenerate = _face
        else:
            _face_to_regenerate = surface

        surf = BRepAdaptor_Surface(_face_to_regenerate, True)
        bsrf = surf.BSpline()
        face = make_face(bsrf, 1e-6)
        return face


if __name__ == "__main__":
    maker = ffd_maker()
    maker.make_ffd_volume(np.array([0,0,0]), np.array([1,1,1]))
    # print(maker.ffd_volumes[0].control_points())






















