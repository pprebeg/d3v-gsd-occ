import numpy as np
import scipy as sci
from OCC.Core.BRep import BRep_Tool_Pnt
from OCC.Core.BRepAlgoAPI import BRepAlgoAPI_Section
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.GProp import GProp_GProps
from OCC.Core.Geom import Geom_Plane
from OCC.Core.TopAbs import TopAbs_VERTEX
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.gp import gp_Vec, gp_Pnt, gp_Ax3, gp_Trsf, gp_Dir
from OCC.Display.SimpleGui import init_display

import occhullmanipulation as occhm

props_orig = GProp_GProps(gp_Pnt(0, 0, 0))
draft = occhm.draft
bbox = occhm.bounding_box


def apply_lackenby(shape, desired_cp, cp):
    sub_faces = occhm.get_sections(shape)
    delta_cp = desired_cp-cp
    print(f'Applying Lackenby, cp variation is: {delta_cp:.3f}')
    lack = occhm.lackenby(shape, delta_cp, 0.00)
    new_faces = []
    for face, lack_dx in zip(sub_faces, lack):
        if lack_dx != 0:
            obj = face
            dx = lack_dx
            dy = 0
            dz = 0
            # Get the origin point of the object
            origin = obj.Location().Transformation().TranslationPart()
            # Compute the displacement vector
            disp = gp_Vec(dx, dy, dz)
            gp_Ax3(gp_Pnt(), gp_Dir(disp))
            trsf = gp_Trsf()
            trsf.SetTranslation(disp)
            transformer = BRepBuilderAPI_Transform(trsf)
            # Apply the transformation to the surface
            transformer.Perform(face)
            # Get the mirrored surface
            new_faces.append(transformer.Shape())
    new_sub = sub_faces[0:len(sub_faces) - len(new_faces)] + new_faces
    if False:
        display, start_display, add_menu, add_function_to_menu = init_display()
        display.DisplayShape(shape, update=True, color='black')
        for w in new_sub:
            display.DisplayShape(w, update=True, color='red')
        start_display()
    return new_sub


def coord_x(stations_list):
    coords_x = []
    sub_faces = stations_list
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        coords_x.append(props_orig.CentreOfMass().X())
    return coords_x


def station_areas(stations_list):
    areas = []
    sub_faces = stations_list
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        areas.append(2 * props_orig.Mass())
    return areas


def stations_bwl(stations_list):
    bwl_list = []
    sub_faces = stations_list
    water_plane = Geom_Plane(gp_Pnt(0, 0, draft), gp_Dir(0, 0, 1))
    coords_x = coord_x(stations_list)
    for face in sub_faces:
        water_plane_section = BRepAlgoAPI_Section(face, water_plane, True).Shape()
        vert_expl = TopExp_Explorer(water_plane_section, TopAbs_VERTEX)
        vert_list = []
        while vert_expl.More():
            vert_list.append(BRep_Tool_Pnt(vert_expl.Current()).Y())
            vert_expl.Next()
        if vert_list:
            bwl = max(vert_list)
            bwl_list.append(2*bwl)
        else:
            pass
    for i in range(len(coords_x) - len(bwl_list)):
        bwl_list.append(0)
    return bwl_list


def water_plane_area(stations_list):
    bwl_list = stations_bwl(stations_list)
    coords_x = coord_x(stations_list)
    return sci.integrate.simps(bwl_list, coords_x)


def volume(stations_list):
    areas = station_areas(stations_list)
    coords_x = coord_x(stations_list)
    vol = np.abs(np.around(sci.integrate.simps(areas, coords_x), 3))
    return vol


def long_center_buoyancy(stations_list):
    areas = station_areas(stations_list)
    coords_x = coord_x(stations_list)
    mom_vol = [area * x for area, x in zip(areas, coords_x)]
    int_mom_vol = sci.integrate.simps(mom_vol, coords_x)
    xcb = int_mom_vol / volume(stations_list)
    return xcb


def stations_vertical_center_buoyancy(stations_list):
    station_zb = []
    sub_faces = stations_list
    for face in sub_faces:
        brepgprop_SurfaceProperties(face, props_orig)
        station_zb.append(props_orig.CentreOfMass().Z())
    return [draft - zb for zb in station_zb]


def vertical_center_buoyancy(stations_list):
    areas = station_areas(stations_list)
    zb_list = stations_vertical_center_buoyancy(stations_list)
    static_mom = [area * z for area, z in zip(areas, zb_list)]
    coords_x = coord_x(stations_list)
    return sci.integrate.simps(static_mom, coords_x) / volume(stations_list)


def metacentric_radius(stations_list):
    bwl_list = [bwl/2 for bwl in stations_bwl(stations_list)]
    coords_x = coord_x(stations_list)
    i_xx = 0  # inertia in around x-axis
    for ind in range(len(bwl_list)):
        if ind > 0:
            i_xx += ((coords_x[ind] - coords_x[ind - 1]) * (bwl_list[ind] + bwl_list[ind - 1]) ** 3) / 12
    return i_xx / volume(stations_list)


def metacentre(stations_list):
    areas = station_areas(stations_list)
    max_area_face_index = areas.index(max(areas))
    max_area_face = stations_list[max_area_face_index]
    bbox_face = bbox(max_area_face)
    keel = bbox_face[4]
    kb = vertical_center_buoyancy(stations_list) - keel
    bm = metacentric_radius(stations_list)
    km = bm + kb
    return km


def long_center_floatation(stations_list):
    bwl_list = stations_bwl(stations_list)
    coords_x = coord_x(stations_list)
    moment_x = 0
    for ind in range(len(bwl_list)):
        if ind > 0:
            moment_x += ((coords_x[ind] - coords_x[ind - 1]) * ((bwl_list[ind] + bwl_list[ind - 1]) / 2)) * (
                    coords_x[ind - 1] + (coords_x[ind] - coords_x[ind - 1]) / 2)
    return moment_x/water_plane_area(stations_list)

