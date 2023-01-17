from OCC.Core.TopoDS import TopoDS_Shape
import openmesh as om
from OCC.Core.Tesselator import ShapeTesselator
from OCC.Core.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeTorus,
    BRepPrimAPI_MakeSphere)
import numpy as np

def get_open_mesh_from_TopoDS_using_shape_tesselator(shape:TopoDS_Shape,mesh_quality=0.2)->om.TriMesh:
    tess = ShapeTesselator(shape)
    tess.Compute(compute_edges=True,mesh_quality=mesh_quality)
    vpos = tess.GetVerticesPositionAsTuple()
    nvert = tess.ObjGetVertexCount()
    ntria = tess.ObjGetTriangleCount()
    points = np.array(vpos).reshape((ntria*3,3))
    fvi = np.arange(ntria * 3).reshape((ntria, 3))
    mesh = om.TriMesh(points, fvi)
    return mesh

if __name__ == "__main__":
    pass