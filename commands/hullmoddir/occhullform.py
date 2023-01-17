from hullformdir.hullform import *
from typing import Set,Dict,List
from OCC.Core.IGESControl import IGESControl_Reader
from OCC.Core.TopAbs import TopAbs_WIRE,TopAbs_FACE
from OCC.Core.TopoDS import TopoDS_Shape
import csv
from hullmoddir.occhelperfun import get_open_mesh_from_TopoDS_using_shape_tesselator
def get_relative_path(full_path,referent_file):
    relative_path=''
    return relative_path

def get_full_path(file_name,referent_file):
    full_path=''
    return full_path

class OCCHullform(HullForm):
    def __init__(self, fileName,name="",translation=np.zeros(3)):
        super().__init__(fileName, name,translation)
        self._shipdata:Dict[str,object]={}
        self._surfaces= []
        self._curves = []
        occ_entities = None
        filename_no_ext, file_extension = os.path.splitext(fileName)
        if file_extension == ".hoc":
            self._shipdata = self.read_hoc_file()
            surface_path = self._shipdata.get('surface_file_relative_path')
            if surface_path is not None:
                full_path = get_full_path(surface_path,self.filename)
                occ_entities = self.read_hoc_file(full_path)
        elif file_extension == ".igs":
            occ_entities = self.read_igs_file(self.filename)
        if occ_entities is not None:
            self._surfaces= occ_entities[1]
            self._curves = occ_entities[0]
        self.regenerateHullHorm()

    def regenerateHullHorm(self):
        if len(self._surfaces) > 0:
            self.mesh = get_open_mesh_from_TopoDS_using_shape_tesselator(self._surfaces[0],mesh_quality=0.05)

    def modify_form(self):
        #do some change with poles
        self.regenerateHullHorm()
        self.emit_geometries_rebuild()
        print('Hull form changed')

    def read_hoc_file(self,file_path:str):
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

    def read_igs_file(self,file_path):
        curves = []
        surfaces=[]
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
                elif shp.ShapeType() ==TopAbs_FACE:
                    surfaces.append(shp)
        return (curves,surfaces)