from optbase import AnalysisExecutor, OptimizationProblem, SimpleInputOutputArrayAnalysisExecutor
from optbase import DesignVariable, DesignConstraint, DesignObjective, NdArrayGetSetConnector, NdArrayGetConnector
from optlib_scipy import ScipyOptimizationAlgorithm
from optbase import AnalysisResultType, BasicGetSetConnector, BasicGetConnector
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCC.Core.gp import gp_Pnt
from itertools import product
import hull_manipulation as hm
import numpy as np

from hullmoddir.occhullform import OCCHullform, CADDeformation


class CallbackPoleXGetSetConnector(BasicGetSetConnector):

    def __init__(self, pole):
        self._pole = pole

    @property
    def value(self):
        x = self._pole.X()
        return x

    @value.setter
    def value(self, value):
        x = value
        self._pole.SetX(x)


class CallbackPoleXGetConnector(BasicGetSetConnector):

    def __init__(self, pole):
        self._pole = pole

    @property
    def value(self):
        x = self._pole.X()
        return x


class OCCHullForm_SurfaceFitCurves_Analysis(SimpleInputOutputArrayAnalysisExecutor):

    def __init__(self, hullform: OCCHullform, hullmanipulation: hm.HullManipulation):
        # super().__init__()
        self.hullform: OCCHullform = hullform
        self.hullmanipulation: hm.HullManipulation(hullform.surface,
                                                   station_starts='fore') = hullmanipulation
        self.poles_pos: np.ndarray = None
        self.new_coords_x = hullmanipulation.lackenby(delta_cp=0.03, delta_p=0.01)
        self.new_stations = hullmanipulation.get_station_curve_gp_points()
        self.nurbs_surface = hullmanipulation.nurbs_surface()
        self.bbox = hullmanipulation.bounding_box()
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        n_poles = n_poles_u * n_poles_v
        for station, coords in self.new_stations.items():
            station_list = list(self.new_stations.keys())
            ind = station_list.index(station)
            new_x = self.new_coords_x[ind]
            for x in coords:
                x.SetX(new_x)
        super().__init__(n_poles, n_poles + 1 - n_poles_v)

    def analyze(self):
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        poles_coordinates = np.zeros(shape=(n_poles_u * n_poles_v, 3))
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            pole = self.nurbs_surface.Pole(u + 1, v + 1)
            poles_coordinates[pole_id, :] = pole.X(), pole.Y(), pole.Z()
            new_pole = gp_Pnt(self.inarray[pole_id], poles_coordinates[pole_id][1], poles_coordinates[pole_id][2])
            self.nurbs_surface.SetPole(u + 1, v + 1, new_pole)
        n_pnts = 10 - 2
        pnt_list = []
        self.outarray[0] = 0
        for lists in self.new_stations.values():
            step = round(len(lists) / n_pnts)
            li = [lists[0]] + [lists[i] for i in list(range(1, len(lists), step))] + [lists[-1]]
            pnt_list += [li]
        for lists in pnt_list:
            for p in lists:
                proj = GeomAPI_ProjectPointOnSurf(p, self.nurbs_surface)
                dist = proj.Distance(1) ** 2
                self.outarray[0] += dist
        return AnalysisResultType.OK


class OCCHullForm_SurfaceFitCurves_OptimizationProblem(OptimizationProblem):

    def __init__(self, name, hullform: OCCHullform, hullmanipulation: hm.HullManipulation):
        super().__init__(name)
        # AnalysisExecutor
        am = OCCHullForm_SurfaceFitCurves_Analysis(hullform, hullmanipulation)
        n_poles_u = am.nurbs_surface.NbUPoles()
        n_poles_v = am.nurbs_surface.NbVPoles()
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            pole = am.nurbs_surface.Pole(u + 1, v + 1)
            if u + 1 > 1 and v + 1 > 1:
                pole_prev = am.nurbs_surface.Pole(u, v)
                self.add_design_variable(DesignVariable('Pole ' + str(pole_id),
                                                        NdArrayGetSetConnector(am.inarray, pole_id),
                                                        pole_prev.X() * 1.1, pole.X() * 1.1))

            else:
                # x starts at fore, positive to aftwards
                self.add_design_variable(DesignVariable('Pole ' + str(pole_id),
                                                        NdArrayGetSetConnector(am.inarray, pole_id),
                                                        pole.X() * 0.9, pole.X() * 1.1))

        self.add_objective(DesignObjective('Sum of sq dist', NdArrayGetConnector(am.outarray, 0)))
        self.add_analysis_executor(am)


def surface_through_curves_fit(hullform: OCCHullform, hullmanipulation: hm.HullManipulation):
    op = OCCHullForm_SurfaceFitCurves_OptimizationProblem(hullform, hullmanipulation)
    opt_ctrl = {'maxiter': 100000}
    op.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000', 'SLSQP', opt_ctrl)
    if True:
        sol = op.optimize()
        op.print_output()
    else:
        sol = op.optimize_and_write(out_folder_path)
