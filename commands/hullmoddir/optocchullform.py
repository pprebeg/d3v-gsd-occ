from optbase import AnalysisExecutor, OptimizationProblem,SingleobjectiveOptimizationOutput
from optbase import DesignVariable, DesignConstraint, DesignObjective, NdArrayGetSetConnector, NdArrayGetConnector
from optlib_scipy import ScipyOptimizationAlgorithm
from optlib_pymoo_proto import PymooOptimizationAlgorithmMulti, PymooOptimizationAlgorithmSingle
from optbase import AnalysisResultType, BasicGetSetConnector, BasicGetConnector
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from itertools import product
from hullmoddir.occhullform import OCCHullform
from geometry_extend import GeometryExtension
import openmesh as om


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


class PoleXDiffGetCallbackConnector(BasicGetConnector):
    """Class used to connect and normalize objective functions and such."""

    def __init__(self, prev_pole, pole):
        pass
        self._prev_pole = prev_pole
        self._pole = pole

    @property
    def value(self):
        return self._pole.X() - self._prev_pole.X()


class SurfDistGetCallbackConnector(BasicGetConnector):
    """Class used to connect and normalize objective functions and such."""

    def __init__(self, dist):
        pass
        self._dist = dist

    @property
    def value(self):
        return self._dist


class PoleXInitGetCallbackConnector(BasicGetConnector):
    """Class used to connect and normalize objective functions and such."""

    def __init__(self, pole):
        pass
        self._pole = pole

    @property
    def value(self):
        return self._pole.X()


class OCCHullForm_SurfaceFitCurves_Analysis(AnalysisExecutor):

    def __init__(self, hullform: OCCHullform):
        super().__init__()
        self.hullform: OCCHullform = hullform
        self.new_coords_x = hullform.new_coords_x()
        self.new_stations = hullform.new_stations()
        self.nurbs_surface = hullform.nurbs_surface()
        self.bbox = hullform.bounding_box()
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        self.poles = []
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        self.dists = 0
        for pole_id, (u, v) in enumerate(pole_ids):
            pole = self.nurbs_surface.Pole(u + 1, v + 1)
            self.poles.append(pole)
        for station, coords in self.new_stations.items():
            station_list = list(self.new_stations.keys())
            ind = station_list.index(station)
            new_x = self.new_coords_x[ind]
            for x in coords:
                x.SetX(new_x)

    def analyze(self):
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            pole = self.poles[pole_id]
            self.nurbs_surface.SetPole(u + 1, v + 1, pole)
        n_pnts = 10 - 2
        pnt_list = []
        for lists in self.new_stations.values():
            step = round(len(lists) / n_pnts)
            li = [lists[0]] + [lists[i] for i in list(range(1, len(lists), step))] + [lists[-1]]
            pnt_list += [li]
        for lists in pnt_list:
            for p in lists:
                proj = GeomAPI_ProjectPointOnSurf(p, self.nurbs_surface)
                dist = proj.Distance(1) ** 2
                self.dists += dist
        return AnalysisResultType.OK


class OCCHullForm_SurfaceFitCurves_OptimizationProblem(OptimizationProblem):

    def __init__(self, name, hullform: OCCHullform):
        super().__init__(name)
        # AnalysisExecutor
        am = OCCHullForm_SurfaceFitCurves_Analysis(hullform)
        n_poles_u = am.nurbs_surface.NbUPoles()
        n_poles_v = am.nurbs_surface.NbVPoles()
        dist = am.dists
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            if pole_id > n_poles_v-1:
                pole = am.poles[pole_id]
                pole_con = CallbackPoleXGetSetConnector(pole)
                lb = pole.X() * 0.9
                ub = pole.X() * 1.1
                self.add_design_variable(DesignVariable(str(f' Pole {pole_id}: {u+1}|{v+1}'),
                                                        pole_con, min(lb, ub), max(lb, ub)))
                pole_prev = am.poles[pole_id - n_poles_v]
                const_con = PoleXDiffGetCallbackConnector(pole_prev, pole)
                self.add_constraint(DesignConstraint(
                    str(f' Pole Diff. Constr. - {pole_id - n_poles_v}|{pole_id}'),
                    const_con, 0))
        dist_conn = SurfDistGetCallbackConnector(dist)
        self.add_objective(DesignObjective(' Sum of sq dist', dist_conn))
        self.add_analysis_executor(am)


F = False
T = True


def surface_through_curves_fit(hullform: OCCHullform):
    init_geo= GeometryExtension('opt_geo_initial')
    hullform.regenerateHullHorm()
    fvi = hullform.mesh.fv_indices()
    points = hullform.mesh.points()
    init_geo.mesh = om.TriMesh(points, fvi)
    init_geo.emit_geometry_built()
    op = OCCHullForm_SurfaceFitCurves_OptimizationProblem('Hullmod', hullform)
    op._init_opt_problem()
    x0 = op.get_initial_design()
    op.evaluate(x0)
    op._opt_output = SingleobjectiveOptimizationOutput('InitialDesign', op.name, 0, 1, op.get_current_sol())
    op.print_output()
    if T:  # nelder-mead
        termination = ('n_eval', 20)
        nm_ctrl = {'termination': termination}
        op.opt_algorithm = PymooOptimizationAlgorithmSingle('nelder-mead_default', 'nelder-mead', alg_ctrl=nm_ctrl)

    if F:  # Sequential Least SQuares Programming
        opt_ctrl = {'maxiter': 1000}
        op.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000', 'SLSQP', opt_ctrl)
    if True:
        sol = op.optimize()
        op.print_output()
    else:
        sol = op.optimize_and_write(out_folder_path)
    final_geo = GeometryExtension('opt_geo_final')
    hullform.regenerateHullHorm()
    fvi = hullform.mesh.fv_indices()
    points = hullform.mesh.points()
    final_geo.mesh = om.TriMesh(points, fvi)
    final_geo.emit_geometry_built()
    return sol
