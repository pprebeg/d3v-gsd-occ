import numpy as np
from optbase import AnalysisExecutor, OptimizationProblem, SingleobjectiveOptimizationOutput,ConstrType
from optbase import NdArrayGetConnector
from optbase import DesignVariable, DesignConstraint, DesignObjective, NdArrayGetSetConnector, NdArrayGetConnector
from optlib_scipy import ScipyOptimizationAlgorithm
from optlib_pymoo_proto import PymooOptimizationAlgorithmMulti, PymooOptimizationAlgorithmSingle
from optbase import AnalysisResultType, BasicGetSetConnector, BasicGetConnector
from OCC.Core.GeomAPI import GeomAPI_ProjectPointOnSurf
from itertools import product
from hullmoddir.occhullform import OCCHullform
from geometry_extend import GeometryExtension
import openmesh as om
from OCC.Display.SimpleGui import init_display


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
        return (self._pole.X() - self._prev_pole.X())/self._pole.X()


class SurfDistGetCallbackConnector(BasicGetConnector):
    """Class used to connect and normalize objective functions and such."""

    def __init__(self, dist):
        pass
        self._dist = dist

    @property
    def value(self):
        return self._dist()


class OCCHullForm_SurfaceFitCurves_Analysis(AnalysisExecutor):

    def __init__(self, hullform: OCCHullform):
        super().__init__()
        self.hullform: OCCHullform = hullform
        #self.new_coords_x = hullform.new_coords_x()
        self.new_stations = hullform.new_stations()
        self.nurbs_surface = hullform.nurbs_surface()
        self.bbox = hullform.bounding_box()
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        self.poles = []
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            pole = self.nurbs_surface.Pole(u + 1, v + 1)
            self.poles.append(pole)
        self._sum_sqr = self.calc_SumSqr_arrays()
        #for station, coords in self.new_stations.items():
        #    station_list = list(self.new_stations.keys())
        #    ind = station_list.index(station)
        #    new_x = self.new_coords_x[ind]
        #    for x in coords:
        #        x.SetX(new_x)

    def analyze(self):
        # display, start_display, add_menu, add_function_to_menu = init_display()
        # display.DisplayShape(self.nurbs_surface, update=True, color='red')
        n_poles_u = self.nurbs_surface.NbUPoles()
        n_poles_v = self.nurbs_surface.NbVPoles()
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            if pole_id > n_poles_v - 1:
                pole = self.poles[pole_id]
                self.nurbs_surface.SetPole(u + 1, v + 1, pole)
                #print(pole_id, pole.X())
        #print('Analysis called, dist =', self.dists())
        # display.DisplayShape(self.nurbs_surface, update=True, color='black')
        # start_display()
        self.calc_SumSqr_arrays(self._sum_sqr)
        print('Analysis called -1, dist =', self._sum_sqr[self._sum_sqr.size-1])

        return AnalysisResultType.OK

    def dists(self):
        distance = 0
        n_pnts = 10 - 2
        pnt_list = []
        for lists in self.new_stations.values():
            step = round(len(lists) / n_pnts)
            li = [lists[0]] + [lists[i] for i in list(range(1, len(lists), step))] + [lists[-1]]
            pnt_list += [li]
        for lists in pnt_list:
            for p in lists:
                proj = GeomAPI_ProjectPointOnSurf(p, self.nurbs_surface)
                dist = proj.Distance(1) ** 4
                distance += dist
        print('Dist called =', distance)
        return distance

    def calc_SumSqr_arrays(self,ss_array:np.ndarray = None):
        distance = 0
        n_pnts = 10 - 2
        pnt_list = []
        for lists in self.new_stations.values():
            step = round(len(lists) / n_pnts)
            li = [lists[0]] + [lists[i] for i in list(range(1, len(lists), step))] + [lists[-1]]
            pnt_list += [li]
        if ss_array is None:
            out_ss_array = np.ndarray(len(pnt_list)+1)
        else:
            out_ss_array = ss_array
        i = 0
        for lists in pnt_list:
            dist_station = 0.0
            for p in lists:
                proj = GeomAPI_ProjectPointOnSurf(p, self.nurbs_surface)
                dist = proj.Distance(1) ** 4
                dist_station+= dist
                distance += dist
            out_ss_array[i] = dist_station
            i+=1
        out_ss_array[i] = distance
        return out_ss_array

    # def surf(self):
    #     surface = self.nurbs_surface
    #     n_poles_u = surface.NbUPoles()
    #     n_poles_v = surface.NbVPoles()
    #     pole_ids = product(range(n_poles_u), range(n_poles_v))
    #     for pole_id, (u, v) in enumerate(pole_ids):
    #         if pole_id > n_poles_v - 1:
    #             pole = self.poles[pole_id]
    #             surface.SetPole(u + 1, v + 1, pole)
    #     return surface


class OCCHullForm_SurfaceFitCurves_OptimizationProblem(OptimizationProblem):

    def __init__(self, name, hullform: OCCHullform):
        super().__init__(name)
        # AnalysisExecutor
        am = OCCHullForm_SurfaceFitCurves_Analysis(hullform)
        n_poles_u = am.nurbs_surface.NbUPoles()
        n_poles_v = am.nurbs_surface.NbVPoles()
        dist = am.dists()
        # surf = am.surf()
        # n_poles_u = surf.NbUPoles()
        # n_poles_v = surf.NbVPoles()
        dx=am.bbox[6]
        dx_bnd= dx/25
        pole_ids = product(range(n_poles_u), range(n_poles_v))
        for pole_id, (u, v) in enumerate(pole_ids):
            if pole_id > n_poles_v - 1:
                pole = am.poles[pole_id]
                if pole.X() > 0.01 or True:
                    pole_con = CallbackPoleXGetSetConnector(pole)
                    #lb = pole.X() * 0.9
                    #ub = pole.X() * 1.1
                    lb = min(pole.X() * 0.9, (pole.X()-dx_bnd))
                    ub = max(pole.X() * 1.1, (pole.X()+dx_bnd))
                    self.add_design_variable(DesignVariable(str(f' Pole {pole_id}: {u + 1}|{v + 1}'),
                                                            pole_con, min(lb, ub), max(lb, ub),True))
                    pole_prev = am.poles[pole_id - n_poles_v]
                    const_con = PoleXDiffGetCallbackConnector(pole_prev, pole)
                    self.add_constraint(DesignConstraint(
                        str(f' Diff: ({pole_id}: {u + 1}|{v + 1}) - ({pole_id - n_poles_v}: {u + 1}|{v})'),
                        const_con, 0))
        if False:
            dist_conn = SurfDistGetCallbackConnector(am.dists)
            self.add_constraint(DesignConstraint('Con Sum of sq dist', dist_conn,0.001,ConstrType.LT))
            self.add_objective(DesignObjective('Obj Sum of sq dist', dist_conn))
        else:
            for i in range(am._sum_sqr.size-1):
                conn= NdArrayGetConnector(am._sum_sqr,i)
                self.add_constraint(DesignConstraint('C_SS_st_{0}'.format(i+1), conn, 100.0, ConstrType.LT))
            conn = NdArrayGetConnector(am._sum_sqr, am._sum_sqr.size-1)
            self.add_objective(DesignObjective('O_SS_tot', conn))
        self.add_analysis_executor(am)


F = False
T = True


def surface_through_curves_fit(hullform: OCCHullform):
    init_geo = GeometryExtension('opt_geo_initial')
    hullform.regenerateHullHorm()
    fvi = hullform.mesh.fv_indices()
    points = hullform.mesh.points()
    init_geo.mesh = om.TriMesh(points, fvi)
    init_geo.emit_geometry_built()
    op = OCCHullForm_SurfaceFitCurves_OptimizationProblem('Hullmod', hullform)
    if False:
        op._init_opt_problem()
        x0 = op.get_initial_design()
        op.evaluate(x0)
        op._opt_output = SingleobjectiveOptimizationOutput('InitialDesign', op.name, 0, 1, op.get_current_sol())
        op.print_output()

        x0 = op.get_initial_design(True)
        op.evaluate(x0)
        op._opt_output = SingleobjectiveOptimizationOutput('RandDesign', op.name, 0, 1, op.get_current_sol())
        op.print_output()
        return
    if F:  # nelder-mead
        termination = ('n_eval', 1000)
        nm_ctrl = {'termination': termination}
        op.opt_algorithm = PymooOptimizationAlgorithmSingle('nelder-mead_default', 'nelder-mead', alg_ctrl=nm_ctrl)
    if F:#ga
        pop_size = 10
        num_iter = 30
        max_evaluations = pop_size * num_iter
        termination = ('n_eval', max_evaluations)
        ga_ctrl = {'pop_size': pop_size,'termination': termination}
        op.opt_algorithm = PymooOptimizationAlgorithmSingle('ga_default', 'ga', alg_ctrl=ga_ctrl)
    if T:  # Sequential Least SQuares Programming
        opt_ctrl = {'disp': True,
                    'maxiter': 1000,
                    'eps': 0.01}
        op.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000', 'SLSQP', opt_ctrl)
    if True:
        x0 = op.get_initial_design()
        x0 = np.array(x0)*1.01
        sol = op.optimize(x0)
        print(sol)
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
