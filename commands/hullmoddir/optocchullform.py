from optbase import AnalysisExecutor,OptimizationProblem
from optbase import DesignVariable, DesignConstraint, DesignObjective
from optlib_scipy import ScipyOptimizationAlgorithm
from optbase import AnalysisResultType, BasicGetSetConnector, BasicGetConnector
import numpy as np

from hullmoddir.occhullform import OCCHullform,CADDeformation

class CallbackPoleXGetSetConnector(BasicGetSetConnector):

    def __init__(self,pole):
        self._pole = pole

    @property
    def value(self):
        xyz = self._pole.XYZ()
        return  xyz[0]

    @value.setter
    def value(self, value):
        xyz = self._pole.XYZ()
        xyz[0]=value
        self._pole.SetXYZ(xyz)

class CallbackPoleXGetConnector(BasicGetSetConnector):

    def __init__(self,pole):
        self._pole = pole

    @property
    def value(self):
        xyz = self._pole.XYZ()
        return  xyz[0]

class OCCHullForm_SurfaceFitCurves_Analysis(AnalysisExecutor):

    def __init__(self, hullform:OCCHullform):
        super().__init__()
        self.hullform:OCCHullform = hullform
        self.poles_pos:np.ndarray = None

    def analyze(self):
        pass
        return AnalysisResultType.OK


class OCCHullForm_SurfaceFitCurves_OptimizationProblem(OptimizationProblem):

    def __init__(self, name, hullform:OCCHullform):
        super().__init__(name)
        # AnalysisExecutor
        am = OCCHullForm_SurfaceFitCurves_Analysis(hullform)
        self.add_analysis_executor(am)

def surface_through_curves_fit(hullform:OCCHullform):
    op= OCCHullForm_SurfaceFitCurves_OptimizationProblem(hullform)
    opt_ctrl = {'maxiter': 100000}
    op.opt_algorithm = ScipyOptimizationAlgorithm('SLSQP_mi=1000','SLSQP',opt_ctrl)
    if True:
        sol = op.optimize()
        op.print_output()
    else:
        sol = op.optimize_and_write(out_folder_path)

