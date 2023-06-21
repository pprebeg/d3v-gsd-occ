from optbase import *
from optlib_scipy import ScipyOptimizationAlgorithm


#varijable: pomicanje 4 tocke po x i y, treba li po z?
class form_analysis_model():
    def __init__(self, Hullform = None):
        self.Hullform = Hullform
        box_center = np.array([L/2,B/4,H/2])   #b/4 jer se gleda samo jedna strana
        box_length = [L,B/2,H]
        self.Hullform.make_ffd_volume(self, box_center = box_center, n_control_points = box_length, box_rotation_angles = [0,0,0], knots_to_add = [0,0,0])  #treba mi box center i box len
        self.ffd = self.Hullform._ffd_volumes[0]

        # indeksi redova na pramcu:
        # row1
        bp1_ind = [4,1,0]
        bp2_ind = [4,1,1]
        self.brow1 = [bp1_ind, bp2_ind]     #lista indeksa rowov-a
        # row2
        bp3_ind = [5,1,0]
        bp4_ind = [5,1,1]
        self.brow2 = [bp3_ind, bp4_ind]

        # indeksi redova na krmi:
        # row1
        ap1_ind  = [1,1,0]
        ap2_ind  = [1,1,1]
        self.arow1 = [ap1_ind, ap2_ind]
        # row2
        ap3_ind  = [2,1,0]
        ap4_ind  = [2,1,1]
        self.arow2 = [ap3_ind, ap4_ind]


    #bow functions:
    def set_b_mo_x1(self, mo_x):    #ffd displacement of bow row1
        self.b_mo_x1 = mo_x
        for ind in self.brow1:
            self.ffd.array_mu_x[ind] = self.b_mo_x1

    def set_b_mo_x2(self, mo_x):    #ffd displacement of brow1
        self.b_mo_x2 = mo_x
        for ind in self.brow2:
            self.ffd.array_mu_x[ind] = self.b_mo_x2

    def get_b_mo_x1(self):
        return self.b_mo_x1

    def get_b_mo_x2(self):
        return self.b_mo_x2

    #aft functions
    def set_a_mo_x1(self, mo_x):  # ffd displacement of brow1
        self.a_mo_x1 = mo_x
        for ind in self.arow1:
            self.ffd.array_mu_x[ind] = self.a_mo_x1

    def set_a_mo_x2(self, mo_x):  # ffd displacement of brow1
        self.a_mo_x2 = mo_x
        for ind in self.arow2:
            self.ffd.array_mu_x[ind] = self.a_mo_x2

    def get_a_mo_x1(self):
        return self.a_mo_x1

    def get_a_mo_x2(self):
        return self.a_mo_x2

    # dodaj za y i z ako treba

    # def calc_resistance(self):
    #     self.resistance = somefunction()

    
    # def get_volume(self)
    #     return volume     #koje su granice volumena?

    #... ostali koeficjenti




    def analyze(self):
        pass
      # self.calc_resistance()
      # self.calc_volume()




class form_OptimizationProblem(OptimizationProblem):
    def __init__(self,name=''):
        if name == '':
            name = 'form'
        super().__init__(name)
        am = form_analysis_model()
        b_mo_x1 = DesignVariable('b_mo_x1', CallbackGetSetConnector(am.get_b_mo_x1, am.set_b_mo_x1),0,5)
        b_mo_x2 = DesignVariable('b_mo_x2', CallbackGetSetConnector(am.get_b_mo_x2, am.set_b_mo_x2),0,5)
        #krma  varijable:
        a_mo_x1 = DesignVariable('a_mo_x1', CallbackGetSetConnector(am.get_a_mo_x1, am.set_a_mo_x1),0,5)
        a_mo_x2 = DesignVariable('a_mo_x2', CallbackGetSetConnector(am.get_a_mo_x2, am.set_a_mo_x2),0,5)

        self.add_design_variable(b_mo_x1)
        self.add_design_variable(b_mo_x1)
        self.add_design_variable(b_mo_x1)
        self.add_design_variable(b_mo_x1)

        self.add_objective(DesignObjective('resistance', CallbackGetConnector(am.get_resistance)))

        self.add_constraint(DesignConstraint('volume', CallbackGetConnector(am.get_volume), 25.0, ConstrType.LT))
        self.add_analysis_executor(am)
