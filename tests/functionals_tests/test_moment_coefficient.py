import unittest
import numpy as np

import openmdao.api as om

from openaerostruct.functionals.moment_coefficient import MomentCoefficient
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        wing_dict = {"name": "wing", "mesh": np.zeros((2, 7)), "symmetry": True}
        tail_dict = {"name": "tail", "mesh": np.zeros((3, 5)), "symmetry": False}

        surfaces = [wing_dict, tail_dict]

        comp = MomentCoefficient(surfaces=surfaces)

        run_test(self, comp, complex_flag=True, method="cs")

    def test2(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = MomentCoefficient(surfaces=surfaces)

        indep_var_comp = om.IndepVarComp()

        indep_var_comp.add_output("S_ref_total", val=1e4, units="m**2")
        indep_var_comp.add_output("cg", val=np.array([-10.0, 10.0, -10.0]), units="m")

        group.add_subsystem("moment_calc", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.S_ref_total", "moment_calc.S_ref_total")
        group.connect("indep_var_comp.cg", "moment_calc.cg")

        run_test(self, group, complex_flag=True, method="cs")


if __name__ == "__main__":
    unittest.main()
