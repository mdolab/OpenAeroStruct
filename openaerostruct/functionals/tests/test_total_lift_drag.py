import unittest

import openmdao.api as om

from openaerostruct.functionals.total_lift_drag import TotalLiftDrag
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        wing_dict = {"name": "wing", "num_y": 7, "num_x": 2, "symmetry": True}
        tail_dict = {"name": "tail", "num_y": 5, "num_x": 3, "symmetry": False}

        surfaces = [wing_dict, tail_dict]

        comp = TotalLiftDrag(surfaces=surfaces)

        run_test(self, comp)

    # This is known to have some issues for sufficiently small values of S_ref_total
    # There is probably a derivative bug somewhere in the moment_coefficient.py calcs
    def test2(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = TotalLiftDrag(surfaces=surfaces)

        indep_var_comp = om.IndepVarComp()

        indep_var_comp.add_output("S_ref_total", val=10.0, units="m**2")

        group.add_subsystem("moment_calc", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.S_ref_total", "moment_calc.S_ref_total")

        run_test(self, group)


if __name__ == "__main__":
    unittest.main()
