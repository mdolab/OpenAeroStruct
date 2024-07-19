import unittest
import numpy as np

import openmdao.api as om

from openaerostruct.aerodynamics.solve_matrix import SolveMatrix
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        group = om.Group()
        comp = SolveMatrix(surfaces=surfaces)

        indep_var_comp = om.IndepVarComp()

        system_size = 0
        for surface in surfaces:
            nx = surface["mesh"].shape[0]
            ny = surface["mesh"].shape[1]
            system_size += (nx - 1) * (ny - 1)

        indep_var_comp.add_output("rhs", val=np.ones((system_size)), units="m/s")
        indep_var_comp.add_output("mtx", val=np.identity((system_size)), units="1/m")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("solve_matrix", comp, promotes=["*"])

        run_test(self, group)


if __name__ == "__main__":
    unittest.main()
