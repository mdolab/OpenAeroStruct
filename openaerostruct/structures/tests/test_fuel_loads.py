import unittest

import openmdao.api as om
from openaerostruct.structures.fuel_loads import FuelLoads
from openaerostruct.utils.testing import run_test, get_default_surfaces
import numpy as np


class Test(unittest.TestCase):
    def test_0(self):
        surface = get_default_surfaces()[0]

        comp = FuelLoads(surface=surface)

        group = om.Group()

        indep_var_comp = om.IndepVarComp()

        nodesval = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]], dtype=complex)

        indep_var_comp.add_output("nodes", val=nodesval, units="m")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("load", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, atol=1e-2, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
