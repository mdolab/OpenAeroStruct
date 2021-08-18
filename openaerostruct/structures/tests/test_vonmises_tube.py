import unittest
import numpy as np

import openmdao.api as om
from openaerostruct.structures.vonmises_tube import VonMisesTube
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        group = om.Group()

        comp = VonMisesTube(surface=surface)

        indep_var_comp = om.IndepVarComp()

        ny = surface["mesh"].shape[1]

        #  define the spar with y out the wing
        nodes = np.zeros((ny, 3))
        nodes[:, 0] = np.linspace(0, 0.01, ny)
        nodes[:, 1] = np.linspace(0, 1, ny)

        radius = 0.01 * np.ones((ny - 1))

        disp = np.zeros((ny, 6))
        for i in range(6):
            disp[:, i] = np.linspace(0, 0.001, ny)

        indep_var_comp.add_output("nodes", val=nodes, units="m")
        indep_var_comp.add_output("radius", val=radius, units="m")
        indep_var_comp.add_output("disp", val=disp, units="m")

        group.add_subsystem("vm_comp", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.nodes", "vm_comp.nodes")
        group.connect("indep_var_comp.radius", "vm_comp.radius")
        group.connect("indep_var_comp.disp", "vm_comp.disp")

        run_test(self, group, complex_flag=True, compact_print=True, method="cs", step=1e-40, atol=2e-4, rtol=1e-8)


if __name__ == "__main__":
    unittest.main()
