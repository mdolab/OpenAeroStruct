import unittest
import numpy as np

import openmdao.api as om
from openaerostruct.structures.vonmises_wingbox import VonMisesWingbox
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        # turn down some of these properties, so the absolute deriv error isn't magnified
        surface["E"] = 7
        surface["G"] = 3
        surface["yield"] = 0.02

        surface["strength_factor_for_upper_skin"] = 1.0

        comp = VonMisesWingbox(surface=surface)

        group = om.Group()

        indep_var_comp = om.IndepVarComp()

        ny = surface["mesh"].shape[1]

        nodesval = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])

        indep_var_comp.add_output("nodes", val=nodesval, units="m")
        indep_var_comp.add_output("disp", val=np.ones((ny, 6)), units="m")
        indep_var_comp.add_output("Qz", val=np.ones((ny - 1)), units="m**3")
        indep_var_comp.add_output("Iz", val=np.ones((ny - 1)))
        indep_var_comp.add_output("J", val=np.ones((ny - 1)), units="m**4")
        indep_var_comp.add_output("A_enc", val=np.ones((ny - 1)), units="m**2")
        indep_var_comp.add_output("spar_thickness", val=np.ones((ny - 1)), units="m")
        indep_var_comp.add_output("skin_thickness", val=np.ones((ny - 1)), units="m")
        indep_var_comp.add_output("htop", val=np.ones((ny - 1)), units="m")
        indep_var_comp.add_output("hbottom", val=np.ones((ny - 1)), units="m")
        indep_var_comp.add_output("hfront", val=np.ones((ny - 1)), units="m")
        indep_var_comp.add_output("hrear", val=np.ones((ny - 1)), units="m")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("vonmises_wingbox", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, step=1e-8, atol=2e-5, compact_print=True)


if __name__ == "__main__":
    unittest.main()
