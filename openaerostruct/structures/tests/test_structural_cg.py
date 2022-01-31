import unittest
import numpy as np

import openmdao.api as om
from openaerostruct.structures.structural_cg import StructuralCG
from openaerostruct.utils.testing import run_test, get_default_surfaces

np.random.seed(1)


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        group = om.Group()

        comp = StructuralCG(surface=surface)

        indep_var_comp = om.IndepVarComp()

        ny = surface["mesh"].shape[1]

        rng = np.random.default_rng(0)
        indep_var_comp.add_output("nodes", val=rng.random((ny, 3)), units="m")
        indep_var_comp.add_output("structural_mass", val=1.0, units="kg")
        indep_var_comp.add_output("element_mass", val=np.ones((ny - 1)), units="kg")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("structural_cg", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, compact_print=False)


if __name__ == "__main__":
    unittest.main()
