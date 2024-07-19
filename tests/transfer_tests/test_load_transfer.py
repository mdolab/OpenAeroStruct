import unittest
import numpy as np

import openmdao.api as om
from openaerostruct.transfer.load_transfer import LoadTransfer
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]
        group = om.Group()

        comp = LoadTransfer(surface=surface)

        indep_var_comp = om.IndepVarComp()

        nx = surface["mesh"].shape[0]
        ny = surface["mesh"].shape[1]

        rng = np.random.default_rng(0)
        indep_var_comp.add_output("def_mesh", val=rng.random((nx, ny, 3)), units="m")
        indep_var_comp.add_output("sec_forces", val=rng.random((nx - 1, ny - 1, 3)), units="N")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("load_transfer", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, compact_print=False)


if __name__ == "__main__":
    unittest.main()
