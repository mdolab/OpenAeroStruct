import unittest
import numpy as np

import openmdao.api as om
from openaerostruct.structures.weight import Weight
from openaerostruct.utils.testing import run_test, get_default_surfaces


np.random.seed(314)


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]
        ny = surface["mesh"].shape[1]

        group = om.Group()

        ivc = om.IndepVarComp()
        rng = np.random.default_rng(0)
        ivc.add_output("nodes", val=rng.random((ny, 3)), units="m")

        comp = Weight(surface=surface)

        group.add_subsystem("ivc", ivc, promotes=["*"])
        group.add_subsystem("comp", comp, promotes=["*"])

        run_test(self, group, compact_print=False, complex_flag=True)


if __name__ == "__main__":
    unittest.main()
