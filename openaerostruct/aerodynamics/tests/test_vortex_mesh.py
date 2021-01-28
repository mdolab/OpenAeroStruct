import unittest

from openaerostruct.aerodynamics.vortex_mesh import VortexMesh
from openaerostruct.utils.testing import run_test, get_default_surfaces, get_ground_effect_surfaces

import numpy as np
np.set_printoptions(linewidth=200)


class Test(unittest.TestCase):

    def test(self):
        surfaces = get_default_surfaces()

        comp = VortexMesh(surfaces=surfaces)

        run_test(self, comp)

    def test_groundplane(self):
        surfaces = get_ground_effect_surfaces()

        comp = VortexMesh(surfaces=surfaces)

        run_test(self, comp, atol=1e6)

if __name__ == '__main__':
    unittest.main()
