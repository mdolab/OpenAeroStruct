import unittest

from openaerostruct.aerodynamics.panel_forces_surf import PanelForcesSurf
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = PanelForcesSurf(surfaces=surfaces)

        run_test(self, comp)


if __name__ == "__main__":
    unittest.main()
