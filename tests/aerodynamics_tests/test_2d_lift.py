import unittest

from openaerostruct.aerodynamics.lift_coeff_2D import LiftCoeff2D
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = LiftCoeff2D(surface=surfaces[0])

        run_test(self, comp, method="cs", complex_flag=True)

    def test2(self):
        surfaces = get_default_surfaces()

        comp = LiftCoeff2D(surface=surfaces[1])

        run_test(self, comp, method="cs", complex_flag=True)


if __name__ == "__main__":
    unittest.main()
