import unittest

from openaerostruct.functionals.sum_areas import SumAreas
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = SumAreas(surfaces=surfaces)

        run_test(self, comp)


if __name__ == "__main__":
    unittest.main()
