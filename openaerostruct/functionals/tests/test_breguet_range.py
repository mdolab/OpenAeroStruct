import unittest

from openaerostruct.functionals.breguet_range import BreguetRange
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = BreguetRange(surfaces=surfaces)

        run_test(self, comp, complex_flag=True)


if __name__ == "__main__":
    unittest.main()
