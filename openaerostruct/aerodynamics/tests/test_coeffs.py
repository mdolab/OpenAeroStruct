import unittest

from openaerostruct.aerodynamics.coeffs import Coeffs
from openaerostruct.utils.testing import run_test


class Test(unittest.TestCase):
    def test(self):
        comp = Coeffs()

        run_test(self, comp)


if __name__ == "__main__":
    unittest.main()
