import unittest

from openaerostruct.structures.failure_ks import FailureKS
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        # turn down some of these properties, so the absolute deriv error isn't magnified
        surface["E"] = 7
        surface["G"] = 3
        surface["yield"] = 0.02

        comp = FailureKS(surface=surface)

        run_test(self, comp, complex_flag=True, step=1e-40, method="cs", compact_print=False)


if __name__ == "__main__":
    unittest.main()
