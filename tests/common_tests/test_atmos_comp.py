import unittest

from openaerostruct.common.atmos_comp import AtmosComp
from openaerostruct.utils.testing import run_test


class Test(unittest.TestCase):
    def test(self):
        comp = AtmosComp()

        run_test(self, comp, method="fd", atol=1e20, rtol=1e-4)


if __name__ == "__main__":
    unittest.main()
