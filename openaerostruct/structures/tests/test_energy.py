import unittest

from openaerostruct.structures.energy import Energy
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surface = get_default_surfaces()[0]

        comp = Energy(surface=surface)

        run_test(self, comp)


if __name__ == "__main__":
    unittest.main()
