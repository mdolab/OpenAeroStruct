import unittest

from openaerostruct.aerodynamics.get_vectors import GetVectors
from openaerostruct.utils.testing import run_test, get_default_surfaces, get_ground_effect_surfaces


class Test(unittest.TestCase):

    def test(self):
        surfaces = get_default_surfaces()

        comp = GetVectors(surfaces=surfaces, num_eval_points=10, eval_name='test_name')

        run_test(self, comp)

    def test_groundplane(self):
        surfaces = get_ground_effect_surfaces()

        comp = GetVectors(surfaces=surfaces, num_eval_points=10, eval_name='test_name')

        run_test(self, comp)


if __name__ == '__main__':
    unittest.main()
