import unittest

from openaerostruct.aerodynamics.eval_velocities import EvalVelocities
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = EvalVelocities(surfaces=surfaces, eval_name="TestEval", num_eval_points=11)

        run_test(self, comp)


if __name__ == "__main__":
    unittest.main()
