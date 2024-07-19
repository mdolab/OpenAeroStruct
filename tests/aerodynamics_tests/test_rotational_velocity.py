import unittest

import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials

from openaerostruct.aerodynamics.rotational_velocity import RotationalVelocity
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = RotationalVelocity(surfaces=surfaces)

        run_test(self, comp)

    def test_rotation_option_derivatives(self):
        surfaces = get_default_surfaces()

        comp = RotationalVelocity(surfaces=surfaces)

        prob = om.Problem()
        prob.model.add_subsystem("comp", comp)
        prob.setup(force_alloc_complex=True)

        prob["comp.omega"] = np.array([0.3, 0.4, -0.1])
        prob["comp.cg"] = np.array([0.1, 0.6, 0.4])
        rng = np.random.default_rng(0)
        prob["comp.coll_pts"] = rng.random(prob["comp.coll_pts"].shape)
        prob.run_model()

        check = prob.check_partials(compact_print=True, method="cs", step=1e-40)

        assert_check_partials(check)


if __name__ == "__main__":
    unittest.main()
