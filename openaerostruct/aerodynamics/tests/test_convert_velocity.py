import unittest

import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials

from openaerostruct.aerodynamics.convert_velocity import ConvertVelocity
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        surfaces = get_default_surfaces()

        comp = ConvertVelocity(surfaces=surfaces)

        run_test(self, comp)

    def test_rotation_option_derivatives(self):
        surfaces = get_default_surfaces()

        comp = ConvertVelocity(surfaces=surfaces, rotational=True)

        prob = om.Problem()
        prob.model.add_subsystem("comp", comp)
        prob.setup(force_alloc_complex=True)

        rng = np.random.default_rng(0)
        prob["comp.rotational_velocities"] = rng.random(prob["comp.rotational_velocities"].shape)
        prob["comp.beta"] = 15.0
        prob.run_model()

        check = prob.check_partials(compact_print=True, method="cs", step=1e-40)

        assert_check_partials(check)


if __name__ == "__main__":
    unittest.main()
