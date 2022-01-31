import unittest
import numpy as np
import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials
from openaerostruct.common.reynolds_comp import ReynoldsComp


class Test(unittest.TestCase):
    def test_reynolds_derivs(self):
        comp = ReynoldsComp()

        prob = om.Problem()
        prob.model.add_subsystem("comp", comp, promotes=["*"])
        prob.setup(force_alloc_complex=True)

        rng = np.random.default_rng(0)
        prob["rho"] = rng.random()
        prob["mu"] = rng.random()
        prob["v"] = rng.random()
        prob.run_model()

        check = prob.check_partials(compact_print=True, method="cs", step=1e-40)

        assert_check_partials(check)


if __name__ == "__main__":
    unittest.main()
