import unittest
import openmdao.api as om
from openmdao.utils.assert_utils import assert_check_partials

from openaerostruct.aerodynamics.eval_mtx import EvalVelMtx
from openaerostruct.utils.testing import run_test, get_default_surfaces, get_ground_effect_surfaces


class Test(unittest.TestCase):

    def test(self):
        surfaces = get_default_surfaces()

        comp = EvalVelMtx(surfaces=surfaces, num_eval_points=2, eval_name='test_name')

        run_test(self, comp, complex_flag=True)

    def test_assembled_jac(self):
        surfaces = get_default_surfaces()

        comp = EvalVelMtx(surfaces=surfaces, num_eval_points=2, eval_name='test_name')

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)

        prob.model.linear_solver = om.DirectSolver(assemble_jac=True)
        prob.model.options['assembled_jac_type'] = 'csc'

        prob.setup(force_alloc_complex=True)

        prob.run_model()

        data = prob.check_partials(compact_print=True, out_stream=None, method='cs', step=1e-40)
        assert_check_partials(data, atol=1e20, rtol=1e-6)

class GroundEffectTest(unittest.TestCase):

    def test(self):
        surfaces = get_ground_effect_surfaces()

        comp = EvalVelMtx(surfaces=surfaces, num_eval_points=2, eval_name='test_name')

        run_test(self, comp, complex_flag=True)

    def test_assembled_jac(self):
        surfaces = get_ground_effect_surfaces()

        comp = EvalVelMtx(surfaces=surfaces, num_eval_points=2, eval_name='test_name')

        prob = om.Problem()
        prob.model.add_subsystem('comp', comp)

        prob.model.linear_solver = om.DirectSolver(assemble_jac=True)
        prob.model.options['assembled_jac_type'] = 'csc'

        prob.setup(force_alloc_complex=True)

        prob.run_model()

        data = prob.check_partials(compact_print=True, out_stream=None, method='cs', step=1e-40)
        assert_check_partials(data, atol=1e20, rtol=1e-6)

if __name__ == '__main__':
    unittest.main()
