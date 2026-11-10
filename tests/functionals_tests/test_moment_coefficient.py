import unittest
import numpy as np

import openmdao.api as om

from openaerostruct.functionals.moment_coefficient import MomentCoefficient
from openmdao.utils.assert_utils import assert_check_partials

from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):
    def test(self):
        wing_dict = {"name": "wing", "mesh": np.zeros((2, 7)), "symmetry": True}
        tail_dict = {"name": "tail", "mesh": np.zeros((3, 5)), "symmetry": False}

        surfaces = [wing_dict, tail_dict]

        comp = MomentCoefficient(surfaces=surfaces)

        run_test(self, comp, complex_flag=True, method="cs")

    def test2(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = MomentCoefficient(surfaces=surfaces)

        indep_var_comp = om.IndepVarComp()

        indep_var_comp.add_output("S_ref_total", val=1e4, units="m**2")
        indep_var_comp.add_output("cg", val=np.array([-10.0, 10.0, -10.0]), units="m")

        group.add_subsystem("moment_calc", comp)
        group.add_subsystem("indep_var_comp", indep_var_comp)

        group.connect("indep_var_comp.S_ref_total", "moment_calc.S_ref_total")
        group.connect("indep_var_comp.cg", "moment_calc.cg")

        run_test(self, group, complex_flag=True, method="cs")

    def check_repeated_partials(self, surfaces, num_calls=3):
        prob = om.Problem(reports=False)
        prob.model.add_subsystem("comp", MomentCoefficient(surfaces=surfaces))
        prob.setup(force_alloc_complex=True)

        # Use random inputs so that every partial is nonzero and distinguishable from stale values
        rng = np.random.default_rng(0)
        for name, meta in prob.model.comp.get_io_metadata(iotypes="input", return_rel_names=False).items():
            prob.set_val(name, rng.uniform(0.5, 1.5, size=meta["shape"]))
        prob.run_model()

        # Each check_partials call triggers a fresh compute_partials call
        for _ in range(num_calls):
            check = prob.check_partials(compact_print=True, method="cs", out_stream=None)
            assert_check_partials(check, atol=1e-8, rtol=1e-8)

    def test_repeated_partials_symmetric_surface(self):
        """
        Given a MomentCoefficient with a single symmetric surface and random inputs,
        when compute_partials is called several times in a row,
        then the partials of M and CM match complex step every time, not only on the first call.
        """
        surfaces = [{"name": "wing", "mesh": np.zeros((2, 7)), "symmetry": True}]
        self.check_repeated_partials(surfaces)

    def test_repeated_partials_nonsymmetric_surface(self):
        """
        Given a MomentCoefficient with a single non-symmetric surface and random inputs,
        when compute_partials is called several times in a row,
        then the partials of M and CM match complex step every time, not only on the first call.
        """
        surfaces = [{"name": "wing", "mesh": np.zeros((3, 5)), "symmetry": False}]
        self.check_repeated_partials(surfaces)

    def test_repeated_partials_multiple_surfaces(self):
        """
        Given a MomentCoefficient with a symmetric wing and a non-symmetric tail and random inputs,
        when compute_partials is called several times in a row,
        then the partials of M and CM match complex step every time, including the cg partials summed over surfaces.
        """
        surfaces = [
            {"name": "wing", "mesh": np.zeros((2, 7)), "symmetry": True},
            {"name": "tail", "mesh": np.zeros((3, 5)), "symmetry": False},
        ]
        self.check_repeated_partials(surfaces)


if __name__ == "__main__":
    unittest.main()
