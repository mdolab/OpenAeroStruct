import unittest
import numpy as np

from openmdao.utils.assert_utils import assert_near_equal
import openmdao.api as om
from openaerostruct.structures.compute_point_mass_loads import ComputePointMassLoads
from openaerostruct.utils.testing import run_test, get_default_surfaces


derivs_added = False


class Test(unittest.TestCase):
    @unittest.skipUnless(derivs_added, "Analytic derivs not added yet")
    def test_derivs(self):
        surface = get_default_surfaces()[0]

        surface["n_point_masses"] = 2

        comp = ComputePointMassLoads(surface=surface)

        group = om.Group()

        indep_var_comp = om.IndepVarComp()

        nodesval = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])

        point_masses = np.array([[2.0, 1.0]])

        point_mass_locations = np.array([[2.1, 0.1, 0.2], [3.2, 1.2, 0.3]])

        indep_var_comp.add_output("nodes", val=nodesval, units="m")
        indep_var_comp.add_output("point_masses", val=point_masses, units="kg")
        indep_var_comp.add_output("point_mass_locations", val=point_mass_locations, units="m")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("compute_point_mass_loads", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, step=1e-8, atol=1e-5, compact_print=True)

    @unittest.skipUnless(derivs_added, "Analytic derivs not added yet")
    def test_simple_values(self):
        surface = get_default_surfaces()[0]

        surface["n_point_masses"] = 1

        comp = ComputePointMassLoads(surface=surface)

        group = om.Group()

        indep_var_comp = om.IndepVarComp()

        nodesval = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 2.0, 0.0], [0.0, 3.0, 0.0]])

        point_masses = np.array([[1 / 9.8]])

        point_mass_locations = np.array([[0.55012, 0.1, 0.0]])

        indep_var_comp.add_output("nodes", val=nodesval, units="m")
        indep_var_comp.add_output("point_masses", val=point_masses, units="kg")
        indep_var_comp.add_output("point_mass_locations", val=point_mass_locations, units="m")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("compute_point_mass_loads", comp, promotes=["*"])

        prob = run_test(self, group, complex_flag=True, step=1e-8, atol=1e-5, compact_print=True)

        truth_array = np.array([0, 0, -1.0, 0.0, 0.55012, 0.0])

        assert_near_equal(prob["comp.loads_from_point_masses"][0, :], truth_array, 1e-6)


if __name__ == "__main__":
    unittest.main()
