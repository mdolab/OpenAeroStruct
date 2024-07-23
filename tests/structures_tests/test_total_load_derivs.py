import unittest

from openaerostruct.structures.wing_weight_loads import StructureWeightLoads
from openaerostruct.structures.total_loads import TotalLoads
from openaerostruct.utils.testing import run_test, get_default_surfaces
import openmdao.api as om
import numpy as np


class Test(unittest.TestCase):
    def test_0(self):
        surface = get_default_surfaces()[0]

        comp = TotalLoads(surface=surface)

        run_test(self, comp, complex_flag=True)

    def test_1(self):
        surface = get_default_surfaces()[0]
        surface["struct_weight_relief"] = True

        comp = TotalLoads(surface=surface)

        run_test(self, comp, complex_flag=True)

    def test_2(self):
        surface = get_default_surfaces()[0]
        surface["distributed_fuel_weight"] = True

        comp = TotalLoads(surface=surface)

        run_test(self, comp, complex_flag=True)

    def test_structural_mass_loads(self):
        surface = get_default_surfaces()[0]

        comp = StructureWeightLoads(surface=surface)

        group = om.Group()

        indep_var_comp = om.IndepVarComp()

        ny = surface["mesh"].shape[1]

        # carefully chosen "random" values that give non-uniform derivatives outputs that are good for testing
        nodesval = np.array([[1.0, 2.0, 4.0], [20.0, 22.0, 7.0], [8.0, 17.0, 14.0], [13.0, 14.0, 16.0]], dtype=complex)
        element_mass_val = np.arange(ny - 1) + 1

        indep_var_comp.add_output("nodes", val=nodesval, units="m")
        indep_var_comp.add_output("element_mass", val=element_mass_val, units="kg")

        group.add_subsystem("indep_var_comp", indep_var_comp, promotes=["*"])
        group.add_subsystem("load", comp, promotes=["*"])

        run_test(self, group, complex_flag=True, compact_print=True)


if __name__ == "__main__":
    unittest.main()
