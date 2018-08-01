import unittest
import numpy as np

from openmdao.api import Group, IndepVarComp
from openaerostruct.aerodynamics.geometry import VLMGeometry
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):

    def test(self):
        surfaces = get_default_surfaces()

        group = Group()

        comp = VLMGeometry(surface=surfaces[0])

        indep_var_comp = IndepVarComp()

        indep_var_comp.add_output('def_mesh', val=surfaces[0]['mesh'], units='m')

        group.add_subsystem('geom', comp)
        group.add_subsystem('indep_var_comp', indep_var_comp)

        group.connect('indep_var_comp.def_mesh', 'geom.def_mesh')

        run_test(self, group)

if __name__ == '__main__':
    unittest.main()
