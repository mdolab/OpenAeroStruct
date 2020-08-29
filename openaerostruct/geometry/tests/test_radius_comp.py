from numpy.testing import assert_almost_equal, assert_equal

import unittest

import openmdao.api as om

from openaerostruct.geometry.radius_comp import RadiusComp
from openaerostruct.utils.testing import run_test, get_default_surfaces

import numpy as np

class Test(unittest.TestCase):

    def test(self):
        surfaces = get_default_surfaces()

        group = om.Group()

        comp = RadiusComp(surface=surfaces[0])
        ny = surfaces[0]['mesh'].shape[1]

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('mesh', val=surfaces[0]['mesh'], units='m')
        indep_var_comp.add_output('t_over_c', val=np.linspace(0.1,0.5,num=ny-1))

        group.add_subsystem('radius', comp)
        group.add_subsystem('indep_var_comp', indep_var_comp)

        group.connect('indep_var_comp.mesh', 'radius.mesh')
        group.connect('indep_var_comp.t_over_c', 'radius.t_over_c')

        run_test(self, group)


if __name__ == '__main__':
    unittest.main()
