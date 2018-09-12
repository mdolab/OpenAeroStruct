import unittest
import numpy as np

from openmdao.api import Group, IndepVarComp
from openaerostruct.transfer.load_transfer import LoadTransfer 
from openaerostruct.utils.testing import run_test, get_default_surfaces


class Test(unittest.TestCase):

    def test(self):
        surface = get_default_surfaces()[0]
        ##surface['num_x'] = 7
        #surface['num_y'] = 2
        group = Group()
        
        comp = LoadTransfer(surface=surface)

        indep_var_comp = IndepVarComp()

        ny = surface['num_y']
        nx = surface['num_x']
        
        
        indep_var_comp.add_output('def_mesh', val=np.random.random((nx, ny, 3)), units='m')
        indep_var_comp.add_output('sec_forces', val=np.random.random((nx-1, ny-1, 3)), units='N')
        
        #indep_var_comp.add_output('loadA', val=np.random.random((ny,3)), units='N')
        #indep_var_comp.add_output('loadB', val=np.random.random((ny,3)), units='m')
        indep_var_comp.add_output('load', val=np.random.random((ny,6)), units='N')
        
        group.add_subsystem('indep_var_comp', indep_var_comp, promotes=['*'])
        group.add_subsystem('load_transfer', comp, promotes=['*'])
        
        run_test(self, group, complex_flag=True, compact_print=False)

if __name__ == '__main__':
    unittest.main()
