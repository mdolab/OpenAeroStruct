from openmdao.utils.assert_utils import assert_rel_error
import unittest
import numpy as np

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.transfer.displacement_transfer import DisplacementTransfer
from openaerostruct.structures.struct_groups import SpatialBeamAlone

import openmdao.api as om


class Test(unittest.TestCase):

    def test(self):

        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 7,
                     'wing_type' : 'CRM',
                     'symmetry' : True,
                     'num_twist_cp' : 5}

        mesh, twist_cp = generate_mesh(mesh_dict)

        surf_dict = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'fem_model_type' : 'tube',

                    'mesh' : mesh,

                    # Structural values are based on aluminum 7075
                    'E' : 70.e9,            # [Pa] Young's modulus of the spar
                    'G' : 30.e9,            # [Pa] shear modulus of the spar
                    'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
                    'mrho' : 3.e3,          # [kg/m^3] material density
                    'fem_origin' : 0.35,    # normalized chordwise location of the spar
                    't_over_c_cp' : np.array([0.15]),      # maximum airfoil thickness
                    'thickness_cp' : np.ones((3)) * .1,
                    'wing_weight_ratio' : 2.,
                    'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                    'distributed_fuel_weight' : False,
                    'exact_failure_constraint' : False,
                    }

        # Create the problem and assign the model group
        prob = om.Problem()

        ny = surf_dict['mesh'].shape[1]

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('loads', val=np.zeros((ny, 6)), units='N')
        indep_var_comp.add_output('load_factor', val=1.)

        struct_group = SpatialBeamAlone(surface=surf_dict)

        # Add indep_vars to the structural group
        struct_group.add_subsystem('indep_vars',
             indep_var_comp,
             promotes=['*'])

        prob.model.add_subsystem(surf_dict['name'], struct_group, promotes=['*'])

        prob.driver = om.ScipyOptimizeDriver()

        # Setup problem and add design variables, constraint, and objective
        prob.model.add_design_var('thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
        prob.model.add_constraint('failure', upper=0.)
        prob.model.add_constraint('thickness_intersects', upper=0.)

        # Add design variables, constraisnt, and objective on the problem
        prob.model.add_objective('structural_mass', scaler=1e-5)

        # Set up the problem
        prob.setup()

        prob.run_driver()

        assert_rel_error(self, prob['structural_mass'][0], 13601.162582, 1e-4)
        assert_rel_error(self, prob['disp'][1, 2], 0., 1e-4)



if __name__ == '__main__':
    unittest.main()
