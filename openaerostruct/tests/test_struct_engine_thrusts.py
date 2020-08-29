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
        mesh_dict = {'num_y' : 31,
                     'wing_type' : 'rect',
                     'span' : 10,
                     'symmetry' : True}

        mesh = generate_mesh(mesh_dict)

        surface = {
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
                    'fem_origin' : 0.5,    # normalized chordwise location of the spar
                    't_over_c_cp' : np.array([0.15]),      # maximum airfoil thickness
                    'thickness_cp' : np.ones((3)) * .1,
                    'wing_weight_ratio' : 2.,
                    'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                    'distributed_fuel_weight' : False,
                    'exact_failure_constraint' : False,
                    }

        # Create the problem and assign the model group
        prob = om.Problem()

        ny = surface['mesh'].shape[1]
        surface['n_point_masses'] = 1

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('loads', val=np.zeros((ny, 6)), units='N')
        indep_var_comp.add_output('load_factor', val=1.)

        engine_thrusts = np.array([[10.]])

        point_mass_locations = np.array([[1., -10., 0.]])

        indep_var_comp.add_output('engine_thrusts', val=engine_thrusts, units='N')
        indep_var_comp.add_output('point_mass_locations', val=point_mass_locations, units='m')

        struct_group = SpatialBeamAlone(surface=surface)

        # Add indep_vars to the structural group
        struct_group.add_subsystem('indep_vars',
             indep_var_comp,
             promotes=['*'])

        prob.model.add_subsystem(surface['name'], struct_group, promotes=['*'])

        # Set up the problem
        prob.setup()

        prob.run_model()

        assert_rel_error(self, prob['vonmises'][-1, 0], 300352.285697, 1e-4)

    def test_multiple_masses(self):

        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 31,
                     'wing_type' : 'rect',
                     'span' : 10,
                     'symmetry' : True}

        mesh = generate_mesh(mesh_dict)

        surface = {
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
                    'fem_origin' : 0.5,    # normalized chordwise location of the spar
                    't_over_c_cp' : np.array([0.15]),      # maximum airfoil thickness
                    'thickness_cp' : np.ones((3)) * .1,
                    'wing_weight_ratio' : 2.,
                    'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                    'distributed_fuel_weight' : False,
                    'exact_failure_constraint' : False,
                    }

        # Create the problem and assign the model group
        prob = om.Problem()

        ny = surface['mesh'].shape[1]
        surface['n_point_masses'] = 2

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('loads', val=np.zeros((ny, 6)), units='N')
        indep_var_comp.add_output('load_factor', val=1.)

        engine_thrusts = np.array([[10., 20.]])

        point_mass_locations = np.array([[1., -1., 0.],
                                         [1., -2., 0.]])

        indep_var_comp.add_output('engine_thrusts', val=engine_thrusts, units='N')
        indep_var_comp.add_output('point_mass_locations', val=point_mass_locations, units='m')

        struct_group = SpatialBeamAlone(surface=surface)

        # Add indep_vars to the structural group
        struct_group.add_subsystem('indep_vars',
             indep_var_comp,
             promotes=['*'])

        prob.model.add_subsystem(surface['name'], struct_group, promotes=['*'])

        # Set up the problem
        prob.setup()

        prob.run_model()

        assert_rel_error(self, prob['vonmises'][-1, 0], 137509.870021, 1e-4)


if __name__ == '__main__':
    unittest.main()
