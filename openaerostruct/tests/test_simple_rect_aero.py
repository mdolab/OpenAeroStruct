import unittest

import numpy as np

import openmdao.api as om
from openmdao.utils.assert_utils import assert_rel_error

from openaerostruct.geometry.utils import generate_mesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.aerodynamics.aero_groups import AeroPoint


class Test(unittest.TestCase):

    def test(self):

        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 5,
                     'num_x' : 2,
                     'wing_type' : 'rect',
                     'symmetry' : True}

        mesh = generate_mesh(mesh_dict)

        surf_dict = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'type' : 'aero',
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted', # how we compute the wing area,
                                             # can be 'wetted' or 'projected'

                    'twist_cp' : np.array([0.]),

                    'mesh' : mesh,

                    # Aerodynamic performance of the lifting surface at
                    # an angle of attack of 0 (alpha=0).
                    # These CL0 and CD0 values are added to the CL and CD
                    # obtained from aerodynamic analysis of the surface to get
                    # the total CL and CD.
                    # These CL0 and CD0 values do not vary wrt alpha.
                    'CL0' : 0.0,            # CL of the surface at alpha=0
                    'CD0' : 0.015,            # CD of the surface at alpha=0

                    # Airfoil properties for viscous drag calculation
                    'k_lam' : 0.05,         # percentage of chord with laminar
                                            # flow, used for viscous drag
                    't_over_c_cp' : np.array([0.15]),      # thickness over chord ratio (NACA0015)
                    'c_max_t' : .303,       # chordwise location of maximum (NACA0015)
                                            # thickness
                    'with_viscous' : True,
                    'with_wave' : False,     # if true, compute wave drag
                    }

        surfaces = [surf_dict]

        # Create the problem and the model group
        prob = om.Problem()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('v', val=248.136, units='m/s')
        indep_var_comp.add_output('alpha', val=5., units='deg')
        indep_var_comp.add_output('Mach_number', val=0.84)
        indep_var_comp.add_output('re', val=1.e6, units='1/m')
        indep_var_comp.add_output('rho', val=0.38, units='kg/m**3')
        indep_var_comp.add_output('cg', val=np.zeros((3)), units='m')

        prob.model.add_subsystem('prob_vars',
            indep_var_comp,
            promotes=['*'])

        # Loop over each surface in the surfaces list
        for surface in surfaces:

            geom_group = Geometry(surface=surface)

            # Add tmp_group to the problem as the name of the surface.
            # Note that is a group and performance group for each
            # individual surface.
            prob.model.add_subsystem(surface['name'], geom_group)

        # Loop through and add a certain number of aero points
        for i in range(1):

            # Create the aero point group and add it to the model
            aero_group = AeroPoint(surfaces=surfaces)
            point_name = 'aero_point_{}'.format(i)
            prob.model.add_subsystem(point_name, aero_group)

            # Connect flow properties to the analysis point
            prob.model.connect('v', point_name + '.v')
            prob.model.connect('alpha', point_name + '.alpha')
            prob.model.connect('Mach_number', point_name + '.Mach_number')
            prob.model.connect('re', point_name + '.re')
            prob.model.connect('rho', point_name + '.rho')
            prob.model.connect('cg', point_name + '.cg')

            # Connect the parameters within the model for each aero point
            for surface in surfaces:

                name = surface['name']

                # Connect the mesh from the geometry component to the analysis point
                prob.model.connect(name + '.mesh', point_name + '.' + name + '.def_mesh')

                # Perform the connections with the modified names within the
                # 'aero_states' group.
                prob.model.connect(name + '.mesh', point_name + '.aero_states.' + name + '_def_mesh')

                prob.model.connect(name + '.t_over_c', point_name + '.' + name + '_perf.' + 't_over_c')

        # Set up the problem
        prob.setup()

        prob.run_model()

        assert_rel_error(self, prob['aero_point_0.wing_perf.CD'][0], 0.03487336411850356, 1e-6)
        assert_rel_error(self, prob['aero_point_0.wing_perf.CL'][0], 0.4615561217697067, 1e-6)
        assert_rel_error(self, prob['aero_point_0.CM'][0], 0.0, 1e-6)
        assert_rel_error(self, prob['aero_point_0.CM'][1], -0.11507021674483686, 1e-6)
        assert_rel_error(self, prob['aero_point_0.CM'][2], 0.0, 1e-6)




if __name__ == '__main__':
    unittest.main()
