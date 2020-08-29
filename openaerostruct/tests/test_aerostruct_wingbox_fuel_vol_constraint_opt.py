from openmdao.utils.assert_utils import assert_rel_error
import unittest
import numpy as np

from openaerostruct.geometry.utils import generate_mesh

from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

import openmdao.api as om
from openaerostruct.structures.wingbox_fuel_vol_delta import WingboxFuelVolDelta

# Provide coordinates for a portion of an airfoil for the wingbox cross-section as an nparray with dtype=complex (to work with the complex-step approximation for derivatives).
# These should be for an airfoil with the chord scaled to 1.
# We use the 10% to 60% portion of the NASA SC2-0612 airfoil for this case
# We use the coordinates available from airfoiltools.com. Using such a large number of coordinates is not necessary.
# The first and last x-coordinates of the upper and lower skins must be the same

upper_x = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6], dtype = 'complex128')
lower_x = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6], dtype = 'complex128')
upper_y = np.array([ 0.0447,  0.046,  0.0472,  0.0484,  0.0495,  0.0505,  0.0514,  0.0523,  0.0531,  0.0538, 0.0545,  0.0551,  0.0557, 0.0563,  0.0568, 0.0573,  0.0577,  0.0581,  0.0585,  0.0588,  0.0591,  0.0593,  0.0595,  0.0597,  0.0599,  0.06,    0.0601,  0.0602,  0.0602,  0.0602,  0.0602,  0.0602,  0.0601,  0.06,    0.0599,  0.0598,  0.0596,  0.0594,  0.0592,  0.0589,  0.0586,  0.0583,  0.058,   0.0576,  0.0572,  0.0568,  0.0563,  0.0558,  0.0553,  0.0547,  0.0541], dtype = 'complex128')
lower_y = np.array([-0.0447, -0.046, -0.0473, -0.0485, -0.0496, -0.0506, -0.0515, -0.0524, -0.0532, -0.054, -0.0547, -0.0554, -0.056, -0.0565, -0.057, -0.0575, -0.0579, -0.0583, -0.0586, -0.0589, -0.0592, -0.0594, -0.0595, -0.0596, -0.0597, -0.0598, -0.0598, -0.0598, -0.0598, -0.0597, -0.0596, -0.0594, -0.0592, -0.0589, -0.0586, -0.0582, -0.0578, -0.0573, -0.0567, -0.0561, -0.0554, -0.0546, -0.0538, -0.0529, -0.0519, -0.0509, -0.0497, -0.0485, -0.0472, -0.0458, -0.0444], dtype = 'complex128')

class Test(unittest.TestCase):

    def test(self):
        """
        This is an opt problem that tests the fuel volume constraint with the wingbox model
        """

        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 7,
                     'num_x' : 2,
                     'wing_type' : 'CRM',
                     'symmetry' : True,
                     'num_twist_cp' : 6,
                     'chord_cos_spacing' : 0,
                     'span_cos_spacing' : 0,
                     }

        mesh, twist_cp = generate_mesh(mesh_dict)

        surf_dict = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted', # how we compute the wing area,
                                             # can be 'wetted' or 'projected'
                    'fem_model_type' : 'wingbox',

                    'spar_thickness_cp' : np.array([0.004, 0.005, 0.005, 0.008, 0.008, 0.01]), # [m]
                    'skin_thickness_cp' : np.array([0.005, 0.01, 0.015, 0.020, 0.025, 0.026]),
                    'twist_cp' : np.array([4., 5., 8., 8., 8., 9.]),
                    'mesh' : mesh,

                    'data_x_upper' : upper_x,
                    'data_x_lower' : lower_x,
                    'data_y_upper' : upper_y,
                    'data_y_lower' : lower_y,
                    'strength_factor_for_upper_skin' : 1.,

                    # Aerodynamic performance of the lifting surface at
                    # an angle of attack of 0 (alpha=0).
                    # These CL0 and CD0 values are added to the CL and CD
                    # obtained from aerodynamic analysis of the surface to get
                    # the total CL and CD.
                    # These CL0 and CD0 values do not vary wrt alpha.
                    'CL0' : 0.0,            # CL of the surface at alpha=0
                    'CD0' : 0.0078,            # CD of the surface at alpha=0

                    # Airfoil properties for viscous drag calculation
                    'k_lam' : 0.05,         # percentage of chord with laminar
                                            # flow, used for viscous drag
                    't_over_c_cp' : np.array([0.08, 0.08, 0.08, 0.10, 0.10, 0.08]),
                    'original_wingbox_airfoil_t_over_c' : 0.12,
                    'c_max_t' : .38,       # chordwise location of maximum thickness
                    'with_viscous' : True,
                    'with_wave' : False,     # if true, compute wave drag

                    # Structural values are based on aluminum 7075
                    'E' : 73.1e9,              # [Pa] Young's modulus
                    'G' : (73.1e9/2/1.33),     # [Pa] shear modulus (calculated using E and the Poisson's ratio here)
                    'yield' : (420.e6 / 1.5),  # [Pa] allowable yield stress
                    'mrho' : 2.78e3,           # [kg/m^3] material density
                    'strength_factor_for_upper_skin' : 1.0, # the yield stress is multiplied by this factor for the upper skin
                    # 'fem_origin' : 0.35,    # normalized chordwise location of the spar
                    'wing_weight_ratio' : 1.25,
                    'struct_weight_relief' : True,
                    'distributed_fuel_weight' : False,
                    # Constraints
                    'exact_failure_constraint' : False, # if false, use KS function
                    'fuel_density' : 803.,      # [kg/m^3] fuel density (only needed if the fuel-in-wing volume constraint is used)
                    'Wf_reserve' :15000.,       # [kg] reserve fuel mass
                    }

        surfaces = [surf_dict]

        # Create the problem and assign the model group
        prob = om.Problem()

        # Add problem information as an independent variables component
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('v', val=.85 * 295.07, units='m/s')
        indep_var_comp.add_output('alpha', val=0., units='deg')
        indep_var_comp.add_output('Mach_number', val=0.85)
        indep_var_comp.add_output('re', val=0.348*295.07*.85*1./(1.43*1e-5), units='1/m')
        indep_var_comp.add_output('rho', val=0.348, units='kg/m**3')
        indep_var_comp.add_output('CT', val=0.53/3600, units='1/s')
        indep_var_comp.add_output('R', val=14.307e6, units='m')
        indep_var_comp.add_output('W0', val=148000 + surf_dict['Wf_reserve'],  units='kg')
        indep_var_comp.add_output('speed_of_sound', val=295.07, units='m/s')
        indep_var_comp.add_output('load_factor', val=1.)
        indep_var_comp.add_output('empty_cg', val=np.zeros((3)), units='m')

        prob.model.add_subsystem('prob_vars',
             indep_var_comp,
             promotes=['*'])

        # Loop over each surface in the surfaces list
        for surface in surfaces:

            # Get the surface name and create a group to contain components
            # only for this surface
            name = surface['name']

            aerostruct_group = AerostructGeometry(surface=surface)

            # Add tmp_group to the problem with the name of the surface.
            prob.model.add_subsystem(name, aerostruct_group)

        # Loop through and add a certain number of aero points
        for i in range(1):

            point_name = 'AS_point_{}'.format(i)
            # Connect the parameters within the model for each aero point

            # Create the aero point group and add it to the model
            AS_point = AerostructPoint(surfaces=surfaces)

            prob.model.add_subsystem(point_name, AS_point)

            # Connect flow properties to the analysis point
            prob.model.connect('v', point_name + '.v')
            prob.model.connect('alpha', point_name + '.alpha')
            prob.model.connect('Mach_number', point_name + '.Mach_number')
            prob.model.connect('re', point_name + '.re')
            prob.model.connect('rho', point_name + '.rho')
            prob.model.connect('CT', point_name + '.CT')
            prob.model.connect('R', point_name + '.R')
            prob.model.connect('W0', point_name + '.W0')
            prob.model.connect('speed_of_sound', point_name + '.speed_of_sound')
            prob.model.connect('empty_cg', point_name + '.empty_cg')
            prob.model.connect('load_factor', point_name + '.load_factor')

            for surface in surfaces:

                com_name = point_name + '.' + name + '_perf.'
                prob.model.connect(name + '.local_stiff_transformed', point_name + '.coupled.' + name + '.local_stiff_transformed')
                prob.model.connect(name + '.nodes', point_name + '.coupled.' + name + '.nodes')

                # Connect aerodyamic mesh to coupled group mesh
                prob.model.connect(name + '.mesh', point_name + '.coupled.' + name + '.mesh')
                prob.model.connect(name + '.element_mass', point_name + '.coupled.' + name + '.element_mass')

                # Connect performance calculation variables
                prob.model.connect(name + '.nodes', com_name + 'nodes')
                prob.model.connect(name + '.cg_location', point_name + '.' + 'total_perf.' + name + '_cg_location')
                prob.model.connect(name + '.structural_mass', point_name + '.' + 'total_perf.' + name + '_structural_mass')

                # Connect wingbox properties to von Mises stress calcs
                prob.model.connect(name + '.Qz', com_name + 'Qz')
                prob.model.connect(name + '.J', com_name + 'J')
                prob.model.connect(name + '.A_enc', com_name + 'A_enc')
                prob.model.connect(name + '.htop', com_name + 'htop')
                prob.model.connect(name + '.hbottom', com_name + 'hbottom')
                prob.model.connect(name + '.hfront', com_name + 'hfront')
                prob.model.connect(name + '.hrear', com_name + 'hrear')

                prob.model.connect(name + '.spar_thickness', com_name + 'spar_thickness')
                prob.model.connect(name + '.t_over_c', com_name + 't_over_c')

            #=======================================================================================
            # Here we add the fuel volume constraint componenet to the model
            #=======================================================================================
            prob.model.add_subsystem('fuel_vol_delta', WingboxFuelVolDelta(surface=surface))
            prob.model.connect('AS_point_0.fuelburn', 'fuel_vol_delta.fuelburn')
            prob.model.connect('wing.struct_setup.fuel_vols', 'fuel_vol_delta.fuel_vols')
            #=======================================================================================
            #=======================================================================================

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['tol'] = 1e-9

        # prob.driver = om.pyOptSparseDriver()
        # # prob.driver.add_recorder(om.SqliteRecorder("cases.sql"))
        # prob.driver.options['optimizer'] = "SNOPT"
        # prob.driver.opt_settings['Major optimality tolerance'] = 1e-6
        # prob.driver.opt_settings['Major feasibility tolerance'] = 1e-8
        # prob.driver.opt_settings['Major iterations limit'] = 200

        prob.model.add_objective('AS_point_0.fuelburn', scaler=1e-5)

        prob.model.add_design_var('wing.twist_cp', lower=-15., upper=15., scaler=0.1)
        prob.model.add_design_var('wing.spar_thickness_cp', lower=0.003, upper=0.1, scaler=1e2)
        prob.model.add_design_var('wing.skin_thickness_cp', lower=0.003, upper=0.1, scaler=1e2)
        prob.model.add_design_var('wing.geometry.t_over_c_cp', lower=0.07, upper=0.2, scaler=10.)

        prob.model.add_constraint('AS_point_0.CL', equals=0.5)
        prob.model.add_constraint('AS_point_0.wing_perf.failure', upper=0.)

        #=======================================================================================
        # Here we add the fuel volume constraint
        #=======================================================================================
        prob.model.add_constraint('fuel_vol_delta.fuel_vol_delta', lower=0.)
        #=======================================================================================
        #=======================================================================================

        # Set up the problem
        prob.setup()

        # om.view_model(prob)

        prob.run_driver()

        # prob.check_partials(form='central', compact_print=True)

        print(prob['AS_point_0.fuelburn'][0])
        print(prob['wing.structural_mass'][0]/1.25)
        print(prob['fuel_vol_delta.fuel_vol_delta'][0])

        assert_rel_error(self, prob['AS_point_0.fuelburn'][0], 82019.934119, 1e-5)
        assert_rel_error(self, prob['wing.structural_mass'][0]/1.25, 12211.382514, 1e-4)
        assert_rel_error(self, prob['fuel_vol_delta.fuel_vol_delta'][0], 35.3020703438527, 1e-4)


if __name__ == '__main__':
    unittest.main()
