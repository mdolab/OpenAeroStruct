from openmdao.utils.assert_utils import assert_rel_error, assert_check_partials
import unittest
from openaerostruct.utils.constants import grav_constant


"""
This is not a physically meaningful case, in that it is not set up to give
the correct results or with the correct coupled physics between the two points.
Instead, this is simply testing the morphing aspect of the code using a
multipoint aerostruct formulation.
"""

class Test(unittest.TestCase):

    def test(self):
        import numpy as np

        from openaerostruct.geometry.utils import generate_mesh
        from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint
        from openaerostruct.utils.constants import grav_constant

        import openmdao.api as om

        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 11,
                     'num_x' : 2,
                     'wing_type' : 'CRM',
                     'symmetry' : True,
                     'num_twist_cp' : 5}

        mesh, twist_cp = generate_mesh(mesh_dict)

        surface = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted', # how we compute the wing area,
                                             # can be 'wetted' or 'projected'
                    'fem_model_type' : 'tube',

                    'thickness_cp' : np.array([.1, .2, .3]),

                    'twist_cp' : twist_cp,
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

                    # Structural values are based on aluminum 7075
                    'E' : 70.e9,            # [Pa] Young's modulus of the spar
                    'G' : 30.e9,            # [Pa] shear modulus of the spar
                    'yield' : 500.e6 / 2.5, # [Pa] yield stress divided by 2.5 for limiting case
                    'mrho' : 3.e3,          # [kg/m^3] material density
                    'fem_origin' : 0.35,    # normalized chordwise location of the spar
                    'wing_weight_ratio' : 2.,
                    'struct_weight_relief' : False,    # True to add the weight of the structure to the loads on the structure
                    'distributed_fuel_weight' : False,
                    # Constraints
                    'exact_failure_constraint' : False, # if false, use KS function
                    }

        # Create the problem and assign the model group
        prob = om.Problem()

        # Add problem information as an independent variables component
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('v', val=[248.136, 0.5 * 340.], units='m/s')
        indep_var_comp.add_output('alpha', val=[5., 10.], units='deg')
        indep_var_comp.add_output('Mach_number', val=[0.84, 0.5])
        indep_var_comp.add_output('re', val=[1.e6, 0.5e6], units='1/m')
        indep_var_comp.add_output('rho', val=[0.38, .764], units='kg/m**3')
        indep_var_comp.add_output('CT', val=grav_constant * 17.e-6, units='1/s')
        indep_var_comp.add_output('R', val=11.165e6, units='m')
        indep_var_comp.add_output('W0', val=0.4 * 3e5,  units='kg')
        indep_var_comp.add_output('speed_of_sound', val=295.4, units='m/s')
        indep_var_comp.add_output('load_factor', val=[1., 2.5])
        indep_var_comp.add_output('empty_cg', val=np.zeros((3)), units='m')

        prob.model.add_subsystem('prob_vars',
             indep_var_comp,
             promotes=['*'])

        # Add morphing variables as an independent variables component
        morphing_vars = om.IndepVarComp()
        morphing_vars.add_output('t_over_c_cp', val=np.array([0.15]))
        morphing_vars.add_output('thickness_cp', val=np.array([0.01, 0.01, 0.01]), units='m')
        morphing_vars.add_output('twist_cp_0', val=np.array([2., 3., 4., 4., 4.]), units='deg')
        morphing_vars.add_output('twist_cp_1', val=np.array([4., 4., 4., 5., 6.]), units='deg')

        prob.model.add_subsystem('morphing_vars',
             morphing_vars,
             promotes=['*'])

        # Connect geometric design variables to each point
        prob.model.connect('t_over_c_cp', 'AS_point_0.wing.geometry.t_over_c_cp')
        prob.model.connect('t_over_c_cp', 'AS_point_1.wing.geometry.t_over_c_cp')

        prob.model.connect('thickness_cp', 'AS_point_0.wing.tube_group.thickness_cp')
        prob.model.connect('thickness_cp', 'AS_point_1.wing.tube_group.thickness_cp')

        prob.model.connect('twist_cp_0', 'AS_point_0.wing.geometry.twist_cp')
        prob.model.connect('twist_cp_1', 'AS_point_1.wing.geometry.twist_cp')

        for point in range(2):

            name = 'wing'

            point_name = 'AS_point_{}'.format(point)

            # Create the aero point group and add it to the model
            AS_point = AerostructPoint(surfaces=[surface])

            prob.model.add_subsystem(point_name, AS_point)

            aerostruct_group = AerostructGeometry(surface=surface, connect_geom_DVs=False)
            AS_point.add_subsystem(name, aerostruct_group)

            # Connect flow properties to the analysis point
            prob.model.connect('alpha', point_name + '.alpha', src_indices=[point])
            prob.model.connect('v', point_name + '.v', src_indices=[point])
            prob.model.connect('Mach_number', point_name + '.Mach_number', src_indices=[point])
            prob.model.connect('re', point_name + '.re', src_indices=[point])
            prob.model.connect('rho', point_name + '.rho', src_indices=[point])
            prob.model.connect('CT', point_name + '.CT')
            prob.model.connect('R', point_name + '.R')
            prob.model.connect('W0', point_name + '.W0')
            prob.model.connect('speed_of_sound', point_name + '.speed_of_sound')
            prob.model.connect('empty_cg', point_name + '.empty_cg')
            prob.model.connect('load_factor', point_name + '.load_factor', src_indices=[point])

            com_name = point_name + '.' + name + '_perf'
            AS_point.connect(name + '.local_stiff_transformed', 'coupled.' + name + '.local_stiff_transformed')
            AS_point.connect(name + '.nodes', 'coupled.' + name + '.nodes')

            # Connect aerodyamic mesh to coupled group mesh
            AS_point.connect(name + '.mesh', 'coupled.' + name + '.mesh')

            # Connect performance calculation variables
            AS_point.connect(name + '.radius', name + '_perf' + '.radius')
            AS_point.connect(name + '.thickness', name + '_perf' + '.thickness')
            AS_point.connect(name + '.nodes', name + '_perf' + '.nodes')
            AS_point.connect(name + '.cg_location', 'total_perf.' + name + '_cg_location')
            AS_point.connect(name + '.structural_mass', 'total_perf.' + name + '_structural_mass')
            AS_point.connect(name + '.geometry.t_over_c', name + '_perf' + '.t_over_c')
            AS_point.connect(name + '.geometry.t_over_c', name + '.t_over_c')

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['tol'] = 1e-9
        prob.driver.options['maxiter'] = 2

        recorder = om.SqliteRecorder("morphing_aerostruct.db")
        prob.driver.add_recorder(recorder)
        prob.driver.recording_options['record_derivatives'] = True
        prob.driver.recording_options['includes'] = ['*']

        # Setup problem and add design variables, constraint, and objective
        prob.model.add_design_var('twist_cp_0', lower=-10., upper=15.)
        prob.model.add_design_var('twist_cp_1', lower=-10., upper=15.)
        prob.model.add_design_var('thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
        prob.model.add_constraint('AS_point_0.wing_perf.failure', upper=0.)
        prob.model.add_constraint('AS_point_0.wing_perf.thickness_intersects', upper=0.)
        prob.model.add_constraint('AS_point_1.wing_perf.failure', upper=0.)
        prob.model.add_constraint('AS_point_1.wing_perf.thickness_intersects', upper=0.)

        # Add design variables, constraisnt, and objective on the problem
        prob.model.add_design_var('alpha', lower=-15., upper=15.)
        prob.model.add_constraint('AS_point_0.L_equals_W', equals=0.)
        prob.model.add_constraint('AS_point_1.L_equals_W', equals=0.)
        prob.model.add_objective('AS_point_0.fuelburn', scaler=1e-5)

        # Set up the problem
        prob.setup(check=True)

        # om.view_model(prob)

        prob.run_model()

        # Check the partials at the initial point in the design space,
        # only care about relative error
        data = prob.check_partials(compact_print=True, out_stream=None, method='cs', step=1e-40)

        assert_check_partials(data, atol=1e20, rtol=1e-6)

        # Run the optimizer for 2 iterations
        prob.run_driver()

        # Check the partials at this point in the design space
        data = prob.check_partials(compact_print=True, out_stream=None, method='cs', step=1e-40)
        assert_check_partials(data, atol=1e20, rtol=1e-6)


if __name__ == '__main__':
    unittest.main()
