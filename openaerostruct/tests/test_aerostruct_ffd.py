from openmdao.utils.assert_utils import assert_rel_error
import unittest
import numpy as np
from openaerostruct.utils.constants import grav_constant

try:
    import pygeo
    pygeo_flag = True
except:
    pygeo_flag = False

@unittest.skipUnless(pygeo_flag, "pyGeo is required.")
class Test(unittest.TestCase):

    def test(self):
        from openaerostruct.geometry.utils import generate_mesh, write_FFD_file

        from openaerostruct.integration.aerostruct_groups import AerostructGeometry, AerostructPoint

        import openmdao.api as om
        from pygeo import DVGeometry


        # Create a dictionary to store options about the surface
        mesh_dict = {'num_y' : 5,
                     'num_x' : 2,
                     'wing_type' : 'CRM',
                     'symmetry' : True,
                     'num_twist_cp' : 5}

        mesh, twist_cp = generate_mesh(mesh_dict)

        surf_dict = {
                    # Wing definition
                    'name' : 'wing',        # name of the surface
                    'symmetry' : True,     # if true, model one half of wing
                                            # reflected across the plane y = 0
                    'S_ref_type' : 'wetted', # how we compute the wing area,
                                             # can be 'wetted' or 'projected'
                    'fem_model_type' : 'tube',

                    'thickness_cp' : np.array([.1, .2, .3]),

                    'mesh' : mesh,

                    'geom_manipulator' : 'FFD',
                    'mx' : 2,
                    'my' : 3,

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

        surfaces = [surf_dict]

        # Create the problem and assign the model group
        prob = om.Problem()

        # Add problem information as an independent variables component
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('v', val=248.136, units='m/s')
        indep_var_comp.add_output('alpha', val=5., units='deg')
        indep_var_comp.add_output('Mach_number', val=0.84)
        indep_var_comp.add_output('re', val=1.e6, units='1/m')
        indep_var_comp.add_output('rho', val=0.38, units='kg/m**3')
        indep_var_comp.add_output('CT', val=grav_constant * 17.e-6, units='1/s')
        indep_var_comp.add_output('R', val=11.165e6, units='m')
        indep_var_comp.add_output('W0', val=0.4 * 3e5,  units='kg')
        indep_var_comp.add_output('speed_of_sound', val=295.4, units='m/s')
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

            filename = write_FFD_file(surface, surface['mx'], surface['my'])
            DVGeo = DVGeometry(filename)
            aerostruct_group = AerostructGeometry(surface=surface, DVGeo=DVGeo)

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

                com_name = point_name + '.' + name + '_perf'
                prob.model.connect(name + '.local_stiff_transformed', point_name + '.coupled.' + name + '.local_stiff_transformed')
                prob.model.connect(name + '.nodes', point_name + '.coupled.' + name + '.nodes')

                # Connect aerodyamic mesh to coupled group mesh
                prob.model.connect(name + '.mesh', point_name + '.coupled.' + name + '.mesh')

                # Connect performance calculation variables
                prob.model.connect(name + '.radius', com_name + '.radius')
                prob.model.connect(name + '.thickness', com_name + '.thickness')
                prob.model.connect(name + '.nodes', com_name + '.nodes')
                prob.model.connect(name + '.cg_location', point_name + '.' + 'total_perf.' + name + '_cg_location')
                prob.model.connect(name + '.structural_mass', point_name + '.' + 'total_perf.' + name + '_structural_mass')
                prob.model.connect(name + '.t_over_c', com_name + '.t_over_c')

        # Import the Scipy Optimizer and set the driver of the problem to use
        # it, which defaults to an SLSQP optimization method
        prob.driver = om.ScipyOptimizeDriver()

        recorder = om.SqliteRecorder("aerostruct_ffd.db")
        prob.driver.add_recorder(recorder)
        prob.driver.recording_options['record_derivatives'] = True
        prob.driver.recording_options['includes'] = ['*']

        # Setup problem and add design variables, constraint, and objective
        prob.model.add_design_var('wing.shape', lower=-3, upper=2)
        prob.model.add_design_var('wing.thickness_cp', lower=0.01, upper=0.5, scaler=1e2)
        prob.model.add_constraint('AS_point_0.wing_perf.failure', upper=0.)
        prob.model.add_constraint('AS_point_0.wing_perf.thickness_intersects', upper=0.)

        # Add design variables, constraisnt, and objective on the problem
        prob.model.add_design_var('alpha', lower=-10., upper=10.)
        prob.model.add_constraint('AS_point_0.L_equals_W', equals=0.)
        prob.model.add_objective('AS_point_0.fuelburn', scaler=1e-5)

        # iprofile.setup()
        # iprofile.start()

        # Set up the problem
        prob.setup()

        # om.view_model(prob, outfile='aerostruct_ffd', show_browser=False)

        # prob.run_model()
        prob.run_driver()

        # prob.check_partials(compact_print=True)

        # print("\nWing CL:", prob['aero_point_0.wing_perf.CL'])
        # print("Wing CD:", prob['aero_point_0.wing_perf.CD'])


        # from helper import plot_3d_points
        #
        # mesh = prob['aero_point_0.wing.def_mesh']
        # plot_3d_points(mesh)
        #
        # filename = mesh_dict['wing_type'] + '_' + str(mesh_dict['num_x']) + '_' + str(mesh_dict['num_y'])
        # filename += '_' + str(surf_dict['mx']) + '_' + str(surf_dict['my']) + '.mesh'
        # np.save(filename, mesh)

        assert_rel_error(self, prob['AS_point_0.fuelburn'][0], 97680.8964568375, 1e-3)


if __name__ == '__main__':
    unittest.main()
