from openmdao.utils.assert_utils import assert_rel_error
import unittest


class Test(unittest.TestCase):
    def test(self):
        import numpy as np
        import openmdao.api as om

        from openaerostruct.geometry.utils import generate_mesh
        from openaerostruct.geometry.geometry_group import Geometry
        from openaerostruct.aerodynamics.aero_groups import AeroPoint
        from openaerostruct.integration.multipoint_comps import MultiCD

        # Create a dictionary to store options about the surface
        mesh_dict = {
            'num_y': 5,
            'num_x': 3,
            'wing_type': 'CRM',
            'symmetry': True,
            'num_twist_cp': 5,
            'span_cos_spacing': 0.0,
        }

        mesh, twist_cp = generate_mesh(mesh_dict)

        surf_dict = {
            # Wing definition
            'name': 'wing',  # name of the surface
            'symmetry': True,  # if true, model one half of wing
            # reflected across the plane y = 0
            'S_ref_type': 'wetted',  # how we compute the wing area,
            # can be 'wetted' or 'projected'
            'fem_model_type': 'tube',
            'mesh': mesh,
            'twist_cp': twist_cp,
            # Aerodynamic performance of the lifting surface at
            # an angle of attack of 0 (alpha=0).
            # These CL0 and CD0 values are added to the CL and CD
            # obtained from aerodynamic analysis of the surface to get
            # the total CL and CD.
            # These CL0 and CD0 values do not vary wrt alpha.
            'CL0': 0.0,  # CL of the surface at alpha=0
            'CD0': 0.015,  # CD of the surface at alpha=0
            # Airfoil properties for viscous drag calculation
            'k_lam': 0.05,  # percentage of chord with laminar
            # flow, used for viscous drag
            't_over_c_cp': np.array([0.15]),  # thickness over chord ratio (NACA0015)
            'c_max_t': 0.303,  # chordwise location of maximum (NACA0015)
            # thickness
            'with_viscous': True,  # if true, compute viscous drag
            'with_wave': False,  # if true, compute wave drag
        }

        surfaces = [surf_dict]

        n_points = 2

        # Create the problem and the model group
        prob = om.Problem()

        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('v', val=248.136, units='m/s')
        indep_var_comp.add_output('alpha', val=np.ones(n_points) * 6.64, units='deg')
        indep_var_comp.add_output('Mach_number', val=0.84)
        indep_var_comp.add_output('re', val=1.0e6, units='1/m')
        indep_var_comp.add_output('rho', val=0.38, units='kg/m**3')
        indep_var_comp.add_output('cg', val=np.zeros((3)), units='m')

        prob.model.add_subsystem('prob_vars', indep_var_comp, promotes=['*'])

        # Loop over each surface and create the geometry groups
        for surface in surfaces:
            # Get the surface name and create a group to contain components only for this surface.
            name = surface['name']
            geom_group = Geometry(surface=surface)

            # Add geom_group to the problem with the name of the surface.
            prob.model.add_subsystem(name + '_geom', geom_group)

        # Loop through and add a certain number of aero points
        for i in range(n_points):
            # Create the aero point group and add it to the model
            aero_group = AeroPoint(surfaces=surfaces)
            point_name = 'aero_point_{}'.format(i)
            prob.model.add_subsystem(point_name, aero_group)

            # Connect flow properties to the analysis point
            prob.model.connect('v', point_name + '.v')
            prob.model.connect('alpha', point_name + '.alpha', src_indices=[i])
            prob.model.connect('Mach_number', point_name + '.Mach_number')
            prob.model.connect('re', point_name + '.re')
            prob.model.connect('rho', point_name + '.rho')
            prob.model.connect('cg', point_name + '.cg')

            # Connect the parameters within the model for each aero point
            for surface in surfaces:
                name = surface['name']

                # Connect the drag coeff at each point to the multi_CD component, which does the summation.
                prob.model.connect(point_name + '.CD', 'multi_CD.' + str(i) + '_CD')

                # Connect the mesh from the geometry component to the analysis point
                prob.model.connect(name + '_geom.mesh', point_name + '.' + name + '.def_mesh')

                # Perform the connections with the modified names within the 'aero_states' group.
                prob.model.connect(name + '_geom.mesh', point_name + '.aero_states.' + name + '_def_mesh')
                prob.model.connect(name + '_geom.t_over_c', point_name + '.' + name + '_perf.' + 't_over_c')

        prob.model.add_subsystem('multi_CD', MultiCD(n_points=n_points), promotes_outputs=['CD'])

        prob.driver = om.ScipyOptimizeDriver()
        prob.driver.options['tol'] = 1e-9

        # Setup problem and add design variables, constraint, and objective
        # design variables are the wing twist and angle-of-attack at each point.
        prob.model.add_design_var('alpha', lower=-15, upper=15)
        prob.model.add_design_var('wing_geom.twist_cp', lower=-5, upper=8)

        # set different target CL value at each point.
        prob.model.add_constraint('aero_point_0.wing_perf.CL', equals=0.45)
        prob.model.add_constraint('aero_point_1.wing_perf.CL', equals=0.5)

        # objective is the sum of CDs at each point.
        prob.model.add_objective('CD', scaler=1e4)

        # Set up the problem and run optimization.
        prob.setup()
        prob.run_driver()

        assert_rel_error(self, prob['aero_point_0.wing_perf.CL'][0], 0.45, 1e-6)
        assert_rel_error(self, prob['aero_point_0.wing_perf.CD'][0], 0.0323165143, 1e-6)
        assert_rel_error(self, prob['aero_point_1.wing_perf.CL'][0], 0.5, 1e-6)
        assert_rel_error(self, prob['aero_point_1.wing_perf.CD'][0], 0.0337665178, 1e-6)


if __name__ == '__main__':
    unittest.main()
