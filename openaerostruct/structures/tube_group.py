import openmdao.api as om
from openaerostruct.structures.section_properties_tube import SectionPropertiesTube
from openaerostruct.geometry.radius_comp import RadiusComp


class TubeGroup(om.Group):
    """ Group that contains everything needed for a structural-only problem. """

    def initialize(self):
        self.options.declare('surface', types=dict)
        self.options.declare('connect_geom_DVs', default=True)

    def setup(self):
        surface = self.options['surface']
        connect_geom_DVs = self.options['connect_geom_DVs']
        mesh = surface['mesh']
        ny = mesh.shape[1]

        if connect_geom_DVs:
            # Add independent variables that do not belong to a specific component
            indep_var_comp = om.IndepVarComp()

            # Add structural components to the surface-specific group
            self.add_subsystem('indep_vars',
                     indep_var_comp,
                     promotes=['*'])

        if 'thickness_cp' in surface.keys():
            n_cp = len(surface['thickness_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('thickness_bsp', om.BsplinesComp(
                in_name='thickness_cp', out_name='thickness', units='m',
                num_control_points=n_cp, num_points=int(ny-1),
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['thickness_cp'], promotes_outputs=['thickness'])
            if connect_geom_DVs:
                indep_var_comp.add_output('thickness_cp', val=surface['thickness_cp'], units='m')

        if 'radius_cp' in surface.keys():
            n_cp = len(surface['radius_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('radius_bsp', om.BsplinesComp(
                in_name='radius_cp', out_name='radius', units='m',
                num_control_points=n_cp, num_points=int(ny-1),
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['radius_cp'], promotes_outputs=['radius'])
            if connect_geom_DVs:
                indep_var_comp.add_output('radius_cp', val=surface['radius_cp'], units='m')
        else:
            self.add_subsystem('radius_comp',
                RadiusComp(surface=surface),
                promotes_inputs=['mesh', 't_over_c'],
                promotes_outputs=['radius'])

        self.add_subsystem('tube',
            SectionPropertiesTube(surface=surface),
            promotes_inputs=['thickness', 'radius'],
            promotes_outputs=['A', 'Iy', 'Iz', 'J'])
