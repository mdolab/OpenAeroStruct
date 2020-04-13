import openmdao.api as om
from openaerostruct.structures.section_properties_wingbox import SectionPropertiesWingbox
from openaerostruct.structures.wingbox_geometry import WingboxGeometry


class WingboxGroup(om.Group):
    """ Group that contains everything needed for a structural-only problem. """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']
        ny = surface['mesh'].shape[1]

        if 'spar_thickness_cp' in surface.keys() or 'skin_thickness_cp' in surface.keys():
            # Add independent variables that do not belong to a specific component
            indep_var_comp = om.IndepVarComp()

            # Add structural components to the surface-specific group
            self.add_subsystem('indep_vars',
                     indep_var_comp,
                     promotes=['*'])

        if 'spar_thickness_cp' in surface.keys():
            n_cp = len(surface['spar_thickness_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('spar_thickness_bsp', om.BsplinesComp(
                in_name='spar_thickness_cp', out_name='spar_thickness',
                num_control_points=n_cp, num_points=int(ny-1), units='m',
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['spar_thickness_cp'], promotes_outputs=['spar_thickness'])
            indep_var_comp.add_output('spar_thickness_cp', val=surface['spar_thickness_cp'], units='m')

        if 'skin_thickness_cp' in surface.keys():
            n_cp = len(surface['skin_thickness_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('skin_thickness_bsp', om.BsplinesComp(
                in_name='skin_thickness_cp', out_name='skin_thickness',
                num_control_points=n_cp, num_points=int(ny-1), units='m',
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['skin_thickness_cp'], promotes_outputs=['skin_thickness'])
            indep_var_comp.add_output('skin_thickness_cp', val=surface['skin_thickness_cp'], units='m')

        self.add_subsystem('wingbox_geometry',
            WingboxGeometry(surface=surface),
            promotes_inputs=['mesh'],
            promotes_outputs=['fem_chords', 'fem_twists', 'streamwise_chords'])

        self.add_subsystem('wingbox',
            SectionPropertiesWingbox(surface=surface),
            promotes_inputs=['spar_thickness', 'skin_thickness', 't_over_c', 'fem_chords', 'fem_twists', 'streamwise_chords'],
            promotes_outputs=['A', 'Iy', 'Qz', 'Iz', 'J', 'A_enc', 'A_int', 'htop', 'hbottom', 'hfront', 'hrear'])
