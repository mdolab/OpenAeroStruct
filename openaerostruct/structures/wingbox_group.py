from openmdao.api import Group, ExplicitComponent, BsplinesComp
from openaerostruct.geometry.geometry_mesh import GeometryMesh
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.structures.spatial_beam_states import SpatialBeamStates
from openaerostruct.structures.spatial_beam_functionals import SpatialBeamFunctionals
from openaerostruct.structures.spatial_beam_setup import SpatialBeamSetup
from openaerostruct.structures.section_properties_wingbox import SectionPropertiesWingbox
from openaerostruct.structures.wingbox_geometry import WingboxGeometry

from openmdao.api import IndepVarComp, Group


class WingboxGroup(Group):
    """ Group that contains everything needed for a structural-only problem. """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']
        ny = surface['num_y']

        if 'spar_thickness_cp' in surface.keys() or 'skin_thickness_cp' in surface.keys():
            # Add independent variables that do not belong to a specific component
            indep_var_comp = IndepVarComp()

            # Add structural components to the surface-specific group
            self.add_subsystem('indep_vars',
                     indep_var_comp,
                     promotes=['*'])

        if 'spar_thickness_cp' in surface.keys():
            n_cp = len(surface['spar_thickness_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('spar_thickness_bsp', BsplinesComp(
                in_name='spar_thickness_cp', out_name='spar_thickness',
                num_control_points=n_cp, num_points=int(ny-1),
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['spar_thickness_cp'], promotes_outputs=['spar_thickness'])
            indep_var_comp.add_output('spar_thickness_cp', val=surface['spar_thickness_cp'], units='m')

        if 'skin_thickness_cp' in surface.keys():
            n_cp = len(surface['skin_thickness_cp'])
            # Add bspline components for active bspline geometric variables.
            self.add_subsystem('skin_thickness_bsp', BsplinesComp(
                in_name='skin_thickness_cp', out_name='skin_thickness',
                num_control_points=n_cp, num_points=int(ny-1),
                bspline_order=min(n_cp, 4), distribution='uniform'),
                promotes_inputs=['skin_thickness_cp'], promotes_outputs=['skin_thickness'])
            indep_var_comp.add_output('skin_thickness_cp', val=surface['skin_thickness_cp'], units='m')

        self.add_subsystem('wingbox_geometry',
            WingboxGeometry(surface=surface),
            promotes_inputs=['mesh'],
            promotes_outputs=['fem_chords', 'fem_twists', 'streamwise_chords'])

        self.add_subsystem('wingbox',
            SectionPropertiesWingbox(surface=surface),
            promotes_inputs=['spar_thickness', 'skin_thickness', 'toverc', 'fem_chords', 'fem_twists', 'streamwise_chords'],
            promotes_outputs=['A', 'Iy', 'Qz', 'Iz', 'J', 'A_enc', 'htop', 'hbottom', 'hfront', 'hrear'])
