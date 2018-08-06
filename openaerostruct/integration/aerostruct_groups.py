from openaerostruct.geometry.geometry_mesh import GeometryMesh
from openaerostruct.aerodynamics.geometry import VLMGeometry
from openaerostruct.geometry.geometry_group import Geometry
from openaerostruct.transfer.displacement_transfer import DisplacementTransfer
from openaerostruct.structures.section_properties_tube import SectionPropertiesTube
from openaerostruct.structures.spatial_beam_setup import SpatialBeamSetup
from openaerostruct.structures.spatial_beam_states import SpatialBeamStates
from openaerostruct.aerodynamics.functionals import VLMFunctionals
from openaerostruct.structures.spatial_beam_functionals import SpatialBeamFunctionals
from openaerostruct.functionals.total_performance import TotalPerformance
from openaerostruct.transfer.load_transfer import LoadTransfer
from openaerostruct.aerodynamics.states import VLMStates
from openaerostruct.structures.tube_group import TubeGroup
from openaerostruct.structures.wingbox_group import WingboxGroup

from openmdao.api import IndepVarComp, Problem, Group, NewtonSolver, ScipyIterativeSolver, LinearBlockGS, NonlinearBlockGS, DirectSolver, LinearBlockGS, LinearRunOnce, ExplicitComponent, PetscKSP


class Aerostruct(Group):

    def initialize(self):
        self.options.declare('surface', types=dict)
        self.options.declare('DVGeo', default=None)

    def setup(self):
        surface = self.options['surface']
        DVGeo = self.options['DVGeo']

        geom_promotes = []

        if 'twist_cp' in surface.keys():
            geom_promotes.append('twist_cp')
        if 'toverc_cp' in surface.keys():
            geom_promotes.append('toverc_cp')
            geom_promotes.append('toverc')
        if 'mx' in surface.keys():
            geom_promotes.append('shape')

        self.add_subsystem('geometry',
            Geometry(surface=surface, DVGeo=DVGeo),
            promotes_inputs=[],
            promotes_outputs=['mesh'] + geom_promotes)

        if surface['fem_model_type'] == 'tube':
            tube_promotes = []
            if 'thickness_cp' in surface.keys():
                tube_promotes.append('thickness_cp')
            self.add_subsystem('tube_group',
                TubeGroup(surface=surface),
                promotes_inputs=['mesh'],
                promotes_outputs=['A', 'Iy', 'Iz', 'J', 'radius', 'thickness'] + tube_promotes)
        elif surface['fem_model_type'] == 'wingbox':
            wingbox_promotes = []
            if 'skin_thickness_cp' in surface.keys() and 'spar_thickness_cp' in surface.keys():
                wingbox_promotes.append('skin_thickness_cp')
                wingbox_promotes.append('spar_thickness_cp')
                wingbox_promotes.append('skin_thickness')
                wingbox_promotes.append('spar_thickness')
            elif 'skin_thickness_cp' in surface.keys() or 'spar_thickness_cp' in surface.keys():
                raise NameError('Please have both skin and spar thickness as design variables, not one or the other.')

            self.add_subsystem('wingbox_group',
                WingboxGroup(surface=surface),
                promotes_inputs=['mesh', 'toverc'],
                promotes_outputs=['A', 'Iy', 'Iz', 'J', 'Qz', 'A_enc', 'htop', 'hbottom', 'hfront', 'hrear'] + wingbox_promotes)
        else:
            raise NameError('Please select a valid `fem_model_type` from either `tube` or `wingbox`.')

        self.add_subsystem('struct_setup',
            SpatialBeamSetup(surface=surface),
            promotes_inputs=['mesh', 'A', 'Iy', 'Iz', 'J', 'load_factor'],
            promotes_outputs=['nodes', 'K', 'structural_weight', 'cg_location', 'element_weights'])

class CoupledAS(Group):

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.add_subsystem('struct_states',
            SpatialBeamStates(surface=surface),
            promotes_inputs=['K', 'forces', 'loads', 'element_weights'], promotes_outputs=['disp'])

        self.add_subsystem('def_mesh',
            DisplacementTransfer(surface=surface),
            promotes_inputs=['mesh', 'disp'], promotes_outputs=['def_mesh'])

        self.add_subsystem('aero_geom',
            VLMGeometry(surface=surface),
            promotes_inputs=['def_mesh'], promotes_outputs=['b_pts', 'c_pts', 'widths', 'cos_sweep', 'lengths', 'chords', 'normals', 'S_ref'])

        self.linear_solver = LinearRunOnce()


class CoupledPerformance(Group):

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.add_subsystem('aero_funcs',
            VLMFunctionals(surface=surface),
            promotes_inputs=['v', 'alpha', 'M', 're', 'rho', 'widths', 'cos_sweep', 'lengths', 'S_ref', 'sec_forces'], promotes_outputs=['CDv', 'L', 'D', 'CL1', 'CDi', 'CD', 'CL'])

        if surface['fem_model_type'] == 'tube':
            self.add_subsystem('struct_funcs',
                SpatialBeamFunctionals(surface=surface),
                promotes_inputs=['thickness', 'radius', 'nodes', 'disp'], promotes_outputs=['thickness_intersects', 'vonmises', 'failure'])

        elif surface['fem_model_type'] == 'wingbox':
            self.add_subsystem('struct_funcs',
                SpatialBeamFunctionals(surface=surface),
                promotes_inputs=['Qz', 'Iz', 'J', 'A_enc', 'spar_thickness', 'skin_thickness', 'htop', 'hbottom', 'hfront', 'hrear', 'nodes', 'disp'], promotes_outputs=['vonmises', 'failure'])
        else:
            raise NameError('Please select a valid `fem_model_type` from either `tube` or `wingbox`.')


class AerostructPoint(Group):

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        surfaces = self.options['surfaces']

        coupled = Group()

        for surface in surfaces:

            name = surface['name']

            # Connect the output of the loads component with the FEM
            # displacement parameter. This links the coupling within the coupled
            # group that necessitates the subgroup solver.
            coupled.connect(name + '_loads.loads', name + '.loads')

            # Perform the connections with the modified names within the
            # 'aero_states' group.
            coupled.connect(name + '.def_mesh', 'aero_states.' + name + '_def_mesh')
            coupled.connect(name + '.b_pts', 'aero_states.' + name + '_b_pts')
            coupled.connect(name + '.c_pts', 'aero_states.' + name + '_c_pts')
            coupled.connect(name + '.normals', 'aero_states.' + name + '_normals')
            coupled.connect(name + '.cos_sweep', 'aero_states.' + name + '_cos_sweep')
            coupled.connect(name + '.widths', 'aero_states.' + name + '_widths')

            # Connect the results from 'coupled' to the performance groups
            coupled.connect(name + '.def_mesh', name + '_loads.def_mesh')
            coupled.connect('aero_states.' + name + '_sec_forces', name + '_loads.sec_forces')

            # Connect the results from 'aero_states' to the performance groups
            self.connect('coupled.aero_states.' + name + '_sec_forces', name + '_perf' + '.sec_forces')

            # Connection performance functional variables
            self.connect(name + '_perf.CL', 'total_perf.' + name + '_CL')
            self.connect(name + '_perf.CD', 'total_perf.' + name + '_CD')
            self.connect('coupled.aero_states.' + name + '_sec_forces', 'total_perf.' + name + '_sec_forces')
            self.connect('coupled.' + name + '.chords', name + '_perf.aero_funcs.chords')

            # Connect parameters from the 'coupled' group to the performance
            # groups for the individual surfaces.
            self.connect('coupled.' + name + '.disp', name + '_perf.disp')
            self.connect('coupled.' + name + '.S_ref', name + '_perf.S_ref')
            self.connect('coupled.' + name + '.widths', name + '_perf.widths')
            # self.connect('coupled.' + name + '.chords', name + '_perf.chords')
            self.connect('coupled.' + name + '.lengths', name + '_perf.lengths')
            self.connect('coupled.' + name + '.cos_sweep', name + '_perf.cos_sweep')

            # Connect parameters from the 'coupled' group to the total performance group.
            self.connect('coupled.' + name + '.S_ref', 'total_perf.' + name + '_S_ref')
            self.connect('coupled.' + name + '.widths', 'total_perf.' + name + '_widths')
            self.connect('coupled.' + name + '.chords', 'total_perf.' + name + '_chords')
            self.connect('coupled.' + name + '.b_pts', 'total_perf.' + name + '_b_pts')

            # Add components to the 'coupled' group for each surface.
            # The 'coupled' group must contain all components and parameters
            # needed to converge the aerostructural system.
            coupled_AS_group = CoupledAS(surface=surface)

            coupled.add_subsystem(name, coupled_AS_group)

            # TODO: add this info to the options
            # prob.model.add_options(surface['name'] + 'yield_stress', surface['yield'])
            # prob.model.add_options(surface['name'] + 'fem_origin', surface['fem_origin'])

        # Add a single 'aero_states' component for the whole system within the
        # coupled group.
        coupled.add_subsystem('aero_states',
            VLMStates(surfaces=surfaces),
            promotes_inputs=['v', 'alpha', 'rho', 'M'])

        # Explicitly connect parameters from each surface's group and the common
        # 'aero_states' group.
        for surface in surfaces:
            name = surface['name']

            # Add a loads component to the coupled group
            coupled.add_subsystem(name + '_loads', LoadTransfer(surface=surface))

        """
        ### Change the solver settings here ###
        """

        # Set solver properties for the coupled group
        # coupled.linear_solver = ScipyIterativeSolver()
        # coupled.linear_solver.precon = LinearRunOnce()

        coupled.nonlinear_solver = NonlinearBlockGS(use_aitken=True)
        coupled.nonlinear_solver.options['maxiter'] = 50
        coupled.nonlinear_solver.options['atol'] = 5e-6
        coupled.nonlinear_solver.options['rtol'] = 1e-12

        # coupled.jacobian = DenseJacobian()
        coupled.linear_solver = DirectSolver()

        # coupled.nonlinear_solver = NewtonSolver(solve_subsystems=True)
        # coupled.nonlinear_solver.options['maxiter'] = 50
        coupled.nonlinear_solver.options['iprint'] = 0

        """
        ### End change of solver settings ###
        """

        # Add the coupled group to the model problem
        self.add_subsystem('coupled', coupled, promotes_inputs=['v', 'alpha', 'rho', 'M'])

        for surface in surfaces:
            name = surface['name']

            # Add a performance group which evaluates the data after solving
            # the coupled system
            perf_group = CoupledPerformance(surface=surface)

            self.add_subsystem(name + '_perf', perf_group, promotes_inputs=["rho", "v", "alpha", "re", "M"])

        # Add functionals to evaluate performance of the system.
        # Note that only the interesting results are promoted here; not all
        # of the parameters.
        self.add_subsystem('total_perf',
                 TotalPerformance(surfaces=surfaces),
                 promotes_inputs=['CL', 'CD', 'v', 'rho', 'empty_cg', 'total_weight', 'CT', 'a', 'R', 'M', 'W0', 'load_factor'],
                 promotes_outputs=['L_equals_W', 'fuelburn', 'CM', 'cg'])
