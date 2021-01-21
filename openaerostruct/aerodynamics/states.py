import openmdao.api as om
from openaerostruct.aerodynamics.get_vectors import GetVectors
from openaerostruct.aerodynamics.collocation_points import CollocationPoints
from openaerostruct.aerodynamics.eval_mtx import EvalVelMtx
from openaerostruct.aerodynamics.convert_velocity import ConvertVelocity
from openaerostruct.aerodynamics.mtx_rhs import VLMMtxRHSComp
from openaerostruct.aerodynamics.solve_matrix import SolveMatrix
from openaerostruct.aerodynamics.horseshoe_circulations import HorseshoeCirculations
from openaerostruct.aerodynamics.eval_velocities import EvalVelocities
from openaerostruct.aerodynamics.rotational_velocity import RotationalVelocity
from openaerostruct.aerodynamics.mesh_point_forces import MeshPointForces
from openaerostruct.aerodynamics.panel_forces import PanelForces
from openaerostruct.aerodynamics.panel_forces_surf import PanelForcesSurf
from openaerostruct.aerodynamics.vortex_mesh import VortexMesh


class VLMStates(om.Group):
    """
    Group that houses all components to compute the aerodynamic states.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('rotational', False, types=bool,
                             desc="Set to True to turn on support for computing angular velocities")

    def setup(self):
        surfaces = self.options['surfaces']
        rotational = self.options['rotational']

        num_collocation_points = 0
        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            num_collocation_points += (ny - 1) * (nx - 1)

        num_force_points = num_collocation_points

        # Get collocation points
        self.add_subsystem('collocation_points',
             CollocationPoints(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['coll_pts', 'force_pts', 'bound_vecs'])

        # Compute the vortex mesh based off the deformed aerodynamic mesh
        self.add_subsystem('vortex_mesh',
            VortexMesh(surfaces=surfaces),
            promotes_inputs=['*'],
            promotes_outputs=['*'])

        # Get vectors from mesh points to collocation points
        self.add_subsystem('get_vectors',
             GetVectors(surfaces=surfaces, num_eval_points=num_collocation_points,
                eval_name='coll_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Construct matrix based on rings, not horseshoes
        self.add_subsystem('mtx_assy',
             EvalVelMtx(surfaces=surfaces, num_eval_points=num_collocation_points,
                eval_name='coll_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Convert freestream velocity to array of velocities
        if rotational:
            self.add_subsystem('rotational_velocity',
                 RotationalVelocity(surfaces=surfaces),
                 promotes_inputs=['*'],
                 promotes_outputs=['*'])

        self.add_subsystem('convert_velocity',
             ConvertVelocity(surfaces=surfaces, rotational=rotational),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Construct RHS and full matrix of system
        self.add_subsystem('mtx_rhs',
             VLMMtxRHSComp(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Solve Mtx RHS to get ring circs
        self.add_subsystem('solve_matrix',
             SolveMatrix(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Convert ring circs to horseshoe circs
        self.add_subsystem('horseshoe_circulations',
             HorseshoeCirculations(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Eval force vectors
        self.add_subsystem('get_vectors_force',
             GetVectors(surfaces=surfaces, num_eval_points=num_force_points,
                eval_name='force_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Set up force mtx
        self.add_subsystem('mtx_assy_forces',
             EvalVelMtx(surfaces=surfaces, num_eval_points=num_force_points,
                eval_name='force_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Multiply by horseshoe circs to get velocities
        self.add_subsystem('eval_velocities',
             EvalVelocities(surfaces=surfaces, num_eval_points=num_force_points,
                eval_name='force_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Get sectional panel forces
        self.add_subsystem('panel_forces',
             PanelForces(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Get panel forces for each lifting surface individually
        self.add_subsystem('panel_forces_surf',
             PanelForcesSurf(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # Get nodal forces for each lifting surface individually
        self.add_subsystem('mesh_point_forces_surf',
             MeshPointForces(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])
