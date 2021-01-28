"""
Class definition for CompressibleVLMStates.
"""
import openmdao.api as om

from openaerostruct.aerodynamics.get_vectors import GetVectors
from openaerostruct.aerodynamics.collocation_points import CollocationPoints
from openaerostruct.aerodynamics.eval_mtx import EvalVelMtx
from openaerostruct.aerodynamics.convert_velocity import ConvertVelocity
from openaerostruct.aerodynamics.mtx_rhs import VLMMtxRHSComp
from openaerostruct.aerodynamics.solve_matrix import SolveMatrix
from openaerostruct.aerodynamics.horseshoe_circulations import HorseshoeCirculations
from openaerostruct.aerodynamics.eval_velocities import EvalVelocities
from openaerostruct.aerodynamics.mesh_point_forces import MeshPointForces
from openaerostruct.aerodynamics.panel_forces import PanelForces
from openaerostruct.aerodynamics.panel_forces_surf import PanelForcesSurf
from openaerostruct.aerodynamics.rotational_velocity import RotationalVelocity
from openaerostruct.aerodynamics.vortex_mesh import VortexMesh

from openaerostruct.aerodynamics.pg_transform import PGTransform, InversePGTransform


class CompressibleVLMStates(om.Group):
    """
    Group that contains the states for a compressible aerodynamic analysis.

    This is done in three steps:

        1. Convert all VLM geometries to Prandtl-Glaert domain.

        2. Solve the VLM problem in the Prandtl-Glaert domain using
        VLMStates group.

        3. Convert the resulting forces back to the physical domain through
        the inverse Prandtl-Glauert transform to recover the compressible
        forces.
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

        #----------------------------------------------------------------
        # Step 0: Need to calculate a few things prior to transformation.
        #----------------------------------------------------------------

        # Get collocation points
        self.add_subsystem('collocation_points',
             CollocationPoints(surfaces=surfaces),
             promotes_inputs=['*'])

        # Convert freestream velocity to array of velocities
        if rotational:
            self.add_subsystem('rotational_velocity',
                 RotationalVelocity(surfaces=surfaces),
                 promotes_inputs=['cg', 'omega'])

            self.connect('collocation_points.coll_pts', 'rotational_velocity.coll_pts')

        #----------------------------------------
        # Step 1: Transform geometry to PG domain
        #----------------------------------------

        self.connect('collocation_points.coll_pts', 'pg_transform.coll_pts')
        self.connect('collocation_points.bound_vecs', 'pg_transform.bound_vecs')
        self.connect('collocation_points.force_pts', 'pg_transform.force_pts')
        if rotational:
            self.connect('rotational_velocity.rotational_velocities', 'pg_transform.rotational_velocities')

        prom_in = ['alpha', 'beta', 'Mach_number']
        for surface in surfaces:
            name = surface['name']
            vname = name + '_def_mesh'
            prom_in.append(vname)
            self.connect('pg_transform.' + vname + '_pg', 'vortex_mesh.' + vname)

            vname = name + '_normals'
            prom_in.append(vname)
            self.connect('pg_transform.' + vname + '_pg', 'mtx_rhs.' + vname)

        self.add_subsystem('pg_transform', PGTransform(surfaces=surfaces, rotational=rotational),
            promotes_inputs=prom_in)

        self.connect('pg_transform.coll_pts_pg', 'coll_pts')
        self.connect('pg_transform.bound_vecs_pg', 'bound_vecs')
        self.connect('pg_transform.force_pts_pg', 'force_pts')
        if rotational:
            self.connect('pg_transform.rotational_velocities_pg', 'rotational_velocities')

        #---------------------------------------------------
        # Step 2: Solve incompressible problem in PG domain
        #---------------------------------------------------

        # Compute the vortex mesh based off the deformed aerodynamic mesh
        self.add_subsystem('vortex_mesh',
            VortexMesh(surfaces=surfaces),
            promotes_outputs=['*'])

        # Get vectors from mesh points to collocation points
        self.add_subsystem('get_vectors',
             GetVectors(surfaces=surfaces, num_eval_points=num_collocation_points,
                eval_name='coll_pts'),
             promotes_inputs=['*'],
             promotes_outputs=['*'])

        # In the PG domain, alpha and beta are zero.
        indep_var_comp = om.IndepVarComp()
        indep_var_comp.add_output('alpha_pg', val=0., units='deg')
        indep_var_comp.add_output('beta_pg', val=0., units='deg')
        self.add_subsystem('pg_frame', indep_var_comp)

        # Construct matrix based on rings, not horseshoes
        self.add_subsystem('mtx_assy',
             EvalVelMtx(surfaces=surfaces, num_eval_points=num_collocation_points,
                eval_name='coll_pts'),
             promotes_inputs=['*_vectors'],
             promotes_outputs=['*'])

        self.connect('pg_frame.alpha_pg', 'mtx_assy.alpha')

        # Convert freestream velocity to array of velocities
        # Note, don't want to promote Alpha or Beta here because we are in the transformed system.
        prom_in = ['v']
        if rotational:
            prom_in.append('rotational_velocities')
        self.add_subsystem('convert_velocity',
             ConvertVelocity(surfaces=surfaces, rotational=rotational),
             promotes_inputs=prom_in,
             promotes_outputs=['*'])

        self.connect('pg_frame.alpha_pg', 'convert_velocity.alpha')
        self.connect('pg_frame.beta_pg', 'convert_velocity.beta')

        # Construct RHS and full matrix of system
        self.add_subsystem('mtx_rhs',
             VLMMtxRHSComp(surfaces=surfaces),
             promotes_inputs=['freestream_velocities', '*coll_pts_vel_mtx'],
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
        # Note, don't want to promote Alpha here because we are in the transformed system.
        self.add_subsystem('mtx_assy_forces',
             EvalVelMtx(surfaces=surfaces, num_eval_points=num_force_points,
                eval_name='force_pts'),
             promotes_inputs=['*_force_pts_vectors'],
             promotes_outputs=['*'])

        self.connect('pg_frame.alpha_pg', 'mtx_assy_forces.alpha')

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
             promotes_inputs=['*'])

        #-----------------------------------------------------------
        # Step 3: Transform forces from PG domain to physical domain
        #-----------------------------------------------------------

        prom_out = []
        for surface in surfaces:
            vname = surface['name'] + '_sec_forces'
            prom_out.append(vname)
            self.connect('panel_forces_surf.' + vname, 'inverse_pg_transform.' + vname + '_pg')

        self.add_subsystem('inverse_pg_transform', InversePGTransform(surfaces=surfaces),
            promotes_inputs=['alpha', 'beta', 'Mach_number'],
            promotes_outputs=prom_out)

        #---------------------------------------------------------------
        # Step 4: Mesh point forces are downstream, already transformed.
        #---------------------------------------------------------------

        # Get nodal forces for each lifting surface individually
        self.add_subsystem('mesh_point_forces_surf',
             MeshPointForces(surfaces=surfaces),
             promotes_inputs=['*'],
             promotes_outputs=['*'])
