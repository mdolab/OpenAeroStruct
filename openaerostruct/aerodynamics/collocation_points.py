import numpy as np

import openmdao.api as om


class CollocationPoints(om.ExplicitComponent):
    """
    Compute the Cartesian locations of the collocation points, the force
    analysis points, and the bound vortex vectors for the VLM analysis.
    These points are 3/4 back of the front of the panel in the
    chordwise direction, and halfway across the panel in the spanwise direction.

    We enforce the flow tangency condition at these collocation points when
    solving for the circulations of the lifting surfaces.

    Parameters
    ----------
    def_mesh[nx, ny, 3] : numpy array
        We have a mesh for each lifting surface in the problem.
        That is, if we have both a wing and a tail surface, we will have both
        `wing_def_mesh` and `tail_def_mesh` as inputs.

    Returns
    -------
    coll_pts[num_eval_points, 3] : numpy array
        The xyz coordinates of the collocation points used in the VLM analysis.
        This array contains points for all lifting surfaces in the problem.
    force_pts[num_eval_points, 3] : numpy array
        The xyz coordinates of the force points used in the VLM analysis.
        We evaluate the velocity of the air at these points to get the sectional
        forces acting on the panel. This includes both the freestream and the
        induced velocity acting at these points.
        This array contains points for all lifting surfaces in the problem.
    bound_vecs[num_eval_points, 3] : numpy array
        The vectors representing the bound vortices for each panel in the
        problem.
        This array contains points for all lifting surfaces in the problem.

    """

    def initialize(self):
        self.options.declare('surfaces', types=list)

    def setup(self):
        num_eval_points = 0

        # Loop through all the surfaces to determine the total number
        # of evaluation points.
        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]

            num_eval_points += (nx - 1) * (ny - 1)

        self.add_output('coll_pts', shape=(num_eval_points, 3), units='m')
        self.add_output('force_pts', shape=(num_eval_points, 3), units='m')
        self.add_output('bound_vecs', shape=(num_eval_points, 3), units='m')

        eval_indices = np.arange(num_eval_points * 3).reshape((num_eval_points, 3))

        ind_eval_points_1 = 0
        ind_eval_points_2 = 0
        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            # Keep track of how many evaluation points come from this surface.
            ind_eval_points_2 += (nx - 1) * (ny - 1)

            # Take in a deformed mesh for each surface.
            mesh_name = name + '_def_mesh'
            self.add_input(mesh_name, shape=(nx, ny, 3), units='m')

            mesh_indices = np.arange(nx * ny * 3).reshape(
                (nx, ny, 3))

            # Compute the Jacobian for `coll_pts` wrt the meshes.
            # These do not change; the Jacobian is linear.
            rows = np.tile(eval_indices[ind_eval_points_1:ind_eval_points_2, :].flatten(), 4)
            cols = np.concatenate([
                mesh_indices[0:-1, 0:-1, :].flatten(),
                mesh_indices[1:  , 0:-1, :].flatten(),
                mesh_indices[0:-1, 1:  , :].flatten(),
                mesh_indices[1:  , 1:  , :].flatten(),
            ])
            data = np.concatenate([
                0.25 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # FR
                0.75 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # BR
                0.25 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # FL
                0.75 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # BL
            ])
            self.declare_partials('coll_pts', mesh_name, val=data, rows=rows, cols=cols)

            # Compute the Jacobian for `force_pts` wrt the meshes.
            # These do not change; the Jacobian is linear.
            data = np.concatenate([
                0.75 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # FR
                0.25 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # BR
                0.75 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # FL
                0.25 * 0.5 * np.ones((nx - 1) * (ny - 1) * 3),  # BL
            ])
            self.declare_partials('force_pts', mesh_name, val=data, rows=rows, cols=cols)

            # Compute the Jacobian for `bound_vecs` wrt the meshes.
            # These do not change; the Jacobian is linear.
            data = np.concatenate([
                 0.75 * np.ones((nx - 1) * (ny - 1) * 3),  # FR
                 0.25 * np.ones((nx - 1) * (ny - 1) * 3),  # BR
                -0.75 * np.ones((nx - 1) * (ny - 1) * 3),  # FL
                -0.25 * np.ones((nx - 1) * (ny - 1) * 3),  # BL
            ])
            self.declare_partials('bound_vecs', mesh_name, val=data, rows=rows, cols=cols)

            ind_eval_points_1 += (nx - 1) * (ny - 1)

    def compute(self, inputs, outputs):
        ind_eval_points_1 = 0
        ind_eval_points_2 = 0

        # Loop through each surface and compute the corresponding outputs,
        # paying special attention to the total number of evaluation points
        # in the system and each surface's place within the final arrays.
        for surface in self.options['surfaces']:
            mesh = surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            ind_eval_points_2 += (nx - 1) * (ny - 1)

            mesh_name = name + '_def_mesh'

            # The collocation points are 3/4 chord down the panel and in the
            # midpoint spanwise.
            outputs['coll_pts'][ind_eval_points_1:ind_eval_points_2, :] = (
                0.25 * 0.5 * inputs[mesh_name][0:-1, 0:-1, :] +
                0.75 * 0.5 * inputs[mesh_name][1:  , 0:-1, :] +
                0.25 * 0.5 * inputs[mesh_name][0:-1, 1:  , :] +
                0.75 * 0.5 * inputs[mesh_name][1:  , 1:  , :]
            ).reshape(((nx - 1) * (ny - 1), 3))

            # The force points are 1/4 chord down the panel and in the midpoint
            # spanwise.
            outputs['force_pts'][ind_eval_points_1:ind_eval_points_2, :] = (
                0.75 * 0.5 * inputs[mesh_name][0:-1, 0:-1, :] +
                0.25 * 0.5 * inputs[mesh_name][1:  , 0:-1, :] +
                0.75 * 0.5 * inputs[mesh_name][0:-1, 1:  , :] +
                0.25 * 0.5 * inputs[mesh_name][1:  , 1:  , :]
            ).reshape(((nx - 1) * (ny - 1), 3))

            # The bound vectors are computed at the 1/4 chord line.
            outputs['bound_vecs'][ind_eval_points_1:ind_eval_points_2, :] = (
                 0.75 * inputs[mesh_name][0:-1, 0:-1, :] +
                 0.25 * inputs[mesh_name][1:  , 0:-1, :] +
                -0.75 * inputs[mesh_name][0:-1, 1:  , :] +
                -0.25 * inputs[mesh_name][1:  , 1:  , :]
            ).reshape(((nx - 1) * (ny - 1), 3))

            # Increment the indices based on the amount contributed by this
            # surface.
            ind_eval_points_1 += (nx - 1) * (ny - 1)
