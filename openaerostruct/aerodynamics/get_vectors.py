import numpy as np

import openmdao.api as om


class GetVectors(om.ExplicitComponent):
    """
    Compute the vectors going from the vortex mesh points to the evaluation
    points, where the evluation points are either the colloation points
    or the force points used in the AIC computations.

    Parameters
    ----------
    vortex_mesh[nx, ny, 3] : numpy array
        The actual aerodynamic mesh used in VLM calculations, where we look
        at the rings of the panels instead of the panels themselves. That is,
        this mesh coincides with the quarter-chord panel line, except for the
        final row, where it lines up with the trailing edge.
    eval_name[num_eval_points, 3] : numpy array
        These are the evaluation points, either collocation or force points.

    Returns
    -------
    vectors[num_eval_points, nx, ny, 3] : numpy array
        The actual velocities experienced at the evaluation points for each
        lifting surface in the system. This is the summation of the freestream
        velocities and the induced velocities caused by the circulations.
    """

    def initialize(self):
        self.options.declare('surfaces', types=list)
        self.options.declare('num_eval_points', types=int)
        self.options.declare('eval_name', types=str)

    def setup(self):
        surfaces = self.options['surfaces']
        num_eval_points = self.options['num_eval_points']
        eval_name = self.options['eval_name']

        # Take in the evaluation points
        self.add_input(eval_name, val=np.zeros((num_eval_points, 3)), units='m')

        for surface in surfaces:
            mesh=surface['mesh']
            nx = mesh.shape[0]
            ny = mesh.shape[1]
            name = surface['name']

            ground_effect = surface.get('groundplane', False)

            vectors_name = '{}_{}_vectors'.format(name, eval_name)

            # This is where we handle the symmetry in the VLM method.
            # If it's symmetric, we need to effectively mirror the mesh by
            # accounting for the ghost mesh. We do this by using an artificially
            # larger mesh here.

            if surface['symmetry']:
                actual_ny_size = ny * 2 - 1
            else:
                actual_ny_size = ny
            if ground_effect:
                actual_nx_size = nx * 2
            else:
                actual_nx_size = nx

            self.add_input(name + '_vortex_mesh', val=np.zeros((actual_nx_size, actual_ny_size, 3)), units='m')
            self.add_output(vectors_name, val=np.ones((num_eval_points, actual_nx_size, actual_ny_size, 3)), units='m')

            # Set up indices so we can get the rows and cols for the delcare
            vector_indices = np.arange(num_eval_points * actual_nx_size * actual_ny_size * 3)
            mesh_indices = np.outer(
                np.ones(num_eval_points, int),
                np.arange(actual_nx_size * actual_ny_size * 3),
            ).flatten()
            eval_indices = np.einsum('il,jk->ijkl',
                np.arange(num_eval_points * 3).reshape((num_eval_points, 3)),
                np.ones((actual_nx_size, actual_ny_size), int),
            ).flatten()

            self.declare_partials(vectors_name, name + '_vortex_mesh', val=-1., rows=vector_indices, cols=mesh_indices)
            self.declare_partials(vectors_name, eval_name, val= 1., rows=vector_indices, cols=eval_indices)

    def compute(self, inputs, outputs):
        surfaces = self.options['surfaces']
        num_eval_points = self.options['num_eval_points']
        eval_name = self.options['eval_name']

        # At the end of the day, all this component is doing is computing
        # vectors that go from the mesh to the evaluation points. We have
        # a lot of vector algebra to make sure everything's lined up okay
        # and in a usable data format.

        for surface in surfaces:
            nx = surface['mesh'].shape[0]
            ny = surface['mesh'].shape[1]
            name = surface['name']

            mesh_name = name + '_vortex_mesh'
            vectors_name = '{}_{}_vectors'.format(name, eval_name)

            ground_effect = surface.get('groundplane', False)

            mesh_reshaped = np.einsum('i,jkl->ijkl', np.ones(num_eval_points), inputs[mesh_name])
            if surface['symmetry'] and ground_effect:
                eval_points_reshaped = np.einsum('il,jk->ijkl',
                    inputs[eval_name],
                    np.ones((2*nx, 2*ny-1)),
                )
            elif surface['symmetry']:
                eval_points_reshaped = np.einsum('il,jk->ijkl',
                    inputs[eval_name],
                    np.ones((nx, 2*ny-1)),
                )
            else:
                eval_points_reshaped = np.einsum('il,jk->ijkl',
                    inputs[eval_name],
                    np.ones((nx, ny)),
                )

            # Actually subtract the vectors.
            outputs[vectors_name] = eval_points_reshaped - mesh_reshaped
