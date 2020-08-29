import numpy as np

import openmdao.api as om

from openaerostruct.utils.vector_algebra import get_array_indices


np.random.seed(314)

class DisplacementTransfer(om.ExplicitComponent):
    """
    Apply the computed FEM displacements and rotations on the aerodynamic mesh
    to obtain the deformed mesh.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Original undeformed aerodynamic mesh.
    disp[ny, 6] : numpy array
        Displacements and rotations acting on the structural spar which come
        from solving the FEM system. Contains displacements for all six degrees
        of freedom, including displacements in the x, y, and z directions, and
        rotations about the x, y, and z axes.
    transformation_matrix[ny, 3, 3] : numpy array
        Array containing the transformation matrices to apply the rotations
        from the FEM results. These are the angles that come from the FEM
        solver that rotate the structural spar. They will rotate the mesh based
        on the both the angle (which comes from `disp`) and from the distance
        of the aerodynamic mesh nodes to the structural mesh nodes (basically
        a moment arm).

    Returns
    -------
    def_mesh[nx, ny, 3] : numpy array
        The final deformed aerodynamic mesh for the lifting surface based on
        the FEM results.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        mesh=surface['mesh']
        self.nx = mesh.shape[0]
        self.ny = mesh.shape[1]

        self.add_input('mesh', val=np.ones((self.nx, self.ny, 3)), units='m')
        self.add_input('disp', val=np.ones((self.ny, 6)), units='m')
        self.add_input('transformation_matrix', val=np.ones((self.ny, 3, 3)))
        self.add_input('nodes', val=np.ones((self.ny, 3)), units='m')

        self.add_output('def_mesh', val=np.random.random_sample((self.nx, self.ny, 3)), units='m')

        # Create index arrays for each relevant input and output.
        # This allows us to set up the rows and cols for the sparse Jacobians.
        disp_indices = get_array_indices(self.ny, 6)
        nodes_indices = get_array_indices(self.ny, 3)
        mesh_indices = get_array_indices(self.nx, self.ny, 3)
        transform_indices = get_array_indices(self.ny, 3, 3)
        mesh_disp_indices = get_array_indices(self.nx, self.ny, 3)

        # Set up the rows and cols for `def_mesh` wrt `disp`
        rows = mesh_disp_indices.flatten()
        cols = np.einsum('i,jk->ijk', np.ones(self.nx), disp_indices[:, :3]).flatten()
        self.declare_partials('def_mesh', 'disp', val=1., rows=rows, cols=cols)

        # Set up the rows and cols for `def_mesh` wrt `nodes`
        rows = np.einsum('ijk,l->ijkl', mesh_disp_indices, np.ones(3, int)).flatten()
        cols = np.einsum('ik,jl->ijkl', np.ones((self.nx, 3), int), nodes_indices).flatten()
        self.declare_partials('def_mesh', 'nodes', rows=rows, cols=cols)

        # Set up the rows and cols for `def_mesh` wrt `mesh`
        rows = np.einsum('ijk,l->ijkl', mesh_disp_indices, np.ones(3, int)).flatten()
        cols = np.einsum('ijl,k->ijkl', mesh_indices, np.ones(3, int)).flatten()
        self.declare_partials('def_mesh', 'mesh', rows=rows, cols=cols)

        # Set up the rows and cols for `def_mesh` wrt `transformation_matrix`
        rows = np.einsum('ijl,k->ijkl', mesh_disp_indices, np.ones(3, int)).flatten()
        cols = np.einsum('jlk,i->ijkl', transform_indices, np.ones(self.nx, int)).flatten()
        self.declare_partials('def_mesh', 'transformation_matrix', rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        # Get the location of the spar
        nodes = inputs['nodes']

        # First set the deformed mesh with the undeformed mesh values
        outputs['def_mesh'] = inputs['mesh'].copy()

        # Add the translational displacements to the deformed mesh.
        # These are simply the x,y,z displacements getting added to all nodal
        # mesh points.
        outputs['def_mesh'] += np.einsum('i,jk->ijk',
            np.ones(self.nx), inputs['disp'][:, :3])

        # Compute the moment arms from the aerodynamic mesh points to the
        # structural mesh points.
        moment_arms = inputs['mesh'] - nodes

        # Apply the transformation matrix to the moment arms to get the
        # rotational displacements from the FEM results transformed to the
        # aerodynamic mesh. Then add these to the deformed mesh.
        outputs['def_mesh'] += np.einsum('lij,klj->kli',
                                         inputs['transformation_matrix'],
                                         moment_arms)

    def compute_partials(self, inputs, partials):
        partials['def_mesh', 'nodes'] = -np.einsum('i,jlk->ijlk',
            np.ones(self.nx), inputs['transformation_matrix']).flatten()

        partials['def_mesh', 'mesh'] = np.einsum('i,jlk->ijlk',
            np.ones(self.nx), inputs['transformation_matrix']).flatten()
        partials['def_mesh', 'mesh'] += np.tile(np.eye(3), self.nx * self.ny).flatten(order='F')

        partials['def_mesh', 'transformation_matrix'] = np.einsum('ijk,l->ijkl',
            inputs['mesh'] - np.einsum('i,jk->ijk', np.ones(self.nx), inputs['nodes']),
            np.ones(3)).flatten()
