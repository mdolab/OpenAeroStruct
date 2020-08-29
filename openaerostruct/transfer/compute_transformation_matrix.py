import numpy as np

import openmdao.api as om

from openaerostruct.utils.vector_algebra import get_array_indices


class ComputeTransformationMatrix(om.ExplicitComponent):
    """
    Compute the transformation matrix used to apply the rotations obtained
    from the FEM system to the aerodynamic mesh.

    Technically the order of rotations (x, y, z) matter, but here we assume
    they are independent. This allows us to construct a single transformation
    matrix and apply it to the moment arms going from the structural mesh
    to the aerodynamic mesh.

    Parameters
    ----------
    disp[ny, 6] : numpy array
        Displacements and rotations acting on the structural spar which come
        from solving the FEM system. Contains displacements for all six degrees
        of freedom, including displacements in the x, y, and z directions, and
        rotations about the x, y, and z axes.

    Returns
    -------
    transformation_matrix[ny, 3, 3] : numpy array
        Array containing the transformation matrices to apply the rotations
        from the FEM results. These are the angles that come from the FEM
        solver that rotate the structural spar. They will rotate the mesh based
        on the both the angle (which comes from `disp`) and from the distance
        of the aerodynamic mesh nodes to the structural mesh nodes (basically
        a moment arm).
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']
        mesh=surface['mesh']
        self.nx = mesh.shape[0]
        self.ny = mesh.shape[1]

        self.add_input('disp', val=np.zeros((self.ny, 6)), units='m')
        self.add_output('transformation_matrix', shape=(self.ny, 3, 3))

        # Create index arrays for each relevant input and output.
        # This allows us to set up the rows and cols for the sparse Jacobians.
        disp_indices = get_array_indices(self.ny, 6)
        transform_indices = get_array_indices(self.ny, 3, 3)

        # Set up the rows and cols for `transformation_matrix` wrt `disp`
        rows = np.einsum('ijk,l->ijkl',
            transform_indices,
            np.ones(3, int)).flatten()
        cols = np.einsum('il,jk->ijkl',
            get_array_indices(self.ny, 6)[:, 3:],
            np.ones((3, 3), int)).flatten()
        self.declare_partials('transformation_matrix', 'disp', rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        # We populate the diagonal of each transformation matrix to account
        # for the equivalent identity matrices we're adding in later steps.
        # We need to do this because we're treating the three rotational
        # matrices as independent and combining them all by adding.
        outputs['transformation_matrix'] = 0.
        for i in range(3):
            outputs['transformation_matrix'][:, i, i] -= 2.

        # These are the rotations obtained from the FEM solution
        rx = inputs['disp'][:, 3]
        ry = inputs['disp'][:, 4]
        rz = inputs['disp'][:, 5]

        # We then apply rotations to each corresponding entry in the
        # transformation matrix, cumulatively adding or subtracting to
        # obtain the final result.

        # T[ 1:,  1:] += [[cos(rx), -sin(rx)], [ sin(rx), cos(rx)]]
        outputs['transformation_matrix'][:, 1, 1] += np.cos(rx)
        outputs['transformation_matrix'][:, 1, 2] -= np.sin(rx)
        outputs['transformation_matrix'][:, 2, 1] += np.sin(rx)
        outputs['transformation_matrix'][:, 2, 2] += np.cos(rx)

        # T[::2, ::2] += [[cos(ry),  sin(ry)], [-sin(ry), cos(ry)]]
        outputs['transformation_matrix'][:, 0, 0] += np.cos(ry)
        outputs['transformation_matrix'][:, 0, 2] += np.sin(ry)
        outputs['transformation_matrix'][:, 2, 0] -= np.sin(ry)
        outputs['transformation_matrix'][:, 2, 2] += np.cos(ry)

        # T[ :2,  :2] += [[cos(rz), -sin(rz)], [ sin(rz), cos(rz)]]
        outputs['transformation_matrix'][:, 0, 0] += np.cos(rz)
        outputs['transformation_matrix'][:, 0, 1] -= np.sin(rz)
        outputs['transformation_matrix'][:, 1, 0] += np.sin(rz)
        outputs['transformation_matrix'][:, 1, 1] += np.cos(rz)

    def compute_partials(self, inputs, partials):
        # Because of the way this is constructed, the derivatives are rather
        # straightforward. We simply take the trigonometric deriv for each
        # entry in the matrix.
        rx = inputs['disp'][:, 3]
        ry = inputs['disp'][:, 4]
        rz = inputs['disp'][:, 5]

        derivs = partials['transformation_matrix', 'disp'].reshape((self.ny, 3, 3, 3))
        derivs[:, :, :, :] = 0.

        derivs[:, 1, 1, 0] -= np.sin(rx)
        derivs[:, 1, 2, 0] -= np.cos(rx)
        derivs[:, 2, 1, 0] += np.cos(rx)
        derivs[:, 2, 2, 0] -= np.sin(rx)

        derivs[:, 0, 0, 1] -= np.sin(ry)
        derivs[:, 0, 2, 1] += np.cos(ry)
        derivs[:, 2, 0, 1] -= np.cos(ry)
        derivs[:, 2, 2, 1] -= np.sin(ry)

        derivs[:, 0, 0, 2] -= np.sin(rz)
        derivs[:, 0, 1, 2] -= np.cos(rz)
        derivs[:, 1, 0, 2] += np.cos(rz)
        derivs[:, 1, 1, 2] -= np.sin(rz)
