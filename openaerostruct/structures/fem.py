"""Define the LinearSystemComp class."""

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu

import openmdao.api as om


class FEM(om.ImplicitComponent):
    """
    Component that solves a linear system, Ax=b.

    Designed to handle small, dense linear systems (Ax=B) that can be efficiently solved with
    sparse lu-decomposition. It can be vectorized to either solve for multiple right hand sides,
    or to solve multiple linear systems.

    A is represented sparsely as a local_stiff_transformed, which is an ny x 12 x 12 array.

    Attributes
    ----------
    _lup : None or list(object)
        matrix factorizations returned from scipy.linag.lu_factor for each A matrix
    k_cols : ndarray
        Cached column indices for sparse representation of stiffness matrix.
    k_rows : ndarray
        Cached row indices for sparse representation of stiffness matrix.
    k_data : ndarray
        Cached values for sparse representation of stiffness matrix.
    """

    def __init__(self, **kwargs):
        """
        Intialize the LinearSystemComp component.

        Parameters
        ----------
        **kwargs : dict of keyword arguments
            Keyword arguments that will be mapped into the Component options.
        """
        super(FEM, self).__init__(**kwargs)
        self._lup = None
        self.k_cols = None
        self.k_rows = None
        self.k_data = None

    def initialize(self):
        """
        Declare options.
        """
        self.options.declare('surface', types=dict)
        self.options.declare('vec_size', types=int, default=1,
                             desc='Number of linear systems to solve.')

    def setup(self):
        """
        Matrix and RHS are inputs, solution vector is the output.
        """
        surface = self.options['surface']
        self.ny = ny = surface['mesh'].shape[1]
        self.size = size = int(6 * ny + 6)

        vec_size = self.options['vec_size']
        full_size = size * vec_size

        self._lup = []
        shape = (vec_size, size) if vec_size > 1 else (size, )

        init_locK = np.tile(np.eye(12).flatten(), ny-1).reshape(ny-1, 12, 12)

        self.add_input('local_stiff_transformed', val=init_locK)
        self.add_input('forces', val=np.ones(shape), units='N')
        self.add_output('disp_aug', shape=shape, val=.1, units='m')

        # Set up the derivatives.
        row_col = np.arange(full_size, dtype="int")

        self.declare_partials('disp_aug', 'forces', val=np.full(full_size, -1.0),
                              rows=row_col, cols=row_col)

        # The derivative of residual wrt displacements is the stiffness matrix K. We can use the
        # sparsity pattern here and when constucting the sparse matrix, so save rows and cols.

        base_row = np.repeat(0, 6)
        base_col = np.arange(6)

        # Upper diagonal blocks
        rows1 = np.tile(base_row, 6*(ny-1)) + np.repeat(np.arange(6*(ny-1)), 6)
        col = np.tile(base_col + 6, 6)
        cols1 = np.tile(col, ny-1) + np.repeat(6*np.arange(ny-1), 36)

        # Lower diagonal blocks
        rows2 = np.tile(base_row + 6, 6*(ny-1)) + np.repeat(np.arange(6*(ny-1)), 6)
        col = np.tile(base_col, 6)
        cols2 = np.tile(col, ny-1) + np.repeat(6*np.arange(ny-1), 36)

        # Main diagonal blocks, root
        rows3 = np.tile(base_row, 6) + np.repeat(np.arange(6), 6)
        cols3 = np.tile(base_col, 6)

        # Main diagonal blocks, tip
        rows4 = np.tile(base_row + (ny-1)*6, 6) + np.repeat(np.arange(6), 6)
        cols4 = np.tile(base_col + (ny-1)*6, 6)

        # Main diagonal blocks, interior
        rows5 = np.tile(base_row + 6, 6*(ny-2)) + np.repeat(np.arange(6*(ny-2)), 6)
        col = np.tile(base_col + 6, 6)
        cols5 = np.tile(col, ny-2) + np.repeat(6*np.arange(ny-2), 36)

        # Find constrained nodes based on closeness to specified cg point
        symmetry = self.options['surface']['symmetry']
        if symmetry:
            idx = self.ny - 1
        else:
            idx = (self.ny - 1) // 2

        index = 6 * idx
        num_dofs = 6 * ny
        arange = np.arange(6)

        # Fixed boundary condition.
        rows6 = index + arange
        cols6 = num_dofs + arange

        self.k_rows = rows = np.concatenate([rows1, rows2, rows3, rows4, rows5, rows6, cols6])
        self.k_cols = cols = np.concatenate([cols1, cols2, cols3, cols4, cols5, cols6, rows6])

        sp_size = len(rows)
        vec_rows = np.tile(rows, vec_size) + np.repeat(sp_size*np.arange(vec_size), sp_size)
        vec_cols = np.tile(cols, vec_size) + np.repeat(sp_size*np.arange(vec_size), sp_size)

        self.declare_partials(of='disp_aug', wrt='disp_aug', rows=vec_rows, cols=vec_cols)

        base_row = np.tile(0, 12)
        base_col = np.arange(12)
        row = np.tile(base_row, 12) + np.repeat(np.arange(12), 12)
        col = np.tile(base_col, 12) + np.repeat(12*np.arange(12), 12)
        rows = np.tile(row, ny-1) + np.repeat(6*np.arange(ny-1), 144)
        cols = np.tile(col, ny-1) + np.repeat(144*np.arange(ny-1), 144)

        self.declare_partials('disp_aug', 'local_stiff_transformed', rows=rows, cols=cols)

    def apply_nonlinear(self, inputs, outputs, residuals):
        """
        R = Ax - b.

        Parameters
        ----------
        inputs : Vector
            unscaled, dimensional input variables read via inputs[key]
        outputs : Vector
            unscaled, dimensional output variables read via outputs[key]
        residuals : Vector
            unscaled, dimensional residuals written to via residuals[key]
        """
        K = self.assemble_CSC_K(inputs)
        residuals['disp_aug'] = K.dot(outputs['disp_aug']) - inputs['forces']

    def solve_nonlinear(self, inputs, outputs):
        """
        Use numpy to solve Ax=b for x.

        Parameters
        ----------
        inputs : Vector
            unscaled, dimensional input variables read via inputs[key]
        outputs : Vector
            unscaled, dimensional output variables read via outputs[key]
        """
        # lu factorization for use with solve_linear
        K = self.assemble_CSC_K(inputs)
        self._lup = splu(K)
        outputs['disp_aug'] = self._lup.solve(inputs['forces'])

    def linearize(self, inputs, outputs, J):
        """
        Compute the non-constant partial derivatives.

        Parameters
        ----------
        inputs : Vector
            unscaled, dimensional input variables read via inputs[key]
        outputs : Vector
            unscaled, dimensional output variables read via outputs[key]
        J : Jacobian
            sub-jac components written to jacobian[output_name, input_name]
        """
        x = outputs['disp_aug']
        vec_size = self.options['vec_size']
        ny = self.ny

        idx = np.tile(np.tile(np.arange(12), 12), ny-1) + np.repeat(6*np.arange(ny-1), 144)
        J['disp_aug', 'local_stiff_transformed'] = np.tile(x[idx], vec_size)

        J['disp_aug', 'disp_aug'] = np.tile(self.k_data, vec_size)

    def solve_linear(self, d_outputs, d_residuals, mode):
        r"""
        Back-substitution to solve the derivatives of the linear system.

        If mode is:
            'fwd': d_residuals \|-> d_outputs

            'rev': d_outputs \|-> d_residuals

        Parameters
        ----------
        d_outputs : Vector
            unscaled, dimensional quantities read via d_outputs[key]
        d_residuals : Vector
            unscaled, dimensional quantities read via d_residuals[key]
        mode : str
            either 'fwd' or 'rev'
        """
        vec_size = self.options['vec_size']

        if mode == 'fwd':
            if vec_size > 1:
                for j in range(vec_size):
                    d_outputs['disp_aug'] = self._lup.solve(d_residuals['disp_aug'][j])
            else:
                d_outputs['disp_aug'] = self._lup.solve(d_residuals['disp_aug'])
        else:
            if vec_size > 1:
                for j in range(vec_size):
                    d_residuals['disp_aug'] = self._lup.solve(d_outputs['disp_aug'][j])
            else:
                d_residuals['disp_aug'] = self._lup.solve(d_outputs['disp_aug'])

    def assemble_CSC_K(self, inputs):
        """
        Assemble the stiffness matrix in sparse CSC format.

        Returns
        -------
        ndarray
            Stiffness matrix as dense ndarray.
        """
        k_loc = inputs['local_stiff_transformed']
        size = self.size

        data1 = k_loc[:, :6, 6:].flatten()
        data2 = k_loc[:, 6:, :6].flatten()
        data3 = k_loc[0, :6, :6].flatten()
        data4 = k_loc[-1, 6:, 6:].flatten()
        data5 = (k_loc[0:-1, 6:, 6:] + k_loc[1:, :6, :6]).flatten()
        data6 = np.full((6, ), 1e9)

        self.k_data = data = np.concatenate([data1, data2, data3, data4, data5, data6, data6])

        return coo_matrix((data, (self.k_rows, self.k_cols)), shape=(size, size)).tocsc()
