import numpy as np

import openmdao.api as om


row_indices = np.arange(12)
col_indices = np.array([
    0, 6,
    3, 9,
    2, 4, 8, 10,
    1, 5, 7, 11,
])
mtx = np.zeros((12, 12))
mtx[row_indices, col_indices] = 1.


class LocalStiffPermuted(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = ny = surface['mesh'].shape[1]

        self.add_input('local_stiff', shape=(ny - 1, 12, 12))
        self.add_output('local_stiff_permuted', shape=(ny - 1, 12, 12))

        indices = np.arange(144 * (ny - 1)).reshape((ny - 1, 12, 12))
        ones = np.ones((12, 12), int)

        # data = np.empty((ny - 1, 12, 12, 12, 12))
        # for row1 in range(12):
        #     for col1 in range(12):
        #         for row2 in range(12):
        #             for col2 in range(12):
        #                 data[:, row1, col1, row2, col2] = mtx[row1, row2] * mtx[col2, col1]
        data = np.einsum('i,jl,mk->ijklm', np.ones(ny - 1), mtx.T, mtx)

        data = data.flatten()
        rows = np.einsum('ijk,lm->ijklm', indices, ones).flatten()
        cols = np.einsum('ilm,jk->ijklm', indices, ones).flatten()
        self.declare_partials('local_stiff_permuted', 'local_stiff', val=data, rows=rows, cols=cols)

    def compute(self, inputs, outputs):
        outputs['local_stiff_permuted'] = np.einsum('jl,ilm,mk->ijk', mtx.T, inputs['local_stiff'], mtx)
