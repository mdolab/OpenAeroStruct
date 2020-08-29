import numpy as np

import openmdao.api as om


class ComputeNodes(om.ExplicitComponent):
    """
    Compute FEM nodes based on aerodynamic mesh.

    The FEM nodes are placed at fem_origin * chord,
    with the default fem_origin = 0.35.

    Parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.

    Returns
    -------
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']
        mesh = surface['mesh']
        nx = mesh.shape[0]
        ny = mesh.shape[1]

        if surface['fem_model_type'] == 'tube':
            self.fem_origin = surface['fem_origin']
        else:
            y_upper = surface['data_y_upper']
            x_upper = surface['data_x_upper']
            y_lower = surface['data_y_lower']

            self.fem_origin = (x_upper[0]  * (y_upper[0]  - y_lower[0]) +
                               x_upper[-1] * (y_upper[-1] - y_lower[-1])) / \
                             ((y_upper[0]  -  y_lower[0]) + (y_upper[-1] - y_lower[-1]))

        self.add_input('mesh', val=np.zeros((nx, ny, 3)), units='m')
        self.add_output('nodes', val=np.zeros((ny, 3)), units='m')

        w = self.fem_origin
        n = ny * 3

        data = np.zeros((2 * n))
        data[:n] = 1 - w
        data[n:] = w

        arange = np.arange(n)

        rows = np.hstack((arange, arange))
        cols = np.hstack((arange, arange + (nx - 1) * n))

        self.declare_partials('nodes', 'mesh', rows=rows, cols=cols, val=data)

    def compute(self, inputs, outputs):
        w = self.fem_origin
        mesh = inputs['mesh']
        outputs['nodes'] = (1-w) * mesh[0, :, :] + w * mesh[-1, :, :]
