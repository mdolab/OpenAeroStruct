import numpy as np

import openmdao.api as om

from openaerostruct.utils.vector_algebra import compute_norm


class Length(om.ExplicitComponent):

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        surface = self.options['surface']

        self.ny = ny = surface['mesh'].shape[1]

        self.add_input('nodes', shape=(ny, 3), units='m')
        self.add_output('element_lengths', shape=ny - 1, units='m')

        mesh_indices = np.arange(3 * ny).reshape((ny, 3))
        length_indices = np.arange(ny - 1)

        rows = np.tile(np.outer(length_indices, np.ones(3, int)).flatten(), 2)
        cols = np.concatenate([
            mesh_indices[:-1, :].flatten(),
            mesh_indices[1: , :].flatten(),
        ])
        self.declare_partials('element_lengths', 'nodes', rows=rows, cols=cols)
        self.set_check_partial_options(wrt='nodes', method='fd', step=1e-6)

    def compute(self, inputs, outputs):
        vec = inputs['nodes'][1:, :] - inputs['nodes'][:-1, :]

        outputs['element_lengths'] = np.linalg.norm(vec, axis=-1)

    def compute_partials(self, inputs, partials):
        ny = self.ny

        vec = inputs['nodes'][1:, :] - inputs['nodes'][:-1, :]

        derivs = partials['element_lengths', 'nodes'].reshape((2, ny - 1, 3))
        derivs[0, :, :] = -vec / compute_norm(vec)
        derivs[1, :, :] =  vec / compute_norm(vec)
