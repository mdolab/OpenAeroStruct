import numpy as np

import openmdao.api as om


def norm(vec):
    return np.sqrt(np.sum(vec**2))

class WingboxFuelVol(om.ExplicitComponent):
    """
    Computes the internal volumes of the wingbox segments.

    parameters
    ----------
    nodes[ny, 3] : numpy array
        Coordinates of FEM nodes.
    A_int[ny-1] : numpy array
        Internal cross-sectional area of each wingbox segment.

    Returns
    -------
    fuel_vols[ny-1] : numpy array
        Internal volume of each wingbox segment.
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['mesh'].shape[1]

        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')
        self.add_input('A_int', val=np.zeros((self.ny-1)), units='m**2')
        self.add_output('fuel_vols', val=np.zeros((self.ny-1)), units='m**3')

        self.declare_partials('*', '*', method='cs')

    def compute(self, inputs, outputs):
        nodes = inputs['nodes']

        element_lengths = np.zeros(self.ny-1, dtype=type(nodes[0, 0]))

        for i in range(self.ny - 1):
            element_lengths[i] = norm(nodes[i+1] - nodes[i])

        # Next we multiply the element lengths with the A_int for the internal volumes of the wingobox segments
        vols = element_lengths * inputs['A_int']

        outputs['fuel_vols'] = vols
