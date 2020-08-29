import numpy as np

import openmdao.api as om
from openaerostruct.structures.utils import radii


class SparWithinWing(om.ExplicitComponent):
    """
    Create a constraint to see if the spar is within the wing.
    This is based on the wing's t/c and the spar radius.

    .. warning::
        This component has not been extensively tested.
        It may require additional coding to work as intended.

    parameters
    ----------
    mesh[nx, ny, 3] : numpy array
        Array defining the nodal points of the lifting surface.
    radius[ny-1] : numpy array
        Radius of each element of the FEM spar.
    t_over_c[ny-1] : numpy array
        The streamwise thickness-to-chord ratio of each VLM panel.

    Returns
    -------
    spar_within_wing[ny-1] : numpy array
        If all the values are negative, each element is within the wing,
        based on the surface's t_over_c value and the current chord.
        Set a constraint with
        `OASProblem.add_constraint('spar_within_wing', upper=0.)`
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['mesh'].shape[1]
        nx = surface['mesh'].shape[0]

        self.add_input('mesh', val=np.zeros((nx, self.ny, 3)), units='m')
        self.add_input('radius', val=np.zeros((self.ny-1)), units='m')
        self.add_input('t_over_c', val=np.zeros((self.ny-1)))
        self.add_output('spar_within_wing', val=np.zeros((self.ny-1)), units='m')

        self.declare_partials('spar_within_wing', 'mesh', method='cs')

        arange = np.arange(self.ny - 1)
        self.declare_partials('spar_within_wing', 'radius', rows=arange, cols=arange, val=1.)

    def compute(self, inputs, outputs):
        mesh = inputs['mesh']
        t_over_c = inputs['t_over_c']
        max_radius = radii(mesh, t_over_c)
        outputs['spar_within_wing'] = inputs['radius'] - max_radius
