import numpy as np

import openmdao.api as om


class NonIntersectingThickness(om.ExplicitComponent):
    """
    Create a constraint so the thickness of the spar does not intersect
    itself in the center of the spar. Basically, the thickness must be less
    than or equal to the radius.

    parameters
    ----------
    thickness[ny-1] : numpy array
        Thickness of each element of the FEM spar.
    radius[ny-1] : numpy array
        Radius of each element of the FEM spar.

    Returns
    -------
    thickness_intersects[ny-1] : numpy array
        If all the values are negative, each element does not intersect itself.
        If a value is positive, then the thickness within the hollow spar
        intersects itself and presents an impossible design.
        Add a constraint as
        `OASProblem.add_constraint('thickness_intersects', upper=0.)`
    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['mesh'].shape[1]

        self.add_input('thickness', val=np.zeros((self.ny-1)), units='m')
        self.add_input('radius', val=np.zeros((self.ny-1)), units='m')
        self.add_output('thickness_intersects', val=np.zeros((self.ny-1)), units='m')

        arange = np.arange(self.ny-1)
        self.declare_partials('thickness_intersects', 'thickness', rows=arange, cols=arange, val=1.)
        self.declare_partials('thickness_intersects', 'radius', rows=arange, cols=arange, val=-1.)

    def compute(self, inputs, outputs):
        outputs['thickness_intersects'] = inputs['thickness'] - inputs['radius']
