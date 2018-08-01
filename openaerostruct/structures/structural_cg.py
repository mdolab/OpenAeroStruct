from __future__ import print_function, division
import numpy as np

from openmdao.api import ExplicitComponent

try:
    from openaerostruct.fortran import OAS_API
    fortran_flag = True
    data_type = float
except:
    fortran_flag = False
    data_type = complex

class StructuralCG(ExplicitComponent):
    """ Compute center-of-gravity location of the spar elements.

    parameters
    ----------
    A[ny-1] : numpy array
        Areas for each FEM element.
    nodes[ny, 3] : numpy array
        Flattened array with coordinates for each FEM node.

    Returns
    -------
    structural_weight : float
        Weight of the structural spar.
    cg_location[3] : numpy array
        Location of the structural spar's cg.

    """

    def initialize(self):
        self.options.declare('surface', types=dict)

    def setup(self):
        self.surface = surface = self.options['surface']

        self.ny = surface['num_y']

        self.add_input('nodes', val=np.zeros((self.ny, 3)), units='m')#, dtype=data_type))
        self.add_input('structural_weight', val=0., units='N')
        self.add_input('element_weights', val=np.zeros((self.ny-1)), units='N')
        self.add_output('cg_location', val=np.zeros((3)), units='m')#, dtype=data_type))

        self.declare_partials('*', '*')

    def compute(self, inputs, outputs):
        nodes = inputs['nodes']
        structural_weight = inputs['structural_weight']
        element_weights = inputs['element_weights']

        # Calculate the center-of-gravity location of the spar elements only
        center_of_elements = (nodes[1:, :] + nodes[:-1, :]) / 2.
        cg_loc = np.sum(center_of_elements.T * element_weights, axis=1) / structural_weight

        # If the tube is symmetric, double the computed weight and set the
        # y-location of the cg to 0, at the symmetry plane
        if self.surface['symmetry']:
            cg_loc[1] = 0.
            cg_loc *= 2.

        outputs['cg_location'] = cg_loc

    def compute_partials(self, inputs, partials):
        nodes = inputs['nodes']
        structural_weight = inputs['structural_weight']
        element_weights = inputs['element_weights']
